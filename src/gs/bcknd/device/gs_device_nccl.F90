! Copyright (c) 2025, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Defines NCCL based gather-scatter communication
module gs_device_nccl
  use num_types, only : rp, c_rp
  use gs_comm, only : gs_comm_t
  use gs_ops
  use stack, only : stack_i4_t
  use comm, only : pe_size, pe_rank
  use htable, only : htable_i4_t
  use device
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding, only : c_sizeof, c_int32_t, &
       c_ptr, C_NULL_PTR, c_size_t, c_associated, c_int
  implicit none
  private

  !> Buffers for non-blocking communication and packing/unpacking
  type, private :: gs_device_nccl_buf_t
     integer, allocatable :: ndofs(:) !< Number of dofs
     integer, allocatable :: offset(:) !< Offset into buf
     integer :: total !< Total number of dofs
     type(c_ptr) :: buf_d = C_NULL_PTR !< Device buffer
     type(c_ptr) :: dof_d = C_NULL_PTR !< Dof mapping for pack/unpack
   contains
     procedure, pass(this) :: init => gs_device_nccl_buf_init
     procedure, pass(this) :: free => gs_device_nccl_buf_free
  end type gs_device_nccl_buf_t

  !> Gather-scatter communication using NCCL
  !! The arrays are indexed per PE like @a send_pe and @ recv_pe.
  type, public, extends(gs_comm_t) :: gs_device_nccl_t
     type(gs_device_nccl_buf_t) :: send_buf
     type(gs_device_nccl_buf_t) :: recv_buf
     type(c_ptr), allocatable :: stream(:)
     type(c_ptr), allocatable :: event(:)
     integer :: nb_strtgy
     type(c_ptr) :: send_event = C_NULL_PTR
   contains
     procedure, pass(this) :: init => gs_device_nccl_init
     procedure, pass(this) :: free => gs_device_nccl_free
     procedure, pass(this) :: nbsend => gs_device_nccl_nbsend
     procedure, pass(this) :: nbrecv => gs_device_nccl_nbrecv
     procedure, pass(this) :: nbwait => gs_device_nccl_nbwait
  end type gs_device_nccl_t

#ifdef HAVE_HIP
  interface
     subroutine hip_gs_pack(u_d, buf_d, dof_d, offset, n, stream) &
          bind(c, name='hip_gs_pack')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: n, offset
       type(c_ptr), value :: u_d, buf_d, dof_d, stream
     end subroutine hip_gs_pack
  end interface

  interface
     subroutine hip_gs_unpack(u_d, op, buf_d, dof_d, offset, n, stream) &
          bind(c, name='hip_gs_unpack')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: op, offset, n
       type(c_ptr), value :: u_d, buf_d, dof_d, stream
     end subroutine hip_gs_unpack
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_gs_pack(u_d, buf_d, dof_d, offset, n, stream) &
          bind(c, name='cuda_gs_pack')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: n, offset
       type(c_ptr), value :: u_d, buf_d, dof_d, stream
     end subroutine cuda_gs_pack
  end interface

  interface
     subroutine cuda_gs_unpack(u_d, op, buf_d, dof_d, offset, n, stream) &
          bind(c, name='cuda_gs_unpack')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: op, offset, n
       type(c_ptr), value :: u_d, buf_d, dof_d, stream
     end subroutine cuda_gs_unpack
  end interface
#endif

  interface
     subroutine device_nccl_sendrecv(sbuf_d, soffset, scount, srank, &
          rbuf_d, roffset, rcount, rrank, nbytes, stream) &
          bind(c, name='device_nccl_sendrecv')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: soffset, scount, roffset, rcount
       integer(c_int), value :: srank, rrank, nbytes
       type(c_ptr), value :: sbuf_d, rbuf_d, stream
     end subroutine device_nccl_sendrecv
  end interface

contains

  subroutine gs_device_nccl_buf_init(this, pe_order, dof_stack, mark_dupes)
    class(gs_device_nccl_buf_t), intent(inout) :: this
    integer, allocatable, intent(inout) :: pe_order(:)
    type(stack_i4_t), allocatable, intent(inout) :: dof_stack(:)
    logical, intent(in) :: mark_dupes
    integer, allocatable :: dofs(:)
    integer :: i, j, total
    integer(c_size_t) :: sz
    type(htable_i4_t) :: doftable
    integer :: dupe, marked, k
    real(c_rp) :: rp_dummy
    integer(c_int32_t) :: i4_dummy


    allocate(this%ndofs(size(pe_order)))
    allocate(this%offset(size(pe_order)))

    total = 0
    do i = 1, size(pe_order)
       this%ndofs(i) = dof_stack(pe_order(i))%size()
       this%offset(i) = total
       total = total + this%ndofs(i)
    end do

    this%total = total

    sz = c_sizeof(rp_dummy) * total
    call device_alloc(this%buf_d, sz)

    sz = c_sizeof(i4_dummy) * total
    call device_alloc(this%dof_d, sz)

    if (mark_dupes) call doftable%init(2*total)
    allocate(dofs(total))

    ! Copy from dof_stack into dofs, optionally marking duplicates with doftable
    marked = 0
    do i = 1, size(pe_order)
       ! %array() breaks on cray
       select type (arr => dof_stack(pe_order(i))%data)
       type is (integer)
          do j = 1, this%ndofs(i)
             k = this%offset(i) + j
             if (mark_dupes) then
                if (doftable%get(arr(j), dupe) .eq. 0) then
                   if (dofs(dupe) .gt. 0) then
                      dofs(dupe) = -dofs(dupe)
                      marked = marked + 1
                   end if
                   dofs(k) = -arr(j)
                   marked = marked + 1
                else
                   call doftable%set(arr(j), k)
                   dofs(k) = arr(j)
                end if
             else
                dofs(k) = arr(j)
             end if
          end do
       end select
    end do

    call device_memcpy(dofs, this%dof_d, total, HOST_TO_DEVICE, sync=.false.)

    deallocate(dofs)
    call doftable%free()

  end subroutine gs_device_nccl_buf_init

  subroutine gs_device_nccl_buf_free(this)
    class(gs_device_nccl_buf_t), intent(inout) :: this

    if (allocated(this%ndofs)) deallocate(this%ndofs)
    if (allocated(this%offset)) deallocate(this%offset)

    if (c_associated(this%buf_d)) call device_free(this%buf_d)
    if (c_associated(this%dof_d)) call device_free(this%dof_d)
  end subroutine gs_device_nccl_buf_free

  !> Initialise NCCL based communication method
  subroutine gs_device_nccl_init(this, send_pe, recv_pe)
    class(gs_device_nccl_t), intent(inout) :: this
    type(stack_i4_t), intent(inout) :: send_pe
    type(stack_i4_t), intent(inout) :: recv_pe
    integer :: i

    call this%init_order(send_pe, recv_pe)

    call this%send_buf%init(this%send_pe, this%send_dof, .false.)
    call this%recv_buf%init(this%recv_pe, this%recv_dof, .true.)

#if defined(HAVE_HIP) || defined(HAVE_CUDA)
    ! Create a set of non-blocking streams
    allocate(this%stream(size(this%recv_pe)))
    do i = 1, size(this%recv_pe)
       call device_stream_create_with_priority(this%stream(i), 1, STRM_HIGH_PRIO)
    end do

    allocate(this%event(size(this%recv_pe)))
    do i = 1, size(this%recv_pe)
       call device_event_create(this%event(i), 2)
    end do
#endif

  end subroutine gs_device_nccl_init

  !> Deallocate MPI based communication method
  subroutine gs_device_nccl_free(this)
    class(gs_device_nccl_t), intent(inout) :: this
    integer :: i

    call this%send_buf%free()
    call this%recv_buf%free()

    call this%free_order()
    call this%free_dofs()

#if defined(HAVE_HIP) || defined(HAVE_CUDA)
    if (allocated(this%stream)) then
       do i = 1, size(this%stream)
          call device_stream_destroy(this%stream(i))
       end do
       deallocate(this%stream)
    end if
#endif

  end subroutine gs_device_nccl_free

  !> Post non-blocking send operations
  subroutine gs_device_nccl_nbsend(this, u, n, deps, strm)
    class(gs_device_nccl_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    type(c_ptr), intent(inout) :: deps
    type(c_ptr), intent(inout) :: strm
    integer :: i
    type(c_ptr) :: u_d

    u_d = device_get_ptr(u)

    do i = 1, size(this%send_pe)
       call device_stream_wait_event(this%stream(i), deps, 0)
#ifdef HAVE_HIP
       call hip_gs_pack(u_d, &
            this%send_buf%buf_d, &
            this%send_buf%dof_d, &
            this%send_buf%offset(i), &
            this%send_buf%ndofs(i), &
            this%stream(i))
#elif HAVE_CUDA
       call cuda_gs_pack(u_d, &
            this%send_buf%buf_d, &
            this%send_buf%dof_d, &
            this%send_buf%offset(i), &
            this%send_buf%ndofs(i), &
            this%stream(i))
#else
       call neko_error('gs_device_nccl: no backend')
#endif
    end do


    ! Everything else is done in the wait routine

  end subroutine gs_device_nccl_nbsend

  !> Post non-blocking receive operations
  subroutine gs_device_nccl_nbrecv(this)
    class(gs_device_nccl_t), intent(inout) :: this
    integer :: i

    ! Everything is done in the wait routine

  end subroutine gs_device_nccl_nbrecv

  !> Wait for non-blocking operations
  subroutine gs_device_nccl_nbwait(this, u, n, op, strm)
    class(gs_device_nccl_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: u
    type(c_ptr), intent(inout) :: strm
    integer :: op, done_req, i
    type(c_ptr) :: u_d
    real(c_rp) :: rp_dummy
    integer(c_int) :: nbytes

    u_d = device_get_ptr(u)
    nbytes = c_sizeof(rp_dummy)

    do i = 1, size(this%send_pe)

       call device_nccl_sendrecv(this%send_buf%buf_d, &
            nbytes*this%send_buf%offset(i), &
            this%send_buf%ndofs(i), &
            this%send_pe(i), &
            this%recv_buf%buf_d, &
            nbytes*this%recv_buf%offset(i), &
            this%recv_buf%ndofs(i), &
            this%recv_pe(i), &
            nbytes, &
            this%stream(i))

#ifdef HAVE_HIP
       call hip_gs_unpack(u_d, op, &
            this%recv_buf%buf_d, &
            this%recv_buf%dof_d, &
            this%recv_buf%offset(i), &
            this%recv_buf%ndofs(i), &
            this%stream(i))
#elif HAVE_CUDA
       call cuda_gs_unpack(u_d, op, &
            this%recv_buf%buf_d, &
            this%recv_buf%dof_d, &
            this%recv_buf%offset(i), &
            this%recv_buf%ndofs(i), &
            this%stream(i))
#else
       call neko_error('gs_device_mpi: no backend')
#endif
       call device_event_record(this%event(i), this%stream(i))
    end do

    ! Sync non-blocking streams
    do done_req = 1, size(this%recv_pe)
       call device_stream_wait_event(strm, &
            this%event(done_req), 0)
    end do

  end subroutine gs_device_nccl_nbwait

end module gs_device_nccl
