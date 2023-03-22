! Copyright (c) 2019-2021, The Neko Authors
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
module mesh_conn
  use num_types
  implicit none
  private

  public :: mesh_conn_t

  ! Connectivity information regarding vertices, faces and edges.
  ! It contains both object global numbering and communication information
  type obj_conn_t
     integer(i4) :: lnum, lown ! number of local and owned objects
     integer(i8) :: goff ! global object offset
     integer(i8) :: gnum ! global number of objects
     integer(i4), allocatable, dimension(:,:) :: lmap ! element
                                       ! vertices/faces/edges to object mapping
     integer(i4) :: nrank, nshare ! number of MPI ranks sharing objects and
                                  ! number of shared objects
     integer(i8), allocatable, dimension(:) :: lgidx ! global indexing of
                                                  ! unique objects of given type
     integer(i4), allocatable, dimension(:) :: lrank ! list of ranks sharing
                                                     ! objects
     integer(i4), allocatable, dimension(:) :: lshare ! list of shared objects
     integer(i4), allocatable, dimension(:) :: loff ! offset in the lshare list
   contains
     procedure, public, pass(this) :: free => obj_conn_free
  end type obj_conn_t

  !> Base type for mesh connectivity information
  type :: mesh_conn_t
     logical, private :: ifnonconf_ !< Nonconforming mesh flag
     integer(i4) :: nel ! local number of elements
     integer(i4), allocatable, dimension(:,:) :: falg ! face alignment
     integer(i4), allocatable, dimension(:,:) :: ealg ! edge alignment
     type(obj_conn_t) :: vert ! vertex info
     type(obj_conn_t) :: face ! face info
     type(obj_conn_t) :: edge ! edges info; 3D mesh only
   contains
     procedure, pass(this) :: nonconf => mesh_conn_nonconf
     procedure, pass(this) :: set_nonconf => mesh_conn_set_nonconf
     procedure, pass(this) :: free => mesh_conn_free
  end type mesh_conn_t

contains

  ! Type bounded routines
  subroutine obj_conn_free(this)
    ! argument list
    class(obj_conn_t), intent(inout) :: this

    ! Deallocate arrays
    if (allocated(this%lmap)) deallocate(this%lmap)
    if (allocated(this%lgidx)) deallocate(this%lgidx)
    if (allocated(this%lrank)) deallocate(this%lrank)
    if (allocated(this%lshare)) deallocate(this%lshare)
    if (allocated(this%loff)) deallocate(this%loff)

    ! Reset registers
    this%lnum = 0
    this%lown = 0
    this%goff = 0
    this%gnum = 0
    this%nrank = 0
    this%nshare = 0

    return
  end subroutine obj_conn_free

  !> Return the mesh nonconformity flag
  pure function mesh_conn_nonconf(this) result(ifnonconf)
    class(mesh_conn_t), intent(in) :: this
    logical :: ifnonconf
    ifnonconf = this%ifnonconf_
  end function mesh_conn_nonconf

  !> Set the mesh nonconformity flag
  subroutine mesh_conn_set_nonconf(this, ifnonconf)
    class(mesh_conn_t), intent(inout) :: this
    logical, intent(in) :: ifnonconf
    this%ifnonconf_ = ifnonconf
  end subroutine mesh_conn_set_nonconf

  !> Free memory
  subroutine mesh_conn_free(this)
    ! argument list
    class(mesh_conn_t), intent(inout) :: this

    ! Deallocate arrays
    if (allocated(this%falg)) deallocate(this%falg)
    if (allocated(this%ealg)) deallocate(this%ealg)

    ! Free types
    call this%vert%free()
    call this%face%free()
    call this%edge%free()

    ! Reset registers
    this%nel = 0

    ! Reset mesh nonconformity flag
    call this%set_nonconf(.false.)

    return
  end subroutine mesh_conn_free

end module mesh_conn
