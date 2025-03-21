! Copyright (c) 2008-2020, UCHICAGO ARGONNE, LLC.
!
! The UChicago Argonne, LLC as Operator of Argonne National
! Laboratory holds copyright in the Software. The copyright holder
! reserves all rights except those expressly granted to licensees,
! and U.S. Government license rights.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
! 1. Redistributions of source code must retain the above copyright
! notice, this list of conditions and the disclaimer below.
!
! 2. Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the disclaimer (as noted below)
! in the documentation and/or other materials provided with the
! distribution.
!
! 3. Neither the name of ANL nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
! UCHICAGO ARGONNE, LLC, THE U.S. DEPARTMENT OF
! ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
! TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! Additional BSD Notice
! ---------------------
! 1. This notice is required to be provided under our contract with
! the U.S. Department of Energy (DOE). This work was produced at
! Argonne National Laboratory under Contract
! No. DE-AC02-06CH11357 with the DOE.
!
! 2. Neither the United States Government nor UCHICAGO ARGONNE,
! LLC nor any of their employees, makes any warranty,
! express or implied, or assumes any liability or responsibility for the
! accuracy, completeness, or usefulness of any information, apparatus,
! product, or process disclosed, or represents that its use would not
! infringe privately-owned rights.
!
! 3. Also, reference herein to any specific commercial products, process,
! or services by trade name, trademark, manufacturer or otherwise does
! not necessarily constitute or imply its endorsement, recommendation,
! or favoring by the United States Government or UCHICAGO ARGONNE LLC.
! The views and opinions of authors expressed
! herein do not necessarily state or reflect those of the United States
! Government or UCHICAGO ARGONNE, LLC, and shall
! not be used for advertising or product endorsement purposes.
!
!> Project x onto X, the space of old solutions and back again
!! @note In this code we assume that the matrix project for the
!! pressure Ax does not vary in time.
module projection
  use num_types, only : rp, c_rp
  use math, only : rzero, glsc3, add2, copy, cmult
  use coefs, only : coef_t
  use ax_product, only : ax_t
  use bc_list, only : bc_list_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use neko_config, only : NEKO_BCKND_DEVICE, NEKO_BLK_SIZE, &
       NEKO_DEVICE_MPI, NEKO_BCKND_OPENCL
  use device, only : device_alloc, HOST_TO_DEVICE, device_memcpy, &
       device_get_ptr, device_free, device_map
  use device_math, only : device_glsc3, device_add2s2, device_cmult, &
       device_rzero, device_copy, device_add2, device_add2s2_many, &
       device_glsc3_many
  use device_projection, only : device_proj_on, device_project_ortho
  use profiler, only : profiler_start_region, profiler_end_region
  use logger, only : LOG_SIZE, neko_log
  use utils, only : neko_warning
  use bc_list, only : bc_list_t
  use time_step_controller, only : time_step_controller_t
  use comm, only : NEKO_COMM, pe_rank, MPI_Allreduce, MPI_IN_PLACE, &
       MPI_SUM, MPI_REAL_PRECISION
  use, intrinsic :: iso_c_binding, only : c_ptr, c_size_t, &
       c_sizeof, C_NULL_PTR, c_loc, c_associated
  implicit none
  private
  public :: proj_ortho

  type, public :: projection_t
     real(kind=rp), allocatable :: xx(:,:)
     real(kind=rp), allocatable :: bb(:,:)
     real(kind=rp), allocatable :: xbar(:)
     type(c_ptr), allocatable :: xx_d(:)
     type(c_ptr), allocatable :: bb_d(:)
     type(c_ptr) :: xbar_d = C_NULL_PTR
     type(c_ptr) :: alpha_d = C_NULL_PTR
     type(c_ptr) :: xx_d_d = C_NULL_PTR
     type(c_ptr) :: bb_d_d = C_NULL_PTR
     integer :: m, L
     real(kind=rp) :: tol = 1e-7_rp
     !logging variables
     real(kind=rp) :: proj_res
     integer :: proj_m = 0
     integer :: activ_step ! steps to activate projection
   contains
     procedure, pass(this) :: clear => bcknd_clear
     procedure, pass(this) :: project_on => bcknd_project_on
     procedure, pass(this) :: project_back => bcknd_project_back
     procedure, pass(this) :: log_info => print_proj_info
     procedure, pass(this) :: init => projection_init
     procedure, pass(this) :: free => projection_free
     procedure, pass(this) :: pre_solving => projection_pre_solving
     procedure, pass(this) :: post_solving => projection_post_solving
  end type projection_t

contains

  subroutine projection_init(this, n, L, activ_step)
    class(projection_t), target, intent(inout) :: this
    integer, intent(in) :: n
    integer, optional, intent(in) :: L, activ_step
    integer :: i
    integer(c_size_t) :: ptr_size
    type(c_ptr) :: ptr
    real(c_rp) :: dummy

    call this%free()

    if (present(L)) then
       this%L = L
    else
       this%L = 20
    end if

    if (present(activ_step)) then
       this%activ_step = activ_step
    else
       this%activ_step = 5
    end if

    this%m = 0

    allocate(this%xx(n, this%L))
    allocate(this%bb(n, this%L))
    allocate(this%xbar(n))
    allocate(this%xx_d(this%L))
    allocate(this%bb_d(this%L))
    call rzero(this%xbar, n)
    do i = 1, this%L
       call rzero(this%xx(1, i), n)
       call rzero(this%bb(1, i), n)
    end do
    if (NEKO_BCKND_DEVICE .eq. 1) then

       call device_map(this%xbar, this%xbar_d, n)
       call device_alloc(this%alpha_d, int(c_sizeof(dummy)*this%L, c_size_t))

       call device_rzero(this%xbar_d, n)
       call device_rzero(this%alpha_d, this%L)

       do i = 1, this%L
          this%xx_d(i) = C_NULL_PTR
          call device_map(this%xx(:, i), this%xx_d(i), n)
          call device_rzero(this%xx_d(i), n)
          this%bb_d(i) = C_NULL_PTR
          call device_map(this%bb(:, i), this%bb_d(i), n)
          call device_rzero(this%bb_d(i), n)
       end do

       ptr_size = c_sizeof(C_NULL_PTR) * this%L
       call device_alloc(this%xx_d_d, ptr_size)
       ptr = c_loc(this%xx_d)
       call device_memcpy(ptr, this%xx_d_d, ptr_size, &
            HOST_TO_DEVICE, sync = .false.)
       call device_alloc(this%bb_d_d, ptr_size)
       ptr = c_loc(this%bb_d)
       call device_memcpy(ptr, this%bb_d_d, ptr_size, &
            HOST_TO_DEVICE, sync = .false.)
    end if


  end subroutine projection_init

  subroutine projection_free(this)
    class(projection_t), intent(inout) :: this
    integer :: i
    if (allocated(this%xx)) then
       deallocate(this%xx)
    end if
    if (allocated(this%bb)) then
       deallocate(this%bb)
    end if
    if (allocated(this%xbar)) then
       deallocate(this%xbar)
    end if
    if (allocated(this%xx_d)) then
       do i = 1, this%L
          if (c_associated(this%xx_d(i))) then
             call device_free(this%xx_d(i))
          end if
       end do
    end if
    if (c_associated(this%xx_d_d)) then
       call device_free(this%xx_d_d)
    end if
    if (c_associated(this%xbar_d)) then
       call device_free(this%xbar_d)
    end if
    if (c_associated(this%alpha_d)) then
       call device_free(this%alpha_d)
    end if
    if (allocated(this%bb_d)) then
       do i = 1, this%L
          if (c_associated(this%bb_d(i))) then
             call device_free(this%bb_d(i))
          end if
       end do
    end if
    if (c_associated(this%bb_d_d)) then
       call device_free(this%bb_d_d)
    end if

  end subroutine projection_free

  subroutine projection_pre_solving(this, b, tstep, coef, n, dt_controller, &
       string)
    class(projection_t), intent(inout) :: this
    integer, intent(inout) :: n
    real(kind=rp), intent(inout), dimension(n) :: b
    integer, intent(in) :: tstep
    class(coef_t), intent(inout) :: coef
    type(time_step_controller_t), intent(in) :: dt_controller
    character(len=*), optional :: string

    if (tstep .gt. this%activ_step .and. this%L .gt. 0) then
       if (dt_controller%if_variable_dt) then
          ! the time step at which dt is changed
          if (dt_controller%dt_last_change .eq. 0) then
             call this%clear(n)
          else if (dt_controller%dt_last_change .gt. this%activ_step - 1) then
             ! activate projection some steps after dt is changed
             ! note that dt_last_change start from 0
             call this%project_on(b, coef, n)
             if (present(string)) then
                call this%log_info(string)
             end if
          end if
       else
          call this%project_on(b, coef, n)
          if (present(string)) then
             call this%log_info(string)
          end if
       end if
    end if

  end subroutine projection_pre_solving

  subroutine projection_post_solving(this, x, Ax, coef, bclst, gs_h, n, tstep, &
       dt_controller)
    class(projection_t), intent(inout) :: this
    integer, intent(inout) :: n
    class(Ax_t), intent(inout) :: Ax
    class(coef_t), intent(inout) :: coef
    class(bc_list_t), intent(inout) :: bclst
    type(gs_t), intent(inout) :: gs_h
    real(kind=rp), intent(inout), dimension(n) :: x
    integer, intent(in) :: tstep
    type(time_step_controller_t), intent(in) :: dt_controller

    if (tstep .gt. this%activ_step .and. this%L .gt. 0) then
       if (.not.(dt_controller%if_variable_dt) .or. &
            (dt_controller%dt_last_change .gt. this%activ_step - 1)) then
          call this%project_back(x, Ax, coef, bclst, gs_h, n)
       end if
    end if

  end subroutine projection_post_solving

  subroutine bcknd_project_on(this, b, coef, n)
    class(projection_t), intent(inout) :: this
    integer, intent(inout) :: n
    class(coef_t), intent(inout) :: coef
    real(kind=rp), intent(inout), dimension(n) :: b
    call profiler_start_region('Project on', 16)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_project_on(this, b, coef, n)
    else
       call cpu_project_on(this, b, coef, n)
    end if
    call profiler_end_region('Project on', 16)
  end subroutine bcknd_project_on

  subroutine bcknd_project_back(this, x, Ax, coef, bclst, gs_h, n)
    class(projection_t), intent(inout) :: this
    integer, intent(inout) :: n
    class(Ax_t), intent(inout) :: Ax
    class(coef_t), intent(inout) :: coef
    class(bc_list_t), intent(inout) :: bclst
    type(gs_t), intent(inout) :: gs_h
    real(kind=rp), intent(inout), dimension(n) :: x
    type(c_ptr) :: x_d

    call profiler_start_region('Project back', 17)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       x_d = device_get_ptr(x)
       ! Restore desired solution
       if (this%m .gt. 0) call device_add2(x_d, this%xbar_d, n)
       if (this%m .eq. this%L) then
          this%m = 1
       else
          this%m = min(this%m+1, this%L)
       end if

       call device_copy(this%xx_d(this%m), x_d, n) ! Update (X,B)

    else
       if (this%m .gt. 0) call add2(x, this%xbar, n) ! Restore desired solution
       if (this%m .eq. this%L) then
          this%m = 1
       else
          this%m = min(this%m+1, this%L)
       end if

       call copy(this%xx(1, this%m), x, n) ! Update (X,B)
    end if

    call Ax%compute(this%bb(1, this%m), x, coef, coef%msh, coef%Xh)
    call gs_h%gs_op_vector(this%bb(1, this%m), n, GS_OP_ADD)
    call bclst%apply_scalar(this%bb(1, this%m), n)

    call proj_ortho(this, coef, n)
    call profiler_end_region('Project back', 17)
  end subroutine bcknd_project_back



  subroutine cpu_project_on(this, b, coef, n)
    class(projection_t), intent(inout) :: this
    integer, intent(inout) :: n
    class(coef_t), intent(inout) :: coef
    real(kind=rp), intent(inout), dimension(n) :: b
    integer :: i, j, k, l, ierr
    real(kind=rp) :: work(this%L), alpha(this%L), s

    associate(xbar => this%xbar, xx => this%xx, &
         bb => this%bb)

      if (this%m .le. 0) return

      !First round of CGS
      call rzero(alpha, this%m)
      this%proj_res = sqrt(glsc3(b, b, coef%mult, n) / coef%volume)
      this%proj_m = this%m
      do i = 1, n, NEKO_BLK_SIZE
         j = min(NEKO_BLK_SIZE, n-i+1)
         do k = 1, this%m
            s = 0.0_rp
            do l = 0, (j-1)
               s = s + xx(i+l, k) * coef%mult(i+l,1,1,1) * b(i+l)
            end do
            alpha(k) = alpha(k) + s
         end do
      end do

      !First one outside loop to avoid zeroing xbar and bbar
      call MPI_Allreduce(MPI_IN_PLACE, alpha, this%m, &
           MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)

      call rzero(work, this%m)

      do i = 1, n, NEKO_BLK_SIZE
         j = min(NEKO_BLK_SIZE, n-i+1)
         do l = 0, (j-1)
            xbar(i+l) = alpha(1) * xx(i+l,1)
            b(i+l) = b(i+l) - alpha(1) * bb(i+l,1)
         end do
         do k = 2, this%m
            do l = 0, (j-1)
               xbar(i+l) = xbar(i+l) + alpha(k) * xx(i+l,k)
               b(i+l) = b(i+l)- alpha(k) * bb(i+l,k)
            end do
         end do
         !Second round of CGS
         do k = 1, this%m
            s = 0.0_rp
            do l = 0, (j-1)
               s = s + xx(i+l,k) * coef%mult(i+l,1,1,1) * b(i+l)
            end do
            work(k) = work(k) + s
         end do
      end do

      call MPI_Allreduce(work, alpha, this%m, &
           MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)

      do i = 1, n, NEKO_BLK_SIZE
         j = min(NEKO_BLK_SIZE, n-i+1)
         do k = 1, this%m
            do l = 0, (j-1)
               xbar(i+l) = xbar(i+l) + alpha(k) * xx(i+l,k)
               b(i+l) = b(i+l) - alpha(k) * bb(i+l,k)
            end do
         end do
      end do
    end associate
  end subroutine cpu_project_on

  subroutine device_project_on(this, b, coef, n)
    class(projection_t), intent(inout) :: this
    integer, intent(inout) :: n
    class(coef_t), intent(inout) :: coef
    real(kind=rp), intent(inout), dimension(n) :: b
    real(kind=rp) :: alpha(this%L)
    type(c_ptr) :: b_d
    integer :: i
    b_d = device_get_ptr(b)

    associate(xbar_d => this%xbar_d, xx_d => this%xx_d, xx_d_d => this%xx_d_d, &
         bb_d => this%bb_d, bb_d_d => this%bb_d_d, alpha_d => this%alpha_d)

      if (this%m .le. 0) return



      this%proj_res = sqrt(device_glsc3(b_d, b_d, coef%mult_d, n)/coef%volume)
      this%proj_m = this%m
      if (NEKO_DEVICE_MPI .and. (NEKO_BCKND_OPENCL .ne. 1)) then
         call device_proj_on(alpha_d, b_d, xx_d_d, bb_d_d, &
              coef%mult_d, xbar_d, this%m, n)
      else
         if (NEKO_BCKND_OPENCL .eq. 1) then
            do i = 1, this%m
               alpha(i) = device_glsc3(b_d, xx_d(i), coef%mult_d, n)
            end do
         else
            call device_glsc3_many(alpha, b_d, xx_d_d, coef%mult_d, this%m, n)
            call device_memcpy(alpha, alpha_d, this%m, &
                 HOST_TO_DEVICE, sync = .false.)
         end if
         call device_rzero(xbar_d, n)
         if (NEKO_BCKND_OPENCL .eq. 1) then
            do i = 1, this%m
               call device_add2s2(xbar_d, xx_d(i), alpha(i), n)
            end do
            call cmult(alpha, -1.0_rp, this%m)
         else
            call device_add2s2_many(xbar_d, xx_d_d, alpha_d, this%m, n)
            call device_cmult(alpha_d, -1.0_rp, this%m)
         end if

         if (NEKO_BCKND_OPENCL .eq. 1) then
            do i = 1, this%m
               call device_add2s2(b_d, bb_d(i), alpha(i), n)
               alpha(i) = device_glsc3(b_d, xx_d(i), coef%mult_d, n)
            end do
         else
            call device_add2s2_many(b_d, bb_d_d, alpha_d, this%m, n)
            call device_glsc3_many(alpha, b_d, xx_d_d, coef%mult_d, this%m, n)
            call device_memcpy(alpha, alpha_d, this%m, &
                 HOST_TO_DEVICE, sync = .false.)
         end if

         if (NEKO_BCKND_OPENCL .eq. 1) then
            do i = 1, this%m
               call device_add2s2(xbar_d, xx_d(i), alpha(i), n)
               call cmult(alpha, -1.0_rp, this%m)
               call device_add2s2(b_d, bb_d(i), alpha(i), n)
            end do
         else
            call device_add2s2_many(xbar_d, xx_d_d, alpha_d, this%m, n)
            call device_cmult(alpha_d, -1.0_rp, this%m)
            call device_add2s2_many(b_d, bb_d_d, alpha_d, this%m, n)
         end if
      end if

    end associate
  end subroutine device_project_on

  !Choose between CPU or device for proj_ortho
  subroutine proj_ortho(this, coef, n)
    class(projection_t), intent(inout) :: this
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: n

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_proj_ortho(this, this%xx_d, this%bb_d, coef%mult_d, n)
    else
       call cpu_proj_ortho (this, this%xx, this%bb, coef%mult, n)
    end if
  end subroutine proj_ortho

  !This is a lot more primitive than on the CPU
  subroutine device_proj_ortho(this, xx_d, bb_d, w_d, n)
    type(projection_t), intent(inout) :: this
    integer, intent(in) :: n
    type(c_ptr), dimension(this%L) :: xx_d, bb_d
    type(c_ptr), intent(in) :: w_d
    real(kind=rp) :: nrm, scl
    real(kind=rp) :: alpha(this%L)
    integer :: i

    associate(m => this%m, xx_d_d => this%xx_d_d, &
         bb_d_d => this%bb_d_d, alpha_d => this%alpha_d)

      if (m .le. 0) return

      if (NEKO_DEVICE_MPI .and. (NEKO_BCKND_OPENCL .ne. 1)) then
         call device_project_ortho(alpha_d, bb_d(m), xx_d_d, bb_d_d, &
              w_d, xx_d(m), this%m, n, nrm)
      else
         if (NEKO_BCKND_OPENCL .eq. 1)then
            do i = 1, m
               alpha(i) = device_glsc3(bb_d(m), xx_d(i), w_d,n)
            end do
         else
            call device_glsc3_many(alpha, bb_d(m), xx_d_d, w_d, m, n)
         end if
         nrm = sqrt(alpha(m))
         call cmult(alpha, -1.0_rp,m)
         if (NEKO_BCKND_OPENCL .eq. 1)then
            do i = 1, m - 1
               call device_add2s2(xx_d(m), xx_d(i), alpha(i), n)
               call device_add2s2(bb_d(m), bb_d(i), alpha(i), n)

               alpha(i) = device_glsc3(bb_d(m), xx_d(i), w_d, n)
            end do
         else
            call device_memcpy(alpha, alpha_d, this%m, &
                 HOST_TO_DEVICE, sync = .false.)
            call device_add2s2_many(xx_d(m), xx_d_d, alpha_d, m-1, n)
            call device_add2s2_many(bb_d(m), bb_d_d, alpha_d, m-1, n)

            call device_glsc3_many(alpha, bb_d(m), xx_d_d, w_d, m, n)
         end if
         call cmult(alpha, -1.0_rp,m)
         if (NEKO_BCKND_OPENCL .eq. 1)then
            do i = 1, m - 1
               call device_add2s2(xx_d(m), xx_d(i), alpha(i), n)
               call device_add2s2(bb_d(m), bb_d(i), alpha(i), n)
               alpha(i) = device_glsc3(bb_d(m), xx_d(i), w_d, n)
            end do
         else
            call device_memcpy(alpha, alpha_d, m, &
                 HOST_TO_DEVICE, sync = .false.)
            call device_add2s2_many(xx_d(m), xx_d_d, alpha_d, m-1, n)
            call device_add2s2_many(bb_d(m), bb_d_d, alpha_d, m-1, n)
            call device_glsc3_many(alpha, bb_d(m), xx_d_d, w_d, m, n)
         end if
      end if

      alpha(m) = device_glsc3(xx_d(m), w_d, bb_d(m), n)
      alpha(m) = sqrt(alpha(m))

      if (alpha(m) .gt. this%tol*nrm) then !New vector is linearly independent
         scl = 1.0_rp / alpha(m)
         call device_cmult(xx_d(m), scl, n)
         call device_cmult(bb_d(m), scl, n)


      else !New vector is not linearly independent, forget about it
         if (pe_rank .eq. 0) then
            call neko_warning('New vector not linearly indepependent!')
         end if
         m = m - 1 !Remove column
      end if

    end associate

  end subroutine device_proj_ortho


  subroutine cpu_proj_ortho(this, xx, bb, w, n)
    type(projection_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n, this%L), intent(inout) :: xx, bb
    real(kind=rp), dimension(n), intent(in) :: w
    real(kind=rp) :: nrm, scl1, scl2, c, s
    real(kind=rp) :: alpha(this%L), beta(this%L)
    integer :: i, j, k, l, h, ierr

    associate(m => this%m)

      if (m .le. 0) return !No vectors to ortho-normalize

      ! AX = B
      ! Calculate dx, db: dx = x-XX^Tb, db=b-BX^Tb
      call rzero(alpha, m)
      do i = 1, n, NEKO_BLK_SIZE
         j = min(NEKO_BLK_SIZE, n-i+1)
         do k = 1, m !First round CGS
            s = 0.0_rp
            c = 0.0_rp
            do l = 0, (j-1)
               s = s + xx(i+l,k) * w(i+l) * bb(i+l,m)
               c = c + bb(i+l,k) * w(i+l) * xx(i+l,m)
            end do
            alpha(k) = alpha(k) + 0.5_rp * (s + c)
         end do
      end do

      call MPI_Allreduce(MPI_IN_PLACE, alpha, this%m, &
           MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)

      nrm = sqrt(alpha(m)) !Calculate A-norm of new vector


      do i = 1, n, NEKO_BLK_SIZE
         j = min(NEKO_BLK_SIZE, n-i+1)
         do k = 1, m-1
            do l = 0, (j-1)
               xx(i+l,m) = xx(i+l,m) - alpha(k) * xx(i+l,k)
               bb(i+l,m) = bb(i+l,m) - alpha(k) * bb(i+l,k)
            end do
         end do
      end do
      call rzero(beta,m)

      do i = 1, n, NEKO_BLK_SIZE
         j = min(NEKO_BLK_SIZE, n-i+1)
         do k = 1, m-1
            s = 0.0_rp
            c = 0.0_rp
            do l = 0, (j-1)
               s = s + xx(i+l,k) * w(i+l) * bb(i+l,m)
               c = c + bb(i+l,k) * w(i+l) * xx(i+l,m)
            end do
            beta(k) = beta(k) + 0.5_rp * (s + c)
         end do
      end do

      call MPI_Allreduce(MPI_IN_PLACE, beta, this%m-1, &
           MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)

      alpha(m) = 0.0_rp

      do i = 1, n, NEKO_BLK_SIZE
         j = min(NEKO_BLK_SIZE,n-i+1)
         do k = 1, m-1
            do l = 0, (j-1)
               xx(i+l,m) = xx(i+l,m) - beta(k) * xx(i+l,k)
               bb(i+l,m) = bb(i+l,m) - beta(k) * bb(i+l,k)
            end do
         end do
         s = 0.0_rp
         do l = 0, (j-1)
            s = s + xx(i+l,m) * w(i+l) * bb(i+l,m)
         end do
         alpha(m) = alpha(m) + s
      end do
      do k = 1, m-1
         alpha(k) = alpha(k) + beta(k)
      end do

      !alpha(m) = glsc3(xx(1,m), w, bb(1,m), n)
      call MPI_Allreduce(MPI_IN_PLACE, alpha(m), 1, &
           MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
      alpha(m) = sqrt(alpha(m))
      !dx and db now stored in last column of xx and bb

      if (alpha(m) .gt. this%tol*nrm) then !New vector is linearly independent
         !Normalize dx and db
         scl1 = 1.0_rp / alpha(m)
         do i = 0, (n - 1)
            xx(1+i,m) = scl1 * xx(1+i,m)
            bb(1+i,m) = scl1 * bb(1+i,m)
         end do

      else !New vector is not linearly independent, forget about it
         k = m !location of rank deficient column
         if (pe_rank .eq. 0) then
            call neko_warning('New vector not linearly indepependent!')
         end if
         m = m - 1 !Remove column
      end if

    end associate

  end subroutine cpu_proj_ortho

  subroutine print_proj_info(this, string)
    class(projection_t), intent(in) :: this
    character(len=*), intent(in) :: string
    character(len=LOG_SIZE) :: log_buf

    if (this%proj_m .gt. 0) then
       write(log_buf, '(A,A)') 'Projection ', string
       call neko_log%message(log_buf)
       write(log_buf, '(A,A)') 'Proj. vec.:', '   Orig. residual:'
       call neko_log%message(log_buf)
       write(log_buf, '(I11,3x, E15.7,5x)') this%proj_m, this%proj_res
       call neko_log%message(log_buf)
    end if

  end subroutine print_proj_info

  subroutine bcknd_clear(this, n)
    class(projection_t), intent(inout) :: this
    integer, intent(in) :: n
    integer :: i, j

    this%m = 0
    this%proj_m = 0

    do i = 1, this%L
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_rzero(this%xx_d(i), n)
          call device_rzero(this%xx_d(i), n)
       else
          do j = 1, n
             this%xx(j,i) = 0.0_rp
             this%bb(j,i) = 0.0_rp
          end do
       end if
    end do

  end subroutine bcknd_clear

end module projection
