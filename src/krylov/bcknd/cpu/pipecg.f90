! Copyright (c) 2021-2024, The Neko Authors
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
!> Defines a pipelined Conjugate Gradient methods
module pipecg
  use krylov, only : ksp_t, ksp_monitor_t, KSP_MAX_ITER
  use precon, only : pc_t
  use ax_product, only : ax_t
  use num_types, only: rp
  use field, only : field_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use bc_list, only : bc_list_t
  use math, only : glsc3, rzero, copy, abscmp
  use comm
  implicit none
  private

  integer, parameter :: PIPECG_P_SPACE = 7

  !> Pipelined preconditioned conjugate gradient method
  type, public, extends(ksp_t) :: pipecg_t
     real(kind=rp), allocatable :: p(:)
     real(kind=rp), allocatable :: q(:)
     real(kind=rp), allocatable :: r(:)
     real(kind=rp), allocatable :: s(:)
     real(kind=rp), allocatable :: u(:,:)
     real(kind=rp), allocatable :: w(:)
     real(kind=rp), allocatable :: z(:)
     real(kind=rp), allocatable :: mi(:)
     real(kind=rp), allocatable :: ni(:)
   contains
     !> Constructor.
     procedure, pass(this) :: init => pipecg_init
     !> Destructor.
     procedure, pass(this) :: free => pipecg_free
     !> Solve the linear system.
     procedure, pass(this) :: solve => pipecg_solve
     !> Solve the coupled linear system.
     procedure, pass(this) :: solve_coupled => pipecg_solve_coupled
  end type pipecg_t

contains

  !> Initialise a pipelined PCG solver
  subroutine pipecg_init(this, n, max_iter, M, rel_tol, abs_tol, monitor)
    class(pipecg_t), target, intent(inout) :: this
    integer, intent(in) :: max_iter
    class(pc_t), optional, intent(in), target :: M
    integer, intent(in) :: n
    real(kind=rp), optional, intent(in) :: rel_tol
    real(kind=rp), optional, intent(in) :: abs_tol
    logical, optional, intent(in) :: monitor

    call this%free()

    allocate(this%p(n))
    allocate(this%q(n))
    allocate(this%r(n))
    allocate(this%s(n))
    allocate(this%u(n,PIPECG_P_SPACE+1))
    allocate(this%w(n))
    allocate(this%z(n))
    allocate(this%mi(n))
    allocate(this%ni(n))
    if (present(M)) then
       this%M => M
    end if

    if (present(rel_tol) .and. present(abs_tol) .and. present(monitor)) then
       call this%ksp_init(max_iter, rel_tol, abs_tol, monitor = monitor)
    else if (present(rel_tol) .and. present(abs_tol)) then
       call this%ksp_init(max_iter, rel_tol, abs_tol)
    else if (present(monitor) .and. present(abs_tol)) then
       call this%ksp_init(max_iter, abs_tol = abs_tol, monitor = monitor)
    else if (present(rel_tol) .and. present(monitor)) then
       call this%ksp_init(max_iter, rel_tol, monitor = monitor)
    else if (present(rel_tol)) then
       call this%ksp_init(max_iter, rel_tol = rel_tol)
    else if (present(abs_tol)) then
       call this%ksp_init(max_iter, abs_tol = abs_tol)
    else if (present(monitor)) then
       call this%ksp_init(max_iter, monitor = monitor)
    else
       call this%ksp_init(max_iter)
    end if

  end subroutine pipecg_init

  !> Deallocate a pipelined PCG solver
  subroutine pipecg_free(this)
    class(pipecg_t), intent(inout) :: this

    call this%ksp_free()

    if (allocated(this%p)) then
       deallocate(this%p)
    end if
    if (allocated(this%q)) then
       deallocate(this%q)
    end if
    if (allocated(this%r)) then
       deallocate(this%r)
    end if
    if (allocated(this%s)) then
       deallocate(this%s)
    end if
    if (allocated(this%u)) then
       deallocate(this%u)
    end if
    if (allocated(this%w)) then
       deallocate(this%w)
    end if
    if (allocated(this%z)) then
       deallocate(this%z)
    end if
    if (allocated(this%mi)) then
       deallocate(this%mi)
    end if
    if (allocated(this%ni)) then
       deallocate(this%ni)
    end if

    nullify(this%M)


  end subroutine pipecg_free

  !> Pipelined PCG solve
  function pipecg_solve(this, Ax, x, f, n, coef, blst, gs_h, niter) result(ksp_results)
    class(pipecg_t), intent(inout) :: this
    class(ax_t), intent(in) :: Ax
    type(field_t), intent(inout) :: x
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: f
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blst
    type(gs_t), intent(inout) :: gs_h
    type(ksp_monitor_t) :: ksp_results
    integer, optional, intent(in) :: niter
    integer :: iter, max_iter, i, j, k, ierr, p_cur, p_prev, u_prev
    real(kind=rp) :: rnorm, rtr, reduction(3), norm_fac
    real(kind=rp) :: alpha(PIPECG_P_SPACE), beta(PIPECG_P_SPACE)
    real(kind=rp) :: gamma1, gamma2, delta
    real(kind=rp) :: tmp1, tmp2, tmp3, x_plus(NEKO_BLK_SIZE)
    type(MPI_Request) :: request
    type(MPI_Status) :: status

    if (present(niter)) then
       max_iter = niter
    else
       max_iter = this%max_iter
    end if
    norm_fac = 1.0_rp / sqrt(coef%volume)

    associate(p => this%p, q => this%q, r => this%r, s => this%s, &
         u => this%u, w => this%w, z => this%z, mi => this%mi, ni => this%ni)

      p_prev = PIPECG_P_SPACE
      u_prev = PIPECG_P_SPACE+1
      p_cur = 1
      call rzero(x%x, n)
      call rzero(z, n)
      call rzero(q, n)
      call rzero(p, n)
      call rzero(s, n)
      call copy(r, f, n)
      call this%M%solve(u(1,u_prev), r, n)
      call Ax%compute(w, u(1,u_prev), coef, x%msh, x%Xh)
      call gs_h%op(w, n, GS_OP_ADD)
      call blst%apply(w, n)

      rtr = glsc3(r, coef%mult, r, n)
      rnorm = sqrt(rtr)*norm_fac
      ksp_results%res_start = rnorm
      ksp_results%res_final = rnorm
      ksp_results%iter = 0
      if(abscmp(rnorm, 0.0_rp)) return

      gamma1 = 0.0_rp
      tmp1 = 0.0_rp
      tmp2 = 0.0_rp
      tmp3 = 0.0_rp
      do i = 1, n
         tmp1 = tmp1 + r(i) * coef%mult(i,1,1,1) * u(i,u_prev)
         tmp2 = tmp2 + w(i) * coef%mult(i,1,1,1) * u(i,u_prev)
         tmp3 = tmp3 + r(i) * coef%mult(i,1,1,1) * r(i)
      end do
      reduction(1) = tmp1
      reduction(2) = tmp2
      reduction(3) = tmp3

      call this%monitor_start('PipeCG')
      do iter = 1, max_iter
         call MPI_Iallreduce(MPI_IN_PLACE, reduction, 3, &
              MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, request, ierr)

         call this%M%solve(mi, w, n)
         call Ax%compute(ni, mi, coef, x%msh, x%Xh)
         call gs_h%op(ni, n, GS_OP_ADD)
         call blst%apply(ni, n)

         call MPI_Wait(request, status, ierr)
         gamma2 = gamma1
         gamma1 = reduction(1)
         delta = reduction(2)
         rtr = reduction(3)

         rnorm = sqrt(rtr)*norm_fac
         call this%monitor_iter(iter, rnorm)
         if (rnorm .lt. this%abs_tol) exit

         if (iter .gt. 1) then
            beta(p_cur) = gamma1 / gamma2
            alpha(p_cur) = gamma1 / (delta - (beta(p_cur) * gamma1/alpha(p_prev)))
         else
            beta(p_cur) = 0.0_rp
            alpha(p_cur) = gamma1/delta
         end if

         tmp1 = 0.0_rp
         tmp2 = 0.0_rp
         tmp3 = 0.0_rp
         do i = 0, n, NEKO_BLK_SIZE
            if (i + NEKO_BLK_SIZE .le. n) then
               do k = 1, NEKO_BLK_SIZE
                  z(i+k) = beta(p_cur) * z(i+k) + ni(i+k)
                  q(i+k) = beta(p_cur) * q(i+k) + mi(i+k)
                  s(i+k) = beta(p_cur) * s(i+k) + w(i+k)
                  r(i+k) = r(i+k) - alpha(p_cur) * s(i+k)
                  u(i+k,p_cur) = u(i+k,u_prev) - alpha(p_cur) * q(i+k)
                  w(i+k) = w(i+k) - alpha(p_cur) * z(i+k)
                  tmp1 = tmp1 + r(i+k) * coef%mult(i+k,1,1,1) * u(i+k,p_cur)
                  tmp2 = tmp2 + w(i+k) * coef%mult(i+k,1,1,1) * u(i+k,p_cur)
                  tmp3 = tmp3 + r(i+k) * coef%mult(i+k,1,1,1) * r(i+k)
               end do
            else
               do k = 1, n-i
                  z(i+k) = beta(p_cur) * z(i+k) + ni(i+k)
                  q(i+k) = beta(p_cur) * q(i+k) + mi(i+k)
                  s(i+k) = beta(p_cur) * s(i+k) + w(i+k)
                  r(i+k) = r(i+k) - alpha(p_cur) * s(i+k)
                  u(i+k,p_cur) = u(i+k,u_prev) - alpha(p_cur) * q(i+k)
                  w(i+k) = w(i+k) - alpha(p_cur) * z(i+k)
                  tmp1 = tmp1 + r(i+k) * coef%mult(i+k,1,1,1) * u(i+k,p_cur)
                  tmp2 = tmp2 + w(i+k) * coef%mult(i+k,1,1,1) * u(i+k,p_cur)
                  tmp3 = tmp3 + r(i+k) * coef%mult(i+k,1,1,1) * r(i+k)
               end do
            end if
         end do

         reduction(1) = tmp1
         reduction(2) = tmp2
         reduction(3) = tmp3

         if (p_cur .eq. PIPECG_P_SPACE) then
            do i = 0, n, NEKO_BLK_SIZE
               if (i + NEKO_BLK_SIZE .le. n) then
                  do k = 1, NEKO_BLK_SIZE
                     x_plus(k) = 0.0_rp
                  end do
                  p_prev = PIPECG_P_SPACE+1
                  do j = 1, p_cur
                     do k = 1, NEKO_BLK_SIZE
                        p(i+k) = beta(j) * p(i+k) + u(i+k,p_prev)
                        x_plus(k) = x_plus(k) + alpha(j) * p(i+k)
                     end do
                     p_prev = j
                  end do
                  do k = 1, NEKO_BLK_SIZE
                     x%x(i+k,1,1,1) = x%x(i+k,1,1,1) + x_plus(k)
                     u(i+k,PIPECG_P_SPACE+1) = u(i+k,PIPECG_P_SPACE)
                  end do
               else
                  do k = 1, n-i
                     x_plus(1) = 0.0_rp
                     p_prev = PIPECG_P_SPACE + 1
                     do j = 1, p_cur
                        p(i+k) = beta(j) * p(i+k) + u(i+k,p_prev)
                        x_plus(1) = x_plus(1) + alpha(j) * p(i+k)
                        p_prev = j
                     end do
                     x%x(i+k,1,1,1) = x%x(i+k,1,1,1) + x_plus(1)
                     u(i+k,PIPECG_P_SPACE+1) = u(i+k,PIPECG_P_SPACE)
                  end do
               end if
            end do
            p_prev = p_cur
            u_prev = PIPECG_P_SPACE+1
            alpha(1) = alpha(p_cur)
            beta(1) = beta(p_cur)
            p_cur = 1
         else
            u_prev = p_cur
            p_prev = p_cur
            p_cur = p_cur + 1
         end if
      end do

      if ( p_cur .ne. 1) then
         do i = 0, n, NEKO_BLK_SIZE
            if (i + NEKO_BLK_SIZE .le. n) then
               do k = 1, NEKO_BLK_SIZE
                  x_plus(k) = 0.0_rp
               end do
               p_prev = PIPECG_P_SPACE+1
               do j = 1, p_cur
                  do k = 1, NEKO_BLK_SIZE
                     p(i+k) = beta(j) * p(i+k) + u(i+k,p_prev)
                     x_plus(k) = x_plus(k) + alpha(j) * p(i+k)
                  end do
                  p_prev = j
               end do
               do k = 1, NEKO_BLK_SIZE
                  x%x(i+k,1,1,1) = x%x(i+k,1,1,1) + x_plus(k)
                  u(i+k,PIPECG_P_SPACE+1) = u(i+k,PIPECG_P_SPACE)
               end do
            else
               do k = 1, n-i
                  x_plus(1) = 0.0_rp
                  p_prev = PIPECG_P_SPACE + 1
                  do j = 1, p_cur
                     p(i+k) = beta(j) * p(i+k) + u(i+k,p_prev)
                     x_plus(1) = x_plus(1) + alpha(j) * p(i+k)
                     p_prev = j
                  end do
                  x%x(i+k,1,1,1) = x%x(i+k,1,1,1) + x_plus(1)
                  u(i+k,PIPECG_P_SPACE+1) = u(i+k,PIPECG_P_SPACE)
               end do
            end if
         end do
      end if
      call this%monitor_stop()
      ksp_results%res_final = rnorm
      ksp_results%iter = iter
      ksp_results%converged = this%is_converged(iter, rnorm)

    end associate

  end function pipecg_solve

  !> Pipelined PCG coupled solve
  function pipecg_solve_coupled(this, Ax, x, y, z, fx, fy, fz, &
       n, coef, blstx, blsty, blstz, gs_h, niter) result(ksp_results)
    class(pipecg_t), intent(inout) :: this
    class(ax_t), intent(in) :: Ax
    type(field_t), intent(inout) :: x
    type(field_t), intent(inout) :: y
    type(field_t), intent(inout) :: z
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: fx
    real(kind=rp), dimension(n), intent(in) :: fy
    real(kind=rp), dimension(n), intent(in) :: fz
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blstx
    type(bc_list_t), intent(inout) :: blsty
    type(bc_list_t), intent(inout) :: blstz
    type(gs_t), intent(inout) :: gs_h
    type(ksp_monitor_t), dimension(3) :: ksp_results
    integer, optional, intent(in) :: niter

    ksp_results(1) = this%solve(Ax, x, fx, n, coef, blstx, gs_h, niter)
    ksp_results(2) = this%solve(Ax, y, fy, n, coef, blsty, gs_h, niter)
    ksp_results(3) = this%solve(Ax, z, fz, n, coef, blstz, gs_h, niter)

  end function pipecg_solve_coupled

end module pipecg


