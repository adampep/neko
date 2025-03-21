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
!> Defines various Bi-Conjugate Gradient Stabilized methods
module bicgstab
  use num_types, only: rp
  use krylov, only : ksp_t, ksp_monitor_t, KSP_MAX_ITER
  use precon, only : pc_t
  use ax_product, only : ax_t
  use field, only : field_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use bc_list, only : bc_list_t
  use math, only : glsc3, rzero, copy, NEKO_EPS, add2s2, x_update, &
       p_update, abscmp
  use utils, only : neko_error
  implicit none
  private

  !> Standard preconditioned Bi-Conjugate Gradient Stabilized method
  type, public, extends(ksp_t) :: bicgstab_t
     real(kind=rp), allocatable :: p(:)
     real(kind=rp), allocatable :: p_hat(:)
     real(kind=rp), allocatable :: r(:)
     real(kind=rp), allocatable :: s(:)
     real(kind=rp), allocatable :: s_hat(:)
     real(kind=rp), allocatable :: t(:)
     real(kind=rp), allocatable :: v(:)
   contains
     !> Constructor.
     procedure, pass(this) :: init => bicgstab_init
     !> Destructor.
     procedure, pass(this) :: free => bicgstab_free
     !> Solve the system.
     procedure, pass(this) :: solve => bicgstab_solve
     !> Solve the coupled system.
     procedure, pass(this) :: solve_coupled => bicgstab_solve_coupled
  end type bicgstab_t

contains

  !> Constructor.
  subroutine bicgstab_init(this, n, max_iter, M, rel_tol, abs_tol, monitor)
    class(bicgstab_t), target, intent(inout) :: this
    class(pc_t), optional, intent(in), target :: M
    integer, intent(in) :: n
    integer, intent(in) :: max_iter
    real(kind=rp), optional, intent(in) :: rel_tol
    real(kind=rp), optional, intent(in) :: abs_tol
    logical, optional, intent(in) :: monitor


    call this%free()

    allocate(this%p(n))
    allocate(this%p_hat(n))
    allocate(this%r(n))
    allocate(this%s(n))
    allocate(this%s_hat(n))
    allocate(this%t(n))
    allocate(this%v(n))
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

  end subroutine bicgstab_init

  !> Deallocate a standard BiCGSTAB solver
  subroutine bicgstab_free(this)
    class(bicgstab_t), intent(inout) :: this

    call this%ksp_free()

    if (allocated(this%v)) then
       deallocate(this%v)
    end if

    if (allocated(this%r)) then
       deallocate(this%r)
    end if

    if (allocated(this%t)) then
       deallocate(this%t)
    end if

    if (allocated(this%p)) then
       deallocate(this%p)
    end if

    if (allocated(this%p_hat)) then
       deallocate(this%p_hat)
    end if

    if (allocated(this%s)) then
       deallocate(this%s)
    end if

    if (allocated(this%s_hat)) then
       deallocate(this%s_hat)
    end if

    nullify(this%M)


  end subroutine bicgstab_free

  !> Bi-Conjugate Gradient Stabilized method solve
  function bicgstab_solve(this, Ax, x, f, n, coef, blst, gs_h, niter) result(ksp_results)
    class(bicgstab_t), intent(inout) :: this
    class(ax_t), intent(in) :: Ax
    type(field_t), intent(inout) :: x
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: f
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blst
    type(gs_t), intent(inout) :: gs_h
    type(ksp_monitor_t) :: ksp_results
    integer, optional, intent(in) :: niter
    integer :: iter, max_iter
    real(kind=rp) :: rnorm, rtr, norm_fac, gamma
    real(kind=rp) :: beta, alpha, omega, rho_1, rho_2

    if (present(niter)) then
       max_iter = niter
    else
       max_iter = this%max_iter
    end if
    norm_fac = 1.0_rp / sqrt(coef%volume)

    associate(r => this%r, t => this%t, s => this%s, v => this%v, p => this%p, &
         s_hat => this%s_hat, p_hat => this%p_hat)

      call rzero(x%x, n)
      call copy(r, f, n)

      rtr = sqrt(glsc3(r, coef%mult, r, n))
      rnorm = rtr * norm_fac
      gamma = rnorm * this%rel_tol
      ksp_results%res_start = rnorm
      ksp_results%res_final = rnorm
      ksp_results%iter = 0
      if(abscmp(rnorm, 0.0_rp)) return
      call this%monitor_start('BiCGStab')
      do iter = 1, max_iter

         rho_1 = glsc3(r, coef%mult, f ,n)

         if (abs(rho_1) .lt. NEKO_EPS) then
            call neko_error('Bi-CGStab rho failure')
         end if

         if (iter .eq. 1) then
            call copy(p, r, n)
         else
            beta = (rho_1 / rho_2) * (alpha / omega)
            call p_update(p, r, v, beta, omega, n)
         end if

         call this%M%solve(p_hat, p, n)
         call Ax%compute(v, p_hat, coef, x%msh, x%Xh)
         call gs_h%op(v, n, GS_OP_ADD)
         call blst%apply(v, n)
         alpha = rho_1 / glsc3(f, coef%mult, v, n)
         call copy(s, r, n)
         call add2s2(s, v, -alpha, n)
         rtr = glsc3(s, coef%mult, s, n)
         rnorm = sqrt(rtr) * norm_fac
         if (rnorm .lt. this%abs_tol .or. rnorm .lt. gamma) then
            call add2s2(x%x, p_hat, alpha,n)
            exit
         end if

         call this%M%solve(s_hat, s, n)
         call Ax%compute(t, s_hat, coef, x%msh, x%Xh)
         call gs_h%op(t, n, GS_OP_ADD)
         call blst%apply(t, n)
         omega = glsc3(t, coef%mult, s, n) &
              / glsc3(t, coef%mult, t, n)
         call x_update(x%x, p_hat, s_hat, alpha, omega, n)
         call copy(r, s, n)
         call add2s2(r, t, -omega, n)

         rtr = glsc3(r, coef%mult, r, n)
         rnorm = sqrt(rtr) * norm_fac
         call this%monitor_iter(iter, rnorm)
         if (rnorm .lt. this%abs_tol .or. rnorm .lt. gamma) then
            exit
         end if

         if (omega .lt. NEKO_EPS) then
            call neko_error('Bi-CGstab omega failure')
         end if
         rho_2 = rho_1

      end do
    end associate
    call this%monitor_stop()
    ksp_results%res_final = rnorm
    ksp_results%iter = iter
    ksp_results%converged = this%is_converged(iter, rnorm)
  end function bicgstab_solve

  !> Standard BiCGSTAB coupled solve
  function bicgstab_solve_coupled(this, Ax, x, y, z, fx, fy, fz, &
       n, coef, blstx, blsty, blstz, gs_h, niter) result(ksp_results)
    class(bicgstab_t), intent(inout) :: this
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

  end function bicgstab_solve_coupled

end module bicgstab


