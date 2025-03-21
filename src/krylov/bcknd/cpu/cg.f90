! Copyright (c) 2020-2024, The Neko Authors
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
!> Defines various Conjugate Gradient methods
module cg
  use num_types, only: rp, xp
  use krylov, only : ksp_t, ksp_monitor_t, KSP_MAX_ITER
  use precon, only : pc_t
  use ax_product, only : ax_t
  use field, only : field_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use bc_list, only : bc_list_t
  use math, only : glsc3, rzero, copy, abscmp
  use comm
  implicit none
  private

  integer, parameter :: CG_P_SPACE = 7

  !> Standard preconditioned conjugate gradient method
  type, public, extends(ksp_t) :: cg_t
     real(kind=rp), allocatable :: w(:)
     real(kind=rp), allocatable :: r(:)
     real(kind=rp), allocatable :: p(:,:)
     real(kind=rp), allocatable :: z(:)
     real(kind=rp), allocatable :: alpha(:)
   contains
     procedure, pass(this) :: init => cg_init
     procedure, pass(this) :: free => cg_free
     procedure, pass(this) :: solve => cg_solve
     procedure, pass(this) :: solve_coupled => cg_solve_coupled
  end type cg_t

contains

  !> Initialise a standard PCG solver
  subroutine cg_init(this, n, max_iter, M, rel_tol, abs_tol, monitor)
    class(cg_t), intent(inout), target :: this
    integer, intent(in) :: max_iter
    class(pc_t), optional, intent(in), target :: M
    integer, intent(in) :: n
    real(kind=rp), optional, intent(in) :: rel_tol
    real(kind=rp), optional, intent(in) :: abs_tol
    logical, optional, intent(in) :: monitor

    call this%free()

    allocate(this%w(n))
    allocate(this%r(n))
    allocate(this%p(n,CG_P_SPACE))
    allocate(this%z(n))
    allocate(this%alpha(CG_P_SPACE))

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

  end subroutine cg_init

  !> Deallocate a standard PCG solver
  subroutine cg_free(this)
    class(cg_t), intent(inout) :: this

    call this%ksp_free()

    if (allocated(this%w)) then
       deallocate(this%w)
    end if

    if (allocated(this%r)) then
       deallocate(this%r)
    end if

    if (allocated(this%p)) then
       deallocate(this%p)
    end if

    if (allocated(this%z)) then
       deallocate(this%z)
    end if

    if (allocated(this%alpha)) then
       deallocate(this%alpha)
    end if

    nullify(this%M)

  end subroutine cg_free

  !> Standard PCG solve
  function cg_solve(this, Ax, x, f, n, coef, blst, gs_h, niter) result(ksp_results)
    class(cg_t), intent(inout) :: this
    class(ax_t), intent(in) :: Ax
    type(field_t), intent(inout) :: x
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(in) :: f
    type(coef_t), intent(inout) :: coef
    type(bc_list_t), intent(inout) :: blst
    type(gs_t), intent(inout) :: gs_h
    type(ksp_monitor_t) :: ksp_results
    integer, optional, intent(in) :: niter
    integer :: iter, max_iter, i, j, k, p_cur, p_prev
    real(kind=rp) :: rnorm, rtr, rtz2, rtz1, x_plus(NEKO_BLK_SIZE)
    real(kind=rp) :: beta, pap, norm_fac

    if (present(niter)) then
       max_iter = niter
    else
       max_iter = this%max_iter
    end if
    norm_fac = 1.0_rp / sqrt(coef%volume)

    associate(w => this%w, r => this%r, p => this%p, &
         z => this%z, alpha => this%alpha)

      rtz1 = 1.0_rp
      call rzero(x%x, n)
      call rzero(p(1,CG_P_SPACE), n)
      call copy(r, f, n)

      rtr = glsc3(r, coef%mult, r, n)
      rnorm = sqrt(rtr) * norm_fac
      ksp_results%res_start = rnorm
      ksp_results%res_final = rnorm
      ksp_results%iter = 0
      p_prev = CG_P_SPACE
      p_cur = 1
      if(abscmp(rnorm, 0.0_rp)) return
      call this%monitor_start('CG')
      do iter = 1, max_iter
         call this%M%solve(z, r, n)
         rtz2 = rtz1
         rtz1 = glsc3(r, coef%mult, z, n)

         beta = rtz1 / rtz2
         if (iter .eq. 1) beta = 0.0_rp
         do i = 1, n
            p(i,p_cur) = z(i) + beta * p(i,p_prev)
         end do

         call Ax%compute(w, p(1,p_cur), coef, x%msh, x%Xh)
         call gs_h%op(w, n, GS_OP_ADD)
         call blst%apply(w, n)

         pap = glsc3(w, coef%mult, p(1,p_cur), n)

         alpha(p_cur) = rtz1 / pap
         call second_cg_part(rtr, r, coef%mult, w, alpha(p_cur), n)
         rnorm = sqrt(rtr) * norm_fac
         call this%monitor_iter(iter, rnorm)

         if ((p_cur .eq. CG_P_SPACE) .or. &
              (rnorm .lt. this%abs_tol) .or. iter .eq. max_iter) then
            do i = 0, n, NEKO_BLK_SIZE
               if (i + NEKO_BLK_SIZE .le. n) then
                  do k = 1, NEKO_BLK_SIZE
                     x_plus(k) = 0.0_rp
                  end do
                  do j = 1, p_cur
                     do k = 1, NEKO_BLK_SIZE
                        x_plus(k) = x_plus(k) + alpha(j) * p(i+k,j)
                     end do
                  end do
                  do k = 1, NEKO_BLK_SIZE
                     x%x(i+k,1,1,1) = x%x(i+k,1,1,1) + x_plus(k)
                  end do
               else
                  do k = 1, n-i
                     x_plus(1) = 0.0_rp
                     do j = 1, p_cur
                        x_plus(1) = x_plus(1) + alpha(j) * p(i+k,j)
                     end do
                     x%x(i+k,1,1,1) = x%x(i+k,1,1,1) + x_plus(1)
                  end do
               end if
            end do
            p_prev = p_cur
            p_cur = 1
            if (rnorm .lt. this%abs_tol) exit
         else
            p_prev = p_cur
            p_cur = p_cur + 1
         end if
      end do
    end associate
    call this%monitor_stop()
    ksp_results%res_final = rnorm
    ksp_results%iter = iter
    ksp_results%converged = this%is_converged(iter, rnorm)
  end function cg_solve

  subroutine second_cg_part(rtr, r, mult, w, alpha, n)
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: r(n), rtr
    real(kind=xp) :: tmp
    real(kind=rp), intent(in) ::mult(n), w(n), alpha
    integer :: i, ierr

    tmp = 0.0_xp
    do i = 1, n
       r(i) = r(i) - alpha*w(i)
       tmp = tmp + r(i) * r(i) * mult(i)
    end do
    call MPI_Allreduce(MPI_IN_PLACE, tmp, 1, &
         MPI_EXTRA_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    rtr = tmp

  end subroutine second_cg_part

  !> Standard PCG coupled solve
  function cg_solve_coupled(this, Ax, x, y, z, fx, fy, fz, &
       n, coef, blstx, blsty, blstz, gs_h, niter) result(ksp_results)
    class(cg_t), intent(inout) :: this
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

  end function cg_solve_coupled

end module cg


