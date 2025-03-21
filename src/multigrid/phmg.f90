! Copyright (c) 2024-2025, The Neko Authors
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
!> Hybrid ph-multigrid preconditioner
module phmg
  use num_types, only : rp
  use precon, only : pc_t
  use gather_scatter, only : gs_t, GS_OP_ADD
  use space, only : space_t, GLL
  use dofmap, only : dofmap_t
  use field, only : field_t
  use coefs, only : coef_t
  use mesh, only : mesh_t
  use bc, only : bc_t
  use bc_list, only : bc_list_t
  use dirichlet, only : dirichlet_t
  use utils, only : neko_error
  use cheby, only : cheby_t
  use cheby_device, only : cheby_device_t
  use jacobi, only : jacobi_t
  use device_jacobi, only : device_jacobi_t
  use ax_product, only : ax_t, ax_helm_factory
  use tree_amg_multigrid, only : tamg_solver_t
  use interpolation, only : interpolator_t
  use math, only : copy, col2, add2
  use device, only : device_get_ptr
  use device_math, only : device_rzero, device_copy, device_add2, device_sub3,&
       device_add2s2, device_invcol2, device_glsc2, device_col2
  use profiler, only : profiler_start_region, profiler_end_region
  use neko_config, only: NEKO_BCKND_DEVICE
  use krylov, only : ksp_t, ksp_monitor_t, KSP_MAX_ITER, &
       krylov_solver_factory
  use logger, only : neko_log, LOG_SIZE
  use, intrinsic :: iso_c_binding
  implicit none
  private


  type, private :: phmg_lvl_t
     integer :: lvl = -1
     type(space_t), pointer :: Xh
     type(dofmap_t), pointer :: dm_Xh
     type(gs_t), pointer :: gs_h
     type(cheby_t) :: cheby
     type(cheby_device_t) :: cheby_device
     type(jacobi_t) :: jacobi
     type(device_jacobi_t) :: device_jacobi
     type(coef_t), pointer :: coef
     type(bc_list_t) :: bclst
     type(dirichlet_t) :: bc
     type(field_t) :: r, w, z
  end type phmg_lvl_t

  type, public :: phmg_hrchy_t
     type(phmg_lvl_t), allocatable :: lvl(:)
  end type phmg_hrchy_t


  type, public, extends(pc_t) :: phmg_t
     type(tamg_solver_t) :: amg_solver
     integer :: nlvls
     type(phmg_hrchy_t) :: phmg_hrchy
     class(ax_t), allocatable :: ax
     type(interpolator_t), allocatable :: intrp(:)
     type(mesh_t), pointer :: msh
   contains
     procedure, pass(this) :: init => phmg_init
     procedure, pass(this) :: free => phmg_free
     procedure, pass(this) :: solve => phmg_solve
     procedure, pass(this) :: update => phmg_update
  end type phmg_t

contains

  subroutine phmg_init(this, msh, Xh, coef, dof, gs_h, bclst)
    class(phmg_t), intent(inout), target :: this
    type(mesh_t), intent(inout), target :: msh
    type(space_t), intent(inout), target :: Xh
    type(coef_t), intent(in), target :: coef
    type(dofmap_t), intent(in), target :: dof
    type(gs_t), intent(inout), target :: gs_h
    type(bc_list_t), intent(inout), target :: bclst
    integer :: lx_crs, lx_mid
    integer, allocatable :: lx_lvls(:)
    integer :: n, i, j
    class(bc_t), pointer :: bc_j
    logical :: use_jacobi, use_cheby
    use_jacobi = .false.
    use_cheby = .true.

    this%msh => msh

    this%nlvls = Xh%lx - 1

    allocate(lx_lvls(0:this%nlvls - 1))
    do i = 1, this%nlvls -1
       lx_lvls(i) = Xh%lx - i
    end do

    allocate(this%phmg_hrchy%lvl(0:this%nlvls - 1))

    this%phmg_hrchy%lvl(0)%lvl = 0
    this%phmg_hrchy%lvl(0)%Xh => Xh
    this%phmg_hrchy%lvl(0)%coef => coef
    this%phmg_hrchy%lvl(0)%dm_Xh => dof
    this%phmg_hrchy%lvl(0)%gs_h => gs_h

    do i = 1, this%nlvls - 1
       allocate(this%phmg_hrchy%lvl(i)%Xh)
       allocate(this%phmg_hrchy%lvl(i)%dm_Xh)
       allocate(this%phmg_hrchy%lvl(i)%gs_h)
       allocate(this%phmg_hrchy%lvl(i)%coef)

       this%phmg_hrchy%lvl(i)%lvl = i
       call this%phmg_hrchy%lvl(i)%Xh%init(GLL, lx_lvls(i), lx_lvls(i), &
            lx_lvls(i))
       call this%phmg_hrchy%lvl(i)%dm_Xh%init(msh, this%phmg_hrchy%lvl(i)%Xh)
       call this%phmg_hrchy%lvl(i)%gs_h%init(this%phmg_hrchy%lvl(i)%dm_Xh)
       call this%phmg_hrchy%lvl(i)%coef%init(this%phmg_hrchy%lvl(i)%gs_h)
    end do

    do i = 0, this%nlvls - 1
       call this%phmg_hrchy%lvl(i)%r%init(this%phmg_hrchy%lvl(i)%dm_Xh)
       call this%phmg_hrchy%lvl(i)%w%init(this%phmg_hrchy%lvl(i)%dm_Xh)
       call this%phmg_hrchy%lvl(i)%z%init(this%phmg_hrchy%lvl(i)%dm_Xh)

       if (use_cheby) then
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call this%phmg_hrchy%lvl(i)%cheby_device%init(this%phmg_hrchy%lvl(i)%dm_Xh%size(), KSP_MAX_ITER)
          else
             call this%phmg_hrchy%lvl(i)%cheby%init(this%phmg_hrchy%lvl(i)%dm_Xh%size(), KSP_MAX_ITER)
          end if
       end if

       if (use_jacobi) then
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call this%phmg_hrchy%lvl(i)%device_jacobi%init(this%phmg_hrchy%lvl(i)%coef, &
                  this%phmg_hrchy%lvl(i)%dm_Xh, &
                  this%phmg_hrchy%lvl(i)%gs_h)
          else
             call this%phmg_hrchy%lvl(i)%jacobi%init(this%phmg_hrchy%lvl(i)%coef, &
                  this%phmg_hrchy%lvl(i)%dm_Xh, &
                  this%phmg_hrchy%lvl(i)%gs_h)
          end if
       end if

       this%phmg_hrchy%lvl(i)%coef%ifh2 = coef%ifh2
       call copy(this%phmg_hrchy%lvl(i)%coef%h1, coef%h1, &
            this%phmg_hrchy%lvl(i)%dm_Xh%size())

       call this%phmg_hrchy%lvl(i)%bc%init_base(this%phmg_hrchy%lvl(i)%coef)
       if (bclst%size() .gt. 0 ) then
          do j = 1, bclst%size()
             bc_j => bclst%get(j)
             call this%phmg_hrchy%lvl(i)%bc%mark_facets(bc_j%marked_facet)
          end do
       end if
       call this%phmg_hrchy%lvl(i)%bc%finalize()
       call this%phmg_hrchy%lvl(i)%bc%set_g(0.0_rp)
       call this%phmg_hrchy%lvl(i)%bclst%init()
       call this%phmg_hrchy%lvl(i)%bclst%append(this%phmg_hrchy%lvl(i)%bc)
    end do

    ! Create backend specific Ax operator
    call ax_helm_factory(this%ax, full_formulation = .false.)

    ! Interpolator Fine + mg levels
    allocate(this%intrp(this%nlvls - 1))
    do i = 1, this%nlvls -1
       call this%intrp(i)%init(this%phmg_hrchy%lvl(i-1)%Xh, &
            this%phmg_hrchy%lvl(i)%Xh)
    end do

    call this%amg_solver%init(this%ax, this%phmg_hrchy%lvl(this%nlvls -1)%Xh, &
         this%phmg_hrchy%lvl(this%nlvls -1)%coef, this%msh, &
         this%phmg_hrchy%lvl(this%nlvls-1)%gs_h, 4, &
         this%phmg_hrchy%lvl(this%nlvls -1)%bclst, 1)

  end subroutine phmg_init

  subroutine phmg_free(this)
    class(phmg_t), intent(inout) :: this
  end subroutine phmg_free

  subroutine phmg_solve(this, z, r, n)
    class(phmg_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: z
    real(kind=rp), dimension(n), intent(inout) :: r
    type(c_ptr) :: z_d, r_d
    type(ksp_monitor_t) :: ksp_results


    associate( mglvl => this%phmg_hrchy%lvl)
      if (NEKO_BCKND_DEVICE .eq. 1) then
         z_d = device_get_ptr(z)
         r_d = device_get_ptr(r)
         !We should not work with the input
         call device_copy(mglvl(0)%r%x_d, r_d, n)
         call device_rzero(mglvl(0)%z%x_d, n)
         call device_rzero(mglvl(0)%w%x_d, n)
         call phmg_mg_cycle(mglvl(0)%z, mglvl(0)%r, mglvl(0)%w, 0, this%nlvls -1, &
              mglvl, this%intrp, this%msh, this%Ax, this%amg_solver)
         call device_copy(z_d, mglvl(0)%z%x_d, n)
      else
         !We should not work with the input
         call copy(mglvl(0)%r%x, r, n)

         mglvl(0)%z%x = 0.0_rp
         mglvl(0)%w%x = 0.0_rp

         call phmg_mg_cycle(mglvl(0)%z, mglvl(0)%r, mglvl(0)%w, 0, this%nlvls -1, &
              mglvl, this%intrp, this%msh, this%Ax, this%amg_solver)

         call copy(z, mglvl(0)%z%x, n)
      end if
    end associate

  end subroutine phmg_solve

  subroutine phmg_update(this)
    class(phmg_t), intent(inout) :: this
  end subroutine phmg_update


  recursive subroutine phmg_mg_cycle(z, r, w, lvl, clvl, &
       mg, intrp, msh, Ax, amg_solver)
    type(ksp_monitor_t) :: ksp_results
    integer :: lvl, clvl
    type(phmg_lvl_t) :: mg(0:clvl)
    type(interpolator_t) :: intrp(1:clvl)
    type(tamg_solver_t), intent(inout) :: amg_solver
    class(ax_t), intent(inout) :: Ax
    type(mesh_t), intent(inout) :: msh
    type(field_t) :: z, r, w
    integer :: i
    real(kind=rp) :: val

    call profiler_start_region('PHMG_cycle', 8)
    !>----------<!
    !> SMOOTH   <!
    !>----------<!
    call profiler_start_region('PHMG_PreSmooth', 9)
    !!!!!!--------------------------------------------------------------------------------------------------------------
    !call phmg_jacobi_smoother(z, r, w, mg(lvl), msh, Ax, mg(lvl)%dm_Xh%size(), lvl)
    !!!!!!--------------------------------------------------------------------------------------------------------------
    if (NEKO_BCKND_DEVICE .eq. 1) then
       ksp_results = mg(lvl)%cheby_device%solve(Ax, z, &
            r%x, mg(lvl)%dm_Xh%size(), &
            mg(lvl)%coef, mg(lvl)%bclst, &
            mg(lvl)%gs_h, niter = 11)
    else
       ksp_results = mg(lvl)%cheby%solve(Ax, z, &
            r%x, mg(lvl)%dm_Xh%size(), &
            mg(lvl)%coef, mg(lvl)%bclst, &
            mg(lvl)%gs_h, niter = 15)
    end if
    call profiler_end_region('PHMG_PreSmooth', 9)

    !>----------<!
    !> Residual <!
    !>----------<!
    call Ax%compute(w%x, z%x, mg(lvl)%coef, msh, mg(lvl)%Xh)
    call mg(lvl)%gs_h%op(w%x, mg(lvl)%dm_Xh%size(), GS_OP_ADD)
    !call mg(lvl)%bclst%apply_scalar(w%x, mg(lvl)%dm_Xh%size())

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_sub3(w%x_d, r%x_d, w%x_d, mg(lvl)%dm_Xh%size())
    else
       w%x = r%x - w%x
    end if

    !>----------<!
    !> Restrict <!
    !>----------<!
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col2(w%x_d, mg(lvl)%coef%mult_d, mg(lvl)%dm_Xh%size())
    else
       call col2(w%x, mg(lvl)%coef%mult, mg(lvl)%dm_Xh%size())
    end if

    call profiler_start_region('PHMG_map_to_coarse', 9)
    call intrp(lvl+1)%map(mg(lvl+1)%r%x, w%x, msh%nelv, mg(lvl+1)%Xh)
    call profiler_end_region('PHMG_map_to_coarse', 9)

    call mg(lvl+1)%gs_h%op(mg(lvl+1)%r%x, mg(lvl+1)%dm_Xh%size(), GS_OP_ADD)

    call mg(lvl+1)%bclst%apply_scalar( &
         mg(lvl+1)%r%x, &
         mg(lvl+1)%dm_Xh%size())
    !>----------<!
    !> SOLVE    <!
    !>----------<!
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_rzero(mg(lvl+1)%z%x_d, mg(lvl+1)%dm_Xh%size())
    else
       mg(lvl+1)%z%x = 0.0_rp
    end if
    if (lvl+1 .eq. clvl) then
       call profiler_start_region('PHMG_tAMG_coarse_grid', 10)

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call amg_solver%device_solve(mg(lvl+1)%z%x, &
               mg(lvl+1)%r%x, &
               mg(lvl+1)%z%x_d, &
               mg(lvl+1)%r%x_d, &
               mg(lvl+1)%dm_Xh%size())
       else
          call amg_solver%solve(mg(lvl+1)%z%x, &
               mg(lvl+1)%r%x, &
               mg(lvl+1)%dm_Xh%size())
       end if
       call profiler_end_region('PHMG_tAMG_coarse_grid', 10)

       call mg(lvl+1)%bclst%apply_scalar( &
            mg(lvl+1)%z%x,&
            mg(lvl+1)%dm_Xh%size())
    else
       call phmg_mg_cycle(mg(lvl+1)%z, mg(lvl+1)%r, mg(lvl+1)%w, lvl+1, &
            clvl, mg, intrp, msh, Ax, amg_solver)
    end if

    !>----------<!
    !> Project  <!
    !>----------<!
    call profiler_start_region('PHMG_map_to_fine', 9)
    call intrp(lvl+1)%map(w%x, mg(lvl+1)%z%x, msh%nelv, mg(lvl)%Xh)
    call profiler_end_region('PHMG_map_to_fine', 9)

    call mg(lvl)%gs_h%op(w%x, mg(lvl)%dm_Xh%size(), GS_OP_ADD)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col2(w%x_d, mg(lvl)%coef%mult_d, mg(lvl)%dm_Xh%size())
    else
       call col2(w%x, mg(lvl)%coef%mult, mg(lvl)%dm_Xh%size())
    end if

    call mg(lvl)%bclst%apply_scalar(w%x, mg(lvl)%dm_Xh%size())

    !>----------<!
    !> Correct  <!
    !>----------<!
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_add2(z%x_d, w%x_d, mg(lvl)%dm_Xh%size())
    else
       z%x = z%x + w%x
    end if

    !>----------<!
    !> SMOOTH   <!
    !>----------<!
    call profiler_start_region('PHMG_PostSmooth', 9)
    !!!!!!--------------------------------------------------------------------------------------------------------------
    !call phmg_jacobi_smoother(z, r, w, mg(lvl), msh, Ax, mg(lvl)%dm_Xh%size(), lvl)
    !!!!!!--------------------------------------------------------------------------------------------------------------
    if (NEKO_BCKND_DEVICE .eq. 1) then
       ksp_results = mg(lvl)%cheby_device%solve(Ax, z, &
            r%x, mg(lvl)%dm_Xh%size(), &
            mg(lvl)%coef, mg(lvl)%bclst, &
            mg(lvl)%gs_h, niter = 11)
    else
       ksp_results = mg(lvl)%cheby%solve(Ax, z, &
            r%x, mg(lvl)%dm_Xh%size(), &
            mg(lvl)%coef, mg(lvl)%bclst, &
            mg(lvl)%gs_h, niter = 15)
    end if
    call profiler_end_region('PHMG_PostSmooth', 9)

    call profiler_end_region('PHMG_cycle', 8)
  end subroutine phmg_mg_cycle

  subroutine phmg_jacobi_smoother(z, r, w, mg, msh, Ax, n, lvl)
    type(phmg_lvl_t) :: mg
    class(ax_t), intent(inout) :: Ax
    type(mesh_t), intent(inout) :: msh
    type(field_t), intent(inout) :: z, r, w
    integer, intent(in) :: n, lvl
    integer :: i, iblk, ni, niblk

    ni = 1
    niblk = 3

    do i = 1, ni
       call Ax%compute(w%x, z%x, mg%coef, msh, mg%Xh)
       call mg%gs_h%op(w%x, n, GS_OP_ADD)
       call mg%bclst%apply_scalar(w%x, n)
       call device_sub3(w%x_d, r%x_d, w%x_d, n)

       call mg%device_jacobi%solve(w%x, w%x, n)

       call device_add2s2(z%x_d, w%x_d, 0.8_rp, n)

       do iblk = 1, niblk
          call Ax%compute(w%x, z%x, mg%coef, msh, mg%Xh)
          call device_invcol2(w%x_d, mg%coef%mult_d, n)
          call device_sub3(w%x_d, r%x_d, w%x_d, n)
          call mg%device_jacobi%solve(w%x, w%x, n)
          call device_add2s2(z%x_d, w%x_d, 0.8_rp, n)
       end do
    end do
  end subroutine phmg_jacobi_smoother


  subroutine phmg_resid_monitor(z, r, w, mg, msh, Ax, lvl, typ)
    integer :: lvl, typ
    type(phmg_lvl_t) :: mg
    class(ax_t), intent(inout) :: Ax
    type(mesh_t), intent(inout) :: msh
    type(field_t) :: z, r, w
    real(kind=rp) :: val
    character(len=LOG_SIZE) :: log_buf
    call Ax%compute(w%x, z%x, mg%coef, msh, mg%Xh)
    call mg%gs_h%op(w%x, mg%dm_Xh%size(), GS_OP_ADD)
    call mg%bclst%apply_scalar(w%x, mg%dm_Xh%size())
    call device_sub3(w%x_d, r%x_d, w%x_d, mg%dm_Xh%size())
    val = device_glsc2(w%x_d, w%x_d, mg%dm_Xh%size())
    if (typ .eq. 1) then
       write(log_buf, '(A15,I4,F12.6)') 'PRESMOO - PRE', lvl, val
    else if (typ .eq. 2) then
       write(log_buf, '(A15,I4,F12.6)') 'PRESMOO -POST', lvl, val
    else if (typ .eq. 3) then
       write(log_buf, '(A15,I4,F12.6)') 'POSTSMOO- PRE', lvl, val
    else if (typ .eq. 4) then
       write(log_buf, '(A15,I4,F12.6)') 'POSTSMOO-POST', lvl, val
    else if (typ .eq. 5) then
       write(log_buf, '(A15,I4,F12.6)') 'TAMG - PRE', lvl, val
    else if (typ .eq. 6) then
       write(log_buf, '(A15,I4,F12.6)') 'TAMG -POST', lvl, val
    else
       write(log_buf, '(A15,I4,F12.6)') 'RESID', lvl, val
    end if
    call neko_log%message(log_buf)
  end subroutine phmg_resid_monitor

end module phmg
