!> Classic NEKTON formulation
!! Compute pressure and velocity using consistent approximation spaces.
module fluid_plan1
  use fluid_scheme_incompressible, only : fluid_scheme_incompressible_t, &
       fluid_scheme_solver_factory, fluid_scheme_precon_factory, &
       space_t, dofmap_t, gs_t, coef_t, gs_init, gs_free, space_free, &
       mesh_t, rp, time_step_controller_t, time_scheme_controller_t
  implicit none

  type, extends(fluid_scheme_incompressible_t) :: fluid_plan1_t
     type(space_t) :: Yh !< Function space for pressure \f$ Xh - 2 \f$
     type(dofmap_t) :: dm_Yh !< Dofmap associated with \f$ Yh \f$
     type(gs_t) :: gs_Yh !< Gather-scatter associated with \f$ Y_h \f$
     type(coef_t) :: c_Yh !< Coefficients associated with \f$ Y_h \f$
     !>@todo Remaning plan1 related data, ax, precon etc
   contains
     procedure, pass(this) :: init => fluid_plan1_init
     procedure, pass(this) :: free => fluid_plan1_free
     procedure, pass(this) :: step => fluid_plan1_step
  end type fluid_plan1_t

contains

  subroutine fluid_plan1_init(this, msh, lx, param)
    class(fluid_plan1_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(in) :: lx
    type(param_t), target, intent(inout) :: param
    character(len=15), parameter :: scheme = 'plan1 (Pn/Pn-2)'
    integer :: lx2

    call this%free()

    !> Setup velocity fields on the space \f$ Xh \f$
    call this%scheme_init(msh, lx, param, kspv_init=.true., scheme=scheme)

    !> Setup pressure field and related space \f$ Yh \f$
    lx2 = lx - 2
    if (msh%gdim .eq. 2) then
       call this%Yh%init(GLL, lx2, lx2)
    else
       call this%Yh%init(GLL, lx2, lx2, lx2)
    end if

    call this%dm_Yh%init(msh, this%Yh)

    call this%p%init(this%dm_Yh)

    call gs_init(this%gs_Yh, this%dm_Yh)

    call this%c_Yh%init(this%gs_Yh)

    call fluid_scheme_solver_factory(this%ksp_prs, this%dm_Yh%size(), &
         param%ksp_prs, param%abstol_prs)
    call fluid_scheme_precon_factory(this%pc_prs, this%ksp_prs, &
         this%c_Yh, this%dm_Yh, this%gs_Yh, this%bclst_prs, param%pc_prs)


  end subroutine fluid_plan1_init

  subroutine fluid_plan1_free(this)
    class(fluid_plan1_t), intent(inout) :: this

    ! Deallocate velocity and pressure fields
    call this%scheme_free()

    ! Deallocate Pressure dofmap and space
    call space_free(this%Yh)

    call gs_free(this%gs_Yh)

    call this%c_Yh%free()

  end subroutine fluid_plan1_free

  subroutine fluid_plan1_step(this, t, tstep, ext_bdf, dt_controller)
    class(fluid_plan1_t), target, intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(time_scheme_controller_t), intent(in) :: ext_bdf
    type(time_step_controller_t), intent(in) :: dt_controller

    if (this%freeze) return

  end subroutine fluid_plan1_step

end module fluid_plan1
