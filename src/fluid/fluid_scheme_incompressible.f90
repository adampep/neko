! Copyright (c) 2020-2025, The Neko Authors
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
!> Fluid formulations
module fluid_scheme_incompressible
  use fluid_scheme_base, only : fluid_scheme_base_t
  use gather_scatter, only : gs_t, GS_OP_MIN, GS_OP_MAX
  use mean_sqr_flow, only : mean_sqr_flow_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use checkpoint, only : chkp_t
  use mean_flow, only : mean_flow_t
  use num_types, only : rp, i8
  use comm
  use fluid_source_term, only: fluid_source_term_t
  use field, only : field_t
  use space, only : space_t, GLL
  use dofmap, only : dofmap_t
  use zero_dirichlet, only : zero_dirichlet_t
  use krylov, only : ksp_t, krylov_solver_factory, KSP_MAX_ITER
  use coefs, only: coef_t
  use usr_inflow, only : usr_inflow_t, usr_inflow_eval
  use dirichlet, only : dirichlet_t
  use field_dirichlet, only : field_dirichlet_t
  use field_dirichlet_vector, only: field_dirichlet_vector_t
  use jacobi, only : jacobi_t
  use sx_jacobi, only : sx_jacobi_t
  use device_jacobi, only : device_jacobi_t
  use hsmg, only : hsmg_t
  use phmg, only : phmg_t
  use precon, only : pc_t, precon_factory, precon_destroy
  use fluid_stats, only : fluid_stats_t
  use bc, only : bc_t
  use bc_list, only : bc_list_t
  use mesh, only : mesh_t, NEKO_MSH_MAX_ZLBLS, NEKO_MSH_MAX_ZLBL_LEN
  use math, only : cfill, add2s2, glsum
  use device_math, only : device_cfill, device_add2s2
  use operators, only : cfl
  use logger, only : neko_log, LOG_SIZE, NEKO_LOG_VERBOSE
  use field_registry, only : neko_field_registry
  use json_utils, only : json_get, json_get_or_default, json_extract_object, &
       json_extract_item
  use json_module, only : json_file, json_core, json_value
  use scratch_registry, only : scratch_registry_t
  use user_intf, only : user_t, dummy_user_material_properties, &
       user_material_properties
  use utils, only : neko_error, neko_warning
  use field_series, only : field_series_t
  use time_step_controller, only : time_step_controller_t
  use field_math, only : field_cfill, field_add2s2
  use shear_stress, only : shear_stress_t
  use gradient_jump_penalty, only : gradient_jump_penalty_t
  use device, only : device_event_sync, device_stream_wait_event, &
       glb_cmd_queue, glb_cmd_event
  implicit none
  private

  !> Base type of all fluid formulations
  type, abstract, extends(fluid_scheme_base_t) :: fluid_scheme_incompressible_t
     !> The source term for the momentum equation.
     type(fluid_source_term_t) :: source_term
     class(ksp_t), allocatable :: ksp_vel !< Krylov solver for velocity
     class(ksp_t), allocatable :: ksp_prs !< Krylov solver for pressure
     class(pc_t), allocatable :: pc_vel !< Velocity Preconditioner
     class(pc_t), allocatable :: pc_prs !< Velocity Preconditioner
     integer :: vel_projection_dim !< Size of the projection space for ksp_vel
     integer :: pr_projection_dim !< Size of the projection space for ksp_pr
     integer :: vel_projection_activ_step !< Steps to activate projection for ksp_vel
     integer :: pr_projection_activ_step !< Steps to activate projection for ksp_pr
     logical :: strict_convergence !< Strict convergence for the velocity solver
     !> Gradient jump panelty
     logical :: if_gradient_jump_penalty
     type(gradient_jump_penalty_t) :: gradient_jump_penalty_u
     type(gradient_jump_penalty_t) :: gradient_jump_penalty_v
     type(gradient_jump_penalty_t) :: gradient_jump_penalty_w
     !> Extrapolation velocity fields for LES
     type(field_t), pointer :: u_e => null() !< Extrapolated x-Velocity
     type(field_t), pointer :: v_e => null() !< Extrapolated y-Velocity
     type(field_t), pointer :: w_e => null() !< Extrapolated z-Velocity

     type(mean_flow_t) :: mean !< Mean flow field
     type(fluid_stats_t) :: stats !< Fluid statistics
     type(mean_sqr_flow_t) :: mean_sqr !< Mean squared flow field
     logical :: forced_flow_rate = .false. !< Is the flow rate forced?

     !> The turbulent kinematic viscosity field name
     character(len=:), allocatable :: nut_field_name

     !> Global number of GLL points for the fluid (not unique)
     integer(kind=i8) :: glb_n_points
     !> Global number of GLL points for the fluid (unique)
     integer(kind=i8) :: glb_unique_points
     type(scratch_registry_t) :: scratch !< Manager for temporary fields
   contains
     !> Constructor for the base type
     procedure, pass(this) :: init_base => fluid_scheme_init_base
     procedure, pass(this) :: scheme_free => fluid_scheme_free
     !> Validate that all components are properly allocated
     procedure, pass(this) :: validate => fluid_scheme_validate
     !> Apply pressure boundary conditions
     procedure, pass(this) :: bc_apply_vel => fluid_scheme_bc_apply_vel
     !> Apply velocity boundary conditions
     procedure, pass(this) :: bc_apply_prs => fluid_scheme_bc_apply_prs
     !> Compute the CFL number
     procedure, pass(this) :: compute_cfl => fluid_compute_cfl
     !> Set rho and mu
     procedure, pass(this) :: set_material_properties => &
          fluid_scheme_set_material_properties

     !> Update variable material properties
     procedure, pass(this) :: update_material_properties => &
          fluid_scheme_update_material_properties
     !> Linear solver factory, wraps a KSP constructor
     procedure, nopass :: solver_factory => fluid_scheme_solver_factory
     !> Preconditioner factory
     procedure, pass(this) :: precon_factory_ => fluid_scheme_precon_factory
  end type fluid_scheme_incompressible_t

  interface
     !> Initialise a fluid scheme
     module subroutine fluid_scheme_factory(object, type_name)
       class(fluid_scheme_base_t), intent(inout), allocatable :: object
       character(len=*) :: type_name
     end subroutine fluid_scheme_factory
  end interface

  public :: fluid_scheme_incompressible_t, fluid_scheme_factory

contains

  !> Initialize common data for the current scheme
  subroutine fluid_scheme_init_base(this, msh, lx, params, scheme, user, &
       kspv_init)
    implicit none
    class(fluid_scheme_incompressible_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(in) :: lx
    character(len=*), intent(in) :: scheme
    type(json_file), target, intent(inout) :: params
    type(user_t), target, intent(in) :: user
    logical, intent(in) :: kspv_init
    type(dirichlet_t) :: bdry_mask
    character(len=LOG_SIZE) :: log_buf
    real(kind=rp), allocatable :: real_vec(:)
    real(kind=rp) :: real_val, kappa, B, z0
    logical :: logical_val
    integer :: integer_val, ierr
    type(json_file) :: wm_json
    character(len=:), allocatable :: string_val1, string_val2
    real(kind=rp) :: GJP_param_a, GJP_param_b

    !
    ! SEM simulation fundamentals
    !

    this%msh => msh

    if (msh%gdim .eq. 2) then
       call this%Xh%init(GLL, lx, lx)
    else
       call this%Xh%init(GLL, lx, lx, lx)
    end if

    call this%dm_Xh%init(msh, this%Xh)

    call this%gs_Xh%init(this%dm_Xh)

    call this%c_Xh%init(this%gs_Xh)

    ! Local scratch registry
    this%scratch = scratch_registry_t(this%dm_Xh, 10, 2)

    ! Assign a name
    call json_get_or_default(params, 'case.fluid.name', this%name, "fluid")

    !
    ! First section of fluid log
    !

    call neko_log%section('Fluid')
    write(log_buf, '(A, A)') 'Type       : ', trim(scheme)
    call neko_log%message(log_buf)
    write(log_buf, '(A, A)') 'Name       : ', trim(this%name)
    call neko_log%message(log_buf)

    !
    ! Material properties
    !
    call this%set_material_properties(params, user)

    !
    ! Turbulence modelling and variable material properties
    !
    if (params%valid_path('case.fluid.nut_field')) then
       call json_get(params, 'case.fluid.nut_field', this%nut_field_name)
       this%variable_material_properties = .true.
    else
       this%nut_field_name = ""
    end if

    ! Fill mu and rho field with the physical value

    call this%mu_field%init(this%dm_Xh, "mu")
    call this%rho_field%init(this%dm_Xh, "rho")
    call field_cfill(this%mu_field, this%mu, this%mu_field%size())
    call field_cfill(this%rho_field, this%rho, this%mu_field%size())

    ! Since mu, rho is a field, and the none-stress simulation fetches
    ! data from the host arrays, we need to mirror the constant
    ! material properties on the host
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call cfill(this%mu_field%x, this%mu, this%mu_field%size())
       call cfill(this%rho_field%x, this%rho, this%rho_field%size())
    end if


    ! Projection spaces
    call json_get_or_default(params, &
         'case.fluid.velocity_solver.projection_space_size', &
         this%vel_projection_dim, 20)
    call json_get_or_default(params, &
         'case.fluid.pressure_solver.projection_space_size', &
         this%pr_projection_dim, 20)
    call json_get_or_default(params, &
         'case.fluid.velocity_solver.projection_hold_steps', &
         this%vel_projection_activ_step, 5)
    call json_get_or_default(params, &
         'case.fluid.pressure_solver.projection_hold_steps', &
         this%pr_projection_activ_step, 5)


    call json_get_or_default(params, 'case.fluid.freeze', this%freeze, .false.)

    if (params%valid_path("case.fluid.flow_rate_force")) then
       this%forced_flow_rate = .true.
    end if


    if (lx .lt. 10) then
       write(log_buf, '(A, I1)') 'Poly order : ', lx-1
    else if (lx .ge. 10) then
       write(log_buf, '(A, I2)') 'Poly order : ', lx-1
    else
       write(log_buf, '(A, I3)') 'Poly order : ', lx-1
    end if
    call neko_log%message(log_buf)
    this%glb_n_points = int(this%msh%glb_nelv, i8)*int(this%Xh%lxyz, i8)
    this%glb_unique_points = int(glsum(this%c_Xh%mult, this%dm_Xh%size()), i8)

    write(log_buf, '(A, I0)') 'GLL points : ', this%glb_n_points
    call neko_log%message(log_buf)
    write(log_buf, '(A, I0)') 'Unique pts.: ', this%glb_unique_points
    call neko_log%message(log_buf)

    write(log_buf, '(A,ES13.6)') 'rho        :', this%rho
    call neko_log%message(log_buf)
    write(log_buf, '(A,ES13.6)') 'mu         :', this%mu
    call neko_log%message(log_buf)

    call json_get(params, 'case.numerics.dealias', logical_val)
    write(log_buf, '(A, L1)') 'Dealias    : ', logical_val
    call neko_log%message(log_buf)

    write(log_buf, '(A, L1)') 'LES        : ', this%variable_material_properties
    call neko_log%message(log_buf)

    call json_get_or_default(params, 'case.output_boundary', logical_val, &
         .false.)
    write(log_buf, '(A, L1)') 'Save bdry  : ', logical_val
    call neko_log%message(log_buf)

    !
    ! Setup right-hand side fields.
    !
    allocate(this%f_x)
    allocate(this%f_y)
    allocate(this%f_z)
    call this%f_x%init(this%dm_Xh, fld_name = "fluid_rhs_x")
    call this%f_y%init(this%dm_Xh, fld_name = "fluid_rhs_y")
    call this%f_z%init(this%dm_Xh, fld_name = "fluid_rhs_z")

    ! Initialize the source term
    call this%source_term%init(this%f_x, this%f_y, this%f_z, this%c_Xh, user)
    call this%source_term%add(params, 'case.fluid.source_terms')

    ! Initialize velocity solver
    if (kspv_init) then
       call neko_log%section("Velocity solver")
       call json_get_or_default(params, &
            'case.fluid.velocity_solver.max_iterations', &
            integer_val, KSP_MAX_ITER)
       call json_get(params, 'case.fluid.velocity_solver.type', string_val1)
       call json_get(params, 'case.fluid.velocity_solver.preconditioner', &
            string_val2)
       call json_get(params, 'case.fluid.velocity_solver.absolute_tolerance', &
            real_val)
       call json_get_or_default(params, &
            'case.fluid.velocity_solver.monitor', &
            logical_val, .false.)

       call neko_log%message('Type       : ('// trim(string_val1) // &
            ', ' // trim(string_val2) // ')')

       write(log_buf, '(A,ES13.6)') 'Abs tol    :', real_val
       call neko_log%message(log_buf)
       call this%solver_factory(this%ksp_vel, this%dm_Xh%size(), &
            string_val1, integer_val, real_val, logical_val)
       call this%precon_factory_(this%pc_vel, this%ksp_vel, &
            this%c_Xh, this%dm_Xh, this%gs_Xh, this%bcs_vel, string_val2)
       call neko_log%end_section()
    end if

    ! Strict convergence for the velocity solver
    call json_get_or_default(params, 'case.fluid.strict_convergence', &
         this%strict_convergence, .false.)

    ! Assign velocity fields
    call neko_field_registry%add_field(this%dm_Xh, 'u')
    call neko_field_registry%add_field(this%dm_Xh, 'v')
    call neko_field_registry%add_field(this%dm_Xh, 'w')
    this%u => neko_field_registry%get_field('u')
    this%v => neko_field_registry%get_field('v')
    this%w => neko_field_registry%get_field('w')

    !! Initialize time-lag fields
    call this%ulag%init(this%u, 2)
    call this%vlag%init(this%v, 2)
    call this%wlag%init(this%w, 2)

    ! Initiate gradient jump penalty
    call json_get_or_default(params, &
         'case.fluid.gradient_jump_penalty.enabled',&
         this%if_gradient_jump_penalty, .false.)

    if (this%if_gradient_jump_penalty .eqv. .true.) then
       if ((this%dm_Xh%xh%lx - 1) .eq. 1) then
          call json_get_or_default(params, &
               'case.fluid.gradient_jump_penalty.tau',&
               GJP_param_a, 0.02_rp)
          GJP_param_b = 0.0_rp
       else
          call json_get_or_default(params, &
               'case.fluid.gradient_jump_penalty.scaling_factor',&
               GJP_param_a, 0.8_rp)
          call json_get_or_default(params, &
               'case.fluid.gradient_jump_penalty.scaling_exponent',&
               GJP_param_b, 4.0_rp)
       end if
       call this%gradient_jump_penalty_u%init(params, this%dm_Xh, this%c_Xh, &
            GJP_param_a, GJP_param_b)
       call this%gradient_jump_penalty_v%init(params, this%dm_Xh, this%c_Xh, &
            GJP_param_a, GJP_param_b)
       call this%gradient_jump_penalty_w%init(params, this%dm_Xh, this%c_Xh, &
            GJP_param_a, GJP_param_b)
    end if

    call neko_field_registry%add_field(this%dm_Xh, 'u_e')
    call neko_field_registry%add_field(this%dm_Xh, 'v_e')
    call neko_field_registry%add_field(this%dm_Xh, 'w_e')
    this%u_e => neko_field_registry%get_field('u_e')
    this%v_e => neko_field_registry%get_field('v_e')
    this%w_e => neko_field_registry%get_field('w_e')

    call neko_log%end_section()

  end subroutine fluid_scheme_init_base

  subroutine fluid_scheme_free(this)
    class(fluid_scheme_incompressible_t), intent(inout) :: this

    !
    ! Free everything related to field_dirichlet BCs
    !

    call this%Xh%free()

    if (allocated(this%ksp_vel)) then
       call this%ksp_vel%free()
       deallocate(this%ksp_vel)
    end if

    if (allocated(this%ksp_prs)) then
       call this%ksp_prs%free()
       deallocate(this%ksp_prs)
    end if

    if (allocated(this%pc_vel)) then
       call precon_destroy(this%pc_vel)
       deallocate(this%pc_vel)
    end if

    if (allocated(this%pc_prs)) then
       call precon_destroy(this%pc_prs)
       deallocate(this%pc_prs)
    end if

    call this%source_term%free()

    call this%gs_Xh%free()

    call this%c_Xh%free()

    call this%scratch%free()

    nullify(this%u)
    nullify(this%v)
    nullify(this%w)
    nullify(this%p)

    if (this%variable_material_properties) then
       nullify(this%u_e)
       nullify(this%v_e)
       nullify(this%w_e)
    end if

    call this%ulag%free()
    call this%vlag%free()
    call this%wlag%free()


    if (associated(this%f_x)) then
       call this%f_x%free()
    end if

    if (associated(this%f_y)) then
       call this%f_y%free()
    end if

    if (associated(this%f_z)) then
       call this%f_z%free()
    end if

    nullify(this%f_x)
    nullify(this%f_y)
    nullify(this%f_z)

    call this%rho_field%free()
    call this%mu_field%free()

    ! Free gradient jump penalty
    if (this%if_gradient_jump_penalty .eqv. .true.) then
       call this%gradient_jump_penalty_u%free()
       call this%gradient_jump_penalty_v%free()
       call this%gradient_jump_penalty_w%free()
    end if

  end subroutine fluid_scheme_free

  !> Validate that all fields, solvers etc necessary for
  !! performing time-stepping are defined
  subroutine fluid_scheme_validate(this)
    class(fluid_scheme_incompressible_t), target, intent(inout) :: this
    ! Variables for retrieving json parameters
    logical :: logical_val

    if ( (.not. associated(this%u)) .or. &
         (.not. associated(this%v)) .or. &
         (.not. associated(this%w)) .or. &
         (.not. associated(this%p))) then
       call neko_error('Fields are not registered')
    end if

    if ( (.not. allocated(this%u%x)) .or. &
         (.not. allocated(this%v%x)) .or. &
         (.not. allocated(this%w%x)) .or. &
         (.not. allocated(this%p%x))) then
       call neko_error('Fields are not allocated')
    end if

    if (.not. allocated(this%ksp_vel)) then
       call neko_error('No Krylov solver for velocity defined')
    end if

    if (.not. allocated(this%ksp_prs)) then
       call neko_error('No Krylov solver for pressure defined')
    end if

    !
    ! Setup checkpoint structure (if everything is fine)
    !
    call this%chkp%init(this%u, this%v, this%w, this%p)

  end subroutine fluid_scheme_validate

  !> Apply all boundary conditions defined for velocity
  !! Here we perform additional gs operations to take care of
  !! shared points between elements that have different BCs, as done in Nek5000.
  !! @todo Why can't we call the interface here?
  subroutine fluid_scheme_bc_apply_vel(this, t, tstep, strong)
    class(fluid_scheme_incompressible_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    logical, intent(in) :: strong

    call this%bcs_vel%apply_vector(&
         this%u%x, this%v%x, this%w%x, this%dm_Xh%size(), t, tstep, strong)
    call this%gs_Xh%op(this%u, GS_OP_MIN, glb_cmd_event)
    call this%gs_Xh%op(this%v, GS_OP_MIN, glb_cmd_event)
    call this%gs_Xh%op(this%w, GS_OP_MIN, glb_cmd_event)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_stream_wait_event(glb_cmd_queue, glb_cmd_event, 0)
    end if

    call this%bcs_vel%apply_vector(&
         this%u%x, this%v%x, this%w%x, this%dm_Xh%size(), t, tstep, strong)
    call this%gs_Xh%op(this%u, GS_OP_MAX, glb_cmd_event)
    call this%gs_Xh%op(this%v, GS_OP_MAX, glb_cmd_event)
    call this%gs_Xh%op(this%w, GS_OP_MAX, glb_cmd_event)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_stream_wait_event(glb_cmd_queue, glb_cmd_event, 0)
    end if

  end subroutine fluid_scheme_bc_apply_vel

  !> Apply all boundary conditions defined for pressure
  !! @todo Why can't we call the interface here?
  subroutine fluid_scheme_bc_apply_prs(this, t, tstep)
    class(fluid_scheme_incompressible_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    call this%bcs_prs%apply(this%p, t, tstep)
    call this%gs_Xh%op(this%p,GS_OP_MIN, glb_cmd_event)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_stream_wait_event(glb_cmd_queue, glb_cmd_event, 0)
    end if

    call this%bcs_prs%apply(this%p, t, tstep)
    call this%gs_Xh%op(this%p,GS_OP_MAX, glb_cmd_event)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_stream_wait_event(glb_cmd_queue, glb_cmd_event, 0)
    end if


  end subroutine fluid_scheme_bc_apply_prs

  !> Initialize a linear solver
  !! @note Currently only supporting Krylov solvers
  subroutine fluid_scheme_solver_factory(ksp, n, solver, &
       max_iter, abstol, monitor)
    class(ksp_t), allocatable, target, intent(inout) :: ksp
    integer, intent(in), value :: n
    character(len=*), intent(in) :: solver
    integer, intent(in) :: max_iter
    real(kind=rp), intent(in) :: abstol
    logical, intent(in) :: monitor

    call krylov_solver_factory(ksp, n, solver, max_iter, abstol, &
         monitor = monitor)

  end subroutine fluid_scheme_solver_factory

  !> Initialize a Krylov preconditioner
  subroutine fluid_scheme_precon_factory(this, pc, ksp, coef, dof, gs, bclst, &
       pctype)
    class(fluid_scheme_incompressible_t), intent(inout) :: this
    class(pc_t), allocatable, target, intent(inout) :: pc
    class(ksp_t), target, intent(inout) :: ksp
    type(coef_t), target, intent(in) :: coef
    type(dofmap_t), target, intent(in) :: dof
    type(gs_t), target, intent(inout) :: gs
    type(bc_list_t), target, intent(inout) :: bclst
    character(len=*) :: pctype

    call precon_factory(pc, pctype)

    select type (pcp => pc)
    type is (jacobi_t)
       call pcp%init(coef, dof, gs)
    type is (sx_jacobi_t)
       call pcp%init(coef, dof, gs)
    type is (device_jacobi_t)
       call pcp%init(coef, dof, gs)
    type is (hsmg_t)
       if (len_trim(pctype) .gt. 4) then
          if (index(pctype, '+') .eq. 5) then
             call pcp%init(dof%msh, dof%Xh, coef, dof, gs, bclst, &
                  trim(pctype(6:)))
          else
             call neko_error('Unknown coarse grid solver')
          end if
       else
          call pcp%init(dof%msh, dof%Xh, coef, dof, gs, bclst)
       end if
    type is (phmg_t)
       call pcp%init(dof%msh, dof%Xh, coef, dof, gs, bclst)
    end select

    call ksp%set_pc(pc)

  end subroutine fluid_scheme_precon_factory

  !> Compute CFL
  function fluid_compute_cfl(this, dt) result(c)
    class(fluid_scheme_incompressible_t), intent(in) :: this
    real(kind=rp), intent(in) :: dt
    real(kind=rp) :: c

    c = cfl(dt, this%u%x, this%v%x, this%w%x, &
         this%Xh, this%c_Xh, this%msh%nelv, this%msh%gdim)

  end function fluid_compute_cfl


  !> Update the values of `mu_field` if necessary.
  subroutine fluid_scheme_update_material_properties(this)
    class(fluid_scheme_incompressible_t), intent(inout) :: this
    type(field_t), pointer :: nut
    integer :: n

    this%mu_field = this%mu
    if (this%variable_material_properties) then
       nut => neko_field_registry%get_field(this%nut_field_name)
       n = nut%size()
       call field_add2s2(this%mu_field, nut, this%rho, n)
    end if
  end subroutine fluid_scheme_update_material_properties

  !> Sets rho and mu
  !! @param params The case paramter file.
  !! @param user The user interface.
  subroutine fluid_scheme_set_material_properties(this, params, user)
    class(fluid_scheme_incompressible_t), intent(inout) :: this
    type(json_file), intent(inout) :: params
    type(user_t), target, intent(in) :: user
    character(len=LOG_SIZE) :: log_buf
    ! A local pointer that is needed to make Intel happy
    procedure(user_material_properties), pointer :: dummy_mp_ptr
    logical :: nondimensional
    real(kind=rp) :: dummy_lambda, dummy_cp

    dummy_mp_ptr => dummy_user_material_properties

    if (.not. associated(user%material_properties, dummy_mp_ptr)) then

       write(log_buf, '(A)') "Material properties must be set in the user&
       & file!"
       call neko_log%message(log_buf)
       call user%material_properties(0.0_rp, 0, this%rho, this%mu, &
            dummy_cp, dummy_lambda, params)
    else
       ! Incorrect user input
       if (params%valid_path('case.fluid.Re') .and. &
            (params%valid_path('case.fluid.mu') .or. &
            params%valid_path('case.fluid.rho'))) then
          call neko_error("To set the material properties for the fluid," // &
               " either provide Re OR mu and rho in the case file.")

          ! Non-dimensional case
       else if (params%valid_path('case.fluid.Re')) then

          write(log_buf, '(A)') 'Non-dimensional fluid material properties &
          & input.'
          call neko_log%message(log_buf, lvl = NEKO_LOG_VERBOSE)
          write(log_buf, '(A)') 'Density will be set to 1, dynamic viscosity to&
          & 1/Re.'
          call neko_log%message(log_buf, lvl = NEKO_LOG_VERBOSE)

          ! Read Re into mu for further manipulation.
          call json_get(params, 'case.fluid.Re', this%mu)
          write(log_buf, '(A)') 'Read non-dimensional material properties'
          call neko_log%message(log_buf)
          write(log_buf, '(A,ES13.6)') 'Re         :', this%mu
          call neko_log%message(log_buf)

          ! Set rho to 1 since the setup is non-dimensional.
          this%rho = 1.0_rp
          ! Invert the Re to get viscosity.
          this%mu = 1.0_rp/this%mu
          ! Dimensional case
       else
          call json_get(params, 'case.fluid.mu', this%mu)
          call json_get(params, 'case.fluid.rho', this%rho)
       end if

    end if
  end subroutine fluid_scheme_set_material_properties

end module fluid_scheme_incompressible
