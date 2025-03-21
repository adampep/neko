! Copyright (c) 2022-2024, The Neko Authors
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
!> Contains the scalar_scheme_t type.

module scalar_scheme
  use gather_scatter, only : gs_t
  use checkpoint, only : chkp_t
  use num_types, only: rp
  use field, only : field_t
  use field_list, only: field_list_t
  use space, only : space_t
  use dofmap, only : dofmap_t
  use krylov, only : ksp_t, krylov_solver_factory, KSP_MAX_ITER
  use coefs, only : coef_t
  use dirichlet, only : dirichlet_t
  use neumann, only : neumann_t
  use jacobi, only : jacobi_t
  use device_jacobi, only : device_jacobi_t
  use sx_jacobi, only : sx_jacobi_t
  use hsmg, only : hsmg_t
  use bc, only : bc_t
  use bc_list, only : bc_list_t
  use precon, only : pc_t, precon_factory, precon_destroy
  use field_dirichlet, only: field_dirichlet_t, field_dirichlet_update
  use mesh, only : mesh_t, NEKO_MSH_MAX_ZLBLS, NEKO_MSH_MAX_ZLBL_LEN
  use facet_zone, only : facet_zone_t
  use time_scheme_controller, only : time_scheme_controller_t
  use logger, only : neko_log, LOG_SIZE, NEKO_LOG_VERBOSE
  use field_registry, only : neko_field_registry
  use usr_scalar, only : usr_scalar_t, usr_scalar_bc_eval
  use json_utils, only : json_get, json_get_or_default, json_extract_item
  use json_module, only : json_file
  use user_intf, only : user_t, dummy_user_material_properties, &
       user_material_properties
  use utils, only : neko_error
  use comm, only: NEKO_COMM, MPI_INTEGER, MPI_SUM
  use scalar_source_term, only : scalar_source_term_t
  use field_series, only : field_series_t
  use math, only : cfill, add2s2
  use field_math, only : field_add2s2
  use device_math, only : device_cfill, device_add2s2
  use neko_config, only : NEKO_BCKND_DEVICE
  use field_series, only : field_series_t
  use time_step_controller, only : time_step_controller_t
  use gradient_jump_penalty, only : gradient_jump_penalty_t
  implicit none

  !> Base type for a scalar advection-diffusion solver.
  type, abstract :: scalar_scheme_t
     !> A name that can be used to distinguish this solver in e.g. user routines
     character(len=:), allocatable :: name
     !> x-component of Velocity
     type(field_t), pointer :: u
     !> y-component of Velocity
     type(field_t), pointer :: v
     !> z-component of Velocity
     type(field_t), pointer :: w
     !> The scalar.
     type(field_t), pointer :: s
     !> Lag arrays, i.e. solutions at previous timesteps.
     type(field_series_t) :: slag
     !> Function space \f$ X_h \f$.
     type(space_t), pointer :: Xh
     !> Dofmap associated with \f$ X_h \f$.
     type(dofmap_t), pointer :: dm_Xh
     !> Gather-scatter associated with \f$ X_h \f$.
     type(gs_t), pointer :: gs_Xh
     !> Coefficients associated with \f$ X_h \f$.
     type(coef_t), pointer :: c_Xh
     !> Right-hand side.
     type(field_t), pointer :: f_Xh => null()
     !> The source term for equation.
     type(scalar_source_term_t) :: source_term
     !> Krylov solver.
     class(ksp_t), allocatable :: ksp
     !> Max iterations in the Krylov solver.
     integer :: ksp_maxiter
     !> Projection space size.
     integer :: projection_dim
     !< Steps to activate projection for ksp
     integer :: projection_activ_step
     !> Preconditioner.
     class(pc_t), allocatable :: pc
     !> List of boundary conditions, including the user one.
     type(bc_list_t) :: bcs
     !> Case paramters.
     type(json_file), pointer :: params
     !> Mesh.
     type(mesh_t), pointer :: msh => null()
     !> Checkpoint for restarts.
     type(chkp_t) :: chkp
     !> Thermal diffusivity.
     real(kind=rp) :: lambda
     !> The variable lambda field
     type(field_t) :: lambda_field
     !> The turbulent kinematic viscosity field name
     character(len=:), allocatable :: nut_field_name
     !> Density.
     real(kind=rp) :: rho
     !> Specific heat capacity.
     real(kind=rp) :: cp
     !> Turbulent Prandtl number.
     real(kind=rp) :: pr_turb
     !> Is lambda varying in time? Currently only due to LES models.
     logical :: variable_material_properties = .false.
     !> Gradient jump panelty
     logical :: if_gradient_jump_penalty
     type(gradient_jump_penalty_t) :: gradient_jump_penalty
   contains
     !> Constructor for the base type.
     procedure, pass(this) :: scheme_init => scalar_scheme_init
     !> Destructor for the base type.
     procedure, pass(this) :: scheme_free => scalar_scheme_free
     !> Validate successful initialization.
     procedure, pass(this) :: validate => scalar_scheme_validate
     !> Set lambda and cp
     procedure, pass(this) :: set_material_properties => &
          scalar_scheme_set_material_properties
     !> Update variable material properties
     procedure, pass(this) :: update_material_properties => &
          scalar_scheme_update_material_properties
     !> Constructor.
     procedure(scalar_scheme_init_intrf), pass(this), deferred :: init
     !> Destructor.
     procedure(scalar_scheme_free_intrf), pass(this), deferred :: free
     !> Solve for the current timestep.
     procedure(scalar_scheme_step_intrf), pass(this), deferred :: step
     !> Restart from a checkpoint.
     procedure(scalar_scheme_restart_intrf), pass(this), deferred :: restart
  end type scalar_scheme_t

  !> Abstract interface to initialize a scalar formulation
  abstract interface
     subroutine scalar_scheme_init_intrf(this, msh, coef, gs, params, &
          numerics_params, user, ulag, vlag, wlag, time_scheme, rho)
       import scalar_scheme_t
       import json_file
       import coef_t
       import gs_t
       import mesh_t
       import user_t
       import field_series_t
       import time_scheme_controller_t
       import rp
       class(scalar_scheme_t), target, intent(inout) :: this
       type(mesh_t), target, intent(in) :: msh
       type(coef_t), target, intent(in) :: coef
       type(gs_t), target, intent(inout) :: gs
       type(json_file), target, intent(inout) :: params
       type(json_file), target, intent(inout) :: numerics_params
       type(user_t), target, intent(in) :: user
       type(field_series_t), target, intent(in) :: ulag, vlag, wlag
       type(time_scheme_controller_t), target, intent(in) :: time_scheme
       real(kind=rp), intent(in) :: rho
     end subroutine scalar_scheme_init_intrf
  end interface

  !> Abstract interface to restart a scalar formulation
  abstract interface
     subroutine scalar_scheme_restart_intrf(this, dtlag, tlag)
       import scalar_scheme_t
       import chkp_t
       import rp
       class(scalar_scheme_t), target, intent(inout) :: this
       real(kind=rp) :: dtlag(10), tlag(10)
     end subroutine scalar_scheme_restart_intrf
  end interface

  !> Abstract interface to dealocate a scalar formulation
  abstract interface
     subroutine scalar_scheme_free_intrf(this)
       import scalar_scheme_t
       class(scalar_scheme_t), intent(inout) :: this
     end subroutine scalar_scheme_free_intrf
  end interface

  !> Abstract interface to compute a time-step
  abstract interface
     subroutine scalar_scheme_step_intrf(this, t, tstep, dt, ext_bdf, &
          dt_controller)
       import scalar_scheme_t
       import time_scheme_controller_t
       import time_step_controller_t
       import rp
       class(scalar_scheme_t), intent(inout) :: this
       real(kind=rp), intent(in) :: t
       integer, intent(in) :: tstep
       real(kind=rp), intent(in) :: dt
       type(time_scheme_controller_t), intent(in) :: ext_bdf
       type(time_step_controller_t), intent(in) :: dt_controller
     end subroutine scalar_scheme_step_intrf
  end interface

contains

  !> Initialize all related components of the current scheme
  !! @param msh The mesh.
  !! @param c_Xh The coefficients.
  !! @param gs_Xh The gather-scatter.
  !! @param params The parameter dictionary in json.
  !! @param scheme The name of the scalar scheme.
  !! @param user Type with user-defined procedures.
  !! @param rho The density of the fluid.
  subroutine scalar_scheme_init(this, msh, c_Xh, gs_Xh, params, scheme, user, &
       rho)
    class(scalar_scheme_t), target, intent(inout) :: this
    type(mesh_t), target, intent(in) :: msh
    type(coef_t), target, intent(in) :: c_Xh
    type(gs_t), target, intent(inout) :: gs_Xh
    type(json_file), target, intent(inout) :: params
    character(len=*), intent(in) :: scheme
    type(user_t), target, intent(in) :: user
    real(kind=rp), intent(in) :: rho
    ! IO buffer for log output
    character(len=LOG_SIZE) :: log_buf
    ! Variables for retrieving json parameters
    logical :: logical_val
    real(kind=rp) :: real_val, solver_abstol
    integer :: integer_val, ierr
    character(len=:), allocatable :: solver_type, solver_precon, field_name
    real(kind=rp) :: GJP_param_a, GJP_param_b

    this%u => neko_field_registry%get_field('u')
    this%v => neko_field_registry%get_field('v')
    this%w => neko_field_registry%get_field('w')

    ! Assign a name
    call json_get_or_default(params, 'name', this%name, 'scalar')

    call neko_log%section('Scalar')
    call json_get(params, 'solver.type', solver_type)
    call json_get(params, 'solver.preconditioner', &
         solver_precon)
    call json_get(params, 'solver.absolute_tolerance', &
         solver_abstol)

    call json_get_or_default(params, &
         'solver.projection_space_size', &
         this%projection_dim, 20)
    call json_get_or_default(params, &
         'solver.projection_hold_steps', &
         this%projection_activ_step, 5)


    write(log_buf, '(A, A)') 'Type       : ', trim(scheme)
    call neko_log%message(log_buf)
    write(log_buf, '(A, A)') 'Name       : ', trim(this%name)
    call neko_log%message(log_buf)
    call neko_log%message('Ksp scalar : ('// trim(solver_type) // &
         ', ' // trim(solver_precon) // ')')
    write(log_buf, '(A,ES13.6)') ' `-abs tol :', solver_abstol
    call neko_log%message(log_buf)

    this%Xh => this%u%Xh
    this%dm_Xh => this%u%dof
    this%params => params
    this%msh => msh

    call json_get_or_default(params, 'field_name', field_name, 's')

    if (.not. neko_field_registry%field_exists(field_name)) then
       call neko_field_registry%add_field(this%dm_Xh, field_name)
    end if

    this%s => neko_field_registry%get_field(field_name)

    call this%slag%init(this%s, 2)

    this%gs_Xh => gs_Xh
    this%c_Xh => c_Xh

    !
    ! Material properties
    !
    this%rho = rho
    call this%set_material_properties(params, user)

    write(log_buf, '(A,ES13.6)') 'rho        :', this%rho
    call neko_log%message(log_buf)
    write(log_buf, '(A,ES13.6)') 'lambda     :', this%lambda
    call neko_log%message(log_buf)
    write(log_buf, '(A,ES13.6)') 'cp         :', this%cp
    call neko_log%message(log_buf)

    !
    ! Turbulence modelling and variable material properties
    !
    if (params%valid_path('nut_field')) then
       call json_get(params, 'Pr_t', this%pr_turb)
       call json_get(params, 'nut_field', this%nut_field_name)
       this%variable_material_properties = .true.
    else
       this%nut_field_name = ""
    end if

    write(log_buf, '(A,L1)') 'LES        : ', this%variable_material_properties
    call neko_log%message(log_buf)

    ! Fill lambda field with the physical value
    call this%lambda_field%init(this%dm_Xh, "lambda")
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cfill(this%lambda_field%x_d, this%lambda, &
            this%lambda_field%size())
    else
       call cfill(this%lambda_field%x, this%lambda, this%lambda_field%size())
    end if

    !
    ! Setup right-hand side field.
    !
    allocate(this%f_Xh)
    call this%f_Xh%init(this%dm_Xh, fld_name = "scalar_rhs")

    ! Initialize the source term
    call this%source_term%init(this%f_Xh, this%c_Xh, user)
    call this%source_term%add(params, 'source_terms')

    ! todo parameter file ksp tol should be added
    call json_get_or_default(params, &
         'solver.max_iterations', &
         integer_val, KSP_MAX_ITER)
    call json_get_or_default(params, &
         'solver.monitor', &
         logical_val, .false.)
    call scalar_scheme_solver_factory(this%ksp, this%dm_Xh%size(), &
         solver_type, integer_val, solver_abstol, logical_val)
    call scalar_scheme_precon_factory(this%pc, this%ksp, &
         this%c_Xh, this%dm_Xh, this%gs_Xh, this%bcs, solver_precon)

    ! Initiate gradient jump penalty
    call json_get_or_default(params, &
         'gradient_jump_penalty.enabled',&
         this%if_gradient_jump_penalty, .false.)

    if (this%if_gradient_jump_penalty .eqv. .true.) then
       if ((this%dm_Xh%xh%lx - 1) .eq. 1) then
          call json_get_or_default(params, &
               'gradient_jump_penalty.tau',&
               GJP_param_a, 0.02_rp)
          GJP_param_b = 0.0_rp
       else
          call json_get_or_default(params, &
               'gradient_jump_penalty.scaling_factor',&
               GJP_param_a, 0.8_rp)
          call json_get_or_default(params, &
               'gradient_jump_penalty.scaling_exponent',&
               GJP_param_b, 4.0_rp)
       end if
       call this%gradient_jump_penalty%init(params, this%dm_Xh, this%c_Xh, &
            GJP_param_a, GJP_param_b)
    end if

    call neko_log%end_section()

  end subroutine scalar_scheme_init


  !> Deallocate a scalar formulation
  subroutine scalar_scheme_free(this)
    class(scalar_scheme_t), intent(inout) :: this

    nullify(this%Xh)
    nullify(this%dm_Xh)
    nullify(this%gs_Xh)
    nullify(this%c_Xh)
    nullify(this%params)

    if (allocated(this%ksp)) then
       call this%ksp%free()
       deallocate(this%ksp)
    end if

    if (allocated(this%pc)) then
       call precon_destroy(this%pc)
       deallocate(this%pc)
    end if

    call this%source_term%free()

    call this%bcs%free()

    call this%lambda_field%free()
    call this%slag%free()

    ! Free gradient jump penalty
    if (this%if_gradient_jump_penalty .eqv. .true.) then
       call this%gradient_jump_penalty%free()
    end if

  end subroutine scalar_scheme_free

  !> Validate that all fields, solvers etc necessary for
  !! performing time-stepping are defined
  subroutine scalar_scheme_validate(this)
    class(scalar_scheme_t), target, intent(inout) :: this

    if ( (.not. allocated(this%u%x)) .or. &
         (.not. allocated(this%v%x)) .or. &
         (.not. allocated(this%w%x)) .or. &
         (.not. allocated(this%s%x))) then
       call neko_error('Fields are not allocated')
    end if

    if (.not. allocated(this%ksp)) then
       call neko_error('No Krylov solver for velocity defined')
    end if

    if (.not. associated(this%Xh)) then
       call neko_error('No function space defined')
    end if

    if (.not. associated(this%dm_Xh)) then
       call neko_error('No dofmap defined')
    end if

    if (.not. associated(this%c_Xh)) then
       call neko_error('No coefficients defined')
    end if

    if (.not. associated(this%f_Xh)) then
       call neko_error('No rhs allocated')
    end if

    if (.not. associated(this%params)) then
       call neko_error('No parameters defined')
    end if

    !
    ! Setup checkpoint structure (if everything is fine)
    !
!    @todo no io for now
!    call this%chkp%init(this%u, this%v, this%w, this%p)

  end subroutine scalar_scheme_validate

  !> Initialize a linear solver
  !! @note Currently only supporting Krylov solvers
  subroutine scalar_scheme_solver_factory(ksp, n, solver, max_iter, &
       abstol, monitor)
    class(ksp_t), allocatable, target, intent(inout) :: ksp
    integer, intent(in), value :: n
    integer, intent(in) :: max_iter
    character(len=*), intent(in) :: solver
    real(kind=rp) :: abstol
    logical, intent(in) :: monitor

    call krylov_solver_factory(ksp, n, solver, max_iter, &
         abstol, monitor = monitor)

  end subroutine scalar_scheme_solver_factory

  !> Initialize a Krylov preconditioner
  subroutine scalar_scheme_precon_factory(pc, ksp, coef, dof, gs, bclst, &
       pctype)
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
    end select

    call ksp%set_pc(pc)

  end subroutine scalar_scheme_precon_factory

  !> Update the values of `lambda_field` if necessary.
  subroutine scalar_scheme_update_material_properties(this)
    class(scalar_scheme_t), intent(inout) :: this
    type(field_t), pointer :: nut
    integer :: n
    ! Factor to transform nu_t to lambda_t
    real(kind=rp) :: lambda_factor

    lambda_factor = this%rho*this%cp/this%pr_turb
    this%lambda_field = this%lambda

    if (this%variable_material_properties) then
       nut => neko_field_registry%get_field(this%nut_field_name)
       n = nut%size()
       call field_add2s2(this%lambda_field, nut, lambda_factor, n)
    end if

  end subroutine scalar_scheme_update_material_properties

  !> Set lamdba and cp.
  !! @param params The case file configuration dictionary.
  !! @param user The user interface.
  subroutine scalar_scheme_set_material_properties(this, params, user)
    class(scalar_scheme_t), intent(inout) :: this
    type(json_file), intent(inout) :: params
    type(user_t), target, intent(in) :: user
    character(len=LOG_SIZE) :: log_buf
    ! A local pointer that is needed to make Intel happy
    procedure(user_material_properties), pointer :: dummy_mp_ptr
    real(kind=rp) :: dummy_mu, dummy_rho

    dummy_mp_ptr => dummy_user_material_properties

    if (.not. associated(user%material_properties, dummy_mp_ptr)) then

       write(log_buf, '(A)') "Material properties must be set in the user " // &
            "file!"
       call neko_log%message(log_buf)
       call user%material_properties(0.0_rp, 0, dummy_rho, dummy_mu, &
            this%cp, this%lambda, params)
    else
       if (params%valid_path('Pe') .and. &
            (params%valid_path('lambda') .or. &
            params%valid_path('cp'))) then
          call neko_error("To set the material properties for the scalar, " // &
               "either provide Pe OR lambda and cp in the case file.")
          ! Non-dimensional case
       else if (params%valid_path('Pe')) then
          write(log_buf, '(A)') 'Non-dimensional scalar material properties' //&
               ' input.'
          call neko_log%message(log_buf, lvl = NEKO_LOG_VERBOSE)
          write(log_buf, '(A)') 'Specific heat capacity will be set to 1,'
          call neko_log%message(log_buf, lvl = NEKO_LOG_VERBOSE)
          write(log_buf, '(A)') 'conductivity to 1/Pe. Assumes density is 1.'
          call neko_log%message(log_buf, lvl = NEKO_LOG_VERBOSE)

          ! Read Pe into lambda for further manipulation.
          call json_get(params, 'Pe', this%lambda)
          write(log_buf, '(A,ES13.6)') 'Pe         :', this%lambda
          call neko_log%message(log_buf)

          ! Set cp and rho to 1 since the setup is non-dimensional.
          this%cp = 1.0_rp
          this%rho = 1.0_rp
          ! Invert the Pe to get conductivity
          this%lambda = 1.0_rp/this%lambda
          ! Dimensional case
       else
          call json_get(params, 'lambda', this%lambda)
          call json_get(params, 'cp', this%cp)
       end if

    end if
  end subroutine scalar_scheme_set_material_properties

end module scalar_scheme
