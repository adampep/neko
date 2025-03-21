# User File {#user-file}

\tableofcontents

The user file is a fortran file where the user can implement their own functions
to extend the capabilities of the default Neko executable. The user file can be
used for setting advanced initial/boundary conditions, source terms, I/O
operations, and interactions with the Neko framework.

## Compiling and running

The user file is a regular Fortran `.f90` file that needs to be compiled with
`makeneko`, located in the `bin` folder of your neko installation. To compile a
user file `user.f90`, run:

```bash
makeneko user.f90
```

If everything goes well, you should observe the following output:

```bash
N E K O build tool, Version 0.7.99
(build: 2024-02-13 on x86_64-pc-linux-gnu using gnu)

Building user NEKO ... done!
```

Compiling your user file with `makeneko` will create a `neko` executable, which
you will need to execute with your case file as an argument. For example, if
your case file is called `user.case`:

```bash
./neko user.case
```

Or in parallel using MPI:

```bash
mpirun -n 8 ./neko user.case
```


## High-level structure

The current high-level structure of the user file is shown below.

```fortran
module user
  use neko
  implicit none

contains

  ! Register user defined functions here (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u


  end subroutine user_setup

end module user

```

The user file implements the `user` module. The `user` modules contains a
subroutine named `user_setup`, which we use to interface the internal procedures
defined in `src/common/user_intf.f90` with the subroutines that you will
implement in your user file. Each user subroutine should be implemented under
the `contains` statement, below `user_setup`.

@note The above code snippet is the most basic code structure for the user file.
Compiling it and running it would be equivalent to running the "vanilla" neko
executable `bin/neko` in your local neko installation folder.

## Default user functions

The following user functions, if defined in the user file, will always be
executed, regardless of what is set in the case file:

- [user_startup](@ref user-file_init-and-final): For early access to the case
  file, its manipulation or initializing simple user parameter variables.
- [user_init_modules](@ref user-file_init-and-final): For initializing user
  variables and objects
- [user_finalize_modules](@ref user-file_init-and-final): For finalizing, e.g
  freeing variables and terminating processes
- [user_check](@ref user-file_user-check): Executed at the end of every time
  step, for e.g. computing and/or outputting user defined quantities.
- [material_properties](@ref user-file_mat-prop): For computing and setting
  material properties such as `rho`, `mu`, `cp` and `lambda`.
- [user_mesh_setup](@ref user-file_user-mesh-setup): For applying a deformation
  to the mesh element nodes, before the simulation time loop.
- [scalar_user_bc](@ref user-file_scalar-bc): For applying boundary conditions
  to the scalar, on all zones that are not already specified with uniform
  dirichlet values e.g. `d=1`. For more information on the scalar, see the
  [relevant section of the case file](@ref case-file_scalar).

### Initializing and finalizing {#user-file_init-and-final}

Three subroutines `user_startup`, `user_init_modules` and
`user_finalize_modules` may be used to initialize/finalize any user defined
variables, external objects, or processes.  The `user_startup` routine is called
immediately after the case file is read in, meaning that no important objects
are initialized yet (e.g. the mesh, fluid, etc.). The `user_init_modules` and
`user_finalize_modules` are respectively executed right before/after the
simulation time loop.

In most cases, one can use `user_init_modules` routine and not the
`user_startup`. The latter is only necessary when you want to either manipulate
the case file programmatically before it is used in the simulation, or define
some variables that will be used in the constructors of some of the types.
An example of the latter is defining some material property constants that
can be used in the user `material_properties` routine, which is run by the
constructor of the `case%fluid` object.


```fortran
  ! Manipulate the case file and define simple user variables
  subroutine startup(params)
    type(json_file), intent(inout) :: params

    ! insert your initialization code here

  end subroutine startup

  ! Initialize user variables or external objects
  subroutine initialize(t, u, v, w, p, coef, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params

    ! insert your initialization code here

  end subroutine initialize

  ! Finalize user variables or external objects
  subroutine finalize(t, params)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: params

    ! insert your code here

  end subroutine initialize
```

In the example above, the subroutines `startup`, `initialize` and `finalize`
contain the actual implementations. They must also be interfaced to the internal
procedures `user_init_modules` and `user_finalize_modules` in `user_setup`:

```fortran

  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u

    u%user_startup => startup
    u%user_init_modules => initialize
    u%user_finalize_modules => finalize

  end subroutine user_setup

```

@note All three routines are independent of each
other. Using one does not require the use of the other.

### Computing at every time step {#user-file_user-check}

The subroutine `user_check` is executed at the end of every time step. It can be
used for computing and/or outputting your own variables/quantities at every time
step.
```fortran
  ! This is called at the end of every time step
  subroutine usercheck(t, tstep, u, v, w, p, coef, param)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: param

    ! insert code below

  end subroutine usercheck

```

In the example above, the subroutine `usercheck` contains the actual
implementation, and needs to be registered by adding:

```fortran
u%user_check => usercheck
```

to our `user_setup`.

### Setting material properties {#user-file_mat-prop}

`material_properties` allows for more complex computations and setting of
various material properties, such as `rho`, `mu` for the fluid and `cp`,
`lambda` for the scalar. The example below is taken from the
[rayleigh_benard_cylinder
example](https://github.com/ExtremeFLOW/neko/blob/564686b127ff75a362a06126c6b23e9b4e21879e/examples/rayleigh_benard_cylinder/rayleigh.f90#L22C1-L38C41).

```fortran

  subroutine set_material_properties(t, tstep, rho, mu, cp, lambda, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    real(kind=rp), intent(inout) :: rho, mu, cp, lambda
    type(json_file), intent(inout) :: params

    ! Re and Pr computed in `user_startup`

    mu = 1.0_rp / Re
    lambda = mu / Pr
    rho = 1.0_rp
    cp = 1.0_rp
  end subroutine set_material_properties

```

And of course not forgetting to register our function in `user_setup` by adding
the following line:

```fortran
u%material_properties => set_material_properties
```

### Runtime mesh deformation {#user-file_user-mesh-setup}

This user function allows for the modification of the mesh at runtime, by acting
on the element nodes of the mesh specified in the case file. This function is
only called once before the simulation time loop. The example below is taken
from the [tgv
example](https://github.com/ExtremeFLOW/neko/blob/a0613606360240e5059e65d6d98f4a57cf73e237/examples/tgv/tgv.f90#L27-L42).

```fortran
  ! Rescale mesh
  subroutine user_mesh_scale(msh)
    type(mesh_t), intent(inout) :: msh
    integer :: i, p, nvert
    real(kind=rp) :: d
    d = 4._rp

    ! original mesh has size 0..8 to be mapped onto -pi..pi
    ! will be updated later to a method giving back the vertices of the mesh
    nvert = size(msh%points)
    do i = 1, nvert
       msh%points(i)%x(1) = (msh%points(i)%x(1) - d) / d * pi
       msh%points(i)%x(2) = (msh%points(i)%x(2) - d) / d * pi
       msh%points(i)%x(3) = (msh%points(i)%x(3) - d) / d * pi
    end do

  end subroutine user_mesh_scale

```

The registering of the above function in `user_setup` should then be done as follows:

```fortran
    u%user_mesh_setup => user_mesh_scale
```

### Scalar boundary conditions {#user-file_scalar-bc}

This user function can be used to specify the scalar boundary values, on all
boundaries of type `user_pointwise`. See [relevant section of the case
file](@ref case-file_scalar). The example below sets the scalar boundary
condition values to be a linear function of the `z` coordinate (taken from the
[rayleigh_benard
example](https://github.com/ExtremeFLOW/neko/blob/aa72ad9bf34cbfbac0ee893c045639fdd095f80a/examples/rayleigh_benard_cylinder/rayleigh.f90#L41-L63)).

```fortran

  subroutine set_scalar_boundary_conditions(s, x, y, z, nx, ny, nz, ix, iy, iz, ie, t, tstep)
    real(kind=rp), intent(inout) :: s
    real(kind=rp), intent(in) :: x
    real(kind=rp), intent(in) :: y
    real(kind=rp), intent(in) :: z
    real(kind=rp), intent(in) :: nx
    real(kind=rp), intent(in) :: ny
    real(kind=rp), intent(in) :: nz
    integer, intent(in) :: ix
    integer, intent(in) :: iy
    integer, intent(in) :: iz
    integer, intent(in) :: ie
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    ! This will be used on all zones without labels
    s = 1.0_rp - z

  end subroutine set_scalar_boundary_conditions

```

This function will be called on all the points on the relevant boundaries. The
registering of the above function in `user_setup` should be done as follows:

```fortran
    u%scalar_user_bc => set_scalar_boundary_conditions
```

### User defined simulation components {#user-file_simcomps}

In addition to the case-specific user functions, the user can also define their
own simulation components. This can be done by writing a new type which extends
the \ref simulation_component_t type, and implementing the necessary functions
for the new type. The user can then specify the component in the list of
simulation components in the case file. The setting `is_user` should be set to
`true` in the JSON object for the new simulation component. The typename is used
to extract the settings for the simulation component from the JSON file.

```json
{
    "type": "my_comp",
    "is_user": true,
    // other settings
}
```
```fortran
  subroutine user_simcomp(params)
    type(json_file), intent(inout) :: params
    type(user_simcomp_t), allocatable :: my_simcomp
    type(json_file) :: simcomp_settings

    ! Allocate a simulation component
    allocate(my_simcomp)
    simcomp_settings = simulation_component_user_settings("my_comp", params)
    call neko_simcomps%add_user_simcomp(my_simcomp, simcomp_settings)

  end subroutine user_simcomp
```

In the example above, the subroutine `user_simcomp` contains the actual
implementation, and needs to be registered by adding:

```fortran
    u%init_user_simcomp => user_simcomp
```

A full example of a user-defined simulation component can be found in the
examples.

## Case-specific user functions

As explained in the [case file](case-file.md) page, certain components of the
simulation can be set to be user defined. These components and their associated
user functions are:

| Description                          | User function                                                   | JSON Object in the case file                                    |
| ------------------------------------ | --------------------------------------------------------------- | --------------------------------------------------------------- |
| Fluid initial condition              | [fluid_user_ic](@ref user-file_user-ic)                         | `case.fluid.initial_condition`                                  |
| Scalar initial condition             | [scalar_user_ic](@ref user-file_user-ic)                        | `case.scalar.initial_condition`                                 |
| Fluid inflow boundary condition      | [fluid_user_if](@ref user-file_fluid-user-if)                   | `case.fluid.inflow_condition`                                   |
| Scalar boundary conditions           | [scalar_user_bc](@ref user-file_scalar-bc)                      | (user function is always called)                                |
| Fluid source term                    | [fluid_user_f_vector or fluid_user_f](@ref user-file_user-f)    | `case.fluid.source_terms`                                       |
| Scalar source term                   | [scalar_user_f_vector or scalar_user_f](@ref user-file_user-f)  | `case.scalar.source_terms`                                      |
| Fluid and Scalar boundary conditions | [field_dirichlet_update](@ref user-file_field-dirichlet-update) | `case.fluid.boundary_types` and/or `case.scalar.boundary_types` |

Note that `scalar_user_bc` is included for completeness but is technically not case-specific.

### Fluid and Scalar initial conditions {#user-file_user-ic}

Enabling user defined initial conditions for the fluid and/or scalar is done by
setting the `initial_condition.type` to `"user"` in the relevant sections of the
case file, `case.fluid` and/or `case.scalar`.

```.json

"case": {
    "fluid": {
        "initial_condition": {
            "type": "user"
        }
    }
}
```

See the relevant sections on the [fluid](@ref case-file_fluid-ic) and
[scalar](@ref case-file_scalar) initial conditions in the
[case file page](@ref case-file) for more details.

The associated user functions for the fluid and/or scalar initial conditions can
then be added to the user file. An example for the fluid taken from the
[advecting cone example](https://github.com/ExtremeFLOW/neko/blob/aa72ad9bf34cbfbac0ee893c045639fdd095f80a/examples/advecting_cone/advecting_cone.f90#L48-L75),
is shown below.

```fortran

  !> Set the advecting velocity field.
  subroutine set_velocity(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params
    integer :: i, e, k, j
    real(kind=rp) :: x, y

    do i = 1, u%dof%size()
       x = u%dof%x(i,1,1,1)
       y = u%dof%y(i,1,1,1)

       ! Angular velocity is pi, giving a full rotation in 2 sec
       u%x(i,1,1,1) = -y*pi
       v%x(i,1,1,1) = x*pi
       w%x(i,1,1,1) = 0
    end do

     if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(u%x, u%x_d, u%dof%size(), &
                          HOST_TO_DEVICE, sync=.false.)
       call device_memcpy(v%x, v%x_d, v%dof%size(), &
                          HOST_TO_DEVICE, sync=.false.)
       call device_memcpy(w%x, w%x_d, w%dof%size(), &
                          HOST_TO_DEVICE, sync=.false.)
    end if

  end subroutine set_velocity

```

@note Notice the use of the `NEKO_BCKND_DEVICE` flag, which will be set to 1 if
running on GPUs, and the calls to `device_memcpy` to transfer data between the
host and the device. See [Running on GPUs](@ref user-file_tips_running-on-gpus) for
more information on how this works.

The same can be done for the scalar, with the example below also inspired from
the
[advecting cone example](https://github.com/ExtremeFLOW/neko/blob/aa72ad9bf34cbfbac0ee893c045639fdd095f80a/examples/advecting_cone/advecting_cone.f90#L14-L45):

```fortran

  !> User initial condition for the scalar
  subroutine set_s_ic(s, params)
    type(field_t), intent(inout) :: s
    type(json_file), intent(inout) :: params
    integer :: i, e, k, j
    real(kind=rp) :: cone_radius, mux, muy, x, y, r, theta

    ! Center of the cone
    mux = 1
    muy = 0

    cone_radius = 0.5

    do i = 1, s%dof%size()
       x = s%dof%x(i,1,1,1) - mux
       y = s%dof%y(i,1,1,1) - muy

       r = sqrt(x**2 + y**2)
       theta = atan2(y, x)

       ! Check if the point is inside the cone's base
       if (r > cone_radius) then
         s%x(i,1,1,1) = 0.0
       else
         s%x(i,1,1,1) = 1.0 - r / cone_radius
       endif
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(s%x, s%x_d, s%dof%size(), &
                          HOST_TO_DEVICE, sync=.false.)
    end if

  end subroutine set_s_ic

```

@note Notice the use of the `NEKO_BCKND_DEVICE` flag, which will be set to 1 if
running on GPUs, and the calls to `device_memcpy` to transfer data between the
host and the device. See [Running on GPUs](@ref user-file_tips_running-on-gpus) for
more information on how this works.

We should also add of the following lines in `user_setup`, registering our user
functions `set_velocity` and `set_s_ic` to be used as the fluid and scalar
initial conditions:

```fortran
u%fluid_user_ic => set_velocity
u%scalar_user_ic => set_s_ic
```

### Fluid pointwise velocity boundary values {#user-file_fluid-user-if}

This routine will be used to compute the velocity on all boundaries of type
`user_velocity_pointwise`.
See the [the relevant section](@ref case-file_fluid-boundary-conditions) in the
[case file page](@ref case-file) for more details. The associated user
function for the fluid velocity can then be added to the user file.
An example inspired from the
[lid-driven cavity example](https://github.com/ExtremeFLOW/neko/blob/aa72ad9bf34cbfbac0ee893c045639fdd095f80a/examples/lid/lid.f90#L29-L53)
is shown below.

```fortran
 ! user-defined boundary condition
  subroutine user_bc(u, v, w, x, y, z, nx, ny, nz, ix, iy, iz, ie, t, tstep)
    real(kind=rp), intent(inout) :: u
    real(kind=rp), intent(inout) :: v
    real(kind=rp), intent(inout) :: w
    real(kind=rp), intent(in) :: x
    real(kind=rp), intent(in) :: y
    real(kind=rp), intent(in) :: z
    real(kind=rp), intent(in) :: nx
    real(kind=rp), intent(in) :: ny
    real(kind=rp), intent(in) :: nz
    integer, intent(in) :: ix
    integer, intent(in) :: iy
    integer, intent(in) :: iz
    integer, intent(in) :: ie
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    real(kind=rp) lsmoothing
    lsmoothing = 0.05_rp    ! length scale of smoothing at the edges

    u = step( x/lsmoothing ) * step( (1._rp-x)/lsmoothing )
    v = 0._rp
    w = 0._rp

  end subroutine user_bc
```

We should also add of the following line in `user_setup`, registering our user
function `user_bc` to be used as the fluid inflow conditions:

```fortran
u%fluid_user_if => user_bc
```

### Fluid and scalar source terms {#user-file_user-f}

Enabling user defined source terms for the fluid and/or scalar is done by adding
JSON Objects to the `case.fluid.source_terms` and/or `case.scalar.source_terms`
lists.

```.json

"case": {
    "fluid": {
        "source_terms":
        [
            {
                "type": "user_vector"
            }
        ]
    }
}
```

See the relevant sections on the [fluid](@ref case-file_fluid-source-term) and
[scalar](@ref case-file_scalar) source terms in the [case file page](@ref case-file) for
more details.

@attention There are two variants of the source term user functions: `_user_f`
and `_user_f_vector`. The former is called when setting `"user_pointwise"` as
the source term type, while the latter requires the use of the `"user_vector"`
keyword in the case file. The pointwise variant, `fluid_user_f` is not supported
on GPUs. In general, `fluid_user_f_vector` is the prefered variant, and is the
one which will be use in our examples below. The same applies for the scalar
source term user functions.

The associated user functions for the fluid and/or scalar source terms can then
be added to the user file. An example for the fluid, taken from the
[rayleigh_benard_cylinder example](https://github.com/ExtremeFLOW/neko/blob/49925b7a04a638259db3b1ddd54349ca57f5d207/examples/rayleigh_benard_cylinder/rayleigh.f90#L101C1-L121C44),
is shown below.

```fortran
  ! Sets the z-component of the fluid forcing term = scalar
  subroutine set_bousinesq_forcing_term(f, t)
    class(fluid_user_source_term_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t

    ! Retrieve u,v,w,s fields from the field registry
    type(field_t), pointer :: u, v, w, s
    u => neko_field_registry%get_field('u')
    v => neko_field_registry%get_field('v')
    w => neko_field_registry%get_field('w')
    s => neko_field_registry%get_field('s')

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_rzero(f%u_d,f%dm%size())
       call device_rzero(f%v_d,f%dm%size())
       call device_copy(f%w_d,s%x_d,f%dm%size())
    else
       call rzero(f%u,f%dm%size())
       call rzero(f%v,f%dm%size())
       call copy(f%w,s%x,f%dm%size())
    end if
  end subroutine set_bousinesq_forcing_term
```

@note Notice the use of the `neko_field_registry` to retrieve the velocity and
scalar fields. See [Registries](@ref user-file_tips_registries) for more information
about registries in neko. @note Notice the use of the `NEKO_BCKND_DEVICE` flag,
which will be set to 1 if running on GPUs, and the use of `device_` functions.
See [Running on GPUs](@ref user-file_tips_running-on-gpus) for more information on
how this works.

The same can be done for the scalar, with the example below also taken from the
[scalar_mms example](https://github.com/ExtremeFLOW/neko/blob/49925b7a04a638259db3b1ddd54349ca57f5d207/examples/scalar_mms/scalar_mms.f90#L28-L47):

```fortran

  !> Set source term
  subroutine set_source(f, t)
    class(scalar_user_source_term_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t
    real(kind=rp) :: x, y
    integer :: i

    do i = 1, f%dm%size()
       x = f%dm%x(i,1,1,1)
       y = f%dm%y(i,1,1,1)

       ! 0.01 is the viscosity
       f%s(i,1,1,1) = cos(x) - 0.01 * sin(x) - 1.0_rp
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(f%s, f%s_d, f%dm%size(), &
                          HOST_TO_DEVICE, sync=.false.)
    end if

  end subroutine set_source

```

@note Notice the use of the `NEKO_BCKND_DEVICE` flag, which will be set to 1 if
running on GPUs, and the call to `device_memcpy` to transfer data between the
host and the device. See [Running on GPUs](@ref user-file_tips_running-on-gpus) for
more information on how this works.

We should also add of the following lines in `user_setup`, registering our user
functions `set_boussinesq_forcing_term` and `set_source` to be used as the fluid
and scalar source terms:

```fortran
u%fluid_user_f_vector => set_boussinesq_forcing_term
u%scalar_user_f_vector => set_source
```

### Complex fluid and/or scalar boundary conditions {#user-file_field-dirichlet-update}

This user function can be used to specify Dirichlet boundary values for velocity
components `u,v,w`, the pressure `p`, and/or the scalar `s`. This type of
boundary condition allows for time-dependent velocity profiles (currently not
possible with the `user_pointiwse` boundary condition) or non-uniform pressure
profiles to e.g. impose an outlet pressure computed from another simulation.

The user routine is called by the `user_velocity` and `user_pressure` boundary
conditions for the fluid, and the `user` boundary condtition for the scalar.
Once the appropriate boundaries have been identified, the user function
`field_dirichlet_update` should be used to compute and apply the desired values
to our velocity/pressure/scalar field(s). The prefix "field" in
`field_dirichlet_update` refers to the fact that a list of entire fields is
passed down for the user to edit.

The fields that are passed down depends on the underlying boundary conditions.
For `user_velocity`, it will be a list of 3 fields with names `u`, `v`, `w`; For
`user_pressure`, a list with a single field `p`; For the scalar, also a single
field with the same name as the solution field for the scalar (`s` by default).

It is crucial to understand that all three boundary condition will call the same
routine! So, if one has, for example, both `user_velocity` for the fluid and
`user` for the scalar, it is necessary to have an `if` statement in the user
routine to distinguish between the two cases. The convenient way to do that is
to check the size of the passed field list and the names of the fields inside.
For example, if there is one field and it is called `s`, one executes the code
setting the boundary values for the scalar `s`.

Note that the fields that one manipulates in the user routine are not the actual
velocity fields, etc. Instead, the code does a masked copy from the dummy fields
manipulated in the user routine to the actual fields of unknowns, with the mask
corresponding to the boundary faces. So, even if you somehow manipulate the
fields elsewhere in the domain inside the user routine, that will not affect the
solution.

In the following example, we indicate in `case.fluid.boundary_conditions` that
we would like to apply a velocity profile on the boundary number 1 (in this
case, the inlet boundary). On boundary number 2 (the outlet boundary), we change
the pressure. In `case.scalar.boundary_conditions`, we indicate the same for the
scalar on boundaries 1 and 2 (inlet and outlet).

The header of the user function is given in the code snippet below.

```fortran
  subroutine dirichlet_update(field_bc_list, bc, coef, t, tstep)
    type(field_list_t), intent(inout) :: field_bc_list
    type(bc_list_t), intent(inout) :: bc_bc_list
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
```

The arguments and their purpose are as follows:

* `field_bc_list` is the list of the field that can be edited. It is a list
of `field_t` objects.
  * The field `i` contained in `field_bc_list` is accessed using
  `field_bc_list%%items(i)%%ptr` and will refer to a `field_t` object.
* `bc` contains a `field_dirichlet_t` object to help access the boundary indices
  through the boundary mask, `msk`.
  * The boundary mask of the `bc `object is accessed via `bc%%msk`. It contains
  the linear indices of each GLL point on the boundary facets. @note
  `msk(0)` contains the size of the array. The first boundary index is `msk(1)`.
* `coef` is a `coef_t` object containing various numerical parameters and
variables, such as the polynomial order `lx`, derivatives, facet normals...
* `t`, `tstep` are self-explanatory.

Links to the documentation to learn more about what the types mentioned above
contain and how to use them: `src/field/field.f90`, `src/bc/bc.f90`,
`src/sem/coef.f90`.

The user function should be registered in `user_setup` with the following line:

```fortran
u%user_dirichlet_update => dirichlet_update
```

A very simple example illustrating the above is shown below, which is taken from the
[cyl_boundary_layer example](https://github.com/ExtremeFLOW/neko/blob/feature/field_bcs/examples/cyl_boundary_layer/cyl_bl.f90)

```fortran
  !! Initial example of using user specified dirichlet bcs
  !! Parameters:
  !! -----------
  !! field_bc_list:     List of fields from which the BC conditions zill be extracted.
  !!                    If called by the fluid, contains (u,v,w,p).
  !!                    If called by the scalar, contains (s).
  !! bc:                The BC containing the boundary mask, etc.
  !! coef:              Coef object.
  !! t:                 Current time.
  !! tstep:             Current time step.
  subroutine dirichlet_update(field_bc_list, bc, coef, t, tstep)
    type(field_list_t), intent(inout) :: field_bc_list
    type(field_dirichlet_t), intent(in) :: bc
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    integer :: i
    real(kind=rp) :: y,z

    ! Only do this at the first time step since our BCs are constant in time.
    if (tstep .ne. 1) return

    ! Check that we are being called by the fluid via the name of the field
    if (field_bc_list%items(1)%ptr%name .eq. "u") then

       associate(u => field_bc_list%items(1)%ptr, &
                 v => field_bc_list%items(2)%ptr, &
                 w => field_bc_list%items(3)%ptr)
         !
         ! Perform operations on u%x, v%x, w%x here
         ! Here we are applying very simple uniform boundaries (u,v,w) = (1,0,0)
         ! Technically the values are put in the interior as well, but this
         ! does not matter, only the boundary values will be copied to the
         ! actual fields
         u = 1.0_rp
         v = 0.0_rp
         w = 0.0_rp

       end associate

    ! Check that we are being called by the user_pressure bc via the name
    ! of the field
    else if (field_bc_list%items(1)%ptr%name .eq. "p") then
       associate( p => field_bc_list%items(1)%ptr)
         !
         ! Perform operations on the pressure field here
         !

         do i = 1, bc%msk(0)
            p%x(bc%msk(i), 1, 1, 1) = -1
         end do

       end associate

    ! Check that we are being called by the scalar via the name of the field
    else if (field_bc_list%items(1)%ptr%name .eq. "s") then

       associate( s => field_bc_list%items(1)%ptr)
         !
         ! Perform operations on the scalar field here
         !

         do i = 1, bc%msk(0)
            y = bc%dof%y(bc%msk(i), 1, 1, 1)
            z = bc%dof%z(bc%msk(i), 1, 1, 1)
            s%x(bc%msk(i), 1, 1, 1) = sin(y)*sin(z)
         end do

       end associate

    end if
  end subroutine dirichlet_update
```

This example is applying constant dirichlet values at the selected boundaries
for the velocity components and presure. The scalar is applied a function
`s(y,z) = sin(y)*sin(z)` to demonstrate the usage of boundary masks.

@attention The notation `u = 1.0_rp` is only possible because of the overloading
of the assignement operator `=` in `field_t`. In general, a field's array should
be accessed and modified with `u%%x`.

Note that we are only applying our boundary values at the first timestep, which
is done simply with the line `if (tstep .ne. 1) return`. This is a trick that
can be used for time-independent boundary profiles that require some kind of
time consuming operation like interpolation or reading from a file, which would
add overhead if executed at every time step.

@attention All the rules for [Running on GPUs](@ref
user-file_tips_running-on-gpus) apply when working on field arrays. Use
`device_memcpy` to make sure the device arrays are also updated.

## Additional remarks and tips

### Running on GPUs {#user-file_tips_running-on-gpus}

When running on GPUs, special care must be taken when using certain user
functions. The short explanation is that the device (GPU) has its own memory and
cannot directly access the memory on the host (CPU). This means that data and
more specifically arrays must be copied manually from the host to the device
(see device::device_memcpy).

@attention In some cases, data transfer via `device_memcpy` is avoidable. Neko
has some device math functions implemented that operate directly on device
arrays. If you can decompose whatever operations you are performing in a user
function into a set of instructions from the `math` module (e.g. `cadd`,
`cfill`, `sub2`, ...), you may use the corresponding `device_math` functions to
[offload work to the GPU](@ref accelerators_offload-work). See the
[fluid forcing code snippet](@ref user-file_user-f) for a simple example. For
more advanced examples, see the
[rayleigh_benard example](https://github.com/ExtremeFLOW/neko/blob/49925b7a04a638259db3b1ddd54349ca57f5d207/examples/rayleigh_benard/rayleigh.f90#L96-119)
or the
[tgv example](https://github.com/ExtremeFLOW/neko/blob/49925b7a04a638259db3b1ddd54349ca57f5d207/examples/tgv/tgv.f90#L146-172).

To illustrate this, let us have a look at the
[fluid initial condition code snippet](@ref user-file_user-ic):

```fortran

  !> Set the advecting velocity field.
  subroutine set_velocity(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params
    integer :: i, e, k, j
    real(kind=rp) :: x, y

    !
    ! 1. Set the initial condition in fields u%x, v%x, w%x
    !
    do i = 1, u%dof%size()
       x = u%dof%x(i,1,1,1)
       y = u%dof%y(i,1,1,1)

       ! Angular velocity is pi, giving a full rotation in 2 sec
       u%x(i,1,1,1) = -y*pi
       v%x(i,1,1,1) = x*pi
       w%x(i,1,1,1) = 0
    end do

    !
    ! 2. Copy the data set in u%x, v%x, w%x to the device arrays
    ! u%x_d, v%x_d, w%x_d.
    !
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(u%x, u%x_d, u%dof%size(), &
                          HOST_TO_DEVICE, sync=.false.)
       call device_memcpy(v%x, v%x_d, v%dof%size(), &
                          HOST_TO_DEVICE, sync=.false.)
       call device_memcpy(w%x, w%x_d, w%dof%size(), &
                          HOST_TO_DEVICE, sync=.false.)
    end if

  end subroutine set_velocity

```

The code above is used to set the fluid initial condition, by specifying the
values of fields `u,v,w` (and `p`) at all points in the domain. Notice that we
have divided the above code into two parts.

In the first part, we set the velocity components `u=-y*pi*`, `v=x*pi*`, and
`w=0`, which updates the velocity field arrays `u%%x, v%%x, w%%x` allocated on
the **host (CPU)**. If we were to run on GPUs, these lines of code would only
act on the velocity arrays on the host (CPU), leaving the device (GPU) arrays
untouched.

We take care of this in the second part, for all three velocity arrays. To
update the device (GPU) arrays, we use `device_memcpy` to copy the data
contained in a host (CPU) array to a device (GPU) array. Looking at the details
of the `device_memcpy` calls, we note the following:
- Device arrays are refered to by appending the suffix `_d` to the host array
  variable name (e.g. `u%%x` and `u%%x_d`). This is the standard in Neko.
- We specify the direction of the data movement with the flag `HOST_TO_DEVICE`.
  Other flags can also be used to move data from device to host
  (`DEVICE_TO_HOST`) or device to device (`DEVICE_TO_DEVICE`). See the
  [accelerators page](@ref accelerators_data-transfer) for more details on this.
- The `sync` argument is a non-optional argument which dictates wether or not to
  perform the data transfer synchronously.

@attention Use asynchronous data transfers at your own risk! If you are unsure,
use `sync = .true.` as a starting point.

Finally, observe that we use the flag `NEKO_BCKND_DEVICE` to check if we are
indeed running on GPUs. In that case, `NEKO_BCKND_DEVICE` would be equal to 1.

### Registries {#user-file_tips_registries}

Neko uses the concept of `registry` as a practical way to retrieve fields and
point zones anywhere in the user file.

The field registry `neko_field_registry` is often used in user functions where
certain fields are not directly accessible as arguments. One can retrieve any
field in the registry by its `name` with `neko_field_registry%%get_field(name)`.
Default fields that are added to the registry are `u,v,w,p` and `s` if running
with the scalar enabled. For a practical example of usage, see the
[rayleigh benard example](https://github.com/ExtremeFLOW/neko/blob/49925b7a04a638259db3b1ddd54349ca57f5d207/examples/rayleigh_benard/rayleigh.f90#L102-L105)

Other fields may be added to the registry by various simulation components. For
example:
- If running with `simulation_components.vorticity` enabled, the fields
  `omega_x, omega_y, omega_z` will be accessible in the registry.
- If running with `simulation_components.lambda2` enabled, the field `lambda2`
  will be accessible in the registry.

@note You can add your own fields to the registry with `neko_field_registry%%add_field` (see field_registry::field_add).

The point zone registry, `neko_point_zone_registry`, can be used to retrieve
pointers to `point_zone_t` objects defined in the case file. See
[using point zones](#point-zones_using-point-zones) for detailed instructions.
