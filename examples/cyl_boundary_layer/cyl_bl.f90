!>User file to set up at rotating cylinder. Martin Karp 30/5-2022

module user
  use neko
  implicit none
  ! Case parameters
  real(kind=rp), parameter :: h = 1.0
  real(kind=rp), parameter :: gam = 20.0
  real(kind=rp), parameter :: rad = h/gam
  real(kind=rp), parameter :: n = 7
  real(kind=rp), parameter :: pw = 1/n
  real(kind=rp), parameter :: ucl = 1
  real(kind=rp), parameter :: alpha = 1
  real(kind=rp), parameter :: u_th2 = ucl*alpha
  real(kind=rp), parameter :: u_rho = 0.0
  real(kind=rp), parameter :: u_axial = 0.0
  real(kind=rp), parameter :: y0 = 0.0
  real(kind=rp), parameter :: y1 = 0.0
  real(kind=rp), parameter :: delta = 0.005*h
contains


  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%fluid_user_ic => user_ic
    u%fluid_user_if => user_inflow_eval
    u%user_dirichlet_update => dirichlet_update
  end subroutine user_setup

  subroutine cylinder_deform(msh)
    type(mesh_t), intent(inout) :: msh
    msh%apply_deform => cylinder_gen_curve
  end subroutine cylinder_deform

  !> Make the cylinder more circular!
  subroutine cylinder_gen_curve(msh, x, y, z, lx, ly, lz)
    class(mesh_t) :: msh
    integer, intent(in) :: lx, ly, lz
    real(kind=rp), intent(inout) :: x(lx, lx, lx, msh%nelv)
    real(kind=rp), intent(inout) :: y(lx, lx, lx, msh%nelv)
    real(kind=rp), intent(inout) :: z(lx, lx, lx, msh%nelv)
    type(tuple_i4_t) :: el_and_facet
    real(kind=rp) :: th
    integer :: e, i, j ,k, l,  facet

    !The cylinders zone number is 7
    do l = 1, msh%labeled_zones(7)%size
       el_and_facet = msh%labeled_zones(7)%facet_el(l)
       facet = el_and_facet%x(1)
       e = el_and_facet%x(2)
       do k = 1, lz
          do j = 1, ly
              do i = 1, lx
                 if (index_is_on_facet(i, j, k, lx, ly, lz, facet)) then
                    th = atan2(z(i,j,k,e), x(i,j,k,e))
                    x(i,j,k,e) = rad * cos(th)
                    z(i,j,k,e) = rad * sin(th)
                 end if
              end do
          end do
       end do
    end do
  end subroutine cylinder_gen_curve

  subroutine user_inflow_eval(u, v, w, x, y, z, nx, ny, nz, ix, iy, iz, ie, t, &
       tstep)
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
    real(kind=rp) ::  u_th, dist, th, yy
    real(kind=rp) ::  arg

!   Two different regions (inflow & cyl) have the label 'v  '
!   Let compute the distance from the (0,0) in the x-y plane
!   to identify the proper one
    dist = sqrt(x**2 + z**2)

! --- INFLOW
    if (dist .gt. 1.1*rad) then
       u =  ucl*y**pw
    end if
! ---

    w = 0.0
    v = 0.0
! --- SPINNING CYLINDER

    if (dist .lt. 1.5*rad .and. y .gt. 0.1) then
       th = atan2(z, x)
       u = cos(th)*u_rho - sin(th)*u_th2
       w = sin(th)*u_rho + cos(th)*u_th2
    end if

! ---


!     Smoothing function for the velocity u_th on the spinning cylinder
!     to avoid gap in the at the bottom wall

!     u_th is smoothed if z0 < z < delta
!     u_th=1 if z >= delta


    yy = y + abs(y0) ! coordinate shift

    if (dist .lt. 1.5*rad) then
       if (yy .lt. delta) then
          arg  = yy/delta
          u_th = u_th2/(1.0_rp+exp(1.0_rp/(arg-1.0_rp)+1.0_rp/arg))
       else
          u_th = u_th2
       end if

       th = atan2(z,x)

       u = cos(th)*u_rho - sin(th)*u_th
       w = sin(th)*u_rho + cos(th)*u_th
    end if
  end subroutine user_inflow_eval

  ! User defined initial condition
  subroutine user_ic(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params
    integer :: i
    real(kind=rp) :: y

    do i = 1, u%dof%size()
       y = u%dof%y(i,1,1,1)
       u%x(i,1,1,1) =  ucl*y**pw
       v%x(i,1,1,1) = 0.0
       w%x(i,1,1,1) = 0.0
    end do
  end subroutine user_ic

  ! Initial example of using user specified dirichlet bcs
  ! Note: This subroutine will be called two times, once in the fluid solver, and once
  ! in the scalar solver (if enabled).
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

    ! Only do this at the first time step since our BCs are constants.
    if (tstep .ne. 1) return

    ! Check that we are being called by the fluid via the name of the field
    if (field_bc_list%items(1)%ptr%name .eq. "u") then

       associate(u => field_bc_list%items(1)%ptr, &
            v => field_bc_list%items(2)%ptr, &
            w => field_bc_list%items(3)%ptr)

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

end module user
