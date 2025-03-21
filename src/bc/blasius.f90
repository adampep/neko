! Copyright (c) 2025, The Neko Authors
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
!> Defines a Blasius profile dirichlet condition
module blasius
  use num_types, only : rp
  use coefs, only : coef_t
  use utils, only : nonlinear_index
  use device, only : HOST_TO_DEVICE, device_memcpy, device_free, device_alloc
  use device_inhom_dirichlet
  use flow_profile
  use, intrinsic :: iso_fortran_env
  use, intrinsic :: iso_c_binding
  use bc, only : bc_t
  use json_module, only : json_file
  use json_utils, only : json_get
  implicit none
  private

  !> Blasius profile for inlet (vector valued).
  !! @warning Works only with axis-aligned sugar-cube elements and assumes
  !! the boundary is alinged with zOy.
  type, public, extends(bc_t) :: blasius_t
     real(kind=rp), dimension(3) :: uinf = [0d0, 0d0, 0d0]
     real(kind=rp) :: delta
     procedure(blasius_profile), nopass, pointer :: bla => null()
     type(c_ptr), private :: blax_d = C_NULL_PTR
     type(c_ptr), private :: blay_d = C_NULL_PTR
     type(c_ptr), private :: blaz_d = C_NULL_PTR
   contains
     procedure, pass(this) :: apply_scalar => blasius_apply_scalar
     procedure, pass(this) :: apply_vector => blasius_apply_vector
     procedure, pass(this) :: apply_scalar_dev => blasius_apply_scalar_dev
     procedure, pass(this) :: apply_vector_dev => blasius_apply_vector_dev
     procedure, pass(this) :: set_params => blasius_set_params
     !> Constructor
     procedure, pass(this) :: init => blasius_init
     !> Constructor from components
     procedure, pass(this) :: init_from_components => &
          blasius_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => blasius_free
     !> Finalize.
     procedure, pass(this) :: finalize => blasius_finalize
  end type blasius_t

contains

  !> Constructor
  !! @param[in] coef The SEM coefficients.
  !! @param[inout] json The JSON object configuring the boundary condition.
  subroutine blasius_init(this, coef, json)
    class(blasius_t), intent(inout), target :: this
    type(coef_t), intent(in) :: coef
    type(json_file), intent(inout) :: json
    real(kind=rp) :: delta
    real(kind=rp), allocatable :: uinf(:)
    character(len=:), allocatable :: approximation

    call this%init_base(coef)

    call json_get(json, 'delta', delta)
    call json_get(json, 'approximation', approximation)
    call json_get(json, 'freestream_velocity', uinf)

    if (size(uinf) .ne. 3) then
       call neko_error("The uinf keyword for the blasius profile should be an &
& array of 3 reals")
    end if

    call this%init_from_components(coef, delta, uinf, approximation)

  end subroutine blasius_init

  !> Constructor from components
  !! @param[in] coef The SEM coefficients.
  !! @param[in] delta The boundary layer thickness.
  !! @param[in] uinf The freestream velocity.
  !! @param[in] approximation The type of approximation for the profile.
  subroutine blasius_init_from_components(this, coef, delta, uinf, &
       approximation)
    class(blasius_t), intent(inout), target :: this
    type(coef_t), intent(in) :: coef
    real(kind=rp) :: delta
    real(kind=rp) :: uinf(3)
    character(len=*) :: approximation

    call this%init_base(coef)

    this%delta = delta
    this%uinf = uinf

    select case (trim(approximation))
    case ('linear')
       this%bla => blasius_linear
    case ('quadratic')
       this%bla => blasius_quadratic
    case ('cubic')
       this%bla => blasius_cubic
    case ('quartic')
       this%bla => blasius_quartic
    case ('sin')
       this%bla => blasius_sin
    case default
       call neko_error('Invalid Blasius approximation')
    end select
  end subroutine blasius_init_from_components

  subroutine blasius_free(this)
    class(blasius_t), target, intent(inout) :: this

    call this%free_base()
    nullify(this%bla)

    if (c_associated(this%blax_d)) then
       call device_free(this%blax_d)
    end if

    if (c_associated(this%blay_d)) then
       call device_free(this%blay_d)
    end if

    if (c_associated(this%blaz_d)) then
       call device_free(this%blaz_d)
    end if

  end subroutine blasius_free

  !> No-op scalar apply
  subroutine blasius_apply_scalar(this, x, n, t, tstep, strong)
    class(blasius_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    logical, intent(in), optional :: strong
  end subroutine blasius_apply_scalar

  !> No-op scalar apply (device version)
  subroutine blasius_apply_scalar_dev(this, x_d, t, tstep, strong)
    class(blasius_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    logical, intent(in), optional :: strong
  end subroutine blasius_apply_scalar_dev

  !> Apply blasius conditions (vector valued)
  subroutine blasius_apply_vector(this, x, y, z, n, t, tstep, strong)
    class(blasius_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(inout),  dimension(n) :: y
    real(kind=rp), intent(inout),  dimension(n) :: z
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    logical, intent(in), optional :: strong
    integer :: i, m, k, idx(4), facet
    logical :: strong_ = .true.

    if (present(strong)) strong_ = strong

    associate(xc => this%coef%dof%x, yc => this%coef%dof%y, &
              zc => this%coef%dof%z, nx => this%coef%nx, ny => this%coef%ny, &
              nz => this%coef%nz, lx => this%coef%Xh%lx)
      m = this%msk(0)
      if (strong_) then
         do i = 1, m
            k = this%msk(i)
            facet = this%facet(i)
            idx = nonlinear_index(k, lx, lx, lx)
            select case (facet)
            case (1, 2)
               x(k) = this%bla(zc(idx(1), idx(2), idx(3), idx(4)), &
                  this%delta, this%uinf(1))
               y(k) = 0.0_rp
               z(k) = 0.0_rp
            case (3, 4)
               x(k) = 0.0_rp
               y(k) = this%bla(xc(idx(1), idx(2), idx(3), idx(4)), &
                  this%delta, this%uinf(2))
               z(k) = 0.0_rp
            case (5, 6)
               x(k) = 0.0_rp
               y(k) = 0.0_rp
               z(k) = this%bla(yc(idx(1), idx(2), idx(3), idx(4)), &
                  this%delta, this%uinf(3))
            end select
         end do
      end if
    end associate
  end subroutine blasius_apply_vector

  !> Apply blasius conditions (vector valued) (device version)
  subroutine blasius_apply_vector_dev(this, x_d, y_d, z_d, t, tstep, strong)
    class(blasius_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    type(c_ptr) :: y_d
    type(c_ptr) :: z_d
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    logical, intent(in), optional :: strong
    integer :: i, m, k, idx(4), facet
    integer(c_size_t) :: s
    real(kind=rp), allocatable :: bla_x(:), bla_y(:), bla_z(:)
    logical :: strong_ = .true.

    if (present(strong)) strong_ = strong

    associate(xc => this%coef%dof%x, yc => this%coef%dof%y, &
              zc => this%coef%dof%z, nx => this%coef%nx, ny => this%coef%ny, &
              nz => this%coef%nz, lx => this%coef%Xh%lx , &
              blax_d => this%blax_d, blay_d => this%blay_d, &
              blaz_d => this%blaz_d)

      m = this%msk(0)


      ! Pretabulate values during first call to apply
      if (.not. c_associated(blax_d) .and. strong_ .and. this%msk(0) .gt. 0) then
         allocate(bla_x(m), bla_y(m), bla_z(m)) ! Temp arrays

         if (rp .eq. REAL32) then
            s = m * 4
         else if (rp .eq. REAL64) then
            s = m * 8
         end if

         call device_alloc(blax_d, s)
         call device_alloc(blay_d, s)
         call device_alloc(blaz_d, s)

         do i = 1, m
            k = this%msk(i)
            facet = this%facet(i)
            idx = nonlinear_index(k, lx, lx, lx)
            select case (facet)
            case (1,2)
               bla_x(i) = this%bla(zc(idx(1), idx(2), idx(3), idx(4)), &
                    this%delta, this%uinf(1))
               bla_y(i) = 0.0_rp
               bla_z(i) = 0.0_rp
            case (3,4)
               bla_x(i) = 0.0_rp
               bla_y(i) = this%bla(xc(idx(1), idx(2), idx(3), idx(4)), &
                    this%delta, this%uinf(2))
               bla_z(i) = 0.0_rp
            case (5,6)
               bla_x(i) = 0.0_rp
               bla_y(i) = 0.0_rp
               bla_z(i) = this%bla(yc(idx(1), idx(2), idx(3), idx(4)), &
                    this%delta, this%uinf(3))
            end select
         end do

         call device_memcpy(bla_x, blax_d, m, HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(bla_y, blay_d, m, HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(bla_z, blaz_d, m, HOST_TO_DEVICE, sync = .true.)

         deallocate(bla_x, bla_y, bla_z)
      end if

      if (strong_ .and. this%msk(0) .gt. 0) then
         call device_inhom_dirichlet_apply_vector(this%msk_d, x_d, y_d, z_d, &
              blax_d, blay_d, blaz_d, m)
      end if

    end associate

  end subroutine blasius_apply_vector_dev

  !> Set Blasius parameters
  subroutine blasius_set_params(this, uinf, delta, type)
    class(blasius_t), intent(inout) :: this
    real(kind=rp), intent(in) :: uinf(3)
    real(kind=rp), intent(in) :: delta
    character(len=*) :: type
    this%delta = delta
    this%uinf = uinf

    select case (trim(type))
    case ('linear')
       this%bla => blasius_linear
    case ('quadratic')
       this%bla => blasius_quadratic
    case ('cubic')
       this%bla => blasius_cubic
    case ('quartic')
       this%bla => blasius_quartic
    case ('sin')
       this%bla => blasius_sin
    case default
       call neko_error('Invalid Blasius approximation')
    end select
  end subroutine blasius_set_params

  !> Finalize
  subroutine blasius_finalize(this)
    class(blasius_t), target, intent(inout) :: this

    call this%finalize_base()
  end subroutine blasius_finalize
end module blasius
