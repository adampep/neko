! Copyright (c) 2020-2023, The Neko Authors
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
!> Defines inflow dirichlet conditions
module usr_inflow
  use num_types, only : rp
  use coefs, only : coef_t
  use inflow, only : inflow_t
  use device
  use device_inhom_dirichlet
  use utils, only : neko_error, nonlinear_index, neko_warning
  use bc, only : bc_t
  use json_module, only : json_file
  implicit none
  private

  !> User defined dirichlet condition for velocity.
  type, public, extends(bc_t) :: usr_inflow_t
     procedure(usr_inflow_eval), nopass, pointer :: eval => null()
     type(c_ptr), private :: usr_x_d = C_NULL_PTR
     type(c_ptr), private :: usr_y_d = C_NULL_PTR
     type(c_ptr), private :: usr_z_d = C_NULL_PTR
   contains
     procedure, pass(this) :: apply_scalar => usr_inflow_apply_scalar
     procedure, pass(this) :: apply_vector => usr_inflow_apply_vector
     procedure, pass(this) :: validate => usr_inflow_validate
     procedure, pass(this) :: set_eval => usr_inflow_set_eval
     procedure, pass(this) :: apply_vector_dev => usr_inflow_apply_vector_dev
     procedure, pass(this) :: apply_scalar_dev => usr_inflow_apply_scalar_dev
     !> Constructor
     procedure, pass(this) :: init => usr_inflow_init
     !> Destructor.
     procedure, pass(this) :: free => usr_inflow_free
     !> Finalize.
     procedure, pass(this) :: finalize => usr_inflow_finalize
  end type usr_inflow_t

  abstract interface

     !> Abstract interface defining a user defined inflow condition (pointwise)
     !! @param u The x componenet of the velocity in this point
     !! @param v The y componenet of the velocity in this point
     !! @param w The w componenet of the velocity in this point
     !! @param x The x coord in this point
     !! @param y The y coord in this point
     !! @param z The z coord in this point
     !! @param nx The x component of the facet normal in this point
     !! @param ny The y component of the facet normal in this point
     !! @param nz The z component of the facet normal in this point
     !! @param ix The r idx of this point
     !! @param iy The s idx of this point
     !! @param iz The t idx of this point
     !! @param ie The element idx of this point
     !! @param t Current time
     !! @param tstep Current time-step
     subroutine usr_inflow_eval(u, v, w, x, y, z, nx, ny, nz, &
                                ix, iy, iz, ie, t, tstep)
       import rp
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
     end subroutine usr_inflow_eval
  end interface

  public :: usr_inflow_eval

contains

  !> Constructor
  !! @param[in] coef The SEM coefficients.
  !! @param[inout] json The JSON object configuring the boundary condition.
  subroutine usr_inflow_init(this, coef, json)
    class(usr_inflow_t), intent(inout), target :: this
    type(coef_t), intent(in) :: coef
    type(json_file), intent(inout) ::json

    call this%init_base(coef)
  end subroutine usr_inflow_init

  subroutine usr_inflow_free(this)
    class(usr_inflow_t), target, intent(inout) :: this

    call this%free_base()

    if (c_associated(this%usr_x_d)) then
       call device_free(this%usr_x_d)
    end if

    if (c_associated(this%usr_y_d)) then
       call device_free(this%usr_y_d)
    end if

    if (c_associated(this%usr_z_d)) then
       call device_free(this%usr_z_d)
    end if

  end subroutine usr_inflow_free

  !> No-op scalar apply
  subroutine usr_inflow_apply_scalar(this, x, n, t, tstep, strong)
    class(usr_inflow_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    logical, intent(in), optional :: strong
  end subroutine usr_inflow_apply_scalar

  !> No-op scalar apply (device version)
  subroutine usr_inflow_apply_scalar_dev(this, x_d, t, tstep, strong)
    class(usr_inflow_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    logical, intent(in), optional :: strong
  end subroutine usr_inflow_apply_scalar_dev

  !> Apply user defined inflow conditions (vector valued)
  subroutine usr_inflow_apply_vector(this, x, y, z, n, t, tstep, strong)
    class(usr_inflow_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout),  dimension(n) :: x
    real(kind=rp), intent(inout),  dimension(n) :: y
    real(kind=rp), intent(inout),  dimension(n) :: z
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    logical, intent(in), optional :: strong
    integer :: i, m, k, idx(4), facet, tstep_
    real(kind=rp) :: t_
    logical :: strong_ = .true.

    if (present(strong)) strong_ = strong

    if (present(t)) then
       t_ = t
    else
       t_ = 0.0_rp
    end if

    if (present(tstep)) then
       tstep_ = tstep
    else
       tstep_ = 1
    end if

    associate(&
         xc => this%coef%dof%x, yc => this%coef%dof%y, zc => this%coef%dof%z, &
         nx => this%coef%nx, ny => this%coef%ny, nz => this%coef%nz, &
         lx => this%coef%Xh%lx)
      m = this%msk(0)
      if (strong_) then
         do i = 1, m
            k = this%msk(i)
            facet = this%facet(i)
            idx = nonlinear_index(k, lx, lx, lx)
            select case (facet)
            case (1,2)
               call this%eval(x(k), y(k), z(k), &
                  xc(idx(1), idx(2), idx(3), idx(4)), &
                  yc(idx(1), idx(2), idx(3), idx(4)), &
                  zc(idx(1), idx(2), idx(3), idx(4)), &
                  nx(idx(2), idx(3), facet, idx(4)), &
                  ny(idx(2), idx(3), facet, idx(4)), &
                  nz(idx(2), idx(3), facet, idx(4)), &
                  idx(1), idx(2), idx(3), idx(4), &
                  t_, tstep_)
            case (3,4)
               call this%eval(x(k), y(k), z(k), &
                  xc(idx(1), idx(2), idx(3), idx(4)), &
                  yc(idx(1), idx(2), idx(3), idx(4)), &
                  zc(idx(1), idx(2), idx(3), idx(4)), &
                  nx(idx(1), idx(3), facet, idx(4)), &
                  ny(idx(1), idx(3), facet, idx(4)), &
                  nz(idx(1), idx(3), facet, idx(4)), &
                  idx(1), idx(2), idx(3), idx(4), &
                  t_, tstep_)
            case (5,6)
               call this%eval(x(k), y(k), z(k), &
                  xc(idx(1), idx(2), idx(3), idx(4)), &
                  yc(idx(1), idx(2), idx(3), idx(4)), &
                  zc(idx(1), idx(2), idx(3), idx(4)), &
                  nx(idx(1), idx(2), facet, idx(4)), &
                  ny(idx(1), idx(2), facet, idx(4)), &
                  nz(idx(1), idx(2), facet, idx(4)), &
                  idx(1), idx(2), idx(3), idx(4), &
                  t_, tstep_)
            end select
         end do
      end if
    end associate

  end subroutine usr_inflow_apply_vector

  subroutine usr_inflow_apply_vector_dev(this, x_d, y_d, z_d, t, tstep, strong)
    class(usr_inflow_t), intent(inout), target :: this
    type(c_ptr) :: x_d
    type(c_ptr) :: y_d
    type(c_ptr) :: z_d
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    logical, intent(in), optional :: strong
    integer :: i, m, k, idx(4), facet, tstep_
    integer(c_size_t) :: s
    real(kind=rp) :: t_
    real(kind=rp), allocatable :: x(:)
    real(kind=rp), allocatable :: y(:)
    real(kind=rp), allocatable :: z(:)
    logical :: strong_ = .true.

    if (present(strong)) strong_ = strong

    if (present(t)) then
       t_ = t
    else
       t_ = 0.0_rp
    end if

    if (present(tstep)) then
       tstep_ = tstep
    else
       tstep_ = 1
    end if

    associate(&
         xc => this%coef%dof%x, yc => this%coef%dof%y, zc => this%coef%dof%z, &
         nx => this%coef%nx, ny => this%coef%ny, nz => this%coef%nz, &
         lx => this%coef%Xh%lx, usr_x_d => this%usr_x_d, &
         usr_y_d => this%usr_y_d, usr_z_d => this%usr_z_d)

      m = this%msk(0)


      ! Pretabulate values during first call to apply
      if (.not. c_associated(usr_x_d) .and. strong_ .and. &
          (this%msk(0) .gt. 0)) then
         allocate(x(m), y(m), z(m)) ! Temp arrays

         s = m*rp

         call device_alloc(usr_x_d, s)
         call device_alloc(usr_y_d, s)
         call device_alloc(usr_z_d, s)

         associate(xc => this%coef%dof%x, yc => this%coef%dof%y, &
              zc => this%coef%dof%z, &
              nx => this%coef%nx, ny => this%coef%ny, nz => this%coef%nz, &
              lx => this%coef%Xh%lx)
           do i = 1, m
              k = this%msk(i)
              facet = this%facet(i)
              idx = nonlinear_index(k, lx, lx, lx)
              select case (facet)
              case (1,2)
                 call this%eval(x(i), y(i), z(i), &
                      xc(idx(1), idx(2), idx(3), idx(4)), &
                      yc(idx(1), idx(2), idx(3), idx(4)), &
                      zc(idx(1), idx(2), idx(3), idx(4)), &
                      nx(idx(2), idx(3), facet, idx(4)), &
                      ny(idx(2), idx(3), facet, idx(4)), &
                      nz(idx(2), idx(3), facet, idx(4)), &
                      idx(1), idx(2), idx(3), idx(4), &
                      t_, tstep_)
              case (3,4)
                 call this%eval(x(i), y(i), z(i), &
                      xc(idx(1), idx(2), idx(3), idx(4)), &
                      yc(idx(1), idx(2), idx(3), idx(4)), &
                      zc(idx(1), idx(2), idx(3), idx(4)), &
                      nx(idx(1), idx(3), facet, idx(4)), &
                      ny(idx(1), idx(3), facet, idx(4)), &
                      nz(idx(1), idx(3), facet, idx(4)), &
                      idx(1), idx(2), idx(3), idx(4), &
                      t_, tstep_)
              case (5,6)
                 call this%eval(x(i), y(i), z(i), &
                      xc(idx(1), idx(2), idx(3), idx(4)), &
                      yc(idx(1), idx(2), idx(3), idx(4)), &
                      zc(idx(1), idx(2), idx(3), idx(4)), &
                      nx(idx(1), idx(2), facet, idx(4)), &
                      ny(idx(1), idx(2), facet, idx(4)), &
                      nz(idx(1), idx(2), facet, idx(4)), &
                      idx(1), idx(2), idx(3), idx(4), &
                      t_, tstep_)
              end select
           end do
         end associate

         call device_memcpy(x, usr_x_d, m, HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(y, usr_y_d, m, HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(z, usr_z_d, m, HOST_TO_DEVICE, sync = .true.)

         deallocate(x, y, z)
      end if

      if (strong_ .and. (this%msk(0) .gt. 0)) then
         call device_inhom_dirichlet_apply_vector(this%msk_d, x_d, y_d, z_d, &
              usr_x_d, usr_y_d, usr_z_d, m)
      end if

    end associate

  end subroutine usr_inflow_apply_vector_dev

  !> Assign user provided eval function
  !! @param user_eval User specified boundary condition for u,v,w (vector)
  subroutine usr_inflow_set_eval(this, usr_eval)
    class(usr_inflow_t), intent(inout) :: this
    procedure(usr_inflow_eval) :: usr_eval
    this%eval => usr_eval
  end subroutine usr_inflow_set_eval

  !> Validate user inflow condition
  subroutine usr_inflow_validate(this)
    class(usr_inflow_t), intent(inout) :: this
    logical :: valid

    valid = .true. ! Assert it's going to be ok...
    if (.not. associated(this%coef)) then
       call neko_warning('Missing coefficients')
       valid = .false.
    end if

    if (.not. associated(this%eval)) then
       call neko_warning('Missing eval function')
       valid = .false.
    end if

    if (.not. valid) then
       call neko_error('Invalid user defined inflow condition')
    end if

  end subroutine usr_inflow_validate

  !> Finalize
  subroutine usr_inflow_finalize(this)
    class(usr_inflow_t), target, intent(inout) :: this

    call this%finalize_base()
  end subroutine usr_inflow_finalize
end module usr_inflow
