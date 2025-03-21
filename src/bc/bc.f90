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
!> Defines a boundary condition
module bc
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use device, only : device_get_ptr, HOST_TO_DEVICE, device_memcpy, &
       device_free, device_map, DEVICE_TO_HOST
  use iso_c_binding, only: c_associated
  use dofmap, only : dofmap_t
  use coefs, only : coef_t
  use space, only : space_t
  use mesh, only : mesh_t, NEKO_MSH_MAX_ZLBLS, NEKO_MSH_MAX_ZLBL_LEN
  use facet_zone, only : facet_zone_t
  use stack, only : stack_i4t2_t
  use tuple, only : tuple_i4_t
  use field, only : field_t
  use gather_scatter
  use math, only : relcmp
  use utils, only : neko_error, linear_index, split_string
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR
  use json_module, only : json_file
  implicit none
  private

  !> Base type for a boundary condition
  type, public, abstract :: bc_t
     !> The linear index of each node in each boundary facet
     integer, allocatable :: msk(:)
     !> A list of facet ids (1 to 6), one for each element in msk
     integer, allocatable :: facet(:)
     !> Map of degrees of freedom
     type(dofmap_t), pointer :: dof
     !> SEM coefficients
     type(coef_t), pointer :: coef
     !> The mesh
     type(mesh_t), pointer :: msh
     !> The function space
     type(space_t), pointer :: Xh
     !> Index tuples (facet, element) marked as part of the boundary condition
     type(stack_i4t2_t) :: marked_facet
     !> Device pointer for msk
     type(c_ptr) :: msk_d = C_NULL_PTR
     !> Device pointer for facet
     type(c_ptr) :: facet_d = C_NULL_PTR
     !> Wether the bc is strongly enforced. Essentially valid for all Dirichlet
     !! types of bcs. These need to be masked out for solvers etc, so that
     !! values are not affected.
     !! Mixed bcs are, by convention, weak.
     logical :: strong = .true.
   contains
     !> Constructor
     procedure, pass(this) :: init_base => bc_init_base
     !> Destructor
     procedure, pass(this) :: free_base => bc_free_base
     !> Mark a facet on an element as part of the boundary condition
     procedure, pass(this) :: mark_facet => bc_mark_facet
     !> Mark all facets from a (facet, element) tuple list
     procedure, pass(this) :: mark_facets => bc_mark_facets
     !> Mark all facets from a zone
     procedure, pass(this) :: mark_zone => bc_mark_zone
     !> Finalize the construction of the bc by populating the msk and facet
     !! arrays
     procedure, pass(this) :: finalize_base => bc_finalize_base

     !> Apply the boundary condition to a scalar field. Dispatches to the CPU
     !! or the device version.
     procedure, pass(this) :: apply_scalar_generic => bc_apply_scalar_generic
     !> Apply the boundary condition to a vector field. Dispatches to the CPU
     !! or the device version.
     procedure, pass(this) :: apply_vector_generic => bc_apply_vector_generic
     !> Write a field showing the mask of the bcs
     procedure, pass(this) :: debug_mask_ => bc_debug_mask
     !> Apply the boundary condition to a scalar field on the CPU.
     procedure(bc_apply_scalar), pass(this), deferred :: apply_scalar
     !> Apply the boundary condition to a vector field on the CPU.
     procedure(bc_apply_vector), pass(this), deferred :: apply_vector
     !> Device version of \ref apply_scalar on the device.
     procedure(bc_apply_scalar_dev), pass(this), deferred :: apply_scalar_dev
     !> Device version of \ref apply_vector on the device.
     procedure(bc_apply_vector_dev), pass(this), deferred :: apply_vector_dev
     !> Deferred destructor.
     procedure(bc_destructor), pass(this), deferred :: free
     !> Deferred constructor.
     procedure(bc_constructor), pass(this), deferred :: init
     !> Deferred finalizer.
     procedure(bc_finalize), pass(this), deferred :: finalize
  end type bc_t

  !> Pointer to a @ref `bc_t`.
  type, public :: bc_ptr_t
     class(bc_t), pointer :: ptr => null()
  end type bc_ptr_t

  ! Helper type to have an array of polymorphic bc_t objects.
  type, public :: bc_alloc_t
     class(bc_t), allocatable :: obj
  end type bc_alloc_t


  abstract interface
     !> Constructor
     subroutine bc_constructor(this, coef, json)
       import :: bc_t, coef_t, json_file
       class(bc_t), intent(inout), target :: this
       type(coef_t), intent(in) :: coef
       type(json_file), intent(inout) :: json
     end subroutine bc_constructor
  end interface

  abstract interface
     !> Destructor
     subroutine bc_destructor(this)
       import :: bc_t
       class(bc_t), intent(inout), target :: this
     end subroutine bc_destructor
  end interface

  abstract interface
     !> Finalize by building the mask and facet arrays.
     subroutine bc_finalize(this)
       import :: bc_t
       class(bc_t), intent(inout), target :: this
     end subroutine bc_finalize
  end interface

  abstract interface
     !> Apply the boundary condition to a scalar field
     !! @param x The field for which to apply the boundary condition.
     !! @param n The size of x.
     !! @param t Current time.
     !! @param tstep Current time-step.
     !! @param strong Whether we are setting a strong or a weak bc.
     subroutine bc_apply_scalar(this, x, n, t, tstep, strong)
       import :: bc_t
       import :: rp
       class(bc_t), intent(inout) :: this
       integer, intent(in) :: n
       real(kind=rp), intent(inout), dimension(n) :: x
       real(kind=rp), intent(in), optional :: t
       integer, intent(in), optional :: tstep
       logical, intent(in), optional :: strong
     end subroutine bc_apply_scalar
  end interface

  abstract interface
     !> Apply the boundary condition to a vector field
     !! @param x The x comp of the field for which to apply the bc.
     !! @param y The y comp of the field for which to apply the bc.
     !! @param z The z comp of the field for which to apply the bc.
     !! @param n The size of x, y, and z.
     !! @param t Current time.
     !! @param tstep Current time-step.
     !! @param strong Whether we are setting a strong or a weak bc.
     subroutine bc_apply_vector(this, x, y, z, n, t, tstep, strong)
       import :: bc_t
       import :: rp
       class(bc_t), intent(inout) :: this
       integer, intent(in) :: n
       real(kind=rp), intent(inout), dimension(n) :: x
       real(kind=rp), intent(inout), dimension(n) :: y
       real(kind=rp), intent(inout), dimension(n) :: z
       real(kind=rp), intent(in), optional :: t
       integer, intent(in), optional :: tstep
       logical, intent(in), optional :: strong
     end subroutine bc_apply_vector
  end interface

  abstract interface
     !> Apply the boundary condition to a scalar field on the device
     !! @param x_d Device pointer to the field.
     !! @param t The time value.
     !! @param tstep The time iteration.
     !! @param strong Whether we are setting a strong or a weak bc.
     subroutine bc_apply_scalar_dev(this, x_d, t, tstep, strong)
       import :: c_ptr
       import :: bc_t
       import :: rp
       class(bc_t), intent(inout), target :: this
       type(c_ptr) :: x_d
       real(kind=rp), intent(in), optional :: t
       integer, intent(in), optional :: tstep
       logical, intent(in), optional :: strong
     end subroutine bc_apply_scalar_dev
  end interface

  abstract interface
     !> Apply the boundary condition to a vector field on the device.
     !! @param x_d Device pointer to the values to be applied for the x comp.
     !! @param y_d Device pointer to the values to be applied for the y comp.
     !! @param z_d Device pointer to the values to be applied for the z comp.
     !! @param t The time value.
     !! @param tstep Current time-step.
     !! @param strong Whether we are setting a strong or a weak bc.
     subroutine bc_apply_vector_dev(this, x_d, y_d, z_d, t, tstep, strong)
       import :: c_ptr
       import :: bc_t
       import :: rp
       class(bc_t), intent(inout), target :: this
       type(c_ptr) :: x_d
       type(c_ptr) :: y_d
       type(c_ptr) :: z_d
       real(kind=rp), intent(in), optional :: t
       integer, intent(in), optional :: tstep
       logical, intent(in), optional :: strong
     end subroutine bc_apply_vector_dev
  end interface

contains

  !> Constructor
  !! @param dof Map of degrees of freedom.
  subroutine bc_init_base(this, coef)
    class(bc_t), intent(inout) :: this
    type(coef_t), target, intent(in) :: coef

    call this%free_base

    this%dof => coef%dof
    this%coef => coef
    this%Xh => this%dof%Xh
    this%msh => this%dof%msh

    call this%marked_facet%init()

  end subroutine bc_init_base

  !> Destructor for the base type, `bc_t`.
  subroutine bc_free_base(this)
    class(bc_t), intent(inout) :: this

    call this%marked_facet%free()

    nullify(this%Xh)
    nullify(this%msh)
    nullify(this%dof)
    nullify(this%coef)

    if (allocated(this%msk)) then
       deallocate(this%msk)
    end if

    if (allocated(this%facet)) then
       deallocate(this%facet)
    end if

    if (c_associated(this%msk_d)) then
       call device_free(this%msk_d)
       this%msk_d = C_NULL_PTR
    end if

    if (c_associated(this%facet_d)) then
       call device_free(this%facet_d)
       this%facet_d = C_NULL_PTR
    end if

  end subroutine bc_free_base

  !> Apply the boundary condition to a vector field. Dispatches to the CPU
  !! or the device version.
  !! @param x The x comp of the field for which to apply the bc.
  !! @param y The y comp of the field for which to apply the bc.
  !! @param z The z comp of the field for which to apply the bc.
  !! @param n The size of x, y, and z.
  !! @param t Current time.
  !! @param tstep The current time iteration.
  subroutine bc_apply_vector_generic(this, x, y, z, n, t, tstep)
    class(bc_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    real(kind=rp), intent(inout), dimension(n) :: y
    real(kind=rp), intent(inout), dimension(n) :: z
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    type(c_ptr) :: x_d
    type(c_ptr) :: y_d
    type(c_ptr) :: z_d
    integer :: i


    if (NEKO_BCKND_DEVICE .eq. 1) then
       x_d = device_get_ptr(x)
       y_d = device_get_ptr(y)
       z_d = device_get_ptr(z)
       if (present(t) .and. present(tstep)) then
          call this%apply_vector_dev(x_d, y_d, z_d, t = t, tstep = tstep)
       else if (present(t)) then
          call this%apply_vector_dev(x_d, y_d, z_d, t = t)
       else if (present(tstep)) then
          call this%apply_vector_dev(x_d, y_d, z_d, tstep = tstep)
       else
          call this%apply_vector_dev(x_d, y_d, z_d)
       end if
    else
       if (present(t) .and. present(tstep)) then
          call this%apply_vector(x, y, z, n, t = t, tstep = tstep)
       else if (present(t)) then
          call this%apply_vector(x, y, z, n, t = t)
       else if (present(tstep)) then
          call this%apply_vector(x, y, z, n, tstep = tstep)
       else
          call this%apply_vector(x, y, z, n)
       end if
    end if

  end subroutine bc_apply_vector_generic

  !> Apply the boundary condition to a scalar field. Dispatches to the CPU
  !! or the device version.
  !! @param x The x comp of the field for which to apply the bc.
  !! @param y The y comp of the field for which to apply the bc.
  !! @param z The z comp of the field for which to apply the bc.
  !! @param n The size of x, y, and z.
  !! @param t Current time.
  !! @param tstep The current time iteration.
  subroutine bc_apply_scalar_generic(this, x, n, t, tstep)
    class(bc_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    real(kind=rp), intent(in), optional :: t
    integer, intent(in), optional :: tstep
    type(c_ptr) :: x_d
    integer :: i


    if (NEKO_BCKND_DEVICE .eq. 1) then
       x_d = device_get_ptr(x)
       if (present(t) .and. present(tstep)) then
          call this%apply_scalar_dev(x_d, t = t, tstep = tstep)
       else if (present(t)) then
          call this%apply_scalar_dev(x_d, t = t)
       else if (present(tstep)) then
          call this%apply_scalar_dev(x_d, tstep = tstep)
       else
          call this%apply_scalar_dev(x_d)
       end if
    else
       if (present(t) .and. present(tstep)) then
          call this%apply_scalar(x, n, t = t, tstep = tstep)
       else if (present(t)) then
          call this%apply_scalar(x, n, t = t)
       else if (present(tstep)) then
          call this%apply_scalar(x, n, tstep = tstep)
       else
          call this%apply_scalar(x, n)
       end if
    end if

  end subroutine bc_apply_scalar_generic

  !> Mark @a facet on element @a el as part of the boundary condition
  !! @param facet The index of the facet.
  !! @param el The index of the element.
  subroutine bc_mark_facet(this, facet, el)
    class(bc_t), intent(inout) :: this
    integer, intent(in) :: facet
    integer, intent(in) :: el
    type(tuple_i4_t) :: t

    t%x = [facet, el]
    call this%marked_facet%push(t)

  end subroutine bc_mark_facet

  !> Mark all facets from a (facet, el) tuple list
  !! @param facet_list The list of tuples.
  subroutine bc_mark_facets(this, facet_list)
    class(bc_t), intent(inout) :: this
    type(stack_i4t2_t), intent(inout) :: facet_list
    type(tuple_i4_t), pointer :: fp(:)
    integer :: i

    fp => facet_list%array()
    do i = 1, facet_list%size()
       call this%marked_facet%push(fp(i))
    end do

  end subroutine bc_mark_facets

  !> Mark all facets from a zone
  !! @param bc_zone Boundary zone to be marked.
  subroutine bc_mark_zone(this, bc_zone)
    class(bc_t), intent(inout) :: this
    class(facet_zone_t), intent(in) :: bc_zone
    integer :: i
    do i = 1, bc_zone%size
       call this%marked_facet%push(bc_zone%facet_el(i))
    end do
  end subroutine bc_mark_zone

  !> Finalize the construction of the bc by populting the `msk` and `facet`
  !! arrays.
  !! @param only_facets, if the bc is only to be applied on facets.
  !! Relevant for bcs where the normal direction is important
  !! and where and shared dofs should not be included.
  !! @details This will linearize the marked facet's indicies in the msk array.
  !!
  subroutine bc_finalize_base(this, only_facets)
    class(bc_t), target, intent(inout) :: this
    logical, optional, intent(in) :: only_facets
    type(tuple_i4_t), pointer :: bfp(:)
    type(tuple_i4_t) :: bc_facet
    type(field_t) :: test_field
    integer :: facet_size, facet, el
    logical :: only_facet = .false.
    integer :: i, j, k, l, msk_c
    integer :: lx, ly, lz, n
    lx = this%Xh%lx
    ly = this%Xh%ly
    lz = this%Xh%lz
    if ( present(only_facets)) then
       only_facet = only_facets
    else
       only_facet = .false.
    end if
    !>@todo add 2D case

    ! Note we assume that lx = ly = lz
    facet_size = lx**2
    allocate(this%msk(0:facet_size * this%marked_facet%size()))
    allocate(this%facet(0:facet_size * this%marked_facet%size()))

    msk_c = 0
    bfp => this%marked_facet%array()

    ! Loop through each (facet, element) id tuple
    ! Then loop over all the nodes of the face and compute their linear index
    ! This index goes into this%msk, whereas the corresponding face id goes into
    ! this%facet
    do i = 1, this%marked_facet%size()
       bc_facet = bfp(i)
       facet = bc_facet%x(1)
       el = bc_facet%x(2)
       select case (facet)
       case (1)
          do l = 1, lz
             do k = 1, ly
                msk_c = msk_c + 1
                this%msk(msk_c) = linear_index(1, k, l, el, lx, ly, lz)
                this%facet(msk_c) = 1
             end do
          end do
       case (2)
          do l = 1, lz
             do k = 1, ly
                msk_c = msk_c + 1
                this%msk(msk_c) = linear_index(lx, k, l, el, lx, ly, lz)
                this%facet(msk_c) = 2
             end do
          end do
       case (3)
          do l = 1, lz
             do j = 1, lx
                msk_c = msk_c + 1
                this%msk(msk_c) = linear_index(j, 1, l, el, lx, ly, lz)
                this%facet(msk_c) = 3
             end do
          end do
       case (4)
          do l = 1, lz
             do j = 1, lx
                msk_c = msk_c + 1
                this%msk(msk_c) = linear_index(j, ly, l, el, lx, ly, lz)
                this%facet(msk_c) = 4
             end do
          end do
       case (5)
          do k = 1, ly
             do j = 1, lx
                msk_c = msk_c + 1
                this%msk(msk_c) = linear_index(j, k, 1, el, lx, ly, lz)
                this%facet(msk_c) = 5
             end do
          end do
       case (6)
          do k = 1, ly
             do j = 1, lx
                msk_c = msk_c + 1
                this%msk(msk_c) = linear_index(j, k, lz, el, lx, ly, lz)
                this%facet(msk_c) = 6
             end do
          end do
       end select
    end do
    if ( .not. only_facet) then
       !Makes check for points not on facet that should have bc applied
       call test_field%init(this%dof)

       n = test_field%size()
       test_field%x = 0.0_rp
       !Apply this bc once
       do i = 1, msk_c
          test_field%x(this%msk(i),1,1,1) = 1.0
       end do
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy(test_field%x, test_field%x_d, n, &
               HOST_TO_DEVICE, sync=.true.)
       end if
       !Check if some point that was not zeroed was zeroed on another element
       call this%coef%gs_h%op(test_field,GS_OP_ADD)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy(test_field%x, test_field%x_d, n, &
               DEVICE_TO_HOST, sync=.true.)
       end if
       msk_c = 0
       do i = 1, this%dof%size()
          if (test_field%x(i,1,1,1) .gt. 0.5) then
             msk_c = msk_c + 1
          end if
       end do
       !Allocate new mask
       deallocate(this%msk)
       allocate(this%msk(0:msk_c))
       j = 1
       do i = 1, this%dof%size()
          if (test_field%x(i,1,1,1) .gt. 0.5) then
             this%msk(j) = i
             j = j + 1
          end if
       end do

       call test_field%free()
    end if

    this%msk(0) = msk_c
    this%facet(0) = msk_c

    if (NEKO_BCKND_DEVICE .eq. 1) then
       n = msk_c + 1
       call device_map(this%msk, this%msk_d, n)
       call device_map(this%facet, this%facet_d, n)

       call device_memcpy(this%msk, this%msk_d, n, &
            HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(this%facet, this%facet_d, n, &
            HOST_TO_DEVICE, sync = .true.)
    end if

  end subroutine bc_finalize_base

  !> Write a field showing the mask of the bc
  !! @details The mask will be marked with 1.
  !! @param file_name The name of the fld file.
  subroutine bc_debug_mask(this, file_name)
    use field, only : field_t
    use file, only : file_t

    class(bc_t), intent(inout) :: this
    character(len=*), intent(in) :: file_name
    type(field_t) :: bdry_field
    integer:: i, m, k
    type(file_t) :: dump_file

    call bdry_field%init(this%coef%dof, 'bdry')
    m = this%msk(0)
    do i = 1, m
       k = this%msk(i)
       bdry_field%x(k,1,1,1) = 1.0_rp
    end do
    dump_file = file_t(file_name)
    call dump_file%write(bdry_field)

  end subroutine bc_debug_mask
end module bc
