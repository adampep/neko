! Copyright (c) 2018-2023, The Neko Authors
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
!> Connectivity vertex type
module vertex_cnn
  use num_types, only : i4
  use utils, only : neko_error
  use polytope_cnn, only : polytope_cnn_t
  implicit none
  private

  public :: vertex_cnn_t, vertex_cnn_ptr, vertex_ncnf_cnn_t, vertex_ncnf_cnn_ptr

  ! object information
  integer(i4), public, parameter :: NEKO_VERTEX_DIM = 0
  integer(i4), public, parameter :: NEKO_VERTEX_NFACET = 0
  integer(i4), public, parameter :: NEKO_VERTEX_NRIDGE = 0
  integer(i4), public, parameter :: NEKO_VERTEX_NPEAK = 0

  !> Vertex type for global communication
  !! @details Vertex as the only realisation of zero-dimensional polytope
  !! (monon) and contains unique global id only. Vertex has no alignment.
  type, extends(polytope_cnn_t) :: vertex_cnn_t
   contains
     !> Initialise vertex
     procedure, pass(this) :: init => vertex_init
     !> vertex equality check
     procedure, pass(this) :: equal => vertex_equal
     generic :: operator(.eq.) => equal
  end type vertex_cnn_t

  !> Pointer to a vertex type
  type ::  vertex_cnn_ptr
     type(vertex_cnn_t), pointer :: obj
  end type vertex_cnn_ptr

  !> Vertex type for nonconforming meshes
  !! @details Vertex can be either independent (located at edge, face, cell
  !! corner of all neighbours; marked 0), facet hanging (located at facet
  !! centre of some neighbours in case of faces or cells; marked 1) or
  !! ridge hanging (located at ridge centre of some neighbours in case of
  !! cells; marked 2). There are no operations on vertex, so no procedure
  !! pointers.
  type :: vertex_ncnf_cnn_t
     ! vertex pointer
     type(vertex_cnn_ptr) :: vertex
     !> hanging information
     integer(i4) :: hanging = 0
   contains
     !> Initialise vertex pointer
     procedure, pass(this) :: init => vertex_hanging_init
     !> Free vertex data
     procedure, pass(this) :: free => vertex_hanging_free
     !> Set hanging information
     procedure, pass(this) :: set_hng => vertex_hanging_set
     !> Get hanging information
     procedure, pass(this) :: hng => vertex_hanging_get
  end type vertex_ncnf_cnn_t

  !> Pointer to a hanging vertex type
  type ::  vertex_ncnf_cnn_ptr
     type(vertex_ncnf_cnn_t), pointer :: obj
  end type vertex_ncnf_cnn_ptr

contains

  !> @brief Initialise vertex with global id
  !! @parameter[in]   id     unique id
  subroutine vertex_init(this, id)
    class(vertex_cnn_t), intent(inout) :: this
    integer(i4), intent(in) :: id

    call this%set_dim(NEKO_VERTEX_DIM)
    call this%set_nelem(NEKO_VERTEX_NFACET, NEKO_VERTEX_NRIDGE,&
         & NEKO_VERTEX_NPEAK)
    call this%set_id(id)

    return
  end subroutine vertex_init

  !> @brief Check if two vertices are the same
  !! @return   equal
  pure function vertex_equal(this, other) result(equal)
    class(vertex_cnn_t), intent(in) :: this
    class(polytope_cnn_t), intent(in) :: other
    logical :: equal

    ! check polygon information
    equal = this%equal_poly(other)
    if (equal) then
       ! check global id
       equal = (this%id() == other%id())
    end if
    return
  end function vertex_equal

  !> @brief Initialise vertex pointer
  !! @parameter[in]   vrt     vertex
  subroutine vertex_hanging_init(this, vrt)
    class(vertex_ncnf_cnn_t), intent(inout) :: this
    type(vertex_cnn_t), target, intent(in) :: vrt
    call this%free()
    this%vertex%obj => vrt
    return
  end subroutine vertex_hanging_init

  !> @brief free vertex pointer and hanging information
  subroutine vertex_hanging_free(this)
    class(vertex_ncnf_cnn_t), intent(inout) :: this
    this%vertex%obj => null()
    this%hanging = 0
    return
  end subroutine vertex_hanging_free

  !> @brief Set hanging information
  !! @parameter[in]   hng     hanging information
  subroutine vertex_hanging_set(this, hng)
    class(vertex_ncnf_cnn_t), intent(inout) :: this
    integer(i4), intent(in) :: hng
    if (hng >= 0 .and. hng <= 2) then
       this%hanging = hng
    else
       call neko_error('Inconsistent vertex hanging information.')
    end if
    return
  end subroutine vertex_hanging_set

  !> @brief Get hanging information
  !! @return   hng
  pure function vertex_hanging_get(this) result(hng)
    class(vertex_ncnf_cnn_t), intent(in) :: this
    integer(i4) :: hng
    hng = this%hanging
    return
  end function vertex_hanging_get

end module vertex_cnn