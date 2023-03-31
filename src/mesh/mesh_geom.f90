! Copyright (c) 2019-2021, The Neko Authors
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
module mesh_geom
  use num_types
  use element
  use point
  use quad
  use hex
  implicit none
  private

  public :: mesh_geom_t

  type, private :: mesh_element_t
     class(element_t), allocatable :: el
  end type mesh_element_t

  !> Base type for mesh geometrical information
  type :: mesh_geom_t
     integer(i4) :: nel ! local number of elements
     integer(i4) :: mpts ! local number of unique points in the mesh
     type(point_t), allocatable :: points(:) ! list of local unique points
     type(mesh_element_t), allocatable :: elements(:) ! list of local elements
     logical, allocatable :: dfrmd_el(:) ! element deformation flag
   contains
     procedure, pass(this) :: free => mesh_geom_free
  end type mesh_geom_t

contains

  !> Free memory
  subroutine mesh_geom_free(this)
    ! argument list
    class(mesh_geom_t), intent(inout) :: this

    integer(i4) :: il

    ! Deallocate arrays
    if (allocated(this%points)) deallocate(this%points)
    if (allocated(this%elements)) then
       do il = 1, this%nel
          call this%elements(il)%el%free()
       end do
       deallocate(this%elements)
    end if
    if (allocated(this%dfrmd_el)) deallocate(this%dfrmd_el)

    ! Reset registers
    this%nel = 0
    this%mpts = 0

    return
  end subroutine mesh_geom_free

end module mesh_geom
