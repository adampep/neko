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
module mesh_cnstr
  use num_types
  use mesh
  implicit none
  private

  public :: mesh_manager_transfer_t, mesh_cnstr_t

  !> Type for data transfer between mesh manager and neko
  type mesh_manager_transfer_t
     integer(i4) :: mm_nel ! local element number on mesh manager side
     integer(i4) :: nk_nel ! local element number on neko side
     ! mesh manager <=> neko element distribution mapping
     ! (global element number, process id)
     integer(i4), allocatable, dimension(:, :) :: elmap_mm2nk, elmap_nk2mm
   contains
     procedure, public, pass(this) :: free => manager_transfer_free
  end type mesh_manager_transfer_t

  !> Base type for mesh construction
  type, abstract :: mesh_cnstr_t
     ! element mapping between mesh manager and neko
     type(mesh_manager_transfer_t) :: msh_trs
   contains
     procedure(mesh_cnstr_extract), pass(this), deferred :: mesh_get
     procedure(mesh_cnstr_free), pass(this), deferred :: free
  end type mesh_cnstr_t

  abstract interface
     subroutine mesh_cnstr_extract(this, msh)
       import :: mesh_cnstr_t
       import :: mesh_t
       class(mesh_cnstr_t), intent(inout) :: this
       type(mesh_t), intent(inout) :: msh
     end subroutine mesh_cnstr_extract
  end interface

  abstract interface
     subroutine mesh_cnstr_free(this)
       import :: mesh_cnstr_t
       class(mesh_cnstr_t), intent(inout) :: this
     end subroutine mesh_cnstr_free
  end interface

contains

  subroutine manager_transfer_free(this)
    ! argument list
    class(mesh_manager_transfer_t), intent(inout) :: this

    ! Deallocate arrays
    if (allocated(this%elmap_mm2nk)) deallocate(this%elmap_mm2nk)
    if (allocated(this%elmap_nk2mm)) deallocate(this%elmap_nk2mm)

    ! Reset registers
    this%mm_nel = 0
    this%nk_nel = 0

    return
  end subroutine manager_transfer_free

end module mesh_cnstr
