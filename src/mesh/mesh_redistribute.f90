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
module mesh_redistribute
  use num_types
  implicit none
  private

  public :: mesh_manager_transfer_t, mesh_reconstruct_transfer_t

  !> Type for data transfer between mesh manager and neko
  type mesh_manager_transfer_t
     ! mesh manager <=> neko element distribution mapping
     ! (global element number, process id)
     integer(i4), allocatable, dimension(:, :) :: elmap_p2n, elmap_n2p
   contains
     procedure, public, pass(this) :: free => manager_transfer_free
  end type mesh_manager_transfer_t

  !> Type for data transfer during mesh reconstruction
  type mesh_reconstruct_transfer_t
     ! arrays to store global element re-mapping to perform refinement on
     ! neko side
     ! map_nr - local number of unchanged elements
     ! rfn_nr - local number of refined elements
     ! crs_nr - local number of coarsened elements
     ! rfn_nr_a - local number of elements after refinement
     ! (including refined ones)
     ! crs_nr_s - local number of elements for coarsening (child owner);
     ! THIS IS MOST PROBABLY NOT NEEDED
     ! crs_nr_b - local number of elements before coarsening
     ! (parent owner; including one that would disappear)
     ! nelvo - old number of local elements
     integer(i4) :: map_nr, rfn_nr, crs_nr, rfn_nr_a, crs_nr_s, crs_nr_b, nelvo
     ! elgl_map - element number/process id mapping data for unchanged elements
     ! (old gl. num., old loc. num., old proc. id)
     ! elgl_rfn - element number/process id mapping data for refined elements
     ! (ch. gl. num., old p. gl. num., old p. loc. num., old p. proc. id, ch. pos.)
     integer(i4), allocatable, dimension(:, :) :: elgl_map, elgl_rfn
     ! elgl_crs - element number/process id mapping data for coarsened elements
     ! (new gl. num., old ch. gl. num., old ch. loc. num., old ch. proc. id)
     integer(i4), allocatable, dimension(:, :, :) :: elgl_crs
   contains
     procedure, public, pass(this) :: free => reconstruct_transfer_free
  end type mesh_reconstruct_transfer_t

contains

  subroutine manager_transfer_free(this)
    ! argument list
    class(mesh_manager_transfer_t), intent(inout) :: this

    ! Deallocate arrays
    if (allocated(this%elmap_p2n)) deallocate(this%elmap_p2n)
    if (allocated(this%elmap_n2p)) deallocate(this%elmap_n2p)

    return
  end subroutine manager_transfer_free

  subroutine reconstruct_transfer_free(this)
    ! argument list
    class(mesh_reconstruct_transfer_t), intent(inout) :: this

    ! Reset registers
    this%map_nr = 0
    this%rfn_nr = 0
    this%crs_nr = 0
    this%rfn_nr_a = 0
    this%crs_nr_s = 0
    this%crs_nr_b = 0
    this%nelvo = 0

    ! Deallocate arrays
    if (allocated(this%elgl_map)) deallocate(this%elgl_map)
    if (allocated(this%elgl_rfn)) deallocate(this%elgl_rfn)
    if (allocated(this%elgl_crs)) deallocate(this%elgl_crs)

    return
  end subroutine reconstruct_transfer_free

end module mesh_redistribute
