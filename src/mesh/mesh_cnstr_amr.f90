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
module mesh_cnstr_amr
  use num_types
  use mesh_redistribute
  use mesh_cnstr
  implicit none
  private

  public :: mesh_cnstr_amr_t

  !> Base type for mesh construction in the refinement/coarsening process
  type, extends(mesh_cnstr_t), abstract :: mesh_cnstr_amr_t
   contains
     procedure(mesh_cnstr_refine), pass(this), deferred :: refine
  end type mesh_cnstr_amr_t

  abstract interface
     subroutine mesh_cnstr_refine(this, ref_mark, el_gidx, msh_trs, level_max,&
       &ifmod, msh_rcn)
       import :: mesh_cnstr_amr_t
       import :: mesh_manager_transfer_t
       import :: mesh_reconstruct_transfer_t
       import :: i4
       class(mesh_cnstr_amr_t), intent(inout) :: this
       integer(i4), dimension(:), intent(in) :: ref_mark, el_gidx
       type(mesh_manager_transfer_t), intent(in) :: msh_trs
       integer(i4), intent(in) :: level_max
       logical, intent(out) :: ifmod
       type(mesh_reconstruct_transfer_t), intent(out) :: msh_rcn
     end subroutine mesh_cnstr_refine
  end interface

contains

end module mesh_cnstr_amr
