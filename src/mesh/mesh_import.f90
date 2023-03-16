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
module mesh_import
  implicit none
  private
  
  !> Base type for mesh data import from mesh managers
  type, public, abstract :: mesh_import_t 
     logical, private :: initialised_ = .false. !< Mesh manager initialisation flag
   contains
     procedure, pass(this) :: initialised => manager_init
     procedure, pass(this) :: set_init => manager_set_init
     procedure, pass(this) :: free => manager_free
  end type mesh_import_t

contains

  !> Return the mesh manager initialisation flag
  pure function manager_init(this) result(initialised)
    class(mesh_import_t), intent(in) :: this
    logical :: initialised
    initialised = this%initialised_
  end function manager_init

  !> Set the mesh manager initialisation flag
  subroutine manager_set_init(this, initialised)
    class(mesh_import_t), intent(inout) :: this
    logical, intent(in) :: initialised
    this%initialised_ = initialised
  end subroutine manager_set_init

  !> Set the mesh manager initialisation flag to false
  subroutine manager_free(this)
    class(mesh_import_t), intent(inout) :: this
    this%initialised_ = .false.
  end subroutine manager_free

end module mesh_import
