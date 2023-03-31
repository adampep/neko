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
module mesh_manager
  use num_types
  use logger
  use utils
  use mesh_cnstr
  use p4est
  implicit none
  private

  public :: mesh_manager_t

  !> Type for mesh manager
  type mesh_manager_t
     logical, private :: initialised_ = .false. ! mesh manager
                                                ! initialisation flag
     class(mesh_cnstr_t), allocatable :: msh_imp ! mesh manager's imported mesh
   contains
     procedure, pass(this) :: initialised => manager_init
     procedure, pass(this) :: set_init => manager_set_init
     procedure, pass(this) :: free => manager_free
  end type mesh_manager_t

  interface mesh_manager_t
     module procedure mesh_manager_init
  end interface mesh_manager_t

contains

  !> Initialise mesh manager
  function mesh_manager_init(mesh_file, log_threshold) result(this)
    character(len=*), intent(in) :: mesh_file
    integer(i4), intent(in), optional :: log_threshold
    type(mesh_manager_t), target :: this
    character(len=80) :: suffix

    call this%free()

    call filename_suffix(mesh_file, suffix)

    if (suffix == "p4est") then
       call neko_log%message('Initialising mesh manager p4est')
       allocate(p4_mesh_cnstr_t :: this%msh_imp)

       if (present(log_threshold)) then
          call p4_manager_init(mesh_file, log_threshold)
       else
          call p4_manager_init(mesh_file)
       end if
    else
       call neko_error('Unknown mesh manager')
    end if

    call this%set_init(.true.)

    return
  end function mesh_manager_init

  !> Return the mesh manager initialisation flag
  pure function manager_init(this) result(initialised)
    class(mesh_manager_t), intent(in) :: this
    logical :: initialised
    initialised = this%initialised_
    return
  end function manager_init

  !> Set the mesh manager initialisation flag
  subroutine manager_set_init(this, initialised)
    class(mesh_manager_t), intent(inout) :: this
    logical, intent(in) :: initialised
    this%initialised_ = initialised
    return
  end subroutine manager_set_init

  !> Finalise mesh manager and free memory
  subroutine manager_free(this, log_threshold)
    class(mesh_manager_t), intent(inout) :: this
    integer(i4), intent(in), optional :: log_threshold

    if (allocated(this%msh_imp)) then
       ! free memory
       call this%msh_imp%free()
       ! finalise mesh manager
       select type (imp => this%msh_imp)
       type is(p4_mesh_cnstr_t)
          call neko_log%message('Finalising mesh manager p4est')
          if (present(log_threshold)) then
             call p4_manager_free(log_threshold)
          else
             call p4_manager_free()
          end if
       end select
       deallocate(this%msh_imp)
    end if

    this%initialised_ = .false.
    return
  end subroutine manager_free

end module mesh_manager
