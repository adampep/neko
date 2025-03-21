! Copyright (c) 2024, The Neko Authors
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
!
!> Defines a factory subroutine for simulation components.
submodule (simulation_component) simulation_component_fctry
  use vorticity, only : vorticity_t
  use force_torque, only : force_torque_t
  use fluid_stats_simcomp, only : fluid_stats_simcomp_t
  use lambda2, only : lambda2_t
  use probes, only : probes_t
  use les_simcomp, only : les_simcomp_t
  use utils, only : concat_string_array, neko_error
  use field_writer, only : field_writer_t
  use weak_grad, only : weak_grad_t
  use derivative, only : derivative_t
  use spectral_error, only: spectral_error_t
  use utils, only : neko_type_error

  ! List of all possible types created by the factory routine
  character(len=20) :: SIMCOMPS_KNOWN_TYPES(9) = [character(len=20) :: &
       "vorticity", &
       "lambda2", &
       "probes", &
       "les_model", &
       "field_writer", &
       "fluid_stats", &
       "weak_grad", &
       "force_torque", &
       "spectral_error"]

contains

  !> Simulation component factory. Both constructs and initializes the object.
  !! @param object The object to be created and initialized.
  !! @param json JSON object initializing the simulation component.
  !! @param case The simulation case.
  module subroutine simulation_component_factory(object, json, case)
    class(simulation_component_t), allocatable, intent(inout) :: object
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=:), allocatable :: type_name
    character(len=:), allocatable :: type_string
    logical :: is_user

    ! Check if this is a user-defined component
    call json_get_or_default(json, "is_user", is_user, .false.)
    if (is_user) return

    ! Get the type name
    call json_get(json, "type", type_name)

    ! Allocate
    call simulation_component_allocator(object, type_name)

    ! Initialize
    call object%init(json, case)

  end subroutine simulation_component_factory

  !> Simulation component allocator.
  !! @param object The object to be allocated.
  !! @param type_name The name of the simcomp type.
  module subroutine simulation_component_allocator(object, type_name)
    class(simulation_component_t), allocatable, intent(inout) :: object
    character(len=*), intent(in):: type_name

    select case (trim(type_name))
    case ("vorticity")
       allocate(vorticity_t::object)
    case ("lambda2")
       allocate(lambda2_t::object)
    case ("probes")
       allocate(probes_t::object)
    case ("les_model")
       allocate(les_simcomp_t::object)
    case ("field_writer")
       allocate(field_writer_t::object)
    case ("weak_grad")
       allocate(weak_grad_t::object)
    case ("derivative")
       allocate(derivative_t::object)
    case ("force_torque")
       allocate(force_torque_t::object)
    case ("fluid_stats")
       allocate(fluid_stats_simcomp_t::object)
    case ("spectral_error")
       allocate(spectral_error_t::object)
    case default
       call neko_type_error("simulation component", trim(type_name), &
            SIMCOMPS_KNOWN_TYPES)
    end select

  end subroutine simulation_component_allocator

end submodule simulation_component_fctry
