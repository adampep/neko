! Copyright (c) 2023, The Neko Authors
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
!> Simulation components are objects that encapsulate functionality that can be
!! fit to a particular compute pattern.
!! @note
!! The canonical way to abbreviate simulation_component is simcomp.
module simulation_component
  use num_types, only : rp
  use json_module, only : json_file
  use case, only : case_t
  use time_based_controller, only : time_based_controller_t
  use json_utils, only : json_get_or_default, json_get
  implicit none
  private

  !> Base abstract class for simulation components.
  type, abstract, public :: simulation_component_t
     !> Pointer to the simulation case.
     type(case_t), pointer :: case
     !> Controller for when to run `preprocess`.
     type(time_based_controller_t) :: preprocess_controller
     !> Controller for when to run `compute`.
     type(time_based_controller_t) :: compute_controller
     !> Controller for when to do output.
     type(time_based_controller_t) :: output_controller
     !> The execution order, lowest excutes first.
     integer :: order
   contains
     !> Constructor for the simulation_component_t (base) class.
     procedure, pass(this) :: init_base => simulation_component_init_base
     !> Destructor for the simulation_component_t (base) class.
     procedure, pass(this) :: free_base => simulation_component_free_base
     !> Wrapper for calling `set_counter` for the time based controllers.
     !! Serves as the public interface.
     procedure, pass(this) :: restart => simulation_component_restart_wrapper
     !> Wrapper for calling `preprocess_` based on the `preprocess_controller`.
     !! Serves as the public interface.
     procedure, pass(this) :: preprocess => &
          simulation_component_preprocess_wrapper
     !> Wrapper for calling `compute_` based on the `compute_controller`.
     !! Serves as the public interface.
     procedure, pass(this) :: compute => simulation_component_compute_wrapper
     !> The common constructor using a JSON dictionary.
     procedure(simulation_component_init), pass(this), deferred :: init
     !> Destructor.
     procedure(simulation_component_free), pass(this), deferred :: free
     !> The preprocessing function to be executed during the run.
     procedure, pass(this) :: preprocess_
     !> The main function to be executed during the run.
     procedure, pass(this) :: compute_
     !> The restart function to be called upon restarting simulation
     procedure, pass(this) :: restart_
  end type simulation_component_t

  !> A helper type that is needed to have an array of polymorphic objects
  type, public :: simulation_component_wrapper_t
     class(simulation_component_t), allocatable :: simcomp
  end type simulation_component_wrapper_t


  abstract interface
     !> The common constructor using a JSON dictionary.
     !! @param json The JSON with properties.
     !! @param case The case_t object.
     subroutine simulation_component_init(this, json, case)
       import simulation_component_t, json_file, case_t
       class(simulation_component_t), intent(inout) :: this
       type(json_file), intent(inout) :: json
       class(case_t), intent(inout), target :: case
     end subroutine simulation_component_init
  end interface

  abstract interface
     !> Destructor.
     subroutine simulation_component_free(this)
       import simulation_component_t
       class(simulation_component_t), intent(inout) :: this
     end subroutine simulation_component_free
  end interface

  interface
     !> Simulation component factory.
     !! Both constructs and initializes the object.
     !! @param object The object to be created and initialized.
     !! @param json JSON object initializing the simulation component.
     !! @param case The simulation case.
     module subroutine simulation_component_factory(object, json, case)
       class(simulation_component_t), allocatable, intent(inout) :: object
       type(json_file), intent(inout) :: json
       class(case_t), intent(inout), target :: case
     end subroutine simulation_component_factory
  end interface

  interface
     !> Simulation component allocator.
     !! @param object The object to be allocated.
     !! @param type_name The name of the simcomp type.
     module subroutine simulation_component_allocator(object, type_name)
       class(simulation_component_t), allocatable, intent(inout) :: object
       character(len=*), intent(in):: type_name
     end subroutine simulation_component_allocator
  end interface

  public :: simulation_component_factory

contains
  !> Constructor for the `simulation_component_t` (base) class.
  subroutine simulation_component_init_base(this, json, case)
    class(simulation_component_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=:), allocatable :: preprocess_control, compute_control, output_control
    real(kind=rp) :: preprocess_value, compute_value, output_value
    integer :: order

    this%case => case

    ! We default to preprocess every time-step
    call json_get_or_default(json, "preprocess_control", preprocess_control, &
         "tsteps")
    call json_get_or_default(json, "preprocess_value", preprocess_value, 1.0_rp)

    ! We default to compute every time-step
    call json_get_or_default(json, "compute_control", compute_control, &
         "tsteps")
    call json_get_or_default(json, "compute_value", compute_value, 1.0_rp)

    if (compute_control .eq. "fluid_output") then
       call json_get(this%case%params, 'case.fluid.output_control', &
            compute_control)
       call json_get(this%case%params, 'case.fluid.output_value', &
            compute_value)
    end if

    ! We default to output whenever we execute
    call json_get_or_default(json, "output_control", output_control, &
         compute_control)
    call json_get_or_default(json, "output_value", output_value, &
         compute_value)

    if (output_control == "global") then
       call json_get(this%case%params, 'case.fluid.output_control', &
            output_control)
       call json_get(this%case%params, 'case.fluid.output_value', &
            output_value)
    end if

    call json_get_or_default(json, "order", order, -1)
    this%order = order

    call this%preprocess_controller%init(case%end_time, preprocess_control, &
         preprocess_value)
    call this%compute_controller%init(case%end_time, compute_control, &
         compute_value)
    call this%output_controller%init(case%end_time, output_control, &
         output_value)

  end subroutine simulation_component_init_base

  !> Destructor for the `simulation_component_t` (base) class.
  subroutine simulation_component_free_base(this)
    class(simulation_component_t), intent(inout) :: this

    nullify(this%case)
  end subroutine simulation_component_free_base

  !> Wrapper for calling `preprocess_` based on the `preprocess_controller`.
  !! Serves as the public interface.
  !! @param t The time value.
  !! @param tstep The current time-step
  subroutine simulation_component_preprocess_wrapper(this, t, tstep)
    class(simulation_component_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    if (this%preprocess_controller%check(t, tstep)) then
       call this%preprocess_(t, tstep)
       call this%preprocess_controller%register_execution()
    end if
  end subroutine simulation_component_preprocess_wrapper

  !> Wrapper for calling `compute_` based on the `compute_controller`.
  !! Serves as the public interface.
  !! @param t The time value.
  !! @param tstep The current time-step
  subroutine simulation_component_compute_wrapper(this, t, tstep)
    class(simulation_component_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    if (this%compute_controller%check(t, tstep)) then
       call this%compute_(t, tstep)
       call this%compute_controller%register_execution()
    end if
  end subroutine simulation_component_compute_wrapper

  !> Wrapper for calling `set_counter_` based for the controllers.
  !! Serves as the public interface.
  !! @param t The time value.
  subroutine simulation_component_restart_wrapper(this, t)
    class(simulation_component_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t

    call this%compute_controller%set_counter(t)
    call this%output_controller%set_counter(t)
    call this%restart_(t)

  end subroutine simulation_component_restart_wrapper

  !> Dummy restart function.
  !! @param t The time value.
  subroutine restart_(this, t)
    class(simulation_component_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t

    ! Do nothing
  end subroutine restart_

  !> Dummy preprocessing function.
  !! @param t The time value.
  !! @param tstep The current time-step
  subroutine preprocess_(this, t, tstep)
    class(simulation_component_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    ! Do nothing
  end subroutine preprocess_

  !> Dummy compute function.
  !! @param t The time value.
  !! @param tstep The current time-step
  subroutine compute_(this, t, tstep)
    class(simulation_component_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    ! Do nothing
  end subroutine compute_
end module simulation_component
