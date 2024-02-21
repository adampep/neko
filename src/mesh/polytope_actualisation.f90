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
 !> Abstract type for polytope actualisation class for elements building blocks
module polytope_actualisation
  use num_types, only : i4
  use polytope, only : polytope_t
  use polytope_aligned, only : polytope_aligned_t
  implicit none
  private

  public :: polytope_actualisation_t

  !> Base type for a polytope actualisation.
  !! @details This is an abstract type building on topology and alignment data
  !! with nonconforming and interpolation information. Note that vertices can be
  !! hanging but have no interpolation operation. This type corresponds to
  !! the realisation of an abstract objects as parts of an existing higher
  !! dimension elements.
  type, extends(polytope_aligned_t), abstract :: polytope_actualisation_t
     !> Is there any interpolation operator active (excluding identity)
     logical, private :: ifinterpolation_ = .false.
     !> Hanging information
     integer(i4), private :: hanging_ = 0
     !> position in the higher dimension polytope
     integer(i4), private :: position_ = -1
   contains
     !> Free aligned polytope and interpolation data
     procedure, pass(this) :: free => polytope_actualisation_free
     !> Initialise general data
     procedure, pass(this) :: init_act => polytope_actualisation_init_act
     !> Return interpolation flag
     procedure, pass(this) :: ifintp => polytope_actualisation_ifintp_get
     !> Return hanging node information
     procedure, pass(this) :: hng => polytope_actualisation_hng_get
     !> Return position in the higher dimension object
     procedure, pass(this) :: pos => polytope_actualisation_pos_get
     !> Initialise a polytope actualisation
     procedure(polytope_actualisation_init), pass(this), deferred :: init
  end type polytope_actualisation_t

  !> Initialise a polytope with hanging information
  !! @parameter[in]   pltp   polytope
  !! @parameter[in]   algn   alignment information
  !! @parameter[in]   ifint  interpolation flag
  !! @parameter[in]   hng    hanging information
  !! @parameter[in]   pos    position in the higher order element
  abstract interface
     subroutine polytope_actualisation_init(this, pltp, algn, ifint, hng, pos)
       import i4
       import polytope_t
       import polytope_actualisation_t
       class(polytope_actualisation_t), intent(inout) :: this
       class(polytope_t), target, intent(in) :: pltp
       integer(i4), intent(in) :: algn, hng, pos
       logical, intent(in) :: ifint
     end subroutine polytope_actualisation_init
  end interface

contains

  !> Free aligned polytope and interpolation data
  subroutine polytope_actualisation_free(this)
    class(polytope_actualisation_t), intent(inout) :: this

    call this%free_base()

    this%ifinterpolation_ = .false.
    this%hanging_ = 0
    this%position_ = -1
  end subroutine polytope_actualisation_free

  !> Initialise general data
  !! @parameter[in]   pltp   polytope
  !! @parameter[in]   ifalgn if non identity alignment
  !! @parameter[in]   ifint  interpolation flag
  !! @parameter[in]   hng    hanging information
  !! @parameter[in]   pos    position in the higher order element
  subroutine polytope_actualisation_init_act(this, pltp, ifalgn, ifint, hng, &
       & pos)
    class(polytope_actualisation_t), intent(inout) :: this
    class(polytope_t), target, intent(in) :: pltp
    logical, intent(in)  :: ifalgn, ifint
    integer(i4), intent(in) :: hng, pos

    call this%init_base(pltp, ifalgn)
    this%ifinterpolation_ = ifint
    this%hanging_ = hng
    this%position_ = pos
  end subroutine polytope_actualisation_init_act

  !> @brief Get interpolation flag
  !! @return   ifintp
  pure function polytope_actualisation_ifintp_get(this) result(ifintp)
    class(polytope_actualisation_t), intent(in) :: this
    logical :: ifintp
    ifintp = this%ifinterpolation_
  end function polytope_actualisation_ifintp_get

  !> @brief Get hanging information
  !! @return   hng
  pure function polytope_actualisation_hng_get(this) result(hng)
    class(polytope_actualisation_t), intent(in) :: this
    integer(i4) :: hng
    hng = this%hanging_
  end function polytope_actualisation_hng_get

  !> @brief Get position information
  !! @return   pos
  pure function polytope_actualisation_pos_get(this) result(pos)
    class(polytope_actualisation_t), intent(in) :: this
    integer(i4) :: pos
    pos = this%position_
  end function polytope_actualisation_pos_get

end module polytope_actualisation