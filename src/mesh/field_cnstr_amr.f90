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
module field_cnstr_amr
  use num_types
  use logger
  use mesh_cnstr_amr
  use mesh
  use space
  use speclib
  use math
  use mxm_wrapper
  implicit none
  private

  public :: field_cnstr_amr_t

  !> Base type for a single field reconstruction during refinement/coarsening
  type :: field_cnstr_amr_t
     ! mesh pointer
     type(mesh_t), pointer :: msh => null()
     ! functional space pointer
     type(space_t), pointer :: Xh => null()
     ! reconstruction transfer data
     type(mesh_reconstruct_transfer_t), pointer :: rcn_trs => null()

     ! interpolation operators
     ! coarse to fine
     real(dp), allocatable, dimension(:,:,:) :: x_cr_to_fn, x_cr_to_fnT
     real(dp), allocatable, dimension(:,:,:) :: y_cr_to_fn, y_cr_to_fnT
     real(dp), allocatable, dimension(:,:,:) :: z_cr_to_fn, z_cr_to_fnT
     ! fine to coarse
     real(dp), allocatable, dimension(:,:,:) :: x_fn_to_cr, x_fn_to_crT
     real(dp), allocatable, dimension(:,:,:) :: y_fn_to_cr, y_fn_to_crT
     real(dp), allocatable, dimension(:,:,:) :: z_fn_to_cr, z_fn_to_crT
     ! point multiplicity; use for coarsening after children face summation
     ! for now I assume lx=ly=lz
     ! in general there should be 3 face and edge arrays
     ! element
     real(dp), allocatable, dimension(:,:,:) :: el_mult
     ! face
     real(dp), allocatable, dimension(:,:) :: fc_mult
     ! edge
     real(dp), allocatable, dimension(:) :: ed_mult

     ! work space
     real(dp), allocatable, dimension(:,:,:,:) :: tmp, ftmp

   contains
     procedure, pass(this) :: free => field_cnstr_amr_free
     procedure, pass(this) :: ftmp_get => fld_allocate_ftmp
     procedure, pass(this) :: ftmp_free => fld_deallocate_ftmp
     procedure, pass(this) :: refine => fld_rcn_refine_coarsen_single
  end type field_cnstr_amr_t

  interface field_cnstr_amr_t
     module procedure field_cnstr_amr_init
  end interface field_cnstr_amr_t

contains

  !> Initialise field reconstruction data
  !! @param[in]   msh      mesh information
  !! @param[in]   Xh       spectral space info
  !! @param[in]   rcn_trs  transfer data for field reconstruction
  !! @return  fld_rcn
  function field_cnstr_amr_init(msh, Xh, rcn_trs) result(fld_rcn)
    ! argument list
    type(mesh_t), target, intent(in) :: msh
    type(space_t), target, intent(in) :: Xh
    type( mesh_reconstruct_transfer_t), target, intent(in) :: rcn_trs
    ! result
    type(field_cnstr_amr_t) :: fld_rcn
    ! local variables
    integer(i4) :: il, jl, kl, nt2
    real(dp), allocatable, dimension(:) :: tmpl

    call fld_rcn%free()

    ! for now GLL only
    if (Xh%t /= 1) call neko_error('Only GLL space supported for AMR.')

    ! which mesh, space and transfer data this workflow relies on
    fld_rcn%msh => msh
    fld_rcn%Xh => Xh
    fld_rcn%rcn_trs => rcn_trs

    ! work array
    allocate(tmpl(max(Xh%lx, Xh%ly, Xh%lz)))

    ! interpolation operators
    ! x-direction
    allocate(fld_rcn%x_cr_to_fn(Xh%lx, Xh%lx, 2), source = 0.0_dp)
    allocate(fld_rcn%x_cr_to_fnT, source = fld_rcn%x_cr_to_fn)
    allocate(fld_rcn%x_fn_to_cr, source = fld_rcn%x_cr_to_fn)
    allocate(fld_rcn%x_fn_to_crT, source = fld_rcn%x_cr_to_fn)

    ! coarse -> fine
    ! negative
    do il = 1, Xh%lx
       tmpl(il) = 0.5_dp*(Xh%zg(il,1) - 1.0_dp)
    end do
    call igllm(fld_rcn%x_cr_to_fn, fld_rcn%x_cr_to_fnT, Xh%zg(:, 1), &
         &tmpl, Xh%lx, Xh%lx, Xh%lx, Xh%lx)
    ! positive; we use symmetry
    do jl = 1, Xh%lx
       do il = 1, Xh%lx
          fld_rcn%x_cr_to_fn(Xh%lx-il+1, Xh%lx-jl+1, 2) =&
               & fld_rcn%x_cr_to_fn(il, jl, 1)
          fld_rcn%x_cr_to_fnT(Xh%lx-il+1, Xh%lx-jl+1, 2) =&
               & fld_rcn%x_cr_to_fnT(il, jl, 1)
       end do
    end do
    ! fine -> coarse
    ! negative
    nt2 = Xh%lx/2 + mod(Xh%lx, 2)
    do il=1, nt2
       tmpl(il) = 2.0_dp*Xh%zg(il, 1) + 1.0_dp
    end do
    call igllm(fld_rcn%x_fn_to_cr, fld_rcn%x_fn_to_crT, Xh%zg(:, 1), &
         &tmpl, Xh%lx, nt2, Xh%lx, Xh%lx)
    ! positive; we use symmetry
    do jl = 1, Xh%lx
       do il = 1, nt2
          fld_rcn%x_fn_to_cr(Xh%lx-il+1, Xh%lx-jl+1, 2) =&
               & fld_rcn%x_fn_to_cr(il, jl, 1)
          fld_rcn%x_fn_to_crT(Xh%lx-il+1, Xh%lx-jl+1, 2) =&
               & fld_rcn%x_fn_to_crT(il, jl, 1)
       end do
    end do

    ! y-direction
    allocate(fld_rcn%y_cr_to_fn(Xh%ly, Xh%ly, 2), source = 0.0_dp)
    allocate(fld_rcn%y_cr_to_fnT, source = fld_rcn%y_cr_to_fn)
    allocate(fld_rcn%y_fn_to_cr, source = fld_rcn%y_cr_to_fn)
    allocate(fld_rcn%y_fn_to_crT, source = fld_rcn%y_cr_to_fn)

    ! coarse -> fine
    ! negative
    do il = 1, Xh%ly
       tmpl(il) = 0.5_dp*(Xh%zg(il, 2) - 1.0_dp)
    end do
    call igllm(fld_rcn%y_cr_to_fn, fld_rcn%y_cr_to_fnT, Xh%zg(:, 2), &
         &tmpl, Xh%ly, Xh%ly, Xh%ly, Xh%ly)
    ! positive; we use symmetry
    do jl = 1, Xh%ly
       do il = 1, Xh%ly
          fld_rcn%y_cr_to_fn(Xh%ly-il+1, Xh%ly-jl+1, 2) =&
               & fld_rcn%y_cr_to_fn(il, jl, 1)
          fld_rcn%y_cr_to_fnT(Xh%ly-il+1, Xh%ly-jl+1, 2) =&
               & fld_rcn%y_cr_to_fnT(il, jl, 1)
       end do
    end do
    ! fine -> coarse
    ! negative
    nt2 = Xh%ly/2 + mod(Xh%ly, 2)
    do il = 1, nt2
       tmpl(il) = 2.0_dp*Xh%zg(il, 2) + 1.0_dp
    end do
    call igllm(fld_rcn%y_fn_to_cr, fld_rcn%y_fn_to_crT, Xh%zg(:, 2), &
         &tmpl, Xh%ly, nt2, Xh%ly, Xh%ly)
    ! positive; we use symmetry
    do jl = 1, Xh%ly
       do il = 1, nt2
          fld_rcn%y_fn_to_cr(Xh%ly-il+1, Xh%ly-jl+1, 2) =&
               & fld_rcn%y_fn_to_cr(il, jl, 1)
          fld_rcn%y_fn_to_crT(Xh%ly-il+1, Xh%ly-jl+1, 2) =&
               & fld_rcn%y_fn_to_crT(il, jl, 1)
       end do
    end do

    ! z-direction
    if (msh%gdim == 3) then
       allocate(fld_rcn%z_cr_to_fn(Xh%lz, Xh%lz, 2), source = 0.0_dp)
       allocate(fld_rcn%z_cr_to_fnT, source = fld_rcn%z_cr_to_fn)
       allocate(fld_rcn%z_fn_to_cr, source = fld_rcn%z_cr_to_fn)
       allocate(fld_rcn%z_fn_to_crT, source = fld_rcn%z_cr_to_fn)

       ! coarse -> fine
       ! negative
       do il = 1, Xh%lz
          tmpl(il) = 0.5_dp*(Xh%zg(il, 3) - 1.0_dp)
       end do
       call igllm(fld_rcn%z_cr_to_fn, fld_rcn%z_cr_to_fnT, Xh%zg(:, 3), &
            &tmpl, Xh%lz, Xh%lz, Xh%lz, Xh%lz)
       ! positive; we use symmetry
       do jl = 1, Xh%lz
          do il = 1, Xh%lz
             fld_rcn%z_cr_to_fn(Xh%lz-il+1, Xh%lz-jl+1, 2) =&
                  & fld_rcn%z_cr_to_fn(il, jl, 1)
             fld_rcn%z_cr_to_fnT(Xh%lz-il+1, Xh%lz-jl+1, 2) =&
                  & fld_rcn%z_cr_to_fnT(il, jl, 1)
          end do
       end do
       ! fine -> coarse
       ! negative
       nt2 = Xh%lz/2 + mod(Xh%lz, 2)
       do il = 1, nt2
          tmpl(il) = 2.0_dp*Xh%zg(il, 3) + 1.0_dp
       end do
       call igllm(fld_rcn%z_fn_to_cr, fld_rcn%z_fn_to_crT, Xh%zg(:, 3), &
            &tmpl, Xh%lz, nt2, Xh%lz, Xh%lz)
       ! positive; we use symmetry
       do jl = 1, Xh%lz
          do il = 1, nt2
             fld_rcn%z_fn_to_cr(Xh%lz-il+1, Xh%lz-jl+1, 2) =&
                  & fld_rcn%z_fn_to_cr(il, jl, 1)
             fld_rcn%z_fn_to_crT(Xh%lz-il+1, Xh%lz-jl+1, 2) =&
                  & fld_rcn%z_fn_to_crT(il, jl, 1)
          end do
       end do
    end if

    ! multiplicity arrays
    if ((Xh%lx /= Xh%ly).or.((Xh%lz /= 1).and.(Xh%lz /= Xh%lx))) &
         & call neko_error('AMR does not support lx /= ly /= lz (for lz /= 1)')
    allocate(fld_rcn%el_mult(Xh%lx, Xh%ly, Xh%lz), source = 0.0_dp)
    allocate(fld_rcn%fc_mult(Xh%lx, Xh%lx), source = 0.0_dp)
    allocate(fld_rcn%ed_mult(Xh%lx), source = 0.0_dp)
    ! X
    if (mod(Xh%lx, 2) == 1) then
       il = Xh%lx/2 + 1
       do kl = 1, Xh%lz
          do jl = 1, Xh%ly
             fld_rcn%el_mult(il, jl, kl) = fld_rcn%el_mult(il, jl, kl) + 1.0_dp
          end do
       end do
    end if
    ! Y
    if (mod(Xh%ly, 2) == 1) then
       jl = Xh%ly/2 + 1
       do kl = 1, Xh%lz
          do il = 1, Xh%lx
             fld_rcn%el_mult(il, jl, kl) = fld_rcn%el_mult(il, jl, kl) + 1.0_dp
          end do
       end do
       if (mod(Xh%lx, 2) == 1) then
          il = Xh%lx/2 + 1
          do kl = 1, Xh%lz
             fld_rcn%el_mult(il, jl, kl) = fld_rcn%el_mult(il, jl, kl) + 1.0_dp
          end do
       end if
    end if

    ! Z
    if (msh%gdim == 3) then
       if (mod(Xh%lz, 2) == 1) then
          kl = Xh%lz/2 + 1
          do jl = 1, Xh%ly
             do il = 1, Xh%lx
                fld_rcn%el_mult(il, jl, kl) =&
                     & fld_rcn%el_mult(il, jl, kl) + 1.0_dp
             end do
          end do
          if (mod(Xh%lx, 2) == 1) then
             il = Xh%lx/2 + 1
             do jl = 1, XH%ly
                fld_rcn%el_mult(il, jl, kl) =&
                     & fld_rcn%el_mult(il, jl, kl) + 1.0_dp
             end do
          end if
          if (mod(Xh%ly, 2) == 1) then
             jl = Xh%ly/2 + 1
             do il = 1, Xh%lx
                fld_rcn%el_mult(il, jl, kl) =&
                     & fld_rcn%el_mult(il, jl, kl) + 1.0_dp
             end do
          end if
          if ((mod(Xh%lx, 2) == 1).and.(mod(Xh%ly,2 ) == 1)) then
             il = Xh%lx/2 + 1
             jl = Xh%ly/2 + 1
             fld_rcn%el_mult(il, jl, kl) = fld_rcn%el_mult(il, jl, kl) + 1.0_dp
          end if
       end if
    end if
    ! calculate inverse
    nt2 = Xh%lx*Xh%ly*Xh%lz
    call invcol1(fld_rcn%el_mult, nt2)

    ! to get proper J-1 on faces and edges for fast diagonalisation method
    ! I assume here LX1=LY1=LZ1, so only one array is needed
    !call ftovecl(fld_rcn%fc_mult,fld_rcn%el_mult,1,Xh%lx,Xh%ly,Xh%lz)
    !call etovec(fld_rcn%ed_mult,1,fld_rcn%el_mult,Xh%lx,Xh%ly,Xh%lz)
    ! THIS SHOLD BE CHECKED, BUT SHOULD BE FINE
    fld_rcn%fc_mult(:,:) = fld_rcn%el_mult(:,:, 1)
    fld_rcn%ed_mult(:) = fld_rcn%el_mult(:, 1, 1)

    ! work space
    allocate(fld_rcn%tmp(Xh%lx, Xh%ly, Xh%lz, 3))

    deallocate(tmpl)

    return
  end function field_cnstr_amr_init

  !> Free the type
  subroutine field_cnstr_amr_free(this)
    ! argument list
    class(field_cnstr_amr_t), intent(inout) :: this

    ! Clean pointers
    NULLIFY(this%msh, this%Xh, this%rcn_trs)

    ! Free memory
    if (allocated(this%x_cr_to_fn)) deallocate(this%x_cr_to_fn)
    if (allocated(this%x_cr_to_fnT)) deallocate(this%x_cr_to_fnT)
    if (allocated(this%x_fn_to_cr)) deallocate(this%x_fn_to_cr)
    if (allocated(this%x_fn_to_crT)) deallocate(this%x_fn_to_crT)

    if (allocated(this%y_cr_to_fn)) deallocate(this%y_cr_to_fn)
    if (allocated(this%y_cr_to_fnT)) deallocate(this%y_cr_to_fnT)
    if (allocated(this%y_fn_to_cr)) deallocate(this%y_fn_to_cr)
    if (allocated(this%y_fn_to_crT)) deallocate(this%y_fn_to_crT)

    if (allocated(this%z_cr_to_fn)) deallocate(this%z_cr_to_fn)
    if (allocated(this%z_cr_to_fnT)) deallocate(this%z_cr_to_fnT)
    if (allocated(this%z_fn_to_cr)) deallocate(this%z_fn_to_cr)
    if (allocated(this%z_fn_to_crT)) deallocate(this%z_fn_to_crT)

    if (allocated(this%el_mult)) deallocate(this%el_mult)
    if (allocated(this%fc_mult)) deallocate(this%fc_mult)
    if (allocated(this%ed_mult)) deallocate(this%ed_mult)

    if (allocated(this%tmp)) deallocate(this%tmp)
    if (allocated(this%ftmp)) deallocate(this%ftmp)

    return
  end subroutine field_cnstr_amr_free

  !> Allocate temporary array to store filed vector
  subroutine fld_allocate_ftmp(this)
    ! argument list
    class(field_cnstr_amr_t), intent(inout) :: this

    ! this requires some calculation
    !allocate(this%ftmp(,,,))
    return
  end subroutine fld_allocate_ftmp

  !> Allocate temporary array to store filed vector
  subroutine fld_deallocate_ftmp(this)
    ! argument list
    class(field_cnstr_amr_t), intent(inout) :: this

    if (allocated(this%ftmp)) deallocate(this%ftmp)
    return
  end subroutine fld_deallocate_ftmp

  !> @brief Perform refinement, transfer and coarsening of a single field
  !! @param[inout] vcf     refined/coarsened vector
  subroutine fld_rcn_refine_coarsen_single(this, vfc)
    ! argument list
    class(field_cnstr_amr_t), intent(inout) :: this
    real(dp), dimension(:,:,:,:), intent(inout) :: vfc
    ! local variable
    integer(i4) :: il
    integer(i4), dimension(4) :: lshape
    real(dp), allocatable, dimension(:,:,:,:) :: svfc

    ! here I should copy vector to vfc to ftmp and deallocate vfc

    ! local refinement
    if (this%rcn_trs%rfn_nr > 0) then
       if (mod(this%rcn_trs%rfn_nr, this%msh%npts) /= 0) &
            & call neko_error('Number of ref elem not multiply of nvert')
       call fld_rcn_refine_single(this, vfc)
    end if

    ! place for transfer/sorting
    ! JUST TESTING VERSION; THIS HAS TO BE REWRITTEN
    lshape(:) = shape(vfc)
    allocate(svfc(lshape(1),lshape(2),lshape(3),lshape(4)))
    do il = 1, this%msh%nelv
       if (this%rcn_trs%elgl_map(3, il) /= pe_rank) &
            & call neko_error('No data transfer between ranks for now')
       svfc(:,:,:,this%rcn_trs%elgl_map(1, il) - this%msh%offset_el) =&
            & vfc(:,:,:,il)
    end do
    vfc(:,:,:,:) = svfc(:,:,:,:)
    deallocate(svfc)

!!$    test : block
!!$      integer(i4) :: il
!!$      do il = 1, this%msh%nelv
!!$         write(*,*) 'TEST mapping', pe_rank, lshape(:), il,&
!!$              & this%rcn_trs%elgl_map(1, il) - this%msh%offset_el,&
!!$                   & this%rcn_trs%elgl_map(3, il)
!!$      end do
!!$    end block test


    ! local coarsening
    if (this%rcn_trs%crs_nr > 0) then
       call fld_rcn_coarsen_single(this, vfc)
    end if

    ! here vfc should be allocated and data should be moved from ftmp to vfc

    return
  end subroutine fld_rcn_refine_coarsen_single

  !> @brief Perform a single field refinement operation
  !! @param[inout] ref     refinement data
  !! @param[inout] vcf     refined vector
  subroutine fld_rcn_refine_single(ref, vfc)
    ! argument list
    type(field_cnstr_amr_t), intent(inout) :: ref
    real(dp), dimension(:,:,:,:), intent(inout) :: vfc
    ! local variables
    integer(i4) :: il, jl, itmp
    integer(i4), dimension(3) :: ch_pos

    ! I assume el_lst(1) gives position of the coarse block
    ! and final ch_pos() = 1,1,1
    ! loop over refined elements
    do il= 1, ref%rcn_trs%rfn_nr, ref%msh%npts ! should be changed to number of childeren?
       ! copy coarse element to a temporary array
       ref%tmp(:,:,:, 3) = vfc(:,:,:, ref%rcn_trs%elgl_rfn(3, il))
!!$       test1 : block
!!$         integer :: kl,ll
!!$         if (pe_rank == 0) then
!!$         do kl = 1, ref%Xh%lx
!!$            do ll = 1, ref%Xh%lx
!!$               write(*,*) "OLD",ref%rcn_trs%elgl_rfn(3, il),ll,kl,
!!$                    ref%tmp(:,ll,kl, 3)
!!$            end do
!!$            write(*,*) ' '
!!$         end do
!!$         write(*,*) '===========================================',pe_rank, &
!!$              & ref%rcn_trs%elgl_rfn(3, il)
!!$         end if
!!$       end block test1
       ! loop over all the children
       do jl= 1, ref%msh%npts
          ! get child position
          ! new position in the array
          itmp = ref%rcn_trs%elgl_rfn(3, il + jl - 1)
          ch_pos(3) = (jl-1)/4 + 1 ! z position
          ch_pos(2) = mod((jl - 1)/2, 2) + 1 ! y position
          ch_pos(1) = mod(jl - 1, 2) +1 ! x position
          ! refine
          call fld_rcn_map_ctof(ref, ch_pos, ref%tmp(:,:,:,3), vfc(:,:,:,itmp))
!!$          test2 : block
!!$            integer :: kl,ll
!!$            if (pe_rank == 0) then
!!$            do kl = 1, ref%Xh%lx
!!$               do ll = 1, ref%Xh%lx
!!$                  write(*,*) "NEW",itmp,ll,kl,vfc(:,ll,kl, itmp)
!!$               end do
!!$               write(*,*) ' '
!!$            end do
!!$            write(*,*) '==============================================',&
!!$            & pe_rank, itmp,ch_pos(:)
!!$            endif
!!$          end block test2
       end do
    end do

    return
  end subroutine fld_rcn_refine_single

  !> @brief Perform a single field coarsening operation
  !! @param[inout] ref     refinement data
  !! @param[inout] vcf     coarsened vector
  subroutine fld_rcn_coarsen_single(ref, vfc)
    ! argument list
    type(field_cnstr_amr_t), intent(inout) :: ref
    real(dp), dimension(:,:,:,:), intent(inout) :: vfc
    ! local variables
    integer(i4) :: il, jl, itmp
    integer(i4), dimension(3) :: ch_pos

    ! I assume el_lst(1) gives position of the coarse block
    ! and final ch_pos() = 1,1,1
    ! loop over coarsened elements
    do il= 1, ref%rcn_trs%crs_nr
       ! loop over all the children
       do jl= 1, ref%msh%npts
          ! get child position
          itmp = ref%rcn_trs%elgl_crs(2, jl, il) ! new position in the array
          ch_pos(3) = (jl-1)/4 + 1 ! z position
          ch_pos(2) = mod((jl - 1)/2, 2) + 1 ! y position
          ch_pos(1) = mod(jl - 1, 2) +1 ! x position
          ! coarsen
          call fld_rcn_map_ftoc(ref, ch_pos, vfc(:,:,:,itmp), ref%tmp(:,:,:,3))
          ! sum contributions
          if (jl == 1) then
             vfc(:,:,:,itmp) = ref%tmp(:,:,:,3)
          else
             vfc(:,:,:,itmp) =  vfc(:,:,:,itmp) + ref%tmp(:,:,:,3)
          end if
       end do
       vfc(:,:,:,itmp) =  vfc(:,:,:,itmp)*ref%el_mult(:,:,:)
    end do

    return
  end subroutine fld_rcn_coarsen_single

  !> @brief Map a single coarse element to a fine one
  !! @param[inout] ref     refinement data
  !! @param[in]    ch_pos  child position
  !! @param[in]    vc      coarse element vector
  !! @param[out]   vf      fine element vector
  subroutine fld_rcn_map_ctof(ref, ch_pos, vc, vf)
    ! argument list
    type(field_cnstr_amr_t), intent(inout) :: ref
    integer(i4), dimension(3), intent(in) :: ch_pos
    real(dp), dimension(:,:,:), intent(in) :: vc
    real(dp), dimension(:,:,:), intent(out) :: vf
    ! local variables
    integer(i4) :: iz

    if (ref%msh%gdim == 3) then ! 3D
       call mxm(ref%x_cr_to_fn(:,:, ch_pos(1)), ref%Xh%lx, vc, &
            & ref%Xh%lx, ref%tmp(:,:,:, 1), ref%Xh%lyz)
       do iz = 1, ref%Xh%lz
          call mxm(ref%tmp(:,:, iz, 1), ref%Xh%lx,&
               & ref%y_cr_to_fnT(:,:, ch_pos(2)),&
               & ref%Xh%ly, ref%tmp(:,:, iz, 2), ref%Xh%ly)
       end do
       call mxm(ref%tmp(:,:,:, 2), ref%Xh%lxy,&
            & ref%z_cr_to_fnT(:,:, ch_pos(3)),&
            & ref%Xh%lz, vf, ref%Xh%lz)
    else ! 2D
       call mxm(ref%x_cr_to_fn(:,:, ch_pos(1)), ref%Xh%lx, vc, ref%Xh%lz,&
            & ref%tmp(:,:,1,1), ref%Xh%lyz)
       call mxm(ref%tmp(:,:,1,1), ref%Xh%lx, ref%y_cr_to_fnT(:,:, ch_pos(2)),&
            & ref%Xh%ly, vf, ref%Xh%ly)
    end if

    return
  end subroutine fld_rcn_map_ctof

  !> @brief Map a single fine element to a coarse one
  !! @param[inout] ref     refinement data
  !! @param[in]    ch_pos  child position
  !! @param[in]    vf      fine element vector
  !! @param[out]   vc      coarse element vector
  subroutine fld_rcn_map_ftoc(ref, ch_pos, vf, vc)
    ! argument list
    type(field_cnstr_amr_t), intent(inout) :: ref
    integer(i4), dimension(3), intent(in) :: ch_pos
    real(dp), dimension(:,:,:), intent(in) :: vf
    real(dp), dimension(:,:,:), intent(out) :: vc
    ! local variables
    integer(i4) :: iz

    if (ref%msh%gdim == 3) then ! 3D
       call mxm(ref%x_fn_to_cr(:,:, ch_pos(1)), ref%Xh%lx, vf, &
            & ref%Xh%lx, ref%tmp(:,:,:, 1), ref%Xh%lyz)
       do iz = 1, ref%Xh%lz
          call mxm(ref%tmp(:,:, iz, 1), ref%Xh%lx,&
               & ref%y_fn_to_crT(:,:, ch_pos(2)),&
               & ref%Xh%ly, ref%tmp(:,:, iz, 2), ref%Xh%ly)
       end do
       call mxm(ref%tmp(:,:,:, 2), ref%Xh%lxy,&
            & ref%z_fn_to_crT(:,:, ch_pos(3)),&
            & ref%Xh%lz, vc, ref%Xh%lz)
    else ! 2D
       call mxm(ref%x_fn_to_cr(:,:, ch_pos(1)), ref%Xh%lx, vf, ref%Xh%lz,&
            & ref%tmp(:,:,1,1), ref%Xh%lyz)
       call mxm(ref%tmp(:,:,1,1), ref%Xh%lx, ref%y_fn_to_crT(:,:, ch_pos(2)),&
            & ref%Xh%ly, vc, ref%Xh%ly)
    end if

    return
  end subroutine fld_rcn_map_ftoc

end module field_cnstr_amr
