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
!> Main interface for h-type amr in neko
module amr
  use mpi_f08
  use num_types
  use comm
  use logger
  use parameters
  use mesh_cnstr_amr
  use mesh_manager
  use field_cnstr_amr
  use p4est
  use mesh
  use fluid_method
  use fluid_plan1
  use fluid_plan4
  use fluid_pnpn
  use user_intf
  implicit none

  private
  public :: amr_t

  ! type for fields reconstruction (refinement/transfer/coarsening)
  type amr_t
     logical, private :: ifamr_ = .false. ! is amr used

     ! mesh manager and imported data
     type(mesh_manager_t) :: msh_mngr

     ! field reconstruction data
     type(field_cnstr_amr_t) :: fld_rcn
   contains
     procedure, pass(this) :: ifamr => amr_is_used
     procedure, pass(this) :: get_mesh => amr_get_mesh
     procedure, pass(this) :: set_fctr => amr_init_field_cnstr
     procedure, pass(this) :: msh_refine => amr_msh_refine
     procedure, pass(this) :: fld_refine => amr_fld_refine
     procedure, pass(this) :: free => amr_free
  end type amr_t

  interface amr_t
     module procedure amr_init
  end interface amr_t

contains

  !> Initialise mesh manager
  !! @return amr
  function amr_init(mesh_file) result(amr)
    ! argument list
    character(len=*), intent(in) :: mesh_file
    ! result
    type(amr_t) :: amr

    call amr%free()

    ! set up mesh manager
    amr%msh_mngr = mesh_manager_t(mesh_file)

    ! I do not initialise field constructor here, as space Xh is not
    ! available yet

    amr%ifamr_ = .true.

    return
  end function amr_init

  !> Return amr usage flag
  pure function amr_is_used(this) result(used)
    class(amr_t), intent(in) :: this
    logical :: used
    used = this%ifamr_

    return
  end function amr_is_used

  !> Provide mesh from the mesh manger
  function amr_get_mesh(this) result(msh)
    ! argument list
    class(amr_t), intent(inout) :: this
    ! result
    type(mesh_t) :: msh

    call this%msh_mngr%msh_imp%mesh_get(msh)

    return
  end function amr_get_mesh

  !> Initialise field constructor
  subroutine amr_init_field_cnstr(this, msh, Xh)
    ! argument list
    class(amr_t), intent(inout) :: this
    type(mesh_t), intent(in) :: msh
    type(space_t), intent(in) :: Xh

    select type(msh_imp => this%msh_mngr%msh_imp)
    class is(mesh_cnstr_amr_t)
       this%fld_rcn = field_cnstr_amr_t(msh, Xh, msh_imp%rcn_trs)
    class default
       call neko_error('Invalid refinement class type')
    end select

    return
  end subroutine amr_init_field_cnstr

  !> Finalise mesh manager and clean memory
  subroutine amr_free(this)
    ! argument list
    class(amr_t), intent(inout) :: this

    ! Free types
    call this%fld_rcn%free()
    call this%msh_mngr%free()

    this%ifamr_ = .false.

    return
  end subroutine amr_free

  !> Perform mesh refinement
  !! @param[inout] msh      mesh information
  !! @param[in]    param    run-time parameters
  !! @param[in]    usr      user interface
  !! @param[in]    time     simulation time
  !! @param[in]    tstep    time step
  subroutine amr_msh_refine(this, msh, param, usr, time, tstep)
    ! argument list
    class(amr_t), intent(inout) :: this
    type(mesh_t), intent(inout) :: msh
    type(param_t), intent(in) :: param
    type(user_t), intent(in) :: usr
    real(kind=rp), intent(in) :: time
    integer, intent(in) :: tstep
    ! local variables
    integer(i4) :: level_max, il, jl, ielo, itmp, ips, ich
    integer(i4) :: ierr, gnum, goff
    logical :: ifmod
    integer(i4), dimension(5,8) :: itmpv
    integer(i4), allocatable, dimension(:) :: ref_mark, el_gidx
    integer(i4), allocatable, dimension(:,:) :: lel_map

    call neko_log%section("Mesh refinement")

    ! get refinement flag
    allocate(ref_mark(msh%nelv), source = 0)
    call usr%amr_refmark(ref_mark, time, tstep, msh, param)

    ! check if any refinement is marked
    itmp = minval(ref_mark)
    call MPI_Allreduce(itmp, il, 1, MPI_INTEGER, MPI_MAX, NEKO_COMM, ierr)
    itmp = maxval(ref_mark)
    call MPI_Allreduce(itmp, jl, 1, MPI_INTEGER, MPI_MIN, NEKO_COMM, ierr)
    if (il == 0 .and. jl == 0) then
       call neko_log%message('No refinemnt mark set')
       call neko_log%end_section()
       return
    end if

    ! save old number of local elements
    !amr_rcn%msh_rcn%nelvo = msh%nelv ! this should be done in other place
    ! perform refinement/coarsening on mesh manager side
    level_max = param%amrlmax
    ! global elements numbers; THIS IS JUST A TEMPORARY HACK
    allocate(el_gidx(msh%nelv))
    do il = 1, msh%nelv
       el_gidx(il) = msh%offset_el + il
    end do
    select type(msh_imp => this%msh_mngr%msh_imp)
    class is (mesh_cnstr_amr_t)
       call msh_imp%refine(ref_mark, el_gidx, level_max, ifmod)
    class default
       call neko_error('Invalid refinement class type')
    end select
    deallocate(el_gidx)

    if (ifmod)  then
       ! Import mesh from mesh manager to neko
       call this%msh_mngr%msh_imp%mesh_get(msh)

!!$       ! PLACE FOR NEW MESH PARTITIONING
!!$
!!$       ! PLACE FOR TRANSFER/SORTING OF MAPPING DATA DATA
!!$       ! THIS SHOULD BE CHANGED AFTER ADDING COMMUNICATION AS IT WOULD INCLUDE COPY
!!$       ! collect information: old gidx, old nid, new gidx, new nid
!!$       ! the max size of tmp arrays could be max(nelv_new,nelv_old)*children_number
!!$       
!!$
!!$       ! Reshuffle data to get propper mapping on the destination rank
!!$       ! THIS SHOULD BE CHANGED AFTER ADDING COMMUNICATION
!!$       allocate(lel_map(3,max(msh%nelv,amr_rcn%msh_rcn%nelvo)*msh%npts)) ! at this point I do not know the proper size, so take a safe value
!!$       lel_map(:,:) = 0
!!$       lel_map(3,:) = -1
!!$       do il = 1, msh%nelv ! loop witn neko distribution
!!$          if (amr_rcn%msh_rcn%elgl_map(1,il) /= 0) then
!!$             lel_map(1, amr_rcn%msh_rcn%elgl_map(2,il)) = msh%offset_el + il ! this is just speciffic to a test case; in general wrong
!!$             lel_map(3, amr_rcn%msh_rcn%elgl_map(2,il)) = pe_rank ! this is just speciffic to a test case; in general wrong
!!$          end if
!!$       end do
!!$       call MOVE_ALLOC(lel_map,amr_rcn%msh_rcn%elgl_map)
!!$       ! PLACE TO UPDATE amr_rcn%msh_rcn%map_nr TO NEKO LOCAL VALUE
!!$
!!$       ! PLACE FOR TRANSFER/SORTING OF REFINEMENT DATA
!!$
!!$       ! Recalculate local position of refined element on the destination rank
!!$       ! THIS SHOULD BE CHANGED AFTER ADDING COMMUNICATION AS IT WOULD INCLUDE COPY
!!$       ! RIGHT NOW IT IS JUST LOCAL POSITION UPDATE
!!$       ! SORTING IS IMPORTANT
!!$       ! loop over refined elements
!!$       itmp = amr_rcn%msh_rcn%nelvo
!!$       do il = 1, amr_rcn%msh_rcn%rfn_nr, amr_rcn%msh%npts ! loop over parent owner
!!$          ips = itmp
!!$          ielo = amr_rcn%msh_rcn%elgl_rfn(3, il)
!!$          ! THIS COPY IS HERE AS THERE IS NO COMMUNICATION
!!$          itmpv(:, 1:amr_rcn%msh%npts) = &
!!$               & amr_rcn%msh_rcn%elgl_rfn(:, il:il + amr_rcn%msh%npts - 1)
!!$          ! loop over all the children
!!$          do jl = 1, amr_rcn%msh%npts
!!$             ! which child
!!$             ich = itmpv(5, jl)
!!$             ! local parent position sanity check
!!$             if (itmpv(3, jl) /= ielo) &
!!$                  & call neko_error('Children do not share parent')
!!$             ! new local position at the end of the array
!!$             if (ich == 0) then
!!$                amr_rcn%msh_rcn%elgl_rfn(3, il + ich) = ielo
!!$             else
!!$                itmp = itmp + 1
!!$                amr_rcn%msh_rcn%elgl_rfn(3, il + ich) = ips + ich
!!$             end if
!!$             ! new global element number
!!$             amr_rcn%msh_rcn%elgl_rfn(1, il + ich) = itmpv(1, jl)
!!$             ! parent local position in the array
!!$             amr_rcn%msh_rcn%elgl_rfn(2, il + ich) = ielo
!!$          end do
!!$       end do
!!$       ! save local number of elements after refinement
!!$       amr_rcn%msh_rcn%rfn_nr_a = itmp
!!$       ! PLACE TO UPDATE amr_rcn%msh_rcn%rfn_nr TO NEKO LOCAL VALUE
!!$
!!$       ! update global rank mapping
!!$       do il = 1, amr_rcn%msh_rcn%rfn_nr ! loop over parent owner
!!$          jl = amr_rcn%msh_rcn%elgl_rfn(3, il)
!!$          if (amr_rcn%msh_rcn%elgl_map(1, jl) == 0) then
!!$             amr_rcn%msh_rcn%elgl_map(1, jl) = amr_rcn%msh_rcn%elgl_rfn(1, il)
!!$             amr_rcn%msh_rcn%elgl_map(3, jl) = pe_rank ! this is just speciffic to a test case; in general wrong
!!$          else
!!$             call neko_error('Refinement transfer index already used; refinement.')
!!$          end if
!!$       end do
!!$
!!$       ! Provide global numbering of elements used for coarsening
!!$       call MPI_Allreduce(amr_rcn%msh_rcn%crs_nr, gnum, 1, MPI_INTEGER, MPI_SUM, & !????????
!!$            & NEKO_COMM, ierr)
!!$       ! get global offset
!!$       call MPI_Scan(amr_rcn%msh_rcn%crs_nr, goff, 1, MPI_INTEGER, MPI_SUM, &
!!$            & NEKO_COMM, ierr)
!!$       goff = goff - amr_rcn%msh_rcn%crs_nr
!!$       if (amr_rcn%msh_rcn%crs_nr > 0) then
!!$          itmp = msh%glb_nelv + (msh%npts - 1)*goff
!!$          do il = 1, amr_rcn%msh_rcn%crs_nr ! loop on mesh manager distribution
!!$             do jl = 2, msh%npts
!!$                itmp = itmp + 1
!!$                amr_rcn%msh_rcn%elgl_crs(1, jl, il) = itmp
!!$             end do
!!$          end do
!!$       end if
!!$
!!$       ! PLACE FOR COARSENING DATA TRANSFER AND SORTING; child owner
!!$
!!$       ! save local number of elements to be coarsened (child owner)
!!$       ! crs_nr_s = ???? ! THIS IS JUST TEST SPECIFFIC in general it should be result of communication
!!$       if (amr_rcn%msh_rcn%crs_nr == 0) then
!!$          amr_rcn%msh_rcn%crs_nr_s = 0
!!$       else
!!$          amr_rcn%msh_rcn%crs_nr_s = msh%nelv
!!$       end if
!!$       ! sanity check (sum of unchaged, refined parents and corsend children has to be equal nelvo)
!!$       if (amr_rcn%msh_rcn%nelvo /= amr_rcn%msh_rcn%crs_nr_s + &
!!$            & amr_rcn%msh_rcn%map_nr + &
!!$            & amr_rcn%msh_rcn%rfn_nr/msh%npts) call neko_error('Inconsiten')
!!$
!!$       ! update global rank mapping
!!$       ! THIS LOOP SHOULD LOOK DIFFERENT WITH COMMUNICATION
!!$       do il = 1, amr_rcn%msh_rcn%crs_nr ! loop over elements on child owner
!!$          do jl = 1, msh%npts
!!$             if (amr_rcn%msh_rcn%elgl_map(1, amr_rcn%msh_rcn%elgl_crs(3, jl, il)) == 0) then
!!$                amr_rcn%msh_rcn%elgl_map(1, jl) = amr_rcn%msh_rcn%elgl_crs(1, jl, il)
!!$                amr_rcn%msh_rcn%elgl_map(3, jl) = pe_rank ! this is just speciffic to a test case; in general wrong
!!$             else
!!$                call neko_error('Refinement transfer index already used; coarsening.')
!!$             end if
!!$          end do
!!$       end do
!!$
!!$       ! PLACE FOR COARSENING DATA TRANSFER AND SORTING; coarsened element owner
!!$
!!$       ! recalculate new childen position
!!$       ! THIS LOOP HAS TO BE CHANGED FOR COMMUNICATION
!!$       itmp = msh%nelv
!!$       do il = 1, amr_rcn%msh_rcn%crs_nr ! loop over children on parent owner
!!$          ips = itmp
!!$          !ielo = PARENT GLOBAL NUMBER
!!$          ! loop over all the children
!!$          do jl = 1, amr_rcn%msh%npts
!!$             ! sanity check; global parent position
!!$             ! NOT NEEDED FOR THIS EXAMPLE
!!$             ! which child
!!$             ! EQUAL TO JL IN THIS EXAMPLE
!!$             ich = jl
!!$             ! new global element number
!!$             ! DOESN'T CHANGE IN THIS EXAMPLE
!!$             ! local child position
!!$             if (ich == 1) then
!!$                amr_rcn%msh_rcn%elgl_crs(1, jl, il) = il ! this is not correct in general
!!$             else
!!$                itmp = itmp + 1
!!$                amr_rcn%msh_rcn%elgl_crs(1, jl, il) = ips + ich - 1
!!$             end if
!!$          end do
!!$       end do
!!$       ! save local number of elements before coarsening
!!$       amr_rcn%msh_rcn%crs_nr_b = itmp
!!$       ! PLACE TO UPDATE amr_rcn%msh_rcn%crs_nr TO NEKO LOCAL VALUE (parent owner)
          
    end if ! ifmod

!!$    testing : block
!!$      integer :: ierr
!!$      write(*,*) 'TESTING0', pe_rank, amr_rcn%msh_rcn%nelvo, amr_rcn%msh%nelv, &
!!$           & amr_rcn%msh_rcn%map_nr, amr_rcn%msh_rcn%rfn_nr, amr_rcn%msh_rcn%crs_nr,&
!!$           & amr_rcn%msh_rcn%rfn_nr_a, amr_rcn%msh_rcn%crs_nr_s, amr_rcn%msh_rcn%crs_nr_b, goff, gnum
!!$      call MPI_Barrier(NEKO_COMM, ierr)
!!$      !call neko_log%end_section()
!!$      return
!!$      !call neko_error('This is not an error.')
!!$    end block testing

    ! free memory
    deallocate(ref_mark)

    call neko_log%end_section()

    return
  end subroutine amr_msh_refine

  !> Perform field refinement
  !! @param[inout]  msh      mesh information
  !! @param[inout]  fld      field
  !! @param[in]     param    run-time parameters
  !! @param[in]     usr      user interface
  !! @param[in]     time     simulation time
  !! @param[in]     tstep    time step
  subroutine amr_fld_refine(this, msh, fld, param, usr, time, tstep)
    ! argument list
    class(amr_t), intent(inout) :: this
    type(mesh_t), intent(inout) :: msh
    class(fluid_scheme_t), intent(inout) :: fld
    type(param_t), intent(in) :: param
    type(user_t), intent(in) :: usr
    real(kind=rp), intent(in) :: time
    integer, intent(in) :: tstep
    ! local variables
    integer(i4) :: itmp
    real(dp), allocatable, dimension(:,:,:,:) :: tmpv

    ! ALL THIS SHOULD BE PART OF TYPE EXTENSION, BUT FOR TESTING IT SHOULD BE FINE
    select type(fld)
    type is(fluid_plan1_t)
       call neko_error('Nothing done for plan1')
    type is(fluid_plan4_t)
       call neko_error('Nothing done for plan4')
    type is(fluid_pnpn_t)
       ! perform mesh refinement
       call this%msh_refine(msh, param, usr, time, tstep)

!       !test refinement on a single field
!       allocate(tmpv(fld%Xh%lx, fld%Xh%ly, fld%Xh%lz,msh%nelv))
!       tmpv(:,:,:,1:amr_rcn%msh_rcn%nelvo) = fld%dm_Xh%x(:,:,:,:)
!       call amr_rcn_refine_coarsen_single(amr_rcn, tmpv)
!       ! reinitialise variables
!       ! degrees of freedom; THIS SHOULD BE POSSIBLY CHANGED TO SUBROUTINE TAKING VARIABLE
!       !call fld%dm_Xh%resize()

  

!       deallocate(tmpv)

       testing : block
         integer :: ierr
         write(*,*) 'TEST size', pe_rank,  msh%nelv, msh%glb_nelv
         call MPI_Barrier(NEKO_COMM, ierr)
         call neko_error('This is not an error.')
       end block testing
    end select

    return
  end subroutine amr_fld_refine

end module amr
