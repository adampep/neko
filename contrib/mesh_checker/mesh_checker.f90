! Copyright (c) 2024-2025, The Neko Authors
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
! Mesh diagnostics tool
program mesh_checker
  use neko
  implicit none

  character(len=NEKO_FNAME_LEN) :: inputchar, mesh_fname
  type(file_t) :: mesh_file
  type(mesh_t) :: msh
  integer :: argc, i, n_labeled, ierr
  character(len=LOG_SIZE) :: log_buf
  integer :: total_size, inlet_size, wall_size, periodic_size
  integer :: outlet_size, symmetry_size, outlet_normal_size
  type(dirichlet_t) :: bdry_mask
  type(field_t) :: bdry_field
  type(file_t) :: bdry_file
  type(space_t) :: Xh
  type(gs_t) :: gs
  type(coef_t) :: coef
  type(dofmap_t) :: dofmap
  logical :: write_zone_ids

  argc = command_argument_count()

  if (argc .lt. 1) then
     if (pe_rank .eq. 0) then
        write(*,*) 'Usage: ./mesh_checker mesh.nmsh [--write_zone_indices]'
        write(*,*) '--write_zone_indices : write a field file with boundaries &
             & marked by zone index.'
     end if
     stop
  end if

  call neko_init

  call get_command_argument(1, inputchar)
  read(inputchar, *) mesh_fname

  write_zone_ids = .false.

  do i = 2, argc
     call get_command_argument(i, inputchar)
     if (trim(inputchar) == "--write_zone_indices") then
        write_zone_ids = .true.
     end if
  end do

  mesh_file = file_t(trim(mesh_fname))

  call mesh_file%read(msh)

  call MPI_Allreduce(msh%inlet%size, inlet_size, 1, &
       MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
  call MPI_Allreduce(msh%wall%size, wall_size, 1, &
       MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
  call MPI_Allreduce(msh%outlet%size, outlet_size, 1, &
       MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
  call MPI_Allreduce(msh%outlet_normal%size, outlet_normal_size, 1, &
       MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
  call MPI_Allreduce(msh%sympln%size, symmetry_size, 1, &
       MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
  call MPI_Allreduce(msh%periodic%size, periodic_size, 1, &
       MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)

  total_size = inlet_size + wall_size + outlet_size + outlet_normal_size &
       + symmetry_size

  if (pe_rank .eq. 0) then
     write(*,*) ''
     write(*,*) '--------------Size-------------'
     write(*,*) 'Number of elements: ', msh%glb_nelv
     write(*,*) 'Number of points:   ', msh%glb_mpts
     write(*,*) 'Number of faces:    ', msh%glb_mfcs
     write(*,*) 'Number of edges:    ', msh%glb_meds
     write(*,*) ''
     write(*,*) '--------------Zones------------'
     write(*,*) 'Number of built-in inlet faces:         ', inlet_size
     write(*,*) 'Number of built-in wall faces:          ', wall_size
     write(*,*) 'Number of built-in outlet faces:        ', outlet_size
     write(*,*) 'Number of built-in outlet-normal faces: ', outlet_normal_size
     write(*,*) 'Number of built-in symmetry faces:      ', symmetry_size
     write(*,*) 'Number of periodic faces:               ', periodic_size
     write(*,*) ''
     if (total_size .gt. 0) then
        write(*,*) 'WARNING: Your mesh contains non-periodic "built-in" zones,&
             & which are deprecated. Your mesh must have been created natively&
             & by Nek5000, e.g. with genbox.'
     end if
     write(*,*) 'Labeled zones: '
  end if

  do i = 1, size(msh%labeled_zones)

     call MPI_Allreduce(msh%labeled_zones(i)%size, total_size, 1, &
          MPI_INTEGER, MPI_SUM, NEKO_COMM, ierr)
     if (total_size .gt. 0 .and. pe_rank .eq. 0) then
        write(*,'(A,I2,A,I0,A)') '    Zone ', i, ': ', total_size, &
             ' faces'

     end if
  end do

  if (write_zone_ids) then
     if (pe_rank .eq. 0) write(*,*) 'Writing zone ids to zone_indices0.f00000'
     call Xh%init(1, 2, 2, 2)
     call dofmap%init(msh, Xh)
     call gs%init(dofmap)
     call coef%init(gs)

     call bdry_field%init(dofmap)

     do i = 1, size(msh%labeled_zones)
        call bdry_mask%init_from_components(coef, real(i, kind=rp))
        call bdry_mask%mark_zone(msh%labeled_zones(i))
        call bdry_mask%finalize()
        call bdry_mask%apply_scalar(bdry_field%x, dofmap%size())
        call bdry_mask%free()
     end do

     bdry_file = file_t('zone_indices.fld')
     call bdry_file%write(bdry_field)

     call Xh%free()
     call gs%free()
     call dofmap%free()
     call bdry_field%free()
     call msh%free()
     call bdry_mask%free()
  end if

  if (pe_rank .eq. 0) write(*,*) 'Done'
  call neko_finalize

end program mesh_checker

