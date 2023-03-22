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
!> Main interface for data exchange and manipulation between p4est and neko
module p4est
  use mpi_f08
  use num_types
  use comm
  use logger
  use utils
  use math
  use tuple
  use htable
  use mesh_cnstr_amr
  use mesh_redistribute
  use distdata
  use point
  use quad
  use hex
  use mesh_conn
  use mesh
  implicit none

  private
  public :: p4_mesh_cnstr_t, p4_manager_init, p4_manager_free

  ! Following node types contain geometrical and connectivity information
  ! for the mesh.
  ! Type for independent nodes.
  ! These are globally unique nodes not sharing physical coordinates and global
  ! id in the communicator. It is a minimal set of points needed to build
  ! a consistent element connectivity map. This set of nodes does not include
  ! e.g. repeated many times nodes on periodic faces sharing the same global
  ! id in the communicator.
  type p4_node_ind_t
     integer(i4) :: lown ! number of owned independent nodes
     integer(i4) :: lshr ! number of independent owned shared nodes
     integer(i4) :: loff ! local offset of owned independent nodes
     integer(i4) :: lnum ! local number of independent nodes
     integer(i8), allocatable, dimension(:) :: gidx ! global indexing of
                                                    ! unique nodes
     integer(i4), allocatable, dimension(:) :: ndown ! node owner (mpi rank)
     ! physical node coordinates
     real(kind=dp), allocatable, dimension(:,:) :: coord
     ! NOTICE! p4est provides only approximate coordinates that are exact
     ! for max refinement level = 0, or for linear elements (for any ref. lev.).
   contains
     procedure, public, pass(this) :: free => p4_node_ind_free
  end type p4_node_ind_t

  ! Type for periodic nodes
  ! These are element vertices that are not included in independent node list
  ! as they are located on a periodic boundary (have unique physical
  ! coordinates) and share the same global id in the communicator. These
  ! nodes do not include hanging nodes. The hanging nodes located on the
  ! periodic bc are marked as hanging nodes only.
  ! Notice! The element vertices from the set of vertices on periodic boundaries
  ! may share the same global id in the communicator, so we can map periodic
  ! nodes to the independent ones.
  type p4_node_per_t
     integer(i4) :: lnum ! local number of periodic nodes
     integer(i4), allocatable, dimension(:) :: lmap ! local periodic to
                                                    ! independent node mapping
     real(kind=dp), allocatable, dimension(:,:) :: coord ! physical node
                                                         ! coordinates
     ! NOTICE! p4est provides only approximate coordinates that are exact
     ! for max refinement level = 0, or for linear elements (for any ref. lev.).
   contains
     procedure, public, pass(this) :: free => p4_node_per_free
  end type p4_node_per_t

  ! Type for hanging nodes
  ! Hanging nodes are neither independent nor periodic. They are located in
  ! the centre of the nonconforming face or edge (have unique physical
  ! coordinates) and can be maped to the independent face/edge vertices. This
  ! mapping is not unique and will depend on the element position in a tree.
  ! There are two types of hanging nodes h2 (dependent on 2 independent nodes;
  ! 2D face and 3D edge) and h4 (dependent on 4 independent nodes; 3D face).
  type p4_node_hng_t
     integer(i4) :: lnum ! local number of given type hanging nodes
     integer(i4), allocatable, dimension(:,:) :: lmap ! local hanging to
                                                      ! independent node mapping
     real(kind=dp), allocatable, dimension(:,:) :: coord ! physical node
                                                         ! coordinates
     ! NOTICE! p4est provides only approximate coordinates that are exact
     ! for max refinement level = 0, or for linear elements (for any ref. lev.).
   contains
     procedure, public, pass(this) :: free => p4_node_hng_free
  end type p4_node_hng_t

  ! Type for connectivity information regarding vertices, faces and edges.
  ! It contains both object global numbering and communication information
  type p4_obj_conn_t
     integer(i4) :: lnum, lown ! number of local and owned objects
     integer(i8) :: goff ! global object offset
     integer(i8) :: gnum ! global number of objects
     integer(i4), allocatable, dimension(:,:) :: lmap ! element
                                                      ! vertices/faces/edges
                                                      ! to object mapping
     integer(i4) :: nrank, nshare ! number of MPI ranks sharing objects and
                                  ! number of shared objects
     integer(i8), allocatable, dimension(:) :: lgidx ! global indexing of unique
                                                     ! objects of given type
     integer(i4), allocatable, dimension(:) :: lrank ! list of ranks sharing
                                                     ! objects
     integer(i4), allocatable, dimension(:) :: lshare ! list of shared objects
     integer(i4), allocatable, dimension(:) :: loff ! offset in the lshare list
   contains
     procedure, public, pass(this) :: free => p4_obj_conn_free
  end type p4_obj_conn_t

  ! This type contains main element information
  type p4_element_t
     integer(i4) :: nelt, nelv ! local number of temperature and velocity
                               ! elements (Nek5000 specific; CHANGE)
     integer(i8) :: nelgt, nelgto ! global element number, global element offset
     integer(i8), allocatable, dimension(:) :: gidx ! global element index
     integer(i4), allocatable, dimension(:) :: level ! element refinement level
     integer(i4), allocatable, dimension(:) :: igrp ! element group
                                                    ! (not used right now)
     integer(i4), allocatable, dimension(:,:) :: crv ! face curvature flag
                                                     ! (not used right now)
     integer(i4), allocatable, dimension(:,:) :: bc ! face boundary condition;
                                                    ! -1- periodic,
                                                    ! 0-internal,
                                                    ! 0< user specified

     ! Local mapping of element vertices to the nodes
     ! for:
     ! 1<= vnmap(iv,iel) <= nin - independent node
     ! nin < vnmap(iv,iel) <= nin + npe - periodic node
     ! nin + npe < vnmap(iv,iel) <= nin + npe + nhf - face hanging node
     ! nin + +npe + nhf < vnmap(iv,iel) <= nin + npe + nhf + nhe - edge hanging node
     ! where
     ! nin = number of local independent nodes
     ! npe = number of local periodic nodes
     ! nhf = number of local face hanging nodes
     ! nhe = number of local edge hanging nodes
     ! vnmap uses symmetric vertex notation with (r,s,t) being a local
     ! counterpart of (x,y,z):
     !             3+--------+4    ^ s
     !             /        /|     |
     !            /        / |     |
     !           /        /  |     |
     !         7+--------+8  +2    +----> r
     !          |        |  /     /
     !          |        | /     /
     !          |        |/     /
     !         5+--------+6    t
     !
     integer(i4), allocatable, dimension(:,:) :: vnmap

     ! Face orientation.
     ! In 2D case permutation array takes 2 values
     !          0 => no permutations
     !          1 => row permutations
     ! In 3D case there are 8 possible transformations for faces
     ! Vertex numbering  (0,1,2,3)
     !                       2--3
     !                       |  |
     !                       0--1
     ! They are numbered according to p4est by 0..7 and stored in falg
     ! { 0, 1, 2, 3 }      0 => identity
     ! { 0, 2, 1, 3 }      1 => T
     ! { 1, 0, 3, 2 }      2 => P_x
     ! { 1, 3, 0, 2 }      3 => P_x T
     ! { 2, 0, 3, 1 }      4 => P_y T
     ! { 2, 3, 0, 1 }      5 => P_y
     ! { 3, 1, 2, 0 }      6 => P_y P_x T
     ! { 3, 2, 1, 0 }      7 => P_y P_x
     ! where   T - transpose
     !               P_x - permutation in x (rows)
     !               P_y - permutation in y (columns)
     integer(i4), allocatable, dimension(:,:) :: falg ! face alignment
     ! Edge orientation is similar to 2D face, however p4est does not build
     ! a global edge orientation and keeps just relative orientation of two
     ! edges. That is why I redo it introducing global orientation based on
     ! vertex global id. THIS HES TO BE CHANGED
     integer(i4), allocatable, dimension(:,:) :: ealg ! edge alignment

     ! Topology information: hanging elements, faces, edges and vertices
     ! The element is marked as hanging if it contains at least one hanging
     ! object.
     ! Hanging face
     ! 2D
     !   = -1 if the face is not hanging,
     !   = 0 if the face is the first half,
     !   = 1 if the face is the second half.
     ! 3D
     !   = -1 if the face is not hanging,
     !   = the corner of the full face that it touches:
     ! Hanging edge
     !   = -1 if the edge is not hanging,
     !   =  0 if the edge is the first half of a full edge,
     !        but neither of the two faces touching the
     !        edge is hanging,
     !   =  1 if the edge is the second half of a full edge,
     !        but neither of the two faces touching the
     !        edge is hanging,
     !   =  2 if the edge is the first half of a full edge
     !        and is on the boundary of a full face,
     !   =  3 if the edge is the second half of a full edge
     !        and is on the boundary of a full face,
     !   =  4 if the edge is in the middle of a full face.
     !        See the diagram below for clarification.
     !  o...............o o...............o +---2---+.......o o.......+---3---+
     !  :               : :               : |       |       : :       |       |
     !  :               : :               : 3   2   4       : :       4   3   3
     !  :               : :               : |       |       : :       |       |
     !  +---4---+       : :       +---4---+ +---4---+       : :       +---4---+
     !  |       |       : :       |       | :               : :               :
     !  2   0   4       : :       4   1   2 :               : :               :
     !  |       |       : :       |       | :               : :               :
     !  +---2---+.......o o.......+---3---+ o...............o o...............o
     !
     !                     o                  +-------+
     !                     :                  |\       \
     !                     :                  1 \       \
     !                     :                  |  +-------+
     !                     +-------+          +  |       |
     !                     |\       \         :\ |       |
     !                     0 \       \        : \|       |
     !                     |  +-------+       :  +-------+
     !                     +  |       |       o
     !                      \ |       |
     !                       \|       |
     !                        +-------+
     integer(i4), allocatable, dimension(:) :: hngel ! hanging element list
                      ! (1 if at least one hanging face or edge otherwise 0)
     integer(i4), allocatable, dimension(:,:) :: hngfc ! hanging face list;
                                    ! position of hanging face; otherwise -1
     integer(i4), allocatable, dimension(:,:) :: hnged ! hanging edge list;
                          ! position of hanging edge; otherwise -1 (3D only)
!!!!!!!!!!!!! vertex hanging information can be directly extracted from vnmap;
!!!!!!!!!!!!! It can bee I will have to add an array

     ! Flag elements that can be coarsened and share the same parent
     ! fmlm(1,LELT) - mark of a parent (not a real element number as parents
     ! do not exist on neko side)
     !                         0 - element that cannot be coarsened
     !                        >0 - family mark
     ! fmlm(2,LELT) - vertex number shared by all the family members
     integer(i8), allocatable, dimension(:,:) :: fmlm

     ! Following data could be recreated from nodes, but for now I stick to
     ! what p4est can provide.
     ! Node info for building communicators; includes global indexing of various
     ! element objects
     ! vertex%lmap uses symmetric vertex notation with (r,s,t) being a local
     ! counterpart of (x,y,z):
     !             3+--------+4    ^ s
     !             /        /|     |
     !            /        / |     |
     !           /        /  |     |
     !         7+--------+8  +2    +----> r
     !          |        |  /     /
     !          |        | /     /
     !          |        |/     /
     !         5+--------+6    t
     !
     type(p4_obj_conn_t) :: vert ! vertex info
     ! face%lmap uses symmetric face notation with (r,s,t) being a local
     ! counterpart of (x,y,z):
     !              +--------+     ^ s
     !             /        /|     |
     !            /    4   / |     |
     !      1--> /        /  |     |
     !          +--------+ 2 +     +----> r
     !          |        |  /     /
     !          |    6   | /     /
     !          |        |/     /
     !          +--------+     t
     !               3
     type(p4_obj_conn_t) :: face ! face info
     ! edge%lmap uses symmetric edge notation with (r,s,t) being a local
     ! counterpart of (x,y,z):
     !              +---2----+     ^ s
     !             /        /|     |
     !           11       12 6     |
     !           /        /  |     |
     !          +---4----+   +     +----> r
     !          |        |  /     /
     !          7        8 10    /
     !          |        |/     /
     !          +---3----+     t
     !
     type(p4_obj_conn_t) :: edge ! edges info; 3D mesh only
   contains
     procedure, public, pass(this) :: free => p4_element_free
  end type p4_element_t

  ! Type combining all the mesh information that can be extracted from p4est
  type, extends(mesh_cnstr_amr_t) :: p4_mesh_cnstr_t
     ! general information
     integer :: dim, maxl, maxg ! dimension, local and global max refinement
                                ! level

     ! unique node information
     type(p4_node_ind_t) :: indn ! independent nodes
     type(p4_node_per_t) :: pern ! periodic nodes
     type(p4_node_hng_t) :: fchn ! face hanging nodes
     type(p4_node_hng_t) :: edhn ! edge hanging nodes

     ! element information
     type(p4_element_t) :: elem
   contains
     procedure, public, pass(this) :: msh_get => p4_msh_get
     procedure, public, pass(this) :: refine => p4_refine
     procedure, public, pass(this) :: free => p4_mesh_cnstr_free
  end type p4_mesh_cnstr_t

  ! connectivity parameter arrays
  ! face vertices
  integer, parameter, dimension(4,6) :: p4_vface = reshape(&
       &(/ 1,3,5,7 , 2,4,6,8 , 1,2,5,6 , 3,4,7,8 , 1,2,3,4 , 5,6,7,8 /),&
       &shape(p4_vface))

  ! edge vertices
  integer, parameter, dimension(2,12) :: p4_vedge  = reshape(&
       &(/ 1,2 , 3,4 , 5,6 , 7,8 , 1,3 , 2,4 , 5,7 , 6,8 , 1,5 , 2,6 , 3,7&
       & , 4,8 /),shape(p4_vedge))

  ! edge related faces
  integer, parameter, dimension(2,12) :: p4_eface  = reshape(&
       &(/ 3,5 , 4,5 , 3,6 , 4,6 , 1,5 , 2,5 , 1,6 , 2,6 , 1,3 , 2,3 , 1,4&
       & , 2,4 /),shape(p4_eface))

  ! corner related faces
  integer, parameter, dimension(3,8) :: p4_cface = reshape(&
       &(/ 1,3,5 , 2,3,5 , 1,4,5 , 2,4,5 , 1,3,6 , 2,3,6 , 1,4,6 , 2,4,6 /),&
       &shape(p4_cface))

  ! corner related edges
  integer, parameter, dimension(3,8) :: p4_cedge = reshape(&
       &(/ 1,5,9 , 1,6,10 , 2,5,11 , 2,6,12 , 3,7,9 , 3,8,10 , 4,7,11 ,&
       & 4,8,12 /),shape(p4_cedge))

  ! corner to face corner
  integer, parameter, dimension(6,8) :: p4_cfcrn = reshape(&
       &(/ 1,-1, 1,-1, 1,-1 , -1, 1, 2,-1, 2,-1 ,  2,-1,-1, 1, 3,-1 &
       &, -1, 2,-1, 2, 4,-1 ,  3,-1, 3,-1,-1, 1 , -1, 3, 4,-1,-1, 2 &
       &,  4,-1,-1, 3,-1, 3 , -1, 4,-1, 4,-1, 4 /),shape(p4_cfcrn))

  ! to calculate neighbour face corner
  integer, parameter, dimension(6,6) :: p4_rt =reshape( &
       &(/ 1,2,2,1,1,2 , 3,1,1,2,2,1 , 3,1,1,2,2,1 , 1,3,3,1,1,2 &
       &, 1,3,3,1,1,2 , 3,1,1,3,3,1 /),shape(p4_rt))
  integer, parameter, dimension(4,3) :: p4_qt = reshape(&
       &(/ 2,3,6,7 , 1,4,5,8 , 1,5,4,8 /),shape(p4_qt))
  integer, parameter, dimension(4,8) :: p4_pt = reshape(&
       &(/ 1,2,3,4 , 1,3,2,4 , 2,1,4,3 , 2,4,1,3 , 3,1,4,2 , 3,4,1,2 &
       &, 4,2,3,1 , 4,3,2,1 /),shape(p4_pt))

  ! default log threshold - production; for more info see sc.h
  integer, parameter :: p4_lp_production = 6

  ! set of interfaces to call p4est functions
  interface
     subroutine wp4est_init(fmpicomm, catch_signals, print_backtrace, &
          & log_threshold) bind(c, name = 'wp4est_init')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int), value :: fmpicomm, catch_signals, &
            & print_backtrace, log_threshold
     end subroutine wp4est_init

     subroutine wp4est_finalize(log_priority) bind(c, name = 'wp4est_finalize')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int), value :: log_priority
     end subroutine wp4est_finalize

     subroutine wp4est_cnn_del() bind(c, name = 'wp4est_cnn_del')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_cnn_del

     subroutine wp4est_cnn_valid(is_valid) bind(c, name = 'wp4est_cnn_valid')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int) :: is_valid
     end subroutine wp4est_cnn_valid

     subroutine wp4est_tree_del() bind(c, name = 'wp4est_tree_del')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_tree_del

     subroutine wp4est_tree_valid(is_valid) bind(c, name = 'wp4est_tree_valid')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int) :: is_valid
     end subroutine wp4est_tree_valid

     subroutine wp4est_tree_save(save_data,filename) &
          & bind(c, name = 'wp4est_tree_save')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int), value :: save_data
       character(kind=c_char), dimension(*) :: filename
     end subroutine wp4est_tree_save

     subroutine wp4est_tree_load(fmpicomm, load_data, filename) &
          & bind(c, name = 'wp4est_tree_load')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int), value :: fmpicomm, load_data
       character(kind=c_char), dimension(*) :: filename
     end subroutine wp4est_tree_load

     subroutine wp4est_ghost_new() bind(c, name = 'wp4est_ghost_new')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_ghost_new

     subroutine wp4est_ghost_del() bind(c, name = 'wp4est_ghost_del')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_ghost_del

     subroutine wp4est_mesh_new() bind(c, name = 'wp4est_mesh_new')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_mesh_new

     subroutine wp4est_mesh_del() bind(c, name = 'wp4est_mesh_del')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_mesh_del

     subroutine wp4est_nodes_new() bind(c, name = 'wp4est_nodes_new')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_nodes_new

     subroutine wp4est_nodes_del() bind(c, name = 'wp4est_nodes_del')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_nodes_del

     subroutine wp4est_lnodes_new(degree) bind(c, name = 'wp4est_lnodes_new')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int), value :: degree
     end subroutine wp4est_lnodes_new

     subroutine wp4est_lnodes_del() bind(c, name = 'wp4est_lnodes_del')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_lnodes_del

     subroutine wp4est_part(partforcoarsen) bind(c, name = 'wp4est_part')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int), value :: partforcoarsen
     end subroutine wp4est_part

     subroutine wp4est_bc_check() bind(c, name = 'wp4est_bc_check')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_bc_check

     subroutine wp4est_msh_get_size(mdim, nelgt, nelgto, nelt, nelv, maxl) &
          & bind(c, name = 'wp4est_msh_get_size')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int) :: mdim, nelv, maxl
       integer(c_int32_t) :: nelt
       integer(c_int64_t) :: nelgt, nelgto
     end subroutine wp4est_msh_get_size

     subroutine wp4est_nds_get_size(nowin, nowsh, oowin, nin, nhf, nhe) &
          & bind(c, name = 'wp4est_nds_get_size')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int) :: nowin, nowsh, oowin, nin, nhf, nhe
     end subroutine wp4est_nds_get_size

     subroutine wp4est_nds_get_ind(nglid, nown, ncoord)&
          & bind(c, name = 'wp4est_nds_get_ind')
       USE, INTRINSIC :: ISO_C_BINDING
       type(c_ptr), value :: nglid, nown, ncoord
     end subroutine wp4est_nds_get_ind

     subroutine wp4est_nds_get_hfc(depend, ncoord)&
          & bind(c, name = 'wp4est_nds_get_hfc')
       USE, INTRINSIC :: ISO_C_BINDING
       type(c_ptr), value :: depend, ncoord
     end subroutine wp4est_nds_get_hfc

     subroutine wp4est_nds_get_hed(depend, ncoord)&
          & bind(c, name = 'wp4est_nds_get_hed')
       USE, INTRINSIC :: ISO_C_BINDING
       type(c_ptr), value :: depend, ncoord
     end subroutine wp4est_nds_get_hed

     subroutine wp4est_nds_get_vmap(vmap)&
          & bind(c, name = 'wp4est_nds_get_vmap')
       USE, INTRINSIC :: ISO_C_BINDING
       type(c_ptr), value :: vmap
     end subroutine wp4est_nds_get_vmap

     subroutine wp4est_elm_get_dat(gidx, level, igrp, crv, bc, coord, falg) &
          & bind(c, name = 'wp4est_elm_get_dat')
       USE, INTRINSIC :: ISO_C_BINDING
       type(c_ptr), value :: gidx, level, igrp, crv, bc, coord, falg
     end subroutine wp4est_elm_get_dat

     subroutine wp4est_elm_get_lnode(lnnum, lnown, lnoff, lnodes) &
          & bind(c, name = 'wp4est_elm_get_lnode')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int) :: lnnum, lnown
       integer(c_int64_t) :: lnoff
       type(c_ptr), value :: lnodes
     end subroutine wp4est_elm_get_lnode

     subroutine wp4est_sharers_get_size(nrank, nshare) &
          & bind(c, name = 'wp4est_sharers_get_size')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int) :: nrank, nshare
     end subroutine wp4est_sharers_get_size

     subroutine wp4est_sharers_get_ind(nglid, lrank, loff, lshare) &
          & bind(c, name = 'wp4est_sharers_get_ind')
       USE, INTRINSIC :: ISO_C_BINDING
       type(c_ptr), value :: nglid, lrank, loff, lshare
     end subroutine wp4est_sharers_get_ind

     subroutine wp4est_hang_get_info(hang_elm, hang_fsc, hang_edg) &
          & bind(c, name = 'wp4est_hang_get_info')
       USE, INTRINSIC :: ISO_C_BINDING
       type(c_ptr), value :: hang_elm, hang_fsc, hang_edg
     end subroutine wp4est_hang_get_info

     subroutine wp4est_fml_get_info(family, nelf) &
          & bind(c, name = 'wp4est_fml_get_info')
       USE, INTRINSIC :: ISO_C_BINDING
       type(c_ptr), value :: family
       integer(c_int) :: nelf
     end subroutine wp4est_fml_get_info

     subroutine wp4est_refine(max_level) &
          & bind(c, name = 'wp4est_refine')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int), value :: max_level
     end subroutine wp4est_refine

     subroutine wp4est_coarsen() &
          & bind(c, name = 'wp4est_coarsen')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_coarsen

     subroutine wp4est_balance() &
          & bind(c, name = 'wp4est_balance')
       USE, INTRINSIC :: ISO_C_BINDING
     end subroutine wp4est_balance

     subroutine wp4est_tree_copy(quad_data) &
          & bind(c, name = 'wp4est_tree_copy')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int), value :: quad_data
     end subroutine wp4est_tree_copy

     subroutine wp4est_tree_check(check, quad_data) &
          & bind(c, name = 'wp4est_tree_check')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int) :: check
       integer(c_int), value :: quad_data
     end subroutine wp4est_tree_check

     subroutine wp4est_refm_put(ref_mark) &
          & bind(c, name = 'wp4est_refm_put')
       USE, INTRINSIC :: ISO_C_BINDING
       type(c_ptr), value :: ref_mark
     end subroutine wp4est_refm_put

     subroutine wp4est_egmap_put(el_gnum, el_lnum, el_nid) &
          & bind(c, name = 'wp4est_egmap_put')
       USE, INTRINSIC :: ISO_C_BINDING
       type(c_ptr), value :: el_gnum, el_lnum, el_nid
     end subroutine wp4est_egmap_put

     subroutine wp4est_msh_get_hst(map_nr, rfn_nr, crs_nr, &
          & elgl_map, elgl_rfn, elgl_crs) &
          & bind(c, name = 'wp4est_msh_get_hst')
       USE, INTRINSIC :: ISO_C_BINDING
       integer(c_int) :: map_nr, rfn_nr, crs_nr
       type(c_ptr), value :: elgl_map, elgl_rfn, elgl_crs
     end subroutine wp4est_msh_get_hst

  end interface

contains

  ! Type bounded routines
  subroutine p4_node_ind_free(this)
    ! argument list
    class(p4_node_ind_t), intent(inout) :: this

    ! Deallocate arrays
    if (allocated(this%gidx)) deallocate(this%gidx)
    if (allocated(this%ndown)) deallocate(this%ndown)
    if (allocated(this%coord)) deallocate(this%coord)

    ! Reset registers
    this%lown = 0
    this%lshr = 0
    this%loff = 0
    this%lnum = 0

    return
  end subroutine p4_node_ind_free

  subroutine p4_node_per_free(this)
    ! argument list
    class(p4_node_per_t), intent(inout) :: this

    ! Deallocate arrays
    if (allocated(this%lmap)) deallocate(this%lmap)
    if (allocated(this%coord)) deallocate(this%coord)

    ! Reset registers
    this%lnum = 0

    return
  end subroutine p4_node_per_free

  subroutine p4_node_hng_free(this)
    ! argument list
    class(p4_node_hng_t), intent(inout) :: this

    ! Deallocate arrays
    if (allocated(this%lmap)) deallocate(this%lmap)
    if (allocated(this%coord)) deallocate(this%coord)

    ! Reset registers
    this%lnum = 0

    return
  end subroutine p4_node_hng_free

  subroutine p4_obj_conn_free(this)
    ! argument list
    class(p4_obj_conn_t), intent(inout) :: this

    ! Deallocate arrays
    if (allocated(this%lmap)) deallocate(this%lmap)
    if (allocated(this%lgidx)) deallocate(this%lgidx)
    if (allocated(this%lrank)) deallocate(this%lrank)
    if (allocated(this%lshare)) deallocate(this%lshare)
    if (allocated(this%loff)) deallocate(this%loff)

    ! Reset registers
    this%lnum = 0
    this%lown = 0
    this%goff = 0
    this%gnum = 0
    this%nrank = 0
    this%nshare = 0

    return
  end subroutine p4_obj_conn_free

  subroutine p4_element_free(this)
    ! argument list
    class(p4_element_t), intent(inout) :: this

    ! Deallocate arrays
    if (allocated(this%gidx)) deallocate(this%gidx)
    if (allocated(this%level)) deallocate(this%level)
    if (allocated(this%igrp)) deallocate(this%igrp)
    if (allocated(this%crv)) deallocate(this%crv)
    if (allocated(this%bc)) deallocate(this%bc)
    if (allocated(this%vnmap)) deallocate(this%vnmap)
    if (allocated(this%falg)) deallocate(this%falg)
    if (allocated(this%ealg)) deallocate(this%ealg)
    if (allocated(this%hngel)) deallocate(this%hngel)
    if (allocated(this%hngfc)) deallocate(this%hngfc)
    if (allocated(this%hnged)) deallocate(this%hnged)
    if (allocated(this%fmlm)) deallocate(this%fmlm)

    ! Free types
    call this%vert%free()
    call this%face%free()
    call this%edge%free()

    ! Reset registers
    this%nelt = 0
    this%nelv = 0
    this%nelgt = 0
    this%nelgto = 0

    return
  end subroutine p4_element_free

  subroutine p4_mesh_cnstr_free(this)
    ! argument list
    class(p4_mesh_cnstr_t), intent(inout) :: this

    ! Free types
    call this%indn%free()
    call this%pern%free()
    call this%fchn%free()
    call this%edhn%free()
    call this%elem%free()

    ! Reset registers
    this%dim = 0
    this%maxl = 0
    this%maxg = 0

    return
  end subroutine p4_mesh_cnstr_free

#ifdef HAVE_P4EST
  ! Initialise p4est as mesh manager
  subroutine p4_manager_init(mesh_file, log_threshold)
    character(len=*), intent(in) :: mesh_file
    integer, intent(in), optional :: log_threshold

    ! sc, p4est initialisation
    integer, parameter :: catch_signals=1, print_backtrace=1, load_data = 1
    integer :: log_thr, slen, ierr, is_valid
    character(18) :: log_cval
    ! p4est file reading
    character(len=80) :: suffix
    character(len=2) :: nlp

    ! get log threshold for p4est and sc libraries
    if (present(log_threshold)) then
       log_thr = log_threshold
    else
       if (pe_rank.eq.0) then
          call getenv('NEKOP4LOGL',log_cval)
          slen = len_trim(log_cval)
          if (slen.gt.0) then
             read(log_cval,'(I2)',iostat=ierr) log_thr
             if (ierr.gt.0) log_thr = p4_lp_production
          else
             log_thr = p4_lp_production ! default
          end if
       end if
    end if
    ! broadcast log threshold
    call MPI_Bcast(log_thr, 1, MPI_INTEGER, 0, NEKO_COMM, ierr)
    write(nlp,'(I2)') log_thr
    call neko_log%message('p4est will use log level '//nlp)
    call neko_log%message(&
         &'This log level can be changed by setting shell variable NEKOP4LOGL')
    call neko_log%message(&
         &"To suppress p4est logs 'export NEKOP4LOGL=9'")
    ! initialise sc and p4est
    call wp4est_init(NEKO_COMM%mpi_val, catch_signals, &
         & print_backtrace, log_thr)

    ! p4est file reading
    call filename_suffix(mesh_file, suffix)
    if (trim(suffix) /= 'p4est') call neko_error('Expected p4est mesh file')
    call neko_log%message('Reading p4est file ' // trim(mesh_file))
    call wp4est_tree_load(NEKO_COMM%mpi_val,load_data, &
         & trim(mesh_file)//c_null_char)
    ! check consistency of a forest and a connectivity
    call wp4est_tree_valid(is_valid)
    if (is_valid == 0) call neko_error('Invalid p4est tree')
    call wp4est_cnn_valid(is_valid)
    if (is_valid == 0) call neko_error('Invalid p4est connectivity')
    ! perform partitioning on p4est side
    call wp4est_part(1)

    return
  end subroutine p4_manager_init

  ! Finalise p4est as mesh manager
  subroutine p4_manager_free(log_priority)
    integer, intent(in), optional :: log_priority

    integer :: log_thr

    if (present(log_priority)) then
       log_thr = log_priority
    else
       log_thr = p4_lp_production ! default
    end if
    ! destroy p4est objects
    call wp4est_tree_del()
    call wp4est_cnn_del()
    ! finalise p4est and sc
    call wp4est_finalize(log_thr)

    return
  end subroutine p4_manager_free

  ! Construct mesh information using p4est imported data
  subroutine p4_msh_get(this, msh)
    ! argument list
    class(p4_mesh_cnstr_t), intent(inout) :: this
    type(mesh_t), intent(inout) :: msh
    ! local variables
    integer(i4) :: il

    call neko_log%section("Mesh")
    ! import data from p4est
    call p4_mesh_import_data(this)

    call neko_log%message('Build the mesh')

    ! What follows is just a hack, as I'm trying to minimise
    ! changes outside this file
    ! There are some differences in a concept, so not everything
    ! can be filled in properly

    ! Initialise mesh type
    call mesh_free(msh)
    ! general mesh info
    msh%gdim = this%dim ! grid dimension
    ! general element info
    msh%nelv = this%elem%nelv ! local number of elements
    msh%glb_nelv = this%elem%nelgt ! global number of elelments
    msh%offset_el = this%elem%nelgto ! global element offset;
    ! NOT NEEDED IN THIS VERSION
    msh%npts = 2**this%dim ! number of element vertices

    ! Point definition differs so I include here independent points only
    msh%mpts = this%elem%vert%lnum ! local number of unique vertices
    ! (independent nodes only)
    msh%mfcs = this%elem%face%lnum ! local number of unique faces
    msh%meds = this%elem%edge%lnum ! local number of unique edges

    msh%glb_mpts = this%elem%vert%gnum ! global number of unique vertices
    msh%glb_mfcs = this%elem%face%gnum ! global number of unique faces
    msh%glb_meds = this%elem%edge%gnum ! global number of unique edges

    ! this doesn't seem to be used in original version and
    ! I do not need it for p4est either
    !msh%max_pts_id =

    ! If there are any elements
    if (this%elem%nelv > 0) then

       ! Fill in point/node information
       call p4_nodes_fill(msh, this)

       ! Fill in element information
       call p4_element_fill(msh, this)

       ! msh%htp - not used as hanging node global id can be
       ! identified for a given element only
       ! msh%htf - not used as faces global id is already provided by p4est
       ! msh%hte - not used as edge global id is already provided by p4est

       ! Fill in neighbour information
       call p4_neighbour_fill(msh, this)

       ! Fill in distdata information
       call p4_distdata_fill(msh, this)

       ! Fill in boundary condition information
       call p4_bc_fill(msh, this)


       ! curved edges
       !! @todo curves for p4est
       call msh%curve%init(msh%nelv)
       ! P4EST DOES NOT PROVIDE CURVATURE DATA DIRECTLY,
       ! SO LET'S LEAVE IT FOR NOW
       ! this data should be loaded from another file??
       !call mesh_mark_curve_element(msh, el_idx, curve_data, type)

       ! fill connectivity information
       call p4_connect_fill(msh, this)

       ! mesh finalisation
       call mesh_generate_flags(msh) ! deformation flag;
       ! THIS IS NOT 100% CORRECT

       call msh%wall%finalize() ! BC finalisation
       call msh%inlet%finalize()
       call msh%outlet%finalize()
       call msh%outlet_normal%finalize()
       call msh%sympln%finalize()
       call msh%periodic%finalize()
       do il = 1, NEKO_MSH_MAX_ZLBLS
          call msh%labeled_zones(il)%finalize()
       end do
       call msh%curve%finalize() ! curved sides finalisation
    end if

    msh%lconn = .true.
    msh%lnumr = .true.

    call neko_log%end_section()

    return
  end subroutine p4_msh_get

  !> perform refinement based on refinement flag
  subroutine p4_refine(this, ref_mark, el_gidx, msh_trs, level_max,&
       &ifmod, msh_rcn)
    ! argument list
    class(p4_mesh_cnstr_t), intent(inout) :: this
    integer(i4), dimension(:), intent(in) :: ref_mark, el_gidx
    type(mesh_manager_transfer_t), intent(in) :: msh_trs
    integer(i4), intent(in) :: level_max
    logical, intent(out) :: ifmod
    type(mesh_reconstruct_transfer_t), intent(out) :: msh_rcn
    ! local variables
    integer(i4) :: il, itmp
    integer(i4), target, allocatable, dimension(:) :: pref_mark, pel_gnum,&
         & pel_lnum, pel_nid
    integer(i4), allocatable, dimension(:) :: nel_gnum, nel_lnum, nel_nid
    ! element restructure data
    integer(i4) :: map_nr, rfn_nr, crs_nr
    integer(i4), target, allocatable, dimension(:, :) :: elgl_map, elgl_rfn
    integer(i4), target, allocatable, dimension(:, :, :) :: elgl_crs

    ! THIS SHOULD BE PROBABLY MOVED FROM HERE
    ! partition strategy for p4est
    ! 0 - equal element count
    ! 1 - octant families prepared for coarsening
    integer(i4), parameter :: p4part = 1

    ! copy parameter for forest comparison in p4est
    ! 0 - no check of quadrant data
    ! 1 - check if quadrant data are identical
    integer(i4), parameter :: p4test = 1

    itmp = size(ref_mark)
    if (itmp /= size(el_gidx)) call neko_error('Inconsistent array size')

    ! put refinement mark in p4est
    allocate(pref_mark(this%elem%nelv))
    ! PLACE FOR DATA DISTRIBUTION (NEKO => P4EST) using
    ! mesh_manager_transfer_t
    ! in general itmp can be different than this%elem%nelv, but for now
    ! they are the same
    pref_mark(1:itmp) = ref_mark(1:itmp)
    ! PLACE FOR LOCAL CONSISTENCY CHECKS
    call wp4est_refm_put(c_loc(pref_mark))

    ! put element distribution info in p4est
    allocate(pel_gnum(this%elem%nelv), pel_lnum(this%elem%nelv),&
         & pel_nid(this%elem%nelv))
    allocate(nel_gnum(itmp), nel_lnum(itmp), nel_nid(itmp))
    !
    do il = 1, itmp
       nel_gnum(il) = el_gidx(il)
       nel_lnum(il) = il
    end do
    nel_nid(:) = pe_rank
    ! PLACE FOR DATA DISTRIBUTION (NEKO=> P4EST) using
    ! mesh_manager_transfer_t
    ! in general itmp can be different than this%elem%nelv, but for now
    ! they are the same
    pel_gnum(1:itmp) = nel_gnum(1:itmp)
    pel_lnum(1:itmp) = nel_lnum(1:itmp)
    pel_nid(1:itmp) = nel_nid(1:itmp)
    ! this possibly could be reconstructed from mesh_manager_transfer_t,
    ! so no communication would be necessary
    call wp4est_egmap_put(c_loc(pel_gnum), c_loc(pel_lnum), c_loc(pel_nid))

    ! perform local refine/coarsen/balance on p4est side
    call wp4est_tree_copy(p4test)
    call wp4est_refine(level_max)
    call wp4est_coarsen()
    call wp4est_balance()
    ! perform partitioning on p4est side
    call wp4est_part(p4part)
    call wp4est_tree_check(itmp, p4test)

    if (itmp == 0 ) then
       ifmod = .true.
    else
       ifmod = .false.
    end if

    ! if mesh was modified extract data to restructure elements
    if (ifmod) then
       itmp = 2**this%dim ! number of vertices
       allocate(elgl_map(3, itmp*this%elem%nelv), &
            & elgl_rfn(5, itmp*this%elem%nelv), &
            & elgl_crs(4, itmp, this%elem%nelv))
       call wp4est_msh_get_hst(map_nr, rfn_nr, crs_nr, c_loc(elgl_map), &
            & c_loc(elgl_rfn), c_loc(elgl_crs))
       ! move data
       msh_rcn%map_nr = map_nr
       msh_rcn%rfn_nr = rfn_nr
       msh_rcn%crs_nr = crs_nr
       call MOVE_ALLOC(elgl_map, msh_rcn%elgl_map)
       call MOVE_ALLOC(elgl_rfn, msh_rcn%elgl_rfn)
       call MOVE_ALLOC(elgl_crs, msh_rcn%elgl_crs)
    end if

    ! free memory
    deallocate(pref_mark, pel_gnum, pel_lnum, pel_nid, nel_gnum,&
         & nel_lnum, nel_nid)

    return
  end subroutine p4_refine

  ! Import data from p4est and store in a local type
  subroutine p4_mesh_import_data(p4)
    ! argument list
    type(p4_mesh_cnstr_t), intent(inout) :: p4
    ! local variables
    character(len=LOG_SIZE) :: log_buf
    integer :: ierr, nvert, nface, nedge
    integer(i8) :: itmp8
    integer(i4), allocatable, target, dimension(:) :: itmp4v1, itmp4v2, itmp4v3
    integer(i4), allocatable, target, dimension(:,:) :: itmp4v21, itmp4v22,&
         & itmp4v23
    integer(i8), allocatable, target, dimension(:) :: itmp8v1
    real(dp), allocatable, target, dimension(:,:) :: rtmpv1
    real(dp), allocatable, target, dimension(:,:,:) :: rtmpv2
    character(len=*), parameter ::&
         & frmt="('gdim = ', i1, ', nelements =', i9,', max ref. lev. = ',i2)"

       call neko_log%message('Import mesh data from p4est')
       ! clean memory
       call p4%free()

       ! create p4est ghost zones, mesh and nodes
       call wp4est_ghost_new()
       call wp4est_mesh_new()
       call wp4est_nodes_new()

       ! check boundary condition mark
       call wp4est_bc_check()

       associate(dim=>p4%dim, nelv=>p4%elem%nelv)
         ! get mesh size and distribution information
         call wp4est_msh_get_size(p4%dim,p4%elem%nelgt, p4%elem%nelgto, &
              & p4%elem%nelt, p4%elem%nelv, p4%maxl)
         ! get max refinement level across all ranks
         call MPI_Allreduce(p4%maxl, p4%maxg, 1, &
              MPI_INTEGER, MPI_MAX, NEKO_COMM, ierr)

         write(log_buf,frmt) dim, p4%elem%nelgt, p4%maxg
         call neko_log%message(log_buf)

         if (nelv > 0) then
            ! get nodes and their coordinates
            nvert = 2**dim
            nface = 2*dim
            call wp4est_nds_get_size(p4%indn%lown, p4%indn%lshr, &
                 & p4%indn%loff, p4%indn%lnum, p4%fchn%lnum, &
                 & p4%edhn%lnum)

            ! get independent node info
            if (p4%indn%lnum > 0) then
               allocate(itmp8v1(p4%indn%lnum), &
                    & itmp4v1(p4%indn%lnum), rtmpv1(dim, &
                    & p4%indn%lnum))
               call wp4est_nds_get_ind(c_loc(itmp8v1), &
                    & c_loc(itmp4v1),c_loc(rtmpv1))
               call MOVE_ALLOC(itmp8v1, p4%indn%gidx)
               call MOVE_ALLOC(itmp4v1, p4%indn%ndown)
               call MOVE_ALLOC(rtmpv1, p4%indn%coord)
            end if

            ! get hanging node info; for refined mesh only
            if (p4%maxl > 0) then
               if (dim == 2) then
                  if (p4%fchn%lnum > 0) then
                     ! each face hanging node is defined by 4 independent nodes
                     allocate(itmp4v21(2, p4%fchn%lnum), &
                          & rtmpv1(dim, p4%fchn%lnum))
                     call wp4est_nds_get_hfc(c_loc(itmp4v21), c_loc(rtmpv1))
                     call MOVE_ALLOC(itmp4v21, p4%fchn%lmap)
                     call MOVE_ALLOC(rtmpv1, p4%fchn%coord)
                  end if
               else if (dim == 3) then
                  if (p4%fchn%lnum > 0) then
                     ! each face hanging node is defined by 4 independent nodes
                     allocate(itmp4v21(4, p4%fchn%lnum), &
                          & rtmpv1(p4%dim, p4%fchn%lnum))
                     call wp4est_nds_get_hfc(c_loc(itmp4v21), c_loc(rtmpv1))
                     call MOVE_ALLOC(itmp4v21, p4%fchn%lmap)
                     call MOVE_ALLOC(rtmpv1, p4%fchn%coord)
                  end if
                  if (p4%edhn%lnum > 0) then
                     ! each edge hanging node is defined by 2 independent nodes
                     allocate(itmp4v21(2, p4%edhn%lnum), &
                          & rtmpv1(p4%dim, p4%edhn%lnum))
                     call wp4est_nds_get_hed(c_loc(itmp4v21), c_loc(rtmpv1))
                     call MOVE_ALLOC(itmp4v21, p4%edhn%lmap)
                     call MOVE_ALLOC(rtmpv1, p4%edhn%coord)
                  end if
               end if
            end if

            ! get element vertex mapping to nodes
            allocate(itmp4v21(nvert, nelv))
            call wp4est_nds_get_vmap(c_loc(itmp4v21))
            call MOVE_ALLOC(itmp4v21, p4%elem%vnmap)

            ! get element info
            allocate(itmp8v1(nelv), itmp4v1(nelv), itmp4v2(nelv), &
                 & itmp4v21(nface, nelv), itmp4v22(nface, nelv), &
                 & rtmpv2(dim, nvert, nelv), itmp4v23(nface, nelv))
            call wp4est_elm_get_dat(c_loc(itmp8v1), c_loc(itmp4v1), &
                 & c_loc(itmp4v2), c_loc(itmp4v21), c_loc(itmp4v22), &
                 & c_loc(rtmpv2), c_loc(itmp4v23))
            call MOVE_ALLOC(itmp8v1, p4%elem%gidx)
            call MOVE_ALLOC(itmp4v1, p4%elem%level)
            call MOVE_ALLOC(itmp4v2, p4%elem%igrp)
            call MOVE_ALLOC(itmp4v21, p4%elem%crv)
            call MOVE_ALLOC(itmp4v22, p4%elem%bc)
            call MOVE_ALLOC(itmp4v23, p4%elem%falg)

            ! get periodic nodes
            call p4_node_periodic_get(p4, rtmpv2, dim, nvert, nelv)
            deallocate(rtmpv2)

            call wp4est_lnodes_new(1)
            ! get hanging object info; based on lnode information
            allocate(itmp4v1(nelv), itmp4v21(nface, nelv))
            if (dim == 3) then
               nedge = 12
            else
               nedge = 0  !!!!! THIS SHOULD BE CHECKED FOR 2D SIMULATION !!!!!
            end if
            allocate(itmp4v22(nedge, nelv))
            call wp4est_hang_get_info(c_loc(itmp4v1), c_loc(itmp4v21), &
                 & c_loc(itmp4v22))
            call MOVE_ALLOC(itmp4v1, p4%elem%hngel)
            call MOVE_ALLOC(itmp4v21, p4%elem%hngfc)
            call MOVE_ALLOC(itmp4v22, p4%elem%hnged) &
                 &!!!!! THIS SHOULD BE CHECKED FOR 2D SIMULATION !!!

            ! get global indexes of element vertices
            allocate(itmp4v21(nvert, nelv))
            call wp4est_elm_get_lnode(p4%elem%vert%lnum, &
                 & p4%elem%vert%lown, p4%elem%vert%goff, &
                 & c_loc(itmp4v21))
            call MOVE_ALLOC(itmp4v21, p4%elem%vert%lmap)
            call wp4est_sharers_get_size(p4%elem%vert%nrank, &
                 & p4%elem%vert%nshare)
            allocate(itmp8v1(p4%elem%vert%lnum), &
                 & itmp4v1(p4%elem%vert%nrank), &
                 & itmp4v2(p4%elem%vert%nrank+1), &
                 & itmp4v3(p4%elem%vert%nshare))
            call wp4est_sharers_get_ind(c_loc(itmp8v1), c_loc(itmp4v1), &
                 & c_loc(itmp4v2), c_loc(itmp4v3))
            call MOVE_ALLOC(itmp8v1, p4%elem%vert%lgidx)
            call MOVE_ALLOC(itmp4v1, p4%elem%vert%lrank)
            call MOVE_ALLOC(itmp4v2, p4%elem%vert%loff)
            call MOVE_ALLOC(itmp4v3, p4%elem%vert%lshare)

            ! get global number of vertices
            itmp8 = p4%elem%vert%lown
            call MPI_Allreduce(itmp8, p4%elem%vert%gnum, 1, &
                 MPI_INTEGER8, MPI_SUM, NEKO_COMM, ierr)
            ! get global indexes of element faces
            call wp4est_lnodes_del()

            ! get global indexes of element faces
            call wp4est_lnodes_new(-1)
            allocate(itmp4v21(nface, nelv))
            call wp4est_elm_get_lnode(p4%elem%face%lnum, p4%elem%face%lown, &
                 & p4%elem%face%goff, c_loc(itmp4v21))
            call MOVE_ALLOC(itmp4v21, p4%elem%face%lmap)
            call wp4est_sharers_get_size(p4%elem%face%nrank, &
                 & p4%elem%face%nshare)
            allocate(itmp8v1(p4%elem%face%lnum), &
                 & itmp4v1(p4%elem%face%nrank), &
                 & itmp4v2(p4%elem%face%nrank+1), &
                 & itmp4v3(p4%elem%face%nshare))
            call wp4est_sharers_get_ind(c_loc(itmp8v1), c_loc(itmp4v1), &
                 & c_loc(itmp4v2), c_loc(itmp4v3))
            call MOVE_ALLOC(itmp8v1, p4%elem%face%lgidx)
            call MOVE_ALLOC(itmp4v1, p4%elem%face%lrank)
            call MOVE_ALLOC(itmp4v2, p4%elem%face%loff)
            call MOVE_ALLOC(itmp4v3, p4%elem%face%lshare)

            ! get global number of faces
            itmp8 = p4%elem%face%lown
            call MPI_Allreduce(itmp8, p4%elem%face%gnum, 1, &
                 MPI_INTEGER8, MPI_SUM, NEKO_COMM, ierr)
            call wp4est_lnodes_del()

            ! get edge information; this is not straightforward, as p4est
            ! does not provide a simple way to extract this information
            if (dim == 3) call p4_edge_get(p4)

            ! get family information
            call p4_family_get(p4)

         end if
       end associate

       ! destroy p4est nodes, mesh and ghost cells
       call wp4est_nodes_del()
       call wp4est_mesh_del()
       call wp4est_ghost_del()

    return
  end subroutine p4_mesh_import_data

  ! this subroutine should be adjusted for 2D; not done yet
  subroutine p4_node_periodic_get(p4, vcoord, dim, nvert, nelv)
    ! argument list
    type(p4_mesh_cnstr_t), intent(inout) :: p4
    integer, intent(in) :: dim, nvert, nelv
    real(dp), dimension(dim, nvert, nelv), intent(in) :: vcoord
    !local variables
    integer :: il, jl, kl, itmp
    integer :: nvper
    integer :: nin, nhf, nhe, nvt
    logical :: ifvequal
    type(point_t) :: ptsv, ptsn
    type(tuple_i4_t) :: idx, tmp
    type(htable_pt_t) :: htpts
    integer(i4), allocatable, dimension(:,:) :: pvmap
    integer(i4), allocatable, dimension(:) :: pnmap
    real(dp), allocatable, dimension(:,:) :: pncoord

    ! find a number of vertices not included in node arrays
    nvper = 0
    ! get shifts in mapping array
    nin = p4%indn%lnum
    nhf = nin + p4%fchn%lnum
    nhe = nhf + p4%edhn%lnum
    ! initialise htable
    call htpts%init(4, idx)
    allocate(pvmap(nvert, nelv), pnmap(nvert*nelv), pncoord(dim, nvert*nelv))
    pvmap = -1
    do il = 1, nelv
       do jl = 1, nvert
          ! get element vertex coordinates
          ptsv = point_t(vcoord(:, jl, il))
          ! get mapped node coordinates
          nvt = p4%elem%vnmap(jl, il)
          if (nvt <= nin) then
             ptsn = point_t(p4%indn%coord(:, nvt))
          else if (nvt <= nhf) then
             ptsn = point_t(p4%fchn%coord(:, nvt - nin))
          else if (nvt <= nhe) then
             ptsn = point_t(p4%edhn%coord(:, nvt - nhf))
          else
             call neko_error('Inconsistent vnmap value.')
          end if
          ifvequal = ptsv.ne.ptsn
          if (ifvequal) then
             ! check if this vertex belongs to periodic bc
             ifvequal = .false.
             do kl = 1, dim
                ifvequal = p4%elem%bc(p4_cface(kl, jl), il) == -1
                if (ifvequal) then
                   ! hanging nodes should not be marked periodic
                   if (nvt > nin) &
                        & call neko_error('Hangign node marked as periodic.')
                   if (htpts%get(ptsv, tmp) > 0) then
                      ! new node
                      nvper = nvper + 1
                      idx = (/ nvper, nvt /)
                      call htpts%set(ptsv, idx)
                      call ptsv%set_id(nvper)
                      pvmap(jl, il) = nvper ! vertex to periodic node mapping
                      pnmap(nvper) = nvt ! periodic node to independent node mapping
                      pncoord(1:dim, nvper) = ptsv%x(1:dim) ! periodic node coordinates
                   else
                      ! exisitng node
                      ! check mapping consistency
                      if (nvt /= tmp%x(2)) &
                           & call neko_error('Inconsistent mapping of per. nodes to indep. nodes.')
                      pvmap(jl, il) = tmp%x(1) ! vertex to periodic node mapping
                   end if
                   exit
                end if
             end do
             if (.not.ifvequal) then
                call neko_error('Periodic node not a vertex of periodic face')
             end if
          end if
       end do
    end do

    ! put data to the type
    p4%pern%lnum = nvper
    allocate(p4%pern%lmap(nvper), p4%pern%coord(dim, nvper))
    p4%pern%lmap(1:nvper) = pnmap(1:nvper)
    p4%pern%coord(:,1:nvper) = pncoord(:, 1:nvper)

    ! update vertex to node mapping array
    do il = 1, nelv
       do jl = 1, nvert
          if (pvmap(jl, il) == -1) then
             if (p4%elem%vnmap(jl, il) > nin) &
                  & p4%elem%vnmap(jl, il) = p4%elem%vnmap(jl, il) + nvper
          else
             p4%elem%vnmap(jl, il) = nin + pvmap(jl, il)
          end if
       end do
    end do

    deallocate(pvmap, pnmap, pncoord)
    call htpts%free()

    return
  end subroutine p4_node_periodic_get

  ! This routine provides global edge numbering and orientation
  subroutine p4_edge_get(p4)
    ! argument list
    type(p4_mesh_cnstr_t), intent(inout) :: p4
    !local variables
    integer, parameter :: nface = 6, nedge = 12
    integer :: il, jl, kl, itmp, ierr, nlnode, last, nlown
    integer(i8) :: glnode, goff, itmp8
    integer(i4) :: felnum, felown, fenrank, fenshare
    integer(i8) :: feloff
    integer(i4), allocatable, target, dimension(:,:) :: felmap ! global index of element faces and edges
    integer(i4), allocatable, target, dimension(:) :: felrank, feloffs, felshare, feind
    integer(i8), allocatable, target, dimension(:) :: felgidx
    integer(i4), allocatable, dimension(:) :: esort, eind, eoffset ! sorting arrays for edges
    integer(i4), allocatable, dimension(:) :: eloffs
    integer(i4), allocatable, dimension(:,:) :: rbuf, sbuf ! send/receive buffers
    type(MPI_Request), allocatable, dimension(:) :: request ! MPI request
    type(MPI_Status), allocatable, dimension(:) :: status ! MPI status

    if (p4%dim == 3) then
       associate(nelv=>p4%elem%nelv, lnum=>p4%elem%edge%lnum, lown=>p4%elem%edge%lown, &
            & nrank=>p4%elem%edge%nrank, nshare=>p4%elem%edge%nshare)
         ! import from p4est combined face and edge information
         allocate(felmap(nface+nedge, nelv))
         call wp4est_lnodes_new(-2)
         call wp4est_elm_get_lnode(felnum, felown, feloff, c_loc(felmap))
         call wp4est_sharers_get_size(fenrank, fenshare)
         allocate(felgidx(felnum), felrank(fenrank), feloffs(fenrank+1), &
              & felshare(fenshare), feind(fenshare))
         call wp4est_sharers_get_ind(c_loc(felgidx), c_loc(felrank), c_loc(feloffs), &
              & c_loc(felshare))
         call wp4est_lnodes_del()

         ! extract local edge mapping
         allocate(esort(nedge*nelv), eind(nedge*nelv), eoffset(nedge*nelv+1))
         do il = 1, nelv
            do jl = 1, nedge
               esort((il-1)*nedge + jl) = felmap(nface + jl, il)
            end do
         end do
         ! sort edge mapping
         itmp = nedge*nelv
         call sorti4(esort, eind, itmp)
         ! compress the list removing multiplicities
         nlnode = 1
         last = esort(1)
         eoffset(1) = 1
         do il = 2, itmp
            if (last /= esort(il)) then
               last = esort(il)
               nlnode = nlnode + 1
               eoffset(nlnode) = il
               esort(nlnode) = last
            end if
         end do
         eoffset(nlnode+1) = itmp + 1
         ! find a number of owned edges
         nlown = 0
         do il = 1, nlnode
            if (esort(il) <= felown) then
               nlown = il
            else
               exit
            end if
         end do
         ! get global number of unique (owned) edges
         itmp8 = nlown
         call MPI_Allreduce(itmp8, glnode, 1, MPI_INTEGER8, MPI_SUM, &
              & NEKO_COMM, ierr)
         ! get global offset
         call MPI_Scan(itmp8, goff, 1, MPI_INTEGER8, MPI_SUM, &
              & NEKO_COMM, ierr)
         goff = goff - itmp8

         ! reduce sharing information excluding all the face information
         ! take advantage of the fact both arrays esort and sections of felshare are sorted
         allocate(eloffs(fenrank+1))
         fenshare = 1
         eloffs(1) = fenshare
         do il=1, fenrank ! mpi rank loop
            ! unfortunately the list of the shared nodes on a local rank seems to be
            ! a concatenated list rahter than a sorted one;
            ! that is why for a local rank perform sorting
            if (felrank(il) == pe_rank) then
               itmp = feloffs(il + 1) - feloffs(il)
               call sorti4(felshare(feloffs(il):feloffs(il + 1) - 1), &
                    & feind(feloffs(il):feloffs(il + 1) - 1), itmp)
            end if
            jl = feloffs(il)
            kl = 1
            do
               if (esort(kl) == felshare(jl)) then
                  felshare(fenshare) = kl ! point to the new local edge position
                  fenshare = fenshare + 1
                  jl = jl + 1
               else if (esort(kl) > felshare(jl)) then
                  jl = jl + 1
               else
                  kl = kl + 1
               end if
               if ((jl == feloffs(il+1)).or.(kl > nlnode)) then
                  eloffs(il+1) = fenshare
                  exit
               end if
            end do
         end do
         ! correct number of entries
         fenshare = fenshare-1

         ! start communication
         allocate(rbuf(2, fenshare), sbuf(2, fenshare), request(fenrank), status(fenrank))
         ! set non-blocking receive counting ranks
         itmp = 0
         kl = 0
         do il=1, fenrank ! mpi rank loop
            if (eloffs(il) < eloffs(il+1)) then
               itmp = itmp + 1 ! count ranks
               if (felrank(il) < pe_rank) then ! receive from ranks with lower id only
                  kl = kl + 1 ! count messages
                  jl = 2*(eloffs(il+1) - eloffs(il))
                  call MPI_IRecv(rbuf(:, eloffs(il):eloffs(il+1)-1), jl, MPI_INTEGER, &
                       & felrank(il), jl, NEKO_COMM, request(kl), ierr)
               end if
            end if
         end do

         ! renumber owned nodes and start collecting information
         p4%elem%edge%lnum = nlnode
         p4%elem%edge%lown = nlown
         p4%elem%edge%gnum = glnode
         p4%elem%edge%goff = goff
         p4%elem%edge%nrank = itmp
         p4%elem%edge%nshare = fenshare

         ! allocate arrays
         allocate(p4%elem%edge%lmap(nedge, nelv), p4%elem%edge%lgidx(lnum), &
              & p4%elem%edge%lrank(nrank), p4%elem%edge%loff(nrank+1), &
              & p4%elem%edge%lshare(nshare))

         ! copy shared rank and offset
         itmp = 0
         p4%elem%edge%loff(1) = 1
         do il=1, fenrank
            if (eloffs(il) < eloffs(il+1)) then
               itmp = itmp + 1 ! count ranks
               p4%elem%edge%lrank(itmp) = felrank(il)
               p4%elem%edge%loff(itmp+1) = eloffs(il+1)
            end if
         end do
         if (p4%elem%edge%loff(nrank+1) /= nshare + 1) &
              & call neko_error('Invalid shared offset for edges')
         ! copy shared nodes
         p4%elem%edge%lshare(1:nshare) = felshare(1:nshare)

         ! renumber owned edges
         do il= 1, lown
            p4%elem%edge%lgidx(il) = p4%elem%edge%goff +  il
         end do
         ! mark not owned edges
         do il= lown + 1, lnum
            p4%elem%edge%lgidx(il) = -1
         end do

         ! redistribute information to get numberring of not owned nodes
         ! for safety reason I store in a buffer an old and a new global number
         itmp = 0
         do il=1, nrank
            if (p4%elem%edge%lrank(il) > pe_rank) then ! send to ranks with higher id only
               itmp = itmp + 1 ! count messages
               do jl = p4%elem%edge%loff(il), p4%elem%edge%loff(il+1) -1
                  sbuf(1,jl) = felgidx(esort(p4%elem%edge%lshare(jl))) ! old global number
                  sbuf(2,jl) = p4%elem%edge%lgidx(p4%elem%edge%lshare(jl))  ! new global number
               end do
               jl = 2*(p4%elem%edge%loff(il+1) - p4%elem%edge%loff(il))
               call MPI_Send(sbuf(:, p4%elem%edge%loff(il):p4%elem%edge%loff(il+1)-1), &
                    & jl, MPI_INTEGER, p4%elem%edge%lrank(il), jl, NEKO_COMM, ierr)
            end if
         end do

         ! reconstruct edge mapping
         deallocate(eloffs)
         allocate(eloffs(nedge*nelv))
         ! inflate sorted mapping
         do il = 1, lnum
            do jl = eoffset(il), eoffset(il+1) - 1
               eloffs(jl) = il
            end do
         end do
         ! sort it back
         itmp = nedge*nelv
         call reordi4(eloffs, eind, itmp)
         ! copy data
         do il = 1, nelv
            do jl = 1, nedge
               p4%elem%edge%lmap(jl, il) = eloffs((il-1)*nedge + jl)
            end do
         end do

         ! finalize communication
         ! this could be more fancy, but most probably the gain would be minimal
         call MPI_Waitall(kl, request, status, ierr)

         ! fill in missing global id for shared edges
         do il=1, nrank
            if (p4%elem%edge%lrank(il) < pe_rank) then ! conside ranks with lower id only
               do jl = p4%elem%edge%loff(il), p4%elem%edge%loff(il+1) -1
                  if (p4%elem%edge%lshare(jl) > lnum) then ! wrong mapping
                     call neko_error('Invalid shared mapping for edges; > lnum')
                  else if (p4%elem%edge%lshare(jl) > lown) then ! shared nodes
                     if (p4%elem%edge%lgidx(p4%elem%edge%lshare(jl)) == -1) then ! renumber
                        if (felgidx(esort(p4%elem%edge%lshare(jl))) == &
                             & rbuf(1, jl)) then ! correct old global id
                           p4%elem%edge%lgidx(p4%elem%edge%lshare(jl)) = &
                                & rbuf(2, jl)
                        else
                           call neko_error('Invalid old global id for edges; 1')
                        end if
                     else ! node already set
                        if (felgidx(esort(p4%elem%edge%lshare(jl))) == &
                             & rbuf(1, jl)) then ! correct old global id
                           if (rbuf(2, jl) /= -1) &
                                & call neko_error('Invalid global id for shared edges; this one shoud not be set yet')
                        else
                           call neko_error('Invalid old global id for edges; 2')
                        end if
                     end if
                  else ! owned nodes; lower mpi rank shoud own the nodes
                     call neko_error('Invalid shared mapping for edges; <= lown')
                  end if
               end do
            end if
         end do

         ! check if all the ids are set properly
         do il = lown + 1, lnum
            if (p4%elem%edge%lgidx(il) == -1) &
                 & call neko_error('Invalid old global id for edges; shared id not set')
            if (p4%elem%edge%lgidx(il) >= p4%elem%edge%lgidx(1)) &
                 & call neko_error('Invalid old global id for edges; shared id bigger than owned')
         end do

         ! get global edge orientation based on vertex global number
         ! this is just an easy (dirty) and temporary hack, and will be changed in the future !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         allocate(p4%elem%ealg(nedge, nelv))
         do il = 1, nelv
            do jl = 1, nedge
               if (p4%elem%vert%lmap(p4_vedge(1, jl), il) < &
                    & p4%elem%vert%lmap(p4_vedge(2, jl), il)) then
                  p4%elem%ealg(jl, il) = 0
               else
                  p4%elem%ealg(jl, il) = 1
               end if
            end do
         end do

         deallocate(felmap, felgidx, felrank, feloffs, felshare, feind, &
              & esort, eind, eoffset, eloffs)
         deallocate(rbuf, sbuf, request, status)
       end associate
    end if

    return
  end subroutine p4_edge_get

  ! This routine provides information about sets of elements that share
  ! the same parent and could be destroyed togeter during coarsening step.
  subroutine p4_family_get(p4)
    ! argument list
    type(p4_mesh_cnstr_t), intent(inout) :: p4
    !local variables
    integer(i4) :: nvert, ierr, il, jl
    integer(i8) :: itmp8
    integer(i4) :: nlfam ! local number of families
    integer(i8) :: nlfam_off ! global family nuimber offset
    integer(i8), allocatable, target, dimension(:) :: itmp4v1

    ! import family mark
    associate(dim=>p4%dim, nelv=>p4%elem%nelv)
      allocate(itmp4v1(nelv))
      call wp4est_fml_get_info(c_loc(itmp4v1), nlfam)

      ! test number of families; it must be multiply of number of vertices
      nvert = 2**dim
      if ((nlfam > nelv).or.(mod(nlfam, nvert) /= 0)) &
           & call neko_error('Invalid number of families.')

      ! get global family offset
      nlfam = nlfam/nvert
      itmp8 = nlfam
      call MPI_Scan(itmp8, nlfam_off, 1, MPI_INTEGER8, MPI_SUM, NEKO_COMM, ierr)
      nlfam_off = nlfam_off - itmp8

      ! mark all local families
      allocate(p4%elem%fmlm(2, nelv))
      p4%elem%fmlm = 0
      do il=1, nlfam
         nlfam_off = nlfam_off + 1
         do jl = 1, nvert
            itmp8 = itmp4v1(jl + (il-1)*nvert) - p4%elem%nelgto
            ! one could test if 0<itmp8<=nelv
            p4%elem%fmlm(1, int(itmp8,i4)) = nlfam_off
            p4%elem%fmlm(1, int(itmp8,i4)) = nvert + 1 - jl
         end do
      end do
    end associate

    deallocate(itmp4v1)

    return
  end subroutine p4_family_get

  ! Fill the mesh type with node information from p4_mesh_cnstr type
  subroutine p4_nodes_fill(msh, p4)
    ! argument list
    type(mesh_t), intent(inout) :: msh
    type(p4_mesh_cnstr_t), intent(in) :: p4
    !local variables
    integer(i4) :: il
    integer(i4) :: itmp
    integer(i4) :: nin, npe, nhf, nhe

    ! Fill in with all unique local (not global) points.
    ! The point id is set to a local id, as at this stage the mapping to the global id is not unique
    ! get shifts in vertex mapping array
    nin = p4%indn%lnum
    npe = p4%pern%lnum
    nhf = p4%fchn%lnum
    nhe = p4%edhn%lnum
    allocate(msh%points(nin+npe+nhf+nhe))
    ! independent nodes; the only ones with direct global mapping
    do il = 1, nin
       itmp = nin + il ! this because point constructor has id with intent(inout) attribute
       msh%points(il) = point_t(p4%indn%coord(:, il), itmp)
    end do
    ! periodic nodes; global id through mapping to independent nodes (independent on element position)
    do il = 1, npe
       itmp = nin + il ! this because point constructor has id with intent(inout) attribute
       msh%points(itmp) = point_t(p4%pern%coord(:, il), itmp)
    end do
    npe = npe + nin
    ! face hanging node; global mapping depend on the element position
    do il = 1, nhf
       itmp = npe + il
       msh%points(itmp) = point_t(p4%fchn%coord(:, il), itmp)
    end do
    nhf = nhf + npe
    if (msh%gdim == 3) then
       ! edge hanging node; global mapping depend on the element position
       do il = 1, nhe
          itmp = nhf + il
          msh%points(itmp) = point_t(p4%edhn%coord(:, il), itmp)
       end do
    end if

    return
  end subroutine p4_nodes_fill

  ! Fill the mesh type with element information from p4_mesh_cnstr type
  subroutine p4_element_fill(msh, p4)
    ! argument list
    type(mesh_t), intent(inout) :: msh
    type(p4_mesh_cnstr_t), intent(in) :: p4
    !local variables
    integer(i4) :: il
    integer(i4) :: itmp
    class(element_t), pointer :: ep

    ! Fill in element array
    ! Mapping of the vertex to local node is sufficient to get unique global,
    ! so I could set point id to the global one here, but entity id is just int4
    ! and hex/quad uses pointers
    allocate(msh%elements(msh%nelv))
    if (msh%gdim == 3) then
       do il = 1, msh%nelv
          allocate(hex_t::msh%elements(il)%e)
          ! get global element id; NOTE there is type casting from int8
          itmp = int(p4%elem%gidx(il),i4)
          ! set hex in the array using vertex to node mapping
          select type(ep => msh%elements(il)%e)
          type is (hex_t)
             call ep%init(itmp, &
                  & msh%points(p4%elem%vnmap(1, il)), &
                  & msh%points(p4%elem%vnmap(2, il)), &
                  & msh%points(p4%elem%vnmap(3, il)), &
                  & msh%points(p4%elem%vnmap(4, il)), &
                  & msh%points(p4%elem%vnmap(5, il)), &
                  & msh%points(p4%elem%vnmap(6, il)), &
                  & msh%points(p4%elem%vnmap(7, il)), &
                  & msh%points(p4%elem%vnmap(8, il)))
          class default
             call neko_error('Invalid element type')
          end select
       end do
    else if (msh%gdim == 2) then
       call neko_error("2D not supported for p4est yet")
    else
       call neko_error("Invalid dimension")
    end if

    ! deformed element flag; It is set in mesh_generate_flags (dependent on vertex ordering)
    allocate(msh%dfrmd_el(msh%nelv))

    return
  end subroutine p4_element_fill

  ! Fill the mesh type with neighbour information from p4_mesh_cnstr type
  subroutine p4_neighbour_fill(msh, p4)
    ! argument list
    type(mesh_t), intent(inout) :: msh
    type(p4_mesh_cnstr_t), intent(in) :: p4

    ! get vertex neighbours; not used
    !call p4_vertex_neighbour_fill(msh, p4)

    ! USED IN gather_scatter.f90; MORE INVESTIGATION NEEDED
    ! get face neighbours; is it reall y needed?
    call p4_face_neighbour_fill(msh, p4)

    ! USED IN gather_scatter.f90; MORE INVESTIGATION NEEDED
    ! get neighbour rank info; based on vertex sharing info
    call p4_rank_neighbour_fill(msh, p4)

    return
  end subroutine p4_neighbour_fill

  ! Fill the mesh type with vertex neighbour information
  subroutine p4_vertex_neighbour_fill(msh, p4)
    ! argument list
    type(mesh_t), intent(inout) :: msh
    type(p4_mesh_cnstr_t), intent(in) :: p4
    !local variables
    integer(i4) :: il, jl, kl, ll, ml
    integer(i4) :: itmp, ierr
    integer(i4), allocatable, dimension(:) :: cmoff ! offest in send/receive buffers
    integer(i4) :: neighn
    integer(i4), pointer, dimension(:) :: neighl
    integer(i4), allocatable, dimension(:) :: rbuf, sbuf ! send/receive buffers
    type(MPI_Request), allocatable, dimension(:) :: request ! MPI request
    type(MPI_Status), allocatable, dimension(:) :: status ! MPI status

    ! Get list of vertex neighbours (elements a given vertex belongs to)
    ! local vertex neighbours
    allocate(msh%point_neigh(msh%mpts))
    do il = 1, msh%mpts
       call msh%point_neigh(il)%init()
    end do
    do il = 1, msh%nelv
       ! get global element id; NOTE there is type casting from int8
       itmp = int(p4%elem%gidx(il), i4)
       do jl = 1, msh%npts
          call msh%point_neigh(p4%elem%vert%lmap(jl, il))%push(itmp)
       end do
    end do

    ! nonlocal node neighbours
    ! count size of send/receive buffers
    allocate(cmoff(p4%elem%vert%nrank), request(p4%elem%vert%nrank-1), &
         & status(p4%elem%vert%nrank-1))
    itmp = 1
    cmoff(itmp) = 1
    do il= 1, p4%elem%vert%nrank ! mpi rank loop
       if (p4%elem%vert%lrank(il) /= pe_rank) then
          itmp = itmp + 1 ! mpi rank loop
          cmoff(itmp) = cmoff(itmp - 1)
          do jl = p4%elem%vert%loff(il), p4%elem%vert%loff(il+1) -1
             ! sum number of local vertex neighbours for a given rank
             neighn = msh%point_neigh(p4%elem%vert%lshare(jl))%size()
             cmoff(itmp) = cmoff(itmp) + 2 + neighn
          end do
       end if
    end do
    allocate(rbuf(cmoff(p4%elem%vert%nrank)), sbuf(cmoff(p4%elem%vert%nrank)))
    ! set non-blocking receive
    itmp = 0
    do il=1, p4%elem%vert%nrank ! mpi rank loop
       if (p4%elem%vert%lrank(il) /= pe_rank) then
          itmp = itmp + 1 ! count messages
          jl = cmoff(itmp+1) - cmoff(itmp)
          call MPI_IRecv(rbuf(cmoff(itmp):cmoff(itmp+1)-1), jl, MPI_INTEGER, &
               & p4%elem%vert%lrank(il), jl, NEKO_COMM, request(itmp), ierr)
       end if
    end do

    ! redistribute vertex information
    itmp = 0
    do il = 1, p4%elem%vert%nrank ! mpi rank loop
       if (p4%elem%vert%lrank(il) /= pe_rank) then
          itmp = itmp + 1 ! count messages
          ll = cmoff(itmp)
          do jl = p4%elem%vert%loff(il), p4%elem%vert%loff(il+1) - 1
             sbuf(ll) = p4%elem%vert%lgidx(p4%elem%vert%lshare(jl)) ! vertex global id
             ll = ll + 1
             ! extract local neighbours
             neighn = msh%point_neigh(p4%elem%vert%lshare(jl))%size()
             neighl => msh%point_neigh(p4%elem%vert%lshare(jl))%array()
             sbuf(ll) = neighn ! local number of neighbours
             ll = ll + 1
             do kl = 1, neighn
                sbuf(ll) = neighl(kl)
                ll = ll + 1
             end do
          end do
          ! sanity check
          if (ll /= cmoff(itmp + 1)) &
               & call neko_error('Inconsistent number of local vertex neighbours; receive.')
          jl = cmoff(itmp+1) - cmoff(itmp)
          call MPI_Send(sbuf(cmoff(itmp):cmoff(itmp+1)-1), jl, MPI_INTEGER, &
               & p4%elem%vert%lrank(il), jl, NEKO_COMM, ierr)
       end if
    end do

    ! finalize communication
    call MPI_Waitall(p4%elem%vert%nrank - 1, request, status, ierr)

    ! extract data
    itmp = 0
    do il = 1, p4%elem%vert%nrank ! mpi rank loop
       if (p4%elem%vert%lrank(il) /= pe_rank) then
          itmp = itmp + 1 ! count messages
          ll = cmoff(itmp)
          do jl = p4%elem%vert%loff(il), p4%elem%vert%loff(il+1) - 1
             ! check vertex global id
             if (p4%elem%vert%lgidx(p4%elem%vert%lshare(jl)) /= rbuf(ll)) &
                  & call neko_error('Inconsistent global vertex number in neighbour exchange.')
             ll = ll + 1
             neighn = rbuf(ll)
             ll = ll + 1
             do kl = 1, neighn
                ml = -rbuf(ll)
                call msh%point_neigh(p4%elem%vert%lshare(jl))%push(ml)
                ll = ll + 1
             end do
          end do
          ! sanity check
          if (ll /= cmoff(itmp + 1)) &
               & call neko_error('Inconsistent number of local vertex neighbours; extract.')
       end if
    end do

    deallocate(cmoff, request, status, rbuf, sbuf)
    return
  end subroutine p4_vertex_neighbour_fill

  ! Fill the mesh type with face neighbour information
  subroutine p4_face_neighbour_fill(msh, p4)
    ! argument list
    type(mesh_t), intent(inout) :: msh
    type(p4_mesh_cnstr_t), intent(in) :: p4
    !local variables
    integer(i4) :: il, jl
    integer(i4) :: itmp, ierr
    type(stack_i4t2_t), allocatable :: tmp_neigh(:) ! to get face neighbours
    type(tuple_i4_t) :: ttmp
    integer(i4) :: neighn
    type(tuple_i4_t), pointer, dimension(:) :: neighl
    integer(i4), allocatable, dimension(:,:) :: rbuf, sbuf ! send/receive buffers
    type(MPI_Request), allocatable, dimension(:) :: request ! MPI request
    type(MPI_Status), allocatable, dimension(:) :: status ! MPI status

    ! Get list of face neighbours
    ! local face neighbours
    if (msh%gdim == 3) then
       ! gether data
       allocate(tmp_neigh(msh%mfcs))
       do il = 1, msh%mfcs
          call tmp_neigh(il)%init()
       end do
       do il = 1, msh%nelv
          do jl = 1, NEKO_HEX_NFCS
             ! I work with local element number
             ttmp%x = (/ jl, il /)
             call tmp_neigh(p4%elem%face%lmap(jl, il))%push(ttmp)
          end do
       end do
       ! fill in type array
       allocate(msh%facet_neigh(NEKO_HEX_NFCS, msh%nelv)) ! face neigbour
       do il = 1, msh%mfcs
          neighn = tmp_neigh(il)%size()
          neighl => tmp_neigh(il)%array()
          select case(neighn)
          case(1)
             ! no local neighbour
             msh%facet_neigh(neighl(1)%x(1), neighl(1)%x(2)) = 0
             ! check if face mapping is correct
             if (p4%elem%face%lmap(neighl(1)%x(1), neighl(1)%x(2)) /= il) &
                  & call neko_error('Inconsistent local face number and neighbour list.')
          case(2)
             ! swap global element numbers
             msh%facet_neigh(neighl(1)%x(1), neighl(1)%x(2)) = &
                  & int(p4%elem%gidx(neighl(2)%x(2)), i4)
             msh%facet_neigh(neighl(2)%x(1), neighl(2)%x(2)) = &
                  & int(p4%elem%gidx(neighl(1)%x(2)), i4)
             if ((p4%elem%face%lmap(neighl(1)%x(1), neighl(1)%x(2)) /= il).or. &
                  & (p4%elem%face%lmap(neighl(2)%x(1), neighl(2)%x(2)) /= il)) &
                  & call neko_error('Inconsistent local face number and neighbour list.')
          case default
             ! none or too many neighbours
             call neko_error('Invalid number of local face neighbours.')
          end select
       end do

       ! non-local face neighbours
       ! start communication
       allocate(rbuf(2, p4%elem%face%nshare), sbuf(2, p4%elem%face%nshare), &
            & request(p4%elem%face%nrank), status(p4%elem%face%nrank))
       ! set non-blocking receive
       itmp = 0
       do il = 1, p4%elem%face%nrank ! mpi rank loop
          if (p4%elem%face%lrank(il) /= pe_rank) then
             itmp = itmp + 1 ! count messages
             jl = 2*(p4%elem%face%loff(il+1) - p4%elem%face%loff(il))
             call MPI_IRecv(rbuf(:, p4%elem%face%loff(il):p4%elem%face%loff(il+1)-1), &
                  & jl, MPI_INTEGER, p4%elem%face%lrank(il), jl, NEKO_COMM, request(itmp), ierr)
          end if
       end do

       ! redistribute edge information
       do il = 1, p4%elem%face%nrank ! mpi rank loop
          if (p4%elem%face%lrank(il) /= pe_rank) then
             do jl = p4%elem%face%loff(il), p4%elem%face%loff(il+1) -1
                ! has the face local local neighbour
                neighn = tmp_neigh(p4%elem%face%lshare(jl))%size()
                neighl => tmp_neigh(p4%elem%face%lshare(jl))%array()
                select case(neighn)
                case(1)
                   ! no local neighbour
                   sbuf(1,jl) = p4%elem%face%lgidx(p4%elem%face%lshare(jl)) ! face global number
                   sbuf(2,jl) = int(p4%elem%gidx(neighl(1)%x(2)), i4) !element global number
                case(2)
                   call neko_error('Face with a local neighbour cannot have nonlocal one.')
                case default
                   ! none or too many neighbours
                   call neko_error('Invalid number of local face neighbours.')
                end select
             end do
             jl = 2*(p4%elem%face%loff(il+1) - p4%elem%face%loff(il))
             call MPI_Send(sbuf(:, p4%elem%face%loff(il):p4%elem%face%loff(il+1)-1), &
                  & jl, MPI_INTEGER, p4%elem%face%lrank(il), jl, NEKO_COMM, ierr)
          end if
       end do

       ! finalize communication
       call MPI_Waitall(p4%elem%face%nrank - 1, request, status, ierr)

       ! extract data
       do il = 1, p4%elem%face%nrank ! mpi rank loop
          if (p4%elem%face%lrank(il) /= pe_rank) then
             do jl = p4%elem%face%loff(il), p4%elem%face%loff(il+1) -1
                ! does the face global id match
                if (p4%elem%face%lgidx(p4%elem%face%lshare(jl)) /= rbuf(1,jl)) &
                     & call neko_error('Global face number does not match.')
                ! has the face local local neighbour
                neighn = tmp_neigh(p4%elem%face%lshare(jl))%size()
                neighl => tmp_neigh(p4%elem%face%lshare(jl))%array()
                select case(neighn)
                case(1)
                   ! no local neighbour
                   ! none of the faces can have more than one nonlocal neighbour
                   if (msh%facet_neigh(neighl(1)%x(1),neighl(1)%x(2)) /= 0) &
                        & call neko_error('Nonlocal face neighbour already set.')
                   msh%facet_neigh(neighl(1)%x(1),neighl(1)%x(2)) = - rbuf(2,jl)
                case(2)
                   call neko_error('Face with a local neighbour cannot have nonlocal one.')
                case default
                   ! none or too many neighbours
                   call neko_error('Invalid number of local face neighbours.')
                end select
             end do
          end if
       end do

       ! free memory
       do il = 1, msh%mfcs
          call tmp_neigh(il)%free()
       end do
       deallocate(tmp_neigh)
       deallocate(rbuf, sbuf, request, status)
    else if (msh%gdim == 2) then
       call neko_error("2D not supported for p4est yet")
    else
       call neko_error("Invalid dimension")
    end if

    return
  end subroutine p4_face_neighbour_fill

  ! Fill neighbour rank using vertex sharing information
  subroutine p4_rank_neighbour_fill(msh, p4)
    ! argument list
    type(mesh_t), intent(inout) :: msh
    type(p4_mesh_cnstr_t), intent(in) :: p4
    !local variables
    integer(i4) :: il, jl, src, dst
    type(stack_i4_t) :: neigh_order


    ! initialize neighbour ranks
    allocate(msh%neigh(0:pe_size-1))
    msh%neigh = .false.

    ! p4%elem%vert%lrank
    if (pe_size > 1) then
       ! get rank list based on vertex information
       do il= 1, p4%elem%vert%nrank ! mpi rank loop
          if (p4%elem%vert%lrank(il) /= pe_rank) & ! not sure I have to exclude myself
               & msh%neigh(p4%elem%vert%lrank(il)) = .true.
       end do

       !
       ! Generate neighbour exchange order
       !
       call neigh_order%init(pe_size)

       do il = 1, pe_size - 1
          src = modulo(pe_rank - il + pe_size, pe_size)
          dst = modulo(pe_rank + il, pe_size)
          if (msh%neigh(src) .or. msh%neigh(dst)) then
             jl = il ! adhere to standards...
             call neigh_order%push(jl)
          end if
       end do

       allocate(msh%neigh_order(neigh_order%size()))
       select type(order => neigh_order%data)
       type is (integer)
          do il = 1, neigh_order%size()
             msh%neigh_order(il) = order(il)
          end do
       end select
       call neigh_order%free()
    else
       allocate(msh%neigh_order(1))
       msh%neigh_order = 1
    end if

    return
  end subroutine p4_rank_neighbour_fill

  ! Fill the mesh type with distdata information from p4_mesh_cnstr type
  subroutine p4_distdata_fill(msh, p4)
    ! argument list
    type(mesh_t), intent(inout) :: msh
    type(p4_mesh_cnstr_t), intent(in) :: p4
    !local variables
    integer(i4) :: il, jl, kl
    integer(i4) :: itmp

    ! Set connectivity information
    call distdata_init(msh%ddata)
    ! Add shared points; local id
    do il = 1, p4%elem%vert%nrank
       if (p4%elem%vert%lrank(il) == pe_rank) then
          do jl = p4%elem%vert%loff(il), p4%elem%vert%loff(il+1) - 1
             call distdata_set_shared_point(msh%ddata, p4%elem%vert%lshare(jl))
          end do
          exit
       end if
    end do
    ! THERE WOULD BE GOOD TO ADD LOCAL TO GOBAL VERTEX NUMBERRING AS WELL

    ! Add shared faces; local id
    do il = 1, p4%elem%face%nrank
       if (p4%elem%face%lrank(il) == pe_rank) then
          itmp = il
          do jl = p4%elem%face%loff(il), p4%elem%face%loff(il+1) - 1
             call distdata_set_shared_facet(msh%ddata, p4%elem%face%lshare(jl))
          end do
          exit
       end if
    end do
    ! Get faces that are shared; this is very simple and silly way of doing, but for now it is enough
    do il = p4%elem%face%loff(il), p4%elem%face%loff(il+1) - 1
       elemf : do jl = 1, msh%nelv
          do kl = 1, NEKO_HEX_NFCS
             if (p4%elem%face%lmap(kl, jl) == p4%elem%face%lshare(il)) then
                call distdata_set_shared_el_facet(msh%ddata, jl, kl)
                ! check consistency of facet_neigh
                if (msh%facet_neigh(kl, jl) >= 0) &
                     & call neko_error('Inconsistent facet_neigh and face shared list.')
                exit elemf
             end if
          end do
       end do elemf
    end do
    ! local to global id mapping
    allocate(msh%ddata%local_to_global_facet(msh%mfcs))
    do il = 1, msh%mfcs
       ! global face id; notice type casting
       itmp = int(p4%elem%face%lgidx(il),i4)
       call distdata_set_local_to_global_facet(msh%ddata, il, itmp)
    end do

    ! Add shared edges; local id
    if (msh%gdim == 3) then
       do il = 1, p4%elem%edge%nrank
          if (p4%elem%edge%lrank(il) == pe_rank) then
             do jl = p4%elem%edge%loff(il), p4%elem%edge%loff(il+1) - 1
                call distdata_set_shared_edge(msh%ddata, p4%elem%edge%lshare(jl))
             end do
             exit
          end if
       end do
       ! local to global id mapping
       allocate(msh%ddata%local_to_global_edge(msh%meds))
       do il = 1, msh%meds
          ! global edge id; notice type casting
          itmp = int(p4%elem%edge%lgidx(il),i4)
          call distdata_set_local_to_global_edge(msh%ddata, il, itmp)
       end do
    end if

    msh%ldist = .true.

    return
  end subroutine p4_distdata_fill

  ! Fill the mesh type with boundarry condition information from p4_mesh_cnstr type
  subroutine p4_bc_fill(msh, p4)
    ! argument list
    type(mesh_t), intent(inout) :: msh
    type(p4_mesh_cnstr_t), intent(in) :: p4
    !local variables
    integer(i4) :: il, jl, ie, if
    integer(i4) :: itmp
    integer(i4) :: p_f, p_e ! dummy variables for setting periodic bc
    integer(i4), dimension(4) :: pt_id ! dummy variable for setting periodic bc

    ! Initialize element boundary condition
    allocate(msh%facet_type(2 * msh%gdim, msh%nelv))
    msh%facet_type = 0
    call msh%wall%init(msh%nelv)
    call msh%inlet%init(msh%nelv)
    call msh%outlet%init(msh%nelv)
    call msh%outlet_normal%init(msh%nelv)
    call msh%sympln%init(msh%nelv)
    call msh%periodic%init(msh%nelv)

    allocate(msh%labeled_zones(NEKO_MSH_MAX_ZLBLS))
    do il = 1, NEKO_MSH_MAX_ZLBLS
       call msh%labeled_zones(il)%init(msh%nelv)
    end do

    ! I consider periodic and labeled bc only
    !! @todo periodic bc for p4est
    ! There is inconsistency here, as I treat bc differently than neko.
    ! PLACE FOR DISCUSSION; I HAVE TO CHECK WHAT GMSH DOES WITH BC
    ! ALL THIS HAVE TO BE CHECKED IN THE FUTURE
    do il = 1, msh%nelv
       ie = il ! once again problem with inout attributre
       do jl = 1, NEKO_HEX_NFCS
          if = jl ! once again problem with inout attributre
          select case(p4%elem%bc(jl,il))
          case(-1) ! this is my periodic bc mark, which is inconsistent with neko (5)
             ! I don't think I have to do here much more than just mark a face,
             ! as periodicity is already taken into account in communicator
             ! IT WOULD BE SUFFICIENT TO HAVE THIS ZONE TO BE LIKE ANY OTHER ONE
             call mesh_mark_periodic_facet(msh, if, ie, p_f, p_e, pt_id)
             ! There is no need to call mesh apply, as I do not stick to
             ! node numbering anyhow.
          case (1:NEKO_MSH_MAX_ZLBLS) ! everything elese marked as labeled bc
             itmp = p4%elem%bc(jl,il) ! once again problem with inout attribute
             call mesh_mark_labeled_facet(msh, if, ie, itmp)
          end select
       end do
    end do

    return
  end subroutine p4_bc_fill

  ! Fill the mesh type with global connectivity information from p4_mesh_cnstr type
  subroutine p4_connect_fill(msh, p4)
    ! argument list
    type(mesh_t), intent(inout) :: msh
    type(p4_mesh_cnstr_t), intent(in) :: p4
    !local variables

    ! allocate connectivity type
    if (allocated(msh%msh_conn)) then
       call msh%msh_conn%free()
       deallocate(msh%msh_conn)
    end if
    allocate(mesh_conn_t :: msh%msh_conn)

    ! Set mesh nonconformity flag
    call msh%msh_conn%set_nonconf(.true.)

    ! copy data
    msh%msh_conn%nel = p4%elem%nelt

    allocate(msh%msh_conn%falg, source = p4%elem%falg)
    allocate(msh%msh_conn%ealg, source = p4%elem%ealg)

    msh%msh_conn%vert%lnum = p4%elem%vert%lnum
    msh%msh_conn%vert%lown = p4%elem%vert%lown
    msh%msh_conn%vert%goff = p4%elem%vert%goff
    msh%msh_conn%vert%gnum = p4%elem%vert%gnum
    allocate(msh%msh_conn%vert%lmap, source = p4%elem%vert%lmap)
    allocate(msh%msh_conn%vert%lgidx, source = p4%elem%vert%lgidx)
    allocate(msh%msh_conn%vert%lrank, source = p4%elem%vert%lrank)
    allocate(msh%msh_conn%vert%lshare, source = p4%elem%vert%lshare)
    allocate(msh%msh_conn%vert%loff, source = p4%elem%vert%loff)

    msh%msh_conn%face%lnum = p4%elem%face%lnum
    msh%msh_conn%face%lown = p4%elem%face%lown
    msh%msh_conn%face%goff = p4%elem%face%goff
    msh%msh_conn%face%gnum = p4%elem%face%gnum
    allocate(msh%msh_conn%face%lmap, source = p4%elem%face%lmap)
    allocate(msh%msh_conn%face%lgidx, source = p4%elem%face%lgidx)
    allocate(msh%msh_conn%face%lrank, source = p4%elem%face%lrank)
    allocate(msh%msh_conn%face%lshare, source = p4%elem%face%lshare)
    allocate(msh%msh_conn%face%loff, source = p4%elem%face%loff)

    if (msh%gdim == 3) then
       msh%msh_conn%edge%lnum = p4%elem%edge%lnum
       msh%msh_conn%edge%lown = p4%elem%edge%lown
       msh%msh_conn%edge%goff = p4%elem%edge%goff
       msh%msh_conn%edge%gnum = p4%elem%edge%gnum
       allocate(msh%msh_conn%edge%lmap, source = p4%elem%edge%lmap)
       allocate(msh%msh_conn%edge%lgidx, source = p4%elem%edge%lgidx)
       allocate(msh%msh_conn%edge%lrank, source = p4%elem%edge%lrank)
       allocate(msh%msh_conn%edge%lshare, source = p4%elem%edge%lshare)
       allocate(msh%msh_conn%edge%loff, source = p4%elem%edge%loff)
    end if

    return
  end subroutine p4_connect_fill
#else

  subroutine p4_manager_init(mesh_file, log_threshold)
    character(len=*), intent(in) :: mesh_file
    integer, intent(in), optional :: log_threshold

    call neko_error('NEKO needs to be built with P4EST support')
  end subroutine p4_manager_init

  subroutine p4_manager_free(log_priority)
    integer, intent(in), optional :: log_priority

    call neko_error('NEKO needs to be built with P4EST support')
  end subroutine p4_manager_free

  subroutine p4_msh_get(this, msh)
    class(p4_mesh_cnstr_t), intent(inout) :: p4
    type(mesh_t), intent(inout) :: msh

    call neko_error('NEKO needs to be built with P4EST support')
    return
  end subroutine p4_msh_get

  subroutine p4_refine(this, ref_mark, el_gidx, msh_trs, level_max,&
       &ifmod, msh_rcn)
    class(p4_mesh_cnstr_t), intent(inout) :: this
    integer(i4), dimension(:), intent(in) :: ref_mark, el_gidx
    type(mesh_manager_transfer_t), intent(in) :: msh_trs
    integer(i4), intent(in) :: level_max
    logical, intent(out) :: ifmod
    type(mesh_reconstruct_transfer_t), intent(out) :: msh_rcn

    call neko_error('NEKO needs to be built with P4EST support')
    return
  end subroutine p4_refine

#endif

  ! Following stuff should be in math, but right now there is no clear
  ! division of routines as
  ! flipv, swap and reord are not in math. For now I leave it here
  !> Use Heap Sort (p 231 Num. Rec., 1st Ed.)
  subroutine sorti4(a, ind, n)
    integer, intent(in) :: n
    integer(i4), intent(inout) :: a(n)
    integer, intent(inout) :: ind(n)
    integer(i4) :: aa
    integer :: j, ir, i, ii, l
    do j = 1, n
       ind(j) = j
    end do

    if (n.le.1) return

    l = n/2+1
    ir = n
    do while (.true.)
       if (l.gt.1) then
          l=l-1
          aa  = a  (l)
          ii  = ind(l)
       else
          aa =   a(ir)
          ii = ind(ir)
          a(ir) =   a( 1)
          ind(ir) = ind( 1)
          ir=ir-1
          if (ir.eq.1) then
             a(1) = aa
             ind(1) = ii
             return
          endif
       endif
       i=l
       j=l+l
       do while (j .le. ir)
          if (j.lt.ir) then
             if ( a(j).lt.a(j+1) ) j=j+1
          endif
          if (aa.lt.a(j)) then
             a(i) = a(j)
             ind(i) = ind(j)
             i=j
             j=j+j
          else
             j=ir+1
          endif
       end do
       a(i) = aa
       ind(i) = ii
    end do
  end subroutine sorti4

  !> sort the array acording to ind vector
  subroutine swapi4(b, ind, n)
    integer, intent(in) :: n
    integer(i4), intent(inout) :: b(n)
    integer, intent(inout) :: ind(n)
    integer(i4) :: temp(n)
    integer :: i, jj

    do i = 1, n
       temp(i)=b(i)
    end do
    do i = 1, n
       jj=ind(i)
       b(i)=temp(jj)
    end do

  end subroutine swapi4

  !> reorder the array - inverse of swap
  subroutine reordi4(b, ind, n)
    integer, intent(in) :: n
    integer(i4), intent(inout) :: b(n)
    integer, intent(inout) :: ind(n)
    integer(i4) :: temp(n)
    integer :: i, jj

    do i = 1, n
       temp(i)=b(i)
    end do
    do i = 1, n
       jj=ind(i)
       b(jj)=temp(i)
    end do

  end subroutine reordi4

end module p4est
