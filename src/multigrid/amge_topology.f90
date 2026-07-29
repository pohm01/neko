! Copyright (c) 2026, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the conditions of the
! BSD 3-Clause License are met (see Neko's COPYING file).
!
!> Prototype: hierarchical macroelement topology extraction for
!! agglomeration-based (AMGe-style) coarse spaces.
!!
!! Given a partition of the (hexahedral) element set into macroelements,
!! this module extracts the macroentity complex needed to build
!! trace-compatible coarse spaces:
!!
!!   exposed facets (with macroelement-pair labels)
!!     -> skeleton edges  (carriers of macroedges)
!!     -> macrovertices   (coarse degrees of freedom)
!!     -> macroedges      (chains of skeleton edges between macrovertices)
!!     -> macrofaces      (same-label facet patches split by the skeleton,
!!                         each recording its bounding macroedges)
!!
!! Skeleton edge rules: an edge is on the skeleton if its incident
!! exposed facets satisfy any of
!!   (a) >= 2 distinct macroelement-pair labels (triple junctions,
!!       interfaces meeting the domain boundary);
!!   (b) incidence count /= 2 (rim / non-manifold);
!!   (c) exactly two same-label facets sharing a common element
!!       (dihedral crease of a macroelement interface).
!! Macrovertex rules: skeleton degree /= 2 or >= 2 distinct incident
!! edge signatures (signature = sorted label set of the edge's exposed
!! facets); closed skeleton components get a loop breakpoint.
!!
!! DESIGN NOTES FOR NEKO INTEGRATION (prototype status):
!!  * The core takes plain arrays (element -> 8 global vertex ids), so it
!!    is framework-agnostic; a thin wrapper for mesh_t is sketched at the
!!    bottom of this file. In Neko the vertex ids would come from
!!    msh%elements(e)%e%pts(k)%p%id().
!!  * The internal open-addressing tuple hash (tuple_map_t below) is a
!!    self-contained stand-in for Neko's htable_i4t2_t / htable_i4t4_t,
!!    which mesh.f90 already uses for exactly this kind of edge/facet
!!    numbering; swap it out when integrating.
!!  * The local face/edge numbering of the reference hex matches the
!!    Octave prototype, NOT Neko's hex_t facet convention; remap when
!!    integrating.
!!  * Serial / rank-local: agglomerates are assumed not to cross rank
!!    boundaries. A facet with a single local element gets label
!!    (part,0), which conflates the physical boundary with rank
!!    interfaces; production code should consult msh%facet_neigh (and
!!    the distributed data in mesh_t) to distinguish them and to
!!    exchange skeleton flags across ranks.
!!  * Natural home: alongside the tree_amg coarse-grid infrastructure,
!!    complementing its aggregation with interface-aware macroentity
!!    structure for trace-compatible interpolation.
!!
!! MULTILEVEL USE (the coarsen_level_3d pattern): the extraction does
!! NOT run on the original hex_verts with a composed partition -- that
!! is a different ("flat") algorithm posed on the wrong space. Instead,
!! each level is described by the generic entity tables macro_mesh_t
!! (faces with vertex/edge lists, edges as vertex pairs, facet->element
!! incidence); the hex path is only the finest-level FRONT-END filling
!! those tables. After extraction, build_next_level() emits the next
!! level's tables: macroedges collapse to coarse edges (their endpoint
!! macrovertex pairs), macrofaces collapse to coarse faces (macrovertex
!! polygons bounded by macroedges), and facet->element incidence comes
!! from the stored labels. The same init_tables() then runs unchanged
!! on the coarse tables. Note the crease rule (c) lifts correctly:
!! coarse face->element incidence stores the incident level-l ELEMENTS
!! (macroelements), so two coplanar coarse interface faces between the
!! same pair of NEW macroelements have disjoint element pairs (no
!! crease, they merge), while two faces of one old macroelement meeting
!! at a dihedral share it (crease fires).
module amge_topology
  use num_types, only : i4, i8
  use utils, only : neko_error
  implicit none
  private
  public :: macro_mesh_init_hex

  !> A macroedge: an open or closed chain of vertices (local indices);
  !! chain(1) == chain(size(chain)) marks a closed loop (breakpoint at
  !! chain(1)).
  type :: macro_edge_t
     integer(i4), allocatable :: chain(:)     !< vertex chain (local ids)
     integer(i4), allocatable :: edge_ids(:)  !< fine edge indices along it
  end type macro_edge_t

  !> A macroface: a connected same-label patch of exposed facets.
  type :: macro_face_t
     integer(i4), allocatable :: facet_ids(:) !< fine facet indices
     integer(i4), allocatable :: verts(:)     !< unique vertex set (local ids)
     integer(i4), allocatable :: bnd_medge(:) !< bounding macroedge indices
     integer(i4) :: label(2)                  !< sorted macro pair, 0=exterior
  end type macro_face_t

  !> Generic per-level entity tables: what init_tables consumes and what
  !! build_next_level emits. The hex front-end fills these for level 1;
  !! coarse levels are polyhedral (variable face arity), hence CSR.
  type, public :: macro_mesh_t
     integer(i4) :: n_verts = 0
     integer(i4) :: n_elem = 0
     integer(i4) :: n_face = 0
     integer(i4) :: n_edge = 0
     integer(i4), allocatable :: face_vtx_ptr(:)  !< CSR: face -> vertices
     integer(i4), allocatable :: face_vtx_idx(:)
     integer(i4), allocatable :: face_edge_ptr(:) !< CSR: face -> edges
     integer(i4), allocatable :: face_edge_idx(:)
     integer(i4), allocatable :: edge_vtx(:,:)    !< (2, n_edge)
     integer(i4), allocatable :: face_el(:,:)     !< (2, n_face) incident elems
     integer(i4), allocatable :: face_nel(:)      !< 1 or 2
     integer(i4), allocatable :: vert_id(:)       !< local -> global vertex id
     !> Is this vertex shared with another MPI rank? Inherited from the
     !! finest-level (Q1) dofmap%shared_dof and propagated by
     !! build_next_level, which is exact because a macrovertex IS a fine
     !! vertex. Needed by the gather-scatter mapping to split dofs into
     !! rank-local and communicated groups.
     logical, allocatable :: shared_vtx(:)
   contains
     procedure, pass(this) :: free => macro_mesh_free
  end type macro_mesh_t

  !> Extracted macroentity complex.
  type, public :: macro_topology_t
     integer(i4) :: n_verts = 0               !< number of local vertices
     integer(i4) :: n_edges = 0               !< number of fine edges
     integer(i4) :: n_facets = 0              !< number of fine facets
     integer(i4) :: n_mv = 0                  !< number of macrovertices
     integer(i4) :: n_medge = 0               !< number of macroedges
     integer(i4) :: n_mface = 0               !< number of macrofaces
     logical, allocatable :: is_mv(:)         !< macrovertex flag per vertex
     integer(i4), allocatable :: vert_id(:)   !< local -> global vertex id
     logical, allocatable :: shared_vtx(:)    !< inherited rank-shared flag
     type(macro_edge_t), allocatable :: medge(:)
     type(macro_face_t), allocatable :: mface(:)
     integer(i4) :: n_macro = 0               !< number of macroelements
   contains
     procedure, pass(this) :: init_tables => macro_topology_init_tables
     procedure, pass(this) :: build_next_level => macro_topology_build_next_level
     procedure, pass(this) :: free => macro_topology_free
  end type macro_topology_t

  !> Minimal open-addressing hash for sorted integer tuples (length <= 4).
  !! Stand-in for htable_i4t2_t / htable_i4t4_t; keys stored in full and
  !! compared exactly, so hashing is only an acceleration.
  type :: tuple_map_t
     integer(i4), allocatable :: keys(:,:)    !< (4, cap); 0 = empty slot
     integer(i4), allocatable :: val(:)
     integer(i4) :: cap = 0
     integer(i4) :: n = 0
  end type tuple_map_t

  ! Reference-hex local facets and edges, in NEKO hex_t vertex ordering
  ! (src/mesh/hex.f90). Vertices are numbered lexicographically within
  ! the cell: x fastest, then y, then z, so in symmetric coordinates
  !   v1=(-,-,-) v2=(+,-,-) v3=(-,+,-) v4=(+,+,-)   [z-]
  !   v5=(-,-,+) v6=(+,-,+) v7=(-,+,+) v8=(+,+,+)   [z+]
  ! The 6 facets are ordered -x,+x,-y,+y,-z,+z. IMPORTANT: each facet's
  ! four vertices are listed in PERIMETER (cycle) order, not lexicographic
  ! order, because the face-edge builder below forms boundary edges from
  ! consecutive entries; a lexicographic listing would form diagonals.
  ! (Vertex-set membership and the sorted tuple keys are orientation-
  ! independent, so only this perimeter property matters.)
  integer(i4), parameter :: loc_face(4,6) = reshape( &
       [1,3,7,5,  2,4,8,6,  1,2,6,5,  3,4,8,7,  1,2,4,3,  5,6,8,7], [4,6])
  integer(i4), parameter :: loc_edge(2,12) = reshape( &
       [1,2, 1,3, 1,5, 2,4, 2,6, 3,4, 3,7, 4,8, 5,6, 5,7, 6,8, 7,8], [2,12])

  integer(i4), parameter :: MAX_SIG_LEN = 16  !< max labels per edge signature

contains

  !> Hex front-end: fill the generic tables (faces as quads, 12 edges per
  !! hex) via sorted-tuple hashing. Prototype uses the internal tuple map;
  !! in Neko this is htable_i4t2_t / htable_i4t4_t territory (mesh.f90).
  subroutine macro_mesh_init_hex(mmsh, nelv, hex_verts)
    type(macro_mesh_t), intent(inout) :: mmsh
    integer(i4), intent(in) :: nelv
    integer(i4), intent(in) :: hex_verts(8, nelv)
    type(tuple_map_t) :: vmap, emap, fmap
    integer(i4), allocatable :: lverts(:,:)
    integer(i4) :: e, t, q, v, g, f, nv, ne, nf, cyc(4), pr(2), key(4)

    call mmsh%free()
    call tmap_init(vmap, 16 * nelv)
    call tmap_init(emap, 32 * nelv)
    call tmap_init(fmap, 16 * nelv)
    allocate(lverts(8, nelv))
    allocate(mmsh%edge_vtx(2, 12 * nelv))
    allocate(mmsh%face_vtx_ptr(6 * nelv + 1), mmsh%face_vtx_idx(4 * 6 * nelv))
    allocate(mmsh%face_edge_ptr(6 * nelv + 1), mmsh%face_edge_idx(4 * 6 * nelv))
    allocate(mmsh%face_el(2, 6 * nelv), mmsh%face_nel(6 * nelv))
    allocate(mmsh%vert_id(8 * nelv))
    mmsh%face_nel = 0
    nv = 0; ne = 0; nf = 0
    mmsh%face_vtx_ptr(1) = 0
    mmsh%face_edge_ptr(1) = 0

    do e = 1, nelv
       do t = 1, 8
          key = [hex_verts(t, e), 0, 0, 0]
          call tmap_find_or_add(vmap, key, nv + 1, v)
          if (v .eq. nv + 1) then
             nv = v
             mmsh%vert_id(v) = hex_verts(t, e)
          end if
          lverts(t, e) = v
       end do
       do t = 1, 12
          pr = lverts(loc_edge(:, t), e)
          call sort2(pr)
          key = [pr(1), pr(2), 0, 0]
          call tmap_find_or_add(emap, key, ne + 1, g)
          if (g .eq. ne + 1) then
             ne = g
             mmsh%edge_vtx(:, g) = pr
          end if
       end do
       do t = 1, 6
          cyc = lverts(loc_face(:, t), e)
          key = cyc
          call sort4(key)
          call tmap_find_or_add(fmap, key, nf + 1, f)
          if (f .eq. nf + 1) then
             nf = f
             mmsh%face_vtx_ptr(f + 1) = mmsh%face_vtx_ptr(f) + 4
             mmsh%face_edge_ptr(f + 1) = mmsh%face_edge_ptr(f) + 4
             do q = 1, 4
                mmsh%face_vtx_idx(mmsh%face_vtx_ptr(f) + q) = cyc(q)
                pr = [cyc(q), cyc(mod(q, 4) + 1)]
                call sort2(pr)
                key = [pr(1), pr(2), 0, 0]
                call tmap_find_or_add(emap, key, -1, g)
                if (g .lt. 0) call neko_error('macro_topology: face edge missing')
                mmsh%face_edge_idx(mmsh%face_edge_ptr(f) + q) = g
             end do
          end if
          mmsh%face_nel(f) = mmsh%face_nel(f) + 1
          if (mmsh%face_nel(f) .gt. 2) &
               call neko_error('macro_topology: facet with > 2 elements')
          mmsh%face_el(mmsh%face_nel(f), f) = e
       end do
    end do
    mmsh%n_verts = nv
    mmsh%n_edge = ne
    mmsh%n_face = nf
    mmsh%n_elem = nelv
    call tmap_free(vmap)
    call tmap_free(emap)
    call tmap_free(fmap)
  end subroutine macro_mesh_init_hex

  !> Level-generic extraction: runs on any macro_mesh_t, fine (from the
  !! hex front-end) or coarse (from build_next_level).
  subroutine macro_topology_init_tables(this, mmsh, part, n_macro)
    class(macro_topology_t), intent(inout) :: this
    type(macro_mesh_t), intent(in) :: mmsh
    integer(i4), intent(in) :: part(mmsh%n_elem), n_macro

    logical, allocatable :: exposed(:)
    integer(i4), allocatable :: face_lab(:,:)
    integer(i8), allocatable :: face_labkey(:)
    integer(i4), allocatable :: ef_xadj(:), ef_adj(:)
    logical, allocatable :: is_skel(:)
    integer(i4), allocatable :: edge_sig(:)
    integer(i8), allocatable :: sig_dat(:,:)
    integer(i4), allocatable :: sig_len(:)
    integer(i4) :: n_sig
    integer(i4), allocatable :: vdeg(:), vsig(:)
    integer(i4), allocatable :: uf(:)
    logical, allocatable :: root_has_mv(:)
    integer(i4), allocatable :: root_minv(:)
    integer(i4), allocatable :: vs_xadj(:), vs_adj(:)
    logical, allocatable :: eused(:)
    integer(i4), allocatable :: edge2me(:)
    integer(i4), allocatable :: f_uf(:), f_grp(:)
    integer(i4), allocatable :: vmark(:), memark(:)
    integer(i4), allocatable :: tmp_chain(:), tmp_eids(:)
    integer(i8) :: labs(MAX_SIG_LEN)
    integer(i4) :: fl(64), nfl
    integer(i4) :: e, f, g, v, a, b, q, nf, ne, nv, deg, nu, r
    integer(i4) :: n_exp, v0, cur, nxt, ge, cnt, grp, ngrp, tclen, telen

    call this%free()
    nf = mmsh%n_face
    ne = mmsh%n_edge
    nv = mmsh%n_verts
    this%n_verts = nv
    this%n_edges = ne
    this%n_facets = nf
    this%n_macro = n_macro
    allocate(this%vert_id(nv), this%shared_vtx(nv))
    this%vert_id = mmsh%vert_id(1:nv)
    if (allocated(mmsh%shared_vtx)) then
       this%shared_vtx = mmsh%shared_vtx(1:nv)
    else
       this%shared_vtx = .false.
    end if

    ! ---------------- phase 2: exposed facets + labels --------------------
    allocate(exposed(nf), face_lab(2, nf), face_labkey(nf))
    exposed = .false.
    face_labkey = -1_i8
    n_exp = 0
    do f = 1, nf
       if (mmsh%face_nel(f) .eq. 1) then
          exposed(f) = .true.
          face_lab(:, f) = [0_i4, part(mmsh%face_el(1, f))]
       else if (part(mmsh%face_el(1, f)) .ne. part(mmsh%face_el(2, f))) then
          exposed(f) = .true.
          a = part(mmsh%face_el(1, f)); b = part(mmsh%face_el(2, f))
          face_lab(:, f) = [min(a, b), max(a, b)]
       end if
       if (exposed(f)) then
          n_exp = n_exp + 1
          face_labkey(f) = int(face_lab(1, f), i8) * int(n_macro + 1, i8) &
                         + int(face_lab(2, f), i8)
       end if
    end do

    ! ---------------- phase 3: edge -> exposed facets (CRS) ---------------
    allocate(ef_xadj(ne + 1), ef_adj(size(mmsh%face_edge_idx)))
    ef_xadj = 0
    do f = 1, nf
       if (.not. exposed(f)) cycle
       do q = mmsh%face_edge_ptr(f) + 1, mmsh%face_edge_ptr(f + 1)
          g = mmsh%face_edge_idx(q)
          ef_xadj(g + 1) = ef_xadj(g + 1) + 1
       end do
    end do
    do g = 1, ne
       ef_xadj(g + 1) = ef_xadj(g + 1) + ef_xadj(g)
    end do
    block
      integer(i4), allocatable :: fill(:)
      allocate(fill(ne)); fill = 0
      do f = 1, nf
         if (.not. exposed(f)) cycle
         do q = mmsh%face_edge_ptr(f) + 1, mmsh%face_edge_ptr(f + 1)
            g = mmsh%face_edge_idx(q)
            fill(g) = fill(g) + 1
            ef_adj(ef_xadj(g) + fill(g)) = f
         end do
      end do
    end block

    ! ---------------- phase 4: skeleton edges + signatures ----------------
    allocate(is_skel(ne), edge_sig(ne))
    allocate(sig_dat(MAX_SIG_LEN, ne), sig_len(ne))
    is_skel = .false.
    edge_sig = 0
    n_sig = 0
    do g = 1, ne
       deg = ef_xadj(g + 1) - ef_xadj(g)
       if (deg .eq. 0) cycle
       if (deg .gt. MAX_SIG_LEN) call neko_error('macro_topology: raise MAX_SIG_LEN')
       do q = 1, deg
          labs(q) = face_labkey(ef_adj(ef_xadj(g) + q))
       end do
       call sort_i8(labs(1:deg))
       nu = 1
       do q = 2, deg
          if (labs(q) .ne. labs(nu)) then
             nu = nu + 1
             labs(nu) = labs(q)
          end if
       end do
       if (nu .ge. 2) then
          is_skel(g) = .true.
       else if (deg .ne. 2) then
          is_skel(g) = .true.
       else
          ! rule (c): two same-label facets sharing an element = crease
          a = ef_adj(ef_xadj(g) + 1); b = ef_adj(ef_xadj(g) + 2)
          if (facets_share_element(mmsh%face_el(:, a), mmsh%face_nel(a), &
               mmsh%face_el(:, b), mmsh%face_nel(b))) is_skel(g) = .true.
       end if
       if (is_skel(g)) edge_sig(g) = sig_lookup(sig_dat, sig_len, n_sig, labs(1:nu))
    end do

    ! ---------------- phase 5: macrovertices ------------------------------
    allocate(vdeg(nv), vsig(nv), this%is_mv(nv))
    vdeg = 0
    vsig = 0
    do g = 1, ne
       if (.not. is_skel(g)) cycle
       do q = 1, 2
          v = mmsh%edge_vtx(q, g)
          vdeg(v) = vdeg(v) + 1
          if (vsig(v) .eq. 0) then
             vsig(v) = edge_sig(g)
          else if (vsig(v) .gt. 0 .and. vsig(v) .ne. edge_sig(g)) then
             vsig(v) = -1
          end if
       end do
    end do
    this%is_mv = .false.
    do v = 1, nv
       if (vdeg(v) .eq. 0) cycle
       if (vdeg(v) .ne. 2 .or. vsig(v) .eq. -1) this%is_mv(v) = .true.
    end do

    ! loop breakpoints
    allocate(uf(nv))
    do v = 1, nv
       uf(v) = v
    end do
    do g = 1, ne
       if (is_skel(g)) call uf_union(uf, mmsh%edge_vtx(1, g), mmsh%edge_vtx(2, g))
    end do
    allocate(root_has_mv(nv), root_minv(nv))
    root_has_mv = .false.
    root_minv = 0
    do v = 1, nv
       if (vdeg(v) .eq. 0) cycle
       r = uf_find(uf, v)
       if (this%is_mv(v)) root_has_mv(r) = .true.
       if (root_minv(r) .eq. 0) root_minv(r) = v
       if (this%vert_id(v) .lt. this%vert_id(root_minv(r))) root_minv(r) = v
    end do
    do v = 1, nv
       if (vdeg(v) .eq. 0) cycle
       r = uf_find(uf, v)
       if (.not. root_has_mv(r)) then
          this%is_mv(root_minv(r)) = .true.
          root_has_mv(r) = .true.
       end if
    end do
    this%n_mv = count(this%is_mv)

    ! ---------------- phase 7: vertex -> skeleton edges (CRS) -------------
    allocate(vs_xadj(nv + 1))
    vs_xadj = 0
    do g = 1, ne
       if (.not. is_skel(g)) cycle
       vs_xadj(mmsh%edge_vtx(1, g) + 1) = vs_xadj(mmsh%edge_vtx(1, g) + 1) + 1
       vs_xadj(mmsh%edge_vtx(2, g) + 1) = vs_xadj(mmsh%edge_vtx(2, g) + 1) + 1
    end do
    do v = 1, nv
       vs_xadj(v + 1) = vs_xadj(v + 1) + vs_xadj(v)
    end do
    allocate(vs_adj(vs_xadj(nv + 1)))
    block
      integer(i4), allocatable :: fill(:)
      allocate(fill(nv)); fill = 0
      do g = 1, ne
         if (.not. is_skel(g)) cycle
         do q = 1, 2
            v = mmsh%edge_vtx(q, g)
            fill(v) = fill(v) + 1
            vs_adj(vs_xadj(v) + fill(v)) = g
         end do
      end do
    end block

    ! ---------------- phase 8: macroedges (chain tracing) -----------------
    allocate(eused(ne), edge2me(ne))
    eused = .false.
    edge2me = 0
    allocate(this%medge(max(count(is_skel), 1)))
    allocate(tmp_chain(count(is_skel) + 1), tmp_eids(max(count(is_skel), 1)))
    this%n_medge = 0
    do v0 = 1, nv
       if (.not. this%is_mv(v0)) cycle
       do q = vs_xadj(v0) + 1, vs_xadj(v0 + 1)
          g = vs_adj(q)
          if (eused(g)) cycle
          cur = v0
          ge = g
          tclen = 1; tmp_chain(1) = v0
          telen = 0
          do
             eused(ge) = .true.
             telen = telen + 1; tmp_eids(telen) = ge
             if (mmsh%edge_vtx(1, ge) .eq. cur) then
                nxt = mmsh%edge_vtx(2, ge)
             else
                nxt = mmsh%edge_vtx(1, ge)
             end if
             tclen = tclen + 1; tmp_chain(tclen) = nxt
             cur = nxt
             if (this%is_mv(cur)) exit
             nxt = 0
             do a = vs_xadj(cur) + 1, vs_xadj(cur + 1)
                if (.not. eused(vs_adj(a))) then
                   nxt = vs_adj(a)
                   exit
                end if
             end do
             if (nxt .eq. 0) exit
             ge = nxt
          end do
          this%n_medge = this%n_medge + 1
          associate (me => this%medge(this%n_medge))
            allocate(me%chain(tclen), me%edge_ids(telen))
            me%chain = tmp_chain(1:tclen)
            me%edge_ids = tmp_eids(1:telen)
          end associate
          edge2me(tmp_eids(1:telen)) = this%n_medge
       end do
    end do
    if (any(is_skel .and. .not. eused)) &
         call neko_error('macro_topology: untraced skeleton edges remain')

    ! ---------------- phase 9: macrofaces ----------------------------------
    allocate(f_uf(nf))
    do f = 1, nf
       f_uf(f) = f
    end do
    do g = 1, ne
       if (is_skel(g)) cycle
       nfl = 0
       do q = ef_xadj(g) + 1, ef_xadj(g + 1)
          nfl = nfl + 1
          if (nfl .gt. 64) call neko_error('macro_topology: raise fl buffer')
          fl(nfl) = ef_adj(q)
       end do
       do a = 1, nfl
          do b = a + 1, nfl
             if (face_labkey(fl(a)) .eq. face_labkey(fl(b))) &
                  call uf_union(f_uf, fl(a), fl(b))
          end do
       end do
    end do
    allocate(f_grp(nf))
    f_grp = 0
    ngrp = 0
    do f = 1, nf
       if (.not. exposed(f)) cycle
       r = uf_find(f_uf, f)
       if (f_grp(r) .eq. 0) then
          ngrp = ngrp + 1
          f_grp(r) = ngrp
       end if
       f_grp(f) = f_grp(r)
    end do
    this%n_mface = ngrp
    allocate(this%mface(max(ngrp, 1)))
    allocate(vmark(nv), memark(max(this%n_medge, 1)))
    vmark = 0
    memark = 0
    do grp = 1, ngrp
       cnt = 0
       do f = 1, nf
          if (exposed(f) .and. f_grp(f) .eq. grp) cnt = cnt + 1
       end do
       associate (mf => this%mface(grp))
         allocate(mf%facet_ids(cnt))
         cnt = 0
         do f = 1, nf
            if (.not. (exposed(f) .and. f_grp(f) .eq. grp)) cycle
            cnt = cnt + 1
            mf%facet_ids(cnt) = f
            mf%label = face_lab(:, f)
         end do
         cnt = 0
         do a = 1, size(mf%facet_ids)
            f = mf%facet_ids(a)
            do q = mmsh%face_vtx_ptr(f) + 1, mmsh%face_vtx_ptr(f + 1)
               v = mmsh%face_vtx_idx(q)
               if (vmark(v) .ne. grp) then
                  vmark(v) = grp
                  cnt = cnt + 1
               end if
            end do
         end do
         allocate(mf%verts(cnt))
         cnt = 0
         do a = 1, size(mf%facet_ids)
            f = mf%facet_ids(a)
            do q = mmsh%face_vtx_ptr(f) + 1, mmsh%face_vtx_ptr(f + 1)
               v = mmsh%face_vtx_idx(q)
               if (vmark(v) .eq. grp) then
                  vmark(v) = -grp
                  cnt = cnt + 1
                  mf%verts(cnt) = v
               end if
            end do
         end do
         cnt = 0
         do a = 1, size(mf%facet_ids)
            f = mf%facet_ids(a)
            do q = mmsh%face_edge_ptr(f) + 1, mmsh%face_edge_ptr(f + 1)
               g = mmsh%face_edge_idx(q)
               if (is_skel(g) .and. edge2me(g) .gt. 0) then
                  if (memark(edge2me(g)) .ne. grp) then
                     memark(edge2me(g)) = grp
                     cnt = cnt + 1
                  end if
               end if
            end do
         end do
         allocate(mf%bnd_medge(cnt))
         cnt = 0
         do a = 1, size(mf%facet_ids)
            f = mf%facet_ids(a)
            do q = mmsh%face_edge_ptr(f) + 1, mmsh%face_edge_ptr(f + 1)
               g = mmsh%face_edge_idx(q)
               if (is_skel(g) .and. edge2me(g) .gt. 0) then
                  if (memark(edge2me(g)) .eq. grp) then
                     memark(edge2me(g)) = -grp
                     cnt = cnt + 1
                     mf%bnd_medge(cnt) = edge2me(g)
                  end if
               end if
            end do
         end do
         if (size(mf%bnd_medge) .eq. 0) then
            cnt = 0
            do a = 1, size(mf%verts)
               if (this%is_mv(mf%verts(a))) cnt = cnt + 1
            end do
            if (cnt .eq. 0) call neko_error( &
                 'macro_topology: closed macroface without boundary (not handled)')
         end if
       end associate
    end do
  end subroutine macro_topology_init_tables

  !> Emit the next level's generic tables from this extraction: macroedges
  !! collapse to coarse edges (endpoint macrovertex pairs), macrofaces to
  !! coarse faces (macrovertex polygons bounded by macroedges), and coarse
  !! facet->element incidence comes from the stored labels. The recursion
  !! closes: init_tables runs unchanged on the result.
  subroutine macro_topology_build_next_level(this, mmshC)
    class(macro_topology_t), intent(in) :: this
    type(macro_mesh_t), intent(inout) :: mmshC
    integer(i4), allocatable :: cidx(:)
    integer(i4) :: v, k, a, cnt

    call mmshC%free()
    allocate(cidx(this%n_verts))
    cidx = 0
    cnt = 0
    do v = 1, this%n_verts
       if (this%is_mv(v)) then
          cnt = cnt + 1
          cidx(v) = cnt
       end if
    end do
    mmshC%n_verts = cnt
    mmshC%n_elem = this%n_macro
    mmshC%n_edge = this%n_medge
    mmshC%n_face = this%n_mface
    allocate(mmshC%vert_id(cnt), mmshC%shared_vtx(cnt))
    mmshC%shared_vtx = .false.
    do v = 1, this%n_verts
       if (cidx(v) .gt. 0) then
          mmshC%vert_id(cidx(v)) = this%vert_id(v)
          ! a macrovertex IS a fine vertex, so its shared status is
          ! inherited verbatim -- no communication needed to determine it
          if (allocated(this%shared_vtx)) &
               mmshC%shared_vtx(cidx(v)) = this%shared_vtx(v)
       end if
    end do
    allocate(mmshC%edge_vtx(2, max(this%n_medge, 1)))
    do k = 1, this%n_medge
       associate (ch => this%medge(k)%chain)
         mmshC%edge_vtx(1, k) = min(cidx(ch(1)), cidx(ch(size(ch))))
         mmshC%edge_vtx(2, k) = max(cidx(ch(1)), cidx(ch(size(ch))))
       end associate
    end do
    allocate(mmshC%face_el(2, max(this%n_mface, 1)), mmshC%face_nel(max(this%n_mface, 1)))
    allocate(mmshC%face_vtx_ptr(this%n_mface + 1), mmshC%face_edge_ptr(this%n_mface + 1))
    mmshC%face_vtx_ptr(1) = 0
    mmshC%face_edge_ptr(1) = 0
    do k = 1, this%n_mface
       cnt = 0
       do a = 1, size(this%mface(k)%verts)
          if (this%is_mv(this%mface(k)%verts(a))) cnt = cnt + 1
       end do
       mmshC%face_vtx_ptr(k + 1) = mmshC%face_vtx_ptr(k) + cnt
       mmshC%face_edge_ptr(k + 1) = mmshC%face_edge_ptr(k) + size(this%mface(k)%bnd_medge)
       mmshC%face_nel(k) = 0
       do a = 1, 2
          if (this%mface(k)%label(a) .gt. 0) then
             mmshC%face_nel(k) = mmshC%face_nel(k) + 1
             mmshC%face_el(mmshC%face_nel(k), k) = this%mface(k)%label(a)
          end if
       end do
    end do
    allocate(mmshC%face_vtx_idx(mmshC%face_vtx_ptr(this%n_mface + 1)))
    allocate(mmshC%face_edge_idx(mmshC%face_edge_ptr(this%n_mface + 1)))
    do k = 1, this%n_mface
       cnt = 0
       do a = 1, size(this%mface(k)%verts)
          v = this%mface(k)%verts(a)
          if (this%is_mv(v)) then
             cnt = cnt + 1
             mmshC%face_vtx_idx(mmshC%face_vtx_ptr(k) + cnt) = cidx(v)
          end if
       end do
       do a = 1, size(this%mface(k)%bnd_medge)
          mmshC%face_edge_idx(mmshC%face_edge_ptr(k) + a) = this%mface(k)%bnd_medge(a)
       end do
    end do
  end subroutine macro_topology_build_next_level

  !> Deallocate a level table set.
  subroutine macro_mesh_free(this)
    class(macro_mesh_t), intent(inout) :: this
    if (allocated(this%face_vtx_ptr)) deallocate(this%face_vtx_ptr)
    if (allocated(this%face_vtx_idx)) deallocate(this%face_vtx_idx)
    if (allocated(this%face_edge_ptr)) deallocate(this%face_edge_ptr)
    if (allocated(this%face_edge_idx)) deallocate(this%face_edge_idx)
    if (allocated(this%edge_vtx)) deallocate(this%edge_vtx)
    if (allocated(this%face_el)) deallocate(this%face_el)
    if (allocated(this%face_nel)) deallocate(this%face_nel)
    if (allocated(this%vert_id)) deallocate(this%vert_id)
    if (allocated(this%shared_vtx)) deallocate(this%shared_vtx)
    this%n_verts = 0; this%n_elem = 0; this%n_face = 0; this%n_edge = 0
  end subroutine macro_mesh_free


  !> Deallocate everything.
  subroutine macro_topology_free(this)
    class(macro_topology_t), intent(inout) :: this
    integer(i4) :: k
    if (allocated(this%is_mv)) deallocate(this%is_mv)
    if (allocated(this%vert_id)) deallocate(this%vert_id)
    if (allocated(this%shared_vtx)) deallocate(this%shared_vtx)
    if (allocated(this%medge)) then
       do k = 1, this%n_medge
          if (allocated(this%medge(k)%chain)) deallocate(this%medge(k)%chain)
          if (allocated(this%medge(k)%edge_ids)) deallocate(this%medge(k)%edge_ids)
       end do
       deallocate(this%medge)
    end if
    if (allocated(this%mface)) then
       do k = 1, this%n_mface
          if (allocated(this%mface(k)%facet_ids)) deallocate(this%mface(k)%facet_ids)
          if (allocated(this%mface(k)%verts)) deallocate(this%mface(k)%verts)
          if (allocated(this%mface(k)%bnd_medge)) deallocate(this%mface(k)%bnd_medge)
       end do
       deallocate(this%mface)
    end if
    this%n_mv = 0; this%n_medge = 0; this%n_mface = 0
    this%n_verts = 0; this%n_edges = 0; this%n_facets = 0
  end subroutine macro_topology_free

  ! ======================= internal helpers ==============================

  !> True if the two facets' element lists intersect.
  pure function facets_share_element(ela, na, elb, nb) result(res)
    integer(i4), intent(in) :: ela(2), elb(2), na, nb
    logical :: res
    integer(i4) :: a, b
    res = .false.
    do a = 1, na
       do b = 1, nb
          if (ela(a) .eq. elb(b)) then
             res = .true.
             return
          end if
       end do
    end do
  end function facets_share_element

  !> Signature dictionary lookup/insert (linear scan; fine at these sizes).
  function sig_lookup(sig_dat, sig_len, n_sig, labs) result(id)
    integer(i8), intent(inout) :: sig_dat(:,:)
    integer(i4), intent(inout) :: sig_len(:), n_sig
    integer(i8), intent(in) :: labs(:)
    integer(i4) :: id, s
    do s = 1, n_sig
       if (sig_len(s) .eq. size(labs)) then
          if (all(sig_dat(1:sig_len(s), s) .eq. labs)) then
             id = s
             return
          end if
       end if
    end do
    n_sig = n_sig + 1
    sig_len(n_sig) = size(labs)
    sig_dat(1:size(labs), n_sig) = labs
    id = n_sig
  end function sig_lookup

  ! ---- union-find ----
  recursive function uf_find(uf, v) result(r)
    integer(i4), intent(inout) :: uf(:)
    integer(i4), intent(in) :: v
    integer(i4) :: r
    if (uf(v) .eq. v) then
       r = v
    else
       r = uf_find(uf, uf(v))
       uf(v) = r
    end if
  end function uf_find

  subroutine uf_union(uf, a, b)
    integer(i4), intent(inout) :: uf(:)
    integer(i4), intent(in) :: a, b
    integer(i4) :: ra, rb
    ra = uf_find(uf, a)
    rb = uf_find(uf, b)
    if (ra .ne. rb) uf(min(ra, rb)) = max(ra, rb)
  end subroutine uf_union

  ! ---- tiny sorted-tuple hash map (swap for neko htable when integrating) ----
  subroutine tmap_init(map, cap_hint)
    type(tuple_map_t), intent(inout) :: map
    integer(i4), intent(in) :: cap_hint
    map%cap = 2 * cap_hint + 17
    allocate(map%keys(4, map%cap), map%val(map%cap))
    map%keys = 0
    map%n = 0
  end subroutine tmap_init

  subroutine tmap_free(map)
    type(tuple_map_t), intent(inout) :: map
    if (allocated(map%keys)) deallocate(map%keys)
    if (allocated(map%val)) deallocate(map%val)
    map%cap = 0; map%n = 0
  end subroutine tmap_free

  !> Find key; if absent and new_val >= 0, insert with new_val.
  !! Returns the stored value, or -1 if absent and new_val < 0.
  subroutine tmap_find_or_add(map, key, new_val, val)
    type(tuple_map_t), intent(inout) :: map
    integer(i4), intent(in) :: key(4), new_val
    integer(i4), intent(out) :: val
    integer(i8) :: h
    integer(i4) :: s, probe
    h = 1469598103934665603_i8
    do s = 1, 4
       h = ieor(h, int(key(s), i8)) * 1099511628211_i8
    end do
    probe = int(modulo(h, int(map%cap, i8)), i4) + 1
    do
       if (map%keys(1, probe) .eq. 0) then
          if (new_val .lt. 0) then
             val = -1
             return
          end if
          if (map%n .ge. map%cap / 2) &
               call neko_error('macro_topology: tuple map over half full')
          map%keys(:, probe) = key
          map%val(probe) = new_val
          map%n = map%n + 1
          val = new_val
          return
       end if
       if (all(map%keys(:, probe) .eq. key)) then
          val = map%val(probe)
          return
       end if
       probe = probe + 1
       if (probe .gt. map%cap) probe = 1
    end do
  end subroutine tmap_find_or_add

  ! ---- small sorts ----
  pure subroutine sort2(a)
    integer(i4), intent(inout) :: a(2)
    integer(i4) :: t
    if (a(1) .gt. a(2)) then
       t = a(1); a(1) = a(2); a(2) = t
    end if
  end subroutine sort2

  pure subroutine sort4(a)
    integer(i4), intent(inout) :: a(4)
    integer(i4) :: i, j, t
    do i = 2, 4
       t = a(i)
       j = i - 1
       do while (j .ge. 1)
          if (a(j) .le. t) exit
          a(j + 1) = a(j)
          j = j - 1
       end do
       a(j + 1) = t
    end do
  end subroutine sort4

  pure subroutine sort_i8(a)
    integer(i8), intent(inout) :: a(:)
    integer(i8) :: t
    integer(i4) :: i, j
    do i = 2, size(a)
       t = a(i)
       j = i - 1
       do while (j .ge. 1)
          if (a(j) .le. t) exit
          a(j + 1) = a(j)
          j = j - 1
       end do
       a(j + 1) = t
    end do
  end subroutine sort_i8

end module amge_topology
