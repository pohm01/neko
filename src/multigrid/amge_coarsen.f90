! Copyright (c) 2026, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the conditions of the
! BSD 3-Clause License are met (see Neko's COPYING file).
!
!> Prototype: one coarsening step of the trace-compatible local Galerkin
!! (AMGe-style) hierarchy -- the Fortran counterpart of the Octave
!! coarsen_level_3d. Builds on the macroentity extraction of
!! macro_topology and populates the shared cell-wise data structures of
!! macroelm_data (melm_amg_level_t with a flat CSR elm_vtx_ptr/idx dof
!! layout, matrix_t local/coarse blocks, and the transfer stored on the
!! coarse level's %tr, with PiTilde: coarse -> finer level). All
!! solve-phase work then happens in the duplicated cell-wise storage via
!! the melm_* operations of macroelm_data.
!!
!! Dense pseudoinverses use LAPACK dsyev (Neko already links LAPACK).
!! Serial / rank-local, same MPI caveats as macro_topology.
module amge_coarsen
  use num_types, only : i4, rp
  use utils, only : neko_error
  use matrix, only : matrix_t
  use comm, only : NEKO_COMM, pe_size
  use amge_topology, only : macro_topology_t, &
       macro_topology_build_next_level_owned, macro_mesh_elem_face_csr, &
       macro_mesh_splice_ghost
  use amge_level, only : amge_level_t
  use amge_ghost, only : amge_ghost_t, amge_ghost_exchange, amge_ghost_part_ext
  implicit none
  private

  real(rp), parameter :: PINV_RTOL = 1e-12_rp

  !> Toggle for the O(1)-per-macroelement debug checks in schur_extend and
  !! the O(1)-per-macroedge/macroface checks in compute_edge_maps/
  !! compute_face_maps (an eigendecomposition plus a few dense matmuls on
  !! small local blocks each). Cheap relative to the rest of setup, but not
  !! free; default on since this module is prototype/debug-stage code, but
  !! easy to flip off for a build that wants to skip the extra work.
  logical, public :: AMGE_DEBUG_CHECKS = .true.

  !> Running worst-case diagnostics accumulated across a level's
  !! macroelements by schur_extend, verifying eq.(48)-(53) of the AMGe
  !! theory note (see [[amge_theory_reference]]): symmetry of S_dM (52),
  !! interior stationarity of the harmonic extension (48), the "smart"
  !! Schur-complement A_c against the naive Pi~^T A Pi~ (51), and SPSD-ness
  !! of A_c (required by eq. 3/6).
  type, public :: schur_check_t
     real(rp) :: sym_err = 0.0_rp
     real(rp) :: harmonic_err = 0.0_rp
     real(rp) :: galerkin_err = 0.0_rp
     real(rp) :: min_eig = huge(1.0_rp)
     integer(i4) :: n_checked = 0
  end type schur_check_t

  !> Small dense-block scratch type used only inside this module for the
  !! ragged per-entity trace maps (Q_E, Q_F).
  type :: dblk_t
     real(rp), allocatable :: x(:,:)
  end type dblk_t

  type :: ivec_t
     integer(i4), allocatable :: x(:)
  end type ivec_t

  public :: coarsen_level_3d, agglomerate_level

contains

  !> Compact-growth agglomeration on the element dual graph (adjacency =
  !! shared facet): absorb the frontier element with the most
  !! connections into the current cluster.
  !! (Integration note: tree_amg's greedy aggregation plays this role.)
  subroutine agglomerate_level(lvl, target_size, part, n_macro)
    type(amge_level_t), intent(in) :: lvl
    integer(i4), intent(in) :: target_size
    integer(i4), allocatable, intent(inout) :: part(:)
    integer(i4), intent(out) :: n_macro
    integer(i4), allocatable :: adj_ptr(:), adj_idx(:), fill(:), frontier(:)
    logical, allocatable :: infr(:)
    integer(i4) :: ne, f, e, a, b, s, cnt, nfr, best, bscore, sc, q, c
    integer(i4) :: i, j, tmp
    integer(i4), allocatable :: rand_order(:)
    real(kind=rp) :: r

    ne = lvl%mmsh%n_elem
    if (allocated(part)) deallocate(part)
    allocate(part(ne))
    allocate(adj_ptr(ne + 1))
    adj_ptr = 0
    do f = 1, lvl%mmsh%n_face
       if (lvl%mmsh%face_nel(f) .eq. 2) then
          adj_ptr(lvl%mmsh%face_el(1, f) + 1) = adj_ptr(lvl%mmsh%face_el(1, f) + 1) + 1
          adj_ptr(lvl%mmsh%face_el(2, f) + 1) = adj_ptr(lvl%mmsh%face_el(2, f) + 1) + 1
       end if
    end do
    do e = 1, ne
       adj_ptr(e + 1) = adj_ptr(e + 1) + adj_ptr(e)
    end do
    allocate(adj_idx(adj_ptr(ne + 1)), fill(ne))
    fill = 0
    do f = 1, lvl%mmsh%n_face
       if (lvl%mmsh%face_nel(f) .eq. 2) then
          a = lvl%mmsh%face_el(1, f); b = lvl%mmsh%face_el(2, f)
          fill(a) = fill(a) + 1; adj_idx(adj_ptr(a) + fill(a)) = b
          fill(b) = fill(b) + 1; adj_idx(adj_ptr(b) + fill(b)) = a
       end if
    end do

    ! Initialize a random permutation
    allocate( rand_order( ne ) )
    do i = 1, ne
       rand_order(i) = i
    end do
    !! Shuffle rand_order using Fisher-Yates algorithm
    !do i = ne, 2, -1
    !   call random_number(r)
    !   j = int(r * real(i, kind=rp)) + 1
    !   tmp = rand_order(i)
    !   rand_order(i) = rand_order(j)
    !   rand_order(j) = tmp
    !end do

    part = 0
    n_macro = 0
    allocate(frontier(ne), infr(ne))
    do i = 1, ne
       s = rand_order(i)
       if (part(s) .ne. 0) cycle
       n_macro = n_macro + 1
       part(s) = n_macro
       cnt = 1
       nfr = 0
       infr = .false.
       do q = adj_ptr(s) + 1, adj_ptr(s + 1)
          e = adj_idx(q)
          if (part(e) .eq. 0 .and. .not. infr(e)) then
             nfr = nfr + 1; frontier(nfr) = e; infr(e) = .true.
          end if
       end do
       do while (cnt .lt. target_size .and. nfr .gt. 0)
          best = 0; bscore = -1
          do a = 1, nfr
             e = frontier(a)
             if (part(e) .ne. 0) cycle
             sc = 0
             do q = adj_ptr(e) + 1, adj_ptr(e + 1)
                if (part(adj_idx(q)) .eq. n_macro) sc = sc + 1
             end do
             if (sc .gt. bscore) then
                bscore = sc; best = a
             end if
          end do
          if (best .eq. 0) exit
          c = frontier(best)
          frontier(best) = frontier(nfr); nfr = nfr - 1
          part(c) = n_macro
          cnt = cnt + 1
          do q = adj_ptr(c) + 1, adj_ptr(c + 1)
             e = adj_idx(q)
             if (part(e) .eq. 0 .and. .not. infr(e)) then
                nfr = nfr + 1; frontier(nfr) = e; infr(e) = .true.
             end if
          end do
       end do
    end do
  end subroutine agglomerate_level

  !> One coarsening step: lvl + part -> coarse level lvlC + transfer tr
  !! (extracted topology returned too). Steps: extract macroentities;
  !! one shared Q_E per macroedge and Q_F per macroface; per
  !! macroelement assemble, fill Q_dM hierarchically (vertices -> edges
  !! -> faces), extend harmonically, form A_c = Q' S Q; emit the coarse
  !! numerical level (tables via build_next_level).
  subroutine coarsen_level_3d(lvl, part, n_macro, topo, lvlC, use_ghost)
    !> intent(inout), not (in): amge_gs_correct_shared_count below drives
    !! lvl%gsh's comm object (transient send/recv buffers), not just reads
    !! from lvl.
    type(amge_level_t), intent(inout) :: lvl
    integer(i4), intent(in) :: part(:), n_macro
    type(macro_topology_t), intent(inout) :: topo
    type(amge_level_t), intent(inout) :: lvlC
    !> Ghost-extend the rank-boundary topology/trace-map extraction for this
    !! step (every level; see amge_ghost.f90). Ignored (treated as
    !! .false.) whenever pe_size == 1, so a single-rank run is unaffected.
    logical, optional, intent(in) :: use_ghost

    type(dblk_t), allocatable :: qe(:), qf(:)
    type(ivec_t), allocatable :: fb(:), fi(:)      ! face boundary/interior verts
    type(ivec_t), allocatable :: melems(:)
    type(ivec_t) :: myfaces, myedges, bset, cg
    integer(i4), allocatable :: g2l(:), cidx(:), bpos_of(:)
    logical, allocatable :: markv(:), marke(:)
    real(rp), allocatable :: am(:,:), qdm(:,:)
    real(rp), allocatable :: w(:)
    integer(i4) :: m, e, k, a, q, v, nd, nb, nc, nmv
    logical :: ghosted
    type(schur_check_t) :: schk
    type(amge_ghost_t) :: gh
    type(amge_level_t) :: lvl_ext
    integer(i4), allocatable :: part_ext(:)
    integer(i4), allocatable :: p_vtx_ptr(:), p_vtx_idx(:)
    integer(i4), allocatable :: p_face_ptr(:), p_fv_ptr(:), p_fv_idx(:)
    integer(i4), allocatable :: p_fe_ptr(:)
    integer(i4), allocatable :: p_fe_vtx(:,:)
    integer(i4), allocatable :: p_amat_ptr(:)
    real(rp), allocatable :: p_amat(:)
    integer(i4), allocatable :: g2l_ext(:)
    integer(i4) :: n_ext, ge, nv_e
    integer(i4) :: moff, n_macro_topo

    ghosted = .false.
    if (present(use_ghost)) ghosted = use_ghost .and. (pe_size > 1)
    moff = 0

    if (ghosted) then
       ! one vertex-layer ghost exchange so this rank's topology/trace-map
       ! extraction sees a complete neighborhood at every rank-boundary
       ! entity -- see amge_ghost.f90 for the rationale.
       call build_ghost_payload_local(lvl, p_vtx_ptr, p_vtx_idx, &
            p_face_ptr, p_fv_ptr, p_fv_idx, p_fe_ptr, p_fe_vtx, &
            p_amat_ptr, p_amat)
       call amge_ghost_exchange(NEKO_COMM%mpi_val, p_vtx_ptr, p_vtx_idx, &
            p_face_ptr, p_fv_ptr, p_fv_idx, p_fe_ptr, p_fe_vtx, &
            p_amat_ptr, p_amat, part, n_macro, gh)
       call amge_ghost_part_ext(gh, part, part_ext)
       moff = gh%macro_offset
       n_macro_topo = gh%n_macro_global

       ! splice the ghost elements' own face/edge/vertex sub-tables (just
       ! received) onto a copy of this level's own mesh -- NOT a rebuild
       ! from a raw vertex list, which only works at level 0 (see
       ! amge_ghost.f90's header and macro_mesh_splice_ghost's doc comment)
       call macro_mesh_splice_ghost(lvl%mmsh, gh%n_ghost, gh%vtx_ptr, &
            gh%vtx_idx, gh%face_ptr, gh%face_vtx_ptr, gh%face_vtx_idx, &
            gh%face_edge_ptr, gh%face_edge_vtx, lvl_ext%mmsh)

       ! lvl_ext's elm_vtx_ptr/idx + AM: local elements keep their own
       ! (unchanged local numbering, since macro_mesh_splice_ghost leaves
       ! vertices 1..lvl%mmsh%n_verts untouched); ghost elements are
       ! appended after them, their vertex ids translated through a fresh
       ! global -> extended-local map built off lvl_ext%mmsh%vert_id.
       n_ext = lvl%mmsh%n_elem + gh%n_ghost
       allocate(lvl_ext%elm_vtx_ptr(n_ext + 1), lvl_ext%AM(n_ext))
       lvl_ext%elm_vtx_ptr(1:lvl%mmsh%n_elem + 1) = lvl%elm_vtx_ptr
       allocate(lvl_ext%elm_vtx_idx(lvl%elm_vtx_ptr(lvl%mmsh%n_elem + 1) + &
            gh%vtx_ptr(gh%n_ghost + 1)))
       lvl_ext%elm_vtx_idx(1:lvl%elm_vtx_ptr(lvl%mmsh%n_elem + 1)) = lvl%elm_vtx_idx
       do e = 1, lvl%mmsh%n_elem
          call lvl_ext%AM(e)%init(lvl%ndof_el(e), lvl%ndof_el(e))
          lvl_ext%AM(e)%x = lvl%AM(e)%x
       end do
       allocate(g2l_ext(maxval(lvl_ext%mmsh%vert_id)))
       g2l_ext = 0
       do v = 1, lvl_ext%mmsh%n_verts
          g2l_ext(lvl_ext%mmsh%vert_id(v)) = v
       end do
       do ge = 1, gh%n_ghost
          e = lvl%mmsh%n_elem + ge
          nv_e = gh%vtx_ptr(ge + 1) - gh%vtx_ptr(ge)
          lvl_ext%elm_vtx_ptr(e + 1) = lvl_ext%elm_vtx_ptr(e) + nv_e
          do k = 1, nv_e
             lvl_ext%elm_vtx_idx(lvl_ext%elm_vtx_ptr(e) + k) = &
                  g2l_ext(gh%vtx_idx(gh%vtx_ptr(ge) + k))
          end do
          call lvl_ext%AM(e)%init(nv_e, nv_e)
          lvl_ext%AM(e)%x = reshape( &
               gh%amat(gh%amat_ptr(ge) + 1 : gh%amat_ptr(ge + 1)), [nv_e, nv_e])
       end do

       call topo%init_tables(lvl_ext%mmsh, part_ext, n_macro_topo)
       call compute_edge_maps(lvl_ext, topo, qe)
       call compute_face_maps(lvl_ext, topo, fb, fi, qf)
    else
       call topo%init_tables(lvl%mmsh, part, n_macro)
       call compute_edge_maps(lvl, topo, qe)
       call compute_face_maps(lvl, topo, fb, fi, qf)
    end if

    ! member lists per macroelement
    allocate(melems(n_macro))
    block
      integer(i4), allocatable :: cnts(:)
      allocate(cnts(n_macro))
      cnts = 0
      do e = 1, lvl%mmsh%n_elem
         cnts(part(e)) = cnts(part(e)) + 1
      end do
      do m = 1, n_macro
         allocate(melems(m)%x(cnts(m)))
      end do
      cnts = 0
      do e = 1, lvl%mmsh%n_elem
         m = part(e)
         cnts(m) = cnts(m) + 1
         melems(m)%x(cnts(m)) = e
      end do
    end block

    ! coarse dof numbering
    allocate(cidx(lvl%mmsh%n_verts))
    cidx = 0
    nmv = 0
    do v = 1, lvl%mmsh%n_verts
       if (topo%is_mv(v)) then
          nmv = nmv + 1
          cidx(v) = nmv
       end if
    end do

    ! transfer lives on the COARSE level, PiTilde: coarse -> finer level
    associate (tr => lvlC%tr)
    call tr%free()
    tr%n_fine = lvl%mmsh%n_verts
    tr%n_coarse = nmv
    tr%n_melm = n_macro
    tr%n_fine_unique = lvl%mmsh%n_verts
    tr%n_coarse_unique = nmv
    allocate(tr%maps(n_macro))
    allocate(tr%winv(lvl%mmsh%n_verts), w(lvl%mmsh%n_verts))
    w = 0.0_rp

    ! coarse level shell: CSR dof lists (elm_vtx_ptr/idx) + matrix_t blocks
    call lvlC%mmsh%free()
    if (ghosted) then
       call macro_topology_build_next_level_owned(topo, lvlC%mmsh, &
            lvl%mmsh%n_verts, moff, moff + n_macro)
    else
       call topo%build_next_level(lvlC%mmsh)
    end if
    if (lvlC%mmsh%n_verts .ne. nmv) then
       call neko_error('coarsen_level_3d: coarse mesh vertex count ' // &
            'disagrees with the transfer''s n_coarse')
    end if
    allocate(lvlC%elm_vtx_ptr(n_macro + 1), lvlC%AM(n_macro))
    lvlC%elm_vtx_ptr(1) = 0
    allocate(g2l(lvl%mmsh%n_verts), markv(lvl%mmsh%n_verts))
    allocate(marke(max(topo%n_medge, 1)))
    g2l = 0
    markv = .false.
    marke = .false.

    ! first pass counts coarse dofs per macroelement to size elm_vtx_ptr
    do m = 1, n_macro
       call macro_dof_list(lvl, melems(m)%x, markv, tr%maps(m)%fdofs)
       nc = 0
       do a = 1, size(tr%maps(m)%fdofs)
          if (topo%is_mv(tr%maps(m)%fdofs(a))) nc = nc + 1
       end do
       lvlC%elm_vtx_ptr(m + 1) = lvlC%elm_vtx_ptr(m) + nc
    end do
    allocate(lvlC%elm_vtx_idx(lvlC%elm_vtx_ptr(n_macro + 1)))

    do m = 1, n_macro
       associate (mp => tr%maps(m))
         nd = size(mp%fdofs)
         do a = 1, nd
            g2l(mp%fdofs(a)) = a
         end do
         ! child indexing (flat chloc_ptr / chloc_idx)
         mp%child_idx = melems(m)%x
         allocate(mp%chloc_ptr(size(melems(m)%x) + 1))
         mp%chloc_ptr(1) = 0
         do q = 1, size(melems(m)%x)
            mp%chloc_ptr(q + 1) = mp%chloc_ptr(q) + lvl%ndof_el(melems(m)%x(q))
         end do
         allocate(mp%chloc_idx(mp%chloc_ptr(size(melems(m)%x) + 1)))
         allocate(am(nd, nd))
         am = 0.0_rp
         do q = 1, size(melems(m)%x)
            e = melems(m)%x(q)
            do a = 1, lvl%ndof_el(e)
               mp%chloc_idx(mp%chloc_ptr(q) + a) = &
                    g2l(lvl%elm_vtx_idx(lvl%elm_vtx_ptr(e) + a))
            end do
            call add_block(lvl, e, g2l, am)
         end do

         ! incident entities, boundary set, coarse dofs. topo labels with
         ! GLOBAL macro ids when ghost-extended (moff==0 otherwise).
         call incident_entities(topo, m + moff, marke, myfaces, myedges)
         call boundary_set(topo, myfaces, myedges, mp%fdofs, markv, bset, cg)
         nb = size(bset%x)
         nc = size(cg%x)
         allocate(mp%cdofs(nc))
         do a = 1, nc
            mp%cdofs(a) = cidx(cg%x(a))
         end do

         ! Q_dM: hierarchical fill (vertices -> edges -> faces)
         allocate(bpos_of(lvl%mmsh%n_verts))
         do a = 1, nb
            bpos_of(bset%x(a)) = a
         end do
         allocate(qdm(nb, nc))
         qdm = 0.0_rp
         call fill_qdm(topo, myfaces, myedges, qe, fb, fi, qf, bpos_of, &
                       cidx, cg%x, qdm)

         ! interior harmonic extension, PiTilde, local Galerkin
         call mp%PiTilde%init(nd, nc)
         call lvlC%AM(m)%init(nc, nc)
         call schur_extend(am, mp%fdofs, bset%x, markv, qdm, &
                           mp%PiTilde%x, lvlC%AM(m)%x, schk)

         do a = 1, nc
            lvlC%elm_vtx_idx(lvlC%elm_vtx_ptr(m) + a) = mp%cdofs(a)
         end do
         do a = 1, nd
            w(mp%fdofs(a)) = w(mp%fdofs(a)) + 1.0_rp
            g2l(mp%fdofs(a)) = 0
         end do
         deallocate(am, qdm, bpos_of)
       end associate
    end do
    ! w(v) so far only counts THIS RANK's own macroelements -- correct it
    ! to the true cross-rank total at fine vertices lvl%mmsh%shared_vtx
    ! flags as shared, using lvl's own (already-initialized) gs handle.
    call lvl%gsh%correct_shared_count(lvl%elm_vtx_idx, w)
    where (w .gt. 0.0_rp)
       tr%winv = 1.0_rp / w
    elsewhere
       tr%winv = 1.0_rp
    end where
    end associate
    if (ghosted) call gh%free()
    if (AMGE_DEBUG_CHECKS .and. schk%n_checked .gt. 0) then
       write(*, '("   [check] per-macroelement (n=", I0, "): S sym = ", ES10.3, &
            & "  harmonic resid = ", ES10.3, "  |PiT^TAPi - Ac| = ", ES10.3, &
            & "  min eig(Ac) = ", ES10.3)') schk%n_checked, schk%sym_err, &
            schk%harmonic_err, schk%galerkin_err, schk%min_eig
    end if
  end subroutine coarsen_level_3d

  !> Extract, for every local element, its own vertex/face/edge CSR
  !! sub-tables from lvl%mmsh (already built, at any level) plus its
  !! dense matrix -- the payload amge_ghost_exchange ships to
  !! neighboring ranks. Unlike the old hex-only build_hv_amat_local, this
  !! makes no arity assumption: it reuses whatever facet/edge structure
  !! lvl%mmsh already has (built by macro_mesh_init_hex at level 0, or
  !! build_next_level[_owned] at level >= 1) rather than re-deriving it
  !! from a fixed template -- which is impossible above level 0 anyway,
  !! see amge_ghost.f90's header. A face shared by two local elements is
  !! reported once per element (harmless duplication: the receiver dedups
  !! purely by vertex-set identity, see macro_mesh_splice_ghost).
  subroutine build_ghost_payload_local(lvl, vtx_ptr, vtx_idx, &
       face_ptr, face_vtx_ptr, face_vtx_idx, face_edge_ptr, face_edge_vtx, &
       amat_ptr, amat)
    type(amge_level_t), intent(in) :: lvl
    integer(i4), allocatable, intent(out) :: vtx_ptr(:), vtx_idx(:)
    integer(i4), allocatable, intent(out) :: face_ptr(:), face_vtx_ptr(:)
    integer(i4), allocatable, intent(out) :: face_vtx_idx(:)
    integer(i4), allocatable, intent(out) :: face_edge_ptr(:)
    integer(i4), allocatable, intent(out) :: face_edge_vtx(:,:)
    integer(i4), allocatable, intent(out) :: amat_ptr(:)
    real(rp), allocatable, intent(out) :: amat(:)
    integer(i4), allocatable :: eface_ptr(:), eface_idx(:)
    integer(i4) :: nelv, e, k, off, a, f, g, nv_e, slot

    nelv = lvl%mmsh%n_elem
    call macro_mesh_elem_face_csr(lvl%mmsh, eface_ptr, eface_idx)

    allocate(vtx_ptr(nelv + 1), face_ptr(nelv + 1), amat_ptr(nelv + 1))
    vtx_ptr(1) = 0; face_ptr(1) = 0; amat_ptr(1) = 0
    do e = 1, nelv
       nv_e = lvl%ndof_el(e)
       vtx_ptr(e + 1) = vtx_ptr(e) + nv_e
       amat_ptr(e + 1) = amat_ptr(e) + nv_e * nv_e
       face_ptr(e + 1) = face_ptr(e) + (eface_ptr(e + 1) - eface_ptr(e))
    end do
    allocate(vtx_idx(vtx_ptr(nelv + 1)), amat(amat_ptr(nelv + 1)))

    allocate(face_vtx_ptr(face_ptr(nelv + 1) + 1))
    allocate(face_edge_ptr(face_ptr(nelv + 1) + 1))
    face_vtx_ptr(1) = 0
    face_edge_ptr(1) = 0
    slot = 0
    do e = 1, nelv
       do a = eface_ptr(e) + 1, eface_ptr(e + 1)
          f = eface_idx(a)
          slot = slot + 1
          face_vtx_ptr(slot + 1) = face_vtx_ptr(slot) + &
               (lvl%mmsh%face_vtx_ptr(f + 1) - lvl%mmsh%face_vtx_ptr(f))
          face_edge_ptr(slot + 1) = face_edge_ptr(slot) + &
               (lvl%mmsh%face_edge_ptr(f + 1) - lvl%mmsh%face_edge_ptr(f))
       end do
    end do
    allocate(face_vtx_idx(face_vtx_ptr(face_ptr(nelv + 1) + 1)))
    allocate(face_edge_vtx(2, face_edge_ptr(face_ptr(nelv + 1) + 1)))

    slot = 0
    do e = 1, nelv
       off = lvl%elm_vtx_ptr(e)
       do k = 1, lvl%ndof_el(e)
          vtx_idx(vtx_ptr(e) + k) = lvl%mmsh%vert_id(lvl%elm_vtx_idx(off + k))
       end do
       amat(amat_ptr(e) + 1 : amat_ptr(e + 1)) = &
            reshape(lvl%AM(e)%x, [lvl%ndof_el(e)**2])
       do a = eface_ptr(e) + 1, eface_ptr(e + 1)
          f = eface_idx(a)
          slot = slot + 1
          do k = 1, lvl%mmsh%face_vtx_ptr(f + 1) - lvl%mmsh%face_vtx_ptr(f)
             face_vtx_idx(face_vtx_ptr(slot) + k) = lvl%mmsh%vert_id( &
                  lvl%mmsh%face_vtx_idx(lvl%mmsh%face_vtx_ptr(f) + k))
          end do
          do k = 1, lvl%mmsh%face_edge_ptr(f + 1) - lvl%mmsh%face_edge_ptr(f)
             g = lvl%mmsh%face_edge_idx(lvl%mmsh%face_edge_ptr(f) + k)
             face_edge_vtx(1, face_edge_ptr(slot) + k) = &
                  lvl%mmsh%vert_id(lvl%mmsh%edge_vtx(1, g))
             face_edge_vtx(2, face_edge_ptr(slot) + k) = &
                  lvl%mmsh%vert_id(lvl%mmsh%edge_vtx(2, g))
          end do
       end do
    end do
  end subroutine build_ghost_payload_local

  ! ================== trace maps ==================

  !> One Q_E per macroedge: energy minimization over the neighborhood of
  !! elements sharing >= 2 chain vertices, Schur onto the chain, harmonic
  !! interior given the endpoint values. Rows in chain order; closed
  !! chains get a single column (loop breakpoint).
  subroutine compute_edge_maps(lvl, topo, qe)
    type(amge_level_t), intent(in) :: lvl
    type(macro_topology_t), intent(in) :: topo
    type(dblk_t), allocatable, intent(out) :: qe(:)
    logical, allocatable :: inchain(:)
    integer(i4), allocatable :: g2l(:)
    type(ivec_t) :: touching, dofs
    real(rp), allocatable :: a(:,:), sm(:,:), piv(:,:)
    integer(i4) :: k, nch, nb, nf, i
    logical :: closed
    real(rp) :: sym_err, min_eig

    allocate(qe(topo%n_medge))
    allocate(inchain(lvl%mmsh%n_verts), g2l(lvl%mmsh%n_verts))
    inchain = .false.
    g2l = 0
    sym_err = 0.0_rp
    min_eig = huge(1.0_rp)
    do k = 1, topo%n_medge
       associate (ch => topo%medge(k)%chain)
         nch = size(ch)
         closed = ch(1) .eq. ch(nch)
         nb = merge(nch - 1, nch, closed)
         do i = 1, nb
            inchain(ch(i)) = .true.
         end do
         call neighborhood_dofs(lvl, inchain, 2, ch(1:nb), dofs, touching, g2l)
         allocate(a(size(dofs%x), size(dofs%x)))
         a = 0.0_rp
         do i = 1, size(touching%x)
            call add_block(lvl, touching%x(i), g2l, a)
         end do
         allocate(sm(nb, nb))
         call schur_onto_leading(a, nb, sm)
         if (AMGE_DEBUG_CHECKS) then
            ! eq.(31): S_E is a Schur complement of a symmetric matrix, so
            ! it must be symmetric and PSD (the neighborhood energy (32)
            ! cannot be negative).
            sym_err = max(sym_err, maxval(abs(sm - transpose(sm))))
            min_eig = min(min_eig, min_eigenvalue(sm))
         end if
         ! harmonic along the chain: identity at breakpoints, interior
         ! rows = -pinv(S_FF) S_FC
         allocate(qe(k)%x(nb, merge(1, 2, closed)))
         qe(k)%x = 0.0_rp
         qe(k)%x(1, 1) = 1.0_rp
         if (.not. closed) qe(k)%x(nb, 2) = 1.0_rp
         nf = merge(nb - 1, nb - 2, closed)
         if (nf .gt. 0) then
            allocate(piv(nf, nf))
            if (closed) then
               call sym_pinv(sm(2:nb, 2:nb), piv)
               qe(k)%x(2:nb, 1:1) = -matmul(piv, sm(2:nb, 1:1))
            else
               call sym_pinv(sm(2:nb-1, 2:nb-1), piv)
               qe(k)%x(2:nb-1, :) = -matmul(piv, sm(2:nb-1, [1, nb]))
            end if
            deallocate(piv)
         end if
         do i = 1, nb
            inchain(ch(i)) = .false.
         end do
         do i = 1, size(dofs%x)
            g2l(dofs%x(i)) = 0
         end do
         deallocate(a, sm)
       end associate
    end do
    if (AMGE_DEBUG_CHECKS .and. topo%n_medge .gt. 0) then
       write(*, '("   [check] Q_E (n=", I0, "): S_E sym = ", ES10.3, &
            & "  min eig(S_E) = ", ES10.3)') topo%n_medge, sym_err, min_eig
    end if
  end subroutine compute_edge_maps

  !> One Q_F per macroface, with a globally fixed sorted split of the
  !! patch into boundary vertices (bounding-chain vertices +
  !! macrovertices) and interior vertices: energy minimization over the
  !! neighborhood of elements sharing >= 3 patch vertices, Schur onto the
  !! patch, then Q_F = -pinv(S_ii) S_ib.
  subroutine compute_face_maps(lvl, topo, fb, fi, qf)
    type(amge_level_t), intent(in) :: lvl
    type(macro_topology_t), intent(in) :: topo
    type(ivec_t), allocatable, intent(out) :: fb(:), fi(:)
    type(dblk_t), allocatable, intent(out) :: qf(:)
    logical, allocatable :: mark(:), inface(:)
    integer(i4), allocatable :: g2l(:), lead(:)
    type(ivec_t) :: touching, dofs
    real(rp), allocatable :: a(:,:), sm(:,:), piv(:,:)
    integer(i4) :: k, nb, ni, i
    real(rp) :: sym_err, min_eig
    integer(i4) :: n_checked

    allocate(fb(topo%n_mface), fi(topo%n_mface), qf(topo%n_mface))
    allocate(mark(lvl%mmsh%n_verts), inface(lvl%mmsh%n_verts))
    allocate(g2l(lvl%mmsh%n_verts))
    mark = .false.
    inface = .false.
    g2l = 0
    sym_err = 0.0_rp
    min_eig = huge(1.0_rp)
    n_checked = 0
    do k = 1, topo%n_mface
       call face_split(topo, k, mark, fb(k), fi(k))
       nb = size(fb(k)%x)
       ni = size(fi(k)%x)
       allocate(qf(k)%x(ni, nb))
       if (ni .eq. 0) cycle
       associate (mf => topo%mface(k))
         do i = 1, size(mf%verts)
            inface(mf%verts(i)) = .true.
         end do
         allocate(lead(nb + ni))
         lead(1:nb) = fb(k)%x
         lead(nb+1:nb+ni) = fi(k)%x
         call neighborhood_dofs(lvl, inface, 3, lead, dofs, touching, g2l)
         allocate(a(size(dofs%x), size(dofs%x)))
         a = 0.0_rp
         do i = 1, size(touching%x)
            call add_block(lvl, touching%x(i), g2l, a)
         end do
         allocate(sm(nb + ni, nb + ni), piv(ni, ni))
         call schur_onto_leading(a, nb + ni, sm)
         if (AMGE_DEBUG_CHECKS) then
            ! same rationale as compute_edge_maps: S_F is a Schur
            ! complement of a symmetric matrix, so symmetric and PSD.
            sym_err = max(sym_err, maxval(abs(sm - transpose(sm))))
            min_eig = min(min_eig, min_eigenvalue(sm))
            n_checked = n_checked + 1
         end if
         call sym_pinv(sm(nb+1:nb+ni, nb+1:nb+ni), piv)
         qf(k)%x = -matmul(piv, sm(nb+1:nb+ni, 1:nb))
         do i = 1, size(mf%verts)
            inface(mf%verts(i)) = .false.
         end do
         do i = 1, size(dofs%x)
            g2l(dofs%x(i)) = 0
         end do
         deallocate(a, sm, piv, lead)
       end associate
    end do
    if (AMGE_DEBUG_CHECKS .and. n_checked .gt. 0) then
       write(*, '("   [check] Q_F (n=", I0, "): S_F sym = ", ES10.3, &
            & "  min eig(S_F) = ", ES10.3)') n_checked, sym_err, min_eig
    end if
  end subroutine compute_face_maps

  !> Sorted boundary/interior vertex split of macroface k: boundary =
  !! bounding-chain vertices + macrovertices in the patch; interior =
  !! the remaining patch vertices.
  subroutine face_split(topo, k, mark, bv, iv)
    type(macro_topology_t), intent(in) :: topo
    integer(i4), intent(in) :: k
    logical, intent(inout) :: mark(:)
    type(ivec_t), intent(inout) :: bv, iv
    integer(i4) :: a, q, v, nb, ni
    associate (mf => topo%mface(k))
      nb = 0
      do a = 1, size(mf%bnd_medge)
         associate (ch => topo%medge(mf%bnd_medge(a))%chain)
           do q = 1, size(ch)
              if (.not. mark(ch(q))) then
                 mark(ch(q)) = .true.; nb = nb + 1
              end if
           end do
         end associate
      end do
      do q = 1, size(mf%verts)
         v = mf%verts(q)
         if (topo%is_mv(v) .and. .not. mark(v)) then
            mark(v) = .true.; nb = nb + 1
         end if
      end do
      ni = 0
      do q = 1, size(mf%verts)
         if (.not. mark(mf%verts(q))) ni = ni + 1
      end do
      if (allocated(bv%x)) deallocate(bv%x)
      if (allocated(iv%x)) deallocate(iv%x)
      allocate(bv%x(nb), iv%x(ni))
      ni = 0
      do q = 1, size(mf%verts)
         v = mf%verts(q)
         if (.not. mark(v)) then
            ni = ni + 1
            iv%x(ni) = v
         end if
      end do
      ! collect + unmark boundary (chains first, then leftover macrovertices)
      nb = 0
      do a = 1, size(mf%bnd_medge)
         associate (ch => topo%medge(mf%bnd_medge(a))%chain)
           do q = 1, size(ch)
              if (mark(ch(q))) then
                 mark(ch(q)) = .false.
                 nb = nb + 1
                 bv%x(nb) = ch(q)
              end if
           end do
         end associate
      end do
      do q = 1, size(mf%verts)
         v = mf%verts(q)
         if (mark(v)) then
            mark(v) = .false.
            nb = nb + 1
            bv%x(nb) = v
         end if
      end do
      call sort_i4(bv%x)
      call sort_i4(iv%x)
    end associate
  end subroutine face_split

  ! ================== per-macroelement pieces ==================

  !> Sorted unique dof list of the given elements (scratch mark array).
  subroutine macro_dof_list(lvl, elems, markv, dofs)
    type(amge_level_t), intent(in) :: lvl
    integer(i4), intent(in) :: elems(:)
    logical, intent(inout) :: markv(:)
    integer(i4), allocatable, intent(inout) :: dofs(:)
    integer(i4) :: q, e, a, v, nd
    nd = 0
    do q = 1, size(elems)
       e = elems(q)
       do a = 1, lvl%ndof_el(e)
          v = lvl%elm_vtx_idx(lvl%elm_vtx_ptr(e) + a)
          if (.not. markv(v)) then
             markv(v) = .true.
             nd = nd + 1
          end if
       end do
    end do
    if (allocated(dofs)) deallocate(dofs)
    allocate(dofs(nd))
    nd = 0
    do q = 1, size(elems)
       e = elems(q)
       do a = 1, lvl%ndof_el(e)
          v = lvl%elm_vtx_idx(lvl%elm_vtx_ptr(e) + a)
          if (markv(v)) then
             markv(v) = .false.
             nd = nd + 1
             dofs(nd) = v
          end if
       end do
    end do
    call sort_i4(dofs)
  end subroutine macro_dof_list

  !> Incident macrofaces (label contains m) and the union of their
  !! bounding macroedges.
  subroutine incident_entities(topo, m, marke, myfaces, myedges)
    type(macro_topology_t), intent(in) :: topo
    integer(i4), intent(in) :: m
    logical, intent(inout) :: marke(:)
    type(ivec_t), intent(inout) :: myfaces, myedges
    integer(i4) :: k, a, cnt
    cnt = 0
    do k = 1, topo%n_mface
       if (topo%mface(k)%label(1) .eq. m .or. topo%mface(k)%label(2) .eq. m) cnt = cnt + 1
    end do
    if (allocated(myfaces%x)) deallocate(myfaces%x)
    allocate(myfaces%x(cnt))
    cnt = 0
    do k = 1, topo%n_mface
       if (topo%mface(k)%label(1) .eq. m .or. topo%mface(k)%label(2) .eq. m) then
          cnt = cnt + 1
          myfaces%x(cnt) = k
       end if
    end do
    cnt = 0
    do a = 1, size(myfaces%x)
       associate (bm => topo%mface(myfaces%x(a))%bnd_medge)
         do k = 1, size(bm)
            if (.not. marke(bm(k))) then
               marke(bm(k)) = .true.
               cnt = cnt + 1
            end if
         end do
       end associate
    end do
    if (allocated(myedges%x)) deallocate(myedges%x)
    allocate(myedges%x(cnt))
    cnt = 0
    do a = 1, size(myfaces%x)
       associate (bm => topo%mface(myfaces%x(a))%bnd_medge)
         do k = 1, size(bm)
            if (marke(bm(k))) then
               marke(bm(k)) = .false.
               cnt = cnt + 1
               myedges%x(cnt) = bm(k)
            end if
         end do
       end associate
    end do
    call sort_i4(myedges%x)
  end subroutine incident_entities

  !> Boundary vertex set B (sorted; asserted to lie inside the dof list)
  !! and coarse dofs cG (macrovertices among dofs, sorted).
  subroutine boundary_set(topo, myfaces, myedges, dofs, markv, bset, cg)
    type(macro_topology_t), intent(in) :: topo
    type(ivec_t), intent(in) :: myfaces, myedges
    integer(i4), intent(in) :: dofs(:)
    logical, intent(inout) :: markv(:)
    type(ivec_t), intent(inout) :: bset, cg
    integer(i4) :: a, k, q, v, nb, nc
    nb = 0
    do a = 1, size(myfaces%x)
       k = myfaces%x(a)
       do q = 1, size(topo%mface(k)%verts)
          v = topo%mface(k)%verts(q)
          if (.not. markv(v)) then
             markv(v) = .true.; nb = nb + 1
          end if
       end do
    end do
    do a = 1, size(myedges%x)
       k = myedges%x(a)
       do q = 1, size(topo%medge(k)%chain)
          v = topo%medge(k)%chain(q)
          if (.not. markv(v)) then
             markv(v) = .true.; nb = nb + 1
          end if
       end do
    end do
    nc = 0
    do a = 1, size(dofs)
       if (topo%is_mv(dofs(a))) then
          nc = nc + 1
          if (.not. markv(dofs(a))) then
             markv(dofs(a)) = .true.; nb = nb + 1
          end if
       end if
    end do
    if (allocated(bset%x)) deallocate(bset%x)
    if (allocated(cg%x)) deallocate(cg%x)
    allocate(bset%x(nb), cg%x(nc))
    nb = 0; nc = 0
    do a = 1, size(dofs)              ! dofs sorted -> B, cG emitted sorted
       v = dofs(a)
       if (markv(v)) then
          markv(v) = .false.
          nb = nb + 1
          bset%x(nb) = v
       end if
       if (topo%is_mv(v)) then
          nc = nc + 1
          cg%x(nc) = v
       end if
    end do
    if (nb .ne. size(bset%x)) call neko_error( &
         'macro_coarsen: boundary verts not contained in macroelement dofs')
  end subroutine boundary_set

  !> Hierarchical fill of Q_dM: macrovertex identity rows, macroedge
  !! chain rows from Q_E, macroface interior rows from Q_F applied to
  !! the already-filled boundary rows.
  subroutine fill_qdm(topo, myfaces, myedges, qe, fb, fi, qf, bpos_of, &
                      cidx, cg, qdm)
    type(macro_topology_t), intent(in) :: topo
    type(ivec_t), intent(in) :: myfaces, myedges
    type(dblk_t), intent(in) :: qe(:), qf(:)
    type(ivec_t), intent(in) :: fb(:), fi(:)
    integer(i4), intent(in) :: bpos_of(:), cidx(:), cg(:)
    real(rp), intent(inout) :: qdm(:,:)
    integer(i4), allocatable :: cpos(:)
    real(rp), allocatable :: brows(:,:)
    integer(i4) :: a, k, q, nch, i, nc
    logical :: closed

    nc = size(cg)
    allocate(cpos(maxval(cidx)))
    do a = 1, nc
       cpos(cidx(cg(a))) = a
    end do
    do a = 1, nc                                   ! tier 1: macrovertices
       qdm(bpos_of(cg(a)), a) = 1.0_rp
    end do
    do a = 1, size(myedges%x)                      ! tier 2: macroedges
       k = myedges%x(a)
       associate (ch => topo%medge(k)%chain)
         nch = size(ch)
         if (nch .le. 2) cycle
         closed = ch(1) .eq. ch(nch)
         do q = 2, nch - 1
            qdm(bpos_of(ch(q)), cpos(cidx(ch(1)))) = qe(k)%x(q, 1)
            if (.not. closed) &
                 qdm(bpos_of(ch(q)), cpos(cidx(ch(nch)))) = qe(k)%x(q, 2)
         end do
       end associate
    end do
    do a = 1, size(myfaces%x)                      ! tier 3: macrofaces
       k = myfaces%x(a)
       if (size(fi(k)%x) .eq. 0) cycle
       allocate(brows(size(fb(k)%x), nc))
       do i = 1, size(fb(k)%x)
          brows(i, :) = qdm(bpos_of(fb(k)%x(i)), :)
       end do
       brows = matmul(qf(k)%x, brows)              ! now (nI x nc)
       do i = 1, size(fi(k)%x)
          qdm(bpos_of(fi(k)%x(i)), :) = brows(i, :)
       end do
       deallocate(brows)
    end do
  end subroutine fill_qdm

  !> Interior harmonic extension and local Galerkin:
  !!   Xi = pinv(A_II) A_IB,  S = A_BB - A_BI Xi,
  !!   PiTilde = [Q ; -Xi Q],  A_c = Q' S Q.
  !! @param chk  running debug-check accumulator (updated in place; see
  !!             schur_check_t and AMGE_DEBUG_CHECKS). Verifies, per
  !!             macroelement, eq.(52) (S symmetry), eq.(48) (the harmonic
  !!             extension's interior stationarity: A_{I_M,:} Pi~_M = 0),
  !!             eq.(51) (the Schur-complement shortcut ac against the
  !!             brute-force Pi~^T A Pi~), and eq.(3)/(6) (A_c SPSD).
  subroutine schur_extend(am, dofs, bset, markv, qdm, pit, ac, chk)
    real(rp), intent(in) :: am(:,:), qdm(:,:)
    integer(i4), intent(in) :: dofs(:), bset(:)
    logical, intent(inout) :: markv(:)
    real(rp), intent(out) :: pit(:,:), ac(:,:)
    type(schur_check_t), intent(inout) :: chk
    integer(i4), allocatable :: bp(:), ip(:)
    real(rp), allocatable :: piv(:,:), xi(:,:), s(:,:)
    real(rp), allocatable :: resid(:,:), ac_bf(:,:)
    integer(i4) :: a, nd, nb, ni, ib, ii

    nd = size(dofs)
    nb = size(bset)
    ni = nd - nb
    do a = 1, nb
       markv(bset(a)) = .true.
    end do
    allocate(bp(nb), ip(max(ni, 1)))
    ib = 0; ii = 0
    do a = 1, nd
       if (markv(dofs(a))) then
          ib = ib + 1; bp(ib) = a
       else
          ii = ii + 1; ip(ii) = a
       end if
    end do
    do a = 1, nb
       markv(bset(a)) = .false.
    end do
    if (ib .ne. nb) call neko_error('macro_coarsen: B/dofs split mismatch')

    if (ni .gt. 0) then
       allocate(piv(ni, ni), xi(ni, nb), s(nb, nb))
       call sym_pinv(am(ip(1:ni), ip(1:ni)), piv)
       xi = matmul(piv, am(ip(1:ni), bp))
       s = am(bp, bp) - matmul(am(bp, ip(1:ni)), xi)
       pit(bp, :) = qdm
       pit(ip(1:ni), :) = -matmul(xi, qdm)
       ac = matmul(transpose(qdm), matmul(s, qdm))

       if (AMGE_DEBUG_CHECKS) then
          ! eq.(52): S_dM must be symmetric (Schur complement of a
          ! symmetric am).
          chk%sym_err = max(chk%sym_err, maxval(abs(s - transpose(s))))
          ! eq.(48): the interior values are the energy-minimizing
          ! (harmonic) extension of the boundary trace, i.e. the discrete
          ! Euler-Lagrange condition A_{I_M,:} Pi~_M = 0 must hold at the
          ! constructed Pi~_M -- this is the optimality condition the
          ! argmin in (48) actually asserts, checked directly rather than
          ! trusted from the algebra.
          allocate(resid(ni, size(pit, 2)))
          resid = matmul(am(ip(1:ni), :), pit)
          chk%harmonic_err = max(chk%harmonic_err, maxval(abs(resid)))
          deallocate(resid)
       end if
    else
       pit(bp, :) = qdm
       ac = matmul(transpose(qdm), matmul(am(bp, bp), qdm))
    end if

    if (AMGE_DEBUG_CHECKS) then
       ! eq.(51): the Schur-complement shortcut above must equal the naive
       ! (defining) Pi~^T A Pi~ -- an independent brute-force cross-check
       ! of schur_extend's own algebra, not just a re-derivation of it.
       allocate(ac_bf(size(ac,1), size(ac,2)))
       ac_bf = matmul(transpose(pit), matmul(am, pit))
       chk%galerkin_err = max(chk%galerkin_err, maxval(abs(ac_bf - ac)))
       deallocate(ac_bf)
       ! eq.(3)/(6): A_c must be SPSD.
       chk%min_eig = min(chk%min_eig, min_eigenvalue(ac))
       chk%n_checked = chk%n_checked + 1
    end if
  end subroutine schur_extend

  ! ================== shared low-level helpers ==================

  !> Neighborhood of elements sharing >= min_shared marked vertices, and
  !! its dof list ordered [lead (as given), others sorted]; fills the
  !! global-to-local scratch map accordingly (caller resets via dofs).
  subroutine neighborhood_dofs(lvl, marked, min_shared, lead, dofs, touching, g2l)
    type(amge_level_t), intent(in) :: lvl
    logical, intent(in) :: marked(:)
    integer(i4), intent(in) :: min_shared, lead(:)
    type(ivec_t), intent(inout) :: dofs, touching
    integer(i4), intent(inout) :: g2l(:)
    integer(i4) :: e, a, v, cnt, nt, noth, nlead
    nlead = size(lead)
    nt = 0
    do e = 1, lvl%mmsh%n_elem
       cnt = 0
       do a = 1, lvl%ndof_el(e)
          if (marked(lvl%elm_vtx_idx(lvl%elm_vtx_ptr(e) + a))) cnt = cnt + 1
       end do
       if (cnt .ge. min_shared) nt = nt + 1
    end do
    if (allocated(touching%x)) deallocate(touching%x)
    allocate(touching%x(nt))
    nt = 0
    do e = 1, lvl%mmsh%n_elem
       cnt = 0
       do a = 1, lvl%ndof_el(e)
          if (marked(lvl%elm_vtx_idx(lvl%elm_vtx_ptr(e) + a))) cnt = cnt + 1
       end do
       if (cnt .ge. min_shared) then
          nt = nt + 1
          touching%x(nt) = e
       end if
    end do
    noth = 0
    do e = 1, nt
       do a = 1, lvl%ndof_el(touching%x(e))
          v = lvl%elm_vtx_idx(lvl%elm_vtx_ptr(touching%x(e)) + a)
          if (.not. marked(v) .and. g2l(v) .eq. 0) then
             g2l(v) = -1
             noth = noth + 1
          end if
       end do
    end do
    if (allocated(dofs%x)) deallocate(dofs%x)
    allocate(dofs%x(nlead + noth))
    dofs%x(1:nlead) = lead
    noth = nlead
    do e = 1, nt
       do a = 1, lvl%ndof_el(touching%x(e))
          v = lvl%elm_vtx_idx(lvl%elm_vtx_ptr(touching%x(e)) + a)
          if (g2l(v) .eq. -1) then
             g2l(v) = -2
             noth = noth + 1
             dofs%x(noth) = v
          end if
       end do
    end do
    call sort_i4(dofs%x(nlead+1:noth))
    do a = 1, noth
       g2l(dofs%x(a)) = a
    end do
  end subroutine neighborhood_dofs

  !> Accumulate element e's local matrix into a dense block via g2l.
  subroutine add_block(lvl, e, g2l, a)
    type(amge_level_t), intent(in) :: lvl
    integer(i4), intent(in) :: e, g2l(:)
    real(rp), intent(inout) :: a(:,:)
    integer(i4) :: n, i, j
    n = lvl%ndof_el(e)
    do j = 1, n
       do i = 1, n
          a(g2l(lvl%elm_vtx_idx(lvl%elm_vtx_ptr(e) + i)), &
            g2l(lvl%elm_vtx_idx(lvl%elm_vtx_ptr(e) + j))) = &
               a(g2l(lvl%elm_vtx_idx(lvl%elm_vtx_ptr(e) + i)), &
                 g2l(lvl%elm_vtx_idx(lvl%elm_vtx_ptr(e) + j))) + &
               lvl%AM(e)%x(i, j)
       end do
    end do
  end subroutine add_block

  !> Schur complement onto the leading nl dofs: S = A_ll - A_lo pinv(A_oo) A_ol.
  subroutine schur_onto_leading(a, nl, s)
    real(rp), intent(in) :: a(:,:)
    integer(i4), intent(in) :: nl
    real(rp), intent(out) :: s(:,:)
    real(rp), allocatable :: piv(:,:)
    integer(i4) :: nd, no
    nd = size(a, 1)
    no = nd - nl
    if (no .eq. 0) then
       s = a(1:nl, 1:nl)
       return
    end if
    allocate(piv(no, no))
    call sym_pinv(a(nl+1:nd, nl+1:nd), piv)
    s = a(1:nl, 1:nl) - matmul(a(1:nl, nl+1:nd), matmul(piv, a(nl+1:nd, 1:nl)))
  end subroutine schur_onto_leading

  !> Moore-Penrose pseudoinverse of a symmetric matrix via LAPACK dsyev,
  !! relative eigenvalue threshold (the translation of pinv in the
  !! Octave/numpy references).
  subroutine sym_pinv(a, p)
    real(rp), intent(in) :: a(:,:)
    real(rp), intent(out) :: p(:,:)
    real(rp), allocatable :: v(:,:), w(:), work(:)
    real(rp) :: tol
    integer(i4) :: n, info, lwork, i, j, k
    n = size(a, 1)
    if (n .eq. 0) return
    allocate(v(n, n), w(n))
    v = a
    lwork = max(64 * n, 1)
    allocate(work(lwork))
    call dsyev('V', 'U', n, v, n, w, work, lwork, info)
    if (info .ne. 0) call neko_error('macro_coarsen: dsyev failed')
    tol = PINV_RTOL * max(maxval(abs(w)), tiny(1.0_rp))
    p = 0.0_rp
    do k = 1, n
       if (abs(w(k)) .gt. tol) then
          do j = 1, n
             do i = 1, n
                p(i, j) = p(i, j) + v(i, k) * v(j, k) / w(k)
             end do
          end do
       end if
    end do
  end subroutine sym_pinv

  !> Smallest eigenvalue of a symmetric matrix (LAPACK dsyev, eigenvectors
  !! not needed so 'N'). Debug-check helper: SPSD requires this to be >=
  !! -tol for every locally-assembled/coarsened matrix.
  function min_eigenvalue(a) result(lam_min)
    real(rp), intent(in) :: a(:,:)
    real(rp) :: lam_min
    real(rp), allocatable :: v(:,:), w(:), work(:)
    integer(i4) :: n, info, lwork
    n = size(a, 1)
    if (n .eq. 0) then
       lam_min = 0.0_rp
       return
    end if
    allocate(v(n, n), w(n))
    v = a
    lwork = max(64 * n, 1)
    allocate(work(lwork))
    call dsyev('N', 'U', n, v, n, w, work, lwork, info)
    if (info .ne. 0) call neko_error('macro_coarsen: dsyev failed (min_eigenvalue)')
    lam_min = minval(w)
  end function min_eigenvalue

  ! ---- housekeeping ----



  pure subroutine sort_i4(a)
    integer(i4), intent(inout) :: a(:)
    integer(i4) :: i, j, t
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
  end subroutine sort_i4

end module amge_coarsen
