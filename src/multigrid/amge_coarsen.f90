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
  use comm, only : NEKO_COMM, pe_size, pe_rank
  use mpi_f08, only : MPI_Barrier
  use amge_topology, only : macro_topology_t, macro_mesh_t, &
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

  !> Toggle for level-0 face-matching frontier growth in
  !! agglomerate_face_growth (see that subroutine's header): grows each
  !! macroelement with the same connectivity-score frontier method as
  !! always, but seeds each new one from a previously-finished
  !! macroelement's WHOLE exposed face (not a single point) and, when
  !! growth touches an already-finished neighbor, snaps its own boundary
  !! to match that neighbor's whole face instead of a partial, jagged
  !! touch. Only fires where lvl%mmsh%elm_loc_face is allocated (level 0
  !! -- coarser levels' macroelements have no fixed face template to
  !! walk); everywhere else this degrades to plain single-element
  !! lexicographic seeding (still correct, just without the face-
  !! matching benefit). Default on; disable for A/B comparison.
  logical, public :: amge_use_face_growth = .true.

  !> target_size is a WEAK bound in the ordinary bscore growth loop
  !! (greedy_grow_from_seedset): once a macroelement reaches target_size,
  !! a further frontier candidate is still absorbed if it touches at
  !! least this many already-included elements -- e.g. a leftover element
  !! wedged against 3 faces of an otherwise-finished macroelement is a
  !! bad fit anywhere else (whatever it doesn't join here can match at
  !! most its remaining faces elsewhere), so it's better folded in even
  !! if that overshoots target_size than left to become an awkward
  !! small/singleton macroelement of its own. Does NOT apply to the
  !! face-matching snap-in (see greedy_grow_from_seedset's header),
  !! which has no upper limit at all -- only to the plain connectivity-
  !! score competition. Growth via this weak bound still hard-stops at
  !! amge_overgrow_factor * target_size (see below).
  integer(i4), public :: amge_nice_shape_min_touch = 3

  !> Hard ceiling on a macroelement's size, as a multiple of target_size,
  !! applied only to the "weak bound" growth amge_nice_shape_min_touch
  !! allows past target_size in the ordinary bscore loop -- a safety
  !! valve against a long chain of well-connected candidates growing one
  !! macroelement unboundedly. The face-matching snap-in is unaffected by
  !! this (see amge_nice_shape_min_touch).
  real(rp), public :: amge_overgrow_factor = 2.0_rp

  !> Hard cap on the number of DISTINCT already-finished neighbors a
  !! single growing macroelement will snap-absorb a face patch from (see
  !! greedy_grow_from_seedset's header) -- not a cap on absorbed element
  !! count, which stays uncapped per neighbor once snapped into. Prevents
  !! an extended seam of thin neighbors from being zippered through one
  !! bite at a time; 3 mirrors the routine's own motivating "wedged
  !! against 3 faces" scenario.
  integer(i4), public :: amge_max_snap_neighbors = 3

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
  !! shared facet). At level 0, growth is driven by
  !! agglomerate_face_growth (face-matching frontier growth -- see its
  !! header); everywhere else (no fixed face template to walk) it
  !! degrades to the same connectivity-score frontier method, seeded one
  !! element at a time in lexicographic order.
  !! (Integration note: tree_amg's greedy aggregation plays this role.)
  subroutine agglomerate_level(lvl, target_size, part, n_macro)
    type(amge_level_t), intent(in) :: lvl
    integer(i4), intent(in) :: target_size
    integer(i4), allocatable, intent(inout) :: part(:)
    integer(i4), intent(out) :: n_macro
    integer(i4), allocatable :: adj_ptr(:), adj_idx(:), fill(:)
    integer(i4) :: ne, f, e, a, b
    integer(i4), allocatable :: sizes(:)   !< debug: final size of each cluster

    ne = lvl%mmsh%n_elem
    if (allocated(part)) deallocate(part)
    allocate(part(ne))

    ! ---- element dual graph: e1 -- e2 iff they share an interior facet
    ! (face_nel==2). Two-pass CSR build: first pass counts each element's
    ! degree into adj_ptr and prefix-sums it into offsets, second pass
    ! (using fill(:) as a per-element running write cursor) scatters the
    ! actual neighbor ids into adj_idx at those offsets. A facet with
    ! face_nel==1 is a mesh/rank boundary and contributes no edge.
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

    part = 0
    n_macro = 0
    allocate(sizes(ne))   ! upper bound on n_macro; only sizes(1:n_macro) used

    call agglomerate_face_growth(lvl, adj_ptr, adj_idx, target_size, part, &
         n_macro, sizes)

    if (AMGE_DEBUG_CHECKS) call report_size_stats('agglomeration: elements/macroelement', &
         sizes(1:n_macro))
  end subroutine agglomerate_level

  !> Level-0 face-matching frontier growth (replaces the old fixed-shape
  !! brick growth): grows each macroelement via greedy_grow_from_seedset
  !! (the same connectivity-score frontier method used throughout this
  !! module), but seeds every macroelement after the first from the
  !! WHOLE exposed-face footprint of a previously-finished macroelement
  !! (via face_patch, other_id=0) rather than a single point, and snaps a
  !! growing macroelement's own boundary to match a neighbor's whole face
  !! when it bumps into one (also via face_patch, handled inside
  !! greedy_grow_from_seedset) instead of a partial, jagged touch.
  !!
  !! `face_q_e`/`face_q_slot` is a queue of (element, slot) pairs still
  !! facing unassigned territory, pushed after each macroelement finishes
  !! (one entry per exposed slot of its own members -- redundant entries
  !! from the same underlying face patch are harmless, just some wasted
  !! re-computation, since face_patch re-checks "still unassigned" and
  !! simply yields nothing the second time). Popping an entry re-verifies
  !! it's still valid (may have been claimed meanwhile by a different
  !! path) before calling face_patch.
  !!
  !! Whenever the queue is empty (including the very first call, and any
  !! time a whole connected region has been fully claimed), the next
  !! macroelement is bootstrapped from the next unassigned element in
  !! lexicographic order instead -- this is what makes the loop cover
  !! every element with no separate fallback pass needed afterward, and
  !! is also the ONLY growth mode at levels >= 1 (elm_loc_face never
  !! allocated there), reproducing the plain single-element frontier
  !! method exactly.
  !!
  !! No-op on the face-matching side (falls straight through to plain
  !! lexicographic bootstrap seeding throughout) unless
  !! lvl%mmsh%elm_loc_face is allocated (level 0 only).
  subroutine agglomerate_face_growth(lvl, adj_ptr, adj_idx, target_size, &
       part, n_macro, sizes)
    type(amge_level_t), intent(in) :: lvl
    integer(i4), intent(in) :: adj_ptr(:), adj_idx(:)
    integer(i4), intent(in) :: target_size
    integer(i4), intent(inout) :: part(:)
    integer(i4), intent(inout) :: n_macro
    integer(i4), intent(inout) :: sizes(:)
    integer(i4) :: ne, s, cnt, e0, dslot, npatch, k, nseed, i, cursor
    integer(i4), allocatable :: stamp(:), queue(:), patch(:), far(:)
    integer(i4), allocatable :: worklist(:), seed_set(:)
    integer(i4), allocatable :: face_q_e(:), face_q_slot(:)
    integer(i4) :: epoch, patch_cap, qh, qt, fe, fslot
    logical :: use_faces

    ne = lvl%mmsh%n_elem
    use_faces = amge_use_face_growth .and. allocated(lvl%mmsh%elm_loc_face)

    patch_cap = max(4 * target_size, 8)
    allocate(stamp(ne)); stamp = 0; epoch = 0
    allocate(queue(ne))
    allocate(patch(patch_cap), far(patch_cap))
    ! Sized to ne, not target_size: the face-matching snap-in in
    ! greedy_grow_from_seedset has no upper limit (see its header), so a
    ! macroelement's size is no longer bounded by target_size or even
    ! amge_overgrow_factor * target_size.
    allocate(worklist(max(ne, 1)))
    allocate(seed_set(patch_cap))
    allocate(face_q_e(6 * ne), face_q_slot(6 * ne))
    qh = 1; qt = 0
    cursor = 1

    do
       nseed = 0
       if (use_faces) then
          do while (qt .ge. qh)
             e0 = face_q_e(qh); dslot = face_q_slot(qh); qh = qh + 1
             if (part(e0) .eq. 0) cycle   ! shouldn't happen, be safe
             call step_neighbor(lvl%mmsh, lvl%mmsh%elm_loc_face, e0, &
                  dslot, fe, fslot)
             if (fe .ne. 0) cycle   ! no longer exposes "unassigned"
             call face_patch(lvl%mmsh, lvl%mmsh%elm_loc_face, part, e0, &
                  dslot, part(e0), 0, stamp, epoch, queue, patch, far, npatch)
             if (npatch .eq. 0) cycle
             nseed = 0
             do k = 1, npatch
                if (far(k) .ne. 0) then
                   if (part(far(k)) .eq. 0) then
                      nseed = nseed + 1
                      seed_set(nseed) = far(k)
                   end if
                end if
             end do
             if (nseed .gt. 0) exit
          end do
       end if

       if (nseed .eq. 0) then
          ! bootstrap: next unassigned element in lexicographic order
          ! (cursor only ever advances, so this is O(ne) amortized over
          ! the whole run, not per call)
          s = 0
          do while (cursor .le. ne)
             if (part(cursor) .eq. 0) then
                s = cursor
                exit
             end if
             cursor = cursor + 1
          end do
          if (s .eq. 0) exit   ! every element assigned -- done
          nseed = 1; seed_set(1) = s
       end if

       n_macro = n_macro + 1
       call greedy_grow_from_seedset(lvl, adj_ptr, adj_idx, target_size, &
            part, n_macro, seed_set, nseed, stamp, epoch, queue, patch, &
            far, worklist, cnt)
       sizes(n_macro) = cnt

       if (use_faces) then
          do i = 1, cnt
             e0 = worklist(i)
             do k = 1, 6
                call step_neighbor(lvl%mmsh, lvl%mmsh%elm_loc_face, e0, k, &
                     fe, fslot)
                if (fe .ne. 0) cycle
                qt = qt + 1
                face_q_e(qt) = e0
                face_q_slot(qt) = k
             end do
          end do
       end if
    end do
  end subroutine agglomerate_face_growth

  !> Grow one macroelement (id `new_id`) from `seed_set` via the SAME
  !! best-connected-neighbor frontier growth used throughout this module,
  !! generalized to (a) start from multiple seed elements at once (all
  !! immediately absorbed, before any frontier scoring), and (b) snap to
  !! a neighboring macroelement's whole face when growth touches one:
  !! whenever a newly-absorbed element's neighbor is found to already
  !! belong to some OTHER macroelement C (not 0, not new_id), look up C's
  !! whole face patch there (face_patch, own=C, other_id=new_id) via
  !! slot_toward, and absorb EVERY still-unassigned far-side element of
  !! that patch, unconditionally -- a partial face match would defeat
  !! the whole point of snapping to it, so this has no target_size cap
  !! at all; the only per-neighbor bound is face_patch's own patch_cap
  !! safety valve -- then scan their own neighbors too, via `worklist`
  !! (so absorption can cascade: a snapped-in element may itself bump
  !! into yet another macroelement). An element that turns out to
  !! already belong to a THIRD macroelement at a genuine 3-way junction
  !! is simply left alone -- an unavoidable corner, not a bug.
  !!
  !! Separately, the NUMBER OF DISTINCT neighbors snapped into (as
  !! opposed to elements absorbed per neighbor) IS capped, at
  !! amge_max_snap_neighbors: without this, a growing macroelement
  !! sitting along an extended seam of thin already-finished neighbors
  !! can zipper through dozens of them, one small bite at a time, wildly
  !! overshooting target_size (observed: 103 elements from 34 distinct
  !! neighbors on one case). A repeat bump into an already-snapped
  !! neighbor doesn't count against this cap -- only encountering a new
  !! one after the cap is reached is skipped, leaving that neighbor's
  !! patch with its current owner.
  !!
  !! The ordinary bscore competition, in contrast, treats target_size as
  !! a WEAK bound: once cnt reaches target_size, the best remaining
  !! frontier candidate is still absorbed if it's well-connected enough
  !! (should_accept, see amge_nice_shape_min_touch's header), up to a
  !! hard ceiling of amge_overgrow_factor * target_size.
  !!
  !! On return, `worklist(1:cnt)` holds every element absorbed into this
  !! macroelement, in absorption order -- the caller reuses this
  !! directly as the finished macroelement's member list, no separate
  !! scan needed. `worklist` must be sized >= lvl%mmsh%n_elem by the
  !! caller: since the snap-in has no upper limit, cnt is no longer
  !! bounded by target_size or even amge_overgrow_factor * target_size.
  subroutine greedy_grow_from_seedset(lvl, adj_ptr, adj_idx, target_size, &
       part, new_id, seed_set, n_seed, stamp, epoch, queue, patch, far, &
       worklist, cnt)
    type(amge_level_t), intent(in) :: lvl
    integer(i4), intent(in) :: adj_ptr(:), adj_idx(:)
    integer(i4), intent(in) :: target_size
    integer(i4), intent(inout) :: part(:)
    integer(i4), intent(in) :: new_id
    integer(i4), intent(in) :: seed_set(:), n_seed
    integer(i4), intent(inout) :: stamp(:), epoch, queue(:), patch(:), far(:)
    integer(i4), intent(inout) :: worklist(:)
    integer(i4), intent(out) :: cnt
    integer(i4), allocatable :: frontier(:)
    logical, allocatable :: infr(:)
    logical, allocatable :: snapped(:)
    integer(i4) :: ne, nfr, best, bscore, sc, a, q, e, c, si, hard_max
    integer(i4) :: wn, wi, v, dslot, npatch, k, n_snapped
    logical :: use_faces

    ne = lvl%mmsh%n_elem
    use_faces = allocated(lvl%mmsh%elm_loc_face)
    hard_max = max(target_size, &
         nint(real(target_size, rp) * amge_overgrow_factor))
    allocate(frontier(ne), infr(ne))
    infr = .false.
    nfr = 0
    cnt = 0
    wn = 0
    allocate(snapped(ne)); snapped = .false.; n_snapped = 0

    do si = 1, n_seed
       part(seed_set(si)) = new_id
       cnt = cnt + 1
       wn = wn + 1; worklist(wn) = seed_set(si)
    end do

    wi = 1
    do
       do while (wi .le. wn)
          v = worklist(wi); wi = wi + 1
          do q = adj_ptr(v) + 1, adj_ptr(v + 1)
             e = adj_idx(q)
             if (part(e) .eq. 0) then
                if (.not. infr(e)) then
                   nfr = nfr + 1; frontier(nfr) = e; infr(e) = .true.
                end if
             else if (use_faces .and. part(e) .ne. new_id) then
                ! cap the number of DISTINCT neighbors we'll snap into
                ! (amge_max_snap_neighbors) -- a repeat bump into a
                ! neighbor already snapped is unaffected (still folds in
                ! its whole patch, uncapped in element count); only a
                ! genuinely new neighbor past the cap is skipped, leaving
                ! its patch with its current owner.
                if (.not. snapped(part(e))) then
                   if (n_snapped .ge. amge_max_snap_neighbors) cycle
                   snapped(part(e)) = .true.
                   n_snapped = n_snapped + 1
                end if
                dslot = slot_toward(lvl%mmsh, lvl%mmsh%elm_loc_face, e, v)
                if (dslot .eq. 0) cycle
                call face_patch(lvl%mmsh, lvl%mmsh%elm_loc_face, part, e, &
                     dslot, part(e), new_id, stamp, epoch, queue, patch, &
                     far, npatch)
                do k = 1, npatch
                   if (far(k) .ne. 0) then
                      if (part(far(k)) .eq. 0) then
                         part(far(k)) = new_id
                         cnt = cnt + 1
                         wn = wn + 1; worklist(wn) = far(k)
                      end if
                   end if
                end do
             end if
          end do
       end do
       if (nfr .eq. 0) exit
       ! grow one element at a time: score every frontier candidate by
       ! how many of ITS neighbors are already in the current
       ! macroelement (sc), absorb the best-connected one (best/bscore)
       ! -- this is what keeps macroelements compact/blob-shaped rather
       ! than growing thin tendrils. target_size is a WEAK bound here:
       ! once reached, the winner is still absorbed if it's well
       ! connected enough (should_accept), up to hard_max.
       best = 0; bscore = -1
       do a = 1, nfr
          e = frontier(a)
          if (part(e) .ne. 0) cycle
          sc = 0
          do q = adj_ptr(e) + 1, adj_ptr(e + 1)
             if (part(adj_idx(q)) .eq. new_id) sc = sc + 1
          end do
          if (sc .gt. bscore) then
             bscore = sc; best = a
          end if
       end do
       if (best .eq. 0) exit
       if (.not. should_accept(cnt, target_size, hard_max, bscore, &
            amge_nice_shape_min_touch)) exit
       c = frontier(best)
       frontier(best) = frontier(nfr); nfr = nfr - 1
       part(c) = new_id
       cnt = cnt + 1
       wn = wn + 1; worklist(wn) = c
    end do
  end subroutine greedy_grow_from_seedset

  !> target_size acceptance rule for the ordinary bscore growth loop
  !! (see greedy_grow_from_seedset, and amge_nice_shape_min_touch/
  !! amge_overgrow_factor's headers for the rationale): always accept
  !! below target_size; past it, accept only a candidate connected
  !! enough to be a clearly better fit here than elsewhere (sc >=
  !! min_touch), and only up to hard_max.
  pure function should_accept(cnt, target_size, hard_max, sc, min_touch) &
       result(accept)
    integer(i4), intent(in) :: cnt, target_size, hard_max, sc, min_touch
    logical :: accept

    accept = (cnt .lt. target_size) .or. &
         (cnt .lt. hard_max .and. sc .ge. min_touch)
  end function should_accept

  !> Find the local slot of element e whose neighbor is v (e, v assumed
  !! face-adjacent, e.g. via the dual graph adj_ptr/adj_idx) -- a small
  !! O(6) scan, used to recover "the direction back toward v" needed to
  !! call face_patch when growth bumps into v's macroelement. Returns 0
  !! if no such slot is found (shouldn't happen given e, v are adjacent,
  !! but callers treat 0 as "skip" defensively).
  function slot_toward(mmsh, elm_loc_face, e, v) result(slot)
    type(macro_mesh_t), intent(in) :: mmsh
    integer(i4), intent(in) :: elm_loc_face(:,:)
    integer(i4), intent(in) :: e, v
    integer(i4) :: slot
    integer(i4) :: f, oe

    do slot = 1, 6
       f = elm_loc_face(slot, e)
       if (mmsh%face_nel(f) .eq. 2) then
          oe = mmsh%face_el(1, f)
          if (oe .eq. e) oe = mmsh%face_el(2, f)
          if (oe .eq. v) return
       end if
    end do
    slot = 0
  end function slot_toward

  !> Flood-fill the flat patch of `own`'s own elements (own is always a
  !! real macroelement id here, never "unassigned") that all expose the
  !! SAME `other_id` (0 = boundary/unassigned; otherwise another
  !! macroelement's id) through the SAME local slot `dslot`, starting
  !! from seed element `e0`. Walks the 4 in-plane slots (the two
  !! axis-pairs other than dslot's own pair) via step_neighbor,
  !! restricted to elements with part==own; stops expanding a direction
  !! when it leaves `own` (that's own's own edge or corner) or its dslot
  !! no longer exposes the same other_id (a genuine change of what's
  !! beyond, not just a pass-through point).
  !!
  !! Reuses the same "trust a fresh element's numeric slot to mean the
  !! same physical direction" assumption the old brick-growth code relied
  !! on at its branch points, valid here without extra bookkeeping
  !! because at level 0 the slot-pairing {1,2}/{3,4}/{5,6} is a FIXED
  !! template invariant (macro_mesh_init_hex), not something that needs
  !! per-element rediscovery.
  !!
  !! `stamp`/`epoch` are a reusable "visited" marker and `queue` a
  !! reusable BFS scratch array (both sized >= mmsh%n_elem, owned by the
  !! top-level caller) so repeated calls -- one per queued face, one per
  !! bump-into event -- don't each pay an O(n_elem) allocation. `patch`/
  !! `far` are caller-sized; a patch that would exceed that size aborts
  !! (n=0) as a cheap safety valve for "this region isn't reliably
  !! locally structured" -- same philosophy as the old brick-growth
  !! duplicate-abort; callers fall back to single-element behavior.
  !!
  !! Two uses: (1) seeding -- own = a just-finished macroelement,
  !! other_id = 0, far (deduplicated, re-checked still-unassigned)
  !! becomes the next macroelement's seed set; (2) snapping -- own = the
  !! macroelement just bumped into, other_id = the growing macroelement's
  !! own id, e0 = the contacted element, dslot = its slot facing back
  !! toward the growing macroelement, far = the growing macroelement's
  !! own candidate elements to pull in directly.
  subroutine face_patch(mmsh, elm_loc_face, part, e0, dslot, own, other_id, &
       stamp, epoch, queue, patch, far, n)
    type(macro_mesh_t), intent(in) :: mmsh
    integer(i4), intent(in) :: elm_loc_face(:,:)
    integer(i4), intent(in) :: part(:)
    integer(i4), intent(in) :: e0, dslot, own, other_id
    integer(i4), intent(inout) :: stamp(:)
    integer(i4), intent(inout) :: epoch
    integer(i4), intent(inout) :: queue(:)
    integer(i4), intent(out) :: patch(:), far(:)
    integer(i4), intent(out) :: n
    integer(i4), parameter :: OPP(6) = [2, 1, 4, 3, 6, 5]
    integer(i4) :: qhead, qtail, v, w, wslot, fe, fslot, p, v_other, cap

    epoch = epoch + 1
    cap = size(patch)
    qhead = 1; qtail = 1
    queue(1) = e0
    stamp(e0) = epoch
    n = 0
    do while (qhead .le. qtail)
       v = queue(qhead); qhead = qhead + 1
       call step_neighbor(mmsh, elm_loc_face, v, dslot, fe, fslot)
       v_other = 0
       if (fe .ne. 0) v_other = part(fe)
       ! Accept a far side that's EITHER already other_id or still
       ! unassigned (0) -- reject only a genuine third macroelement (the
       ! 3-way-junction case). For the seeding call (other_id=0) this is
       ! unchanged (both conditions collapse to the same test); for the
       ! snapping call (other_id=new_id) this is what actually lets
       ! still-unassigned far-side elements into the patch, rather than
       ! only ones already claimed by new_id (which face_patch would
       ! then be handing back as "new" to absorb -- a no-op).
       if (v_other .ne. other_id .and. v_other .ne. 0) cycle
       n = n + 1
       if (n .gt. cap) then
          n = 0
          return
       end if
       patch(n) = v
       far(n) = fe
       do p = 1, 6
          if (p .eq. dslot .or. p .eq. OPP(dslot)) cycle
          call step_neighbor(mmsh, elm_loc_face, v, p, w, wslot)
          if (w .eq. 0) cycle
          if (part(w) .ne. own) cycle
          if (stamp(w) .eq. epoch) cycle
          stamp(w) = epoch
          qtail = qtail + 1
          queue(qtail) = w
       end do
    end do
  end subroutine face_patch

  !> Step from element e via local face-slot `slot` (1..6, loc_face
  !! ordering {1,2}=-x/+x, {3,4}=-y/+y, {5,6}=-z/+z, see
  !! amge_topology.f90's macro_mesh_init_hex header) to its neighbor,
  !! re-discovering the neighbor's own slot for the shared face so the
  !! caller can continue walking the SAME physical direction via that
  !! slot's opposite. This is what makes the walk robust to per-element
  !! local-numbering differences: only each element's own opposite-face
  !! pairing is assumed (always true of any hex), never that two
  !! face-adjacent elements number their local slots the same way. Used
  !! by face_patch's flood-fill and by slot_toward/
  !! greedy_grow_from_seedset's bump-into handling.
  !! next_e=0 at a physical/rank boundary (face_nel==1) OR when the
  !! neighbor landed on has no matching slot in its own elm_loc_face row
  !! -- only possible if elm_loc_face is partially populated, which
  !! doesn't currently happen (level 0's elm_loc_face, from
  !! macro_mesh_init_hex, is always fully populated); kept as a graceful
  !! "treat like a boundary" fallback rather than an error regardless.
  subroutine step_neighbor(mmsh, elm_loc_face, e, slot, next_e, next_slot)
    type(macro_mesh_t), intent(in) :: mmsh
    integer(i4), intent(in) :: elm_loc_face(:,:)
    integer(i4), intent(in) :: e, slot
    integer(i4), intent(out) :: next_e, next_slot
    integer(i4), parameter :: OPP(6) = [2, 1, 4, 3, 6, 5]
    integer(i4) :: f, kk

    f = elm_loc_face(slot, e)
    if (mmsh%face_nel(f) .eq. 1) then
       next_e = 0; next_slot = 0
       return
    end if
    next_e = mmsh%face_el(1, f)
    if (next_e .eq. e) next_e = mmsh%face_el(2, f)
    do kk = 1, 6
       if (elm_loc_face(kk, next_e) .eq. f) then
          next_slot = OPP(kk)
          return
       end if
    end do
    next_e = 0; next_slot = 0
  end subroutine step_neighbor

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
    integer(i4), allocatable :: dsizes(:)   !< debug: coarse dofs per macroelement
    integer :: ierr   !< debug: MPI_Barrier status after [check] prints
    logical :: ghosted
    type(schur_check_t) :: schk
    type(amge_ghost_t) :: gh
    type(amge_level_t) :: lvl_ext
    integer(i4), allocatable :: part_ext(:)
    integer(i4), allocatable :: p_vtx_ptr(:), p_vtx_idx(:)
    integer(i4), allocatable :: p_face_ptr(:), p_fv_ptr(:), p_fv_idx(:)
    integer(i4), allocatable :: p_face_anchor(:,:)
    integer(i4), allocatable :: p_fe_ptr(:)
    integer(i4), allocatable :: p_fe_vtx(:,:)
    integer(i4), allocatable :: p_fe_anchor(:,:)
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
            p_face_ptr, p_fv_ptr, p_fv_idx, p_face_anchor, p_fe_ptr, &
            p_fe_vtx, p_fe_anchor, p_amat_ptr, p_amat)
       call amge_ghost_exchange(NEKO_COMM%mpi_val, p_vtx_ptr, p_vtx_idx, &
            p_face_ptr, p_fv_ptr, p_fv_idx, p_face_anchor, p_fe_ptr, &
            p_fe_vtx, p_fe_anchor, p_amat_ptr, p_amat, part, n_macro, gh)
       call amge_ghost_part_ext(gh, part, part_ext)
       moff = gh%macro_offset
       n_macro_topo = gh%n_macro_global

       ! splice the ghost elements' own face/edge/vertex sub-tables (just
       ! received) onto a copy of this level's own mesh -- NOT a rebuild
       ! from a raw vertex list, which only works at level 0 (see
       ! amge_ghost.f90's header and macro_mesh_splice_ghost's doc comment)
       call macro_mesh_splice_ghost(lvl%mmsh, gh%n_ghost, gh%vtx_ptr, &
            gh%vtx_idx, gh%face_ptr, gh%face_vtx_ptr, gh%face_vtx_idx, &
            gh%face_anchor, gh%face_edge_ptr, gh%face_edge_vtx, &
            gh%face_edge_anchor, lvl_ext%mmsh)

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
    allocate(dsizes(n_macro))
    do m = 1, n_macro
       call macro_dof_list(lvl, melems(m)%x, markv, tr%maps(m)%fdofs)
       nc = 0
       do a = 1, size(tr%maps(m)%fdofs)
          if (topo%is_mv(tr%maps(m)%fdofs(a))) nc = nc + 1
       end do
       lvlC%elm_vtx_ptr(m + 1) = lvlC%elm_vtx_ptr(m) + nc
       dsizes(m) = nc
    end do
    allocate(lvlC%elm_vtx_idx(lvlC%elm_vtx_ptr(n_macro + 1)))
    if (AMGE_DEBUG_CHECKS) call report_size_stats('coarse dofs/macroelement', dsizes)

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
       write(*, '("   [check] rank ", I0, ": per-macroelement (n=", I0, "): S sym = ", ES10.3, &
            & "  harmonic resid = ", ES10.3, "  |PiT^TAPi - Ac| = ", ES10.3, &
            & "  min eig(Ac) = ", ES10.3)') pe_rank, schk%n_checked, schk%sym_err, &
            schk%harmonic_err, schk%galerkin_err, schk%min_eig
    end if
    ! gated on AMGE_DEBUG_CHECKS alone (uniform across ranks), NOT on
    ! schk%n_checked (which can vary per rank) -- every rank must reach
    ! this call or the barrier deadlocks
    if (AMGE_DEBUG_CHECKS) call MPI_Barrier(NEKO_COMM, ierr)
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
       face_ptr, face_vtx_ptr, face_vtx_idx, face_anchor, face_edge_ptr, &
       face_edge_vtx, face_edge_anchor, amat_ptr, amat)
    type(amge_level_t), intent(in) :: lvl
    integer(i4), allocatable, intent(out) :: vtx_ptr(:), vtx_idx(:)
    integer(i4), allocatable, intent(out) :: face_ptr(:), face_vtx_ptr(:)
    integer(i4), allocatable, intent(out) :: face_vtx_idx(:)
    integer(i4), allocatable, intent(out) :: face_anchor(:,:)
    integer(i4), allocatable, intent(out) :: face_edge_ptr(:)
    integer(i4), allocatable, intent(out) :: face_edge_vtx(:,:)
    integer(i4), allocatable, intent(out) :: face_edge_anchor(:,:)
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
    allocate(face_anchor(4, face_ptr(nelv + 1)))
    allocate(face_edge_vtx(2, face_edge_ptr(face_ptr(nelv + 1) + 1)))
    allocate(face_edge_anchor(2, face_edge_ptr(face_ptr(nelv + 1) + 1)))

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
          face_anchor(:, slot) = lvl%mmsh%face_anchor(:, f)
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
             face_edge_anchor(:, face_edge_ptr(slot) + k) = lvl%mmsh%edge_anchor(:, g)
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
    integer :: ierr
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
       write(*, '("   [check] rank ", I0, ": Q_E (n=", I0, "): S_E sym = ", ES10.3, &
            & "  min eig(S_E) = ", ES10.3)') pe_rank, topo%n_medge, sym_err, min_eig
    end if
    ! gated on AMGE_DEBUG_CHECKS alone (uniform across ranks), NOT on
    ! topo%n_medge (which can vary per rank)
    if (AMGE_DEBUG_CHECKS) call MPI_Barrier(NEKO_COMM, ierr)
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
    integer :: ierr
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
       write(*, '("   [check] rank ", I0, ": Q_F (n=", I0, "): S_F sym = ", ES10.3, &
            & "  min eig(S_F) = ", ES10.3)') pe_rank, n_checked, sym_err, min_eig
    end if
    ! gated on AMGE_DEBUG_CHECKS alone (uniform across ranks), NOT on
    ! n_checked (which can vary per rank)
    if (AMGE_DEBUG_CHECKS) call MPI_Barrier(NEKO_COMM, ierr)
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

  !> Debug helper: print min/avg/median/max of an integer size array
  !! (aggregate sizes, coarse dofs per macroelement, ...) in one common
  !! format, gated by the caller on AMGE_DEBUG_CHECKS. Both call sites
  !! gate on that same (rank-uniform) flag alone, so it's safe to place
  !! the debug MPI_Barrier here rather than at each call site.
  subroutine report_size_stats(label, sizes)
    character(len=*), intent(in) :: label
    integer(i4), intent(in) :: sizes(:)
    integer(i4), allocatable :: sorted(:)
    real(rp) :: avg, med
    integer(i4) :: n
    integer :: ierr
    n = size(sizes)
    sorted = sizes
    call sort_i4(sorted)
    avg = real(sum(sizes), rp) / real(n, rp)
    if (mod(n, 2) == 1) then
       med = real(sorted((n + 1) / 2), rp)
    else
       med = 0.5_rp * real(sorted(n / 2) + sorted(n / 2 + 1), rp)
    end if
    write(*, '("   [check] rank ", I0, ": ", A, " (n=", I0, "): min/avg/median/max = ", &
         & I0, "/", F6.2, "/", F6.2, "/", I0)') pe_rank, label, n, sorted(1), avg, med, sorted(n)
    call MPI_Barrier(NEKO_COMM, ierr)
  end subroutine report_size_stats

end module amge_coarsen
