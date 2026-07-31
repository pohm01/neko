! Copyright (c) 2026, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the conditions of the
! BSD 3-Clause License are met (see Neko's COPYING file).
!
!> One-layer ghost (halo) exchange that makes the rank-local AMGe
!! macroentity extraction and trace maps agree with a serial run, at
!! EVERY coarsening level (not just the finest).
!!
!! WHY A HALO IS NEEDED AT ALL. With agglomeration kept rank-local, no
!! macroelement crosses a rank boundary, but two things on the interface
!! still go wrong without communication:
!!   (1) a facet whose only local element is mine looks like a physical
!!       boundary, so it is labelled (0,m) instead of (m,m_remote). The
!!       skeleton/macrovertex rules read those labels, so the two sides
!!       can disagree about which fine vertices are coarse dofs;
!!   (2) the trace maps Q_E / Q_F minimize energy over the elements
!!       touching the entity, which straddle the interface. Computed from
!!       half the neighborhood on each side, the two sides produce
!!       DIFFERENT maps for the SAME macroedge/macroface, so the
!!       prolongation is discontinuous across the rank boundary and
!!       trace compatibility -- the property the whole construction rests
!!       on -- is lost.
!!
!! WHY ONE VERTEX-LAYER IS EXACTLY ENOUGH. Every element incident to an
!! interface vertex, edge or facet shares a vertex with some local
!! element, hence is a ghost under a vertex-based halo. The trace-map
!! neighborhoods are "elements sharing >= 2 chain vertices" (Q_E) and
!! ">= 3 patch vertices" (Q_F): both consist of elements sharing vertices
!! with the entity, so they are fully contained in one vertex layer.
!! Nothing deeper is ever consulted.
!!
!! WHY THE PAYLOAD IS A CSR SUB-TABLE, NOT A RAW VERTEX LIST. At the
!! finest level a "ghost element" is a hex: its facet/edge structure is
!! derivable from its 8 vertices via a fixed template. Above level 0, a
!! macroelement's faces/edges are NOT derivable from its vertex set at
!! all -- they exist only because they were inherited from the level
!! below (macro_topology_t%build_next_level[_owned]). So instead of
!! shipping a ghost element's raw vertex list and re-deriving its
!! topology, each element being ghosted ships its OWN already-computed
!! face/edge sub-tables (see build_ghost_payload_local in
!! amge_coarsen.f90, which extracts them from that level's own mmsh) and
!! the receiver grafts them onto its own mesh by vertex-set/vertex-pair
!! IDENTITY (see macro_mesh_splice_ghost in amge_topology.f90), never by
!! re-deriving structure from a template. This subsumes the finest-level
!! hex case too (nc=8, template already baked into lvl(0)%mmsh) -- there
!! is only one code path, for every level.
!!
!! HOW THIS MODULE IS USED. It returns ghost element data (CSR vertex,
!! face and edge sub-tables, local matrices, global macro ids) to graft
!! onto the caller's own mesh via macro_mesh_splice_ghost. The only
!! changes needed downstream are: agglomerate only the LOCAL elements,
!! and emit macroelements only for locally OWNED macro ids (see
!! amge_ghost_owned_mask and macro_topology_build_next_level_owned).
!!
!! Communication is a standard rendezvous: vertex ids are hashed to a
!! distributed directory which reports co-owners, then each pair of
!! neighbouring ranks exchanges the elements touching their shared
!! vertices, together with those elements' own face/edge sub-tables and
!! local matrices.
module amge_ghost
  use mpi
  use num_types, only : i4, rp
  use utils, only : neko_error
  implicit none
  private

  !> Ghost element data (CSR, variable arity), to be spliced onto the
  !! caller's own mesh by macro_mesh_splice_ghost.
  type, public :: amge_ghost_t
     integer(i4) :: n_ghost = 0
     integer(i4) :: n_local = 0
     integer(i4) :: n_macro_global = 0            !< total macros over all ranks
     integer(i4) :: macro_offset = 0              !< my local id 1 -> global id offset+1
     !> CSR (n_ghost+1): global vertex ids per ghost element
     integer(i4), allocatable :: vtx_ptr(:)
     integer(i4), allocatable :: vtx_idx(:)
     !> CSR (n_ghost+1): faces per ghost element (flat face-slot list)
     integer(i4), allocatable :: face_ptr(:)
     !> CSR (n_faces+1), over the flat face-slot list: global vertex ids
     !! per face
     integer(i4), allocatable :: face_vtx_ptr(:)
     integer(i4), allocatable :: face_vtx_idx(:)
     !> CSR (n_faces+1), over the flat face-slot list: bounding edges per
     !! face, given as global vertex-id endpoint PAIRS (edges have no
     !! separate identity of their own on the wire -- identity is a
     !! vertex pair, established on receipt)
     integer(i4), allocatable :: face_edge_ptr(:)
     integer(i4), allocatable :: face_edge_vtx(:,:)  !< (2, *)
     integer(i4), allocatable :: gmacro(:)         !< (n_ghost) global macro id
     !> CSR (n_ghost+1): flattened local matrix per ghost element
     !! (row-major/column-major matches matrix_t%x's own storage, size
     !! nv_e**2 where nv_e is that element's own vertex count)
     integer(i4), allocatable :: amat_ptr(:)
     real(rp), allocatable :: amat(:)
   contains
     procedure, pass(this) :: free => amge_ghost_free
  end type amge_ghost_t

  public :: amge_ghost_exchange, amge_ghost_part_ext, amge_ghost_owned_mask

contains

  !> Build the ghost layer. All ids in/out are GLOBAL mesh vertex ids.
  !! @param comm            MPI communicator (Neko: NEKO_COMM)
  !! @param vtx_ptr/vtx_idx CSR (nelv+1): global vertex ids per local element
  !! @param face_ptr        CSR (nelv+1): faces per local element (flat
  !!                        face-slot list)
  !! @param face_vtx_ptr/idx  CSR, over the flat face-slot list: global
  !!                        vertex ids per face
  !! @param face_edge_ptr/vtx  CSR, over the flat face-slot list:
  !!                        bounding-edge global vertex-id pairs per face
  !! @param amat_ptr/amat   CSR (nelv+1): flattened local matrix per
  !!                        local element
  !! @param part_local      (nelv) my LOCAL macro ids, 1..n_macro_local
  !! @param n_macro_local   number of macroelements I own
  !! @param gh              result
  subroutine amge_ghost_exchange(comm, vtx_ptr, vtx_idx, &
       face_ptr, face_vtx_ptr, face_vtx_idx, face_edge_ptr, face_edge_vtx, &
       amat_ptr, amat, part_local, n_macro_local, gh)
    integer, intent(in) :: comm
    integer(i4), intent(in) :: vtx_ptr(:), vtx_idx(:)
    integer(i4), intent(in) :: face_ptr(:), face_vtx_ptr(:), face_vtx_idx(:)
    integer(i4), intent(in) :: face_edge_ptr(:)
    integer(i4), intent(in) :: face_edge_vtx(:,:)
    integer(i4), intent(in) :: amat_ptr(:)
    real(rp), intent(in) :: amat(:)
    integer(i4), intent(in) :: part_local(:), n_macro_local
    type(amge_ghost_t), intent(inout) :: gh

    integer :: ierr, myrank, nrank, p
    integer(i4) :: nelv, i, j, k, e, nv, cnt, g, src, ndir, npair
    integer(i4), allocatable :: myv(:)                     ! my unique vertex gids
    integer, allocatable :: scnt(:), sdsp(:), rcnt(:), rdsp(:)
    integer, allocatable :: scnt2(:), sdsp2(:), rcnt2(:), rdsp2(:)
    integer(i4), allocatable :: sbuf(:), rbuf(:)
    integer(i4), allocatable :: dgid(:), dsrc(:), ord(:)
    logical, allocatable :: shared_with(:,:)               ! (nv, 0:nrank-1)
    logical, allocatable :: send_elem(:)
    integer(i4) :: offset
    ! ---- phase 4: variable-length per-element payload ----
    integer, allocatable :: ecnt(:), edsp(:), recnt(:), redsp(:)
    integer, allocatable :: icnt(:), idsp(:), ircnt(:), irdsp(:)
    integer, allocatable :: rlen(:), rdsp2r(:), rrcnt(:), rrdsp(:)
    integer(i4), allocatable :: eints(:), rints(:)
    real(rp), allocatable :: ereal(:), rreal(:)
    integer(i4) :: a, f, nvf, nef, nv_e, pos_i, pos_r
    integer(i4) :: n_gf_tot, n_gvref_tot, n_gfvref_tot, n_geref_tot
    integer(i4) :: idx, pi, pr_, gidx, gf, gfv, gfe

    call MPI_Comm_rank(comm, myrank, ierr)
    call MPI_Comm_size(comm, nrank, ierr)
    call gh%free()
    nelv = size(vtx_ptr) - 1
    gh%n_local = nelv

    ! ---------------- globally unique macro ids ----------------
    offset = 0
    call MPI_Exscan(n_macro_local, offset, 1, MPI_INTEGER, MPI_SUM, comm, ierr)
    if (myrank == 0) offset = 0
    gh%macro_offset = offset
    call MPI_Allreduce(n_macro_local, gh%n_macro_global, 1, MPI_INTEGER, &
         MPI_SUM, comm, ierr)

    ! ---------------- my unique vertex gids ----------------
    cnt = vtx_ptr(nelv + 1)
    allocate(myv(max(cnt, 1)))
    myv(1:cnt) = vtx_idx(1:cnt)
    call sort_i4(myv(1:cnt))
    nv = 0
    do i = 1, cnt
       if (i == 1 .or. myv(i) /= myv(max(i-1,1))) then
          nv = nv + 1
          myv(nv) = myv(i)
       end if
    end do

    ! ---------------- phase 1: vertices -> distributed directory ----------------
    allocate(scnt(nrank), sdsp(nrank), rcnt(nrank), rdsp(nrank))
    scnt = 0
    do i = 1, nv
       p = int(mod(int(myv(i) - 1, i4), int(nrank, i4)))
       scnt(p + 1) = scnt(p + 1) + 1
    end do
    call MPI_Alltoall(scnt, 1, MPI_INTEGER, rcnt, 1, MPI_INTEGER, comm, ierr)
    call make_displs(scnt, sdsp)
    call make_displs(rcnt, rdsp)
    allocate(sbuf(max(sum(scnt), 1)), rbuf(max(sum(rcnt), 1)))
    block
      integer, allocatable :: fill(:)
      allocate(fill(nrank)); fill = 0
      do i = 1, nv
         p = int(mod(int(myv(i) - 1, i4), int(nrank, i4)))
         fill(p + 1) = fill(p + 1) + 1
         sbuf(sdsp(p + 1) + fill(p + 1)) = myv(i)
      end do
    end block
    call MPI_Alltoallv(sbuf, scnt, sdsp, MPI_INTEGER, &
                       rbuf, rcnt, rdsp, MPI_INTEGER, comm, ierr)

    ! ---------------- phase 2: directory finds co-owners ----------------
    ndir = sum(rcnt)
    allocate(dgid(max(ndir,1)), dsrc(max(ndir,1)))
    do p = 1, nrank
       do i = 1, rcnt(p)
          dgid(rdsp(p) + i) = rbuf(rdsp(p) + i)
          dsrc(rdsp(p) + i) = p - 1
       end do
    end do
    allocate(ord(max(ndir,1)))
    call sort_index_i4(dgid, ndir, ord)
    allocate(scnt2(nrank), sdsp2(nrank), rcnt2(nrank), rdsp2(nrank))
    scnt2 = 0
    i = 1
    do while (i <= ndir)
       j = i
       do while (j < ndir)
          if (dgid(ord(j+1)) /= dgid(ord(i))) exit
          j = j + 1
       end do
       if (j > i) then
          do k = i, j
             scnt2(dsrc(ord(k)) + 1) = scnt2(dsrc(ord(k)) + 1) + 2 * (j - i)
          end do
       end if
       i = j + 1
    end do
    call MPI_Alltoall(scnt2, 1, MPI_INTEGER, rcnt2, 1, MPI_INTEGER, comm, ierr)
    call make_displs(scnt2, sdsp2)
    call make_displs(rcnt2, rdsp2)
    deallocate(sbuf, rbuf)
    allocate(sbuf(max(sum(scnt2),1)), rbuf(max(sum(rcnt2),1)))
    block
      integer, allocatable :: fill(:)
      integer(i4) :: aa, bb
      allocate(fill(nrank)); fill = 0
      i = 1
      do while (i <= ndir)
         j = i
         do while (j < ndir)
            if (dgid(ord(j+1)) /= dgid(ord(i))) exit
            j = j + 1
         end do
         do aa = i, j
            do bb = i, j
               if (aa == bb) cycle
               p = dsrc(ord(aa))
               fill(p + 1) = fill(p + 1) + 1
               sbuf(sdsp2(p + 1) + fill(p + 1)) = dgid(ord(aa))
               fill(p + 1) = fill(p + 1) + 1
               sbuf(sdsp2(p + 1) + fill(p + 1)) = dsrc(ord(bb))
            end do
         end do
         i = j + 1
      end do
    end block
    call MPI_Alltoallv(sbuf, scnt2, sdsp2, MPI_INTEGER, &
                       rbuf, rcnt2, rdsp2, MPI_INTEGER, comm, ierr)

    ! ---------------- phase 3: which of my vertices are shared, with whom ----
    npair = sum(rcnt2) / 2
    allocate(shared_with(nv, 0:nrank-1))
    shared_with = .false.
    do i = 1, npair
       g = rbuf(2*i - 1)
       src = rbuf(2*i)
       k = bsearch_i4(myv, nv, g)
       if (k > 0) shared_with(k, src) = .true.
    end do

    ! ---------------- phase 4: which elements to send, and their sizes ----
    allocate(send_elem(nelv))
    allocate(ecnt(nrank), icnt(nrank), rlen(nrank))
    ecnt = 0; icnt = 0; rlen = 0
    do p = 0, nrank - 1
       if (p == myrank) cycle
       call mark_touching(p, send_elem)
       do e = 1, nelv
          if (.not. send_elem(e)) cycle
          nv_e = vtx_ptr(e + 1) - vtx_ptr(e)
          ecnt(p + 1) = ecnt(p + 1) + 1
          icnt(p + 1) = icnt(p + 1) + elem_int_len(e)
          rlen(p + 1) = rlen(p + 1) + nv_e * nv_e
       end do
    end do
    allocate(recnt(nrank), ircnt(nrank), rrcnt(nrank))
    call MPI_Alltoall(ecnt, 1, MPI_INTEGER, recnt, 1, MPI_INTEGER, comm, ierr)
    call MPI_Alltoall(icnt, 1, MPI_INTEGER, ircnt, 1, MPI_INTEGER, comm, ierr)
    call MPI_Alltoall(rlen, 1, MPI_INTEGER, rrcnt, 1, MPI_INTEGER, comm, ierr)
    gh%n_ghost = sum(recnt)

    allocate(edsp(nrank), idsp(nrank), rdsp2r(nrank))
    allocate(redsp(nrank), irdsp(nrank), rrdsp(nrank))
    call make_displs(ecnt, edsp)
    call make_displs(icnt, idsp)
    call make_displs(rlen, rdsp2r)
    call make_displs(recnt, redsp)
    call make_displs(ircnt, irdsp)
    call make_displs(rrcnt, rrdsp)

    allocate(eints(max(sum(icnt), 1)), ereal(max(sum(rlen), 1)))
    allocate(rints(max(sum(ircnt), 1)), rreal(max(sum(rrcnt), 1)))

    ! ---- pack ----
    do p = 0, nrank - 1
       if (p == myrank) cycle
       if (ecnt(p + 1) == 0) cycle
       call mark_touching(p, send_elem)
       pos_i = 0
       pos_r = 0
       do e = 1, nelv
          if (.not. send_elem(e)) cycle
          nv_e = vtx_ptr(e + 1) - vtx_ptr(e)
          eints(idsp(p+1) + pos_i + 1) = nv_e
          eints(idsp(p+1) + pos_i + 2 : idsp(p+1) + pos_i + 1 + nv_e) = &
               vtx_idx(vtx_ptr(e) + 1 : vtx_ptr(e + 1))
          pos_i = pos_i + 1 + nv_e
          eints(idsp(p+1) + pos_i + 1) = face_ptr(e + 1) - face_ptr(e)
          pos_i = pos_i + 1
          do a = face_ptr(e) + 1, face_ptr(e + 1)
             nvf = face_vtx_ptr(a + 1) - face_vtx_ptr(a)
             eints(idsp(p+1) + pos_i + 1) = nvf
             eints(idsp(p+1) + pos_i + 2 : idsp(p+1) + pos_i + 1 + nvf) = &
                  face_vtx_idx(face_vtx_ptr(a) + 1 : face_vtx_ptr(a + 1))
             pos_i = pos_i + 1 + nvf
             nef = face_edge_ptr(a + 1) - face_edge_ptr(a)
             eints(idsp(p+1) + pos_i + 1) = nef
             pos_i = pos_i + 1
             do k = 1, nef
                eints(idsp(p+1) + pos_i + 1) = face_edge_vtx(1, face_edge_ptr(a) + k)
                eints(idsp(p+1) + pos_i + 2) = face_edge_vtx(2, face_edge_ptr(a) + k)
                pos_i = pos_i + 2
             end do
          end do
          eints(idsp(p+1) + pos_i + 1) = part_local(e) + offset
          pos_i = pos_i + 1
          ereal(rdsp2r(p+1) + pos_r + 1 : rdsp2r(p+1) + pos_r + nv_e * nv_e) = &
               amat(amat_ptr(e) + 1 : amat_ptr(e + 1))
          pos_r = pos_r + nv_e * nv_e
       end do
    end do

    call MPI_Alltoallv(eints, icnt, idsp, MPI_INTEGER, &
                       rints, ircnt, irdsp, MPI_INTEGER, comm, ierr)
    call MPI_Alltoallv(ereal, rlen, rdsp2r, MPI_DOUBLE_PRECISION, &
                       rreal, rrcnt, rrdsp, MPI_DOUBLE_PRECISION, comm, ierr)

    ! ---- parse: pass 1 (sizes), pass 2 (fill) ----
    n_gvref_tot = 0; n_gf_tot = 0; n_gfvref_tot = 0; n_geref_tot = 0
    pi = 0
    do idx = 1, gh%n_ghost
       nv_e = rints(pi + 1)
       n_gvref_tot = n_gvref_tot + nv_e
       pi = pi + 1 + nv_e
       gf = rints(pi + 1)
       n_gf_tot = n_gf_tot + gf
       pi = pi + 1
       do a = 1, gf
          nvf = rints(pi + 1)
          n_gfvref_tot = n_gfvref_tot + nvf
          pi = pi + 1 + nvf
          nef = rints(pi + 1)
          n_geref_tot = n_geref_tot + nef
          pi = pi + 1 + 2 * nef
       end do
       pi = pi + 1  ! gmacro
    end do

    allocate(gh%vtx_ptr(gh%n_ghost + 1), gh%face_ptr(gh%n_ghost + 1))
    allocate(gh%amat_ptr(gh%n_ghost + 1), gh%gmacro(max(gh%n_ghost, 1)))
    allocate(gh%vtx_idx(max(n_gvref_tot, 1)))
    allocate(gh%face_vtx_ptr(n_gf_tot + 1), gh%face_edge_ptr(n_gf_tot + 1))
    allocate(gh%face_vtx_idx(max(n_gfvref_tot, 1)))
    allocate(gh%face_edge_vtx(2, max(n_geref_tot, 1)))
    allocate(gh%amat(max(sum(rrcnt), 1)))
    gh%vtx_ptr(1) = 0; gh%face_ptr(1) = 0; gh%amat_ptr(1) = 0
    gh%face_vtx_ptr(1) = 0; gh%face_edge_ptr(1) = 0

    pi = 0; pr_ = 0; gfv = 0; gfe = 0
    do idx = 1, gh%n_ghost
       nv_e = rints(pi + 1)
       gh%vtx_ptr(idx + 1) = gh%vtx_ptr(idx) + nv_e
       gh%vtx_idx(gh%vtx_ptr(idx) + 1 : gh%vtx_ptr(idx + 1)) = &
            rints(pi + 2 : pi + 1 + nv_e)
       pi = pi + 1 + nv_e
       gf = rints(pi + 1)
       gh%face_ptr(idx + 1) = gh%face_ptr(idx) + gf
       pi = pi + 1
       do a = 1, gf
          gidx = gh%face_ptr(idx) + a
          nvf = rints(pi + 1)
          gh%face_vtx_ptr(gidx + 1) = gh%face_vtx_ptr(gidx) + nvf
          gh%face_vtx_idx(gh%face_vtx_ptr(gidx) + 1 : gh%face_vtx_ptr(gidx + 1)) = &
               rints(pi + 2 : pi + 1 + nvf)
          pi = pi + 1 + nvf
          nef = rints(pi + 1)
          gh%face_edge_ptr(gidx + 1) = gh%face_edge_ptr(gidx) + nef
          pi = pi + 1
          do k = 1, nef
             gh%face_edge_vtx(1, gh%face_edge_ptr(gidx) + k) = rints(pi + 1)
             gh%face_edge_vtx(2, gh%face_edge_ptr(gidx) + k) = rints(pi + 2)
             pi = pi + 2
          end do
       end do
       gh%gmacro(idx) = rints(pi + 1)
       pi = pi + 1
       nv_e = gh%vtx_ptr(idx + 1) - gh%vtx_ptr(idx)
       gh%amat_ptr(idx + 1) = gh%amat_ptr(idx) + nv_e * nv_e
       gh%amat(gh%amat_ptr(idx) + 1 : gh%amat_ptr(idx + 1)) = &
            rreal(pr_ + 1 : pr_ + nv_e * nv_e)
       pr_ = pr_ + nv_e * nv_e
    end do

  contains

    !> Mark which of my local elements touch a vertex shared with rank p.
    subroutine mark_touching(p, mark)
      integer(i4), intent(in) :: p
      logical, intent(inout) :: mark(:)
      integer(i4) :: ee, kk, ii
      mark = .false.
      do ee = 1, nelv
         do kk = vtx_ptr(ee) + 1, vtx_ptr(ee + 1)
            ii = bsearch_i4(myv, nv, vtx_idx(kk))
            if (ii > 0) then
               if (shared_with(ii, p)) then
                  mark(ee) = .true.
                  exit
               end if
            end if
         end do
      end do
    end subroutine mark_touching

    !> Length (in ints) of element e's self-describing payload record.
    function elem_int_len(e) result(len_)
      integer(i4), intent(in) :: e
      integer(i4) :: len_
      integer(i4) :: aa
      len_ = 2 + (vtx_ptr(e + 1) - vtx_ptr(e)) + 1  ! nv_e, verts, nf_e, gmacro
      do aa = face_ptr(e) + 1, face_ptr(e + 1)
         len_ = len_ + 2 + (face_vtx_ptr(aa + 1) - face_vtx_ptr(aa)) + &
              2 * (face_edge_ptr(aa + 1) - face_edge_ptr(aa))
      end do
    end function elem_int_len

  end subroutine amge_ghost_exchange

  !> Build the extended macro-id list (part_ext) that
  !! topo%init_tables(mmsh_ext, part_ext, gh%n_macro_global) consumes:
  !! this rank's own macro ids (offset to be globally unique) followed by
  !! the ghost elements' already-global macro ids, in the SAME order
  !! macro_mesh_splice_ghost appends ghost elements
  !! (mmsh%n_elem+1 .. mmsh%n_elem+n_ghost).
  subroutine amge_ghost_part_ext(gh, part_local, part_ext)
    type(amge_ghost_t), intent(in) :: gh
    integer(i4), intent(in) :: part_local(:)
    integer(i4), allocatable, intent(out) :: part_ext(:)
    integer(i4) :: n
    n = gh%n_local + gh%n_ghost
    allocate(part_ext(n))
    part_ext(1:gh%n_local) = part_local + gh%macro_offset
    if (gh%n_ghost > 0) part_ext(gh%n_local+1:n) = gh%gmacro(1:gh%n_ghost)
  end subroutine amge_ghost_part_ext

  !> Mask of macro ids this rank OWNS (must emit macroelements for).
  !! Ghost macros appear only to supply labels and trace-map
  !! neighborhoods; their macroelements belong to another rank.
  subroutine amge_ghost_owned_mask(gh, n_macro_local, owned)
    type(amge_ghost_t), intent(in) :: gh
    integer(i4), intent(in) :: n_macro_local
    logical, allocatable, intent(out) :: owned(:)
    integer(i4) :: m
    allocate(owned(gh%n_macro_global))
    owned = .false.
    do m = 1, n_macro_local
       owned(gh%macro_offset + m) = .true.
    end do
  end subroutine amge_ghost_owned_mask

  subroutine amge_ghost_free(this)
    class(amge_ghost_t), intent(inout) :: this
    if (allocated(this%vtx_ptr)) deallocate(this%vtx_ptr)
    if (allocated(this%vtx_idx)) deallocate(this%vtx_idx)
    if (allocated(this%face_ptr)) deallocate(this%face_ptr)
    if (allocated(this%face_vtx_ptr)) deallocate(this%face_vtx_ptr)
    if (allocated(this%face_vtx_idx)) deallocate(this%face_vtx_idx)
    if (allocated(this%face_edge_ptr)) deallocate(this%face_edge_ptr)
    if (allocated(this%face_edge_vtx)) deallocate(this%face_edge_vtx)
    if (allocated(this%gmacro)) deallocate(this%gmacro)
    if (allocated(this%amat_ptr)) deallocate(this%amat_ptr)
    if (allocated(this%amat)) deallocate(this%amat)
    this%n_ghost = 0; this%n_local = 0
  end subroutine amge_ghost_free

  ! ---------------- small helpers ----------------

  subroutine make_displs(cnt, dsp)
    integer, intent(in) :: cnt(:)
    integer, intent(out) :: dsp(:)
    integer :: i
    dsp(1) = 0
    do i = 2, size(cnt)
       dsp(i) = dsp(i-1) + cnt(i-1)
    end do
  end subroutine make_displs

  pure subroutine sort_i4(a)
    integer(i4), intent(inout) :: a(:)
    call qsort_i4(a, 1, size(a))
  end subroutine sort_i4

  pure recursive subroutine qsort_i4(a, lo, hi)
    integer(i4), intent(inout) :: a(:)
    integer(i4), intent(in) :: lo, hi
    integer(i4) :: i, j, p, t
    if (lo >= hi) return
    p = a((lo + hi) / 2)
    i = lo; j = hi
    do
       ! nested guards: Fortran does not short-circuit .and.
       do
          if (i > hi) exit
          if (a(i) >= p) exit
          i = i + 1
       end do
       do
          if (j < lo) exit
          if (a(j) <= p) exit
          j = j - 1
       end do
       if (i > j) exit
       t = a(i); a(i) = a(j); a(j) = t
       i = i + 1; j = j - 1
    end do
    call qsort_i4(a, lo, j)
    call qsort_i4(a, i, hi)
  end subroutine qsort_i4

  !> ord(1:n) = permutation sorting key(1:n) ascending
  subroutine sort_index_i4(key, n, ord)
    integer(i4), intent(in) :: key(:)
    integer(i4), intent(in) :: n
    integer(i4), intent(out) :: ord(:)
    integer(i4) :: i
    do i = 1, n
       ord(i) = i
    end do
    if (n > 1) call qsort_ord(key, ord, 1, n)
  end subroutine sort_index_i4

  recursive subroutine qsort_ord(key, ord, lo, hi)
    integer(i4), intent(in) :: key(:)
    integer(i4), intent(inout) :: ord(:)
    integer(i4), intent(in) :: lo, hi
    integer(i4) :: i, j, p, t
    if (lo >= hi) return
    p = key(ord((lo + hi) / 2))
    i = lo; j = hi
    do
       do
          if (i > hi) exit
          if (key(ord(i)) >= p) exit
          i = i + 1
       end do
       do
          if (j < lo) exit
          if (key(ord(j)) <= p) exit
          j = j - 1
       end do
       if (i > j) exit
       t = ord(i); ord(i) = ord(j); ord(j) = t
       i = i + 1; j = j - 1
    end do
    call qsort_ord(key, ord, lo, j)
    call qsort_ord(key, ord, i, hi)
  end subroutine qsort_ord

  pure function bsearch_i4(a, n, v) result(k)
    integer(i4), intent(in) :: a(:), n, v
    integer(i4) :: k, lo, hi, mid
    lo = 1; hi = n; k = 0
    do while (lo <= hi)
       mid = (lo + hi) / 2
       if (a(mid) == v) then
          k = mid
          return
       else if (a(mid) < v) then
          lo = mid + 1
       else
          hi = mid - 1
       end if
    end do
  end function bsearch_i4

end module amge_ghost
