! Copyright (c) 2026, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the conditions of the
! BSD 3-Clause License are met (see Neko's COPYING file).
!
!> One-layer ghost (halo) exchange that makes the rank-local AMGe
!! macroentity extraction and trace maps agree with a serial run.
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
!! HOW THIS MODULE IS USED. It returns an EXTENDED element list
!! (local elements followed by ghosts) with globally unique macroelement
!! labels. Feed that list to the ordinary serial path:
!!   * macro_mesh_init_hex sees two incident elements on interface facets
!!     and labels them (m_mine, m_remote) -- issue (1) disappears;
!!   * the trace-map neighborhoods pick up the ghost elements and their
!!     local matrices -- issue (2) disappears.
!! The only changes needed downstream are: agglomerate only the LOCAL
!! elements, and emit macroelements only for locally OWNED macro ids
!! (see amge_ghost_owned_mask).
!!
!! Communication is a standard rendezvous: vertex ids are hashed to a
!! distributed directory which reports co-owners, then each pair of
!! neighbouring ranks exchanges the elements touching their shared
!! vertices, together with those elements' local matrices.
module amge_ghost
  use mpi
  use num_types, only : i4, rp
  use utils, only : neko_error
  implicit none
  private

  !> Ghost element data, to be appended after the local elements.
  type, public :: amge_ghost_t
     integer(i4) :: n_ghost = 0
     integer(i4) :: n_local = 0
     integer(i4) :: n_macro_global = 0            !< total macros over all ranks
     integer(i4) :: macro_offset = 0              !< my local id 1 -> global id offset+1
     integer(i4), allocatable :: vtx(:,:)         !< (8, n_ghost) global vertex ids
     integer(i4), allocatable :: gmacro(:)        !< (n_ghost) global macro id
     real(rp), allocatable :: amat(:,:)           !< (64, n_ghost) col-major 8x8
   contains
     procedure, pass(this) :: free => amge_ghost_free
  end type amge_ghost_t

  public :: amge_ghost_exchange, amge_ghost_extended_list, amge_ghost_owned_mask

contains

  !> Build the ghost layer. All ids in/out are GLOBAL mesh vertex ids.
  !! @param comm       MPI communicator (Neko: NEKO_COMM)
  !! @param hv_local   (8, nelv) global vertex ids of my elements
  !! @param amat_local (64, nelv) my element matrices, column-major 8x8
  !! @param part_local (nelv) my LOCAL macro ids, 1..n_macro_local
  !! @param n_macro_local number of macroelements I own
  !! @param gh         result
  subroutine amge_ghost_exchange(comm, hv_local, amat_local, part_local, &
                                 n_macro_local, gh)
    integer, intent(in) :: comm
    integer(i4), intent(in) :: hv_local(:,:)
    real(rp), intent(in) :: amat_local(:,:)
    integer(i4), intent(in) :: part_local(:), n_macro_local
    type(amge_ghost_t), intent(inout) :: gh

    integer :: ierr, myrank, nrank, p
    integer(i4) :: nelv, i, j, k, e, nv, cnt, g, src, ndir, npair
    integer(i4), allocatable :: myv(:)                     ! my unique vertex gids
    integer, allocatable :: scnt(:), sdsp(:), rcnt(:), rdsp(:)
    integer, allocatable :: scnt2(:), sdsp2(:), rcnt2(:), rdsp2(:)
    integer(i4), allocatable :: sbuf(:), rbuf(:)
    integer(i4), allocatable :: dgid(:), dsrc(:), ord(:)
    integer(i4), allocatable :: pgid(:), prnk(:)           ! (gid, other rank) pairs
    logical, allocatable :: shared_with(:,:)               ! (nv, 0:nrank-1)
    logical, allocatable :: send_elem(:)
    integer(i4), allocatable :: eints(:), rints(:)
    real(rp), allocatable :: ereal(:), rreal(:)
    integer(i4) :: offset

    call MPI_Comm_rank(comm, myrank, ierr)
    call MPI_Comm_size(comm, nrank, ierr)
    call gh%free()
    nelv = size(hv_local, 2)
    gh%n_local = nelv

    ! ---------------- globally unique macro ids ----------------
    ! exclusive prefix sum: my local macro m becomes  m + offset
    offset = 0
    call MPI_Exscan(n_macro_local, offset, 1, MPI_INTEGER, MPI_SUM, comm, ierr)
    if (myrank == 0) offset = 0
    gh%macro_offset = offset
    call MPI_Allreduce(n_macro_local, gh%n_macro_global, 1, MPI_INTEGER, &
         MPI_SUM, comm, ierr)

    ! ---------------- my unique vertex gids ----------------
    allocate(myv(8 * nelv))
    cnt = 0
    do e = 1, nelv
       do k = 1, 8
          cnt = cnt + 1
          myv(cnt) = hv_local(k, e)
       end do
    end do
    call sort_i4(myv)
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
    ! for each gid group, every member learns about every other member
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
      integer(i4) :: a, b
      allocate(fill(nrank)); fill = 0
      i = 1
      do while (i <= ndir)
         j = i
         do while (j < ndir)
            if (dgid(ord(j+1)) /= dgid(ord(i))) exit
            j = j + 1
         end do
         do a = i, j
            do b = i, j
               if (a == b) cycle
               p = dsrc(ord(a))
               fill(p + 1) = fill(p + 1) + 1
               sbuf(sdsp2(p + 1) + fill(p + 1)) = dgid(ord(a))
               fill(p + 1) = fill(p + 1) + 1
               sbuf(sdsp2(p + 1) + fill(p + 1)) = dsrc(ord(b))
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

    ! ---------------- phase 4: exchange elements touching shared vertices ----
    ! payload per element: 9 ints (8 vertex gids + global macro id) + 64 reals
    allocate(send_elem(nelv))
    scnt = 0
    do p = 0, nrank - 1
       if (p == myrank) cycle
       send_elem = .false.
       do e = 1, nelv
          do k = 1, 8
             i = bsearch_i4(myv, nv, hv_local(k, e))
             if (i > 0) then
                if (shared_with(i, p)) then
                   send_elem(e) = .true.
                   exit
                end if
             end if
          end do
       end do
       scnt(p + 1) = count(send_elem)
    end do
    call MPI_Alltoall(scnt, 1, MPI_INTEGER, rcnt, 1, MPI_INTEGER, comm, ierr)
    call make_displs(scnt, sdsp)
    call make_displs(rcnt, rdsp)
    gh%n_ghost = sum(rcnt)

    allocate(eints(max(9 * sum(scnt), 1)), ereal(max(64 * sum(scnt), 1)))
    allocate(rints(max(9 * gh%n_ghost, 1)), rreal(max(64 * gh%n_ghost, 1)))
    block
      integer(i4) :: pos
      do p = 0, nrank - 1
         if (scnt(p + 1) == 0) cycle
         send_elem = .false.
         do e = 1, nelv
            do k = 1, 8
               i = bsearch_i4(myv, nv, hv_local(k, e))
               if (i > 0) then
                  if (shared_with(i, p)) then
                     send_elem(e) = .true.
                     exit
                  end if
               end if
            end do
         end do
         pos = 0
         do e = 1, nelv
            if (.not. send_elem(e)) cycle
            eints(9 * sdsp(p + 1) + 9 * pos + 1 : 9 * sdsp(p + 1) + 9 * pos + 8) = &
                 hv_local(:, e)
            eints(9 * sdsp(p + 1) + 9 * pos + 9) = part_local(e) + offset
            ereal(64 * sdsp(p + 1) + 64 * pos + 1 : 64 * sdsp(p + 1) + 64 * pos + 64) = &
                 amat_local(:, e)
            pos = pos + 1
         end do
      end do
    end block
    block
      integer, allocatable :: sc9(:), sd9(:), rc9(:), rd9(:)
      integer, allocatable :: sc64(:), sd64(:), rc64(:), rd64(:)
      allocate(sc9(nrank), sd9(nrank), rc9(nrank), rd9(nrank))
      allocate(sc64(nrank), sd64(nrank), rc64(nrank), rd64(nrank))
      sc9 = 9 * scnt;  rc9 = 9 * rcnt
      sc64 = 64 * scnt; rc64 = 64 * rcnt
      sd9 = 9 * sdsp;  rd9 = 9 * rdsp
      sd64 = 64 * sdsp; rd64 = 64 * rdsp
      call MPI_Alltoallv(eints, sc9, sd9, MPI_INTEGER, &
                         rints, rc9, rd9, MPI_INTEGER, comm, ierr)
      call MPI_Alltoallv(ereal, sc64, sd64, MPI_DOUBLE_PRECISION, &
                         rreal, rc64, rd64, MPI_DOUBLE_PRECISION, comm, ierr)
    end block

    allocate(gh%vtx(8, max(gh%n_ghost,1)), gh%gmacro(max(gh%n_ghost,1)))
    allocate(gh%amat(64, max(gh%n_ghost,1)))
    do e = 1, gh%n_ghost
       gh%vtx(:, e) = rints(9*(e-1) + 1 : 9*(e-1) + 8)
       gh%gmacro(e) = rints(9*(e-1) + 9)
       gh%amat(:, e) = rreal(64*(e-1) + 1 : 64*(e-1) + 64)
    end do
  end subroutine amge_ghost_exchange

  !> Concatenate local + ghost into the extended element list that the
  !! ordinary (serial) extraction and coarsening consume.
  subroutine amge_ghost_extended_list(gh, hv_local, amat_local, part_local, &
                                      hv_ext, amat_ext, part_ext)
    type(amge_ghost_t), intent(in) :: gh
    integer(i4), intent(in) :: hv_local(:,:), part_local(:)
    real(rp), intent(in) :: amat_local(:,:)
    integer(i4), allocatable, intent(out) :: hv_ext(:,:), part_ext(:)
    real(rp), allocatable, intent(out) :: amat_ext(:,:)
    integer(i4) :: n
    n = gh%n_local + gh%n_ghost
    allocate(hv_ext(8, n), amat_ext(64, n), part_ext(n))
    hv_ext(:, 1:gh%n_local) = hv_local
    amat_ext(:, 1:gh%n_local) = amat_local
    part_ext(1:gh%n_local) = part_local + gh%macro_offset
    if (gh%n_ghost > 0) then
       hv_ext(:, gh%n_local+1:n) = gh%vtx(:, 1:gh%n_ghost)
       amat_ext(:, gh%n_local+1:n) = gh%amat(:, 1:gh%n_ghost)
       part_ext(gh%n_local+1:n) = gh%gmacro(1:gh%n_ghost)
    end if
  end subroutine amge_ghost_extended_list

  !> Mask of macro ids this rank OWNS (must emit macroelements for).
  !! Ghost macros appear in the extended list only to supply labels and
  !! trace-map neighborhoods; their macroelements belong to another rank.
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
    if (allocated(this%vtx)) deallocate(this%vtx)
    if (allocated(this%gmacro)) deallocate(this%gmacro)
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
