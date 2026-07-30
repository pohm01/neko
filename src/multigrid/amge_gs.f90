! Copyright (c) 2026, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the conditions of the
! BSD 3-Clause License are met (see Neko's COPYING file).
!
!> Gather-scatter for AMGe levels, without ever building a dofmap.
!!
!! WHY THIS WORKS. v := G G^T v -- sum every copy of a dof, across
!! elements and MPI ranks, and broadcast the sum back -- only needs, per
!! storage slot, (a) the GLOBAL id of the dof stored there and (b) whether
!! that dof is shared with another rank. Neko's gs_t gets both by walking a
!! dofmap_t as a tensor-product (lx,ly,lz,nelv) grid, which coarse AMGe
!! levels do not have (a macroelement is a polyhedron with a variable
!! number of vertex dofs). AMGe levels DO have both directly, though:
!!      (a) mmsh%vert_id  -- global ids, valid on every level because a
!!          macrovertex IS a fine vertex and build_next_level inherits the
!!          id verbatim from msh%elements(e)%e%pts(k)%p%id();
!!      (b) mmsh%shared_vtx -- set on the finest (Q1) level from
!!          dofmap%shared_dof (amge_mesh_set_shared_from_dofmap) and
!!          inherited the same way.
!!
!! THE DESIGN. amge_gs_map_t builds the local/shared slot <-> global-id
!! mapping directly from the level's CSR element-dof lists (see its own
!! doc comment below). amge_gs_map_init_comm then wires the shared group
!! into a real gs_t's MPI schedule by populating gs%shared_dofs and
!! calling gs_schedule -- which never reads gs%dofmap or gs%Xh, only
!! gs%shared_dofs and gs%comm, so it is reusable verbatim. amge_gs_t is
!! the per-level handle a caller actually owns: it wraps one
!! amge_gs_map_t plus one gs_t (used only for its comm), and its op(u)
!! reproduces v := G G^T v with a true cross-rank reduction for shared
!! dofs, no packing/padding required since the map already operates
!! directly on the level's native (elm_vtx_ptr, elm_vtx_idx) layout.
!!
!! FUTURE OPTIMIZATION (not implemented): level 0's vector layout
!! (8 corner dofs, nelv) in Neko hex_t vertex order is byte-identical to a
!! Q1 (lx=ly=lz=2) Neko field, so level 0 could instead reuse the
!! coef/dofmap/gs_h the solver already owns when a Q1 one exists, rather
!! than building its own map -- a zero-copy shortcut orthogonal to the
!! map-based design used here for every level (including 0).
!!
!! MPI CAVEATS THAT REMAIN:
!!  * macro_topology is still rank-local: a facet with one local element
!!    is labelled (part,0), conflating the physical boundary with a rank
!!    interface. Consult msh%facet_neigh before trusting label 0, or the
!!    extracted macroentities differ across ranks.
!!  * agglomerates must not cross ranks unless the agglomeration is made
!!    global; the gs above then correctly stitches the coarse dofs that
!!    lie ON rank interfaces, which is the common case.
!!  * the coarsest level should eventually be gathered to one rank (or
!!    solved redundantly) rather than gs'd, once it is small.
module amge_gs
  use num_types, only : i4, i8, rp
  use dofmap, only : dofmap_t
  use gather_scatter, only : gs_t, GS_OP_ADD, gs_schedule, gs_comm_alloc_init, &
       GS_COMM_MPI, GS_COMM_MPIGPU
  use neko_config, only : NEKO_DEVICE_MPI
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR
  implicit none
  private

  !> Gather-scatter mapping for one AMGe level.
  type, public :: amge_gs_map_t
     integer(i4) :: ndofs = 0        !< storage slots on this level
     integer(i4) :: nlocal = 0       !< slots in the rank-local group
     integer(i4) :: nshared = 0      !< slots in the shared group
     integer(i4) :: nlocal_blks = 0  !< distinct local gather targets
     integer(i4) :: nshared_blks = 0 !< distinct shared gather targets
     integer(i4), allocatable :: local_dof_gs(:)   !< slot index
     integer(i4), allocatable :: local_gs_dof(:)   !< gather-buffer index
     integer(i4), allocatable :: local_blk_len(:)  !< multiplicities
     integer(i4), allocatable :: local_blk_off(:)  !< offset of each block
     integer(i4), allocatable :: shared_dof_gs(:)
     integer(i4), allocatable :: shared_gs_dof(:)
     integer(i4), allocatable :: shared_blk_len(:)
     integer(i4), allocatable :: shared_blk_off(:) !< offset of each block
     integer(i4), allocatable :: shared_gid(:)     !< global id per shared slot
   contains
     procedure, pass(this) :: init => amge_gs_map_init
     procedure, pass(this) :: init_comm => amge_gs_map_init_comm
     procedure, pass(this) :: apply => amge_gs_map_apply
     procedure, pass(this) :: free => amge_gs_map_free
  end type amge_gs_map_t

  !> Per-level gather-scatter handle: a local/shared mapping plus the gs_t
  !! that owns its MPI comm schedule.
  type, public :: amge_gs_t
     type(amge_gs_map_t) :: map !< local/shared slot <-> global id mapping
     type(gs_t) :: gs_h         !< owns only the MPI comm (gs_h%comm)
   contains
     procedure, pass(this) :: init => amge_gs_init
     procedure, pass(this) :: op => amge_gs_op
     procedure, pass(this) :: free => amge_gs_free
     procedure, pass(this) :: correct_shared_count => amge_gs_correct_shared_count
  end type amge_gs_t

  public :: amge_gs_map_shared_ids, amge_mesh_set_shared_from_dofmap

contains

  !> Build the gather-scatter for one AMGe level: the local/shared mapping
  !! plus (if any dof is shared) the MPI schedule to exchange it.
  !! Takes the level's raw index arrays rather than the level type, so
  !! this module sits BELOW amge_level and a level can own a handle.
  !! @param nelm         elements on this level
  !! @param elm_vtx_ptr  (nelm+1) CSR offsets into elm_vtx_idx
  !! @param elm_vtx_idx  (ndofs) local vertex index of each storage slot
  !! @param vert_id      (n_verts) local vertex -> GLOBAL dof id
  !! @param shared_vtx   (n_verts) local vertex -> shared with another rank?
  !! @param comm_bcknd   optional comm backend override
  subroutine amge_gs_init(this, nelm, elm_vtx_ptr, elm_vtx_idx, vert_id, &
                          shared_vtx, comm_bcknd)
    class(amge_gs_t), intent(inout) :: this
    integer(i4), intent(in) :: nelm
    integer(i4), intent(in) :: elm_vtx_ptr(:), elm_vtx_idx(:), vert_id(:)
    logical, intent(in) :: shared_vtx(:)
    integer, optional, intent(in) :: comm_bcknd

    call this%map%init(nelm, elm_vtx_ptr, elm_vtx_idx, vert_id, shared_vtx)
    call this%map%init_comm(this%gs_h, comm_bcknd)
  end subroutine amge_gs_init

  !> v := G G^T v. Replaces amge_gs_placeholder verbatim at every call
  !! site: call lvl%gsh%op(v%x) (sites: amge_smooth_l1 x2, calc_resid, and
  !! the restriction assembly in amge_flat_vcycle). No elm_vtx_ptr/n_dofs
  !! needed -- the map already operates directly on the level's native
  !! (elm_vtx_ptr, elm_vtx_idx) layout, so there is nothing to pack/unpack.
  !!
  !! Once wired up, delete amge_gs_placeholder, amge_gather and
  !! amge_scatter_add (they survive only as the placeholder's internals;
  !! mult/valence is unaffected and can also be computed with one gs on a
  !! vector of ones: v = 1 ; gs ; mult = 1/v).
  subroutine amge_gs_op(this, u)
    class(amge_gs_t), intent(inout) :: this
    real(rp), intent(inout) :: u(:)
    call this%map%apply(u, this%gs_h)
  end subroutine amge_gs_op

  subroutine amge_gs_free(this)
    class(amge_gs_t), intent(inout) :: this
    call this%map%free()
    call this%gs_h%free()
  end subroutine amge_gs_free

  !> Correct a per-UNIQUE-vertex scalar (e.g. a locally-computed
  !! duplication/membership count) to its TRUE cross-rank total, at
  !! exactly the vertices this level's gs handle already knows are shared
  !! (this%map%shared_gid) -- values at every other vertex are already
  !! exact from local information alone and are left untouched.
  !!
  !! This is a DIFFERENT granularity than op()/apply(): op() sums a
  !! per-DUPLICATED-SLOT field (one entry per (element, corner) copy),
  !! which is the wrong thing to reach for here -- broadcasting a
  !! per-vertex count out to every local duplicate slot and summing via
  !! op() would multiply in the local slot count, not just add the
  !! cross-rank contribution. Here w is already exactly one value per
  !! unique vertex, so only the compact shared-block exchange (the same
  !! nbsend/nbrecv/nbwait primitives apply() uses for its shared group) is
  !! needed, with no local gather/scatter step at all.
  !! @param elm_vtx_idx  the level's own CSR slot -> unique-vertex map
  !!                      (used to translate this%map's slot-indexed
  !!                      shared_dof_gs into the vertex-indexed w)
  !! @param w  length = number of unique vertices on this level; corrected
  !!           in place at shared vertices
  subroutine amge_gs_correct_shared_count(this, elm_vtx_idx, w)
    class(amge_gs_t), intent(inout) :: this
    integer(i4), intent(in) :: elm_vtx_idx(:)
    real(rp), intent(inout) :: w(:)
    integer(i4), allocatable :: vert2blk(:)
    real(rp), allocatable :: sbuf(:)
    integer(i4) :: k, v, tid
    type(c_ptr) :: deps, strm

    if (this%map%nshared_blks .eq. 0) return

    allocate(vert2blk(size(w)))
    vert2blk = 0
    do k = 1, this%map%nshared
       v = elm_vtx_idx(this%map%shared_dof_gs(k))
       vert2blk(v) = this%map%shared_gs_dof(k)
    end do

    allocate(sbuf(this%map%nshared_blks))
    sbuf = 0.0_rp
    do v = 1, size(w)
       if (vert2blk(v) .gt. 0) sbuf(vert2blk(v)) = w(v)
    end do

    tid = 0
    deps = C_NULL_PTR
    strm = C_NULL_PTR
    call this%gs_h%comm%nbsend(sbuf, this%map%nshared_blks, tid, deps, strm)
    call this%gs_h%comm%nbrecv(tid)
    call this%gs_h%comm%nbwait(sbuf, this%map%nshared_blks, GS_OP_ADD, strm)

    do v = 1, size(w)
       if (vert2blk(v) .gt. 0) w(v) = sbuf(vert2blk(v))
    end do
  end subroutine amge_gs_correct_shared_count

  !> Gather-scatter MAPPING for AMGe levels -- the analogue of Neko's
  !! gs_init_mapping, built from macro_mesh_t instead of dofmap_t.
  !!
  !! WHY A SEPARATE ROUTINE IS NEEDED. Neko's gs_init_mapping walks a
  !! dofmap_t as a (lx,ly,lz,nelv) tensor grid and uses that geometric
  !! structure to enumerate dofs and to split them into interior / facet
  !! groups for the blocked kernels. AMGe coarse levels have no such
  !! structure: a macroelement is a polyhedron with a variable number of
  !! vertex dofs, so lx/ly/lz do not exist. What gs_init_mapping actually
  !! REQUIRES, however, is only two things per storage slot:
  !!      (a) the GLOBAL id of the dof stored there,
  !!      (b) whether that dof is shared with another MPI rank,
  !! everything else being enumeration order and performance grouping.
  !! AMGe has both:
  !!      (a) mmsh%vert_id  -- global ids, valid on every level because a
  !!          macrovertex IS a fine vertex and build_next_level inherits
  !!          the id verbatim;
  !!      (b) mmsh%shared_vtx -- set on the finest (Q1) level from
  !!          dofmap%shared_dof and inherited the same way.
  !! So the mapping can be built for every level from the CSR element-dof
  !! lists alone.
  !!
  !! WHAT THIS PRODUCES. The same three-array pattern Neko's gs uses, in
  !! two groups (rank-local and shared):
  !!      dof_gs(k)  : index into the level vector (the storage slot)
  !!      gs_dof(k)  : index into the gather buffer (one entry per distinct
  !!                   global dof in that group)
  !!      blk_len(b) : run length of consecutive k with equal gs_dof, i.e.
  !!                   the multiplicity of gather target b -- entries are
  !!                   ordered so each target's slots are contiguous.
  !! With these,
  !!      gather :  buf = 0 ;  buf(gs_dof(k)) += u(dof_gs(k))
  !!      scatter:  u(dof_gs(k)) = buf(gs_dof(k))
  !! reproduces v := G G^T v exactly (see amge_gs_map_apply, and the test).
  !!
  !! Slots whose dof has multiplicity 1 AND is not shared are omitted
  !! entirely: summing a single copy is the identity, so they need neither
  !! gather nor communication. This matches Neko's behaviour.
  !!
  !! DISTRIBUTED USE. The shared group is the part that needs MPI: after
  !! the local gather into shared_buf, the entries must be summed across
  !! ranks by global id before scattering back. amge_gs_map_shared_ids
  !! returns the global id of each shared gather slot, in gather-buffer
  !! order, which is exactly what is needed to initialise a gs_comm_t (or
  !! the rendezvous exchange in amge_ghost).

  !> Build the mapping for one level. Inputs are exactly the AMGe level
  !! data: the CSR element-dof lists plus the per-vertex global id and
  !! shared flag carried by macro_mesh_t.
  !! @param nelm         elements on this level
  !! @param elm_vtx_ptr  (nelm+1) CSR offsets into elm_vtx_idx
  !! @param elm_vtx_idx  (ndofs) local vertex index of each storage slot
  !! @param vert_id      (n_verts) local vertex -> GLOBAL dof id
  !! @param shared_vtx   (n_verts) local vertex -> shared with another rank?
  subroutine amge_gs_map_init(this, nelm, elm_vtx_ptr, elm_vtx_idx, &
                              vert_id, shared_vtx)
    class(amge_gs_map_t), intent(inout) :: this
    integer(i4), intent(in) :: nelm
    integer(i4), intent(in) :: elm_vtx_ptr(:), elm_vtx_idx(:), vert_id(:)
    logical, intent(in) :: shared_vtx(:)
    integer(i4), allocatable :: gid(:), ord(:)
    logical, allocatable :: shr(:)
    integer(i4) :: n, s, i, j, run, nl, ns, nlb, nsb, tgt

    call this%free()
    n = elm_vtx_ptr(nelm + 1)
    this%ndofs = n
    if (n == 0) return

    ! slot -> (global id, shared flag)
    allocate(gid(n), shr(n), ord(n))
    do s = 1, n
       gid(s) = vert_id(elm_vtx_idx(s))
       shr(s) = shared_vtx(elm_vtx_idx(s))
    end do

    ! order slots by global id so that all copies of a dof are contiguous.
    ! (Neko gets this contiguity from the tensor-grid walk; here we sort.
    !  Sorting by GLOBAL id also makes the mapping independent of local
    !  element numbering, hence reproducible across rank counts.)
    call sort_index(gid, n, ord)

    ! ---- pass 1: size the two groups ----
    nl = 0; ns = 0; nlb = 0; nsb = 0
    i = 1
    do while (i <= n)
       j = i
       do while (j < n)
          if (gid(ord(j + 1)) /= gid(ord(i))) exit
          j = j + 1
       end do
       run = j - i + 1
       if (shr(ord(i))) then
          ns = ns + run
          nsb = nsb + 1
       else if (run > 1) then
          nl = nl + run
          nlb = nlb + 1
       end if
       ! run == 1 and not shared: identity under gather-scatter, omitted
       i = j + 1
    end do

    this%nlocal = nl; this%nshared = ns
    this%nlocal_blks = nlb; this%nshared_blks = nsb
    allocate(this%local_dof_gs(max(nl,1)), this%local_gs_dof(max(nl,1)))
    allocate(this%local_blk_len(max(nlb,1)))
    allocate(this%shared_dof_gs(max(ns,1)), this%shared_gs_dof(max(ns,1)))
    allocate(this%shared_blk_len(max(nsb,1)), this%shared_gid(max(nsb,1)))

    ! ---- pass 2: fill ----
    nl = 0; ns = 0; nlb = 0; nsb = 0
    i = 1
    do while (i <= n)
       j = i
       do while (j < n)
          if (gid(ord(j + 1)) /= gid(ord(i))) exit
          j = j + 1
       end do
       run = j - i + 1
       if (shr(ord(i))) then
          nsb = nsb + 1
          tgt = nsb
          this%shared_blk_len(nsb) = run
          this%shared_gid(nsb) = gid(ord(i))
          do s = i, j
             ns = ns + 1
             this%shared_dof_gs(ns) = ord(s)
             this%shared_gs_dof(ns) = tgt
          end do
       else if (run > 1) then
          nlb = nlb + 1
          tgt = nlb
          this%local_blk_len(nlb) = run
          do s = i, j
             nl = nl + 1
             this%local_dof_gs(nl) = ord(s)
             this%local_gs_dof(nl) = tgt
          end do
       end if
       i = j + 1
    end do

    ! block offsets: 0-based, blk_off(b) = sum of blk_len(1:b-1). Same
    ! convention as Neko's gs_find_blks, kept in case a future caller wants
    ! to drive gs_bcknd_t's blocked kernels directly instead of the plain
    ! loops in amge_gs_map_apply.
    allocate(this%local_blk_off(max(nlb,1)))
    if (nlb > 0) then
       this%local_blk_off(1) = 0
       do i = 2, nlb
          this%local_blk_off(i) = this%local_blk_off(i-1) &
               + this%local_blk_len(i-1)
       end do
    end if
    allocate(this%shared_blk_off(max(nsb,1)))
    if (nsb > 0) then
       this%shared_blk_off(1) = 0
       do i = 2, nsb
          this%shared_blk_off(i) = this%shared_blk_off(i-1) &
               + this%shared_blk_len(i-1)
       end do
    end if
  end subroutine amge_gs_map_init

  !> Wire this level's shared dofs into a real gs_t's MPI schedule, without
  !! ever building a dofmap. gs_schedule (Neko's crystal-router owner
  !! election that turns shared global ids into per-peer send/recv lists)
  !! only reads gs%shared_dofs and gs%comm -- never gs%dofmap or gs%Xh --
  !! so it is reusable verbatim once those two are populated from this map.
  !! @param gs         a gs_t owned by the caller; freed here before use, so
  !!                    it may be freshly declared or a stale one being
  !!                    re-initialized.
  !! @param comm_bcknd  optional comm backend override, else the same
  !!                     NEKO_DEVICE_MPI default gs_init uses.
  subroutine amge_gs_map_init_comm(this, gs, comm_bcknd)
    class(amge_gs_map_t), intent(in) :: this
    type(gs_t), intent(inout) :: gs
    integer, optional, intent(in) :: comm_bcknd
    integer :: comm_bcknd_, k, idx, idummy
    integer(i8) :: gid8

    call gs%free()

    if (present(comm_bcknd)) then
       comm_bcknd_ = comm_bcknd
    else if (NEKO_DEVICE_MPI) then
       comm_bcknd_ = GS_COMM_MPIGPU
    else
       comm_bcknd_ = GS_COMM_MPI
    end if

    ! global id -> compact shared-buffer index, exactly what gs_schedule
    ! expects in gs%shared_dofs (built there from dofmap%dof instead).
    call gs%shared_dofs%init(max(this%nshared_blks, 1), idummy)
    do k = 1, this%nshared_blks
       gid8 = int(this%shared_gid(k), i8)
       idx = k
       call gs%shared_dofs%set(gid8, idx)
    end do

    call gs_comm_alloc_init(gs%comm, comm_bcknd_)
    call gs_schedule(gs)
  end subroutine amge_gs_map_init_comm

  !> Gather-scatter using the mapping: u := G G^T u. The local group is
  !! always a same-rank reduction (plain loops below). The shared group is
  !! a true cross-rank reduction when @a gs is present -- built by
  !! init_comm, so gs%comm%send_dof/recv_dof already key on the same
  !! compact 1..nshared_blks ids as shared_gs_dof/shared_gid here -- via
  !! Neko's non-blocking gs_comm_t exchange (nbsend/nbrecv/nbwait, the same
  !! primitives gs_t itself calls). Omitting @a gs (or nshared_blks == 0,
  !! e.g. a serial run) performs the local-sum-only reduction, which is
  !! correct whenever no dof in this group is actually shared cross-rank.
  subroutine amge_gs_map_apply(this, u, gs)
    class(amge_gs_map_t), intent(in) :: this
    real(rp), intent(inout) :: u(:)
    type(gs_t), optional, intent(inout) :: gs
    real(rp), allocatable :: lbuf(:), sbuf(:)
    integer(i4) :: k
    integer :: tid
    type(c_ptr) :: deps, strm

    if (this%nlocal_blks > 0) then
       allocate(lbuf(this%nlocal_blks))
       lbuf = 0.0_rp
       do k = 1, this%nlocal
          lbuf(this%local_gs_dof(k)) = lbuf(this%local_gs_dof(k)) &
               + u(this%local_dof_gs(k))
       end do
       do k = 1, this%nlocal
          u(this%local_dof_gs(k)) = lbuf(this%local_gs_dof(k))
       end do
    end if

    if (this%nshared_blks > 0) then
       allocate(sbuf(this%nshared_blks))
       sbuf = 0.0_rp
       do k = 1, this%nshared
          sbuf(this%shared_gs_dof(k)) = sbuf(this%shared_gs_dof(k)) &
               + u(this%shared_dof_gs(k))
       end do

       if (present(gs)) then
          tid = 0
          deps = C_NULL_PTR
          strm = C_NULL_PTR
          call gs%comm%nbsend(sbuf, this%nshared_blks, tid, deps, strm)
          call gs%comm%nbrecv(tid)
          call gs%comm%nbwait(sbuf, this%nshared_blks, GS_OP_ADD, strm)
       end if

       do k = 1, this%nshared
          u(this%shared_dof_gs(k)) = sbuf(this%shared_gs_dof(k))
       end do
    end if
  end subroutine amge_gs_map_apply

  !> Global id of each shared gather slot, in gather-buffer order --
  !! the key list a gs_comm_t / rendezvous exchange is built from.
  subroutine amge_gs_map_shared_ids(this, ids, n)
    type(amge_gs_map_t), intent(in) :: this
    integer(i4), allocatable, intent(out) :: ids(:)
    integer(i4), intent(out) :: n
    n = this%nshared_blks
    allocate(ids(max(n,1)))
    if (n > 0) ids(1:n) = this%shared_gid(1:n)
  end subroutine amge_gs_map_shared_ids

  !> Seed the finest level's shared flags from the solver's own dofmap.
  !! Call once on level 0; build_next_level then inherits the flags
  !! downward. @a dof can be at ANY polynomial order (not only a literal
  !! Q1 dofmap): only its 8 hex-corner entries per element are read, at
  !! the same tensor-grid corners gs_init_mapping reads off any dofmap and
  !! macro_mesh_init_hex assigns pts to (Neko hex_t vertex order: x
  !! fastest, then y, then z), which is exactly AMGe's level-0 CSR slot
  !! order -- so no reshape of the full (lx,ly,lz,nelv) array is needed or
  !! correct (that would only work by coincidence at lx=ly=lz=2).
  !! @param dof  the solver's dofmap for this mesh
  subroutine amge_mesh_set_shared_from_dofmap(n_verts, nelm, elm_vtx_ptr, &
                                              elm_vtx_idx, dof, shared_vtx)
    integer(i4), intent(in) :: n_verts, nelm
    integer(i4), intent(in) :: elm_vtx_ptr(:), elm_vtx_idx(:)
    type(dofmap_t), intent(in) :: dof
    logical, allocatable, intent(inout) :: shared_vtx(:)
    integer(i4) :: e, t, off, lx, ly, lz
    logical :: corner(8)

    if (allocated(shared_vtx)) deallocate(shared_vtx)
    allocate(shared_vtx(n_verts))
    shared_vtx = .false.

    lx = dof%Xh%lx; ly = dof%Xh%ly; lz = dof%Xh%lz
    do e = 1, nelm
       corner = [ dof%shared_dof(1, 1, 1, e), dof%shared_dof(lx, 1, 1, e), &
            dof%shared_dof(1, ly, 1, e), dof%shared_dof(lx, ly, 1, e), &
            dof%shared_dof(1, 1, lz, e), dof%shared_dof(lx, 1, lz, e), &
            dof%shared_dof(1, ly, lz, e), dof%shared_dof(lx, ly, lz, e) ]
       off = elm_vtx_ptr(e)
       do t = 1, 8
          if (corner(t)) shared_vtx(elm_vtx_idx(off + t)) = .true.
       end do
    end do
  end subroutine amge_mesh_set_shared_from_dofmap

  subroutine amge_gs_map_free(this)
    class(amge_gs_map_t), intent(inout) :: this
    if (allocated(this%local_dof_gs)) deallocate(this%local_dof_gs)
    if (allocated(this%local_gs_dof)) deallocate(this%local_gs_dof)
    if (allocated(this%local_blk_len)) deallocate(this%local_blk_len)
    if (allocated(this%local_blk_off)) deallocate(this%local_blk_off)
    if (allocated(this%shared_dof_gs)) deallocate(this%shared_dof_gs)
    if (allocated(this%shared_gs_dof)) deallocate(this%shared_gs_dof)
    if (allocated(this%shared_blk_len)) deallocate(this%shared_blk_len)
    if (allocated(this%shared_blk_off)) deallocate(this%shared_blk_off)
    if (allocated(this%shared_gid)) deallocate(this%shared_gid)
    this%ndofs = 0; this%nlocal = 0; this%nshared = 0
    this%nlocal_blks = 0; this%nshared_blks = 0
  end subroutine amge_gs_map_free

  !> ord(1:n) = permutation sorting key ascending (stable enough here)
  subroutine sort_index(key, n, ord)
    integer(i4), intent(in) :: key(:)
    integer(i4), intent(in) :: n
    integer(i4), intent(out) :: ord(:)
    integer(i4) :: i
    do i = 1, n
       ord(i) = i
    end do
    if (n > 1) call qsort_ord(key, ord, 1, n)
  end subroutine sort_index

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

end module amge_gs
