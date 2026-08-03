module amge_utils
  use num_types, only : i4, rp, dp
  use utils, only : neko_error
  use math, only : copy, col2, rzero, add2s1
  use mesh, only : mesh_t
  use space, only : space_t
  use coefs, only : coef_t
  use ax_product, only : ax_t
  use bc_list, only: bc_list_t
  use matrix, only : matrix_t
  use amge_topology, only : macro_topology_t
  use amge_level, only : amge_level_t, amge_vec_t, amge_level_transfer_t
  use comm, only : NEKO_COMM, pe_rank, MPI_REAL_PRECISION
  use mpi_f08, only : MPI_Allreduce, MPI_SUM, MPI_MAX, MPI_IN_PLACE, MPI_INTEGER, &
       MPI_Barrier
  implicit none
  private

  ! module-level stash for probed columns (sketch convenience; a real
  ! version fuses this into the probe loop and avoids the O(nxyz^2 nelv)
  ! buffer). (nxyz, nelv, nxyz): column j of element e's local matrix.
  real(rp), allocatable :: g_cols(:,:,:)

  public :: amge_axpy, scale_by_valence
  public :: amge_gather, amge_scatter_add, amge_gs_placeholder
  public :: q1_hex
  public :: assemble_dense, build_p_dense, check_transition, check_invariants
  public :: check_spsd, check_constant_reproduction
  public :: assemble_dense_global, build_p_dense_global

  !> Above this many global dofs, the cross-rank debug checks below skip
  !! the dense O(glb_n^2) reconstruction (memory and Allreduce cost) rather
  !! than silently attempting it on a real production-sized level.
  integer(i4), parameter :: MAX_GLB_DEBUG_SIZE = 5000
  public :: amge_fill_AM_from_ax

contains

  !> y := y + alpha x (cell-wise, purely local).
  subroutine amge_axpy(alpha, x, y)
    real(rp), intent(in) :: alpha
    type(amge_vec_t), intent(in) :: x
    type(amge_vec_t), intent(inout) :: y
    y%x = y%x + alpha * x%x
  end subroutine amge_axpy

  subroutine scale_by_valence(lvl, u)
    type(amge_level_t), intent(inout) :: lvl
    type(amge_vec_t), intent(inout) :: u
    call col2(u%x, lvl%mult, u%n_dofs)
  end subroutine scale_by_valence

  !> Gather G: assembled (unique) vector -> duplicated cell-wise copies.
  subroutine amge_gather(lvl, u, v)
    type(amge_level_t), intent(in) :: lvl
    real(rp), intent(in) :: u(:)              !< length n_verts
    type(amge_vec_t), intent(inout) :: v
    integer(i4) :: e, p, off
    do e = 1, lvl%nelm()
       off = lvl%elm_vtx_ptr(e)
       do p = 1, lvl%ndof_el(e)
          v%x(off + p) = u(lvl%elm_vtx_idx(off + p))
       end do
    end do
  end subroutine amge_gather

  !> Scatter-add G^T (direct stiffness summation): sum duplicated copies
  !! of a dual vector into the unique assembled vector. The sole
  !! inter-element communication primitive (Neko: gs_t with GS_OP_ADD).
  subroutine amge_scatter_add(lvl, v, u)
    type(amge_level_t), intent(in) :: lvl
    type(amge_vec_t), intent(in) :: v
    real(rp), intent(out) :: u(:)             !< length n_verts
    integer(i4) :: e, p, off
    u = 0.0_rp
    do e = 1, lvl%nelm()
       off = lvl%elm_vtx_ptr(e)
       do p = 1, lvl%ndof_el(e)
          u(lvl%elm_vtx_idx(off + p)) = u(lvl%elm_vtx_idx(off + p)) + v%x(off + p)
       end do
    end do
  end subroutine amge_scatter_add

  !TODO: This GS operator is a placeholder.
  subroutine amge_gs_placeholder(lvl, v)
    type(amge_level_t), intent(inout) :: lvl
    type(amge_vec_t), intent(inout) :: v
    real(kind=rp), allocatable :: u(:)
    allocate(u(lvl%mmsh%n_verts))
    call amge_scatter_add(lvl, v, u)
    call amge_gather(lvl, u, v)
  end subroutine amge_gs_placeholder

  !TODO: FOR TESTING ONLY.
  !TODO: REPLACE WITH NEKO HEX ASSEMBLY FUNCTION
  !> Trilinear hex stiffness, 2x2x2 Gauss (port of the Octave element).
  subroutine q1_hex(coords, ak)
    real(rp), intent(in) :: coords(8, 3)
    real(rp), intent(out) :: ak(8, 8)
    real(rp), parameter :: xs(8) = [-1, 1,-1, 1,-1, 1,-1, 1] * 1.0_rp
    real(rp), parameter :: ys(8) = [-1,-1, 1, 1,-1,-1, 1, 1] * 1.0_rp
    real(rp), parameter :: zs(8) = [-1,-1,-1,-1, 1, 1, 1, 1] * 1.0_rp
    real(rp) :: gp(2), xi, eta, zeta, dn(3, 8), jm(3, 3), jinv(3, 3), dj, g(3, 8)
    integer(i4) :: a, b, c, i
    gp = [-1.0_rp, 1.0_rp] / sqrt(3.0_rp)
    ak = 0.0_rp
    do a = 1, 2
       do b = 1, 2
          do c = 1, 2
             xi = gp(a); eta = gp(b); zeta = gp(c)
             do i = 1, 8
                dn(1, i) = 0.125_rp * xs(i) * (1 + eta * ys(i)) * (1 + zeta * zs(i))
                dn(2, i) = 0.125_rp * ys(i) * (1 + xi * xs(i)) * (1 + zeta * zs(i))
                dn(3, i) = 0.125_rp * zs(i) * (1 + xi * xs(i)) * (1 + eta * ys(i))
             end do
             jm = matmul(dn, coords)
             dj = jm(1,1)*(jm(2,2)*jm(3,3) - jm(2,3)*jm(3,2)) &
                - jm(1,2)*(jm(2,1)*jm(3,3) - jm(2,3)*jm(3,1)) &
                + jm(1,3)*(jm(2,1)*jm(3,2) - jm(2,2)*jm(3,1))
             if (dj .le. 0.0_rp) error stop 'non-positive hex Jacobian'
             jinv(1,1) =  (jm(2,2)*jm(3,3) - jm(2,3)*jm(3,2)) / dj
             jinv(1,2) = -(jm(1,2)*jm(3,3) - jm(1,3)*jm(3,2)) / dj
             jinv(1,3) =  (jm(1,2)*jm(2,3) - jm(1,3)*jm(2,2)) / dj
             jinv(2,1) = -(jm(2,1)*jm(3,3) - jm(2,3)*jm(3,1)) / dj
             jinv(2,2) =  (jm(1,1)*jm(3,3) - jm(1,3)*jm(3,1)) / dj
             jinv(2,3) = -(jm(1,1)*jm(2,3) - jm(1,3)*jm(2,1)) / dj
             jinv(3,1) =  (jm(2,1)*jm(3,2) - jm(2,2)*jm(3,1)) / dj
             jinv(3,2) = -(jm(1,1)*jm(3,2) - jm(1,2)*jm(3,1)) / dj
             jinv(3,3) =  (jm(1,1)*jm(2,2) - jm(1,2)*jm(2,1)) / dj
             g = matmul(jinv, dn)
             ak = ak + matmul(transpose(g), g) * dj
          end do
       end do
    end do
  end subroutine q1_hex

  !> Dense assembly of a level's operator.
  subroutine assemble_dense(lvl, a)
    type(amge_level_t), intent(in) :: lvl
    real(rp), allocatable, intent(out) :: a(:,:)
    integer(i4) :: e, n, i, j, li, lj
    allocate(a(lvl%mmsh%n_verts, lvl%mmsh%n_verts))
    a = 0.0_rp
    do e = 1, lvl%mmsh%n_elem
       n = lvl%ndof_el(e)
       do j = 1, n
          lj = lvl%elm_vtx_idx(lvl%elm_vtx_ptr(e) + j)
          do i = 1, n
             li = lvl%elm_vtx_idx(lvl%elm_vtx_ptr(e) + i)
             a(li, lj) = a(li, lj) + lvl%AM(e)%x(i, j)
          end do
       end do
    end do
  end subroutine assemble_dense

  !> Dense conforming prolongation from the transfer blocks (eq 62).
  subroutine build_p_dense(tr, p)
    type(amge_level_transfer_t), intent(in) :: tr
    real(rp), allocatable, intent(out) :: p(:,:)
    integer(i4) :: m, a, q
    allocate(p(tr%n_fine, tr%n_coarse))
    p = 0.0_rp
    do m = 1, tr%n_melm
       associate (mp => tr%maps(m))
         do q = 1, size(mp%cdofs)
            do a = 1, size(mp%fdofs)
               p(mp%fdofs(a), mp%cdofs(q)) = &
                    p(mp%fdofs(a), mp%cdofs(q)) + mp%PiTilde%x(a, q)
            end do
         end do
       end associate
    end do
    do a = 1, tr%n_fine
       p(a, :) = p(a, :) * tr%winv(a)
    end do
  end subroutine build_p_dense

  !> Upper bound on the range of GLOBAL ids appearing in local_ids, across
  !! every rank -- large enough to safely index a dense array by global id
  !! directly (ids need not be dense/contiguous globally; unused rows/cols
  !! are simply always zero).
  function glb_size(local_ids) result(glb_n)
    integer(i4), intent(in) :: local_ids(:)
    integer(i4) :: glb_n
    glb_n = maxval(local_ids)
    call MPI_Allreduce(MPI_IN_PLACE, glb_n, 1, MPI_INTEGER, MPI_MAX, NEKO_COMM)
  end function glb_size

  !> Cross-rank dense assembly of a level's operator, indexed by GLOBAL
  !! vertex id: each rank fills only the entries its own local elements
  !! touch (at their global-id positions), then one Allreduce(SUM) totals
  !! every rank's partial contribution at a shared (rank-boundary) id --
  !! exactly the piece assemble_dense (rank-local only) is missing. Debug
  !! only: O(glb_n^2) memory/communication, skipped above
  !! MAX_GLB_DEBUG_SIZE.
  !! @param ok  .false. if skipped for size; caller should not trust a_glb
  subroutine assemble_dense_global(lvl, a_glb, ok)
    type(amge_level_t), intent(in) :: lvl
    real(rp), allocatable, intent(out) :: a_glb(:,:)
    logical, intent(out) :: ok
    integer(i4) :: e, n, i, j, li, lj, glb_n

    glb_n = glb_size(lvl%mmsh%vert_id)
    ok = glb_n .le. MAX_GLB_DEBUG_SIZE
    if (.not. ok) return

    allocate(a_glb(glb_n, glb_n))
    a_glb = 0.0_rp
    do e = 1, lvl%mmsh%n_elem
       n = lvl%ndof_el(e)
       do j = 1, n
          lj = lvl%mmsh%vert_id(lvl%elm_vtx_idx(lvl%elm_vtx_ptr(e) + j))
          do i = 1, n
             li = lvl%mmsh%vert_id(lvl%elm_vtx_idx(lvl%elm_vtx_ptr(e) + i))
             a_glb(li, lj) = a_glb(li, lj) + lvl%AM(e)%x(i, j)
          end do
       end do
    end do
    call MPI_Allreduce(MPI_IN_PLACE, a_glb, glb_n * glb_n, MPI_REAL_PRECISION, &
         MPI_SUM, NEKO_COMM)
  end subroutine assemble_dense_global

  !> Cross-rank dense conforming prolongation (eq 62), indexed by GLOBAL
  !! ids on both the fine (row) and coarse (column) side. Each rank's
  !! per-macroelement sum-then-scale-by-winv is already the correctly
  !! weighted PARTIAL contribution at a shared fine dof (winv itself is
  !! now cross-rank-correct, see amge_gs_correct_shared_count); summing
  !! those partials across ranks -- rather than using only one rank's
  !! partial, as build_p_dense does -- recovers the true global P. Debug
  !! only: see assemble_dense_global's caveats.
  !! @param ok  .false. if skipped for size; caller should not trust p_glb
  subroutine build_p_dense_global(lvf, lvc, p_glb, ok)
    type(amge_level_t), intent(in) :: lvf, lvc
    real(rp), allocatable, intent(out) :: p_glb(:,:)
    logical, intent(out) :: ok
    integer(i4) :: fv_loc, cv_loc, glb_nf, glb_nc, fv, cv
    real(rp), allocatable :: p_loc(:,:)

    glb_nf = glb_size(lvf%mmsh%vert_id)
    glb_nc = glb_size(lvc%mmsh%vert_id)
    ok = (glb_nf .le. MAX_GLB_DEBUG_SIZE) .and. (glb_nc .le. MAX_GLB_DEBUG_SIZE)
    if (.not. ok) return

    ! p_loc is already this rank's COMPLETE local partial (build_p_dense
    ! itself sums every local macroelement touching a given (fine, coarse)
    ! pair, then scales by winv) -- placing each of ITS entries once at
    ! their global position is the correct per-rank contribution; looping
    ! per-macroelement again here would double count any (fine, coarse)
    ! pair touched by more than one local macroelement.
    call build_p_dense(lvc%tr, p_loc)
    allocate(p_glb(glb_nf, glb_nc))
    p_glb = 0.0_rp
    do cv_loc = 1, lvc%mmsh%n_verts
       cv = lvc%mmsh%vert_id(cv_loc)
       do fv_loc = 1, lvf%mmsh%n_verts
          fv = lvf%mmsh%vert_id(fv_loc)
          p_glb(fv, cv) = p_glb(fv, cv) + p_loc(fv_loc, cv_loc)
       end do
    end do
    call MPI_Allreduce(MPI_IN_PLACE, p_glb, glb_nf * glb_nc, MPI_REAL_PRECISION, &
         MPI_SUM, NEKO_COMM)
  end subroutine build_p_dense_global

  !> Internal machine-precision checks of one transition.
  subroutine check_transition(lvf, lvc)
    type(amge_level_t), intent(in) :: lvf, lvc
    real(rp), allocatable :: af(:,:), ac(:,:), p(:,:), pap(:,:)
    real(rp) :: gerr, cerr, rerr, serr
    integer(i4) :: m, a, q, fv, cv, glb_nf, glb_nc
    integer :: ierr
    logical :: ok_a, ok_p

    ! glb_size is collective (MPI_Allreduce): every rank must call it, even
    ! though only rank 0 prints -- gating the CALL (not just the print) to
    ! one rank would deadlock the others waiting in the reduction.
    glb_nf = glb_size(lvf%mmsh%vert_id)
    glb_nc = glb_size(lvc%mmsh%vert_id)
    call assemble_dense_global(lvc, ac, ok_a)
    call build_p_dense_global(lvf, lvc, p, ok_p)
    if (.not. (ok_a .and. ok_p)) then
       if (pe_rank .eq. 0) write(*, '("   [check] rank ", I0, ": skipped: glb_size(fine) = ", &
            & I0, "  glb_size(coarse) = ", I0, "  (cap = ", I0, ")")') &
            pe_rank, glb_nf, glb_nc, MAX_GLB_DEBUG_SIZE
       ! ok_a/ok_p are already globally agreed (derived from all-reduced
       ! glb_nf/glb_nc), so every rank takes this branch together -- safe
       ! to call the barrier here unconditionally
       call MPI_Barrier(NEKO_COMM, ierr)
       return
    end if
    call assemble_dense_global(lvf, af, ok_a)

    allocate(pap(size(ac,1), size(ac,2)))
    pap = matmul(transpose(p), matmul(af, p))
    gerr = sqrt(sum((pap - ac)**2))
    ! conformity: each rank's own PiTilde block must equal the TRUE global
    ! conforming P restricted to its (translated) global positions -- with
    ! p now cross-rank-summed, this is a genuine check that every rank's
    ! trace maps agree with every other rank's on a shared macroedge/face
    ! (eq. 56), not merely self-consistent with this rank's own data.
    cerr = 0.0_rp
    block
      integer(i4) :: worst_m, worst_a, worst_q, mm, aa, qq
      worst_m = 0; worst_a = 0; worst_q = 0
      do m = 1, lvc%tr%n_melm
         associate (mp => lvc%tr%maps(m))
           do q = 1, size(mp%cdofs)
              cv = lvc%mmsh%vert_id(mp%cdofs(q))
              do a = 1, size(mp%fdofs)
                 fv = lvf%mmsh%vert_id(mp%fdofs(a))
                 if (abs(mp%PiTilde%x(a, q) - p(fv, cv)) .gt. cerr) then
                    cerr = abs(mp%PiTilde%x(a, q) - p(fv, cv))
                    worst_m = m; worst_a = a; worst_q = q
                 end if
              end do
           end do
         end associate
      end do
      ! TEMP DEBUG: dump every macroelement's own contribution at the
      ! worst-offending (fv, cv) pair, to see whether they genuinely
      ! disagree (real trace-incompatibility) or whether only one
      ! macroelement contributes non-zero winv-scaled mass.
      if (cerr .gt. 1.0e-8_rp .and. pe_rank .eq. 0 .and. worst_m .gt. 0) then
         associate (mp => lvc%tr%maps(worst_m))
           cv = lvc%mmsh%vert_id(mp%cdofs(worst_q))
           fv = lvf%mmsh%vert_id(mp%fdofs(worst_a))
         end associate
         write(*, '("   [debug] worst conformity: fv=", I0, " cv=", I0, &
              & " p_glb=", ES14.6, " w(fv)=", ES14.6)') fv, cv, p(fv, cv), &
              1.0_rp / lvc%tr%winv(lvc%tr%maps(worst_m)%fdofs(worst_a))
         do mm = 1, lvc%tr%n_melm
            associate (mp2 => lvc%tr%maps(mm))
              ! is fv itself a coarse dof (macrovertex) of this macroelement?
              do qq = 1, size(mp2%cdofs)
                 if (lvc%mmsh%vert_id(mp2%cdofs(qq)) .eq. fv) &
                      write(*, '("   [debug]   macroelm ", I0, " has fv AS ITS OWN cdof q=", &
                           & I0)') mm, qq
              end do
              do qq = 1, size(mp2%cdofs)
                 if (lvc%mmsh%vert_id(mp2%cdofs(qq)) .ne. cv) cycle
                 do aa = 1, size(mp2%fdofs)
                    if (lvf%mmsh%vert_id(mp2%fdofs(aa)) .ne. fv) cycle
                    write(*, '("   [debug]   macroelm ", I0, " Pi(a=", I0, &
                         & ",q=", I0, ") = ", ES14.6)') mm, aa, qq, mp2%PiTilde%x(aa, qq)
                 end do
              end do
            end associate
         end do
      end if
    end block
    rerr = maxval(abs(sum(ac, dim=2)))
    serr = maxval(abs(ac - transpose(ac)))
    if (pe_rank .eq. 0) then
       write(*, '("   [check] rank ", I0, ": ||P^T A P - A_c||_F = ", ES10.3, &
            & "  conformity = ", ES10.3, "  rowsum = ", ES10.3, "  sym = ", ES10.3)') &
            pe_rank, gerr, cerr, rerr, serr
    end if
    ! called by every rank unconditionally, not just rank 0, or the
    ! barrier would deadlock
    call MPI_Barrier(NEKO_COMM, ierr)
  end subroutine check_transition

  !> Invariants: every macroedge chain terminates at macrovertices (or is
  !! a closed loop through its breakpoint), and every macroface's
  !! bounding macroedge chains lie inside its vertex set.
  subroutine check_invariants(topo)
    type(macro_topology_t), intent(in) :: topo
    integer(i4) :: k, a, b, v
    integer :: ierr
    logical :: ok1, ok2, found
    ok1 = .true.
    do k = 1, topo%n_medge
       associate (ch => topo%medge(k)%chain)
         ok1 = ok1 .and. topo%is_mv(ch(1)) .and. topo%is_mv(ch(size(ch)))
       end associate
    end do
    ok2 = .true.
    do k = 1, topo%n_mface
       do a = 1, size(topo%mface(k)%bnd_medge)
          associate (ch => topo%medge(topo%mface(k)%bnd_medge(a))%chain)
            do b = 1, size(ch)
               v = ch(b)
               found = any(topo%mface(k)%verts .eq. v)
               ok2 = ok2 .and. found
            end do
          end associate
       end do
    end do
    write(*, '("   [check] rank ", I0, ": chain endpoints are mv: ", L1, ' // &
         '";  rim chains inside face verts: ", L1)') pe_rank, ok1, ok2
    call MPI_Barrier(NEKO_COMM, ierr)
  end subroutine check_invariants

  !> SPSD check: eq.(3)/(6) of the AMGe theory note requires A^ell_M >= 0
  !! at every level, for every macroelement. Reports the worst (most
  !! negative) eigenvalue across all of a level's local matrices; should
  !! be >= -tol (0, up to roundoff).
  subroutine check_spsd(lvl)
    type(amge_level_t), intent(in) :: lvl
    real(rp) :: worst
    integer(i4) :: e
    integer :: ierr
    worst = huge(1.0_rp)
    do e = 1, lvl%nelm()
       worst = min(worst, min_eigenvalue(lvl%AM(e)%x))
    end do
    write(*, '("   [check] rank ", I0, ": SPSD: min eig over ", I0, &
         & " local matrices = ", ES10.3)') pe_rank, lvl%nelm(), worst
    call MPI_Barrier(NEKO_COMM, ierr)
  end subroutine check_spsd

  !> Constant reproduction: Pi~_M must map the all-ones coarse vector to
  !! the all-ones fine vector on every macroelement (eq. 50), i.e. every
  !! row of PiTilde sums to 1. This is what lets the fine-level nullspace
  !! (constant functions, A_K 1 = 0) propagate correctly up the hierarchy;
  !! together with the rowsum check already in check_transition (A_c 1 = 0,
  !! checked at the assembled/global level), a failure here versus there
  !! isolates whether a discrepancy is in the interpolation itself or in
  !! the coarse operator.
  subroutine check_constant_reproduction(lvc)
    type(amge_level_t), intent(in) :: lvc
    real(rp) :: worst
    integer(i4) :: m, a
    integer :: ierr
    worst = 0.0_rp
    do m = 1, lvc%tr%n_melm
       associate (mp => lvc%tr%maps(m))
         do a = 1, size(mp%fdofs)
            worst = max(worst, abs(sum(mp%PiTilde%x(a, :)) - 1.0_rp))
         end do
       end associate
    end do
    write(*, '("   [check] rank ", I0, ": constant reproduction: max |rowsum(Pi~) - 1| = ", &
         & ES10.3)') pe_rank, worst
    call MPI_Barrier(NEKO_COMM, ierr)
  end subroutine check_constant_reproduction

  !> Smallest eigenvalue of a symmetric matrix (LAPACK dsyev). Local copy
  !! of amge_coarsen's own helper -- not worth a cross-module dependency
  !! for one ~15-line function.
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
    if (info .ne. 0) call neko_error('amge_utils: dsyev failed (min_eigenvalue)')
    lam_min = minval(w)
  end function min_eigenvalue

  ! ================== fill from ax_t ==================
  !> SKETCH: filling the AMGe fine-level element matrices AM(e) from a real
  !! Neko matrix-free operator (ax_t), instead of the placeholder q1_hex.
  !!
  !! The AMGe hierarchy is VERTEX-based: each hex contributes an 8x8 local
  !! matrix on its corner dofs. A Neko ax_t, however, acts on the full GLL
  !! tensor grid (lxyz = lx*ly*lz nodes per element). Two coherent routes:
  !!
  !!  (A) PROBE + STATIC CONDENSATION (this sketch, amge_fill_AM_from_ax):
  !!      materialize each element's full lxyz x lxyz matrix by probing the
  !!      matrix-free operator with unit local vectors, then Schur-complement
  !!      the non-corner nodes onto the 8 corners. Gives an 8x8 AM(e) that is
  !!      spectrally faithful to the high-order operator -- the right choice
  !!      when AMGe preconditions the actual SEM system.
  !!
  !!  (B) DIRECT Q1 REBUILD (amge_fill_AM_q1_from_coef, sketched in comments):
  !!      skip ax_t entirely and assemble the 8x8 trilinear-vertex operator
  !!      from coef's geometry. Cheaper, but it preconditions a Q1 proxy, not
  !!      the real operator. This is what the placeholder q1_hex approximates.
  !!
  !! Probing is embarrassingly parallel over (element, local column) and uses
  !! only the element-local action, so the ax_t is called with a single-element
  !! or full-field unit vector and NO gather-scatter (we want the UNASSEMBLED
  !! element blocks, exactly the A_K of the AMGe splitting).

  !> Fill lvl%AM(e) (8x8 corner-vertex blocks) from a matrix-free ax_t.
  !! @param lvl   AMGe fine level (mmsh + elm_vtx already built; AM filled here)
  !! @param ax    matrix-free operator (ax_t)
  !! @param coef  geometry/metric coefficients
  !! @param Xh    function space (gives lxyz, and the corner-node indices)
  !! @param msh   Neko mesh
  !! DIRICHLET BOUNDARY CONDITIONS. Periodic dofs are already identified
  !! (via coef%dof, same as amge_gs.f90's shared_vtx) and Neumann needs no
  !! special treatment (a natural condition in the weak form ax_t already
  !! encodes). Homogeneous Dirichlet (the residual/correction this
  !! hierarchy solves for is exactly 0 there) is applied HERE, once, at
  !! the finest level, by zeroing the row AND column of every Dirichlet
  !! local dof in each element's full lxyz x lxyz block before it is ever
  !! condensed or coarsened -- the standard FEM Dirichlet-elimination
  !! trick. Everything downstream (agglomeration, macrovertex/macroedge/
  !! macroface extraction, schur_extend, the smoother) is generic in
  !! AM(e)/mesh connectivity and needs no BC-awareness of its own: once a
  !! dof's row/col is isolated here, schur_extend's pseudoinverse of the
  !! (now block-diagonal-at-that-dof) interior block automatically gives
  !! it a zero row of PiTilde at every subsequent level, so no coarse
  !! correction is ever prolonged onto it.
  subroutine amge_fill_AM_from_ax(lvl, ax, coef, Xh, msh, blst)
    type(amge_level_t), intent(inout) :: lvl
    class(ax_t), intent(inout) :: ax
    type(coef_t), intent(inout) :: coef
    type(space_t), intent(inout) :: Xh
    type(mesh_t), intent(inout) :: msh
    type(bc_list_t), target, intent(inout) :: blst
    real(rp), allocatable, target :: u(:,:), w(:,:)
    real(rp), allocatable, target :: sentinel(:,:)
    real(rp), pointer :: u_flat(:), w_flat(:), sentinel_flat(:)
    real(rp), allocatable :: Ae(:,:)            ! full lxyz x lxyz element block
    integer(i4), allocatable :: corner(:)       ! lxyz-index of each of the 8 corners
    logical, allocatable :: is_dirichlet(:,:)   ! (nxyz, nelv) local dof mask
    integer(i4) :: nxyz, e, col, i
    integer :: n

    nxyz = Xh%lxyz
    n = nxyz * msh%nelv
    call corner_node_indices(Xh, corner)        ! the 8 tensor-grid corners

    allocate(u(nxyz, msh%nelv), w(nxyz, msh%nelv), Ae(nxyz, nxyz))
    ! blst%apply takes a flat (n) array; bc_list_t's dummy is explicit-
    ! shape, so a rank-2 actual doesn't sequence-associate through the
    ! GENERIC apply binding (gfortran rejects it) -- rank-remap a pointer
    ! onto the SAME storage instead of copying.
    u_flat(1:n) => u
    w_flat(1:n) => w

    ! ---- which (local node, element) pairs are strongly Dirichlet -----
    ! detect via a sentinel value: blst%apply overwrites a Dirichlet entry
    ! with its target g (0 for the residual/correction field this
    ! hierarchy solves); any entry that changed from the untouched
    ! sentinel was written by some bc in blst. Works regardless of g's
    ! actual value as long as it isn't 1 (true for a homogeneous target).
    allocate(sentinel(nxyz, msh%nelv), is_dirichlet(nxyz, msh%nelv))
    sentinel_flat(1:n) => sentinel
    sentinel = 1.0_rp
    call blst%apply(sentinel_flat, n)
    is_dirichlet = (sentinel .ne. 1.0_rp)
    deallocate(sentinel)

    ! ---- probe: one ax_t apply per local column materializes column j of
    !      every element's block simultaneously (all elements share the
    !      same local node numbering, so a single unit-vector field probes
    !      column j for all e at once). nxyz applies total, not nxyz*nelv. ----
    do e = 1, msh%nelv
       if (allocated(lvl%AM(e)%x)) call lvl%AM(e)%free()
       call lvl%AM(e)%init(8, 8)
    end do

    do col = 1, nxyz
       u = 0.0_rp
       u(col, :) = 1.0_rp                         ! unit vector at local node col, all elements
       ! never perturb a Dirichlet dof (zeros that column)
       call blst%apply(u_flat, n)
       call ax%compute(w, u, coef, msh, Xh)       ! UNASSEMBLED element action
       ! and ignore whatever raw coupling the operator reports AT a
       ! Dirichlet dof (zeros that row) -- together, the standard
       ! zero-row-and-column Dirichlet elimination, applied per fine
       ! element before any assembly/condensation happens
       call blst%apply(w_flat, n)
       ! w(:,e) is now column `col` of element e's full local matrix
       call stash_column(col, w, msh%nelv, nxyz)
       !! Note, we could handle element-by-element to keep memory small,
       !! but for now we just stash into a big buffer.
    end do

    ! ---- condense each element's full block onto its 8 corners ----
    do e = 1, msh%nelv
       call get_element_block(e, nxyz, Ae)        ! full lxyz x lxyz for element e
       ! the row/col zeroing above leaves a Dirichlet dof's diagonal zero
       ! too (an all-zero row+col+diagonal), which is singular -- and
       ! condense_to_corners' interior solve needs an invertible A_II.
       ! Restore a positive diagonal, decoupled from every other dof: the
       ! standard convention, and what makes the elimination self-
       ! propagate correctly through Schur condensation/coarsening.
       do i = 1, nxyz
          if (is_dirichlet(i, e)) Ae(i, i) = 1.0_rp
       end do
       call condense_to_corners(Ae, nxyz, corner, lvl%AM(e)%x)
    end do

    call clear_stash()
  end subroutine amge_fill_AM_from_ax

  !> Static condensation of a full element matrix onto its 8 corner dofs:
  !!   partition nodes into C (the 8 corners) and I (all the rest);
  !!   A_cc^schur = A_CC - A_CI A_II^{-1} A_IC.
  !! For a Q1 element (lx=2) I is empty and this returns A_CC unchanged.
  subroutine condense_to_corners(Ae, nxyz, corner, A8)
    real(rp), intent(in) :: Ae(:,:)
    integer(i4), intent(in) :: nxyz, corner(:)
    real(rp), intent(out) :: A8(8,8)
    logical, allocatable :: is_c(:)
    integer(i4), allocatable :: ci(:), ii(:)
    real(rp), allocatable :: Aii(:,:), Aic(:,:), Aci(:,:), tmp(:,:)
    integer(i4) :: i, nI, k
    allocate(is_c(nxyz)); is_c = .false.
    do i = 1, 8
       is_c(corner(i)) = .true.
    end do
    nI = nxyz - 8
    allocate(ci(8), ii(max(nI,1)))
    k = 0
    do i = 1, nxyz
       if (is_c(i)) then
          k = k + 1; ci(k) = i
       end if
    end do
    if (nI == 0) then
       A8 = Ae(ci, ci)
       return
    end if
    k = 0
    do i = 1, nxyz
       if (.not. is_c(i)) then
          k = k + 1; ii(k) = i
       end if
    end do
    allocate(Aii(nI,nI), Aic(nI,8), Aci(8,nI), tmp(nI,8))
    Aii = Ae(ii, ii); Aic = Ae(ii, ci); Aci = Ae(ci, ii)
    ! tmp = Aii^{-1} Aic  (SPD solve; in Neko use a Cholesky/LU or the
    ! existing dense solvers -- here a plain LAPACK solve)
    call spd_solve(Aii, Aic, tmp, nI, 8)
    A8 = Ae(ci, ci) - matmul(Aci, tmp)
  end subroutine condense_to_corners

  !> The 8 corner nodes of the (lx,ly,lz) tensor grid, in the SAME order
  !! as macro_mesh_init_hex expects the hex pts (Neko hex_t vertex order:
  !! x fastest, then y, then z). Local node index = i + (j-1)lx + (k-1)lx*ly.
  subroutine corner_node_indices(Xh, corner)
    type(space_t), intent(in) :: Xh
    integer(i4), allocatable, intent(out) :: corner(:)
    integer(i4) :: lx, ly, lz
    lx = Xh%lx; ly = Xh%ly; lz = Xh%lz
    allocate(corner(8))
    corner(1) = idx(1,  1,  1,  lx, ly)
    corner(2) = idx(lx, 1,  1,  lx, ly)
    corner(3) = idx(1,  ly, 1,  lx, ly)
    corner(4) = idx(lx, ly, 1,  lx, ly)
    corner(5) = idx(1,  1,  lz, lx, ly)
    corner(6) = idx(lx, 1,  lz, lx, ly)
    corner(7) = idx(1,  ly, lz, lx, ly)
    corner(8) = idx(lx, ly, lz, lx, ly)
  contains
    pure integer(i4) function idx(i,j,k,lx,ly)
      integer(i4), intent(in) :: i,j,k,lx,ly
      idx = i + (j-1)*lx + (k-1)*lx*ly
    end function idx
  end subroutine corner_node_indices

  ! ---- SPD solve X = A^{-1} B via LAPACK (dposv); Neko has dense solvers too ----
  subroutine spd_solve(A, B, X, n, nrhs)
    real(rp), intent(in) :: A(:,:), B(:,:)
    real(rp), intent(out) :: X(:,:)
    integer(i4), intent(in) :: n, nrhs
    real(rp), allocatable :: Af(:,:)
    integer(i4) :: info
    allocate(Af(n,n)); Af = A(1:n,1:n)
    X = B(1:n,1:nrhs)
    call dposv('U', n, nrhs, Af, n, X, n, info)
    if (info /= 0) call neko_error('amge_ax_extract: dposv failed (interior not SPD?)')
  end subroutine spd_solve

  ! ---- stash accessors (see g_cols declared at module scope) ----
  subroutine stash_column(col, w, nelv, nxyz)
    integer(i4), intent(in) :: col, nelv, nxyz
    real(rp), intent(in) :: w(:,:)
    if (.not. allocated(g_cols)) allocate(g_cols(nxyz, nelv, nxyz))
    g_cols(:, :, col) = w
  end subroutine stash_column
  subroutine get_element_block(e, nxyz, Ae)
    integer(i4), intent(in) :: e, nxyz
    real(rp), intent(out) :: Ae(:,:)
    integer(i4) :: col
    do col = 1, nxyz
       Ae(:, col) = g_cols(:, e, col)
    end do
  end subroutine get_element_block
  subroutine clear_stash()
    if (allocated(g_cols)) deallocate(g_cols)
  end subroutine clear_stash

end module amge_utils
