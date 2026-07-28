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

  !> Internal machine-precision checks of one transition.
  subroutine check_transition(lvf, lvc)
    type(amge_level_t), intent(in) :: lvf, lvc
    real(rp), allocatable :: af(:,:), ac(:,:), p(:,:), pap(:,:)
    real(rp) :: gerr, cerr, rerr, serr
    integer(i4) :: m, a, q
    call assemble_dense(lvf, af)
    call assemble_dense(lvc, ac)
    call build_p_dense(lvc%tr, p)
    allocate(pap(lvc%tr%n_coarse, lvc%tr%n_coarse))
    pap = matmul(transpose(p), matmul(af, p))
    gerr = sqrt(sum((pap - ac)**2))
    ! conformity: each PiTilde block equals the conforming P restricted
    cerr = 0.0_rp
    do m = 1, lvc%tr%n_melm
       associate (mp => lvc%tr%maps(m))
         do q = 1, size(mp%cdofs)
            do a = 1, size(mp%fdofs)
               cerr = max(cerr, abs(mp%PiTilde%x(a, q) - p(mp%fdofs(a), mp%cdofs(q))))
            end do
         end do
       end associate
    end do
    rerr = maxval(abs(sum(ac, dim=2)))
    serr = maxval(abs(ac - transpose(ac)))
    write(*, '("   [check] ||P^T A P - A_c||_F = ", ES10.3, "  conformity = ", ES10.3, &
         & "  rowsum = ", ES10.3, "  sym = ", ES10.3)') gerr, cerr, rerr, serr
  end subroutine check_transition

  !> Invariants: every macroedge chain terminates at macrovertices (or is
  !! a closed loop through its breakpoint), and every macroface's
  !! bounding macroedge chains lie inside its vertex set.
  subroutine check_invariants(topo)
    type(macro_topology_t), intent(in) :: topo
    integer(i4) :: k, a, b, v
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
    write(*, '("   [check] chain endpoints are mv: ", L1, ' // &
         '";  rim chains inside face verts: ", L1)') ok1, ok2
  end subroutine check_invariants

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
  subroutine amge_fill_AM_from_ax(lvl, ax, coef, Xh, msh, blst)
    type(amge_level_t), intent(inout) :: lvl
    class(ax_t), intent(inout) :: ax
    type(coef_t), intent(inout) :: coef
    type(space_t), intent(inout) :: Xh
    type(mesh_t), intent(inout) :: msh
    type(bc_list_t), target, intent(inout) :: blst
    real(rp), allocatable :: u(:,:), w(:,:)
    real(rp), allocatable :: Ae(:,:)            ! full lxyz x lxyz element block
    integer(i4), allocatable :: corner(:)       ! lxyz-index of each of the 8 corners
    integer(i4) :: nxyz, e, col, i
    integer :: n

    nxyz = Xh%lxyz
    n = nxyz * msh%nelv
    call corner_node_indices(Xh, corner)        ! the 8 tensor-grid corners

    allocate(u(nxyz, msh%nelv), w(nxyz, msh%nelv), Ae(nxyz, nxyz))

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
       call ax%compute(w, u, coef, msh, Xh)       ! UNASSEMBLED element action
       ! w(:,e) is now column `col` of element e's full local matrix
       ! Would be nice if we could locally apply the BC here
       ! and have it track through the rest of the problem
       !call blst%apply_scalar(w, n)
       do e = 1, msh%nelv
          ! (handled per-element below to keep memory small; here we could
          !  stash into a big buffer, but we condense element-by-element)
       end do
       ! stash this column into a per-element accumulator
       call stash_column(col, w, msh%nelv, nxyz)
    end do

    ! ---- condense each element's full block onto its 8 corners ----
    do e = 1, msh%nelv
       call get_element_block(e, nxyz, Ae)        ! full lxyz x lxyz for element e
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
