module amge
  use num_types, only : i4, rp, dp
  use comm
  use mpi_f08, only: MPI_Allreduce, MPI_MIN, MPI_IN_PLACE, MPI_INTEGER
  use utils, only : neko_error
  use math, only : copy, col2, rzero, add2s1
  use mesh, only : mesh_t
  use space, only : space_t
  use coefs, only : coef_t
  use ax_product, only : ax_t
  use bc_list, only: bc_list_t
  use matrix, only : matrix_t
  use amge_topology, only : macro_topology_t, macro_mesh_t, macro_mesh_init_hex
  use amge_coarsen, only : coarsen_level_3d, agglomerate_level
  use amge_level, only : amge_level_t, amge_vec_t, amge_level_transfer_t
  use amge_utils, only : amge_axpy, scale_by_valence, &
       amge_gather, amge_scatter_add, amge_gs_placeholder, &
       assemble_dense, build_p_dense, check_transition, check_invariants, &
       check_spsd, check_constant_reproduction, &
       amge_fill_AM_from_ax, q1_hex
  use amge_gs, only : amge_mesh_set_shared_from_dofmap
  implicit none
  private

  !> A hierarchy for AMGe
  type, public :: amge_hierarchy_t
    integer :: nlvls = 2 !< number of levels in the hierarchy
    type(amge_level_t), allocatable :: lvl(:) !< amg levels in the hierarchy
    integer :: target_agg_size = 8 !< target agg size
    integer :: min_grid_elm = 1 !< minimum elements on coarse grid
  contains
    procedure, pass(this) :: init => amge_hierarchy_init
    procedure, pass(this) :: vcycle => amge_flat_vcycle
  end type amge_hierarchy_t

  type, public :: amge_solver_t
    type(amge_hierarchy_t) :: amg
  contains
    procedure, pass(this) :: solve => amge_solve
  end type amge_solver_t

contains

  ! ================== matrix operator ==================

  !> Cell-local operator: primal in, UNASSEMBLED dual out. Each element
  !! applies its own local matrix to its own copies; no scatter.
  subroutine amge_apply(lvl, xin, yout)
    type(amge_level_t), intent(in) :: lvl
    type(amge_vec_t), intent(in) :: xin
    type(amge_vec_t), intent(inout) :: yout
    integer(i4) :: e, off, n
    do e = 1, lvl%nelm()
       off = lvl%elm_vtx_ptr(e)
       n = lvl%ndof_el(e)
       yout%x(off + 1 : off + n) = matmul(lvl%AM(e)%x, xin%x(off + 1 : off + n))
    end do
  end subroutine amge_apply

  ! ================== grid transfers ==================

  !> Cell-local prolongation, coarse level (l+1) -> fine level l. Each
  !! coarse macroelement applies its PiTilde to its own coarse copies and
  !! writes its children's fine blocks (one parent per child: no
  !! accumulation, no weights). Output is continuous across macroelement
  !! boundaries by trace compatibility.
  subroutine amge_prolong(lvls, l, xm_c, xm_f)
    type(amge_level_t), intent(inout) :: lvls(0:)
    integer(i4), intent(in) :: l
    type(amge_vec_t), intent(in) :: xm_c
    type(amge_vec_t), intent(inout) :: xm_f
    real(rp), allocatable :: um(:)
    integer(i4) :: m, k, i, chld, coff, foff, nc, cptr
    associate (lc => lvls(l + 1))
      do m = 1, lc%tr%n_melm
         coff = lc%elm_vtx_ptr(m)
         nc = lc%ndof_el(m)
         um = matmul(lc%tr%maps(m)%PiTilde%x, xm_c%x(coff + 1 : coff + nc))
         do k = 1, size(lc%tr%maps(m)%child_idx)
            chld = lc%tr%maps(m)%child_idx(k)
            cptr = lc%tr%maps(m)%chloc_ptr(k)
            foff = lvls(l)%elm_vtx_ptr(chld)
            do i = 1, lvls(l)%ndof_el(chld)
               xm_f%x(foff + i) = um(lc%tr%maps(m)%chloc_idx(cptr + i))
            end do
         end do
      end do
    end associate
  end subroutine amge_prolong

  !> Cell-local restriction of an assembled vector, fine level l ->
  !! coarse level (l+1). Each coarse macroelement collects its children's
  !! partial duals into its own frame and applies PiTilde^T.
  subroutine amge_restrict(lvls, l, rm_f, rm_c)
    type(amge_level_t), intent(inout) :: lvls(0:)
    integer(i4), intent(in) :: l
    type(amge_vec_t), intent(in) :: rm_f
    type(amge_vec_t), intent(inout) :: rm_c
    real(rp), allocatable :: z(:)
    integer(i4) :: m, k, i, chld, coff, foff, nc, nf, cptr
    rm_c%x = 0.0_rp
    associate (lc => lvls(l + 1))
      do m = 1, lc%tr%n_melm
         nf = lc%tr%maps(m)%PiTilde%get_nrows()
         nc = lc%ndof_el(m)
         if (allocated(z)) deallocate(z)
         allocate(z(nf))
         z = 0.0_rp
         do k = 1, size(lc%tr%maps(m)%child_idx)
            chld = lc%tr%maps(m)%child_idx(k)
            cptr = lc%tr%maps(m)%chloc_ptr(k)
            foff = lvls(l)%elm_vtx_ptr(chld)
            do i = 1, lvls(l)%ndof_el(chld)
               z(lc%tr%maps(m)%chloc_idx(cptr + i)) = rm_f%x(foff + i)
            end do
         end do
         ! eq-(62) duplication weight, then PiTilde^T
         do i = 1, nf
           z(i) = z(i) * lc%tr%winv(lc%tr%maps(m)%fdofs(i))
         end do
         coff = lc%elm_vtx_ptr(m)
         rm_c%x(coff + 1 : coff + nc) = matmul(transpose(lc%tr%maps(m)%PiTilde%x), z)
      end do
    end associate
  end subroutine amge_restrict

  !> Cell-local restriction of an UNASSEMBLED dual, fine level l ->
  !! coarse level (l+1). Each coarse macroelement collects its children's
  !! partial duals into its own frame and applies PiTilde^T. Reproduces
  !! the weighted conforming restriction without any weights (the
  !! partial sums are already partitioned among the children).
  subroutine amge_restrict_ua(lvls, l, rm_f, rm_c)
    type(amge_level_t), intent(inout) :: lvls(0:)
    integer(i4), intent(in) :: l
    type(amge_vec_t), intent(in) :: rm_f
    type(amge_vec_t), intent(inout) :: rm_c
    real(rp), allocatable :: z(:)
    integer(i4) :: m, k, i, chld, coff, foff, nc, nf, cptr
    rm_c%x = 0.0_rp
    associate (lc => lvls(l + 1))
      do m = 1, lc%tr%n_melm
         nf = lc%tr%maps(m)%PiTilde%get_nrows()
         nc = lc%ndof_el(m)
         if (allocated(z)) deallocate(z)
         allocate(z(nf))
         z = 0.0_rp
         do k = 1, size(lc%tr%maps(m)%child_idx)
            chld = lc%tr%maps(m)%child_idx(k)
            cptr = lc%tr%maps(m)%chloc_ptr(k)
            foff = lvls(l)%elm_vtx_ptr(chld)
            do i = 1, lvls(l)%ndof_el(chld)
               z(lc%tr%maps(m)%chloc_idx(cptr + i)) = &
                    z(lc%tr%maps(m)%chloc_idx(cptr + i)) + rm_f%x(foff + i)
            end do
         end do
         coff = lc%elm_vtx_ptr(m)
         rm_c%x(coff + 1 : coff + nc) = matmul(transpose(lc%tr%maps(m)%PiTilde%x), z)
      end do
    end associate
  end subroutine amge_restrict_ua

  ! ================== smoother ==================

  !> Literal l1-Jacobi sweep in cell-wise storage
  !! Works on assembled rhs
  subroutine amge_smooth_l1(lvl, dl1, u, bt, nu)
    type(amge_level_t), intent(inout) :: lvl
    real(rp), intent(in) :: dl1(:)            !< assembled l1 diagonal
    type(amge_vec_t), intent(inout) :: u      !< primal iterate (continuous)
    type(amge_vec_t), intent(in) :: bt        !< assembled cell-wise rhs
    integer(i4), intent(in) :: nu
    type(amge_vec_t) :: rt
    integer(i4) :: it, e, p, off
    call lvl%new_vec(rt)
    do it = 1, nu
       call amge_apply(lvl, u, rt)
       call lvl%gsh%op(rt%x)
       rt%x = bt%x - rt%x
       do e = 1, lvl%nelm()
          off = lvl%elm_vtx_ptr(e)
          do p = 1, lvl%ndof_el(e)
             if (abs(dl1(lvl%elm_vtx_idx(off + p))) > 0) then
                u%x(off + p) = u%x(off + p) + rt%x(off+p) / dl1(lvl%elm_vtx_idx(off + p))
             end if
          end do
       end do
       call lvl%gsh%op(u%x)
       call scale_by_valence(lvl, u)
    end do
    call rt%free()
  end subroutine amge_smooth_l1

  ! ================== hierarchy ==================

  !> Initialize a macroelement amg hierarchy from a neko mesh_t
  subroutine amge_hierarchy_init(this, ax, Xh, coef, msh, blst, nlvls)
    class(amge_hierarchy_t), intent(inout) :: this
    class(ax_t), intent(inout) :: ax
    type(coef_t), intent(inout) :: coef
    type(space_t), intent(inout) :: Xh
    type(mesh_t), target, intent(inout) :: msh
    type(bc_list_t), target, intent(inout) :: blst
    integer, intent(in) :: nlvls
    type(macro_topology_t) :: topo
    integer, allocatable :: part(:)
    integer :: l, nm, glb_min_elm
    this%nlvls = nlvls
    allocate(this%lvl(0:this%nlvls-1))

    ! Create finest level
    call amge_level_init_from_mesh(this%lvl(0), msh)
    call amge_mesh_set_shared_from_dofmap( this%lvl(0)%mmsh%n_verts, &
         this%lvl(0)%mmsh%n_elem, this%lvl(0)%elm_vtx_ptr, &
         this%lvl(0)%elm_vtx_idx, coef%dof, &
         this%lvl(0)%mmsh%shared_vtx)
    call amge_fill_AM_from_ax(this%lvl(0), ax, coef, Xh, msh, blst)
    call this%lvl(0)%data_init(0)
    call check_spsd(this%lvl(0))

    ! Build the rest of the hierarchy
    do l = 1, this%nlvls-1
       call agglomerate_level(this%lvl(l-1), this%target_agg_size, part, nm)
       ! Ghost-extend the rank-boundary topology at every coarsening step
       ! (amge_ghost.f90's CSR wire format and macro_mesh_splice_ghost are
       ! arity-generic, not tied to level 0's hex shape) -- see
       ! amge_ghost.f90's header for why a halo is needed at all, and
       ! amge_gs.f90 for the invariant this restores (a shared coarse dof
       ! must be a live, identically-classified unknown on every rank that
       ! shares it). coarsen_level_3d ignores this on single-rank runs.
       call coarsen_level_3d(this%lvl(l-1), part, nm, topo, this%lvl(l), &
                             use_ghost=.true.)
       call check_invariants(topo)
       call this%lvl(l)%data_init(l)
       write(*, '("level ", I0, "->", I0, " elms: ", I0, " -> ", I0, " dofs: ", I0, " -> ", I0)') &
         (l-1), l, this%lvl(l-1)%nelm(), this%lvl(l)%nelm(), &
         this%lvl(l)%tr%n_fine, this%lvl(l)%tr%n_coarse
       call check_transition(this%lvl(l-1), this%lvl(l))
       call check_spsd(this%lvl(l))
       call check_constant_reproduction(this%lvl(l))
       deallocate(part)
       ! Check if further coarsening is possible/wanted
       glb_min_elm = this%lvl(l)%mmsh%n_elem
       call MPI_Allreduce(MPI_IN_PLACE, glb_min_elm, 1, &
            MPI_INTEGER, MPI_MIN, NEKO_COMM)
       if (glb_min_elm <= this%min_grid_elm) exit
    end do
    if (l < this%nlvls-1) then
       this%nlvls = l+1
       !TODO: reduce allocation size of this%lvl
       !this%lvl = this%lvl(0:l)
    end if
  end subroutine amge_hierarchy_init

  !> Initialize a macroelement amg level from a neko mesh_t
  !! @param lvl amge_level to be filled
  !! @param msh neko mesh_t object
  subroutine amge_level_init_from_mesh(lvl, msh)
    type(amge_level_t), intent(inout) :: lvl
    type(mesh_t), intent(in) :: msh
    real(kind=dp) :: coords(8, 3)
    integer(i4), allocatable :: hv(:,:), g2l(:)
    integer(i4) :: e, k, v, npts
    call lvl%free()
    allocate(hv(8, msh%nelv))
    do e = 1, msh%nelv
       do k = 1, 8
          hv(k, e) = msh%elements(e)%e%pts(k)%p%id()
       end do
    end do
    call macro_mesh_init_hex(lvl%mmsh, msh%nelv, hv)

    npts = msh%mpts
    if (.not.(npts == lvl%mmsh%n_verts)) then
      print *, "VERT COUNTING ISSUE!"
    end if
    ! g2l is indexed by GLOBAL vertex id, which ranges up to msh%glb_mpts
    ! (the mesh-wide unique point count) -- NOT msh%mpts (this rank's own
    ! local unique point count, "npts" above). Sizing to npts was an
    ! out-of-bounds write on every rank whose elements reference global
    ! ids beyond its own small local count, i.e. on any real MPI partition;
    ! it only happened to be safe on a single rank, where glb_mpts == mpts.
    allocate(g2l(msh%glb_mpts))
    g2l = 0
    do v = 1, lvl%mmsh%n_verts
       g2l(lvl%mmsh%vert_id(v)) = v
    end do
    allocate(lvl%elm_vtx_ptr(msh%nelv + 1), lvl%elm_vtx_idx(8 * msh%nelv), lvl%AM(msh%nelv))
    lvl%elm_vtx_ptr(1) = 0
    do e = 1, msh%nelv
       lvl%elm_vtx_ptr(e + 1) = lvl%elm_vtx_ptr(e) + 8
       call lvl%AM(e)%init(8, 8)
       do k = 1, 8
          lvl%elm_vtx_idx(lvl%elm_vtx_ptr(e) + k) = g2l(hv(k, e))
          coords(k, :) = msh%elements(e)%e%pts(k)%p%x(:)
       end do
       call q1_hex(coords, lvl%AM(e)%x)
    end do
  end subroutine amge_level_init_from_mesh

  !> V-cycle on hierarchy, non-recursive version
  subroutine amge_flat_vcycle(this)
    class(amge_hierarchy_t), intent(inout) :: this
    integer :: l, lmax
    lmax = this%nlvls-1
    ! Traverse down hierarchy to coarse grid
    do l = 0, lmax-1
       associate( lvl => this%lvl(l), lc => this%lvl(l+1) )
         call amge_smooth_l1(lvl, lvl%dl1, lvl%x, lvl%b, 1)
         call calc_resid(lvl, lvl%r, lvl%x, lvl%b)
         call amge_restrict(this%lvl, l, lvl%r, lc%b)
         ! Assemble the per-macroelement contributions
         call lc%gsh%op(lc%b%x)
         call rzero(lc%x%x, lc%x%n_dofs)
       end associate
    end do
    ! Coarse grid solve
    associate( lvl => this%lvl(lmax) )
      call amge_smooth_l1(lvl, lvl%dl1, lvl%x, lvl%b, 1)
    end associate
    ! Traverse up hierarchy to fine grid
    do l = lmax-1, 0, -1
       associate( lvl => this%lvl(l), lc => this%lvl(l+1) )
         call amge_prolong(this%lvl, l, lc%x, lvl%r) !r as a tmp workspace
         call amge_axpy(1.0_rp, lvl%r, lvl%x) !r as a tmp workspace
         call amge_smooth_l1(lvl, lvl%dl1, lvl%x, lvl%b, 1)
       end associate
    end do
  end subroutine amge_flat_vcycle

  !> Solve function for AMGe. Matching interface for pc_t
  !! @param z The solution to be returned
  !! @param r The right-hand side
  !! @param n Number of dofs
  subroutine amge_solve(this, z, r, n)
    class(amge_solver_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: z
    real(kind=rp), dimension(n), intent(inout) :: r
    call copy(this%amg%lvl(0)%b%x, r, n) ! Need to convert to cell-local RHS
    call rzero(this%amg%lvl(0)%x%x, n)
    call this%amg%vcycle()
    call copy(z, this%amg%lvl(0)%x%x, n)
  end subroutine amge_solve

  ! ================== math operations ==================

  !> Wrapper function for calculating residual r = b - Ax
  !! @param lvl AMG level structure on which to compute residual
  !! @param r Output. Assembled residual
  !! @param x Input. Assembled element-local vector
  !! @param b Input. Assembled RHS
  subroutine calc_resid(lvl, r, x, b)
    type(amge_level_t), intent(inout) :: lvl
    type(amge_vec_t), intent(inout) :: r
    type(amge_vec_t), intent(in) :: x
    type(amge_vec_t), intent(in) :: b
    real(kind=rp), allocatable :: u(:)
    call amge_apply(lvl, x, r)
    call lvl%gsh%op(r%x)
    call add2s1(r%x, b%x, -1.0_rp, r%n_dofs)
  end subroutine calc_resid

end module amge
