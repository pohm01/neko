module amge
  use num_types, only : i4, rp, dp
  use comm
  use mpi_f08, only: MPI_Allreduce, MPI_MIN, MPI_IN_PLACE, MPI_INTEGER, MPI_Barrier
  use utils, only : neko_error
  use math, only : copy, col2, rzero, add2s1, glsc3
  use mesh, only : mesh_t
  use space, only : space_t
  use coefs, only : coef_t
  use ax_product, only : ax_t
  use bc_list, only: bc_list_t
  use matrix, only : matrix_t
  use amge_topology, only : macro_topology_t, macro_mesh_t, macro_mesh_init_hex
  use amge_coarsen, only : coarsen_level_3d, agglomerate_level
  use amge_level, only : amge_level_t, amge_vec_t, amge_level_transfer_t, &
       amge_apply
  use amge_utils, only : amge_axpy, &
       amge_gather, amge_scatter_add, amge_gs_placeholder, &
       assemble_dense, build_p_dense, check_transition, check_invariants, &
       check_spsd, check_constant_reproduction, &
       amge_fill_AM_from_ax, q1_hex
  use amge_gs, only : amge_mesh_set_shared_from_dofmap
  use amge_smoother, only : amge_smoother_wrapper_t, amge_smoother_alloc, &
       AMGE_SMOOTHER_CHEBY
  use profiler, only : profiler_start_region, profiler_end_region
  implicit none
  private

  !> A hierarchy for AMGe
  type, public :: amge_hierarchy_t
    integer :: nlvls = 2 !< number of levels in the hierarchy
    type(amge_level_t), allocatable :: lvl(:) !< amg levels in the hierarchy
    !> Smoother per level (a heterogeneous, per-level choice of concrete
    !! type -- see amge_smoother.f90's amge_smoother_wrapper_t for why
    !! this needs a wrapper rather than a plain polymorphic array).
    type(amge_smoother_wrapper_t), allocatable :: smoo(:)
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

  ! ================== hierarchy ==================

  !> Initialize a macroelement amg hierarchy from a neko mesh_t
  subroutine amge_hierarchy_init(this, ax, Xh, coef, msh, blst, nlvls, sm_itr, &
       smoother_kind)
    class(amge_hierarchy_t), intent(inout) :: this
    class(ax_t), intent(inout) :: ax
    type(coef_t), intent(inout) :: coef
    type(space_t), intent(inout) :: Xh
    type(mesh_t), target, intent(inout) :: msh
    type(bc_list_t), target, intent(inout) :: blst
    integer, intent(in) :: nlvls
    integer, intent(in) :: sm_itr
    !> Per-level smoother choice (AMGE_SMOOTHER_CHEBY/AMGE_SMOOTHER_JACOBI
    !! from amge_smoother.f90), indexed 0:nlvls-1. Defaults to Chebyshev
    !! on every level when not given.
    integer, intent(in), optional :: smoother_kind(0:)
    type(macro_topology_t) :: topo
    integer, allocatable :: part(:)
    integer :: l, nm, glb_min_elm, ierr
    this%nlvls = nlvls
    allocate(this%lvl(0:this%nlvls-1))

    ! Create finest level
    call amge_level_init_from_mesh(this%lvl(0), msh)
    call amge_mesh_set_shared_from_dofmap( this%lvl(0)%mmsh%n_verts, &
         this%lvl(0)%mmsh%n_elem, this%lvl(0)%elm_vtx_ptr, &
         this%lvl(0)%elm_vtx_idx, coef%dof, &
         this%lvl(0)%mmsh%shared_vtx)
    call amge_fill_AM_from_ax(this%lvl(0), ax, coef, Xh, msh, blst)
    call this%lvl(0)%data_init(0, sm_itr)
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
       call this%lvl(l)%data_init(l, sm_itr)
       write(*, '("   [check] rank ", I0, ": level ", I0, "->", I0, " elms: ", I0, " -> ", I0, &
            & " dofs: ", I0, " -> ", I0)') &
         pe_rank, (l-1), l, this%lvl(l-1)%nelm(), this%lvl(l)%nelm(), &
         this%lvl(l)%tr%n_fine, this%lvl(l)%tr%n_coarse
       call MPI_Barrier(NEKO_COMM, ierr)
       if (this%lvl(l)%tr%n_coarse == 0) then
          call neko_error("AMGe: coarse grid has no dofs. This is known topology bug.")
       end if
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

    ! Initialize the smoother on every level now that nlvls is final (an
    ! early exit above can shrink it from what was requested)
    allocate(this%smoo(0:this%nlvls-1))
    do l = 0, this%nlvls-1
       if (present(smoother_kind)) then
          call amge_smoother_alloc(this%smoo(l)%obj, smoother_kind(l))
       else
          call amge_smoother_alloc(this%smoo(l)%obj, AMGE_SMOOTHER_CHEBY)
       end if
       call this%smoo(l)%obj%init(this%lvl(l), l, sm_itr)
    end do
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

    if (.not.(msh%mpts == lvl%mmsh%n_verts)) then
      print *, "WARNING: vertex counts: neko msh", msh%mpts, "AMGe mmsh", lvl%mmsh%n_verts
    end if
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
    integer :: l, lmax, si
    character(len=2) :: lvl_name
    lmax = this%nlvls-1
    ! Traverse down hierarchy to coarse grid
    do l = 0, lmax-1
       write(lvl_name, '(I0)') l
       call profiler_start_region( "AMGe_level_" // trim(lvl_name))
       associate( lvl => this%lvl(l), lc => this%lvl(l+1) )
         call this%smoo(l)%obj%solve(lvl, lvl%x, lvl%b, zero_init=.true.)
         call calc_resid(lvl, lvl%r, lvl%x, lvl%b)
         call amge_restrict(this%lvl, l, lvl%r, lc%b)
         ! Assemble the per-macroelement contributions
         call lc%gsh%op(lc%b%x)
         call rzero(lc%x%x, lc%x%n_dofs)
       end associate
       call profiler_end_region( "AMGe_level_" // trim(lvl_name))
    end do
    ! Coarse grid solve
    write(lvl_name, '(I0)') lmax
    call profiler_start_region( "AMGe_level_" // trim(lvl_name))
    associate( lvl => this%lvl(lmax) )
      call this%smoo(lmax)%obj%solve(lvl, lvl%x, lvl%b, zero_init=.true.)
    end associate
    call profiler_end_region( "AMGe_level_" // trim(lvl_name))
    ! Traverse up hierarchy to fine grid
    do l = lmax-1, 0, -1
       write(lvl_name, '(I0)') l
       call profiler_start_region( "AMGe_level_" // trim(lvl_name))
       associate( lvl => this%lvl(l), lc => this%lvl(l+1) )
         call amge_prolong(this%lvl, l, lc%x, lvl%r) !r as a tmp workspace
         call amge_axpy(1.0_rp, lvl%r, lvl%x) !r as a tmp workspace
         call this%smoo(l)%obj%solve(lvl, lvl%x, lvl%b)
       end associate
       call profiler_end_region( "AMGe_level_" // trim(lvl_name))
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
    real(kind=rp) :: rnorm_pre, rnorm_post
    associate( lvl0 => this%amg%lvl(0) )
      call copy(lvl0%b%x, r, n) ! Need to convert to cell-local RHS
      call rzero(lvl0%x%x, n)

      ! weighted by lvl0%mult (1/duplication count), same convention as
      ! every other global reduction in AMGe (e.g. amge_cheby_power),
      ! since lvl0%b%x is cell-wise (duplicated) storage
      rnorm_pre = sqrt(glsc3(lvl0%b%x, lvl0%mult, lvl0%b%x, n))

      call this%amg%vcycle()

      ! true residual of the correction just computed, b - A*z
      call calc_resid(lvl0, lvl0%r, lvl0%x, lvl0%b)
      rnorm_post = sqrt(glsc3(lvl0%r%x, lvl0%mult, lvl0%r%x, n))

      write(*, '("   [check] rank ", I0, ": AMGe vcycle residual: before = ", &
           & ES12.5, "  after = ", ES12.5)') pe_rank, rnorm_pre, rnorm_post

      call copy(z, lvl0%x%x, n)
    end associate
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
