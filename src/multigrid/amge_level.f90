module amge_level
  use num_types, only : i4, rp
  use utils, only : neko_error
  use matrix, only : matrix_t
  use amge_gs, only : amge_gs_t
  use amge_topology, only : macro_mesh_t
  implicit none
  private

  !> Cell-wise (macroelement-local, duplicated) vector. x(:) is indexed
  !! by the owning level's elm_vtx_ptr: element e occupies
  !! x(elm_vtx_ptr(e)+1 : elm_vtx_ptr(e+1)), and position p within that
  !! block corresponds to global dof elm_vtx_idx(elm_vtx_ptr(e)+p).
  type, public :: amge_vec_t
     integer(i4) :: n_melm = 0        !< number of elements
     integer(i4) :: n_dofs = 0        !< total dofs counting duplicates (size of x)
     integer(i4) :: n_dofs_unique = 0 !< number of unique dofs
     logical :: assembled = .false. !< track if vector is assembled or not
     real(rp), allocatable :: x(:)
   contains
     procedure, pass(this) :: free => amge_vec_free
  end type amge_vec_t

  !> Per-macroelement transfer block (PiTilde: coarse level -> finer
  !! level, restricted to this macroelement and its children).
  type :: tr_map_t
     integer(i4), allocatable :: fdofs(:)   !< sorted fine dofs on this macroelm
     integer(i4), allocatable :: cdofs(:)   !< coarse dof ids on this macroelm
     type(matrix_t) :: PiTilde              !< (nF x nC), fine <- coarse
     integer(i4), allocatable :: child_idx(:)             !< member fine elements
     integer(i4), allocatable :: chloc_ptr(:), chloc_idx(:) !< child dof positions in fdofs
   contains
     procedure, pass(this) :: free => tr_map_free
  end type tr_map_t

  !> Transfer information between a level and the next finer level.
  type, public :: amge_level_transfer_t
     integer(i4) :: n_fine = 0, n_coarse = 0, n_melm = 0
     integer(i4) :: n_fine_unique = 0, n_coarse_unique = 0
     real(rp), allocatable :: winv(:)       !< 1 / duplication count (eq 62)
     type(tr_map_t), allocatable :: maps(:) !< per-macroelement transfer
   contains
     procedure, pass(this) :: free => amge_level_transfer_free
  end type amge_level_transfer_t

  !> A numerical level of the AMGe hierarchy.
  type, public :: amge_level_t
     integer :: level                           !< level index in hierarchy
     type(macro_mesh_t) :: mmsh                 !< level mesh (entity tables)
     integer(i4), allocatable :: elm_vtx_ptr(:) !< CSR dof list per element
     integer(i4), allocatable :: elm_vtx_idx(:)
     type(matrix_t), allocatable :: AM(:) !< local matrix per element
     logical :: topo_done = .true.        !< mmsh and topo info filled
     type(amge_level_transfer_t) :: tr    !< transfer to the finer level
     type(amge_gs_t) :: gsh               !< gather-scatter on level
     logical :: gsh_ready = .false.
     !! Workspace vectors on level
     type(amge_vec_t) :: x !< Element-local solution vector (Assembled)
     type(amge_vec_t) :: b !< Element-local RHS (Assembled)
     type(amge_vec_t) :: r !< Workspace (undetermined).
     !! Workspace vectors on level
     real(rp), allocatable :: mult(:) !< multiplicity of dof !TODO: this duplicates winv in tr (but after gather)
     !! Smoother info
     real(kind=rp), allocatable :: dl1(:) !< Storage for l1 diagonal (Assembled?)
     integer :: sm_itr !< Smoother iterations on level
   contains
     procedure, pass(this) :: data_init => amge_level_data_init
     procedure, pass(this) :: free => amge_level_free
     procedure, pass(this) :: ndof_el => amge_level_ndof_el
     procedure, pass(this) :: nelm => amge_level_nelm
     procedure, pass(this) :: new_vec => amge_level_new_vec
  end type amge_level_t

contains

  ! ================== level / vector helpers ==================

  !> Initialize workspace data on level
  !! Needs mmsh to already be filled
  !! Needs elm_vtx_ptr elm_vtx_idx to already be filled
  !! Needs AM to already be filled
  subroutine amge_level_data_init(this, level, sm_itr)
    class(amge_level_t), intent(inout) :: this
    integer, intent(in) :: level
    integer, intent(in) :: sm_itr
    logical, allocatable :: shared_vtx(:)
    this%level = level
    this%sm_itr = sm_itr
    if (this%topo_done) then
       ! Initialize gsh on level. shared_vtx is only seeded on level 0 once
       ! amge_mesh_set_shared_from_dofmap is wired into the hierarchy build;
       ! until then, fall back to "nothing is shared" (correct for a
       ! single-rank run, and for any level whose mesh doesn't cross ranks).
       if (allocated(this%mmsh%shared_vtx)) then
          call this%gsh%init(this%nelm(), &
               this%elm_vtx_ptr, this%elm_vtx_idx, &
               this%mmsh%vert_id, this%mmsh%shared_vtx)
       else
          allocate(shared_vtx(this%mmsh%n_verts))
          shared_vtx = .false.
          call this%gsh%init(this%nelm(), &
               this%elm_vtx_ptr, this%elm_vtx_idx, &
               this%mmsh%vert_id, shared_vtx)
       end if
       this%gsh_ready = .true.
       ! Allocate workspace vectors
       call this%new_vec(this%x)
       call this%new_vec(this%b)
       call this%new_vec(this%r)
       ! Allocate and fill multiplicity
       call amge_valence(this)
       ! Fill l1 diagonal
       allocate(this%dl1(this%mmsh%n_verts))
       call amge_setup_l1_diag(this, this%dl1)
    else
       call neko_error("AMGe level needs mesh info before init")
    end if
  end subroutine amge_level_data_init

  !> Calculate the number of dofs on element e
  pure function amge_level_ndof_el(this, e) result(n)
    class(amge_level_t), intent(in) :: this
    integer(i4), intent(in) :: e
    integer(i4) :: n
    n = this%elm_vtx_ptr(e + 1) - this%elm_vtx_ptr(e)
  end function amge_level_ndof_el

  !> Return number of elements on level
  pure function amge_level_nelm(this) result(n)
    class(amge_level_t), intent(in) :: this
    integer(i4) :: n
    n = this%mmsh%n_elem
  end function amge_level_nelm

  !> Allocate a cell-wise vector shaped for this level (x sized to the
  !! total duplicated dof count = elm_vtx_ptr(nelm+1)).
  subroutine amge_level_new_vec(this, v)
    class(amge_level_t), intent(in) :: this
    type(amge_vec_t), intent(inout) :: v
    call v%free()
    v%n_melm = this%nelm()
    v%n_dofs = this%elm_vtx_ptr(this%nelm() + 1)
    v%n_dofs_unique = this%mmsh%n_verts
    allocate(v%x(v%n_dofs))
    v%x = 0.0_rp
  end subroutine amge_level_new_vec

  !> Valence (duplication count) per unique dof: diag(G^T G). v = 1 ; gs ;
  !! mult = 1/v -- gsh%op already sums every copy of a dof, same-rank AND
  !! cross-rank, so a single gs of an all-ones vector gives the exact
  !! global valence at every duplicated slot directly (see amge_gs.f90's
  !! header comment for the rationale).
  !!
  !! winv (on lvl%tr, the transfer FROM this level's own finer neighbor)
  !! is a DIFFERENT quantity -- macroelement-boundary-membership count of
  !! a FINE (pre-coarsening) vertex, not this level's own element-duplication
  !! count -- and is filled in coarsen_level_3d, which is the only place
  !! that has both the fine level's fdofs lists and its own gs handle to
  !! correct them cross-rank; see amge_gs_correct_shared_count.
  subroutine amge_valence(lvl)
    type(amge_level_t), intent(inout) :: lvl
    integer(i4) :: n_dofs
    n_dofs = lvl%elm_vtx_ptr(lvl%nelm()+1)
    if (allocated(lvl%mult)) deallocate(lvl%mult)
    allocate(lvl%mult(n_dofs))
    lvl%mult = 1.0_rp
    call lvl%gsh%op(lvl%mult)
    lvl%mult = 1.0_rp / lvl%mult
  end subroutine amge_valence

  ! ================== smoother ==================

  !> Element-local l1 diagonal d_i = sum_{e ni i} sum_j |(A_e)_ij|,
  !! returned as a unique (assembled) vector. Dominates the exact
  !! assembled l1 diagonal entrywise, preserving the l1-Jacobi guarantee.
  subroutine amge_setup_l1_diag(lvl, d)
    type(amge_level_t), intent(inout) :: lvl
    real(rp), intent(inout) :: d(:)             !< length n_verts
    integer(i4) :: e, i, off, n
    d = 0.0_rp
    do e = 1, lvl%nelm()
       off = lvl%elm_vtx_ptr(e)
       n = lvl%ndof_el(e)
       do i = 1, n
          d(lvl%elm_vtx_idx(off + i)) = d(lvl%elm_vtx_idx(off + i)) &
               + sum(abs(lvl%AM(e)%x(i, :)))
       end do
    end do
  end subroutine amge_setup_l1_diag

  ! ================== housekeeping ==================

  subroutine amge_vec_free(this)
    class(amge_vec_t), intent(inout) :: this
    if (allocated(this%x)) deallocate(this%x)
    this%n_melm = 0; this%n_dofs = 0; this%n_dofs_unique = 0
  end subroutine amge_vec_free

  subroutine tr_map_free(this)
    class(tr_map_t), intent(inout) :: this
    if (allocated(this%fdofs)) deallocate(this%fdofs)
    if (allocated(this%cdofs)) deallocate(this%cdofs)
    call this%PiTilde%free()
    if (allocated(this%child_idx)) deallocate(this%child_idx)
    if (allocated(this%chloc_ptr)) deallocate(this%chloc_ptr)
    if (allocated(this%chloc_idx)) deallocate(this%chloc_idx)
  end subroutine tr_map_free

  subroutine amge_level_transfer_free(this)
    class(amge_level_transfer_t), intent(inout) :: this
    integer(i4) :: m
    if (allocated(this%maps)) then
       do m = 1, size(this%maps)
          call this%maps(m)%free()
       end do
       deallocate(this%maps)
    end if
    if (allocated(this%winv)) deallocate(this%winv)
    this%n_fine = 0; this%n_coarse = 0; this%n_melm = 0
    this%n_fine_unique = 0; this%n_coarse_unique = 0
  end subroutine amge_level_transfer_free

  subroutine amge_level_free(this)
    class(amge_level_t), intent(inout) :: this
    integer(i4) :: e
    call this%mmsh%free()
    if (allocated(this%elm_vtx_ptr)) deallocate(this%elm_vtx_ptr)
    if (allocated(this%elm_vtx_idx)) deallocate(this%elm_vtx_idx)
    if (allocated(this%AM)) then
       do e = 1, size(this%AM)
          call this%AM(e)%free()
       end do
       deallocate(this%AM)
    end if
    call this%tr%free()
  end subroutine amge_level_free

end module amge_level
