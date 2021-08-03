!> Fluid formulations
module fluid_method
  use gather_scatter
  use neko_config
  use parameters
  use num_types
  use source
  use field
  use space
  use dofmap
  use krylov
  use coefs
  use wall
  use inflow
  use usr_inflow
  use dirichlet
  use symmetry
  use cg
  use bicgstab
  use bc
  use jacobi
  use sx_jacobi
  use gmres
  use mesh
  use math
  use abbdf
  use mathops
  use operators
  use hsmg
  use log
  implicit none
  
  !> Base type of all fluid formulations
  type, abstract :: fluid_scheme_t
     type(field_t) :: u         !< x-component of Velocity
     type(field_t) :: v         !< y-component of Velocity
     type(field_t) :: w         !< z-component of Velocity
     type(field_t) :: p         !< Pressure
     type(space_t) :: Xh        !< Function space \f$ X_h \f$
     type(dofmap_t) :: dm_Xh    !< Dofmap associated with \f$ X_h \f$
     type(gs_t) :: gs_Xh        !< Gather-scatter associated with \f$ X_h \f$
     type(coef_t) :: c_Xh       !< Coefficients associated with \f$ X_h \f$
     type(source_t) :: f_Xh     !< Source term associated with \f$ X_h \f$
     class(ksp_t), allocatable  :: ksp_vel     !< Krylov solver for velocity
     class(ksp_t), allocatable  :: ksp_prs     !< Krylov solver for pressure
     class(pc_t), allocatable :: pc_vel        !< Velocity Preconditioner
     class(pc_t), allocatable :: pc_prs        !< Velocity Preconditioner
     type(no_slip_wall_t) :: bc_wall           !< No-slip wall for velocity
     class(inflow_t), allocatable :: bc_inflow !< Dirichlet inflow for velocity
     type(dirichlet_t) :: bc_prs               !< Dirichlet pressure condition
     type(symmetry_t) :: bc_sym                !< Symmetry plane for velocity
     type(bc_list_t) :: bclst_vel              !< List of velocity conditions
     type(bc_list_t) :: bclst_prs              !< List of pressure conditions
     type(field_t) :: bdry                     !< Boundary markings
     type(param_t), pointer :: params          !< Parameters          
     type(mesh_t), pointer :: msh => null()    !< Mesh
   contains
     procedure, pass(this) :: fluid_scheme_init_all
     procedure, pass(this) :: fluid_scheme_init_uvw
     procedure, pass(this) :: scheme_free => fluid_scheme_free
     procedure, pass(this) :: validate => fluid_scheme_validate
     procedure, pass(this) :: bc_apply_vel => fluid_scheme_bc_apply_vel
     procedure, pass(this) :: bc_apply_prs => fluid_scheme_bc_apply_prs
     procedure, pass(this) :: set_usr_inflow => fluid_scheme_set_usr_inflow
     procedure, pass(this) :: compute_cfl => fluid_compute_cfl
     procedure(fluid_method_init), pass(this), deferred :: init
     procedure(fluid_method_free), pass(this), deferred :: free
     procedure(fluid_method_step), pass(this), deferred :: step
     generic :: scheme_init => fluid_scheme_init_all, fluid_scheme_init_uvw
  end type fluid_scheme_t

  !> Abstract interface to initialize a fluid formulation
  abstract interface
     subroutine fluid_method_init(this, msh, lx, param)
       import fluid_scheme_t
       import param_t
       import mesh_t
       class(fluid_scheme_t), intent(inout) :: this
       type(mesh_t), intent(inout) :: msh       
       integer, intent(inout) :: lx
       type(param_t), intent(inout) :: param              
     end subroutine fluid_method_init
  end interface

  !> Abstract interface to dealocate a fluid formulation
  abstract interface
     subroutine fluid_method_free(this)
       import fluid_scheme_t
       class(fluid_scheme_t), intent(inout) :: this
     end subroutine fluid_method_free
  end interface
  
  !> Abstract interface to compute a time-step
  abstract interface
     subroutine fluid_method_step(this, t, tstep, ab_bdf)
       import fluid_scheme_t
       import abbdf_t
       import rp
       class(fluid_scheme_t), intent(inout) :: this
       real(kind=rp), intent(inout) :: t
       integer, intent(inout) :: tstep
       type(abbdf_t), intent(inout) :: ab_bdf
     end subroutine fluid_method_step
  end interface

contains

  !> Initialize common data for the current scheme
  subroutine fluid_scheme_init_common(this, msh, lx, params)
    class(fluid_scheme_t), intent(inout) :: this
    type(mesh_t), intent(inout), target :: msh
    integer, intent(inout) :: lx
    type(param_t), intent(inout), target :: params
    type(dirichlet_t) :: bdry_mask
    
    call neko_log%section('Fluid')
    call neko_log%message('Ksp vel. : ('// trim(params%ksp_vel) // &
         ', ' // trim(params%pc_vel) // ')')
    call neko_log%message('Ksp prs. : ('// trim(params%ksp_prs) // &
         ', ' // trim(params%pc_prs) // ')')
    
    if (msh%gdim .eq. 2) then
       call space_init(this%Xh, GLL, lx, lx)
    else
       call space_init(this%Xh, GLL, lx, lx, lx)
    end if

    this%dm_Xh = dofmap_t(msh, this%Xh)

    this%params => params

    this%msh => msh

    call gs_init(this%gs_Xh, this%dm_Xh)

    call coef_init(this%c_Xh, this%gs_Xh)

    call source_init(this%f_Xh, this%dm_Xh)

    !
    ! Setup velocity boundary conditions
    !
    call bc_list_init(this%bclst_vel)

    if (msh%sympln%size .gt. 0) then
       call this%bc_sym%init(this%dm_Xh)
       call this%bc_sym%mark_zone(msh%sympln)
       call this%bc_sym%finalize()
       call this%bc_sym%init_msk(this%c_Xh)    
       call bc_list_add(this%bclst_vel, this%bc_sym)
    end if

    if (msh%inlet%size .gt. 0) then

       if (trim(params%fluid_inflow) .eq. "default") then
          allocate(inflow_t::this%bc_inflow)
       else if (trim(params%fluid_inflow) .eq. "user") then
          allocate(usr_inflow_t::this%bc_inflow)
       else
          call neko_error('Invalid Inflow condition')
       end if
       
       call this%bc_inflow%init(this%dm_Xh)
       call this%bc_inflow%mark_zone(msh%inlet)
       call this%bc_inflow%finalize()
       call this%bc_inflow%set_inflow(params%uinf)
       call bc_list_add(this%bclst_vel, this%bc_inflow)

       if (trim(params%fluid_inflow) .eq. "user") then
          select type(bc_if => this%bc_inflow)
          type is(usr_inflow_t)
             call bc_if%set_coef(this%C_Xh)
          end select
       end if
    end if
    
    if (msh%wall%size .gt. 0 ) then
       call this%bc_wall%init(this%dm_Xh)
       call this%bc_wall%mark_zone(msh%wall)
       call this%bc_wall%finalize()
       call bc_list_add(this%bclst_vel, this%bc_wall)
    end if
       
    if (params%output_bdry) then

       call neko_log%message('Saving boundary markings')
       
       call field_init(this%bdry, this%dm_Xh, 'bdry')
       this%bdry = real(0d0,rp)
       
       call bdry_mask%init(this%dm_Xh)
       call bdry_mask%mark_zone(msh%wall)
       call bdry_mask%finalize()
       call bdry_mask%set_g(real(1d0,rp))
       call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%n_dofs)
       call bdry_mask%free()

       call bdry_mask%init(this%dm_Xh)
       call bdry_mask%mark_zone(msh%inlet)
       call bdry_mask%finalize()
       call bdry_mask%set_g(real(2d0,rp))
       call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%n_dofs)
       call bdry_mask%free()

       call bdry_mask%init(this%dm_Xh)
       call bdry_mask%mark_zone(msh%outlet)
       call bdry_mask%finalize()
       call bdry_mask%set_g(real(3d0,rp))
       call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%n_dofs)
       call bdry_mask%free()

       call bdry_mask%init(this%dm_Xh)
       call bdry_mask%mark_zone(msh%sympln)
       call bdry_mask%finalize()
       call bdry_mask%set_g(real(4d0,rp))
       call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%n_dofs)
       call bdry_mask%free()

       call bdry_mask%init(this%dm_Xh)
       call bdry_mask%mark_zone(msh%periodic)
       call bdry_mask%finalize()
       call bdry_mask%set_g(real(5d0,rp))
       call bdry_mask%apply_scalar(this%bdry%x, this%dm_Xh%n_dofs)
       call bdry_mask%free()
    end if

  end subroutine fluid_scheme_init_common

  !> Initialize all velocity related components of the current scheme
  subroutine fluid_scheme_init_uvw(this, msh, lx, params, kspv_init)
    class(fluid_scheme_t), intent(inout) :: this
    type(mesh_t), intent(inout) :: msh
    integer, intent(inout) :: lx
    type(param_t), intent(inout) :: params
    logical :: kspv_init

    call fluid_scheme_init_common(this, msh, lx, params)
    
    call field_init(this%u, this%dm_Xh, 'u')
    call field_init(this%v, this%dm_Xh, 'v')
    call field_init(this%w, this%dm_Xh, 'w')

    if (kspv_init) then
       call fluid_scheme_solver_factory(this%ksp_vel, this%dm_Xh%size(), &
            params%ksp_vel, params%abstol_vel)
       call fluid_scheme_precon_factory(this%pc_vel, this%ksp_vel, &
            this%c_Xh, this%dm_Xh, this%gs_Xh, this%bclst_vel, params%pc_vel)
    end if

    call neko_log%end_section()
  end subroutine fluid_scheme_init_uvw

  !> Initialize all components of the current scheme
  subroutine fluid_scheme_init_all(this, msh, lx, params, kspv_init, kspp_init)
    class(fluid_scheme_t), intent(inout) :: this
    type(mesh_t), intent(inout) :: msh
    integer, intent(inout) :: lx
    type(param_t), intent(inout) :: params
    logical :: kspv_init
    logical :: kspp_init

    call fluid_scheme_init_common(this, msh, lx, params)
    
    call field_init(this%u, this%dm_Xh, 'u')
    call field_init(this%v, this%dm_Xh, 'v')
    call field_init(this%w, this%dm_Xh, 'w')
    call field_init(this%p, this%dm_Xh, 'p')

    !
    ! Setup pressure boundary conditions
    !
    call bc_list_init(this%bclst_prs)
    if (msh%outlet%size .gt. 0) then
       call this%bc_prs%init(this%dm_Xh)
       call this%bc_prs%mark_zone(msh%outlet)
       call this%bc_prs%finalize()
       call this%bc_prs%set_g(real(0d0,rp))
       call bc_list_add(this%bclst_prs, this%bc_prs)
    end if

    if (kspv_init) then
       call fluid_scheme_solver_factory(this%ksp_vel, this%dm_Xh%size(), &
            params%ksp_vel, params%abstol_vel)
       call fluid_scheme_precon_factory(this%pc_vel, this%ksp_vel, &
            this%c_Xh, this%dm_Xh, this%gs_Xh, this%bclst_vel, params%pc_vel)
    end if

    if (kspp_init) then
       call fluid_scheme_solver_factory(this%ksp_prs, this%dm_Xh%size(), &
            params%ksp_prs, params%abstol_prs)
       call fluid_scheme_precon_factory(this%pc_prs, this%ksp_prs, &
            this%c_Xh, this%dm_Xh, this%gs_Xh, this%bclst_prs, params%pc_prs)
    end if


    call neko_log%end_section()
    
  end subroutine fluid_scheme_init_all

  !> Deallocate a fluid formulation
  subroutine fluid_scheme_free(this)
    class(fluid_scheme_t), intent(inout) :: this

    call field_free(this%u)
    call field_free(this%v)
    call field_free(this%w)
    call field_free(this%p)
    call field_free(this%bdry)

    if (allocated(this%bc_inflow)) then
       call this%bc_inflow%free()
    end if

    call this%bc_wall%free()
    call this%bc_sym%free()

    call space_free(this%Xh)    

    if (allocated(this%ksp_vel)) then
       call this%ksp_vel%free()
       deallocate(this%ksp_vel)
    end if

    if (allocated(this%ksp_prs)) then
       call this%ksp_prs%free()
       deallocate(this%ksp_prs)
    end if

    call gs_free(this%gs_Xh)

    call coef_free(this%c_Xh)

    call source_free(this%f_Xh)

    call bc_list_free(this%bclst_vel)

    nullify(this%params)
    
  end subroutine fluid_scheme_free

  !> Validate that all fields, solvers etc necessary for
  !! performing time-stepping are defined
  subroutine fluid_scheme_validate(this)
    class(fluid_scheme_t), intent(inout) :: this

    if ( (.not. allocated(this%u%x)) .or. &
         (.not. allocated(this%v%x)) .or. &
         (.not. allocated(this%w%x)) .or. &
         (.not. allocated(this%p%x))) then
       call neko_error('Fields are not allocated')
    end if

    if (.not. allocated(this%ksp_vel)) then
       call neko_error('No Krylov solver for velocity defined')
    end if
    
    if (.not. allocated(this%ksp_prs)) then
       call neko_error('No Krylov solver for pressure defined')
    end if

    if (.not. associated(this%f_Xh%eval)) then
       call neko_error('No source term defined')
    end if

    if (.not. associated(this%params)) then
       call neko_error('No parameters defined')
    end if

    select type(ip => this%bc_inflow)
    type is(usr_inflow_t)
       call ip%validate
    end select

  end subroutine fluid_scheme_validate

  !> Apply all boundary conditions defined for velocity
  !! @todo Why can't we call the interface here?
  subroutine fluid_scheme_bc_apply_vel(this)
    class(fluid_scheme_t), intent(inout) :: this
    call bc_list_apply_vector(this%bclst_vel,&
         this%u%x, this%v%x, this%w%x, this%dm_Xh%n_dofs)
  end subroutine fluid_scheme_bc_apply_vel
  
  !> Apply all boundary conditions defined for pressure
  !! @todo Why can't we call the interface here?
  subroutine fluid_scheme_bc_apply_prs(this)
    class(fluid_scheme_t), intent(inout) :: this
    call bc_list_apply_scalar(this%bclst_prs, this%p%x, this%p%dof%n_dofs)
  end subroutine fluid_scheme_bc_apply_prs
  
  !> Initialize a linear solver
  !! @note Currently only supporting Krylov solvers
  subroutine fluid_scheme_solver_factory(ksp, n, solver, abstol)
    class(ksp_t), allocatable, intent(inout) :: ksp
    integer, intent(in), value :: n
    character(len=20), intent(inout) :: solver
    real(kind=rp) :: abstol
    if (trim(solver) .eq. 'cg') then
       allocate(cg_t::ksp)
    else if (trim(solver) .eq. 'gmres') then
       allocate(gmres_t::ksp)
    else if (trim(solver) .eq. 'bicgstab') then
       allocate(bicgstab_t::ksp)
    else
       call neko_error('Unknown linear solver')
    end if

    select type(kp => ksp)
    type is(cg_t)
       call kp%init(n, abs_tol = abstol)
    type is(gmres_t)
       call kp%init(n, abs_tol = abstol)
    type is(bicgstab_t)
       call kp%init(n, abs_tol = abstol)
    end select
    
  end subroutine fluid_scheme_solver_factory

  !> Initialize a Krylov preconditioner
  subroutine fluid_scheme_precon_factory(pc, ksp, coef, dof, gs, bclst, pctype)
    class(pc_t), allocatable, intent(inout), target :: pc
    class(ksp_t), allocatable, intent(inout) :: ksp
    type(coef_t), intent(inout) :: coef
    type(dofmap_t), intent(inout) :: dof
    type(gs_t), intent(inout) :: gs
    type(bc_list_t), intent(inout) :: bclst
    character(len=20) :: pctype
    
    if (trim(pctype) .eq. 'jacobi') then
       if (NEKO_BCKND_SX .eq. 1) then
          allocate(sx_jacobi_t::pc)
       else
          allocate(jacobi_t::pc)
       end if
    else if (trim(pctype) .eq. 'hsmg') then
       allocate(hsmg_t::pc)
    else
       call neko_error('Unknown preconditioner')
    end if

    select type(pcp => pc)
    type is(jacobi_t)
       call pcp%init(coef, dof, gs)
    type is(hsmg_t)
       call pcp%init(dof%msh, dof%Xh, coef, dof, gs, bclst)
    end select

    call ksp%set_pc(pc)
    
  end subroutine fluid_scheme_precon_factory

  !> Initialize a user defined inflow condition
  subroutine fluid_scheme_set_usr_inflow(this, usr_eval)
    class(fluid_scheme_t), intent(inout) :: this
    procedure(usr_inflow_eval) :: usr_eval

    if (this%msh%inlet%size .gt. 0) then
       select type(bc_if => this%bc_inflow)
       type is(usr_inflow_t)
          call bc_if%set_eval(usr_eval)
       class default
          call neko_error("Not a user defined inflow condition")
       end select
    end if
    
  end subroutine fluid_scheme_set_usr_inflow

  !> Compute CFL
  function fluid_compute_cfl(this, dt) result(c)
    class(fluid_scheme_t), intent(in) :: this
    real(kind=rp), intent(in) :: dt
    real(kind=rp) :: c

    c = cfl(dt, this%u%x, this%v%x, this%w%x, &
         this%Xh, this%c_Xh, this%msh%nelv, this%msh%gdim)
    
  end function fluid_compute_cfl
     
end module fluid_method
