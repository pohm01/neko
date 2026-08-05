! Copyright (c) 2026, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the conditions of the
! BSD 3-Clause License are met (see Neko's COPYING file).
!
!> Chebyshev smoother for one AMGe level, mirroring tree_amg_smoother's
!! amg_cheby_t (same power-iteration eigenvalue estimate, same 3-term
!! Chebyshev recurrence from Saad's textbook).
!!
!! The one structural difference from tree_amg's version: EVERY AMGe
!! level -- not just the finest -- uses cell-wise (duplicated) per-
!! macroelement storage (amge_vec_t), so there is no "already unique,
!! unweighted" branch the way tree_amg has for its aggregated coarse
!! levels. The "matvec" here is always amge_apply (unassembled local
!! action) followed by a gather-scatter (lvl%gsh%op) to sum shared
!! copies, and inner products are always weighted by lvl%mult (1/
!! duplication count), the AMGe analogue of Neko's coef%mult.
module amge_smoother
  use amge_level, only : amge_level_t, amge_vec_t, amge_apply
  use num_types, only : i4, rp
  use math, only : glsc3, sub2, cmult2, copy
  use logger, only : neko_log, LOG_SIZE
  implicit none
  private

  !> Chebyshev iteration smoother for one AMGe level.
  type, public :: amge_cheby_t
     type(amge_vec_t) :: d, w, r
     real(kind=rp) :: tha, dlt
     integer :: lvl
     integer :: power_its = 250
     integer :: max_iter = 10
     logical :: recompute_eigs = .true.
   contains
     procedure, pass(this) :: init => amge_cheby_init
     procedure, pass(this) :: solve => amge_cheby_solve
     procedure, pass(this) :: comp_eig => amge_cheby_power
     procedure, pass(this) :: free => amge_cheby_free
  end type amge_cheby_t

contains

  !> Initialize the Chebyshev smoother for one AMGe level.
  !! @param lvl_obj the AMGe level this smoother will act on (used only
  !!                 to size the work vectors via new_vec; not stored)
  !! @param lvl      level index, for logging only
  !! @param max_iter Chebyshev degree (number of iterations per smooth)
  subroutine amge_cheby_init(this, lvl_obj, lvl, max_iter)
    class(amge_cheby_t), intent(inout) :: this
    type(amge_level_t), intent(in) :: lvl_obj
    integer, intent(in) :: lvl
    integer, intent(in) :: max_iter

    call lvl_obj%new_vec(this%d)
    call lvl_obj%new_vec(this%w)
    call lvl_obj%new_vec(this%r)
    this%lvl = lvl
    this%max_iter = max_iter
    this%recompute_eigs = .true.

    call amge_smoo_monitor(lvl, this)
  end subroutine amge_cheby_init

  !> Free cheby data
  subroutine amge_cheby_free(this)
    class(amge_cheby_t), intent(inout) :: this
    call this%d%free()
    call this%w%free()
    call this%r%free()
  end subroutine amge_cheby_free

  !> The assembled action of this level's operator: amge_apply's
  !! unassembled local action, then a gather-scatter to sum every
  !! duplicated (same-macroelement AND cross-rank) copy of a dof -- the
  !! AMGe analogue of tamg's amg%matvec.
  subroutine amge_cheby_matvec(lvl_obj, xin, yout)
    type(amge_level_t), intent(inout) :: lvl_obj
    type(amge_vec_t), intent(in) :: xin
    type(amge_vec_t), intent(inout) :: yout
    call amge_apply(lvl_obj, xin, yout)
    call lvl_obj%gsh%op(yout%x)
  end subroutine amge_cheby_matvec

  !> Power method to approximate the largest eigenvalue of this level's
  !! assembled operator, mirroring tree_amg_smoother's amg_cheby_power.
  !! @param lvl_obj the AMGe level to estimate eigenvalues for
  subroutine amge_cheby_power(this, lvl_obj)
    class(amge_cheby_t), intent(inout) :: this
    type(amge_level_t), intent(inout) :: lvl_obj
    real(kind=rp) :: lam, b, a, rn
    real(kind=rp), parameter :: boost = 1.1_rp
    real(kind=rp), parameter :: lam_factor = 30.0_rp
    real(kind=rp) :: wtw, dtw, dtd
    integer, allocatable :: fixed_seed(:), saved_seed(:)
    integer :: i, rnd_n, n

    associate(w => this%w, d => this%d)
      n = d%n_dofs

      ! Save current random seed and set a fixed seed (reproducible
      ! eigenvalue estimate across runs)
      call random_seed(size = rnd_n)
      allocate(saved_seed(rnd_n), fixed_seed(rnd_n))
      fixed_seed = 3901
      call random_seed(get = saved_seed)
      call random_seed(put = fixed_seed)

      do i = 1, n
         call random_number(rn)
         d%x(i) = rn + 10.0_rp
      end do

      ! Restore saved random seed
      call random_seed(put = saved_seed)

      ! make the random initial vector consistent across duplicated
      ! copies of shared dofs before iterating (same role as tamg's
      ! gs_h%op on its level-0 initial vector)
      call lvl_obj%gsh%op(d%x)

      ! Power method to get lambda max
      do i = 1, this%power_its
         call amge_cheby_matvec(lvl_obj, d, w)
         wtw = glsc3(w%x, lvl_obj%mult, w%x, n)
         call cmult2(d%x, w%x, 1.0_rp / sqrt(wtw), n)
      end do

      call amge_cheby_matvec(lvl_obj, d, w)
      dtw = glsc3(d%x, lvl_obj%mult, w%x, n)
      dtd = glsc3(d%x, lvl_obj%mult, d%x, n)
      lam = dtw / dtd
      b = lam * boost
      a = lam / lam_factor
      this%tha = (b + a) / 2.0_rp
      this%dlt = (b - a) / 2.0_rp

      this%recompute_eigs = .false.
      call amge_cheby_monitor(this%lvl, lam)
    end associate
  end subroutine amge_cheby_power

  !> Chebyshev smoother, from Saad's iterative methods textbook.
  !! @param lvl_obj  the AMGe level to smooth on
  !! @param x        the iterate to be updated in place
  !! @param f        the (already gather-scattered) cell-wise rhs
  subroutine amge_cheby_solve(this, lvl_obj, x, f, zero_init)
    class(amge_cheby_t), intent(inout) :: this
    type(amge_level_t), intent(inout) :: lvl_obj
    type(amge_vec_t), intent(inout) :: x
    type(amge_vec_t), intent(in) :: f
    logical, optional, intent(in) :: zero_init
    integer :: iter, max_iter, i, n
    real(kind=rp) :: rhok, rhokp1, s1, thet, delt, tmp1, tmp2
    logical :: zero_initial_guess

    if (this%recompute_eigs) call this%comp_eig(lvl_obj)

    if (present(zero_init)) then
       zero_initial_guess = zero_init
    else
       zero_initial_guess = .false.
    end if
    max_iter = this%max_iter
    n = x%n_dofs

    associate(w => this%w, r => this%r, d => this%d)
      call copy(r%x, f%x, n)
      if (.not. zero_initial_guess) then
         call amge_cheby_matvec(lvl_obj, x, w)
         call sub2(r%x, w%x, n)
      end if

      thet = this%tha
      delt = this%dlt
      s1 = thet / delt
      rhok = 1.0_rp / s1

      ! First iteration
      !OCL NORECURRENCE, NOVREC, NOALIAS
      !DIR$ CONCURRENT
      !DIR$ IVDEP
      !GCC$ ivdep
      !$omp parallel do
      do i = 1, n
         d%x(i) = 1.0_rp / thet * r%x(i)
         x%x(i) = x%x(i) + d%x(i)
      end do
      !$omp end parallel do

      ! Rest of iterations
      do iter = 2, max_iter
         call amge_cheby_matvec(lvl_obj, d, w)

         rhokp1 = 1.0_rp / (2.0_rp * s1 - rhok)
         tmp1 = rhokp1 * rhok
         tmp2 = 2.0_rp * rhokp1 / delt
         rhok = rhokp1

         !$omp parallel private(i)
         !OCL NORECURRENCE, NOVREC, NOALIAS
         !DIR$ CONCURRENT
         !DIR$ IVDEP
         !GCC$ ivdep
         !$omp do
         do i = 1, n
            r%x(i) = r%x(i) - w%x(i)
            d%x(i) = tmp1 * d%x(i) + tmp2 * r%x(i)
            x%x(i) = x%x(i) + d%x(i)
         end do
         !$omp end do
         !$omp end parallel
      end do
    end associate
  end subroutine amge_cheby_solve

  subroutine amge_smoo_monitor(lvl, smoo)
    integer, intent(in) :: lvl
    class(amge_cheby_t), intent(in) :: smoo
    character(len=LOG_SIZE) :: log_buf

    write(log_buf, '(A8,I2,A28)') '-- level', lvl, '-- init smoother: Chebyshev'
    call neko_log%message(log_buf)
    write(log_buf, '(A22,I6)') 'Iterations:', smoo%max_iter
    call neko_log%message(log_buf)
  end subroutine amge_smoo_monitor

  subroutine amge_cheby_monitor(lvl, lam)
    integer, intent(in) :: lvl
    real(kind=rp), intent(in) :: lam
    character(len=LOG_SIZE) :: log_buf

    write(log_buf, '(A13,I2,A29,F12.3)') '-- AMGe level', lvl, &
         '-- Chebyshev approx. max eig', lam
    call neko_log%message(log_buf)
  end subroutine amge_cheby_monitor

end module amge_smoother
