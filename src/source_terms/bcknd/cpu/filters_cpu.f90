! Copyright (c) 2024, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> CPU implementations of the filter functions.
module filters_cpu
  use num_types, only: rp
  implicit none
  private
  public :: smooth_step_cpu, step_function_cpu, permeability_cpu

contains

  !> @brief Apply a smooth step function to a scalar.
  subroutine smooth_step_cpu(x, edge0, edge1, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: x
    real(kind=rp), intent(in) :: edge0, edge1

    x = smooth_step_kernel(x, edge0, edge1)

  end subroutine smooth_step_cpu

  !> @brief Apply a step function to a scalar.
  subroutine step_function_cpu(x, x_step, value0, value1, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: x
    real(kind=rp), intent(in) :: x_step, value0, value1

    x = step_function_kernel(x, x_step, value0, value1)

  end subroutine step_function_cpu

  !> @brief Apply a permeability function to a scalar.
  subroutine permeability_cpu(x, k_0, k_1, q, n)
    integer, intent(in) :: n
    real(kind=rp), dimension(n), intent(inout) :: x
    real(kind=rp), intent(in) :: k_0, k_1, q

    x = permeability_kernel(x, k_0, k_1, q)

  end subroutine permeability_cpu


  ! ========================================================================== !
  ! Internal functions and subroutines
  ! ========================================================================== !

  !> @brief Apply a smooth step function to a scalar.
  elemental function smooth_step_kernel(x, edge0, edge1) result(res)
    real(kind=rp), intent(in) :: x
    real(kind=rp), intent(in) :: edge0, edge1
    real(kind=rp) :: res, t

    t = clamp_kernel((x - edge0) / (edge1 - edge0), 0.0_rp, 1.0_rp)

    res = t**3 * (t * (6.0_rp * t - 15.0_rp) + 10.0_rp)

  end function smooth_step_kernel

  !> @brief Clamp a value between two limits.
  elemental function clamp_kernel(x, lowerlimit, upperlimit) result(res)
    real(kind=rp), intent(in) :: x
    real(kind=rp), intent(in) :: lowerlimit, upperlimit
    real(kind=rp) :: res

    res = max(lowerlimit, min(upperlimit, x))
  end function clamp_kernel

  !> @brief Apply a step function to a scalar.
  elemental function step_function_kernel(x, x_step, value0, value1) result(res)
    real(kind=rp), intent(in) :: x, x_step, value0, value1
    real(kind=rp) :: res

    res = merge(value0, value1, x > x_step)

  end function step_function_kernel

  !> @brief Apply a permeability function to a scalar.
  elemental function permeability_kernel(x, k_0, k_1, q) result(perm)
    real(kind=rp), intent(in) :: x, k_0, k_1, q
    real(kind=rp) :: perm

    perm = k_0 + (k_1 - k_0) * x * (q + 1.0_rp) / (q + x)

  end function permeability_kernel


end module filters_cpu
