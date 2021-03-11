module operators
  use num_types  
  use opr_cpu
  use opr_sx
  use space  
  use coefs
  use field
  use math
  implicit none
  
contains
  
  subroutine dudxyz (du,u,dr,ds,dt,coef)
!--------------------------------------------------------------
!
!     du  - dU/dx or dU/dy or dU/dz
!     u   - a field variable defined on mesh 1
!     dr  - dr/dx or dr/dy or dr/dz  
!     ds  - ds/dx or ds/dy or ds/dz
!     dt  - dt/dx or dt/dy or dt/dz
    type(coef_t), intent(in), target :: coef
    real(kind=dp), dimension(coef%Xh%lx,coef%Xh%ly,coef%Xh%lz,coef%msh%nelv), intent(inout) ::  du
    real(kind=dp), dimension(coef%Xh%lx,coef%Xh%ly,coef%Xh%lz,coef%msh%nelv), intent(inout) ::  u, dr, ds, dt

    if (1 .eq. 1) then !! TODO check this in neko_config
       call opr_cpu_dudxyz(du, u, dr, ds, dt, coef)
    else
       call opr_sx_dudxyz(du, u, dr, ds, dt, coef)
    end if
    
  end subroutine dudxyz

  !> Equals wgradm1 in nek5000. Gradient of velocity vectors.
  subroutine opgrad(ux,uy,uz,u,coef) ! weak form of grad 

  !Compute gradient of T -- mesh 1 to mesh 1 (vel. to vel.)

    type(coef_t), intent(in) :: coef  
    real(kind=dp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: ux
    real(kind=dp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: uy
    real(kind=dp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: uz
    real(kind=dp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: u

    if (1 .eq. 1) then !! TODO check this in neko_config
       call opr_cpu_opgrad(ux, uy, uz, u, coef)
    else
       call opr_sx_opgrad(ux, uy, uz, u, coef)
    end if
    
  end subroutine opgrad
  
  !> Othogonalize with regard to vector (1,1,1,1,1,1...,1)^T.
  subroutine ortho(x,n ,glb_n)
    integer, intent(in) :: n
    integer, intent(in) :: glb_n
    real(kind=dp), dimension(n), intent(inout) :: x
    real(kind=dp) :: rlam

    rlam = glsum(x,n)/glb_n
    call cadd(x,-rlam,n)

  end subroutine ortho
  
  !> Compute DT*X (entire field)
  !> This needs to be revised... the loop over n1,n2 is probably unesccssary
  subroutine cdtp (dtx,x,dr,ds,dt, coef)
    type(coef_t), intent(in) :: coef
    real(kind=dp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: dtx
    real(kind=dp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: x
    real(kind=dp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: dr
    real(kind=dp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: ds
    real(kind=dp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: dt

    if (1 .eq. 1) then !! TODO check this in neko_config
       call opr_cpu_cdtp(dtx, x, dr, ds, dt, coef)
    else
       call opr_sx_cdtp(dtx, x, dr, ds, dt, coef)
    end if
    
  end subroutine cdtp
   
  subroutine conv1(du,u, vx, vy, vz, Xh, coef, nelv, gdim)  ! used to be conv1n
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    integer, intent(in) :: nelv, gdim
    real(kind=dp), intent(inout) ::  du(Xh%lxyz,nelv)
    real(kind=dp), intent(inout), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  u
    real(kind=dp), intent(inout), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  vx
    real(kind=dp), intent(inout), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  vy
    real(kind=dp), intent(inout), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  vz

    call opr_cpu_conv1(du, u, vx, vy, vz, Xh, coef, nelv, gdim)

  end subroutine conv1

  subroutine curl(w1, w2, w3, u1, u2, u3, work1, work2, c_Xh)
    type(field_t), intent(inout) :: w1
    type(field_t), intent(inout) :: w2
    type(field_t), intent(inout) :: w3
    type(field_t), intent(inout) :: u1
    type(field_t), intent(inout) :: u2
    type(field_t), intent(inout) :: u3
    type(field_t), intent(inout) :: work1
    type(field_t), intent(inout) :: work2
    type(coef_t), intent(inout)  :: c_Xh

    call opr_cpu_curl(w1, w2, w3, u1, u2, u3, work1, work2, c_Xh)

  end subroutine curl

  function cfl(dt, u, v, w, Xh, coef, nelv, gdim)
    type(space_t) :: Xh
    type(coef_t) :: coef
    integer :: nelv, gdim
    real(kind=dp) :: dt
    real(kind=dp) ::  du(Xh%lxyz,nelv)
    real(kind=dp), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  u, v, w
    real(kind=dp) :: cflr, cfls, cflt, cflm, cfl_temp(1)
    real(kind=dp) :: ur, us, ut
    real(kind=dp) :: cfl
    integer :: i,j,k,e
    cfl_temp(1) = 0d0
    if (gdim .eq. 3) then
       do e=1,nelv
          do k=1,Xh%lz
          do j=1,Xh%ly
          do i=1,Xh%lx
             ur = ( u(i,j,k,e)*coef%drdx(i,j,k,e) &
                +   v(i,j,k,e)*coef%drdy(i,j,k,e) &
                +   w(i,j,k,e)*coef%drdz(i,j,k,e) ) * coef%jacinv(i,j,k,e)
             us = ( u(i,j,k,e)*coef%dsdx(i,j,k,e) &
                +   v(i,j,k,e)*coef%dsdy(i,j,k,e) &
                +   w(i,j,k,e)*coef%dsdz(i,j,k,e) ) * coef%jacinv(i,j,k,e)
             ut = ( u(i,j,k,e)*coef%dtdx(i,j,k,e) &
                +   v(i,j,k,e)*coef%dtdy(i,j,k,e) &
                +   w(i,j,k,e)*coef%dtdz(i,j,k,e) ) * coef%jacinv(i,j,k,e)
 
             cflr = abs(dt*ur*Xh%dr_inv(i))
             cfls = abs(dt*us*Xh%ds_inv(j))
             cflt = abs(dt*ut*Xh%dt_inv(k))
 
             cflm = cflr + cfls + cflt
             cfl_temp(1)  = max(cfl_temp(1),cflm)
          enddo
          enddo
          enddo
       enddo
    else
       do e=1,nelv
          do j=1,Xh%ly
          do i=1,Xh%lx
             ur = ( u(i,j,1,e)*coef%drdx(i,j,1,e) &
                +   v(i,j,1,e)*coef%drdy(i,j,1,e) ) * coef%jacinv(i,j,1,e)
             us = ( u(i,j,1,e)*coef%dsdx(i,j,1,e) &
                +   v(i,j,1,e)*coef%dsdy(i,j,1,e) ) * coef%jacinv(i,j,1,e)

             cflr = abs(dt*ur*Xh%dr_inv(i))
             cfls = abs(dt*us*Xh%ds_inv(j))

             cflm = cflr + cfls
             cfl_temp(1)  = max(cfl_temp(1),cflm)

          enddo
          enddo
       enddo
    endif
    cfl = glmax(cfl_temp,1)
  end function cfl
  
end module operators
