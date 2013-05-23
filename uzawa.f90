!=================================================================================================
subroutine uzawa(it)
  use module_grid
  use module_BC
  use module_flow
  use module_tmpvar
  implicit none

  integer it
  du = 0d0; dv = 0d0; dw = 0d0
  fx = 0d0; fy = 0d0; fz = 0d0;

  do it=1,maxit
     call volumeForce(rho,rho1,rho2,dpdx,dpdy,dpdz,BuoyancyCase,fx,fy,fz,gx,gy,gz,du,dv,dw, &
          rho_ave)
     if(dosolids) then
        du = du*umask; dv = dv*vmask; dw = dw*wmask
     endif
     do k=ks,ke;  do j=js,je; do i=is,ieu;    ! pressure term in x-momentum balance
        du(i,j,k) = du(i,j,k)-(2.0/dxh(i))*(p(i+1,j,k)-p(i,j,k))/(rho(i+1,j,k)+rho(i,j,k))
     enddo; enddo; enddo
     
     do k=ks,ke;  do j=js,jev; do i=is,ie;    ! pressure term in y-momentum balance
        dv(i,j,k) = dv(i,j,k)-(2.0/dyh(j))*(p(i,j+1,k)-p(i,j,k))/(rho(i,j+1,k)+rho(i,j,k))
     enddo; enddo; enddo
      
     do k=ks,kew;  do j=js,je; do i=is,ie;   ! pressure term in z-momentum balance
        dw(i,j,k) = dw(i,j,k)-(2.0/dzh(k))*(p(i,j,k+1)-p(i,j,k))/(rho(i,j,k+1)+rho(i,j,k))


     call SetupUvel_st(u,du,rho,mu,rho1,mu1,A)
     call Linearsolver_one_it(A,u) 

     call SetupVvel_st(v,dv,rho,mu,rho1,mu1,A)
     call Linearsolver_one_it(A,v) 

     call SetupWvel_st(w,dw,rho,mu,rho1,mu1,A)
     call Linearsolver_one_it(A,w) 

     call SetVelocityBC(u,v,w,umask,vmask,wmask)
     !-----------------------------------------PROJECTION STEP-----------------------------------------
     
     call SetupPoisson_st(u,v,w,umask,vmask,wmask,rho,A)
     tmp=0.d0
     call LinearSolver_one_it(A,tmp)
           
     do k=ks,ke;  do j=js,je; do i=is,ieu;    ! CORRECT THE u-velocity 
        u(i,j,k)=u(i,j,k)-dt*(2.0/dxh(i))*(p(i+1,j,k)-p(i,j,k))/(rho(i+1,j,k)+rho(i,j,k))
     enddo; enddo; enddo
     
     do k=ks,ke;  do j=js,jev; do i=is,ie;    ! CORRECT THE v-velocity
        v(i,j,k)=v(i,j,k)-dt*(2.0/dyh(j))*(p(i,j+1,k)-p(i,j,k))/(rho(i,j+1,k)+rho(i,j,k))
     enddo; enddo; enddo
      
     do k=ks,kew;  do j=js,je; do i=is,ie;   ! CORRECT THE w-velocity
        w(i,j,k)=w(i,j,k)-dt*(2.0/dzh(k))*(p(i,j,k+1)-p(i,j,k))/(rho(i,j,k+1)+rho(i,j,k))
     enddo; enddo; enddo
     call calcresidual2(u,v,w,umask,vmask,wmask,mu,rho,gx,gy,gz)
  enddo
  return
  end subroutine uzawa
!=================================================================================================
!=================================================================================================
! Solves the following linear equiation:
! A7*Uijk = umask*(A1*Ui-1jk + A2*Ui+1jk + A3*Uij-1k + 
!           A4*Uij+1k + A5*Uijk-1 + A6*Uijk+1 + A8)
!-------------------------------------------------------------------------------------------------
subroutine LinearSolver_one_it(A,u,umask,beta)
  use module_grid
  use module_BC
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: u
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: umask
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(in) :: A
  real(8), intent(in) :: beta
  integer :: i,j,k, ierr
  integer :: req(12),sta(MPI_STATUS_SIZE,12)
!-------------------------------------- ONE ITERATION of IMPLICIT DIFFUSION --------------------------

    do k=ks,ke; do j=js,je; do i=is,ie
      u(i,j,k)=umask(i,j,k)*((1d0-beta)*u(i,j,k)+beta* 1d0/A(i,j,k,7)*(&
        A(i,j,k,1) * u(i-1,j,k) + A(i,j,k,2) * u(i+1,j,k) +            &
        A(i,j,k,3) * u(i,j-1,k) + A(i,j,k,4) * u(i,j+1,k) +            &
        A(i,j,k,5) * u(i,j,k-1) + A(i,j,k,6) * u(i,j,k+1) + A(i,j,k,8)))
    enddo; enddo; enddo

    call ghost_x(u,1,req( 1: 4)); call ghost_y(u,1,req( 5: 8)); call ghost_z(u,1,req( 9:12))
    call MPI_WAITALL(12,req,sta,ierr)
end subroutine LinearSolver_one_it
!=================================================================================================
