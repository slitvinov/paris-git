!=================================================================================================
!=================================================================================================
! Paris-0.1
!
! authors Stephane Zaleski (zaleski@dalembert.upmc.fr) 
!         Ruben Scardovelli 
!         Daniel Fuster
!
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation; either version 2 of the
! License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
! 02111-1307, USA.  
!=================================================================================================
! module_VOF: Contains definition of variables for the Volume of Fluid interface tracking.
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------------------------
subroutine swp(us,c,f,d,vof1,vof2,vof3)
  use module_vof
  implicit none
  real (8)  , dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: us
  integer, intent(in) :: d
  real (8)  , dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: c,vof1,vof2,vof3
  integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: f
  if (VOF_advect=='Dick_Yue') then  ! Yue-Weymouth = Eulerian Implicit + central cell stuff
     call swpr(us,c,f,d,vof1,vof2,vof3)  
  elseif (VOF_advect=='CIAM') then  ! CIAM == Lagrangian Explicit
     call swpz(us,c,f,d,vof1,vof2,vof3)
  else
     call pariserror("*** unknown vof scheme")
  endif

  call get_flags_and_clip()
  call do_all_ghost(c)
  call do_all_ighost(f)

end subroutine swp
!
subroutine swpmom(us,c,f,d,mom1,mom2,mom3,mom)
  use module_vof
  implicit none
  integer, intent(in) :: d
  real(8)  , dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: us,mom
  real(8)  , dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: c
  real(8)  , dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: mom1,mom2,mom3
  integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: f

  if (VOF_advect=='Dick_Yue') then  ! Yue-Weymouth = Eulerian Implicit + central cell stuff
     call pariserror("***Dick_Yue + conserving momentum not implemented yet")
  elseif (VOF_advect=='CIAM') then  ! CIAM == Lagrangian Explicit
     call swpzmom(us,c,f,d,mom1,mom2,mom3,mom)
  else
     call pariserror("*** unknown vof scheme")
  endif

  call get_flags_and_clip()
  call do_all_ghost(mom)

end subroutine swpmom
!  Implements the CIAM (Lagrangian Explicit, onto square)
!  advection method of Jie Li. 
! 
! ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
! split advection of the interface along the x (d=1), y (d=2) and z (d=3)
! directions
! ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
subroutine swpz(us,c,f,d,vof1,vof2,vof3)
  !***
  use module_grid
  use module_flow
  use module_vof
  implicit none
  logical error
  integer i,j,k
  integer i0,j0,k0
  integer i1,j1,k1
  integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: f
  integer, intent(in) :: d
  real (8)  , dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: us
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: c,vof1,vof2,vof3
  real(8) mm1,mm2
  real(8) a1,a2,alpha,fl3d, stencil3x3(-1:1,-1:1,-1:1)
  real(8) dm(3)
  real(8) fl3dnew,nr(3),deltax(3),x0(3)
  intrinsic dmax1,dmin1
  !***
  call init_i0j0k0 (d,i0,j0,k0)

  if(ng.lt.2) call pariserror("wrong ng")
  do k=ks-1,ke+1
     do j=js-1,je+1
        do i=is-1,ie+1
           if (d.eq.1) then
              a2 = us(i,j,k)*dt/dxh(i)
              a1 = us(i-1,j,k)*dt/dxh(i-1)
           elseif (d.eq.2) then
              a2 = us(i,j,k)*dt/dyh(j)
              a1 = us(i,j-1,k)*dt/dyh(j-1)
           elseif (d.eq.3) then
              a2 = us(i,j,k)*dt/dzh(k)
              a1 = us(i,j,k-1)*dt/dzh(k-1)
           endif
           !***
           !     3 cases: 1: default (c=0. and fluxes=0.); 2: c=1.; 3:c>0.
           !***
           vof1(i,j,k) = 0.0d0
           vof2(i,j,k) = 0.0d0
           vof3(i,j,k) = 0.0d0

           ! we need to introduce full/empty flags

           if (c(i,j,k) .eq. 1.0d0) then
              vof1(i,j,k) = dmax1(-a1,0.d0)
              vof2(i,j,k) = 1.d0 - dmax1(a1,0.d0) + dmin1(a2,0.d0)
              vof3(i,j,k) = dmax1(a2,0.d0)

           else if (c(i,j,k) .gt. 0.d0) then
              do i1=-1,1; do j1=-1,1; do k1=-1,1
                stencil3x3(i1,j1,k1) = c(i+i1,j+j1,k+k1)
              enddo;enddo;enddo

              if(oldvof) then
                 call fit_plane(c(i,j,k),d,a1,a2,stencil3x3,dm,alpha,error)
                 if (error) cycle
              else
                 call fit_plane_new(c(i,j,k),d,a1,a2,stencil3x3,nr,alpha,error)
                 if(error) cycle
              endif

              mm1 = dmax1(a1,0.0d0)
              mm2 = 1.d0 - mm1 + dmin1(0.d0,a2)

              if(oldvof) then
                 if (a1 .lt. 0.d0) &
                      vof1(i,j,k) = fl3d(dm(1),dm(2),dm(3),alpha,a1  ,-a1)
                 if (a2 .gt. 0.d0) &
                      vof3(i,j,k) = fl3d(dm(1),dm(2),dm(3),alpha,1.d0,a2)
                 vof2(i,j,k) = fl3d(dm(1),dm(2),dm(3),alpha,mm1,mm2)
              else
                 x0=0d0
                 deltax=1d0
                 if(a1.lt.0d0) then
                    x0(d)=a1
                    deltax(d)=-a1
                    vof1(i,j,k) = fl3dnew(nr,alpha,x0,deltax)
                 endif
                 if(a2.gt.0d0) then
                    x0(d)=1d0
                    deltax(d)=a2
                    vof3(i,j,k) = fl3dnew(nr,alpha,x0,deltax)
                 endif
                 x0(d)=mm1
                 deltax(d)=mm2
                 vof2(i,j,k) = fl3dnew(nr,alpha,x0,deltax)
              endif

           endif
        enddo
     enddo
  enddo
  !
  ! assume that ghost layers take care of the boundary conditions. 
  ! so i-1, i+1 needs to be computed. 
  ! at least the ghost layers is-2, is-1, ie+1,ie+2  need to be there
  ! at the beginning of the subroutine, so that fluxes vof1,vof3 are computed
  ! for is-1, ie+1. 
  !    (1) new values of c and  clip it: 0. <= c <= 1.
  !    (2) apply proper boundary conditions to c
  !*(1)* 
  do k=ks,ke
     do j=js,je
        do i=is,ie
           c(i,j,k) = vof3(i-i0,j-j0,k-k0) + vof2(i,j,k) + vof1(i+i0,j+j0,k+k0)
           c(i,j,k) = dmax1(0.0d0,dmin1(1.0d0,c(i,j,k)))
        enddo
     enddo
  enddo
  !*(2)*
  call setvofbc(c,f)
  !***
end subroutine swpz

subroutine swpzmom(us,c,f,d,mom1,mom2,mom3,mom)
!  !***
  use module_grid
  use module_flow
  use module_vof
  use module_BC
  implicit none
  logical error
  integer i,j,k
  integer i0,j0,k0
  integer i1,j1,k1
  integer inv(3)
  integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: f
  integer, intent(in) :: d
  real(8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: us
  real(8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: mom
  real(8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: c
  real(8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: mom1,mom2,mom3
  real(8) mm1,mm2,vof
  real(8) a1,a2,alpha,fl3d,uavg,rhoavg
  real(8) mxyz(3), dm(3), stencil3x3(-1:1,-1:1,-1:1)
  intrinsic dmax1,dmin1
  !***
  call init_i0j0k0 (d,i0,j0,k0)

  if(ng.lt.2) call pariserror("wrong ng")
  do k=ks-1,ke+1
     do j=js-1,je+1
        do i=is-1,ie+1
           if (d.eq.1) then
              a2 = us(i,j,k)*dt/dxh(i)
              a1 = us(i-1,j,k)*dt/dxh(i-1)
           elseif (d.eq.2) then
              a2 = us(i,j,k)*dt/dyh(j)
              a1 = us(i,j-1,k)*dt/dyh(j-1)
           elseif (d.eq.3) then
              a2 = us(i,j,k)*dt/dzh(k)
              a1 = us(i,j,k-1)*dt/dzh(k-1)
           endif
            
           !momentum for full cells
           mm1 = dmax1(a1,0.0d0)
           mm2 = 1.d0 - mm1 + dmin1(0.d0,a2)

           mom1(i,j,k)  = dmax1(-a1,0.d0)*mom(i,j,k)
           mom3(i,j,k)  = dmax1(a2,0.d0) *mom(i,j,k)
           mom2(i,j,k)  = mom(i,j,k) - mom1(i,j,k) - mom3(i,j,k)
           uavg = mom(i,j,k)/(rho2*c(i,j,k)+rho1*(1.d0-c(i,j,k)))

           if ((c(i,j,k) .gt. 0.d0).and.(c(i,j,k) .lt. 1.d0)) then

              do i1=-1,1; do j1=-1,1; do k1=-1,1
                stencil3x3(i1,j1,k1) = c(i+i1,j+j1,k+k1)
              enddo;enddo;enddo
              call fit_plane(c(i,j,k),d,a1,a2,stencil3x3,dm,alpha,error)
              if (error) cycle

              if (a1 .lt. 0.d0) then
                     vof = fl3d(dm(1),dm(2),dm(3),alpha,a1  ,-a1)
                     mom1(i,j,k) = (rho2*vof + rho1*(-a1 - vof))*uavg
              endif
              if (a2 .gt. 0.d0) then
                     vof = fl3d(dm(1),dm(2),dm(3),alpha,1.d0,a2)
                     mom3(i,j,k) = (rho2*vof + rho1*(a2 - vof))*uavg
              endif
              
              vof = fl3d(dm(1),dm(2),dm(3),alpha,mm1,mm2)
              mom2(i,j,k) = (rho2*vof + rho1*(mm2 - vof))*uavg

           endif
        enddo
     enddo
  enddo

  do k=ks,ke
    do j=js,je
      do i=is,ie
        mom(i,j,k)  = mom1(i+i0,j+j0,k+k0)+mom2(i,j,k)+mom3(i-i0,j-j0,k-k0)
      enddo
    enddo
  enddo

  call SetMomentumBC(us,c,mom,d,umask,rho1,rho2) 
  !***
end subroutine swpzmom

subroutine fit_plane(vof,d,a1,a2,stencil3x3,dm,alpha,error)
  use module_grid
  !***
  !     (1) normal vector: dmx,dmy,dmz, and |dmx|+|dmy|+|dmz| = 1.
  !     (2) dmx,dmy,dmz>0.
  !     (3) get alpha;               (4) back to original plane;
  !     (5) lagrangian advection;    (6) get fluxes
  !*(1)*
  integer i0,j0,k0
  integer inv(3)
  integer, intent(in) :: d
  real(8) a1,a2,alpha,al3d,vof
  real(8) mxyz(3), dm(3), stencil3x3(-1:1,-1:1,-1:1)
  logical error
  intrinsic dmax1,dmin1

  error=.TRUE.
  call mycs(stencil3x3,mxyz)
  !*(2)*  
  inv(:) = 1
  do i0=1,3
    if (mxyz(i0) .lt. 0.0d0) then
      mxyz(i0) = -mxyz(i0)
      inv(i0) = -1
    endif
  enddo
  !*(3)*  
! TEMPORARY - Note: avoid calculating alpha when mxyz is zero, which occurs when cvof is 
!                   a very small non-zero number in an isolated cell  
  if ( mxyz(1) == 0.d0 .and. mxyz(2) == 0.d0 .and. mxyz(3) == 0.d0 ) return
! END TEMPORARY  

  alpha = al3d(mxyz(1),mxyz(2),mxyz(3),vof)
  !*(4)*  
  mxyz(:) = inv*mxyz
  alpha = alpha + dmin1(0.d0,mxyz(1)) + dmin1(0.d0,mxyz(2)) + &
        dmin1(0.d0,mxyz(3))
  !*(5)*  
  mxyz(d) = mxyz(d)/(1.0d0 - a1 + a2)
  alpha = alpha + mxyz(d)*a1

  if (d.eq.1) then
    dm = mxyz
    elseif (d.eq.2) then
    dm(1) = mxyz(2)
    dm(2) = mxyz(3)
    dm(3) = mxyz(1)
    elseif (d.eq.3) then
    dm(1) = mxyz(3)
    dm(2) = mxyz(1)
    dm(3) = mxyz(2)
  endif
  
  error=.FALSE.

end subroutine fit_plane

subroutine fit_plane_new(vof,d,a1,a2,stencil3x3,mxyz,alpha,error)
  use module_grid
  integer i0,j0,k0
  integer inv(3)
  integer, intent(in) :: d
  real(8) a1,a2,alpha,al3dnew,vof
  real(8) mxyz(3), stencil3x3(-1:1,-1:1,-1:1)
  logical error
  intrinsic dmax1,dmin1

  error=.TRUE.
  call mycs(stencil3x3,mxyz)
! TEMPORARY - Note: avoid calculating alpha when mxyz is zero, which occurs when cvof is 
!                   a very small non-zero number in an isolated cell  
  if ( mxyz(1) == 0.d0 .and. mxyz(2) == 0.d0 .and. mxyz(3) == 0.d0 ) return
! END TEMPORARY  

  alpha = al3dnew(mxyz,vof)
  mxyz(d) = mxyz(d)/(1.0d0 - a1 + a2)
  alpha = alpha + mxyz(d)*a1
  
  error=.FALSE.

end subroutine fit_plane_new

!
!=================================================================================================
! split 1D advection of the interface along the x,y,z (d=1,2,3) directions
!
! Following the advection method of Weymouth & Yue : Weymouth, G D, and Dick K P Yue,
! "Conservative Volume-of-Fluid Method for Free-Surface Simulations on Cartesian-Grids."
! Journal of Computational Physics 229, no. 8 (April 2010): 2853-2865. doi:10.1016/j.jcp.2009.12.018.
!=================================================================================================
!=================================================================================================
SUBROUTINE swpr(us,c,f,dir,vof1,cg,vof3)
!***
    USE module_grid
    USE module_flow
    USE module_vof
    
    IMPLICIT NONE
    include 'mpif.h'
    INTEGER :: i,j,k
    INTEGER :: invx,invy,invz,ii,jj,kk,i0,j0,k0
    INTEGER, INTENT(IN) :: dir
    REAL (8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), INTENT(IN) :: us
    REAL (8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), INTENT(INOUT) :: c,vof1,cg,vof3
    integer, dimension(imin:imax,jmin:jmax,kmin:kmax),  intent(inout) :: f
    REAL(8), TARGET :: dmx,dmy,dmz,dxyz
    REAL(8), POINTER :: dm1,dm2,dm3
    REAL(8) :: a1,a2,alpha,AL3D,FL3D
    REAL(8) :: AL3DNEW, FL3DNEW, x0(3), deltax(3)
    real(8) :: mxyz(3),stencil3x3(-1:1,-1:1,-1:1)
    INTRINSIC DMAX1,DMIN1
!
  if(ng < 2) call pariserror("wrong ng")
  ii=0; jj=0; kk=0
  if (dir == 1) then
     ii=1; dm1 => dmx;  dm2 => dmy;  dm3 => dmz 
  else if (dir == 2) then
     jj=1; dm1 => dmy;  dm2 => dmz;  dm3 => dmx 
  else if (dir == 3) then
     kk=1; dm1 => dmz;  dm2 => dmx;  dm3 => dmy 
  endif
  dxyz = dxh(is)
  call test_cell_size()

  do k=ks-1,ke+1
     do j=js-1,je+1
        do i=is-1,ie+1
           a2 = us(i,j,k)*dt/dxyz
           a1 = us(i-ii,j-jj,k-kk)*dt/dxyz
           !  default: fluxes=0. (good also for c=0.)
           vof1(i,j,k) = 0.d0
           vof3(i,j,k) = 0.d0
           !  c = 1.
           if (c(i,j,k) == 1.0d0) then
              vof1(i,j,k) = DMAX1(-a1,0.d0)
              vof3(i,j,k) = DMAX1(a2,0.d0)
           ! 0. < c < 1.
           else if (c(i,j,k) > 0.d0) then
              ! local stencil and normal vector: |dmx|+|dmy|+|dmz| = 1.
              do i0=-1,1; do j0=-1,1; do k0=-1,1
                 stencil3x3(i0,j0,k0) = c(i+i0,j+j0,k+k0)
              enddo;enddo;enddo
              call mycs(stencil3x3,mxyz)
              if(oldvof) then
                 dmx = mxyz(1); dmy = mxyz(2); dmz = mxyz(3)
                 ! positive dmx,dmy,dmz
                 invx = 1; invy = 1; invz = 1
                 if (dmx < 0.d0) then
                    dmx = -dmx; invx = -1
                 endif
                 if (dmy < 0.d0) then
                    dmy = -dmy; invy = -1
                 endif
                 if (dmz < 0.0d0) then
                    dmz = -dmz; invz = -1
                 endif
                 ! get alpha
                 alpha = AL3D(dmx,dmy,dmz,c(i,j,k))
                 ! back to the original plane
                 dmx = invx*dmx
                 dmy = invy*dmy
                 dmz = invz*dmz
                 alpha = alpha + DMIN1(0.d0,dmx) + DMIN1(0.d0,dmy) + DMIN1(0.d0,dmz)
              else
                 alpha = AL3DNEW(mxyz,c(i,j,k))
              endif

              ! Eulerian advection
              if(oldvof) then
                 if (a1 < 0.d0) &
                      vof1(i,j,k) = FL3D(dm1,dm2,dm3,alpha,0.d0,-a1)
                 if (a2 > 0.d0) &
                      vof3(i,j,k) = FL3D(dm1,dm2,dm3,alpha,1.d0-a2,a2)
              else
                 x0=0d0
                 deltax=1d0
                 if(a1<0d0) then
                    deltax(dir)=-a1
                    vof1(i,j,k) = FL3DNEW(mxyz,alpha,x0,deltax)
                 endif
                 if(a2>0d0) then
                    x0(dir)=1d0-a2
                    deltax(dir)=a2
                    vof3(i,j,k) = FL3DNEW(mxyz,alpha,x0,deltax)
                 endif
              endif
              
           endif
        enddo
     enddo
  enddo
  ! assume that ghost layers take care of the boundary conditions, then 
  ! fluxes vof1,vof3 must be computed for is-1, ie+1 
  ! new clipped values of c (0. <= c <= 1)
  do k=ks,ke
     do j=js,je
        do i=is,ie
           a2 = us(i,j,k)*dt/dxyz
           a1 = us(i-ii,j-jj,k-kk)*dt/dxyz
           c(i,j,k) = c(i,j,k) - (vof3(i,j,k) - vof1(i+ii,j+jj,k+kk)) + & 
                      (vof3(i-ii,j-jj,k-kk) - vof1(i,j,k)) + cg(i,j,k)*(a2-a1);
!!$           c(i,j,k) = DMAX1(0.d0,DMIN1(1.d0,c(i,j,k)))
           if (c(i,j,k) < EPSC) then
              c(i,j,k) = 0.d0
           elseif (c(i,j,k) >  (1.d0 - EPSC)) then
              c(i,j,k) = 1.d0
           endif
        enddo
     enddo
  enddo
  ! apply proper boundary conditions to c
  call setvofbc(c,f)
end subroutine swpr

subroutine init_i0j0k0 (d,i0,j0,k0)
  integer d,i0,j0,k0

  i0=0;j0=0;k0=0

  if (d.eq.1) then
    i0=1
  elseif (d.eq.2) then
    j0=1
  elseif (d.eq.3) then
    k0=1
  endif

end subroutine init_i0j0k0
!=================================================================================================
!=================================================================================================
! ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
! PROGRAM TO FIND alpha IN: m1 x1 + m2 x2 + m3 x3 = alpha,
! GIVEN m1+m2+m3=1 (all > 0) AND THE VOLUMETRIC FRACTION cc
! ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
function al3d(b1,b2,b3,cc)
  !***
  implicit none
  real(8) m1,m2,m3,cc,b1,b2,b3,tmp,pr,ch,mm,m12
  real(8) p,p12,q,teta,cs,al3d
  real(8) untier,v1,v2,v3
  parameter (untier=1.d0/3.d0)
  intrinsic dmax1,dmin1,dsqrt,dacos,dcos
  !***  
  !     (1) order coefficients: m1<m2<m3; (2) get ranges: v1<v2<v3;
  !     (3) limit ch (0.d0 < ch < 0.5d0); (4) calculate alpha
  !*(1)* 
  m1 = dmin1(b1,b2)
  m3 = dmax1(b1,b2)
  m2 = b3
  if (m2 .lt. m1) then
     tmp = m1
     m1 = m2
     m2 = tmp
  else if (m2 .gt. m3) then
     tmp = m3
     m3 = m2
     m2 = tmp
  endif
  !*(2)*
  m12 = m1 + m2 
  pr  = DMAX1(6.d0*m1*m2*m3,1.d-50)
  V1  = m1*m1*m1/pr
  V2  = V1 + 0.5d0*(m2-m1)/m3
  if (m3 .LT. m12) then
     mm = m3
     V3 = (m3*m3*(3.d0*m12-m3) + m1*m1*(m1-3.d0*m3) +&
          m2*m2*(m2-3.d0*m3))/pr
  else
     mm = m12
     V3 = 0.5d0*mm/m3
  endif
  !*(3)*
  ch = DMIN1(cc,1.d0-cc)
  !*(4)*      
  if (ch .LT. V1) then
     !***         AL3D = cbrt(pr*ch)
     AL3D = (pr*ch)**UNTIER
  else if (ch .LT. V2) then
     AL3D = 0.5d0*(m1 + DSQRT(m1*m1 + 8.d0*m2*m3*(ch-V1)))
  else if (ch .LT. V3) then
     p = 2.d0*m1*m2
     q = 1.5d0*m1*m2*(m12 - 2.d0*m3*ch)
     p12 = DSQRT(p)
     teta = DACOS(q/(p*p12))/3.d0
     cs = DCOS(teta)
     AL3D = p12*(DSQRT(3.d0*(1.d0-cs*cs)) - cs) + m12
  else if (m12 .LT. m3) then
     AL3D = m3*ch + 0.5d0*mm
  else 
     p = m1*(m2+m3) + m2*m3 - 0.25d0
     q = 1.5d0*m1*m2*m3*(0.5d0-ch)
     p12 = DSQRT(p)
     teta = DACOS(q/(p*p12))/3.0
     cs = DCOS(teta)
     AL3D = p12*(DSQRT(3.d0*(1.d0-cs*cs)) - cs) + 0.5d0
  endif

  if (cc .GT. 0.5d0)  AL3D = 1.d0 - AL3D
  !***
  return
end function al3d
! ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
! PROGRAM TO FIND THE "CUT VOLUME" V0 GIVEN r0, dr0 AND
! m1 x1 + m2 x2 + m3 x3 = alpha
! ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
function fl3d(m1,m2,m3,alpha,r0,dr0)
  !***
  implicit none
  real(8) m1,m2,m3,alpha,r0,dr0,fl3D
  real(8) al,al0,n1,n2,n3,b1,b2,b3,b12,bm,tmp,pr
  INTRINSIC DMAX1,DMIN1,DABS
  !***
  !     (1) move origin to r0 along r ;  (2) reflect parallelepiped;
  !     (3) limit alpha (0<= al0 <=0.5); (4) order coefficients: b1<b2<b3;
  !     (5) calculate volume (NOTE: it is assumed:s0=t0=0; ds0=dt0=1.)
  !*(1)*
  al = alpha - m1*r0
  !*(2)*
  al = al + DMAX1(0.d0,-m1*dr0)+DMAX1(0.d0,-m2)+DMAX1(0.d0,-m3)
  tmp = DABS(m1)*dr0 + DABS(m2) + DABS(m3)
  n1 = DABS(m1)/tmp
  n2 = DABS(m2)/tmp
  n3 = DABS(m3)/tmp
  al = DMAX1(0.d0,DMIN1(1.d0,al/tmp))
  !*(3)*
  al0 = DMIN1(al,1.d0-al)
  !*(4)*
  b1 = DMIN1(n1*dr0,n2)
  b3 = DMAX1(n1*dr0,n2)
  b2 = n3
  if (b2 .LT. b1) then
     tmp = b1
     b1 = b2
     b2 = tmp
  else if (b2 .GT. b3) then
     tmp = b3
     b3 = b2
     b2 = tmp
  endif
  b12 = b1 + b2
  bm = DMIN1(b12,b3)
  pr = DMAX1(6.d0*b1*b2*b3,1.0d-50)
  !*5*     
  if (al0 .LT. b1) then
     tmp = al0*al0*al0/pr
  else if (al0 .LT. b2) then
     tmp = 0.5d0*al0*(al0-b1)/(b2*b3) +  b1*b1*b1/pr
  else if (al0 .LT. bm) then
     tmp = (al0*al0*(3.d0*b12-al0) + b1*b1*(b1-3.d0*al0) +&
          b2*b2*(b2-3.d0*al0))/pr
  else if (b12 .LT. b3) then
     tmp = (al0 - 0.5d0*bm)/b3
  else
     tmp = (al0*al0*(3.d0-2.d0*al0) + b1*b1*(b1-3.d0*al0) +&
          b2*b2*(b2-3.d0*al0) + b3*b3*(b3-3.d0*al0))/pr
  endif

  if (al .LE. 0.5d0) then
     FL3D = tmp*dr0
  else
     FL3D = (1.d0-tmp)*dr0
  endif
  !***  
  return
end function fl3d

subroutine ls2vof_in_cell(stencil3x3,c,nflag)
  implicit none
  real(8), intent(out):: c
  integer, intent(out):: nflag
  real(8) :: zero, one, norml1
  real(8) :: mx,my,mz,alpha
  real(8) :: fl3d
  real(8) :: mxyz(3),stencil3x3(-1:1,-1:1,-1:1)

  zero=0.d0
  one=1.d0
  !***
  !     (1) gradient*32: mx,my,mz; (2) mx,my,mz>0. and mx+my+mz = 1.;
  !     (3) normalize alpha = level set at center. Cell units. 
  !     (4) shift alpha to origin=vertex;   (5) get volume from alpha.  
  !
  !     *(1)*  
  !***
  call fd32(stencil3x3,mxyz)
  !***
  !     *(2)*  
  !***
  mx = dabs(mxyz(1)); my = dabs(mxyz(2)); mz = dabs(mxyz(3))
  norml1 = mx+my+mz
  mx = mx/norml1;     my = my/norml1;     mz = mz/norml1
  !***
  !     *(3)*  
  !***
  ! the factor is 32 because grad ls=(1,0,0) gives mx=32.
  alpha = 32.d0*stencil3x3(0,0,0)/norml1   
  !***
  !     *(4)*  
  !***
  alpha = alpha + 0.5d0
  !***
  !     *(5)*  
  !***
  if(alpha.ge.1.d0) then 
     c = 1.d0
     nflag = 1
  else if (alpha.le.0.d0) then
     c = 0.d0
     nflag = 0 
  else 
     c = fl3d(mx,my,mz,alpha,zero,one)
     nflag = 2
  end if
  return
end subroutine ls2vof_in_cell
!
! *-----------------------------------------------------* 
! *  MYC - Mixed Youngs and Central Scheme              *
! * returns normal normalized so that |mx|+|my|+|mz| = 1* 
! *-----------------------------------------------------*
! 
!
!Known problems: the index (1,1,1), i.e. the central cell
!in the block, never occurs: neither in the central scheme
!nor in Youngs' method. Therefore an isolated droplet will have
!a normal with all components to zero. I took care of the
!division-by-zero issue, but not of this one.
!
!Ruben
!
!
! Translated into f90 by Stephane Z.
!
subroutine mycs(c,mxyz)
  !***
  implicit none
  real(8) c(0:2,0:2,0:2)
  real(8) mxyz(0:2)
  real(8) m1,m2,m(0:3,0:2),t0,t1,t2
  integer cn
  real(8), parameter  :: NOT_ZERO=1.e-30

  ! write the plane as: sgn(mx) X =  my Y +  mz Z + alpha 
  !                           m00 X = m01 Y + m02 Z + alpha 

  m1 = c(0,1,0) + c(0,1,2) + c(0,0,1) + c(0,2,1) + &
       c(0,1,1)
  m2 = c(2,1,0) + c(2,1,2) + c(2,0,1) + c(2,2,1) + &
       c(2,1,1)

  if(m1>m2) then
     m(0,0) = 1.
  else
     m(0,0) = -1.
  end if

  m1 = c(0,0,1)+ c(2,0,1)+ c(1,0,1)
  m2 = c(0,2,1)+ c(2,2,1)+ c(1,2,1)
  m(0,1) = 0.5*(m1-m2)

  m1 = c(0,1,0)+ c(2,1,0)+ c(1,1,0)
  m2 = c(0,1,2)+ c(2,1,2)+ c(1,1,2)
  m(0,2) = 0.5*(m1-m2)

  ! write the plane as: sgn(my) Y =  mx X +  mz Z + alpha, 
  !                          m11 Y = m10 X + m12 Z + alpha.

  m1 = c(0,0,1) + c(0,2,1) + c(0,1,1)
  m2 = c(2,0,1) + c(2,2,1) + c(2,1,1)
  m(1,0) = 0.5*(m1-m2)

  m1 = c(1,0,0) + c(1,0,2) + c(2,0,1) + c(0,0,1) +&
       c(1,0,1)
  m2 = c(1,2,0) + c(1,2,2) + c(2,2,1) + c(0,2,1) +&
       c(1,2,1)


  if(m1>m2) then
     m(1,1) = 1.
  else
     m(1,1) = -1.
  end if

  m1 = c(1,0,0)+ c(1,1,0)+ c(1,2,0)
  m2 = c(1,0,2)+ c(1,1,2)+ c(1,2,2)
  m(1,2) = 0.5*(m1-m2)

  ! write the plane as: sgn(mz) Z =  mx X +  my Y + alpha 
  !                          m22 Z = m20 X + m21 Y + alpha

  m1 = c(0,1,0)+ c(0,1,2)+ c(0,1,1)
  m2 = c(2,1,0)+ c(2,1,2)+ c(2,1,1)
  m(2,0) = 0.5*(m1-m2)

  m1 = c(1,0,0)+ c(1,0,2)+ c(1,0,1)
  m2 = c(1,2,0)+ c(1,2,2)+ c(1,2,1)
  m(2,1) = 0.5*(m1-m2)

  m1 = c(0,1,0) + c(2,1,0) + c(1,0,0) + c(1,2,0) +&
       c(1,1,0)
  m2 = c(0,1,2) + c(2,1,2) + c(1,0,2) + c(1,2,2) +&
       c(1,1,2)

  if(m1>m2) then
     m(2,2) = 1.
  else
     m(2,2) = -1.
  end if

  ! normalize each set (mx,my,mz): |mx|+|my|+|mz| = 1

  t0 = DABS(m(0,0)) + DABS(m(0,1)) + DABS(m(0,2))
  m(0,0) = m(0,0)/t0
  m(0,1) = m(0,1)/t0
  m(0,2) = m(0,2)/t0

  t0 = DABS(m(1,0)) + DABS(m(1,1)) + DABS(m(1,2))
  m(1,0) = m(1,0)/t0
  m(1,1) = m(1,1)/t0
  m(1,2) = m(1,2)/t0

  t0 = DABS(m(2,0)) + DABS(m(2,1)) + DABS(m(2,2))
  m(2,0) = m(2,0)/t0
  m(2,1) = m(2,1)/t0
  m(2,2) = m(2,2)/t0

  ! choose among the three central schemes */ 
  t0 = DABS(m(0,0))
  t1 = DABS(m(1,1))
  t2 = DABS(m(2,2))

  cn = 0
  if (t1 > t0) then
    t0 = t1
    cn = 1
  endif

  if (t2 > t0) cn = 2

  ! Youngs-CIAM scheme */  
  
  call fd32(c,m(3,0:2))

  ! normalize the set (mx,my,mz): |mx|+|my|+|mz| = 1 

  t0 = DABS(m(3,0)) + DABS(m(3,1)) + DABS(m(3,2)) + NOT_ZERO
  m(3,0) = m(3,0)/t0
  m(3,1) = m(3,1)/t0
  m(3,2) = m(3,2)/t0

  ! choose between the previous choice and Youngs-CIAM 
  t0 = DABS (m(3,0))
  t1 = DABS (m(3,1))
  t2 = DABS (m(3,2))
  if (t1 > t0)  t0 = t1
  if (t2 > t0)  t0 = t2

  if (DABS(m(cn,cn)) > t0)  cn = 3

  ! components of the normal vector */
  mxyz(0) = m(cn,0)
  mxyz(1) = m(cn,1)
  mxyz(2) = m(cn,2)

  return 
  end subroutine mycs

! *----------------------------------------------------------------* 
! *  FD32 - Youngs Finite Difference Gradient Scheme               *
! *  the gradient is computed with a multiplicative factor of -32: *
! *  mm = - 32 * grad (c)                                          *
! *----------------------------------------------------------------*
! 
!Known problems: the index (1,1,1), i.e. the central cell
!in the block, never occurs:
!Therefore an isolated droplet will have
!a normal with all components to zero. 
!
!Ruben
!
!
! Translated into f90 by Stephane Z.
!
subroutine fd32(c,mm)
  !***
  implicit none
  real(8), intent(in) :: c(0:2,0:2,0:2)
  real(8), intent(out) :: mm(0:2)
  real(8) :: m1,m2

  m1 = c(0,0,0) + c(0,2,0) + c(0,0,2) + c(0,2,2) +&
       2.d0*(c(0,0,1) + c(0,2,1) + c(0,1,0) + c(0,1,2)) +&
       4.d0*c(0,1,1)
  m2 = c(2,0,0) + c(2,2,0) + c(2,0,2) + c(2,2,2) +&
       2.d0*(c(2,0,1) + c(2,2,1) + c(2,1,0) + c(2,1,2)) +&
       4.d0*c(2,1,1)
  mm(0) = m1-m2

  m1 = c(0,0,0) + c(0,0,2) + c(2,0,0) + c(2,0,2) +&
       2.d0*(c(0,0,1) + c(2,0,1) + c(1,0,0) + c(1,0,2)) +&
       4.d0*c(1,0,1)
  m2 = c(0,2,0) + c(0,2,2) + c(2,2,0) + c(2,2,2) +&
       2.d0*(c(0,2,1) + c(2,2,1) + c(1,2,0) + c(1,2,2)) +&
       4.d0*c(1,2,1)
  mm(1) = m1-m2

  m1 = c(0,0,0) + c(0,2,0) + c(2,0,0) + c(2,2,0) +&
       2.d0*(c(0,1,0) + c(2,1,0) + c(1,0,0) + c(1,2,0)) +&
       4.d0*c(1,1,0)
  m2 = c(0,0,2) + c(0,2,2) + c(2,0,2) + c(2,2,2) +&
       2.d0*(c(0,1,2) + c(2,1,2) + c(1,0,2) + c(1,2,2)) +&
       4.d0*c(1,1,2)
  mm(2) = m1-m2

  return 
  end subroutine fd32

!
! *----------------------------------------------------------------* 
! *  youngs - Youngs Finite Difference Gradient Scheme             *
! *  the gradient is normed so that |mx|+|my|+|mz| = 1             *
! *----------------------------------------------------------------*
! 
! Known problems: the index (1,1,1), i.e. the central cell
! in the block, never occurs:
! Therefore an isolated droplet will have
! a normal with all components to zero. 
!
!
subroutine youngs(c,mm)
  !***
  implicit none
  real(8), intent(in) :: c(0:2,0:2,0:2)
  real(8), intent(out) :: mm(0:2)
  integer :: i
  real(8) :: norm

  call fd32(c,mm)
  norm = abs(mm(0))+abs(mm(1))+abs(mm(2)) ! + TINY
  do i=0,2
     mm(i) = mm(i)/norm
  enddo
  return 
end subroutine youngs

! NEW VOF FUNCTIONS

!===============================================================================
! DESCRIPTION OF FUNCTION AL3DNEW:
! compute in the unit cube the coefficient alpha that satisfies
! [nr]*[x] = nr1*x1 + nr2*x2 + nr3*x3 = alpha
! INPUT: normal coefficients (nr1,nr2,nr3) and volume fraction cc 
! NOTE : the normal coefficients MUST satisfy the relation |nr1|+|nr2|+|nr3|=1
!-------------------------------------------------------------------------------
FUNCTION AL3DNEW(nr,cc)

  IMPLICIT NONE
  REAL(8), INTENT(IN):: nr(3),cc
  REAL(8) :: AL3DNEW
  REAL(8) :: cch,c01,c02,c03,np1,np2,np3
  REAL(8) :: m1,m2,m3,m12,numer,denom,p,pst,q,arc,csarc
  REAL(8), PARAMETER :: athird=1.d0/3.d0 
  INTRINSIC DABS,DMIN1,DMAX1,DSQRT,DACOS,DCOS

  np1 = DABS(nr(1))                                 ! need positive coefficients
  np2 = DABS(nr(2))
  np3 = DABS(nr(3))
  m1 = DMIN1(np1,np2)                              ! order positive coefficients
  m3 = DMAX1(np1,np2)
  if (np3 < m1) then
     m2 = m1
     m1 = np3
  else if (np3 >= m3) then
     m2 = m3
     m3 = np3
   else
     m2 = np3
  endif
  cch = DMIN1(cc,1.d0-cc)                              ! symmetry: 0 < cch < 1/2
  denom = DMAX1(6.d0*m1*m2*m3,1.d-50)                           ! get cch ranges
  c01 = m1*m1*m1/denom
  c02  = c01 + 0.5d0*(m2-m1)/m3
  m12 = m1 + m2
  if (m12 <= m3) then
     c03 = 0.5d0*m12/m3
  else
     numer = m3*m3*(3.d0*m12-m3) + m1*m1*(m1-3.d0*m3) + m2*m2*(m2-3.d0*m3)
     c03 = numer/denom
  endif

! 1: C<=C1; 2: C1<=C<=C2; 3: C2<=C<=C3; 4: C3<=C<=1/2 (a: m12<=m3; b: m3<m12)) 
  if (cch <= c01) then
     AL3DNEW = (denom*cch)**athird                                       ! case (1)
  else if (cch <= c02) then 
     AL3DNEW = 0.5d0*(m1 + DSQRT(m1*m1 + 8.d0*m2*m3*(cch-c01)))          ! case (2)
  else if (cch <= c03) then
     p = 2.d0*m1*m2
     q = 1.5d0*m1*m2*(m12 - 2.d0*m3*cch)
     pst = DSQRT(p)
     arc = athird*DACOS(q/(p*pst))
     csarc = DCOS(arc)
     AL3DNEW = pst*(DSQRT(3.d0*(1.d0-csarc*csarc)) - csarc) + m12        ! case (3)
  else if (m12 <= m3) then
     AL3DNEW = m3*cch + 0.5d0*m12                                       ! case (4a)
  else                                  
     p = m12*m3 + m1*m2 - 0.25d0                                     
     q = 1.5d0*m1*m2*m3*(0.5d0-cch)
     pst = DSQRT(p)
     arc = athird*DACOS(q/(p*pst))
     csarc = DCOS(arc)
     AL3DNEW = pst*(DSQRT(3.d0*(1.d0-csarc*csarc)) - csarc) + 0.5d0     ! case (4b)
  endif

  if (cc > 0.5d0)  AL3DNEW = 1.d0 - AL3DNEW
  
! compute alpha for the real coefficients  
  AL3DNEW = AL3DNEW + DMIN1(0.d0,nr(1)) + DMIN1(0.d0,nr(2)) + DMIN1(0.d0,nr(3))
  
END FUNCTION AL3DNEW

!===============================================================================
! DESCRIPTION OF FUNCTION CC3DNEW:
! compute in the unit cube the volume cut by the plane
! [nr]*[x] = nr1*x1 + nr2*x2 + nr3*x3 = alpha
! where the normal is pointing outward from the reference phase
! INPUT: normal coefficients (nr1,nr2,nr3) and plane constant alpha 
! NOTE : the normal coefficients MUST satisfy the relation |nr1|+|nr2|+|nr3|=1
!-------------------------------------------------------------------------------
FUNCTION CC3DNEW(nr,alpha)

  IMPLICIT NONE
  REAL(8), INTENT(IN):: nr(3),alpha
  REAL(8) :: CC3DNEW
  REAL(8) :: al,alh,np1,np2,np3,m1,m2,m3,m12,mm,denom,top
  INTRINSIC DMAX1,DMIN1,DABS

  np1 = DABS(nr(1))                                 ! need positive coefficients
  np2 = DABS(nr(2))
  np3 = DABS(nr(3))
  m1 = DMIN1(np1,np2)                              ! order positive coefficients
  m3 = DMAX1(np1,np2)
  if (np3 < m1) then
     m2 = m1
     m1 = np3
  else if (np3 >= m3) then
     m2 = m3
     m3 = np3
   else
     m2 = np3
  endif
  al = alpha + DMAX1(0.d0,-nr(1))+DMAX1(0.d0,-nr(2))+DMAX1(0.d0,-nr(3))
  al = DMAX1(0.d0,DMIN1(1.d0,al))                                  ! safe limits
  alh = DMIN1(al,1.d0-al)                              ! symmetry: 0 < alh < 1/2
  m12 = m1 + m2
  mm = DMIN1(m12,m3)
  denom = DMAX1(6.d0*m1*m2*m3,1.0d-50)

! 1: al<=m1; 2: m1<=al<=m2; 3: m2<=al<=mm; 4: mm<=al<=1/2 (a:m12<=m3; b:m3<m12)) 
  if (alh <= m1) then
    CC3DNEW = alh*alh*alh/denom                                          ! case (1)
  else if (alh <= m2) then
    CC3DNEW = 0.5d0*alh*(alh-m1)/(m2*m3) +  m1*m1*m1/denom               ! case (2)
  else if (alh <= mm) then
    top = alh*alh*(3.d0*m12-alh) + m1*m1*(m1-3.d0*alh)              
    CC3DNEW = (top + m2*m2*(m2-3.d0*alh))/denom                          ! case (3)
  else if (m12 <= m3) then
    CC3DNEW = (alh - 0.5d0*m12)/m3                                      ! case (4a)
  else
    top = alh*alh*(3.d0-2.d0*alh) + m1*m1*(m1-3.d0*alh)             
    CC3DNEW = (top + m2*m2*(m2-3.d0*alh) + m3*m3*(m3-3.d0*alh))/denom   ! case (4b)
  endif
  
  if (al > 0.5d0) CC3DNEW = 1.d0 - CC3DNEW

END FUNCTION CC3DNEW

!===============================================================================
! DESCRIPTION OF FUNCTION FL3DNEW:
! compute in the right hexahedron starting at (x01,x02,x03) and of sides 
! (dx1,dx2,dx3) the volume cut by the plane
! [nr]*[x] = nr1*x1 + nr2*x2 + nr3*x3 = alpha
! INPUT: normal coefficients (nr1,nr2,nr3), plane constant alpha, starting
!        point (x01,x02,x03), sides (dx1,dx2,dx3)
! NOTE : the normal coefficients MUST satisfy the relation |nr1|+|nr2|+|nr3|=1
!-------------------------------------------------------------------------------
FUNCTION FL3DNEW(nr,alpha,x0,dx)

  IMPLICIT NONE
  REAL(8), INTENT(IN):: nr(3),x0(3),dx(3),alpha
  REAL(8) :: FL3DNEW
  REAL(8) :: al,almax,alh,np1,np2,np3,m1,m2,m3,m12,mm,denom,frac,top
  INTRINSIC DMAX1,DMIN1,DABS

! move origin to x0 
  al = alpha - nr(1)*x0(1) - nr(2)*x0(2) - nr(3)*x0(3)
! reflect the figure when negative coefficients
  al = al + DMAX1(0.d0,-nr(1)*dx(1)) + DMAX1(0.d0,-nr(2)*dx(2)) &
          + DMAX1(0.d0,-nr(3)*dx(3)) 
  np1 = DABS(nr(1))                                 ! need positive coefficients
  np2 = DABS(nr(2))
  np3 = DABS(nr(3))
  almax = np1*dx(1) + np2*dx(2) + np3*dx(3)                       
  al = DMAX1(0.d0,DMIN1(1.d0,al/almax))           !get new al within safe limits
  alh = DMIN1(al,1.d0-al)                              ! symmetry: 0 < alh < 1/2


! normalized equation: m1*y1 + m2*y2 + m3*y3 = alh, with m1 <= m2 <= m3
  np1 = np1/almax;
  np2 = np2/almax;
  np3 = np3/almax;
  m1 = DMIN1(np1*dx(1),np2*dx(2))                           ! order coefficients
  m3 = DMAX1(np1*dx(1),np2*dx(2))
  top = np3*dx(3)
  if (top < m1) then
     m2 = m1
     m1 = top
  else if (top >= m3) then
     m2 = m3
     m3 = top
   else
     m2 = top
  endif
  m12 = m1 + m2
  mm = DMIN1(m12,m3)
  denom = DMAX1(6.d0*m1*m2*m3,1.0d-50)

! 1: al<=m1; 2: m1<=al<=m2; 3: m2<=al<=mm; 4: mm<=al<=1/2 (a:m12<=m3; b:m3<m12)) 
  if (alh <= m1) then
     frac = alh*alh*alh/denom                                         ! case (1)
  else if (alh <= m2) then
     frac = 0.5d0*alh*(alh-m1)/(m2*m3) +  m1*m1*m1/denom              ! case (2)
  else if (alh <= mm) then
     top = alh*alh*(3.d0*m12-alh) + m1*m1*(m1-3.d0*alh)              
     frac = (top + m2*m2*(m2-3.d0*alh))/denom                         ! case (3)
  else if (m12 <= m3) then
     frac = (alh - 0.5d0*m12)/m3                                     ! case (4a)
  else
     top = alh*alh*(3.d0-2.d0*alh) + m1*m1*(m1-3.d0*alh)             
     frac = (top + m2*m2*(m2-3.d0*alh) + m3*m3*(m3-3.d0*alh))/denom  ! case (4b)
  endif

  top = dx(1)*dx(2)*dx(3)
  if (al <= 0.5d0) then
     FL3DNEW = frac*top
  else
     FL3DNEW = (1.d0-frac)*top
  endif

END FUNCTION FL3DNEW

!===============================================================================
! DESCRIPTION OF FUNCTION AREA3D:
! compute in the unit cube the area cut by the plane
! [nr]*[x] = nr1*x1 + nr2*x2 + nr3*x3 = alpha
! INPUT: normal coefficients (nr1,nr2,nr3) and volume fraction cc 
! NOTE : the normal coefficients MUST satisfy the relation |nr1|+|nr2|+|nr3|=1;
!        the cut area A is invariant with respects to reflections, ordering  
!        of the coefficients and midpoint (i.e. A(alpha) = A(1-alpha))
!-------------------------------------------------------------------------------
FUNCTION AREA3D(nr,cc)

  IMPLICIT NONE
  REAL(8), INTENT(IN):: nr(3),cc
  REAL(8) :: AREA3D
  REAL(8) :: cch,al,c00,c01,c02,c03,np1,np2,np3
  REAL(8) :: m1,m2,m3,m12,numer,denom,p,pst,q,arc,csarc
  REAL(8), PARAMETER :: athird=1.d0/3.d0 
  INTRINSIC DABS,DMIN1,DMAX1,DSQRT,DACOS,DCOS

  np1 = DABS(nr(1))                                 ! need positive coefficients
  np2 = DABS(nr(2))
  np3 = DABS(nr(3))
  m1 = DMIN1(np1,np2)                              ! order positive coefficients
  m3 = DMAX1(np1,np2)
  if (np3 < m1) then
     m2 = m1
     m1 = np3
  else if (np3 >= m3) then
     m2 = m3
     m3 = np3
   else
     m2 = np3
  endif
  cch = DMIN1(cc,1.d0-cc)                              ! symmetry: 0 < cch < 1/2
  denom = DMAX1(6.d0*m1*m2*m3,1.d-50)                           ! get cch ranges
  c01 = m1*m1*m1/denom
  c02  = c01 + 0.5d0*(m2-m1)/m3
  m12 = m1 + m2
  if (m12 <= m3) then
     c03 = 0.5d0*m12/m3
  else
     numer = m3*m3*(3.d0*m12-m3) + m1*m1*(m1-3.d0*m3) + m2*m2*(m2-3.d0*m3)
     c03 = numer/denom
  endif
  c00 = 3.d0*DSQRT(m1*m1 + m2*m2 + m3*m3)/denom

! 1: C<=C1; 2: C1<=C<=C2; 3: C2<=C<=C3; 4: C3<=C<=1/2 (a: m12<=m3; b: m3<m12)) 
  if (cch <= c01) then
     al = (denom*cch)**athird                                         
     AREA3D = c00*al*al                                               ! case (1)
  else if (cch <= c02) then 
    al = 0.5d0*(m1 + DSQRT(m1*m1 + 8.d0*m2*m3*(cch - c01)))           
    AREA3D = c00*m1*(2.d0*al-m1)                                      ! case (2)
  else if (cch <= c03) then
     p = 2.d0*m1*m2                                                   
     q = 1.5d0*m1*m2*(m12 - 2.d0*m3*cch)
     pst = DSQRT(p)
     arc = athird*DACOS(q/(p*pst))
     csarc = DCOS(arc)
     al = pst*(DSQRT(3.d0*(1.d0-csarc*csarc)) - csarc) + m12
     AREA3D = c00*(-al*al + m1*(2.d0*al-m1) + m2*(2.d0*al-m2))        ! case (3)
  else if (m12 <= m3) then                                  
     al = m3*cch + 0.5d0*m12                                         
     AREA3D = c00*2.d0*m1*m2                                         ! case (4a) 
  else                                  
     p = m12*m3 + m1*m2 - 0.25d0                                     
     q = 1.5d0*m1*m2*m3*(0.5d0-cch)
     pst = DSQRT(p)
     arc = athird*DACOS(q/(p*pst))
     csarc = DCOS(arc)
     al = pst*(DSQRT(3.d0*(1.d0-csarc*csarc)) - csarc) + 0.5d0 
     AREA3D = c00*(2.d0*al*(1.d0-al) - (m1*m1 + m2*m2 + m3*m3))      ! case (4b)
  endif

END FUNCTION AREA3D

!===============================================================================
! DESCRIPTION OF FUNCTION CENT3D:
! compute in the unit cube the coordinates of the centroid of the area 
! cut by the plane
! [nr]*[x] = nr1*x1 + nr2*x2 + nr3*x3 = alpha
! INPUT: normal coefficients nr1, nr2, nr3 and volume fraction cc 
! NOTE : the normal coefficients MUST satisfy the relation |nr1|+|nr2|+|nr3|=1;
!        the centroid coordinates change with the reflections required for 
!        negative normal coefficients, the midpoint symmetry and the ordering 
!        of the coefficients    
!-------------------------------------------------------------------------------
FUNCTION CENT3D(nr,cc)

  IMPLICIT NONE
  REAL(8), INTENT(IN):: nr(3),cc
  REAL(8), DIMENSION(3) :: CENT3D,ctd0
  REAL(8) :: cch,al,c01,c02,c03,np1,np2,np3
  REAL(8) :: m1,m2,m3,m12,numer,denom,p,pst,q,arc,csarc,top,bot
  REAL(8), PARAMETER :: athird=1.d0/3.d0 
  INTEGER, DIMENSION(3) :: index
  INTRINSIC DMAX1,DMIN1,DSQRT,DACOS,DCOS

  np1 = DABS(nr(1))                                 ! need positive coefficients
  np2 = DABS(nr(2))
  np3 = DABS(nr(3))
! order positive coefficients and get index ordering
  if (np1 <= np2) then                             
    m1 = np1
    m3 = np2
    index(1) = 1                                            
    index(3) = 2
  else
    m1 = np2
    m3 = np1
    index(1) = 2
    index(3) = 1
  endif
  
  if (np3 < m1) then
    m2 = m1
    m1 = np3
    index(2) = index(1)
    index(1) = 3
  else if (np3 >= m3) then
    m2 = m3
    m3 = np3
    index(2) = index(3)
    index(3) = 3
  else
    m2 = np3
    index(2) = 3
  endif
  cch = DMIN1(cc,1.d0-cc)                              ! symmetry: 0 < cch < 1/2
  denom = DMAX1(6.d0*m1*m2*m3,1.d-50)                           ! get cch ranges
  c01 = m1*m1*m1/denom
  c02  = c01 + 0.5d0*(m2-m1)/m3
  m12 = m1 + m2
  if (m12 <= m3) then
    c03 = 0.5d0*m12/m3
  else
    numer = m3*m3*(3.d0*m12-m3) + m1*m1*(m1-3.d0*m3) + m2*m2*(m2-3.d0*m3)
    c03 = numer/denom
  endif
  
! 1: C<=C1; 2: C1<=C<=C2; 3: C2<=C<=C3; 4: C3<=C<=1/2 (a: m12<=m3; b: m3<m12)) 
  if (cch <= c01) then
    al = (denom*cch)**athird 
    ctd0(1) = athird*al/m1
    ctd0(2) = athird*al/m2
    ctd0(3) = athird*al/m3                                            ! case (1)
  else if (cch <= c02) then 
    al = 0.5d0*(m1 + DSQRT(m1*m1 + 8.d0*m2*m3*(cch - c01)))    
    top = m1*m1 + 3.*al*(al-m1)
    bot = 3.*(2.*al-m1)
    ctd0(1) = (3.*al-2.*m1)/bot
    ctd0(2) = top/(m2*bot)
    ctd0(3) = top/(m3*bot)                                            ! case (2)
  else if (cch <= c03) then
    p = 2.d0*m1*m2                                                   
    q = 1.5d0*m1*m2*(m12 - 2.d0*m3*cch)
    pst = DSQRT(p)
    arc = athird*DACOS(q/(p*pst))
    csarc = DCOS(arc)
    al = pst*(DSQRT(3.d0*(1.d0-csarc*csarc)) - csarc) + m12
    bot = 3.*((al-m2)*(al-m2) + (al-m1)*(al-m1) - al*al)
    top = (al-m2)*(al-m2)*(al-m2) + m1*m1*(2.*m1-3.*al)
    ctd0(1) = top/(m1*bot)
    top = (al-m1)*(al-m1)*(al-m1) + m2*m2*(2.*m2-3.*al)
    ctd0(2) = top/(m2*bot)
!   top = (al-m1)*(al-m1)*(al-m1) + (al-m2)*(al-m2)*(al-m2) - al*al*al
!   ctd0(3) = top/(m3*bot)                                           
    ctd0(3) = (al - m1*ctd0(1) - m2*ctd0(2))/m3                       ! case (3)
  else if (m12 <= m3) then                                  
    al = m3*cch + 0.5d0*m12                                         
    ctd0(1) = 0.5d0
    ctd0(2) = 0.5d0
    ctd0(3) = (al-0.5d0*m12)/m3                                      ! case (4a) 
  else                                  
    p = m12*m3 + m1*m2 - 0.25d0                                     
    q = 1.5d0*m1*m2*m3*(0.5d0-cch)
    pst = DSQRT(p)
    arc = athird*DACOS(q/(p*pst))
    csarc = DCOS(arc)
    al = pst*(DSQRT(3.d0*(1.d0-csarc*csarc)) - csarc) + 0.5d0 
    bot = 2.d0*al*(1.d0-al) - (m1*m1 + m2*m2 + m3*m3)
    top = (m1*m1*m1 + m2*m2*m2 + m3*m3*m3) -3.*al*(m1*m1 + m2*m2 + m3*m3) 
    top = athird*(top + 3.*al*al - 2.*al*al*al)
    ctd0(1) = (top/m1 - (al-m1)*(al-m1))/bot
    ctd0(2) = (top/m2 - (al-m2)*(al-m2))/bot
!   ctd0(3) = (top/m3 - (al-m3)*(al-m3))/bot                         
    ctd0(3) = (al - m1*ctd0(1) - m2*ctd0(2))/m3                      ! case (4b)
  endif

! back from symmetry
  if (cc > 0.5d0) then
    ctd0(1) = 1.d0 - ctd0(1)
    ctd0(2) = 1.d0 - ctd0(2)
    ctd0(3) = 1.d0 - ctd0(3)
  endif

! get correct indexing
  CENT3D(index(1)) = ctd0(1)                              
  CENT3D(index(2)) = ctd0(2)
  CENT3D(index(3)) = ctd0(3)

! take care of negative coefficients
  if (nr(1) < 0.d0) CENT3D(1) = 1.d0 - CENT3D(1) 
  if (nr(2) < 0.d0) CENT3D(2) = 1.d0 - CENT3D(2)
  if (nr(3) < 0.d0) CENT3D(3) = 1.d0 - CENT3D(3)

END FUNCTION CENT3D
