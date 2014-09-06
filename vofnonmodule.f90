!=================================================================================================
!=================================================================================================
! PARIS  Parallel Robust Interface Simulator 
! OLD VOF FUNCTIONS
!
! authors Stephane Zaleski   zaleski@dalembert.upmc.fr
!         Ruben Scardovelli  ruben.scardovelli@unibo.it
!         Daniel Fuster      dfuster@gmail.com
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
! 
!-------------------------------------------------------------------------------------------------
subroutine swp(us,c,f,d,vof1,vof2,vof3)
  use module_vof
  implicit none
  integer, intent(in) :: d
  real (8)  , dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: us
  real (8)  , dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: c,vof1,vof2,vof3
  integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: f
  if (VOF_advect=='Dick_Yue') then  ! Yue-Weymouth = Eulerian Implicit + central cell stuff
     call swpr(us,c,f,d,vof1,vof2,vof3)  
  elseif (VOF_advect=='CIAM') then  ! CIAM == Lagrangian Explicit
     call swpz(us,c,f,d,vof1,vof2,vof3)
  else
     call pariserror("*** unknown vof scheme")
  endif

  call get_flags_and_clip(c,f)
  call do_all_ghost(c)
  call do_all_ighost(f)
  call get_vof_phase(c)

end subroutine swp

subroutine swp_stg(us,c,f,d,vof1,vof2,vof3,dir)
  use module_vof
  implicit none
  integer, intent(in) :: d,dir
  real (8)  , dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: us
  real (8)  , dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: c,vof1,vof2,vof3
  integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: f
  if (VOF_advect=='Dick_Yue') then  ! Yue-Weymouth = Eulerian Implicit + central cell stuff
     call pariserror("*** not implemented yet")
  elseif (VOF_advect=='CIAM') then  ! CIAM == Lagrangian Explicit
     call swpz_stg(us,c,f,d,vof1,vof2,vof3,dir)
  else
     call pariserror("*** unknown vof scheme")
  endif

  call get_flags_and_clip(c,f)
  call do_all_ghost(c)
  call do_all_ighost(f)
  call get_vof_phase(c)

end subroutine swp_stg


subroutine swpmom(us,c,d,mom1,mom2,mom3,mom)
  use module_vof
  implicit none
  integer, intent(in) :: d
  real(8)  , dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: us,mom
  real(8)  , dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: c
  real(8)  , dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: mom1,mom2,mom3

  if (VOF_advect=='Dick_Yue') then  ! Yue-Weymouth = Eulerian Implicit + central cell stuff
     call swprmom(us,c,d,mom1,mom2,mom3,mom)
  elseif (VOF_advect=='CIAM') then  ! CIAM == Lagrangian Explicit
     call swpzmom(us,c,d,mom1,mom2,mom3,mom)
  else
     call pariserror("*** unknown vof scheme")
  endif

  call do_all_ghost(mom)

end subroutine swpmom

subroutine swpmom_stg(us,c,d,mom1,mom2,mom3,mom,dir)
  use module_vof
  implicit none
  integer, intent(in) :: d,dir
  real(8)  , dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: us,mom
  real(8)  , dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: c
  real(8)  , dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: mom1,mom2,mom3

  if (VOF_advect=='Dick_Yue') then  ! Yue-Weymouth = Eulerian Implicit + central cell stuff
     call pariserror("*** not implemented yet")
  elseif (VOF_advect=='CIAM') then  ! CIAM == Lagrangian Explicit
     call swpzmom_stg(us,c,d,mom1,mom2,mom3,mom,dir)
  else
     call pariserror("*** unknown vof scheme")
  endif

  call do_all_ghost(mom)

end subroutine swpmom_stg

!
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
  real(8) a1,a2,alpha,fl3dnew, stencil3x3(-1:1,-1:1,-1:1)
  real(8) nr(3),deltax(3),x0(3)
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

           ! @@@ we need to introduce full/empty flags

           if (c(i,j,k) .eq. 1.0d0) then
              vof1(i,j,k) = dmax1(-a1,0.d0)
              vof2(i,j,k) = 1.d0 - dmax1(a1,0.d0) + dmin1(a2,0.d0)
              vof3(i,j,k) = dmax1(a2,0.d0)
           else if (c(i,j,k) .gt. 0.d0) then
              do i1=-1,1; do j1=-1,1; do k1=-1,1
                stencil3x3(i1,j1,k1) = c(i+i1,j+j1,k+k1)
              enddo;enddo;enddo
              call fit_plane_new(c(i,j,k),d,a1,a2,stencil3x3,nr,alpha,error)
              if(error) cycle
              mm1 = dmax1(a1,0.0d0)
              mm2 = 1.d0 - mm1 + dmin1(0.d0,a2)
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

subroutine swpz_stg(us,c,f,d,vof1,vof2,vof3,dir)
  !***
  use module_grid
  use module_flow
  use module_vof
  implicit none
  logical error
  integer i,j,k
  integer i0,j0,k0
  integer i1,j1,k1
  integer i2,j2,k2,dir
  integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: f
  integer, intent(in) :: d
  real (8)  , dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: us
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: c,vof1,vof2,vof3
  real(8) mm1,mm2
  real(8) a1,a2,alpha,fl3dnew, stencil3x3(-1:1,-1:1,-1:1)
  real(8) nr(3),deltax(3),x0(3)
  intrinsic dmax1,dmin1
  !***
  call init_i0j0k0 (d,i0,j0,k0)
  call init_i0j0k0 (dir,i2,j2,k2)

  if(ng.lt.2) call pariserror("wrong ng")
  do k=ks-1,ke+1
     do j=js-1,je+1
        do i=is-1,ie+1
         if (d.eq.1) then
              a2 = 0.5d0*(us(i,j,k)+us(i+i2,j+j2,k+k2))*dt/dxh(i)
              a1 = 0.5d0*(us(i-1,j,k)+us(i-1+i2,j+j2,k+k2))*dt/dxh(i-1)
           elseif (d.eq.2) then
              a2 = 0.5d0*(us(i,j,k)+us(i+i2,j+j2,k+k2))*dt/dyh(j)
              a1 = 0.5d0*(us(i,j-1,k)+us(i+i2,j-1+j2,k+k2))*dt/dyh(j-1)
           elseif (d.eq.3) then
              a2 = 0.5d0*(us(i,j,k)+us(i+i2,j+j2,k+k2))*dt/dzh(k)
              a1 = 0.5d0*(us(i,j,k-1)+us(i+i2,j+j2,k-1+k2))*dt/dzh(k-1)
           endif
           !***
           !     3 cases: 1: default (c=0. and fluxes=0.); 2: c=1.; 3:c>0.
           !***
           vof1(i,j,k) = 0.0d0
           vof2(i,j,k) = 0.0d0
           vof3(i,j,k) = 0.0d0

           ! @@@ we need to introduce full/empty flags

           if (c(i,j,k) .eq. 1.0d0) then
              vof1(i,j,k) = dmax1(-a1,0.d0)
              vof2(i,j,k) = 1.d0 - dmax1(a1,0.d0) + dmin1(a2,0.d0)
              vof3(i,j,k) = dmax1(a2,0.d0)
           else if (c(i,j,k) .gt. 0.d0) then
              do i1=-1,1; do j1=-1,1; do k1=-1,1
                stencil3x3(i1,j1,k1) = c(i+i1,j+j1,k+k1)
              enddo;enddo;enddo
              call fit_plane_new(c(i,j,k),d,a1,a2,stencil3x3,nr,alpha,error)
              if(error) cycle
              mm1 = dmax1(a1,0.0d0)
              mm2 = 1.d0 - mm1 + dmin1(0.d0,a2)
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
end subroutine swpz_stg

subroutine swpzmom(us,c,d,mom1,mom2,mom3,mom)
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
  integer, intent(in) :: d
  real(8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: us
  real(8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: mom
  real(8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: c
  real(8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: mom1,mom2,mom3
  real(8) mm1,mm2,vof
  real(8) a1,a2,alpha,uavg,rhoavg
  REAL(8) deltax(3),x0(3),fl3dnew
  real(8) mxyz(3), nr(3), stencil3x3(-1:1,-1:1,-1:1)
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
!           if ((us(i,j,k)+us(i-i0,j-j0,k-k0)).gt.0.d0) then
!               uavg = us(i-i0,j-j0,k-k0)
!           else
!               uavg = us(i,j,k)
!           endif

           if ((c(i,j,k) .gt. 0.d0).and.(c(i,j,k) .lt. 1.d0)) then
              do i1=-1,1; do j1=-1,1; do k1=-1,1
                 stencil3x3(i1,j1,k1) = c(i+i1,j+j1,k+k1)
              enddo;enddo;enddo
              call fit_plane_new(c(i,j,k),d,a1,a2,stencil3x3,nr,alpha,error)
              if(error) cycle
              x0=0d0
              deltax=1d0
              if(a1.lt.0d0) then
                 x0(d)=a1
                 deltax(d)=-a1
                 vof = fl3dnew(nr,alpha,x0,deltax)
                 mom1(i,j,k) = (rho2*vof + rho1*(-a1 - vof))*uavg
              endif
              if(a2.gt.0d0) then
                 x0(d)=1d0
                 deltax(d)=a2
                 vof = fl3dnew(nr,alpha,x0,deltax)
                 mom3(i,j,k) = (rho2*vof + rho1*(a2 - vof))*uavg
              endif
              x0(d)=mm1
              deltax(d)=mm2
              vof = fl3dnew(nr,alpha,x0,deltax)
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

subroutine swpzmom_stg(us,c,d,mom1,mom2,mom3,mom,dir)
!  !***
  use module_grid
  use module_flow
  use module_vof
  use module_BC
  implicit none
  logical error
  integer i,j,k,dir
  integer i0,j0,k0
  integer i1,j1,k1
  integer i2,j2,k2
  integer inv(3)
  integer, intent(in) :: d
  real(8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: us
  real(8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: mom
  real(8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: c
  real(8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: mom1,mom2,mom3
  real(8) mm1,mm2,vof
  real(8) a1,a2,alpha,uavg,rhoavg
  REAL(8) deltax(3),x0(3),fl3dnew
  real(8) mxyz(3), nr(3), stencil3x3(-1:1,-1:1,-1:1)
  intrinsic dmax1,dmin1
  !***
  call init_i0j0k0 (d,i0,j0,k0)
  call init_i0j0k0 (dir,i2,j2,k2)

  if(ng.lt.2) call pariserror("wrong ng")
  do k=ks-1,ke+1
     do j=js-1,je+1
        do i=is-1,ie+1
           if (d.eq.1) then
              a2 = 0.5d0*(us(i,j,k)+us(i+i2,j+j2,k+k2))*dt/dxh(i)
              a1 = 0.5d0*(us(i-1,j,k)+us(i-1+i2,j+j2,k+k2))*dt/dxh(i-1)
           elseif (d.eq.2) then
              a2 = 0.5d0*(us(i,j,k)+us(i+i2,j+j2,k+k2))*dt/dyh(j)
              a1 = 0.5d0*(us(i,j-1,k)+us(i+i2,j-1+j2,k+k2))*dt/dyh(j-1)
           elseif (d.eq.3) then
              a2 = 0.5d0*(us(i,j,k)+us(i+i2,j+j2,k+k2))*dt/dzh(k)
              a1 = 0.5d0*(us(i,j,k-1)+us(i+i2,j+j2,k+k2))*dt/dzh(k-1)
           endif
            
           !momentum for full cells
           mm1 = dmax1(a1,0.0d0)
           mm2 = 1.d0 - mm1 + dmin1(0.d0,a2)

           mom1(i,j,k)  = dmax1(-a1,0.d0)*mom(i,j,k)
           mom3(i,j,k)  = dmax1(a2,0.d0) *mom(i,j,k)
           mom2(i,j,k)  = mom(i,j,k) - mom1(i,j,k) - mom3(i,j,k)
           uavg = mom(i,j,k)/(rho2*c(i,j,k)+rho1*(1.d0-c(i,j,k)))
!           if ((us(i,j,k)+us(i-i0,j-j0,k-k0)).gt.0.d0) then
!               uavg = us(i-i0,j-j0,k-k0)
!           else
!               uavg = us(i,j,k)
!           endif

           if ((c(i,j,k) .gt. 0.d0).and.(c(i,j,k) .lt. 1.d0)) then
              do i1=-1,1; do j1=-1,1; do k1=-1,1
                 stencil3x3(i1,j1,k1) = c(i+i1,j+j1,k+k1)
              enddo;enddo;enddo
              call fit_plane_new(c(i,j,k),d,a1,a2,stencil3x3,nr,alpha,error)
              if(error) cycle
              x0=0d0
              deltax=1d0
              if(a1.lt.0d0) then
                 x0(d)=a1
                 deltax(d)=-a1
                 vof = fl3dnew(nr,alpha,x0,deltax)
                 mom1(i,j,k) = (rho2*vof + rho1*(-a1 - vof))*uavg
              endif
              if(a2.gt.0d0) then
                 x0(d)=1d0
                 deltax(d)=a2
                 vof = fl3dnew(nr,alpha,x0,deltax)
                 mom3(i,j,k) = (rho2*vof + rho1*(a2 - vof))*uavg
              endif
              x0(d)=mm1
              deltax(d)=mm2
              vof = fl3dnew(nr,alpha,x0,deltax)
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
end subroutine swpzmom_stg

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
    REAL(8) :: a1,a2,alpha
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
  ! call test_cell_size()

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
              alpha = AL3DNEW(mxyz,c(i,j,k))
              ! Eulerian advection
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

SUBROUTINE swprmom(us,c,dir,mom1,cg,mom3,mom)
!***
    USE module_grid
    USE module_flow
    USE module_vof
    use module_BC
    
    IMPLICIT NONE
    include 'mpif.h'
    INTEGER :: i,j,k
    INTEGER :: invx,invy,invz,ii,jj,kk,i0,j0,k0
    INTEGER, INTENT(IN) :: dir
    REAL (8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), INTENT(INOUT) :: us
    real(8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: mom
    REAL (8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), INTENT(INOUT) :: c,mom1,cg,mom3
    REAL(8), TARGET :: dmx,dmy,dmz,dxyz
    REAL(8), POINTER :: dm1,dm2,dm3
    REAL(8) :: a1,a2,alpha,vof,uavg
    REAL(8) :: AL3DNEW, FL3DNEW, x0(3), deltax(3)
    real(8) :: mxyz(3),stencil3x3(-1:1,-1:1,-1:1)
    INTRINSIC DMAX1,DMIN1

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
           mom1(i,j,k)  = dmax1(-a1,0.d0)*mom(i,j,k)
           mom3(i,j,k)  = dmax1(a2,0.d0) *mom(i,j,k)
           uavg = mom(i,j,k)/(rho2*c(i,j,k)+rho1*(1.d0-c(i,j,k)))
           !  c = 1.
           if ((c(i,j,k) .gt. 0.d0).and.(c(i,j,k) .lt. 1.d0)) then
              ! local stencil and normal vector: |dmx|+|dmy|+|dmz| = 1.
              do i0=-1,1; do j0=-1,1; do k0=-1,1
                 stencil3x3(i0,j0,k0) = c(i+i0,j+j0,k+k0)
              enddo;enddo;enddo
              call mycs(stencil3x3,mxyz)
              alpha = AL3DNEW(mxyz,c(i,j,k))
              ! Eulerian advection
              x0=0d0
              deltax=1d0
              if(a1<0d0) then
                 deltax(dir)=-a1
                 vof = FL3DNEW(mxyz,alpha,x0,deltax)
                 mom1(i,j,k) = (rho2*vof + rho1*(-a1 - vof))*uavg
              endif
              if(a2>0d0) then
                 x0(dir)=1d0-a2
                 deltax(dir)=a2
                 vof = FL3DNEW(mxyz,alpha,x0,deltax)
                 mom3(i,j,k) = (rho2*vof + rho1*(a2 - vof))*uavg
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
           uavg = mom(i,j,k)/(rho2*c(i,j,k)+rho1*(1.d0-c(i,j,k)))
           mom(i,j,k)  = mom(i,j,k) -  (mom3(i,j,k) - mom1(i+ii,j+jj,k+kk)) + & 
                      (mom3(i-ii,j-jj,k-kk) - mom1(i,j,k)) &
                      + (rho2*cg(i,j,k)+rho1*(1.d0-cg(i,j,k)))*uavg*(a2-a1);
        enddo
     enddo
  enddo
  ! apply proper boundary conditions 
  call SetMomentumBC(us,c,mom,dir,umask,rho1,rho2) 
end subroutine swprmom

!=================================================================================================

subroutine ls2vof_in_cell(stencil3x3,c,nflag)
  implicit none
  real(8), intent(out):: c
  integer, intent(out):: nflag
  real(8) :: zero(3), one(3), norml1
  real(8) :: mx,my,mz,alpha
  real(8) :: fl3dnew
  real(8) :: mxyz(3),stencil3x3(-1:1,-1:1,-1:1)

  zero=0d0
  one=1d0
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
  mxyz = dabs(mxyz)
  norml1 = sum(mxyz)
  mxyz = mxyz/norml1
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
     c = fl3dnew(mxyz,alpha,zero,one)
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
subroutine youngs(c,mm)
  !***
  implicit none
  real(8), intent(in) :: c(0:2,0:2,0:2)
  real(8), intent(out) :: mm(0:2)
  integer :: i
  real(8) :: norm

  call fd32(c,mm)
  norm = sum(abs(mm)) ! + TINY
  mm = mm/norm

  return 
end subroutine youngs

