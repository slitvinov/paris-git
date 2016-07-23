!================================================================================================
!=================================================================================================
! Paris-0.1
!
! Boiling extensions
! written by Leon Malan
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
!=================================================================================================
Module module_boil
use module_grid
implicit none

!variables
real(8), allocatable, dimension(:,:,:) :: Te, dTe, Te_old, s_v, kc, kc_old, Cp, Cp_old
real(8), allocatable, dimension(:,:,:) :: energy
real(8), allocatable, dimension(:,:,:,:) :: dPg, dug
real(8) :: Cp1, Cp2, h_fg, kc1, kc2, T_sat, T0_1, T0_2
real(8) :: mdot, qdot
real(8) :: ec1, ec2 ! Thermal expansion coeffs for Boussinesq
real(8) :: BC_T(1:6) 
integer :: BDRY_T(1:6)

contains

  !Read parameters & initialize
  subroutine ReadHeatParameters
    use module_VOF
    use module_flow
    implicit none
    integer :: i,j,k,ierr,BC
    namelist /heatparameters/ Cp1, Cp2, h_fg, kc1, kc2, T_sat, T0_1, T0_2,&
         BDRY_T, BC_T, ec1, ec2

    Cp1 = 1.0d0; Cp2 = 10.0; h_fg=1.0d30; 
    kc1=1.0d0; kc2=1.0d3; T_sat=1.0d0
    T0_1=1.0d1; T0_2=1.0d1
    BDRY_T = 0
    BC_T = 1.0d0

    open(unit=10, file='inputheat', status='old', action='read', iostat=ierr)
    if (ierr .ne. 0) call err_no_out_dir("ReadHeatParameters: error opening 'inputheat' file --- perhaps it does not exist ?")
    read(10,NML=heatparameters)
    close(10)

    do BC=1,6
       if (.not.(BDRY_T(BC)==0 .or. BDRY_T(BC)==1)) then
          write(*,'("Error, Temp BC must be either 0 or 1")')
          call pariserror("Invalid temperature BC")
       endif
    enddo

    ! Allocate variables
    allocate (Te(imin:imax,jmin:jmax,kmin:kmax),dTe(imin:imax,jmin:jmax,kmin:kmax),&
         s_v(imin:imax,jmin:jmax,kmin:kmax),kc(imin:imax,jmin:jmax,kmin:kmax),&
         Cp(imin:imax,jmin:jmax,kmin:kmax), energy(imin:imax,jmin:jmax,kmin:kmax),&
         dPg(imin:imax,jmin:jmax,kmin:kmax,1:3),dug(imin:imax,jmin:jmax,kmin:kmax,1:3) )

    if (itime_scheme==2) then
       allocate( Te_old(imin:imax,jmin:jmax,kmin:kmax), kc_old(imin:imax,jmin:jmax,kmin:kmax),&
       Cp_old(imin:imax,jmin:jmax,kmin:kmax))
    endif

    !Initcondition
    if (TwoPhase) then
       if (DoVof) then
          do i=is,ie; do j=js,je;  do k=ks,ke
             Te(i,j,k) = cvof(i,j,k)*T0_2 + (1.0d0-cvof(i,j,k))*T0_1
          enddo; enddo; enddo
       endif
    else
       Te = T0_1
    endif

  end subroutine ReadHeatParameters

!get volume source at interface
!get heat source at interface
!solve energy, conjugate heat transfer

  !Heat diffusion
  subroutine TempDiffusion
    use module_flow
    implicit none
    integer :: i,j,k
    do i=is,ie; do j=js,je;  do k=ks,ke
       dTe(i,j,k) = dTe(i,j,k) + 1.d0/(rho(i,j,k)*cp(i,j,k))*&
            ( (kc(i+1,j,k)+kc(i,j,k))*(Te(i+1,j,k)-Te(i,j,k)) -&
            (kc(i,j,k)+kc(i-1,j,k))*(Te(i,j,k)-Te(i-1,j,k)) )/( dxh(i)*dx(i)*2.0d0 ) + &
            ( (kc(i,j+1,k)+kc(i,j,k))*(Te(i,j+1,k)-Te(i,j,k)) -&
            (kc(i,j,k)+kc(i,j-1,k))*(Te(i,j,k)-Te(i,j-1,k)) )/( dyh(j)*dy(j)*2.0d0 ) + &
            ( (kc(i,j,k+1)+kc(i,j,k))*(Te(i,j,k+1)-Te(i,j,k)) -&
            (kc(i,j,k)+kc(i,j,k-1))*(Te(i,j,k)-Te(i,j,k-1)) )/( dzh(k)*dz(k)*2.0d0 )
    enddo; enddo; enddo
  end subroutine TempDiffusion

  subroutine get_heat_source
    implicit none
    !Calculate temp gradients at interface in normal direction
    qdot=1.0d2
    mdot=qdot/h_fg
  end subroutine get_heat_source

  subroutine vofandenergysweeps(tswap)
    use module_VOF
    use module_flow
    use module_tmpvar
    implicit none
    integer, dimension(imin:imax,jmin:jmax,kmin:kmax) :: dummy_flag 
    integer, intent(in) :: tswap
    integer :: i,j,k

    !get cell energies, energies at s-1 and e+1 needed for BC's
    do i=is-1,ie+1; do j=js-1,je+1;  do k=ks-1,ke+1
       energy(i,j,k) = rho(i,j,k)*cp(i,j,k)*Te(i,j,k)
    enddo; enddo; enddo

    tmp = cvof
    dummy_flag = vof_flag

    if (MOD(tswap,3).eq.0) then  ! do z then x then y 
       call swpz_energy(w,tmp,3)
       call swp(w,tmp,dummy_flag,3,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))

       call swpz_energy(u,tmp,1)
       call swp(u,tmp,dummy_flag,1,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))

       call swpz_energy(v,tmp,2)
       call swp(v,tmp,dummy_flag,2,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))
    elseif (MOD(tswap,3).eq.1) then ! do y z x
       call swpz_energy(v,tmp,2)
       call swp(v,tmp,dummy_flag,2,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))

       call swpz_energy(w,tmp,3)
       call swp(w,tmp,dummy_flag,3,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))

       call swpz_energy(u,tmp,1)
       call swp(u,tmp,dummy_flag,1,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))
    else ! do x y z
       call swpz_energy(u,tmp,1)
       call swpz(u,tmp,dummy_flag,1,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))

       call swpz_energy(v,tmp,2)
       call swpz(v,tmp,dummy_flag,2,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))

       call swpz_energy(w,tmp,3)
       call swpz(w,tmp,dummy_flag,3,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))
    endif

    call get_dT_from_energy(tmp)
  end subroutine vofandenergysweeps

  subroutine swpz_energy(us,c,d)
    use module_VOF
    use module_flow
    use module_BC
    implicit none
    logical error
    integer i,j,k
    integer i0,j0,k0
    integer i1,j1,k1
    integer, intent(in) :: d
    real(8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: us
    real(8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: c
    real(8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax) :: en1,en2,en3
    real(8) mm1,mm2,vof
    real(8) a1,a2,alpha,T_adv
    real(8) rhocp_ll, rhocp_l, rhocp_r, rhocp_rr, rhocp_c
    real(8) T_l, T_r,d_int
    REAL(8) deltax(3),x0(3),fl3d
    real(8) nr(3), stencil3x3(-1:1,-1:1,-1:1)
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
                d_int=dx(i)/2.d0
             elseif (d.eq.2) then
                a2 = us(i,j,k)*dt/dyh(j)
                a1 = us(i,j-1,k)*dt/dyh(j-1)
                d_int=dy(j)/2.d0
             elseif (d.eq.3) then
                a2 = us(i,j,k)*dt/dzh(k)
                a1 = us(i,j,k-1)*dt/dzh(k-1)
                d_int=dz(k)/2.d0
             endif
             
             ! Get new density and heat capacity after vofsweep
             rhocp_ll=(rho2*cp2*c(i-i0*2,j-j0*2,k-k0*2)+rho1*cp1*(1.d0-c(i-i0*2,j-j0*2,k-k0*2)))
             rhocp_l=(rho2*cp2*c(i-i0,j-j0,k-k0)+rho1*cp1*(1.d0-c(i-i0,j-j0,k-k0)))
             rhocp_c=(rho2*cp2*c(i   ,j   ,k   )+rho1*cp1*(1.d0-c(i   ,j   ,k   )))
             rhocp_r=(rho2*cp2*c(i+i0,j+j0,k+k0)+rho1*cp1*(1.d0-c(i+i0,j+j0,k+k0)))
             rhocp_rr=(rho2*cp2*c(i+i0*2,j+j0*2,k+k0*2)+rho1*cp1*(1.d0-c(i+i0*2,j+j0*2,k+k0*2)))
             
             ! Use new cell temperature after previous energy and vof sweep
             T_adv = energy(i,j,k)/(rhocp_c)

             ! Interpolate for T-1/2 (T_l) and T+1/2 (T_r) using chosen AdvectionScheme
             ! For UpWind, T_l=T_r=T_adv
             if (.not.(i==is-1 .or. i==ie+1 .or. j==js-1 .or. j==je+1 .or. k==ks-1 .or. k==ke+1)) then
                if (us(i-i0,j-j0,k-k0)>0.d0) then
                   T_l=interpole3(energy(i-i0*2,j-j0*2,k-k0*2)/(rhocp_ll),energy(i-i0,j-j0,k-k0)/(rhocp_l),&
                        T_adv,AdvectionScheme,d_int)
                else
                   T_l=interpole3(energy(i+i0,j+j0,k+k0)/(rhocp_r),T_adv,energy(i-i0,j-j0,k-k0)/(rhocp_l),&
                        AdvectionScheme,d_int)
                endif

                if (us(i,j,k)>0.d0) then
                   T_r=interpole3(energy(i-i0,j-j0,k-k0)/(rhocp_l),T_adv,energy(i+i0,j+j0,k+k0)/(rhocp_r),&
                        AdvectionScheme,d_int)
                else
                   T_r=interpole3(energy(i+i0*2,j+j0*2,k+k0*2)/(rhocp_rr),energy(i+i0,j+j0,k+k0)/(rhocp_r),&
                        T_adv,AdvectionScheme,d_int)
                endif
             else
                T_l=T_adv
                T_r=T_adv
             endif

             ! Energy for full cells
             mm1 = dmax1(a1,0.0d0)
             mm2 = 1.d0 - mm1 + dmin1(0.d0,a2)

             if (T_adv > 0.d0) then
                en1(i,j,k) = dmax1(-a1,0.d0)*energy(i,j,k)/T_adv*T_l
                en3(i,j,k) = dmax1( a2,0.d0)*energy(i,j,k)/T_adv*T_r
                en2(i,j,k) = energy(i,j,k) - en1(i,j,k) - en3(i,j,k)
             else
                en1(i,j,k) = dmax1(-a1,0.d0)*energy(i,j,k)
                en3(i,j,k) = dmax1( a2,0.d0)*energy(i,j,k)
                en2(i,j,k) = energy(i,j,k) - en1(i,j,k) - en3(i,j,k)
             endif

             if ((c(i,j,k) .gt. 0.d0).and.(c(i,j,k) .lt. 1.d0)) then
                do i1=-1,1; do j1=-1,1; do k1=-1,1
                   stencil3x3(i1,j1,k1) = c(i+i1,j+j1,k+k1)
                enddo;enddo;enddo
                call fit_plane_new(c(i,j,k),d,a1,a2,stencil3x3,nr,alpha,error)
                if(error) then
                   !write(*,*)"WARNING: new plane error!"
                   cycle
                endif
                x0=0.d0
                deltax=1.d0
                if(a1<0.d0) then
                   x0(d)=a1
                   deltax(d)=-a1
                   vof = fl3d(nr,alpha,x0,deltax)
                   en1(i,j,k) = ( rho2*cp2*vof+rho1*cp1*(-a1-vof) )*T_l
                   !if (k==3 .and. i==20) write(*,'("Cut cell en1: ",e14.5)')en1(i,j,k)
                endif
                if(a2>0.d0) then
                   x0(d)=1.d0
                   deltax(d)=a2
                   vof = fl3d(nr,alpha,x0,deltax)
                   en3(i,j,k) = ( rho2*cp2*vof+rho1*cp1*(a2-vof) )*T_r
                   !if (k==3 .and. i==20) then
                   !   write(*,'("Cut cell en3, vof, a2: ",3e14.5)')en3(i,j,k),vof,a2
                   !   write(*,'("en3 alpha, normals",4e14.5)')alpha,nr(1:3)
                   !endif
                endif
                x0(d)=mm1
                deltax(d)=mm2
                vof = fl3d(nr,alpha,x0,deltax)
                en2(i,j,k) = ( rho2*cp2*vof+rho1*cp1*(mm2-vof) )*T_adv
                !if (k==3 .and. i==20) then
                !   write(*,'("Cut cell en2, vof, a1,a2: ",4e14.5)')en2(i,j,k),vof,a1,a2
                !   write(*,'("en2 alpha, normals",4e14.5)')alpha,nr(1:3)
                !endif
             endif
          enddo
       enddo
    enddo

    do k=ks,ke
       do j=js,je
          do i=is,ie
             energy(i,j,k)  = en1(i+i0,j+j0,k+k0)+en2(i,j,k)+en3(i-i0,j-j0,k-k0)
             !if (k==3 .and. i==20 .and. c(i,j,k)<1.d0 .and. c(i,j,k)>0.d0 ) then
             !   write(*,'("Cut cell energy, en1 right, en2, en3 left: ",4e14.5)')&
             !        energy(i,j,k),en1(i+i0,j+j0,k+k0),en2(i,j,k),en3(i-i0,j-j0,k-k0)
             !endif
          enddo
       enddo
    enddo

    call do_all_ghost(energy)

  end subroutine swpz_energy

  subroutine get_dT_from_energy(c)
    use module_flow
    implicit none
    real(8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: c
    integer :: i,j,k

    do k=ks,ke
       do j=js,je
          do i=is,ie
             dTe(i,j,k) = dTe(i,j,k) + 1.d0/dt*(energy(i,j,k)/&
                  ((rho2*c(i,j,k)+rho1*(1.d0-c(i,j,k)))*(cp2*c(i,j,k)+cp1*(1.d0-c(i,j,k)))) - Te(i,j,k))
          enddo
       enddo
    enddo

  end subroutine get_dT_from_energy

  !Boundary values for temperature
  subroutine SetTempBc
    implicit none
    integer :: i,j,k

    !BC on x-
    if (coords(1)==0) then
       do j=js,je; do k=ks,ke
          if (BDRY_T(1)==0) then !Dirichlet
             Te(is-1,j,k)=2.0d0*BC_T(1)-Te(is,j,k)
          else !Neumann
             Te(is-1,j,k)=Te(is,j,k)-BC_T(1)*dxh(is-1)
          endif
       enddo; enddo
    endif

    !BC on x+
    if (coords(1)==nPx-1) then
       do j=js,je; do k=ks,ke
          if (BDRY_T(4)==0) then 
             Te(ie+1,j,k)=2.0d0*BC_T(4)-Te(ie,j,k)
          else
             Te(ie+1,j,k)=Te(ie,j,k)+BC_T(4)*dxh(ie)
          endif
       enddo; enddo
    endif

    !BC on y-
    if (coords(2)==0) then
       do i=is,ie; do k=ks,ke
          if (BDRY_T(2)==0) then !Dirichlet
             Te(i,js-1,k)=2.0d0*BC_T(2)-Te(i,js,k)
          else !Neumann
             Te(i,js-1,k)=Te(i,js,k)-BC_T(2)*dyh(js-1)
          endif
       enddo; enddo
    endif

    !BC on y+
    if (coords(2)==nPy-1) then
       do i=is,ie; do k=ks,ke
          if (BDRY_T(5)==0) then 
             Te(i,je+1,k)=2.0d0*BC_T(5)-Te(i,je,k)
          else
             Te(i,je+1,k)=Te(i,je,k)+BC_T(5)*dyh(je)
          endif
       enddo; enddo
    endif

    !BC on z-
    if (coords(3)==0) then
       do i=is,ie; do j=js,je
          if (BDRY_T(3)==0) then !Dirichlet
             Te(i,j,ks-1)=2.0d0*BC_T(3)-Te(i,j,ks)
          else !Neumann
             Te(i,j,ks-1)=Te(i,j,ks)-BC_T(3)*dzh(ks-1)
          endif
       enddo; enddo
    endif

    !BC on z+
    if (coords(3)==nPz-1) then
       do i=is,ie; do j=js,je
          if (BDRY_T(6)==0) then 
             Te(i,j,ke+1)=2.0d0*BC_T(6)-Te(i,j,ke)
          else
             Te(i,j,ke+1)=Te(i,j,ke)+BC_T(6)*dzh(ke)
          endif
       enddo; enddo
    endif 

  end subroutine SetTempBC

  !microlayer

end Module module_boil
