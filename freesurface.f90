!================================================================================================
!=================================================================================================
! Paris-0.1
!
! Free surface extensions
! written by Leon Malan and Stephane Zaleski
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
!-------------------------------------------------------------------------------------------------
 subroutine set_topology(vof_phase,iout)
  use module_grid
  use module_freesurface
  use module_IO
  implicit none
  integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: vof_phase
  integer :: i,j,k,level,iout
  logical, dimension(imin:imax,jmin:jmax,kmin:kmax) :: u_assigned, v_assigned, w_assigned
  !initialize pmask to 1d0 and all masks to -1 and unassigned
  !pmask = 1d0
  u_cmask = -1; v_cmask = -1; w_cmask = -1
  u_assigned = .false.; v_assigned=.false.; w_assigned=.false.
  debug = .true.
  if (debug) then
     Open(unit=20,FILE=TRIM(out_path)//'/Top_0-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout,padding))//'.txt')
     Open(unit=21,FILE=TRIM(out_path)//'/Top_1-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout,padding))//'.txt')
     Open(unit=22,FILE=TRIM(out_path)//'/Top_2-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout,padding))//'.txt')
     Open(unit=19,FILE=TRIM(out_path)//'/Pmask-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout,padding))//'.txt')
  endif
  !First loop to set pmask and level 0 velocities in liq-liq and liq-gas cells
  do k=ks,ke; do j=js,je; do i=is,ie
     if (vof_phase(i,j,k) == 1) then
        !pmask(i,j,k) = 0d0
        if (debug) write(19,13)x(i),y(j),z(k)
        if (vof_phase(i+1,j,k) == 0) then
           u_cmask(i,j,k,0) = 1; u_assigned(i,j,k) = .true.
        endif
        if (vof_phase(i,j+1,k) == 0) then
           v_cmask(i,j,k,0) = 1; v_assigned(i,j,k) = .true.
        endif
        if (vof_phase(i,j,k+1) == 0) then 
           w_cmask(i,j,k,0) = 1; w_assigned(i,j,k) = .true.
        endif
     endif
     if (vof_phase(i,j,k) == 0) then
        if (vof_phase(i+1,j,k) == 1) then
           u_cmask(i,j,k,0) = 1; u_assigned(i,j,k) = .true.
        endif
        if (vof_phase(i,j+1,k) == 1) then
           v_cmask(i,j,k,0) = 1; v_assigned(i,j,k) = .true.
        endif
        if (vof_phase(i,j,k+1) == 1) then 
           w_cmask(i,j,k,0) = 1; w_assigned(i,j,k) = .true.
        endif
        if (vof_phase(i+1,j,k) == 0) then
           u_cmask(i,j,k,0) = 1; u_assigned(i,j,k)=.true.
        endif
        if (vof_phase(i,j+1,k) == 0) then
           v_cmask(i,j,k,0) = 1; v_assigned(i,j,k)=.true.
        endif
        if (vof_phase(i,j,k+1) == 0) then 
           w_cmask(i,j,k,0) = 1; w_assigned(i,j,k)=.true.
        endif
     endif
  enddo; enddo; enddo
  !Set levels 1 to X_level
  do level=1,X_level
  do k=ks,ke; do j=js,je; do i=is,ie
     !Tests: in between gas nodes, neighbour level -1, unassigned
     !u-neighbours
     if (u_cmask(i,j,k,level-1)==1) then
        if (.not.u_assigned(i,j+1,k) .and. vof_phase(i,j+1,k)==1 .and. vof_phase(i+1,j+1,k)==1) then
           u_cmask(i,j+1,k,level) = 1; u_assigned(i,j+1,k)=.true.
        endif
        if (.not.u_assigned(i+1,j,k) .and. vof_phase(i+1,j,k)==1 .and. vof_phase(i+2,j,k)==1) then
           u_cmask(i+1,j,k,level) = 1; u_assigned(i+1,j,k) = .true.
        endif
        if (.not.u_assigned(i,j-1,k) .and. vof_phase(i,j-1,k)==1 .and. vof_phase(i+1,j-1,k)==1) then
           u_cmask(i,j-1,k,level) = 1; u_assigned(i,j-1,k)=.true.
        endif
        if (.not.u_assigned(i-1,j,k) .and. vof_phase(i-1,j,k)==1 .and. vof_phase(i,j,k)==1) then
           u_cmask(i-1,j,k,level) = 1; u_assigned(i-1,j,k) = .true.
        endif
        if (.not.u_assigned(i,j,k+1) .and. vof_phase(i,j,k+1)==1 .and. vof_phase(i+1,j,k+1)==1) then
           u_cmask(i,j,k+1,level) = 1; u_assigned(i,j,k+1)=.true.
        endif
        if (.not.u_assigned(i,j,k-1) .and. vof_phase(i,j,k-1)==1 .and. vof_phase(i+1,j,k-1)==1) then
           u_cmask(i,j,k-1,level) = 1; u_assigned(i,j,k-1) = .true.
        endif
     endif
     !v-neighbours
     if (v_cmask(i,j,k,level-1)==1) then
        if (.not.v_assigned(i,j+1,k) .and. vof_phase(i,j+1,k)==1 .and. vof_phase(i,j+2,k)==1) then
           v_cmask(i,j+1,k,level) = 1; v_assigned(i,j+1,k)=.true.
        endif
        if (.not.v_assigned(i+1,j,k) .and. vof_phase(i+1,j,k)==1 .and. vof_phase(i+1,j+1,k)==1) then
           v_cmask(i+1,j,k,level) = 1; v_assigned(i+1,j,k) = .true.
        endif
        if (.not.v_assigned(i,j-1,k) .and. vof_phase(i,j-1,k)==1 .and. vof_phase(i,j,k)==1) then
           v_cmask(i,j-1,k,level) = 1; v_assigned(i,j-1,k)=.true.
        endif
        if (.not.v_assigned(i-1,j,k) .and. vof_phase(i-1,j,k)==1 .and. vof_phase(i-1,j+1,k)==1) then
           v_cmask(i-1,j,k,level) = 1; v_assigned(i-1,j,k) = .true.
        endif
        if (.not.v_assigned(i,j,k+1) .and. vof_phase(i,j,k+1)==1 .and. vof_phase(i,j+1,k+1)==1) then
           v_cmask(i,j,k+1,level) = 1; v_assigned(i,j,k+1)=.true.
        endif
        if (.not.v_assigned(i,j,k-1) .and. vof_phase(i,j,k-1)==1 .and. vof_phase(i,j+1,k-1)==1) then
           v_cmask(i,j,k-1,level) = 1; v_assigned(i,j,k-1) = .true.
        endif
     endif
     !w-neighbours
     if (w_cmask(i,j,k,level-1)==1) then
        if (.not.w_assigned(i,j+1,k) .and. vof_phase(i,j+1,k)==1 .and. vof_phase(i,j+1,k+1)==1) then
           w_cmask(i,j+1,k,level) = 1; w_assigned(i,j+1,k)=.true.
        endif
        if (.not.w_assigned(i+1,j,k) .and. vof_phase(i+1,j,k)==1 .and. vof_phase(i+1,j,k+1)==1) then
           w_cmask(i+1,j,k,level) = 1; w_assigned(i+1,j,k) = .true.
        endif
        if (.not.w_assigned(i,j-1,k) .and. vof_phase(i,j-1,k)==1 .and. vof_phase(i,j-1,k+1)==1) then
           w_cmask(i,j-1,k,level) = 1; w_assigned(i,j-1,k)=.true.
        endif
        if (.not.w_assigned(i-1,j,k) .and. vof_phase(i-1,j,k)==1 .and. vof_phase(i-1,j,k+1)==1) then
           w_cmask(i-1,j,k,level) = 1; w_assigned(i-1,j,k) = .true.
        endif
        if (.not.w_assigned(i,j,k+1) .and. vof_phase(i,j,k+1)==1 .and. vof_phase(i,j,k+2)==1) then
           w_cmask(i,j,k+1,level) = 1; w_assigned(i,j,k+1)=.true.
        endif
        if (.not.w_assigned(i,j,k-1) .and. vof_phase(i,j,k-1)==1 .and. vof_phase(i,j,k)==1) then
           w_cmask(i,j,k-1,level) = 1; w_assigned(i,j,k-1) = .true.
        endif
     endif
  enddo; enddo; enddo
  enddo
  if (debug) then
     do k=ks,ke; do j=js,je; do i=is,ie
        !write outputs for gnuplot 2D
        if (j==(js+je)/2) then
           if (u_cmask(i,j,k,0)==1) write(20,13)xh(i),y(j),z(k)
           !if (v_cmask(i,j,k,0)==1) write(20,13)x(i),yh(j),z(k)
           if (w_cmask(i,j,k,0)==1) write(20,13)x(i),y(j),zh(k)
           if (u_cmask(i,j,k,1)==1) write(21,13)xh(i),y(j),z(k)
           !if (v_cmask(i,j,k,1)==1) write(21,13)x(i),yh(j),z(k)
           if (w_cmask(i,j,k,1)==1) write(21,13)x(i),y(j),zh(k)
           if (u_cmask(i,j,k,2)==1) write(22,13)xh(i),y(j),z(k)
           !if (v_cmask(i,j,k,2)==1) write(22,13)x(i),yh(j),z(k)
           if (w_cmask(i,j,k,2)==1) write(22,13)x(i),y(j),zh(k)
        endif
     enddo; enddo; enddo
  endif
close(unit=20); close(unit=21); close(unit=22)
13 format(3e14.5)
end subroutine set_topology
!-------------------------------------------------------------------------------------------------
subroutine setuppoisson_fs(umask,vmask,wmask,vof_phase,rhot,dt,A,pmask,cvof,n1,n2,n3,kap,istep)    
  use module_grid
  use module_BC
  use module_2phase
  use module_freesurface
  use module_IO
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: rhot
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: umask,vmask,wmask
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: kap
  real(8), dimension(is:ie,js:je,ks:ke), intent(inout) :: pmask
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: cvof,n1,n2,n3
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(out) :: A
  real(8), dimension(is:ie,js:je,ks:ke) :: x_int, y_int, z_int
  integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: vof_phase
  real(8) :: alpha, x_test, y_test, z_test, n_x, n_y, n_z
  real(8) :: x_l, x_r, y_b, y_t, z_r, z_f
  real(8) :: nr(3),al3dnew,al3dold
  real(8) :: dt
  integer :: i,j,k,l,istep
  !new variables
  real(8) :: mod0, mod1, count

  x_int = 0d0; y_int = 0d0; z_int = 0d0
  x_mod=dxh((is+ie)/2); y_mod=dyh((js+je)/2); z_mod=dzh((ks+ke)/2) !assumes an unstretched grid
  P_g = 0d0
  pmask = 1d0
  debug = .false.
  if (debug) then
     Open(unit=50,file='C_int-'//TRIM(int2text(istep,padding))//'.txt') !remove, debugging
     Open(unit=51,file='COEFFS-'//TRIM(int2text(istep,padding))//'.txt') !remove, debugging
     !Open(unit=52,file='P_g-'//TRIM(int2text(istep,padding))//'.txt') !remove, debugging
     !Open(unit=53,file='Phase-'//TRIM(int2text(istep,padding))//'.txt') !remove, debugging 
  endif 

  !do k=ks,ke; do j=js,je; do i=is,ie
  do k=ks,ke-1; do j=js,je-1; do i=is,ie-1 ! For new method, set upper bound to end-1, since we always check next cell in each direction
     !----------------------------------------------------------------------------------------------
     ! New setup: use topology criteria. If Pressure node neighbours are in different phases, A branch is modified
     ! Check for gas-liquid neighbours
     if(vof_phase(i,j,k)==1) then
        if (cvof(i,j,k)<0.5d0) write(*,'("Vof phase error. Phase test 1, cvof: ",e14.5)')cvof(i,j,k) !debugging
        pmask(i,j,k) = 0d0 !pmask local, have to set
        !Check x-neighbour
        if (vof_phase(i+1,j,k) == 0) then
           if (cvof(i+1,j,k)>=0.5d0) write(*,'("Vof phase error. Phase test 0, cvof: ",e14.5)')cvof(i+1,j,k) !debugging
           count = 0d0; mod0 = 0d0; mod1 = 0d0
           !get_intersection for liq cell
           nr(1) = n1(i+1,j,k); nr(2) = n2(i+1,j,k); nr(3) = n3(i+1,j,k)
           n_x = ABS(nr(1)); n_y = ABS(nr(2)); n_z = ABS(nr(3))
           alpha = al3dold(n_x,n_y,n_z,cvof(i+1,j,k))
           if (n_x > 1d-14) then
              x_test = (alpha - (n_y+n_z)/2d0)/n_x
              if (x_test>-0.5d0 .and. n1(i+1,j,k)>0d0) then
                 count = count + 1d0
                 mod0 = (0.5d0-x_test)*dxh(i)
                 if (debug) write(50,312)x(i+1)-mod0,z(k)
              endif
           endif
           !get_intersection for gas cell
           nr(1) = n1(i,j,k); nr(2) = n2(i,j,k); nr(3) = n3(i,j,k)
           n_x = ABS(n1(i,j,k)); n_y = ABS(n2(i,j,k)); n_z = ABS(n3(i,j,k))
           alpha = al3dold(n_x,n_y,n_z,cvof(i,j,k))
           if (n_x > 1d-14) then
              x_test = (alpha - (n_y+n_z)/2d0)/n_x
              if (x_test<1.5d0 .and. n1(i,j,k)>0d0) then
                 count = count + 1d0
                 mod1 = (1.5d0-x_test)*dxh(i)
                 if (debug) write(50,312)x(i+1)-mod1,z(k)
              endif
           endif
           !average if interfaces in both cells intersect in between nodes
           !adapt A coeff
           if (count > 0d0) then
              x_mod(i,j,k) = (mod0 + mod1)/count
           else
              x_mod(i,j,k) = dxh(i)/2d0
              !write(*,'("WARNING: gas-liq x-pair has no interface intercepts between nodes",3I8)')i,j,k
           endif
           A(i+1,j,k,1) = 2d0*dt*umask(i,j,k)/(dx(i+1)*x_mod(i,j,k)*(rhot(i,j,k)+rhot(i+1,j,k)))
           if (debug) write(51,314)x(i+1),z(k),-x_mod(i,j,k),0d0
           if (A(i+1,j,k,1) /= A(i+1,j,k,1)) write(*,'("A1 NaN? ",2e14.4)')A(i+1,j,k,1), x_mod(i,j,k) !debugging
           P_g(i,j,k,1) = sigma*kap(i,j,k)/dx(i) !curvature taken in gas cell.
        endif
        !Check y-neighbour
        if (vof_phase(i,j+1,k) == 0) then
           if (cvof(i,j+1,k)>=0.5d0) write(*,'("Vof phase error. Phase test 0, cvof: ",e14.5)')cvof(i,j+1,k)
           count = 0d0; mod0 = 0d0; mod1 = 0d0
           !get_intersection for liq cell
           nr(1) = n1(i,j+1,k); nr(2) = n2(i,j+1,k); nr(3) = n3(i,j+1,k)
           n_x = ABS(nr(1)); n_y = ABS(nr(2)); n_z = ABS(nr(3))
           alpha = al3dold(n_x,n_y,n_z,cvof(i,j+1,k))
           if (n_y > 1d-14) then
              y_test = (alpha - (n_x+n_z)/2d0)/n_y
              if (y_test>-0.5d0 .and. n2(i,j+1,k)>0d0) then
                 count = count + 1d0
                 mod0 = (0.5d0-y_test)*dyh(j)
                 !if (debug) write(50,312)x(i),y(j+1)-mod0
              endif
           endif
           !get_intersection for gas cell.
           nr(1) = n1(i,j,k); nr(2) = n2(i,j,k); nr(3) = n3(i,j,k)
           n_x = ABS(n1(i,j,k)); n_y = ABS(n2(i,j,k)); n_z = ABS(n3(i,j,k))
           alpha = al3dold(n_x,n_y,n_z,cvof(i,j,k))
           if (n_y > 1d-14) then
              y_test = (alpha - (n_x+n_z)/2d0)/n_y
              if (y_test<1.5d0 .and. n2(i,j,k)>0d0) then
                 count = count + 1d0
                 mod1 = (1.5d0-y_test)*dyh(j)
                 !if (debug) write(50,312)x(i),y(j+1)-mod1
              endif
           endif
           !average if interfaces in both cells intersect in between nodes
           !adapt A coeff
           if (count > 0d0) then
              y_mod(i,j,k) = (mod0 + mod1)/count
           else
              y_mod(i,j,k) = dyh(j)/2d0
              !write(*,'("WARNING: gas-liq y-pair has no interface intercepts between nodes",3I8)')i,j,k
           endif
           A(i,j+1,k,3) = 2d0*dt*vmask(i,j,k)/(dy(j+1)*y_mod(i,j,k)*(rhot(i,j,k)+rhot(i,j+1,k)))
           !if (debug) write(51,314)x(i),y(j+1),0d0,-y_mod(i,j,k)
           if (A(i,j+1,k,3) /= A(i,j+1,k,3)) write(*,'("A3 NaN? ",2e14.4)')A(i,j+1,k,3), y_mod(i,j,k) !debugging
           P_g(i,j,k,2) = sigma*kap(i,j,k)/dy(j) !curvature taken in gas cell.
        endif
        !Check z-neighbour
        if (vof_phase(i,j,k+1) == 0) then
           if (cvof(i,j,k+1)>=0.5d0) write(*,'("Vof phase error. Phase test 0, cvof: ",e14.5)')cvof(i,j,k+1)
           count = 0d0; mod0 = 0d0; mod1 = 0d0
           !get_intersection for liq cell
           nr(1) = n1(i,j,k+1); nr(2) = n2(i,j,k+1); nr(3) = n3(i,j,k+1)
           n_x = ABS(nr(1)); n_y = ABS(nr(2)); n_z = ABS(nr(3))
           alpha = al3dold(n_x,n_y,n_z,cvof(i,j,k+1))
           if (n_z > 1d-14) then
              z_test = (alpha - (n_x+n_y)/2d0)/n_z
              if (z_test>-0.5d0 .and. n3(i,j,k+1)>0d0) then
                 count = count + 1d0
                 mod0 = (0.5d0-z_test)*dzh(k)
                 if (debug) write(50,312)x(i),z(k+1)-mod0
              endif
           endif
           !get_intersection for gas cell
           nr(1) = n1(i,j,k); nr(2) = n2(i,j,k); nr(3) = n3(i,j,k)
           n_x = ABS(n1(i,j,k)); n_y = ABS(n2(i,j,k)); n_z = ABS(n3(i,j,k))
           alpha = al3dold(n_x,n_y,n_z,cvof(i,j,k))
           if (n_z > 1d-14) then
              z_test = (alpha - (n_x+n_y)/2d0)/n_z
              if (z_test<1.5d0 .and. n3(i,j,k)>0d0) then
                 count = count + 1d0
                 mod1 = (1.5d0-z_test)*dzh(k)
                 if (debug) write(50,312)x(i),z(k+1)-mod0
              endif
           endif
           !average if interfaces in both cells intersect in between nodes
           !adapt A coeff
           if (count > 0d0) then
              z_mod(i,j,k) = (mod0 + mod1)/count
           else
              z_mod(i,j,k) = dzh(k)/2d0
              !write(*,'("WARNING: gas-liq z-pair has no interface intercepts between nodes",3I8)')i,j,k
           endif
           A(i,j,k+1,5) = 2d0*dt*wmask(i,j,k)/(dz(k+1)*z_mod(i,j,k)*(rhot(i,j,k)+rhot(i,j,k+1)))
           if (debug) write(51,314)x(i),z(k+1),0d0,-z_mod(i,j,k)
           if (A(i,j,k+1,5) /= A(i,j,k+1,5)) write(*,'("A5 NaN? ",2e14.4)')A(i,j,k+1,5), z_mod(i,j,k) !debugging
           P_g(i,j,k,3) = sigma*kap(i,j,k)/dz(k) !curvature taken in gas cell.
        endif
     endif
!-------------------------Liquid-gas neighbours
     if (vof_phase(i,j,k)==0) then
        if (cvof(i,j,k)>=0.5d0) write(*,'("Vof phase error. Phase test 0, cvof: ",e14.5)')cvof(i,j,k) 
        !Check x-neighbour
        if (vof_phase(i+1,j,k) == 1) then
           if (cvof(i+1,j,k)<0.5d0) write(*,'("Vof phase error. Phase test 1, cvof: ",e14.5)')cvof(i+1,j,k) !debugging
           count = 0d0; mod0 = 0d0; mod1 = 0d0
           pmask(i+1,j,k) = 1d0
           !get_intersection for gas cell
           nr(1) = n1(i+1,j,k); nr(2) = n2(i+1,j,k); nr(3) = n3(i+1,j,k)
           n_x = ABS(nr(1)); n_y = ABS(nr(2)); n_z = ABS(nr(3))
           alpha = al3dold(n_x,n_y,n_z,cvof(i+1,j,k))
           if (n_x > 1d-14) then
              x_test = (alpha - (n_y+n_z)/2d0)/n_x
              if (x_test<1.5d0 .and. n1(i+1,j,k)<0d0) then
                 count = count + 1d0
                 mod1 = (1.5d0-x_test)*dxh(i)
                 if (debug) write(50,312)x(i)+mod1,z(k)
              endif
           endif
           !get_intersection for liq cell
           nr(1) = n1(i,j,k); nr(2) = n2(i,j,k); nr(3) = n3(i,j,k)
           n_x = ABS(n1(i,j,k)); n_y = ABS(n2(i,j,k)); n_z = ABS(n3(i,j,k))
           alpha = al3dold(n_x,n_y,n_z,cvof(i,j,k))
           if (n_x > 1d-14) then
              x_test = (alpha - (n_y+n_z)/2d0)/n_x
              if (x_test>-0.5d0 .and. n1(i,j,k)<0d0) then
                 count = count + 1d0
                 mod0 = (0.5d0-x_test)*dxh(i)
                 if (debug) write(50,312)x(i)+mod0,z(k)
              endif
           endif
           !average if interfaces in both cells intersect in between nodes
           !adapt A coeff
           if (count > 0d0) then
              x_mod(i,j,k) = (mod0 + mod1)/count
           else
              x_mod(i,j,k) = dxh(i)/2d0
              !write(*,'("WARNING: liq-gas x-pair has no interface intercepts between nodes",3I8)')i,j,k
           endif
           A(i,j,k,2) = 2d0*dt*umask(i,j,k)/(dx(i)*x_mod(i,j,k)*(rhot(i,j,k)+rhot(i+1,j,k)))
           if (A(i,j,k,2) /= A(i,j,k,2)) write(*,'("A2 NaN? ",2e14.4)')A(i,j,k,2), x_mod(i,j,k) !debugging
           if (debug) write(51,314)x(i),z(k),x_mod(i,j,k),0d0
           P_g(i+1,j,k,1) = sigma*kap(i,j,k)/dx(i) !check, fix
        endif
        !Check y-neighbour
        if (vof_phase(i,j+1,k) == 1) then
           if (cvof(i,j+1,k)<0.5d0) write(*,'("Vof phase error. Phase test 1, cvof: ",e14.5)')cvof(i,j+1,k) !debugging
           count = 0d0; mod0 = 0d0; mod1 = 0d0
           pmask(i,j+1,k) = 1d0
           !get_intersection for gas cell
           nr(1) = n1(i,j+1,k); nr(2) = n2(i,j+1,k); nr(3) = n3(i,j+1,k)
           n_x = ABS(nr(1)); n_y = ABS(nr(2)); n_z = ABS(nr(3))
           alpha = al3dold(n_x,n_y,n_z,cvof(i,j+1,k))
           if (n_y > 1d-14) then
              y_test = (alpha - (n_x+n_z)/2d0)/n_y
              if (y_test<1.5d0 .and. n2(i,j+1,k)<0d0) then
                 count = count + 1d0
                 mod1 = (1.5d0-y_test)*dyh(j)
                 !if (debug) write(50,312)x(i),y(j)+mod1
              endif
           endif
           !get_intersection for liq cell
           nr(1) = n1(i,j,k); nr(2) = n2(i,j,k); nr(3) = n3(i,j,k)
           n_x = ABS(n1(i,j,k)); n_y = ABS(n2(i,j,k)); n_z = ABS(n3(i,j,k))
           alpha = al3dold(n_x,n_y,n_z,cvof(i,j,k))
           if (n_y > 1d-14) then
              y_test = (alpha - (n_x+n_z)/2d0)/n_y
              if (y_test>-0.5d0 .and. n2(i,j,k)<0d0) then
                 count = count + 1d0
                 mod0 = (0.5d0-y_test)*dyh(j)
                 !if (debug) write(50,312)x(i),y(j)+mod0
              endif
           endif
           !average if interfaces in both cells intersect in between nodes
           !adapt A coeff
           if (count > 0d0) then
              y_mod(i,j,k) = (mod0 + mod1)/count
           else
              y_mod(i,j,k) = dyh(j)/2d0
              !write(*,'("WARNING: liq-gas y-pair has no interface intercepts between nodes",3I8)')i,j,k
           endif
           A(i,j,k,4) = 2d0*dt*vmask(i,j,k)/(dy(j)*y_mod(i,j,k)*(rhot(i,j,k)+rhot(i,j+1,k)))
           !if (debug) write(51,314)x(i),y(j),0d0,y_mod(i,j,k)
           if (A(i,j,k,4) /= A(i,j,k,4)) write(*,'("A4 NaN? ",2e14.4)')A(i,j,k,4), y_mod(i,j,k) !debugging
           P_g(i,j+1,k,2) = sigma*kap(i,j,k)/dy(j) !check, fix
        endif
        !Check z-neighbour
        if (vof_phase(i,j,k+1) == 1) then
           if (cvof(i,j,k+1)<0.5d0) write(*,'("Vof phase error. Phase test 1, cvof: ",e14.5)')cvof(i,j,k+1) !debugging
           count = 0d0; mod0 = 0d0; mod1 = 0d0
           pmask(i,j,k+1) = 1d0
           !get_intersection for gas cell
           nr(1) = n1(i,j,k+1); nr(2) = n2(i,j,k+1); nr(3) = n3(i,j,k+1)
           n_x = ABS(nr(1)); n_y = ABS(nr(2)); n_z = ABS(nr(3))
           alpha = al3dold(n_x,n_y,n_z,cvof(i,j,k+1))
           if (n_z > 1d-14) then
              z_test = (alpha - (n_x+n_y)/2d0)/n_z
              if (z_test<1.5d0 .and. n3(i,j,k+1)<0d0) then
                 count = count + 1d0
                 mod1 = (1.5d0-z_test)*dzh(k)
                 if (debug) write(50,312)x(i),z(k)+mod1
              endif
           endif
           !get_intersection for liq cell
           nr(1) = n1(i,j,k); nr(2) = n2(i,j,k); nr(3) = n3(i,j,k)
           n_x = ABS(n1(i,j,k)); n_y = ABS(n2(i,j,k)); n_z = ABS(n3(i,j,k))
           alpha = al3dold(n_x,n_y,n_z,cvof(i,j,k))
           
           if (n_z > 1d-14) then
              z_test = (alpha - (n_x+n_y)/2d0)/n_z
              if (z_test>-0.5d0 .and. n3(i,j,k)<0d0) then
                 count = count + 1d0
                 mod0 = (0.5d0-z_test)*dzh(k)
                 if (debug) write(50,312)x(i),z(k)+mod0
              endif
           endif
           !average if interfaces in both cells intersect in between nodes
           !adapt A coeff
           if (count > 0d0) then
              z_mod(i,j,k) = (mod0 + mod1)/count
           else
              z_mod(i,j,k) = dzh(k)/2d0
              !write(*,'("WARNING: liq-gas z-pair has no interface intercepts between nodes",3I8)')i,j,k
           endif
           A(i,j,k,6) = 2d0*dt*wmask(i,j,k)/(dz(k)*z_mod(i,j,k)*(rhot(i,j,k)+rhot(i,j,k+1)))
           if (debug) write(51,314)x(i),z(k),0d0,z_mod(i,j,k)
           if (A(i,j,k,6) /= A(i,j,k,6)) write(*,'("A6 NaN? ",2e14.4)')A(i,j,k,6), z_mod(i,j,k) !debugging
           P_g(i,j,k+1,3) = sigma*kap(i,j,k)/dz(k) !check, fix
        endif
     endif
  enddo;enddo;enddo
close(unit=50); close(unit=51); close(unit=52); close(unit=53)
312 format(2e14.5) !remove, debugging
313 format(3e14.5) !remove, debugging
323 format(2e14.5,e14.4) !remove, debugging
314 format(4e14.5) !remove, debugging
414 format(3e14.5,I8)
  do k=ks,ke; do j=js,je; do i=is,ie
     do l=1,6
        A(i,j,k,l) = pmask(i,j,k)*A(i,j,k,l)
        if (A(i,j,k,l) /= A(i,j,k,l)) write(*,*)'A NaN, error imminent'
     enddo
     A(i,j,k,7) = sum(A(i,j,k,1:6)) + (1d0-pmask(i,j,k))
     A(i,j,k,8) = pmask(i,j,k)*(A(i,j,k,8) + dt/(rhot(i,j,k))*&
          (P_g(i+1,j,k,1)/(dx(i)*x_mod(i,j,k))+P_g(i-1,j,k,1)/(dx(i)*x_mod(i-1,j,k))&
          +P_g(i,j+1,k,2)/(dy(j)*y_mod(i,j,k))+P_g(i,j-1,k,2)/(dy(j)*y_mod(i,j-1,k))&
          +P_g(i,j,k+1,3)/(dz(k)*z_mod(i,j,k))+P_g(i,j,k-1,3)/(dz(k)*z_mod(i,j,k-1))))
!!$     if (A(i,j,k,8) /= A(i,j,k,8)) then
!!$        write(*,'("A8 NaN, error imminent. Neigbours mods :",6e14.5)')x_mod(i-1,j,k),x_mod(i,j,k),&
!!$             y_mod(i,j-1,k),y_mod(i,j,k),z_mod(i,j,k-1),z_mod(i,j,k)
!!$        write(*,'("A8 NaN, error imminent. P_g :",6e14.5)')P_g(i-1,j,k,1),P_g(i+1,j,k,1),&
!!$             P_g(i,j-1,k,2),P_g(i,j+1,k,2),P_g(i,j,k-1,3),P_g(i,j,k+1,3)
!!$        write(*,'("pmask :",e14.5)')pmask(i,j,k)
!!$        write(*,'("sigma :",e14.5)')sigma
!!$        write(*,'("A8 NaN, kap :",6e14.5)')kap(i-1,j,k),kap(i,j,k),&
!!$             kap(i,j-1,k),kap(i,j,k),kap(i,j,k-1),kap(i,j,k)
!!$     endif
  enddo;enddo;enddo
end subroutine setuppoisson_fs
