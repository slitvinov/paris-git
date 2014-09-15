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
  use module_BC
  implicit none
  include 'mpif.h' 
  integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: vof_phase
  integer :: req(16),sta(MPI_STATUS_SIZE,16)
  integer :: i,j,k,level,iout,nbr,ierr
  !initialize all masks to 3
  if (.not.initialize_fs) call pariserror("Error: Free surface variables not initialized")
  u_cmask = 3; v_cmask = 3; w_cmask = 3; pcmask = 3
  debug = .false.
  if (debug) then
     Open(unit=19,FILE=TRIM(out_path)//'/Pmask-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout,padding))//'.txt')
     Open(unit=20,FILE=TRIM(out_path)//'/Top_0-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout,padding))//'.txt')
     Open(unit=21,FILE=TRIM(out_path)//'/Top_1-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout,padding))//'.txt')
     Open(unit=22,FILE=TRIM(out_path)//'/Top_2-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout,padding))//'.txt')
     Open(unit=23,FILE=TRIM(out_path)//'/u_ass-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout,padding))//'.txt')
     Open(unit=24,FILE=TRIM(out_path)//'/v_ass-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout,padding))//'.txt')
     Open(unit=25,FILE=TRIM(out_path)//'/w_ass-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout,padding))//'.txt')
     Open(unit=26,FILE=TRIM(out_path)//'/P1-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout,padding))//'.txt')
     Open(unit=27,FILE=TRIM(out_path)//'/P2-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout,padding))//'.txt')
     Open(unit=28,FILE=TRIM(out_path)//'/P3-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout,padding))//'.txt') 
  endif
  !First loop to set level 0 velocities in liq-liq and liq-gas cells
  do k=ks,ke; do j=js,je; do i=is,ie !loop relies on vof_phase in e+1, but phase is updated up to max.
     if (vof_phase(i,j,k) == 1) then 
        if (debug) write(19,13)x(i),y(j),z(k)
        if (vof_phase(i+1,j,k) == 0) then
           u_cmask(i,j,k) = 0
        endif
        if (vof_phase(i,j+1,k) == 0) then
           v_cmask(i,j,k) = 0
        endif
        if (vof_phase(i,j,k+1) == 0) then 
           w_cmask(i,j,k) = 0
        endif
     endif
     if (vof_phase(i,j,k) == 0) then
        pcmask(i,j,k)=0
        if ((vof_phase(i+1,j,k) == 1) .or. (vof_phase(i+1,j,k) == 0)) then
           u_cmask(i,j,k) = 0
        endif
        if ((vof_phase(i,j+1,k) == 1) .or. (vof_phase(i,j+1,k) == 0)) then
           v_cmask(i,j,k) = 0
        endif
        if ((vof_phase(i,j,k+1) == 1) .or.(vof_phase(i,j,k+1) == 0)) then 
           w_cmask(i,j,k) = 0
        endif
     endif
  enddo; enddo; enddo
  !need to set bc for assigned, masks
  !fill ghost layers
  call ighost_x(u_cmask,2,req(1:4)); call ighost_x(v_cmask,2,req(5:8)); call ighost_x(w_cmask,2,req(9:12))
  call ighost_x(pcmask,2,req(13:16)); call MPI_WAITALL(16,req(1:16),sta(:,1:16),ierr)
  call ighost_y(u_cmask,2,req(1:4)); call ighost_y(v_cmask,2,req(5:8)); call ighost_y(w_cmask,2,req(9:12))
  call ighost_y(pcmask,2,req(13:16)); call MPI_WAITALL(16,req(1:16),sta(:,1:16),ierr)
  call ighost_z(u_cmask,2,req(1:4)); call ighost_z(v_cmask,2,req(5:8)); call ighost_z(w_cmask,2,req(9:12))
  call ighost_z(pcmask,2,req(13:16)); call MPI_WAITALL(16,req(1:16),sta(:,1:16),ierr)
  !Set levels 1 to X_level
  do level=1,X_level
     do k=kmin+1,kmax-1; do j=jmin+1,jmax-1; do i=imin+1,imax-1 !uses indexes s-1 and e+1 to check assigned and set masks. 2 ghost layers
        !u-neighbours
        if (u_cmask(i,j,k)==level-1) then
           do nbr=-1,1,2
              if (u_cmask(i+nbr,j,k)==3) then ! .and. vof_phase(i+nbr,j,k)==1 .and. vof_phase(i+nbr+1,j,k)==1) then
                 u_cmask(i+nbr,j,k) = level
              endif
              if (u_cmask(i,j+nbr,k)==3) then !.and. vof_phase(i,j+nbr,k)==1 .and. vof_phase(i+1,j+nbr,k)==1) then
                 u_cmask(i,j+nbr,k) = level
              endif
              if (u_cmask(i,j,k+nbr)==3) then ! .and. vof_phase(i,j,k+nbr)==1 .and. vof_phase(i+1,j,k+nbr)==1) then
                 u_cmask(i,j,k+nbr) = level
              endif
           enddo
        endif
        !v-neighbours
        if (v_cmask(i,j,k)==level-1) then
           do nbr=-1,1,2
              if (v_cmask(i+nbr,j,k)==3) then ! .and. vof_phase(i+nbr,j,k)==1 .and. vof_phase(i+nbr,j+1,k)==1) then
                 v_cmask(i+nbr,j,k) = level
              endif
              if (v_cmask(i,j+nbr,k)==3) then ! .and. vof_phase(i,j+nbr,k)==1 .and. vof_phase(i,j+nbr+1,k)==1) then
                 v_cmask(i,j+nbr,k) = level
              endif
              if (v_cmask(i,j,k+nbr)==3) then ! .and. vof_phase(i,j,k+nbr)==1 .and. vof_phase(i,j+1,k+nbr)==1) then
                 v_cmask(i,j,k+nbr) = level
              endif
           enddo
        endif
        !w-neighbours
        if (w_cmask(i,j,k)==level-1) then
           do nbr=-1,1,2
              if (w_cmask(i+nbr,j,k)==3) then ! .and. vof_phase(i+nbr,j,k)==1 .and. vof_phase(i+nbr,j,k+1)==1) then
                 w_cmask(i+nbr,j,k) = level
              endif
              if (w_cmask(i,j+nbr,k)==3) then ! .and. vof_phase(i,j+nbr,k)==1 .and. vof_phase(i,j+nbr,k+1)==1) then
                 w_cmask(i,j+nbr,k) = level
              endif
              if (w_cmask(i,j,k+nbr)==3) then ! .and. vof_phase(i,j,k+nbr)==1 .and. vof_phase(i,j,k+nbr+1)==1) then
                 w_cmask(i,j,k+nbr) = level
              endif
           enddo
        endif
        if (pcmask(i,j,k)==level-1) then
           do nbr=-1,1,2
              if (pcmask(i+nbr,j,k)==3) then
                 pcmask(i+nbr,j,k) = level
              endif
              if (pcmask(i,j+nbr,k)==3) then
                 pcmask(i,j+nbr,k) = level
              endif
              if (pcmask(i,j,k+nbr)==3) then
                 pcmask(i,j,k+nbr) = level
              endif
           enddo
        endif
     enddo; enddo; enddo
!!$     call ighost_x(u_cmask(:,:,:),2,req(1:4)); call ighost_x(v_cmask(:,:,:),2,req(5:8))
!!$     call ighost_x(w_cmask(:,:,:),2,req(9:12))
!!$     call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
!!$     call ighost_y(u_cmask(:,:,:),2,req(1:4)); call ighost_y(v_cmask(:,:,:),2,req(5:8))
!!$     call ighost_y(w_cmask(:,:,:),2,req(9:12))
!!$     call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
!!$     call ighost_z(u_cmask(:,:,:),2,req(1:4)); call ighost_z(v_cmask(:,:,:),2,req(5:8))
!!$     call ighost_z(w_cmask(:,:,:),2,req(9:12))
!!$     call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
  enddo

  if (debug) then
     do j=jmin,jmax; do i=imin,imax
        !write outputs for gnuplot 2D
        k=(ks+ke)/2
        if (u_cmask(i,j,k)==0) write(20,13)xh(i),y(j),z(k)
        if (v_cmask(i,j,k)==0) write(20,13)x(i),yh(j),z(k)
        !if (w_cmask(i,j,k)==0) write(20,13)x(i),y(j),zh(k)
        if (u_cmask(i,j,k)==1) write(21,13)xh(i),y(j),z(k)
        if (v_cmask(i,j,k)==1) write(21,13)x(i),yh(j),z(k)
        !if (w_cmask(i,j,k)==1) write(21,13)x(i),y(j),zh(k)
        if (pcmask(i,j,k)==1) write(26,13)x(i),y(j),z(k)
        if (u_cmask(i,j,k)==2) write(22,13)xh(i),y(j),z(k)
        if (v_cmask(i,j,k)==2) write(22,13)x(i),yh(j),z(k)
        !if (w_cmask(i,j,k)==2) write(22,13)x(i),y(j),zh(k)
        if (pcmask(i,j,k)==2) write(27,13)x(i),y(j),z(k)
        if (u_cmask(i,j,k)==3) write(23,13)xh(i),y(j),z(k)
        if (v_cmask(i,j,k)==3) write(24,13)x(i),yh(j),z(k)
        if (w_cmask(i,j,k)==3) write(25,13)x(i),y(j),zh(k)
        if (pcmask(i,j,k)==3) write(28,13)x(i),y(j),z(k)
     enddo; enddo
     close(unit=19); close(unit=20); close(unit=21); close(unit=22); close(unit=23); close(unit=24); close(unit=25)
  endif
13 format(3e14.5)
end subroutine set_topology
!-------------------------------------------------------------------------------------------------
subroutine setuppoisson_fs(vof_phase,rhot,dt,A,pmask,cvof,n1,n2,n3,kap,istep)    
  use module_grid
  use module_BC
  use module_2phase
  use module_freesurface
  use module_IO
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: rhot
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: kap
  real(8), dimension(is:ie,js:je,ks:ke), intent(inout) :: pmask
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: cvof,n1,n2,n3
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(out) :: A
  integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: vof_phase
  integer :: req(24),sta(MPI_STATUS_SIZE,24)
  real(8) :: alpha, x_test, y_test, z_test, n_x, n_y, n_z
  real(8) :: nr(3),al3dnew,al3dold
  real(8) :: dt
  integer :: i,j,k,l,istep,nbr,ierr
  real(8) :: mod0, mod1, count

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
  
  do k=ks,ke; do j=js,je; do i=is,ie
     ! Check for gas-liquid neighbours, set pmask and P_g
     if(vof_phase(i,j,k)==1) then
        if (cvof(i,j,k)<0.5d0) write(*,'("Vof phase error. Phase test 1, cvof: ",e14.5)')cvof(i,j,k) !debugging
        pmask(i,j,k)=0d0
        do nbr=-1,1,2
           if (vof_phase(i+nbr,j,k)==0) P_g(i,j,k,1) = sigma*kap(i,j,k)/dx(i) !!filaments and droplets of one cell will be a problem
           if (vof_phase(i,j+nbr,k)==0) P_g(i,j,k,2) = sigma*kap(i,j,k)/dy(j)
           if (vof_phase(i,j,k+nbr)==0) P_g(i,j,k,3) = sigma*kap(i,j,k)/dz(k)
        enddo
     endif
     ! Liquid-gas neighbours
     if (vof_phase(i,j,k)==0) then
        if (cvof(i,j,k)>=0.5d0) write(*,'("Vof phase error. Phase test 0, cvof: ",e14.5)')cvof(i,j,k) !debugging check
        !Check x-neighbours in 2 directions 
        do nbr=-1,1,2
           if (vof_phase(i+nbr,j,k) == 1) then
              if (cvof(i+nbr,j,k)<0.5d0) write(*,'("Vof phase error. Phase test 1, cvof: ",e14.5)')cvof(i+nbr,j,k) !debugging
              count=0d0; mod0=0d0; mod1=0d0
              !get_intersection for gas cell
              nr(1) = n1(i+nbr,j,k); nr(2) = n2(i+nbr,j,k); nr(3) = n3(i+nbr,j,k)
              n_x = ABS(nr(1)); n_y = ABS(nr(2)); n_z = ABS(nr(3))
              alpha = al3dold(n_x,n_y,n_z,cvof(i+nbr,j,k))
              if (n_x > 1d-14) then
                 x_test = (alpha - (n_y+n_z)/2d0)/n_x
                 if ( (x_test<1.5d0) .and. (n1(i+nbr,j,k)*nbr)<0d0 ) then !warning, mult real by integer
                    count = count + 1d0
                    mod1 = (1.5d0-x_test)*dxh(i+(nbr-1)/2)
                    if (debug) write(50,312)x(i)+nbr*mod1,y(j)!z(k) !warning, multiplication of real with integer
                 endif
              endif
              !get_intersection for liq cell
              nr(1) = n1(i,j,k); nr(2) = n2(i,j,k); nr(3) = n3(i,j,k)
              n_x = ABS(n1(i,j,k)); n_y = ABS(n2(i,j,k)); n_z = ABS(n3(i,j,k))
              alpha = al3dold(n_x,n_y,n_z,cvof(i,j,k))
              if (n_x > 1d-14) then
                 x_test = (alpha - (n_y+n_z)/2d0)/n_x
                 if (x_test>-0.5d0 .and. n1(i,j,k)*nbr<0d0) then
                    count = count + 1d0
                    mod0 = (0.5d0-x_test)*dxh(i+(nbr-1)/2)
                    if (debug) write(50,312)x(i)+nbr*mod0,y(j)
                 endif
              endif
              !enddo
              !average if interfaces in both cells intersect in between nodes
              !adapt A coeff
              if (count > 0d0) then
                 x_mod(i+(nbr-1)/2,j,k) = (mod0 + mod1)/count
              else
                 x_mod(i+(nbr-1)/2,j,k) = dxh(i+(nbr-1)/2)/2d0
                 !write(*,'("WARNING: liq-gas x-pair has no interface intercepts between nodes",3I8)')i,j,k
              endif
              A(i,j,k,2+(nbr-1)/2) = dt/(dx(i)*x_mod(i+(nbr-1)/2,j,k)*(rhot(i,j,k)))
              if (A(i,j,k,2+(nbr-1)/2) /= A(i,j,k,2+(nbr-1)/2)) then
                 write(*,'("A1 or A2 NaN? ",2e14.4)')A(i,j,k,2+(nbr-1)/2),x_mod(i+(nbr-1)/2,j,k) !debugging
              endif
              if (debug) write(51,314)x(i),y(j),nbr*x_mod(i+(nbr-1)/2,j,k),0d0
              !P_g(i+nbr,j,k,1) = sigma*kap(i,j,k)/dx(i) !check, fix
           endif
        enddo
        !Check y-neighbours in both directions
        do nbr=-1,1,2
           if (vof_phase(i,j+nbr,k) == 1) then
              if (cvof(i,j+nbr,k)<0.5d0) write(*,'("Vof phase error. Phase test 1, cvof: ",e14.5)')cvof(i,j+nbr,k) !debugging
              count = 0d0; mod0 = 0d0; mod1 = 0d0
              !get_intersection for gas cell
              nr(1) = n1(i,j+nbr,k); nr(2) = n2(i,j+nbr,k); nr(3) = n3(i,j+nbr,k)
              n_x = ABS(nr(1)); n_y = ABS(nr(2)); n_z = ABS(nr(3))
              alpha = al3dold(n_x,n_y,n_z,cvof(i,j+nbr,k))
              if (n_y > 1d-14) then
                 y_test = (alpha - (n_x+n_z)/2d0)/n_y
                 if ((y_test<1.5d0) .and. (n2(i,j+nbr,k)*nbr)<0d0) then
                    count = count + 1d0
                    mod1 = (1.5d0-y_test)*dyh(j+(nbr-1)/2)
                    if (debug) write(50,312)x(i),y(j)+nbr*mod1
                 endif
              endif
              !get_intersection for liq cell
              nr(1) = n1(i,j,k); nr(2) = n2(i,j,k); nr(3) = n3(i,j,k)
              n_x = ABS(n1(i,j,k)); n_y = ABS(n2(i,j,k)); n_z = ABS(n3(i,j,k))
              alpha = al3dold(n_x,n_y,n_z,cvof(i,j,k))
              if (n_y > 1d-14) then
                 y_test = (alpha - (n_x+n_z)/2d0)/n_y
                 if (y_test>-0.5d0 .and. n2(i,j,k)*nbr<0d0) then
                    count = count + 1d0
                    mod0 = (0.5d0-y_test)*dyh(j+(nbr-1)/2)
                    if (debug) write(50,312)x(i),y(j)+nbr*mod0
                 endif
              endif
              !average if interfaces in both cells intersect in between nodes
              !adapt A coeff
              if (count > 0d0) then
                 y_mod(i,j+(nbr-1)/2,k) = (mod0 + mod1)/count
              else
                 y_mod(i,j+(nbr-1)/2,k) = dyh(j+(nbr-1)/2)/2d0
                 !write(*,'("WARNING: liq-gas y-pair has no interface intercepts between nodes",3I8)')i,j,k
              endif
              A(i,j,k,4+(nbr-1)/2) = dt/(dy(j)*y_mod(i,j+(nbr-1)/2,k)*(rhot(i,j,k)))
              if (debug) write(51,314)x(i),y(j),0d0,nbr*y_mod(i,j+(nbr-1)/2,k)
              if (A(i,j,k,4+(nbr-1)/2) /= A(i,j,k,4+(nbr-1)/2)) then
                 write(*,'("A3 or A4 NaN? ",2e14.4)')A(i,j,k,4+(nbr-1)/2), y_mod(i,j+(nbr-1)/2,k) !debugging
              endif
              !P_g(i,j+1,k,2) = sigma*kap(i,j,k)/dy(j) !check, fix
           endif
        enddo
        !Check z-neighbours
        do nbr=-1,1,2
           if (vof_phase(i,j,k+nbr) == 1) then
              if (cvof(i,j,k+nbr)<0.5d0) write(*,'("Vof phase error. Phase test 1, cvof: ",e14.5)')cvof(i,j,k+nbr) !debugging
              count = 0d0; mod0 = 0d0; mod1 = 0d0
              !get_intersection for gas cell
              nr(1) = n1(i,j,k+nbr); nr(2) = n2(i,j,k+nbr); nr(3) = n3(i,j,k+nbr)
              n_x = ABS(nr(1)); n_y = ABS(nr(2)); n_z = ABS(nr(3))
              alpha = al3dold(n_x,n_y,n_z,cvof(i,j,k+nbr))
              if (n_z > 1d-14) then
                 z_test = (alpha - (n_x+n_y)/2d0)/n_z
                 if (z_test<1.5d0 .and. (n3(i,j,k+nbr)*nbr)<0d0) then
                    count = count + 1d0
                    mod1 = (1.5d0-z_test)*dzh(k+(nbr-1)/2)
                    !if (debug) write(50,312)x(i),z(k)+mod1
                 endif
              endif
              !get_intersection for liq cell
              nr(1) = n1(i,j,k); nr(2) = n2(i,j,k); nr(3) = n3(i,j,k)
              n_x = ABS(n1(i,j,k)); n_y = ABS(n2(i,j,k)); n_z = ABS(n3(i,j,k))
              alpha = al3dold(n_x,n_y,n_z,cvof(i,j,k))
              if (n_z > 1d-14) then
                 z_test = (alpha - (n_x+n_y)/2d0)/n_z
                 if (z_test>-0.5d0 .and. n3(i,j,k)*nbr<0d0) then
                    count = count + 1d0
                    mod0 = (0.5d0-z_test)*dzh(k)
                    !if (debug) write(50,312)x(i),z(k)+mod0
                 endif
              endif
              !average if interfaces in both cells intersect in between nodes
              !adapt A coeff
              if (count > 0d0) then
                 z_mod(i,j,k+(nbr-1)/2) = (mod0 + mod1)/count
              else
                 z_mod(i,j,k+(nbr-1)/2) = dzh(k+(nbr-1)/2)/2d0
                 !write(*,'("WARNING: liq-gas z-pair has no interface intercepts between nodes",3I8)')i,j,k
              endif
              A(i,j,k,6+(nbr-1)/2) = dt/(dz(k)*z_mod(i,j,k+(nbr-1)/2)*(rhot(i,j,k)))
              !if (debug) write(51,314)x(i),z(k),0d0,z_mod(i,j,k)
              if (A(i,j,k,6+(nbr-1)/2) /= A(i,j,k,6+(nbr-1)/2)) then
                 write(*,'("A5 or A6 NaN? ",2e14.4)')A(i,j,k,6), z_mod(i,j,k+(nbr-1)/2) !debugging
              endif
              !P_g(i,j,k+1,3) = sigma*kap(i,j,k)/dz(k) !check, fix
           endif
        enddo
     endif
  enddo;enddo;enddo
  call ghost_x(P_g(:,:,:,1),1,req( 1: 4)); call ghost_y(P_g(:,:,:,2),1,req( 5: 8)); call ghost_z(P_g(:,:,:,3),1,req( 9:12))
  call ghost_x(x_mod,1,req(13:16)); call ghost_y(y_mod,1,req(17:20)); call ghost_z(z_mod,1,req(21:24))
  call MPI_WAITALL(24,req,sta,ierr) 
close(unit=50); close(unit=51); close(unit=52); close(unit=53)
312 format(2e14.5) !remove, debugging
313 format(3e14.5) !remove, debugging
323 format(2e14.5,e14.4) !remove, debugging
314 format(4e14.5) !remove, debugging
414 format(3e14.5,I8)
  do k=ks,ke; do j=js,je; do i=is,ie
     !do l=1,6
     !   A(i,j,k,l) = pmask(i,j,k)*A(i,j,k,l)
     !   if (A(i,j,k,l) /= A(i,j,k,l)) write(*,'("A*pmask is NaN. pmask, A :",2e14.5)')pmask(i,j,k),A(i,j,k,l)
     !enddo
     A(i,j,k,7) = sum(A(i,j,k,1:6)) !+ (1d0-pmask(i,j,k))
     A(i,j,k,8) = A(i,j,k,8) + dt/rhot(i,j,k)*&
          (P_g(i+1,j,k,1)/(dx(i)*x_mod(i,j,k))+P_g(i-1,j,k,1)/(dx(i)*x_mod(i-1,j,k))&
          +P_g(i,j+1,k,2)/(dy(j)*y_mod(i,j,k))+P_g(i,j-1,k,2)/(dy(j)*y_mod(i,j-1,k))&
          +P_g(i,j,k+1,3)/(dz(k)*z_mod(i,j,k))+P_g(i,j,k-1,3)/(dz(k)*z_mod(i,j,k-1)))
     if (A(i,j,k,8) /= A(i,j,k,8)) then
        write(*,'("A8 NaN, error imminent. Neigbours mods :",6e14.5)')x_mod(i-1,j,k),x_mod(i,j,k),&
             y_mod(i,j-1,k),y_mod(i,j,k),z_mod(i,j,k-1),z_mod(i,j,k)
        write(*,'("A8 NaN, error imminent. P_g :",6e14.5)')P_g(i-1,j,k,1),P_g(i+1,j,k,1),&
             P_g(i,j-1,k,2),P_g(i,j+1,k,2),P_g(i,j,k-1,3),P_g(i,j,k+1,3)
        write(*,'("rho :",e14.5)')rhot(i,j,k)
        write(*,'("pmask :",e14.5)')pmask(i,j,k)
!!$        write(*,'("A8 NaN, kap :",6e14.5)')kap(i-1,j,k),kap(i,j,k),&
!!$             kap(i,j-1,k),kap(i,j,k),kap(i,j,k-1),kap(i,j,k)
     endif
  enddo;enddo;enddo
end subroutine setuppoisson_fs
!--------------------------------------------------------------------------------------------------------------------
!!$subroutine setuppoisson_fs2(dt,A,cvof,istep)
!!$  use module_grid
!!$  use module_BC
!!$  use module_2phase
!!$  use module_IO
!!$  implicit none
!!$  !real(8), dimension(is:ie,js:je,ks:ke), intent(inout) :: pmask2
!!$  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: cvof
!!$  real(8), dimension(is:ie,js:je,ks:ke,8), intent(out) :: A
!!$  real(8) :: dt
!!$  integer :: i,j,k
!!$  integer :: istep
!!$
!!$  do k=ks,ke; do j=js,je; do i=is,ie;
!!$     if(cvof(i,j,k) >= 0.5d0) then 
!!$        !pmask2(i,j,k) = 1d0 
!!$        A(i,j,k,1) = 2d0*dt*umask(i-1,j,k)/(dx(i)*dxh(i-1)*(rhot(i-1,j,k)+rhot(i,j,k)))
!!$        A(i,j,k,2) = 2d0*dt*umask(i,j,k)/(dx(i)*dxh(i  )*(rhot(i+1,j,k)+rhot(i,j,k)))
!!$        A(i,j,k,3) = 2d0*dt*vmask(i,j-1,k)/(dy(j)*dyh(j-1)*(rhot(i,j-1,k)+rhot(i,j,k)))
!!$        A(i,j,k,4) = 2d0*dt*vmask(i,j,k)/(dy(j)*dyh(j  )*(rhot(i,j+1,k)+rhot(i,j,k)))
!!$        A(i,j,k,5) = 2d0*dt*wmask(i,j,k-1)/(dz(k)*dzh(k-1)*(rhot(i,j,k-1)+rhot(i,j,k)))
!!$        A(i,j,k,6) = 2d0*dt*wmask(i,j,k)/(dz(k)*dzh(k  )*(rhot(i,j,k+1)+rhot(i,j,k)))
!!$        A(i,j,k,7) = sum(A(i,j,k,1:6))
!!$        A(i,j,k,8) =  -((utmp(i,j,k)-utmp(i-1,j,k))/dx(i) &
!!$             +  (vtmp(i,j,k)-vtmp(i,j-1,k))/dy(j) &
!!$             +  (wtmp(i,j,k)-wtmp(i,j,k-1))/dz(k) )
!!$     else
!!$        !pmask2(i,j,k)=0d0
!!$        do l=1,6
!!$           A(i,j,k,l) = 0d0
!!$        enddo
!!$        A(i,j,k,7) = 1d0
!!$        A(i,j,k,8) = 0d0
!!$     endif
!!$  enddo; enddo; enddo
!!$  !do k=ks,ke; do j=js,je; do i=is,ie
!!$  !   do l=1,6
!!$  !      A(i,j,k,l) = pmask(i,j,k)*A(i,j,k,l)
!!$  !   enddo
!!$  !   A(i,j,k,7) = sum(A(i,j,k,1:6)) + (1d0-pmask(i,j,k))
!!$  !   A(i,j,k,8) = pmask(i,j,k)*A(i,j,k,8)
!!$  !enddo;enddo;enddo
!!$end subroutine setuppoisson_fs2
!--------------------------------------------------------------------------------------------------------------------
