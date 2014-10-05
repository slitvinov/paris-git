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
  integer :: i,j,k,level,iout,nbr,ierr,no

  !initialize all masks to 3
  if (.not.initialize_fs) call pariserror("Error: Free surface variables not initialized")
  u_cmask = 3; v_cmask = 3; w_cmask = 3; pcmask = 3
  !First loop to set level 0 velocities in liq-liq and liq-gas cells
  do k=ks,ke; do j=js,je; do i=is,ie
     if (vof_phase(i,j,k) == 1) then
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
        u_cmask(i,j,k) = 0
        v_cmask(i,j,k) = 0
        w_cmask(i,j,k) = 0
     endif
  enddo; enddo; enddo
  !fill ghost layers for zero masks
  call ighost_x(u_cmask,2,req(1:4)); call ighost_x(v_cmask,2,req(5:8)); call ighost_x(w_cmask,2,req(9:12))
  call ighost_x(pcmask,2,req(13:16)); call MPI_WAITALL(16,req(1:16),sta(:,1:16),ierr)
  call ighost_y(u_cmask,2,req(1:4)); call ighost_y(v_cmask,2,req(5:8)); call ighost_y(w_cmask,2,req(9:12))
  call ighost_y(pcmask,2,req(13:16)); call MPI_WAITALL(16,req(1:16),sta(:,1:16),ierr)
  call ighost_z(u_cmask,2,req(1:4)); call ighost_z(v_cmask,2,req(5:8)); call ighost_z(w_cmask,2,req(9:12))
  call ighost_z(pcmask,2,req(13:16)); call MPI_WAITALL(16,req(1:16),sta(:,1:16),ierr)
  !Set levels 1 to X_level
  do level=1,X_level
     do k=ks,ke; do j=js,je; do i=is,ie
        !u-neighbours
        if (u_cmask(i,j,k)==3) then
           if ((u_cmask(i+1,j,k)==level-1).or.(u_cmask(i-1,j,k)==level-1)&
                .or.(u_cmask(i,j+1,k)==level-1).or.(u_cmask(i,j-1,k)==level-1)&
                .or.(u_cmask(i,j,k+1)==level-1).or.(u_cmask(i,j,k-1)==level-1)) then
           u_cmask(i,j,k)=level
           endif
        endif
        !v-neighbours
        if (v_cmask(i,j,k)==3) then
           if ((v_cmask(i+1,j,k)==level-1).or.(v_cmask(i-1,j,k)==level-1)&
                .or.(v_cmask(i,j+1,k)==level-1).or.(v_cmask(i,j-1,k)==level-1)&
                .or.(v_cmask(i,j,k+1)==level-1).or.(v_cmask(i,j,k-1)==level-1)) then
           v_cmask(i,j,k)=level
           endif
        endif
        !w-neighbours
        if (w_cmask(i,j,k)==3) then
           if ((w_cmask(i+1,j,k)==level-1).or.(w_cmask(i-1,j,k)==level-1)&
                .or.(w_cmask(i,j+1,k)==level-1).or.(w_cmask(i,j-1,k)==level-1)&
                .or.(w_cmask(i,j,k+1)==level-1).or.(w_cmask(i,j,k-1)==level-1)) then
           w_cmask(i,j,k)=level
           endif
        endif
        !p-neighbours
        if (pcmask(i,j,k)==3) then
           if ((pcmask(i+1,j,k)==level-1).or.(pcmask(i-1,j,k)==level-1)&
                .or.(pcmask(i,j+1,k)==level-1).or.(pcmask(i,j-1,k)==level-1)&
                .or.(pcmask(i,j,k+1)==level-1).or.(pcmask(i,j,k-1)==level-1)) then
           pcmask(i,j,k)=level
           endif
        endif
     enddo; enddo; enddo
     call ighost_x(u_cmask,2,req(1:4)); call ighost_x(v_cmask,2,req(5:8)); call ighost_x(w_cmask,2,req(9:12))
     call ighost_x(pcmask,2,req(13:16)); call MPI_WAITALL(16,req(1:16),sta(:,1:16),ierr)
     call ighost_y(u_cmask,2,req(1:4)); call ighost_y(v_cmask,2,req(5:8)); call ighost_y(w_cmask,2,req(9:12))
     call ighost_y(pcmask,2,req(13:16)); call MPI_WAITALL(16,req(1:16),sta(:,1:16),ierr)
     call ighost_z(u_cmask,2,req(1:4)); call ighost_z(v_cmask,2,req(5:8)); call ighost_z(w_cmask,2,req(9:12))
     call ighost_z(pcmask,2,req(13:16)); call MPI_WAITALL(16,req(1:16),sta(:,1:16),ierr)
  enddo
!-Debugging
  debug = .false.
  no = 1
  if (debug .and. mod(iout,no)==0) then
     Open(unit=20,FILE=TRIM(out_path)//'/Top_0-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout/no,padding))//'.txt')
     Open(unit=21,FILE=TRIM(out_path)//'/Top_1-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout/no,padding))//'.txt')
     Open(unit=22,FILE=TRIM(out_path)//'/Top_2-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout/no,padding))//'.txt')
     Open(unit=23,FILE=TRIM(out_path)//'/u_ass-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout/no,padding))//'.txt')
     Open(unit=24,FILE=TRIM(out_path)//'/v_ass-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout/no,padding))//'.txt')
     Open(unit=25,FILE=TRIM(out_path)//'/w_ass-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout/no,padding))//'.txt')
     Open(unit=26,FILE=TRIM(out_path)//'/P1-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout/no,padding))//'.txt')
     Open(unit=27,FILE=TRIM(out_path)//'/P2-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout/no,padding))//'.txt')
     Open(unit=28,FILE=TRIM(out_path)//'/P3-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout/no,padding))//'.txt') 
     !j=(js+je)/2
     k=(ks+ke)/2
     !do k=kmin,kmax; do i=imin,imax
     do j=jmin,jmax; do i=imin,imax
        !write outputs for gnuplot 2D
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
        !if (w_cmask(i,j,k)==3) write(25,13)x(i),y(j),zh(k)
        if (pcmask(i,j,k)==3) write(28,13)x(i),y(j),z(k)
     enddo; enddo
     close(unit=20); close(unit=21); close(unit=22); close(unit=23); close(unit=24); close(unit=25)
     close(unit=26); close(unit=27); close(unit=28)
  endif
13 format(3e14.5)
end subroutine set_topology
!-------------------------------------------------------------------------------------------------
subroutine setuppoisson_fs(vof_phase,rhot,dt,A,cvof,n1,n2,n3,kap,iout)    
  use module_grid
  use module_BC
  use module_2phase
  use module_freesurface
  use module_IO
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: rhot
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: kap
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: cvof,n1,n2,n3
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(out) :: A
  integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: vof_phase
  integer :: req(24),sta(MPI_STATUS_SIZE,24)
  real(8) :: alpha2, x_test2, y_test2, z_test2
  real(8) :: nr(3),al3dnew,x0(3),dc(3),FL3DNEW,n_avg(3)
  real(8) :: dt, limit
  integer :: i,j,k,l,iout,nbr,ierr,no
  real(8) :: c1, c0, c_stag

  x_mod=dxh((is+ie)/2); y_mod=dyh((js+je)/2); z_mod=dzh((ks+ke)/2) !assumes an unstretched grid
  P_gx = 0d0; P_gy = 0d0; P_gz = 0d0
  limit = 1d-11
!Debugging
debug = .false.
  no=1
  if (debug .and. mod(iout,no)==0) then
     Open(unit=50,file=TRIM(out_path)//'/C_int-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout/no,padding))//'.txt')
     !Open(unit=51,file=TRIM(out_path)//'/COEFFS-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout/no,padding))//'.txt')
     Open(unit=55,file=TRIM(out_path)//'/C0-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout/no,padding))//'.txt')
     Open(unit=56,file=TRIM(out_path)//'/C1-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout/no,padding))//'.txt')
  endif
  do k=ks,ke; do j=js,je; do i=is,ie
     !----Cav-liquid neighbours, set P_g in cavity cells
     if(vof_phase(i,j,k)==1) then
        !if (cvof(i,j,k)<0.49d0) write(*,'("Vof phase error. Phase test 1, cvof: ",e14.5)')cvof(i,j,k) !debugging
        do nbr=-1,1,2
           if (vof_phase(i+nbr,j,k)==0) then
              P_gx(i,j,k) = sigma*kap(i,j,k)/dx(i) !!filaments and droplets of one cell will be an issue here
              if (P_gx(i,j,k) /= P_gx(i,j,k)) write(*,*)'WARNING, P_g NaN x-dir'
           endif
           if (vof_phase(i,j+nbr,k)==0) then
              P_gy(i,j,k) = sigma*kap(i,j,k)/dy(j)
              if (P_gy(i,j,k) /= P_gy(i,j,k)) write(*,*)'WARNING, P_g NaN y-dir'
           endif
           if (vof_phase(i,j,k+nbr)==0) then
              P_gz(i,j,k) = sigma*kap(i,j,k)/dz(k)
              if (P_gz(i,j,k) /= P_gz(i,j,k)) write(*,*)'WARNING, P_g NaN z-dir'
           endif
        enddo
        !Check x-neighbour         
        if (vof_phase(i+1,j,k) == 0) then
           !if (cvof(i+1,j,k)>0.499d0) write(*,'("Vof phase error. Phase test 0, cvof: ",e14.5)')cvof(i+1,j,k) !debugging
           !get_intersection for liq cell
           nr(1)=n1(i+1,j,k); nr(2)=n2(i+1,j,k); nr(3)=n3(i+1,j,k)
           alpha2 = al3dnew(nr,cvof(i+1,j,k))
           x0(1) = 0d0; 
           x0(2) = 0d0; x0(3) = 0d0
           dc(1) = 0.5d0 
           dc(2) = 1d0; dc(3) = 1d0
           c0 = FL3DNEW(nr,alpha2,x0,dc)
           if (debug .and. k==(ks+ke)/2 .and. mod(iout,no)==0) write(56,314)x(i+1),y(j),c0,cvof(i+1,j,k)
           !get_intersection for cav cell
           nr(1) = n1(i,j,k); nr(2) = n2(i,j,k); nr(3) = n3(i,j,k)
           alpha2 = al3dnew(nr,cvof(i,j,k))
           x0(1) = 0.5d0; 
           x0(2) = 0d0; x0(3) = 0d0
           dc(1) = 0.5d0; 
           dc(2) = 1d0; dc(3) = 1d0
           c1 = FL3DNEW(nr,alpha2,x0,dc)
           if (debug .and. k==(ks+ke)/2 .and. mod(iout,no)==0) write(55,314)x(i),y(j),c1,cvof(i,j,k)
           c_stag = c1+c0
           if (c0>1d-6) then ! use weighted average if liq VOF is not small, otherwise use cav normals
              n_avg(1)=n1(i,j,k)*c1/c_stag + n1(i+1,j,k)*c0/c_stag
              n_avg(2)=n2(i,j,k)*c1/c_stag + n2(i+1,j,k)*c0/c_stag
              n_avg(3)=n3(i,j,k)*c1/c_stag + n3(i+1,j,k)*c0/c_stag
           else
              n_avg(1)= n1(i,j,k)
              n_avg(2)= n2(i,j,k)
              n_avg(3)= n3(i,j,k)
           endif
           alpha2=al3dnew(n_avg,c_stag)
           if (ABS(n_avg(1))>1d-12) then
              x_test2 = (alpha2 - (n_avg(2)+n_avg(3))/2d0)/n_avg(1)
              if (x_test2>limit .and. x_test2<(1d0-limit)) then
                 x_mod(i,j,k) = dxh(i)*(1d0-x_test2)
              else
                 if (x_test2>=(1d0-limit)) x_mod(i,j,k) = limit*dxh(i)
                 if (x_test2<=(limit)) x_mod(i,j,k) = (1d0-limit)*dxh(i)
              endif
           else
              write(*,'("WARNING: x-branch tiny normal",5e14.5,2I8)')&
                   c_stag,x_mod(i,j,k),cvof(i,j,k),cvof(i+1,j,k),n_avg(1),vof_phase(i,j,k),vof_phase(i+1,j,k)
           endif
           if (debug .and. j==(js+je)/2 .and. mod(iout,no)==0) then
              write(50,314)x(i+1)-x_mod(i,j,k),z(k) 
           endif
        endif
        !-------Check y-neighbour
        if (vof_phase(i,j+1,k) == 0) then
           !if (cvof(i,j+1,k)>0.499d0) write(*,'("Vof phase error. Phase test 0, cvof: ",e14.5)')cvof(i,j+1,k) !debugging
           !get_intersection for liq cell
           nr(1) = n1(i,j+1,k); nr(2) = n2(i,j+1,k); nr(3) = n3(i,j+1,k)
           alpha2 = al3dnew(nr,cvof(i,j+1,k))
           x0(1) = 0d0; x0(3) = 0d0 
           x0(2) = 0d0
           dc(1) = 1d0 
           dc(2) = 0.5d0; dc(3) = 1d0
           c0 = FL3DNEW(nr,alpha2,x0,dc)
           if (debug .and. k==(ks+ke)/2 .and. mod(iout,no)==0) write(56,314)x(i),y(j+1),c0,cvof(i,j+1,k)
           !get_intersection for cav cell----------------------------------------
           nr(1) = n1(i,j,k); nr(2) = n2(i,j,k); nr(3) = n3(i,j,k)
           alpha2 = al3dnew(nr,cvof(i,j,k))
           x0(2) = 0.5d0; 
           x0(1) = 0d0; x0(3) = 0d0
           dc(2) = 0.5d0; 
           dc(1) = 1d0; dc(3) = 1d0
           c1 = FL3DNEW(nr,alpha2,x0,dc)
           if (debug .and. k==(ks+ke)/2 .and. mod(iout,no)==0) write(55,314)x(i),y(j),c1,cvof(i,j,k)
           c_stag = c1+c0
           if (c0>1d-6) then 
              n_avg(1)=n1(i,j,k)*c1/c_stag + n1(i,j+1,k)*c0/c_stag
              n_avg(2)=n2(i,j,k)*c1/c_stag + n2(i,j+1,k)*c0/c_stag
              n_avg(3)=n3(i,j,k)*c1/c_stag + n3(i,j+1,k)*c0/c_stag 
           else
              n_avg(1)= n1(i,j,k)
              n_avg(2)= n2(i,j,k)
              n_avg(3)= n3(i,j,k)
           endif
           alpha2=al3dnew(n_avg,c_stag)
           if (ABS(n_avg(2))>1d-12) then
              y_test2 = (alpha2 - (n_avg(1)+n_avg(3))/2d0)/n_avg(2)
              if (y_test2>limit .and. y_test2<(1d0-limit)) then
                 y_mod(i,j,k) = dyh(j)*(1d0-y_test2)
              else
                 if (y_test2>=(1d0-limit)) y_mod(i,j,k) = limit*dyh(j)
                 if (y_test2<=(limit)) y_mod(i,j,k) = (1d0-limit)*dyh(j)
              endif
           else
              write(*,'("WARNING: y-branch tiny normal",5e14.5,2I8)')&
                   c_stag,y_mod(i,j,k),cvof(i,j,k),cvof(i,j+1,k),n_avg(2),vof_phase(i,j,k),vof_phase(i,j+1,k)
           endif
!!$              if (debug .and. k==(ks+ke)/2 .and. mod(iout,no)==0) then
!!$                 write(50,314)x(i),y(j)+nrl*y_mod(i,j+(1-1)/2,k) 
!!$                 write(51,314)x(i),y(j),0d0,nrl*y_mod(i,j+(1-1)/2,k)
!!$              endif
        endif
        !-------Check z-neighbours
        if (vof_phase(i,j,k+1) == 0) then
           !if (cvof(i,j,k+1)>0.499d0) write(*,'("Vof phase error. Phase test 0, cvof: ",e14.5)')cvof(i,j,k+1) !debugging
           !vof fraction in half of liq cell
           nr(1) = n1(i,j,k+1); nr(2) = n2(i,j,k+1); nr(3) = n3(i,j,k+1)
           alpha2 = al3dnew(nr,cvof(i,j,k+1))
           x0(1) = 0d0; x0(2) = 0d0 
           x0(3) = 0d0
           dc(3) = 0.5d0 
           dc(1) = 1d0; dc(2) = 1d0
           c0 = FL3DNEW(nr,alpha2,x0,dc)
           !vof fraction in half of cav cell
           nr(1) = n1(i,j,k); nr(2) = n2(i,j,k); nr(3) = n3(i,j,k)
           alpha2 = al3dnew(nr,cvof(i,j,k))
           x0(1) = 0d0; x0(2) = 0d0 
           x0(3) = 0.5d0
           dc(3) = 0.5d0 
           dc(1) = 1d0; dc(2) = 1d0
           c1 = FL3DNEW(nr,alpha2,x0,dc)
           c_stag = c1+c0
           if (c0>1d-6) then 
              n_avg(1)=n1(i,j,k)*c1/c_stag + n1(i,j,k+1)*c0/c_stag
              n_avg(2)=n2(i,j,k)*c1/c_stag + n2(i,j,k+1)*c0/c_stag
              n_avg(3)=n3(i,j,k)*c1/c_stag + n3(i,j,k+1)*c0/c_stag
           else
              n_avg(1)= n1(i,j,k)
              n_avg(2)= n2(i,j,k)
              n_avg(3)= n3(i,j,k)
           endif
           alpha2=al3dnew(n_avg,c_stag)
           if (ABS(n_avg(3))>1d-12) then
              z_test2 = (alpha2 - (n_avg(1)+n_avg(2))/2d0)/n_avg(3)
              if (z_test2>limit .and. z_test2<(1d0-limit)) then
                 z_mod(i,j,k) = dzh(k)*(1d0-z_test2)
              else
                 if (z_test2>=(1d0-limit)) z_mod(i,j,k+(1-1)/2) = limit*dzh(k+(1-1)/2)
                 if (z_test2<=(limit)) z_mod(i,j,k+(1-1)/2) = (1d0-limit)*dzh(k+(1-1)/2)
              endif
           else
              write(*,'("WARNING: z-branch tiny normal",5e14.5,2I8)')&
                   c_stag,z_mod(i,j,k),cvof(i,j,k),cvof(i,j,k+1),n_avg(3),vof_phase(i,j,k),vof_phase(i,j,k+1)
           endif
           if (debug .and. j==(js+je)/2 .and. mod(iout,no)==0) then
              write(50,314)x(i),z(k+1)-z_mod(i,j,k) 
              !write(51,314)x(i),z(k+1),0d0,-z_mod(i,j,k)
           endif
        endif
     endif
!----Liquid-cavity neighbours
     if (vof_phase(i,j,k)==0) then
        !if (cvof(i,j,k)>=0.499d0) write(*,'("Vof phase error. Phase test 0, cvof: ",e14.5)')cvof(i,j,k) !debugging check
        !Check x-neighbour
        if (vof_phase(i+1,j,k) == 1) then
           !if (cvof(i+1,j,k)<0.499d0) write(*,'("Vof phase error. Phase test 1, cvof: ",e14.5)')cvof(i+1,j,k) !debugging
           !vof fraction in half of cav cell
           nr(1)=n1(i+1,j,k); nr(2)=n2(i+1,j,k); nr(3)=n3(i+1,j,k)
           alpha2 = al3dnew(nr,cvof(i+1,j,k))
           x0(1) = 0d0; 
           x0(2) = 0d0; x0(3) = 0d0
           dc(1) = 0.5d0 
           dc(2) = 1d0; dc(3) = 1d0
           c1 = FL3DNEW(nr,alpha2,x0,dc)
           if (debug .and. k==(ks+ke)/2 .and. mod(iout,no)==0) write(56,314)x(i+1),y(j),c1,cvof(i+1,j,k)
           !vof fraction in half of liq cell
           nr(1) = n1(i,j,k); nr(2) = n2(i,j,k); nr(3) = n3(i,j,k)
           alpha2 = al3dnew(nr,cvof(i,j,k))
           x0(1) = 0.5d0; 
           x0(2) = 0d0; x0(3) = 0d0
           dc(1) = 0.5d0; 
           dc(2) = 1d0; dc(3) = 1d0
           c0 = FL3DNEW(nr,alpha2,x0,dc)
           if (debug .and. k==(ks+ke)/2 .and. mod(iout,no)==0) write(55,314)x(i),y(j),c0,cvof(i,j,k)
           c_stag = c1+c0
           if (c0>1d-6) then ! use weighted average if liq VOF is not small, otherwise use cav normals
              n_avg(1)=n1(i,j,k)*c0/c_stag + n1(i+1,j,k)*c1/c_stag
              n_avg(2)=n2(i,j,k)*c0/c_stag + n2(i+1,j,k)*c1/c_stag
              n_avg(3)=n3(i,j,k)*c0/c_stag + n3(i+1,j,k)*c1/c_stag
           else
              n_avg(1)= n1(i+1,j,k)
              n_avg(2)= n2(i+1,j,k)
              n_avg(3)= n3(i+1,j,k)
           endif
           alpha2=al3dnew(n_avg,c_stag)
           if (ABS(n_avg(1))>1d-12) then
              x_test2 = (alpha2 - (n_avg(2)+n_avg(3))/2d0)/n_avg(1)
              if (x_test2>limit .and. x_test2<(1d0-limit)) then
                 x_mod(i,j,k) = dxh(i)*x_test2
              else
                 if (x_test2>=(1d0-limit)) x_mod(i,j,k) = (1d0-limit)*dxh(i)
                 if (x_test2<=(limit)) x_mod(i,j,k) = limit*dxh(i)
              endif
           else
              write(*,'("WARNING: x-branch tiny normal",5e14.5,2I8)')&
                   c_stag,x_mod(i,j,k),cvof(i,j,k),cvof(i+1,j,k),n_avg(1),vof_phase(i,j,k),vof_phase(i+1,j,k)
           endif
           if (debug .and. j==(js+je)/2 .and. mod(iout,no)==0) then
              write(50,314)x(i)+x_mod(i,j,k),z(k) 
              !write(51,314)x(i),z(k),x_mod(i,j,k),0d0
           endif
        endif
        !-------Check y-neighbours in both directions 
        if (vof_phase(i,j+1,k) == 1) then
           !if (cvof(i,j+1,k)<0.499d0) write(*,'("Vof phase error. Phase test 1, cvof: ",e14.5)')cvof(i,j+1,k) !debugging
           !vof fraction in half of cav cell
           nr(1) = n1(i,j+1,k); nr(2) = n2(i,j+1,k); nr(3) = n3(i,j+1,k)
           alpha2 = al3dnew(nr,cvof(i,j+1,k))
           x0(1) = 0d0; x0(3) = 0d0 
           x0(2) = 0d0
           dc(1) = 1d0 
           dc(2) = 0.5d0; dc(3) = 1d0
           c1 = FL3DNEW(nr,alpha2,x0,dc)
           if (debug .and. k==(ks+ke)/2 .and. mod(iout,no)==0) write(56,314)x(i),y(j+1),c1,cvof(i,j+1,k)
           !vof fraction in half of liq cell----------------------------------------
           nr(1) = n1(i,j,k); nr(2) = n2(i,j,k); nr(3) = n3(i,j,k)
           alpha2 = al3dnew(nr,cvof(i,j,k))
           x0(2) = 0.5d0; 
           x0(1) = 0d0; x0(3) = 0d0
           dc(2) = 0.5d0; 
           dc(1) = 1d0; dc(3) = 1d0
           c0 = FL3DNEW(nr,alpha2,x0,dc)
           if (debug .and. k==(ks+ke)/2 .and. mod(iout,no)==0) write(55,314)x(i),y(j),c0,cvof(i,j,k)
           c_stag = c1+c0
           if (c0>1d-6) then 
              n_avg(1)=n1(i,j,k)*c0/c_stag + n1(i,j+1,k)*c1/c_stag
              n_avg(2)=n2(i,j,k)*c0/c_stag + n2(i,j+1,k)*c1/c_stag
              n_avg(3)=n3(i,j,k)*c0/c_stag + n3(i,j+1,k)*c1/c_stag 
           else
              n_avg(1)= n1(i,j+1,k)
              n_avg(2)= n2(i,j+1,k)
              n_avg(3)= n3(i,j+1,k)
           endif
           alpha2=al3dnew(n_avg,c_stag)
           if (ABS(n_avg(2))>1d-12) then
              y_test2 = (alpha2 - (n_avg(1)+n_avg(3))/2d0)/n_avg(2)
              if (y_test2>limit .and. y_test2<(1d0-limit)) then
                 y_mod(i,j,k) = dyh(j)*y_test2
              else
                 if (y_test2>=(1d0-limit)) y_mod(i,j,k) = (1d0-limit)*dyh(j)
                 if (y_test2<=(limit)) y_mod(i,j,k) = limit*dyh(j)
              endif
           else 
              write(*,'("WARNING: y-branch tiny normal",5e14.5,2I8)')&
                   c_stag,y_mod(i,j,k),cvof(i,j,k),cvof(i,j+1,k),n_avg(2),vof_phase(i,j,k),vof_phase(i,j+1,k)
           endif
!!$              if (debug .and. k==(ks+ke)/2 .and. mod(iout,no)==0) then
!!$                 write(50,314)x(i),y(j)+nrl*y_mod(i,j+(1-1)/2,k) 
!!$                 write(51,314)x(i),y(j),0d0,nrl*y_mod(i,j+(1-1)/2,k)
!!$              endif
        endif
        !-------Check z-neighbours
        if (vof_phase(i,j,k+1) == 1) then
           !if (cvof(i,j,k+1)<0.499d0) write(*,'("Vof phase error. Phase test 1, cvof: ",e14.5)')cvof(i,j,k+1) !debugging
           !vof fraction in half of cav cell
           nr(1) = n1(i,j,k+1); nr(2) = n2(i,j,k+1); nr(3) = n3(i,j,k+1)
           alpha2 = al3dnew(nr,cvof(i,j,k+1))
           x0(1) = 0d0; x0(2) = 0d0 
           x0(3) = 0d0
           dc(3) = 0.5d0 
           dc(1) = 1d0; dc(2) = 1d0
           c1 = FL3DNEW(nr,alpha2,x0,dc)
           !vof fraction in half of liq cell
           nr(1) = n1(i,j,k); nr(2) = n2(i,j,k); nr(3) = n3(i,j,k)
           alpha2 = al3dnew(nr,cvof(i,j,k))
           x0(1) = 0d0; x0(2) = 0d0 
           x0(3) = 0.5d0
           dc(3) = 0.5d0 
           dc(1) = 1d0; dc(2) = 1d0
           c0 = FL3DNEW(nr,alpha2,x0,dc)
           c_stag = c1+c0
           if (c0>1d-6) then 
              n_avg(1)=n1(i,j,k)*c0/c_stag + n1(i,j,k+1)*c1/c_stag
              n_avg(2)=n2(i,j,k)*c0/c_stag + n2(i,j,k+1)*c1/c_stag
              n_avg(3)=n3(i,j,k)*c0/c_stag + n3(i,j,k+1)*c1/c_stag
           else
              n_avg(1)= n1(i,j,k+1)
              n_avg(2)= n2(i,j,k+1)
              n_avg(3)= n3(i,j,k+1)
           endif
           alpha2=al3dnew(n_avg,c_stag)
           if (ABS(n_avg(3))>1d-12) then
              z_test2 = (alpha2 - (n_avg(1)+n_avg(2))/2d0)/n_avg(3)
              if (z_test2>limit .and. z_test2<(1d0-limit)) then
                 z_mod(i,j,k) = dzh(k)*z_test2
              else
                 if (z_test2>=(1d0-limit)) z_mod(i,j,k) = (1d0-limit)*dzh(k)
                 if (z_test2<=(limit)) z_mod(i,j,k) = limit*dzh(k)
              endif
           else
              write(*,'("WARNING: z-branch tiny normal",5e14.5,2I8)')&
                   c_stag,z_mod(i,j,k),cvof(i,j,k),cvof(i,j,k+1),n_avg(3),vof_phase(i,j,k),vof_phase(i,j,k+1)
           endif
           if (debug .and. j==(js+je)/2 .and. mod(iout,no)==0) then
              write(50,314)x(i),z(k)+z_mod(i,j,k) 
              !write(51,314)x(i),z(k),0d0,z_mod(i,j,k)
           endif
        endif
     endif
  enddo;enddo;enddo
!--Debugging
  call ghost_x(P_gx,1,req(1:4)); call ghost_y(P_gy,1,req(5:8)); call ghost_z(P_gz,1,req(9:12)) 
  call ghost_x(x_mod,1,req(13:16)); call ghost_y(y_mod,1,req(17:20)); call ghost_z(z_mod,1,req(21:24)) 
  call MPI_WAITALL(24,req(1:24),sta(:,1:24),ierr)
!!$  call ghost_x(P_gy,1,req(1:4)); call ghost_y(P_gy,1,req(5:8)); call ghost_z(P_gy,1,req(9:12)) 
!!$  call ghost_x(y_mod,1,req(13:16)); call ghost_y(y_mod,1,req(17:20)); call ghost_z(y_mod,1,req(21:24)) 
!!$  call MPI_WAITALL(24,req(1:24),sta(:,1:24),ierr)
!!$  call ghost_x(P_gz,1,req(1:4)); call ghost_y(P_gz,1,req(5:8)); call ghost_z(P_gz,1,req(9:12)) 
!!$  call ghost_x(z_mod,1,req(13:16)); call ghost_y(z_mod,1,req(17:20)); call ghost_z(z_mod,1,req(21:24)) 
!!$  call MPI_WAITALL(24,req(1:24),sta(:,1:24),ierr)
  if (debug .and. mod(iout,no)==0) then
     Open(unit=52,file=TRIM(out_path)//'/P_int1-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout/no,padding))//'.txt')
     Open(unit=53,file=TRIM(out_path)//'/P_int2-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout/no,padding))//'.txt')
     Open(unit=54,file=TRIM(out_path)//'/P_int3-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout/no,padding))//'.txt') 
!!$     Open(unit=57,file=TRIM(out_path)//'/x_mod-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout/no,padding))//'.txt')
!!$     Open(unit=58,file=TRIM(out_path)//'/y_mod-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout/no,padding))//'.txt')
!!$     Open(unit=59,file=TRIM(out_path)//'/z_mod-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout/no,padding))//'.txt')
     Open(unit=60,file=TRIM(out_path)//'/phase-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout/no,padding))//'.txt')
     Open(unit=61,file=TRIM(out_path)//'/kappa-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout/no,padding))//'.txt')
     Open(unit=62,file=TRIM(out_path)//'/n1-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout/no,padding))//'.txt')
     Open(unit=63,file=TRIM(out_path)//'/n2-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout/no,padding))//'.txt')
     Open(unit=64,file=TRIM(out_path)//'/n3-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout/no,padding))//'.txt') 
     j=(js+je)/2
     !k=(ks+ke)/2
     do k=ks-1,ke+1; do i=is-1,ie+1
        write(52,313)x(i),z(k),P_gx(i,j,k)
        !write(53,313)x(i),y(j),P_gy(i,j,k)
        write(54,313)x(i),z(k),P_gz(i,j,k)
        !write(57,313)xh(i),y(j),x_mod(i,j,k)
        !write(58,313)x(i),y(j),y_mod(i,j,k)
        !write(59,313)x(i),zh(k),z_mod(i,j,k)
        write(60,212)x(i),z(k),vof_phase(i,j,k)
        write(61,313)x(i),z(k),kap(i,j,k)
        Write(62,313)x(i),z(k),n1(i,j,k)
        Write(63,313)x(i),z(k),n2(i,j,k)
        Write(64,313)x(i),z(k),n3(i,j,k)
     enddo; enddo   
     close(unit=50);
     close(unit=52); close(unit=53); close(unit=54)
!!$     close(unit=55); close(unit=56); 
!!$     close(unit=57); close(unit=58); close(unit=59)
     close(unit=60); close(unit=61); close(unit=62); close(unit=63)
     close(unit=64)
  endif
212 format(2e14.5,I8)
312 format(2e14.5) !remove, debugging
313 format(3e14.5) !remove, debugging
323 format(2e14.5,e14.4) !remove, debugging
314 format(4e14.5) !remove, debugging
414 format(3e14.5,I8)
!--------------------------------------------------------------------------------------------------------
  do k=ks,ke; do j=js,je; do i=is,ie
     if (vof_phase(i,j,k)==0) then
        A(i,j,k,1) = dt/(dx(i)*x_mod(i-1,j,k)*(rhot(i,j,k)))
        if (A(i,j,k,1) /= A(i,j,k,1)) write(*,'("ERROR: A1 NaN :",2e14.5)')A(i,j,k,1),x_mod(i-1,j,k) 
        if (A(i,j,k,1)>1d8) then
           write(*,'("Large A1",2e14.5)')A(i,j,k,1),x_mod(i-1,j,k)
           !A(i,j,k,1) = 1d8
        endif
        A(i,j,k,2) = dt/(dx(i)*x_mod(i,j,k)*(rhot(i,j,k)))
        if (A(i,j,k,2) /= A(i,j,k,2)) write(*,'("ERROR: A2 NaN :",2e14.5)')A(i,j,k,2),x_mod(i,j,k)
        if (A(i,j,k,2)>1d8) then
           write(*,'("Large A2",2e14.5)')A(i,j,k,2),x_mod(i,j,k)
           !A(i,j,k,2) = 1d8
        endif
        A(i,j,k,3) = dt/(dy(j)*y_mod(i,j-1,k)*(rhot(i,j,k)))
        if (A(i,j,k,3) /= A(i,j,k,3)) write(*,'("ERROR: A3 NaN :",2e14.5)')A(i,j,k,3),y_mod(i,j-1,k)
        if (A(i,j,k,3)>1d8) then
           write(*,'("Large A3",2e14.5)')A(i,j,k,3),y_mod(i,j-1,k)
           !A(i,j,k,3) = 1d8
        endif
        A(i,j,k,4) = dt/(dy(j)*y_mod(i,j,k)*(rhot(i,j,k)))
        if (A(i,j,k,4) /= A(i,j,k,4)) write(*,'("ERROR: A4 NaN :",2e14.5)')A(i,j,k,4),y_mod(i,j,k)
        if (A(i,j,k,4)>1d8) then
           write(*,'("Large A4",2e14.5)')A(i,j,k,4),y_mod(i,j,k)
           !A(i,j,k,4) = 1d8
        endif
        A(i,j,k,5) = dt/(dz(k)*z_mod(i,j,k-1)*(rhot(i,j,k)))
        if (A(i,j,k,5) /= A(i,j,k,5)) write(*,'("ERROR: A5 NaN :",2e14.5)')A(i,j,k,5),z_mod(i,j,k-1)
        if (A(i,j,k,5)>1d8) then
           write(*,'("Large A5",2e14.5)')A(i,j,k,5),z_mod(i,j,k-1)
           !A(i,j,k,5) = 1d8
        endif
        A(i,j,k,6) = dt/(dz(k)*z_mod(i,j,k)*(rhot(i,j,k)))
        if (A(i,j,k,6) /= A(i,j,k,6)) write(*,'("ERROR: A6 NaN :",2e14.5)')A(i,j,k,6),z_mod(i,j,k)
        if (A(i,j,k,6)>1d8) then
           write(*,'("Large A6",2e14.5)')A(i,j,k,6),z_mod(i,j,k)
           !A(i,j,k,6) = 1d8
        endif
        A(i,j,k,7) = sum(A(i,j,k,1:6))
        A(i,j,k,8) = A(i,j,k,8) + dt/rhot(i,j,k)*&
             (P_gx(i+1,j,k)/(dx(i)*x_mod(i,j,k))+P_gx(i-1,j,k)/(dx(i)*x_mod(i-1,j,k))&
             +P_gy(i,j+1,k)/(dy(j)*y_mod(i,j,k))+P_gy(i,j-1,k)/(dy(j)*y_mod(i,j-1,k))&
             +P_gz(i,j,k+1)/(dz(k)*z_mod(i,j,k))+P_gz(i,j,k-1)/(dz(k)*z_mod(i,j,k-1)))
        if (A(i,j,k,8) /= A(i,j,k,8)) then
           write(*,'("A8 NaN, error imminent. Neigbours mods :",6e14.5)')x_mod(i-1,j,k),x_mod(i,j,k),&
                y_mod(i,j-1,k),y_mod(i,j,k),z_mod(i,j,k-1),z_mod(i,j,k)
           write(*,'("A8 NaN, error imminent. P_g :",6e14.5,2I8)')P_gx(i-1,j,k),P_gx(i+1,j,k),&
                P_gy(i,j-1,k),P_gy(i,j+1,k),P_gz(i,j,k-1),P_gz(i,j,k+1),i,k
           write(*,'("rho :",e14.5)')rhot(i,j,k)
        endif
        if (debug .and. mod(iout,no)==0 .and. j==(js+je)/2) then
           Open(unit=51,file=TRIM(out_path)//'/COEFFS-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout/no,padding))//'.txt')
           write(51,314)x(i),z(k),-x_mod(i-1,j,k),0d0
           write(51,314)x(i),z(k),x_mod(i,j,k),0d0
           write(51,314)x(i),z(k),0d0,-z_mod(i,j,k-1)
           write(51,314)x(i),z(k),0d0,z_mod(i,j,k)
        endif
     endif
  enddo;enddo;enddo
end subroutine setuppoisson_fs
!--------------------------------------------------------------------------------------------------------------------
subroutine setuppoisson_fs2(utmp,vtmp,wtmp,dt,A,vof_phase,istep)
  use module_grid
  use module_BC
  use module_2phase
  use module_freesurface
  use module_IO
  implicit none
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(out) :: A
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: utmp,vtmp,wtmp
  integer(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: vof_phase
  real(8) :: dt,h_mod
  integer :: i,j,k,nbr
  integer :: istep

  do k=ks,ke; do j=js,je; do i=is,ie
     if (pcmask(i,j,k)==1 .or. pcmask(i,j,k)==2) then !rho is 1d0
        A(i,j,k,1) = dt/(dx(i)*dxh(i-1))
        if (A(i,j,k,1) /= A(i,j,k,1)) write(*,'("ERROR: A1 NaN in fs2:",e14.5)')A(i,j,k,1)
        A(i,j,k,2) = dt/(dx(i)*dxh(i))
        if (A(i,j,k,2) /= A(i,j,k,2)) write(*,'("ERROR: A2 NaN in fs2:",e14.5)')A(i,j,k,2)
        A(i,j,k,3) = dt/(dy(j)*dyh(j-1))
        if (A(i,j,k,3) /= A(i,j,k,3)) write(*,'("ERROR: A3 NaN in fs2:",e14.5)')A(i,j,k,3)
        A(i,j,k,4) = dt/(dy(j)*dyh(j))
        if (A(i,j,k,4) /= A(i,j,k,4)) write(*,'("ERROR: A4 NaN in fs2:",e14.5)')A(i,j,k,4)
        A(i,j,k,5) = dt/(dz(k)*dzh(k-1))
        if (A(i,j,k,5) /= A(i,j,k,5)) write(*,'("ERROR: A5 NaN in fs2:",e14.5)')A(i,j,k,5)
        A(i,j,k,6) = dt/(dz(k)*dzh(k))
        if (A(i,j,k,6) /= A(i,j,k,6)) write(*,'("ERROR: A6 NaN in fs2:",e14.5)')A(i,j,k,6)

!!$     if (vof_phase(i,j,k)==1) then
!!$        do nbr=-1,1,2
!!$           if (vof_phase(i+nbr,j,k) == 0) then
!!$              h_mod = dxh(i+(nbr-1)/2)-x_mod(i+(nbr-1)/2,j,k)
!!$              if (h_mod > (1d-10)) A(i,j,k,2+(nbr-1)/2) = dt/(dx(i)*h_mod)
!!$           endif
!!$           if (vof_phase(i,j+nbr,k) == 0) then
!!$              h_mod = dyh(j+(nbr-1)/2)-y_mod(i,j+(nbr-1)/2,k)
!!$              if (h_mod > (1d-10)) A(i,j,k,4+(nbr-1)/2) = dt/(dy(j)*h_mod)
!!$           endif
!!$           if (vof_phase(i,j,k+nbr) == 0) then
!!$              h_mod = dzh(k+(nbr-1)/2)-z_mod(i,j,k+(nbr-1)/2)
!!$              if (h_mod > (1d-10)) A(i,j,k,6+(nbr-1)/2) = dt/(dz(k)*h_mod)
!!$           endif
!!$        enddo
!!$     endif
        A(i,j,k,7) = sum(A(i,j,k,1:6))
        A(i,j,k,8) =  -((utmp(i,j,k)-utmp(i-1,j,k))/dx(i) &
             +  (vtmp(i,j,k)-vtmp(i,j-1,k))/dy(j) &
             +  (wtmp(i,j,k)-wtmp(i,j,k-1))/dz(k))
        if (A(i,j,k,8) /= A(i,j,k,8)) write(*,'("ERROR: A8 NaN in fs2:",e14.5)')A(i,j,k,8) 
     endif
  enddo; enddo; enddo
end subroutine setuppoisson_fs2
!--------------------------------------------------------------------------------------------------------------------
subroutine discrete_divergence(u,v,w,iout)
  use module_grid
  use module_freesurface
  use module_IO
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: u,v,w 
  real(8), dimension(is:ie,js:je,ks:ke) :: div, t
  real(8), dimension(0:3) :: divtot, domain
  integer :: i,j,k,l,iout,ierr

!OPEN(unit=20,file=TRIM(out_path)//'/divergence-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout,padding))//'.txt')
OPEN(unit=21,file='div_type.txt',access='append')

divtot = 0d0

do k=ks,ke; do j=js,je; do i=is,ie
   div(i,j,k)=(u(i-1,j,k)-u(i,j,k))*dy(j)*dz(k)+(v(i,j-1,k)-v(i,j,k))*dx(i)*dz(k)+(w(i,j,k-1)-w(i,j,k))*dx(i)*dy(j)
   !write(20,14)x(i),y(j),z(k),div(i,j,k)
   do l=0,3
      if (pcmask(i,j,k)==l) divtot(l)=divtot(l)+ABS(div(i,j,k))
   enddo
enddo; enddo; enddo
divtot=divtot/(Nx*Ny*Nz)
call MPI_ALLREDUCE(divtot,domain,4, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_Cart, ierr)
write(21,15)iout,domain(0),domain(1),domain(2),domain(3)
!close(unit=20)
close(unit=21)
14 format(4e14.5)
15 format(I8,4e14.5)
end subroutine discrete_divergence
!--------------------------------------------------------------------------------------------------------------------
