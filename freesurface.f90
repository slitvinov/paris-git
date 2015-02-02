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
  !call set_topology_bc
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

end subroutine set_topology
!-------------------------------------------------------------------------------------------------
subroutine setuppoisson_fs(utmp,vtmp,wtmp,vof_phase,rhot,dt,A,cvof,n1,n2,n3,kap,iout)    
  use module_grid
  use module_BC
  use module_2phase
  use module_freesurface
  use module_IO
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: utmp,vtmp,wtmp,rhot
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: kap
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: cvof,n1,n2,n3
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(inout) :: A
  integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: vof_phase
  integer :: req(24),sta(MPI_STATUS_SIZE,24)
  real(8) :: alpha2, x_test2, y_test2, z_test2
  real(8) :: nr(3),al3dnew,x0(3),dc(3),FL3DNEW,n_avg(3)
  real(8) :: dt, limit
  integer :: i,j,k,l,iout,nbr,ierr,no
  real(8) :: c1, c0, c_stag, c_min
  real(8), dimension(4) :: P_bc

  x_mod=dxh((is+ie)/2); y_mod=dyh((js+je)/2); z_mod=dzh((ks+ke)/2) !assumes an unstretched grid
  P_gx = 0d0; P_gy = 0d0; P_gz = 0d0
  limit = 1d-3/dx((is+ie)/2)
  c_min = 1d-2

  do k=ks,ke; do j=js,je; do i=is,ie
     !----Cav-liquid neighbours, set P_g in cavity cells
     if(vof_phase(i,j,k)==1) then
        if (vof_phase(i+1,j,k)==0) then
           P_gx(i,j,k) = sigma*kap(i,j,k)/dx(i) !!filaments and droplets of one cell will be an issue here
           if (P_gx(i,j,k) /= P_gx(i,j,k)) write(*,*)'WARNING, P_g NaN x-dir'
        endif
        if (vof_phase(i,j+1,k)==0) then
           P_gy(i,j,k) = sigma*kap(i,j,k)/dy(j)
           if (P_gy(i,j,k) /= P_gy(i,j,k)) write(*,*)'WARNING, P_g NaN y-dir'
        endif
        if (vof_phase(i,j,k+1)==0) then
           P_gz(i,j,k) = sigma*kap(i,j,k)/dz(k)
           if (P_gz(i,j,k) /= P_gz(i,j,k)) write(*,*)'WARNING, P_g NaN z-dir'
        endif

        !Check x-neighbour         
        if (vof_phase(i+1,j,k) == 0) then
           !get_intersection for liq cell
           nr(1)=n1(i+1,j,k); nr(2)=n2(i+1,j,k); nr(3)=n3(i+1,j,k)
           alpha2 = al3dnew(nr,cvof(i+1,j,k))
           x0(1) = 0d0; 
           x0(2) = 0d0; x0(3) = 0d0
           dc(1) = 0.5d0 
           dc(2) = 1d0; dc(3) = 1d0
           c0 = FL3DNEW(nr,alpha2,x0,dc)
           !get_intersection for cav cell
           nr(1) = n1(i,j,k); nr(2) = n2(i,j,k); nr(3) = n3(i,j,k)
           alpha2 = al3dnew(nr,cvof(i,j,k))
           x0(1) = 0.5d0; 
           x0(2) = 0d0; x0(3) = 0d0
           dc(1) = 0.5d0; 
           dc(2) = 1d0; dc(3) = 1d0
           c1 = FL3DNEW(nr,alpha2,x0,dc)
           c_stag = c1+c0
           if (c0>c_min) then
              nr(1)=n1(i+1,j,k)*c0/c_stag + n1(i,j,k)*c1/c_stag
              nr(2)=n2(i+1,j,k)*c0/c_stag + n2(i,j,k)*c1/c_stag
              nr(3)=n3(i+1,j,k)*c0/c_stag + n3(i,j,k)*c1/c_stag
              n_avg(1) = nr(1)/(ABS(nr(1))+ABS(nr(2))+ABS(nr(3)))
              n_avg(2) = nr(2)/(ABS(nr(1))+ABS(nr(2))+ABS(nr(3)))
              n_avg(3) = nr(3)/(ABS(nr(1))+ABS(nr(2))+ABS(nr(3)))
           else
              n_avg(1)= n1(i,j,k)
              n_avg(2)= n2(i,j,k)
              n_avg(3)= n3(i,j,k)
           endif
           if (ABS((ABS(n_avg(1))+ABS(n_avg(2))+ABS(n_avg(3)))-1d0) > 1d-12) then
              write(*,*)'Normals not normalised'
              call pariserror('Normals not normalised')
           endif
           alpha2=al3dnew(n_avg,c_stag)
           if (ABS(n_avg(1))>1d-12) then
              x_test2 = (alpha2 - (n_avg(2)+n_avg(3))/2d0)/n_avg(1)
              if (x_test2 < 0.5d0) P_gx(i,j,k) = sigma*kap(i+1,j,k)/dx(i+1)
              x_mod(i,j,k) = dxh(i)*(1d0-x_test2)
              if (x_mod(i,j,k)>dxh(i)) x_mod(i,j,k) = dxh(i)
              if (x_mod(i,j,k)<limit*dxh(i)) x_mod(i,j,k) = limit*dxh(i)
              if (x_mod(i,j,k) /= x_mod(i,j,k)) then
                 write(*,'("x_mod NaN. c_st, n, vofs, phases:",6e14.5,5I8)')&
                      c_stag,n_avg(1),n_avg(2),n_avg(3),cvof(i,j,k),cvof(i+1,j,k),vof_phase(i,j,k),vof_phase(i+1,j,k),i,j,k
              endif
           endif
        endif
        !-------Check y-neighbour
        if (vof_phase(i,j+1,k) == 0) then
           !get_intersection for liq cell
           nr(1) = n1(i,j+1,k); nr(2) = n2(i,j+1,k); nr(3) = n3(i,j+1,k)
           alpha2 = al3dnew(nr,cvof(i,j+1,k))
           x0(1) = 0d0; x0(3) = 0d0 
           x0(2) = 0d0
           dc(1) = 1d0 
           dc(2) = 0.5d0; dc(3) = 1d0
           c0 = FL3DNEW(nr,alpha2,x0,dc)
           !get_intersection for cav cell----------------------------------------
           nr(1) = n1(i,j,k); nr(2) = n2(i,j,k); nr(3) = n3(i,j,k)
           alpha2 = al3dnew(nr,cvof(i,j,k))
           x0(2) = 0.5d0; 
           x0(1) = 0d0; x0(3) = 0d0
           dc(2) = 0.5d0; 
           dc(1) = 1d0; dc(3) = 1d0
           c1 = FL3DNEW(nr,alpha2,x0,dc)
           c_stag = c1+c0
           if (c0>c_min) then ! use weighted average if liq VOF is not small, otherwise use cav normals
              nr(1)=n1(i,j+1,k)*c0/c_stag + n1(i,j,k)*c1/c_stag
              nr(2)=n2(i,j+1,k)*c0/c_stag + n2(i,j,k)*c1/c_stag
              nr(3)=n3(i,j+1,k)*c0/c_stag + n3(i,j,k)*c1/c_stag
              n_avg(1) = nr(1)/(ABS(nr(1))+ABS(nr(2))+ABS(nr(3)))
              n_avg(2) = nr(2)/(ABS(nr(1))+ABS(nr(2))+ABS(nr(3)))
              n_avg(3) = nr(3)/(ABS(nr(1))+ABS(nr(2))+ABS(nr(3)))
           else
              n_avg(1)= n1(i,j,k)
              n_avg(2)= n2(i,j,k)
              n_avg(3)= n3(i,j,k)
           endif
           if (ABS((ABS(n_avg(1))+ABS(n_avg(2))+ABS(n_avg(3)))-1d0) > 1d-12) then
              write(*,*)'Normals not normalised'
              call pariserror('Normals not normalised')
           endif
           alpha2=al3dnew(n_avg,c_stag)
           if (ABS(n_avg(2))>1d-12) then
              y_test2 = (alpha2 - (n_avg(1)+n_avg(3))/2d0)/n_avg(2)
              if (y_test2 < 0.5d0) P_gy(i,j,k) = sigma*kap(i,j+1,k)/dy(j+1)
              y_mod(i,j,k) = dyh(j)*(1d0-y_test2)
              if (y_mod(i,j,k)>dyh(j)) y_mod(i,j,k) = dyh(j)
              if (y_mod(i,j,k)<limit*dyh(j)) y_mod(i,j,k) = limit*dyh(j)
              if (y_mod(i,j,k) /= y_mod(i,j,k)) write(*,'("y_mod NaN. c_st, n, vofs, phases:",4e14.5,2I8)')&
                   c_stag,n_avg(2),cvof(i,j,k),cvof(i,j+1,k),vof_phase(i,j,k),vof_phase(i,j+1,k)
           endif
        endif
        !-------Check z-neighbours
        if (vof_phase(i,j,k+1) == 0) then
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
           if (c0>c_min) then ! use weighted average if liq VOF is not small, otherwise use cav normals
              nr(1)=n1(i,j,k+1)*c0/c_stag + n1(i,j,k)*c1/c_stag
              nr(2)=n2(i,j,k+1)*c0/c_stag + n2(i,j,k)*c1/c_stag
              nr(3)=n3(i,j,k+1)*c0/c_stag + n3(i,j,k)*c1/c_stag
              n_avg(1) = nr(1)/(ABS(nr(1))+ABS(nr(2))+ABS(nr(3)))
              n_avg(2) = nr(2)/(ABS(nr(1))+ABS(nr(2))+ABS(nr(3)))
              n_avg(3) = nr(3)/(ABS(nr(1))+ABS(nr(2))+ABS(nr(3)))
           else
              n_avg(1)= n1(i,j,k)
              n_avg(2)= n2(i,j,k)
              n_avg(3)= n3(i,j,k)
           endif
           if (ABS((ABS(n_avg(1))+ABS(n_avg(2))+ABS(n_avg(3)))-1d0) > 1d-12) then
              write(*,*)'Normals not normalised'
              call pariserror('Normals not normalised')
           endif
           alpha2=al3dnew(n_avg,c_stag)
           if (ABS(n_avg(3))>1d-12) then
              z_test2 = (alpha2 - (n_avg(1)+n_avg(2))/2d0)/n_avg(3)
              if (z_test2 < 0.5d0) P_gz(i,j,k) = sigma*kap(i,j,k+1)/dz(k+1)
              z_mod(i,j,k) = dzh(k)*(1d0-z_test2)
              if (z_mod(i,j,k)>dzh(k)) z_mod(i,j,k) = dzh(k)
              if (z_mod(i,j,k)<limit*dzh(k)) z_mod(i,j,k) = limit*dzh(k)
              if (z_mod(i,j,k) /= z_mod(i,j,k)) write(*,'("z_mod NaN. c_st, n, vofs, phases:",4e14.5,2I8)')&
                   c_stag,n_avg(3),cvof(i,j,k),cvof(i,j,k+1),vof_phase(i,j,k),vof_phase(i,j,k+1)
           endif
        endif
     endif
!----Liquid-cavity neighbours-----------------------------------------------------
     if (vof_phase(i,j,k)==0) then
        A(i,j,k,1) = 2d0*dt/(dx(i)*dxh(i-1)*(rhot(i-1,j,k)+rhot(i,j,k)))
        A(i,j,k,2) = 2d0*dt/(dx(i)*dxh(i  )*(rhot(i+1,j,k)+rhot(i,j,k)))
        A(i,j,k,3) = 2d0*dt/(dy(j)*dyh(j-1)*(rhot(i,j-1,k)+rhot(i,j,k)))
        A(i,j,k,4) = 2d0*dt/(dy(j)*dyh(j  )*(rhot(i,j+1,k)+rhot(i,j,k)))
        A(i,j,k,5) = 2d0*dt/(dz(k)*dzh(k-1)*(rhot(i,j,k-1)+rhot(i,j,k)))
        A(i,j,k,6) = 2d0*dt/(dz(k)*dzh(k  )*(rhot(i,j,k+1)+rhot(i,j,k)))
        A(i,j,k,7) = sum(A(i,j,k,1:6))
        A(i,j,k,8) =  -( (utmp(i,j,k)-utmp(i-1,j,k))/dx(i) &
             +  (vtmp(i,j,k)-vtmp(i,j-1,k))/dy(j) &
             +  (wtmp(i,j,k)-wtmp(i,j,k-1))/dz(k) )
        ! set Laplace pressure jumps
        if (vof_phase(i+1,j,k)==1) then
           P_gx(i+1,j,k) = sigma*kap(i+1,j,k)/dx(i) !!filaments and droplets of one cell will be an issue here
           if (P_gx(i+1,j,k) /= P_gx(i+1,j,k)) write(*,*)'WARNING, P_g NaN x-dir'
        endif
        if (vof_phase(i,j+1,k)==1) then
           P_gy(i,j+1,k) = sigma*kap(i,j+1,k)/dy(j)
           if (P_gy(i,j+1,k) /= P_gy(i,j+1,k)) write(*,*)'WARNING, P_g NaN y-dir'
        endif
        if (vof_phase(i,j,k+1)==1) then
           P_gz(i,j,k+1) = sigma*kap(i,j,k+1)/dz(k)
           if (P_gz(i,j,k+1) /= P_gz(i,j,k+1)) write(*,*)'WARNING, P_g NaN z-dir'
        endif
        !Check x-neighbour
        if (vof_phase(i+1,j,k) == 1) then
           !vof fraction in half of cav cell
           nr(1)=n1(i+1,j,k); nr(2)=n2(i+1,j,k); nr(3)=n3(i+1,j,k)
           alpha2 = al3dnew(nr,cvof(i+1,j,k))
           x0(1) = 0d0; 
           x0(2) = 0d0; x0(3) = 0d0
           dc(1) = 0.5d0 
           dc(2) = 1d0; dc(3) = 1d0
           c1 = FL3DNEW(nr,alpha2,x0,dc)
           !vof fraction in half of liq cell
           nr(1) = n1(i,j,k); nr(2) = n2(i,j,k); nr(3) = n3(i,j,k)
           alpha2 = al3dnew(nr,cvof(i,j,k))
           x0(1) = 0.5d0; 
           x0(2) = 0d0; x0(3) = 0d0
           dc(1) = 0.5d0; 
           dc(2) = 1d0; dc(3) = 1d0
           c0 = FL3DNEW(nr,alpha2,x0,dc)
           c_stag = c1+c0
           if (c0>c_min) then ! use weighted average if liq VOF is not small, otherwise use cav normals
              nr(1)=n1(i,j,k)*c0/c_stag + n1(i+1,j,k)*c1/c_stag
              nr(2)=n2(i,j,k)*c0/c_stag + n2(i+1,j,k)*c1/c_stag
              nr(3)=n3(i,j,k)*c0/c_stag + n3(i+1,j,k)*c1/c_stag
              n_avg(1) = nr(1)/(ABS(nr(1))+ABS(nr(2))+ABS(nr(3)))
              n_avg(2) = nr(2)/(ABS(nr(1))+ABS(nr(2))+ABS(nr(3)))
              n_avg(3) = nr(3)/(ABS(nr(1))+ABS(nr(2))+ABS(nr(3)))
           else
              n_avg(1)= n1(i+1,j,k)
              n_avg(2)= n2(i+1,j,k)
              n_avg(3)= n3(i+1,j,k)
           endif
           if (ABS((ABS(n_avg(1))+ABS(n_avg(2))+ABS(n_avg(3)))-1d0) > 1d-12) then
              write(*,*)'Normals not normalised'
              call pariserror('Normals not normalised')
           endif
           alpha2=al3dnew(n_avg,c_stag)
           if (ABS(n_avg(1))>1d-12) then
              x_test2 = (alpha2 - (n_avg(2)+n_avg(3))/2d0)/n_avg(1)
              if (x_test2 < 0.5d0) P_gx(i+1,j,k) = sigma*kap(i,j,k)/dx(i)
              x_mod(i,j,k) = dxh(i)*x_test2
              if (x_mod(i,j,k)>dxh(i)) x_mod(i,j,k) = dxh(i)
              if (x_mod(i,j,k)<limit*dxh(i)) x_mod(i,j,k) = limit*dxh(i)
              if (x_mod(i,j,k) /= x_mod(i,j,k)) write(*,'("x_mod NaN. c_st, n, vofs, phases:",4e14.5,2I8)')&
                   c_stag,n_avg(1),cvof(i,j,k),cvof(i+1,j,k),vof_phase(i,j,k),vof_phase(i+1,j,k)
           endif
        endif
        !-------Check y-neighbours in both directions 
        if (vof_phase(i,j+1,k) == 1) then
           !vof fraction in half of cav cell
           nr(1) = n1(i,j+1,k); nr(2) = n2(i,j+1,k); nr(3) = n3(i,j+1,k)
           alpha2 = al3dnew(nr,cvof(i,j+1,k))
           x0(1) = 0d0; x0(3) = 0d0 
           x0(2) = 0d0
           dc(1) = 1d0 
           dc(2) = 0.5d0; dc(3) = 1d0
           c1 = FL3DNEW(nr,alpha2,x0,dc)
           !vof fraction in half of liq cell----------------------------------------
           nr(1) = n1(i,j,k); nr(2) = n2(i,j,k); nr(3) = n3(i,j,k)
           alpha2 = al3dnew(nr,cvof(i,j,k))
           x0(2) = 0.5d0; 
           x0(1) = 0d0; x0(3) = 0d0
           dc(2) = 0.5d0; 
           dc(1) = 1d0; dc(3) = 1d0
           c0 = FL3DNEW(nr,alpha2,x0,dc)
           c_stag = c1+c0
           if (c0>c_min) then ! use weighted average if liq VOF is not small, otherwise use cav normals
              nr(1)=n1(i,j,k)*c0/c_stag + n1(i,j+1,k)*c1/c_stag
              nr(2)=n2(i,j,k)*c0/c_stag + n2(i,j+1,k)*c1/c_stag
              nr(3)=n3(i,j,k)*c0/c_stag + n3(i,j+1,k)*c1/c_stag
              n_avg(1) = nr(1)/(ABS(nr(1))+ABS(nr(2))+ABS(nr(3)))
              n_avg(2) = nr(2)/(ABS(nr(1))+ABS(nr(2))+ABS(nr(3)))
              n_avg(3) = nr(3)/(ABS(nr(1))+ABS(nr(2))+ABS(nr(3)))
           else
              n_avg(1)= n1(i,j+1,k)
              n_avg(2)= n2(i,j+1,k)
              n_avg(3)= n3(i,j+1,k)
           endif
           if (ABS((ABS(n_avg(1))+ABS(n_avg(2))+ABS(n_avg(3)))-1d0) > 1d-12) then
              write(*,*)'Normals not normalised'
              call pariserror('Normals not normalised')
           endif
           alpha2=al3dnew(n_avg,c_stag)
           if (ABS(n_avg(2))>1d-12) then
              y_test2 = (alpha2 - (n_avg(1)+n_avg(3))/2d0)/n_avg(2)
              if (y_test2 < 0.5d0) P_gy(i,j+1,k) = sigma*kap(i,j,k)/dy(j)
              y_mod(i,j,k) = dyh(j)*y_test2
              if (y_mod(i,j,k)>dyh(j)) y_mod(i,j,k) = dyh(j)
              if (y_mod(i,j,k)<limit*dyh(j)) y_mod(i,j,k) = limit*dyh(j)
              if (y_mod(i,j,k) /= y_mod(i,j,k)) write(*,'("y_mod NaN. c_st, n, vofs, phases:",4e14.5,2I8)')&
                   c_stag,n_avg(2),cvof(i,j,k),cvof(i,j+1,k),vof_phase(i,j,k),vof_phase(i,j+1,k)
           endif
        endif
        !-------Check z-neighbours
        if (vof_phase(i,j,k+1) == 1) then
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
           if (c0>c_min) then ! use weighted average if liq VOF is not small, otherwise use cav normals
              nr(1)=n1(i,j,k)*c0/c_stag + n1(i,j,k+1)*c1/c_stag
              nr(2)=n2(i,j,k)*c0/c_stag + n2(i,j,k+1)*c1/c_stag
              nr(3)=n3(i,j,k)*c0/c_stag + n3(i,j,k+1)*c1/c_stag
              n_avg(1) = nr(1)/(ABS(nr(1))+ABS(nr(2))+ABS(nr(3)))
              n_avg(2) = nr(2)/(ABS(nr(1))+ABS(nr(2))+ABS(nr(3)))
              n_avg(3) = nr(3)/(ABS(nr(1))+ABS(nr(2))+ABS(nr(3)))
           else
              n_avg(1)= n1(i,j,k+1)
              n_avg(2)= n2(i,j,k+1)
              n_avg(3)= n3(i,j,k+1)
           endif
           if (ABS((ABS(n_avg(1))+ABS(n_avg(2))+ABS(n_avg(3)))-1d0) > 1d-12) then
              write(*,*)'Normals not normalised'
              call pariserror('Normals not normalised')
           endif
           alpha2=al3dnew(n_avg,c_stag)
           if (ABS(n_avg(3))>1d-12) then
              z_test2 = (alpha2 - (n_avg(1)+n_avg(2))/2d0)/n_avg(3)
              if (z_test2 < 0.5d0) P_gz(i,j,k+1) = sigma*kap(i,j,k)/dz(k)
              z_mod(i,j,k) = dzh(k)*z_test2
              if (z_mod(i,j,k)>dzh(k)) z_mod(i,j,k) = dzh(k)
              if (z_mod(i,j,k)<limit*dzh(k)) z_mod(i,j,k) = limit*dzh(k)
              if (z_mod(i,j,k) /= z_mod(i,j,k)) write(*,'("z_mod NaN. c_st, n, vofs, phases:",4e14.5,2I8)')&
                   c_stag,n_avg(3),cvof(i,j,k),cvof(i,j,k+1),vof_phase(i,j,k),vof_phase(i,j,k+1)
           endif
        endif
     endif
  enddo;enddo;enddo
  call ghost_x(P_gx,1,req(1:4)); call ghost_y(P_gy,1,req(5:8)); call ghost_z(P_gz,1,req(9:12)) 
  call ghost_x(x_mod,1,req(13:16)); call ghost_y(y_mod,1,req(17:20)); call ghost_z(z_mod,1,req(21:24)) 
  call MPI_WAITALL(24,req(1:24),sta(:,1:24),ierr)
!!$  call ghost_y(P_gx,1,req(1:4)); call ghost_y(P_gy,1,req(5:8)); call ghost_y(P_gz,1,req(9:12)) 
!!$  call ghost_y(x_mod,1,req(13:16)); call ghost_y(y_mod,1,req(17:20)); call ghost_y(z_mod,1,req(21:24)) 
!!$  call MPI_WAITALL(24,req(1:24),sta(:,1:24),ierr)
!!$  call ghost_z(P_gx,1,req(1:4)); call ghost_z(P_gy,1,req(5:8)); call ghost_z(P_gz,1,req(9:12)) 
!!$  call ghost_z(x_mod,1,req(13:16)); call ghost_z(y_mod,1,req(17:20)); call ghost_z(z_mod,1,req(21:24)) 
!!$  call MPI_WAITALL(24,req(1:24),sta(:,1:24),ierr)
!--------------------------------------------------------------------------------------------------------
  do k=ks,ke; do j=js,je; do i=is,ie
     if (vof_phase(i,j,k)==0) then
        A(i,j,k,1) = 2d0*dt/(dx(i)*x_mod(i-1,j,k)*(rhot(i-1,j,k)+rhot(i,j,k)))
        if (A(i,j,k,1) /= A(i,j,k,1)) write(*,'("ERROR: A1 NaN :",2e14.5)')A(i,j,k,1),x_mod(i-1,j,k) 
        A(i,j,k,2) = 2d0*dt/(dx(i)*x_mod(i,j,k)*(rhot(i+1,j,k)+rhot(i,j,k)))
        if (A(i,j,k,2) /= A(i,j,k,2)) write(*,'("ERROR: A2 NaN :",2e14.5)')A(i,j,k,2),x_mod(i,j,k)
        A(i,j,k,3) = 2d0*dt/(dy(j)*y_mod(i,j-1,k)*(rhot(i,j-1,k)+rhot(i,j,k)))
        if (A(i,j,k,3) /= A(i,j,k,3)) write(*,'("ERROR: A3 NaN :",2e14.5)')A(i,j,k,3),y_mod(i,j-1,k)
        A(i,j,k,4) = 2d0*dt/(dy(j)*y_mod(i,j,k)*(rhot(i,j+1,k)+rhot(i,j,k)))
        if (A(i,j,k,4) /= A(i,j,k,4)) write(*,'("ERROR: A4 NaN :",2e14.5)')A(i,j,k,4),y_mod(i,j,k)
        A(i,j,k,5) = 2d0*dt/(dz(k)*z_mod(i,j,k-1)*(rhot(i,j,k-1)+rhot(i,j,k)))
        if (A(i,j,k,5) /= A(i,j,k,5)) write(*,'("ERROR: A5 NaN :",2e14.5)')A(i,j,k,5),z_mod(i,j,k-1)
        A(i,j,k,6) = 2d0*dt/(dz(k)*z_mod(i,j,k)*(rhot(i,j,k+1)+rhot(i,j,k)))
        if (A(i,j,k,6) /= A(i,j,k,6)) write(*,'("ERROR: A6 NaN :",2e14.5)')A(i,j,k,6),z_mod(i,j,k)
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
        if (pcmask(i,j,k).ne.0) write(*,'("Error topology, phase 0, pcmask :",i8)')pcmask(i,j,k)
     endif
  enddo;enddo;enddo

  P_bc = 0d0
  ! dp/dn = 0 for inflow bc on face 1 == x- : do not correct u(is-1)
  ! inflow bc on other faces not implemented yet.  !@@ generalize this ! 
  if(coords(1)==0) then
     if(bdry_cond(1)==3) then  
        A(is,:,:,7) = A(is,:,:,7) - A(is,:,:,1)
        A(is,:,:,1) = 0d0
        ! pressure boundary condition
     else if(bdry_cond(1)==5 .and. (.not.RP_test)) then 
        A(is,:,:,8) = (2d0/3d0)*BoundaryPressure(1)  ! P_0 =  1/3 (Pinner - P_b) + P_b
        A(is,:,:,7) = 1d0                      ! P_0  - 1/3 Pinner =  2/3 P_b
        A(is,:,:,1:6) = 0d0                    ! A7 P_is + A2 P_is+1 = A8 
        A(is,:,:,2) = 1d0/3d0 !sign due to definition in Poisson solver
        P_bc(1) = 1d0
     endif
  endif
  ! dp/dn = 0 for outflow/fixed velocity bc on face 4 == x+
  ! outflow/fixed velocity bc on other faces not implemented yet.  
  if(coords(1)==Npx-1) then
     if(bdry_cond(4)==4) then
        A(ie,:,:,7) = A(ie,:,:,7) - A(ie,:,:,2)
        A(ie,:,:,2) = 0d0
        ! pressure boundary condition
     else if(bdry_cond(4)==5 .and. (.not.RP_test)) then
        A(ie,:,:,8) = (2d0/3d0)*BoundaryPressure(2)
        A(ie,:,:,7) = 1d0  
        A(ie,:,:,2:6) = 0d0
        A(ie,:,:,1) = 1d0/3d0
        P_bc(2) = 1d0
     endif
  endif
 
  if (.not. RP_test) then
     ! Pressure BC for y-
     if(coords(2)==0 .and. (bdry_cond(2)==5)) then
        A(:,js,:,8) = (2d0/3d0)*BoundaryPressure(3)
        A(is,js,:,8) = (2d0/3d0)*(BoundaryPressure(3)+BoundaryPressure(1)*P_bc(1))
        A(ie,js,:,8) = (2d0/3d0)*(BoundaryPressure(3)+BoundaryPressure(2)*P_bc(2))
        A(:,js,:,7) = 1d0
        A(is,js,:,7) = 1d0 + P_bc(1)
        A(ie,js,:,7) = 1d0 + P_bc(2)
        A(:,js,:,1:6) = 0d0      
        A(is,js,:,2) = 1d0/3d0*P_bc(1)
        A(ie,js,:,1) = 1d0/3d0*P_bc(2)
        A(:,js,:,4) = 1d0/3d0
        P_bc(3) = 1d0
     endif
     ! Pressure BC for y+
     if(coords(2)==Npy-1 .and. (bdry_cond(5)==5) ) then
        A(:,je,:,8) = (2d0/3d0)*BoundaryPressure(4)
        A(is,je,:,8) = (2d0/3d0)*(BoundaryPressure(4)+BoundaryPressure(1)*P_bc(1))
        A(ie,je,:,8) = (2d0/3d0)*(BoundaryPressure(4)+BoundaryPressure(2)*P_bc(2))
        A(:,je,:,7) = 1d0  
        A(is,je,:,7) = 1d0 + P_bc(1)
        A(ie,je,:,7) = 1d0 + P_bc(2)
        A(:,je,:,1:6) = 0d0
        A(is,je,:,2) = 1d0/3d0*P_bc(1)
        A(ie,je,:,1) = 1d0/3d0*P_bc(2)
        A(:,je,:,3) = 1d0/3d0
        P_bc(4) = 1d0
     endif

     ! Pressure BC for z-
     if(coords(3)==0 .and. (bdry_cond(3)==5)) then
        A(:,:,ks,8) = (2d0/3d0)*BoundaryPressure(5)
        A(is,js,ks,8) = (2d0/3d0)*(BoundaryPressure(5)+BoundaryPressure(1)*P_bc(1)+BoundaryPressure(3)*P_bc(3))
        A(ie,js,ks,8) = (2d0/3d0)*(BoundaryPressure(5)+BoundaryPressure(2)*P_bc(2)+BoundaryPressure(3)*P_bc(3))
        A(is,je,ks,8) = (2d0/3d0)*(BoundaryPressure(5)+BoundaryPressure(1)*P_bc(1)+BoundaryPressure(4)*P_bc(4))
        A(ie,je,ks,8) = (2d0/3d0)*(BoundaryPressure(5)+BoundaryPressure(2)*P_bc(2)+BoundaryPressure(4)*P_bc(4))
        A(:,:,ks,7) = 1d0
        A(is,js,ks,7) = 1d0 + P_bc(1) + P_bc(3)
        A(ie,js,ks,7) = 1d0 + P_bc(2) + P_bc(3)
        A(is,je,ks,7) = 1d0 + P_bc(1) + P_bc(4)
        A(ie,je,ks,7) = 1d0 + P_bc(2) + P_bc(4)
        A(:,:,ks,1:6) = 0d0   
        A(is,js,ks,2) = 1d0/3d0*P_bc(1); A(is,js,ks,4) = 1d0/3d0*P_bc(3)
        A(ie,js,ks,1) = 1d0/3d0*P_bc(2); A(ie,js,ks,4) = 1d0/3d0*P_bc(3)
        A(is,je,ks,2) = 1d0/3d0*P_bc(1); A(is,je,ks,3) = 1d0/3d0*P_bc(4)
        A(ie,je,ks,1) = 1d0/3d0*P_bc(2); A(ie,je,ks,3) = 1d0/3d0*P_bc(4)
        A(:,:,ks,6) = 1d0/3d0
     endif
     ! Pressure BC for z+
     if(coords(3)==Npz-1 .and. (bdry_cond(6)==5) ) then
        A(:,:,ke,8) = (2d0/3d0)*BoundaryPressure(6)
        A(is,js,ke,8) = (2d0/3d0)*(BoundaryPressure(6)+BoundaryPressure(1)*P_bc(1)+BoundaryPressure(3)*P_bc(3))
        A(ie,js,ke,8) = (2d0/3d0)*(BoundaryPressure(6)+BoundaryPressure(2)*P_bc(2)+BoundaryPressure(3)*P_bc(3))
        A(is,je,ke,8) = (2d0/3d0)*(BoundaryPressure(6)+BoundaryPressure(1)*P_bc(1)+BoundaryPressure(4)*P_bc(4))
        A(ie,je,ke,8) = (2d0/3d0)*(BoundaryPressure(6)+BoundaryPressure(2)*P_bc(2)+BoundaryPressure(4)*P_bc(4))
        A(:,:,ke,7) = 1d0  
        A(is,js,ke,7) = 1d0 + P_bc(1) + P_bc(3)
        A(ie,js,ke,7) = 1d0 + P_bc(2) + P_bc(3)
        A(is,je,ke,7) = 1d0 + P_bc(1) + P_bc(4)
        A(ie,je,ke,7) = 1d0 + P_bc(2) + P_bc(4)
        A(:,:,ke,1:6) = 0d0
        A(is,js,ke,2) = 1d0/3d0*P_bc(1); A(is,js,ks,4) = 1d0/3d0*P_bc(3)
        A(ie,js,ke,1) = 1d0/3d0*P_bc(2); A(ie,js,ks,4) = 1d0/3d0*P_bc(3)
        A(is,je,ke,2) = 1d0/3d0*P_bc(1); A(is,je,ks,3) = 1d0/3d0*P_bc(4)
        A(ie,je,ke,1) = 1d0/3d0*P_bc(2); A(ie,je,ks,3) = 1d0/3d0*P_bc(4)
        A(:,:,ke,5) = 1d0/3d0
     endif
  endif
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
  integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: vof_phase
  real(8) :: dt,h_mod
  integer :: i,j,k,nbr
  integer :: istep

  do k=ks,ke; do j=js,je; do i=is,ie
     if (pcmask(i,j,k)==1 .or. pcmask(i,j,k)==2) then !rho is 1d0
        A(i,j,k,1) = 1d0/(dx(i)*dxh(i-1))
        if (A(i,j,k,1) /= A(i,j,k,1)) write(*,'("ERROR: A1 NaN in fs2:",e14.5)')A(i,j,k,1)
        A(i,j,k,2) = 1d0/(dx(i)*dxh(i))
        if (A(i,j,k,2) /= A(i,j,k,2)) write(*,'("ERROR: A2 NaN in fs2:",e14.5)')A(i,j,k,2)
        A(i,j,k,3) = 1d0/(dy(j)*dyh(j-1))
        if (A(i,j,k,3) /= A(i,j,k,3)) write(*,'("ERROR: A3 NaN in fs2:",e14.5)')A(i,j,k,3)
        A(i,j,k,4) = 1d0/(dy(j)*dyh(j))
        if (A(i,j,k,4) /= A(i,j,k,4)) write(*,'("ERROR: A4 NaN in fs2:",e14.5)')A(i,j,k,4)
        A(i,j,k,5) = 1d0/(dz(k)*dzh(k-1))
        if (A(i,j,k,5) /= A(i,j,k,5)) write(*,'("ERROR: A5 NaN in fs2:",e14.5)')A(i,j,k,5)
        A(i,j,k,6) = 1d0/(dz(k)*dzh(k))
        if (A(i,j,k,6) /= A(i,j,k,6)) write(*,'("ERROR: A6 NaN in fs2:",e14.5)')A(i,j,k,6)

        if (pcmask(i,j,k)==2 .and. vof_phase(i,j,k)/=1) call pariserror("Fatal topology error in FS 2nd projection")

        if (pcmask(i,j,k)==1) then
           if (vof_phase(i,j,k)/=1) call pariserror("Fatal topology error in FS 2nd projection")
           do nbr=-1,1,2
              if (pcmask(i+nbr,j,k) == 0) then
                 A(i,j,k,2+(nbr-1)/2) = 0d0
              endif
              if (pcmask(i,j+nbr,k) == 0) then
                 A(i,j,k,4+(nbr-1)/2) = 0d0
              endif
              if (pcmask(i,j,k+nbr) == 0) then
                 A(i,j,k,6+(nbr-1)/2) = 0d0
              endif
           enddo
        endif
        A(i,j,k,7) = sum(A(i,j,k,1:6))
        A(i,j,k,8) =  -1d0*((utmp(i,j,k)-utmp(i-1,j,k))/dx(i) &
             +  (vtmp(i,j,k)-vtmp(i,j-1,k))/dy(j) &
             +  (wtmp(i,j,k)-wtmp(i,j,k-1))/dz(k))
        if (A(i,j,k,8) /= A(i,j,k,8)) write(*,'("ERROR: A8 NaN in fs2:",e14.5)')A(i,j,k,8) 
     endif
  enddo; enddo; enddo
end subroutine setuppoisson_fs2
!--------------------------------------------------------------------------------------------------------------------
!=================================================================================================
!
!  Extrapolation of velocities for free surface
!
!=================================================================================================
subroutine extrapolate_velocities()
  use module_BC
  use module_grid
  use module_flow
  use module_2phase
  use module_freesurface
  implicit none
  integer :: i,j,k,level,nbr
  real(8) :: nr(3)
  real(8) :: x_vel, xcount

  !New simple extrapolation: Velocities of neighbours at lower topological level averaged
  do level = 1, X_level
     do k=ks,ke; do j=js,je; do i=is,ie
        if (u_cmask(i,j,k) == level) then
           xcount = 0d0; x_vel = 0d0
           do nbr=-1,1,2
              if (u_cmask(i+nbr,j,k)==level-1) then
                 xcount = xcount+1d0
                 x_vel = x_vel + u(i+nbr,j,k)
              endif
           enddo
           do nbr=-1,1,2
              if (u_cmask(i,j+nbr,k)==level-1) then
                 xcount = xcount+1d0
                 x_vel = x_vel + u(i,j+nbr,k)
              endif
           enddo
           do nbr=-1,1,2
              if (u_cmask(i,j,k+nbr)==level-1) then
                 xcount = xcount+1d0
                 x_vel = x_vel + u(i,j,k+nbr)
              endif
           enddo
           if (xcount>0d0) then
              u(i,j,k) = x_vel/xcount
           endif
        endif
        if (v_cmask(i,j,k) == level) then
           xcount = 0d0; x_vel = 0d0
           do nbr = -1,1,2
              if (v_cmask(i+nbr,j,k)==level-1) then
                 xcount = xcount+1d0
                 x_vel = x_vel + v(i+nbr,j,k)
              endif
           enddo
           do nbr = -1,1,2
              if (v_cmask(i,j+nbr,k)==level-1) then
                 xcount = xcount+1d0
                 x_vel = x_vel + v(i,j+nbr,k)
              endif
           enddo
           do nbr = -1,1,2
              if (v_cmask(i,j,k+nbr)==level-1) then
                 xcount = xcount+1d0
                 x_vel = x_vel + v(i,j,k+nbr)
              endif
           enddo
           if (xcount>0d0) then
              v(i,j,k) = x_vel/xcount
           endif
        endif
        if (w_cmask(i,j,k) == level) then
           xcount = 0d0; x_vel = 0d0
           do nbr = -1,1,2
              if (w_cmask(i+nbr,j,k)==level-1) then
                 xcount = xcount+1d0
                 x_vel = x_vel + w(i+nbr,j,k)
              endif
           enddo
           do nbr = -1,1,2
              if (w_cmask(i,j+nbr,k)==level-1) then
                 xcount = xcount+1d0
                 x_vel = x_vel + w(i,j+nbr,k)
              endif
           enddo
           do nbr = -1,1,2
              if (w_cmask(i,j,k+nbr)==level-1) then
                 xcount = xcount+1d0
                 x_vel = x_vel + w(i,j,k+nbr)
              endif
           enddo
           if (xcount>0d0) then
              w(i,j,k) = x_vel/xcount
           endif
        endif
     enddo; enddo; enddo
     call do_ghost_vector(u,v,w)
  enddo
end subroutine extrapolate_velocities
!=============================================================================================================================================
subroutine discrete_divergence(u,v,w,iout)
  use module_grid
  use module_freesurface
  use module_IO
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: u,v,w 
  real(8), dimension(is:ie+1,js:je+1,ks:ke+1) :: div
  real(8), dimension(0:3) :: divtot, domain, n_level, n_total
  real(8) :: avg2
  integer :: i,j,k,l,iout,ierr,prank
  character(len=30) :: filename

  divtot = 0d0; n_level = 0d0

  do k=ks,ke+1; do j=js,je+1; do i=is,ie+1
     div(i,j,k)=ABS((u(i-1,j,k)-u(i,j,k))/dx(i)+(v(i,j-1,k)-v(i,j,k))/dy(j)+(w(i,j,k-1)-w(i,j,k))/dz(k)) !divergence
     !div(i,j,k)=ABS((u(i-1,j,k)-u(i,j,k))*dy(j)*dz(k)+(v(i,j-1,k)-v(i,j,k))*dx(i)*dz(k)+(w(i,j,k-1)-w(i,j,k))*dx(i)*dy(j)) !volume error
     !if (pcmask(i,j,k) == 1 .or. pcmask(i,j,k) == 2) then
     !   div(i,j,k) = div(i,j,k)*dt
     !endif
     do l=0,3
        if (pcmask(i,j,k)==l) then
           divtot(l)=divtot(l)+div(i,j,k)
           n_level(l) = n_level(l) + 1d0
        endif
     enddo
  enddo; enddo; enddo
  call MPI_ALLREDUCE(divtot,domain,4, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_Cart, ierr)
  call MPI_ALLREDUCE(n_level,n_total,4, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_Cart, ierr)
  avg2 = (domain(1)+domain(2))/(n_total(1)+n_total(2))
  do l=0,3
     if (n_total(l) > 1d-10) domain(l)=domain(l)/n_total(l)
  enddo
  if (rank==0) then
     OPEN(unit=21,file='div_type.txt',position='append')
     write(21,115)iout,domain(0),domain(1),domain(2),domain(3),avg2
     close(unit=21)
  endif
115 format(I8,5e14.5)
end subroutine discrete_divergence
!--------------------------------------------------------------------------------------------------------------------
subroutine get_ref_volume
  use module_2phase
  use module_BC
  use module_freesurface
  implicit none
  real(8), parameter :: pi=3.1415926535897932384626433833

  if (NumBubble /= 1) call pariserror('For the Rayleigh-Plesset test, only a single bubble is allowed')
  V_0 = 4d0/3d0*pi*(R_ref**3d0)
  R_RK = rad(1)
  dR_RK = 0d0 
  P_inf = BoundaryPressure(1)
  ddR_RK = (P_ref*(R_ref/R_RK)**(3d0*gamma)-P_inf)/R_RK
  !if (rank==0) then
  !   write(*,'("RP test initial bubble volume:",e14.5)')V_0
  !   write(*,'("Free surface ddR: ",e14.5)')ddR_RK
  !endif
end subroutine get_ref_volume
!--------------------------------------------------------------------------------------------------------------------
! This is a straight copy of NewSolver. The idea is to not clutter NewSolver with all the FreeSurface flags and tests, 
! therefore it was copied here and named FreeSolver.
!Solves the following linear equiation:
! A7*Pijk = A1*Pi-1jk + A2*Pi+1jk + A3*Pij-1k + 
!           A4*Pij+1k + A5*Pijk-1 + A6*Pijk+1 + A8
!-------------------------------------------------------------------------------------------------
subroutine FreeSolver(A,p,maxError,beta,maxit,it,ierr,iout,time,tres2)
  use module_grid
  use module_BC
  use module_IO
  use module_freesurface
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: p
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(in) :: A
  real(8), intent(in) :: beta, maxError
  integer, intent(in) :: maxit
  integer, intent(out) :: it, ierr
  real(8) :: res1,res2,resinf,intvol
  real(8) :: tres2, cells, p_c, Vol, time, tcells
  integer :: i,j,k, iout
  integer :: req(12),sta(MPI_STATUS_SIZE,12)
  logical :: mask(imin:imax,jmin:jmax,kmin:kmax)
  integer, parameter :: norm=2, relaxtype=1
  logical, parameter :: recordconvergence=.false.
  integer, save :: itime=0
! Open file for convergence history
  if(rank==0.and.recordconvergence) then
     OPEN(UNIT=89,FILE=TRIM(out_path)//'/convergence_history-'//TRIM(int2text(itime,padding))//'.txt')
  endif
  itime=itime+1
  if (RP_test) then
     call get_vol(Vol)
     p_c = P_ref*(V_0/Vol)**gamma
     if(mod(iout,10)==0 .and. rank==0 .and. solver_flag==2) then
        OPEN(unit=11,file='RP_volume',position='append')
        write(11,2)time,Vol
     endif
  else
     p_c = 0d0
  endif
2 format(2e14.5)
  if (solver_flag == 0) call pariserror("Free Surface solver flag needs to be 1 or 2")
  do k=ks,ke; do j=js,je; do i=is,ie
     if (solver_flag == 1 .and. pcmask(i,j,k) /= 0) p(i,j,k) = p_c
     if (solver_flag == 2 .and. pcmask(i,j,k)==3) p(i,j,k) = 0d0
  enddo; enddo; enddo
  !--------------------------------------ITERATION LOOP--------------------------------------------  
  do it=1,maxit
     if(relaxtype==2) then 
        call LineRelax_fs(A,p,beta)
     elseif(relaxtype==1) then
        call RedBlackRelax_fs(A,p,beta)
     endif
!---------------------------------CHECK FOR CONVERGENCE-------------------------------------------
    res1 = 0d0; res2=0.d0; resinf=0.d0; intvol=0.d0; cells = 0d0
    call ghost_x(p,1,req( 1: 4)); call ghost_y(p,1,req( 5: 8)); call ghost_z(p,1,req( 9:12))
    do k=ks+1,ke-1; do j=js+1,je-1; do i=is+1,ie-1
       if ((pcmask(i,j,k)==0 .and. solver_flag==1)&
            .or.((pcmask(i,j,k)==1 .or. pcmask(i,j,k)==2) .and. solver_flag==2)) then
          res2=res2+abs(-p(i,j,k) * A(i,j,k,7) +                           &
               A(i,j,k,1) * p(i-1,j,k) + A(i,j,k,2) * p(i+1,j,k) +            &
               A(i,j,k,3) * p(i,j-1,k) + A(i,j,k,4) * p(i,j+1,k) +            &
               A(i,j,k,5) * p(i,j,k-1) + A(i,j,k,6) * p(i,j,k+1) + A(i,j,k,8) )**norm
          cells = cells + 1d0
       endif
    enddo; enddo; enddo
    call MPI_WAITALL(12,req,sta,ierr)
    mask=.true.
    mask(is+1:ie-1,js+1:je-1,ks+1:ke-1)=.false.
    do k=ks,ke; do j=js,je; do i=is,ie
       if(mask(i,j,k) .and. ((pcmask(i,j,k)==0 .and. solver_flag==1) .or. &
            ((pcmask(i,j,k)==1 .or. pcmask(i,j,k)==2) .and. solver_flag==2))) then
          res2=res2+abs(-p(i,j,k) * A(i,j,k,7) +&
               A(i,j,k,1) * p(i-1,j,k) + A(i,j,k,2) * p(i+1,j,k) +            &
               A(i,j,k,3) * p(i,j-1,k) + A(i,j,k,4) * p(i,j+1,k) +            &
               A(i,j,k,5) * p(i,j,k-1) + A(i,j,k,6) * p(i,j,k+1) + A(i,j,k,8) )**norm
          cells = cells + 1d0
       endif
    enddo; enddo; enddo
    call catch_divergence_fs(res2,cells,ierr)
    call MPI_ALLREDUCE(res2, tres2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr)
    call MPI_ALLREDUCE(cells, tcells, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr)
    if(norm==2) tres2=sqrt(tres2)
    tres2 = tres2/tcells
    if(rank==0.and.mod(it,10) == 0.and.recordconvergence) write(89,310) it, solver_flag, tres2
310 format(2I6,'  ',(e14.5))
    if (tres2<maxError) then 
       if(rank==0.and.recordconvergence) close(89)
       exit
    endif
  enddo
  if(rank==0.and.recordconvergence) close(89)
  if(it==maxit+1 .and. rank==0) then
     write(*,*) 'Warning: LinearSolver reached maxit: ||res||: ',tres2
     write(*,'("Solver flag:",I8)')solver_flag
  endif
contains
  subroutine catch_divergence_fs(res2,cells,ierr)
    real(8), intent(in) :: res2, cells
    integer, intent(out) :: ierr
    logical :: diverged=.false.
    logical :: extended=.true.
    integer :: l
    if(extended) then
       do k=ks,ke; do j=js,je; do i=is,ie
          do l=1,8
             if ((pcmask(i,j,k)==0 .and. solver_flag==1) &
                  .or.((pcmask(i,j,k)==1 .or. pcmask(i,j,k)==2) .and. solver_flag==2)) then
                if(A(i,j,k,l)/=A(i,j,k,l).or.p(i,j,k)/=p(i,j,k)) then
                   diverged=.true.
                endif
             endif
          enddo
          if(diverged) then
             OPEN(UNIT=88,FILE=TRIM(out_path)//'/message-rank-'//TRIM(int2text(rank,padding))//'.txt')
             write(88,*) "ijk rank",i,j,k,rank
             write(88,*) "A",  A(i,j,k,:)
             write(88,*) "p",  p(i,j,k)
             write(88,'("P neighbours: ",6e14.5)')p(i-1,j,k),p(i+1,j,k),p(i,j-1,k),p(i,j+1,k),&
                  p(i,j,k-1),p(i,j,k+1)
             write(88,*) 'A or p is NaN after',it,'iterations at rank ',rank
             write(88,'("Res2, cells, solver_flag :",2e14.5,I8)')res2,cells,solver_flag
             write(88,'("Pcmask 1-7: ",7I8)')pcmask(i-1,j,k),pcmask(i+1,j,k),pcmask(i,j-1,k),pcmask(i,j+1,k),&
                  pcmask(i,j,k-1),pcmask(i,j,k+1),pcmask(i,j,k)
             close(88)
             if(rank<=30) print*,'A or p is NaN after',it,'iterations at rank ',rank
             call pariserror("A or p is NaN")
             exit
          endif
       end do; end do; end do
    endif
    if ((res2*npx*npy*npz)>1.d16 ) then
       if(rank<=30) then
          print*,'Pressure solver diverged after',it,'iterations at rank ',rank
          write(*,'("Res2, cells, solver_flag :",2e14.5,I8)')res2,cells,solver_flag
       endif
       call pariserror("freesolver error")
    else if (res2 .ne. res2) then 
       if(rank<=30) then
          print*, 'it:',it,'Pressure residual value is invalid at rank', rank
          write(*,'("Res2, cells, solver_flag :",2e14.5,I8)')res2,cells,solver_flag
       endif
       call pariserror("freesolver error")
    else
       ierr=0
    endif
  end subroutine catch_divergence_fs
  subroutine get_vol(Vol_RP)
    !use module_grid
    use module_VOF
    !use module_2phase
    !implicit none
    !include 'mpif.h'
    real(8) :: Vol_RP, bub_vol, vol_tot
    integer :: i,j,k,ierr
    bub_vol = 0d0
    do k=ks,ke; do j=js,je; do i=is,ie
       bub_vol = bub_vol + cvof(i,j,k)*dx(i)*dy(j)*dz(k)
    enddo; enddo; enddo
    call MPI_ALLREDUCE(bub_vol,vol_tot,1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr)
    Vol_RP = vol_tot/(xLength*yLength*zLength)
  end subroutine get_vol
end subroutine FreeSolver
!--------------------------------------ONE RELAXATION ITERATION (SMOOTHER)----------------------------------
subroutine RedBlackRelax_fs(A,p,beta)
  use module_grid
  use module_BC
  use module_IO
  use module_freesurface
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: p
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(in) :: A
  real(8), intent(in) :: beta
  integer :: req(12),sta(MPI_STATUS_SIZE,12)
  integer :: i,j,k,ierr
  integer :: isw,jsw,ksw,ipass
  ksw=1
  do ipass=1,2
     jsw=ksw
     do k=ks,ke
        isw=jsw
        do j=js,je
           do i=isw+is-1,ie,2
              if ((pcmask(i,j,k)==0 .and. solver_flag==1) &
                   .or.((pcmask(i,j,k)==1 .or. pcmask(i,j,k)==2) .and. solver_flag==2)) then
                 p(i,j,k)=(1d0-beta)*p(i,j,k) + (beta/A(i,j,k,7))*(              &
                      A(i,j,k,1) * p(i-1,j,k) + A(i,j,k,2) * p(i+1,j,k) +        &
                      A(i,j,k,3) * p(i,j-1,k) + A(i,j,k,4) * p(i,j+1,k) +        &
                      A(i,j,k,5) * p(i,j,k-1) + A(i,j,k,6) * p(i,j,k+1) + A(i,j,k,8))   
              endif
           enddo
           isw=3-isw
        enddo
        jsw=3-jsw
     enddo
     ksw=3-ksw
     if(ipass==1) then
        call ghost_x(p,1,req( 1: 4)); call ghost_y(p,1,req( 5: 8)); call ghost_z(p,1,req( 9:12))
        call MPI_WAITALL(12,req,sta,ierr)
     endif
  enddo
end subroutine RedBlackRelax_fs
!--------------------------------------ONE RELAXATION ITERATION (SMOOTHER)----------------------------------
subroutine LineRelax_fs(A,p,beta)
  use module_grid
  use module_freesurface
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: p
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(in) :: A
  real(8), intent(in) :: beta
  integer :: i,j,k
!--------------------------------------ITERATION LOOP--------------------------------------------  
  do k=ks,ke; do j=js,je; do i=is,ie
     if ((pcmask(i,j,k)==0 .and. solver_flag==1).or.((pcmask(i,j,k)==1 .or. pcmask(i,j,k)==2) .and. solver_flag==2)) then
        p(i,j,k)=(1d0-beta)*p(i,j,k) + (beta/A(i,j,k,7))*(              &
             A(i,j,k,1) * p(i-1,j,k) + A(i,j,k,2) * p(i+1,j,k) +        &
             A(i,j,k,3) * p(i,j-1,k) + A(i,j,k,4) * p(i,j+1,k) +        &
             A(i,j,k,5) * p(i,j,k-1) + A(i,j,k,6) * p(i,j,k+1) + A(i,j,k,8))
     endif
  enddo; enddo; enddo
end subroutine LineRelax_fs
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!------------VTK output routines------------------------------------------------------------------
!
! These routines output data used to debug and/or display variables in the free surface part of the
! code.
!
!-------------------------------------------------------------------------------------------------
subroutine append_visit_fs(index,iout)
  use module_flow
  use module_grid
  use module_IO
  use module_freesurface
  implicit none
  character(len=10) :: file
  character(len=40) :: file_root
  integer index, iout, prank
  if(rank.ne.0) call pariserror("rank.ne.0 in append")
  file = TRIM(visit_file(index))
  if(.not.vtk_open(index)) then
     OPEN(UNIT=100,FILE=TRIM(file)//'.visit')
     write(100,10) nPdomain
10   format('!NBLOCKS ',I4)
     vtk_open(index) = .true.
  else
     OPEN(UNIT=100,FILE=TRIM(file)//'.visit',position='append')
  endif

  file_root = TRIM(out_path)//'/VTK/'//TRIM(file_short(index))
  do prank=0,NpDomain-1
     write(100,11) TRIM(file_root)//TRIM(int2text(iout,padding))//'-'//TRIM(int2text(prank,padding))//'.vtk'
11   format(A)
  enddo
  close(100)
end subroutine  append_visit_fs
!-------------------------------------------------------------------------------------------------
subroutine VTK_scalar_struct(index,iout,var)
  use module_flow
  use module_grid
  use module_IO
  use module_freesurface
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax) :: var 
  character(len=40) :: file_root
  integer index, iout, i,j,k
  file_root = TRIM(out_path)//'/VTK/'//TRIM(file_short(index))//TRIM(int2text(iout,padding))//'-'
  !Write to VTK file
  OPEN(UNIT=8,FILE=TRIM(file_root)//TRIM(int2text(rank,padding))//'.vtk')
  write(8,10)
  write(8,11) time
  write(8,12)
  write(8,13)
  write(8,14)ie+1-is+1,je+1-js+1,ke+1-ks+1
  write(8,15) x(is),y(js),z(ks)
  write(8,16) x(is+1)-x(is),y(js+1)-y(js),z(ks+1)-z(ks)
10 format('# vtk DataFile Version 2.0')
11 format('grid, time ',F16.8)
12 format('ASCII')
13 format('DATASET STRUCTURED_POINTS')
14 format('DIMENSIONS ',I5,I5,I5)
15 format('ORIGIN ',F16.8,F16.8,F16.8)
16 format('SPACING ',F16.8,F16.8,F16.8)

  write(8,19)(ie+1-is+1)*(je+1-js+1)*(ke+1-ks+1)
  write(8,17)file_short(index)
  write(8,18)
19 format('POINT_DATA ',I17)
17 format('SCALARS ',A20,' float 1')
18 format('LOOKUP_TABLE default')

  do k=ks,ke+1; do j=js,je+1; do i=is,ie+1
     write(8,210) var(i,j,k)
  enddo; enddo; enddo
210 format(e14.5)
  close(8)
end subroutine VTK_scalar_struct
!==============================================================================================================================
! Routine to numerically integrate the Rayleigh-Plesset equation
subroutine Integrate_RP(dt,t)
  !use module_2phase
  use module_freesurface
  implicit none 
  integer :: j
  real(8) :: dt, t, Volume
  real(8), parameter :: pi=3.141592653589793238462643383
  integer, parameter :: nvar = 2
  real(8) :: y(nvar), ytmp(nvar)
  real(8) :: dydt1(nvar), dydt2(nvar),dydt3(nvar),dydt0(nvar),dydt4(nvar)

    ! start at previous RP values
    y(1) = R_RK
    y(2) = dR_RK
    ! RK4
    !1st step
    ytmp(:)  = y(:)
    call func(ytmp,dydt0)
    !2nd step
    do j=1,nvar
       ytmp(j) = y(j) + dydt0(j)*dt/2.d0
    enddo
    call func(ytmp,dydt1)
    !3rd step
    do j=1,nvar
       ytmp(j) = y(j) + dydt1(j)*dt/2.d0
    enddo
    call func(ytmp,dydt2)
    !4th step
    do j=1,nvar
       ytmp(j) = y(j) + dydt2(j)*dt
    enddo
    call func(ytmp,dydt3)

    do j=1,nvar
       y(j) = y(j) + dt*(dydt0(j)/6d0 + dydt1(j)/3d0 + dydt2(j)/3d0 + dydt3(j)/6d0)
    enddo
    Volume = 4d0/3d0*pi*y(1)**3d0
    call func(y,dydt4)
    R_RK = y(1); dR_RK = y(2); ddR_RK = dydt4(2)
    
contains
subroutine func(y,dydt)
  use module_2phase
  use module_freesurface
    real(8) :: P_c
    real(8) :: y(2), dydt(2)
    P_c = P_ref*(R_ref/y(1))**(3d0*gamma)-2d0*sigma/y(1)

    dydt(2) = -3d0*(y(2)**2d0)/(2d0*y(1)) + (P_c - P_inf)/y(1)
    dydt(1) = y(2)
  end subroutine func
end subroutine Integrate_RP
!=======================================================================================================================================
subroutine write_RP_test(t) 
  use module_freesurface
  implicit none
  real(8) :: t, vol, p_mid, p_corner
  real(8), parameter :: pi=3.141592653589793238462643383
  vol = 4d0/3d0*pi*R_RK**3d0
  p_mid = pressure(0.5d0)
  p_corner = pressure(sqrt(2d0)/2d0)
  OPEN(UNIT=20,FILE='RK_int_RP.txt',position='append')
  WRITE(20,2) t, R_RK, dR_RK, ddR_RK, vol, p_mid, p_corner
  CLOSE(unit=20)
2 format(7e14.5)
contains
  function pressure(r)
    use module_2phase
    use module_freesurface
    real(8) :: r, P_l, pressure
    P_l = P_ref*(R_ref/R_RK)**(3d0*gamma)-2d0*sigma/R_RK
    pressure = P_l - (dR_RK**2d0 * R_RK**4d0/(2d0*r**4d0) - (ddR_RK*R_RK**2d0 + 2d0*R_RK*dR_RK**2d0)/r +&
         ddR_RK*R_RK + 3d0/2d0*dR_RK**2d0)
  end function pressure
end subroutine write_RP_test
!======================================================================================================================================
! Idea is to set free outflow in the radial direction for the RP test case. In combination with variable Pressure bx's, we can have more
! accurate results on a truncated domain.
!
subroutine set_RP_radial_velocity(u,v,w)
  use module_grid
  use module_2phase
  use module_freesurface
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: u, v, w
  real(8) :: di1, di2
  integer :: ih1, ih2
  integer :: i,j,k

  !x- face
  if (coords(1)==0) then
     !write(*,*)'x- face loop'
     do j=js,je; do k=ks,ke
        ! Components to be done independently, starting xith u on the boundary
        if (y(j).ge.yc(1)) then 
           ih1 = -1
        else
           ih1 = 1
        endif
        if (z(k).ge.zc(1)) then 
           ih2 = -1
        else
           ih2 = 1
        endif
        !get cuts, will be non-dimensional if dx=dy=dz
        di1 = ABS(yc(1)-y(j))/ABS(xc(1)-xh(is-1))
        di2 = ABS(zc(1)-z(k))/ABS(xc(1)-xh(is-1))
        !set
        u(is-1,j,k)=u(is,j,k)*(1d0-di1)*(1d0-di2)+u(is,j+ih1,k)*di1*(1d0-di2)&
             +u(is,j,k+ih2)*(1d0-di1)*di2+u(is,j+ih1,k+ih2)*di1*di2
        !same for v
        if (yh(j).ge.yc(1)) then 
           ih1 = -1
        else
           ih1 = 1
        endif
        if (zh(k).ge.zc(1)) then 
           ih2 = -1
        else
           ih2 = 1
        endif
        !get cuts, will be non-dimensional if dx=dy=dz
        di1 = ABS(yc(1)-yh(j))/ABS(xc(1)-x(is-1))
        di2 = ABS(zc(1)-zh(k))/ABS(xc(1)-x(is-1))
        !set
        v(is-1,j,k)=v(is,j,k)*(1d0-di1)*(1d0-di2)+v(is,j+ih1,k)*di1*(1d0-di2)&
             +v(is,j,k+ih2)*(1d0-di1)*di2+v(is,j+ih1,k+ih2)*di1*di2
        !same for w
        w(is-1,j,k)=w(is,j,k)*(1d0-di1)*(1d0-di2)+w(is,j+ih1,k)*di1*(1d0-di2)&
             +w(is,j,k+ih2)*(1d0-di1)*di2+w(is,j+ih1,k+ih2)*di1*di2
        !set u(is-2) to allow cell is-1 to be divergence free
     enddo; enddo
  endif
  !x+ face
  if (coords(1)==nPx-1) then
    !write(*,*)'x+ face loop' 
     do j=js,je; do k=ks,ke
        ! Components to be done independently, starting xith u on the boundary
        if (y(j).ge.yc(1)) then 
           ih1 = -1
        else
           ih1 = 1
        endif
        if (z(k).ge.zc(1)) then 
           ih2 = -1
        else
           ih2 = 1
        endif
        !get cuts, will be non-dimensional if dx=dy=dz
        di1 = ABS(yc(1)-y(j))/ABS(xc(1)-xh(ie))
        di2 = ABS(zc(1)-z(k))/ABS(xc(1)-xh(ie))
        !set
        u(ie,j,k)=u(ie-1,j,k)*(1d0-di1)*(1d0-di2)+u(ie-1,j+ih1,k)*di1*(1d0-di2)&
             +u(ie-1,j,k+ih2)*(1d0-di1)*di2+u(ie-1,j+ih1,k+ih2)*di1*di2
        !same for v
        if (yh(j).ge.yc(1)) then 
           ih1 = -1
        else
           ih1 = 1
        endif
        if (zh(k).ge.zc(1)) then 
           ih2 = -1
        else
           ih2 = 1
        endif
        !get cuts, will be non-dimensional if dx=dy=dz
        di1 = ABS(yc(1)-yh(j))/ABS(xc(1)-x(ie+1))
        di2 = ABS(zc(1)-zh(k))/ABS(xc(1)-x(ie+1))
        !set
        v(ie,j,k)=v(ie-1,j,k)*(1d0-di1)*(1d0-di2)+v(ie-1,j+ih1,k)*di1*(1d0-di2)&
             +v(ie-1,j,k+ih2)*(1d0-di1)*di2+v(ie-1,j+ih1,k+ih2)*di1*di2
        !same for w
        w(ie,j,k)=w(ie-1,j,k)*(1d0-di1)*(1d0-di2)+w(ie-1,j+ih1,k)*di1*(1d0-di2)&
             +w(ie-1,j,k+ih2)*(1d0-di1)*di2+w(ie-1,j+ih1,k+ih2)*di1*di2
        !set u(ie+2) to allow cell to be divergence free
     enddo; enddo
  endif

  !y-face
  if (coords(2)==0) then
     !write(*,*)'y- face loop'
     do i=is,ie; do k=ks,ke
        if (x(i).ge.xc(1)) then 
           ih1 = -1
        else
           ih1 = 1
        endif
        if (z(k).ge.zc(1)) then 
           ih2 = -1
        else
           ih2 = 1
        endif
        !get cuts, will be non-dimensional if dx=dy=dz
        di1 = ABS(xc(1)-x(i))/ABS(yc(1)-yh(is-1))
        di2 = ABS(zc(1)-z(k))/ABS(yc(1)-yh(is-1))
        !set
        v(i,js-1,k)=v(i,js,k)*(1d0-di1)*(1d0-di2)+v(is,js+ih1,k)*di1*(1d0-di2)&
             +v(i,js,k+ih2)*(1d0-di1)*di2+v(i,js+ih1,k+ih2)*di1*di2
        !same for v
        if (xh(i).ge.xc(1)) then 
           ih1 = -1
        else
           ih1 = 1
        endif
        if (zh(k).ge.zc(1)) then 
           ih2 = -1
        else
           ih2 = 1
        endif
        !get cuts, will be non-dimensional if dx=dy=dz
        di1 = ABS(xc(1)-xh(i))/ABS(yc(1)-y(is-1))
        di2 = ABS(zc(1)-zh(k))/ABS(yc(1)-y(is-1))
        !set
        u(i,js-1,k)=u(i,js,k)*(1d0-di1)*(1d0-di2)+u(i+ih1,js,k)*di1*(1d0-di2)&
             +u(i,js,k+ih2)*(1d0-di1)*di2+u(i+ih1,js,k+ih2)*di1*di2
        !same for w
        w(i,js-1,k)=w(i,js,k)*(1d0-di1)*(1d0-di2)+w(i+ih1,js,k)*di1*(1d0-di2)&
             +w(i,js,k+ih2)*(1d0-di1)*di2+w(i+ih1,js,k+ih2)*di1*di2
        !set u(is-2) to allow cell is-1 to be divergence free
     enddo; enddo
  endif
  !y+ face
  if (coords(2)==nPy-1) then
     !write(*,*)'y+ face loop'
     do i=is,ie; do k=ks,ke
        if (x(i).ge.xc(1)) then 
           ih1 = -1
        else
           ih1 = 1
        endif
        if (z(k).ge.zc(1)) then 
           ih2 = -1
        else
           ih2 = 1
        endif
        !get cuts, will be non-dimensional if dx=dy=dz
        di1 = ABS(xc(1)-x(i))/ABS(yc(1)-yh(je))
        di2 = ABS(zc(1)-z(k))/ABS(yc(1)-yh(je))
        !set
        v(i,je,k)=v(i,je-1,k)*(1d0-di1)*(1d0-di2)+v(i+ih1,je-1,k)*di1*(1d0-di2)&
             +v(i,je-1,k+ih2)*(1d0-di1)*di2+v(i+ih1,je-1,k+ih2)*di1*di2
        !same for u
        if (xh(i).ge.xc(1)) then 
           ih1 = -1
        else
           ih1 = 1
        endif
        if (zh(k).ge.zc(1)) then 
           ih2 = -1
        else
           ih2 = 1
        endif
        !get cuts, will be non-dimensional if dx=dy=dz
        di1 = ABS(xc(1)-xh(i))/ABS(yc(1)-y(je+1))
        di2 = ABS(zc(1)-zh(k))/ABS(yc(1)-y(je+1))
        !set
        u(i,je+1,k)=u(i,je,k)*(1d0-di1)*(1d0-di2)+u(i+ih1,je,k)*di1*(1d0-di2)&
             +u(i,je,k+ih2)*(1d0-di1)*di2+u(i+ih1,je,k+ih2)*di1*di2
        !same for w
        w(i,je+1,k)=w(i,je,k)*(1d0-di1)*(1d0-di2)+w(i+ih1,je,k)*di1*(1d0-di2)&
             +w(i,je,k+ih2)*(1d0-di1)*di2+w(i+ih1,je,k+ih2)*di1*di2
        !set u(ie+2) to allow cell to be divergence free
     enddo; enddo
  endif
end subroutine set_RP_radial_velocity
!====================================================================================================================================================
subroutine set_RP_pressure(p)
  use module_grid
  use module_2phase
  use module_freesurface
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: p
  real(8) :: r, P_l
  integer :: i,j,k
  P_l = P_ref*(R_ref/R_RK)**(3d0*gamma)-2d0*sigma/R_RK
  if (coords(1)==0) then 
     do k=ks-1,ke+1; do j=js-1,je+1
        r = sqrt((x(is-1)-xc(1))**2d0 + (y(j)-yc(1))**2d0 + (z(k)-zc(1))**2d0)
        p(is-1,j,k) = P_l - (dR_RK**2d0 * R_RK**4d0/(2d0*r**4d0) - (ddR_RK*R_RK**2d0 + 2d0*R_RK*dR_RK**2d0)/r +&
             ddR_RK*R_RK + 3d0/2d0*dR_RK**2d0)
     enddo; enddo
  endif
  if(coords(1)==Npx-1) then
     do k=ks-1,ke+1; do j=js-1,je+1
        r = sqrt((x(ie+1)-xc(1))**2d0 + (y(j)-yc(1))**2d0 + (z(k)-zc(1))**2d0)
        p(ie+1,j,k) = P_l - (dR_RK**2d0 * R_RK**4d0/(2d0*r**4d0) - (ddR_RK*R_RK**2d0 + 2d0*R_RK*dR_RK**2d0)/r +&
             ddR_RK*R_RK + 3d0/2d0*dR_RK**2d0)
     enddo; enddo 
  endif
  if(coords(2)==0) then
     do k=ks-1,ke+1; do i=is-1,ie+1
        r = sqrt((x(i)-xc(1))**2d0 + (y(js-1)-yc(1))**2d0 + (z(k)-zc(1))**2d0)
        p(i,js-1,k) = P_l - (dR_RK**2d0 * R_RK**4d0/(2d0*r**4d0) - (ddR_RK*R_RK**2d0 + 2d0*R_RK*dR_RK**2d0)/r +&
             ddR_RK*R_RK + 3d0/2d0*dR_RK**2d0)
     enddo; enddo 
  endif
  if(coords(2)==Npy-1) then
    do k=ks-1,ke+1; do i=is-1,ie+1
        r = sqrt((x(i)-xc(1))**2d0 + (y(je+1)-yc(1))**2d0 + (z(k)-zc(1))**2d0)
        p(i,je+1,k) = P_l - (dR_RK**2d0 * R_RK**4d0/(2d0*r**4d0) - (ddR_RK*R_RK**2d0 + 2d0*R_RK*dR_RK**2d0)/r +&
             ddR_RK*R_RK + 3d0/2d0*dR_RK**2d0)
     enddo; enddo  
  endif

  ! Pressure BC for z-
  if(coords(3)==0) then
     do j=js-1,je+1; do i=is-1,ie+1
        r = sqrt((x(i)-xc(1))**2d0 + (y(j)-yc(1))**2d0 + (z(ks-1)-zc(1))**2d0)
        p(i,j,ks-1) = P_l - (dR_RK**2d0 * R_RK**4d0/(2d0*r**4d0) - (ddR_RK*R_RK**2d0 + 2d0*R_RK*dR_RK**2d0)/r +&
             ddR_RK*R_RK + 3d0/2d0*dR_RK**2d0)
     enddo; enddo 
  endif
  ! Pressure BC for z+
  if(coords(3)==Npz-1) then
     do j=js-1,je+1; do i=is-1,ie+1
        r = sqrt((x(i)-xc(1))**2d0 + (y(j)-yc(1))**2d0 + (z(ke+1)-zc(1))**2d0)
        p(i,j,ke+1) = P_l - (dR_RK**2d0 * R_RK**4d0/(2d0*r**4d0) - (ddR_RK*R_RK**2d0 + 2d0*R_RK*dR_RK**2d0)/r +&
             ddR_RK*R_RK + 3d0/2d0*dR_RK**2d0)
     enddo; enddo 
  endif
end subroutine set_RP_pressure
!====================================================================================================================================================
subroutine initialize_P_RP(p)
  use module_grid
  use module_2phase
  use module_freesurface
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: p
  real(8) :: r, P_l
  integer :: i,j,k,ierr
  P_l = P_ref*(R_ref/R_RK)**(3d0*gamma)-2d0*sigma/R_RK
  do k=ks,ke; do j=js,je; do i=is,ie
     r = sqrt((x(i)-xc(1))**2d0 + (y(j)-yc(1))**2d0 + (z(k)-zc(1))**2d0)
     if (r .ge. rad(1)) then
        p(i,j,k) = P_l - (dR_RK**2d0 * R_RK**4d0/(2d0*r**4d0) - (ddR_RK*R_RK**2d0 + 2d0*R_RK*dR_RK**2d0)/r +&
             ddR_RK*R_RK + 3d0/2d0*dR_RK**2d0)
     endif
  enddo; enddo; enddo
end subroutine initialize_P_RP
!====================================================================================================================================================
