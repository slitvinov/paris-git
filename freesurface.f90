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
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
! 02111-1307, USA.
!=================================================================================================
!=================================================================================================
!!$subroutine fs_sweep(iout,t)
!!$  use module_grid
!!$  use module_vof
!!$  use module_surface_tension
!!$  use module_freesurface
!!$  use module_IO
!!$  use module_BC
!!$  use module_lag_part
!!$  implicit none
!!$  include 'mpif.h'
!!$  real(8) :: t
!!$  integer :: iout
!!$
!!$  if (check_stray_liquid .and. mod(itimestep,n_stray_liquid)==0) then
!!$     call tag_bubbles(0,iout,t)
!!$     call check_topology(.false.)
!!$     call ReleaseTag2DropTable
!!$     if (fill_ghost) then
!!$        call clean_debris
!!$        call do_all_ghost(cvof)
!!$        call get_flags_and_clip(cvof,vof_flag)
!!$        call get_vof_phase(cvof) !cvof updated above from min to max
!!$     endif
!!$  endif
!!$  if (.not. (test_capwave .or. test_plane)) then
!!$     call tag_bubbles(1,iout,t)
!!$     if (.not.(RP_test) .and. (mod(itimestep,n_check_topo)==0)) then
!!$        call check_topology(.true.)
!!$        if (fill_ghost) then
!!$           call clean_debris
!!$           call do_all_ghost(cvof)
!!$           call get_flags_and_clip(cvof,vof_flag)
!!$           call get_vof_phase(cvof) !cvof updated above from min to max
!!$           call ReleaseTag2DropTable
!!$           call tag_bubbles(1,iout,t)
!!$           call get_all_heights(iout)
!!$        endif
!!$     endif
!!$  endif
!!$  call set_topology(vof_phase,itimestep) !vof_phase updated in vofsweeps
!!$
!!$endsubroutine fs_sweep
!-------------------------------------------------------------------------------------------------
!This subroutine sets up neighbour levels inside the gas (void) phase for bubble calculations.
!The breakdown is like this: level "0" in the liquid (or interface?), "1" first cell layer inside
!the gas, and "2" the second layer. Note that only "1" cells are needed for pressure calculation,
!"2" for velocity.
subroutine set_topology(phase,iout)
  use module_grid
  use module_freesurface
  use module_IO
  use module_BC
  use module_lag_part
  implicit none
  include 'mpif.h'
  integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: phase
  integer :: req(4),sta(MPI_STATUS_SIZE,4)
  integer :: i,j,k,level,iout,ierr
  integer :: level3, l3sum

  !initialize all masks to 3
  if (.not.initialize_fs) call pariserror("Error: Free surface variables not initialized")
  u_cmask = 3; v_cmask = 3; w_cmask = 3; pcmask = 3
  !First loop to set level 0 velocities in liq-liq and liq-gas cells
  do k=ks,ke; do j=js,je; do i=is,ie
     if (phase(i,j,k) == 1) then
        if (phase(i+1,j,k) == 0) then
           u_cmask(i,j,k) = 0
        endif
        if (phase(i,j+1,k) == 0) then
           v_cmask(i,j,k) = 0
        endif
        if (phase(i,j,k+1) == 0) then
           w_cmask(i,j,k) = 0
        endif
     endif

!!$     if (phase(i,j,k) == 0 .or. (implode_flag(tag_id(i,j,k))) ) then
     if (phase(i,j,k) == 0) then
        pcmask(i,j,k)=0
        u_cmask(i,j,k)=0
        v_cmask(i,j,k)=0
        w_cmask(i,j,k)=0
     endif
  enddo; enddo; enddo
  call ighost_x(u_cmask,2,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
  call ighost_y(u_cmask,2,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
  call ighost_z(u_cmask,2,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
  call ighost_x(v_cmask,2,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
  call ighost_y(v_cmask,2,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
  call ighost_z(v_cmask,2,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
  call ighost_x(w_cmask,2,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
  call ighost_y(w_cmask,2,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
  call ighost_z(w_cmask,2,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
  call ighost_x(pcmask,2,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
  call ighost_y(pcmask,2,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
  call ighost_z(pcmask,2,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
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
     call ighost_x(u_cmask,2,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
     call ighost_y(u_cmask,2,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
     call ighost_z(u_cmask,2,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
     call ighost_x(v_cmask,2,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
     call ighost_y(v_cmask,2,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
     call ighost_z(v_cmask,2,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
     call ighost_x(w_cmask,2,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
     call ighost_y(w_cmask,2,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
     call ighost_z(w_cmask,2,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
     call ighost_x(pcmask,2,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
     call ighost_y(pcmask,2,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
     call ighost_z(pcmask,2,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
  enddo
end subroutine set_topology
!-------------------------------------------------------------------------------------------------
!checking for: imploding bubbles, and for the liquid-inside-bubble entrapment situation
subroutine check_topology(is_gas)
  use module_grid
  use module_2phase
  use module_freesurface
  use module_IO
  use module_flow
  use module_BC
  use module_Lag_part
  implicit none
  include 'mpif.h'
  real(8) :: volume
  integer :: dropid
  integer :: i,j,k,bub,ierr,rankid
  logical :: remove
  logical, intent(in) :: is_gas

  implode_flag = .false.
  fill_ghost = .false.
  remove = .false.
  v_source=0.d0

  if (num_drop(rank)>0) then
     do bub=1,num_drop(rank)
        volume = drops(bub)%element%vol
        dropid = drops(bub)%element%id
        !write(*,'("Volume of drop ",I5," in rank ",I5", : ",2e14.5)')bub,rank,volume,volume/(dx(is)**3.d0)
        if (volume > 1.d-9*dx(is)**3.d0) then
           if (is_gas) then
              if (volume < 125.0*dx(is)**3.d0) then
                 !write(*,'("Volume <  125 cells of drop ",I5," in rank ",I5," with id ",I5, " : ",2e14.5)')&
                 !     bub,rank,dropid,volume,volume/(dx(is)**3.d0)
                 call bub_implode(dropid,.false.)
              endif
           else
              !if (volume<0.4*(xLength*yLength*zLength))&
              !     write(*,'("Single element ",I5," found in rank ",I5,"with volume: "e14.5)')bub,rank,volume
              if (volume < 50.0*dx(is)**3.d0) then
                 write(*,'("Removing detached liquid in rank: ",I5)')rank
                 call clear_stray_liquid(dropid)
                 remove = .true.
              endif
           endif
        else
           write(*,'("Bubble volume error in topology check. Vol from table: ",e14.5)')volume
        endif
     enddo
  endif

  if ( num_drop_merge(rank) > 0 ) then
     do bub=1,num_drop_merge(rank)
        volume = drops_merge(bub)%element%vol
        dropid = drops_merge(bub)%element%id
        !write(*,'("Volume of drop_merge ",I5," in rank ",I5", : ",e14.5)')bub,rank,volume
        if (volume > 1.d-9*dx(is)**3.d0) then
           if (is_gas) then
              if (volume < 125.0*dx(is)**3.d0) then
                 !write(*,'("Volume < of drop_merge ",I5," in rank ",I5", : ",e14.5)')bub,rank,volume
                 call bub_implode(dropid,.true.)
              endif
           else
              if (volume<0.4*(xLength*yLength*zLength))&
                   write(*,'("Merged element ",I5," found in rank ",I5,"with volume: "e14.5)')bub,rank,volume
              if (volume < 50.0*dx(is)**3.d0) then
                 call clear_stray_liquid(dropid) !<-clears liquid inside bubbles
                 write(*,'("Removing detached liquid (merged) in rank: ",I5)')rank
                 remove = .true.
              endif
           endif
        else
           write(*,'("Bubble volume error in topology check. Vol from table: ",e14.5)')volume
        endif
     enddo
  endif
  call MPI_ALLREDUCE(remove,fill_ghost, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_Active, ierr)
contains
  subroutine bub_implode(bub_id,merged) !uses bub_id from lpp-type procedures
    use module_vof
    use module_flow
    implicit none
    integer :: i,j,k,bub_id,inbr,jnbr,knbr
    integer :: req(4),sta(MPI_STATUS_SIZE,4),ierr
    real(8) :: d_clean
    integer :: max_implode
    logical :: merged
    max_implode = 0
    remove =.false.
!!$    do k=ks,ke; do j=js,je; do i=is,ie
!!$
!!$       if (tag_id(i,j,k)==bub_id) then !check if we are in the correct bubble
!!$          if (.not.implode_flag(tag_id(i,j,k))) then
!!$             implode_flag(tag_id(i,j,k)) = .true.
!!$             write(*,'("COLLAPSING BUBBLE DETECTED IN RANK: ",I4)')rank
!!$          endif
!!$          v_source(i,j,k) = (u(i,j,k)-u(i-1,j,k))/dx(i)+&
!!$               (v(i,j,k)-v(i,j-1,k))/dy(j)+(w(i,j,k)-w(i,j,k-1))/dz(k)
!!$       endif
!!$    enddo; enddo; enddo

    if (volume <= 64.0*dx(is)**3.d0) then
       remove =.true.
       write(*,'("Collapsing bubble to be removed in rank ",I4", with rank ID: ",I4)')rank,bub_id
    endif
    if (remove) then
       do k=ks,ke; do j=js,je; do i=is,ie
          if (tag_id(i,j,k)==bub_id) then
             cvof(i,j,k)=0.0d0
             v_source(i,j,k) = 0.0d0
             if (itime_scheme==2) cvofold(i,j,k)=0.0d0
             implode_flag(tag_id(i,j,k))=.false.
          endif
       enddo; enddo; enddo
    endif

  end subroutine bub_implode
  !-----------------------------------------------------------------------
  subroutine clear_stray_liquid(bub_id)
    use module_vof
    implicit none
    integer :: i,j,k,bub_id
    do k=ks,ke; do j=js,je; do i=is,ie
       if (tag_id(i,j,k)==bub_id) then
          cvof(i,j,k)=1.0d0
       endif
    enddo; enddo; enddo
  end subroutine clear_stray_liquid
end subroutine check_topology
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
  integer :: inbr, jnbr, knbr
  real(8), dimension(-2:2,-2:2,-2:2) :: vel_val
  integer, dimension(-2:2,-2:2,-2:2) :: velmask
  logical :: x_success
  real(8) :: varx

  !Simple 1st order extrapolation: Velocities of neighbours at lower topological level averaged
  !Applies the linear equation from FreeSurface documentation (report)
  if (order_extrap == 1) then
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
        call do_ghost_vector(u,v,w) !<mpi
     enddo
  else if (order_extrap ==2) then
     do level = 1, X_level !Loop extrapolation levels
        do k=ks,ke; do j=js,je; do i=is,ie !Loop i,j,k
           !setup matrix
           ! Zero sums
           if (u_cmask(i,j,k) == level) then
              !get 5x5x5 mask
              do knbr=-2,2; do jnbr=-2,2; do inbr=-2,2
                 velmask(inbr,jnbr,knbr)=u_cmask(i+inbr,j+jnbr,k+knbr)
                 vel_val(inbr,jnbr,knbr)=u(i+inbr,j+jnbr,k+knbr)
              enddo; enddo; enddo
              call setupmatrix(velmask,vel_val,level,varx,x_success)
              if (x_success) u(i,j,k)=varx
           endif
           if (v_cmask(i,j,k) == level) then
              !get 5x5x5 mask
              do knbr=-2,2; do jnbr=-2,2; do inbr=-2,2
                 velmask(inbr,jnbr,knbr)=v_cmask(i+inbr,j+jnbr,k+knbr)
                 vel_val(inbr,jnbr,knbr)=v(i+inbr,j+jnbr,k+knbr)
              enddo; enddo; enddo
              call setupmatrix(velmask,vel_val,level,varx,x_success)
              if (x_success) v(i,j,k)=varx
           endif
           if (w_cmask(i,j,k) == level) then
              !get 5x5x5 mask
              do inbr=-2,2; do jnbr=-2,2; do knbr=-2,2
                 velmask(inbr,jnbr,knbr)=w_cmask(i+inbr,j+jnbr,k+knbr)
                 vel_val(inbr,jnbr,knbr)=w(i+inbr,j+jnbr,k+knbr)
              enddo; enddo; enddo
              call setupmatrix(velmask,vel_val,level,varx,x_success)
              if (x_success) w(i,j,k)=varx
           endif
        enddo; enddo; enddo !end loop i,j,k
        call do_ghost_vector(u,v,w) !mpi again
     enddo !end loop levels
  endif
contains
  subroutine setupmatrix(topmask,vel,lvl,var,extra_vel_found)
    use module_surface_tension
    implicit none
    real(8) :: var
    real(8), dimension(-2:2,-2:2,-2:2), intent(in) :: vel
    integer, dimension(-2:2,-2:2,-2:2), intent(in) :: topmask
    integer :: lvl,row,col,mult
    integer :: inr,jnr,knr
    real(8), dimension(1:4,1:4) :: x_m, x_im, test
    real(8), dimension(1:4) :: dh, rhs
    real(8) :: den, psi, maxc
    integer :: fit_cells
    logical :: inverse_success, extra_vel_found
    var = 0.0d0
    fit_cells=0
    rhs=0.0d0
    x_m =0.0d0
    !maxc=0.0d0
    !check eligible neighbours
    extra_vel_found=.false.
    do knr=-2,2; do jnr=-2,2; do inr=-2,2
       if (topmask(inr,jnr,knr)<=lvl-1) then
          fit_cells=fit_cells+1
          den = inr**2.0d0 + jnr**2.0d0 + knr**2.0d0
          if (den >= 1.0d0) then
             psi=1.0d0/den
          else
             psi=1.0d0
          endif
          dh(1)=1.0d0
          dh(2)=inr*dx(is) !only unstretched regular grids
          dh(3)=jnr*dy(js)
          dh(4)=knr*dz(ks)
          do col=1,4; do row=1,4
             x_m(row,col) = x_m(row,col) + dh(row)*dh(col)*psi
          enddo; enddo
          do row=1,4
             rhs(row) = rhs(row) + dh(row)*vel(inr,jnr,knr)*psi
          enddo
       endif
    enddo; enddo; enddo
    do row=1,4
       maxc = 0.0d0
       do col=1,4
          maxc=MAX(maxc,x_m(row,col))
       enddo
       if (ABS(maxc)>1d-10) then
          do col=1,4
             x_m(row,col)=x_m(row,col)/maxc
          enddo
          rhs(row) = rhs(row)/maxc
       endif
    enddo
    !solve matrix locally
    if (fit_cells>3) then
       !calling Stanley's FindInverseMatrix to solve for the extrapolated
       !velocity, since it is a matrix equation
       call FindInverseMatrix(x_m,x_im,4,inverse_success)
       if (inverse_success) then
!!$          !check product matrix*inverse
!!$          test = 0.0d0
!!$          do row=1,4; do col=1,4
!!$             test(row,col) = 0.0d0
!!$             do mult=1,4
!!$                test(row,col) = test(row,col) + x_m(row,mult)*x_im(mult,col)
!!$             enddo; enddo
!!$          enddo
!!$          write(*,'("Product of  matrix and inverted matrix: ")')
!!$          do row=1,4
!!$             write(*,'(4e16.5)')test(row,1:4)
!!$          enddo
!!$          write(*,*)' '
          !get extrapolated velocity
          !write(*,'("Matrix inverted successfully using ",I4, "points.")')fit_cells
          do col=1,4
             !write(*,'("Inverted row, col",I4,e14.5)')col,x_im(4,col)
             var=var+x_im(1,col)*rhs(col)
          enddo
          !write(*,'("Extrapolated velocity: ",e14.5)')var
          !debug
          if (var /= var) then
             open(unit=70,file="Extr_v_NaN.txt",position="append")
             write(70,'("Extrapolated vel NaN after fit success at time step: ",I10)')itimestep
             write(70,'("ijk rank",4I4)')i,j,k,rank
             write(70,'("x, y, z: ",3e14.5)')x(i),y(j),z(k)
             write(70,'("Cvof 1-7: ",7e14.5)')cvof(i-1,j,k),cvof(i+1,j,k),cvof(i,j-1,k),cvof(i,j+1,k),&
                  cvof(i,j,k-1),cvof(i,j,k+1),cvof(i,j,k)
             write(70,'("  ")')
             write(70,'("Topology mask: ")')
             do jnr=-2,2; do knr=-2,2
                write(70,'(" ",3I4," :",e14.5)')topmask(-2:2,jnr,knr)
             enddo; enddo
             write(70,'("  ")')
             write(70,'("Velocity mask: ")')
             do jnr=-2,2; do knr=-2,2
                write(70,'(" ",3I4," :",e14.5)')vel(-2:2,jnr,knr)
             enddo; enddo
             write(70,'("  ")')
             write(70,'("Phase 1-7: ",7I8)')vof_phase(i-1,j,k),vof_phase(i+1,j,k),vof_phase(i,j-1,k),vof_phase(i,j+1,k),&
                  vof_phase(i,j,k-1),vof_phase(i,j,k+1),vof_phase(i,j,k)
             write(70,'("Pcmask 1-7: ",7I8)')pcmask(i-1,j,k),pcmask(i+1,j,k),pcmask(i,j-1,k),pcmask(i,j+1,k),&
                  pcmask(i,j,k-1),pcmask(i,j,k+1),pcmask(i,j,k)
             write(70,'("  ")')
             close(70)
          else
             extra_vel_found=.true.
          endif
       else
          !write(*,*)'Error: Matrix could not be inverted for velocity extrapolation, reverting to 1st order'
          !write(*,'("LS fit failed using ",I5," points")')fit_cells
!!$          do row=1,4; do col=1,4
!!$             write(*,'("A",2I4,e14.5)')row,col,x_m(row,col)
!!$          enddo; enddo
          !revert to 1st order
          var=rhs(1)/(1.0d0*fit_cells)
          if (var /= var) then
             open(unit=70,file="Extr_v_NaN.txt",position="append")
             write(70,'("Extrapolated vel NaN reverting 1st, > 3, at time step: ",I10)')itimestep
             write(70,'("ijk rank",4I4)')i,j,k,rank
             write(70,'("x, y, z: ",3e14.5)')x(i),y(j),z(k)
             write(70,'("Cvof 1-7: ",7e14.5)')cvof(i-1,j,k),cvof(i+1,j,k),cvof(i,j-1,k),cvof(i,j+1,k),&
                  cvof(i,j,k-1),cvof(i,j,k+1),cvof(i,j,k)
             write(70,'("  ")')
             write(70,'("Topology mask: ")')
             do jnr=-2,2; do knr=-2,2
                write(70,'(" ",3I4," :",e14.5)')topmask(-2:2,jnr,knr)
             enddo; enddo
             write(70,'("  ")')
             write(70,'("Velocity mask: ")')
             do jnr=-2,2; do knr=-2,2
                write(70,'(" ",3I4," :",e14.5)')vel(-2:2,jnr,knr)
             enddo; enddo
             write(70,'("  ")')
             write(70,'("Phase 1-7: ",7I8)')vof_phase(i-1,j,k),vof_phase(i+1,j,k),vof_phase(i,j-1,k),vof_phase(i,j+1,k),&
                  vof_phase(i,j,k-1),vof_phase(i,j,k+1),vof_phase(i,j,k)
             write(70,'("Pcmask 1-7: ",7I8)')pcmask(i-1,j,k),pcmask(i+1,j,k),pcmask(i,j-1,k),pcmask(i,j+1,k),&
                  pcmask(i,j,k-1),pcmask(i,j,k+1),pcmask(i,j,k)
             write(70,'("  ")')
             close(70)
          else
             extra_vel_found=.true. !1st order
          endif
       endif
    else
       if (fit_cells >= 1) then
          var=rhs(1)/(1.0d0*fit_cells)
          if (var /= var) then
             open(unit=70,file="Extr_v_NaN.txt",position="append")
             write(70,'("Extrapolated vel NaN reverting 1st, < 3, at time step: ",I10)')itimestep
             write(70,'("ijk rank",4I4)')i,j,k,rank
             write(70,'("x, y, z: ",3e14.5)')x(i),y(j),z(k)
             write(70,'("Cvof 1-7: ",7e14.5)')cvof(i-1,j,k),cvof(i+1,j,k),cvof(i,j-1,k),cvof(i,j+1,k),&
                  cvof(i,j,k-1),cvof(i,j,k+1),cvof(i,j,k)
             write(70,'("  ")')
             write(70,'("Topology mask: ")')
             do jnr=-2,2; do knr=-2,2
                write(70,'(" ",3I4," :",e14.5)')topmask(-2:2,jnr,knr)
             enddo; enddo
             write(70,'("  ")')
             write(70,'("Velocity mask: ")')
             do jnr=-2,2; do knr=-2,2
                write(70,'(" ",3I4," :",e14.5)')vel(-2:2,jnr,knr)
             enddo; enddo
             write(70,'("  ")')
             write(70,'("Phase 1-7: ",7I8)')vof_phase(i-1,j,k),vof_phase(i+1,j,k),vof_phase(i,j-1,k),vof_phase(i,j+1,k),&
                  vof_phase(i,j,k-1),vof_phase(i,j,k+1),vof_phase(i,j,k)
             write(70,'("Pcmask 1-7: ",7I8)')pcmask(i-1,j,k),pcmask(i+1,j,k),pcmask(i,j-1,k),pcmask(i,j+1,k),&
                  pcmask(i,j,k-1),pcmask(i,j,k+1),pcmask(i,j,k)
             write(70,'("  ")')
             close(70)
          else
             extra_vel_found=.true.
          endif
       else
          write(*,'("WARNING, ISOLATED TOPOLOGICAL STRUCTURE, NO NEIGHBOURS! Level: ",I4)')lvl
          write(*,'("Neighbouring level values:")')
          do inr=-2,2; do jnr=-2,2; do knr=-2,2
             write(*,'("Value at ",3I4," :",e14.5)')inr,jnr,knr,topmask(inr,jnr,knr)
          enddo; enddo; enddo
       endif
    endif
  end subroutine setupmatrix
end subroutine extrapolate_velocities
!=============================================================================================================================================
!for a divergence correction in the 2nd projection
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

  do k=ks,ke; do j=js,je; do i=is,ie
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
!Older (obsoleted) procedure to get reference pressure in single bubble Rayleigh-Plesset test
!"don't worry about it too much:)"
subroutine get_ref_volume(vof)
  use module_grid
  use module_2phase
  use module_BC
  use module_freesurface
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: vof
  real(8), parameter :: pi=3.1415926535897932384626433833
  real(8) :: V_loc
  integer :: i,j,k,ierr

  !V_0 = 4d0/3d0*pi*(R_ref**3d0)
  V_loc = 0.0d0
  do k=ks,ke; do j=js,je; do i=is,ie
     V_loc=V_loc+vof(i,j,k)*dx(i)*dy(j)*dz(k)
  enddo; enddo; enddo
  call MPI_ALLREDUCE(V_loc, V_0, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_Domain, ierr)
  R_RK = rad(1)
  dR_RK = 0d0
  P_inf = BoundaryPressure(1)
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
  real(8) :: res_local, resinf_all
  integer :: i,j,k, iout
  integer :: req(12),sta(MPI_STATUS_SIZE,12)
  logical :: mask(imin:imax,jmin:jmax,kmin:kmax)
  integer, parameter :: norm=2, relaxtype=1
  integer, save :: itime=0
  logical :: use_L_inf = .true.
  ! Open file for convergence history
  if(rank==0.and.recordconvergence) then
     OPEN(UNIT=89,FILE=TRIM(out_path)//'/convergence_history-'//TRIM(int2text(itime,padding))//'.txt')
  endif
  itime=itime+1
2 format(2e14.5)
  if (solver_flag == 0) call pariserror("Free Surface solver flag needs to be 1 or 2")
  do k=ks,ke; do j=js,je; do i=is,ie
     if (solver_flag==1 .and. pcmask(i,j,k)/=0) p(i,j,k) = P_gas(i,j,k)
     if (solver_flag==2 .and. pcmask(i,j,k)==3) p(i,j,k) = 0.d0
  enddo; enddo; enddo
  call ghost_x(p,1,req( 1: 4)); call ghost_y(p,1,req( 5: 8)); call ghost_z(p,1,req( 9:12))
  call MPI_WAITALL(12,req,sta,ierr)
  !--------------------------------------ITERATION LOOP--------------------------------------------
  do it=1,maxit
     if(relaxtype==2) then
        call LineRelax_fs(A,p,beta)
     elseif(relaxtype==1) then
        call RedBlackRelax_fs(A,p,beta)
     endif
     !---------------------------------CHECK FOR CONVERGENCE-------------------------------------------
     res1 = 0d0; res2=0.d0; resinf=0.d0; intvol=0.d0; cells = 0d0; res_local=0.d0
     call ghost_x(p,1,req( 1: 4)); call ghost_y(p,1,req( 5: 8)); call ghost_z(p,1,req( 9:12))
     do k=ks+1,ke-1; do j=js+1,je-1; do i=is+1,ie-1
        if ((pcmask(i,j,k)==0 .and. solver_flag==1)&
             .or.((pcmask(i,j,k)==1 .or. pcmask(i,j,k)==2) .and. solver_flag==2)) then
           res_local = abs(-p(i,j,k) * A(i,j,k,7) +                           &
                A(i,j,k,1) * p(i-1,j,k) + A(i,j,k,2) * p(i+1,j,k) +            &
                A(i,j,k,3) * p(i,j-1,k) + A(i,j,k,4) * p(i,j+1,k) +            &
                A(i,j,k,5) * p(i,j,k-1) + A(i,j,k,6) * p(i,j,k+1) + A(i,j,k,8) )
           res2=res2+res_local**norm
           resinf=MAX(resinf,res_local)
           cells = cells + 1.d0
           if (res2/=res2) then
              write(*,*)'ERROR: RES2 NaN'
              write(*,'("Res_local, res2: ",2e14.5)')res_local,res2
              call debug_details(i,j,k,A)
              call pariserror('FreeSolver Res NaN')
           endif
        endif
     enddo; enddo; enddo
     call MPI_WAITALL(12,req,sta,ierr)
     mask=.true.
     mask(is+1:ie-1,js+1:je-1,ks+1:ke-1)=.false.
     do k=ks,ke; do j=js,je; do i=is,ie
        if(mask(i,j,k) .and. ((pcmask(i,j,k)==0 .and. solver_flag==1) .or. &
             ((pcmask(i,j,k)==1 .or. pcmask(i,j,k)==2) .and. solver_flag==2))) then
           res_local=abs(-p(i,j,k) * A(i,j,k,7) +&
                A(i,j,k,1) * p(i-1,j,k) + A(i,j,k,2) * p(i+1,j,k) +            &
                A(i,j,k,3) * p(i,j-1,k) + A(i,j,k,4) * p(i,j+1,k) +            &
                A(i,j,k,5) * p(i,j,k-1) + A(i,j,k,6) * p(i,j,k+1) + A(i,j,k,8) )
           res2=res2+res_local**norm
           resinf=MAX(resinf,res_local)
           cells = cells + 1d0
           if (res2/=res2) then
              write(*,*)'ERROR: RES2 NaN, proc border'
              write(*,'("Res_local, res2: ",e14.5)')res_local,res2
              call debug_details(i,j,k,A)
              call pariserror('FreeSolver Res NaN in proc border')
           endif
        endif
     enddo; enddo; enddo
     call catch_divergence_fs(res2,cells,ierr)
     call MPI_ALLREDUCE(res2, tres2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr)
     call MPI_ALLREDUCE(resinf, resinf_all, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_Comm_Cart, ierr)
     call MPI_ALLREDUCE(cells, tcells, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr)
     if(norm==2) tres2=sqrt(tres2)
     if (tcells >= 0) then
        tres2 = tres2/tcells
     else
        tres2 = 0
     endif
     if(rank==0.and.mod(it,10) == 0.and.recordconvergence) write(89,310) it, solver_flag, tres2
310  format(2I6,'  ',(e14.5))
     if (.not. use_L_inf) then
        if (tres2<maxError) then
           if(rank==0.and.recordconvergence) close(89)
           exit
        endif
     else
        if (resinf_all<maxError) then
           if(rank==0.and.recordconvergence) close(89)
           exit
        endif
     endif
  enddo

  if(rank==0.and.recordconvergence) close(89)
  if(it==maxit+1 .and. rank==0 .and. solver_flag==1) then
     write(*,*) 'Warning: LinearSolver reached maxit: ||res||: ',tres2
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
    if ((res2*npx*npy*npz)>1.d18 ) then
       if(rank<=100) then
          print*,'Pressure solver diverged after',it,'iterations at rank ',rank
          write(*,'("Res2, cells, solver_flag :",2e14.5,I8)')res2,cells,solver_flag
       endif
       call pariserror("freesolver error")
    else if (res2 .ne. res2) then
       if(rank<=100) then
          print*, 'it:',it,'Pressure residual value is invalid at rank', rank
          write(*,'("Res2, cells, solver_flag :",2e14.5,I8)')res2,cells,solver_flag
       endif
       call pariserror("freesolver error")
    else
       ierr=0
    endif
  end subroutine catch_divergence_fs
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
                 if (p(i,j,k)/=p(i,j,k)) then
                    write(*,*)'ERROR: P NaN'
                    call debug_details(i,j,k,A)
                    call pariserror('P NaN in FreeSolver')
                 endif
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
! Routine to numerically integrate the Rayleigh-Plesset equation (-> for the test suite)
subroutine Integrate_RP(dt,t,rho)
  use module_freesurface
  implicit none
  integer :: j
  real(8) :: dt, t, Volume
  real(8), parameter :: pi=3.141592653589793238462643383
  integer, parameter :: nvar = 2
  real(8) :: y(nvar), ytmp(nvar),rho
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
    dydt(2) = -3d0*(y(2)**2d0)/(2d0*y(1)) + (P_c - P_inf)/(y(1)*rho)
    dydt(1) = y(2)
  end subroutine func
end subroutine Integrate_RP
!=======================================================================================================================================
subroutine write_RP_test(t,rho)
  use module_freesurface
  implicit none
  real(8) :: t, vol, p_mid, p_corner, rho
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
    pressure = P_l - rho*(dR_RK**2d0 * R_RK**4d0/(2d0*r**4d0) - (ddR_RK*R_RK**2d0 + 2d0*R_RK*dR_RK**2d0)/r +&
         ddR_RK*R_RK + 3d0/2d0*dR_RK**2d0)
  end function pressure
end subroutine write_RP_test
!====================================================================================================================================================
!In the RP test where the p_inf is give, this sets the boundary pressure (on domain boundary)
!to a level that would agree with that p_inf.
subroutine set_RP_pressure(p,rho)
  use module_grid
  use module_2phase
  use module_freesurface
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: p
  real(8) :: r, rho
  integer :: i,j,k

  if (coords(1)==0) then
     do k=ks-1,ke+1; do j=js-1,je+1
        r = sqrt((x(is-1)-xc(1))**2d0 + (y(j)-yc(1))**2d0 + (z(k)-zc(1))**2d0)
        p(is-1,j,k) = pressure(r)
     enddo; enddo
  endif
  if(coords(1)==Npx-1) then
     do k=ks-1,ke+1; do j=js-1,je+1
        r = sqrt((x(ie+1)-xc(1))**2d0 + (y(j)-yc(1))**2d0 + (z(k)-zc(1))**2d0)
        p(ie+1,j,k) = pressure(r)
     enddo; enddo
  endif
  if(coords(2)==0) then
     do k=ks-1,ke+1; do i=is-1,ie+1
        r = sqrt((x(i)-xc(1))**2d0 + (y(js-1)-yc(1))**2d0 + (z(k)-zc(1))**2d0)
        p(i,js-1,k) = pressure(r)
     enddo; enddo
  endif
  if(coords(2)==Npy-1) then
     do k=ks-1,ke+1; do i=is-1,ie+1
        r = sqrt((x(i)-xc(1))**2d0 + (y(je+1)-yc(1))**2d0 + (z(k)-zc(1))**2d0)
        p(i,je+1,k) = pressure(r)
     enddo; enddo
  endif

  ! Pressure BC for z-
  if(coords(3)==0) then
     do j=js-1,je+1; do i=is-1,ie+1
        r = sqrt((x(i)-xc(1))**2d0 + (y(j)-yc(1))**2d0 + (z(ks-1)-zc(1))**2d0)
        p(i,j,ks-1) = pressure(r)
     enddo; enddo
  endif
  ! Pressure BC for z+
  if(coords(3)==Npz-1) then
     do j=js-1,je+1; do i=is-1,ie+1
        r = sqrt((x(i)-xc(1))**2d0 + (y(j)-yc(1))**2d0 + (z(ke+1)-zc(1))**2d0)
        p(i,j,ke+1) = pressure(r)
     enddo; enddo
  endif
contains
  function pressure(r)
    real(8) :: P_l, pressure, r
    P_l = P_ref*(R_ref/R_RK)**(3d0*gamma)-2d0*sigma/R_RK
    pressure = P_l - rho*(dR_RK**2d0 * R_RK**4d0/(2d0*r**4d0) - (ddR_RK*R_RK**2d0 + 2d0*R_RK*dR_RK**2d0)/r +&
         ddR_RK*R_RK + 3d0/2d0*dR_RK**2d0)
  end function pressure
end subroutine set_RP_pressure
!====================================================================================================================================================
!initial gas pressure for R-P
subroutine initialize_P_RP(p,rho)
  use module_grid
  use module_2phase
  use module_freesurface
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: p
  real(8) :: r, P_l, rho
  integer :: i,j,k,ierr
  P_l = P_ref*(R_ref/R_RK)**(3d0*gamma)-2d0*sigma/R_RK
  if (rank==0) write(*,'("RP test, P_l, R_RK, dR_RK, ddR_RK: ",4e14.5)')P_l,R_RK, dR_RK, ddR_RK
  do k=ks,ke; do j=js,je; do i=is,ie
     r = sqrt((x(i)-xc(1))**2d0 + (y(j)-yc(1))**2d0 + (z(k)-zc(1))**2d0)
     if (r .ge. rad(1)) then
        p(i,j,k) = P_l - rho*(dR_RK**2d0 * R_RK**4d0/(2d0*r**4d0) - (ddR_RK*R_RK**2d0 + 2d0*R_RK*dR_RK**2d0)/r +&
             ddR_RK*R_RK + 3d0/2d0*dR_RK**2d0)
     else
        p(i,j,k) = P_ref*(R_ref/rad(1))**(3d0*gamma)
     endif
  enddo; enddo; enddo
end subroutine initialize_P_RP
!====================================================================================================================================================
subroutine debug_details(i,j,k,A)
  use module_grid
  use module_flow
  use module_freesurface
  implicit none
  integer :: i,j,k
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(in) :: A
  open(unit=40,file="details.txt",position="append")
  write(40,'("Solver flag:",I5)')solver_flag
  write(40,*) "p",  p(i,j,k)
  write(40,'("A branches: ",8e14.5)')A(i,j,k,:)
  write(40,'("P neighbours ",6e14.5)')p(i-1,j,k),p(i+1,j,k),p(i,j-1,k),p(i,j+1,k),&
       p(i,j,k-1),p(i,j,k+1)
  write(40,'("Mods ",6e14.5)')x_mod(i-1,j,k),x_mod(i,j,k),y_mod(i,j-1,k),y_mod(i,j,k),&
       z_mod(i,j,k-1),z_mod(i,j,k)
  write(40,'("Limits: ",6I8)')is,ie,js,je,ks,ke
  write(40,*) "ijk rank",i,j,k,rank
  write(40,'("Pcmask 1-7: ",7I8)')pcmask(i-1,j,k),pcmask(i+1,j,k),pcmask(i,j-1,k),pcmask(i,j+1,k),&
       pcmask(i,j,k-1),pcmask(i,j,k+1),pcmask(i,j,k)
  close(40)
end subroutine debug_details
!==================================================================================================================
subroutine tag_bubbles(phase_ref,iout,time_stats)
  use module_grid
  use module_Lag_part
  use module_freesurface
  implicit none
  include 'mpif.h'
  real(8), intent(in)  :: time_stats
  integer :: iout, flag
  integer, intent(in) :: phase_ref
  integer :: ierr
  logical, allocatable :: implode_global(:)

  tracked_phase = phase_ref
  call tag_drop()
  if ( nPdomain > 1 ) call tag_drop_all
  call CreateTag2DropTable
  if ( nPdomain > 1 ) call merge_drop_pieces
  if ( MOD(iout,nstats) == 0 ) then
     call drop_statistics(iout,time_stats)
     !! DEBUGGING
!!$     call output_tag(iout/nstats,is,ie+1,js,je+1,ks,ke+1)
!!$     allocate( implode_global(0:NumBubble*28) )
!!$     call MPI_ALLREDUCE(implode_flag,implode_global, (NumBubble*28+1), MPI_LOGICAL, MPI_LOR, MPI_COMM_Active, ierr)
!!$     if (rank == 0) then
!!$        open(unit=77,file='implode_flags.txt',position='append')
!!$        write(77,'("Implode flags at time step: ",I8)')iout
!!$        write(77,'("Time                      : ",e14.5)')time
!!$        do flag = 0, total_num_tag
!!$           write(77,'(I8,": ",L3)')flag,implode_global(flag)
!!$        enddo
!!$        write(77,'(" ")')
!!$        close(77)
!!$     endif
!!$     deallocate( implode_global )
  endif
end subroutine tag_bubbles
!==================================================================================================================
!This allows to set a non-zero P_gas inside the bubbles
subroutine set_bubble_pressure
  use module_grid
  use module_Lag_part
  use module_freesurface
  use module_IO
  use module_2phase
  implicit none
  !real(8), dimension(imin:imax,jmin:jmax,kmin:kmax) :: P_g
  real(8) :: volume
  integer :: i,j,k,dropid

  P_gas = 0.0d0
  if ((NumBubble>0) .and. (P_ref > 1.d-14) .and. .not. (test_capwave .or. test_plane)) then
     do k=ks,ke; do j=js,je; do i=is,ie
        if (pcmask(i,j,k) /= 0) then
           dropid = tag_dropid(tag_id(i,j,k))
           if (.not.(tag_mergeflag(tag_id(i,j,k)) == 1 .or. tag_mergeflag(tag_id(i,j,k)) == 0)) then
              call pariserror('Merge tag should be 0 or 1')
           endif
           if (tag_mergeflag(tag_id(i,j,k)) == 1) then
              volume = drops_merge(dropid)%element%vol
           else
              volume = drops(dropid)%element%vol
           endif
           if (volume > 1.0d-1*dx(is)**3.0d0) then
              P_gas(i,j,k) = P_ref*(V_0/volume)**gamma
           else
              write(*,'("Bubble volume error in FreeSolver. Vol from table: ",e14.5)')volume
           endif
        endif
     enddo;enddo;enddo
  else
     P_gas = P_ref
  endif
end subroutine set_bubble_pressure
!==================================================================================================================
!Smoothens (in time, iters) the boundary velocity condition (check out the itimestep condition).
!could theoretically be used to manipulate this condition in time
subroutine inflow_accelerate
  use module_grid
  use module_BC
  use module_freesurface
  use module_flow
  implicit none
  real(8) :: factor,ix=0.d0
  real(8), parameter :: pi=4.d0*atan(1.d0)
  integer :: i,j,k
  !20160628 new variable fs_refactor added [WA]

  if (itimestep < 2) then
     factor = 0.0d0
  else
     if(fs_refactor) then
        !WA: factor varies sinusoidally, which
        !is stretched over entire time of simulation
        ix=time/endtime
        factor = sin(2.d0*pi*ix)
     else
        !LCM: factor varies smoothly from 0 to 1:
        !which takes step_max steps
        ix=DBLE((itimestep)*1.0)/DBLE(step_max*1.0)
        factor = 0.5d0*(1.d0 + sin(-pi/2.d0 +ix*pi))
     endif

  endif
  !if (rank==0) Write(*,'("Factor: ",e14.5)')factor
  ! inflow boundary condition x- with injection
  if(bdry_cond(1)==3 .and. coords(1)==0    ) then
     do j=jmin,jmax
        do k=kmin,kmax
           u(is-1,j,k)=u(is-1,j,k)*factor
           u(is-2,j,k)=u(is-2,j,k)*factor
        enddo
     enddo
  endif
  ! inflow boundary condition y-
  if(bdry_cond(2)==3 .and. coords(2)==0    ) then
     do i=imin,imax
        do k=kmin,kmax
           v(i,js-1,k)=v(i,js-1,k)*factor
           v(i,js-2,k)=v(i,js-2,k)*factor
        enddo
     enddo
  endif
  ! inflow on z-
  if(bdry_cond(3)==3 .and. coords(3)==0   ) then
     do i=imin,imax
        do j=jmin,jmax
           w(i,j,ks-1)= w(i,j,ks-1)*factor
           w(i,j,ks-2)= w(i,j,ks-2)*factor
        enddo
     enddo
  endif
  ! inflow on x+
  if(bdry_cond(4)==3 .and. coords(1)==nPx-1   ) then
     do j=jmin,jmax
        do k=kmin,kmax
           u(ie,j,k)= u(ie,j,k)*factor
           u(ie+1,j,k)= u(ie+1,j,k)*factor
        enddo
     enddo
  endif
  ! inflow on y+
  if(bdry_cond(5)==3 .and. coords(2)==nPy-1   ) then
     do i=imin,imax
        do k=kmin,kmax
           v(i,je,k)=v(i,je,k)*factor
           v(i,je+1,k)=v(i,je+1,k)*factor
        enddo
     enddo
  endif
  ! inflow on z+
  if(bdry_cond(6)==3 .and. coords(3)==nPz-1   ) then
     do i=imin,imax
        do j=jmin,jmax
           w(i,j,ke)=w(i,j,ke)*factor
           w(i,j,ke+1)=w(i,j,ke+1)*factor
        enddo
     enddo
  endif

  if(rank==0) print *,'step:',itimestep,'t=',time,'tend=',endtime,'x=',ix,'factor:',factor
end subroutine inflow_accelerate
!===================================================================================================================================
subroutine setuppoisson_fs_heights(utmp,vtmp,wtmp,rho,dt,coeff,bub_id)
  use module_grid
  use module_freesurface
  use module_IO
  use module_VOF
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: utmp,vtmp,wtmp
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(out) :: coeff
  integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: bub_id
  real(8), intent(in) :: dt, rho
  integer :: i,j,k,nbr
  integer :: reqd(24),stat(MPI_STATUS_SIZE,24)

  if (.not.(solver_flag==1 .or. solver_flag==2)) call pariserror('Solver_flag FS needs to be either 1 or 2')
  if (solver_flag == 1) then
     call liq_gas2()
  else !(solver_flag == 2)
     if (.not.(solver_flag==2)) call pariserror('ERROR: Solver flag for FS must be set to either 1 or 2')
     do k=ks,ke; do j=js,je; do i=is,ie
        if (pcmask(i,j,k)==1 .or. pcmask(i,j,k)==2) then !rho is 1d0

           coeff(i,j,k,1) = 1.d0/(dx(i)*dxh(i-1))
           coeff(i,j,k,2) = 1.d0/(dx(i)*dxh(i))
           coeff(i,j,k,3) = 1.d0/(dy(j)*dyh(j-1))
           coeff(i,j,k,4) = 1.d0/(dy(j)*dyh(j))
           coeff(i,j,k,5) = 1.d0/(dz(k)*dzh(k-1))
           coeff(i,j,k,6) = 1.d0/(dz(k)*dzh(k))

           if (pcmask(i,j,k)==2 .and. vof_phase(i,j,k)/=1) call pariserror("Fatal topology error in FS 2nd projection")

           if (pcmask(i,j,k)==1) then
              if (vof_phase(i,j,k)/=1) call pariserror("Fatal topology error in FS 2nd projection")
              do nbr=-1,1,2
                 if (pcmask(i+nbr,j,k) == 0) then
                    coeff(i,j,k,2+(nbr-1)/2) = 0d0
                 endif
                 if (pcmask(i,j+nbr,k) == 0) then
                    coeff(i,j,k,4+(nbr-1)/2) = 0d0
                 endif
                 if (pcmask(i,j,k+nbr) == 0) then
                    coeff(i,j,k,6+(nbr-1)/2) = 0d0
                 endif
              enddo
           endif
           coeff(i,j,k,7) = sum(coeff(i,j,k,1:6))
           if (coeff(i,j,k,7)<1.d-12) pcmask(i,j,k) = 0 !Fix for isolated cav cells.
           coeff(i,j,k,8) =  -1.d0*((utmp(i,j,k)-utmp(i-1,j,k))/dx(i) &
                +  (vtmp(i,j,k)-vtmp(i,j-1,k))/dy(j) &
                +  (wtmp(i,j,k)-wtmp(i,j,k-1))/dz(k))
        endif
     enddo; enddo; enddo
  endif
contains
  subroutine liq_gas2()
    use module_BC
    use module_surface_tension
    use module_2phase
    implicit none
    real(8) :: kap(imin:imax,jmin:jmax,kmin:kmax)
    real(8) :: Source
    integer :: i,j,k,l,ierr
    real(8) :: avg_kap, n_kap

    !OPEN(unit=121,file='mods.txt',access='append')
    call get_all_curvatures(kap,2)
    x_mod=dxh((is+ie)/2); y_mod=dyh((js+je)/2); z_mod=dzh((ks+ke)/2) !assumes an unstretched grid
    P_gx = 0d0; P_gy = 0d0; P_gz = 0d0
    do k=ks,ke; do j=js,je; do i=is,ie
       Source = 0.d0
       if ( implode_flag(bub_id(i,j,k)) .and. vof_phase(i,j,k) == 1) Source = v_source(i,j,k)
       coeff(i,j,k,1) = dt/(dx(i)*dx(i)*rho)
       coeff(i,j,k,2) = dt/(dx(i)*dx(i)*rho)
       coeff(i,j,k,3) = dt/(dy(j)*dy(j)*rho)
       coeff(i,j,k,4) = dt/(dy(j)*dy(j)*rho)
       coeff(i,j,k,5) = dt/(dz(k)*dz(k)*rho)
       coeff(i,j,k,6) = dt/(dz(k)*dz(k)*rho)
       coeff(i,j,k,7) = sum(coeff(i,j,k,1:6))
       coeff(i,j,k,8) =  -(Source + (utmp(i,j,k)-utmp(i-1,j,k))/dx(i) &
            +  (vtmp(i,j,k)-vtmp(i,j-1,k))/dy(j) &
            +  (wtmp(i,j,k)-wtmp(i,j,k-1))/dz(k) )
       ! Debugging A8 NaN
       if (coeff(i,j,k,8) /= coeff(i,j,k,8)) then
          write(*,'("A8 NaN, liq_gas at x y z: ",3e14.5)')x(i),y(j),z(k)
          write(*,'("Limits: ",6I8)')is,ie,js,je,ks,ke
          write(*,'("ijk rank",4I4)')i,j,k,rank
          write(*,'("A branches: ",8e14.5)')coeff(i,j,k,:)
          write(*,'("Cvof 1-7: ",7e14.5)')cvof(i-1,j,k),cvof(i+1,j,k),cvof(i,j-1,k),cvof(i,j+1,k),&
               cvof(i,j,k-1),cvof(i,j,k+1),cvof(i,j,k)
          write(*,'("Phase 1-7: ",7I8)')vof_phase(i-1,j,k),vof_phase(i+1,j,k),vof_phase(i,j-1,k),vof_phase(i,j+1,k),&
               vof_phase(i,j,k-1),vof_phase(i,j,k+1),vof_phase(i,j,k)
          write(*,'("Pcmask 1-7: ",7I8)')pcmask(i-1,j,k),pcmask(i+1,j,k),pcmask(i,j-1,k),pcmask(i,j+1,k),&
               pcmask(i,j,k-1),pcmask(i,j,k+1),pcmask(i,j,k)
          write(*,'("S_v 1-7: ",7e14.5)')v_source(i-1,j,k),v_source(i+1,j,k),v_source(i,j-1,k),v_source(i,j+1,k),&
               v_source(i,j,k-1),v_source(i,j,k+1),v_source(i,j,k)
          write(*,'("Velocities in div u: ",6e14.5)')utmp(i,j,k),utmp(i-1,j,k),vtmp(i,j,k),vtmp(i,j-1,k),&
               wtmp(i,j,k),wtmp(i,j,k-1)
       endif
       !===============================================================================================================
       !----Cav-liquid neighbours, set P_g in cavity cells
       if (.not.implode_flag(bub_id(i,j,k))) then
          ! Set Laplace jumps for surface tension
          if(vof_phase(i,j,k)==1) then
             if (vof_phase(i+1,j,k)==0) then
                x_mod(i,j,k)=-1.d0*height(i+1,j,k,1)*dx(i)
                !if (x_mod(i,j,k) /= x_mod(i,j,k)) call mod_details(x_mod(i,j,k),i,j,k,.true.,1)
                if (x_mod(i,j,k)>dx(i) .or. x_mod(i,j,k)<limit*dxh(i))&
                     call staggered_cut(cvof,i,j,k,x_mod(i,j,k),vof_phase(i,j,k),1)
                if (x_mod(i,j,k) /= x_mod(i,j,k)) call mod_details(x_mod(i,j,k),i,j,k,.false.,1)
                !write(121,10)x(i+1),y(j),z(k),-x_mod(i,j,k),0d0,0d0
             endif
             if (vof_phase(i,j+1,k)==0) then
                y_mod(i,j,k)=-1.d0*height(i,j+1,k,3)*dy(j)
                !if (y_mod(i,j,k) /= y_mod(i,j,k)) call mod_details(y_mod(i,j,k),i,j,k,.true.,3)
                if (y_mod(i,j,k)>dy(j) .or. y_mod(i,j,k)<limit*dxh(i))&
                     call staggered_cut(cvof,i,j,k,y_mod(i,j,k),vof_phase(i,j,k),2)
                if (y_mod(i,j,k) /= y_mod(i,j,k)) call mod_details(y_mod(i,j,k),i,j,k,.false.,3)
                !write(121,10)x(i),y(j+1),z(k),0d0,-y_mod(i,j,k),0d0
             endif
             if (vof_phase(i,j,k+1)==0) then
                z_mod(i,j,k)=-1.d0*height(i,j,k+1,5)*dz(k)
                !if (z_mod(i,j,k) /= z_mod(i,j,k)) call mod_details(z_mod(i,j,k),i,j,k,.true.,5)
                if (z_mod(i,j,k)>dz(k) .or. z_mod(i,j,k)<limit*dxh(i))&
                     call staggered_cut(cvof,i,j,k,z_mod(i,j,k),vof_phase(i,j,k),3)
                if (z_mod(i,j,k) /= z_mod(i,j,k)) call mod_details(z_mod(i,j,k),i,j,k,.false.,5)
                !write(121,10)x(i),y(j),z(k+1),0d0,0d0,-z_mod(i,j,k)
             endif
          endif ! Cavity cell

          if(vof_phase(i,j,k)==0) then
             if (vof_phase(i+1,j,k)==1) then
                x_mod(i,j,k)=height(i,j,k,2)*dx(i)
                !if (x_mod(i,j,k) /= x_mod(i,j,k)) call mod_details(x_mod(i,j,k),i,j,k,.true.,2)
                if (x_mod(i,j,k)>dx(i) .or. x_mod(i,j,k)<limit*dxh(i))&
                     call staggered_cut(cvof,i,j,k,x_mod(i,j,k),vof_phase(i,j,k),1)
                if (x_mod(i,j,k) /= x_mod(i,j,k)) call mod_details(x_mod(i,j,k),i,j,k,.false.,2)
                !write(121,10)x(i),y(j),z(k),x_mod(i,j,k),0d0,0d0
             endif
             if (vof_phase(i,j+1,k)==1) then
                y_mod(i,j,k)=height(i,j,k,4)*dy(j)
                !if (y_mod(i,j,k) /= y_mod(i,j,k)) call mod_details(y_mod(i,j,k),i,j,k,.true.,4)
                if (y_mod(i,j,k)>dy(j) .or. y_mod(i,j,k)<limit*dxh(i))&
                     call staggered_cut(cvof,i,j,k,y_mod(i,j,k),vof_phase(i,j,k),2)
                if (y_mod(i,j,k) /= y_mod(i,j,k)) call mod_details(y_mod(i,j,k),i,j,k,.false.,4)
                !write(121,10)x(i),y(j),z(k),0d0,y_mod(i,j,k),0d0
             endif
             if (vof_phase(i,j,k+1)==1) then
                z_mod(i,j,k)=height(i,j,k,6)*dz(k)
                !if (z_mod(i,j,k) /= z_mod(i,j,k)) call mod_details(z_mod(i,j,k),i,j,k,.true.,6)
                if (z_mod(i,j,k)>dz(k) .or. z_mod(i,j,k)<limit*dxh(i))&
                     call staggered_cut(cvof,i,j,k,z_mod(i,j,k),vof_phase(i,j,k),3)
                if (z_mod(i,j,k) /= z_mod(i,j,k)) call mod_details(z_mod(i,j,k),i,j,k,.false.,6)
                !write(121,10)x(i),y(j),z(k),0d0,0d0,z_mod(i,j,k)
             endif
          endif ! Liquid cell
       endif ! we are not imploding
    enddo; enddo; enddo
10  format(6e14.5)
    !close(121)
    call ghost_x(x_mod,1,reqd(1:4)); call ghost_y(y_mod,1,reqd(5:8)); call ghost_z(z_mod,1,reqd(9:12))
    call MPI_WAITALL(12,reqd(1:12),stat(:,1:12),ierr)

    do k=ks,ke; do j=js,je; do i=is,ie
       if ( (.not.implode_flag(bub_id(i,j,k))) .and. vof_phase(i,j,k)==1) then
          do l=-1,1,2
             avg_kap=0.0d0; n_kap=0.d0
             if (vof_phase(i+l,j,k)==0) then
                if (vof_flag(i+l,j,k)==2 .and. kap(i+l,j,k)<1.d6) then
                   n_kap=n_kap+1.d0
                   avg_kap = avg_kap + kap(i+l,j,k)
                endif
                if (vof_flag(i,j,k)==2 .and. kap(i,j,k)<1.d6 ) then
                   n_kap=n_kap+1.d0
                   avg_kap=avg_kap + kap(i,j,k)
                endif
                if (n_kap>=9.99d-1) then
                   avg_kap=avg_kap/n_kap
                else
                   call pariserror('No curvature found in liq-gas pair for FS bubble')
                endif
                P_gx(i,j,k) = sigma*avg_kap/dx(i) !!filaments and droplets of one cell will be an issue here
             endif

             avg_kap=0.0d0; n_kap=0.d0
             if (vof_phase(i,j+l,k)==0) then
                if (vof_flag(i,j+l,k)==2 .and. kap(i,j+l,k)<1.d6) then
                   n_kap=n_kap+1.d0
                   avg_kap = avg_kap + kap(i,j+l,k)
                endif
                if (vof_flag(i,j,k)==2 .and. kap(i,j,k)<1.d6 ) then
                   n_kap=n_kap+1.d0
                   avg_kap=avg_kap + kap(i,j,k)
                endif
                if (n_kap>=9.99d-1) then
                   avg_kap=avg_kap/n_kap
                else
                   call pariserror('No curvature found in liq-gas pair for FS bubble')
                endif
                P_gy(i,j,k) = sigma*avg_kap/dy(j)
             endif

             if (vof_phase(i,j,k+l)==0) then
                if (vof_flag(i,j,k+l)==2 .and. kap(i,j,k+l)<1.d6) then
                   n_kap=n_kap+1.d0
                   avg_kap = avg_kap + kap(i,j,k+l)
                endif
                if (vof_flag(i,j,k)==2 .and. kap(i,j,k)<1.d6 ) then
                   n_kap=n_kap+1.d0
                   avg_kap=avg_kap + kap(i,j,k)
                endif
                if (n_kap>=9.99d-1) then
                   avg_kap=avg_kap/n_kap
                else
                   call pariserror('No curvature found in liq-gas pair for FS bubble')
                endif
                P_gz(i,j,k) = sigma*avg_kap/dz(k)
             endif
          enddo
       endif
    enddo; enddo; enddo
    call ghost_x(P_gx,1,reqd(1:4)); call ghost_y(P_gy,1,reqd(5:8)); call ghost_z(P_gz,1,reqd(9:12))
    call MPI_WAITALL(12,reqd(1:12),stat(:,1:12),ierr)
    !--------------------------------------------------------------------------------------------------------
    do k=ks,ke; do j=js,je; do i=is,ie
       if (vof_phase(i,j,k)==0 .and. (.not.implode_flag(bub_id(i,j,k))) ) then

          coeff(i,j,k,1) = 2.d0*dt/((dx(i)+x_mod(i-1,j,k))*x_mod(i-1,j,k)*rho)
          coeff(i,j,k,2) = 2.d0*dt/((dx(i)+x_mod(i  ,j,k))*x_mod(i  ,j,k)*rho)
          coeff(i,j,k,3) = 2.d0*dt/((dy(j)+y_mod(i,j-1,k))*y_mod(i,j-1,k)*rho)
          coeff(i,j,k,4) = 2.d0*dt/((dy(j)+y_mod(i,j  ,k))*y_mod(i,j  ,k)*rho)
          coeff(i,j,k,5) = 2.d0*dt/((dz(k)+z_mod(i,j,k-1))*z_mod(i,j,k-1)*rho)
          coeff(i,j,k,6) = 2.d0*dt/((dz(k)+z_mod(i,j,k  ))*z_mod(i,j,k  )*rho)
          coeff(i,j,k,7) = sum(coeff(i,j,k,1:6))
          coeff(i,j,k,8) = coeff(i,j,k,8) + coeff(i,j,k,1)*P_gx(i-1,j,k) + coeff(i,j,k,2)*P_gx(i+1,j,k)&
               +coeff(i,j,k,3)*P_gy(i,j-1,k)+coeff(i,j,k,4)*P_gy(i,j+1,k)&
               +coeff(i,j,k,5)*P_gz(i,j,k-1)+coeff(i,j,k,6)*P_gz(i,j,k+1)
          if (coeff(i,j,k,8) /= coeff(i,j,k,8)) then
             write(*,'("A8 NaN, error imminent. Neigbours mods :",6e14.5)')x_mod(i-1,j,k),x_mod(i,j,k),&
                  y_mod(i,j-1,k),y_mod(i,j,k),z_mod(i,j,k-1),z_mod(i,j,k)
             write(*,'("A8 NaN, error imminent. P_g :",6e14.5,2I8)')P_gx(i-1,j,k),P_gx(i+1,j,k),&
                  P_gy(i,j-1,k),P_gy(i,j+1,k),P_gz(i,j,k-1),P_gz(i,j,k+1),i,k
          endif
          if (pcmask(i,j,k).ne.0) write(*,'("Error topology, phase 0, pcmask :",i8)')pcmask(i,j,k)
       endif
    enddo;enddo;enddo

    if (.not. RP_test) call Poisson_BCs(coeff)
  end subroutine liq_gas2
end subroutine setuppoisson_fs_heights
!==============================================================================================================
!calculates theta (height) when HF is not available (by using a staggered cell and fl3d).
subroutine staggered_cut(cvof,i,j,k,theta,phase,d)
  use module_grid
  use module_freesurface
  implicit none
  real(8) :: theta
  integer :: i,j,k,phase,d
  integer :: i0,j0,k0,loc(3)
  real(8), dimension(-1:1,-1:1,-1:1) :: stencil3x3
  real(8) :: n_ref(3), n_nbr(3), x0(3), dc(3), n_stag(3), nr(3)
  real(8) :: alpha2, c_ref, c_nbr, c_stag, c_min=1.0d-1, test
  real(8) :: al3d, fl3d
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: cvof

  if (.not.(phase==1 .or. phase==0)) call pariserror('cell phase has to be 0 or 1 in staggered_cut call')

  loc = 0
  loc(d) = 1
  theta = dxh(is)
  !liq_gas, make general
  !get normals
  do i0=-1,1; do j0=-1,1; do k0=-1,1
     stencil3x3(i0,j0,k0) = cvof(i+i0,j+j0,k+k0)
  enddo;enddo;enddo
  call mycs(stencil3x3,n_ref)
  do i0=-1,1; do j0=-1,1; do k0=-1,1
     stencil3x3(i0,j0,k0) = cvof(i+loc(1)+i0,j+loc(2)+j0,k+loc(3)+k0)
  enddo;enddo;enddo
  call mycs(stencil3x3,n_nbr)
  !=====================
  alpha2 = al3d(n_ref,cvof(i,j,k))
  x0 = 0.0d0; x0(d) = 0.5d0
  dc = 1.0d0; dc(d) = 0.5d0
  c_ref = fl3d(n_ref,alpha2,x0,dc)

  alpha2 = al3d(n_nbr,cvof(i+loc(1),j+loc(2),k+loc(3)))
  x0=0.0d0
  dc = 1.0d0; dc(d) = 0.5d0
  c_nbr = fl3d(n_nbr,alpha2,x0,dc)
  c_stag = c_ref+c_nbr

  if (min(c_ref,c_nbr)>c_min) then
     nr = 0.5d0*(n_ref+n_nbr)
     n_stag(1) = nr(1)/(ABS(nr(1))+ABS(nr(2))+ABS(nr(3)))
     n_stag(2) = nr(2)/(ABS(nr(1))+ABS(nr(2))+ABS(nr(3)))
     n_stag(3) = nr(3)/(ABS(nr(1))+ABS(nr(2))+ABS(nr(3)))
  else
     if (phase==1) then
        n_stag=n_ref
     else
        n_stag=n_nbr
     endif
  endif

  if (ABS((ABS(n_stag(1))+ABS(n_stag(2))+ABS(n_stag(3)))-1.0d0) > 1.0d-12) then
     write(*,*)'Normals not normalised'
     write(*,'("Normals: ",3e14.5)')n_stag(1:3)
     call mod_details(-1.0d2,i,j,k,.false.,d)
     call pariserror('Normals not normalised')
  endif

  alpha2=al3d(n_stag,c_stag)
  if (ABS(n_stag(d))>1.0d-12) then
     test = (alpha2 - (n_stag(1)+n_stag(2)+n_stag(3)-n_stag(d))/2.0d0)/n_stag(d)
     if (phase==1) then
        theta = dxh(i)*(1.0d0-test)
     else
        theta = dxh(i)*test
     endif
     if (theta>dxh(i)) theta = dxh(i)
     if (theta<limit*dxh(i)) theta = limit*dxh(i)
  else
     theta = dxh(i) !set to standard length if no cut in staggered cell can be found
  endif
end subroutine staggered_cut
!================================================================================================================
subroutine mod_details(h,i,j,k,height,d)
  use module_grid
  use module_VOF
  use module_freesurface
  use module_flow
  implicit none
  real(8) :: h
  integer :: i,j,k,d
  logical :: height
  integer :: ind(1:3)

  ind = 0
  open(unit=70,file="mod_details.txt",position="append")
  write(70,'("Mod in direction ",I4,": ",e14.5)')d,h
  write(70,'("Time step ",I10)')itimestep
  if (height) then
     write(70,'("NaN using Heights")')
  else
     write(70,'("NaN using staggered cut")')
  endif
  write(70,'("Limits: ",6I8)')is,ie,js,je,ks,ke
  write(70,'("ijk rank",4I4)')i,j,k,rank
  write(70,'("x, y, z: ",3e14.5)')x(i),y(j),z(k)
  write(70,'("Cvof 1-7: ",7e14.5)')cvof(i-1,j,k),cvof(i+1,j,k),cvof(i,j-1,k),cvof(i,j+1,k),&
       cvof(i,j,k-1),cvof(i,j,k+1),cvof(i,j,k)
  write(70,'("Phase 1-7: ",7I8)')vof_phase(i-1,j,k),vof_phase(i+1,j,k),vof_phase(i,j-1,k),vof_phase(i,j+1,k),&
       vof_phase(i,j,k-1),vof_phase(i,j,k+1),vof_phase(i,j,k)
  write(70,'("Pcmask 1-7: ",7I8)')pcmask(i-1,j,k),pcmask(i+1,j,k),pcmask(i,j-1,k),pcmask(i,j+1,k),&
       pcmask(i,j,k-1),pcmask(i,j,k+1),pcmask(i,j,k)
  write(70,'("S_v 1-7: ",7e14.5)')v_source(i-1,j,k),v_source(i+1,j,k),v_source(i,j-1,k),v_source(i,j+1,k),&
       v_source(i,j,k-1),v_source(i,j,k+1),v_source(i,j,k)
  write(70,'("Velocities in div u: ",6e14.5)')u(i,j,k),u(i-1,j,k),v(i,j,k),v(i,j-1,k),&
       w(i,j,k),w(i,j,k-1)
  write(70,'("  ")')
  close(70)
end subroutine mod_details
!=============================================================================================================================
subroutine curvature_sphere(t)
  use module_VOF
  use module_surface_tension
  use module_grid
  use module_freesurface
  implicit none
  include 'mpif.h'
  real(8) :: kap(imin:imax,jmin:jmax,kmin:kmax)
  integer :: i,j,k,l,ierr
  real(8), dimension(1:2) :: errnorm(2), err_glob(2)
  real(8) :: kap_sphere, max_loc, err_inf, gh_cells, t_cells, t
  real(8) :: V_loc, V_t
  logical :: check
  ! initialize errors
  errnorm = 0.0d0; err_glob=0.0d0; gh_cells = 0.0d0; max_loc = 0.0d0

  OPEN(unit=122,file='curve_stats',position='append')
  ! Get bubble volume
  V_loc = 0.0d0
  do k=ks,ke; do j=js,je; do i=is,ie
     V_loc=V_loc+cvof(i,j,k)*dx(i)*dy(j)*dz(k)
  enddo; enddo; enddo
  call MPI_ALLREDUCE(V_loc, V_t, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_Domain, ierr)
  ! Curvature for perfect sphere
  kap_sphere = -1.0d0*((32.0d0*pi)/(V_t*3.0d0))**(1.0/3.0)
  !write(*,'("Curvature sphere: ",e14.5)')kap_sphere
  ! Get calculated curvature
  call get_all_curvatures(kap,2)
  do k=ks,ke; do j=js,je; do i=is,ie
     if(vof_phase(i,j,k)==1) then
        check = .false.
        do l=-1,1,2
           if (vof_phase(i+l,j,k)==0) check = .true.
           if (vof_phase(i,j+l,k)==0) check = .true.
           if (vof_phase(i,j,k+l)==0) check = .true.
        enddo
        if (check) then
           call err_curve(kap(i,j,k)/dx(i),kap_sphere,errnorm,max_loc)
           gh_cells=gh_cells+1.0d0
        endif
     endif
  enddo; enddo; enddo
  call MPI_ALLREDUCE(errnorm(1), err_glob(1), 2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_Domain, ierr)
  call MPI_ALLREDUCE(max_loc, err_inf, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_Comm_Cart, ierr)
  call MPI_ALLREDUCE(gh_cells, t_cells, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr)
  if (t_cells<1.0d0) then
     write(*,*)"WARNING: curvature calc without cut cells!"
  else
     err_glob(1)=err_glob(1)/t_cells
     err_glob(2)=sqrt(err_glob(2))/t_cells
  endif
  if (rank==0) write(122,13)t,err_glob(1),err_glob(2),err_inf

13 format(4e14.5)
  CLOSE(122)
contains
  subroutine err_curve(kapc,ref,err_out,localmax)
    implicit none
    real(8) :: kapc, ref, localmax
    real(8) :: err_out(1:2)
    err_out(1)=err_out(1)+ABS(kapc-ref)
    err_out(2)=err_out(2)+((kapc-ref)**2.0d0)
    localmax=MAX(ABS(kapc-ref),localmax)
  end subroutine err_curve
end subroutine curvature_sphere
!=======================================================================================================
subroutine remove_drops(phase,vof,fill_all) !some cleaning
  use module_grid
  use module_freesurface
  implicit none
  include 'mpif.h'
  integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: phase
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: vof
  logical :: clearcells, fill, fill_all
  integer :: i,j,k,ierr
  integer :: inr,jnr,knr,phase_ref

  fill = .false.
  do i=is,ie; do j=js,je; do k=ks,ke
     clearcells = .true.
     !Check "shell" at 5 cubed
     phase_ref=phase(i-2,j-2,k-2)
     do inr=-2,2,4; do jnr=-2,2; do knr=-2,2
        if (.not.(phase(i+inr,j+jnr,k+knr)==phase_ref)) then
           clearcells = .false.
           exit
        endif
     enddo; enddo; enddo
     do jnr=-2,2,4; do inr=-2,2; do knr=-2,2
        if (.not.(phase(i+inr,j+jnr,k+knr)==phase_ref)) then
           clearcells = .false.
           exit
        endif
     enddo; enddo; enddo
     do knr=-2,2,4; do inr=-2,2; do jnr=-2,2
        if (.not.(phase(i+inr,j+jnr,k+knr)==phase_ref)) then
           clearcells = .false.
           exit
        endif
     enddo; enddo; enddo
     !Set "core" of 3 cubed to ref phase
     if (clearcells) then
        do inr=-2,2; do jnr=-2,2; do knr=-2,2
           if (phase_ref==0) then
              vof(i+inr,j+jnr,k+knr)=0.00d0
           else
              vof(i+inr,j+jnr,k+knr)=1.00d0
           endif
        enddo; enddo; enddo
        fill =.true.
     endif
  enddo; enddo; enddo
  call MPI_ALLREDUCE(fill,fill_all, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_Active, ierr)
end subroutine remove_drops

subroutine setuppoisson_fs_hypre(utmp,vtmp,wtmp,rho,dt,coeff,height,kap,bub_id)
  use module_grid
  use module_BC
  use module_freesurface
  use module_IO
  use module_2phase
  use module_VOF
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax,6), intent(in) :: height
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: utmp,vtmp,wtmp
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: kap
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(out) :: coeff
  integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: bub_id
  real(8), intent(in) :: dt, rho
  integer :: i,j,k,nbr,ierr,l
  integer :: reqd(12),stat(MPI_STATUS_SIZE,12)
  real(8) :: Source
  real(8) :: avg_kap, n_kap
  !!Debugging for curvature errors
!!$  real(8) :: err, err_global, L2_err, Linf_err, kappa_theory, L2_glob, Linf_glob
!!$  integer :: uncomp_curv, uncomp_curv_global, n_avg_kap, n_glob

  !OPEN(unit=121,file='mods.txt',access='append')
  x_mod=dxh((is+ie)/2); y_mod=dyh((js+je)/2); z_mod=dzh((ks+ke)/2) !assumes an unstretched grid
  P_gx = 0d0; P_gy = 0d0; P_gz = 0d0
  do k=ks,ke; do j=js,je; do i=is,ie
     Source = 0.d0
     if (implode_flag(bub_id(i,j,k)) .and. vof_phase(i,j,k) == 1) Source = v_source(i,j,k)
     coeff(i,j,k,1) = dt/(dx(i)*dx(i)*rho)
     coeff(i,j,k,2) = dt/(dx(i)*dx(i)*rho)
     coeff(i,j,k,3) = dt/(dy(j)*dy(j)*rho)
     coeff(i,j,k,4) = dt/(dy(j)*dy(j)*rho)
     coeff(i,j,k,5) = dt/(dz(k)*dz(k)*rho)
     coeff(i,j,k,6) = dt/(dz(k)*dz(k)*rho)
     coeff(i,j,k,7) = sum(coeff(i,j,k,1:6))
     coeff(i,j,k,8) =  -(Source + (utmp(i,j,k)-utmp(i-1,j,k))/dx(i) &
          +  (vtmp(i,j,k)-vtmp(i,j-1,k))/dy(j) &
          +  (wtmp(i,j,k)-wtmp(i,j,k-1))/dz(k) )
     ! Debugging A8 NaN
     if (coeff(i,j,k,8) /= coeff(i,j,k,8)) then
        write(*,'("A8 NaN, liq_gas at x y z: ",3e14.5)')x(i),y(j),z(k)
        write(*,'("Limits: ",6I8)')is,ie,js,je,ks,ke
        write(*,'("ijk rank",4I4)')i,j,k,rank
        write(*,'("A branches: ",8e14.5)')coeff(i,j,k,:)
        write(*,'("Cvof 1-7: ",7e14.5)')cvof(i-1,j,k),cvof(i+1,j,k),cvof(i,j-1,k),cvof(i,j+1,k),&
             cvof(i,j,k-1),cvof(i,j,k+1),cvof(i,j,k)
        write(*,'("Phase 1-7: ",7I8)')vof_phase(i-1,j,k),vof_phase(i+1,j,k),vof_phase(i,j-1,k),vof_phase(i,j+1,k),&
             vof_phase(i,j,k-1),vof_phase(i,j,k+1),vof_phase(i,j,k)
        write(*,'("Pcmask 1-7: ",7I8)')pcmask(i-1,j,k),pcmask(i+1,j,k),pcmask(i,j-1,k),pcmask(i,j+1,k),&
             pcmask(i,j,k-1),pcmask(i,j,k+1),pcmask(i,j,k)
        write(*,'("S_v 1-7: ",7e14.5)')v_source(i-1,j,k),v_source(i+1,j,k),v_source(i,j-1,k),v_source(i,j+1,k),&
             v_source(i,j,k-1),v_source(i,j,k+1),v_source(i,j,k)
        write(*,'("Velocities in div u: ",6e14.5)')utmp(i,j,k),utmp(i-1,j,k),vtmp(i,j,k),vtmp(i,j-1,k),&
             wtmp(i,j,k),wtmp(i,j,k-1)
     endif
     !===============================================================================================================
     !----Cav-liquid neighbours, set P_g in cavity cells
     !if (.not.implode_flag(bub_id(i,j,k))) then
        ! Set Laplace jumps for surface tension
        if(vof_phase(i,j,k)==1) then
           if (vof_phase(i+1,j,k)==0) then
              x_mod(i,j,k)=-1.d0*height(i+1,j,k,1)*dx(i)
              !if (x_mod(i,j,k) /= x_mod(i,j,k)) call mod_details(x_mod(i,j,k),i,j,k,.true.,1)
              if (x_mod(i,j,k)>dx(i) .or. x_mod(i,j,k)<limit*dxh(i))&
                   call staggered_cut(cvof,i,j,k,x_mod(i,j,k),vof_phase(i,j,k),1)
              if (x_mod(i,j,k) /= x_mod(i,j,k)) call mod_details(x_mod(i,j,k),i,j,k,.false.,1)
              !write(121,10)x(i+1),y(j),z(k),-x_mod(i,j,k),0d0,0d0
           endif
           if (vof_phase(i,j+1,k)==0) then
              y_mod(i,j,k)=-1.d0*height(i,j+1,k,3)*dy(j)
              !if (y_mod(i,j,k) /= y_mod(i,j,k)) call mod_details(y_mod(i,j,k),i,j,k,.true.,3)
              if (y_mod(i,j,k)>dy(j) .or. y_mod(i,j,k)<limit*dxh(i))&
                   call staggered_cut(cvof,i,j,k,y_mod(i,j,k),vof_phase(i,j,k),2)
              if (y_mod(i,j,k) /= y_mod(i,j,k)) call mod_details(y_mod(i,j,k),i,j,k,.false.,3)
              !write(121,10)x(i),y(j+1),z(k),0d0,-y_mod(i,j,k),0d0
           endif
           if (vof_phase(i,j,k+1)==0) then
              z_mod(i,j,k)=-1.d0*height(i,j,k+1,5)*dz(k)
              !if (z_mod(i,j,k) /= z_mod(i,j,k)) call mod_details(z_mod(i,j,k),i,j,k,.true.,5)
              if (z_mod(i,j,k)>dz(k) .or. z_mod(i,j,k)<limit*dxh(i))&
                   call staggered_cut(cvof,i,j,k,z_mod(i,j,k),vof_phase(i,j,k),3)
              if (z_mod(i,j,k) /= z_mod(i,j,k)) call mod_details(z_mod(i,j,k),i,j,k,.false.,5)
              !write(121,10)x(i),y(j),z(k+1),0d0,0d0,-z_mod(i,j,k)
           endif
        endif ! Cavity cell

        if(vof_phase(i,j,k)==0) then
           if (vof_phase(i+1,j,k)==1) then
              x_mod(i,j,k)=height(i,j,k,2)*dx(i)
              !if (x_mod(i,j,k) /= x_mod(i,j,k)) call mod_details(x_mod(i,j,k),i,j,k,.true.,2)
              if (x_mod(i,j,k)>dx(i) .or. x_mod(i,j,k)<limit*dxh(i))&
                   call staggered_cut(cvof,i,j,k,x_mod(i,j,k),vof_phase(i,j,k),1)
              if (x_mod(i,j,k) /= x_mod(i,j,k)) call mod_details(x_mod(i,j,k),i,j,k,.false.,2)
              !write(121,10)x(i),y(j),z(k),x_mod(i,j,k),0d0,0d0
           endif
           if (vof_phase(i,j+1,k)==1) then
              y_mod(i,j,k)=height(i,j,k,4)*dy(j)
              !if (y_mod(i,j,k) /= y_mod(i,j,k)) call mod_details(y_mod(i,j,k),i,j,k,.true.,4)
              if (y_mod(i,j,k)>dy(j) .or. y_mod(i,j,k)<limit*dxh(i))&
                   call staggered_cut(cvof,i,j,k,y_mod(i,j,k),vof_phase(i,j,k),2)
              if (y_mod(i,j,k) /= y_mod(i,j,k)) call mod_details(y_mod(i,j,k),i,j,k,.false.,4)
              !write(121,10)x(i),y(j),z(k),0d0,y_mod(i,j,k),0d0
           endif
           if (vof_phase(i,j,k+1)==1) then
              z_mod(i,j,k)=height(i,j,k,6)*dz(k)
              !if (z_mod(i,j,k) /= z_mod(i,j,k)) call mod_details(z_mod(i,j,k),i,j,k,.true.,6)
              if (z_mod(i,j,k)>dz(k) .or. z_mod(i,j,k)<limit*dxh(i))&
                   call staggered_cut(cvof,i,j,k,z_mod(i,j,k),vof_phase(i,j,k),3)
              if (z_mod(i,j,k) /= z_mod(i,j,k)) call mod_details(z_mod(i,j,k),i,j,k,.false.,6)
              !write(121,10)x(i),y(j),z(k),0d0,0d0,z_mod(i,j,k)
           endif
        endif ! Liquid cell
     !endif ! we are not imploding
  enddo; enddo; enddo
10 format(6e14.5)
  !close(121)
  call ghost_x(x_mod,1,reqd(1:4)); call ghost_y(y_mod,1,reqd(5:8)); call ghost_z(z_mod,1,reqd(9:12))
  call MPI_WAITALL(12,reqd(1:12),stat(:,1:12),ierr)

!!$  uncomp_curv=0; n_avg_kap=0; L2_err=0.0d0; Linf_err=0.0d0
!!$  kappa_theory = -2.0d0*dx(is)/r_min
  do k=ks,ke; do j=js,je; do i=is,ie
     !if( .not.implode_flag(bub_id(i,j,k)) .and. vof_phase(i,j,k)==1) then
     if(vof_phase(i,j,k)==1) then
        do l=-1,1,2
           avg_kap=0.0d0; n_kap=0.d0
           if (vof_phase(i+l,j,k)==0) then
              if (vof_flag(i+l,j,k)==2 .and. kap(i+l,j,k)<1.d6) then
                 n_kap=n_kap+1.d0
                 avg_kap = avg_kap + kap(i+l,j,k)
              endif
              if (vof_flag(i,j,k)==2 .and. kap(i,j,k)<1.d6 ) then
                 n_kap=n_kap+1.d0
                 avg_kap=avg_kap + kap(i,j,k)
              endif
              if (n_kap>=9.99d-1) then
                 avg_kap=avg_kap/n_kap
              else
                 call pariserror('No curvature found in liq-gas pair for FS bubble')
              endif
!!$              n_avg_kap=n_avg_kap+1
!!$              err=ABS(avg_kap-kappa_theory)
!!$              L2_err=L2_err+err**2.0d0
!!$              Linf_err=MAX(err,Linf_err)
              P_gx(i,j,k) = sigma*avg_kap/dx(i) !!filaments and droplets of one cell will be an issue here
           endif

           avg_kap=0.0d0; n_kap=0.d0
           if (vof_phase(i,j+l,k)==0) then
              if (vof_flag(i,j+l,k)==2 .and. kap(i,j+l,k)<1.d6) then
                 n_kap=n_kap+1.d0
                 avg_kap = avg_kap + kap(i,j+l,k)
              endif
              if (vof_flag(i,j,k)==2 .and. kap(i,j,k)<1.d6 ) then
                 n_kap=n_kap+1.d0
                 avg_kap=avg_kap + kap(i,j,k)
              endif
              if (n_kap>=9.99d-1) then
                 avg_kap=avg_kap/n_kap
              else
                 call pariserror('No curvature found in liq-gas pair for FS bubble')
              endif
!!$              n_avg_kap=n_avg_kap+1
!!$              err=ABS(avg_kap-kappa_theory)
!!$              L2_err=L2_err+err**2.0d0
!!$              Linf_err=MAX(err,Linf_err)
              P_gy(i,j,k) = sigma*avg_kap/dy(j)
           endif

           if (vof_phase(i,j,k+l)==0) then
              if (vof_flag(i,j,k+l)==2 .and. kap(i,j,k+l)<1.d6) then
                 n_kap=n_kap+1.d0
                 avg_kap = avg_kap + kap(i,j,k+l)
              endif
              if (vof_flag(i,j,k)==2 .and. kap(i,j,k)<1.d6 ) then
                 n_kap=n_kap+1.d0
                 avg_kap=avg_kap + kap(i,j,k)
              endif
              if (n_kap>=9.99d-1) then
                 avg_kap=avg_kap/n_kap
              else
                 call pariserror('No curvature found in liq-gas pair for FS bubble')
              endif
!!$              n_avg_kap=n_avg_kap+1
!!$              err=ABS(avg_kap-kappa_theory)
!!$              L2_err=L2_err+err**2.0d0
!!$              Linf_err=MAX(err,Linf_err)
              P_gz(i,j,k) = sigma*avg_kap/dz(k)
           endif
        enddo
     endif
  enddo; enddo; enddo
  call ghost_x(P_gx,1,reqd(1:4)); call ghost_y(P_gy,1,reqd(5:8)); call ghost_z(P_gz,1,reqd(9:12))
  call MPI_WAITALL(12,reqd(1:12),stat(:,1:12),ierr)
  !! Debugging uncomputed curvatures
!!$  call MPI_ALLREDUCE(uncomp_curv,uncomp_curv_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_Active, ierr)
!!$  call MPI_ALLREDUCE(n_avg_kap,n_glob, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_Active, ierr)
!!$  call MPI_ALLREDUCE(L2_err,L2_glob, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_Active, ierr)
!!$  call MPI_ALLREDUCE(Linf_err,Linf_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_Active, ierr)
!!$  L2_glob = sqrt(L2_glob)/n_glob
!!$  if (rank==0) then
!!$     write(*,'("Number of average curvatures reqd , uncomputed: ",2I8)')n_glob, uncomp_curv_global
!!$     write(*,'("Curvature error norms L2: ",e14.5," L_inf: ",e14.5)')L2_glob,Linf_glob
!!$     write(*,'(" ")')
!!$  endif
  !--------------------------------------------------------------------------------------------------------
  do k=ks,ke; do j=js,je; do i=is,ie
     !if (vof_phase(i,j,k)==0 .and. .not.implode_flag(bub_id(i,j,k))) then
     if (vof_phase(i,j,k)==0) then
        coeff(i,j,k,1) = 2.d0*dt/((dx(i)+x_mod(i-1,j,k))*x_mod(i-1,j,k)*rho)
        coeff(i,j,k,2) = 2.d0*dt/((dx(i)+x_mod(i  ,j,k))*x_mod(i  ,j,k)*rho)
        coeff(i,j,k,3) = 2.d0*dt/((dy(j)+y_mod(i,j-1,k))*y_mod(i,j-1,k)*rho)
        coeff(i,j,k,4) = 2.d0*dt/((dy(j)+y_mod(i,j  ,k))*y_mod(i,j  ,k)*rho)
        coeff(i,j,k,5) = 2.d0*dt/((dz(k)+z_mod(i,j,k-1))*z_mod(i,j,k-1)*rho)
        coeff(i,j,k,6) = 2.d0*dt/((dz(k)+z_mod(i,j,k  ))*z_mod(i,j,k  )*rho)
        coeff(i,j,k,7) = sum(coeff(i,j,k,1:6))
        coeff(i,j,k,8) = coeff(i,j,k,8) + coeff(i,j,k,1)*P_gx(i-1,j,k) + coeff(i,j,k,2)*P_gx(i+1,j,k)&
             +coeff(i,j,k,3)*P_gy(i,j-1,k)+coeff(i,j,k,4)*P_gy(i,j+1,k)&
             +coeff(i,j,k,5)*P_gz(i,j,k-1)+coeff(i,j,k,6)*P_gz(i,j,k+1)
        if (coeff(i,j,k,8) /= coeff(i,j,k,8)) then
           write(*,'("A8 NaN, error imminent. Neigbours mods :",6e14.5)')x_mod(i-1,j,k),x_mod(i,j,k),&
                y_mod(i,j-1,k),y_mod(i,j,k),z_mod(i,j,k-1),z_mod(i,j,k)
           write(*,'("A8 NaN, error imminent. P_g :",6e14.5,2I8)')P_gx(i-1,j,k),P_gx(i+1,j,k),&
                P_gy(i,j-1,k),P_gy(i,j+1,k),P_gz(i,j,k-1),P_gz(i,j,k+1),i,k
        endif
        if (pcmask(i,j,k).ne.0) write(*,'("Error topology, phase 0, pcmask :",i8)')pcmask(i,j,k)
     else if (vof_phase(i,j,k)==1) then
        coeff(i,j,k,1:6) = 0.0d0
        coeff(i,j,k,7) = 1.0d0
        coeff(i,j,k,8) = P_gas(i,j,k) !Should now work for polytropic law and vacuum
     endif
  enddo;enddo;enddo

  call Poisson_BCs(coeff)
end subroutine setuppoisson_fs_hypre
