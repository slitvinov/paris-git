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
  integer :: reqs(48),stat(MPI_STATUS_SIZE,48)
  integer :: i,j,k,level,iout,ierr1
  integer :: level3, l3sum

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
     if (vof_phase(i,j,k) == 0 .or. implode(i,j,k)>0) then
        pcmask(i,j,k)=0
        u_cmask(i,j,k)=0
        v_cmask(i,j,k)=0
        w_cmask(i,j,k)=0
     endif
  enddo; enddo; enddo
  !fill ghost layers for zero masks
  call ighost_x(u_cmask,2,reqs(1:4)); call ighost_x(v_cmask,2,reqs(5:8))
  call ighost_x(w_cmask,2,reqs(9:12)); call ighost_x(pcmask,2,reqs(13:16))
  call ighost_y(u_cmask,2,reqs(17:20)); call ighost_y(v_cmask,2,reqs(21:24))
  call ighost_y(w_cmask,2,reqs(25:28)); call ighost_y(pcmask,2,reqs(29:32))
  call ighost_z(u_cmask,2,reqs(33:36)); call ighost_z(v_cmask,2,reqs(37:40))
  call ighost_z(w_cmask,2,reqs(41:44)); call ighost_z(pcmask,2,reqs(45:48))
  call MPI_WAITALL(48,reqs(1:48),stat,ierr1)
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
     call ighost_x(u_cmask,2,reqs(1:4)); call ighost_x(v_cmask,2,reqs(5:8))
     call ighost_x(w_cmask,2,reqs(9:12)); call ighost_x(pcmask,2,reqs(13:16))
     call ighost_y(u_cmask,2,reqs(17:20)); call ighost_y(v_cmask,2,reqs(21:24))
     call ighost_y(w_cmask,2,reqs(25:28)); call ighost_y(pcmask,2,reqs(29:32))
     call ighost_z(u_cmask,2,reqs(33:36)); call ighost_z(v_cmask,2,reqs(37:40))
     call ighost_z(w_cmask,2,reqs(41:44)); call ighost_z(pcmask,2,reqs(45:48))
     call MPI_WAITALL(48,reqs(1:48),stat,ierr1)
  enddo
end subroutine set_topology
!-------------------------------------------------------------------------------------------------
subroutine check_topology()
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
  integer :: i,j,k,bub,ierr
  logical :: remove,signal
  !cycle through bubs first maybe? Use lagrangian approach rather than Eulerian, since there are much less bubble 
  ! parts which may be imploding than discrete points!
  signal = .false.
  fill_ghost = .false.
  remove = .false.
  v_source=0.d0
  if (num_drop(rank)>0) then
     do bub=1,num_drop(rank)
        dropid = drops(bub)%element%id
        if (.not.(rank==tag_rank(dropid))) write(*,*)"Full bub loop issue"
        if (.not.(tag_mergeflag(dropid) == 1 .or. tag_mergeflag(dropid) == 0)) then
           write(*,'("ERROR. Bub index, rank",3I8)')bub,rank
           call pariserror('Error in Lagrangian loop')
        endif
        volume = drops(bub)%element%vol
        if (volume > 1d-20) then
           if (volume < 125.0*dx(is)**3.d0) then              
              call bub_implode(bub)
           endif
        else
           write(*,'("Bubble volume error in topology check. Vol from table: ",e14.5)')volume
        endif
     enddo
  endif

  if ( num_drop_merge(rank) > 0 ) then 
     do bub=1,num_drop_merge(rank)
        dropid = drops_merge(bub)%element%id
        if (.not.(rank==tag_rank(dropid))) write(*,*)"Merged bub loop issue"
        if (.not.(tag_mergeflag(dropid) == 1 .or. tag_mergeflag(dropid) == 0)) then
           write(*,'("ERROR. Bub index, rank",3I8)')bub,rank
           call pariserror('Error in Lagrangian loop')
        endif
        volume = drops_merge(bub)%element%vol
        if (volume > 1d-20) then
           if (volume < 125.0*dx(is)**3.d0) then
              call bub_implode(bub)
           endif
        else
           write(*,'("Bubble volume error in topology check. Vol from table: ",e14.5)')volume
        endif
     enddo
  endif
  if (signal) write(*,'("COLLAPSING BUBBLE DETECTED IN RANK: ",I4)')rank
  call MPI_ALLREDUCE(remove,fill_ghost, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_Active, ierr)
  !write(*,'("Remove    : ",L8," in rank "I4)')remove,rank
  !write(*,'("Fill_ghost: ",L8," in rank "I4)')fill_ghost,rank
  
contains 
  subroutine bub_implode(bub_id)
    use module_vof
    use module_flow
    implicit none
    
    integer :: i,j,k,bub_id,inbr,jnbr,knbr
    integer :: req(4),sta(MPI_STATUS_SIZE,4),ierr
    do k=ks,ke; do j=js,je; do i=is,ie
       if (vof_phase(i,j,k)==0) then
          v_source(i,j,k) = 0.d0 ! since source is Eulerian, the moving interface may cause source to move from gas to liq
          implode(i,j,k) = 0
       endif
       if (vof_phase(i,j,k)==1) then !check if gas phase
          if (tag_dropid(tag_id(i,j,k))==bub_id) then !check if we are in the correct bubble
             implode(i,j,k) = implode(i,j,k)+1
             v_source(i,j,k) = 0.95*(u(i-1,j,k)-u(i,j,k))/dx(i)+&
                  (v(i,j-1,k)-v(i,j,k))/dy(j)+(w(i,j,k-1)-w(i,j,k))/dz(k)
             if (implode(i,j,k)==1) signal = .true.
             if ((volume < 64.0*dx(is)**3.d0) .or. implode(i,j,k)==20) remove =.true.
          endif
       endif
    enddo; enddo; enddo
    if (remove) then
       do k=ks,ke; do j=js,je; do i=is,ie
          if (tag_dropid(tag_id(i,j,k))==bub_id) then
             cvof(i,j,k) =0.d0
             v_source(i,j,k) = 0.d0
             if (itime_scheme==2) cvofold(i,j,k)=0.d0
             implode(i,j,k) = 0
             do inbr=-2,2; do jnbr=-2,2; do knbr=-2,2
                if (cvof(i+inbr,j+jnbr,k+knbr) > 0.d0) then !bubbles in very close proximity will be issue here
                   cvof(i+inbr,j+jnbr,k+knbr) = 0.d0
                   if (itime_scheme==2) cvofold(i+inbr,j+jnbr,k+knbr)=0.d0
                   implode(i+inbr,j+jnbr,k+knbr)=0
                   v_source(i+inbr,j+jnbr,k+knbr) = 0.d0 !should not be necessary
                endif
             enddo;enddo;enddo
          endif
       enddo; enddo; enddo
       write(*,'("BUBBLE TRACES REMOVED IN RANK: ",I4)')rank
    endif
  end subroutine bub_implode
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
subroutine get_ref_volume
  use module_2phase
  use module_BC
  use module_freesurface
  implicit none
  real(8), parameter :: pi=3.1415926535897932384626433833

  V_0 = 4d0/3d0*pi*(R_ref**3d0)
  R_RK = rad(1)
  dR_RK = 0d0 
  P_inf = BoundaryPressure(1)
  ddR_RK = (P_ref*(R_ref/R_RK)**(3d0*gamma)-P_inf)/R_RK
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
  !use module_Lag_part
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: p
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(in) :: A
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax) :: P_gas
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
  logical, parameter :: recordconvergence=.false.
  integer, save :: itime=0
  logical :: use_L_inf = .true.
! Open file for convergence history
  if(rank==0.and.recordconvergence) then
     OPEN(UNIT=89,FILE=TRIM(out_path)//'/convergence_history-'//TRIM(int2text(itime,padding))//'.txt')
  endif
  itime=itime+1
  if (solver_flag==1) call get_bubble_pressure(P_gas)
2 format(2e14.5)
  if (solver_flag == 0) call pariserror("Free Surface solver flag needs to be 1 or 2")
  do k=ks,ke; do j=js,je; do i=is,ie
     if (solver_flag==1 .and. pcmask(i,j,k)/=0) p(i,j,k) = P_gas(i,j,k)
     if (solver_flag==2 .and. pcmask(i,j,k)==3) p(i,j,k) = 0.d0
  enddo; enddo; enddo
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
310 format(2I6,'  ',(e14.5))
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
! Routine to numerically integrate the Rayleigh-Plesset equation
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
subroutine initialize_P_RP(p,rho)
  use module_grid
  use module_2phase
  use module_freesurface
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: p
  real(8) :: r, P_l, rho
  integer :: i,j,k,ierr
  P_l = P_ref*(R_ref/R_RK)**(3d0*gamma)-2d0*sigma/R_RK
  do k=ks,ke; do j=js,je; do i=is,ie
     r = sqrt((x(i)-xc(1))**2d0 + (y(j)-yc(1))**2d0 + (z(k)-zc(1))**2d0)
     if (r .ge. rad(1)) then
        p(i,j,k) = P_l - rho*(dR_RK**2d0 * R_RK**4d0/(2d0*r**4d0) - (ddR_RK*R_RK**2d0 + 2d0*R_RK*dR_RK**2d0)/r +&
             ddR_RK*R_RK + 3d0/2d0*dR_RK**2d0)
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
  open(unit=40,file="details.txt",access="append")
  write(40,'("Solver flag:",I5)')solver_flag
  write(40,*) "p",  p(i,j,k)
  write(40,'("A branches: ",8e14.5)')A(i,j,k,:)
  write(40,'("P neighbours ",6e14.5)')p(i-1,j,k),p(i+1,j,k),p(i,j-1,k),p(i,j+1,k),&
                  p(i,j,k-1),p(i,j,k+1)
  write(40,'("Limits: ",6I8)')is,ie,js,je,ks,ke
  write(40,*) "ijk rank",i,j,k,rank
  write(40,'("Pcmask 1-7: ",7I8)')pcmask(i-1,j,k),pcmask(i+1,j,k),pcmask(i,j-1,k),pcmask(i,j+1,k),&
       pcmask(i,j,k-1),pcmask(i,j,k+1),pcmask(i,j,k)
  close(40)
end subroutine debug_details
!====================================================================================================================================================
subroutine setuppoisson_fs_new(utmp,vtmp,wtmp,vof_phase,rho,dt,A,cvof,n1,n2,n3,kap)
  use module_grid
  use module_freesurface
  use module_IO
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: utmp,vtmp,wtmp
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: kap
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: cvof,n1,n2,n3
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(inout) :: A
  integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: vof_phase
  real(8) :: dt, rho
  integer :: i,j,k,nbr
  integer :: req(24),sta(MPI_STATUS_SIZE,24)
  
  if (.not.(solver_flag==1 .or. solver_flag==2)) call pariserror('Solver_flag FS needs to be either 1 or 2')

  if (solver_flag == 1) then
     call liq_gas()
  else !(solver_flag == 2)   
     if (.not.(solver_flag==2)) call pariserror('ERROR: Solver flag for FS must be set to either 1 or 2')
     do k=ks,ke; do j=js,je; do i=is,ie
        if (pcmask(i,j,k)==1 .or. pcmask(i,j,k)==2) then !rho is 1d0
              
           A(i,j,k,1) = 1.d0/(dx(i)*dxh(i-1))
           A(i,j,k,2) = 1.d0/(dx(i)*dxh(i))
           A(i,j,k,3) = 1.d0/(dy(j)*dyh(j-1))
           A(i,j,k,4) = 1.d0/(dy(j)*dyh(j))
           A(i,j,k,5) = 1.d0/(dz(k)*dzh(k-1))
           A(i,j,k,6) = 1.d0/(dz(k)*dzh(k))

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
           if (A(i,j,k,7)<1.d-12) pcmask(i,j,k) = 0 !Fix for isolated cav cells.
           A(i,j,k,8) =  -1.d0*((utmp(i,j,k)-utmp(i-1,j,k))/dx(i) &
                +  (vtmp(i,j,k)-vtmp(i,j-1,k))/dy(j) &
                +  (wtmp(i,j,k)-wtmp(i,j,k-1))/dz(k))
        endif
     enddo; enddo; enddo
  endif
  
contains
  subroutine liq_gas()
    use module_BC
    use module_2phase
    implicit none
    real(8) :: limit, c_min
    real(8) :: alpha2, x_test2, y_test2, z_test2
    real(8) :: nr(3),al3dnew,x0(3),dc(3),FL3DNEW,n_avg(3)
    real(8) :: c1, c0, c_stag, Source
    integer :: i,j,k,l,ierr

    x_mod=dxh((is+ie)/2); y_mod=dyh((js+je)/2); z_mod=dzh((ks+ke)/2) !assumes an unstretched grid
    P_gx = 0d0; P_gy = 0d0; P_gz = 0d0

    limit = 1d-4/dx((is+ie)/2)
    c_min = 1d-2

    do k=ks,ke; do j=js,je; do i=is,ie
       Source = 0.d0
       if (implode(i,j,k)>0 .and. vof_phase(i,j,k) == 1) Source = v_source(i,j,k)
       A(i,j,k,1) = dt/(dx(i)*dx(i)*rho)
       A(i,j,k,2) = dt/(dx(i)*dx(i)*rho)
       A(i,j,k,3) = dt/(dy(j)*dy(j)*rho)
       A(i,j,k,4) = dt/(dy(j)*dy(j)*rho)
       A(i,j,k,5) = dt/(dz(k)*dz(k)*rho)
       A(i,j,k,6) = dt/(dz(k)*dz(k)*rho)
       A(i,j,k,7) = sum(A(i,j,k,1:6))
       A(i,j,k,8) =  -(Source + (utmp(i,j,k)-utmp(i-1,j,k))/dx(i) &
            +  (vtmp(i,j,k)-vtmp(i,j-1,k))/dy(j) &
            +  (wtmp(i,j,k)-wtmp(i,j,k-1))/dz(k) )
       if (A(i,j,k,8) /= A(i,j,k,8)) then
          write(*,'("A8 NaN, liq_gas at x y z: ",3e14.5)')x(i),y(j),z(k)
          write(*,'("A branches: ",8e14.5)')A(i,j,k,:)
          write(*,'("Pcmask 1-7: ",7I8)')pcmask(i-1,j,k),pcmask(i+1,j,k),pcmask(i,j-1,k),pcmask(i,j+1,k),&
               pcmask(i,j,k-1),pcmask(i,j,k+1),pcmask(i,j,k)
          write(*,'("Implode_flags 1-7: ",7I8)')implode(i-1,j,k),implode(i+1,j,k),implode(i,j-1,k),implode(i,j+1,k),&
               implode(i,j,k-1),implode(i,j,k+1),implode(i,j,k)
          write(*,'("S_v 1-7: ",7e14.5)')v_source(i-1,j,k),v_source(i+1,j,k),v_source(i,j-1,k),v_source(i,j+1,k),&
               v_source(i,j,k-1),v_source(i,j,k+1),v_source(i,j,k)
          write(*,'("Velocities in div u: ",6e14.5)')utmp(i,j,k),utmp(i-1,j,k),vtmp(i,j,k),vtmp(i,j-1,k),&
               wtmp(i,j,k),wtmp(i,j,k-1)
       endif

       !----Cav-liquid neighbours, set P_g in cavity cells
       if (implode(i,j,k)==0) then
          if(vof_phase(i,j,k)==1) then
             do l=-1,1,2
                if (vof_phase(i+l,j,k)==0) then
                   P_gx(i,j,k) = sigma*kap(i,j,k)/dx(i) !!filaments and droplets of one cell will be an issue here
                endif
                if (vof_phase(i,j+l,k)==0) then
                   P_gy(i,j,k) = sigma*kap(i,j,k)/dy(j)
                endif
                if (vof_phase(i,j,k+l)==0) then
                   P_gz(i,j,k) = sigma*kap(i,j,k)/dz(k)
                endif
             enddo
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
                   !if (z_test2 < 0.5d0) P_gz(i,j,k) = sigma*kap(i,j,k+1)/dz(k+1)
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
                   z_mod(i,j,k) = dzh(k)*z_test2
                   if (z_mod(i,j,k)>dzh(k)) z_mod(i,j,k) = dzh(k)
                   if (z_mod(i,j,k)<limit*dzh(k)) z_mod(i,j,k) = limit*dzh(k)
                   if (z_mod(i,j,k) /= z_mod(i,j,k)) write(*,'("z_mod NaN. c_st, n, vofs, phases:",4e14.5,2I8)')&
                        c_stag,n_avg(3),cvof(i,j,k),cvof(i,j,k+1),vof_phase(i,j,k),vof_phase(i,j,k+1)
                endif !finite normal
             endif ! vof_phase test liq has gas neighbour
          endif ! vof_phase is liq
       endif ! we are not imploding
    enddo; enddo; enddo
    call ghost_x(P_gx,1,req(1:4)); call ghost_y(P_gy,1,req(5:8)); call ghost_z(P_gz,1,req(9:12)) 
    call ghost_x(x_mod,1,req(13:16)); call ghost_y(y_mod,1,req(17:20)); call ghost_z(z_mod,1,req(21:24)) 
    call MPI_WAITALL(24,req(1:24),sta(:,1:24),ierr)
    !--------------------------------------------------------------------------------------------------------
    do k=ks,ke; do j=js,je; do i=is,ie
       if (vof_phase(i,j,k)==0 .and. implode(i,j,k)==0) then
          A(i,j,k,1) = dt/(dx(i)*x_mod(i-1,j,k)*rho)
          A(i,j,k,2) = dt/(dx(i)*x_mod(i  ,j,k)*rho)
          A(i,j,k,3) = dt/(dy(j)*y_mod(i,j-1,k)*rho)
          A(i,j,k,4) = dt/(dy(j)*y_mod(i,j  ,k)*rho)
          A(i,j,k,5) = dt/(dz(k)*z_mod(i,j,k-1)*rho)
          A(i,j,k,6) = dt/(dz(k)*z_mod(i,j,k  )*rho)
          A(i,j,k,7) = sum(A(i,j,k,1:6))
          A(i,j,k,8) = A(i,j,k,8) + dt/rho*&
               (P_gx(i+1,j,k)/(dx(i)*x_mod(i,j,k))+P_gx(i-1,j,k)/(dx(i)*x_mod(i-1,j,k))&
               +P_gy(i,j+1,k)/(dy(j)*y_mod(i,j,k))+P_gy(i,j-1,k)/(dy(j)*y_mod(i,j-1,k))&
               +P_gz(i,j,k+1)/(dz(k)*z_mod(i,j,k))+P_gz(i,j,k-1)/(dz(k)*z_mod(i,j,k-1)))
          if (A(i,j,k,8) /= A(i,j,k,8)) then
             write(*,'("A8 NaN, error imminent. Neigbours mods :",6e14.5)')x_mod(i-1,j,k),x_mod(i,j,k),&
                  y_mod(i,j-1,k),y_mod(i,j,k),z_mod(i,j,k-1),z_mod(i,j,k)
             write(*,'("A8 NaN, error imminent. P_g :",6e14.5,2I8)')P_gx(i-1,j,k),P_gx(i+1,j,k),&
                  P_gy(i,j-1,k),P_gy(i,j+1,k),P_gz(i,j,k-1),P_gz(i,j,k+1),i,k
          endif
          if (pcmask(i,j,k).ne.0) write(*,'("Error topology, phase 0, pcmask :",i8)')pcmask(i,j,k)
       endif
    enddo;enddo;enddo

    if (.not. RP_test) call Poisson_BCs(A)
  end subroutine liq_gas
end subroutine setuppoisson_fs_new
!==================================================================================================================
subroutine tag_bubbles(iout,time_stats)
  use module_grid
  use module_Lag_part
  implicit none
  real(8) :: time_stats
  integer :: iout
  
  call tag_drop()
  if ( nPdomain > 1 ) call tag_drop_all
  call CreateTag2DropTable
  if ( nPdomain > 1 ) call merge_drop_pieces
  !call output_tag(iout,is,ie+1,js,je+1,ks,ke+1)
  if ( MOD(iout,nstats) == 0 ) call drop_statistics(iout,time_stats)
end subroutine tag_bubbles
!==================================================================================================================
subroutine get_bubble_pressure(P_g)
  use module_grid
  use module_Lag_part
  use module_freesurface
  use module_IO
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax) :: P_g
  real(8) :: volume
  integer :: i,j,k,dropid
  
  if (P_ref > 1.d-12) then
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
           if (volume > 1d-12) then
              P_g(i,j,k) = P_ref*(V_0/volume)**gamma
           else
              write(*,'("Bubble volume error in FreeSolver. Vol from table: ",e14.5)')volume
           endif
        endif
     enddo;enddo;enddo
  else 
     P_g = 0.d0
  endif
  call ReleaseTag2DropTable
end subroutine get_bubble_pressure
!==================================================================================================================
subroutine inflow_accelerate
  use module_grid
  use module_BC
  use module_freesurface
  use module_flow
  implicit none
  real(8) :: factor
  real(8), parameter :: pi=3.141592653589793238462643383
  integer :: i,j,k
  
  factor = (1.d0 + sin(-pi/2.d0 + DBLE(itimestep*1.0)/DBLE(step_max*1.0)*pi))/2.
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
  
end subroutine inflow_accelerate
