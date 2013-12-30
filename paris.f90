!=================================================================================================
!=================================================================================================
! Paris-0.2
! Extended from Code: FTC3D2011 (Front Tracking Code for 3D simulations)
! and Surfer. 
! 
! Authors: Sadegh Dabiri, Gretar Tryggvason.
! Author for VOF extensions Stephane Zaleski (zaleski@dalembert.upmc.fr).
! Contact: sdabiri@gmail.com .
! A three dimensional Navier-Stokes flow solver with front tracking for modeling of multiphase 
! flows. Flow can be driven by wall motion, density difference or pressure gradient.
! Boundary conditions supported: wall and periodic
!
!  HISTORY
!
! Version 1.0   1/21/2011   The 3D flow solver for variable density/viscosity is written. 
!                           The density is advected by an ENO scheme.
! Version 2.0   2/25/2011   Parallel implementation.
!                           
! FTC Version 52_4 Implicit momentum diffusion
!                  The density is advected by a QUICK scheme.
!
! Paris version 0.2 Implicit momentum diffusion, VOF and Front-Tracking independent of each other. 
!                   The density is not advected but deduced from VOF or Front Tracking. 
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
!=================================================================================================
! Main program to solve the 3D NS equations for multiphase flows
!-------------------------------------------------------------------------------------------------
Program paris
  use module_flow
  use module_grid
  use module_timer
  use module_BC
  use module_tmpvar
  use module_2phase
  use module_front

  use module_poisson
  use module_IO
  use module_solid
  use module_vof
  use module_output_vof
  use module_hello
  use module_surface_tension
  use module_st_testing
  use module_lag_part

  implicit none
  include 'mpif.h'
  integer :: ierr, icolor
  ! Locals for marching and timing are now in module_timer
  ! real(8) :: start_time, end_time=0.d0
  integer :: req(48),sta(MPI_STATUS_SIZE,48)
  INTEGER :: irank, ii, i, j, k
  real(8) :: residual,cflmax,get_cfl


!---------------------------------------INITIALIZATION--------------------------------------------
  ! Initialize MPI
  call MPI_INIT(ierr)
  if (ierr /= 0) STOP '*** Main: unsuccessful MPI-initialization'
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numProcess, ierr)
  start_time = MPI_Wtime(ierr)

  call ReadParameters
  if(rank==0) write(out,*)'Parameters read successfully'

  !Check consistency of options
  if(rank==0) then
     if(TwoPhase) then
        if((.not.DoFront).and.(.not.DoVOF)) stop 'need a phase tracking for two-phase'
        if(GetPropertiesFromFront.and.(.not.DoFront)) stop 'need Front to get properties from'
        if((.not.GetPropertiesFromFront).and.(.not.DoVOF)) stop 'need VOF to get properties from'
     endif
  endif

  ! check number of processors
  if ((NumProcess < nPdomain+1).and.DoFront) STOP '*** Main: Error with number of processes - Front!'
  if (NumProcess < nPdomain) STOP '*** Main: Error with number of processes!'

  icolor = 0                                  ! Processes 0 to nPdomain-1 are used to solve domain
  if(rank>=nPdomain) icolor = MPI_UNDEFINED
  call MPI_COMM_SPLIT(MPI_COMM_WORLD, icolor, 0, MPI_Comm_Domain, ierr)
  If (ierr /= 0) STOP '*** Main: unsuccessful MPI split'

  icolor = 0                                  ! Process nPdomain is used to solve front
  if(rank>nPdomain) icolor = MPI_UNDEFINED
  if((rank==nPdomain).and.(.not.DoFront)) icolor = MPI_UNDEFINED
  call MPI_COMM_SPLIT(MPI_COMM_WORLD, icolor, 0, MPI_Comm_Active, ierr)
  If (ierr /= 0) STOP '*** Main: unsuccessful MPI split'

  if((rank>nPdomain).or.((rank==nPdomain).and.(.not.DoFront)))then
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_finalize(ierr)
    stop
  endif

  call initialize
  call check_sanity_in_depth()
  call check_stability

  if(DoVOF.and.rank<nPdomain) call initialize_VOF
  if(DoVOF.and.DoLPP) call initialize_LPP

  if(rank<nPdomain) call initialize_solids

  if(DoFront) call InitFront
  if(rank==0) write(out,*)'initialized'
  if(rank==0) write(*  ,*)'initialized'

  if(HYPRE .and. rank<nPdomain) call poi_initialize(mpi_comm_Cart, &
                                                    is,ie,js,je,ks,ke,Nx,Ny,Nz,bdry_cond)
  if(HYPRE .and. rank==0) write(out,*)'hypre initialized'
  if(HYPRE .and. rank==0) write(*  ,*)'hypre initialized'
  
  call InitCondition
  if(rank<nPdomain) then
!-------------------------------------------------------------------------------------------------
!------------------------------------------Begin domain-------------------------------------------
!-------------------------------------------------------------------------------------------------

     ! output initial condition
     ! if(rank==0) start_time = MPI_WTIME()
     if(ICOut .and. rank<nPdomain) then
        call output(0,is,ie+1,js,je+1,ks,ke+1)
!        if(DoVOF) call output_VOF(0,imin,imax,jmin,jmax,kmin,kmax)
        call setvelocityBC(u,v,w,umask,vmask,wmask,time)
        call write_vec_gnuplot(u,v,itimestep)
        call calcstats

        if(rank==0) then
           end_time =  MPI_WTIME()
           open(unit=121,file='stats',access='append')
           write(121,'(20es14.6e2)')time,stats(1:12),dpdx,(stats(8)-stats(9))/dt,end_time-start_time
           close(121)
           write(out,'("Step:",I9," Iterations:",I9," cpu(s):",f10.2)')-1,0,end_time-start_time
           write(*,'("Step:",I6," dt=",es16.5e2," time=",es16.5e2," cpu(s):",f11.3)')   &
                0,0.d0,0.d0,end_time-start_time
           itimestep=0; ii=0
        endif
     endif
     if(test_HF.or.test_LP) then
        ! Exit MPI gracefully
        close(out)
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_finalize(ierr)
        stop
     endif
 

!-----------------------------------------MAIN TIME LOOP------------------------------------------
     call initialize_timer()
     do while(time<EndTime .and. itimestep<nstep)
        if(dtFlag==2)call TimeStepSize(dt)
        time=time+dt
        itimestep=itimestep+1

        if(mod(itimestep,termout)==0) then
           end_time =  MPI_WTIME()
           cflmax = get_cfl(dt)
           if(rank==0) &
                write(out,'("Step: ",I10," dt=",es16.5e2," time=",es16.5e2," cfl="   ,es16.5e2)                 ') &
                itimestep,dt,time                    ,cflmax
           if(rank==0) &
                write(*,  '("Step: ", I6," dt=",es16.5e2," time=",es16.5e2," cpu(s):",f11.3   ," cfl=",es16.5e2)') &
                itimestep,dt,time,end_time-start_time,cflmax
        endif

        if(itime_scheme==2) then
           uold = u
           vold = v
           wold = w
           rhoo = rho
           muold  = mu
           if ( DoVOF .and. DoLPP ) call StoreOldPartSol()
        endif
 !------------------------------------ADVECTION & DIFFUSION----------------------------------------
        do ii=1, itime_scheme
           if(TwoPhase.and.(.not.GetPropertiesFromFront)) then
             call linfunc(rho,rho1,rho2)
             call linfunc(mu,mu1,mu2)
           endif
           call my_timer(2,itimestep,ii)

           ! Receive front from master of front
           if(DoFront) call GetFront('recv')
           call my_timer(13,itimestep,ii)
           if(Implicit) then
              if(Twophase) then 
                 call momentumDiffusion(u,v,w,rho,mu,du,dv,dw)  
              else
                 du = 0d0; dv = 0d0; dw = 0d0
              endif
           else
              call explicitMomDiff(u,v,w,rho,mu,du,dv,dw)
           endif
           call my_timer(3,itimestep,ii)

           if(.not.ZeroReynolds) call momentumConvection(u,v,w,du,dv,dw)
           call my_timer(9,itimestep,ii)

           ! reset the surface tension force on the fixed grid (when surface tension from front)
           fx = 0d0;    dIdx=0d0
           fy = 0d0;    dIdy=0d0
           fz = 0d0;    dIdz=0d0
           ! Wait to finish receiving front
           if(DoFront) then
              call GetFront('wait')
              call Front2GridVector(fx, fy, fz, dIdx, dIdy, dIdz)
              call AdvanceFront2(u, v, w, color, dt)
              
              ! Send the updated front back
              call GetFront('send')
           endif
           call my_timer(13,itimestep,ii)

!------------------------------------VOF STUFF ---------------------------------------------------
           if(DoVOF) then
              if (DoLPP) then
                 call lppsweeps(itimestep)   ! XXX Note: need to update particles somewhere
                 call my_timer(12,itimestep,ii)
              end if ! DoLPP
              call vofsweeps(itimestep)
              call my_timer(4,itimestep,ii)
              call get_all_heights()
              call my_timer(5,itimestep,ii)
              call linfunc(rho,rho1,rho2)
              call surfaceForce(du,dv,dw,rho)
              call my_timer(8,itimestep,ii)
           endif

!------------------------------------END VOF STUFF------------------------------------------------ 
           call volumeForce(rho,rho1,rho2,dpdx,dpdy,dpdz,BuoyancyCase,fx,fy,fz,gx,gy,gz,du,dv,dw, &
                rho_ave)
           if(dosolids) then
              du = du*umask; dv = dv*vmask; dw = dw*wmask
           endif
           call my_timer(2,itimestep,ii)
           if(Implicit) then   
              call SetupUvel(u,du,rho,mu,rho1,mu1,dt,A)
              if(hypre)then
                 call poi_solve(A,u(is:ie,js:je,ks:ke),maxError,maxit,it)
              else
                 call LinearSolver1(A,u,umask,maxError,beta,maxit,it,ierr)
!                 if(rank==0)write(*,'("U implicit momentum diffusion iterations:",I9)')it
             endif
!             call calcresidual1(A,u,umask,residual)
!             if(rank==0)write(*,'("U implicit momentum diffusion residual:",e8.2)') residual
           
              call SetupVvel(v,dv,rho,mu,rho1,mu1,dt,A)
              if(hypre)then
                 call poi_solve(A,v(is:ie,js:je,ks:ke),maxError,maxit,it)
              else
                 call LinearSolver1(A,v,vmask,maxError,beta,maxit,it,ierr)
!                 if(rank==0)write(*  ,'("V implicit momentum diffusion  iterations:",I9)')it
              endif
!             call calcresidual1(A,v,vmask,residual)
!             if(rank==0)write(*,'("V implicit momentum diffusion residual:",e8.2)') residual

              call SetupWvel(w,dw,rho,mu,rho1,mu1,dt,A)
              if(hypre)then
                 call poi_solve(A,w(is:ie,js:je,ks:ke),maxError,maxit,it)
              else
                 call LinearSolver1(A,w,wmask,maxError,beta,maxit,it,ierr)
!                 if(rank==0)write(*  ,'("W implicit momentum diffusion  iterations:",I9)')it
              endif
!             call calcresidual1(A,w,wmask,residual)
!             if(rank==0)write(*,'("W implicit momentum diffusion residual:",e8.2)') residual
           else
              u = u + dt * du
              v = v + dt * dv
              w = w + dt * dw
           endif
           call my_timer(3,itimestep,ii)
           call SetVelocityBC(u,v,w,umask,vmask,wmask,time)

           call ghost_x(u  ,2,req( 1: 4));  call ghost_x(v,2,req( 5: 8)); call ghost_x(w,2,req( 9:12)) 
           call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
           call ghost_y(u  ,2,req( 1: 4));  call ghost_y(v,2,req( 5: 8)); call ghost_y(w,2,req( 9:12)) 
           call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
           call ghost_z(u  ,2,req( 1: 4));  call ghost_z(v,2,req( 5: 8)); call ghost_z(w,2,req( 9:12))
           call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
           call my_timer(1,itimestep,ii)

!-----------------------------------------PROJECTION STEP-----------------------------------------
           call SetPressureBC(umask,vmask,wmask,tmp(is:ie,js:je,ks:ke))
           call SetupPoisson(u,v,w,umask,vmask,wmask,rho,dt,A,tmp(is:ie,js:je,ks:ke))
           ! (div u)*dt < epsilon => div u < epsilon/dt => maxresidual : maxerror/dt 
           if(HYPRE)then
              call poi_solve(A,p(is:ie,js:je,ks:ke),maxError/dt,maxit,it)
              call ghost_x(p,1,req(1:4 ))
              call ghost_y(p,1,req(5:8 ))
              call ghost_z(p,1,req(9:12)) 
              call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
           else
              call LinearSolver(A,p,maxError/dt,beta,maxit,it,ierr)
           endif
           if(mod(itimestep,termout)==0) then
              call calcresidual(A,p,residual)
              if(rank==0)          write(*  ,    '("              pressure residual*dt:   ",e7.1,&
                   &" maxerror: ",e7.1)') residual*dt,maxerror
              if(rank==0.and..not.hypre) write(*,'("              pressure iterations :",I9)')it
           endif
           
           do k=ks,ke;  do j=js,je; do i=is,ieu;    ! CORRECT THE u-velocity 
              u(i,j,k)=u(i,j,k)-dt*(2.0/dxh(i))*(p(i+1,j,k)-p(i,j,k))/(rho(i+1,j,k)+rho(i,j,k))
           enddo; enddo; enddo
           
           do k=ks,ke;  do j=js,jev; do i=is,ie;    ! CORRECT THE v-velocity
              v(i,j,k)=v(i,j,k)-dt*(2.0/dyh(j))*(p(i,j+1,k)-p(i,j,k))/(rho(i,j+1,k)+rho(i,j,k))
           enddo; enddo; enddo
      
           do k=ks,kew;  do j=js,je; do i=is,ie;   ! CORRECT THE w-velocity
              w(i,j,k)=w(i,j,k)-dt*(2.0/dzh(k))*(p(i,j,k+1)-p(i,j,k))/(rho(i,j,k+1)+rho(i,j,k))
           enddo; enddo; enddo
           call my_timer(10,itimestep,ii)
           !--------------------------------------UPDATE COLOR---------------------------------------------
           if (DoFront) then
                 call SetupDensity(dIdx,dIdy,dIdz,A,color)
                 if(hypre)then
                    call poi_solve(A,color(is:ie,js:je,ks:ke),maxError,maxit,it)
                 else
                    call LinearSolver(A,color,maxError,beta,maxit,it,ierr)
                    if(mod(itimestep,termout)==0) then
                       if(rank==0.and..not.hypre)write(*  ,'("              density  iterations:",I9)')it
                    endif
                 endif
                 !adjust color function to 0-1 range
                 do k=ks,ke;  do j=js,je; do i=is,ie
                    color(i,j,k)=min(color(i,j,k),1d0)
                    color(i,j,k)=max(color(i,j,k),0d0)
                 enddo; enddo; enddo
           endif
           call my_timer(13,itimestep,ii)
           call SetVelocityBC(u,v,w,umask,vmask,wmask,time)
           call ghost_x(u  ,2,req( 1: 4));  call ghost_x(v,2,req( 5: 8)); call ghost_x(w,2,req( 9:12)); 
           call ghost_x(color,1,req(13:16));  call MPI_WAITALL(16,req(1:16),sta(:,1:16),ierr)
           call ghost_y(u  ,2,req( 1: 4));  call ghost_y(v,2,req( 5: 8)); call ghost_y(w,2,req( 9:12)); 
           call ghost_y(color,1,req(13:16));  call MPI_WAITALL(16,req(1:16),sta(:,1:16),ierr)
           call ghost_z(u  ,2,req( 1: 4));  call ghost_z(v,2,req( 5: 8)); call ghost_z(w,2,req( 9:12)); 
           call ghost_z(color,1,req(13:16));  call MPI_WAITALL(16,req(1:16),sta(:,1:16),ierr)
           call my_timer(1,itimestep,ii)


!--------------------------------------UPDATE DENSITY/VISCOSITY------------------------------------
           if(TwoPhase) then
              if(GetPropertiesFromFront) then
                 rho = rho2 + (rho1-rho2)*color
                 mu  = mu2  + (mu1 -mu2 )*color
              else
!------------------------------------deduce rho, mu from cvof-------------------------------------
                 call linfunc(rho,rho1,rho2)
                 call linfunc(mu,mu1,mu2)
!------------------------------------END VOF STUFF------------------------------------------------
              endif
           endif
           call my_timer(2,itimestep,ii)

           ! Wait for front to be sent back
           if(DoFront)call GetFront('wait')
           call my_timer(13,itimestep,ii)
        enddo !itime_scheme
        if(itime_scheme==2) then
           u = 0.5*(u+uold)
           v = 0.5*(v+vold)
           w = 0.5*(w+wold)
           rho = 0.5*(rho+rhoo)
           mu  = 0.5*(mu +muold)
           if ( DoVOF .and. DoLPP ) call AveragePartSol()
        endif
!--------------------------------------------OUTPUT-----------------------------------------------
        call calcStats
        call my_timer(2,itimestep,ii)
        if(mod(itimestep,nbackup)==0)call backup_write
        if(mod(itimestep,nout)==0) then 
           call write_vec_gnuplot(u,v,itimestep)
           call output(ITIMESTEP/nout,is,ie+1,js,je+1,ks,ke+1)
           if(DoVOF) call output_VOF(ITIMESTEP/nout,imin,imax,jmin,jmax,kmin,kmax)
           if(rank==0)then
              end_time =  MPI_WTIME()
              write(out,'("Step:",I9," Iterations:",I9," cpu(s):",f10.2)')itimestep,it,end_time-start_time
           endif
        endif
        if(nstats==0) STOP " *** Main: nstats = 0"
        if(mod(itimestep,nstats)==0.and.rank==0)then
              !        open(unit=121,file='track')
              !        write(121,'("Step:",I10," dt=",es16.5e2," time=",es16.5e2)')itimestep,dt,time
              !        write(121,'("            Iterations:",I7," cpu(s):",f10.2)')it,end_time-start_time
              !        close(121)
           open(unit=121,file='stats',access='append')
           write(121,'(20es14.6e2)')time,stats(1:12),dpdx,(stats(8)-stats(9))/dt,end_time-start_time
           close(121)
        endif
        call my_timer(11,itimestep,ii)
     enddo
     !-------------------------------------------------------------------------------------------------
     !--------------------------------------------End domain-------------------------------------------
     !-------------------------------------------------------------------------------------------------
     elseif((rank==nPdomain).and.DoFront) then !front tracking process (rank=nPdomain)
     !-------------------------------------------------------------------------------------------------
     !--------------------------------------------Begin front------------------------------------------
     !-------------------------------------------------------------------------------------------------
     if(ICout) call print_fronts(0,time)
     !---------------------------------------MAIN TIME LOOP--------------------------------------------
     do while(time<EndTime .and. iTimeStep<nstep)
        if(dtFlag==2)call TimeStepSize(dt)
        time=time+dt
        itimestep=itimestep+1
        
        if(itime_scheme==2) call StoreOldFront
        
        do ii=1, itime_scheme
           call CalcSurfaceTension
           do irank=0,nPdomain-1
              call DistributeFront(irank,'send') !,request(1:2,irank))
           enddo
           do irank=0,nPdomain-1
              call DistributeFront(irank,'recv') !,request(1:2,irank))
           enddo
        enddo
   
        if(itime_scheme==2) call AverageFront
        
         if(mod(itimestep,nregrid)==0) then
            print*,'Starting regrid'
            call RegridFront
            print*,'Finished regrid'
         endif
        if(smooth.and.mod(itimestep,nsmooth)==0) call smoothFront
        call CalcVolume
        if(mod(itimestep,nsmooth)==0) then 
           call CorrectVolume
           print*,'Finished volume correction'
        endif
!--------------------------------------------OUTPUT-----------------------------------------------
        if(mod(itimestep,nout)==0)call print_fronts(ITIMESTEP/nout,time)
        if(mod(itimestep,nbackup)==0)call backup_front_write(time,iTimeStep)
        if(mod(itimestep,nstats)==0)then
           open(unit=121,file='statsbub',access='append')
           do i=1, NumBubble
              write(121,'(6es16.8e2)')time,FrontProps(1,i),FrontProps(5:7,i), FrontProps(14,i)
           enddo
           close(121)
        endif
     enddo
!-------------------------------------------------------------------------------------------------
!--------------------------------------------End front--------------------------------------------
  endif
!-------------------------------------------------------------------------------------------------
!--------------- END OF MAIN TIME LOOP ----------------------------------------------------------
  call wrap_up_timer(itimestep)
  if(rank==0) then 
     if(output_format==2) call close_visit_file()
  endif

  if(rank<nPdomain)  call output_at_location()
  if(rank==0)  call final_output(stats(2))
  if(HYPRE) call poi_finalize
  call print_st_stats()
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call MPI_FINALIZE(ierr)
  stop
end program paris
!=================================================================================================
!=================================================================================================
! subroutine TimeStepSize
!-------------------------------------------------------------------------------------------------
subroutine TimeStepSize(deltaT)
  use module_grid
  use module_flow
  implicit none
  include "mpif.h"
  real(8) :: deltaT, h, vmax, dtadv, mydt
  integer :: ierr

  if(rank<nPdomain)then
    h  = minval(dx)
    vmax = maxval(sqrt(u(is:ie,js:je,ks:ke)**2+v(is:ie,js:je,ks:ke)**2+w(is:ie,js:je,ks:ke)**2))
    vmax = max(vmax,1e-3)
  !  dtadv  = min(h/vmax,2d0*nu_min/vmax**2)
    dtadv  = h/vmax
    mydt = CFL*dtadv
    mydt = min(mydt,MaxDt)
  else
    mydt=1.0e10
  endif
  call MPI_ALLREDUCE(mydt, deltat, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_Active, ierr)

end subroutine TimeStepSize
!=================================================================================================
function get_cfl(deltaT)
  use module_grid
  use module_flow
  implicit none
  include 'mpif.h'
  integer :: ierr
  real(8) :: get_cfl,vmax,inbox_cfl,deltaT,h
  vmax = maxval(sqrt(u(is:ie,js:je,ks:ke)**2 + v(is:ie,js:je,ks:ke)**2 + w(is:ie,js:je,ks:ke)**2))
  h  = minval(dx)
  inbox_cfl=vmax*deltaT/h
  call MPI_ALLREDUCE(inbox_cfl, get_cfl, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_Cart, ierr)  
end function get_cfl

!=================================================================================================
! subroutine calcStats
!-------------------------------------------------------------------------------------------------
subroutine calcStats
  use module_grid
  use module_flow
  implicit none
  include "mpif.h"
  integer :: i,j,k,ierr
  real(8) :: vol
  real(8), save :: W_int=-0.02066
  mystats(1:9)=0d0
  do k=ks,ke;  do j=js,je;  do i=is,ie
    vol = dx(i)*dy(j)*dz(k)
! Average u component
    mystats(2)=mystats(2)+u(i,j,k)*vol
    mystats(3)=mystats(3)+fx(i,j,k)*dxh(i)*dy(j)*dz(k)
    mystats(4)=mystats(4)+fy(i,j,k)*dx(i)*dyh(j)*dz(k)
    mystats(5)=mystats(5)+fz(i,j,k)*dx(i)*dy(j)*dzh(k)
    mystats(6)=mystats(6)+rho(i,j,k)*vol
    mystats(7)=mystats(7)+p(i,j,k)*vol
    mystats(8)=mystats(8)+0.5*(rho(i,j,k)+rho(i+1,j,k))*u(i,j,k)*vol
    mystats(9)=mystats(9)+0.5*(rhoo(i,j,k)+rhoo(i+1,j,k))*uold(i,j,k)*vol
  enddo;  enddo;  enddo
  if(js==Ng+1)then
    do k=ks,ke;  do i=is,ie
      mystats(1)=mystats(1)+mu(i,js,k)*u(i,js,k)*dx(i)*dz(k)/(dy(js)/2.0)
    enddo; enddo
  endif
  if(je==Ng+Ny)then
    do k=ks,ke;  do i=is,ie
      mystats(1)=mystats(1)+mu(i,je,k)*u(i,je,k)*dx(i)*dz(k)/(dy(je)/2.0)
    enddo; enddo
  endif
  mystats(2:16) = mystats(2:16)/(xLength*yLength*zLength)
  mystats(1) = mystats(1)/(xLength*zLength*2.0)
  
  call MPI_ALLREDUCE(mystats(1), stats(1), 16, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_Domain, ierr)
  rho_ave = stats(6)
  p_ave = stats(7)
! This stops the code in case the average velocity (stats(2)) becomes NaN.
  if(stats(2).ne.stats(2)) stop "********** Invalid flow rate **********"
  W_int = W_int + dt*(stats(2)-1d0)
  dpdx_stat = 1.0d0*(stats(2)-1d0) + 0.2d0*W_int
  !dpdz_stat = min(max(dpdz_stat,-2),2)
  ! time averages
  averages(1,:)=averages(1,:)+dt
  do k=ks,ke;  do j=js,je;  do i=is,ie
    Vdt = dx(i)*dy(j)*dz(k)*dt
    averages(2,j)=averages(2,j)+Vdt
    averages(3,j)=averages(3,j)+Vdt*(1.0-color(i,j,k))
    averages(4,j)=averages(4,j)+Vdt*u(i,j,k)
    averages(5,j)=averages(5,j)+Vdt*color(i,j,k)*u(i,j,k)
    averages(6,j)=averages(6,j)+Vdt*(u(i,j,k)+u(i-1,j,k))*(v(i,j,k)+v(i,j-1,k))/4.0
    averages(7,j)=averages(7,j)+Vdt*color(i,j,k)*(u(i,j,k)+u(i-1,j,k))*(v(i,j,k)+v(i,j-1,k))/4.0
    averages(8,j)=averages(8,j)+Vdt*(u(i,j,k)+u(i-1,j,k))*(w(i,j,k)+w(i,j,k-1))/4.0
    averages(9,j)=averages(9,j)+Vdt*color(i,j,k)*(u(i,j,k)+u(i-1,j,k))*(w(i,j,k)+w(i,j,k-1))/4.0
  enddo;  enddo;  enddo

  call MPI_REDUCE(averages, allaverages, 10*(Ny+2), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_Domain, ierr)

end subroutine calcStats
!=================================================================================================
!=================================================================================================
! subroutine momentumConvection
! calculates the convection terms in the momentum equations using a QUICK scheme
! and returns them in du, dv, dw
!-------------------------------------------------------------------------------------------------
subroutine momentumConvection(u,v,w,du,dv,dw)
  use module_grid
  use module_tmpvar
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: u, v, w
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: du, dv, dw
  real(8), external :: minabs
  integer :: i,j,k
!  real(8), external :: inter
!  real(8) :: y00,y11,y22
!-----------------------------------QUICK interpolation u-velocity--------------------------------
  do k=ks,ke; do j=js,je; do i=is,ieu+1
    if (u(i-1,j,k)+u(i,j,k)>0.0) then
      work(i,j,k,1) = inter(xh(i-1),xh(i),xh(i-2),x(i),u(i-1,j,k),u(i,j,k),u(i-2,j,k))  !
    else
      work(i,j,k,1) = inter(xh(i),xh(i-1),xh(i+1),x(i),u(i,j,k),u(i-1,j,k),u(i+1,j,k))  !
    endif
  enddo; enddo; enddo
  do k=ks,ke; do j=js-1,je; do i=is-1,ie
    if(v(i,j,k)+v(i+1,j,k)>0.0) then
      work(i,j,k,2) = inter(y(j),y(j+1),y(j-1),yh(j),u(i,j,k),u(i,j+1,k),u(i,j-1,k))  !
    else
      work(i,j,k,2) = inter(y(j+1),y(j),y(j+2),yh(j),u(i,j+1,k),u(i,j,k),u(i,j+2,k))  !
    endif
  enddo; enddo; enddo
  do k=ks-1,ke; do j=js,je; do i=is-1,ie
    if(w(i,j,k)+w(i+1,j,k)>0.0) then
      work(i,j,k,3) = inter(z(k),z(k+1),z(k-1),zh(k),u(i,j,k),u(i,j,k+1),u(i,j,k-1))  !
    else
      work(i,j,k,3) = inter(z(k+1),z(k),z(k+2),zh(k),u(i,j,k+1),u(i,j,k),u(i,j,k+2))  !
    endif
  enddo; enddo; enddo
  do k=ks,ke;  do j=js,je; do i=is,ieu
    du(i,j,k)=du(i,j,k) &
              -    ( work(i+1,j  ,k  ,1)**2 - work(i  ,j  ,k  ,1)**2 )/dxh(i) &
              -0.5*((v(i,j  ,k  )+v(i+1,j  ,k  ))*work(i  ,j  ,k  ,2)-        &
                    (v(i,j-1,k  )+v(i+1,j-1,k  ))*work(i  ,j-1,k  ,2))/dy(j)  &
              -0.5*((w(i,j  ,k  )+w(i+1,j  ,k  ))*work(i  ,j  ,k  ,3)-        &
                    (w(i,j  ,k-1)+w(i+1,j  ,k-1))*work(i  ,j  ,k-1,3))/dz(k)
  enddo; enddo; enddo
!-----------------------------------QUICK interpolation v-velocity--------------------------------
  do k=ks,ke; do j=js,jev+1; do i=is,ie
    if (v(i,j-1,k)+v(i,j,k)>0.0) then
      work(i,j,k,2) = inter(yh(j-1),yh(j),yh(j-2),y(j),v(i,j-1,k),v(i,j,k),v(i,j-2,k))  !
    else
      work(i,j,k,2) = inter(yh(j),yh(j-1),yh(j+1),y(j),v(i,j,k),v(i,j-1,k),v(i,j+1,k))  !
    endif
  enddo; enddo; enddo
  do k=ks,ke; do j=js-1,je; do i=is-1,ie
    if(u(i,j,k)+u(i,j+1,k)>0.0) then
      work(i,j,k,1) = inter(x(i),x(i+1),x(i-1),xh(i),v(i,j,k),v(i+1,j,k),v(i-1,j,k))  !
    else
      work(i,j,k,1) = inter(x(i+1),x(i),x(i+2),xh(i),v(i+1,j,k),v(i,j,k),v(i+2,j,k))  !
    endif
  enddo; enddo; enddo
  do k=ks-1,ke; do j=js-1,je; do i=is,ie
    if(w(i,j,k)+w(i,j+1,k)>0.0) then
      work(i,j,k,3) = inter(z(k),z(k+1),z(k-1),zh(k),v(i,j,k),v(i,j,k+1),v(i,j,k-1))  !
    else
      work(i,j,k,3) = inter(z(k+1),z(k),z(k+2),zh(k),v(i,j,k+1),v(i,j,k),v(i,j,k+2))  !
    endif
  enddo; enddo; enddo
  do k=ks,ke;  do j=js,jev; do i=is,ie
    dv(i,j,k)=dv(i,j,k) &
              -0.5*((u(i  ,j,k  )+u(i  ,j+1,k  ))*work(i  ,j  ,k  ,1)-        &
                    (u(i-1,j,k  )+u(i-1,j+1,k  ))*work(i-1,j  ,k  ,1))/dx(i)  &
              -    ( work(i  ,j+1,k  ,2)**2 - work(i  ,j  ,k  ,2)**2 )/dyh(j) &
              -0.5*((w(i  ,j,k  )+w(i  ,j+1,k  ))*work(i  ,j  ,k  ,3)-        &
                    (w(i  ,j,k-1)+w(i  ,j+1,k-1))*work(i  ,j  ,k-1,3))/dz(k)
  enddo; enddo; enddo
!-----------------------------------QUICK interpolation w-velocity--------------------------------
  do k=ks,kew+1; do j=js,je; do i=is,ie
    if (w(i,j,k-1)+w(i,j,k)>0.0) then
      work(i,j,k,3) = inter(zh(k-1),zh(k),zh(k-2),z(k),w(i,j,k-1),w(i,j,k),w(i,j,k-2))  !
    else
      work(i,j,k,3) = inter(zh(k),zh(k-1),zh(k+1),z(k),w(i,j,k),w(i,j,k-1),w(i,j,k+1))  !
    endif
  enddo; enddo; enddo
  do k=ks-1,ke; do j=js,je; do i=is-1,ie
    if(u(i,j,k)+u(i,j,k+1)>0.0) then
      work(i,j,k,1) = inter(x(i),x(i+1),x(i-1),xh(i),w(i,j,k),w(i+1,j,k),w(i-1,j,k))  !
    else
      work(i,j,k,1) = inter(x(i+1),x(i),x(i+2),xh(i),w(i+1,j,k),w(i,j,k),w(i+2,j,k))  !
    endif
  enddo; enddo; enddo
  do k=ks-1,ke; do j=js-1,je; do i=is,ie
    if(v(i,j,k)+v(i,j,k+1)>0.0) then
      work(i,j,k,2) = inter(y(j),y(j+1),y(j-1),yh(j),w(i,j,k),w(i,j+1,k),w(i,j-1,k))  !
    else
      work(i,j,k,2) = inter(y(j+1),y(j),y(j+2),yh(j),w(i,j+1,k),w(i,j,k),w(i,j+2,k))  !
    endif
  enddo; enddo; enddo
  do k=ks,kew;  do j=js,je; do i=is,ie
    dw(i,j,k)=dw(i,j,k) &
              -0.5*((u(i  ,j  ,k)+u(i  ,j  ,k+1))*work(i  ,j  ,k  ,1)-        &
                    (u(i-1,j  ,k)+u(i-1,j  ,k+1))*work(i-1,j  ,k  ,1))/dx(i)  &
              -0.5*((v(i  ,j  ,k)+v(i  ,j  ,k+1))*work(i  ,j  ,k  ,2)-        &
                    (v(i  ,j-1,k)+v(i  ,j-1,k+1))*work(i  ,j-1,k  ,2))/dy(j)  &
              -    ( work(i  ,j  ,k+1,3)**2 - work(i  ,j  ,k  ,3)**2 )/dzh(k)
  enddo; enddo; enddo
contains
!-------------------------------------------------------------------------------------------------
real(8) function inter(x0,x1,x2,x,y0,y1,y2)
  implicit none
  real(8) x0,x1,x2,x,y0,y1,y2,k,xi,a,b
  ! Interpolation at the cell face (x)
  !     |       |   u>  |
  !     x2      x0  x   x1

  	!QUICK
    xi = (x-x0)/(x1-x0)
    k = (x2-x0)/(x1-x0)
    a = (y2-k**2*y1+(k**2-1.0)*y0)/k/(1.0-k)
    b = (y2-k   *y1+(k   -1.0)*y0)/k/(k-1.0)
    inter = y0+a*xi+b*xi**2
  !  QUICK uniform
  !  inter = 0.75*y0 +0.375*y1 -0.125*y2

    !CD
  !  inter = 0.5*(y0+y1) !+0.125*(2.0*y0-y1-y2)

    !ENO 2nd order
  !  inter = y0 + 0.5*minabs(y1-y0,y0-y2)
end function inter
!-------------------------------------------------------------------------------------------------
end subroutine momentumConvection
!=================================================================================================
!=================================================================================================
! Calculates the diffusion terms in the momentum equation and adds them to du,dv,dw
!-------------------------------------------------------------------------------------------------
subroutine momentumDiffusion(u,v,w,rho,mu,du,dv,dw)
  use module_grid
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: u, v, w, rho, mu
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(out) :: du, dv, dw
  integer :: i,j,k
!-------------------------------------PREDICTED u-velocity-DIFFUSION------------------------------
  do k=ks,ke; do j=js,je; do i=is,ieu
    du(i,j,k)= ( & ! du(i,j,k)+(   &
     +(1./dy(j))*( 0.25*(mu(i,j,k)+mu(i+1,j,k)+mu(i+1,j+1,k)+mu(i,j+1,k))*               &
         ( (1./dxh(i))*(v(i+1,j,k)-v(i,j,k)) ) -    &
                   0.25*(mu(i,j,k)+mu(i+1,j,k)+mu(i+1,j-1,k)+mu(i,j-1,k))*               &
         ( (1./dxh(i))*(v(i+1,j-1,k)-v(i,j-1,k))))  &
     +(1./dz(k))*( 0.25*(mu(i,j,k)+mu(i+1,j,k)+mu(i+1,j,k+1)+mu(i,j,k+1))*               &
         ( (1./dxh(i))*(w(i+1,j,k)-w(i,j,k)) ) -    &
                   0.25*(mu(i,j,k)+mu(i+1,j,k)+mu(i+1,j,k-1)+mu(i,j,k-1))*               &
         ( (1./dxh(i))*(w(i+1,j,k-1)-w(i,j,k-1))))  &
                          )/(0.5*(rho(i+1,j,k)+rho(i,j,k))   )
  enddo; enddo; enddo
!-------------------------------------PREDICTED v-velocity-DIFFUSION------------------------------
  do k=ks,ke; do j=js,jev; do i=is,ie
    dv(i,j,k)= ( & ! dv(i,j,k)+(   &
      (1./dx(i))*( 0.25*(mu(i,j,k)+mu(i+1,j,k)+mu(i+1,j+1,k)+mu(i,j+1,k))*               &
         ( (1./dyh(j))*(u(i,j+1,k)-u(i,j,k)) ) -      &
                   0.25*(mu(i,j,k)+mu(i,j+1,k)+mu(i-1,j+1,k)+mu(i-1,j,k))*               &
         ( (1./dyh(j))*(u(i-1,j+1,k)-u(i-1,j,k)) ) )  &
     +(1./dz(k))*( 0.25*(mu(i,j,k)+mu(i,j+1,k)+mu(i,j+1,k+1)+mu(i,j,k+1))*               &
         ( (1./dyh(j))*(w(i,j+1,k)-w(i,j,k)) ) -    &
                   0.25*(mu(i,j,k)+mu(i,j+1,k)+mu(i,j+1,k-1)+mu(i,j,k-1))*               &
         ( (1./dyh(j))*(w(i,j+1,k-1)-w(i,j,k-1)) ) )  &
                           )/(0.5*(rho(i,j+1,k)+rho(i,j,k))   )
  enddo; enddo; enddo
!-------------------------------------PREDICTED w-velocity-DIFFUSION------------------------------
  do k=ks,kew; do j=js,je; do i=is,ie
    dw(i,j,k)= ( & ! dw(i,j,k)+(   &
      (1./dx(i))*( 0.25*(mu(i,j,k)+mu(i+1,j,k)+mu(i+1,j,k+1)+mu(i,j,k+1))*               &
         ((1./dzh(k))*(u(i,j,k+1)-u(i,j,k)) ) -      &
                   0.25*(mu(i,j,k)+mu(i-1,j,k)+mu(i-1,j,k+1)+mu(i,j,k+1))*               &
         ((1./dzh(k))*(u(i-1,j,k+1)-u(i-1,j,k)) ) )  &
     +(1./dy(j))*( 0.25*(mu(i,j,k)+mu(i,j+1,k)+mu(i,j+1,k+1)+mu(i,j,k+1))*               &
         ((1./dzh(k))*(v(i,j,k+1)-v(i,j,k)) ) -      &
                   0.25*(mu(i,j,k)+mu(i,j-1,k)+mu(i,j-1,k+1)+mu(i,j,k+1))*               &
         ((1./dzh(k))*(v(i,j-1,k+1)-v(i,j-1,k)) ) )  &
                          )/(0.5*(rho(i,j,k+1)+rho(i,j,k))   )
  enddo; enddo; enddo
end subroutine momentumDiffusion
!=================================================================================================
!=================================================================================================
! Calculates the diffusion terms explicitly in the momentum equation and adds them to du,dv,dw
!-------------------------------------------------------------------------------------------------
subroutine explicitMomDiff(u,v,w,rho,mu,du,dv,dw)
  use module_grid
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: u, v, w, rho, mu
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(out) :: du, dv, dw
  integer :: i,j,k
!-------------------------------------PREDICTED u-velocity-DIFFUSION------------------------------
  do k=ks,ke; do j=js,je; do i=is,ieu
    du(i,j,k)= ( & ! du(i,j,k)+(   &
      (1./dxh(i))*2.*(   mu(i+1,j,k)*(1./dx(i+1))*(u(i+1,j,k)-u(i,j,k)) -                &
                         mu(i,j,k)  *(1./dx(i))  *(u(i,j,k)-u(i-1,j,k)) )                &
     +(1./dy(j))*( 0.25*(mu(i,j,k)+mu(i+1,j,k)+mu(i+1,j+1,k)+mu(i,j+1,k))*               &
         ((1./dyh(j))*  (u(i,j+1,k)-u(i,j,k)) + (1./dxh(i))*(v(i+1,j,k)-v(i,j,k)) ) -    &
                   0.25*(mu(i,j,k)+mu(i+1,j,k)+mu(i+1,j-1,k)+mu(i,j-1,k))*               &
         ((1./dyh(j-1))*(u(i,j,k)-u(i,j-1,k)) + (1./dxh(i))*(v(i+1,j-1,k)-v(i,j-1,k))))  &
     +(1./dz(k))*( 0.25*(mu(i,j,k)+mu(i+1,j,k)+mu(i+1,j,k+1)+mu(i,j,k+1))*               &
         ((1./dzh(k))*  (u(i,j,k+1)-u(i,j,k)) + (1./dxh(i))*(w(i+1,j,k)-w(i,j,k)) ) -    &
                   0.25*(mu(i,j,k)+mu(i+1,j,k)+mu(i+1,j,k-1)+mu(i,j,k-1))*               &
         ((1./dzh(k-1))*(u(i,j,k)-u(i,j,k-1)) + (1./dxh(i))*(w(i+1,j,k-1)-w(i,j,k-1))))  &
                          )/(0.5*(rho(i+1,j,k)+rho(i,j,k))   )
  enddo; enddo; enddo
!-------------------------------------PREDICTED v-velocity-DIFFUSION------------------------------
  do k=ks,ke; do j=js,jev; do i=is,ie
    dv(i,j,k)= ( & ! dv(i,j,k)+(   &
      (1./dx(i))*( 0.25*(mu(i,j,k)+mu(i+1,j,k)+mu(i+1,j+1,k)+mu(i,j+1,k))*               &
         ((1./dyh(j))*(u(i,j+1,k)-u(i,j,k)) + (1./dxh(i))*(v(i+1,j,k)-v(i,j,k)) ) -      &
                   0.25*(mu(i,j,k)+mu(i,j+1,k)+mu(i-1,j+1,k)+mu(i-1,j,k))*               &
         ((1./dyh(j))*(u(i-1,j+1,k)-u(i-1,j,k))+ (1./dxh(i-1))*(v(i,j,k)- v(i-1,j,k))))  &
     +(1./dyh(j))*2.*(   mu(i,j+1,k)*(1./dy(j+1))*(v(i,j+1,k)-v(i,j,k)) -                &
                         mu(i,j,k)  *(1./dy(j))*  (v(i,j,k)-v(i,j-1,k)) )                &
     +(1./dz(k))*( 0.25*(mu(i,j,k)+mu(i,j+1,k)+mu(i,j+1,k+1)+mu(i,j,k+1))*               &
         ((1./dzh(k))*  (v(i,j,k+1)-v(i,j,k)) + (1./dyh(j))*(w(i,j+1,k)-w(i,j,k)) ) -    &
                   0.25*(mu(i,j,k)+mu(i,j+1,k)+mu(i,j+1,k-1)+mu(i,j,k-1))*               &
         ((1./dzh(k-1))*(v(i,j,k)-v(i,j,k-1)) + (1./dyh(j))*(w(i,j+1,k-1)-w(i,j,k-1))))  &
                           )/(0.5*(rho(i,j+1,k)+rho(i,j,k))   )
  enddo; enddo; enddo
!-------------------------------------PREDICTED w-velocity-DIFFUSION------------------------------
  do k=ks,kew; do j=js,je; do i=is,ie
    dw(i,j,k)= ( & ! dw(i,j,k)+(   &
      (1./dx(i))*( 0.25*(mu(i,j,k)+mu(i+1,j,k)+mu(i+1,j,k+1)+mu(i,j,k+1))*               &
         ((1./dzh(k))*(u(i,j,k+1)-u(i,j,k)) + (1./dxh(i))*(w(i+1,j,k)-w(i,j,k)) ) -      &
                   0.25*(mu(i,j,k)+mu(i-1,j,k)+mu(i-1,j,k+1)+mu(i,j,k+1))*               &
         ((1./dzh(k))*(u(i-1,j,k+1)-u(i-1,j,k))+ (1./dxh(i-1))*(w(i,j,k)-w(i-1,j,k))) )  &
     +(1./dy(j))*( 0.25*(mu(i,j,k)+mu(i,j+1,k)+mu(i,j+1,k+1)+mu(i,j,k+1))*               &
         ((1./dzh(k))*(v(i,j,k+1)-v(i,j,k)) + (1./dyh(j))*(w(i,j+1,k)-w(i,j,k)) ) -      &
                   0.25*(mu(i,j,k)+mu(i,j-1,k)+mu(i,j-1,k+1)+mu(i,j,k+1))*               &
         ((1./dzh(k))*(v(i,j-1,k+1)-v(i,j-1,k))+ (1./dyh(j-1))*(w(i,j,k)-w(i,j-1,k))) )  &
     +(1./dzh(k))*2.*(   mu(i,j,k+1)*(1./dz(k+1))*(w(i,j,k+1)-w(i,j,k)) -                &
                         mu(i,j,k)  *(1./dz(k))*  (w(i,j,k)-w(i,j,k-1)) )                &
                          )/(0.5*(rho(i,j,k+1)+rho(i,j,k))   )
  enddo; enddo; enddo
end subroutine explicitMomDiff
!=================================================================================================
!=================================================================================================
! Calculates the volume force in the momentum equations and adds them to du,dv,dw
!-------------------------------------------------------------------------------------------------
subroutine volumeForce(rho,rho1,rho2,dpdx,dpdy,dpdz,BuoyancyCase,fx,fy,fz,gx,gy,gz,du,dv,dw, &
                       rho_ave)
!  use module_solid
  use module_grid
  implicit none
  integer, intent(in) :: BuoyancyCase
  real(8), intent(in) :: gx,gy,gz,rho1,rho2, dpdx, dpdy, dpdz, rho_ave
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: rho
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: du, dv, dw, fx,fy,fz
  real(8) :: rro
  integer :: i,j,k
  if(BuoyancyCase==0) then
    rro = 0.0
  elseif(BuoyancyCase==1) then
    rro=rho1
  elseif(BuoyancyCase==2) then
    rro=rho2
  elseif(BuoyancyCase==3) then
    rro=rho_ave
  else
    stop 'volumeForce: invalid buoyancy option'
  endif

  do k=ks,ke;  do j=js,je; do i=is,ieu
    fx(i,j,k)=fx(i,j,k)-dpdx+(0.5*(rho(i+1,j,k)+rho(i,j,k))-rro)*gx
    du(i,j,k)=du(i,j,k) + fx(i,j,k)/(0.5*(rho(i+1,j,k)+rho(i,j,k)))
  enddo; enddo; enddo
  
  do k=ks,ke;  do j=js,jev; do i=is,ie
    fy(i,j,k)=fy(i,j,k)-dpdy+(0.5*(rho(i,j+1,k)+rho(i,j,k))-rro)*gy
    dv(i,j,k)=dv(i,j,k) + fy(i,j,k)/(0.5*(rho(i,j+1,k)+rho(i,j,k)))
  enddo; enddo; enddo

  do k=ks,kew;  do j=js,je; do i=is,ie
    fz(i,j,k)=fz(i,j,k)-dpdz+(0.5*(rho(i,j,k+1)+rho(i,j,k))-rro)*gz
    dw(i,j,k)=dw(i,j,k) + fz(i,j,k)/(0.5*(rho(i,j,k+1)+rho(i,j,k)))
  enddo; enddo; enddo

end subroutine volumeForce
!=================================================================================================
!=================================================================================================
!=================================================================================================
! Calculates the surface force in the momentum equations and adds them to du,dv,dw
!-------------------------------------------------------------------------------------------------
subroutine surfaceForce(du,dv,dw,rho)
!  use module_solid
  use module_grid
  use module_vof
  use module_2phase
  use module_surface_tension
  use module_tmpvar
  use module_timer
  implicit none
  real(8) :: kappa,deltax
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: du, dv, dw, rho
  integer :: i,j,k,l,m,n,ii,itimestep
  deltax=dx(nx)
  call get_all_curvatures(tmp)
  call my_timer(7,itimestep,ii)
  do k=ks,ke;  do j=js,je; do i=is,ieu
     if(abs(cvof(i+1,j,k)-cvof(i,j,k))>EPSC/1d1) then  ! there is a non-zero grad H (H Heaviside function) 
        n=0
        kappa=0d0
        if(tmp(i+1,j,k).lt.1e6) then
           kappa=kappa+tmp(i+1,j,k)
           n=n+1
        endif
        if(tmp(i,j,k).lt.1e6) then
           kappa=kappa+tmp(i,j,k)
           n=n+1
        endif
        if(n==0) then 
           geom_case_count(18) = geom_case_count(18) + 1
           ! call print_cvof_3x3x3(i,j,k)
        else
           kappa=kappa/(deltax*n)
           du(i,j,k) = du(i,j,k) - kappa*sigma*(2.0/dxh(i))*(cvof(i+1,j,k)-cvof(i,j,k))/(rho(i+1,j,k)+rho(i,j,k))
        endif
     endif
  enddo; enddo; enddo
  
  do k=ks,ke;  do j=js,jev; do i=is,ie
     if(abs(cvof(i,j+1,k)-cvof(i,j,k))>EPSC/1d1) then  ! there is a non-zero grad H (H Heaviside function) 
        n=0
        kappa=0d0
        if(tmp(i,j+1,k).lt.1e6) then
           kappa=kappa+tmp(i,j+1,k)
           n=n+1
        endif
        if(tmp(i,j,k).lt.1e6) then
           kappa=kappa+tmp(i,j,k)
           n=n+1
        endif
        if(n==0) then 
           geom_case_count(19) = geom_case_count(19) + 1
           ! call print_cvof_3x3x3(i,j,k)
        else
           kappa=kappa/(deltax*n)
           dv(i,j,k)=dv(i,j,k) - kappa*sigma*(2.0/dyh(j))*(cvof(i,j+1,k)-cvof(i,j,k))/(rho(i,j+1,k)+rho(i,j,k))
        endif
     endif
  enddo; enddo; enddo

  do k=ks,kew;  do j=js,je; do i=is,ie
     if(abs(cvof(i,j,k+1) - cvof(i,j,k))>EPSC/1d1) then  ! there is a non-zero grad H (H Heaviside function) 
        n=0
        kappa=0d0
        if(tmp(i,j,k+1).lt.1e6) then
           kappa=kappa+tmp(i,j,k+1)
           n=n+1
        endif
        if(tmp(i,j,k).lt.1e6) then
           kappa=kappa+tmp(i,j,k)
           n=n+1
        endif
        if(n==0) then  
           geom_case_count(20) = geom_case_count(20) + 1
           ! call print_cvof_3x3x3(i,j,k)
        else
           kappa=kappa/(deltax*n)
           dw(i,j,k)=dw(i,j,k) - kappa*sigma*(2.0/dzh(k))*(cvof(i,j,k+1)-cvof(i,j,k))/(rho(i,j,k+1)+rho(i,j,k))
        endif
     endif
  enddo; enddo; enddo

end subroutine surfaceForce
!=================================================================================================
!=================================================================================================
! subroutine initialize
!   Initializes the flow solver
!   called in:    program paris
!-------------------------------------------------------------------------------------------------
subroutine initialize
  use module_grid
  use module_flow
  use module_BC
  use module_IO
  use module_tmpvar
  use module_hello
  implicit none
  include 'mpif.h'
  integer :: ierr, i,j,k
  Nxt=Nx+2*Ng; Nyt=Ny+2*Ng; Nzt=Nz+2*Ng ! total number of cells

 
  if(rank<nPdomain)then
    allocate(dims(ndim),periodic(ndim),reorder(ndim),coords(ndim),STAT=ierr)
    dims(1) = nPx; dims(2) = nPy; dims(3) = nPz
    
    reorder = 1
    periodic = 0; 
    do i=1,ndim 
       if (bdry_cond(i) == 1) then
          periodic(i) = 1
          if (bdry_cond(i+3) /= 1) then
             if(rank==0) print*,  "inconsistent boundary conditions"
             call MPI_BARRIER(MPI_COMM_WORLD, ierr)
             call MPI_finalize(ierr)
             stop
          endif
       endif
    enddo

    call MPI_Cart_Create(MPI_Comm_Domain,ndim,dims,periodic,reorder,MPI_Comm_Cart,ierr)
    if (ierr /= 0) STOP '*** Grid: unsuccessful Cartesian MPI-initialization'

    call MPI_Cart_Coords(MPI_Comm_Cart,rank,ndim,coords,ierr)
    if (ierr /= 0) STOP '*** Grid: unsuccessful in getting topology coordinates'

!     print *, "rank",rank,"coords",coords
!     stop

    is=coords(1)*Mx+1+Ng; ie=coords(1)*Mx+Mx+Ng; imin=is-Ng; imax=ie+Ng
    js=coords(2)*My+1+Ng; je=coords(2)*My+My+Ng; jmin=js-Ng; jmax=je+Ng
    ks=coords(3)*Mz+1+Ng; ke=coords(3)*Mz+Mz+Ng; kmin=ks-Ng; kmax=ke+Ng
    ieu=ie; if(bdry_cond(1)/=1 .and. coords(1)==nPx-1) ieu=ie-1
    jev=je; if(bdry_cond(2)/=1 .and. coords(2)==nPy-1) jev=je-1
    kew=ke; if(bdry_cond(3)/=1 .and. coords(3)==nPz-1) kew=ke-1

    allocate(  u(imin:imax,jmin:jmax,kmin:kmax), uold(imin:imax,jmin:jmax,kmin:kmax), &
              du(imin:imax,jmin:jmax,kmin:kmax),   fx(imin:imax,jmin:jmax,kmin:kmax), &
               v(imin:imax,jmin:jmax,kmin:kmax), vold(imin:imax,jmin:jmax,kmin:kmax), &
              dv(imin:imax,jmin:jmax,kmin:kmax),   fy(imin:imax,jmin:jmax,kmin:kmax), &
               w(imin:imax,jmin:jmax,kmin:kmax), wold(imin:imax,jmin:jmax,kmin:kmax), &
              dw(imin:imax,jmin:jmax,kmin:kmax),   fz(imin:imax,jmin:jmax,kmin:kmax), &
             rho(imin:imax,jmin:jmax,kmin:kmax), rhoo(imin:imax,jmin:jmax,kmin:kmax), &
            drho(imin:imax,jmin:jmax,kmin:kmax),    p(imin:imax,jmin:jmax,kmin:kmax), &
              mu(imin:imax,jmin:jmax,kmin:kmax),muold(imin:imax,jmin:jmax,kmin:kmax), &
           color(imin:imax,jmin:jmax,kmin:kmax), dIdx(imin:imax,jmin:jmax,kmin:kmax), &
            dIdy(imin:imax,jmin:jmax,kmin:kmax), dIdz(imin:imax,jmin:jmax,kmin:kmax), &
           umask(imin:imax,jmin:jmax,kmin:kmax),vmask(imin:imax,jmin:jmax,kmin:kmax), &
           wmask(imin:imax,jmin:jmax,kmin:kmax))

    allocate(tmp(imin:imax,jmin:jmax,kmin:kmax), work(imin:imax,jmin:jmax,kmin:kmax,3), &
               A(is:ie,js:je,ks:ke,1:8), averages(10,Ng:Ny+Ng+1), oldaverages(10,Ng:Ny+Ng+1), &
               allaverages(10,Ng:Ny+Ng+1))

    allocate(mask(imin:imax,jmin:jmax,kmin:kmax)) 

    du=0.0;dv=0.0;dw=0.0
    u=0.0;v=0.0;w=0.0;p=0.0;tmp=0.0;fx=0.0;fy=0.0;fz=0.0;drho=0.0;rho=0.0;mu=0.0;work=0.0;A=0.0
    averages=0.0; oldaverages=0.0; allaverages=0d0
    umask = 1d0; vmask = 1d0; wmask = 1d0

  else  !   if(rank<nPdomain)then
    is = Ng+1;  ie = Nx+Ng;  imin = 1;  imax = Nxt
    js = Ng+1;  je = Ny+Ng;  jmin = 1;  jmax = Nyt
    ks = Ng+1;  ke = Nz+Ng;  kmin = 1;  kmax = Nzt
    ieu = ie;  if(bdry_cond(1)/=1) ieu = ie-1      ! only if not periodic
    jev = je;  if(bdry_cond(2)/=1) jev = je-1
    kew = ke;  if(bdry_cond(3)/=1) kew = ke-1
  endif !   if(rank<nPdomain)then

  allocate(x(Nxt), xh(Nxt), dx(Nxt), dxh(Nxt), &
           y(Nyt), yh(Nyt), dy(Nyt), dyh(Nyt), &
           z(Nzt), zh(Nzt), dz(Nzt), dzh(Nzt)  )

! Set the stretched grid
  do i=1,Nxt; s=dfloat(i-Ng)/dfloat(Nx); xh(i)=xLength*(xform*s*(0.5d0-s)*(1d0-s)+s); enddo
  do j=1,Nyt; s=dfloat(j-Ng)/dfloat(Ny); yh(j)=yLength*(yform*s*(0.5d0-s)*(1d0-s)+s); enddo
  do k=1,Nzt; s=dfloat(k-Ng)/dfloat(Nz); zh(k)=zLength*(zform*s*(0.5d0-s)*(1d0-s)+s); enddo

  if(read_x)then
    open(unit=12,file=trim(x_file),status='old',action='read')
    read(12,*) (xh(i), i=Ng,Nx+Ng)
    close(12)
    xLength = xh(Nx+Ng)-xh(Ng)
    i=Nx+Ng+1
    xh(i) = 2d0*xh(i-1)-xh(i-2)
    if(Ng==2)xh(1) = 2d0*xh(2)-xh(3)
    if(Ng==2)xh(i+1) = 2d0*xh(i)-xh(i-1)
  endif
  if(read_y)then
    open(unit=12,file=trim(y_file),status='old',action='read')
    read(12,*) (yh(i), i=Ng,Ny+Ng)
    close(12)
    yLength = yh(Ny+Ng)-yh(Ng)
    i=Ny+Ng+1
    yh(i) = 2d0*yh(i-1)-yh(i-2)
    if(Ng==2)yh(1) = 2d0*yh(2)-yh(3)
    if(Ng==2)yh(i+1) = 2d0*yh(i)-yh(i-1)
  endif
  if(read_z)then
    open(unit=12,file=trim(z_file),status='old',action='read')
    read(12,*) (zh(i), i=Ng,Nz+Ng)
    close(12)
    zLength = zh(Nz+Ng)-zh(Ng)
    i=Nz+Ng+1
    zh(i) = 2d0*zh(i-1)-zh(i-2)
    if(Ng==2)zh(1) = 2d0*zh(2)-zh(3)
    if(Ng==2)zh(i+1) = 2d0*zh(i)-zh(i-1)
  endif

  do i=2,Nxt; x(i)=0.5d0*(xh(i)+xh(i-1)); enddo; x(1)=2d0*xh(1)-x(2)
  do j=2,Nyt; y(j)=0.5d0*(yh(j)+yh(j-1)); enddo; y(1)=2d0*yh(1)-y(2)
  do k=2,Nzt; z(k)=0.5d0*(zh(k)+zh(k-1)); enddo; z(1)=2d0*zh(1)-z(2)

  do i=1,Nxt-1; dxh(i)=x(i+1)-x(i); enddo; dxh(Nxt)=dxh(Nxt-1)
  do j=1,Nyt-1; dyh(j)=y(j+1)-y(j); enddo; dyh(Nyt)=dyh(Nyt-1)
  do k=1,Nzt-1; dzh(k)=z(k+1)-z(k); enddo; dzh(Nzt)=dzh(Nzt-1)
  
  do i=2,Nxt; dx(i)=xh(i)-xh(i-1); enddo; dx(1)=dx(2);
  do j=2,Nyt; dy(j)=yh(j)-yh(j-1); enddo; dy(1)=dy(2);
  do k=2,Nzt; dz(k)=zh(k)-zh(k-1); enddo; dz(1)=dz(2);

end subroutine initialize


!=================================================================================================
! subroutine InitCondition
!   Sets the initial conditions
!   called in:    program paris
!-------------------------------------------------------------------------------------------------
subroutine InitCondition
  use module_hello
  use module_grid
  use module_flow
  use module_BC
  use module_front
  use module_tmpvar
  use module_poisson
  use module_IO
  use module_vof
  use module_surface_tension
  use module_st_testing
  use module_lag_part
  implicit none
  include 'mpif.h'
  integer :: i,j,k, ierr, irank, req(12),sta(MPI_STATUS_SIZE,12)
  real(8) :: my_ave
  !---------------------------------------------Domain----------------------------------------------
  if(rank<nPdomain)then
     if(restart)then
        call backup_read
        call SetVelocityBC(u,v,w,umask,vmask,wmask,time)
        call ghost_x(u,2,req( 1: 4)); call ghost_x(v,2,req( 5: 8)); call ghost_x(w,2,req( 9:12))
        call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
        call ghost_y(u,2,req( 1: 4)); call ghost_y(v,2,req( 5: 8)); call ghost_y(w,2,req( 9:12))
        call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
        call ghost_z(u,2,req( 1: 4)); call ghost_z(v,2,req( 5: 8)); call ghost_z(w,2,req( 9:12))
        call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
     else
        ! Set velocities and the color function. 
        ! The color function is used for density and viscosity in the domain 
        ! when set by Front-Tracking.
        color = 0.;  v = 0;  w = 0.
        u = 0.d0;
        if((y(j) - ylength*0.5d0)**2 + (z(k) - zlength*0.5d0)**2.lt.jetradius**2)  u = U_init

        if(DoVOF) then
           call initconditions_VOF()
           call get_all_heights()
        endif
        du = 0d0
     endif
     
     if(test_HF) then
        call test_VOF_HF()
     else if (test_LP) then 
        call test_Lag_part(itimestep)
     endif
     if(DoFront) then
        call GetFront('recv')
        call GetFront('wait')
        call Front2GridVector(fx, fy, fz, dIdx, dIdy, dIdz)
        call SetupDensity(dIdx,dIdy,dIdz,A,color)
        if(hypre) then
           call poi_solve(A,color(is:ie,js:je,ks:ke),maxError,maxit,it)
           call ghost_x(color,1,req( 1: 4))
           call ghost_y(color,1,req( 5: 8))
           call ghost_z(color,1,req( 9:12))
           call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
        else
           call LinearSolver(A,color,maxError,beta,maxit,it,ierr)
           if(rank==0)print*,it,'iterations for initial density.'
        endif
        do k=ks,ke;  do j=js,je; do i=is,ie
           color(i,j,k)=min(color(i,j,k),1d0)
           color(i,j,k)=max(color(i,j,k),0d0)
        enddo; enddo; enddo
     endif

     if(TwoPhase) then
        if(GetPropertiesFromFront) then
           rho = rho2 + (rho1-rho2)*color
           mu  = mu2  + (mu1 -mu2 )*color
        else
           call linfunc(rho,rho1,rho2)
           call linfunc(mu,mu1,mu2)
        endif
     else
        rho=rho1
        mu=mu1
     endif
     my_ave=0.0
     do k=ks,ke;  do j=js,je;  do i=is,ie
        my_ave=my_ave+rho(i,j,k)*dx(i)*dy(j)*dz(k)
     enddo;  enddo;  enddo
     my_ave = my_ave/(xLength*yLength*zLength)
     call MPI_ALLREDUCE(my_ave, rho_ave, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_Domain, ierr)
!----------------------------------------------Front----------------------------------------------
  elseif((rank==nPdomain).and.DoFront)then
     if(restartFront)then
        call backup_front_read(time,iTimeStep)
        call RegridFront
     else
        call InitConditionFront
        write(*,'(3f20.10)') FrontProps(:,1)
        call CalcVolume
        write(*,'(3f20.10)') FrontProps(:,1)
        FrontProps(2,1:NumBubble) = FrontProps(1,1:NumBubble)
     endif

     call CalcSurfaceTension
     do irank=0,nPdomain-1
        call DistributeFront(irank,'send') !,request(1:2,irank))
     enddo

  endif
  !---------------------------------------------End Front----------------------------------------------
  if((.not. restart).or.(.not. restartFront)) then
     time = 0d0
     iTimeStep = 0
  endif
end subroutine InitCondition
!=================================================================================================
! function minabs
!   used for ENO interpolations
!   called in:    subroutine momentumConvection
!                 subroutine density
!-------------------------------------------------------------------------------------------------
function minabs(a,b)
implicit none
real(8) :: minabs, a, b
  if(abs(a)<abs(b)) then
    minabs=a
  else
    minabs=b
  endif
!  minabs = 0.5*(sign(1.0,abs(b)-abs(a))*(a-b)+a+b)
end function minabs
!=================================================================================================
!=================================================================================================
! function text
!   Returns 'number' as a string with length of 'length'
!   called in:    function output
!-------------------------------------------------------------------------------------------------
function int2text(number,length)
  integer :: number, length, i
  character(len=length) :: int2text
  character, dimension(0:9) :: num = (/'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/)
  if(number>=10**length)stop 'int2text error: string is not large enough'
  do i=1,length
    int2text(length+1-i:length+1-i) = num(mod(number/(10**(i-1)),10))
  enddo
end function
!=================================================================================================
function sphere(x,y,z)
  real(8) :: sphere
  real(8) :: x,y,z
  sphere=0.d0
  if(x*x + y*y + z*z < 0.09) sphere=1.d0
end function
!=================================================================================================
! subroutine ReadParameters
!   Reads the parameters in "input" file for the calculation
!   called in:    program paris
!-------------------------------------------------------------------------------------------------
subroutine ReadParameters
  use module_grid
  use module_flow
  use module_2phase
  use module_BC
  use module_IO
  use module_front
  use module_hello
  implicit none
  include 'mpif.h'
  integer :: in, ierr
  real(8) :: xyzrad(4,10000)
  namelist /parameters/ out_path,      Nx,            Ny,            Nz,            Ng,          &
                        xLength,       yLength,       zLength,       gx,            gy,          &
                        gz,            bdry_cond,     dPdx,          dPdy,          dPdz,        &
                        itime_scheme,  nstep,         maxit,         maxError,      &
                        beta,          nout,          TwoPhase,      rho1,          mu1,         &
                        rho2,          mu2,           sigma,         BuoyancyCase,  nPx,         &
                        nPy,           nPz,           amin,          amax,          aspmax,      &
                        MaxPoint,      MaxElem,       MaxFront,      xform,         yform,       &
                        zform,         dt,            nregrid,       GetPropertiesFromFront,     &
                        DoVOF,         DoFront,       Implicit,      U_init,                     &
                        CFL,           EndTime,       MaxDt,         smooth,        nsmooth,     &
                        output_format, read_x,        read_y,        read_z,        x_file,      &
                        y_file,        z_file,        restart,       nBackup,       NumBubble,   &
                        xyzrad,        hypre,         dtFlag,        ICOut,         WallVel,     &
                        maxErrorVol,   restartFront,  nstats,        WallShear,     ZeroReynolds,&
                        restartAverages, termout, excentricity
  in=1
  out=2

  open(unit=in, file='input', status='old', action='read', iostat=ierr)
  if (ierr .ne. 0) stop 'ReadParameters: error opening input file'
  read(UNIT=in,NML=parameters)
  close(in)
  call check_sanity()
  bdry_read=.true.
  if(MaxFront>10000) stop 'Error: ReadParameters: increase size of xyzrad array'

  if(numBubble>MaxFront) stop 'Error: ReadParameters: increase size of xyzrad array (MaxFront)'

  allocate(xc(MaxFront), yc(MaxFront), zc(MaxFront))
  allocate(FrontProps(1:14,MaxFront),rad(MaxFront))

  xc (1:NumBubble) = xyzrad(1,1:NumBubble)
  yc (1:NumBubble) = xyzrad(2,1:NumBubble)
  zc (1:NumBubble) = xyzrad(3,1:NumBubble)
  rad(1:NumBubble) = xyzrad(4,1:NumBubble)


  FrontProps(5,1:NumBubble) = xyzrad(1,1:NumBubble)
  FrontProps(6,1:NumBubble) = xyzrad(2,1:NumBubble)
  FrontProps(7,1:NumBubble) = xyzrad(3,1:NumBubble)
  rad(1:NumBubble) = xyzrad(4,1:NumBubble)

  if(rank==0)then
     call system('mkdir            '//trim(out_path))
     call system('mkdir            '//trim(out_path)//'/VTK')
     call system('cp input         '//trim(out_path))
     !call system('cp paris.f90 '//trim(out_path))
     open(unit=out, file=trim(out_path)//'/output', action='write', iostat=ierr)
     if (ierr .ne. 0) stop 'ReadParameters: error opening output file'
     write(UNIT=out,NML=parameters)
  endif

  ! Number of grid points in streamwise and transverse directions must
  ! be integral multiples of total number of processors
  if(mod(Nx,nPx) /= 0) Stop 'ReadParameters: Nx not divisible by nPx!'
  if(mod(Ny,nPy) /= 0) Stop 'ReadParameters: Ny not divisible by nPy!'
  if(mod(Nz,nPz) /= 0) Stop 'ReadParameters: Nz not divisible by nPz!'
  Mx = Nx/nPx; My = Ny/nPy; Mz = Nz/nPz
  nPdomain = nPx*nPy*nPz

!--- output frequency

  if(termout==0) then
     termout=nout
  endif

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
end subroutine ReadParameters
!=================================================================================================
!=================================================================================================
!=================================================================================================
subroutine pariserror(message) 
  use module_IO
  use module_grid
  implicit none
  include 'mpif.h'
!  integer ierr
  character(*) :: message
  print *, "rank = ",rank
  write(*,*) "ERROR *** ",message, " *** STOP "
  ! Exit MPI gracefully
  close(out)
!   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!   call MPI_finalize(ierr)
  if(rank==0) write(*,*) "Step: last message . . . . ParisExecutionError"
!  if(rank==0) write(*,*) "ERROR *** ",message, " *** STOP "
  stop 
end subroutine pariserror

!=================================================================================================
subroutine check_stability() 
  use module_grid
  use module_flow
  use module_BC
  use module_IO
  use module_front
  use module_hello
  use module_2phase
  implicit none
  include 'mpif.h'
  real(8) von_neumann

  if(dt*mu1 /= 0) then 
     von_neumann = dx(ng)**2*rho1/(dt*mu1)
     if(rank==0) print *, "dx**2*rho/(dt*mu) = ", von_neumann
     if(von_neumann < 6d0.and..not.Implicit) call pariserror("time step too large for viscous terms")
  endif

  if(dt*sigma /= 0) then 
     von_neumann = dx(ng)**3*rho1/(dt**2*sigma)
     if(rank==0) print *, "dx**3*rho/(dt**2*sigma) = ", von_neumann
     if(von_neumann < 4d0) call pariserror("time step too large for ST terms")
  endif

end subroutine check_stability
