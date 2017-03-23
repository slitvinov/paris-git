!=================================================================================================
!
!
! PARIS (Parallel Robust Interface Software) Simulator 
! Main program 
!
! Contact: Stephane Zaleski zaleski@dalembert.upmc.fr
! 
! Authors (in alphabetical order):
!         Tomas Arrufat Jackson
!         Sadegh Dabiri      
!         Daniel Fuster      
! 	  Yue "Stanley" Ling 
!         Jacai Lu
!         Leon Malan
!         Ruben Scardovelli  
!         Gretar Tryggvason  
!         Phil Yecko         
!         Stephane Zaleski   
!
! Extended from or inspired by Codes: 
!      - FTC3D2011 (Front Tracking Code for 3D simulations) by Gretar Tryggvason and Sadegh Dabiri      
!      - Surfer VOF code by Stephane Zaleski, Ruben Scardovelli and others
!      - Gerris Flow Solver by Stephane Popinet
!
! GPL Licence
!
!     This file is part of PARIS.
!
!     PARIS is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
! 
!     PARIS is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with PARIS.  If not, see <http://www.gnu.org/licenses/>.
!
! 
!=================================================================================================
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
Program paris
  use module_flow
  use module_grid
  use module_timer
  use module_BC
  use module_tmpvar
  use module_2phase
  use module_freesurface
  use module_front

  use module_poisson
  use module_mgsolver
!  use module_averages
  use module_IO
  use module_solid
  use module_vof
  use module_output_vof
  
  use module_surface_tension
  use module_st_testing
  use module_lag_part
  use module_output_LPP

#ifdef PHASE_CHANGE
  use module_boil
#endif


  implicit none
  include 'mpif.h'
  integer :: ierr, icolor

  INTEGER :: irank, ii, i, j, k, out_fs
  real(8) :: residual,cflmax,get_cfl_and_check,residualu,residualv,residualw,res_tab(3)
  integer :: itu,itv,itw

!---------------------------------------INITIALIZATION--------------------------------------------
  ! Initialize MPI
  call MPI_INIT(ierr)
  if (ierr /= 0) call err_no_out_dir("*** Main: unsuccessful MPI-initialization")
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numProcess, ierr)
  start_time = MPI_Wtime(ierr)

  call ReadParameters
  if(rank==0) write(out,*)'Parameters read successfully'

  !Check consistency of options
  if(rank==0) then
     if(TwoPhase.or.FreeSurface) then
        if((.not.DoFront).and.(.not.DoVOF)) call pariserror("need a phase tracking for two-phase")
        if(GetPropertiesFromFront.and.(.not.DoFront)) call pariserror("need Front to get properties from")
        if((.not.GetPropertiesFromFront).and.(.not.DoVOF)) call pariserror("need VOF to get properties from")
     endif
     if(TwoPhase.and.FreeSurface) call pariserror("cannnot be both TwoPhase and FreeSurface")
  endif

  ! check number of processors
  if ((NumProcess < nPdomain+1).and.DoFront) call pariserror("*** Main: Error with number of processes - Front!")
  if (NumProcess < nPdomain) call pariserror("*** Main: Error with number of processes!")

  icolor = 0                                  ! Processes 0 to nPdomain-1 are used to solve domain
  if(rank>=nPdomain) icolor = MPI_UNDEFINED
  call MPI_COMM_SPLIT(MPI_COMM_WORLD, icolor, 0, MPI_Comm_Domain, ierr)
  If (ierr /= 0) call pariserror("*** Main: unsuccessful MPI split")

  icolor = 0                                  ! Process nPdomain is used to solve front
  if(rank>nPdomain) icolor = MPI_UNDEFINED
  if((rank==nPdomain).and.(.not.DoFront)) icolor = MPI_UNDEFINED
  call MPI_COMM_SPLIT(MPI_COMM_WORLD, icolor, 0, MPI_Comm_Active, ierr)
  If (ierr /= 0) call pariserror("*** Main: unsuccessful MPI split")

  if((rank>nPdomain).or.((rank==nPdomain).and.(.not.DoFront)))then
     call pariserror("rank>nPdomain).or.((rank==nPdomain).and.(.not.DoFront))")
  endif

  call initialize
  call check_sanity_in_depth()

  if(DoVOF.and.rank<nPdomain) call initialize_VOF
  if(DoLPP) call initialize_LPP

  if(rank<nPdomain) call initialize_solids

  if(DoFront) call InitFront
  if(rank==0) write(out,*)'initialized'
  if(rank==0) write(*  ,*)'initialized'

  if(HYPRE .and. rank<nPdomain) call poi_initialize(mpi_comm_Cart, &
                                                    is,ie,js,je,ks,ke,Nx,Ny,Nz,bdry_cond)
  if(HYPRE .and. rank==0) write(out,*)'hypre initialized'
  if(HYPRE .and. rank==0) write(*  ,*)'hypre initialized'

  call InitCondition

  if (MultiGrid) call init_MG(1)

#ifdef PHASE_CHANGE
  call ReadHeatParameters
  if (rank==0) then
     write(*,'("Running simulation with heat transfer: Parameters read successfully & initialized")')
  endif
  if (TwoPhase.and.(.not.GetPropertiesFromFront)) then
     call linfunc(kc,kc1,kc2,HarmMean)
  else
     kc=kc1
     Cp=Cp1
  endif
  call SetTempBC
  call do_all_ghost(Te)
#endif

    if(rank<nPdomain) then
!-------------------------------------------------------------------------------------------------
!------------------------------------------Begin domain-------------------------------------------
!-------------------------------------------------------------------------------------------------

     ! output initial condition
     ! if(rank==0) start_time = MPI_WTIME()
     if(ICOut .and. rank<nPdomain) then
        if (.not.restart) call output(0,is,ie+1,js,je+1,ks,ke+1)
        if(DoVOF .and. .not.restart) then
           call output_VOF(0,is,ie+1,js,je+1,ks,ke+1)
           if(out_centroid) call output_centroids(0)
           call output_ALL(0,is,ie+1,js,je+1,ks,ke+1,itimestep,5)
           if (FreeSurface) then
              do out_fs = 1,3
                 if (VTK_OUT(out_fs)) then
                    if (rank==0) call append_visit_fs(out_fs,0)
                    SELECT CASE (out_fs)
                    case (1) 
                       tmp = 0d0
                       do k=ks,ke+1; do j=js,je+1; do i=is,ie+1
                          tmp(i,j,k)=ABS((u(i-1,j,k)-u(i,j,k))/dx(i)+(v(i,j-1,k)-v(i,j,k))/dy(j)+(w(i,j,k-1)-w(i,j,k))/dz(k))
                       enddo; enddo; enddo
                       call VTK_scalar_struct(out_fs,0,tmp)
                    case(2)
                       call VTK_scalar_struct(out_fs,0,P_gx)
                    case(3)
                       call VTK_scalar_struct(out_fs,0,v_source)
                    end SELECT
                 endif
              enddo
           endif
        endif
        if(DoLPP .and. DoOutputLPP .and. .not.restart) call output_LPP(0,is,ie+1,js,je+1,ks,ke+1)
        if(test_control_droplet) call do_droplet_test(itimestep,time,REAL(nstats*dt,8))
        if(test_frdroplet.or.test_droplet) call output_droplet(w,time)
        call SetVelocityBC(u,v,w,umask,vmask,wmask,time,dt,0)
        if (FreeSurface) then
           if (inflow .and. (itimestep<=step_max .or. fs_refactor)) call inflow_accelerate  
        endif
        !call write_vec_gnuplot(u,v,cvof,p,itimestep,DoVOF)
        call calcstats

        if(dtFlag==2)call TimeStepSize(dt,vof_phase)
        cflmax = get_cfl_and_check(dt)
        ! stability check positionned after time step chage
        call check_stability

        if(rank==0) then
           end_time =  MPI_WTIME()
           if (time.eq.0.d0) then
              open(unit=121,file='stats')
              ![WA] adding labels 2017/02/07
              write(121,*) '# 1: time 2: shear_y 3: <u> 4: f_x 5: f_y 6: f_z 7: <rho> 8: <p> 9: <mom>^n &
                   & 10: <mom>^(n+1) 11: <C> 12: x_center_of_mass &
                   & 13: E_k1 14: E_k2 15-16: y_mom 17: undefined 18: <p_infl> 19:<p_outfl> 20:Area 21:Ep1 22:Ep2 23: Enstrophy1 & 
                   & 24: Enstrophy2 25: R2_frac 26: Vof_flag2 27: dpdx 28: du/dt 29: time_since_restart'
           else
             open(unit=121,file='stats',position='append')
           endif
           write(121,'(30es14.6e2)')time,stats(1:nstatarray),dpdx,(stats(8)-stats(9))/dt,end_time-start_time
           close(121)
           write(out,'("Step:",I9," Iterations:",I9," cpu(s):",f10.2)')-1,0,end_time-start_time
           write(*,  '("START:", I6," dt=",es16.5e2," time=",es16.5e2," cpu(s):",f11.3   ," cfl=",es16.5e2)') &
                itimestep,dt,time,end_time-start_time,cflmax
!           itimestep=0; ii=0
        endif
     endif

     if (test_MG) then
       call linfunc(rho,rho1,rho2,DensMean)
       p=0.d0; call Setup_testMG(rho,A)
       call NewSolver(A,p,maxError,beta,maxit,it,ierr,ResNormOrderPressure)
       call get_MGtest_err(p)
       end_time =  MPI_WTIME()
       cflmax = get_cfl_and_check(dt)
       if(rank==0) print *, 'cpu =', end_time-start_time
     endif

     if(test_HF.or.test_LP.or.test_MG) then
        ! Exit MPI gracefully
        close(out)
        call print_st_stats(0)
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        if(rank==0) write(*,'("Paris exits succesfully after HF, curvature or LP test")')
        call MPI_finalize(ierr)
        stop
     endif
 !-----------------------------------------MAIN TIME LOOP------------------------------------------
     call initialize_timer()

     do while(time<EndTime .and. itimestep<nstep)
        if (test_capwave.and.mod(itimestep,nstats)==0) call interface_max_min
        if(dtFlag==2)call TimeStepSize(dt,vof_phase)
        time=time+dt
        itimestep=itimestep+1
        if(mod(itimestep,termout)==0) then
           end_time =  MPI_WTIME()
           cflmax = get_cfl_and_check(dt)
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
           if(DoVOF) cvofold  = cvof
           if ( DoLPP ) call StoreOldPartSol()
#ifdef PHASE_CHANGE
           Te_old = Te
           kc_old = kc
           Cp_old = Cp
#endif
        endif
 !------------------------------------ADVECTION & DIFFUSION----------------------------------------
        du = 0d0; dv = 0d0; dw = 0d0
        do ii=1, itime_scheme
           if(TwoPhase.and.(.not.GetPropertiesFromFront)) then
             call linfunc(rho,rho1,rho2,DensMean)
             call linfunc(mu,mu1,mu2,ViscMean) ! Note: harmonic mean matches shear stress better
           endif
           call my_timer(2)

#ifdef PHASE_CHANGE
           dTe = 0.0d0
           call TempDiffusion
#endif
           if( DoLPP ) call StoreBeforeConvectionTerms()
           if(.not.ZeroReynolds) call momentumConvection()
           if( DoLPP ) call StoreAfterConvectionTerms()
           call my_timer(9)
           if(DoVOF) then
#ifdef PHASE_CHANGE
              call vofandenergysweeps(itimestep) !Advects energy using mass flux from vof advection
              call get_heat_source
#endif
              if (DoMOMCONS) then
                 call vofandmomsweepsstaggered(itimestep,time)
                 call vofsweeps(itimestep)
              else
                 call vofsweeps(itimestep)
              endif
              
              if ( MOD(itimestep,nsteps_clean_debris) == 0 .and. do_clean_debris ) then 
                 call clean_debris
                 if ( nPdomain > 1 ) & 
                      call VOFCommGhost    ! Communicate vof field after cleaning
              end if ! nsteps_clean_debris
              
              call my_timer(4)
              call get_all_heights(itimestep)
              call my_timer(5)
              call linfunc(rho,rho1,rho2,DensMean)
              if ((.not.STGhost).and.(.not.FreeSurface).and.(sigma.gt.TINY_DOUBLE)) &
                   call surfaceForce(du,dv,dw,rho)
              call my_timer(8)
              
              if (FreeSurface) then
                 !call fs_sweep(itimestep,time) !Will include all calls below, to be implemented
                 if (check_stray_liquid .and. mod(itimestep,n_stray_liquid)==0) then
                    !calls below are just for cleaning liquid which got into gas
                    !check structures
                    call tag_bubbles(0,itimestep,time)
                    call check_topology(.false.)
                    !for dealloc
                    call ReleaseTag2DropTable
                    !that's valid if you have to fill the ghost cells because u remved sth
                    if (fill_ghost) then
                       call clean_debris
                       call do_all_ghost(cvof)
                       call get_flags_and_clip(cvof,vof_flag)
                       call get_vof_phase(cvof) !cvof updated above from min to max
                    endif
                 endif
                 !removing small gas pockets
                 if (.not. (test_capwave .or. test_plane)) then
                    call tag_bubbles(1,itimestep,time)
                    if (.not.RP_test) then
                       call check_topology(.true.)
                       if (fill_ghost) then
                          call clean_debris
                          call do_all_ghost(cvof)
                          call get_flags_and_clip(cvof,vof_flag)
                          call get_vof_phase(cvof) !cvof updated above from min to max
                          call ReleaseTag2DropTable
                          call tag_bubbles(1,itimestep,time)
                          call get_all_heights(itimestep)
                       endif
                    endif
                 endif
                 !so above code is basically removing stuff and setting tags
                 !after possible removal above
                 call set_topology(vof_phase,itimestep) !vof_phase updated in vofsweeps
                 call set_bubble_pressure
                 if (.not. (test_capwave .or. test_plane)) call ReleaseTag2DropTable
                 call my_timer(15)
              endif ! FreeSurface
           endif
#ifdef PHASE_CHANGE
           if (TwoPhase.and.(.not.GetPropertiesFromFront)) then
              call linfunc(kc,kc1,kc2,HarmMean)
              call linfunc(Cp,Cp1,Cp2,ArithMean)
           endif
           do i=is,ie; do j=js,je; do k=ks,ke
              Te(i,j,k) = Te(i,j,k) + dt*dTe(i,j,k) !Simple explicit treatment of advection/diffusion
           enddo; enddo; enddo
           call SetTempBC
           call do_all_ghost(Te)
#endif
           if (DoLPP) then
                call lppsweeps(itimestep,time,ii)  
                call my_timer(12)
           end if ! DoLPP
!------------------------------------FRONT TRACKING  ---------------------------------------------
           ! Receive front from master of front
           if (.not. FreeSurface) then
              if(DoFront) call GetFront('recv')
              call my_timer(13)
              if(Implicit) then
                 if(Twophase) then 
                    call momentumDiffusion(u,v,w,rho,mu,du,dv,dw)  
                 endif
              else
                 call explicitMomDiff(u,v,w,rho,mu,du,dv,dw)
              endif
              call my_timer(3)

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
              call my_timer(13)

              !------------------------------------END VOF STUFF------------------------------------------------ 
              call volumeForce(rho,rho1,rho2,dpdx,dpdy,dpdz,BuoyancyCase,fx,fy,fz,gx,gy,gz,du,dv,dw, &
                   rho_ave)
              if(dosolids) then
                 du = du*umask; dv = dv*vmask; dw = dw*wmask
              endif
              call my_timer(2)



              if(Implicit) then    !! BOILVEL: Need to get/set here
                 call SetupUvel(u,du,rho,mu,rho1,mu1,dt,A)
                 if(hypre)then
                    call poi_solve(A,u,maxError,maxit,it,HYPRESolverType)
                 else
                    call LinearSolver1(A,u,umask,maxError,beta,maxit,itu,ierr)
                 endif
                 if(mod(itimestep,termout)==0) call calcresidual1(A,u,umask,residualu)
                 call SetupVvel(v,dv,rho,mu,rho1,mu1,dt,A)
                 if(hypre)then
                    call poi_solve(A,v,maxError,maxit,it,HYPRESolverType)
                 else
                    call LinearSolver1(A,v,vmask,maxError,beta,maxit,itv,ierr)
                 endif
                 if(mod(itimestep,termout)==0) call calcresidual1(A,v,vmask,residualv)
                 call SetupWvel(w,dw,rho,mu,rho1,mu1,dt,A)
                 if(hypre)then
                    call poi_solve(A,w,maxError,maxit,it,HYPRESolverType)
                 else
                    call LinearSolver1(A,w,wmask,maxError,beta,maxit,itw,ierr)
                 endif
                 if(mod(itimestep,termout)==0) call calcresidual1(A,w,wmask,residualw)
              else
                 u = u + dt * du
                 v = v + dt * dv
                 w = w + dt * dw
              endif
           else
              do i=is,ie; do j=js,je; do k=ks,ke
                 if (u_cmask(i,j,k)==0) u(i,j,k) = u(i,j,k) + dt*du(i,j,k) 
                 if (v_cmask(i,j,k)==0) v(i,j,k) = v(i,j,k) + dt*dv(i,j,k) 
                 if (w_cmask(i,j,k)==0) w(i,j,k) = w(i,j,k) + dt*dw(i,j,k) 
              enddo; enddo; enddo
           endif
           if(out_sub .and. ii==1 .and. mod(itimestep-itimestepRestart,nout)==0) then
              call output5(itimestep/nout,itimestep,1)
           endif

           call my_timer(3)
           call SetVelocityBC(u,v,w,umask,vmask,wmask,time,dt,0)
           if (FreeSurface) then
              if (inflow .and. (itimestep<=step_max .or. fs_refactor)) call inflow_accelerate  
              call my_timer(15)
           endif
           call do_ghost_vector(u,v,w)
           call my_timer(1)
           if (DoVof .and. debug_par) then
              call write_par_var("U-vel1    ",itimestep,u)
              call write_par_var("V-vel1    ",itimestep,v)
              call write_par_var("W-vel1    ",itimestep,w)
           endif
!-----------------------------------------PROJECTION STEP-----------------------------------------
           call SetPressureBC(umask,vmask,wmask)
           if (.not.FreeSurface) then
             if (STGhost) then
              call SetupPoissonGhost(u,v,w,umask,vmask,wmask,dt,A,tmp,VolumeSource)
             else
              call SetupPoisson(u,v,w,umask,vmask,wmask,rho,dt,A,tmp,VolumeSource)
            endif
           else
              if (FS_Hypre) then
                 call get_all_curvatures(tmp,itimestep)
                 call setuppoisson_fs_hypre(u,v,w,rho1,dt,A,height,tmp,tag_id)
              else
                 solver_flag = 1
                 call setuppoisson_fs_heights(u,v,w,rho1,dt,A,tag_id)
              endif
           endif
           ! (div u)*dt < epsilon => div u < epsilon/dt => maxresidual : maxerror/dt 
           if(HYPRE)then
              tmp = p
              call poi_solve(A,p,maxError/MaxDt*ErrorScaleHYPRE,maxit,it,HYPRESolverType)
              call do_all_ghost(p)
           else
              if (FreeSurface) then
                 if (.not.FS_Hypre) then
                    if (RP_test) call set_RP_pressure(p,rho1)
                    call FreeSolver(A,p,maxError/MaxDt,beta,maxit,it,ierr,itimestep,time,residual)
                 endif
              else
                 call NewSolver(A,p,maxError/MaxDt,beta,maxit,it,ierr,ResNormOrderPressure)
              endif
           endif
           ! If maximum iteration is reached, check residual
           if(it==maxit) then 
              call calcResiduals(A,p,res_tab)
              residual=res_tab(max(ResNormOrderPressure,3))
              if (residual/(maxError/MaxDt) > DivergeTol) then 
                 if ( SwitchHYPRESolver .and. HYPRESolverType == 2 ) then
                    ! if HYPRE-PFMG is diverged, retry with SMG
                    if (rank == 0) write(*,*) "PFMG failed, switch to SMG...",itimestep,ii,residual/(maxError/MaxDt),DivergeTol
                    HYPRESolverType = 1
                    if (residual/(maxError/MaxDt) > DivergeTol) p = tmp !PFMG may diverge,then use old p as IC for iteration
                    call poi_solve(A,p,maxError/MaxDt*ErrorScaleHYPRE,maxit,it,HYPRESolverType)
                    call do_all_ghost(p)
                    call calcResidual(A,p,ResNormOrderPressure,residual)
                    if (residual/(maxError/MaxDt) > DivergeTol ) then  
                       if (rank == 0) write(*,*) "PFMG/SMG pressure solver diverge!",& 
                           itimestep,ii,residual/(maxError/MaxDt),DivergeTol,HYPRESolverType
                       call pariserror("PFMG/SMG pressure solver diverge!")
                    end if ! residual
                 else 
                    call pariserror("Pressure solver diverge! (Adjust DivergeTol or Reduce ErrorScaleHYPRE if HYPRE is used!")
                 end if ! SwitchHYPRESolver
              end if !residual
           end if !it
           if(mod(itimestep,termout)==0 .and. ii==1) then
              !if (.not.FreeSurface) call calcResidual(A,p,ResNormOrderPressure,residual)
              call calcResiduals(A,p,res_tab)
              residual=res_tab(max(ResNormOrderPressure,3))
              if(rank==0) then
                 write(*,'("  pressure residuals*dt L1:",e8.1,"         L2:",e8.1,"       Linf:",e8.1)') &
                      res_tab*MaxDt
                 write(*,'("  pressure iterations     :",I8,  " tolerance :",e8.1," norm order:",I4)')   &
                      it,maxerror,ResNormOrderPressure              
              end if
              if ( SwitchHYPRESolver .and. HYPRESolverType == 1 ) then
                 HYPRESolverType = 2 
                 if (rank == 0) write(*,*) "Switch from HYPRE-SMG to PFMG"
              end if ! SwitchSHYPRESolver
              if ( DynamicAdjustPoiTol ) then 
                 if (HYPRE .and. residual/(maxError/MaxDt) > 2.d0 ) then  
                    ErrorScaleHYPRE = ErrorScaleHYPRE*0.5d0
                    if (rank == 0) write(*,*) "ErrorScaleHYPRE is decreased.",ErrorScaleHYPRE
! TEMPORARY 
                    call poi_solve(A,p,maxError/MaxDt*ErrorScaleHYPRE,maxit,it,HYPRESolverType)
                    call do_all_ghost(p)
                    call calcResidual(A,p,ResNormOrderPressure,residual)
                    if(rank==0) then
                       write(*,'("Solve again:  pressure residual*dt:   ",e8.1,&
                         &" maxerror: ",e8.1)') residual*MaxDt,maxerror
                    end if
! END TEMPORARY
                 else if (HYPRE .and. residual/(maxError/MaxDt) < 0.5d0 ) then  
                    ErrorScaleHYPRE = ErrorScaleHYPRE*2.0d0
                    if (rank == 0) write(*,*) "ErrorScaleHYPRE is increased.",ErrorScaleHYPRE
                 end if ! HYPRE & residual
              end if ! DynamicAdjustPoiTol  
           endif

           if(out_sub .and. ii==1 .and. mod(itimestep-itimestepRestart,nout)==0) then
              call output5(itimestep/nout,itimestep,2)
           endif

           call project_velocity()

           if(out_sub .and. ii==1 .and. mod(itimestep-itimestepRestart,nout)==0) then
              call output5(itimestep/nout,itimestep,3)
           endif
           if(mod(itimestep,nout)==0) call check_corrected_vel(u,umask,itimestep)
           if( DoLPP ) call ComputeSubDerivativeVel()
           call my_timer(10)
           !--------------------------------------UPDATE COLOR---------------------------------------------
           if (DoFront) then
                 call SetupDensity(dIdx,dIdy,dIdz,A,color)
                 if(hypre)then
                    call poi_solve(A,color,maxError,maxit,it,HYPRESolverType)
                 else
                    call NewSolver(A,color,maxError,beta,maxit,it,ierr,ResNormOrderPressure)
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
           
           call my_timer(13)
           call SetVelocityBC(u,v,w,umask,vmask,wmask,time,dt,1)
           call do_ghost_vector(u,v,w)
           call do_all_ghost(color)
           call my_timer(1)        

!----------------------------------EXTRAPOLATION FOR FREE SURFACE---------------------------------
           if (DoVOF .and. FreeSurface) then
              call extrapolate_velocities()
              
              if (do_2nd_projection) then

                 if(out_sub .and. ii==1 .and. mod(itimestep-itimestepRestart,nout)==0) then
                    call output5(itimestep/nout,itimestep,4)
                 endif
                 solver_flag = 2
                 call setuppoisson_fs_heights(u,v,w,rho1,dt,A,tag_id)
                 call FreeSolver(A,p_ext,5.0d-2,beta,10,it,ierr,itimestep,time,residual)
                 if(mod(itimestep,termout)==0) then
                    if(rank==0) then
                       write(*,'("FS2:   pressure residual, L_inf:   ",e8.1,&
                            &" maxerror: ",e8.1)') residual,5.0d-2
                       write(*,'("              pressure iterations :",I9)')it
                    endif
                 endif
                 ! Correct ONLY masked gas velocities at level 1 and 2
                 do k=ks,ke;  do j=js,je; do i=is,ieu    ! CORRECT THE u-velocity 
                    if (u_cmask(i,j,k)==1 .or. u_cmask(i,j,k)==2) then
                       u(i,j,k)=u(i,j,k)-(p_ext(i+1,j,k)-p_ext(i,j,k))/dxh(i)
                    else if (u_cmask(i,j,k)==3) then
                       u(i,j,k) = 0d0
                    endif
                 enddo; enddo; enddo

                 do k=ks,ke;  do j=js,jev; do i=is,ie    ! CORRECT THE v-velocity
                    if (v_cmask(i,j,k)==1 .or. v_cmask(i,j,k)==2) then
                       v(i,j,k)=v(i,j,k)-(p_ext(i,j+1,k)-p_ext(i,j,k))/dyh(j)
                    else if (v_cmask(i,j,k)==3) then
                       v(i,j,k) = 0d0
                    endif
                 enddo; enddo; enddo

                 do k=ks,kew;  do j=js,je; do i=is,ie   ! CORRECT THE w-velocity
                    if (w_cmask(i,j,k)==1 .or. w_cmask(i,j,k)==2) then
                       w(i,j,k)=w(i,j,k)-(p_ext(i,j,k+1)-p_ext(i,j,k))/dzh(k)
                    else if (w_cmask(i,j,k)==3) then
                       w(i,j,k) = 0d0
                    endif
                 enddo; enddo; enddo
                 call SetVelocityBC(u,v,w,umask,vmask,wmask,time,dt,0) !check this
                 call do_ghost_vector(u,v,w)
              endif
              if (mod(itimestep,nstats)==0 .and. mod(ii,itime_scheme)==0) call discrete_divergence(u,v,w,itimestep/nstats)
              call my_timer(15) 
           endif !Extrapolation

!------------------------------------------------------------------------------------------------


!--------------------------------------UPDATE DENSITY/VISCOSITY------------------------------------
           if(TwoPhase) then
              if(GetPropertiesFromFront) then
                 rho = rho2 + (rho1-rho2)*color
                 mu  = mu2  + (mu1 -mu2 )*color
              else
!------------------------------------deduce rho, mu from cvof-------------------------------------
                 call linfunc(rho,rho1,rho2,DensMean)
                 call linfunc(mu,mu1,mu2,ViscMean)
!------------------------------------END VOF STUFF------------------------------------------------
              endif
           endif
           call my_timer(2)

           ! Wait for front to be sent back
           if(DoFront)call GetFront('wait')
           call my_timer(13)
        enddo !itime_scheme
        if(itime_scheme==2) then
           u = 0.5*(u+uold)
           v = 0.5*(v+vold)
           w = 0.5*(w+wold)
           rho = 0.5*(rho+rhoo)
           mu  = 0.5*(mu +muold)
#ifdef PHASE_CHANGE
           Te = 0.5*(Te+Te_old)
           Cp = 0.5*(Cp+Cp_old)
           kc = 0.5*(kc+kc_old)
#endif
           call my_timer(2)
           if(DoVOF) cvof  = 0.5*(cvof +cvofold)
           if(DoVOF) then
              call get_flags_and_clip(cvof,vof_flag)
              call get_vof_phase(cvof) !cvof updated above from min to max
           endif
           call my_timer(4)
           if ( DoLPP ) call AveragePartSol()
           call my_timer(12)
        endif
        if (FreeSurface) then
           if (RP_test) call Integrate_RP(dt,time,rho1)
           if (mod(itimestep,nstats)==0 .and. curve_stats) then
              call curvature_sphere(time)
           endif
        endif
        if (DoLPP) then
           call PartBCWrapper
           call lppvofsweeps(itimestep,time)  
           call SeedParticles
           call my_timer(14)
        end if ! DoLPP
!--------------------------------------------OUTPUT-----------------------------------------------
        if(mod(itimestep,nstats)==0) then
           call calcStats
           if(test_control_droplet) call do_droplet_test(itimestep,time,REAL(nstats*dt,8))
           if ( DoTurbStats .and. time > timeStartTurbStats ) call calcTurbStats(itimestep)
           if( test_KHI2D .or. test_HF ) call h_of_KHI2D(itimestep,time)
           call pressure_stats(time)
        endif
        if(mod(itimestep,nsteps_probe)==0) then
           call probes
        endif
        call my_timer(2)
        if(mod(itimestep,nbackup)==0) then 
           if ( DoFront ) then 
              call backup_write
           else if ( DoVOF ) then 
              call backup_VOF_write
              if ( DoLPP ) call backup_LPP_write 
              call wrap_up_timer(itimestep,iTimeStepRestart)
           end if ! DoFront, DoVOF
        end if ! itimestep
        !        if(mod(itimestep,noutuv)==0) then
        if ( tout > 0.d0 ) then  ! output based on time 
           if ( (tout - (time-timeLastOutput)) < 0.9999*dt ) then 
              nfile = NINT(time/tout)
! Standard outputs
              if(test_frdroplet.or.test_droplet) call output_droplet(w,time)
              call output(nfile,is,ie+1,js,je+1,ks,ke+1)
              if (DoVOF) then
                 call output_VOF(nfile,is,ie+1,js,je+1,ks,ke+1)
                 call output_ALL(nfile,is,ie+1,js,je+1,ks,ke+1,itimestep,5)
                 if(out_centroid) call output_centroids(nfile)
              endif
              if(DoLPP .and. DoOutputLPP) call output_LPP(nfile,is,ie+1,js,je+1,ks,ke+1)
              if(rank==0)then
                 end_time =  MPI_WTIME()
                 write(out,'("Output data at Step:",I9," and time:",f10.2)')itimestep,time
              endif
! End standard outputs ? 
              timeLastOutput = time
           end if! timeAfterOuput
        else                     ! output based on timestep
           if(mod(itimestep-itimestepRestart,nout)==0) then 
              nfile = ITIMESTEP/nout
              if ( tout > 0.d0 .and. dtFlag == 1 ) nfile = NINT(time/tout)
! Standard outputs
              if(test_frdroplet.or.test_droplet) call output_droplet(w,time)
              call write_vec_gnuplot(u,v,cvof,p,itimestep,DoVOF)
              call output(nfile,is,ie+1,js,je+1,ks,ke+1)
              if (DoVOF) then
                 call output_VOF(nfile,is,ie+1,js,je+1,ks,ke+1)
                 call output_ALL(nfile,is,ie+1,js,je+1,ks,ke+1,itimestep,5)
                 if(out_centroid) call output_centroids(nfile)
              endif
              call print_st_stats(nfile)
              if(DoLPP .and. DoOutputLPP) call output_LPP(nfile,is,ie+1,js,je+1,ks,ke+1)
              if(rank==0)then
                 end_time =  MPI_WTIME()
                 write(out,'("Step:",I9," Iterations:",I9," cpu(s):",f10.2)')itimestep,it,end_time-start_time
              endif
! End standard outputs ? 
              if (DoVOF .and. debug_par) then
                 !call get_all_curvatures(tmp,nfile)
                 !call get_all_heights(nfile)
                 call write_par_var("VOF       ",nfile,cvof)           
                 if (FreeSurface) then
                    call write_par_var("P_gx      ",nfile,P_gx)
                    call write_par_var("P_gy      ",nfile,P_gy)
                    call write_par_var("P_gz      ",nfile,P_gz)
                 endif
              endif
           endif
        end if ! tout
        if(nstats==0) call pariserror(" *** Main: nstats = 0")
        if(mod(itimestep,nstats)==0.and.rank==0)then
           !        open(unit=121,file='track')
           !        write(121,'("Step:",I10," dt=",es16.5e2," time=",es16.5e2)')itimestep,dt,time
           !        write(121,'("            Iterations:",I8," cpu(s):",f10.2)')it,end_time-start_time
              !        close(121)
           open(unit=121,file='stats',position='append')
           write(121,'(30es14.6e2)')time,stats(1:nstatarray),dpdx,(stats(8)-stats(9))/dt,end_time-start_time
           close(121)
        endif
        !output for scalar variables used in free surface
        if (FreeSurface) then
           if (RP_test .and. (mod(itimestep,nstats)==0) .and. rank==0) call write_RP_test(time,rho1)
           do out_fs = 1,3
              if (VTK_OUT(out_fs)) then
                 if (mod(itimestep,NOUT_VTK(out_fs))==0) then
                    if (rank==0) call append_visit_fs(out_fs,itimestep/NOUT_VTK(out_fs))
                    SELECT CASE (out_fs)
                    case (1) 
                       tmp = 0d0
                       do k=ks,ke+1; do j=js,je+1; do i=is,ie+1
                          tmp(i,j,k)=ABS((u(i-1,j,k)-u(i,j,k))/dx(i)+(v(i,j-1,k)-v(i,j,k))/dy(j)+(w(i,j,k-1)-w(i,j,k))/dz(k))
                       enddo; enddo; enddo
                       call VTK_scalar_struct(out_fs,itimestep/NOUT_VTK(out_fs),tmp)
                    case(2)
                       call get_all_curvatures(tmp,itimestep)
                       call VTK_scalar_struct(out_fs,itimestep/NOUT_VTK(out_fs),tmp)
                    case(3)
                       call VTK_scalar_struct(out_fs,itimestep/NOUT_VTK(out_fs),v_source)
                    end SELECT
                 endif !output interval
              endif ! VTK_OUT(fs)
           enddo
        endif !FreeSurface
        call my_timer(11)
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
        if(dtFlag==2)call TimeStepSize(dt,vof_phase)
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
           open(unit=121,file='statsbub',position='append')
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
  call wrap_up_timer(itimestep,iTimeStepRestart)
  if(rank==0) then 
     if(output_format==2) call close_visit_file()
  endif

  if(rank<nPdomain)  call output_at_location()
  if(rank==0)  call final_output(stats(2))
  if(HYPRE) call poi_finalize
  if (MultiGrid) call finalize_MG
  if(rank==0) write(*,'("Paris exits succesfully")')
#ifdef PHASE_CHANGE
  call temp_stats
#endif
  call print_st_stats(itimestep/nout+1)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call MPI_FINALIZE(ierr)
  stop
end program paris
!=================================================================================================
subroutine project_velocity ()

  use module_flow
  use module_grid
  use module_freesurface
  use module_vof
  use module_surface_tension
  implicit none
  integer :: i,j,k

  if (.not.FreeSurface) then
     if (STGhost) then
        call project_velocity_staggered(u,umask,dt,p,1)
        call project_velocity_staggered(v,vmask,dt,p,2)
        call project_velocity_staggered(w,wmask,dt,p,3)
     else
        do k=ks,ke;  do j=js,je; do i=is,ieu    ! CORRECT THE u-velocity 
           u(i,j,k)=u(i,j,k)-dt*(2.0*umask(i,j,k)/dxh(i))*(p(i+1,j,k)-p(i,j,k))/(rho(i+1,j,k)+rho(i,j,k))
        enddo; enddo; enddo

        do k=ks,ke;  do j=js,jev; do i=is,ie    ! CORRECT THE v-velocity
           v(i,j,k)=v(i,j,k)-dt*(2.0*vmask(i,j,k)/dyh(j))*(p(i,j+1,k)-p(i,j,k))/(rho(i,j+1,k)+rho(i,j,k))
        enddo; enddo; enddo

        do k=ks,kew;  do j=js,je; do i=is,ie   ! CORRECT THE w-velocity
           w(i,j,k)=w(i,j,k)-dt*(2.0*wmask(i,j,k)/dzh(k))*(p(i,j,k+1)-p(i,j,k))/(rho(i,j,k+1)+rho(i,j,k))
        enddo; enddo; enddo
     endif
  else
     do k=ks,ke;  do j=js,je; do i=is,ieu    ! CORRECT THE u-velocity 
        if (u_cmask(i,j,k)==0) then
           u(i,j,k)=u(i,j,k)-dt/rho(i,j,k)*(p(i+1,j,k)+P_gx(i+1,j,k)-p(i,j,k)-P_gx(i,j,k))/x_mod(i,j,k)
           if (u(i,j,k) /= u(i,j,k)) write(*,'("WARNING u NaN :",2e14.5)')u(i,j,k), x_mod(i,j,k)
        endif
     enddo; enddo; enddo

     do k=ks,ke;  do j=js,jev; do i=is,ie    ! CORRECT THE v-velocity
        if (v_cmask(i,j,k)==0) then
           v(i,j,k)=v(i,j,k)-dt/rho(i,j,k)*(p(i,j+1,k)+P_gy(i,j+1,k)-p(i,j,k)-P_gy(i,j,k))/y_mod(i,j,k)
           if (v(i,j,k) /= v(i,j,k)) write(*,'("WARNING v NaN :",2e14.5)')v(i,j,k), y_mod(i,j,k)
        endif
     enddo; enddo; enddo

     do k=ks,kew;  do j=js,je; do i=is,ie   ! CORRECT THE w-velocity
        if (w_cmask(i,j,k)==0) then
           w(i,j,k)=w(i,j,k)-dt/rho(i,j,k)*(p(i,j,k+1)+P_gz(i,j,k+1)-p(i,j,k)-P_gz(i,j,k))/z_mod(i,j,k)
           if (w(i,j,k) /= w(i,j,k)) write(*,'("WARNING w NaN :",2e14.5)')w(i,j,k), z_mod(i,j,k)
        endif
     enddo; enddo; enddo
     if (RP_test) then
        if (coords(1)==0) then
           do k=ks,ke;  do j=js,je
              u(is-1,j,k)=u(is-1,j,k)-dt/rho(is-1,j,k)*(p(is,j,k)-p(is-1,j,k))/dx(is)
           enddo; enddo
        endif
        if (coords(1)==nPx-1) then
           do k=ks,ke;  do j=js,je
              u(ie,j,k)=u(ie,j,k)-dt/rho(ie,j,k)*(p(ie+1,j,k)-p(ie,j,k))/dx(ie)
           enddo; enddo
        endif
        if (coords(2)==0) then
           do k=ks,ke;  do i=is,ie
              v(i,js-1,k)=v(i,js-1,k)-dt/rho(i,js-1,k)*(p(i,js,k)-p(i,js-1,k))/dy(js)
           enddo; enddo
        endif
        if (coords(2)==nPy-1) then
           do k=ks,ke;  do i=is,ie
              v(i,je,k)=v(i,je,k)-dt/rho(i,je,k)*(p(i,je+1,k)-p(i,je,k))/dy(je)
           enddo; enddo
        endif
        if (coords(3)==0) then
           do j=js,je;  do i=is,ie
              w(i,j,ks-1)=w(i,j,ks-1)-dt/rho(i,j,ks-1)*(p(i,j,ks)-p(i,j,ks-1))/dz(ks)
           enddo; enddo
        endif
        if (coords(3)==nPz-1) then
           do j=js,je;  do i=is,ie
              w(i,j,ke)=w(i,j,ke)-dt/rho(i,j,ke)*(p(i,j,ke+1)-p(i,j,ke))/dz(ke)
           enddo; enddo
        endif
     endif
  endif
end subroutine project_velocity
!=================================================================================================
! subroutine TimeStepSize
!-------------------------------------------------------------------------------------------------
subroutine TimeStepSize(deltaT,vof_phase1)
  use module_grid
  use module_flow
  use module_2phase
  use module_IO
  use module_BC
  use module_freesurface
  use module_VOF
#ifdef PHASE_CHANGE
  use module_boil
#endif
  implicit none
  include "mpif.h"
  integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: vof_phase1
  real(8) :: deltaT, h, vmax, dtadv, mydt
  real(8) :: vmax_phys
  real(8) :: v_local
  integer :: ierr
  integer :: i,j,k
#ifdef PHASE_CHANGE
  real(8) :: dtTdiff
#endif

  vmax_phys = max(1.d-3,max(ugas_inject,uliq_inject))
  if(rank<nPdomain)then
     h  = minval(dx)
     if (.not.FreeSurface) then
        if ( inject_type == 3 ) then 
           vmax_phys = max(1.d-3,max(ugas_inject,uliq_inject))
           vmax      = 0.d0
           do k=ks,ke; do j=js,je; do i=is,ie
              v_local = sqrt(u(i,j,k)**2 + v(i,j,k)**2 + w(i,j,k)**2)
              vmax = MAX(vmax,v_local)
              if ( vmax > vmax_phys*1.d2 ) then
                  OPEN(UNIT=88,FILE=TRIM(out_path)//'/message-rank-'//TRIM(int2text(rank,padding))//'.txt')
                  write(88,*) "Error:Max velocity 100 times larger than physical value",itimestep,rank, vmax,vmax_phys
                  write(88,*) rank, i,j,k,itimestep 
                  write(88,*) u(i-1:i+1,j,k)
                  write(88,*) v(i,j-1:j+1,k)
                  write(88,*) w(i,j,k-1:k+1)
                  write(88,*) cvof(i-1:i+1,j,k)
                  write(88,*) cvof(i,j-1:j+1,k)
                  write(88,*) cvof(i,j,k-1:k+1)
                  CLOSE(88)
                  call pariserror("Max velocity 100 times larger than physical value, something wrong!") 
              end if ! vmax
           enddo; enddo; enddo ! i,j,k 
         else 
           vmax = maxval(sqrt(u(is:ie,js:je,ks:ke)**2+v(is:ie,js:je,ks:ke)**2+w(is:ie,js:je,ks:ke)**2))
        end if ! inject_type 
              
     else
        vmax_phys = 0.d0
        vmax = 0.d0
        do k=ks,ke; do j=js,je; do i=is,ie
           if (vof_phase1(i,j,k)==0) then
              v_local = sqrt(u(i,j,k)**2 + v(i,j,k)**2 + w(i,j,k)**2)
              vmax = MAX(vmax,v_local)
           endif
        enddo; enddo; enddo 
     endif ! FreeSurface
     !print *,"h=",h,"maxv=",vmax,vmax_phys,'freesurf=',FreeSurface     
     dtadv = h/(max(vmax_phys,vmax))
     mydt  = CFL*dtadv
     mydt = min(mydt,MaxDt)
#ifdef PHASE_CHANGE
     dtTdiff = h**2.d0*min(rho1*Cp1/kc1,rho2*Cp2/kc2)/2.d0
     mydt = min(mydt,dtTdiff)
#endif
  else
     mydt=1.0e10
  endif ! rank
  call MPI_ALLREDUCE(mydt, deltat, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_Active, ierr)


end subroutine TimeStepSize
!=================================================================================================
function get_cfl_and_check(deltaT)
  use module_grid
  use module_flow
  use module_VOF
  use module_freesurface
  implicit none
  include 'mpif.h'
  integer :: ierr
  integer :: i,j,k
  real(8) :: get_cfl_and_check,vmax,inbox_cfl,deltaT,h,glogcfl
  real(8) :: v_local
  
  if (.not. Freesurface) then
     vmax = maxval(sqrt(u(is:ie,js:je,ks:ke)**2 + v(is:ie,js:je,ks:ke)**2 + w(is:ie,js:je,ks:ke)**2))
  else
     vmax = 0.d0
     do k=ks,ke; do j=js,je; do i=is,ie
        if (vof_phase(i,j,k)==0) then
           v_local = sqrt(u(i,j,k)**2 + v(i,j,k)**2 + w(i,j,k)**2)
           vmax = MAX(vmax,v_local)
        endif
     enddo; enddo; enddo 
  endif
  h  = minval(dx)
  inbox_cfl=vmax*deltaT/h
  call MPI_ALLREDUCE(inbox_cfl, glogcfl, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_Cart,ierr) 
  get_cfl_and_check=glogcfl
  if(get_cfl_and_check.gt.cflmax_allowed) call pariserror("CFL too large")
end function get_cfl_and_check


subroutine interface_max_min

  use module_grid
  use module_flow
  use module_VOF
  use module_BC
  implicit none
  include "mpif.h"
  integer :: i,j,k,ierr
  real(8) :: int_max_y
  
  int_max_y = -10000000

  call do_all_ghost(cvof)

  !checkme: why does it get the border if I search until je?
  do k=ks,ke;  do j=js,je-10;  do i=is,ie
    if (abs(cvof(i,j,k)-cvof(i,j+1,k)).gt.1d-10) then
      int_max_y = max(int_max_y, y(j)+dy(j)*(cvof(i,j,k) - 0.5d0))
    endif 
  enddo; enddo; enddo;

  call MPI_ALLREDUCE(int_max_y, int_max_y, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_Domain, ierr)

  if (rank==0) then
    if (time.eq.0) then
      open(unit=122,file='interface.dat')
    else
      open(unit=122,file='interface.dat',position='append')
    endif

    write(122,*) time, int_max_y
    close(122)
  endif

end subroutine interface_max_min

!=================================================================================================
! subroutine calcStats
!-------------------------------------------------------------------------------------------------
subroutine calcStats
  use module_grid
  use module_flow
  use module_VOF
#ifdef PHASE_CHANGE
  use module_boil
#endif
  implicit none
  include "mpif.h"
  integer :: i,j,k,ierr
  real(8) :: vol,CC=0d0
  real(8) :: kenergy,vort2,enstrophy
  real(8), save :: W_int=-0.02066
  real(8), save :: height4Stats=0.875d0
  integer :: i0,j0,k0
  real(8) :: stencil3x3(-1:1,-1:1,-1:1),mxyz(3),AREA3D

  nstatarray=25
  if(nstatarray > 100) call pariserror("nstatarray too large")
  mystats(1:nstatarray)=0d0
  do k=ks,ke;  do j=js,je;  do i=is,ie
     vol = dx(i)*dy(j)*dz(k)
     ! Average u component
     mystats(2)=mystats(2)+u(i,j,k)*vol
     mystats(3)=mystats(3)+fx(i,j,k)*dxh(i)*dy(j)*dz(k)
     mystats(4)=mystats(4)+fy(i,j,k)*dx(i)*dyh(j)*dz(k)
     mystats(5)=mystats(5)+fz(i,j,k)*dx(i)*dy(j)*dzh(k)
     mystats(6)=mystats(6)+rho(i,j,k)*vol
     mystats(7)=mystats(7)+p(i,j,k)*vol
     ! average momentum
     mystats(8)=mystats(8)+0.5*(rho(i,j,k)+rho(i+1,j,k))*u(i,j,k)*vol
     mystats(9)=mystats(9)+0.5*(rho(i,j,k)+rho(i+1,j,k))*uold(i,j,k)*vol
     
     if(DoVOF) then
        CC=cvof(i,j,k)
        ! Phase C=1 volume
        mystats(10)=mystats(10)+CC*vol
        ! Phase C=1 center of mass x-coordinate
        mystats(11)=mystats(11)+CC*vol*x(i)
        ! Rise velocity of bubble with gravity in x-direction
        if (test_risingbubble) then
           mystats(24) = mystats(24)+0.50d0*(u(i-1,j,k)+u(i,j,k))*CC*vol
        endif
     endif
     ! kinetic energy
     kenergy = 0.5d0*( (0.5d0*(u(i,j,k)+u(i+1,j,k)))**2.d0 & 
          +(0.5d0*(v(i,j,k)+v(i,j+1,k)))**2.d0 & 
          +(0.5d0*(w(i,j,k)+w(i,j,k+1)))**2.d0)
     if(DoVOF) then
        mystats(12)=mystats(12)+rho(i,j,k)*kenergy*vol*      cvof(i,j,k)
        mystats(13)=mystats(13)+rho(i,j,k)*kenergy*vol*(1.d0-cvof(i,j,k))
     else 
        mystats(13)=mystats(13)+rho(i,j,k)*kenergy*vol
     end if ! DoVOF
     ! y-momentum 
     if(DoVOF) then
      if (test_turb) then
        mystats(14)=mystats(14)+((-u(i,j,k)+u(i+1,j,k))/dx(i))**2*vol
        mystats(15)=mystats(15)+((-v(i,j,k)+v(i,j+1,k))/dy(j))**2*vol
      else
        mystats(14)=mystats(14)+rho(i,j,k)*0.5d0*(v(i,j,k)+v(i,j+1,k))*vol*(     cvof(i,j,k))
        mystats(15)=mystats(15)+rho(i,j,k)*0.5d0*(v(i,j,k)+v(i,j+1,k))*vol*(1.d0-cvof(i,j,k))
      endif
     end if ! (DoVOF)
     ! interfacial area
     if (DoVOF ) then
         ! MODEMI version of interface area 
!        if (     max((cvof(i+1,j,k)-0.5d0)/abs(cvof(i+1,j,k)-0.5d0),0.d0) & 
!             - max((cvof(i-1,j,k)-0.5d0)/abs(cvof(i-1,j,k)-0.5d0),0.d0) /= 0.d0 &
!             .or. max((cvof(i,j+1,k)-0.5d0)/abs(cvof(i,j+1,k)-0.5d0),0.d0) & 
!             - max((cvof(i,j-1,k)-0.5d0)/abs(cvof(i,j-1,k)-0.5d0),0.d0) /= 0.d0 &
!             .or. max((cvof(i,j,k+1)-0.5d0)/abs(cvof(i,j,k+1)-0.5d0),0.d0) & 
!             - max((cvof(i,j,k-1)-0.5d0)/abs(cvof(i,j,k-1)-0.5d0),0.d0) /= 0.d0 & 
!             ) then
!           if ( vof_flag(i,j,k) == 2 ) & 
!                mystats(19) = mystats(19) + dx(i)*dy(j)
!        end if ! cvof
         
         ! PLIC interface area 
         if ( cvof(i,j,k) > 0.d0 .and. cvof(i,j,k) < 1.d0 ) then 
            do i0=-1,1; do j0=-1,1; do k0=-1,1
               stencil3x3(i0,j0,k0) =cvof(i+i0,j+j0,k+k0)
            enddo;enddo;enddo
            call mycs(stencil3x3,mxyz)
            if ( mxyz(1)/=0.d0 .or. mxyz(2)/=0.d0 .or. mxyz(3)/=0.d0 ) &  
               mystats(19) = mystats(19) + AREA3D(mxyz,cvof(i,j,k))
         end if ! cvof(i,j,k)
     end if ! DoVOF
     ! potential energy (considering gravity in y direction)  
     if(DoVOF .and. test_PhaseInversion) then
        mystats(20)=mystats(20)+rho(i,j,k)*Gy*y(j)*vol*      cvof(i,j,k)
        mystats(21)=mystats(21)+rho(i,j,k)*Gy*y(j)*vol*(1.d0-cvof(i,j,k))
     end if ! DoVOF
     ! enstrophy
     if(DoVOF ) then
        vort2 = (0.5d0*(w(i,j+1,k)+w(i,j+1,k+1)-w(i,j-1,k)-w(i,j-1,k+1))/(y(j+1)-y(j-1)) & 
             -0.5d0*(v(i,j,k+1)+v(i,j+1,k+1)-v(i,j,k-1)-v(i,j+1,k-1))/(z(k+1)-z(k-1)))**2.d0 & 
             + (0.5d0*(u(i,j,k+1)+u(i+1,j,k+1)-u(i,j,k-1)-u(i+1,j,k-1))/(z(k+1)-z(k-1)) & 
             -0.5d0*(w(i+1,j,k)+w(i+1,j,k+1)-w(i-1,j,k)-w(i-1,j,k+1))/(x(i+1)-x(i-1)))**2.d0 & 
             + (0.5d0*(v(i+1,j,k)+v(i+1,j+1,k)-v(i-1,j,k)-v(i-1,j+1,k))/(x(i+1)-x(i-1)) & 
             -0.5d0*(u(i,j+1,k)+u(i+1,j+1,k)-u(i,j-1,k)-u(i+1,j-1,k))/(y(j+1)-y(j-1)))**2.d0 
        enstrophy = 0.5*vort2
        mystats(22)=mystats(22)+enstrophy*vol*      cvof(i,j,k)
        mystats(23)=mystats(23)+enstrophy*vol*(1.d0-cvof(i,j,k))
        if(vof_flag(i,j,k)==2) mystats(25)=mystats(25)+1
     end if ! DoVOF
     ! volume fraction of top layer
     if(DoVOF .and. test_PhaseInversion) then
        if ( y(j) > height4Stats ) & 
             mystats(24)=mystats(24)+vol*cvof(i,j,k)/(xLength*zLength*(yLength-height4Stats))
     end if ! DoVOF
  enddo;  enddo;  enddo

! Shear stress on y=0,Ly
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
! Get pressure at entry and exit
  if (test_point_in(ng+1,ny/2+ng,nz/2+ng)) then 
     mystats(17)=p(ng+1,ny/2+ng,nz/2+ng)
  else
     mystats(17)=0d0
  endif
  if (test_point_in(ng+nx,ny/2+ng,nz/2+ng)) then 
     mystats(18)=p(ng+nx,ny/2+ng,nz/2+ng)
  else
     mystats(18)=0d0
  endif
  mystats(2:11) = mystats(2:11)/(xLength*yLength*zLength)
  mystats(1) = mystats(1)/(xLength*zLength*2.0)
  mystats(11) = mystats(11) ! /mystats(10)
  call MPI_ALLREDUCE(mystats(1), stats(1), nstatarray, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_Domain, ierr)
  if (DoVOF .and. test_risingbubble) then
     stats(10:11) = stats(10:11)*(xLength*yLength*zLength)
     if (stats(10) > 1d-15) then
        stats(11) = stats(11)/stats(10)
        stats(24) = stats(24)/stats(10)
     endif
  endif
  rho_ave = stats(6)
  p_ave = stats(7)
! This stops the code in case the average velocity (stats(2)) becomes NaN.
  if(stats(2).ne.stats(2)) call pariserror("********** Invalid flow rate **********")
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
! subroutine probes
!-------------------------------------------------------------------------------------------------
   subroutine probes
      use module_flow
      use module_VOF

      integer :: iprobe,i,j,k
      character(len=30) :: filename 

      do iprobe = 1,num_probes
         i = ijk_probe(iprobe,1)
         j = ijk_probe(iprobe,2)
         k = ijk_probe(iprobe,3)
         if (  i >= is .and. i <= ie .and. & 
               j >= js .and. j <= je .and. & 
               k >= ks .and. k <= ke ) then 
            dat_probe(iprobe,1) = 0.5d0*(u(i,j,k)+u(i+1,j,k))
            dat_probe(iprobe,2) = 0.5d0*(v(i,j,k)+v(i,j+1,k))
            dat_probe(iprobe,3) = 0.5d0*(w(i,j,k)+w(i,j,k+1))
            dat_probe(iprobe,4) = cvof(i,j,k)
            dat_probe(iprobe,5) = p(i,j,k)
            OPEN(UNIT=300+iprobe,FILE=TRIM(out_path)//'/probe-'//TRIM(int2text(iprobe,padding))//'.dat',POSITION='append')
            write(300+iprobe,'(6(E23.16,1X))') time,dat_probe(iprobe,:)
            CLOSE(300+iprobe)
         end if       
      end do ! iprobe

      do iprobe = 1,num_probes_cvof
         i = ijk_probe_cvof(iprobe,1)
         j = ijk_probe_cvof(iprobe,2)
         k = ijk_probe_cvof(iprobe,3)
         if (  i >= is .and. i <= ie .and. & 
               j >= js .and. j <= je .and. & 
               k >= ks .and. k <= ke ) then 
            dat_probe_cvof(iprobe) = cvof(i,j,k)
            OPEN(UNIT=400+iprobe,FILE=TRIM(out_path)//'/probe_cvof-'//TRIM(int2text(iprobe,padding))//'.dat',POSITION='append')
            write(400+iprobe,'(2(E23.16,1X))') time,dat_probe_cvof(iprobe)
            CLOSE(400+iprobe)
         end if       
      end do ! iprobe

      do iprobe = 1,num_probes_linez ! in z direction only
         i = ij_probe_linez(iprobe,1)
         j = ij_probe_linez(iprobe,2)
         if (  i >= is .and. i <= ie .and. & 
               j >= js .and. j <= je ) then
            filename=TRIM(out_path)//'/probe_linez-'//TRIM(int2text(iprobe,padding))//'-'
            OPEN(UNIT=500+iprobe,FILE=TRIM(filename)//TRIM(int2text(rank,padding))//'.dat',POSITION='append')
            do k = ks,ke
               write(500+iprobe,'(7(E23.16,1X))') time,z(k),& 
                  u(i,j,k),v(i,j,k),w(i,j,k),cvof(i,j,k),p(i,j,k)
            end do !k 
            CLOSE(500+iprobe)
         end if ! i,j
      end do ! iprobe
   end subroutine probes

!=================================================================================================
!  subroutine calcTurbStats
!-------------------------------------------------------------------------------------------------
   subroutine calcTurbStats(istep)
      use module_flow
      use module_VOF
      use module_IO
      implicit none
      include "mpif.h"
      
      integer, intent(in) :: istep
      real(8) :: uc(is:ie,js:je,ks:ke),  vc(is:ie,js:je,ks:ke), & 
                 wc(is:ie,js:je,ks:ke),area(is:ie,js:je,ks:ke), & 
                 n1(is:ie,js:je,ks:ke),  n2(is:ie,js:je,ks:ke), & 
                 n3(is:ie,js:je,ks:ke) 
                 
      integer :: i0,j0,k0
      real(8) :: stencil3x3(-1:1,-1:1,-1:1),mxyz(3),AREA3D,nm
      real(8) :: dudx,dvdy,dwdz,dvdx,dwdx,dudy,dwdy,dudz,dvdz
      integer :: counts,cx,cy,i,j,k
      integer :: req(2),sta(MPI_STATUS_SIZE,2),ierr,irank,irank_sum
      character(len=30) :: filename 
      logical :: file_exist
      real(8), dimension(:,:,:,:), allocatable :: turb_vars_all !i,j,ivar,rank 
      real(8), dimension(:,:,:), allocatable :: turb_vars_map !i,j,ivar

      integer :: rootout 

      ! initialize
      if ( .not. calcTurbStats_initialized ) then
         calcTurbStats_initialized = .true.
         if ( TurbStatsOrder == 1 ) then 
            num_turb_vars = 6 
         else if ( TurbStatsOrder == 2 ) then 
            num_turb_vars = 27 
         else if ( TurbStatsOrder == 3 ) then
            num_turb_vars = 72
         else if ( TurbStatsOrder == 4 ) then
            num_turb_vars = 90
         else 
            call pariserror("Statistics of turbulence over 3rd order!")
         end if ! TurbStatsOrder
         allocate(turb_vars(is:ie,js:je,num_turb_vars))
         
         if ( restart ) then  
            call backup_turb_read
         else
            turb_vars = 0.d0 
            iSumTurbStats = 0 
         end if ! restart
      end if ! calcTurbStats_initialized

      ! simulation variables at cell centers
      uc=0.d0;vc=0.d0;wc=0.d0;area=0.d0;n1=0.d0;n2=0.d0;n3=0.d0
      do i= is,ie; do j=js,je; do k=ks,ke
         uc(i,j,k) = 0.5d0*(u(i,j,k) + u(i+1,j,k))
         vc(i,j,k) = 0.5d0*(v(i,j,k) + v(i,j+1,k))
         wc(i,j,k) = 0.5d0*(w(i,j,k) + w(i,j,k+1))
              
         ! calculate interface area
         if ( cvof(i,j,k) > 0.d0 .and. cvof(i,j,k) < 1.d0 ) then 
            do i0=-1,1; do j0=-1,1; do k0=-1,1
               stencil3x3(i0,j0,k0) =cvof(i+i0,j+j0,k+k0)
            enddo;enddo;enddo
            call mycs(stencil3x3,mxyz)
            if ( mxyz(1)/=0.d0 .or. mxyz(2)/=0.d0 .or. mxyz(3)/=0.d0 ) then  
               area(i,j,k) = AREA3D(mxyz,cvof(i,j,k))*dx(i)*dy(j)
               nm = mxyz(1)*mxyz(1) + mxyz(2)*mxyz(2) + mxyz(3)*mxyz(3)
               nm = sqrt(nm)
               n1(i,j,k) = abs(mxyz(1)/nm)
               n2(i,j,k) = abs(mxyz(2)/nm)
               n3(i,j,k) = abs(mxyz(3)/nm)
            end if ! mxyz
         end if ! cvof
      end do; end do; end do
      
      ! calculate turbulent variables and sum over z direction
      iSumTurbStats = iSumTurbStats + 1  
      do k = ks,ke
         ! 1st order stats (1-6)
         if (TurbStatsOrder >= 1 ) then 
         turb_vars(is:ie,js:je,1) = turb_vars(is:ie,js:je,1) + uc(is:ie,js:je,k)
         turb_vars(is:ie,js:je,2) = turb_vars(is:ie,js:je,2) + vc(is:ie,js:je,k)
         turb_vars(is:ie,js:je,3) = turb_vars(is:ie,js:je,3) + wc(is:ie,js:je,k)
         turb_vars(is:ie,js:je,4) = turb_vars(is:ie,js:je,4) + p(is:ie,js:je,k)
         turb_vars(is:ie,js:je,5) = turb_vars(is:ie,js:je,5) + cvof(is:ie,js:je,k)
         turb_vars(is:ie,js:je,6) = turb_vars(is:ie,js:je,6) + area(is:ie,js:je,k)
         end if ! TurStatsOrder 

         ! 2nd order stats (7-27)
         if ( TurbStatsOrder >= 2 ) then 
         turb_vars(is:ie,js:je, 7) = turb_vars(is:ie,js:je,7) & 
                                   + uc(is:ie,js:je,k)*uc(is:ie,js:je,k)
         turb_vars(is:ie,js:je, 8) = turb_vars(is:ie,js:je,8) & 
                                   + vc(is:ie,js:je,k)*vc(is:ie,js:je,k)
         turb_vars(is:ie,js:je, 9) = turb_vars(is:ie,js:je,9) & 
                                   + wc(is:ie,js:je,k)*wc(is:ie,js:je,k)
         turb_vars(is:ie,js:je,10) = turb_vars(is:ie,js:je,10) & 
                                   + uc(is:ie,js:je,k)*vc(is:ie,js:je,k)
         turb_vars(is:ie,js:je,11) = turb_vars(is:ie,js:je,11) & 
                                   + uc(is:ie,js:je,k)*wc(is:ie,js:je,k)
         turb_vars(is:ie,js:je,12) = turb_vars(is:ie,js:je,12) & 
                                   + vc(is:ie,js:je,k)*wc(is:ie,js:je,k)
         turb_vars(is:ie,js:je,13) = turb_vars(is:ie,js:je,13) & 
                                   +  p(is:ie,js:je,k)* p(is:ie,js:je,k)
         turb_vars(is:ie,js:je,14) = turb_vars(is:ie,js:je,14) & 
                                   +  p(is:ie,js:je,k)*uc(is:ie,js:je,k)
         turb_vars(is:ie,js:je,15) = turb_vars(is:ie,js:je,15) & 
                                   +  p(is:ie,js:je,k)*vc(is:ie,js:je,k)
         turb_vars(is:ie,js:je,16) = turb_vars(is:ie,js:je,16) & 
                                   +  p(is:ie,js:je,k)*wc(is:ie,js:je,k)
         turb_vars(is:ie,js:je,17) = turb_vars(is:ie,js:je,17) & 
                                   + cvof(is:ie,js:je,k)*cvof(is:ie,js:je,k)
         turb_vars(is:ie,js:je,18) = turb_vars(is:ie,js:je,18) & 
                                   + cvof(is:ie,js:je,k)*uc(is:ie,js:je,k)
         turb_vars(is:ie,js:je,19) = turb_vars(is:ie,js:je,19) & 
                                   + cvof(is:ie,js:je,k)*vc(is:ie,js:je,k)
         turb_vars(is:ie,js:je,20) = turb_vars(is:ie,js:je,20) & 
                                   + cvof(is:ie,js:je,k)*wc(is:ie,js:je,k)
         turb_vars(is:ie,js:je,21) = turb_vars(is:ie,js:je,21) & 
                                   + area(is:ie,js:je,k)*area(is:ie,js:je,k)
         turb_vars(is:ie,js:je,22) = turb_vars(is:ie,js:je,22) & 
                                   + area(is:ie,js:je,k)*uc(is:ie,js:je,k)
         turb_vars(is:ie,js:je,23) = turb_vars(is:ie,js:je,23) & 
                                   + area(is:ie,js:je,k)*vc(is:ie,js:je,k)
         turb_vars(is:ie,js:je,24) = turb_vars(is:ie,js:je,24) & 
                                   + area(is:ie,js:je,k)*wc(is:ie,js:je,k)
         turb_vars(is:ie,js:je,25) = turb_vars(is:ie,js:je,25) & 
                                   +  p(is:ie,js:je,k)*cvof(is:ie,js:je,k)
         turb_vars(is:ie,js:je,26) = turb_vars(is:ie,js:je,26) & 
                                   +  p(is:ie,js:je,k)*area(is:ie,js:je,k)
         turb_vars(is:ie,js:je,27) = turb_vars(is:ie,js:je,27) & 
                                   + cvof(is:ie,js:je,k)*area(is:ie,js:je,k)
         end if ! TurbStatsOrder

         ! 3rd order stats
         if ( TurbStatsOrder >= 3 ) then 
         turb_vars(is:ie,js:je,28) = turb_vars(is:ie,js:je,28) & 
                                   + uc(is:ie,js:je,k)*uc(is:ie,js:je,k)*cvof(is:ie,js:je,k)
         turb_vars(is:ie,js:je,29) = turb_vars(is:ie,js:je,29) & 
                                   + vc(is:ie,js:je,k)*vc(is:ie,js:je,k)*cvof(is:ie,js:je,k)
         turb_vars(is:ie,js:je,30) = turb_vars(is:ie,js:je,30) & 
                                   + wc(is:ie,js:je,k)*wc(is:ie,js:je,k)*cvof(is:ie,js:je,k)
         turb_vars(is:ie,js:je,31) = turb_vars(is:ie,js:je,31) & 
                                   + uc(is:ie,js:je,k)*vc(is:ie,js:je,k)*cvof(is:ie,js:je,k)
         turb_vars(is:ie,js:je,32) = turb_vars(is:ie,js:je,32) & 
                                   + uc(is:ie,js:je,k)*wc(is:ie,js:je,k)*cvof(is:ie,js:je,k)
         turb_vars(is:ie,js:je,33) = turb_vars(is:ie,js:je,33) & 
                                   + vc(is:ie,js:je,k)*wc(is:ie,js:je,k)*cvof(is:ie,js:je,k)

         do i=is,ie; do j=js,je
            dudx = (u(i,j,k) - u(i-1,j,k))/dx(i)
            dvdy = (v(i,j,k) - v(i,j-1,k))/dy(j)
            dwdz = (w(i,j,k) - w(i,j,k-1))/dz(k)

            dvdx = (v(i+1,j,k)+v(i+1,j-1,k)-v(i-1,j,k)-v(i-1,j-1,k))*0.25d0/dx(i)
            dwdx = (w(i+1,j,k)+w(i+1,j,k-1)-w(i-1,j,k)-w(i-1,j,k-1))*0.25d0/dx(i)
            dudy = (u(i,j+1,k)+u(i-1,j+1,k)-u(i,j-1,k)-u(i-1,j-1,k))*0.25d0/dy(j)
            dwdy = (w(i,j+1,k)+w(i,j+1,k-1)-w(i,j-1,k)-w(i,j-1,k-1))*0.25d0/dy(j)
            dudz = (u(i,j,k+1)+u(i-1,j,k+1)-u(i,j,k-1)-u(i-1,j,k-1))*0.25d0/dz(k)
            dvdz = (v(i,j,k+1)+v(i,j-1,k+1)-v(i,j,k-1)-v(i,j-1,k-1))*0.25d0/dz(k)
         
            turb_vars(i,j,34) = turb_vars(i,j,34) + dudx
            turb_vars(i,j,35) = turb_vars(i,j,35) + dvdy
            turb_vars(i,j,36) = turb_vars(i,j,36) + dwdz
            turb_vars(i,j,37) = turb_vars(i,j,37) + dvdx
            turb_vars(i,j,38) = turb_vars(i,j,38) + dwdx
            turb_vars(i,j,39) = turb_vars(i,j,39) + dudy
            turb_vars(i,j,40) = turb_vars(i,j,40) + dwdy
            turb_vars(i,j,41) = turb_vars(i,j,41) + dudz
            turb_vars(i,j,42) = turb_vars(i,j,42) + dvdz
            turb_vars(i,j,43) = turb_vars(i,j,43) + dudx*dudx
            turb_vars(i,j,44) = turb_vars(i,j,44) + dvdy*dvdy
            turb_vars(i,j,45) = turb_vars(i,j,45) + dwdz*dwdz
            turb_vars(i,j,46) = turb_vars(i,j,46) + dvdx*dvdx
            turb_vars(i,j,47) = turb_vars(i,j,47) + dwdx*dwdx
            turb_vars(i,j,48) = turb_vars(i,j,48) + dudy*dudy
            turb_vars(i,j,49) = turb_vars(i,j,49) + dwdy*dwdy
            turb_vars(i,j,50) = turb_vars(i,j,50) + dudz*dudz
            turb_vars(i,j,51) = turb_vars(i,j,51) + dvdz*dvdz
            turb_vars(i,j,52) = turb_vars(i,j,52) + dvdx*dudy
            turb_vars(i,j,53) = turb_vars(i,j,53) + dwdy*dvdz
            turb_vars(i,j,54) = turb_vars(i,j,54) + dudz*dwdx
         end do; end do !i,j
         turb_vars(is:ie,js:je,55) = turb_vars(is:ie,js:je,55) &
                                   + uc(is:ie,js:je,k)*uc(is:ie,js:je,k)*uc(is:ie,js:je,k)
         turb_vars(is:ie,js:je,56) = turb_vars(is:ie,js:je,56) &
                                   + uc(is:ie,js:je,k)*vc(is:ie,js:je,k)*vc(is:ie,js:je,k)
         turb_vars(is:ie,js:je,57) = turb_vars(is:ie,js:je,57) &
                                   + uc(is:ie,js:je,k)*wc(is:ie,js:je,k)*wc(is:ie,js:je,k)
         turb_vars(is:ie,js:je,58) = turb_vars(is:ie,js:je,58) &
                                   + vc(is:ie,js:je,k)*uc(is:ie,js:je,k)*uc(is:ie,js:je,k)
         turb_vars(is:ie,js:je,59) = turb_vars(is:ie,js:je,59) &
                                   + vc(is:ie,js:je,k)*vc(is:ie,js:je,k)*vc(is:ie,js:je,k)
         turb_vars(is:ie,js:je,60) = turb_vars(is:ie,js:je,60) &
                                   + vc(is:ie,js:je,k)*wc(is:ie,js:je,k)*wc(is:ie,js:je,k)
         turb_vars(is:ie,js:je,61) = turb_vars(is:ie,js:je,61) &
                                   + wc(is:ie,js:je,k)*uc(is:ie,js:je,k)*uc(is:ie,js:je,k)
         turb_vars(is:ie,js:je,62) = turb_vars(is:ie,js:je,62) &
                                   + wc(is:ie,js:je,k)*vc(is:ie,js:je,k)*vc(is:ie,js:je,k)
         turb_vars(is:ie,js:je,63) = turb_vars(is:ie,js:je,63) &
                                   + wc(is:ie,js:je,k)*wc(is:ie,js:je,k)*wc(is:ie,js:je,k)
         turb_vars(is:ie,js:je,64) = turb_vars(is:ie,js:je,64) &
                                   + area(is:ie,js:je,k)*n1(is:ie,js:je,k)
         turb_vars(is:ie,js:je,65) = turb_vars(is:ie,js:je,65) &
                                   + area(is:ie,js:je,k)*n2(is:ie,js:je,k)
         turb_vars(is:ie,js:je,66) = turb_vars(is:ie,js:je,66) &
                                   + area(is:ie,js:je,k)*n3(is:ie,js:je,k)
         turb_vars(is:ie,js:je,67) = turb_vars(is:ie,js:je,67) &
                                   + area(is:ie,js:je,k)**2.d0*n1(is:ie,js:je,k)*n1(is:ie,js:je,k)
         turb_vars(is:ie,js:je,68) = turb_vars(is:ie,js:je,68) &
                                   + area(is:ie,js:je,k)**2.d0*n2(is:ie,js:je,k)*n2(is:ie,js:je,k)
         turb_vars(is:ie,js:je,69) = turb_vars(is:ie,js:je,69) &
                                   + area(is:ie,js:je,k)**2.d0*n3(is:ie,js:je,k)*n3(is:ie,js:je,k)
         turb_vars(is:ie,js:je,70) = turb_vars(is:ie,js:je,70) &
                                   + area(is:ie,js:je,k)**2.d0*n1(is:ie,js:je,k)*n2(is:ie,js:je,k)
         turb_vars(is:ie,js:je,71) = turb_vars(is:ie,js:je,71) &
                                   + area(is:ie,js:je,k)**2.d0*n1(is:ie,js:je,k)*n3(is:ie,js:je,k)
         turb_vars(is:ie,js:je,72) = turb_vars(is:ie,js:je,72) &
                                   + area(is:ie,js:je,k)**2.d0*n2(is:ie,js:je,k)*n3(is:ie,js:je,k)
         end if ! TurbStatsOrder

         ! 4th order stats
         if ( TurbStatsOrder >= 4 ) then

         turb_vars(is:ie,js:je,73) = turb_vars(is:ie,js:je,73) &
                                   + uc(is:ie,js:je,k)*uc(is:ie,js:je,k)*uc(is:ie,js:je,k)*cvof(is:ie,js:je,k)
         turb_vars(is:ie,js:je,74) = turb_vars(is:ie,js:je,74) &
                                   + uc(is:ie,js:je,k)*vc(is:ie,js:je,k)*vc(is:ie,js:je,k)*cvof(is:ie,js:je,k)
         turb_vars(is:ie,js:je,75) = turb_vars(is:ie,js:je,75) &
                                   + uc(is:ie,js:je,k)*wc(is:ie,js:je,k)*wc(is:ie,js:je,k)*cvof(is:ie,js:je,k)
         turb_vars(is:ie,js:je,76) = turb_vars(is:ie,js:je,76) &
                                   + vc(is:ie,js:je,k)*uc(is:ie,js:je,k)*uc(is:ie,js:je,k)*cvof(is:ie,js:je,k)
         turb_vars(is:ie,js:je,77) = turb_vars(is:ie,js:je,77) &
                                   + vc(is:ie,js:je,k)*vc(is:ie,js:je,k)*vc(is:ie,js:je,k)*cvof(is:ie,js:je,k)
         turb_vars(is:ie,js:je,78) = turb_vars(is:ie,js:je,78) &
                                   + vc(is:ie,js:je,k)*wc(is:ie,js:je,k)*wc(is:ie,js:je,k)*cvof(is:ie,js:je,k)
         turb_vars(is:ie,js:je,79) = turb_vars(is:ie,js:je,79) &
                                   + wc(is:ie,js:je,k)*uc(is:ie,js:je,k)*uc(is:ie,js:je,k)*cvof(is:ie,js:je,k)
         turb_vars(is:ie,js:je,80) = turb_vars(is:ie,js:je,80) &
                                   + wc(is:ie,js:je,k)*vc(is:ie,js:je,k)*vc(is:ie,js:je,k)*cvof(is:ie,js:je,k)
         turb_vars(is:ie,js:je,81) = turb_vars(is:ie,js:je,81) &
                                   + wc(is:ie,js:je,k)*wc(is:ie,js:je,k)*wc(is:ie,js:je,k)*cvof(is:ie,js:je,k)

         do i=is,ie; do j=js,je
            dudx = (u(i,j,k) - u(i-1,j,k))/dx(i)
            dvdy = (v(i,j,k) - v(i,j-1,k))/dy(j)
            dwdz = (w(i,j,k) - w(i,j,k-1))/dz(k)

            dvdx = (v(i+1,j,k)+v(i+1,j-1,k)-v(i-1,j,k)-v(i-1,j-1,k))*0.25d0/dx(i)
            dwdx = (w(i+1,j,k)+w(i+1,j,k-1)-w(i-1,j,k)-w(i-1,j,k-1))*0.25d0/dx(i)
            dudy = (u(i,j+1,k)+u(i-1,j+1,k)-u(i,j-1,k)-u(i-1,j-1,k))*0.25d0/dy(j)
            dwdy = (w(i,j+1,k)+w(i,j+1,k-1)-w(i,j-1,k)-w(i,j-1,k-1))*0.25d0/dy(j)
            dudz = (u(i,j,k+1)+u(i-1,j,k+1)-u(i,j,k-1)-u(i-1,j,k-1))*0.25d0/dz(k)
            dvdz = (v(i,j,k+1)+v(i,j-1,k+1)-v(i,j,k-1)-v(i,j-1,k-1))*0.25d0/dz(k)

            turb_vars(i,j,82) = turb_vars(i,j,82) + dudx*dudx*cvof(i,j,k)
            turb_vars(i,j,83) = turb_vars(i,j,83) + dvdy*dvdy*cvof(i,j,k)
            turb_vars(i,j,84) = turb_vars(i,j,84) + dwdz*dwdz*cvof(i,j,k)
            turb_vars(i,j,85) = turb_vars(i,j,85) + dvdx*dvdx*cvof(i,j,k)
            turb_vars(i,j,86) = turb_vars(i,j,86) + dwdx*dwdx*cvof(i,j,k)
            turb_vars(i,j,87) = turb_vars(i,j,87) + dudy*dudy*cvof(i,j,k)
            turb_vars(i,j,88) = turb_vars(i,j,88) + dwdy*dwdy*cvof(i,j,k)
            turb_vars(i,j,89) = turb_vars(i,j,89) + dudz*dudz*cvof(i,j,k)
            turb_vars(i,j,90) = turb_vars(i,j,90) + dvdz*dvdz*cvof(i,j,k)
         end do; end do !i,j

         end if ! TurbStatsOrder
      end do ! 

      ! backup turbulence statistics
      if(mod(iStep,nbackup)==0) then 
         call backup_turb_write
      end if ! iStep

      ! Output turbulence statistics
     if ( mod(iStep,nStepOutputTurbStats) == 0 ) then 

      rootout = nPdomain-1
      if ( rank == rootout ) then 
         filename=TRIM(out_path)//'/turb-vars.dat'
         inquire(FILE=filename,EXIST=file_exist)
         if ( file_exist ) then 
            OPEN(UNIT=101,FILE=filename, STATUS='replace')
         else 
            OPEN(UNIT=101,FILE=filename, STATUS='new')
         end if ! file_exist
         write(101,*) "#",iSumTurbStats
      end if ! rank 

      if ( nPdomain > 1 ) then 
         if ( rank == rootout ) then 
            allocate(turb_vars_all(Mx,My,num_turb_vars,0:nPdomain-1))
            allocate(turb_vars_map(Nx,Ny,num_turb_vars))
            turb_vars_all = 0.d0 
            turb_vars_map = 0.d0 
         end if ! rank 
         ! collect all data to root
         counts = Mx*My*num_turb_vars   ! assuming same number of cells in each block
         if ( rank /= rootout ) then
            call MPI_ISEND(turb_vars(is,js,1),counts, & 
                           MPI_DOUBLE_PRECISION, rootout, 14, MPI_COMM_WORLD, req(1), ierr)
            call MPI_WAIT(req(1),sta(:,1),ierr)
         else
            do irank = 0,nPdomain-1
               if ( irank /= rootout ) then 
                  call MPI_IRECV(turb_vars_all(1,1,1,irank),counts, & 
                                 MPI_DOUBLE_PRECISION, irank, 14, MPI_COMM_WORLD, req(2), ierr)
                  call MPI_WAIT(req(2),sta(:,2),ierr)
               end if ! irank
            end do ! irank
         end if ! rank

         ! sum over z direction and map to global 
         if ( rank == rootout ) then
            turb_vars_all(:,:,:,rootout) = turb_vars(:,:,:)
            do irank = 0,nPdomain-1,nPz
               do irank_sum = 0,npz-1
                  cx = irank/(nPy*nPz)
                  cy = (irank-cx*nPy*nPz)/nPz 
                  i = cx*Mx + 1
                  j = cy*My + 1
                  turb_vars_map(i:i+Mx-1,j:j+My-1,:) & 
                     = turb_vars_map(i:i+Mx-1,j:j+My-1,:) & 
                     + turb_vars_all(1:  Mx  ,1:  My  ,:,irank+irank_sum)
               end do ! irank_sum
            end do ! irank
         end if ! nPz
            
         ! output
         if ( rank == rootout ) then 
            do i=1,Nx
               do j=1,Ny
                  write(101,'(2(E15.8,1X),83(E25.16,1X))') & 
                  (dble(i)-0.5d0)*dx(is),& 
                  (dble(j)-0.5d0)*dy(js),& 
                  turb_vars_map(i,j,1:num_turb_vars)
               end do !j 
               write(101,*)
            end do ! i
         end if ! rank
      else 
         do i=is,ie
            do j=js,je
               write(101,'(2(E15.8,1X),83(E25.16,1X))') & 
               x(i),y(j),turb_vars(i,j,1:num_turb_vars)
            end do !j 
            write(101,*)
         end do ! i
      end if ! nPdomain
      if ( rank == rootout ) then  
         CLOSE(UNIT=101)
         deallocate(turb_vars_all)
         deallocate(turb_vars_map)
      end if ! rank
     end if !iStep

     contains 
      subroutine backup_turb_write
         integer ::i,j
         character(len=100) :: filename
         if (OrganizeOutFolder) then
           filename = trim(out_path)//'/BACKUP_TURB/backup_turb_'//int2text(rank,padding)
         else 
           filename = trim(out_path)//'/backup_turb_'//int2text(rank,padding)
         end if ! OrganizeOutFolder 
         call system('touch '//trim(filename)//'; mv '//trim(filename)//' '//trim(filename)//'.old')
         OPEN(UNIT=101,FILE=trim(filename),status='REPLACE')
         write(101,*) iSumTurbStats
         do i=is,ie; do j=js,je
            write(101,'(83(E25.16,1X))') turb_vars(i,j,1:num_turb_vars)
         enddo; enddo
         CLOSE(UNIT=101)
      end subroutine backup_turb_write

      subroutine backup_turb_read
         integer ::i,j
         character(len=100) :: filename
         logical :: file_exist
         if (OrganizeOutFolder) then
           filename = trim(out_path)//'/BACKUP_TURB/backup_turb_'//int2text(rank,padding)
         else 
           filename = trim(out_path)//'/backup_turb_'//int2text(rank,padding)
         end if ! OrganizeOutFolder 
         inquire(FILE=filename,EXIST=file_exist)
         if ( file_exist ) then 
            OPEN(UNIT=101,FILE=trim(filename),status='unknown',action='read')
            read(101,*) iSumTurbStats
            do i=is,ie; do j=js,je
               read(101,*) turb_vars(i,j,1:num_turb_vars)
            enddo; enddo
            CLOSE(UNIT=101)
         else
            turb_vars = 0.d0 
            iSumTurbStats = 0 
         end if ! file_exit
      end subroutine backup_turb_read

   end subroutine calcTurbStats
!=================================================================================================
subroutine pressure_stats(t)
  use module_grid
  use module_flow
  use module_VOF
  use module_freesurface
  implicit none
  include "mpif.h"
  real(8) :: t
  integer :: i,j,k,ierr
  real(8) :: p_face(1:6), vol_face(1:6), p_face_global(1:6), vol_face_global(1:6), p_min, p_max
  real(8) :: p_liq, p_avg_liq, p_min_global, p_max_global, vol_liq, vol_liq_tot, vol
  real(8) :: p_dyn, p_dyn_avg, ke_box, ke_box_avg

  p_face = 0.0d0
  vol_face = 0.0d0
  !sum pressures on x+ face
  if(p_ind(1)>=is .and. p_ind(1)<=ie) then
     do j=js,je
        do k=ks,ke
           if (j>=p_ind(3) .and. j<=p_ind(4) .and.&
                k>=p_ind(5) .and. k<=p_ind(6) ) then
              i=p_ind(1)
              vol=dx(i)*dy(j)*dz(k)
              p_face(1) = p_face(1) + vol*(1.0d0-cvof(i,j,k))*(p(i,j,k)+&
                   0.5d0*rho(i,j,k)*( (0.5d0*(u(i,j,k)+u(i+1,j,k)))**2.d0 & 
                   +(0.5d0*(v(i,j,k)+v(i,j+1,k)))**2.d0+(0.5d0*(w(i,j,k)+w(i,j,k+1)))**2.d0))
              vol_face(1) = vol_face(1) + vol*(1.0d0-cvof(i,j,k))
           endif
        enddo
     enddo
  endif
  !sum pressures on x- face
  if(p_ind(2)>=is .and. p_ind(2)<=ie) then
     do j=js,je
        do k=ks,ke
           if (j>=p_ind(3) .and. j<=p_ind(4) .and.&
                k>=p_ind(5) .and. k<=p_ind(6) ) then
              i=p_ind(2)
              vol=dx(i)*dy(j)*dz(k)
              p_face(2) = p_face(2) + vol*(1.0d0-cvof(i,j,k))*(p(i,j,k)+&
                   0.5d0*rho(i,j,k)*( (0.5d0*(u(i,j,k)+u(i+1,j,k)))**2.d0 & 
                   +(0.5d0*(v(i,j,k)+v(i,j+1,k)))**2.d0+(0.5d0*(w(i,j,k)+w(i,j,k+1)))**2.d0))
              vol_face(2) = vol_face(2) + vol*(1.0d0-cvof(i,j,k))
           endif
        enddo
     enddo
  endif
  !sum pressures on y- face
  if(p_ind(3)>=js .and. p_ind(3)<=je) then
     do i=is,ie
        do k=ks,ke
           if (i>=p_ind(1) .and. i<=p_ind(2) .and.&
                k>=p_ind(5) .and. k<=p_ind(6) ) then
              j=p_ind(3)
              vol=dx(i)*dy(j)*dz(k)
              p_face(3) = p_face(3) + vol*(1.0d0-cvof(i,j,k))*(p(i,j,k)+&
                   0.5d0*rho(i,j,k)*( (0.5d0*(u(i,j,k)+u(i+1,j,k)))**2.d0 & 
                   +(0.5d0*(v(i,j,k)+v(i,j+1,k)))**2.d0+(0.5d0*(w(i,j,k)+w(i,j,k+1)))**2.d0))  
              vol_face(3) = vol_face(3) + vol*(1.0d0-cvof(i,j,k))
           endif
        enddo
     enddo
  endif
  !sum pressures on y+ face
  if(p_ind(4)>=js .and. p_ind(4)<=je) then
     do i=is,ie
        do k=ks,ke
           if (i>=p_ind(1) .and. i<=p_ind(2) .and.&
                k>=p_ind(5) .and. k<=p_ind(6) ) then
              j=p_ind(4)
              vol=dx(i)*dy(j)*dz(k)
              p_face(4) = p_face(4) + vol*(1.0d0-cvof(i,j,k))*(p(i,j,k)+&
                   0.5d0*rho(i,j,k)*( (0.5d0*(u(i,j,k)+u(i+1,j,k)))**2.d0 & 
                   +(0.5d0*(v(i,j,k)+v(i,j+1,k)))**2.d0+(0.5d0*(w(i,j,k)+w(i,j,k+1)))**2.d0))
              vol_face(4) = vol_face(4) + vol*(1.0d0-cvof(i,j,k))
           endif
        enddo
     enddo
  endif
  !sum pressures on z- face
  if(p_ind(5)>=ks .and. p_ind(5)<=ke) then
     do i=is,ie
        do j=js,je
           if (i>=p_ind(1) .and. i<=p_ind(2) .and.&
                j>=p_ind(3) .and. j<=p_ind(4) ) then
              k=p_ind(5)
              vol=dx(i)*dy(j)*dz(k)
              p_face(5) = p_face(5) + vol*(1.0d0-cvof(i,j,k))*(p(i,j,k)+&
                   0.5d0*rho(i,j,k)*( (0.5d0*(u(i,j,k)+u(i+1,j,k)))**2.d0 & 
                   +(0.5d0*(v(i,j,k)+v(i,j+1,k)))**2.d0+(0.5d0*(w(i,j,k)+w(i,j,k+1)))**2.d0))
              vol_face(5) = vol_face(5) + vol*(1.0d0-cvof(i,j,k))
           endif
        enddo
     enddo
  endif
  ! sum pressures on z+ face
  if(p_ind(6)>=ks .and. p_ind(6)<=ke) then
     do i=is,ie
        do j=js,je
           if (i>=p_ind(1) .and. i<=p_ind(2) .and.&
                j>=p_ind(3) .and. j<=p_ind(4) ) then
              k=p_ind(6)
              vol=dx(i)*dy(j)*dz(k)
              p_face(6) = p_face(6) + vol*(1.0d0-cvof(i,j,k))*(p(i,j,k)+&
                   0.5d0*rho(i,j,k)*( (0.5d0*(u(i,j,k)+u(i+1,j,k)))**2.d0 & 
                   +(0.5d0*(v(i,j,k)+v(i,j+1,k)))**2.d0+(0.5d0*(w(i,j,k)+w(i,j,k+1)))**2.d0))
              vol_face(6) = vol_face(6) + vol*(1.0d0-cvof(i,j,k))
           endif
        enddo
     enddo
  endif
  call MPI_ALLREDUCE(p_face, p_face_global, 6, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_Domain, ierr)
  call MPI_ALLREDUCE(vol_face, vol_face_global, 6, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_Domain, ierr)
  !Calculate face averages
  vol_face_global = vol_face_global + TINY_DOUBLE
  p_face_global(1)=p_face_global(1)/vol_face_global(1)
  p_face_global(2)=p_face_global(2)/vol_face_global(2)
  p_face_global(3)=p_face_global(3)/vol_face_global(3)
  p_face_global(4)=p_face_global(4)/vol_face_global(4)
  p_face_global(5)=p_face_global(5)/vol_face_global(5)
  p_face_global(6)=p_face_global(6)/vol_face_global(6)
  p_max = -1.0d50
  p_min = 1.0d50
  p_liq = 0.0d0
  vol_liq = 0.0d0
  p_dyn = 0.0d0
  ke_box = 0.0d0
  do k=ks,ke; do j=js,je; do i=is,ie
     if ( (i>=p_ind(1)).and.(i<=p_ind(2)) .and.&
        (j>=p_ind(3)).and.(j<=p_ind(4)) .and.&
        (k>=p_ind(5)).and.(k<=p_ind(6)) ) then
        vol=dx(i)*dy(j)*dz(k)
        ke_box = vol*(1.0d0-cvof(i,j,k))*(0.5d0*rho(i,j,k)*( (0.5d0*(u(i,j,k)+u(i+1,j,k)))**2.d0 & 
             +(0.5d0*(v(i,j,k)+v(i,j+1,k)))**2.d0+(0.5d0*(w(i,j,k)+w(i,j,k+1)))**2.d0) )
        p_liq=p_liq+vol*(1.0d0-cvof(i,j,k))*p(i,j,k) + ke_box 
        vol_liq = vol_liq + vol*(1.0d0-cvof(i,j,k))
        p_min = MIN(p_min,p(i,j,k))
        p_max = MAX(p_max,p(i,j,k))
        p_dyn =p_dyn+vol*(1.0d0-cvof(i,j,k))*p(i,j,k) 
     endif
  enddo;enddo;enddo
  call MPI_ALLREDUCE(p_liq, p_avg_liq, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_Domain, ierr)
  call MPI_ALLREDUCE(p_dyn, p_dyn_avg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_Domain, ierr)
  call MPI_ALLREDUCE(ke_box, ke_box_avg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_Domain, ierr)
  call MPI_ALLREDUCE(vol_liq, vol_liq_tot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_Domain, ierr)
  call MPI_ALLREDUCE(p_min, p_min_global, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_Domain, ierr)
  call MPI_ALLREDUCE(p_max, p_max_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_Domain, ierr)
  vol_liq_tot = vol_liq_tot + TINY_DOUBLE
  p_avg_liq = p_avg_liq/vol_liq_tot
  p_dyn_avg = p_dyn_avg/vol_liq_tot
  ke_box_avg = ke_box_avg/vol_liq_tot

  if (rank==0 .and. itimestep==1) then
     open(unit=20,file='pressure_stats.txt',position='append')
     write(20,*) "# 1:time 2:p_face(1-6) 8:p_min 9:p_max 10:pavg_liq 11:p_dyn_avg 12:KE_box_avg 13:Vol_liq"
     close(20)
  endif
  
  if (rank==0) then
     OPEN(UNIT=20,FILE='pressure_stats.txt',position='append')
     write(20,201)t,p_face_global(1:6),p_min_global,p_max_global,p_avg_liq,p_dyn_avg,ke_box_avg,vol_liq_tot
     CLOSE(20)
  endif
201  format(13es14.6e3)

end subroutine pressure_stats
!=================================================================================================
#ifdef PHASE_CHANGE
subroutine temp_stats
  use module_boil
  use module_grid
  implicit none
  integer :: i,j,k

  !if (rank==0) then
  OPEN(UNIT=20,FILE='temp_stats.txt',position='append')
  j=(js+je)/2
  k=(ks+ke)/2
  do i=is,ie
     write(20,201)x(i),Te(i,j,k)
  enddo
  CLOSE(20)
  !endif
201  format(2es14.6e2)

end subroutine temp_stats
#endif
!=================================================================================================
subroutine press_indices
  use module_grid
  use module_freesurface
  use module_2phase
  implicit none
  integer :: q

  p_ind(1) = floor(Nx*coord_min/xLength) + Ng
  p_ind(2) = floor(Nx*(coord_min+var_coord)/xLength) + Ng

  p_ind(3) =  floor(Ny*coord_min/yLength) + Ng
  p_ind(4) =  floor(Ny*(coord_min+var_coord)/yLength) + Ng

  p_ind(5) =  floor(Nz*coord_min/zLength) + Ng
  p_ind(6) =  floor(Nz*(coord_min+var_coord)/zLength) + Ng

  if (rank==0) then
     write(*,'("Grid indices at coord limits for pressure stats: ")')
     do q=1,6    
        write(*,'("Face ",I4,": ",I4)')q,p_ind(q)
     enddo
  endif
end subroutine press_indices
!=========================================================================================================================
subroutine momentumConvection()
  use module_flow
  use module_grid
  use module_tmpvar
  use module_IO

  if (AdvectionScheme=='QUICK') then
    call momentumConvectionQUICK(u,v,w,du,dv,dw)
  elseif (  AdvectionScheme=='ENO' .or. & 
            AdvectionScheme=='Superbee' .or. & 
            AdvectionScheme=='WENO' ) then
!   call momentumConvectionENO(u,v,w,du,dv,dw)
    call momentumConvection_onedim(u,v,w,u,du,1,AdvectionScheme)
    call momentumConvection_onedim(u,v,w,v,dv,2,AdvectionScheme)
    call momentumConvection_onedim(u,v,w,w,dw,3,AdvectionScheme)
  elseif (AdvectionScheme=='UpWind') then
    call momentumConvectionUpWind(u,v,w,du,dv,dw)
  elseif (AdvectionScheme=='Verstappen') then
    call momentumConvectionVerstappen(u,v,w,du,dv,dw)
  elseif (AdvectionScheme=='BCG') then
    call momentumConvectionBCG()
  else
     call pariserror("*** unknown convection scheme")
  endif

  if (out_mom .and. mod(itimestep-itimestepRestart,nout)==0) call write_mom_gnuplot(du,dv,dw,itimestep)

end subroutine momentumConvection

!=================================================================================================
! subroutine momentumConvectionENO
! calculates the convection terms in momentum equation using ENO scheme
! and returns them in du, dv, dw
!-------------------------------------------------------------------------------------------------
subroutine momentumConvection_onedim(u,v,w,phi,dphi,d,AdvScheme)
  use module_grid
  use module_tmpvar
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: u, v, w, phi
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: dphi
  character(20), intent(in) :: AdvScheme 
  integer :: i,j,k,d
  integer, dimension(3) :: is0,is1,ie0,ie1
  integer :: i0,j0,k0
!-------------------------------------ENO interpolation u-velocity--------------------------------
  call init_i0j0k0 (d,i0,j0,k0)

  is0(1)=is-1; is0(2)=js-1; is0(3)=ks-1
  ie0(1)=ie; ie0(2)=je; ie0(3)=ke

  is1 = is0
  ie1 = ie0
  if (d.eq.1) then
    is1(:) = is1(:) + 1
    ie1(1) = ieu+1
  endif

  do k=is1(3),ie1(3); do j=is1(2),ie1(2); do i=is1(1),ie1(1)
    if (u(i-i0,j,k)+u(i,j+j0,k+k0)>0.0) then
      work(i,j,k,1)=phi(i-i0,j,k)+0.5*slope_lim(phi(i-1-i0,j,k),phi(i-i0,j,k),phi(i+1-i0,j,k),AdvScheme)
    else
      work(i,j,k,1)=phi(i+1-i0,j,k)-0.5*slope_lim(phi(i-i0,j,k),phi(i+1-i0,j,k),phi(i+2-i0,j,k),AdvScheme)
    endif
  enddo; enddo; enddo

  is1 = is0
  ie1 = ie0
  if (d.eq.2) then
    is1(:) = is1(:) + 1
    ie1(2) = jev+1
  endif

  do k=is1(3),ie1(3); do j=is1(2),ie1(2); do i=is1(1),ie1(1)
    if(v(i,j-j0,k)+v(i+i0,j,k+k0)>0.0) then
      work(i,j,k,2)=phi(i,j-j0,k)+0.5*slope_lim(phi(i,j-1-j0,k),phi(i,j-j0,k),phi(i,j+1-j0,k),AdvScheme)
    else
      work(i,j,k,2)=phi(i,j+1-j0,k)-0.5*slope_lim(phi(i,j-j0,k),phi(i,j+1-j0,k),phi(i,j+2-j0,k),AdvScheme)
    endif
  enddo; enddo; enddo

  is1 = is0
  ie1 = ie0
  if (d.eq.3) then
    is1(:) = is1(:) + 1
    ie1(3) = kew + 1
  endif

  do k=is1(3),ie1(3); do j=is1(2),ie1(2); do i=is1(1),ie1(1)
    if(w(i,j,k-k0)+w(i+i0,j+j0,k)>0.0) then
      work(i,j,k,3)=phi(i,j,k-k0)+0.5*slope_lim(phi(i,j,k-1-k0),phi(i,j,k-k0),phi(i,j,k+1-k0),AdvScheme)
    else
      work(i,j,k,3)=phi(i,j,k+1-k0)-0.5*slope_lim(phi(i,j,k-k0),phi(i,j,k+1-k0),phi(i,j,k+2-k0),AdvScheme)
    endif
  enddo; enddo; enddo

  is1 = is0 +1
  ie1 = ie0
  if (d.eq.1) then
     ie1(1) = ieu
  elseif (d.eq.2) then
     ie1(2) = jev
  else
     ie1(3) = kew
  endif

  do k=is1(3),ie1(3); do j=is1(2),ie1(2); do i=is1(1),ie1(1)
     dphi(i,j,k)= -0.5*((u(i,j  ,k  )+u(i+i0,j+j0,k+k0))*work(i+i0,j  ,k  ,1)- &
          (u(i-1+i0,j  ,k  )+u(i-1,j+j0,k+k0 ))*work(i-1+i0 ,j  ,k  ,1))/dx(i) &
          -0.5*((v(i,j  ,k  )+v(i+i0,j+j0,k+k0))*work(i  ,j+j0 ,k  ,2)-&
          (v(i,j-1+j0,k  )+v(i+i0,j-1,k+k0 ))*work(i  ,j-1+j0,k  ,2))/dy(j)  &
          -0.5*((w(i,j  ,k  )+w(i+i0,j+j0,k+k0))*work(i  ,j  ,k+k0  ,3)-&
          (w(i,j  ,k-1+k0)+w(i+i0,j+j0,k-1))*work(i  ,j  ,k-1+k0,3))/dz(k)
  enddo; enddo; enddo
 
end subroutine momentumConvection_onedim

subroutine momentumConvectionENO(u,v,w,du,dv,dw)
  use module_grid
  use module_tmpvar
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: u, v, w
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: du, dv, dw
  real(8), external :: minabs
  integer :: i,j,k
!-------------------------------------ENO interpolation u-velocity--------------------------------
  do k=ks,ke; do j=js,je; do i=is,ieu+1
    if (u(i-1,j,k)+u(i,j,k)>0.0) then
      work(i,j,k,1)=u(i-1,j,k)+0.5*minabs((u(i,j,k)-u(i-1,j,k)),(u(i-1,j,k)-u(i-2,j,k)))
    else
      work(i,j,k,1)=u(i,j,k)-0.5*minabs((u(i+1,j,k)-u(i,j,k)),(u(i,j,k)-u(i-1,j,k)))
    endif
  enddo; enddo; enddo
  do k=ks,ke; do j=js-1,je; do i=is-1,ie
    if(v(i,j,k)+v(i+1,j,k)>0.0) then
      work(i,j,k,2)=u(i,j,k)+0.5*minabs((u(i,j+1,k)-u(i,j,k)),(u(i,j,k)-u(i,j-1,k)))
    else
      work(i,j,k,2)=u(i,j+1,k)-0.5*minabs((u(i,j+2,k)-u(i,j+1,k)),(u(i,j+1,k)-u(i,j,k)))
    endif
  enddo; enddo; enddo
  do k=ks-1,ke; do j=js,je; do i=is-1,ie
    if(w(i,j,k)+w(i+1,j,k)>0.0) then
      work(i,j,k,3)=u(i,j,k)+0.5*minabs((u(i,j,k+1)-u(i,j,k)),(u(i,j,k)-u(i,j,k-1)))
    else
      work(i,j,k,3)=u(i,j,k+1)-0.5*minabs((u(i,j,k+2)-u(i,j,k+1)),(u(i,j,k+1)-u(i,j,k)))
    endif
  enddo; enddo; enddo
  do k=ks,ke;  do j=js,je; do i=is,ieu
    du(i,j,k)= -0.5*((u(i,j  ,k  )+u(i+1,j  ,k  ))*work(i+1,j  ,k  ,1)- &
                    (u(i,j  ,k  )+u(i-1,j  ,k  ))*work(i  ,j  ,k  ,1))/dxh(i) &
              -0.5*((v(i,j  ,k  )+v(i+1,j  ,k  ))*work(i  ,j  ,k  ,2)-&
                    (v(i,j-1,k  )+v(i+1,j-1,k  ))*work(i  ,j-1,k  ,2))/dy(j)  &
              -0.5*((w(i,j  ,k  )+w(i+1,j  ,k  ))*work(i  ,j  ,k  ,3)-&
                    (w(i,j  ,k-1)+w(i+1,j  ,k-1))*work(i  ,j  ,k-1,3))/dz(k)
  enddo; enddo; enddo

!-------------------------------------ENO interpolation v-velocity--------------------------------
  do k=ks,ke; do j=js,jev+1; do i=is,ie
    if (v(i,j-1,k)+v(i,j,k)>0.0) then
      work(i,j,k,2)=v(i,j-1,k)+0.5*minabs((v(i,j,k)-v(i,j-1,k)),(v(i,j-1,k)-v(i,j-2,k)))
    else
      work(i,j,k,2)=v(i,j,k)-0.5*minabs((v(i,j+1,k)-v(i,j,k)),(v(i,j,k)-v(i,j-1,k)))
    endif
  enddo; enddo; enddo
  do k=ks,ke; do j=js-1,je; do i=is-1,ie
    if(u(i,j,k)+u(i,j+1,k)>0.0) then
      work(i,j,k,1)=v(i,j,k)+0.5*minabs((v(i+1,j,k)-v(i,j,k)),(v(i,j,k)-v(i-1,j,k)))
    else
      work(i,j,k,1)=v(i+1,j,k)-0.5*minabs((v(i+2,j,k)-v(i+1,j,k)),(v(i+1,j,k)-v(i,j,k)))
    endif
  enddo; enddo; enddo
  do k=ks-1,ke; do j=js-1,je; do i=is,ie
    if(w(i,j,k)+w(i,j+1,k)>0.0) then
      work(i,j,k,3)=v(i,j,k)+0.5*minabs((v(i,j,k+1)-v(i,j,k)),(v(i,j,k)-v(i,j,k-1)))
    else
      work(i,j,k,3)=v(i,j,k+1)-0.5*minabs((v(i,j,k+2)-v(i,j,k+1)),(v(i,j,k+1)-v(i,j,k)))
    endif
  enddo; enddo; enddo
  do k=ks,ke;  do j=js,jev; do i=is,ie
    dv(i,j,k)= &
              -0.5*((u(i  ,j,k  )+u(i  ,j+1,k  ))*work(i  ,j  ,k  ,1)-&
                    (u(i-1,j,k  )+u(i-1,j+1,k  ))*work(i-1,j  ,k  ,1))/dx(i)  &
              -0.5*((v(i  ,j,k  )+v(i  ,j+1,k  ))*work(i  ,j+1,k  ,2)-&
                    (v(i  ,j,k  )+v(i  ,j-1,k  ))*work(i  ,j  ,k  ,2))/dyh(j) &
              -0.5*((w(i  ,j,k  )+w(i  ,j+1,k  ))*work(i  ,j  ,k  ,3)-&
                    (w(i  ,j,k-1)+w(i  ,j+1,k-1))*work(i  ,j  ,k-1,3))/dz(k)
  enddo; enddo; enddo
!-------------------------------------ENO interpolation w-velocity--------------------------------
  do k=ks,kew+1; do j=js,je; do i=is,ie
    if (w(i,j,k-1)+w(i,j,k)>0.0) then
      work(i,j,k,3)=w(i,j,k-1)+0.5*minabs((w(i,j,k)-w(i,j,k-1)),(w(i,j,k-1)-w(i,j,k-2)))
    else
      work(i,j,k,3)=w(i,j,k)-0.5*minabs((w(i,j,k+1)-w(i,j,k)),(w(i,j,k)-w(i,j,k-1)))
    endif
  enddo; enddo; enddo
  do k=ks-1,ke; do j=js,je; do i=is-1,ie
    if(u(i,j,k)+u(i,j,k+1)>0.0) then
      work(i,j,k,1)=w(i,j,k)+0.5*minabs((w(i+1,j,k)-w(i,j,k)),(w(i,j,k)-w(i-1,j,k)))
    else
      work(i,j,k,1)=w(i+1,j,k)-0.5*minabs((w(i+2,j,k)-w(i+1,j,k)),(w(i+1,j,k)-w(i,j,k)))
    endif
  enddo; enddo; enddo
  do k=ks-1,ke; do j=js-1,je; do i=is,ie
    if(v(i,j,k)+v(i,j,k+1)>0.0) then
      work(i,j,k,2)=w(i,j,k)+0.5*minabs((w(i,j+1,k)-w(i,j,k)),(w(i,j,k)-w(i,j-1,k)))
    else
      work(i,j,k,2)=w(i,j+1,k)-0.5*minabs((w(i,j+2,k)-w(i,j+1,k)),(w(i,j+1,k)-w(i,j,k)))
    endif
  enddo; enddo; enddo
  do k=ks,kew;  do j=js,je; do i=is,ie
    dw(i,j,k)= & 
              -0.5*((u(i  ,j  ,k)+u(i  ,j  ,k+1))*work(i  ,j  ,k  ,1)-&
                    (u(i-1,j  ,k)+u(i-1,j  ,k+1))*work(i-1,j  ,k  ,1))/dx(i)  &
              -0.5*((v(i  ,j  ,k)+v(i  ,j  ,k+1))*work(i  ,j  ,k  ,2)-&
                    (v(i  ,j-1,k)+v(i  ,j-1,k+1))*work(i  ,j-1,k  ,2))/dy(j)  &
              -0.5*((w(i  ,j  ,k)+w(i  ,j  ,k+1))*work(i  ,j  ,k+1,3)-&
                    (w(i  ,j  ,k)+w(i  ,j  ,k-1))*work(i  ,j  ,k  ,3))/dzh(k) 
  enddo; enddo; enddo

end subroutine momentumConvectionENO
!=================================================================================================
! subroutine momentumConvectionUpWind
! calculates the convection terms in mumentum equation using UpWind scheme
! and returns them in du, dv, dw
!-------------------------------------------------------------------------------------------------
subroutine momentumConvectionUpWind(u,v,w,du,dv,dw)
  use module_grid
  use module_tmpvar
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: u, v, w
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: du, dv, dw
  real(8), external :: minabs
  integer :: i,j,k
!-------------------------------------UpWind interpolation u-velocity--------------------------------
  do k=ks,ke; do j=js,je; do i=is,ieu+1
    if (u(i-1,j,k)+u(i,j,k)>0.0) then
      work(i,j,k,1)=u(i-1,j,k)
    else
      work(i,j,k,1)=u(i,j,k)
    endif
  enddo; enddo; enddo
  do k=ks,ke; do j=js-1,je; do i=is-1,ie
    if(v(i,j,k)+v(i+1,j,k)>0.0) then
      work(i,j,k,2)=u(i,j,k)
    else
      work(i,j,k,2)=u(i,j+1,k)
    endif
  enddo; enddo; enddo
  do k=ks-1,ke; do j=js,je; do i=is-1,ie
    if(w(i,j,k)+w(i+1,j,k)>0.0) then
      work(i,j,k,3)=u(i,j,k)
    else
      work(i,j,k,3)=u(i,j,k+1)
    endif
  enddo; enddo; enddo
  do k=ks,ke;  do j=js,je; do i=is,ieu
    du(i,j,k)= -0.5*((u(i,j ,k  )+u(i+1,j  ,k  ))*work(i+1,j  ,k  ,1)- &
                    (u(i,j  ,k  )+u(i-1,j  ,k  ))*work(i  ,j  ,k  ,1))/dxh(i) &
              -0.5*((v(i,j  ,k  )+v(i+1,j  ,k  ))*work(i  ,j  ,k  ,2)-&
                    (v(i,j-1,k  )+v(i+1,j-1,k  ))*work(i  ,j-1,k  ,2))/dy(j)  &
              -0.5*((w(i,j  ,k  )+w(i+1,j  ,k  ))*work(i  ,j  ,k  ,3)-&
                    (w(i,j  ,k-1)+w(i+1,j  ,k-1))*work(i  ,j  ,k-1,3))/dz(k)
  enddo; enddo; enddo

!-------------------------------------UpWind interpolation v-velocity--------------------------------
  do k=ks,ke; do j=js,jev+1; do i=is,ie
    if (v(i,j-1,k)+v(i,j,k)>0.0) then
      work(i,j,k,2)=v(i,j-1,k)
    else
      work(i,j,k,2)=v(i,j,k)
    endif
  enddo; enddo; enddo
  do k=ks,ke; do j=js-1,je; do i=is-1,ie
    if(u(i,j,k)+u(i,j+1,k)>0.0) then
      work(i,j,k,1)=v(i,j,k)
    else
      work(i,j,k,1)=v(i+1,j,k)
    endif
  enddo; enddo; enddo
  do k=ks-1,ke; do j=js-1,je; do i=is,ie
    if(w(i,j,k)+w(i,j+1,k)>0.0) then
      work(i,j,k,3)=v(i,j,k)
    else
      work(i,j,k,3)=v(i,j,k+1)
    endif
  enddo; enddo; enddo
  do k=ks,ke;  do j=js,jev; do i=is,ie
    dv(i,j,k)= &
              -0.5*((u(i  ,j,k  )+u(i  ,j+1,k  ))*work(i  ,j  ,k  ,1)-&
                    (u(i-1,j,k  )+u(i-1,j+1,k  ))*work(i-1,j  ,k  ,1))/dx(i)  &
              -0.5*((v(i  ,j,k  )+v(i  ,j+1,k  ))*work(i  ,j+1,k  ,2)-&
                    (v(i  ,j,k  )+v(i  ,j-1,k  ))*work(i  ,j  ,k  ,2))/dyh(j) &
              -0.5*((w(i  ,j,k  )+w(i  ,j+1,k  ))*work(i  ,j  ,k  ,3)-&
                    (w(i  ,j,k-1)+w(i  ,j+1,k-1))*work(i  ,j  ,k-1,3))/dz(k)
  enddo; enddo; enddo
!-------------------------------------UpWind interpolation w-velocity--------------------------------
  do k=ks,kew+1; do j=js,je; do i=is,ie
    if (w(i,j,k-1)+w(i,j,k)>0.0) then
      work(i,j,k,3)=w(i,j,k-1)
    else
      work(i,j,k,3)=w(i,j,k)
    endif
  enddo; enddo; enddo
  do k=ks-1,ke; do j=js,je; do i=is-1,ie
    if(u(i,j,k)+u(i,j,k+1)>0.0) then
      work(i,j,k,1)=w(i,j,k)
    else
      work(i,j,k,1)=w(i+1,j,k)
    endif
  enddo; enddo; enddo
  do k=ks-1,ke; do j=js-1,je; do i=is,ie
    if(v(i,j,k)+v(i,j,k+1)>0.0) then
      work(i,j,k,2)=w(i,j,k)
    else
      work(i,j,k,2)=w(i,j+1,k)
    endif
  enddo; enddo; enddo
  do k=ks,kew;  do j=js,je; do i=is,ie
    dw(i,j,k)= & 
              -0.5*((u(i  ,j  ,k)+u(i  ,j  ,k+1))*work(i  ,j  ,k  ,1)-&
                    (u(i-1,j  ,k)+u(i-1,j  ,k+1))*work(i-1,j  ,k  ,1))/dx(i)  &
              -0.5*((v(i  ,j  ,k)+v(i  ,j  ,k+1))*work(i  ,j  ,k  ,2)-&
                    (v(i  ,j-1,k)+v(i  ,j-1,k+1))*work(i  ,j-1,k  ,2))/dy(j)  &
              -0.5*((w(i  ,j  ,k)+w(i  ,j  ,k+1))*work(i  ,j  ,k+1,3)-&
                    (w(i  ,j  ,k)+w(i  ,j  ,k-1))*work(i  ,j  ,k  ,3))/dzh(k) 
  enddo; enddo; enddo
end subroutine momentumConvectionUpWind
!=================================================================================================
! subroutine momentumConvection
! calculates the convection terms in the momentum equations using a QUICK scheme
! and returns them in du, dv, dw
!-------------------------------------------------------------------------------------------------
subroutine momentumConvectionQUICK(u,v,w,du,dv,dw)
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
      du(i,j,k)= &
      - ( work(i+1,j  ,k  ,1)**2 - work(i  ,j  ,k  ,1)**2 )/dxh(i) &
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
     dv(i,j,k)= &
      -0.5*((u(i  ,j,k  )+u(i  ,j+1,k  ))*work(i  ,j  ,k  ,1)-        &
      (u(i-1,j,k  )+u(i-1,j+1,k  ))*work(i-1,j  ,k  ,1))/dx(i)  &
      -   ( work(i  ,j+1,k  ,2)**2 - work(i  ,j  ,k  ,2)**2 )/dyh(j) &
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
      dw(i,j,k)= &
      -0.5*((u(i  ,j  ,k)+u(i  ,j  ,k+1))*work(i  ,j  ,k  ,1)-        &
      (u(i-1,j  ,k)+u(i-1,j  ,k+1))*work(i-1,j  ,k  ,1))/dx(i)  &
      -0.5*((v(i  ,j  ,k)+v(i  ,j  ,k+1))*work(i  ,j  ,k  ,2)-        &
      (v(i  ,j-1,k)+v(i  ,j-1,k+1))*work(i  ,j-1,k  ,2))/dy(j)  &
      -( work(i  ,j  ,k+1,3)**2 - work(i  ,j  ,k  ,3)**2 )/dzh(k)
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
end subroutine momentumConvectionQUICK

!=================================================================================================
! subroutine momentumConvectionVerstappen
! calculates the convection terms in the momentum equations using the Verstappen
! symmetric scheme
!-------------------------------------------------------------------------------------------------
subroutine momentumConvectionVerstappen (u,v,w,du,dv,dw)
  use module_grid
  use module_tmpvar
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: u, v, w
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: du, dv, dw
  real(8), external :: minabs
  integer :: i,j,k

  !fixme: This only works for regular uniform grids
  do k=ks,ke;  do j=js,je; do i=is,ieu
      work(i,j,k,1) = (u(i,j,k) + u(i+1,j,k))*u(i+1,j,k) &
                     -(u(i,j,k) + u(i-1,j,k))*u(i-1,j,k) &
                     +(v(i,j,k) + v(i+1,j,k))*u(i,j+1,k) &
                     -(v(i,j-1,k) + v(i+1,j-1,k))*u(i,j-1,k) &
                     +(w(i,j,k) + w(i+1,j,k))*u(i,j,k+1) &
                     -(w(i,j,k-1) + w(i+1,j,k-1))*u(i,j,k-1)
      work(i,j,k,1) = work(i,j,k,1)*0.25d0
  enddo; enddo; enddo

  do k=ks,ke;  do j=js,je; do i=is,ieu
      du(i,j,k)= - work(i,j,k,1)/dxh(i)
  enddo; enddo; enddo

  do k=ks,ke;  do j=js,jev; do i=is,ie
       work(i,j,k,1) = (v(i,j,k) + v(i,j+1,k))*v(i,j+1,k) &
                     -(v(i,j,k) + v(i,j-1,k))*v(i,j-1,k) &
                     +(u(i,j,k) + u(i,j+1,k))*v(i+1,j,k) &
                     -(u(i-1,j,k) + u(i-1,j+1,k))*v(i-1,j,k) &
                     +(w(i,j,k) + w(i,j+1,k))*v(i,j,k+1) &
                     -(w(i,j,k-1) + w(i,j+1,k-1))*v(i,j,k-1)
      work(i,j,k,1) = work(i,j,k,1)*0.25d0
  enddo; enddo; enddo

  do k=ks,ke;  do j=js,jev; do i=is,ie
     dv(i,j,k)= - work(i,j,k,1)/dxh(i)
  enddo; enddo; enddo

  do k=ks,kew;  do j=js,je; do i=is,ie
      work(i,j,k,1) = (w(i,j,k) + w(i,j,k+1))*w(i,j,k+1)     &
                     -(w(i,j,k) + w(i,j,k-1))*w(i,j,k-1)     &
                     +(u(i,j,k) + u(i,j,k+1))*w(i+1,j,k)     &
                     -(u(i-1,j,k) + u(i-1,j,k+1))*w(i-1,j,k) &
                     +(v(i,j,k) + v(i,j,k+1))*w(i,j+1,k)     &
                     -(v(i,j-1,k) + v(i,j-1,k+1))*w(i,j-1,k)
      work(i,j,k,1) = work(i,j,k,1)*0.25d0
  enddo; enddo; enddo

  do k=ks,kew;  do j=js,je; do i=is,ie
      dw(i,j,k)= - work(i,j,k,1)/dxh(i)
  enddo; enddo; enddo

end subroutine momentumConvectionVerstappen
!=======================================================================================
! subroutine momentumConvectionBCG
! calculates the convection terms in the momentum equations using a QUICK scheme
! and returns them in du, dv, dw
!-------------------------------------------------------------------------------------------------
subroutine momentumConvectionBCG()
  use module_grid
  use module_tmpvar
  use module_BC
  use module_poisson
  use module_flow
  use module_VOF
  use module_2phase
  use module_surface_tension

  implicit none
  include 'mpif.h'
  real(8), external :: minabs

  integer :: i,j,k,ierr
  integer :: req(48),sta(MPI_STATUS_SIZE,48)
 
  !prediction stage
  call predict_velocity(u,du,1,dt,u,v,w,work(:,:,:,1))
  call predict_velocity(v,dv,2,dt,v,u,w,work(:,:,:,2))
  call predict_velocity(w,dw,3,dt,w,u,v,work(:,:,:,3))

  call SetVelocityBC(work(:,:,:,1),work(:,:,:,2),work(:,:,:,3),umask,vmask,wmask,time,dt,0)
  call do_ghost_vector(work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))
  !-----------------------------------------PROJECTION STEP-----------------------------------------
  tmp = p
  call SetPressureBC(umask,vmask,wmask)
  call SetupPoisson(work(:,:,:,1),work(:,:,:,2),work(:,:,:,3), &
                    umask,vmask,wmask,rho,dt,A,tmp,VolumeSource)
  ! (div u)*dt < epsilon => div u < epsilon/dt => maxresidual : maxerror/dt 
  if(HYPRE)then
    call poi_solve(A,tmp,maxError/MaxDt,maxit,it,HYPRESolverType)
    call ghost_x(tmp,1,req(1:4 ))
    call ghost_y(tmp,1,req(5:8 ))
    call ghost_z(tmp,1,req(9:12)) 
    call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
  else
    call NewSolver(A,tmp,maxError/MaxDt,beta,maxit,it,ierr,ResNormOrderPressure)
  endif
  if (.not.FreeSurface) then
    do k=ks,ke;  do j=js,je; do i=is,ieu    ! CORRECT THE u-velocity 
      work(i,j,k,1)=work(i,j,k,1)-dt*(2.0*umask(i,j,k)/dxh(i))*(tmp(i+1,j,k)-tmp(i,j,k))/ &
                                     (rho(i+1,j,k)+rho(i,j,k))
    enddo; enddo; enddo

    do k=ks,ke;  do j=js,jev; do i=is,ie    ! CORRECT THE v-velocity
      work(i,j,k,2)=work(i,j,k,2)-dt*(2.0*vmask(i,j,k)/dyh(j))*(tmp(i,j+1,k)-tmp(i,j,k))/ &
                                      (rho(i,j+1,k)+rho(i,j,k))
    enddo; enddo; enddo

    do k=ks,kew;  do j=js,je; do i=is,ie   ! CORRECT THE w-velocity
      work(i,j,k,3)=work(i,j,k,3)-dt*(2.0*wmask(i,j,k)/dzh(k))*(tmp(i,j,k+1)-tmp(i,j,k))/ &
                                     (rho(i,j,k+1)+rho(i,j,k))
    enddo; enddo; enddo
  else
    do k=ks,ke;  do j=js,je; do i=is,ieu    ! CORRECT THE u-velocity 
      work(i,j,k,1)=work(i,j,k,1)-dt*(2.0*umask(i,j,k)/(dxh(i)-x_mod(i,j,k)))*(tmp(i+1,j,k)-tmp(i,j,k))/ &
                                      (rho(i+1,j,k)+rho(i,j,k))
    enddo; enddo; enddo

    do k=ks,ke;  do j=js,jev; do i=is,ie    ! CORRECT THE v-velocity
      work(i,j,k,2)=work(i,j,k,2)-dt*(2.0*vmask(i,j,k)/(dyh(j)-y_mod(i,j,k)))*(tmp(i,j+1,k)-tmp(i,j,k))/ &
                                      (rho(i,j+1,k)+rho(i,j,k))
    enddo; enddo; enddo

    do k=ks,kew;  do j=js,je; do i=is,ie   ! CORRECT THE w-velocity
      work(i,j,k,3)=work(i,j,k,3)-dt*(2.0*wmask(i,j,k)/(dzh(k)-z_mod(i,j,k)))*(tmp(i,j,k+1)-tmp(i,j,k))/ &
                                      (rho(i,j,k+1)+rho(i,j,k))
    enddo; enddo; enddo
  endif 

  call do_ghost_vector(work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))

  !advection step
  call predict_velocity(u,du,1,dt,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3),tmp)
  du = tmp
  call predict_velocity(v,dv,2,dt,work(:,:,:,2),work(:,:,:,1),work(:,:,:,3),tmp)
  dv = tmp
  call predict_velocity(w,dw,3,dt,work(:,:,:,3),work(:,:,:,1),work(:,:,:,2),tmp)
  dw = tmp

  !check: required?
  call do_ghost_vector(du,dv,dw)

  work(:,:,:,1) = du
  work(:,:,:,2) = dv
  work(:,:,:,3) = dw

  !cell centered fluxes
  do k=ks-1,ke+1
    do j=js-1,je+1
      do i=is-1,ie+1
        du(i,j,k) = -0.5d0*(work(i,j,k,1)*work(i,j,k,1) - &
                    work(i-1,j,k,1)*work(i-1,j,k,1))/dxh(i) &
                  - 0.5d0*(work(i,j,k,2)+work(i,j-1,k,2))*  &
                    (work(i,j,k,1)-work(i,j-1,k,1))/dyh(j)  &
                  - 0.5d0*(work(i,j,k,3)+work(i,j,k-1,3))*  &
                    (work(i,j,k,1)-work(i,j,k-1,1))/dzh(j)  
        dv(i,j,k) = -0.5d0*(work(i,j,k,2)*work(i,j,k,2) - &
                    work(i-1,j,k,2)*work(i-1,j,k,2))/dyh(i) &
                  - 0.5d0*(work(i,j,k,1)+work(i,j-1,k,1))*  &
                    (work(i,j,k,2)-work(i,j-1,k,2))/dxh(j)  &
                  - 0.5d0*(work(i,j,k,3)+work(i,j,k-1,3))*  &
                    (work(i,j,k,2)-work(i,j,k-1,2))/dzh(j)  
        dw(i,j,k) = -0.5d0*(work(i,j,k,3)*work(i,j,k,3) - &
                    work(i-1,j,k,3)*work(i-1,j,k,3))/dzh(i) &
                  - 0.5d0*(work(i,j,k,2)+work(i,j-1,k,2))*  &
                    (work(i,j,k,3)-work(i,j-1,k,3))/dyh(j)  &
                  - 0.5d0*(work(i,j,k,1)+work(i,j,k-1,1))*  &
                    (work(i,j,k,3)-work(i,j,k-1,3))/dxh(j)  
      enddo
    enddo
  enddo

  call do_ghost_vector(du,dv,dw)

  !reinterpolating at face centers
  work(:,:,:,1) = du
  work(:,:,:,2) = dv
  work(:,:,:,3) = dw

  do k=ks-1,ke+1
    do j=js-1,je+1
      do i=is-1,ie+1
        du(i,j,k) = 0.5d0*(work(i,j,k,1) + work(i+1,j,k,1)) 
        dv(i,j,k) = 0.5d0*(work(i,j,k,2) + work(i+1,j,k,2)) 
        dw(i,j,k) = 0.5d0*(work(i,j,k,3) + work(i+1,j,k,3)) 
      enddo
    enddo
  enddo

  call do_ghost_vector(du,dv,dw)

end subroutine momentumConvectionBCG

subroutine predict_velocity(u,du,d,dt,upred,vpred,wpred,unew)

  use module_grid
  implicit none
  
  integer i,j,k,d,iaux,jaux,kaux
  integer i0,j0,k0
  integer onex, oney, onez, usign

  REAL (8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), INTENT(IN) :: u, du,upred, vpred, wpred
  REAL (8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), INTENT(INOUT) :: unew
  real(8) :: dt, uc, unorm, grad, dxcell,dxcell2
  real(8) :: fyy, fzz
  
  call init_i0j0k0 (d,i0,j0,k0)
  iaux = 0; jaux = 0; kaux = 0
  onex = 0; oney = 0; onez = 0

  do k=ks-1,ke+1
    do j=js-1,je+1
      do i=is-1,ie+1

        uc    = (u(i,j,k) + u(i-i0,j-j0,k-k0))*0.5d0
        usign = nint(sign(1.d0,uc))

        if (d.eq.1) then
          iaux  = -(usign + 1)/2
          onex  = 1
          dxcell = dxh(i)
          dxcell2 = dxh(i+iaux)
        elseif (d.eq.2) then
          jaux  = -(usign + 1)/2
          oney  = 1
          dxcell = dyh(j)
          dxcell2 = dyh(j+jaux)
        else 
          kaux  = -(usign + 1)/2
          onez  = 1
          dxcell = dzh(k)
          dxcell2 = dzh(k+kaux)
        endif
        
        unorm = uc*dt/dxcell

        grad  = (upred(i+iaux,j+jaux,k+kaux) &
                 -upred(i+iaux-onex,j+jaux-oney,k+kaux-onez))/dxcell2

        if (vpred(i+iaux,j+jaux,k+kaux).lt.0.d0) then 
          fyy = upred(i+iaux+onex,j+jaux+oney,k+kaux+onez) &
              - upred(i+iaux,j+jaux,k+kaux)
        else
          fyy = upred(i+iaux,j+jaux,k+kaux) &
              - upred(i+iaux-onex,j+jaux-oney,k+kaux-onez)
        endif

        if (wpred(i+iaux,j+jaux,k+kaux).lt.0.d0) then 
          fzz = upred(i+iaux+onex,j+jaux+oney,k+kaux+onez) &
              - upred(i+iaux,j+jaux,k+kaux)
        else
          fzz = upred(i+iaux,j+jaux,k+kaux) &
              - upred(i+iaux-onex,j+jaux-oney,k+kaux-onez)
        endif
!fixme: check the dxcell values in the traverse terms
        unew(i,j,k) = (u(i+iaux,j+jaux,k+kaux) &
                     + u(i+iaux-onex,j+jaux-oney,k+kaux-onez))*0.5d0 &
        + du(i+iaux,j+jaux,k+kaux)*dt/2.d0 &
        + usign*min(1.d0,1.d0-usign*unorm)*grad*dxcell/2.d0 &
        - dt*vpred(i+iaux,j+jaux,k+kaux)*fyy/(2.d0*dxcell) &
        - dt*wpred(i+iaux,j+kaux,k+kaux)*fzz/(2.d0*dxcell) 

      enddo
    enddo
  enddo

end subroutine predict_velocity

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
    du(i,j,k)= du(i,j,k)+(   &
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
    dv(i,j,k)= dv(i,j,k)+(   &
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
    dw(i,j,k)= dw(i,j,k)+(   &
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
    du(i,j,k)= du(i,j,k)+(   &
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
    dv(i,j,k)= dv(i,j,k)+(   &
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
    dw(i,j,k)= dw(i,j,k)+(   &
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
#ifdef PHASE_CHANGE
  use module_boil
  use module_VOF
#endif

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
    call pariserror("volumeForce: invalid buoyancy option")
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
! Calculates the surface tension force in the momentum equations and adds them to du,dv,dw
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
  include 'mpif.h'
  real(8) :: kappa,deltax
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: du, dv, dw, rho
  integer :: i,j,k,n,ierr
  deltax=dx(nx)
  call get_all_curvatures(tmp,0)
  call my_timer(7)
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
!          if(.not.((vof_flag(i,j,k+1)==1.and.vof_flag(i,j,k)==0) &
!                .or.(vof_flag(i,j,k+1)==0.and.vof_flag(i,j,k)==1))) then
             geom_case_count(18) = geom_case_count(18) + 1
!           endif
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
 !          if(.not.((vof_flag(i,j+1,k)==1.and.vof_flag(i,j,k)==0) &
 !               .or.(vof_flag(i,j+1,k)==0.and.vof_flag(i,j,k)==1))) then
              geom_case_count(19) = geom_case_count(19) + 1
 !          endif
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
 !          if(.not.((vof_flag(i,j,k+1)==1.and.vof_flag(i,j,k)==0) &
 !               .or.(vof_flag(i,j,k+1)==0.and.vof_flag(i,j,k)==1))) then
              geom_case_count(20) = geom_case_count(20) + 1
 !          endif
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
  
  use module_timer
  implicit none
  include 'mpif.h'
  integer :: ierr, i,j,k
  integer :: nd(3)
  Nxt=Nx+2*Ng; Nyt=Ny+2*Ng; Nzt=Nz+2*Ng ! total number of cells

  if(rank==0) call check_integers()
 
  if(rank<nPdomain)then
    allocate(dims(ndim),periodic(ndim),reorder(ndim),coords(ndim),STAT=ierr)
    dims(1) = nPx; dims(2) = nPy; dims(3) = nPz
    
    reorder = 1

! Initialise boundary condition and MPI_Comm_cart periodicity
    periodic = 0 
    do i=1,ndim 
       if (bdry_cond(i) == 1) then
          periodic(i) = 1
          if (bdry_cond(i+3) /= 1) call pariserror("inconsistent boundary conditions")
       endif
    enddo

! determine if check of bc setup is possible
    do i=1,2*ndim
       if (bdry_cond(i) >= 5) check_setup=.false.
    enddo

    call MPI_Cart_Create(MPI_Comm_Domain,ndim,dims,periodic,reorder,MPI_Comm_Cart,ierr)
    if (ierr /= 0) call pariserror("*** Grid: unsuccessful Cartesian MPI-initialization")

    call MPI_Cart_Coords(MPI_Comm_Cart,rank,ndim,coords,ierr)
    if (ierr /= 0) call pariserror("*** Grid: unsuccessful in getting topology coordinates")

!     print *, "rank",rank,"coords",coords
!     stop

    nd(1)=Mx; nd(2)=My; nd(3)=Mz
    call update_bounds(nd)

!     ! For pressure correction with inflow
!     ! dp/dn = 0 for inflow bc on face 1 == x- 
!     ! inflow bc on other faces not implemented yet.  
!     if(bdry_cond(1)==3 .and. coords(1)==0) then
!        istpc=is+1
!     else
!        istpc=is
!     endif

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
           wmask(imin:imax,jmin:jmax,kmin:kmax))  ! 13.5*2=27

    allocate(tmp(imin:imax,jmin:jmax,kmin:kmax), work(imin:imax,jmin:jmax,kmin:kmax,3), &
               A(is:ie,js:je,ks:ke,1:8), averages(10,Ng:Ny+Ng+1), oldaverages(10,Ng:Ny+Ng+1), &
               allaverages(10,Ng:Ny+Ng+1))  ! 39

    allocate(mask(imin:imax,jmin:jmax,kmin:kmax)) ! 40
    call add_2_my_sizer(40,8)

    du=0.0;dv=0.0;dw=0.0
    u=0.0;v=0.0;w=0.0;p=0.0;tmp=0.0;fx=0.0;fy=0.0;fz=0.0;drho=0.0;rho=0.0;mu=0.0;work=0.0;A=0.0
    averages=0.0; oldaverages=0.0; allaverages=0d0
    uold=0.d0;vold=0.d0;wold=0.d0
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
subroutine read_cgd_file(filename, var)

  use module_grid
  use module_flow
  use module_BC

  implicit none
  include 'mpif.h'

  character(len=20) :: filename
  real (8)  , dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: var
  real (8) :: val
  integer  :: i,j,k

  open(unit=12,file=trim(filename),status='old',action='read')
  do i=1,4
    read(12,*); 
  enddo
  do i=1,Nx; do j=1,Ny; do k=1,Nz
    read(12,*) val
    if ((i.ge.is-Ng).and.(i.le.ie)) then 
      if ((j.ge.js-Ng).and.(j.le.je)) then 
        if ((k.ge.ks-Ng).and.(k.le.ke)) then 
          var(i+Ng,j+Ng,k+Ng) = val
        endif
      endif
    endif
  enddo; enddo; enddo
  close(12); 
  call do_all_ghost(var);

end subroutine read_cgd_file
!-------------------------------------------------------------------------------------------------

subroutine InitCondition
  
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
  use module_output_lpp
  use module_output_vof
  use module_freesurface
  use module_mgsolver
  implicit none
  include 'mpif.h'
  integer :: i,j,k, ierr, irank, req(12),sta(MPI_STATUS_SIZE,12)
  real(8) :: my_ave
  real :: erf
  !   real(8) :: NozzleThickness
  !---------------------------------------------Domain----------------------------------------------
  if(rank<nPdomain)then
     if(restart)then
        if ( DoFront ) then 
           call backup_read
        else if ( DoVOF ) then 
           call backup_VOF_read
           call setVOFBC(cvof,vof_flag)
           call do_all_ghost(cvof)
           call do_all_ighost(vof_flag)
           call setVOFBC(cvof,vof_flag)
           call do_ghost_vector(u,v,w)
           if ( do_clean_debris ) call clean_debris
           if ( DoLPP ) call lppvofsweeps(itimestep,time)  
           if (FreeSurface) then
              call press_indices
              if (rank==0) then
                 write(*,'("Restarted simulation with free surface.")') 
                 write(*,'("P_ref, R_ref, R_RK, dR_RK, ddR_RK, V_0: ",6e14.5)')P_ref, R_ref, R_RK, dR_RK, ddR_RK, V_0
              endif
           endif
        end if ! DoFront, DoVOF
        if ( DoLPP ) then 
           call backup_LPP_read
           call SeedParticles
           if ( DoLPP .and. test_injectdrop ) call backup_LPP_write 
        end if ! DoLPP
        call SetVelocityBC(u,v,w,umask,vmask,wmask,time,dt,0)
        call ghost_x(u,2,req( 1: 4)); call ghost_x(v,2,req( 5: 8)); call ghost_x(w,2,req( 9:12))
        call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
        call ghost_y(u,2,req( 1: 4)); call ghost_y(v,2,req( 5: 8)); call ghost_y(w,2,req( 9:12))
        call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
        call ghost_z(u,2,req( 1: 4)); call ghost_z(v,2,req( 5: 8)); call ghost_z(w,2,req( 9:12))
        call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
        write(*,*) "Read backup files done!",rank,time,itimestep
     else
        ! Set velocities and the color function. 
        ! The color function is used for density and viscosity in the domain 
        ! when set by Front-Tracking.
        color = 0.;  v = 0;  w = 0.
        u = U_init;

        !init from cgd file if given in the input file
        if (read_u) then
          call read_cgd_file(u_file, tmp)
          do i=is-1,ie; do j=js,je; do k=ks,ke
            u(i,j,k) = 0.5d0*(tmp(i,j,k) + tmp(i+1,j,k))
          enddo; enddo; enddo
        endif

        if (read_v) then
          call read_cgd_file(v_file, tmp)
          do i=is,ie; do j=js-1,je; do k=ks,ke
            v(i,j,k) = 0.5d0*(tmp(i,j,k) + tmp(i,j+1,k))
          enddo; enddo; enddo
        endif

        if (read_w) then
          call read_cgd_file(w_file, tmp)
          do i=is,ie; do j=js,je; do k=ks-1,ke
            w(i,j,k) = 0.5d0*(tmp(i,j,k) + tmp(i,j,k+1))
          enddo; enddo; enddo
        endif
        
        if(DoLPP) then 
           call SeedParticles
        end if ! DoLPP
        if(DoVOF) then
           call initconditions_VOF()
           call get_all_heights(0)
           if (FreeSurface) call init_FS()
        endif
        du = 0d0

        ! -------------------------------------------------------------------
        ! Initial conditions for different tests
        ! -------------------------------------------------------------------
        ! Test: inject a drop
         if ( test_injectdrop ) then
            do i=imin,imax-1; do j=jmin,jmax-1; do k=kmin,kmax-1
               if((cvof(i,j,k) + cvof(i+1,j,k)) > 0.0d0) u(i,j,k) = 1.5d-1
               if((cvof(i,j,k) + cvof(i,j+1,k)) > 0.0d0) v(i,j,k) =-1.d-1
               if((cvof(i,j,k) + cvof(i,j,k+1)) > 0.0d0) w(i,j,k) =-5.d-1
            enddo; enddo; enddo
         end if
        
        ! -------------------------------------------------------------------
        ! Test: advect a cylinder
         if ( test_cylinder_advection ) then
            do i=imin,imax-1; do j=jmin,jmax-1; do k=kmin,kmax-1
!            if((cvof(i,j,k) + cvof(i+1,j,k)) > 0.0d0) then
            u(i,j,k) = 1.d-2!*cvof(i,j,k)
            v(i,j,k) = 1.d-2!*cvof(i,j,k)
            w(i,j,k) = 0.d0
!           endif
            enddo; enddo; enddo
         endif
     
         if ( test_shear_multiphase ) then
            do i=imin,imax-1; do j=jmin,jmax-1; do k=kmin,kmax-1
               u(i,j,k) = 1.d0*cvof(i,j,k)+15.d0*(1.-cvof(i,j,k))
               v(i,j,k) = 1.d-2*sin(2.d0*PI*x(i))*exp(-(2.d0*(y(i)-0.5d0)/0.1d0)**2)
               w(i,j,k) = 0.d0
            enddo; enddo; enddo
         endif
         !-------------------------------------------------------------------     
         !  Perform tests
         ! -------------------------------------------------------------------

         ! Test: Test height function
         if(test_HF) then
            call test_VOF_HF()
         end if

        ! -------------------------------------------------------------------
        ! Test: Test Lagrangian particle module
         if (test_LP) then 
            call test_Lag_part(itimestep)
         end if ! test_LP
        ! -------------------------------------------------------------------
        ! Test: Test 2d (Quasi-2d) Kelvin-Helmoltz Instability
         if ( test_KHI2D ) then 
            do i=imin,imax; do j=jmin,jmax-1; do k=kmin,kmax-1
               !if( y(j) > yLength*0.5d0 + 0.005*yLength*sin(2.d0*PI*x(i)/xLength)) then 
               !   u(i,j,k) = 1.d1 !*erf( (y(j)-0.5d0*yLength)/(0.1d0*yLength*2.d0) )
               !else
               !   u(i,j,k) = 0.5d1
               !endif
               u(i,j,k) = -0.5d1+(0.5d1+0.5d1)*(1+erf((y(j)-0.5d0*yLength+0.002*xLength* &
                    sin(2.d0*PI*xh(i)/xLength))/(0.008d0*xLength*2.d0)))/2
               ! 2D, only perturb x direction 
            enddo; enddo; enddo
            do i=imin,imax-1; do j=jmin,jmax; do k=kmin,kmax-1
               v(i,j,k) = 1.d0*sin(2.d0*PI*x(i)/xLength)& 
                          *exp(-((yh(j)-yLength*0.5d0)/(0.1d0*xLength))**2.d0)
               ! Nz
               ! perturbation thickness = 0.1*xLength, wavenum(x) = 1
            enddo; enddo; enddo
         end if ! test_KHI_2D 
        ! Test: Test height function (TOMAS)
        !if ( test_KHI2D .or. test_HF) then 
        !   call h_of_KHI2D(itimestep,time)
        !endif
        ! -------------------------------------------------------------------
         if (test_jet ) then 
            ! Test: planar or cylindrical jet with finite length nozzle
            !if ( inject_type == 3 .or. inject_type == 4 ) then
            !   do i = is,ie; do j=js,je; do k = ks,ke
            !      if((cvof(i,j,k) + cvof(i+1,j,k)) > 0.0d0) u(i,j,k) = uliq_inject 
            !   end do; end do; end do
            !end if ! inject_type

            !if ( inject_type == 3 ) then
            !   NozzleThickness = NozzleThick2Cell*dx(is)
            !   do i=imin,imax; do j=jmin,jmax; do k=kmin,kmax
            !      if ( y(j) > radius_liq_inject+NozzleThickness & 
            !         .and. y(j) <= radius_gas_inject ) then
            !         u(i,j,k) = ugas_inject & 
            !            *erf( (y(j) -   radius_liq_inject - NozzleThickness)/blayer_gas_inject ) & 
            !            *erf( (radius_gas_inject - y(j))/blayer_gas_inject )  
            !      end if ! 
            !   end do; end do; end do 
            !end if ! inject_type

            !if ( inject_type == 3 ) then
            !   do i=imin,imax; do j=jmin,jmax; do k=kmin,kmax
            !      if ( y(j) < radius_liq_inject & 
            !         .and. x(i) < radius_gas_inject*10.d0 ) then 
            !         u(i,j,k) = ugas_inject*y(j)/radius_gas_inject
            !      else if ( y(j) > radius_gas_inject & 
            !         .and.  y(j) > radius_liq_inject & 
            !         .and.  x(i) < radius_liq_inject*10.d0 )then
            !         u(i,j,k) = uliq_inject
            !      end if ! y(j)
            !   end do; end do; end do
            !end if ! 
         end if ! test_jet
         ! Test: Shear Droplet
         if (test_shear) then
            do i=imin,imax; do j=jmin,jmax-1; do k=kmin,kmax-1
               u(i,j,k) = 2*(Y(j) - yLength/2)*WallVel(1,1)/yLength
            enddo; enddo; enddo
         endif
         
     
     endif

     if(DoFront) then
        call GetFront('recv')
        call GetFront('wait')
        call Front2GridVector(fx, fy, fz, dIdx, dIdy, dIdz)
        call SetupDensity(dIdx,dIdy,dIdz,A,color)
        if(hypre) then
           call poi_solve(A,color,maxError,maxit,it,HYPRESolverType)
           call ghost_x(color,1,req( 1: 4))
           call ghost_y(color,1,req( 5: 8))
           call ghost_z(color,1,req( 9:12))
           call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
        else
           call NewSolver(A,color,maxError,beta,maxit,it,ierr,ResNormOrderPressure)
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
           call linfunc(rho,rho1,rho2,DensMean)
           call linfunc(mu,mu1,mu2,ViscMean)
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
!  if((.not. restart).or.(.not. restartFront)) then
  if(.not. restart) then
     time = 0d0
     iTimeStep = 0
  endif
  iTimeStepRestart = iTimeStep
  timeLastOutput = DBLE(INT(time/tout))*tout
contains
  subroutine init_FS()
    implicit none
    if (.not. LPP_initialized) then 
      call initialize_LPP()
      LPP_initialized = .true.
      tag_id = 0
    end if !
    call set_topology(vof_phase,itimestep) !vof_phases are updated in initconditions_VOF called above
    call get_ref_volume(cvof)
    call press_indices
    if (RP_test) then
       !call get_ref_volume(cvof)
       call Integrate_RP(dt,time,rho1)
       call initialize_P_RP(p,rho1)  !initialize P field for RP test
       call ghost_x(p,1,req( 1: 4))
       call ghost_y(p,1,req( 5: 8))
       call ghost_z(p,1,req( 9:12))
       call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr) 
    endif
  end subroutine init_FS
end subroutine InitCondition
!=================================================================================================
! function maxabs
!   used for ENO interpolations
!   called in:    subroutine momentumConvection
!                 subroutine density
!-------------------------------------------------------------------------------------------------
function maxabs(a,b)
implicit none
real(8) :: maxabs, a, b
  if(abs(a)>abs(b)) then
    maxabs=a
  else
    maxabs=b
  endif
!  minabs = 0.5*(sign(1.0,abs(b)-abs(a))*(a-b)+a+b)
end function maxabs
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
  if(number>=10**length)call pariserror("int2text error: string is not large enough")
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
  use module_lag_part
  
  implicit none
  include 'mpif.h'
  integer :: in, ierr
  real(8) :: xyzrad(4,10000)
  namelist /parameters/ out_path,      Nx,            Ny,            Nz,            Ng,          &
                        xLength,       yLength,       zLength,       gx,            gy,          &
                        gz,            bdry_cond,     dPdx,          dPdy,          dPdz,        &
                        itime_scheme,  nstep,         maxit,         maxError,                   &
                        beta,          nout,          TwoPhase,      rho1,          mu1,         &
                        rho2,          mu2,           sigma,         BuoyancyCase,  nPx,         &
                        nPy,           nPz,           amin,          amax,          aspmax,      &
                        MaxPoint,      MaxElem,       MaxFront,      xform,         yform,       &
                        zform,         dt,            nregrid,       GetPropertiesFromFront,     &
                        DoVOF,         DoFront,       Implicit,      U_init,        VolumeSource,&
                        CFL,           EndTime,       MaxDt,         smooth,        nsmooth,     &
                        output_format, read_x,        read_y,        read_z,        x_file,      &
                        y_file,        z_file,        restart,       nBackup,       NumBubble,   &
                        xyzrad,        hypre,         dtFlag,        ICOut,         WallVel,     &
                        inject_type,   maxErrorVol,   restartFront,  nstats,        WallShear,   &
                        BoundaryPressure,             ZeroReynolds,  restartAverages,termout,    &  
                        excentricity,  tout,          zip_data,      ugas_inject,   uliq_inject, &  
                        blayer_gas_inject,            tdelay_gas_inject,            padding,     &
                        radius_gas_inject,            radius_liq_inject,     radius_gap_liqgas,  &
                        jetcenter_yc2yLength,         jetcenter_zc2zLength,                      &
                        NozzleThick2Cell,             NozzleLength,                              &
                        cflmax_allowed,               AdvectionScheme, out_mom,   output_fields, & 
                        nsteps_probe,  num_probes,    ijk_probe,                                 &
                        num_probes_cvof,  ijk_probe_cvof,   num_probes_linez, ij_probe_linez,    & 
                        DoTurbStats,   nStepOutputTurbStats, TurbStatsOrder,  timeStartTurbStats,&
                        ResNormOrderPressure,         ErrorScaleHYPRE, DynamicAdjustPoiTol,      & 
                        OutVelSpecified,  MaxFluxRatioPresBC, LateralBdry,                       & 
                        HYPRESolverType, SwitchHYPRESolver, DivergeTol,  uinjectPertAmp,         &
                        ExcludeBoundCellCalcRes, numCellExclude, OrganizeOutFolder,             &  
                        plane, n_p, out_sub, test_MG, MultiGrid, nrelax,u_file, v_file, w_file,  &
                        read_u, read_v, read_w, sym_MG, recordconvergence
 
  Nx = 0; Ny = 4; Nz = 4 ! cause absurd input file that lack nx value to fail. 
  Ng=2;xLength=1d0;yLength=1d0;zLength=1d0
  gx = 0d0; gy=0d0; gz=0d0; bdry_cond = 0
  dPdx = 0d0;  dPdy = 0d0; dPdz = 0d0
  itime_scheme = 1;  nstep = 0; maxit = 50; maxError = 1d-3  
  beta = 1.2; nout = 1 
  TwoPhase = .false.; rho1 = 1d0; mu1 = 0d0
  rho2 = 1d0; mu2 = 0d0; sigma = 0d0; BuoyancyCase = 0; nPx = 1
  nPy = 1; nPz = 1; amin = 0.32; amax = 0.96; aspmax = 1.54
  MaxPoint = 1000000; MaxElem  = 2000000; MaxFront = 100; xform=0d0; yform=0d0; zform=0d0
         dt=0.1d0; nregrid=10; GetPropertiesFromFront = .false.
  DoVOF = .true.;   DoFront = .false.;   Implicit=.false.;   U_init=0d0;   VolumeSource=0d0  
  CFL = 0.5;   EndTime = 0.;  MaxDt = 5d-2;   smooth = .true.;   nsmooth = 20
  output_format = 2;   read_x=.false.;   read_y=.false.;   read_z=.false.
  read_u=.false.;   read_v=.false.;   read_w=.false.
  restart = .false.;  nBackup = 2000;  NumBubble=0
  xyzrad = 0.d0;  hypre=.false.;  dtFlag = 2;  ICout=.false.;   WallVel = 0d0
  inject_type=2 ! redundant
  maxErrorVol=1d-4;   restartfront=.false.;  nstats=10;  WallShear=0d0
  BoundaryPressure=0d0;   ZeroReynolds=.false.;   restartAverages=.false.;   termout=0
  excentricity=0d0;   tout = -1.d0;          zip_data=.false.
  ugas_inject=0.d0;   uliq_inject=0.d0;   uinjectPertAmp = 0.d-2
  blayer_gas_inject=8.d-2; tdelay_gas_inject=1.d-2
  radius_gas_inject=0.1d0; radius_liq_inject=0.1d0 ; radius_gap_liqgas=0d0
  jetcenter_yc2yLength=0.5d0;jetcenter_zc2zLength=0.5d0
  NozzleThick2Cell=2.d0
  NozzleLength = 4.d-3
  padding=5
  cflmax_allowed=0.5d0
  AdvectionScheme = 'QUICK'
  out_mom = .false.
  output_fields = [ .true. , .true. , .true., .true., .true. ]
  nsteps_probe =1; num_probes = 0; ijk_probe = 1; num_probes_cvof = 0; ijk_probe_cvof = 1 
  DoTurbStats = .false.; nStepOutputTurbStats = 1000; TurbStatsOrder = 2
  timeStartTurbStats = 0.d0
  ResNormOrderPressure = 2; ErrorScaleHYPRE = 1.d-2; DynamicAdjustPoiTol=.true. 
  OutVelSpecified = .false.
  MaxFluxRatioPresBC = 0.7d0
  LateralBdry = .false.
  HYPRESolverType = 1; SwitchHYPRESolver = .false.
  DivergeTol = 1.d2
  plane = 0.5d0; n_p = 0.0d0
  out_sub = .false.
  ExcludeBoundCellCalcRes = .false.; numCellExclude = 1
  OrganizeOutFolder = .false.

  in=1
  out=2

  open(unit=in, file='input', status='old', action='read', iostat=ierr)
  if (ierr .ne. 0) call err_no_out_dir("ReadParameters: error opening 'input' file --- perhaps it does not exist ?")
  read(UNIT=in,NML=parameters)
  close(in)
  call check_sanity()
  bdry_read=.true.
  if(MaxFront>10000) call err_no_out_dir("Error: ReadParameters: increase size of xyzrad array")

  if(numBubble>MaxFront) call err_no_out_dir("Error: ReadParameters: increase size of xyzrad array (MaxFront)")

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
     if (OrganizeOutFolder) then 
         call system('mkdir            '//trim(out_path)//'/BACKUP')
         if (DoTurbStats) call system('mkdir  '//trim(out_path)//'/BACKUP_TURB')
         if (DoLPP) then
            call system('mkdir            '//trim(out_path)//'/BACKUPLPP')
            call system('mkdir            '//trim(out_path)//'/ELEMENT_STATS')
         end if ! DoLPP
     end if ! OrganizeOutFloder
     !call system('cp paris.f90 '//trim(out_path))
     open(unit=out, file=trim(out_path)//'/output', action='write', iostat=ierr)
     if (ierr .ne. 0) call pariserror("ReadParameters: error opening output file")
     write(UNIT=out,NML=parameters)
  endif
  call mpi_barrier(MPI_COMM_WORLD, ierr)
  ! Number of grid points in streamwise and transverse directions must
  ! be integral multiples of total number of processors
  if(mod(Nx,nPx) /= 0) call pariserror("ReadParameters: Nx not divisible by nPx!")
  if(mod(Ny,nPy) /= 0) call pariserror("ReadParameters: Ny not divisible by nPy!")
  if(mod(Nz,nPz) /= 0) call pariserror("ReadParameters: Nz not divisible by nPz!")
  Mx = Nx/nPx; My = Ny/nPy; Mz = Nz/nPz
  nPdomain = nPx*nPy*nPz

  ! Check if mesh is uniform in the VOF case
  if(DoVOF.and.rank==0) then
     if ( abs(xLength*dble(Ny)/yLength/dble(Nx) - 1.d0) > 1.d-8 .or. &  
          abs(xLength*dble(Nz)/zLength/dble(Nx) - 1.d0) > 1.d-8 ) & 
          !call pariserror("Mesh is not cubic!")
          print *, "*** WARNING: Mesh is not cubic, and non-uniform VOF is not ready! ***"
  endif

  ! For jet setup
  jetcenter_yc = jetcenter_yc2yLength*yLength
  jetcenter_zc = jetcenter_zc2zLength*zLength

!--- output frequency
  if ( tout > 0.d0 .and. dtFlag == 1) then 
     nout = NINT(tout/dt)
  end if !tout

  if (test_MG) recordconvergence = .TRUE.

  if(termout==0) then
     termout=nout
  endif

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
end subroutine ReadParameters
!=================================================================================================
!=================================================================================================
subroutine parismessage(message) 
  use module_IO
  use module_grid
  implicit none
  include 'mpif.h'
  character(*) :: message
  OPEN(UNIT=88,FILE=TRIM(out_path)//'/message-rank-'//TRIM(int2text(rank,padding))//'.txt')
  write(88,*) message
  if(rank==0) print*,message
  close(88)
end subroutine parismessage
!=================================================================================================
!=================================================================================================
subroutine err_no_out_dir(message) 
! same as pariserror but when out_path directory is not yet created
  use module_IO
  use module_grid
  implicit none
  include 'mpif.h'
  integer :: ierr
  integer :: MPI_errorcode=1
  character(*) :: message
  OPEN(UNIT=89,FILE='error-rank-'//TRIM(int2text(rank,padding))//'.txt')
  write(89,*) message
  if(rank==0) print*,message
  ! Exit MPI gracefully
  close(89)
  close(out)
  call MPI_ABORT(MPI_COMM_WORLD, MPI_errorcode, ierr)
  call MPI_finalize(ierr)
  stop 
end subroutine err_no_out_dir
!=================================================================================================
!=================================================================================================
subroutine pariserror(message) 
  use module_IO
  use module_grid
  implicit none
  include 'mpif.h'
  integer ierr, MPI_errorcode
  character(*) :: message
  OPEN(UNIT=89,FILE=TRIM(out_path)//'/error-rank-'//TRIM(int2text(rank,padding))//'.txt')
  write(89,*) message
  if(rank==0) print*,message
  ! Exit MPI gracefully
  close(89)
  close(out)
  MPI_errorcode=1
  call MPI_ABORT(MPI_COMM_WORLD, MPI_errorcode, ierr)
  call MPI_finalize(ierr)
  stop 
end subroutine pariserror
!=================================================================================================
!=================================================================================================
subroutine check_stability() 
  use module_grid
  use module_flow
  use module_BC
  use module_IO
  use module_front
  
  use module_2phase
  implicit none
  include 'mpif.h'
  real(8) :: von_neumann=10,ststab=10,von_neumann1=10,ststab1=10
  logical :: error=.false.

  if(dt*mu1 /= 0) then 
     von_neumann = dx(ng)**2*rho1/(dt*mu1)
     if(rank==0) print *, "dx**2*rho/(dt*mu)       = ", von_neumann, "dt=",dt
   endif

  if(dt*sigma /= 0) then 
     ststab = dx(ng)**3*rho1/(dt**2*sigma)
     if(rank==0) print *, "dx**3*rho/(dt**2*sigma) = ",ststab , "dt=",dt
  endif

  if(MaxDt*mu1 /= 0) then 
     von_neumann1 = dx(ng)**2*rho1/(MaxDt*mu1)
     if(rank==0) print *, "dx**2*rho/(MaxDt*mu)       = ", von_neumann1, "MaxDt=",MaxDt
  endif

  if(MaxDt*sigma /= 0) then 
     ststab1 = dx(ng)**3*rho1/(MaxDt**2*sigma)
     if(rank==0) print *, "dx**3*rho/(MaxDt**2*sigma) = ", ststab1, "MaxDt=",MaxDt
  endif

  if(von_neumann < 6d0.and..not.Implicit.and.dtFlag==1) then
     call parismessage("time step too large for viscous terms")
     error=.true.
  endif
  if(von_neumann1 < 6d0.and..not.Implicit.and.dtFlag==2) then
     call parismessage("Maximum time step too large for viscous term") 
     error=.true.
  endif
  if(ststab < 2d0.and.dtFlag==1) then
     call parismessage("time step too large for ST term") 
     error=.true.
  endif
  if(ststab1 < 2d0.and.dtFlag==2) then
     call parismessage("Maximum time step too large for ST term") 
     error=.true.
  endif
  if (error) call pariserror("time step too large, aborting")
  return
end subroutine check_stability

subroutine check_integers()
  use module_grid
  implicit none
  integer :: n,big,nbytes,intype
  intype = 4
  print *, "-------------------------"
  print *, "array index integer check"
  print *, "-------------------------"
  nbytes=(mylog2(max(mx,my,mz)))/8
  print *, "log2(box size)", mylog2(max(mx,my,mz))
  nbytes=(3*mylog2(mx)+3)/8
  print *, "nbytes max in box array index", nbytes
  if(nbytes>=intype*8) then
     n=0
     big=1
     print *, "------------"
     print *, "integer check"
     print *, "------------"
     do while (n<=(64/8))
        n=n+1
        big = big*256
        print *,n,big
     end do
     print *, "------------"
     call pariserror("box too large for integer type")
  endif
  contains
    function mylog2(n)
      implicit none
      integer :: mylog2
      integer, intent(in) :: n
      integer :: k
      mylog2=0
      k=n
      if(n<1) call pariserror("mylog2: nonpositive number")
      do while (k>=1) 
         mylog2=mylog2+1
         k=k/2
      end do
      mylog2=mylog2-1
    end function mylog2
end subroutine check_integers

subroutine hello_coucou
  use module_grid
  integer, parameter  :: debug=1
  integer, save :: hello_count=1
  if(debug == 1) then 
  if(rank==0) write(6,*) 'coucou ',hello_count, "Process0"
  if(rank==nPdomain) write(6,*) 'coucou ',hello_count, "Front"
  hello_count = hello_count + 1
  end if
end subroutine hello_coucou

subroutine hello_proc(n)
  use module_grid
  integer, parameter  :: debug=1
  integer, save :: hello_count=1
  integer :: n
  if(debug == 1) then 
  if(rank==n) write(6,*) 'coucou ',hello_count, "Process ", n
  if(rank==nPdomain) write(6,*) 'coucou ',hello_count, "Front"
  hello_count = hello_count + 1
  end if
end subroutine hello_proc


subroutine hello_all
  use module_grid
  integer, parameter  :: debug=1
  integer, save :: hello_count=1
  if(debug == 1) then 
  write(6,*) 'coucou ',hello_count, "Process " , rank
  if(rank==nPdomain) write(6,*) 'coucou ',hello_count, "Front"
  hello_count = hello_count + 1
  end if
end subroutine hello_all

subroutine write_par_var(varname,iout,var)
  use module_grid
  use module_IO
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax) :: var
  integer :: i,j,k,prank,iout
  character(len=10) :: varname

  if (rank==0) then
     OPEN(UNIT=21,FILE=TRIM(varname)//'-holder-'//TRIM(int2text(iout,padding))//'.txt')
     write(21,12)NpDomain
     do prank=0,NpDomain-1
        write(21,13)TRIM(out_path)//'/'//TRIM(varname)//TRIM(int2text(prank,padding))//'-'//TRIM(int2text(iout,padding))//'.3D'
     enddo
     close(21)
  endif
  OPEN(UNIT=20,FILE=TRIM(out_path)//'/'//TRIM(varname)//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout,padding))//'.3D')
  !write(20,13)'X Y Z Var'
  do k=ks,ke; do j=js,je; do i=is,ie
     write(20,14)x(i),y(j),z(k),var(i,j,k)
  enddo;enddo;enddo
  close(20)

12 format('NPROCS',I4)
13 format(A)
14 format(4e17.8)

end subroutine write_par_var
