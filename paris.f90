!=================================================================================================
!=================================================================================================
! Paris-0.1
! Extended from Code: FTC3D2011 (Front Tracking Code for 3D simulations)
! 
! Authors: Sadegh Dabiri, Gretar Tryggvason
! author for minor extenstions Stephane Zaleski(zaleski@dalembert.upmc.fr) 
! Contact: sdabiri@gmail.com
! A three dimensional Navier-Stokes flow solver for modeling of multiphase 
! flows. Flow can be driven by wall motion, density difference or pressure gradient.
! Boundary conditions supported: wall and periodic
! Version 1.0   1/21/2011   The 3D flow solver for variable density/viscosity is written. 
!                           The density is advected by ENO scheme.
! Version 2.0   2/25/2011   Parallel implementation.
!=================================================================================================
!=================================================================================================
!=================================================================================================
! Main program to solve the 3D NS equations for multiphase flows
!-------------------------------------------------------------------------------------------------
Program paris
  use module_flow
  use module_grid
  use module_BC
  use module_tmpvar
  use module_poisson
  use module_IO
  use module_solids
  use module_output_solids
  implicit none
  include 'mpif.h'
  integer :: ierr, i,j,k
  ! Locals for marching and timing
  real(8) :: start_time, end_time !, IC_time, wall_time
  integer :: req(48),sta(MPI_STATUS_SIZE,48)
  logical, allocatable, dimension(:,:,:) :: mask
!---------------------------------------INITIALIZATION--------------------------------------------
  ! Initialize MPI
  call MPI_INIT(ierr)
  If (ierr /= 0) STOP '*** Main: unsuccessful MPI-initialization'
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numProcess, ierr)
  start_time = MPI_Wtime(ierr)

  call ReadParameters
  call ReadSolidParameters
  if(rank==0) write(out,*)'Parameters read successfully'
  ! check number of processors
  If (numProcess .NE. nPx*nPy*nPz) STOP '*** Main: Problem with nPx!'

  call initialize
  if(dosolids) then
     call initsolids
     call output_solids(0,imin,imax,jmin,jmax,kmin,kmax)
     if(rank==0) write(out,*)'solids initialized'
     if(rank==0) write(6,*)'solids initialized'
  endif
  call InitCondition
  if(rank==0) write(out,*)'initialized'
  if(rank==0) write(6,*)'initialized'
  if(HYPRE) call poi_initialize
  if(HYPRE .and. rank==0) write(out,*)'hypre initialized'

  allocate(mask(imin:imax,jmin:jmax,kmin:kmax))
  ! output initial condition
  if(ICOut)call output(0,is,ie,js,je,ks,ke)
  if(rank==0)start_time = MPI_WTIME()
!---------------------------------------MAIN TIME LOOP--------------------------------------------
  do while(time<EndTime .and. itimestep<nstep)
    if(dtFlag==2)call TimeStepSize(dt)
    time=time+dt
    itimestep=itimestep+1
    if(rank==0)write(out,'("Step: ",I6," dt=",es16.5e2," time=",es16.5e2)')itimestep,dt,time
    if(itime_scheme==2) then
      uold = u
      vold = v
      wold = w
      rhoo = rho
      muold  = mu
    endif
    if(dosolids) then
       u = u*(1.d0 -solids) 
       v = v*(1.d0 -solids) 
       w = w*(1.d0 -solids) 
    endif
!------------------------------------ADVECTION & DIFFUSION----------------------------------------
    do ii=1, itime_scheme
      call momentumDiffusion(u,v,w,rho,mu,du,dv,dw)
      call momentumConvection(u,v,w,du,dv,dw)
      call volumeForce(rho,rho1,rho2,dpdx,dpdy,dpdz,BuoyancyCase,fx,fy,fz,gx,gy,gz,du,dv,dw)
      if(TwoPhase) call density(rho,u,v,w,drho)
      u = u + dt * du
      v = v + dt * dv
      w = w + dt * dw
      if(TwoPhase) rho = rho + dt * drho

      call SetVelocityBC(u,v,w)
      
      call ghost_x(u  ,2,req( 1: 4));  call ghost_x(v,2,req( 5: 8)); call ghost_x(w,2,req( 9:12)); 
      call ghost_x(rho,2,req(13:16));  call MPI_WAITALL(16,req(1:16),sta(:,1:16),ierr)
      call ghost_y(u  ,2,req( 1: 4));  call ghost_y(v,2,req( 5: 8)); call ghost_y(w,2,req( 9:12)); 
      call ghost_y(rho,2,req(13:16));  call MPI_WAITALL(16,req(1:16),sta(:,1:16),ierr)
      call ghost_z(u  ,2,req( 1: 4));  call ghost_z(v,2,req( 5: 8)); call ghost_z(w,2,req( 9:12)); 
      call ghost_z(rho,2,req(13:16));  call MPI_WAITALL(16,req(1:16),sta(:,1:16),ierr)

      if(TwoPhase) mu  = (rho-rho1)/(rho2-rho1)*(mu2-mu1)+mu1
    enddo

    if(itime_scheme==2) then
      u = 0.5*(u+uold)
      v = 0.5*(v+vold)
      w = 0.5*(w+wold)
      if(TwoPhase) then
        rho = 0.5*(rho+rhoo)
        mu  = 0.5*(mu +muold )
      endif
    endif
!--------------------------------------PROJECTION STEP--------------------------------------------
    ! Compute source term and the coefficient do p(i,j,k)
!    tmp(is:ie,js:je,ks:ke)=rho(is:ie,js:je,ks:ke)
!    tmp = rho
!    call SetPressureBC(density=tmp,pressure=p)
!    call SetupPoisson(u,v,w,tmp,dt,work)

    tmp = rho
    call SetPressureBC(tmp)
    if(HYPRE)then
      call poi_solve(u,v,w,tmp,dt,p,maxError,maxit)
      call ghost_x(p,1,req(1:4 ))
      call ghost_y(p,1,req(5:8 ))
      call ghost_z(p,1,req(9:12))
      call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
    else
      mask=.true.
      call SetupPoisson(u,v,w,tmp,dt,work,mask)
      call SOR_Solver(work,p,maxError,beta,maxit,it)
    endif
    do k=ks,ke;  do j=js,je; do i=is,ieu;    ! CORRECT THE u-velocity 
      u(i,j,k)=u(i,j,k)-dt*(2.0/dxh(i))*(p(i+1,j,k)-p(i,j,k))/(rho(i+1,j,k)+rho(i,j,k))
    enddo; enddo; enddo
      
    do k=ks,ke;  do j=js,jev; do i=is,ie;    ! CORRECT THE v-velocity
      v(i,j,k)=v(i,j,k)-dt*(2.0/dyh(j))*(p(i,j+1,k)-p(i,j,k))/(rho(i,j+1,k)+rho(i,j,k))
    enddo; enddo; enddo
    
    do k=ks,kew;  do j=js,je; do i=is,ie;   ! CORRECT THE v-velocity
      w(i,j,k)=w(i,j,k)-dt*(2.0/dzh(k))*(p(i,j,k+1)-p(i,j,k))/(rho(i,j,k+1)+rho(i,j,k))
    enddo; enddo; enddo

    call SetVelocityBC(u,v,w)

    call ghost_x(u,2,req( 1: 4)); call ghost_x(v,2,req( 5: 8)); call ghost_x(w,2,req( 9:12))
    call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
    call ghost_y(u,2,req( 1: 4)); call ghost_y(v,2,req( 5: 8)); call ghost_y(w,2,req( 9:12))
    call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
    call ghost_z(u,2,req( 1: 4)); call ghost_z(v,2,req( 5: 8)); call ghost_z(w,2,req( 9:12))
    call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)

!--------------------------------------------OUTPUT-----------------------------------------------
    ! calculate the volumetric flow rate in x-direction
    i=is
    myfr=0d0
    do k=ks,ke;  do j=js,je
      myfr=myfr+u(i,j,k)*dy(j)*dz(k)
    enddo; enddo
    call MPI_ALLREDUCE(myfr, flowrate, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    flowrate = flowrate/dfloat(nPx)
    if(mod(itimestep,nbackup)==0)call backup_write
    if(mod(itimestep,nout)==0) then
       call output(ITIMESTEP/nout,is,ie,js,je,ks,ke)
    endif
    if(rank==0)then
      end_time =  MPI_WTIME()
      write(out,'("Step:",I6," Iterations:",I5," cpu(s):",f10.2)')itimestep,it,end_time-start_time
      open(unit=121,file='track')
      write(121,'("Step:",I6," dt=",es16.5e2," time=",es16.5e2)')itimestep,dt,time
      write(121,'("            Iterations:",I5," cpu(s):",f10.2)')it,end_time-start_time
      close(121)
      write(*,'("Step:",I6," dt=",es16.5e2," time=",es16.5e2)')itimestep,dt,time
      write(*,'("            Iterations:",I5," cpu(s):",f10.2)')it,end_time-start_time
      open(unit=121,file='stats',access='append')
      write(121,'(5es16.5e2)')time,flowrate,end_time-start_time
      close(121)
    endif
  enddo
!--------------- END OF MAIN TIME LOOP ----------------------------------------------------------
  if(rank==0) call close_visit_file()
  call output_at_location()
  if(HYPRE) call poi_finalize
  call MPI_FINALIZE(ierr)
  stop
end program paris
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
subroutine TimeStepSize(deltaT)
  use module_grid
  use module_flow
  implicit none
  include "mpif.h"
  real(8) :: deltaT, h, vmax, dtadv, mydt !, vmaxglob
  real(8), save :: dtvisc, nu_min, nu_max !, hglob
!  logical, save :: first_time=.true.
  integer :: ierr
  h  = min(minval(dx),minval(dy),minval(dz))
  nu_min = mu1/rho1;   if(TwoPhase) nu_min = min(nu_min,mu2/rho2)
  nu_max = mu1/rho1;   if(TwoPhase) nu_max = max(nu_max,mu2/rho2)
  dtvisc=h**2/(6.0d0*nu_max)
  vmax = maxval(sqrt(u(is:ie,js:je,ks:ke)**2+v(is:ie,js:je,ks:ke)**2+w(is:ie,js:je,ks:ke)**2))
  vmax = max(vmax,1e-3)
!  dtadv  = min(h/vmax,2d0*nu_min/vmax**2)
  dtadv  = h/vmax
  mydt = CFL*min(dtvisc,dtadv)
  mydt = min(mydt,MaxDt)

  call MPI_ALLREDUCE(mydt, deltat, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
end subroutine TimeStepSize
!=================================================================================================
!=================================================================================================
! subroutine momentumConvection
! calculates the convection terms in mumentum equation using ENO scheme
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
    du(i,j,k)=du(i,j,k)-0.5*((u(i,j  ,k  )+u(i+1,j  ,k  ))*work(i+1,j  ,k  ,1)-&
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
    dv(i,j,k)=dv(i,j,k)-0.5*((u(i  ,j,k  )+u(i  ,j+1,k  ))*work(i  ,j  ,k  ,1)-&
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
    dw(i,j,k)=dw(i,j,k)-0.5*((u(i  ,j  ,k)+u(i  ,j  ,k+1))*work(i  ,j  ,k  ,1)-&
                    (u(i-1,j  ,k)+u(i-1,j  ,k+1))*work(i-1,j  ,k  ,1))/dx(i)  &
              -0.5*((v(i  ,j  ,k)+v(i  ,j  ,k+1))*work(i  ,j  ,k  ,2)-&
                    (v(i  ,j-1,k)+v(i  ,j-1,k+1))*work(i  ,j-1,k  ,2))/dy(j)  &
              -0.5*((w(i  ,j  ,k)+w(i  ,j  ,k+1))*work(i  ,j  ,k+1,3)-&
                    (w(i  ,j  ,k)+w(i  ,j  ,k-1))*work(i  ,j  ,k  ,3))/dzh(k)
  enddo; enddo; enddo

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
end subroutine momentumDiffusion
!=================================================================================================
!=================================================================================================
! Calculates the volume force in the momentum equations and adds them to du,dv,dw
!-------------------------------------------------------------------------------------------------
subroutine volumeForce(rho,rho1,rho2,dpdx,dpdy,dpdz,BuoyancyCase,fx,fy,fz,gx,gy,gz,du,dv,dw)
  use module_grid
  implicit none
  include "mpif.h"
  integer, intent(in) :: BuoyancyCase
  real(8), intent(in) :: gx,gy,gz,rho1,rho2, dpdx, dpdy, dpdz
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: rho, fx,fy,fz
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: du, dv, dw
  real(8) :: vol, rro, ro
  integer :: i,j,k, ierr
  if(BuoyancyCase==0) then
    rro = 0.0
  elseif(BuoyancyCase==1) then
    rro=rho1
  elseif(BuoyancyCase==2) then
    rro=rho2
  elseif(BuoyancyCase==3) then
    ! calculate average density
    vol=sum(dx(is:ie))*sum(dy(js:je))*sum(dz(ks:ke))
    ro=0.0
    do k=ks,ke; do j=js,je; do i=is,ie
      ro = ro+dx(i)*dy(j)*dz(k)*rho(i,j,k)
    enddo; enddo; enddo
    ro=ro/(xLength*yLength*zLength)
    call MPI_ALLREDUCE(ro, rro, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_World, ierr)
  else
    stop 'volumeForce: invalid buoyancy option'
  endif
   
  do k=ks,ke;  do j=js,je; do i=is,ieu
    du(i,j,k)=du(i,j,k)+(fx(i,j,k)-dpdx)/(0.5*(rho(i+1,j,k)+rho(i,j,k))) + &
        (1.0 -rro/(0.5*(rho(i+1,j,k)+rho(i,j,k))) )*gx
  enddo; enddo; enddo

  do k=ks,ke;  do j=js,jev; do i=is,ie
    dv(i,j,k)=dv(i,j,k)+(fy(i,j,k)-dpdy)/(0.5*(rho(i,j+1,k)+rho(i,j,k))) + &
        (1.0 -rro/(0.5*(rho(i,j+1,k)+rho(i,j,k))) )*gy
  enddo; enddo; enddo

  do k=ks,kew;  do j=js,je; do i=is,ie
    dw(i,j,k)=dw(i,j,k)+(fz(i,j,k)-dpdz)/(0.5*(rho(i,j,k+1)+rho(i,j,k))) + &
        (1.0 -rro/(0.5*(rho(i,j,k+1)+rho(i,j,k))) )*gz
  enddo; enddo; enddo
end subroutine volumeForce
!=================================================================================================
!=================================================================================================
! The Poisson equation for the pressure is setup with matrix A as
! A7*Pijk = A1*Pi-1jk + A2*Pi+1jk + A3*Pij-1k + 
!           A4*Pij+1k + A5*Pijk-1 + A6*Pijk+1 + A8
!-------------------------------------------------------------------------------------------------
subroutine SetupPoisson(utmp,vtmp,wtmp,rhot,dt,A,mask)
  use module_grid
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: utmp,vtmp,wtmp,rhot
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax,8), intent(out) :: A
  logical, intent(in) :: mask(imin:imax,jmin:jmax,kmin:kmax)
  real(8) :: dt
  integer :: i,j,k
  do k=ks,ke; do j=js,je; do i=is,ie;
    if(mask(i,j,k))then
      A(i,j,k,1) = (1./dx(i))*( 1./(dxh(i-1)*(rhot(i-1,j,k)+rhot(i,j,k))) )
      A(i,j,k,2) = (1./dx(i))*( 1./(dxh(i  )*(rhot(i+1,j,k)+rhot(i,j,k))) )
      A(i,j,k,3) = (1./dy(j))*( 1./(dyh(j-1)*(rhot(i,j-1,k)+rhot(i,j,k))) )
      A(i,j,k,4) = (1./dy(j))*( 1./(dyh(j  )*(rhot(i,j+1,k)+rhot(i,j,k))) )
      A(i,j,k,5) = (1./dz(k))*( 1./(dzh(k-1)*(rhot(i,j,k-1)+rhot(i,j,k))) )
      A(i,j,k,6) = (1./dz(k))*( 1./(dzh(k  )*(rhot(i,j,k+1)+rhot(i,j,k))) )
      A(i,j,k,7) = sum(A(i,j,k,1:6))
      A(i,j,k,8) = -0.5/dt *( (utmp(i,j,k)-utmp(i-1,j,k))/dx(i) &
                             +(vtmp(i,j,k)-vtmp(i,j-1,k))/dy(j) &
                             +(wtmp(i,j,k)-wtmp(i,j,k-1))/dz(k)   )
    endif
  enddo; enddo; enddo
end subroutine SetupPoisson
!=================================================================================================
!=================================================================================================
! SOR sover to solve the following linear equiation:
! A7*Pijk = A1*Pi-1jk + A2*Pi+1jk + A3*Pij-1k + 
!           A4*Pij+1k + A5*Pijk-1 + A6*Pijk+1 + A8
!-------------------------------------------------------------------------------------------------
subroutine SOR_Solver(A,p,maxError,beta,maxit,it)
  use module_grid
  use module_BC
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: p
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax,8), intent(in) :: A
  real(8), intent(in) :: beta, maxError
  integer, intent(in) :: maxit
  integer, intent(out) :: it
  real(8) :: res, totalres
  integer :: i,j,k, ierr
  integer :: req(12),sta(MPI_STATUS_SIZE,12)
  logical :: mask(imin:imax,jmin:jmax,kmin:kmax)
!--------------------------------------ITTERATION LOOP--------------------------------------------  
  do it=1,maxit
    do k=ks,ke; do j=js,je; do i=is,ie
      p(i,j,k)=(1.0-beta)*p(i,j,k)+beta* 1.0/A(i,j,k,7)*(              &
        A(i,j,k,1) * p(i-1,j,k) + A(i,j,k,2) * p(i+1,j,k) +            &
        A(i,j,k,3) * p(i,j-1,k) + A(i,j,k,4) * p(i,j+1,k) +            &
        A(i,j,k,5) * p(i,j,k-1) + A(i,j,k,6) * p(i,j,k+1) + A(i,j,k,8))
    enddo; enddo; enddo
!---------------------------------CHECK FOR CONVERGENCE-------------------------------------------
    res = 0.0
      call ghost_x(p,1,req( 1: 4)); call ghost_y(p,1,req( 5: 8)); call ghost_z(p,1,req( 9:12))
      do k=ks+1,ke-1; do j=js+1,je-1; do i=is+1,ie-1
        res=res+abs(-p(i,j,k) * A(i,j,k,7) +                             &
          A(i,j,k,1) * p(i-1,j,k) + A(i,j,k,2) * p(i+1,j,k) +            &
          A(i,j,k,3) * p(i,j-1,k) + A(i,j,k,4) * p(i,j+1,k) +            &
          A(i,j,k,5) * p(i,j,k-1) + A(i,j,k,6) * p(i,j,k+1) + A(i,j,k,8) )
      enddo; enddo; enddo
      call MPI_WAITALL(12,req,sta,ierr)
      mask=.true.
      mask(is+1:ie-1,js+1:je-1,ks+1:ke-1)=.false.
      do k=ks,ke; do j=js,je; do i=is,ie
        if(mask(i,j,k)) res=res+abs(-p(i,j,k) *  A(i,j,k,7) +            &
          A(i,j,k,1) * p(i-1,j,k) + A(i,j,k,2) * p(i+1,j,k) +            &
          A(i,j,k,3) * p(i,j-1,k) + A(i,j,k,4) * p(i,j+1,k) +            &
          A(i,j,k,5) * p(i,j,k-1) + A(i,j,k,6) * p(i,j,k+1) + A(i,j,k,8) )
      enddo; enddo; enddo
    call MPI_ALLREDUCE(res, totalres, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    totalres = totalres/float(Nx*Ny*Nz)
    if (.not.(totalres<1e10)) stop '***** solution has diverged *****'
    if (totalres<maxError) exit
  enddo
  if(it==maxit+1 .and. rank==0) write(*,*) 'Warning: SOR_Solver reached maxit.'
end subroutine SOR_Solver
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

  allocate(dims(ndim),periodic(ndim),reorder(ndim),coords(ndim),STAT=ierr)
  dims(1) = nPx; dims(2) = nPy; dims(3) = nPz

  periodic = 0; do i=1,3; if (bdry_cond(i) =='periodic') periodic(i) = 1; enddo
  reorder = 1
  call MPI_Cart_Create(MPI_Comm_World,ndim,dims,periodic,reorder,MPI_Comm_Cart,ierr)
  if (ierr /= 0) STOP '*** Grid: unsuccessful Cartesian MPI-initialization'

  call MPI_Cart_Coords(MPI_Comm_Cart,rank,ndim,coords,ierr)
  if (ierr /= 0) STOP '*** Grid: unsuccessful in getting topology coordinates'

  is=coords(1)*Mx+1+Ng; ie=coords(1)*Mx+Mx+Ng; imin=is-Ng; imax=ie+Ng
  js=coords(2)*My+1+Ng; je=coords(2)*My+My+Ng; jmin=js-Ng; jmax=je+Ng
  ks=coords(3)*Mz+1+Ng; ke=coords(3)*Mz+Mz+Ng; kmin=ks-Ng; kmax=ke+Ng
  ieu=ie; if(bdry_cond(1)/='periodic' .and. coords(1)==nPx-1) ieu=ie-1
  jev=je; if(bdry_cond(2)/='periodic' .and. coords(2)==nPy-1) jev=je-1
  kew=ke; if(bdry_cond(3)/='periodic' .and. coords(3)==nPz-1) kew=ke-1

  Nxt=Nx+2*Ng; Nyt=Ny+2*Ng; Nzt=Nz+2*Ng ! total number of cells

!  allocate( x(imin:imax),dx(imin:imax),dxh(imin:imax),&
!            y(jmin:jmax),dy(jmin:jmax),dyh(jmin:jmax),&
!            z(kmin:kmax),dz(kmin:kmax),dzh(kmin:kmax) )
  allocate(x(Nxt), xh(Nxt), dx(Nxt), dxh(Nxt), &
           y(Nyt), yh(Nyt), dy(Nyt), dyh(Nyt), &
           z(Nzt), zh(Nzt), dz(Nzt), dzh(Nzt)  )

  allocate(  u(imin:imax,jmin:jmax,kmin:kmax), uold(imin:imax,jmin:jmax,kmin:kmax), &
            du(imin:imax,jmin:jmax,kmin:kmax),   fx(imin:imax,jmin:jmax,kmin:kmax), &
             v(imin:imax,jmin:jmax,kmin:kmax), vold(imin:imax,jmin:jmax,kmin:kmax), &
            dv(imin:imax,jmin:jmax,kmin:kmax),   fy(imin:imax,jmin:jmax,kmin:kmax), &
             w(imin:imax,jmin:jmax,kmin:kmax), wold(imin:imax,jmin:jmax,kmin:kmax), &
            dw(imin:imax,jmin:jmax,kmin:kmax),   fz(imin:imax,jmin:jmax,kmin:kmax), &
           rho(imin:imax,jmin:jmax,kmin:kmax), rhoo(imin:imax,jmin:jmax,kmin:kmax), &
          drho(imin:imax,jmin:jmax,kmin:kmax),    p(imin:imax,jmin:jmax,kmin:kmax), &
            mu(imin:imax,jmin:jmax,kmin:kmax),muold(imin:imax,jmin:jmax,kmin:kmax)  )

  allocate(tmp(imin:imax,jmin:jmax,kmin:kmax), work(imin:imax,jmin:jmax,kmin:kmax,8))

  x=0.0;y=0.0;z=0.0;dx=0.0;dy=0.0;dz=0.0;dxh=0.0;dyh=0.0;dzh=0.0;du=0.0;dv=0.0;dw=0.0
  u=0.0;v=0.0;w=0.0;p=0.0;tmp=0.0;fx=0.0;fy=0.0;fz=0.0;drho=0.0;rho=0.0;mu=0.0;work=00.

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
!=================================================================================================
! subroutine InitCondition
!   Sets the initial conditions
!   called in:    program paris
!-------------------------------------------------------------------------------------------------
subroutine InitCondition
  use module_grid
  use module_flow
  use module_BC
  use module_IO
  use module_tmpvar

  implicit none
  include 'mpif.h'
  integer :: ierr, i,j,k, ib
  integer :: req(48),sta(MPI_STATUS_SIZE,48)

  ! Set the wall velocities
  WallVel=0.0

  if(restart)then
    call backup_read
    call ghost_x(u  ,2,req( 1: 4));  call ghost_x(v,2,req( 5: 8)); call ghost_x(w,2,req( 9:12)); 
    call ghost_x(rho,2,req(13:16));  call MPI_WAITALL(16,req(1:16),sta(:,1:16),ierr)
    call ghost_y(u  ,2,req( 1: 4));  call ghost_y(v,2,req( 5: 8)); call ghost_y(w,2,req( 9:12)); 
    call ghost_y(rho,2,req(13:16));  call MPI_WAITALL(16,req(1:16),sta(:,1:16),ierr)
    call ghost_z(u  ,2,req( 1: 4));  call ghost_z(v,2,req( 5: 8)); call ghost_z(w,2,req( 9:12)); 
    call ghost_z(rho,2,req(13:16));  call MPI_WAITALL(16,req(1:16),sta(:,1:16),ierr)
    mu = mu1
    if(TwoPhase) mu  = (rho-rho1)/(rho2-rho1)*(mu2-mu1)+mu1
  else
    time = 0d0
    itimestep = 0
    rho=rho1; mu=mu1
  ! Set density and viscosity in the domain and the drop
    if(TwoPhase)then
      do ib=1,numBubble
        do i=imin,imax; do j=jmin,jmax; do k=kmin,kmax; 
          if ( (x(i)-xc(ib))**2+(y(j)-yc(ib))**2+(z(k)-zc(ib))**2 < rad(ib)**2) then 
             rho(i,j,k)=rho2; mu(i,j,k)=mu2; 
          endif 
        enddo; enddo; enddo
      enddo
    endif
  endif

end subroutine InitCondition
!=================================================================================================
!=================================================================================================
! subroutine density
!   calculates drho/dt
!   called in:    program paris
!-------------------------------------------------------------------------------------------------
subroutine density(rho,u,v,w,drho)
  use module_grid
  use module_tmpvar
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: u, v, w,rho !,D
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(out) :: drho
  real(8), external :: minabs
  integer :: i,j,k

  do k=ks,ke; do j=js,je; do i=is-1,ie
    if (u(i,j,k)>0.0) then
      work(i,j,k,1)=rho(i  ,j,k)+0.5*minabs((rho(i+1,j,k)-rho(i  ,j,k)),(rho(i  ,j,k)-rho(i-1,j,k)))
    else
      work(i,j,k,1)=rho(i+1,j,k)-0.5*minabs((rho(i+2,j,k)-rho(i+1,j,k)),(rho(i+1,j,k)-rho(i  ,j,k)))
    endif
  enddo; enddo; enddo
  do k=ks,ke; do j=js-1,je; do i=is,ie
    if(v(i,j,k)>0.0) then
      work(i,j,k,2)=rho(i,j  ,k)+0.5*minabs((rho(i,j+1,k)-rho(i,j  ,k)),(rho(i,j  ,k)-rho(i,j-1,k)))
    else
      work(i,j,k,2)=rho(i,j+1,k)-0.5*minabs((rho(i,j+2,k)-rho(i,j+1,k)),(rho(i,j+1,k)-rho(i,j  ,k)))
    endif
  enddo; enddo; enddo
  do k=ks-1,ke; do j=js,je; do i=is,ie
    if(w(i,j,k)>0.0) then
      work(i,j,k,3)=rho(i,j,k  )+0.5*minabs((rho(i,j,k+1)-rho(i,j,k  )),(rho(i,j,k  )-rho(i,j,k-1)))
    else
      work(i,j,k,3)=rho(i,j,k+1)-0.5*minabs((rho(i,j,k+2)-rho(i,j,k+1)),(rho(i,j,k+1)-rho(i,j,k  )))
    endif
  enddo; enddo; enddo
  do k=ks,ke;  do j=js,je; do i=is,ie
    drho(i,j,k) = -( u(i  ,j  ,k  )*work(i  ,j  ,k  ,1)          &
                    -u(i-1,j  ,k  )*work(i-1,j  ,k  ,1) )/dx(i)  &
                  -( v(i  ,j  ,k  )*work(i  ,j  ,k  ,2)          &
                    -v(i  ,j-1,k  )*work(i  ,j-1,k  ,2) )/dy(j)  &
                  -( w(i  ,j  ,k  )*work(i  ,j  ,k  ,3)          &
                    -w(i  ,j  ,k-1)*work(i  ,j  ,k-1,3) )/dz(k)  !&
!          + 0.5*( ( (D(i,j,k)+D(i+1,j,k))*(rho(i+1,j,k)-rho(i,j,k))/dxh(i  )         &
!                  - (D(i,j,k)+D(i-1,j,k))*(rho(i,j,k)-rho(i-1,j,k))/dxh(i-1) )/dx(i) &
!                 +( (D(i,j,k)+D(i,j+1,k))*(rho(i,j+1,k)-rho(i,j,k))/dyh(j  )         &
!                  - (D(i,j,k)+D(i,j-1,k))*(rho(i,j,k)-rho(i,j-1,k))/dyh(j-1) )/dy(j) &
!                 +( (D(i,j,k)+D(i,j,k+1))*(rho(i,j,k+1)-rho(i,j,k))/dzh(k  )         &
!                  - (D(i,j,k)+D(i,j,k-1))*(rho(i,j,k)-rho(i,j,k-1))/dzh(k-1) )/dz(k) )
  enddo; enddo; enddo

end subroutine density
!=================================================================================================
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
! subroutine ReadParameters
!   Reads the parameters in "input" file for the calculation
!   called in:    program paris
!-------------------------------------------------------------------------------------------------
subroutine ReadParameters
  use module_grid
  use module_flow
  use module_BC
  use module_IO
  implicit none
  include 'mpif.h'
  integer :: in, ierr
  integer, parameter :: MaxFront=1000
  real(8) :: xyzrad(4,MaxFront)
  namelist /parameters/ out_path , Nx, Ny, Nz, Ng, xLength, yLength, zLength, gx, gy, gz, hypre, &
                        bdry_cond, dpdx, dpdy, dpdz, itime_scheme, nstep, maxit, maxError,  &
                        beta, nout, TwoPhase, rho1, mu1, rho2, mu2, sigma,BuoyancyCase, &
                        nPx,nPy,nPz, EndTime, MaxDt, CFL, xform, yform, zform, output_format, &
                        nbackup, x_file, y_file, z_file, read_x, read_y, read_z, restart, dt, &
                        dtFlag, ICOut, NumBubble, xyzrad

  in=1
  out=2

  open(unit=in, file='input', status='old', action='read', iostat=ierr)
  if (ierr .ne. 0) stop 'ReadParameters: error opening input file'
  read(UNIT=in,NML=parameters)
  close(in)

  if(numBubble>MaxFront) stop 'Error: ReadParameters: increase size of xyzrad array (MaxFront)'
  allocate(rad(MaxFront), xc(MaxFront), yc(MaxFront), zc(MaxFront))
  xc (1:NumBubble) = xyzrad(1,1:NumBubble)
  yc (1:NumBubble) = xyzrad(2,1:NumBubble)
  zc (1:NumBubble) = xyzrad(3,1:NumBubble)
  rad(1:NumBubble) = xyzrad(4,1:NumBubble)

  if(rank==0)then
    call system('mkdir            '//trim(out_path))
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
  if(numProcess /= nPx*nPy*nPz)then
    if(rank==0) write(out,*) 'Error: incorrect number of processors.'
    call MPI_Finalize(ierr)
    stop
  endif
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
end subroutine ReadParameters
!=================================================================================================
!=================================================================================================
