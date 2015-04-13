!=================================================================================================
!=================================================================================================
!
! PARIS (Parallel Robust Interface Software) Simulator 
! Most common modules file
!
! Contact: Stephane Zaleski zaleski@dalembert.upmc.fr
! 
! Authors (in alphabetical order):
!         Thomas Arrufat Jackson
!         Sadegh Dabiri      
!         Daniel Fuster      
! 	  Yue "Stanley" Ling 
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
! module_grid: Contains definition of variables for the grid.
!-------------------------------------------------------------------------------------------------
module module_grid

#ifdef __GFORTRAN__
#elif __INTEL_COMPILER
#else 
  ! If you know the compiler version macro for your compiler add it here and 
  ! notify the authors
#endif

  implicit none
  save
  integer :: Nx, Ny, Nz, Ng ! Ng is the number of ghost cells
  integer :: Nxt, Nyt, Nzt ! total number of cells
  integer :: is, ie, js, je, ks, ke
  integer :: ieu, jev, kew
  real(8), dimension(:), allocatable :: x, xh, dx, dxh
  real(8), dimension(:), allocatable :: y, yh, dy, dyh
  real(8), dimension(:), allocatable :: z, zh, dz, dzh
  real(8) :: xLength, yLength, zLength, xform, yform, zform !non-uniformity of grid
  logical :: TwoPhase
  integer :: nPx, nPy, nPz, Mx, My, Mz, rank, ndim=3, nPdomain, NumProcess
  integer, dimension(:), allocatable :: dims, coords, periodic, reorder
  integer :: MPI_Comm_Cart, MPI_Comm_Domain, MPI_Comm_Active
  integer :: imin, imax, jmin, jmax, kmin, kmax
  real(8), parameter :: TINY_DOUBLE=1d-300 
! added by SZ
contains
!
  function test_point_in(i,j,k)
    integer, intent(in) :: i,j,k
    logical :: test_point_in
    test_point_in = (is <= i).and.(i <= ie).and.(js <= j).and.(j <= je).and.(ks <= k).and.(k <= ke)
  end function test_point_in
!
  function test_point_in_wide(i,j,k)
    integer, intent(in) :: i,j,k
    logical :: test_point_in_wide
    test_point_in_wide = (imin < i).and.(i < imax).and.(jmin < j).and.(j < jmax).and.(kmin < k).and.(k < kmax)
  end function test_point_in_wide
!
  subroutine check_sanity_in_depth()
    integer :: v12345678912345a=1
    integer :: v12345678912345b=2
    call check_sanity()
    if(is < 1) call pariserror( "wrong is")
    if(js < 1) call pariserror( "wrong js")
    if(ks < 1) call pariserror( "wrong ks")
! paranoid programming
    if(v12345678912345a.eq.v12345678912345b) call pariserror('no long variables')
  end subroutine check_sanity_in_depth
!
  subroutine check_sanity()
    if(nx < 1) call err_no_out_dir("wrong nx")
    if(npx < 1) call err_no_out_dir("wrong npx")
    if(ny < 1) call err_no_out_dir("wrong ny")
    if(npy < 1) call err_no_out_dir("wrong npy")
    if(nz < 1) call err_no_out_dir("wrong nz")
    if(npz < 1) call err_no_out_dir("wrong npz")
    if(nx > 32767) call err_no_out_dir("nx too large")  ! why ? INTEGER*4 goes to 2**32-1 but NX*NY*NZ then too large ? 
  end subroutine check_sanity
!
  function EndProc(d)
    integer, intent(in) :: d
    integer EndProc
    if      (d==1) then
       EndProc = nPx-1
    else if (d==2)then
       EndProc = nPy-1
    else if (d==3) then
       EndProc = nPz-1
    else
       endproc=-1
       call pariserror("EndProc: wrong d.")
    endif
  end function EndProc
  function proclimit(d,sign)
    integer, intent(in) :: d,sign
    integer :: proclimit
    if(sign==-1) then
       proclimit=0
    else if(sign==1) then
       proclimit = EndProc(d)
    else
       proclimit=-1
       call pariserror("proclimit: wrong sign")
    endif
  end function proclimit
  function coordstart(d)
    integer, intent(in) :: d
    integer :: coordstart
    if      (d==1) then
       coordstart = is
    else if (d==2) then
       coordstart = js
    else if (d==3) then
       coordstart = ks
    else
       coordstart=-1
       call pariserror("coordstart: wrong d.")
    endif
  end function coordstart
  function coordend(d)
    integer, intent(in) :: d
    integer :: coordend
    if      (d==1) then
       coordend = ie
    else if (d==2) then
       coordend = je
    else if (d==3) then
       coordend = ke
    else
       coordend=-1
       call pariserror("coordend: wrong d.")
    endif
  end function coordend
  function coordlimit(d,sign)
    integer, intent(in) :: d,sign
    integer :: coordlimit
    if(sign==1) then
       coordlimit = coordend(d)
    else if(sign==-1) then
       coordlimit = coordstart(d)
    else
       coordlimit=-1
       call pariserror("coordlimit: wrong sign")
    endif
  end function coordlimit
  
  function slope_lim (val1,val2,val3,AdvScheme)
  implicit none
  real(8), external :: minabs, maxabs
  real (8) :: slope_lim,val1,val2,val3, alpha1, alpha2,alpha
  character(20), intent(in) :: AdvScheme 
  if ( AdvScheme == 'Superbee' ) then
     alpha1 = minabs(val3-val2,2*(val2-val1))
     alpha2 = minabs(2*(val3-val2),val2-val1)
     slope_lim = maxabs(alpha1, alpha2)
  else if ( AdvScheme == 'ENO' ) then 
     slope_lim = minabs((val3-val2),(val2-val1))
  else if ( AdvScheme == 'WENO' ) then 
     alpha1 = 1.0d0/((val2-val1)**2.d0 + 1.d-16)
     alpha2 = 1.0d0/((val3-val2)**2.d0 + 1.d-16) 
     alpha  = alpha1 + alpha2
     slope_lim = (alpha1*(val2-val1)+alpha2*(val3-val2))/alpha 
  endif
  end function slope_lim 
  
  function interpole3 (val1,val2,val3,AdvScheme,delta_x)
  implicit none
  real (8) :: interpole3, val1, val2, val3, delta_x
  character(20), intent(in) :: AdvScheme 
  
  interpole3 = val2 + slope_lim(val1,val2,val3,AdvScheme)*delta_x 
  
  end function interpole3 
  
  
  
end module module_grid
!=================================================================================================
! module_timer: 
!-------------------------------------------------------------------------------------------------
module module_timer
  implicit none
  save
  integer, parameter :: components=14, steps=1000
  real(8) :: times(components), percentage(components), tmp_time, this_time
  real(8) :: start_time, end_time=0.d0
  real(8) :: start_loop, end_loop
  integer :: ierr2
  character(25) :: timer_component(components)
  logical :: timer_initialized = .false.
  logical :: sizer_initialized = .false.
  real(8) :: alloc_size=0
contains
  subroutine add_2_my_sizer_2(nvalues,nbytes)
    use module_grid
    implicit none
    integer :: nvalues,nbytes
    alloc_size=alloc_size + dble(nvalues)*dble(nbytes)*1d-9
  end subroutine add_2_my_sizer_2
  subroutine remove_from_my_sizer_2(nvalues,nbytes)
    use module_grid
    implicit none
    integer :: nvalues,nbytes
    alloc_size=alloc_size - dble(nvalues)*dble(nbytes)*1d-9
  end subroutine remove_from_my_sizer_2

  subroutine add_2_my_sizer(n,nbytes)
    use module_grid
    implicit none
    integer n,nbytes
!    if(.not.sizer_initialized) call initialize_sizer()
! Alloc size in Gb
    alloc_size=alloc_size + dble(n*(nx/npx)*(ny/npy)*(nz/npz))*nbytes*1d-9
  end subroutine add_2_my_sizer
!  subroutine initialize_sizer
!    implicit none
!    alloc_size
!  end subroutine initialize_sizer
  subroutine initialize_timer
    use module_grid
    implicit none
    include 'mpif.h'
    !                     1234567890123456789012345
    timer_component(1) = 'Advection communications'
    timer_component(2) = 'Miscellaneous'
    timer_component(3) = 'Viscous terms'
    timer_component(4) = 'VOF advection'
    timer_component(5) = 'Heights'
    !                     1234567890123456789012345
    timer_component(6) = 'Heights communications'
    timer_component(7) = 'Curvature'
    timer_component(8) = 'Surface tension'
    timer_component(9) = 'Advection'
    timer_component(10) = 'Poisson for pressure'
    timer_component(11) = 'Output to file'
    timer_component(12) = 'Lag particles advection'
    timer_component(13) = 'Front Tracking'
    timer_component(14) = 'Tag, convert, drop stats'
    times=0d0
    this_time =  MPI_WTIME(ierr2)
    tmp_time = this_time
    start_loop = this_time
    timer_initialized=.true.
!    if(rank>0) return
!    open(unit=122,file='timer_stats')
!    write(122,'("Timer initialised at time",es16.2e2)') this_time - start_time
!    close(122)
  end subroutine initialize_timer
  subroutine my_timer(n)
    use module_grid
    implicit none
    include 'mpif.h'
    integer, intent(in) :: n
    real(8) :: elapsed_time
    if(rank>0) return
    if(n>components) call pariserror("n>components")
    if(.not.timer_initialized) call initialize_timer()
    this_time = MPI_WTIME(ierr2)
    elapsed_time =  this_time - tmp_time
    tmp_time = this_time
    times(n) = times(n) + elapsed_time
  end subroutine my_timer
  subroutine wrap_up_timer(itimestep,iTimeStepRestart)
    use module_grid
    implicit none
    include 'mpif.h'
    integer :: n,ierr2,itimestep,iTimeStepRestart
    real(8) :: totaltime, ZZ
    if(rank>0) return
    end_loop = MPI_WTIME(ierr2)
    ZZ = nx*ny*nz/npx/npy/npz*(itimestep-iTimeStepRestart)/(end_loop-start_loop)
    totaltime=0d0
    do n=1,components
       totaltime = totaltime + times(n)
    end do
    percentage = 1d2*times/totaltime
    open(unit=123,file='aggregate_timer_stats')
    do n=1,components
       write(123,'((A),T26," ",f5.1," %")') TRIM(timer_component(n)),percentage(n)
    enddo
    write(123,*) "   "
    write(123,'("Overall speed Z/np",T24,es16.5e2)') ZZ
    write(123,*) "   "
    write(123,'("Allocated/proc ",T24,es16.5e2,"Gbytes")') alloc_size
    write(123,*) "   "
    write(123,'("Current and initial time steps:",2(I7,1X))') itimestep,iTimeStepRestart
    write(123,*) "   "
    write(123,'("Current and initial time:",2(E15.8,1X))') end_loop,start_loop
    close(123)
  end subroutine wrap_up_timer
end module module_timer
!=================================================================================================
!=================================================================================================
! module_flow: Contains definition of variables for the flow solver.
!-------------------------------------------------------------------------------------------------
module module_flow
  implicit none
  save
  real(8), dimension(:,:,:), allocatable :: u, v, w, uold, vold, wold, fx, fy, fz, color
  real(8), dimension(:,:,:), allocatable :: momentum
  real(8), dimension(:,:,:), allocatable :: p, rho, rhoo, muold, mu, dIdx, dIdy, dIdz
  real(8), dimension(:,:,:), allocatable :: umask,vmask,wmask
  real(8), dimension(:,:,:), allocatable :: du,dv,dw,drho,du_c,dv_c,dw_c
  real(8), dimension(:,:), allocatable :: averages,oldaverages, allaverages
  logical, allocatable, dimension(:,:,:) :: mask

  real(8) :: gx, gy, gz, mu1, mu2, r_avg, dt, dtFlag, rho_ave, p_ave, vdt
  real(8) :: max_velocity, maxTime, Time, EndTime, MaxDt, CFL, mystats(100), stats(100)
  logical :: ZeroReynolds,DoVOF, DoFront, Implicit, hypre, GetPropertiesFromFront
  logical :: dosolids = .false.
  real(8) :: rho1, rho2, s
  real(8) :: U_init, VolumeSource, cflmax_allowed
  real(8) :: dpdx, dpdy, dpdz, W_ave  !pressure gradients in case of pressure driven channel flow
  real(8) :: dpdx_stat, dpdy_stat, dpdz_stat
  real(8) :: beta, MaxError,ErrorScaleHYPRE
  integer :: ResNormOrderPressure
  integer :: HYPRESolverType
  integer :: maxit, it, itime_scheme, BuoyancyCase, drive
  integer :: sbx, sby, Nstep
  integer :: maxStep, itmax, iTimeStep, iTimeStepRestart
  integer :: nstatarray
  character(20) :: AdvectionScheme
  
  integer :: num_probes,num_probes_cvof
  integer, parameter :: max_num_probes = 20 
  real(8) :: dat_probe(max_num_probes,5)  !u,v,z,c,p
  integer :: ijk_probe(max_num_probes,3)  !i,j,k
  real(8) :: dat_probe_cvof(max_num_probes)  
  integer :: ijk_probe_cvof(max_num_probes,3)  !i,j,k

  logical :: calcTurbStats_initialized = .false.
  logical :: DoTurbStats
  integer :: iSumTurbStats,nStepOutputTurbStats,TurbStatsOrder,num_turb_vars
  real(8) :: timeStartTurbStats
  real(8), dimension(:,:,:), allocatable :: turb_vars !i,j,ivar
  

end module module_flow
!=================================================================================================
!=================================================================================================
! Temporary variables
!-------------------------------------------------------------------------------------------------
module module_tmpvar
  real(8), dimension(:,:,:,:), allocatable, target :: work, A
  real(8), dimension(:,:,:), allocatable, target :: tmp
  real(8) :: tcpu(100),t0
end module module_tmpvar

!=================================================================================================
!=================================================================================================
! module_2phase: Contains variables for two-phase flow
!-------------------------------------------------------------------------------------------------
module module_2phase
  real(8), dimension( : ), allocatable :: rad, xc, yc, zc, vol
  real(8) :: excentricity(3)
  real(8) :: ugas_inject,uliq_inject
  real(8) :: blayer_gas_inject, tdelay_gas_inject 
  real(8) :: radius_gas_inject, radius_liq_inject, radius_gap_liqgas
  real(8) :: jetcenter_yc2yLength, jetcenter_zc2zLength 
  real(8) :: jetcenter_yc,         jetcenter_zc
  real(8) :: sigma
  integer :: NumBubble
  real(8) :: NozzleThick2Cell,NozzleLength
end module module_2phase
!=================================================================================================
!=================================================================================================
! module_freesurface: Contains variables for the free surface interface condition
!-------------------------------------------------------------------------------------------------
module module_freesurface
  real(8), dimension(:,:,:), allocatable :: x_mod, y_mod, z_mod, p_ext
  real(8), dimension(:,:,:), allocatable :: P_gx, P_gy, P_gz
  real(8), dimension(:,:,:), allocatable :: v_source
  integer, dimension(:,:,:), allocatable :: u_cmask,v_cmask,w_cmask,pcmask
  integer, dimension(1:4) :: NOUT_VTK
  real(8) :: P_ref, gamma, R_ref, V_0, P_inf !eq pressure and polytropic gas exponent
  real(8) :: R_RK, dR_RK, ddR_RK
  integer :: X_level, solver_flag=0
  logical :: FreeSurface, debug=.false., initialize_fs = .false.
  logical :: RP_test
  logical :: imploding
  logical, dimension(1:4) :: VTK_OUT, vtk_open
  character(len=10) :: visit_file(1:4) = (/ "divergence", "curvature ", "deltaPfs  ", "vol_Source" /)
  character(len=3) :: file_short(1:4) = (/ "DIV", "KAP", "Pfs", "S_v" /)
  real(8) :: v_total
end module module_freesurface
!=================================================================================================
!=================================================================================================
! module_IO: Contains input/output variables and procedures
!-------------------------------------------------------------------------------------------------
module module_IO
  implicit none
  save
  integer :: padding
  integer :: opened=0, opened_p=0
  integer :: nout, out, output_format, nbackup, nstats, termout, nfile
  integer :: nsteps_probe
  character(len=20) :: out_path, x_file, y_file, z_file
  logical :: read_x, read_y, read_z, restart, ICOut, restartFront, restartAverages, zip_data
  logical :: out_mom, output_fields(5)
  real(8) :: tout,timeLastOutput
contains
    !=================================================================================================
  subroutine write_vec_gnuplot(u,v,cvof,p,iout,DoVOF)
    use module_grid
    use module_tmpvar
    implicit none
    include 'mpif.h'
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: u, v, cvof, p
    integer, intent(in) :: iout
    integer :: i,j,k,ierr,indices(3)
    real(8) :: norm=0.d0, coeff, vmax
    logical :: DoVOF
    intrinsic dsqrt
    !
    ! writing u,v
    !
    OPEN(UNIT=89,FILE=TRIM(out_path)//'/UV-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout,padding))//'.txt')
    norm=0.d0
    k=(ks+ke)/2
    vmax = maxval(sqrt(u(is:ie,js:je,ks:ke)**2 + v(is:ie,js:je,ks:ke)**2))
    if(vmax.gt.1D50) then 
       indices = maxloc(sqrt(u(is:ie,js:je,ks:ke)**2 + v(is:ie,js:je,ks:ke)**2))
       if(rank==0) print*,"velocity diverged at",indices
       call pariserror("vmax diverged")
    endif
    if(vmax.ne.vmax) call pariserror("invalid vmax")
    call MPI_ALLREDUCE(vmax, norm, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_Cart, ierr)  
    coeff = 0.8d0/(norm + TINY_DOUBLE)
    do i=is,ie; do j=js,je
       write(89,310) x(i),y(j),coeff*dx(i)*0.5d0*(u(i,j,k)+u(i-1,j,k)), &
                  coeff*dy(j)*0.5d0*(v(i,j,k)+v(i,j-1,k))
    enddo; enddo
    close(unit=89)
    !
    if(DoVOF) then
       !
       ! writing cvof
       !
       OPEN(UNIT=89,FILE=TRIM(out_path)//'/CVoF-'//TRIM(int2text(rank,padding))//'-'// &
            TRIM(int2text(iout,padding))//'.txt')
       k=(ks+ke)/2
       do i=is,ie
          do j=js,je
             write(89,310) cvof(i,j,k)
          enddo
          WRITE(89,*) " "
       enddo
       close(unit=89)
    endif
    !
    ! writing p
    !
    OPEN(UNIT=89,FILE=TRIM(out_path)//'/P-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout,padding))//'.txt')
    k=(ks+ke)/2

    do j=js,je
       do i=is,ie
          write(89,310) x(i),y(j),p(i,j,k)
!          write(89,310) p(i,j,k)
       enddo 
       WRITE(89,*) " "
    enddo
    close(unit=89)
310 format(4e14.5)
  end subroutine  write_vec_gnuplot
  !=================================================================================================
  subroutine write_mom_gnuplot(du,dv,dw,iout)
    use module_grid
    implicit none
    include 'mpif.h'
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: du, dv, dw
    integer, intent(in) :: iout
    integer :: i,j,k,ierr
    
    OPEN(UNIT=20,FILE=TRIM(out_path)//'/Mag_accel-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout,padding))//'.txt')
    k=(ks+ke)/2
    do i=is,ie; do j=js,je
       write(20,13) x(i),y(j),sqrt(du(i,j,k)**2d0+dv(i,j,k)**2d0+dw(i,j,k)**2d0)
    enddo; enddo
    close(unit=20)
  
13 format(3e14.5)
  end subroutine  write_mom_gnuplot
  !=================================================================================================
  subroutine output_droplet(w,time)
    use module_grid
    implicit none
    include 'mpif.h'
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: w
    real(8), intent(in) :: time
    integer :: i,j,k
    i=nx/2
    j=ny/2
    k=3*nz/4
    if(test_point_in(i,j,k)) then
       OPEN(UNIT=89,FILE=TRIM(out_path)//'/droplet-test-vel.txt',status='unknown', &
            action='write',position='append')
       write(89,310) time, w(i,j,k)
       close(unit=89)
    endif
310 format(2e14.5)
  end subroutine  output_droplet
!=================================================================================================
  SUBROUTINE append_visit_file(rootname)
    use module_flow
    use module_grid
    implicit none
    character(*) :: rootname
    integer prank
    if(rank.ne.0) call pariserror("rank.ne.0 in append")
    
    if(opened==0) then
       OPEN(UNIT=90,FILE='velocity.visit')
       write(90,10) nPdomain
10     format('!NBLOCKS ',I4)
       opened=1
    else
    	OPEN(UNIT=90,FILE='velocity.visit',position='append')
    endif
    
    do prank=0,NpDomain-1
       write(90,11) rootname//TRIM(int2text(prank,padding))//'.vtk'
11     format(A)
    enddo
  end subroutine  append_visit_file
  
  subroutine close_visit_file()
    close(90)
  end subroutine close_visit_file
!=================================================================================================
! function int2text
!   Returns 'number' as a string with length of 'length'
!   called in:    function output
!-------------------------------------------------------------------------------------------------
function int2text(number,length)
  integer :: number, length, i
  character(len=length) :: int2text
  character, dimension(0:9) :: num = (/'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/)
  if(number>=10**length)print*, 'Warning: int2text: large input number. Increase "length".'
  do i=1,length
    int2text(length+1-i:length+1-i) = num(mod(number/(10**(i-1)),10))
  enddo
end function
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
subroutine backup_write
  use module_flow
  use module_grid
  implicit none
  integer ::i,j,k
  character(len=100) :: filename
  filename = trim(out_path)//'/backup_'//int2text(rank,padding)
  !call system('touch '//trim(filename)//'; mv '//trim(filename)//' '//trim(filename)//'.old')
  OPEN(UNIT=7,FILE=trim(filename),status='REPLACE')
  write(7,1100)time,itimestep,is,ie,js,je,ks,ke
  do k=ks,ke; do j=js,je; do i=is,ie
    write(7,1200) u(i,j,k), v(i,j,k), w(i,j,k), p(i,j,k), color(i,j,k)
  enddo; enddo; enddo
  do i=1,10; do j=js,je
    write(7,1200) averages(i,j)
  enddo; enddo
  close(7)
  if(rank==0)print*,'Backup written at t=',time
  1100 FORMAT(es17.8e3,7I10)
  1200 FORMAT(5es17.8e3)
end subroutine backup_write
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
subroutine backup_read
  use module_flow
  use module_grid
  
  implicit none
  integer ::i,j,k,i1,i2,j1,j2,k1,k2,ierr
  OPEN(UNIT=7,FILE=trim(out_path)//'/backup_'//int2text(rank,padding),status='old',action='read')
  read(7,*)time,itimestep,i1,i2,j1,j2,k1,k2
  if(i1/=is .or. i2/=ie .or. j1/=js .or. j2/=je .or. k1/=ks .or. k2/=ke) &
    call pariserror("Error: backup_read")
  do k=ks,ke; do j=js,je; do i=is,ie
    read(7,*) u(i,j,k), v(i,j,k), w(i,j,k), p(i,j,k), color(i,j,k)
  enddo; enddo; enddo
  if(restartAverages)then
    do i=1,20; do j=js,je
      read(7,*,iostat=ierr) averages(i,j)
      if(ierr/=0)return !no average data, return
    enddo; enddo
    !oldaverages=averages
  endif
  close(7)
!  1100 FORMAT(es25.16e3,7I10)
!  1200 FORMAT(5es25.16e3)
end subroutine backup_read
!=================================================================================================
!=================================================================================================
! subroutine output
!   Writes the solution in the 'out_path' directory in the tecplot or vtk format
!   called in:    program main
!-------------------------------------------------------------------------------------------------
subroutine output(nf,i1,i2,j1,j2,k1,k2)
  implicit none
  integer :: nf,i1,i2,j1,j2,k1,k2
  if(output_format==1) call output1(nf,i1,i2,j1,j2,k1,k2)
  if(output_format==2) call output2(nf,i1,i2,j1,j2,k1,k2)
  if(output_format==3) call output3(nf,i1,i2,j1,j2,k1,k2)
end subroutine output
!-------------------------------------------------------------------------------------------------
subroutine output1(nf,i1,i2,j1,j2,k1,k2)
  use module_flow
  use module_grid
  
!  use module_tmpvar
  !use IO_mod
  implicit none
  integer :: nf,i1,i2,j1,j2,k1,k2,i,j,k
  logical, save :: first_time=.true.
  
  if(first_time)then
    first_time = .false.
    OPEN(UNIT=7,FILE=trim(out_path)//'/grid_'//int2text(rank,3)//'.dat')
    write(7,*) 'FILETYPE = GRID'
    write(7,*) 'variables = "x","y","z"'
    write(7,2100) i2-i1+1, j2-j1+1, k2-k1+1
    do k=k1,k2; do j=j1,j2; do i=i1,i2;
      write(7,1200) x(i),y(j),z(k)
    enddo; enddo; enddo
    close(7)
  endif

  OPEN(UNIT=7,FILE=trim(out_path)//'/plot'//int2text(nf,3)//'_'//int2text(rank,3)//'.dat')
  write(7,*) 'FILETYPE = SOLUTION'
  write(7,*) 'variables = "u","v","w","p","color"'
  write(7,1100) time, i2-i1+1, j2-j1+1, k2-k1+1
  do k=k1,k2; do j=j1,j2; do i=i1,i2;
    write(7,1200) 0.5d0*(u(i,j,k)+u(i-1,j,k)), &
                  0.5d0*(v(i,j,k)+v(i,j-1,k)), &
                  0.5d0*(w(i,j,k)+w(i,j,k-1)), &
                  p(i,j,k) , color(i,j,k)
  enddo; enddo; enddo
  close(7)
  if(rank==0)then
    OPEN(UNIT=7,FILE=trim(out_path)//'/averages.dat',position='append')
    write(7,1110) time
    do j=Ng,Ng+Ny+1
      write(7,1200) y(j),real(allaverages(:,j)-oldaverages(:,j))
    enddo
    close(7)
    oldaverages=allaverages
  endif
!  1000 FORMAT('variables = "x","y","z","u","v","w","p","color"')
  1100 FORMAT('ZONE solutiontime=', 1PG15.7e2, ', I=',I4, 2X, ', J=',I4, 2X, ', K=',I4, 2X)
  1110 FORMAT('ZONE solutiontime=', 1PG15.7e2)
  1200 FORMAT(21es14.6e2)
  2100 FORMAT('ZONE I=',I4, 2X, ', J=',I4, 2X, ', K=',I4, 2X)
end subroutine output1
!-------------------------------------------------------------------------------------------------
subroutine output2(nf,i1,i2,j1,j2,k1,k2)
  use module_flow
  use module_grid

  !use IO_mod
  implicit none
  integer ::nf,i1,i2,j1,j2,k1,k2,i,j,k, itype=5
  !  logical, save :: first_time=.true.
  character(len=30) :: rootname,filename
  rootname=TRIM(out_path)//'/VTK/plot'//TRIM(int2text(nf,padding))//'-'

  if(rank==0) call append_visit_file(TRIM(rootname))

  OPEN(UNIT=8,FILE=TRIM(rootname)//TRIM(int2text(rank,padding))//'.vtk')
  write(8,10)
  write(8,11) time
  write(8,12)
  write(8,13)
  write(8,14)i2-i1+1,j2-j1+1,k2-k1+1
  write(8,15)(i2-i1+1)*(j2-j1+1)*(k2-k1+1)

  do k=k1,k2; do j=j1,j2; do i=i1,i2;
     write(8,320) x(i),y(j),z(k)
  enddo; enddo; enddo

  write(8,19)(i2-i1+1)*(j2-j1+1)*(k2-k1+1)
  if (itype .le. 4)then
     write(8,17)'density'
     write(8,18)
  else
     write(8,20)
  endif
20 format('VECTORS uv double')
  do k=k1,k2; do j=j1,j2; do i=i1,i2;
     if (itype .eq. 1)write(8,210)rho(i,j,k)
     if (itype .eq. 5)write(8,310)0.5*(u(i,j,k)+u(i-1,j,k)), &
          0.5*(v(i,j,k)+v(i,j-1,k)),0.5*(w(i,j,k)+w(i,j,k-1))
  enddo; enddo; enddo
310 format(e14.5,e14.5,e14.5)

  close(8)

  ! TEMPORARY
  if ( zip_data ) then 
     filename = TRIM(rootname)//TRIM(int2text(rank,padding))//'.vtk'
     call system('gzip '//trim(filename))
  end if ! zip_data
  ! END TEMPORARY

10 format('# vtk DataFile Version 2.0')
11 format('grid, time ',F16.8)
12 format('ASCII')
13 format('DATASET STRUCTURED_GRID')
14 format('DIMENSIONS ',I5,I5,I5)
15 format('POINTS ',I17,' double' )
  !16  format('SPACING ',F16.8,F16.8,F16.8)
19 format('POINT_DATA ',I17)
17 format('SCALARS ',A20,' double 1')
18 format('LOOKUP_TABLE default')

210 format(e14.5)
320 format(e14.5,e14.5,e14.5)
end subroutine output2
!-------------------------------------------------------------------------------------------------
subroutine output3(nf,i1,i2,j1,j2,k1,k2)
  use module_flow
  use module_grid
  !use IO_mod
  implicit none
  integer ::nf,i1,i2,j1,j2,k1,k2,i,j,k, itype=5
!  logical, save :: first_time=.true.
  character(len=30) :: rootname,filename
  rootname=TRIM(out_path)//'/VTK/plot'//TRIM(int2text(nf,padding))//'-'

  if(rank==0) call append_visit_file(TRIM(rootname))

  OPEN(UNIT=8,FILE=TRIM(rootname)//TRIM(int2text(rank,padding))//'.vtk')
    write(8,10)
    write(8,11)time
    write(8,12)
    write(8,13)
    write(8,14)i2-i1+1,j2-j1+1,k2-k1+1
    write(8,15) x(i1),y(j1),z(k1)
    write(8,16) x(i1+1)-x(i1),y(j1+1)-y(j1),z(k1+1)-z(k1)
10  format('# vtk DataFile Version 2.0')
11  format('grid, time ',F16.8)
12  format('ASCII')
13  format('DATASET STRUCTURED_POINTS')
14  format('DIMENSIONS ',I5,I5,I5)
15  format('ORIGIN ',F16.8,F16.8,F16.8)
16  format('SPACING ',F16.8,F16.8,F16.8)

    write(8,19)(i2-i1+1)*(j2-j1+1)*(k2-k1+1)
    if (itype .le. 4)then
      write(8,17)'density'
      write(8,18)
    else
      write(8,20)
    endif
19  format('POINT_DATA ',I17)
17  format('SCALARS ',A20,' double 1')
20  format('VECTORS uv double')
18  format('LOOKUP_TABLE default')

    do k=k1,k2; do j=j1,j2; do i=i1,i2;
      if (itype .eq. 1)write(8,210)rho(i,j,k)
      if (itype .eq. 5)write(8,310)0.5*(u(i,j,k)+u(i-1,j,k)), &
         0.5*(v(i,j,k)+v(i,j-1,k)),0.5*(w(i,j,k)+w(i,j,k-1))
    enddo; enddo; enddo
210 format(e14.5)
310 format(e14.5,e14.5,e14.5)

    close(8)

! TEMPORARY
    if ( zip_data ) then 
      filename = TRIM(rootname)//TRIM(int2text(rank,padding))//'.vtk'
      call system('gzip '//trim(filename))
    end if ! zip_data
! END TEMPORARY 
end subroutine output3
!-------------------------------------------------------------------------------------------------
end module module_IO
!=================================================================================================
!=================================================================================================
! module module_BC: Sets the boundary conditions of the problem.
!-------------------------------------------------------------------------------------------------
module module_BC
  use module_grid
  
  implicit none
  integer :: bdry_cond(6), inject_type
  logical :: bdry_read=.false.
  ! bdry_cond(i) = is the type if boundary condition in i'th direction
  ! explicits the boundary condition codes
  !                                   12345678    12345678    12345678    12345678    12345678
   character(len=8) :: expl(0:6) = (/ "wall    ", "periodic", "shear   ", "vel-in  ", "vel-out ", "pressure", "pres-out" /)
  real(8) :: WallVel(6,3), WallShear(6,3), BoundaryPressure(6)
  ! Tangential velocities on the surfaces of domain. First index represent the 
  ! side on which the velocity in the direction of the second index is specified.
  ! The sides are in this order: -x,+x,-y,+y,-z,+z.
  ! Example: WallVel(4,3) represent the W velocity on +y side of the domain.
  ! LM: The same convention is used for BoundaryPressure
  ! SZ: alternately may contain the velocity of the flow for inflow boundary conditions on x+
  logical :: check_setup=.true.
  logical :: OutVelSpecified
  logical :: LateralBdry(6) ! A LateralBdry is defined as parallel to the main stream  
  real(8) :: MaxFluxRatioPresBC  
  contains
!=================================================================================================
!=================================================================================================
! subroutine SetPressureBC: Sets the pressure boundary condition
!-------------------------------------------------------------------------------------------------
    subroutine SetPressureBC(umask,vmask,wmask)
    use module_grid
    implicit none
    include 'mpif.h'
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: umask,vmask,wmask

  ! for walls set the mask to zero  ! @@@ Aijk coefficients should be changed too. 
    if(bdry_cond(1)==0 .or. bdry_cond(1)==2) then
      if(coords(1)==0    ) umask(is-1,js-1:je+1,ks-1:ke+1)=0d0
    endif
    
    if(bdry_cond(4)==0 .or. bdry_cond(4)==2) then
      if(coords(1)==nPx-1) umask(ie,js-1:je+1,ks-1:ke+1)=0d0
    endif

    if(bdry_cond(2)==0 .or. bdry_cond(2)==2) then
      if(coords(2)==0    ) vmask(is-1:ie+1,js-1,ks-1:ke+1)=0d0
    endif 
    
    if(bdry_cond(5)==0 .or. bdry_cond(5)==2) then
      if(coords(2)==nPy-1) vmask(is-1:ie+1,je,ks-1:ke+1)=0d0
    endif 

    if(bdry_cond(3)==0 .or. bdry_cond(3)==2) then
      if(coords(3)==0    ) wmask(is-1:ie+1,js-1:je+1,ks-1)=0d0
    endif
    
    if(bdry_cond(6)==0 .or. bdry_cond(6)==2) then
      if(coords(3)==nPz-1) wmask(is-1:ie+1,js-1:je+1,ke)=0d0
    endif
  end subroutine SetPressureBC
!=================================================================================================
!=================================================================================================
! subroutine SetVelocityBC: Sets the velocity boundary condition
!-------------------------------------------------------------------------------------------------
  subroutine SetVelocityBC(u,v,w,umask,vmask,wmask,t,dt,AfterProjection)
    use module_grid
    use module_2phase
    use module_freesurface
    implicit none
    include 'mpif.h'
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: u, v, w
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: umask,vmask,wmask
    integer, intent (in) :: AfterProjection
    real(8) :: t,dt,fluxin,tfluxin,uaverage
    real(8) :: fluxout(6),tfluxout(6),tfluxout_all,fluxratio
    integer :: i,j,k,ierr
    ! Note: local and global divergence free cannot be perfectly satisfied 
    ! at the mean time for pressure BC (p=p0,du/dn=0), in BC=6, the global 
    ! divergence free is guaranteed by matching the inflow condiction. 
    ! Both local divergence and du/dn=0 are relatixed to achieve that. 
    ! Warning: Be CAUTIOUS when Pressure BC is used !!!  
    
    ! solid obstacles
    u = u*umask
    v = v*vmask
    w = w*wmask

    ! --------------------------------------------------------------------------------------------
    ! Inflow BC
    ! --------------------------------------------------------------------------------------------
    
    ! inflow boundary condition x- with injection
    fluxin=0
    if(bdry_cond(1)==3 .and. coords(1)==0    ) then
       do j=jmin,jmax
          do k=kmin,kmax
             u(is-1,j,k)=WallVel(1,1)*uinject(j,k,t)*umask(is-1,j,k)
             u(is-2,j,k)=WallVel(1,1)*uinject(j,k,t)*umask(is-2,j,k)
             v(is-1,j,k)=0d0
             w(is-1,j,k)=0d0
          enddo
       enddo
       do j=js,je
         do k=ks,ke
            fluxin = fluxin + u(is-1,j,k)
         enddo
       enddo
    endif
    call MPI_ALLREDUCE(fluxin, tfluxin, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr)
    uaverage=tfluxin/(ny*nz)

    ! inflow boundary condition y-
    if(bdry_cond(2)==3 .and. coords(2)==0    ) then
       do i=imin,imax
          do k=kmin,kmax
             v(i,js-1,k)=WallVel(3,2)
             v(i,js-2,k)=WallVel(3,2)
             u(i,js-1,k)=0d0
             w(i,js-1,k)=0d0
          enddo
       enddo
    endif
    ! inflow on z-
    if(bdry_cond(3)==3 .and. coords(3)==0   ) then
       do i=imin,imax
          do j=jmin,jmax
             w(i,j,ks-1)= WallVel(5,3)
             w(i,j,ks-2)= WallVel(5,3)
             u(i,j,ks-1)=0d0
             v(i,j,ks-1)=0d0
          enddo
       enddo
    endif
    ! inflow on x+
    if(bdry_cond(4)==3 .and. coords(1)==nPx-1   ) then
       do j=jmin,jmax
          do k=kmin,kmax
             u(ie,j,k)= WallVel(2,1)
             u(ie+1,j,k)= WallVel(2,1)
             v(ie+1,j,k)=0d0
             w(ie+1,j,k)=0d0
          enddo
       enddo
    endif
     ! inflow on y+
    if(bdry_cond(5)==3 .and. coords(2)==nPy-1   ) then
       do i=imin,imax
          do k=kmin,kmax
             v(i,je,k)= WallVel(4,2)
             v(i,je+1,k)= WallVel(4,2)
             u(i,je+1,k)=0d0
             w(i,je+1,k)=0d0
          enddo
       enddo
    endif
    ! inflow on z+
    if(bdry_cond(6)==3 .and. coords(3)==nPz-1   ) then
       do i=imin,imax
          do j=jmin,jmax
             w(i,j,ke)= WallVel(6,3)
             w(i,j,ke+1)= WallVel(6,3)
             u(i,j,ke+1)=0d0
             v(i,j,ke+1)=0d0
          enddo
       enddo
    endif    

    ! --------------------------------------------------------------------------------------------
    ! Wall BC with velocity specified
    ! --------------------------------------------------------------------------------------------
    ! wall boundary condition x-
    if(bdry_cond(1)==0 .and. coords(1)==0    ) then
        u(is-1,:,:)=0d0
        u(is-2,:,:)=-u(is,:,:)
        v(is-1,:,:)=2*WallVel(1,2)-v(is,:,:)
        w(is-1,:,:)=2*WallVel(1,3)-w(is,:,:)
    endif
    
    ! wall boundary condition x+
    if(bdry_cond(4)==0 .and. coords(1)==nPx-1) then  
        u(ie  ,:,:)=0d0
        u(ie+1,:,:)=-u(ie-1,:,:)                  ! second order extrapolation
        v(ie+1,:,:)=2*WallVel(2,2)-v(ie,:,:)      ! second order extrapolation
        w(ie+1,:,:)=2*WallVel(2,3)-w(ie,:,:)      ! second order extrapolation
    endif

    ! wall boundary condition y-
    if(bdry_cond(2)==0 .and. coords(2)==0    ) then
        v(:,js-1,:)=0d0
        v(:,js-2,:)=-v(:,js,:)
        u(:,js-1,:)=2*WallVel(3,1)-u(:,js,:)
        w(:,js-1,:)=2*WallVel(3,3)-w(:,js,:)
    endif
    ! wall boundary condition y+
    if(bdry_cond(5)==0 .and. coords(2)==nPy-1) then
        v(:,je  ,:)=0d0
        v(:,je+1,:)=-v(:,je-1,:)
        u(:,je+1,:)=2*WallVel(4,1)-u(:,je,:)
        w(:,je+1,:)=2*WallVel(4,3)-w(:,je,:)
    endif
    ! wall boundary condition z-
    if(bdry_cond(3)==0 .and. coords(3)==0    ) then
        w(:,:,ks-1)=0d0
        w(:,:,ks-2)=-w(:,:,ks)
        u(:,:,ks-1)=2*WallVel(5,1)-u(:,:,ks)
        v(:,:,ks-1)=2*WallVel(5,2)-v(:,:,ks)
    endif
    ! wall boundary condition z+
    if(bdry_cond(6)==0 .and. coords(3)==nPz-1) then
        w(:,:,ke  )=0d0
        w(:,:,ke+1)=-w(:,:,ke-1)
        u(:,:,ke+1)=2*WallVel(6,1)-u(:,:,ke)
        v(:,:,ke+1)=2*WallVel(6,2)-v(:,:,ke)
    endif
    
    ! --------------------------------------------------------------------------------------------
    ! Wall BC with shear specified (Shear is by default zero then is eqv to slip-wall BC)
    ! --------------------------------------------------------------------------------------------
    ! wall boundary condition: shear
    if(bdry_cond(1)==2 .and. coords(1)==0    ) then
        u(is-1,:,:)=0d0
        u(is-2,:,:)=-u(is,:,:)
        v(is-1,:,:) = -dxh(is-1)*WallShear(1,2)+v(is,:,:)
        w(is-1,:,:) = -dxh(is-1)*WallShear(1,3)+w(is,:,:)
    endif
    if(bdry_cond(4)==2 .and. coords(1)==nPx-1) then
        u(ie  ,:,:)=0d0
        u(ie+1,:,:)=-u(ie-1,:,:)
        v(ie+1,:,:) = dxh(ie)*WallShear(2,2)+v(ie,:,:)
        w(ie+1,:,:) = dxh(ie)*WallShear(2,3)+w(ie,:,:)
    endif
    if(bdry_cond(2)==2 .and. coords(2)==0    ) then
        v(:,js-1,:)=0d0
        v(:,js-2,:)=-v(:,js,:)
        u(:,js-1,:) = -dyh(js-1)*WallShear(3,1)+u(:,js,:)
        w(:,js-1,:) = -dyh(js-1)*WallShear(3,3)+w(:,js,:)
    endif
    if(bdry_cond(5)==2 .and. coords(2)==nPy-1) then
        v(:,je  ,:)=0d0
        v(:,je+1,:)=-v(:,je-1,:)
        u(:,je+1,:) = dyh(je)*WallShear(4,1)+u(:,je,:)
        w(:,je+1,:) = dyh(je)*WallShear(4,3)+w(:,je,:)
    endif
    if(bdry_cond(3)==2 .and. coords(3)==0    ) then
        w(:,:,ks-1)=0d0
        w(:,:,ks-2)=-w(:,:,ks)
        u(:,:,ks-1) = -dzh(ks-1)*WallShear(5,1)+u(:,:,ks)
        v(:,:,ks-1) = -dzh(ks-1)*WallShear(5,2)+v(:,:,ks)
    endif
    if(bdry_cond(6)==2 .and. coords(3)==nPz-1) then
        w(:,:,ke  )=0d0
        w(:,:,ke+1)=-w(:,:,ke-1)
        u(:,:,ke+1) = dzh(ke)*WallShear(6,1)+u(:,:,ke)
        v(:,:,ke+1) = dzh(ke)*WallShear(6,2)+v(:,:,ke)
    endif

    ! --------------------------------------------------------------------------------------------
    ! Outflow BC with pressure specified (by default zero)
    ! Set zero normal velocity gradient 
    ! --------------------------------------------------------------------------------------------
    if (bdry_cond(1)==5 .and. coords(1)==0)then
       u(is-1,:,:)=u(is,:,:)
#ifndef OLD_BDRY_COND
       u(is-2,:,:)=u(is,:,:)
#endif
       v(is-1,:,:)=v(is,:,:)
       w(is-1,:,:)=w(is,:,:)
    endif
    
    if (bdry_cond(4)==5 .and. coords(1)==nPx-1)then
       u(ie,:,:)=u(ie-1,:,:)
#ifndef OLD_BDRY_COND
       u(ie+1,:,:)=u(ie-1,:,:)
#endif
       v(ie+1,:,:)=v(ie,:,:)
       w(ie+1,:,:)=w(ie,:,:)
    endif
    
    if (bdry_cond(2)==5 .and. coords(2)==0)then
       v(:,js-1,:)=v(:,js,:)
#ifndef OLD_BDRY_COND
       v(:,js-2,:)=v(:,js,:)
#endif
       u(:,js-1,:)=u(:,js,:)
       w(:,js-1,:)=w(:,js,:)
    endif
    
    if (bdry_cond(5)==5 .and. coords(2)==nPy-1)then
       v(:,je,:)=v(:,je-1,:)
#ifndef OLD_BDRY_COND
       v(:,je+1,:)=v(:,je-1,:)
#endif
       u(:,je+1,:)=u(:,je,:)
       w(:,je+1,:)=w(:,je,:)
    endif
    
    if (bdry_cond(3)==5 .and. coords(3)==0)then
       w(:,:,ks-1)=w(:,:,ks)
#ifndef OLD_BDRY_COND
       w(:,:,ks-2)=w(:,:,ks)
#endif
       u(:,:,ks-1)=u(:,:,ks)
       v(:,:,ks-1)=v(:,:,ks)
    endif
    
    if (bdry_cond(6)==5 .and. coords(3)==nPz-1)then
       w(:,:,ke)=w(:,:,ke-1)
#ifndef OLD_BDRY_COND
       w(:,:,ke+1)=w(:,:,ke-1)
#endif
       u(:,:,ke+1)=u(:,:,ke)
       v(:,:,ke+1)=v(:,:,ke)    
    endif 
    
    ! --------------------------------------------------------------------------------------------
    ! Outflow BC with velocity specified
    ! Pressure gradient is set to zero in Poisson_BC
    ! same velocity as opposing inflow. ! @generalize this !!
    ! --------------------------------------------------------------------------------------------
    if(bdry_cond(4)==4 .and. coords(1)==nPx-1) then
        if ( OutVelSpecified ) then 
           do j=js,je; do k=ks,ke
              u(ie,j,k)=uout(j,k,t)
           end do; end do
        else 
           u(ie  ,:,:)=uaverage
#ifndef OLD_BDRY_COND
           u(ie+1,:,:)=uaverage
#endif
        end if ! OutVelSpecified
        v(ie+1,:,:)=v(ie,:,:)
        w(ie+1,:,:)=w(ie,:,:)
    endif

    ! --------------------------------------------------------------------------------------------
    ! Convective BC  (du/dt+Un*du/dx =0) 
    ! --------------------------------------------------------------------------------------------
    !flux = 0.d0 
    !if(bdry_cond(4)==6 .and. coords(1)==nPx-1) then
    !   do j=js,je
    !     do k=ks,ke
    !        flux = flux + u(ie-1,j,k)
    !     enddo
    !   enddo
    !end if ! bdry_cond(4)==6

    !call MPI_ALLREDUCE(flux, tflux, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr)
    !uaverage=tflux/(ny*nz)

    !if(bdry_cond(4)==6 .and. coords(1)==nPx-1) then
    !   u(ie  ,:,:)=u(ie  ,:,:)-dt*uaverage*(u(ie-1,:,:)-u(ie-2,:,:))/dx(ie-1)
    !   v(ie+1,:,:)=v(ie+1,:,:)-dt*uaverage*(v(ie  ,:,:)-v(ie-1,:,:))/dx(ie-1)
    !   w(ie+1,:,:)=w(ie+1,:,:)-dt*uaverage*(w(ie  ,:,:)-w(ie-1,:,:))/dx(ie-1)
    !end if ! bdry_cond(4)==6
     
    ! --------------------------------------------------------------------------------------------
    ! New pressure BC  (du/dn=0 & p=p0) 
    ! --------------------------------------------------------------------------------------------
    ! Note: for pressure BC, vel-BC apply after projection 
    fluxout(:) = 0.d0 
    if      (bdry_cond(1)==6 .and. coords(1)==0     & 
            .and. AfterProjection == 1 .and. .not.LateralBdry(1)) then
       do j=js,je
         do k=ks,ke
            fluxout(1) = fluxout(1) + u(is  ,j,k)
         enddo
       enddo
    else if (bdry_cond(4)==6 .and. coords(1)==nPx-1 & 
            .and. AfterProjection == 1 .and. .not.LateralBdry(4)) then
       do j=js,je
         do k=ks,ke
            fluxout(4) = fluxout(4) + u(ie-1,j,k)
         enddo
       enddo
    else if (bdry_cond(2)==6 .and. coords(2)==0     & 
            .and. AfterProjection == 1 .and. .not.LateralBdry(2)) then
       do i=is,ie
         do k=ks,ke
            fluxout(2) = fluxout(2) + v(i,js  ,k)
         enddo
       enddo
    else if (bdry_cond(5)==6 .and. coords(2)==nPy-1 & 
            .and. AfterProjection == 1.and. .not.LateralBdry(5)) then
       do i=is,ie
         do k=ks,ke
            fluxout(5) = fluxout(5) + v(i,je-1,k)
         enddo
       enddo
    else if (bdry_cond(3)==6 .and. coords(3)==0     & 
            .and. AfterProjection == 1 .and. .not.LateralBdry(3)) then
       do i=is,ie
         do j=js,je
            fluxout(3) = fluxout(3) + w(i,j,ks  )
         enddo
       enddo
    else if (bdry_cond(6)==6 .and. coords(3)==nPz-1 & 
            .and. AfterProjection == 1 .and. .not.LateralBdry(6)) then
       do i=is,ie
         do j=js,je
            fluxout(6) = fluxout(6) + w(i,j,ke-1)
         enddo
       enddo
    endif
    call MPI_ALLREDUCE(fluxout, tfluxout, 6, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr)
    tfluxout_all = sum(tfluxout)
    fluxratio = min(ABS(tfluxin/(tfluxout_all+1.d-16)),MaxFluxRatioPresBC)  ! fluxratio is capped 
    if      (bdry_cond(1)==6 .and. coords(1)==0     & 
      .and. AfterProjection == 1) then
       u(is-1,:,:)=u(is,:,:)*fluxratio
       v(is-1,:,:)=v(is,:,:)
       w(is-1,:,:)=w(is,:,:)
    else if (bdry_cond(4)==6 .and. coords(1)==nPx-1 .and. AfterProjection == 1) then
       u(ie,:,:) = u(ie-1,:,:)*fluxratio 
       v(ie+1,:,:)=v(ie,:,:)
       w(ie+1,:,:)=w(ie,:,:)
    else if (bdry_cond(2)==6 .and. coords(2)==0     .and. AfterProjection == 1) then
       v(:,js-1,:)=v(:,js,:)*fluxratio
       u(:,js-1,:)=u(:,js,:)
       w(:,js-1,:)=w(:,js,:)
    else if (bdry_cond(5)==6 .and. coords(2)==nPy-1 .and. AfterProjection == 1) then
       v(:,je,:)=v(:,je-1,:)*fluxratio
       u(:,je+1,:)=u(:,je,:)
       w(:,je+1,:)=w(:,je,:)
    else if (bdry_cond(3)==6 .and. coords(3)==0     .and. AfterProjection == 1) then
       w(:,:,ks-1)=w(:,:,ks)*fluxratio
       u(:,:,ks-1)=u(:,:,ks)
       v(:,:,ks-1)=v(:,:,ks)
    else if (bdry_cond(6)==6 .and. coords(3)==nPz-1 .and. AfterProjection == 1) then
       w(:,:,ke)=w(:,:,ke-1)*fluxratio
       u(:,:,ke+1)=u(:,:,ke)
       v(:,:,ke+1)=v(:,:,ke)    
    end if ! bdry_cond()

  end subroutine SetVelocityBC

  !=================================================================================================
! subroutine SetMomentumBC: Sets the momentum boundary condition
!-------------------------------------------------------------------------------------------------
  subroutine SetMomentumBC(u,c,mom,d,mask,rho1,rho2,t)
    use module_grid
    use module_2phase
    
    implicit none
    include 'mpif.h'
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: u, c, mom
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: mask
    real(8), intent(in) :: rho1,rho2
    integer, intent(in) :: d
    real(8) :: t,flux,tflux,uaverage
    integer :: i,j,k,ierr
    ! solid obstacles
    u = u*mask
    mom = mom*mask

    ! wall boundary condition
    if(bdry_cond(1)==0 .and. coords(1)==0    ) then
        if (d.eq.1) then
            mom(is-1,:,:)=0d0
            mom(is-2,:,:)=-mom(is,:,:)
        else !y,z
            mom(is-1,:,:)=(2*WallVel(1,2)-u(is,:,:))*(rho2*c(is,:,:) + rho1*(1.d0 - c(is,:,:))) !CHECK!!
        endif
    endif
    ! inflow boundary condition x- with injection

    if(bdry_cond(1)==3 .and. coords(1)==0    ) then
        if (d.eq.1) then
        flux=0
        do j=jmin,jmax
          do k=kmin,kmax
             mom(is-1,j,k)=WallVel(1,1)*uinject(j,k,t)*(rho2*c(is-1,j,k) + rho1*(1.d0 - c(is-1,j,k)))
#ifndef OLD_BDRY_COND
             mom(is-2,j,k)=WallVel(1,1)*uinject(j,k,t)*(rho2*c(is-1,j,k) + rho1*(1.d0 - c(is-1,j,k)))
             if(j<=je.and.j>=js.and.k<=ke.and.k>=ks) flux=flux+WallVel(1,1)*uinject(j,k,t)
#endif
          enddo
        enddo
        else
            mom(is-1,:,:) = 0.d0
        endif
    endif
    call MPI_ALLREDUCE(flux, tflux, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr)
    uaverage=tflux/(ny*nz)


    ! inflow boundary condition y-
    if(bdry_cond(2)==3 .and. coords(2)==0    ) then
       if (d.eq.2) then
        do i=imin,imax
          do k=kmin,kmax
             mom(i,js-1,k)=WallVel(3,2)*(rho2*c(i,js-1,k) + rho1*(1.d0 - c(i,js-1,k)))
             mom(i,js-2,k)=WallVel(3,2)*(rho2*c(i,js-2,k) + rho1*(1.d0 - c(i,js-2,k)))
          enddo
        enddo
       else
           mom(:,js-1,:) = 0.d0
       endif
    endif
    ! inflow on z-
    if(bdry_cond(3)==3 .and. coords(3)==0   ) then
       if (d.eq.3) then
         do i=imin,imax
          do j=jmin,jmax
             mom(i,j,ks-1)= WallVel(5,3)*(rho2*c(i,j,ks-1) + rho1*(1.d0 - c(i,j,ks-1)))
             mom(i,j,ks-2)= WallVel(5,3)*(rho2*c(i,j,ks-2) + rho1*(1.d0 - c(i,j,ks-2)))
          enddo
         enddo
     else
         mom(:,:,ks-1) = 0.d0
     endif
    endif
    ! inflow on x+
    if(bdry_cond(4)==3 .and. coords(1)==nPx-1   ) then
        if (d.eq.1) then
         do j=jmin,jmax
          do k=kmin,kmax
             mom(ie,j,k)  = WallVel(2,1)*(rho2*c(ie,j,k) + rho1*(1.d0 - c(ie,j,k)))
             mom(ie+1,j,k)= WallVel(2,1)*(rho2*c(ie+1,j,k) + rho1*(1.d0 - c(ie+1,j,k)))
          enddo
         enddo
       else
           mom(ie,:,:) = 0.d0
       endif
    endif
     ! inflow on y+
    if(bdry_cond(5)==3 .and. coords(2)==nPy-1   ) then
        if (d.eq.2) then
           do i=imin,imax
              do k=kmin,kmax
                mom(i,je,k)  = WallVel(4,2)*(rho2*c(i,je,k) + rho1*(1.d0 - c(i,je,k)))
                mom(i,je+1,k)= WallVel(4,2)*(rho2*c(i,je+1,k) + rho1*(1.d0 - c(i,je+1,k)))
              enddo
            enddo
        else
            mom(:,je,:) = 0.d0
        endif
    endif
    ! inflow on z+
    if(bdry_cond(6)==3 .and. coords(3)==nPz-1   ) then
       if (d.eq.3) then
        do i=imin,imax
          do j=jmin,jmax
             mom(i,j,ke)  = WallVel(6,3)*(rho2*c(i,j,ke) + rho1*(1.d0 - c(i,j,ke)))
             mom(i,j,ke+1)= WallVel(6,3)*(rho2*c(i,j,ke+1) + rho1*(1.d0 - c(i,j,ke+1)))
          enddo
        enddo
       else
           mom(:,:,ke) = 0.d0
       endif
    endif    

    if(bdry_cond(4)==0 .and. coords(1)==nPx-1) then 
        if (d.eq.1) then
            mom(ie  ,:,:)=0d0
            mom(ie+1,:,:)=-mom(ie-1,:,:)
        else
            mom(ie+1,:,:)=(2*WallVel(2,2)-u(ie,:,:))**(rho2*c(ie,:,:) + rho1*(1.d0 - c(ie,:,:)))
        endif
    endif
    
    ! outflow boundary condition
    if(bdry_cond(4)==4 .and. coords(1)==nPx-1) then  
        if (d.eq.1) then
#ifndef OLD_BDRY_COND
            mom(ie  ,:,:) = uaverage*(rho2*c(ie,:,:) + rho1*(1.d0 - c(ie,:,:)))
            mom(ie+1,:,:) = uaverage*(rho2*c(ie+1,:,:) + rho1*(1.d0 - c(ie+1,:,:)))
#else
            mom(ie  ,:,:)= mom(ie-1,:,:)
            mom(ie+1,:,:)=-mom(ie,:,:)
#endif
        else
            mom(ie+1,:,:)=mom(ie,:,:)
        endif
    endif


    if(bdry_cond(2)==0 .and. coords(2)==0    ) then
        if (d.eq.2) then
            mom(:,js-1,:)=0d0
            mom(:,js-2,:)=-mom(:,js,:)
        else
            mom(:,js-1,:)=(2*WallVel(3,1)-u(:,js,:))*(rho2*c(:,js,:) + rho1*(1.d0 - c(:,js,:)))
        endif
    endif
    if(bdry_cond(5)==0 .and. coords(2)==nPy-1) then
        if (d.eq.2) then
            mom(:,je  ,:)=0d0
            mom(:,je+1,:)=-mom(:,je-1,:)
        else
            mom(:,je+1,:)=(2*WallVel(4,1)-u(:,je,:))*(rho2*c(:,je,:) + rho1*(1.d0 - c(:,je,:)))
        endif
    endif
    if(bdry_cond(3)==0 .and. coords(3)==0    ) then
        if (d.eq.3) then
            mom(:,:,ks-1)=0d0
            mom(:,:,ks-2)=-mom(:,:,ks)
        else
            mom(:,:,ks-1)=(2*WallVel(5,1)-u(:,:,ks))*(rho2*c(:,:,ks) + rho1*(1.d0 - c(:,:,ks)))
        endif
    endif
    if(bdry_cond(6)==0 .and. coords(3)==nPz-1) then
        if (d.eq.3) then
            mom(:,:,ke  )=0d0
            mom(:,:,ke+1)=-mom(:,:,ke-1)
        else
            mom(:,:,ke+1)=(2*WallVel(6,1)-u(:,:,ke))*(rho2*c(:,:,ke) + rho1*(1.d0 - c(:,:,ke)))
        endif
    endif
    ! wall boundary condition: shear
    if(bdry_cond(1)==2 .and. coords(1)==0    ) then
        if (d.eq.2) then
            mom(is-1,:,:) = (-dxh(is-1)*WallShear(1,2)+u(is,:,:)) * &
                            (rho2*c(is,:,:) + rho1*(1.d0 - c(is,:,:)))
        elseif (d.eq.3) then
            mom(is-1,:,:) = (-dxh(is-1)*WallShear(1,3)+u(is,:,:)) * &
                            (rho2*c(is,:,:) + rho1*(1.d0 - c(is,:,:)))
        endif
    endif
    if(bdry_cond(4)==2 .and. coords(1)==nPx-1) then
        if (d.eq.2) then
            mom(ie+1,:,:) = (dxh(ie)*WallShear(2,2)+u(ie,:,:)) * &
                            (rho2*c(ie,:,:) + rho1*(1.d0 - c(ie,:,:)))
        elseif (d.eq.3) then
            mom(ie+1,:,:) = (dxh(ie)*WallShear(2,3)+u(ie,:,:)) * &
                            (rho2*c(ie,:,:) + rho1*(1.d0 - c(ie,:,:)))
        endif
    endif
    if(bdry_cond(2)==2 .and. coords(2)==0    ) then
        if (d.eq.1) then
            mom(:,js-1,:) = (-dyh(js-1)*WallShear(3,1)+u(:,js,:)) * &
                            (rho2*c(:,js,:) + rho1*(1.d0 - c(:,js,:)))
        elseif (d.eq.3) then
            mom(:,js-1,:) = (-dyh(js-1)*WallShear(3,3)+u(:,js,:)) * &
                            (rho2*c(:,js,:) + rho1*(1.d0 - c(:,js,:)))
        endif
    endif
    if(bdry_cond(5)==2 .and. coords(2)==nPy-1) then
        if (d.eq.1) then
            mom(:,je+1,:) = (dyh(je)*WallShear(4,1)+u(:,je,:)) * &
                            (rho2*c(:,je,:) + rho1*(1.d0 - c(:,je,:)))
        elseif (d.eq.3) then
            mom(:,je+1,:) = (dyh(je)*WallShear(4,3)+u(:,je,:)) * &
                            (rho2*c(:,je,:) + rho1*(1.d0 - c(:,je,:)))
        endif
    endif
    if(bdry_cond(3)==2 .and. coords(3)==0    ) then
        if (d.eq.1) then
            mom(:,:,ks-1) = (-dzh(ks-1)*WallShear(5,1)+u(:,:,ks)) * &
                            (rho2*c(:,:,ks) + rho1*(1.d0 - c(:,:,ks)))
        elseif (d.eq.2) then
            mom(:,:,ks-1) = (-dzh(ks-1)*WallShear(5,2)+u(:,:,ks)) * &
                            (rho2*c(:,:,ks) + rho1*(1.d0 - c(:,:,ks)))
        endif
    endif
    if(bdry_cond(6)==2 .and. coords(3)==nPz-1) then
        if (d.eq.1) then
            mom(:,:,ke+1) = (dzh(ke)*WallShear(6,1)+u(:,:,ke)) * &
                            (rho2*c(:,:,ke) + rho1*(1.d0 - c(:,:,ke)))
        elseif (d.eq.2) then
            mom(:,:,ke+1) = (dzh(ke)*WallShear(6,2)+u(:,:,ke)) * &
                            (rho2*c(:,:,ke) + rho1*(1.d0 - c(:,:,ke)))
        endif
    endif
    
    !Set zero normal velocity gradient for pressure boundary condition
    if (bdry_cond(1)==5 .and. coords(1)==0) then
        if (d.eq.1) then
           mom(is-1,:,:) = mom(is,:,:)
#ifndef OLD_BDRY_COND  
           mom(is-2,:,:) = mom(is,:,:)
#else
           mom(is-2,:,:)=-mom(is,:,:)
#endif
       else
           mom(is-1,:,:)=mom(is,:,:)
       endif
    endif
    
    if (bdry_cond(4)==5 .and. coords(1)==nPx-1) then
       if (d.eq.1) then
           mom(ie,:,:)  = mom(ie-1,:,:)
#ifndef OLD_BDRY_COND  
           mom(ie+1,:,:)= mom(ie-1,:,:)
#else
           mom(ie+1,:,:)=-mom(ie-1,:,:)
#endif
       else
           mom(ie+1,:,:)=mom(ie,:,:)
       endif
    endif

! DROP THE IFNDEF OLD_BDRY_COND IN WHAT FOLLOWS

    if (bdry_cond(2)==5 .and. coords(2)==0) then
        if (d.eq.2) then
           mom(:,js-1,:)= mom(:,js,:)
           mom(:,js-2,:)= mom(:,js,:)
        else
           mom(:,js-1,:)=mom(:,js,:)
        endif
    endif
    
    if (bdry_cond(5)==5 .and. coords(2)==nPy-1) then
       if (d.eq.2) then
           mom(:,je,:)  = mom(:,je-1,:)
           mom(:,je+1,:)= mom(:,je-1,:)
       else
           mom(:,je+1,:)=mom(:,je,:)
       endif
    endif
    
    if (bdry_cond(3)==5 .and. coords(3)==0) then
        if (d.eq.3) then
           mom(:,:,ks-1)= mom(:,:,ks)
           mom(:,:,ks-2)= mom(:,:,ks)
        else
           mom(:,:,ks-1)=u(:,:,ks)
        endif
    endif
    
    if (bdry_cond(6)==5 .and. coords(3)==nPz-1) then
        if (d.eq.3) then
            mom(:,:,ke)  = mom(:,:,ke-1)
            mom(:,:,ke+1)= mom(:,:,ke-1)
        else
            mom(:,:,ke+1)=mom(:,:,ke)
        endif
    endif
    
  end subroutine SetMomentumBC

!=================================================================================================
  subroutine do_ghost_vector(us1,us2,us3)
    implicit none
    real(8), dimension(:,:,:) :: us1,us2,us3
    include 'mpif.h'
    integer :: req(48),sta(MPI_STATUS_SIZE,48)
    integer :: ierr

    call ghost_x(us1  ,2,req( 1: 4));  call ghost_x(us2,2,req( 5: 8)); call ghost_x(us3,2,req( 9:12)) 
    call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
    call ghost_y(us1  ,2,req( 1: 4));  call ghost_y(us2,2,req( 5: 8)); call ghost_y(us3,2,req( 9:12)) 
    call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
    call ghost_z(us1  ,2,req( 1: 4));  call ghost_z(us2,2,req( 5: 8)); call ghost_z(us3,2,req( 9:12))
    call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)

  end subroutine do_ghost_vector


    function uinject(j,k,t)
      use module_2phase
      use module_grid
      use module_flow
      implicit none
      integer :: j,k
      real(8) :: t
      real(8) :: uinject
      real(8) :: ryz, low_gas_radius, NozzleThickness
      real(8), parameter :: PI = 3.14159265359d0
      real :: erf
      uinject=0d0

      if(radius_gap_liqgas==0d0) then
         low_gas_radius = radius_liq_inject
      else
         low_gas_radius = radius_gap_liqgas
      endif

      if (inject_type==1) then      ! uniform inflow
         uinject = 1.d0
      elseif( inject_type==2 ) then ! pulsed round jet
         !tdelay_gas_inject = 0.01d0
         if( (y(j) - jetcenter_yc)**2.d0 + (z(k) - jetcenter_zc)**2.d0 .lt. radius_liq_inject**2.d0 ) then 
            uinject=uliq_inject*(1.d0+0.05d0*SIN(10.d0*2.d0*PI*t))
         end if ! y(j)
      elseif( inject_type==5 ) then ! round jet 
         if( (y(j) - jetcenter_yc)**2.d0 + (z(k) - jetcenter_zc)**2.d0 .lt. radius_liq_inject**2.d0 ) then 
            uinject=uliq_inject
         end if ! y(j)
      else if ( inject_type == 3 ) then ! 2d coflowing jet
         NozzleThickness = NozzleThick2Cell*dx(is) 
         if ( y(j) <= radius_liq_inject ) then 
            uinject = uliq_inject & 
                     *erf( (radius_liq_inject - y(j))/blayer_gas_inject ) &
                     *(1.d0 + erf((time-tdelay_gas_inject*0.5d0)/(tdelay_gas_inject*0.25d0)) )*0.5d0
         else if ( y(j) > radius_liq_inject+NozzleThickness .and. y(j) <= radius_gas_inject ) then
            uinject = ugas_inject & 
                     *erf( (y(j) -   radius_liq_inject - NozzleThickness)/blayer_gas_inject ) & 
                     *erf( (radius_gas_inject - y(j))/blayer_gas_inject ) & 
                     *(1.d0 + erf((time-tdelay_gas_inject*0.5d0)/(tdelay_gas_inject*0.25d0)) )*0.5d0
         else 
            uinject = 0.d0 
         end if  !
      else if ( inject_type == 4 ) then ! 3d coaxial jet
         !tdelay_gas_inject = 0.d-2
         ryz = sqrt( (yh(j) - jetcenter_yc)**2.d0 + (zh(k) - jetcenter_zc)**2.d0 )
         if ( ryz < radius_liq_inject ) then 
            uinject = uliq_inject & 
                     *erf( (radius_liq_inject - ryz)/blayer_gas_inject ) &
                     *(1.d0 + erf((time-tdelay_gas_inject*0.5d0)/(tdelay_gas_inject*0.25d0)) )*0.5d0
         else if ( ryz > low_gas_radius .and. ryz < radius_gas_inject ) then
            uinject = ugas_inject & 
                     *erf( (ryz - low_gas_radius)/blayer_gas_inject ) & 
                     *erf( (radius_gas_inject - ryz)/blayer_gas_inject ) & 
                     *(1.d0 + erf((time-tdelay_gas_inject*0.5d0)/(tdelay_gas_inject*0.25d0)) )*0.5d0
         else 
            uinject = 0.d0 
         end if  !
      end if ! y(j), z(k)
    end function uinject

    function uout(j,k,t)
      use module_2phase
      use module_grid
      use module_flow
      implicit none
      integer :: j,k
      real(8) :: t
      real(8) :: u0,y0,y1,ycjet,uout
      real(8), parameter :: SpreadRate = 0.1d0 
      real(8), parameter :: alpha = 0.88d0 
   
      if ( inject_type == 3 ) then ! 2d jet
         ! planar jet similar velocity profile u=u0*f(y1), Pope 2000
         ! int(f^2(y1))=4/(3*alpha)
         y0    = SpreadRate*xLength
         ycjet = radius_liq_inject+0.5d0*radius_gas_inject
         ! liquid jet is negelected when u0 is computed
         u0    = sqrt( ugas_inject*ugas_inject*(radius_gas_inject  & 
                                              - radius_liq_inject) &
                      /(y0*4.d0/3.d0/alpha) )
         if ( y(j) > ycjet ) then 
            y1    = (y(j)-ycjet)/y0
            uout = u0/COSH(alpha*y1)/COSH(alpha*y1)
         else 
            uout = u0
         end if ! y(j)
      end if ! inject_type
    end function uout
!=================================================================================================
!=================================================================================================
! subroutine SetVectorBC: Sets the boundary condition for vectors (called in Front routines)
!-------------------------------------------------------------------------------------------------
  subroutine SetVectorBC(fx,fy,fz)
    use module_grid
    
    implicit none
    include 'mpif.h'
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: fx, fy, fz
    ! Add the vector outside the domain to the neighboring cell inside the domain
    ! This is used for color gradient vector and surface tension forces
    if(bdry_cond(1)==0 .and. coords(1)==0    ) then
        fx(is-1,:,:)=fx(is-1,:,:)+fx(is-2,:,:)
        fy(is  ,:,:)=fy(is  ,:,:)+fy(is-1,:,:)+fy(is-2,:,:)
        fz(is  ,:,:)=fz(is  ,:,:)+fz(is-1,:,:)+fz(is-2,:,:)
    endif
    if(bdry_cond(4)==0 .and. coords(1)==nPx-1) then
        fx(ie,:,:)=fx(ie,:,:)+fx(ie+1,:,:)
        fy(ie,:,:)=fy(ie,:,:)+fy(ie+1,:,:)+fy(ie+2,:,:)
        fz(ie,:,:)=fz(ie,:,:)+fz(ie+1,:,:)+fz(ie+2,:,:)
    endif
    if(bdry_cond(2)==0 .and. coords(2)==0    ) then
        fy(:,js-1,:)=fy(:,js-1,:)+fy(:,js-2,:)
        fx(:,js  ,:)=fx(:,js  ,:)+fx(:,js-1,:)+fx(:,js-2,:)
        fz(:,js  ,:)=fz(:,js  ,:)+fz(:,js-1,:)+fz(:,js-2,:)
    endif
    if(bdry_cond(5)==0 .and. coords(2)==nPy-1) then
        fy(:,je,:)=fy(:,je,:)+fy(:,je+1,:)
        fx(:,je,:)=fx(:,je,:)+fx(:,je+1,:)+fx(:,je+2,:)
        fz(:,je,:)=fz(:,je,:)+fz(:,je+1,:)+fz(:,je+2,:)
    endif
    if(bdry_cond(3)==0 .and. coords(3)==0    ) then
        fz(:,:,ks-1)=fz(:,:,ks-1)+fz(:,:,ks-2)
        fx(:,:,ks  )=fx(:,:,ks  )+fx(:,:,ks-1)+fx(:,:,ks-2)
        fy(:,:,ks  )=fy(:,:,ks  )+fy(:,:,ks-1)+fy(:,:,ks-2)
    endif
    if(bdry_cond(6)==0 .and. coords(3)==nPz-1) then
        fz(:,:,ke)=fz(:,:,ke)+fz(:,:,ke+1)
        fx(:,:,ke)=fx(:,:,ke)+fx(:,:,ke+1)+fx(:,:,ke+2)
        fy(:,:,ke)=fy(:,:,ke)+fy(:,:,ke+1)+fy(:,:,ke+2)
    endif
  end subroutine SetVectorBC
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
  subroutine ghost_x(F,ngh,req)
    use module_grid
    
    implicit none
    include 'mpif.h'
    integer, intent(in) :: ngh ! number of ghost cell layers to fill
    integer, intent(out) :: req(4)
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: F
    integer :: jlen, klen, ierr !,sta(MPI_STATUS_SIZE,4)
    integer, save :: srcL, srcR, destL, destR, face(2)
    logical, save :: first_time=.true.

    if(ngh>Ng) call pariserror("ghost error: not enough ghost layers to fill")
    if(first_time) then
      first_time=.false.
      jlen=jmax-jmin+1; klen=kmax-kmin+1; !ilen=ngh
      call para_type_block3a(imin, imax, jmin, jmax, 1, jlen, klen, MPI_DOUBLE_PRECISION, face(1))
      call para_type_block3a(imin, imax, jmin, jmax, 2, jlen, klen, MPI_DOUBLE_PRECISION, face(2))
      call MPI_CART_SHIFT(MPI_COMM_CART, 0, 1, srcR, destR, ierr)
      call MPI_CART_SHIFT(MPI_COMM_CART, 0,-1, srcL, destL, ierr)
    endif

    call MPI_IRECV(F(is-ngh  ,jmin,kmin),1,face(ngh),srcR ,0,MPI_COMM_Cart,req(1),ierr)
    call MPI_ISEND(F(ie-ngh+1,jmin,kmin),1,face(ngh),destR,0,MPI_COMM_Cart,req(2),ierr)
    call MPI_IRECV(F(ie+1    ,jmin,kmin),1,face(ngh),srcL ,0,MPI_COMM_Cart,req(3),ierr)
    call MPI_ISEND(F(is      ,jmin,kmin),1,face(ngh),destL,0,MPI_COMM_Cart,req(4),ierr)
!    call MPI_WAITALL(4,req,sta,ierr)
  end subroutine ghost_x
!-------------------------------------------------------------------------------------------------
  subroutine ghost_y(F,ngh,req)
    use module_grid
    
    implicit none
    include 'mpif.h'
    integer, intent(in) :: ngh
    integer, intent(out) :: req(4)
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: F
    integer :: ilen, klen, ierr !,sta(MPI_STATUS_SIZE,4)
    integer, save :: srcL, srcR, destL, destR, face(2)
    logical, save :: first_time=.true.

    if(ngh>Ng) call pariserror("ghost error: not enough ghost layers to fill")
    if(first_time)then
      first_time=.false.
      klen=kmax-kmin+1; ilen=imax-imin+1; !jlen=ngh
      call para_type_block3a(imin, imax, jmin, jmax, ilen, 1, klen, MPI_DOUBLE_PRECISION, face(1))
      call para_type_block3a(imin, imax, jmin, jmax, ilen, 2, klen, MPI_DOUBLE_PRECISION, face(2))
      call MPI_CART_SHIFT(MPI_COMM_CART, 1, 1, srcR, destR, ierr)
      call MPI_CART_SHIFT(MPI_COMM_CART, 1,-1, srcL, destL, ierr)
    endif

    call MPI_IRECV(F(imin,js-ngh  ,kmin),1,face(ngh),srcR ,0,MPI_COMM_Cart,req(1),ierr)
    call MPI_ISEND(F(imin,je-ngh+1,kmin),1,face(ngh),destR,0,MPI_COMM_Cart,req(2),ierr)
    call MPI_IRECV(F(imin,je+1    ,kmin),1,face(ngh),srcL ,0,MPI_COMM_Cart,req(3),ierr)
    call MPI_ISEND(F(imin,js      ,kmin),1,face(ngh),destL,0,MPI_COMM_Cart,req(4),ierr)
!    call MPI_WAITALL(4,req,sta,ierr)
  end subroutine ghost_y
!-------------------------------------------------------------------------------------------------
  subroutine ghost_z(F,ngh,req)
    use module_grid
    
    implicit none
    include 'mpif.h'
    integer, intent(in) :: ngh
    integer, intent(out) :: req(4)
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: F
    integer :: ilen, jlen, ierr !,sta(MPI_STATUS_SIZE,4)
    integer, save :: srcL, srcR, destL, destR, face(2)
    logical, save :: first_time=.true.

    if(ngh>Ng) call pariserror("ghost error: not enough ghost layers to fill")
    if(first_time)then
      first_time=.false.
      ilen=imax-imin+1; jlen=jmax-jmin+1; !klen=ngh
      call para_type_block3a(imin, imax, jmin, jmax, ilen, jlen, 1, MPI_DOUBLE_PRECISION, face(1))
      call para_type_block3a(imin, imax, jmin, jmax, ilen, jlen, 2, MPI_DOUBLE_PRECISION, face(2))
      call MPI_CART_SHIFT(MPI_COMM_CART, 2, 1, srcR, destR, ierr)
      call MPI_CART_SHIFT(MPI_COMM_CART, 2,-1, srcL, destL, ierr)
    endif

    call MPI_IRECV(F(imin,jmin,ks-ngh  ),1,face(ngh),srcR ,0,MPI_COMM_Cart,req(1),ierr)
    call MPI_ISEND(F(imin,jmin,ke-ngh+1),1,face(ngh),destR,0,MPI_COMM_Cart,req(2),ierr)
    call MPI_IRECV(F(imin,jmin,ke+1    ),1,face(ngh),srcL ,0,MPI_COMM_Cart,req(3),ierr)
    call MPI_ISEND(F(imin,jmin,ks      ),1,face(ngh),destL,0,MPI_COMM_Cart,req(4),ierr)
!    call MPI_WAITALL(4,req,sta,ierr)
  end subroutine ghost_z
!=================================================================================================
!=================================================================================================
  subroutine ighost_x(F,ngh,req)
    use module_grid
    
    implicit none
    include 'mpif.h'
    integer, intent(in) :: ngh ! number of ghost cell layers to fill
    integer, intent(out) :: req(4)
    integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: F
    integer :: jlen, klen, ierr !,sta(MPI_STATUS_SIZE,4)
    integer, save :: srcL, srcR, destL, destR, face(2)
    logical, save :: first_time=.true.

    if(ngh>Ng) call pariserror("ghost error: not enough ghost layers to fill")
    if(first_time) then
      first_time=.false.
      jlen=jmax-jmin+1; klen=kmax-kmin+1; !ilen=ngh
      call para_type_block3a(imin, imax, jmin, jmax, 1, jlen, klen, MPI_INTEGER, face(1))
      call para_type_block3a(imin, imax, jmin, jmax, 2, jlen, klen, MPI_INTEGER, face(2))
      call MPI_CART_SHIFT(MPI_COMM_CART, 0, 1, srcR, destR, ierr)
      call MPI_CART_SHIFT(MPI_COMM_CART, 0,-1, srcL, destL, ierr)
    endif

    call MPI_IRECV(F(is-ngh  ,jmin,kmin),1,face(ngh),srcR ,0,MPI_COMM_Cart,req(1),ierr)
    call MPI_ISEND(F(ie-ngh+1,jmin,kmin),1,face(ngh),destR,0,MPI_COMM_Cart,req(2),ierr)
    call MPI_IRECV(F(ie+1    ,jmin,kmin),1,face(ngh),srcL ,0,MPI_COMM_Cart,req(3),ierr)
    call MPI_ISEND(F(is      ,jmin,kmin),1,face(ngh),destL,0,MPI_COMM_Cart,req(4),ierr)
!    call MPI_WAITALL(4,req,sta,ierr)
  end subroutine ighost_x
!-------------------------------------------------------------------------------------------------
  subroutine ighost_y(F,ngh,req)
    use module_grid
    
    implicit none
    include 'mpif.h'
    integer, intent(in) :: ngh
    integer, intent(out) :: req(4)
    integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: F
    integer :: ilen, klen, ierr !,sta(MPI_STATUS_SIZE,4)
    integer, save :: srcL, srcR, destL, destR, face(2)
    logical, save :: first_time=.true.

    if(ngh>Ng) call pariserror("ghost error: not enough ghost layers to fill")
    if(first_time)then
      first_time=.false.
      klen=kmax-kmin+1; ilen=imax-imin+1; !jlen=ngh
      call para_type_block3a(imin, imax, jmin, jmax, ilen, 1, klen, MPI_INTEGER, face(1))
      call para_type_block3a(imin, imax, jmin, jmax, ilen, 2, klen, MPI_INTEGER, face(2))
      call MPI_CART_SHIFT(MPI_COMM_CART, 1, 1, srcR, destR, ierr)
      call MPI_CART_SHIFT(MPI_COMM_CART, 1,-1, srcL, destL, ierr)
    endif

    call MPI_IRECV(F(imin,js-ngh  ,kmin),1,face(ngh),srcR ,0,MPI_COMM_Cart,req(1),ierr)
    call MPI_ISEND(F(imin,je-ngh+1,kmin),1,face(ngh),destR,0,MPI_COMM_Cart,req(2),ierr)
    call MPI_IRECV(F(imin,je+1    ,kmin),1,face(ngh),srcL ,0,MPI_COMM_Cart,req(3),ierr)
    call MPI_ISEND(F(imin,js      ,kmin),1,face(ngh),destL,0,MPI_COMM_Cart,req(4),ierr)
!    call MPI_WAITALL(4,req,sta,ierr)
  end subroutine ighost_y
!-------------------------------------------------------------------------------------------------
  subroutine ighost_z(F,ngh,req)
    use module_grid
    
    implicit none
    include 'mpif.h'
    integer, intent(in) :: ngh
    integer, intent(out) :: req(4)
    integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: F
    integer :: ilen, jlen, ierr !,sta(MPI_STATUS_SIZE,4)
    integer, save :: srcL, srcR, destL, destR, face(2)
    logical, save :: first_time=.true.

    if(ngh>Ng) call pariserror("ghost error: not enough ghost layers to fill")
    if(first_time)then
      first_time=.false.
      ilen=imax-imin+1; jlen=jmax-jmin+1; !klen=ngh
      call para_type_block3a(imin, imax, jmin, jmax, ilen, jlen, 1, MPI_INTEGER, face(1))
      call para_type_block3a(imin, imax, jmin, jmax, ilen, jlen, 2, MPI_INTEGER, face(2))
      call MPI_CART_SHIFT(MPI_COMM_CART, 2, 1, srcR, destR, ierr)
      call MPI_CART_SHIFT(MPI_COMM_CART, 2,-1, srcL, destL, ierr)
    endif

    call MPI_IRECV(F(imin,jmin,ks-ngh  ),1,face(ngh),srcR ,0,MPI_COMM_Cart,req(1),ierr)
    call MPI_ISEND(F(imin,jmin,ke-ngh+1),1,face(ngh),destR,0,MPI_COMM_Cart,req(2),ierr)
    call MPI_IRECV(F(imin,jmin,ke+1    ),1,face(ngh),srcL ,0,MPI_COMM_Cart,req(3),ierr)
    call MPI_ISEND(F(imin,jmin,ks      ),1,face(ngh),destL,0,MPI_COMM_Cart,req(4),ierr)
!    call MPI_WAITALL(4,req,sta,ierr)
  end subroutine ighost_z
!=================================================================================================
subroutine lghost_x(F,ngh,req)
    use module_grid    
    implicit none
    include 'mpif.h'
    integer, intent(in) :: ngh ! number of ghost cell layers to fill
    integer, intent(out) :: req(4)
    logical, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: F
    integer :: jlen, klen, ierr !,sta(MPI_STATUS_SIZE,4)
    integer, save :: srcL, srcR, destL, destR, face(2)
    logical, save :: first_time=.true.

    if(ngh>Ng) call pariserror("ghost error: not enough ghost layers to fill")
    if(first_time) then
      first_time=.false.
      jlen=jmax-jmin+1; klen=kmax-kmin+1; !ilen=ngh
      call para_type_block3a(imin, imax, jmin, jmax, 1, jlen, klen, MPI_LOGICAL, face(1))
      call para_type_block3a(imin, imax, jmin, jmax, 2, jlen, klen, MPI_LOGICAL, face(2))
      call MPI_CART_SHIFT(MPI_COMM_CART, 0, 1, srcR, destR, ierr)
      call MPI_CART_SHIFT(MPI_COMM_CART, 0,-1, srcL, destL, ierr)
    endif

    call MPI_IRECV(F(is-ngh  ,jmin,kmin),1,face(ngh),srcR ,0,MPI_COMM_Cart,req(1),ierr)
    call MPI_ISEND(F(ie-ngh+1,jmin,kmin),1,face(ngh),destR,0,MPI_COMM_Cart,req(2),ierr)
    call MPI_IRECV(F(ie+1    ,jmin,kmin),1,face(ngh),srcL ,0,MPI_COMM_Cart,req(3),ierr)
    call MPI_ISEND(F(is      ,jmin,kmin),1,face(ngh),destL,0,MPI_COMM_Cart,req(4),ierr)
!    call MPI_WAITALL(4,req,sta,ierr)
  end subroutine lghost_x
!-------------------------------------------------------------------------------------------------
  subroutine lghost_y(F,ngh,req)
    use module_grid
    
    implicit none
    include 'mpif.h'
    integer, intent(in) :: ngh
    integer, intent(out) :: req(4)
    logical, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: F
    integer :: ilen, klen, ierr !,sta(MPI_STATUS_SIZE,4)
    integer, save :: srcL, srcR, destL, destR, face(2)
    logical, save :: first_time=.true.

    if(ngh>Ng) call pariserror("ghost error: not enough ghost layers to fill")
    if(first_time)then
      first_time=.false.
      klen=kmax-kmin+1; ilen=imax-imin+1; !jlen=ngh
      call para_type_block3a(imin, imax, jmin, jmax, ilen, 1, klen, MPI_LOGICAL, face(1))
      call para_type_block3a(imin, imax, jmin, jmax, ilen, 2, klen, MPI_LOGICAL, face(2))
      call MPI_CART_SHIFT(MPI_COMM_CART, 1, 1, srcR, destR, ierr)
      call MPI_CART_SHIFT(MPI_COMM_CART, 1,-1, srcL, destL, ierr)
    endif

    call MPI_IRECV(F(imin,js-ngh  ,kmin),1,face(ngh),srcR ,0,MPI_COMM_Cart,req(1),ierr)
    call MPI_ISEND(F(imin,je-ngh+1,kmin),1,face(ngh),destR,0,MPI_COMM_Cart,req(2),ierr)
    call MPI_IRECV(F(imin,je+1    ,kmin),1,face(ngh),srcL ,0,MPI_COMM_Cart,req(3),ierr)
    call MPI_ISEND(F(imin,js      ,kmin),1,face(ngh),destL,0,MPI_COMM_Cart,req(4),ierr)
!    call MPI_WAITALL(4,req,sta,ierr)
  end subroutine lghost_y
!-------------------------------------------------------------------------------------------------
  subroutine lghost_z(F,ngh,req)
    use module_grid
    
    implicit none
    include 'mpif.h'
    integer, intent(in) :: ngh
    integer, intent(out) :: req(4)
    logical, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: F
    integer :: ilen, jlen, ierr !,sta(MPI_STATUS_SIZE,4)
    integer, save :: srcL, srcR, destL, destR, face(2)
    logical, save :: first_time=.true.

    if(ngh>Ng) call pariserror("ghost error: not enough ghost layers to fill")
    if(first_time)then
      first_time=.false.
      ilen=imax-imin+1; jlen=jmax-jmin+1; !klen=ngh
      call para_type_block3a(imin, imax, jmin, jmax, ilen, jlen, 1, MPI_LOGICAL, face(1))
      call para_type_block3a(imin, imax, jmin, jmax, ilen, jlen, 2, MPI_LOGICAL, face(2))
      call MPI_CART_SHIFT(MPI_COMM_CART, 2, 1, srcR, destR, ierr)
      call MPI_CART_SHIFT(MPI_COMM_CART, 2,-1, srcL, destL, ierr)
    endif

    call MPI_IRECV(F(imin,jmin,ks-ngh  ),1,face(ngh),srcR ,0,MPI_COMM_Cart,req(1),ierr)
    call MPI_ISEND(F(imin,jmin,ke-ngh+1),1,face(ngh),destR,0,MPI_COMM_Cart,req(2),ierr)
    call MPI_IRECV(F(imin,jmin,ke+1    ),1,face(ngh),srcL ,0,MPI_COMM_Cart,req(3),ierr)
    call MPI_ISEND(F(imin,jmin,ks      ),1,face(ngh),destL,0,MPI_COMM_Cart,req(4),ierr)
!    call MPI_WAITALL(4,req,sta,ierr)
  end subroutine lghost_z
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
  subroutine ghost_xAdd(F,ir1,is1,iwork,req)
    use module_grid
    
    use module_tmpvar
    implicit none
    include 'mpif.h'
    integer, intent(in) :: ir1, is1, iwork
    integer, intent(out) :: req(4)
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: F
    integer :: jlen, klen, ierr !,sta(MPI_STATUS_SIZE,4)
    integer, save :: srcL, srcR, destL, destR, face
    logical, save :: first_time=.true.

    if(first_time)then
      first_time=.false.
      jlen=jmax-jmin+1; klen=kmax-kmin+1; !ilen=ngh
      call para_type_block3a(imin, imax, jmin, jmax, 2, jlen, klen, MPI_DOUBLE_PRECISION, face)
      call MPI_CART_SHIFT(MPI_COMM_CART, 0, 1, srcR, destR, ierr)
      call MPI_CART_SHIFT(MPI_COMM_CART, 0,-1, srcL, destL, ierr)
    endif

    work(ir1:ir1+1,jmin:jmax,kmin:kmax,iwork) = 0d0
    work(ie-1:ie  ,jmin:jmax,kmin:kmax,iwork) = 0d0
    call MPI_IRECV(work(ir1 ,jmin,kmin,iwork),1,face,srcR ,0,MPI_Comm_Cart,req(1),ierr)
    call MPI_ISEND(F   (is1 ,jmin,kmin      ),1,face,destR,0,MPI_Comm_Cart,req(2),ierr)
    call MPI_IRECV(work(ie-1,jmin,kmin,iwork),1,face,srcL ,0,MPI_Comm_Cart,req(3),ierr)
    call MPI_ISEND(F   (is-2,jmin,kmin      ),1,face,destL,0,MPI_Comm_Cart,req(4),ierr)
!    call MPI_WAITALL(4,req,sta,ierr)

  end subroutine ghost_xAdd
!-------------------------------------------------------------------------------------------------
  subroutine ghost_yAdd(F,jr1,js1,iwork,req)
    use module_grid
    
    use module_tmpvar
    implicit none
    include 'mpif.h'
    integer, intent(in) :: jr1, js1, iwork
    integer, intent(out) :: req(4)
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: F
    integer :: ilen, klen, ierr !,sta(MPI_STATUS_SIZE,4)
    integer, save :: srcL, srcR, destL, destR, face
    logical, save :: first_time=.true.

    if(first_time)then
      first_time=.false.
      ilen=imax-imin+1; klen=kmax-kmin+1; !ilen=ngh
      call para_type_block3a(imin, imax, jmin, jmax, ilen, 2, klen, MPI_DOUBLE_PRECISION, face)
      call MPI_CART_SHIFT(MPI_COMM_CART, 1, 1, srcR, destR, ierr)
      call MPI_CART_SHIFT(MPI_COMM_CART, 1,-1, srcL, destL, ierr)
    endif

    work(imin:imax,jr1:jr1+1,kmin:kmax,iwork) = 0d0
    work(imin:imax,je-1:je  ,kmin:kmax,iwork) = 0d0
    call MPI_IRECV(work(imin,jr1 ,kmin,iwork),1,face,srcR ,0,MPI_Comm_Cart,req(1),ierr)
    call MPI_ISEND(F   (imin,js1 ,kmin      ),1,face,destR,0,MPI_Comm_Cart,req(2),ierr)
    call MPI_IRECV(work(imin,je-1,kmin,iwork),1,face,srcL ,0,MPI_Comm_Cart,req(3),ierr)
    call MPI_ISEND(F   (imin,js-2,kmin      ),1,face,destL,0,MPI_Comm_Cart,req(4),ierr)
!    call MPI_WAITALL(4,req,sta,ierr)
  end subroutine ghost_yAdd
!-------------------------------------------------------------------------------------------------
  subroutine ghost_zAdd(F,kr1,ks1,iwork,req)
    use module_grid

      
    use module_tmpvar
    implicit none
    include 'mpif.h'
    integer, intent(in) :: kr1, ks1, iwork
    integer, intent(out) :: req(4)
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: F
    integer :: ilen, jlen, ierr !,sta(MPI_STATUS_SIZE,4)
    integer, save :: srcL, srcR, destL, destR, face
    logical, save :: first_time=.true.

    if(first_time)then
      first_time=.false.
      ilen=imax-imin+1; jlen=jmax-jmin+1; !ilen=ngh
      call para_type_block3a(imin, imax, jmin, jmax, ilen, jlen, 2, MPI_DOUBLE_PRECISION, face)
      call MPI_CART_SHIFT(MPI_COMM_CART, 2, 1, srcR, destR, ierr)
      call MPI_CART_SHIFT(MPI_COMM_CART, 2,-1, srcL, destL, ierr)
    endif

    work(imin:imax,jmin:jmax,kr1:kr1+1,iwork) = 0d0
    work(imin:imax,jmin:jmax,ke-1:ke  ,iwork) = 0d0
    call MPI_IRECV(work(imin,jmin,kr1 ,iwork),1,face,srcR ,0,MPI_Comm_Cart,req(1),ierr)
    call MPI_ISEND(F   (imin,jmin,ks1       ),1,face,destR,0,MPI_Comm_Cart,req(2),ierr)
    call MPI_IRECV(work(imin,jmin,ke-1,iwork),1,face,srcL ,0,MPI_Comm_Cart,req(3),ierr)
    call MPI_ISEND(F   (imin,jmin,ks-2      ),1,face,destL,0,MPI_Comm_Cart,req(4),ierr)
!    call MPI_WAITALL(4,req,sta,ierr)
  end subroutine ghost_zAdd
!-------------------------------------------------------------------------------------------------
end module module_BC
!=================================================================================================
!=================================================================================================
! module_poisson:
! Using hypre library, solves linear equations with matrix A and solution P as
! A7*Pijk = A1*Pi-1jk + A2*Pi+1jk + A3*Pij-1k + 
!           A4*Pij+1k + A5*Pijk-1 + A6*Pijk+1 + A8
! 
! Syntax: call poi_initialize(mpi_comm_in,is,ie,js,je,ks,ke,Nx,Ny,Nz,bdry_cond)
!
!   Input:  integer :: mpi_comm_in               MPI communicator
!                      is,ie,js,je,ks,ke         bounds of each subdomain
!                      Nx,Ny,Nz                  total size of computational domain
!                      bdry_cond(3)              boundary conditions: use 1 for periodic
!
! Syntax: call poi_solve(A,p,maxError,maxIteration,numIteration)
!
!   Input:  real(8) :: A(is:ie,js:je,ks:ke,1:8)  coefficients and source term
!                      maxError                  criterion of convergence
!                      p(is:ie,js:je,ks:ke)      initial guess
!           integer :: maxIteration              maximum number of iterations
!
!   Output: real(8) :: p(is:ie,js:je,ks:ke)      solution
!           integer :: numIteration              number of iterations performed
!
! Syntax: call poi_finalize
!
! written by Sadegh Dabiri sdabiri@nd.edu 10/2011
!-------------------------------------------------------------------------------------------------
module module_poisson
  implicit none
  save
  private
  public :: poi_initialize, poi_solve, poi_finalize
  interface
     FUNCTION GetNumIterationsSMG (isolver,inum) &
          bind(C, name="HYPRE_StructSMGGetNumIterations")
       integer :: GetNumIterationsSMG
       integer :: inum
       integer (kind = 8) , VALUE :: isolver
     END FUNCTION GetNumIterationsSMG
  end interface
  interface
     FUNCTION GetNumIterationsPFMG (isolver,inum) &
          bind(C, name="HYPRE_StructPFMGGetNumIterations")
       integer :: GetNumIterationsPFMG
       integer :: inum
       integer (kind = 8) , VALUE :: isolver
     END FUNCTION GetNumIterationsPFMG
  end interface
  interface
     FUNCTION GetFinalRelative (isolver,rnum) &
          bind(C, name="HYPRE_StructSMGGetFinalRelativeResidualNorm")
       integer :: GetFinalRelative
       real(8) :: rnum
       integer (kind = 8) , VALUE :: isolver
     END FUNCTION GetFinalRelative
  end interface

  integer :: nstencil
  integer (kind = 8) :: grid_obj, stencil, Amat, Bvec, Xvec, solver
  integer, dimension(:), allocatable :: stencil_indices, ilower, iupper
  integer :: mpi_comm_poi,is,ie,js,je,ks,ke,Mx,My,Mz
  integer, parameter :: ndim=3
  contains

!=================================================================================================
!=================================================================================================
  subroutine poi_initialize(mpi_comm_in,iis,iie,jjs,jje,kks,kke,Nx,Ny,Nz,bdry_cond)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: mpi_comm_in,iis,iie,jjs,jje,kks,kke,Nx,Ny,Nz,bdry_cond(3)
    integer :: ierr, periodic_array(3), i
    integer, dimension(:,:), allocatable :: offsets
    
    mpi_comm_poi = mpi_comm_in
    is = iis;   ie = iie;   Mx = ie-is+1
    js = jjs;   je = jje;   My = je-js+1
    ks = kks;   ke = kke;   Mz = ke-ks+1

    nstencil = 2 * ndim + 1
    allocate(ilower(ndim), iupper(ndim), stencil_indices(nstencil) )
    ilower = (/is, js, ks/);  iupper = (/ie, je, ke/)
    periodic_array = 0
    if( bdry_cond(1)==1) periodic_array(1) = Nx
    if( bdry_cond(2)==1) periodic_array(2) = Ny
    if( bdry_cond(3)==1) periodic_array(3) = Nz

    call HYPRE_StructGridCreate(mpi_comm_poi, ndim, grid_obj, ierr)    ! create a 3d grid object
    call HYPRE_StructGridSetExtents(grid_obj, ilower, iupper, ierr)    ! add a box to the grid
    call HYPRE_StructGridSetPeriodic(grid_obj, periodic_array, ierr)   ! set periodic
    call HYPRE_StructGridAssemble(grid_obj, ierr)                      ! assemble the grid
    call HYPRE_StructStencilCreate(ndim, nstencil, stencil, ierr)      ! create a 3d 7-pt stencil

    allocate(offsets(ndim,nstencil), stat=ierr) ! define the geometry of the stencil
    offsets(:,1) = (/ 0, 0, 0/)
    offsets(:,2) = (/-1, 0, 0/)
    offsets(:,3) = (/ 1, 0, 0/)
    offsets(:,4) = (/ 0,-1, 0/)
    offsets(:,5) = (/ 0, 1, 0/)
    offsets(:,6) = (/ 0, 0,-1/)
    offsets(:,7) = (/ 0, 0, 1/)
    do i = 1, nstencil
      stencil_indices(i) = i-1
      call HYPRE_StructStencilSetElement(stencil, stencil_indices(i), offsets(:,i), ierr)
    enddo

    call HYPRE_StructMatrixCreate(mpi_comm_poi, grid_obj, stencil, Amat, ierr)! set up matrix Amat
    call HYPRE_StructMatrixInitialize(Amat, ierr)

    call HYPRE_StructVectorCreate(mpi_comm_poi, grid_obj, Bvec, ierr)! create vector object Bvec
    call HYPRE_StructVectorInitialize(Bvec, ierr)

    call HYPRE_StructVectorCreate(mpi_comm_poi, grid_obj, Xvec, ierr)! create vector object Xvec
    call HYPRE_StructVectorInitialize(Xvec, ierr)

  end subroutine poi_initialize
!=================================================================================================
! Solve Poisson equation. p is initial guess. 
!=================================================================================================
  subroutine poi_solve(A,p,maxError,maxit,num_iterations,HYPRESolverType)
    use module_timer
    use module_grid
    use iso_c_binding, only: c_int,c_int8_t
    implicit none
    include 'mpif.h'
    integer :: ierr, nvalues, ijk, i,j,k
    real(8), dimension(:), allocatable :: values
    real(8), dimension(is:ie,js:je,ks:ke,8), intent(in) :: A
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: p
!    real(8), dimension(is-2:ie+2,js-2:je+2,ks-2:ke+2), intent(inout) :: p
    real(8), intent(in) :: maxError
    integer, intent(in) :: maxit
    integer, intent(out):: num_iterations
    integer, intent(in) :: HYPRESolverType
!     real(8) :: final_res_norm
    integer, parameter :: HYPRESolverSMG  = 1 
    integer, parameter :: HYPRESolverPFMG = 2
#ifdef DEBUG_HYPRE
    real(8) :: timeConstruct,timeSetup,timeSolve,timeCleanup,timeTotal
#endif
!----------------------------------------Fill in matrix Amat--------------------------------------
#ifdef DEBUG_HYPRE
   if ( rank == 0 ) timeConstruct=MPI_WTIME(ierr)
#endif
    num_iterations = 0
    nvalues = mx * my * mz * nstencil
    allocate(values(nvalues), stat=ierr)
    call add_2_my_sizer_2(nvalues,8)

    if(ierr/=0)call pariserror("**** poi_solve: allocation error ****")

    ijk = 1
    do k=ks,ke;  do j=js,je;  do i=is,ie
      values(ijk+1) = -A(i,j,k,1)
      values(ijk+2) = -A(i,j,k,2)
      values(ijk+3) = -A(i,j,k,3)
      values(ijk+4) = -A(i,j,k,4)
      values(ijk+5) = -A(i,j,k,5)
      values(ijk+6) = -A(i,j,k,6)
      values(ijk  ) =  A(i,j,k,7)
      ijk = ijk + 7
    enddo;  enddo;  enddo
    call HYPRE_StructMatrixSetBoxValues(Amat, ilower, iupper, nstencil, stencil_indices, &
                                        values, ierr)
    deallocate(values, stat=ierr)
    call remove_from_my_sizer_2(nvalues,8)

!------------------------------Fill in source term and initial guess------------------------------
    nvalues = mx * my * mz
    allocate(values(nvalues), stat=ierr)

    ijk = 1
    do k=ks,ke;  do j=js,je;  do i=is,ie
      values(ijk) = A(i,j,k,8)
      ijk = ijk + 1
    enddo;  enddo;  enddo
    call HYPRE_StructVectorSetBoxValues(Bvec, ilower, iupper, values, ierr)

    ijk = 1
    do k=ks,ke;  do j=js,je;  do i=is,ie
      values(ijk) = p(i,j,k)
      ijk = ijk + 1
    enddo;  enddo;  enddo
    call HYPRE_StructVectorSetBoxValues(Xvec, ilower, iupper, values, ierr)
    deallocate(values, stat=ierr)
!---------------------------Assemble matrix Amat and vectors Bvec,Xvec----------------------------
    call HYPRE_StructMatrixAssemble(Amat, ierr)
    call HYPRE_StructVectorAssemble(Bvec, ierr)
    call HYPRE_StructVectorAssemble(Xvec, ierr)
#ifdef DEBUG_HYPRE
   if (rank == 0 ) then
      timeConstruct = MPI_WTIME(ierr)-timeConstruct
      write(*,*) "timeConstruct: ", timeConstruct
   end if ! rank
#endif
!---------------------------------------Solve the equations---------------------------------------
#ifdef DEBUG_HYPRE
   if ( rank == 0 ) timeSetup=MPI_WTIME(ierr)
#endif
   if ( HYPRESolverType == HYPRESolverSMG ) then 
    call HYPRE_StructSMGCreate(mpi_comm_poi, solver, ierr)
    call HYPRE_StructSMGSetMaxIter(solver, maxit, ierr)
    call HYPRE_StructSMGSetTol(solver, MaxError, ierr)
    call hypre_structSMGsetLogging(solver, 1, ierr)
    call HYPRE_StructSMGSetPrintLevel(solver,1,ierr) 
    call HYPRE_StructSMGSetup(solver, Amat, Bvec, Xvec, ierr)
   else if ( HYPRESolverType == HYPRESolverPFMG ) then  
    call HYPRE_StructPFMGCreate(mpi_comm_poi, solver, ierr)
    call HYPRE_StructPFMGSetMaxIter(solver, maxit, ierr)
    call HYPRE_StructPFMGSetTol(solver, MaxError, ierr)
    call hypre_structPFMGsetLogging(solver, 1, ierr)
    call HYPRE_StructPFMGSetPrintLevel(solver,1,ierr) 
    call HYPRE_StructPFMGSetRelChange(solver, 1, ierr) 
    call HYPRE_StructPFMGSetRelaxType(solver, 3, ierr) 
    !Red/Black Gauss-Seidel (nonsymmetric: RB pre- and post-relaxation)
    call HYPRE_StructPFMGSetup(solver, Amat, Bvec, Xvec, ierr)
   end if ! HYPRESolverType
#ifdef DEBUG_HYPRE
   if (rank == 0 ) then
      timeSetup = MPI_WTIME(ierr)-timeSetup
      write(*,*) "timeSetup: ", timeSetup
   end if ! rank
#endif

#ifdef DEBUG_HYPRE
   if ( rank == 0 ) timeSolve=MPI_WTIME(ierr)
#endif
   if ( HYPRESolverType == HYPRESolverSMG ) then 
      call HYPRE_StructSMGSolve(solver, Amat, Bvec, Xvec, ierr)
      ierr = GetNumIterationsSMG(solver, num_iterations)
   else if ( HYPRESolverType == HYPRESolverPFMG ) then  
      call HYPRE_StructPFMGSolve(solver, Amat, Bvec, Xvec, ierr)
      ierr = GetNumIterationsPFMG(solver, num_iterations)
   end if ! HYPRESolverType
#ifdef DEBUG_HYPRE
   if (rank == 0 ) then
      timeSolve = MPI_WTIME(ierr)-timeSolve
      write(*,*) "timeSolve: ", timeSolve
   end if ! rank
#endif
!    ierr = Getfinalrelative(solver, final_res_norm)
!    print *, "relative error",final_res_norm
!--------------------------------------Retrieve the solution--------------------------------------
#ifdef DEBUG_HYPRE
   if ( rank == 0 ) timeCleanup=MPI_WTIME(ierr)
#endif
    nvalues = mx * my * mz
    allocate(values(nvalues), stat=ierr)
    call HYPRE_StructVectorGetBoxValues(Xvec, ilower, iupper, values, ierr)
    if ( HYPRESolverType == HYPRESolverSMG ) then 
      call HYPRE_StructSMGDestroy(solver, ierr)
    else if ( HYPRESolverType == HYPRESolverPFMG ) then  
      call HYPRE_StructPFMGDestroy(solver, ierr)
    end if ! HYPRESolverType
#ifdef DEBUG_HYPRE
   if (rank == 0 ) then
      timeCleanup = MPI_WTIME(ierr)-timeCleanup
      write(*,*) "timeCleanup: ", timeCleanup
      write(*,*) " ********************************** "
      timeTotal = timeConstruct + timeSetup + timeSolve + timeCleanup
      write(*,*) "timeConstruct: ", timeConstruct/timeTotal
      write(*,*) "timeSetup    : ", timeSetup/timeTotal
      write(*,*) "timeSolve    : ", timeSolve/timeTotal
      write(*,*) "timeCleanup  : ", timeCleanup/timeTotal
   end if ! rank
#endif

    ijk = 1
    do k=ks,ke;  do j=js,je;  do i=is,ie
      p(i,j,k) = values(ijk)
      ijk = ijk + 1
    enddo;  enddo;  enddo
    deallocate(values, stat=ierr)
  end subroutine poi_solve
!=================================================================================================
!=================================================================================================
  subroutine poi_finalize
    implicit none
    include 'mpif.h'
    integer :: ierr
    call HYPRE_StructGridDestroy(grid_obj, ierr)
    call HYPRE_StructStencilDestroy(stencil, ierr)
    call HYPRE_StructMatrixDestroy(Amat, ierr)
    call HYPRE_StructVectorDestroy(Bvec, ierr)
    call HYPRE_StructVectorDestroy(Xvec, ierr)
    deallocate(ilower, iupper, stencil_indices, stat=ierr)
  end subroutine poi_finalize
!-------------------------------------------------------------------------------------------------
end module module_poisson
!=================================================================================================
!=================================================================================================
!=================================================================================================
! creates new data type for passing non-contiguous 
!-------------------------------------------------------------------------------------------------
SUBROUTINE para_type_block3a(imin, imax, jmin, jmax, ilen, jlen, klen, ioldtype,inewtype)
  implicit none
  INCLUDE 'mpif.h'
  integer :: imin, imax, jmin, jmax, ilen, jlen, klen, ioldtype,inewtype,isize, ierr, itemp, idist
  CALL MPI_TYPE_EXTENT(ioldtype, isize, ierr)
  CALL MPI_TYPE_VECTOR(jlen, ilen, imax - imin + 1, ioldtype, itemp, ierr)
  idist = (imax - imin + 1) * (jmax - jmin + 1) * isize
  CALL MPI_TYPE_HVECTOR(klen, 1, idist, itemp, inewtype, ierr)
  CALL MPI_TYPE_COMMIT(inewtype, ierr)
END
!=================================================================================================
!=================================================================================================
! The Poisson equation for the density is setup with matrix A as
! A7*Pijk = A1*Pi-1jk + A2*Pi+1jk + A3*Pij-1k + 
!           A4*Pij+1k + A5*Pijk-1 + A6*Pijk+1 + A8
!-------------------------------------------------------------------------------------------------
subroutine SetupDensity(dIdx,dIdy,dIdz,A,color) !,mask)
  use module_grid
  use module_BC
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: dIdx,dIdy,dIdz, color
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(out) :: A
  integer :: i,j,k
  logical, save :: first=.true.
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax) :: dI
  do k=ks,ke; do j=js,je; do i=is,ie
      A(i,j,k,1) = 1d0/dx(i)/dxh(i-1)
      A(i,j,k,2) = 1d0/dx(i)/dxh(i  )
      A(i,j,k,3) = 1d0/dy(j)/dyh(j-1)
      A(i,j,k,4) = 1d0/dy(j)/dyh(j  )
      A(i,j,k,5) = 1d0/dz(k)/dzh(k-1)
      A(i,j,k,6) = 1d0/dz(k)/dzh(k  )
      A(i,j,k,7) = sum(A(i,j,k,1:6))
      A(i,j,k,8) = -(dIdx(i,j,k)-dIdx(i-1,j,k))/dx(i) &
                   -(dIdy(i,j,k)-dIdy(i,j-1,k))/dy(j) &
                   -(dIdz(i,j,k)-dIdz(i,j,k-1))/dz(k)
  enddo; enddo; enddo

  if(bdry_cond(1)==0)then
    if(coords(1)==0    ) A(is,:,:,8)=A(is,:,:,8)+A(is,:,:,1)
    if(coords(1)==nPx-1) A(ie,:,:,8)=A(ie,:,:,8)+A(ie,:,:,2)
    if(coords(1)==0    ) A(is,:,:,1) = 0d0
    if(coords(1)==nPx-1) A(ie,:,:,2) = 0d0
  endif
  if(bdry_cond(2)==0)then
    if(coords(2)==0    ) A(:,js,:,8)=A(:,js,:,8)+A(:,js,:,3)
    if(coords(2)==nPy-1) A(:,je,:,8)=A(:,je,:,8)+A(:,je,:,4)
    if(coords(2)==0    ) A(:,js,:,3) = 0d0
    if(coords(2)==nPy-1) A(:,je,:,4) = 0d0
  endif
  if(bdry_cond(3)==0)then
    if(coords(3)==0    ) A(:,:,ks,8)=A(:,:,ks,8)+A(:,:,ks,5)
    if(coords(3)==nPz-1) A(:,:,ke,8)=A(:,:,ke,8)+A(:,:,ke,6)
    if(coords(3)==0    ) A(:,:,ks,5) = 0d0
    if(coords(3)==nPz-1) A(:,:,ke,6) = 0d0
  endif

  if(.not. first)then
     first=.false.
     do k=ks,ke; do j=js,je; do i=is,ie
        dI(i,j,k) = ( (dIdx(i,j,k)+dIdx(i-1,j,k))**2 + &
             (dIdx(i,j,k)+dIdx(i-1,j,k))**2 + &
             (dIdx(i,j,k)+dIdx(i-1,j,k))**2 )
        if( dI(i,j,k)<1e-6 )then
           A(i,j,k,1:6) = 0d0
           A(i,j,k,7) = 1d0
           A(i,j,k,8) = float(floor(color(i,j,k)+0.5))
        endif
     enddo; enddo; enddo
  endif

!  do k=ks,ke; do j=js,je; do i=is,ie
!      A(i,j,k,7) = sum(A(i,j,k,1:6))
!  enddo; enddo; enddo
  ! Anchor a point to 1
!  if(coords(1)==0 .and. coords(2)==0 .and. coords(3)==0) then
!    A(is,js,ks,:)=0d0
  !    A(is,js,ks,7:8)=1d0
!  endif
end subroutine SetupDensity
!=================================================================================================
!=================================================================================================
! The Poisson equation for the pressure is setup with matrix A as
! A7*Pijk = A1*Pi-1jk + A2*Pi+1jk + A3*Pij-1k + 
!           A4*Pij+1k + A5*Pijk-1 + A6*Pijk+1 + A8
!-------------------------------------------------------------------------------------------------
subroutine SetupPoisson(utmp,vtmp,wtmp,umask,vmask,wmask,rhot,dt,A,pmask,VolumeSource) 
  use module_grid
  use module_BC
  use module_2phase

  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: utmp,vtmp,wtmp,rhot
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: umask,vmask,wmask
  real(8), dimension(is:ie,js:je,ks:ke), intent(inout) :: pmask
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(out) :: A
  real(8), intent(in) :: dt, VolumeSource
  integer :: i,j,k,l
  
  do k=ks,ke; do j=js,je; do i=is,ie
    !    if(mask(i,j,k))then
    A(i,j,k,1) = 2d0*dt*umask(i-1,j,k)/(dx(i)*dxh(i-1)*(rhot(i-1,j,k)+rhot(i,j,k)))
    A(i,j,k,2) = 2d0*dt*umask(i,j,k)/(dx(i)*dxh(i  )*(rhot(i+1,j,k)+rhot(i,j,k)))
    A(i,j,k,3) = 2d0*dt*vmask(i,j-1,k)/(dy(j)*dyh(j-1)*(rhot(i,j-1,k)+rhot(i,j,k)))
    A(i,j,k,4) = 2d0*dt*vmask(i,j,k)/(dy(j)*dyh(j  )*(rhot(i,j+1,k)+rhot(i,j,k)))
    A(i,j,k,5) = 2d0*dt*wmask(i,j,k-1)/(dz(k)*dzh(k-1)*(rhot(i,j,k-1)+rhot(i,j,k)))
    A(i,j,k,6) = 2d0*dt*wmask(i,j,k)/(dz(k)*dzh(k  )*(rhot(i,j,k+1)+rhot(i,j,k)))
    A(i,j,k,7) = sum(A(i,j,k,1:6))
    A(i,j,k,8) =  -(  VolumeSource +(utmp(i,j,k)-utmp(i-1,j,k))/dx(i) &
    +  (vtmp(i,j,k)-vtmp(i,j-1,k))/dy(j) &
    +  (wtmp(i,j,k)-wtmp(i,j,k-1))/dz(k) )
    !    endif
  enddo; enddo; enddo

  call Poisson_BCs (A)
  call check_and_debug_Poisson(A,umask,vmask,wmask,rhot,pmask,dt,VolumeSource)

end subroutine SetupPoisson

subroutine Poisson_BCs(A)

  use module_grid
  use module_BC
  use module_flow
  
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(out) :: A
  real(8), dimension(4) :: P_bc
  integer :: i,j,k,l

  P_bc = 0d0
  ! dp/dn = 0 for inflow bc on face 1 == x- : do not correct u(is-1)
  ! inflow bc on other faces not implemented yet.  !@@ generalize this ! 
  if(coords(1)==0) then
     if(bdry_cond(1)==3) then  
        A(is,:,:,7) = A(is,:,:,7) - A(is,:,:,1)
        A(is,:,:,1) = 0d0
        ! pressure boundary condition
     else if(bdry_cond(1)==5) then 
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
     if(bdry_cond(4)==4 .or. bdry_cond(4)==3) then
        A(ie,:,:,7) = A(ie,:,:,7) - A(ie,:,:,2)
        A(ie,:,:,2) = 0d0
        ! pressure boundary condition
     else if(bdry_cond(4)==5) then
        A(ie,:,:,8) = (2d0/3d0)*BoundaryPressure(2)
        A(ie,:,:,7) = 1d0  
        A(ie,:,:,2:6) = 0d0
        A(ie,:,:,1) = 1d0/3d0
        P_bc(2) = 1d0
     else if(bdry_cond(4)==6) then
        ! p0 and du/dn right at the boundary
        A(ie,:,:,7) = A(ie,:,:,7) - A(ie,:,:,1)*2.d0/3.d0 + A(ie,:,:,2)*5.d0/3.d0
        A(ie,:,:,8) = A(ie,:,:,8) + 8.d0/3.d0*A(ie,:,:,2)*BoundaryPressure(2) & 
                    + (u(ie,js:je,ks:ke)-u(ie-1,js:je,ks:ke))/dx(ie) ! remove dudx in source term 
        A(ie,:,:,1) = A(ie,:,:,1)*1.d0/3.d0 
        A(ie,:,:,2) = 0.d0
     endif
  endif

  ! Pressure BC for y-
  if(coords(2)==0) then
     if(bdry_cond(2)==3) then
        A(:,js,:,7) = A(:,js,:,7) - A(:,js,:,3)
        A(:,js,:,3) = 0d0
        ! pressure boundary condition
     else if(bdry_cond(2)==5) then
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
  endif
  ! Pressure BC for y+
  if(coords(2)==Npy-1) then
     if(bdry_cond(5)==3) then
        A(:,je,:,7) = A(:,je,:,7) - A(:,je,:,4)
        A(:,je,:,4) = 0d0
        ! pressure boundary condition
     else if(bdry_cond(5)==5) then
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
  endif
  ! Pressure BC for z-
  if(coords(3)==0) then
     if (bdry_cond(3)==3) then
        A(:,:,ks,7) = A(:,:,ks,7) - A(:,:,ks,5)
        A(:,:,ks,5) = 0d0
        ! pressure boundary condition
     else if(bdry_cond(3)==5) then
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
  endif
  ! Pressure BC for z+
  if(coords(3)==Npz-1) then
     if (bdry_cond(6)==3) then
        A(:,:,ke,7) = A(:,:,ke,7) - A(:,:,ke,6)
        A(:,:,ke,6) = 0d0  
     else if(bdry_cond(6)==5) then
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
  end subroutine Poisson_BCs

  subroutine check_and_debug_Poisson(A,umask,vmask,wmask,rhot,pmask,dt,VolumeSource)

  use module_grid
  use module_BC
  use module_IO

  implicit none
  real(8), intent(in) :: VolumeSource,dt
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(inout) :: A
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: pmask
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: umask,vmask,wmask,rhot
  integer :: i,j,k,l

  ! What follows is a lot of debugging for small A7 values and checking the matrix 
  do k=ks,ke; do j=js,je; do i=is,ie
    if(A(i,j,k,7) .lt. 1d-50)  then
      ! check that we are in solid. Remember that we cannot have an isolated fluid cell exactly on the entrance. 
      if(umask(i-1,j,k).lt.0.5d0.and.umask(i,j,k).lt.0.5d0.and.     &
      vmask(i,j-1,k).lt.0.5d0.and.vmask(i,j,k).lt.0.5d0.and.   &
      wmask(i,j,k-1).lt.0.5d0.and.wmask(i,j,k).lt.0.5d0 ) then ! we are in solid
      if(A(i,j,k,8).gt.1d-50) then ! check A8 for debugging
        OPEN(UNIT=88,FILE=TRIM(out_path)//'/message-rank-'//TRIM(int2text(rank,padding))//'.txt')
        write(88,*) "A8 non zero in solid at ijk + minmax = ",i,j,k,imin,imax,jmin,jmax,kmin,kmax
        write(88,*) "VolumeSource",VolumeSource
        write(88,*) "umask(i-1,j,k),umask(i,j,k),vmask(i,j-1,k)",&
        "vmask(i,j,k),wmask(i,j,k-1),wmask(i,j,k)",         &
        umask(i-1,j,k),umask(i,j,k),vmask(i,j-1,k),         &
        vmask(i,j,k),wmask(i,j,k-1),wmask(i,j,k)
        write(88,*) "dx(i),dy(j),dz(k)",dx(i),dy(j),dz(k)
        close(88)
        call pariserror("A8 non zero in solid") 
      endif
      if(maxval(A(i,j,k,1:6)).gt.1d-50.or.minval(A(i,j,k,1:6)).lt.0d0) then
        call pariserror("inconsistency in A1-6")
      endif
      A(i,j,k,7) = 1d0
    else ! we are not in solid: error.
      OPEN(UNIT=88,FILE=TRIM(out_path)//'/message-rank-'//TRIM(int2text(rank,padding))//'.txt')
      write(88,*) "A7 tiny outside of solid at ijk minmax = ",i,j,k,imin,imax,jmin,jmax,kmin,kmax
      write(88,*) "dt",dt
      write(88,*) "umask(i-1,j,k),umask(i,j,k),vmask(i,j-1,k)",&
      "vmask(i,j,k),wmask(i,j,k-1),wmask(i,j,k)",         &
      umask(i-1,j,k),umask(i,j,k),vmask(i,j-1,k),         &
      vmask(i,j,k),wmask(i,j,k-1),wmask(i,j,k)
      write(88,*) "dx(i),dxh(i-1),dxh(i),rhot(i-1,j,k),rhot(i,j,k),rhot(i+1,j,k)",&
      dx(i),dxh(i-1),dxh(i),rhot(i-1,j,k),rhot(i,j,k),rhot(i+1,j,k)
      write(88,*) "dy(j),dyh(j-1),dyh(j),rhot(i,j-1,k),rhot(i,j,k),rhot(i,j+1,k)",&
      dy(j),dyh(j-1),dyh(j),rhot(i,j-1,k),rhot(i,j,k),rhot(i,j+1,k)
      write(88,*) "dz(k),dzh(k-1),dzh(k),rhot(i,j,k-1),rhot(i,j,k),rhot(i,j,k+1)",&
      dz(k),dzh(k-1),dzh(k),rhot(i,j,k-1),rhot(i,j,k),rhot(i,j,k+1)
      close(88)
      call pariserror("A7 tiny outside of solid. Debug me")
    endif
  endif
enddo; enddo; enddo

if(check_setup) call check_poisson_setup(A,pmask,umask,vmask,wmask)

! End debugging and checking
end subroutine check_and_debug_Poisson

!=================================================================================================
!=================================================================================================
! The equation for the U velocity is setup with matrix A as
! A7*Uijk = A1*Ui-1jk + A2*Ui+1jk + A3*Uij-1k + 
!           A4*Uij+1k + A5*Uijk-1 + A6*Uijk+1 + A8
!-------------------------------------------------------------------------------------------------
subroutine SetupUvel(u,du,rho,mu,rho1,mu1,dt,A)
  use module_grid
  
  use module_BC

  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: u,du,rho,mu
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(out) :: A
  real(8), intent(in) :: dt,rho1,mu1
  real(8) :: rhom
  integer :: i,j,k
  if(TwoPhase) then
     do k=ks,ke; do j=js,je; do i=is,ie;
        rhom = 0.5d0*(rho(i+1,j,k)+rho(i,j,k))
        A(i,j,k,1) = dt/(dx(i  )*dxh(i)*rhom)*2d0*mu(i  ,j,k)
        A(i,j,k,2) = dt/(dx(i+1)*dxh(i)*rhom)*2d0*mu(i+1,j,k)
        A(i,j,k,3) = dt/(dy(j)*dyh(j-1)*rhom)*0.25d0*(mu(i,j,k)+mu(i+1,j,k)+mu(i,j-1,k)+mu(i+1,j-1,k))
        A(i,j,k,4) = dt/(dy(j)*dyh(j  )*rhom)*0.25d0*(mu(i,j,k)+mu(i+1,j,k)+mu(i,j+1,k)+mu(i+1,j+1,k))
        A(i,j,k,5) = dt/(dz(k)*dzh(k-1)*rhom)*0.25d0*(mu(i,j,k)+mu(i+1,j,k)+mu(i,j,k-1)+mu(i+1,j,k-1))
        A(i,j,k,6) = dt/(dz(k)*dzh(k  )*rhom)*0.25d0*(mu(i,j,k)+mu(i+1,j,k)+mu(i,j,k+1)+mu(i+1,j,k+1))
        A(i,j,k,7) = 1d0+sum(A(i,j,k,1:6))
        A(i,j,k,8) = u(i,j,k) + dt*du(i,j,k)
     enddo; enddo; enddo
  else
     do k=ks,ke; do j=js,je; do i=is,ie;
        rhom = rho1
        A(i,j,k,1) = dt/(dx(i  )*dxh(i)*rhom)*mu1
        A(i,j,k,2) = dt/(dx(i+1)*dxh(i)*rhom)*mu1
        A(i,j,k,3) = dt/(dy(j)*dyh(j-1)*rhom)*mu1
        A(i,j,k,4) = dt/(dy(j)*dyh(j  )*rhom)*mu1
        A(i,j,k,5) = dt/(dz(k)*dzh(k-1)*rhom)*mu1
        A(i,j,k,6) = dt/(dz(k)*dzh(k  )*rhom)*mu1
        A(i,j,k,7) = 1d0+sum(A(i,j,k,1:6))
        A(i,j,k,8) = u(i,j,k) + dt*du(i,j,k)
     enddo; enddo; enddo
  endif

!-------------------------------------------------------------------------------------------------
  !wall boundary conditions
  if(bdry_cond(2)==0 .and. coords(2)==0    ) then
    do k=ks,ke; do i=is,ie;
      A(i,js,k,7) = A(i,js,k,7) + A(i,js,k,3)
      A(i,js,k,3) = 0d0
    enddo; enddo
  endif
  if(bdry_cond(5)==0 .and. coords(2)==nPy-1) then
    do k=ks,ke; do i=is,ie;
      A(i,je,k,7) = A(i,je,k,7) + A(i,je,k,4)
      A(i,je,k,4) = 0d0
    enddo; enddo
  endif

  if(bdry_cond(3)==0 .and. coords(3)==0    ) then
    do j=js,je; do i=is,ie;
      A(i,j,ks,7) = A(i,j,ks,7) + A(i,j,ks,5)
      A(i,j,ks,5) = 0d0
    enddo; enddo
  endif
  if(bdry_cond(6)==0 .and. coords(3)==nPz-1) then
    do j=js,je; do i=is,ie;
      A(i,j,ke,7) = A(i,j,ke,7) + A(i,j,ke,6)
      A(i,j,ke,6) = 0d0
    enddo; enddo
  endif
end subroutine SetupUvel
!=================================================================================================
!=================================================================================================
subroutine SetupVvel(v,dv,rho,mu,rho1,mu1,dt,A) 
  use module_grid
  
  use module_BC
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: v,dv,rho,mu
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(out) :: A
  real(8), intent(in) :: dt,rho1,mu1
  real(8) :: rhom
  integer :: i,j,k
  if(TwoPhase) then
     do k=ks,ke; do j=js,je; do i=is,ie;
        rhom = 0.5d0*(rho(i,j+1,k)+rho(i,j,k))
        if((dy(j  )*dyh(j)*rhom)==0d0)then
           write(*,'(10I4,3f7.3)') rank,i,j,k,is,ie,js,je,ks,ke,dy(j),dyh(j),rhom
        endif
        A(i,j,k,3) = dt/(dy(j  )*dyh(j)*rhom)*2d0*mu(i,j  ,k)
        A(i,j,k,4) = dt/(dy(j+1)*dyh(j)*rhom)*2d0*mu(i,j+1,k)
        A(i,j,k,5) = dt/(dz(k)*dzh(k-1)*rhom)*0.25d0*(mu(i,j,k)+mu(i,j+1,k)+mu(i,j,k-1)+mu(i,j+1,k-1))
        A(i,j,k,6) = dt/(dz(k)*dzh(k  )*rhom)*0.25d0*(mu(i,j,k)+mu(i,j+1,k)+mu(i,j,k+1)+mu(i,j+1,k+1))
        A(i,j,k,1) = dt/(dx(i)*dxh(i-1)*rhom)*0.25d0*(mu(i,j,k)+mu(i,j+1,k)+mu(i-1,j,k)+mu(i-1,j+1,k))
        A(i,j,k,2) = dt/(dx(i)*dxh(i  )*rhom)*0.25d0*(mu(i,j,k)+mu(i,j+1,k)+mu(i+1,j,k)+mu(i+1,j+1,k))
        A(i,j,k,7) = 1d0+sum(A(i,j,k,1:6))
        A(i,j,k,8) = v(i,j,k) + dt*dv(i,j,k)
     enddo; enddo; enddo
  else
     do k=ks,ke; do j=js,je; do i=is,ie;
        rhom = rho1
        A(i,j,k,3) = dt/(dy(j  )*dyh(j)*rhom)*mu1
        A(i,j,k,4) = dt/(dy(j+1)*dyh(j)*rhom)*mu1
        A(i,j,k,5) = dt/(dz(k)*dzh(k-1)*rhom)*mu1
        A(i,j,k,6) = dt/(dz(k)*dzh(k  )*rhom)*mu1
        A(i,j,k,1) = dt/(dx(i-1)*dxh(i)*rhom)*mu1
        A(i,j,k,2) = dt/(dx(i)*dxh(i)*rhom)*mu1
        A(i,j,k,7) = 1d0+sum(A(i,j,k,1:6))
        A(i,j,k,8) = v(i,j,k) + dt*dv(i,j,k)
     enddo; enddo; enddo
  endif
  !-------------------------------------------------------------------------------------------------
  !wall boundary conditions
  if(bdry_cond(1)==0 .and. coords(1)==0    ) then
     do k=ks,ke; do j=js,je;
        A(is,j,k,7) = A(is,j,k,7) + A(is,j,k,1)
      A(is,j,k,1) = 0d0
    enddo; enddo
  endif
  if(bdry_cond(4)==0 .and. coords(1)==nPx-1) then
    do k=ks,ke; do j=js,je;
      A(ie,j,k,7) = A(ie,j,k,7) + A(ie,j,k,2)
      A(ie,j,k,2) = 0d0
    enddo; enddo
  endif

  if(bdry_cond(3)==0 .and. coords(3)==0    ) then
    do j=js,je; do i=is,ie;
      A(i,j,ks,7) = A(i,j,ks,7) + A(i,j,ks,5)
      A(i,j,ks,5) = 0d0
    enddo; enddo
  endif
  if(bdry_cond(6)==0 .and. coords(3)==nPz-1) then
    do j=js,je; do i=is,ie;
      A(i,j,ke,7) = A(i,j,ke,7) + A(i,j,ke,6)
      A(i,j,ke,6) = 0d0
    enddo; enddo
  endif
end subroutine SetupVvel
!=================================================================================================
!=================================================================================================
subroutine SetupWvel(w,dw,rho,mu,rho1,mu1,dt,A)
  use module_grid
  
  use module_BC
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: w,dw,rho,mu
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(out) :: A
  real(8), intent(in) :: dt,rho1,mu1
  real(8) :: rhom
  integer :: i,j,k
  if(TwoPhase) then
     do k=ks,ke; do j=js,je; do i=is,ie;
        rhom = 0.5d0*(rho(i,j,k+1)+rho(i,j,k))
        A(i,j,k,5) = dt/(dz(k  )*dzh(k)*rhom)*2d0*mu(i,j,k  )
        A(i,j,k,6) = dt/(dz(k+1)*dzh(k)*rhom)*2d0*mu(i,j,k+1)
        A(i,j,k,1) = dt/(dx(i)*dxh(i-1)*rhom)*0.25d0*(mu(i,j,k)+mu(i,j,k+1)+mu(i-1,j,k)+mu(i-1,j,k+1))
        A(i,j,k,2) = dt/(dx(i)*dxh(i  )*rhom)*0.25d0*(mu(i,j,k)+mu(i,j,k+1)+mu(i+1,j,k)+mu(i+1,j,k+1))
        A(i,j,k,3) = dt/(dy(j)*dyh(j-1)*rhom)*0.25d0*(mu(i,j,k)+mu(i,j,k+1)+mu(i,j-1,k)+mu(i,j-1,k+1))
        A(i,j,k,4) = dt/(dy(j)*dyh(j  )*rhom)*0.25d0*(mu(i,j,k)+mu(i,j,k+1)+mu(i,j+1,k)+mu(i,j+1,k+1))
        A(i,j,k,7) = 1d0+sum(A(i,j,k,1:6))
        A(i,j,k,8) = w(i,j,k) + dt*dw(i,j,k)
     enddo; enddo; enddo
  else
      do k=ks,ke; do j=js,je; do i=is,ie;
        rhom = rho1
        A(i,j,k,5) = dt/(dz(k  )*dzh(k)*rhom)*mu1
        A(i,j,k,6) = dt/(dz(k+1)*dzh(k)*rhom)*mu1
        A(i,j,k,1) = dt/(dx(i)*dxh(i-1)*rhom)*mu1
        A(i,j,k,2) = dt/(dx(i)*dxh(i  )*rhom)*mu1
        A(i,j,k,3) = dt/(dy(j)*dyh(j-1)*rhom)*mu1
        A(i,j,k,4) = dt/(dy(j)*dyh(j  )*rhom)*mu1
        A(i,j,k,7) = 1d0+sum(A(i,j,k,1:6))
        A(i,j,k,8) = w(i,j,k) + dt*dw(i,j,k)
     enddo; enddo; enddo
  endif
  !-------------------------------------------------------------------------------------------------
  !wall boundary conditions
  if(bdry_cond(1)==0 .and. coords(1)==0    ) then
     do k=ks,ke; do j=js,je;
        A(is,j,k,7) = A(is,j,k,7) + A(is,j,k,1)
        A(is,j,k,1) = 0d0
     enddo; enddo
  endif
  if(bdry_cond(4)==0 .and. coords(1)==nPx-1) then
    do k=ks,ke; do j=js,je;
      A(ie,j,k,7) = A(ie,j,k,7) + A(ie,j,k,2)
      A(ie,j,k,2) = 0d0
    enddo; enddo
  endif

  if(bdry_cond(2)==0 .and. coords(2)==0    ) then
    do k=ks,ke; do i=is,ie;
      A(i,js,k,7) = A(i,js,k,7) + A(i,js,k,3)
      A(i,js,k,3) = 0d0
    enddo; enddo
  endif
  if(bdry_cond(5)==0 .and. coords(2)==nPy-1) then
    do k=ks,ke; do i=is,ie;
      A(i,je,k,7) = A(i,je,k,7) + A(i,je,k,4)
      A(i,je,k,4) = 0d0
    enddo; enddo
  endif
end subroutine SetupWvel
!=================================================================================================
!=================================================================================================
! Solves the following linear equiation:
! A7*Pijk = A1*Pi-1jk + A2*Pi+1jk + A3*Pij-1k + 
!           A4*Pij+1k + A5*Pijk-1 + A6*Pijk+1 + A8
!-------------------------------------------------------------------------------------------------
subroutine LinearSolver(A,p,maxError,beta,maxit,it,ierr)
  use module_grid
  
  use module_BC
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: p
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(in) :: A
  real(8), intent(in) :: beta, maxError
  integer, intent(in) :: maxit
  integer, intent(out) :: it, ierr
  real(8) :: res, totalres
  integer :: i,j,k
  integer :: req(12),sta(MPI_STATUS_SIZE,12)
  logical :: mask(imin:imax,jmin:jmax,kmin:kmax)
!--------------------------------------ITERATION LOOP--------------------------------------------  
  do it=1,maxit
    do k=ks,ke; do j=js,je; do i=is,ie
      p(i,j,k)=(1d0-beta)*p(i,j,k)+beta* 1d0/A(i,j,k,7)*(              &
        A(i,j,k,1) * p(i-1,j,k) + A(i,j,k,2) * p(i+1,j,k) +            &
        A(i,j,k,3) * p(i,j-1,k) + A(i,j,k,4) * p(i,j+1,k) +            &
        A(i,j,k,5) * p(i,j,k-1) + A(i,j,k,6) * p(i,j,k+1) + A(i,j,k,8))
    enddo; enddo; enddo
!---------------------------------CHECK FOR CONVERGENCE-------------------------------------------
    res = 0d0
    call ghost_x(p,1,req( 1: 4)); call ghost_y(p,1,req( 5: 8)); call ghost_z(p,1,req( 9:12))
    do k=ks+1,ke-1; do j=js+1,je-1; do i=is+1,ie-1
      res=res+abs(-p(i,j,k) * A(i,j,k,7) +                             &
        A(i,j,k,1) * p(i-1,j,k) + A(i,j,k,2) * p(i+1,j,k) +            &
        A(i,j,k,3) * p(i,j-1,k) + A(i,j,k,4) * p(i,j+1,k) +            &
        A(i,j,k,5) * p(i,j,k-1) + A(i,j,k,6) * p(i,j,k+1) + A(i,j,k,8) )**2
    enddo; enddo; enddo
    call MPI_WAITALL(12,req,sta,ierr)
    mask=.true.
    mask(is+1:ie-1,js+1:je-1,ks+1:ke-1)=.false.
    do k=ks,ke; do j=js,je; do i=is,ie
      if(mask(i,j,k))res=res+abs(-p(i,j,k) * A(i,j,k,7) +              &
        A(i,j,k,1) * p(i-1,j,k) + A(i,j,k,2) * p(i+1,j,k) +            &
        A(i,j,k,3) * p(i,j-1,k) + A(i,j,k,4) * p(i,j+1,k) +            &
        A(i,j,k,5) * p(i,j,k-1) + A(i,j,k,6) * p(i,j,k+1) + A(i,j,k,8) )**2
    enddo; enddo; enddo
    res = res/dble(Nx*Ny*Nz)
    if ( (res*npx*npy*npz>1.d16) .or. (res.ne.res) ) then
      print*,'Pressure residual is too large or NaN, res= ',res 
      print*,'Pressure solver diverged after',it,'iterations at rank ',rank
      stop  !return
    else 
      call MPI_ALLREDUCE(res, totalres, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr)
      totalres=sqrt(totalres)
      if (totalres<maxError) exit
    end if !res
  enddo
  if(it==maxit+1 .and. rank==0) write(*,*) 'Warning: Pressure Solver reached maxit: totalres',totalres,maxError
end subroutine LinearSolver
!=================================================================================================
!=================================================================================================
! Solves the following linear equiation:
! A7*Uijk = umask*(A1*Ui-1jk + A2*Ui+1jk + A3*Uij-1k + 
!           A4*Uij+1k + A5*Uijk-1 + A6*Uijk+1 + A8)
!-------------------------------------------------------------------------------------------------
subroutine LinearSolver1(A,u,umask,maxError,beta,maxit,it,ierr)
  use module_grid
  
  use module_BC
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: u
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: umask
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(in) :: A
  real(8), intent(in) :: beta, maxError
  integer, intent(in) :: maxit
  integer, intent(out) :: it, ierr
  real(8) :: res, totalres
  integer :: i,j,k
  integer :: req(12),sta(MPI_STATUS_SIZE,12)
  logical :: mask(imin:imax,jmin:jmax,kmin:kmax)
!--------------------------------------ITERATION LOOP--------------------------------------------  
  do it=1,maxit
    do k=ks,ke; do j=js,je; do i=is,ie
      u(i,j,k)=umask(i,j,k)*((1d0-beta)*u(i,j,k)+beta* 1d0/A(i,j,k,7)*(              &
        A(i,j,k,1) * u(i-1,j,k) + A(i,j,k,2) * u(i+1,j,k) +            &
        A(i,j,k,3) * u(i,j-1,k) + A(i,j,k,4) * u(i,j+1,k) +            &
        A(i,j,k,5) * u(i,j,k-1) + A(i,j,k,6) * u(i,j,k+1) + A(i,j,k,8)))
    enddo; enddo; enddo
!---------------------------------CHECK FOR CONVERGENCE-------------------------------------------
    res = 0d0
    call ghost_x(u,1,req( 1: 4)); call ghost_y(u,1,req( 5: 8)); call ghost_z(u,1,req( 9:12))
    do k=ks+1,ke-1; do j=js+1,je-1; do i=is+1,ie-1
      res=res+ umask(i,j,k)*abs(-u(i,j,k) * A(i,j,k,7) +                             &
        A(i,j,k,1) * u(i-1,j,k) + A(i,j,k,2) * u(i+1,j,k) +            &
        A(i,j,k,3) * u(i,j-1,k) + A(i,j,k,4) * u(i,j+1,k) +            &
        A(i,j,k,5) * u(i,j,k-1) + A(i,j,k,6) * u(i,j,k+1) + A(i,j,k,8) )**2
    enddo; enddo; enddo
    call MPI_WAITALL(12,req,sta,ierr)
    mask=.true.
    mask(is+1:ie-1,js+1:je-1,ks+1:ke-1)=.false.
    do k=ks,ke; do j=js,je; do i=is,ie
      if(mask(i,j,k))res=res+umask(i,j,k)*abs(-u(i,j,k) * A(i,j,k,7) +              &
        A(i,j,k,1) * u(i-1,j,k) + A(i,j,k,2) * u(i+1,j,k) +            &
        A(i,j,k,3) * u(i,j-1,k) + A(i,j,k,4) * u(i,j+1,k) +            &
        A(i,j,k,5) * u(i,j,k-1) + A(i,j,k,6) * u(i,j,k+1) + A(i,j,k,8) )**2
    enddo; enddo; enddo
    res = res/dble(Nx*Ny*Nz)
    if ( (res*npx*npy*npz>1.d16) .or. (res.ne.res) ) then
      print*,'Viscous term residual is too large or NaN, res= ',res 
      print*,'Viscous term solver diverged after',it,'iterations at rank ',rank
      stop  !return
    else 
      call MPI_ALLREDUCE(res, totalres, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr)
      totalres=sqrt(totalres)
      if (totalres<maxError) exit
    end if !res
  enddo
  if(it==maxit+1 .and. rank==0) write(*,*) 'Warning: Viscous Solver reached maxit: totalres',totalres,maxError
end subroutine LinearSolver1
!=================================================================================================
!=================================================================================================
! Returns the residual. L1 norm of the new divergence. 
!-------------------------------------------------------------------------------------------------
subroutine calcResidual(A,p,NormOrder, residual)
  use module_grid
  use module_BC
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: p
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(in) :: A
  integer, intent(in) :: NormOrder
  real(8) :: res, Residual,locres
  integer :: i,j,k, ierr
  res = 0d0
  do k=ks,ke; do j=js,je; do i=is,ie
      locres = abs(-p(i,j,k) * A(i,j,k,7) +                             &
      A(i,j,k,1) * p(i-1,j,k) + A(i,j,k,2) * p(i+1,j,k) +            &
      A(i,j,k,3) * p(i,j-1,k) + A(i,j,k,4) * p(i,j+1,k) +            &
      A(i,j,k,5) * p(i,j,k-1) + A(i,j,k,6) * p(i,j,k+1) + A(i,j,k,8) )
      if ( NormOrder == 1 ) then 
         res = res + locres
      else if ( NormOrder == 2 ) then 
         res = res + locres*locres
      else
         res = max(res,locres)
      end if ! NormOrder
  enddo; enddo; enddo
  if ( NormOrder == 1 ) then 
      call MPI_ALLREDUCE(res, residual, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr)
      residual = residual/dble(Nx*Ny*Nz)
  else if ( NormOrder == 2 ) then 
      call MPI_ALLREDUCE(res, residual, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr)
      Residual = sqrt(residual)/dble(Nx*Ny*Nz)
  else 
      call MPI_ALLREDUCE(res, residual, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_Comm_Cart, ierr)
  end if ! NormOrder
end subroutine calcResidual
!=================================================================================================
!=================================================================================================
! Returns the residual including pmask  ! @@fixme: dangerous pmask not used that way 
!-------------------------------------------------------------------------------------------------
subroutine calcResidual1(A,p,pmask,Residual)
  use module_grid
  use module_BC
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: p
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: pmask
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(in) :: A
  real(8) :: res, totalres, Residual
  integer :: i,j,k, ierr
  res = 0d0
  do k=ks,ke; do j=js,je; do i=is,ie
    res=res+pmask(i,j,k)*abs(-p(i,j,k) * A(i,j,k,7) +                             &
      A(i,j,k,1) * p(i-1,j,k) + A(i,j,k,2) * p(i+1,j,k) +            &
      A(i,j,k,3) * p(i,j-1,k) + A(i,j,k,4) * p(i,j+1,k) +            &
      A(i,j,k,5) * p(i,j,k-1) + A(i,j,k,6) * p(i,j,k+1) + A(i,j,k,8) )**2
  enddo; enddo; enddo
  res = res/float(Nx*Ny*Nz)
  call MPI_ALLREDUCE(res, totalres, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr)
  Residual = sqrt(totalres)
end subroutine calcResidual1
!=================================================================================================
!=================================================================================================
! Returns the flow rate
!-------------------------------------------------------------------------------------------------
subroutine calcsum(p)
  use module_grid
  use module_BC
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: p
  real(8) :: flow
  call calcsum_shift(p,flow,0,0,0)
end subroutine calcsum

subroutine calcsum2(p,flowrate)
  use module_grid
  use module_BC
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: p
  real(8) :: flow,flowrate
  call calcsum_shift(p,flow,1,0,0)
  call calcsum_shift(p,flow,0,1,0)
  call calcsum_shift(p,flow,0,0,1)
  call calcsum_shift(p,flow,-1,0,0)
  call calcsum_shift(p,flow,0,-1,0)
  call calcsum_shift(p,flow,0,0,-1)
  call calcsum_shift(p,flow,0,0,0)
  flowrate = -1.
end subroutine calcsum2

subroutine calcsum_shift(p,porosity,sx,sy,sz)
  use module_grid
  use module_BC
  use module_IO
  implicit none
  include 'mpif.h'
  integer, intent(in) :: sx,sy,sz
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: p
  real(8) :: flow, totalflow
  real(8), intent(out) :: porosity
  real(8) volume
  integer :: i,j,k, ierr
  flow=0.d0
  volume=Nx*Ny*Nz
  do k=ks,ke; do j=js,je; do i=is,ie;
     flow=flow+p(i+sx,j+sy,k+sz)
  enddo;enddo;enddo
  call MPI_REDUCE(flow, totalflow, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_Domain, ierr)
  porosity=totalflow/volume
  if(rank==0) then 
     open(UNIT=89,file=trim(out_path)//'/porosity.txt')
     write(89,310)  porosity
     close(unit=89)
  endif
  310 format(F17.11)
end subroutine calcsum_shift
!-------------------------------------------------------------------------------------------------

function calc_imax(f)
  use module_grid
  implicit none
  include 'mpif.h'
  integer calc_imax,ierr,imax_loc,imax_all
  integer, dimension(imin:imax,jmin:jmax,kmin:kmax),  intent(in) :: f
  integer :: i,j,k
  imax_loc=-2147483647
  do k=kmin,kmax
     do j=jmin,jmax
        do i=imin,imax
           imax_loc = max(imax_loc,f(i,j,k))
        enddo
     enddo
  enddo
  imax_all = 0
  call MPI_ALLREDUCE(imax_loc, imax_all, 1, MPI_INTEGER, MPI_MAX, MPI_Comm_Cart, ierr)
  calc_imax=imax_all
  end function calc_imax
 
  subroutine THRESHOLD(z)
    implicit none
    real(8),intent(inout) :: z
    if ( z < 0.d0 ) then
       z = 0.d0
    else if ( z > 1.d0 ) then
       z = 1.d0
    end if
  end subroutine THRESHOLD
