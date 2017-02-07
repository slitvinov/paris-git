!=================================================================================================
!=================================================================================================
! PARIS  Parallel Robust Interface Simulator 
!=================================================================================================
! module_VOF: Contains definition of variables for the Volume of Fluid interface tracking.
!-------------------------------------------------------------------------------------------------
!
! Contact: Stephane Zaleski zaleski@dalembert.upmc.fr
! 
! Authors (in alphabetical order):
!         Daniel Fuster      
!         Ruben Scardovelli  
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
!=================================================================================================
module module_VOF 
  use module_grid
  use module_IO
  use module_tmpvar
  use module_BC
  implicit none
  real(8), dimension(:,:,:), allocatable, target :: cvof ! VOF tracer variable
  real(8), dimension(:,:,:), allocatable :: cvofold 
  integer, dimension(:,:,:), allocatable :: vof_flag, tmp_flag ! 
  !   0 empty
  !   1 full
  !   2 fractional
  !   3 unknown
  integer, dimension(:,:,:), allocatable :: vof_phase
  !   0 cvof <= 0.5 liquid phase in free surface
  !   1 cvof > 0.5 gas phase in free surface

  real(8), parameter  :: A_h = 2d0  ! For initialisation of height test
  real(8), parameter  :: TINY = 1d-50

  real(8) :: EPSC  = 1.d-12  ! for clipping vof and setting flags

  character(20) :: vofbdry_cond(6),test_type,vof_advect,cond
  integer :: parameters_read=0, refinement=-1 
  integer :: cylinder_dir = 0
  logical :: test_heights = .false.  
  logical :: test_capwave = .false.
  logical :: test_droplet = .false.
  logical :: test_shear = .false.  
  logical :: test_control_droplet = .false.
  logical :: normal_up = .true.    ! used for the h
  logical :: test_curvature = .false.  
  logical :: test_curvature_2D = .false.  
  logical :: test_HF = .false.
  logical :: test_LP = .false.
  logical :: test_tag = .false.
  logical :: test_D2P = .false.
  logical :: test_bubbles = .false.
  logical :: test_injectdrop = .false.
  logical :: test_jet = .false.
  logical :: test_cylinder_advection = .false.
  logical :: test_parasitic_currents = .false.
  logical :: test_shear_multiphase = .false.
  logical :: test_KHI2D = .false.
  logical :: test_PhaseInversion = .false.
  logical :: test_fsdrop = .false.
  logical :: test_randombubs = .false.
  logical :: test_plane = .false.
  logical :: test_risingbubble = .false.
  logical :: test_lattice_bubs = .false.
  logical :: test_fcc_lattice = .false.
  logical :: test_turb = .false.
  logical :: linfunc_initialized = .false.
  logical :: DoMOMCONS = .false.
  logical :: LinInterp = .false.
  logical :: STGhost = .false.
  logical :: use_Vofi
!   logical :: oldvof
  logical :: do_rotation
  logical :: debug_curvature
  logical :: mixed_heights
  logical :: use_full_heights
  logical :: debug_par
  logical :: staggeredmom =.true.

  real(8) :: b1,b2,b3,b4,hshift
  integer :: nfilter

  logical :: DoLPP = .false.
  logical :: output_filtered_VOF = .false.

  integer, parameter :: ArithMean = 101
  integer, parameter :: HarmMean  = 102
  integer :: ViscMean,DensMean

  real(8), parameter :: PI = 3.14159265359d0

  real(8) :: tot_clean, cleaned

  real(8) :: cwg_threshold
  logical :: do_clean_debris
  integer :: clean_debris_method, nsteps_clean_debris, clean_debris_neighbours

  logical :: filter_random_seeds

contains
!=================================================================================================
!=================================================================================================
!__|__|__|__|
!__|__|_2.5_|    Maximum possible filter value
!__|__|_1.5_|
!__|__|_x|__|
!__|1.9__|__|
!2.9__|__|__|
!
  subroutine initialize_linfunc()
    implicit none
! Gerris method
! Each vertex averages 8 cells  --> 1/8
! Each Cell averages 8 vertices --> 1/8
! 8 vertex cells belong to only one vertex ->1/64
! 12 edge cells belong to 2 vertices -> 2/64
! 6 face cells belong to 4 vertices -> 4/64
! 8/4 + 6*4/64 + 12*2/64 + 8/64 = 1
    b1=8D0/64
    b2=4D0/64
    b3=2D0/64
    b4=1D0/64
    linfunc_initialized = .true.
  end subroutine initialize_linfunc
!------------------------------------------------------------------------
  subroutine filter(field)
    use module_tmpvar
    implicit none
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(out) :: field
    integer :: i,j,k,n
    field=cvof
    do n=1,nfilter
       call do_all_ghost(field)
       do k=ks-1,ke+1; do j=js-1,je+1; do i=is-1,ie+1
          tmp(i,j,k) = b1*field(i,j,k) + & 
               b2*( field(i-1,j,k) + field(i,j-1,k) + field(i,j,k-1) + &
                    field(i+1,j,k) + field(i,j+1,k) + field(i,j,k+1) ) + &
               b3*( field(i+1,j+1,k) + field(i+1,j-1,k) + field(i-1,j+1,k) + field(i-1,j-1,k) + &
                    field(i+1,j,k+1) + field(i+1,j,k-1) + field(i-1,j,k+1) + field(i-1,j,k-1) + &
                    field(i,j+1,k+1) + field(i,j+1,k-1) + field(i,j-1,k+1) + field(i,j-1,k-1) ) + &
               b4*( field(i+1,j+1,k+1) + field(i+1,j+1,k-1) + field(i+1,j-1,k+1) + field(i+1,j-1,k-1) +  &
                    field(i-1,j+1,k+1) + field(i-1,j+1,k-1) + field(i-1,j-1,k+1) + field(i-1,j-1,k-1) )
       enddo; enddo; enddo
       field=tmp
    enddo
  end subroutine filter
!------------------------------------------------------------------------
  subroutine linfunc(field,a1,a2,MeanFlag)
    implicit none
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(out) :: field
    real(8) :: cfiltered
    real(8), intent(in) :: a1,a2
    integer, intent(in) :: MeanFlag
    integer :: i,j,k
    real(8) :: inva1=0d0,inva2=0d0

    if(.not.linfunc_initialized) call initialize_linfunc

    if ( MeanFlag == HarmMean ) then 
      inva1 = 1.d0/(a1+1.d-50)
      inva2 = 1.d0/(a2+1.d-50)
    end if ! MeanFlag 

    if(nfilter==0) then
       if ( MeanFlag == ArithMean ) then  
         field = cvof*(a2-a1)+a1
       else if ( MeanFlag == HarmMean ) then 
         field = 1.d0/(cvof*(inva2-inva1)+inva1)
       end if ! MeanFlag
    else if (nfilter==1) then
       do k=kmin,kmax; do j=jmin,jmax; do i=imin,imax
         if ( k>kmin.and.k<kmax.and.j>jmin.and.j<jmax.and.i>imin.and.i<imax) then 
            cfiltered = b1*cvof(i,j,k) + & 
               b2*( cvof(i-1,j,k) + cvof(i,j-1,k) + cvof(i,j,k-1) + &
                    cvof(i+1,j,k) + cvof(i,j+1,k) + cvof(i,j,k+1) ) + &
               b3*( cvof(i+1,j+1,k) + cvof(i+1,j-1,k) + cvof(i-1,j+1,k) + cvof(i-1,j-1,k) + &
                    cvof(i+1,j,k+1) + cvof(i+1,j,k-1) + cvof(i-1,j,k+1) + cvof(i-1,j,k-1) + &
                    cvof(i,j+1,k+1) + cvof(i,j+1,k-1) + cvof(i,j-1,k+1) + cvof(i,j-1,k-1) ) + &
               b4*( cvof(i+1,j+1,k+1) + cvof(i+1,j+1,k-1) + cvof(i+1,j-1,k+1) + cvof(i+1,j-1,k-1) +  &
                    cvof(i-1,j+1,k+1) + cvof(i-1,j+1,k-1) + cvof(i-1,j-1,k+1) + cvof(i-1,j-1,k-1) )
         else 
            cfiltered = cvof(i,j,k)
         end if ! i,j,k

         if ( MeanFlag == ArithMean ) then  
            field(i,j,k) = cfiltered*(a2-a1)+a1
         else if ( MeanFlag == HarmMean ) then 
            field(i,j,k) = 1.d0/(cfiltered*(inva2-inva1)+inva1)
         end if ! MeanFlag
       enddo; enddo; enddo
    endif
    
    call do_all_ghost(field)

  end subroutine linfunc

  subroutine linfunclocal(field,a1,a2,MeanFlag,i,j,k)
    implicit none
    real(8), intent(out) :: field
    real(8) :: cfiltered
    real(8), intent(in) :: a1,a2
    integer, intent(in) :: MeanFlag
    integer :: i,j,k
    real(8) :: inva1=0d0,inva2=0d0

    if(.not.linfunc_initialized) call initialize_linfunc

    if ( MeanFlag == HarmMean ) then 
      inva1 = 1.d0/(a1+1.d-50)
      inva2 = 1.d0/(a2+1.d-50)
    end if ! MeanFlag 

    if(nfilter==0) then
       if ( MeanFlag == ArithMean ) then  
         field = cvof(i,j,k)*(a2-a1)+a1
       else if ( MeanFlag == HarmMean ) then 
         field = 1.d0/(cvof(i,j,k)*(inva2-inva1)+inva1)
       end if ! MeanFlag
    else if (nfilter==1) then
          cfiltered = b1*cvof(i,j,k) + & 
               b2*( cvof(i-1,j,k) + cvof(i,j-1,k) + cvof(i,j,k-1) + &
                    cvof(i+1,j,k) + cvof(i,j+1,k) + cvof(i,j,k+1) ) + &
               b3*( cvof(i+1,j+1,k) + cvof(i+1,j-1,k) + cvof(i-1,j+1,k) + cvof(i-1,j-1,k) + &
                    cvof(i+1,j,k+1) + cvof(i+1,j,k-1) + cvof(i-1,j,k+1) + cvof(i-1,j,k-1) + &
                    cvof(i,j+1,k+1) + cvof(i,j+1,k-1) + cvof(i,j-1,k+1) + cvof(i,j-1,k-1) ) + &
               b4*( cvof(i+1,j+1,k+1) + cvof(i+1,j+1,k-1) + cvof(i+1,j-1,k+1) + cvof(i+1,j-1,k-1) +  &
                    cvof(i-1,j+1,k+1) + cvof(i-1,j+1,k-1) + cvof(i-1,j-1,k+1) + cvof(i-1,j-1,k-1) )
         if ( MeanFlag == ArithMean ) then  
            field = cfiltered*(a2-a1)+a1
         else if ( MeanFlag == HarmMean ) then 
            field = 1.d0/(cfiltered*(inva2-inva1)+inva1)
         end if ! MeanFlag
    endif
  end subroutine linfunclocal

!=================================================================================================

  subroutine ReadVOFParameters
    use module_flow
    use module_BC
    use module_2phase
    use module_freesurface
    use module_IO
    implicit none
    include 'mpif.h'
    integer ierr,in,i
    logical file_is_there, fs_file_is_there
    logical ViscMeanIsArith, DensMeanIsArith
    namelist /vofparameters/ vofbdry_cond,test_type,VOF_advect,refinement, &
         cylinder_dir, normal_up, DoLPP, &
         FreeSurface, ViscMeanIsArith, DensMeanIsArith, &
         output_filtered_VOF, DoMOMCONS, use_vofi,nfilter, &
         hshift, do_rotation, debug_curvature, mixed_heights, &
         use_full_heights, debug_par, STGhost, &
         r_min, var_r, coord_min, var_coord, var_pos,&
         out_centroid, do_clean_debris,&
         clean_debris_method,nsteps_clean_debris,clean_debris_neighbours,&
         filter_random_seeds, nb, n_cell,&
         cwg_threshold,EPSC,staggeredmom, LinInterp
      
    namelist /FSparameters/ X_level, RP_test, gamma, R_ref, P_ref,&
         VTK_OUT, NOUT_VTK, step_max, limit, curve_stats, order_extrap,&
         do_2nd_projection, check_stray_liquid, n_stray_liquid, &
         FS_HYPRE, fs_refactor, fs_restrictor
    
    vofbdry_cond=['undefined','undefined','undefined','undefined','undefined','undefined']
    test_type='droplet'
    VOF_advect='CIAM'
    refinement=-1 ! redundant
    cylinder_dir=0 ! redundant
    normal_up=.true. ! redundant
    DoLPP=.false.
    out_centroid=.false.
    FreeSurface=.false.
    STGhost=.false.
    ViscMeanIsArith=.true.; DensMeanIsArith=.true.
    output_filtered_VOF=.false. ! redundant
    DoMOMCONS = .false.
    use_vofi = .false.
    nfilter = 0
    hshift =0d0
    do_rotation = .false.
    debug_curvature=.false.
    ! oldvof = .false.  ! .true.
    mixed_heights = .true.
    use_full_heights = .true. 
    debug_par = .false.
    r_min=0.02; var_r=0.03; coord_min=0.15; var_coord=0.7; var_pos=0.0d0
    do_clean_debris = .false.
    clean_debris_method = 1 
    nsteps_clean_debris = 10 
    clean_debris_neighbours = 2
    filter_random_seeds = .false.
    nb = 3; n_cell=3
    cwg_threshold=0d0
    EPSC = 1.d-12
    staggeredmom = .true.

    in=31

    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    inquire(file='inputvof',exist=file_is_there)
    open(unit=in, file='inputvof', status='old', action='read', iostat=ierr)
    if (file_is_there) then
       if(ierr == 0) then
          read(UNIT=in,NML=vofparameters)
          if(rank==0) write(out,*)'VOF parameters read successfully'
       else
          print *, 'rank=',rank,' has error ',ierr,' opening file inputsolids'
       endif
    else
       if (rank == 0) call pariserror("ReadVOFParameters: no 'inputvof' file.")
    endif
    close(in)
     
    do i=1,3
       if(vofbdry_cond(i) == 'undefined') call pariserror("vofbdry_cond undefined")
       if(vofbdry_cond(i+3) == 'undefined') vofbdry_cond(i+3) = vofbdry_cond(i) 
    enddo
    if(refinement==-1) then
       refinement=8
       if(rank==0) write(*,*) "Using default value for refinement"
    else if(refinement==0) then
       if(rank==0) write(*,*) "Warning: no refinement of VOF initial condition"
    endif

    if(ViscMeanIsArith) then
       ViscMean = ArithMean
       if(rank==0) print *, "Using Arithmetic Mean for Viscosity"
    else
       ViscMean = HarmMean
       if(rank==0) print *, "Using Harmonic Mean for Viscosity"
    endif
 
    if(DensMeanIsArith) then
       DensMean = ArithMean
       if(rank==0) print *, "Using Arithmetic Mean for Density"
    else
       DensMean = HarmMean
       if(rank==0) print *, "Using Harmonic Mean for Density"
    endif

 
    if (rank == 0) then 
     !open(unit=out, file=trim(out_path)//'/output', action='write', iostat=ierr)
     !if (ierr .ne. 0) call pariserror("ReadParameters: error opening output file")
     write(UNIT=out,NML=vofparameters)
    end if ! rank

    ! Free Surface parameters to be read from a parameter file called "inputFS"
    if (FreeSurface) then
       X_level = 2; RP_test = .false.; gamma = 1.4d0; R_ref = 1.d0; P_ref = 1.d0; 
       VTK_OUT = .false.; NOUT_VTK = 100000; vtk_open = .false.
       step_max = 1
       fs_refactor = .false.
       fs_restrictor = .false.
       limit = 1.0d-2
       curve_stats = .false.
       order_extrap = 1
       do_2nd_projection = .false.
       check_stray_liquid = .false.
       n_stray_liquid = 50
       FS_Hypre = .false.
       p_ind = -2
       !call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
       inquire(file='inputFS',exist=fs_file_is_there)
       in = 77
       open(unit=in, file='inputFS', status='old', action='read', iostat=ierr)
       if (file_is_there) then
          if(ierr == 0) then
             read(UNIT=in,NML=FSparameters)
             if(rank==0) write(out,*)'Free Surface parameters read successfully'
          else
             print *, 'rank=',rank,' has error ',ierr,' opening file inputFS'
          endif
       else
          if (rank == 0) call pariserror("ReadVOFParameters: no 'inputFS' file and FreeSurface is on.")
       endif
       close(in)
       if (RP_test .and. FS_HYPRE) then
          if (rank==0) print *,'For FS RP test, FS_HYPRE is not available. Switching to standard Gauss-Seidel'
          FS_HYPRE=.false.
       endif
    endif

  end subroutine ReadVOFParameters
!
  subroutine initialize_VOF
    use module_grid
    use module_BC
    use module_flow
    use module_2phase
    use module_freesurface
    implicit none
    include 'mpif.h'
    integer :: ierr,dir,orientation,bdry
    character :: dc(3) = (/ "x","y","z" /)
    call ReadVOFParameters
! Check grid
    if(read_x.or.read_y.or.read_z) then
       if((xform.ne.0.d0).or.(yform.ne.0.d0).or.(zform.ne.0.d0)) then
          if (rank == 0) print *, "VOF does not yet work with variable grids"
          call MPI_Finalize(ierr)
          stop
       endif
    endif
    allocate(cvof(imin:imax,jmin:jmax,kmin:kmax),vof_flag(imin:imax,jmin:jmax,kmin:kmax),&
         vof_phase(imin:imax,jmin:jmax,kmin:kmax))
    if (DoMOMCONS) then
      allocate(momentum(imin:imax,jmin:jmax,kmin:kmax))
      allocate(tmp_flag(imin:imax,jmin:jmax,kmin:kmax))
    endif
    cvof = 0.d0
    vof_flag = 3
    vof_phase = 2
    if (out_centroid) opened_cent=.false.
    !allocate matrices for Free Surface
    if(FreeSurface) then
       allocate(u_cmask(imin:imax,jmin:jmax,kmin:kmax), v_cmask(imin:imax,jmin:jmax,kmin:kmax), &
            w_cmask(imin:imax,jmin:jmax,kmin:kmax), x_mod(imin:imax,jmin:jmax,kmin:kmax), &
            y_mod(imin:imax,jmin:jmax,kmin:kmax), z_mod(imin:imax,jmin:jmax,kmin:kmax), &
            P_gx(imin:imax,jmin:jmax,kmin:kmax), P_gy(imin:imax,jmin:jmax,kmin:kmax),&
            P_gz(imin:imax,jmin:jmax,kmin:kmax), pcmask(imin:imax,jmin:jmax,kmin:kmax), &
            P_gas(imin:imax,jmin:jmax,kmin:kmax),p_ext(imin:imax,jmin:jmax,kmin:kmax), &
            v_source(imin:imax,jmin:jmax,kmin:kmax), implode_flag(0:NumBubble*28))
       u_cmask = 3; v_cmask = 3; w_cmask = 3
       v_source = 0.d0
       x_mod = dxh((is+ie)/2); y_mod = dyh((js+je)/2); z_mod = dzh((ks+ke)/2)
       P_gx = 0.d0; P_gy = 0.d0; P_gz = 0.d0
       pcmask=3; p_ext=0.d0
       initialize_fs = .true.
       vtk_open = .false.
       implode_flag = .false.

       inflow=.true.
       do bdry=1,6
          if (.not.bdry_cond(bdry) == 3) inflow = .false.
       enddo
    endif
    if ( itime_scheme == 2 ) then  
      allocate(cvofold(imin:imax,jmin:jmax,kmin:kmax))
      cvofold = 0.d0
    end if ! itime_scheme
    if(test_type=='droplet') then
       test_droplet = .true.
    else if(test_type=='shear') then
       test_droplet = .true.
       test_shear = .true.
    else if(test_type=='control_droplet') then
       test_control_droplet = .true.
    else if(test_type=='height_test') then
       test_heights = .true.
    else if(test_type=='Curvature_test') then
       test_curvature = .true.
    else if(test_type=='Curvature2D') then
       test_curvature_2D = .true.
    else if(test_type=='Capillary_wave') then
       test_capwave = .true.    
    else if(test_type=='tag_test') then
       test_tag = .true.
    else if(test_type=='D2P_test') then
       test_D2P = .true.
    else if(test_type=='bubbles_test') then
       test_bubbles = .true.
    else if(test_type=='injectdrop_test') then
       test_injectdrop = .true.
    else if(test_type=='jet') then
       test_jet = .true.
    else if(test_type=='cylinder_advection') then
       test_cylinder_advection = .true.
    else if(test_type=='parasitic_currents') then
       test_parasitic_currents = .true.
    else if(test_type=='shear_multiphase') then
       test_shear_multiphase = .true.
    else if(test_type=='KHI2D') then
       test_KHI2D = .true.
    else if(test_type=='PhaseInversion') then
       test_PhaseInversion = .true.
    else if(test_type=='droplet_fs') then
       test_fsdrop = .true.
    else if(test_type=='random_bubs') then
       test_randombubs = .true.
    else if(test_type=='plane') then
       test_plane = .true.   
    else if(test_type=='rising_bubble') then
       test_risingbubble = .true.
    else if(test_type=='cubic_lat_bubs') then
       test_lattice_bubs = .true.
    else if(test_type=='fcc_lat_bubs') then
       test_fcc_lattice = .true.
    else if(test_type=='Iso_Turbulence') then 
       test_turb = .true.
    else
       write(*,*) test_type, rank
       call pariserror("unknown initialization")
    endif
    test_HF = test_heights .or. test_curvature .or. test_curvature_2D
    test_LP = test_tag .or. test_D2P 

! Post-read initialization of boundary conditions
    ! orientation order:    
    ! x- y- z- x+ y+ z+

    if(.not.bdry_read) call pariserror("bdry not read")
    do orientation=4,6
       dir = orientation-3
       if(vofbdry_cond(orientation)=='periodic'.or.vofbdry_cond(dir)=='periodic') then
          vofbdry_cond(orientation) = 'periodic'
          vofbdry_cond(dir) = 'periodic'
          if(bdry_cond(dir) /= 1) &
               call pariserror("periodic for both bdry_cond and vofbdry_cond or none at all")
       endif
    enddo

    if(rank==0) then
       do dir=1,3
          write(6,'(A1,"-: ",(A)," / ",(A),T32," ",A1,"+: ",(A)," / ",(A))') &
               dc(dir),TRIM(expl(bdry_cond(dir))),TRIM(vofbdry_cond(dir)), &
               dc(dir),TRIM(expl(bdry_cond(dir+3))),TRIM(vofbdry_cond(dir+3)) ! ,rank
       enddo
    endif
  end subroutine initialize_VOF
!=================================================================================================
!   a hack to get the flags quickly (temporary)
!=================================================================================================
  subroutine get_flags_and_clip(c,cflag)
    integer :: i,j,k
    real (8)  , dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: c
    integer   , dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: cflag
    if(ng.lt.2) call pariserror("wrong ng")
    do k=kmin,kmax
       do j=jmin,jmax
          do i=imin,imax
             if(c(i,j,k).le.EPSC) then
                cflag(i,j,k) = 0
                c(i,j,k)=0.d0
             else if(c(i,j,k).ge.1.d0-EPSC) then
                cflag(i,j,k) = 1
                c(i,j,k) = 1.d0
             else
                cflag(i,j,k) = 2
             endif
          enddo
       enddo
    enddo
  end subroutine get_flags_and_clip
!=================================================================================================
!  get the majority phase flag
!=================================================================================================
  subroutine get_vof_phase(c)
    use module_IO
    use module_flow
    integer :: i,j,k
    real (8)  , dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: c
    do k=kmin,kmax
       do j=jmin,jmax
          do i=imin,imax
             if (c(i,j,k)<(0.5d0-EPSC)) then
                vof_phase(i,j,k) = 0
             else
                vof_phase(i,j,k) = 1
             endif
          enddo
       enddo
    enddo
  end subroutine get_vof_phase
  !=================================================================================================
  !=================================================================================================
  !  Initialize vof field and flags
  !=================================================================================================
  subroutine initconditions_VOF()
    use module_grid
    use module_flow
    use module_BC
    use module_freesurface
    use module_2phase
    use module_solid

    implicit none
    include 'mpif.h'
    integer :: ierr
    integer , parameter :: ngh=2
    integer :: ipar
    integer, parameter :: root_rank = 0
    integer :: i,j,k
    real(8) :: ryz, sine, NozzleThickness
    
    if( test_D2P .or. ( test_tag .and. NumBubble >= 50 ) ) then 
       if ( rank == root_rank ) then
          call random_bubbles
       endif
       if (filter_random_seeds) then
          call MPI_BCAST(NumBubble, 1, MPI_INTEGER, &
               root_rank, MPI_Comm_Cart, ierr)
       endif
       call MPI_BCAST(rad, NumBubble, MPI_REAL8, &
                      root_rank, MPI_Comm_Cart, ierr)
       call MPI_BCAST(xc , NumBubble, MPI_REAL8, &
                      root_rank, MPI_Comm_Cart, ierr)
       call MPI_BCAST(yc , NumBubble, MPI_REAL8, &
                      root_rank, MPI_Comm_Cart, ierr)
       call MPI_BCAST(zc , NumBubble, MPI_REAL8, &
                      root_rank, MPI_Comm_Cart, ierr)
    end if ! test_D2P
    
    if (test_randombubs) then
       if ( rank == root_rank ) then
          if (NumBubble <= 2) call pariserror('For random bubble test there has to be more than 2 bubbles')
          call random_bubbles
       endif
       if (filter_random_seeds) then
          call MPI_BCAST(NumBubble, 1, MPI_INTEGER, &
               root_rank, MPI_Comm_Cart, ierr)
       endif
       call MPI_BCAST(rad, NumBubble, MPI_REAL8, &
                      root_rank, MPI_Comm_Cart, ierr)
       call MPI_BCAST(xc , NumBubble, MPI_REAL8, &
                      root_rank, MPI_Comm_Cart, ierr)
       call MPI_BCAST(yc , NumBubble, MPI_REAL8, &
                      root_rank, MPI_Comm_Cart, ierr)
       call MPI_BCAST(zc , NumBubble, MPI_REAL8, &
                      root_rank, MPI_Comm_Cart, ierr)
    endif

    if (test_lattice_bubs) then
       if ( rank == root_rank ) then
          if (nb < 2) &
               call pariserror('For cubic lattice bubble test there has to be at least 2 bubbles per coord direction')
          call cubic_lattice_bubbles
          if(fs_restrictor) call restrict_bubbles("sphere")
       endif
       call MPI_BCAST(NumBubble, 1, MPI_INTEGER, root_rank, MPI_Comm_Cart, ierr)
       call MPI_BCAST(rad, NumBubble, MPI_REAL8, &
            root_rank, MPI_Comm_Cart, ierr)
       call MPI_BCAST(xc , NumBubble, MPI_REAL8, &
            root_rank, MPI_Comm_Cart, ierr)
       call MPI_BCAST(yc , NumBubble, MPI_REAL8, &
            root_rank, MPI_Comm_Cart, ierr)
       call MPI_BCAST(zc , NumBubble, MPI_REAL8, &
            root_rank, MPI_Comm_Cart, ierr)
    endif

    if (test_fcc_lattice) then
       if ( rank == root_rank ) then
          if (n_cell < 2) &
               call pariserror('For cubic lattice bubble test there has to be at least 2 cells per coord dir')
          call fcc_lattice_bubbles
       endif
       call MPI_BCAST(NumBubble, 1, MPI_INTEGER, root_rank, MPI_Comm_Cart, ierr)
       call MPI_BCAST(rad, NumBubble, MPI_REAL8, &
            root_rank, MPI_Comm_Cart, ierr)
       call MPI_BCAST(xc , NumBubble, MPI_REAL8, &
            root_rank, MPI_Comm_Cart, ierr)
       call MPI_BCAST(yc , NumBubble, MPI_REAL8, &
            root_rank, MPI_Comm_Cart, ierr)
       call MPI_BCAST(zc , NumBubble, MPI_REAL8, &
            root_rank, MPI_Comm_Cart, ierr)
    endif

    if(test_heights.or.test_capwave) then 
       if(cylinder_dir==0) then
          write(*,*) "IVOF: Warning: cylinder_dir=0 set to 2"
          cylinder_dir=2
       endif
       ipar=cylinder_dir
       call levelset2vof(wave2ls,ipar)
     else if (test_turb) then
       do i = is,ie; do j=js,je; do k = ks,ke
         if (y(j).LT.0.5) then 
           cvof(i,j,k) = 1.d0 
           vof_flag(i,j,k) = 1
         endif
       enddo; enddo; end do
    else if(NumBubble>0) then
       ! one cylinder in -ipar>0 direction otherwise spheres
       ipar=-cylinder_dir
       if (test_cylinder_advection.OR.test_parasitic_currents)  ipar=-3    ! replaces shapes2lscylinder function argument
       call levelset2vof(shapes2ls,ipar)
    else 
       if(bdry_cond(1) /= 3 .and. rank==0 ) write(*,*) &
            "IVOF: Warning: initializes domain entirely to phase 1 (cvof=0)"
       cvof=0.d0
       vof_flag=0
    endif

    ! hard code for initialized a short jet inside the nozzle
    if (test_jet ) then 
      if ( inject_type == 3 ) then
         NozzleThickness = NozzleThick2Cell*dx(is)
         do i = is,ie; do j=js,je; do k = ks,ke
            if ( x(i) < NozzleLength .and. & 
                 (y(j)-jetcenter_yc) < radius_liq_inject ) then 
               cvof(i,j,k) = 1.d0
               vof_flag(i,j,k) = 1
            end if ! 
         end do; end do; end do
      else if ( inject_type == 4 ) then
         do i = is,ie; do j=js,je; do k = ks,ke
            ryz = sqrt((y(j) - jetcenter_yc)**2.d0 + (z(k) - jetcenter_zc)**2.d0 )
            if ( ryz < radius_liq_inject .and. x(i) < NozzleLength) then 
               cvof(i,j,k) = 1.d0
               vof_flag(i,j,k) = 1
            end if ! 
         end do; end do; end do
      end if ! inject_type
    end if ! test_jet

    if ( test_KHI2D ) then 
      do i = imin,imax; do j=jmin,jmax; do k = kmin,kmax
         sine = yLength*0.5d0 + 0.002*xLength*sin(2.d0*PI*x(i)/xLength) 
         !if ( y(j) > 0.5d0*yLength*(1.d0 + 0.1d0*sin(2.d0*PI*x(i)/xLength)) ) then 
         if ( y(j) - 0.5d0*dy(j) > sine) then 
            cvof(i,j,k) = 1.d0
            !vof_flag(i,j,k) = 1 
         else if ((y(j) + 0.5d0*dy(j) > sine) .and. (y(j) - 0.5d0*dy(j) < sine)) then
            sine=sine*Ny/yLength
            cvof(i,j,k) = 1d0 - (sine-REAL(INT(sine)))
            !write(*,*) 'Intermediate value', cvof(i,j,k)
            !vof_flag(i,j,k) = 1
         endif 
         ! Note: initial pertrubation height = 0.05*yLength 
         ! (when boundary layer thickness delta =0.1*yLength)
      enddo; enddo; enddo
    end if ! test_KHI2D

    if ( test_PhaseInversion ) then 
      do i = is,ie; do j=js,je; do k = ks,ke
         if ( bdry_cond(3) == 1 ) then       ! Quasi-2D, periodic in z direction 
            if (x(i) < 0.5d0*xLength .and. y(j) < 0.5d0*yLength ) then 
               cvof(i,j,k) = 1.d0
               vof_flag(i,j,k) = 1
            end if
         else                                ! 3D 
            if (x(i) < 0.5d0*xLength .and. y(j) < 0.5d0*yLength .and. z(k) < 0.5d0*zLength  ) then 
               cvof(i,j,k) = 1.d0
               vof_flag(i,j,k) = 1
            end if
         end if ! bdry_cond(3)
      end do; end do; end do
    end if ! test_PhaseInversion

    if (test_shear_multiphase) then
      do i = is,ie; do j=js,je; do k = ks,ke
        if ( (y(j)-0.5d0)**2 < 0.01d0 ) then
          cvof(i,j,k) = 1.d0
        endif
      end do; end do; end do
    end if

    if (test_fsdrop) then
       if (.not. FreeSurface) call pariserror('FS has to be switched on for the FS plane test')
       do i=is,ie; do j=js,je; do k=ks,ke
          cvof(i,j,k) = 1.d0-cvof(i,j,k)
       enddo; enddo; enddo 
       call get_flags_and_clip(cvof,vof_flag)
       !call get_vof_phase
    endif

    if (test_plane) then
       call levelset2vof(plane2ls,1)
       call get_flags_and_clip(cvof,vof_flag)
    endif

    call do_all_ghost(cvof)
    call do_all_ighost(vof_flag)
    call setVOFBC(cvof,vof_flag)
    call get_vof_phase(cvof)   
    !call do_all_ighost(vof_phase) ! not needed, and BC for VOF phase are not set. 
    return
  end subroutine initconditions_VOF
  !=================================================================================================
  !   Generate random bubbles 
  !=================================================================================================
  subroutine random_bubbles()
    use module_2phase
#ifdef __INTEL_COMPILER
    use IFPORT
#endif
    implicit none
    integer :: ib, rand_seed, shift, check, compare
    real(8) :: d_min
#ifndef __INTEL_COMPILER
    real :: rand
#endif
    rand_seed = ABS(TIME())
    call srand(rand_seed)
    if(NumBubble>2) then 
      do ib=1,NumBubble
         rad(ib) = r_min + rand()*var_r
         xc(ib)  = (coord_min + rand()*var_coord)
         yc(ib)  = (coord_min + rand()*var_coord)
         zc(ib)  = (coord_min + rand()*var_coord)
         !write(*,'("Bubble ",I4," generated at ",3e14.5," with radius: ",e14.5)')ib,xc(ib),yc(ib),zc(ib),rad(ib)
      end do
    end if ! NumBubble
    if (filter_random_seeds) then
       !calculate min separation distance
       d_min = xLength/(1.5d0*(NumBubble*1.0)**(1.0/3.0))
       write(*,'("Minimum separation distance: ",e14.5)')d_min
       do check=1,NumBubble-1
          do compare=check+1,NumBubble
             if (sqrt((xc(compare)-xc(check))**2.d0+(yc(compare)-yc(check))**2.d0+(zc(compare)-zc(check))**2.d0)<d_min) then
                !delete bub
                do shift=compare,NumBubble-1
                   xc(shift)=xc(shift+1)
                   yc(shift)=yc(shift+1)
                   zc(shift)=zc(shift+1)
                   rad(shift)=rad(shift+1)
                enddo
                NumBubble = NumBubble-1
             endif
          enddo
       enddo
       do compare=1,NumBubble
          write(*,'("Bubble ",I8," generated at ",3e14.5," with radius: ",e14.5)')&
               compare,xc(compare),yc(compare),zc(compare),rad(compare)
       enddo
    endif
  end subroutine random_bubbles 
  !=================================================================================================
  !   Generate bubbles in a primitive cubic lattice 
  !=================================================================================================
  subroutine cubic_lattice_bubbles()
    use module_2phase
#ifdef __INTEL_COMPILER
    use IFPORT
#endif
    implicit none
    integer :: ib, jb, kb, rand_seed, index
#ifndef __INTEL_COMPILER
    real :: rand
#endif
    rand_seed = ABS(TIME())
    call srand(rand_seed)
    NumBubble = nb*nb*nb
    write(*,'("Total number of bubbles: ",I10)')NumBubble
    do ib=0,nb-1; do jb=0,nb-1; do kb=0,nb-1
       index = 1 + kb + jb*nb + ib*nb*nb 
       rad(index) = r_min + rand()*var_r
       xc(index)  = coord_min + var_coord/(nb-1)*(ib+var_pos/(3.0**(0.5d0))*(2.0d0*rand()-1.0d0))
       yc(index)  = coord_min + var_coord/(nb-1)*(jb+var_pos/(3.0**(0.5d0))*(2.0d0*rand()-1.0d0))
       zc(index)  = coord_min + var_coord/(nb-1)*(kb+var_pos/(3.0**(0.5d0))*(2.0d0*rand()-1.0d0))
!!$       write(*,'("Bubble ",I10," generated at ",3e14.5," with radius: ",e14.5)')&
!!$            index,xc(index),yc(index),zc(index),rad(index)
    enddo; enddo; enddo
  end subroutine cubic_lattice_bubbles

    !===================================================================================[primitive restrictor WA]
    subroutine restrict_bubbles(shape)
    use module_2phase
#ifdef __INTEL_COMPILER
    use IFPORT
#endif
    implicit none
    character(len=6), intent(in) :: shape
    integer :: ib,cib,tokill,numbub2
    real(8) :: dist,rad_max,scx,scy,scz
    real(8), dimension(:,:), allocatable :: listcopy
    !call deb_inf("|restrict bubbles ")

    tokill=0
    
    !select desired shape, then start deletion    
    select case(shape) !for now, centered with fixed radius
    case("sphere")
       !sphere radius
       rad_max=0.46*xlength
       !sphere center
       scx=xlength/2.d0;scy=ylength/2.d0;scz=zlength/2.d0
       do ib=1,numbubble;
          dist=sqrt((xc(ib)-scx)**2+(yc(ib)-scy)**2+(zc(ib)-scz)**2)
          if(dist>rad_max) then
             print *,'killing',ib,'with distance',dist,'radmax is',rad_max
             tokill=tokill+1
             rad(ib)=0.d0
          endif          
       enddo
       !decrease number of bubbles:
       numbub2=numbubble-tokill
       write(*,*) "restrictor: will kill",tokill,"of",numbubble,"bubbles leaving",numbub2

       !The block below is a primitive attempt to (1) copy the nonzero bubbles, (2) zero the bubble
       !vectors xc,yc,zc and rad and (3) copy the stuff back. This way, all bubbles are at the beginning
       !of the list.
       allocate(listcopy(1:numbub2,1:4))
       !make another pass with rewrite
       cib=1
       do ib=1,numbubble
          if(rad(ib)>0) then
             listcopy(cib,1)=xc(ib);listcopy(cib,2)=yc(ib);listcopy(cib,3)=zc(ib)
             listcopy(cib,4)=rad(ib)             
             cib=cib+1
          endif
       enddo
       numbubble=numbub2
       xc=0;yc=0;zc=0;rad=0
       do ib=1,numbubble
          xc(ib)=listcopy(ib,1)
          yc(ib)=listcopy(ib,2)
          zc(ib)=listcopy(ib,3)
          rad(ib)=listcopy(ib,4)
       enddo


    case default
       if(rank==0 ) call pariserror("No such shape in bubble restrictor")
       stop
    end select
    deallocate(listcopy)
    return
  end subroutine restrict_bubbles
  
  !=================================================================================================
  !   Generate bubbles in a face centered cubic lattice 
  !   Total bubbles should equal:
  !   NumBubble = 4n^3 + 6n^2 + 3n + 1 with n=n_cell 
  !=================================================================================================
  subroutine fcc_lattice_bubbles()
    use module_2phase
#ifdef __INTEL_COMPILER
    use IFPORT
#endif
    implicit none
    integer :: ib, jb, kb, rand_seed, index
    real(8) :: L_d
#ifndef __INTEL_COMPILER
    real :: rand
#endif
    rand_seed = ABS(TIME())
    call srand(rand_seed)
    index=2
    NumBubble=0
    !Calculate mean separation distance in bulk
    L_d=var_coord/(real(n_cell)*sqrt(2.d0))
    do ib=0,2*n_cell; do jb=0,2*n_cell; do kb=0,2*n_cell
       if (mod(index,2)==0) then
          NumBubble=NumBubble+1
          rad(NumBubble) = r_min + rand()*var_r
          xc(NumBubble)= coord_min + ib*var_coord/(n_cell*2.0d0)+var_pos*L_d/(sqrt(3.d0))*(2.d0*rand()-1.d0)
          yc(NumBubble)= coord_min + jb*var_coord/(n_cell*2.0d0)+var_pos*L_d/(sqrt(3.d0))*(2.d0*rand()-1.d0)
          zc(NumBubble)= coord_min + kb*var_coord/(n_cell*2.0d0)+var_pos*L_d/(sqrt(3.d0))*(2.d0*rand()-1.d0)  
       endif
       index=index+1
    enddo; enddo; enddo
    write(*,'("Total number of bubbles: ",I10)')NumBubble
  end subroutine fcc_lattice_bubbles
  !=================================================================================================
  !   Spheres and cylinders
  !=================================================================================================
  function shapes2ls(xx,yy,zz,ipar)
    use module_2phase
    implicit none
    real(8), intent(in) :: xx,zz,yy
    integer, intent(in) :: ipar
    real(8) :: a, cdir(0:3), shapes2ls
    integer ib

    if(.not.(-3<=ipar.and.ipar<=0)) call pariserror("S: invalid ipar")
    cdir(1:3) = 1.d0/(1.d0 + excentricity(1:3))
    cdir(-ipar) = 0.d0
    shapes2ls = -2.d6
    ! ipar < 0 cylinder case
    ! ipar = 0 spheres
    if(ipar < 0.and.NumBubble/=1) call pariserror("S: invalid NumBubbles")
    do ib=1,NumBubble
! x-axis span: 2*rad/(sqrt cdir) = 2*rad*sqrt(1 + e(1)) etc... 
       a = rad(ib)**2 - (cdir(1)*(xx-xc(ib))**2+cdir(2)*(yy-yc(ib))**2+cdir(3)*(zz-zc(ib))**2)
       shapes2ls = MAX(shapes2ls,a)
    end do
  end function shapes2ls
  !=================================================================================================
  !  sine-wave interface
  !=================================================================================================
  function wave2ls(xx,yy,zz,ipar)
    use module_2phase
    implicit none
    real(8) wave2ls
    real(8), intent(in) :: xx,zz,yy
    integer, intent(in) :: ipar
    integer :: hpar,vpar,normalsign
    real(8) ::  vdir(13),hdir(3)

    if(normal_up) then
       normalsign=1
    else
       normalsign=-1
    endif

    if(.not.(1<=ipar.and.ipar<=3)) call pariserror("invalid ipar")
    if (ipar == 1) then
       hpar = 2
       vpar = 3
    else if (ipar == 2) then
       hpar = 1
       vpar = 3
    else
       hpar = 1
       vpar = 2
    endif
       
    vdir = 0.d0
    vdir(vpar) = 1.d0
    hdir = 0.d0
    hdir(hpar) = 1.d0

    wave2ls = vdir(3)*(- zz + zLength/2.d0) + vdir(2)*(-yy + yLength/2.d0) &
         + vdir(1)*(-xx + xLength/2.d0) &
         + 0.5*A_h*dx(nx/2+2)*cos(2.*PI*(hdir(1)*xx/xLength &
         + hdir(2)*yy/yLength + hdir(3)*zz/zLength))
    wave2ls = wave2ls*dble(normalsign) + hshift*dx(nx/2+2)
  end function wave2ls
  !=================================================================================================
  !  Planar interface
  !=================================================================================================
  function plane2ls(xx,yy,zz)
    use module_2phase
    implicit none
    real(8) plane2ls
    real(8), intent(in) :: xx,zz,yy
    plane2ls = plane -xx*n_p(1) -yy*n_p(2) -zz*n_p(3)
  end function plane2ls
  !=================================================================================================
  !   Converts a level-set field into a VOF field
  !=================================================================================================
  subroutine levelset2vof(lsfunction,ipar)
    use module_BC
    implicit none
    real(8), external :: lsfunction
    integer, intent(in) :: ipar
    include 'mpif.h'
    integer :: ierr, req(12),sta(MPI_STATUS_SIZE,12)
    integer , parameter :: ngh=2
 
#ifdef HAVE_VOFI
    integer :: i,j,k,itrue
    integer, parameter :: ndim0=3
    real(8), dimension(3) :: xv,xloc
    real(8) :: h0,fh,cc
    real(8), external :: get_cc,get_fh
#else
    if(use_vofi) then
       if(rank==0) print*, "vofi not available, falling back on refined levelset2vof"
       use_vofi=.false.
    endif
#endif

   

    if(use_vofi) then
#ifdef HAVE_VOFI
       itrue = 1
       ! starting point to get fh
       h0 = dx(3)
       xv(1) = x(1); xv(2) = y(1); xv(3) = z(1);
       fh = get_fh(vofi_lsfunction,xv,h0,ndim0,itrue)
       
       ! xloc: minor vertex of each cell of the grid 
       do k=kmin,kmax; do j=jmin,jmax; do i=imin,imax
          xloc(1) = x(i) - 0.5d0*dx(i); xloc(2) = y(j) - 0.5*dy(j); xloc(3) = z(k) - 0.5d0*dz(k)
          cc = get_cc(vofi_lsfunction,xloc,h0,fh,ndim0)
          cvof(i,j,k) = cc
          if(cc.gt.EPSC.and.cc.lt.1d0-EPSC) then
             vof_flag(i,j,k) = 2
          else if(cc.gt.-EPSC.and.cc.le.EPSC) then
             vof_flag(i,j,k) = 0
          else if(cc.ge.1d0-EPSC.and.cc.lt.1d0+EPSC) then
             vof_flag(i,j,k) = 1
          else
             print *, i,j,k,cc
             call pariserror("Vofi failed to initialize")
          endif
       enddo; enddo; enddo
#endif
    else
       call ls2vof_refined(lsfunction,ipar,1)
       call ighost_x(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ighost_y(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ighost_z(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       
       call ghost_x(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ghost_y(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ghost_z(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call setVOFBC(cvof,vof_flag)
       call ls2vof_refined(lsfunction,ipar,refinement)
    endif
  end subroutine levelset2vof
  
  function vofi_lsfunction(xyz)
    use module_2phase
    real(8) vofi_lsfunction
    real(8), intent(in) :: xyz(3)
    real(8) :: x,y,z
    integer :: ipar

    x = xyz(1); y = xyz(2); z = xyz(3)
    ! repeats the code in init cond
    if(test_heights) then 
       if(cylinder_dir==0) then
          if(rank==0) print *, "IVOF: Warning: cylinder_dir=0 impossible for test_height, set to 2"
          cylinder_dir=2
       endif
       ipar=cylinder_dir
       vofi_lsfunction = wave2ls(x,y,z,ipar)
    else if(NumBubble>0) then
       ! A cylinder with axis in direction -ipar>0 otherwise a set of spheres
       ipar=-cylinder_dir
       if (test_cylinder_advection.OR.test_parasitic_currents) ipar=-3
       vofi_lsfunction = shapes2ls(x,y,z,ipar)
    else
       call pariserror("unexpected case in vofi_ls")
       vofi_lsfunction=0d0 ! to avoid warnings
    endif
    vofi_lsfunction = - vofi_lsfunction
  end function vofi_lsfunction


  !=================================================================================================
  !   Finds cells surrounded by all empty or all full cells
  !=================================================================================================
  function notisolated(istencil3x3,is2D)
    implicit none
    logical :: notisolated
    integer, intent(in) :: is2D
    integer :: istencil3x3(-1:1,-1:1,-1:1)
    integer :: nfrac
    integer :: nfull
    integer :: i0,j0,k0
    integer :: isnot2D, ncells

    nfrac=0; nfull=0; ncells=0
    isnot2D = 1 - is2D
!     write(*,*) "isnot2D = ", isnot2D
!     if(is2D==1) write(*,'(9I1)',advance='no') istencil3x3(:,:,0)
    do i0=-1,1; 
       do j0=-1,1; 
          do k0=-isnot2D,isnot2D
             nfrac = nfrac + istencil3x3(i0,j0,k0)/2
             nfull = nfull + mod(istencil3x3(i0,j0,k0),2)
             ncells = ncells + 1
          enddo; enddo; enddo
    notisolated = .not.(nfrac==0.and.(nfull==0.or.nfull==ncells))
    ! write (*,*) " nfrac,nfull,ncells ",  nfrac,nfull,ncells
  end function notisolated
!=================================================================================================
!
! mapping
!
!=================================================================================================
  subroutine map3x3in2x2(i1,j1,k1,i0,j0,k0)
    implicit none
    integer, intent(in) :: i0,j0,k0
    integer, intent(out) :: i1(-1:1,-1:1,3), j1(-1:1,-1:1,3), k1(-1:1,-1:1,3)
    integer m,n
    do m=-1,1
       do n=-1,1
          !  d=1
          i1(m,n,1) = i0
          j1(m,n,1) = m + j0
          k1(m,n,1) = n + k0
          ! d=2
          i1(m,n,2) = m + i0
          j1(m,n,2) = j0
          k1(m,n,2) = n + k0
          ! d=3
          i1(m,n,3) = m + i0
          j1(m,n,3) = n + j0
          k1(m,n,3) = k0 
       enddo
    enddo
  end subroutine map3x3in2x2

  subroutine ls2vof_refined(lsfunction,ipar,n1)
    implicit none
    real(8), external :: lsfunction
    integer, intent(in) :: ipar,n1
    real(8) :: stencil3x3(-1:1,-1:1,-1:1),dx1,dy1,dz1,x0,y0,z0,x1,y1,z1,a,b
    integer :: i,j,k,i0,j0,k0,l,m,n,s
    integer :: nfrac,nflag,nfull
    integer :: istencil3x3(-1:1,-1:1,-1:1)
    integer :: i1(-1:1,-1:1,3), j1(-1:1,-1:1,3), k1(-1:1,-1:1,3)
    logical :: refinethis 
    real(8) :: count
    integer :: calc_imax
    integer :: dirselect(0:3), d, is2D,max_flag

! initialization
    count=0.d0
    refinethis = .false.

! Initialize 2D/3D switch
    dirselect = 1  ! spheres: all directions selected: default. 
    d = max(ipar,-ipar)
    dirselect(d)=0

!    Some error checking
    if(d>3) call pariserror("wrong ipar")
    if(min(min(nx,ny),nz)<2) call pariserror("minimum dimension nx ny nz too small")
    max_flag=calc_imax(vof_flag)

    if(n1>1.and.max_flag/=2) then
      if(rank==0) then
          if(max_flag==0) then
             write(*,*) "ls2vof_refined: error: single phase flow ? Nothing to initialize !?"
          else
            write(*,*) "ls2vof_refined: maximum vof_flag = ",max_flag,"but expecting maximum flag = 2"
          endif
       endif
       call pariserror("bad flag")
    endif
 
! main loop
    do k=ks,ke; do j=js,je; do i=is,ie
       is2D=0  ! Default
       if(n1>1) then  ! refinement on second pass only
          if(d==0) then ! check for isolated cells
             do i0=-1,1; do j0=-1,1; do k0=-1,1
                istencil3x3(i0,j0,k0) = vof_flag(i+i0,j+j0,k+k0)
             enddo; enddo; enddo
          else if(d>0) then 
             call map3x3in2x2(i1,j1,k1,i,j,k)
             do m=-1,1; do n=-1,1
                istencil3x3(m,n,0) = vof_flag(i1(m,n,d),j1(m,n,d),k1(m,n,d))
             enddo; enddo
             is2D=1
          else
             call pariserror("bad d")
          endif
          refinethis = notisolated(istencil3x3,is2D)
          if(refinethis) count = count + 1.
       endif
! refine and initialize subcells
       if(n1==1.or.refinethis) then  ! if n1>1 and no refine leave cell as is. 
          dx1 = dx(i)/n1; dy1 = dy(j)/n1; dz1 = dz(k)/n1
          nfrac=0; nfull=0
          b=0.d0
          s=n1/2
          do l=0,n1-1; do m=0,n1-1; do n=0,n1-1
             x0 = x(i) - 0.5d0*dx(i) + 0.5d0*dx1 + dx1*l
             y0 = y(j) - 0.5d0*dy(j) + 0.5d0*dy1 + dy1*m
             z0 = z(k) - 0.5d0*dz(k) + 0.5d0*dz1 + dz1*n
             do i0=-1,1; do j0=-1,1; do k0=-1,1
                x1 = x0 + i0*dx1; y1 = y0 + j0*dy1; z1 = z0 + k0*dz1 
                stencil3x3(i0,j0,k0) = lsfunction(x1,y1,z1,ipar)
             enddo; enddo; enddo
             call ls2vof_in_cell(stencil3x3,a,nflag)
             if(nflag==2) then 
                nfrac = nfrac + 1
             else if(nflag==1) then
                nfull = nfull + 1
             endif
             b=b+a   ! *(1)*
          enddo; enddo; enddo
          cvof(i,j,k) = b/(n1**3)
          if(nfrac > 0) then
             vof_flag(i,j,k) = 2
             ! now either all full, all empty, or mix full/empty : 
          else if(nfull==n1**3) then ! all full
             vof_flag(i,j,k) = 1
             cvof(i,j,k) = 1.d0  ! because arithmetic at (1) may cause round off errors
          else if(nfull==0) then ! all empty
             vof_flag(i,j,k) = 0
             cvof(i,j,k) = 0.d0  ! paranoid programming.
          else ! mix of full and empty
             vof_flag(i,j,k) = 2
          end if
       endif
    enddo; enddo; enddo
!    IF(N1>1) write(*,*) "proportion refined ",100.*count/(nx*ny*nz),"%"
    return
  end subroutine ls2vof_refined
  !=================================================================================================
  subroutine c_mask(cvof, cbinary)
    implicit none
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(out) :: cbinary
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: cvof

    where (cvof > 0.5d0)
       cbinary = 1.d0
    elsewhere
       cbinary = 0.d0
    end where
  end subroutine c_mask
  !=================================================================================================
  subroutine get_half_fractions(c,d,i,j,k,vofh1,vofh2)

    implicit none
    logical error
    integer :: i,j,k,d
    integer :: i1,j1,k1
    real(8)  , dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: c
    real(8) vofh1, vofh2
    real(8) alpha,fl3d,stencil3x3(-1:1,-1:1,-1:1)
    real(8) dm(3),x0(3),deltax(3)

    vofh1 = c(i,j,k)
    vofh2 = c(i,j,k)
    if ((c(i,j,k).gt.0.d0).and.(c(i,j,k).lt.1.d0)) then
       do i1=-1,1; do j1=-1,1; do k1=-1,1
          stencil3x3(i1,j1,k1) = c(i+i1,j+j1,k+k1)
       enddo;enddo;enddo
       call fit_plane_new(c(i,j,k),d,0.d0,0.d0,stencil3x3,dm,alpha,error)
       if (error) then
          vofh1 = c(i,j,k)
          vofh2 = c(i,j,k)
       else
          x0=0d0
          deltax=1d0
          deltax(d)=0.5d0
          vofh1 = 2.0d0*fl3d(dm,alpha,x0,deltax)
          x0(d)=0.5d0
          vofh2 = 2.0d0*fl3d(dm,alpha,x0,deltax)
       endif
    endif

  end subroutine get_half_fractions

  subroutine vofsweeps(tswap)
    use module_BC
    use module_flow
    use module_tmpvar
    implicit none
    integer, intent(in) :: tswap

    if (VOF_advect=='Dick_Yue') call c_mask(cvof, work(:,:,:,2))
    if (MOD(tswap,3).eq.0) then  ! do z then x then y 
       call swp(w,cvof,vof_flag,3,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))
       call swp(u,cvof,vof_flag,1,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))
       call swp(v,cvof,vof_flag,2,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))
    elseif (MOD(tswap,3).eq.1) then ! do y z x
       call swp(v,cvof,vof_flag,2,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))
       call swp(w,cvof,vof_flag,3,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))
       call swp(u,cvof,vof_flag,1,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))
    else ! do x y z
       call swp(u,cvof,vof_flag,1,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))
       call swp(v,cvof,vof_flag,2,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))
       call swp(w,cvof,vof_flag,3,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))
   endif
  end subroutine vofsweeps
!=================================================================================================
  subroutine get_momentum_centered(us,d,mom)

  use module_grid
  use module_flow
  use module_tmpvar

  integer i,j,k,d
  integer i0,j0,k0
  real(8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), intent(in)  :: us
  real(8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), intent(out) :: mom
  real(8) :: rhoavg

  tmp = cvof
  work(:,:,:,1) = 0.d0
  work(:,:,:,2) = 0.d0
  call init_i0j0k0 (d,i0,j0,k0)
  
  call do_all_ghost(us)

  do k=ks-1,ke+1
     do j=js-1,je+1
        do i=is-1,ie+1
           rhoavg = tmp(i,j,k)*rho2 + (1.d0-tmp(i,j,k))*rho1
           mom(i,j,k) = 0.5d0*(us(i,j,k)+us(i-i0,j-j0,k-k0))*rhoavg
        enddo
     enddo
  enddo

  call do_all_ghost(tmp)
  call do_all_ghost(mom)

  end subroutine get_momentum_centered

  subroutine get_momentum_staggered(us,d,mom)

  use module_grid
  use module_flow
  use module_tmpvar

  integer i,j,k,d
  integer i0,j0,k0
  real(8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), intent(inout)  :: us
  real(8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), intent(out) :: mom
  real(8) :: rhoavg

  tmp = cvof
  work(:,:,:,1) = 0.d0
  work(:,:,:,2) = 0.d0
  call init_i0j0k0 (d,i0,j0,k0)
  
  do k=ks-1,ke+1
     do j=js-1,je+1
        do i=is-1,ie+1
            call get_half_fractions(cvof,d,i,j,k,work(i,j,k,1),work(i,j,k,2))
        enddo
     enddo
  enddo

  call do_all_ghost(work(:,:,:,1))
  call do_all_ghost(work(:,:,:,2))
  call do_all_ghost(us)

  do k=ks-1,ke+1
     do j=js-1,je+1
        do i=is-1,ie+1
           tmp(i,j,k) = 0.5d0*(work(i,j,k,2)+work(i+i0,j+j0,k+k0,1))
           rhoavg = tmp(i,j,k)*rho2 + (1.d0-tmp(i,j,k))*rho1
           mom(i,j,k) = us(i,j,k)*rhoavg
        enddo
     enddo
  enddo

!  call SetMomentumBC(us,cvof,mom,d,umask,rho1,rho2,0.d0) 
! TEMPORARY 
   if ( d == 1 ) then 
      call SetMomentumBC(us,cvof,mom,d,umask,rho1,rho2,0.d0)
   else if ( d == 2 ) then 
      call SetMomentumBC(us,cvof,mom,d,vmask,rho1,rho2,0.d0)
   else if ( d == 3 ) then 
      call SetMomentumBC(us,cvof,mom,d,wmask,rho1,rho2,0.d0)
   end if !d
! END TEMPORARY 

  call do_all_ghost(tmp)
  call do_all_ghost(mom)

  end subroutine get_momentum_staggered

  subroutine get_velocity_from_momentum_centered (mom,us,der,d)

  use module_grid
  use module_flow
  use module_BC
  use module_tmpvar

  implicit none
  integer :: i,j,k,d
  integer :: i0,j0,k0
  real(8)  , dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: mom
  real(8)  , dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: us,der
  real(8) :: cflag0,rhoavg0,cflag1, rhoavg1

  call init_i0j0k0 (d,i0,j0,k0)

  do k=ks-1,ke+1
    do j=js-1,je+1
      do i=is-1,ie+1
          cflag0 = tmp(i,j,k)
          rhoavg0 = rho2*cflag0 + (1.d0-cflag0)*rho1
          cflag1 = tmp(i+i0,j+j0,k+k0)
          rhoavg1 = rho2*cflag1 + (1.d0-cflag1)*rho1
          if (((cflag0.gt.0.d0).and.(cflag0.lt.1.d0)).OR. &
              ((cflag1.gt.0.d0).and.(cflag1.lt.1.d0))) then
            der(i,j,k) = ((mom(i,j,k)/rhoavg0+mom(i+i0,j+j0,k+k0)/rhoavg1)*0.5d0 - us(i,j,k))/dt
          endif
      enddo
    enddo
  enddo

  call do_all_ghost(der)

end subroutine get_velocity_from_momentum_centered

subroutine get_velocity_from_momentum_staggered (mom,us,der,d)
  use module_grid
  use module_flow
  use module_BC
  use module_tmpvar
  implicit none
  integer :: i,j,k
  integer :: i0,j0,k0,d
  integer :: iend,jend,kend
  integer :: iini,jini,kini
  real(8)  , dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: mom
  real(8)  , dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: us,der
  real(8) :: cflag,rhoavg

  call init_i0j0k0 (d,i0,j0,k0)

  kini=ks
  kend=ke
  jini=js
  jend=je
  iini=is
  iend=ie 

  if (d.eq.1) then
     iend=ieu 
  elseif (d.eq.2) then
     jend=jev 
  elseif (d.eq.3) then
     kend=kew 
  endif

  do k=kini,kend
     do j=jini,jend
        do i=iini,iend
           ! if interface rewrite interface velocity
           !      cflag = (tmp(i,j,k) + cvof(i,j,k) + cvof(i+i0,j+j0,k+k0))/3.d0
           cflag = tmp(i,j,k)

           if (((cflag.gt.0.d0).and.(cflag.lt.1.d0))) then
              rhoavg = rho2*cflag + (1.d0-cflag)*rho1
              der(i,j,k) = (mom(i,j,k)/rhoavg - us(i,j,k))/dt
           endif
        enddo
     enddo
  enddo

  call do_all_ghost(der)

end subroutine get_velocity_from_momentum_staggered

subroutine vofandmomsweepsstaggered(tswap,t)
  use module_BC
  use module_flow
  use module_tmpvar
  implicit none
  integer dir
  integer, intent(in) :: tswap
  real(8) :: t

  if (staggeredmom) then
     do dir=1,3
        if (dir.eq.1) then
           call get_momentum_staggered(u,1,momentum)
        elseif (dir.eq.2) then
           call get_momentum_staggered(v,2,momentum)
        elseif (dir.eq.3) then
           call get_momentum_staggered(w,3,momentum)
        endif

        tmp_flag = vof_flag

        if (VOF_advect=='Dick_Yue') call c_mask(tmp, work(:,:,:,2))
        if (MOD(tswap,3).eq.0) then  ! do z then x then y 
           call swpmom_stg(w,tmp,3,work(:,:,:,1),work(:,:,:,2), &
                work(:,:,:,3),momentum,dir,t)
           call swp_stg(w,tmp,tmp_flag,3,work(:,:,:,1),work(:,:,:,2),&
                work(:,:,:,3),dir)

           call swpmom_stg(u,tmp,1,work(:,:,:,1),work(:,:,:,2), &
                work(:,:,:,3),momentum,dir,t)
           call swp_stg(u,tmp,tmp_flag,1,work(:,:,:,1),work(:,:,:,2),&
                work(:,:,:,3),dir)

           call swpmom_stg(v,tmp,2,work(:,:,:,1),work(:,:,:,2), &
                work(:,:,:,3),momentum,dir,t)
           call swp_stg(v,tmp,tmp_flag,2,work(:,:,:,1),work(:,:,:,2),&
                work(:,:,:,3),dir)

        elseif (MOD(tswap,3).eq.1) then ! do y z x

           call swpmom_stg(v,tmp,2,work(:,:,:,1),work(:,:,:,2), &
                           work(:,:,:,3),momentum,dir,t)
           call swp_stg(v,tmp,tmp_flag,2,work(:,:,:,1),work(:,:,:,2),&
                        work(:,:,:,3),dir)

           call swpmom_stg(w,tmp,3,work(:,:,:,1),work(:,:,:,2), &
                work(:,:,:,3),momentum,dir,t)
           call swp_stg(w,tmp,tmp_flag,3,work(:,:,:,1),work(:,:,:,2),&
                work(:,:,:,3),dir)

           call swpmom_stg(u,tmp,1,work(:,:,:,1),work(:,:,:,2), &
                work(:,:,:,3),momentum,dir,t)
           call swp_stg(u,tmp,tmp_flag,1,work(:,:,:,1),work(:,:,:,2),&
                work(:,:,:,3),dir)

        else ! do x y z

           call swpmom_stg(u,tmp,1,work(:,:,:,1),work(:,:,:,2), &
                work(:,:,:,3),momentum,dir,t)
           call swp_stg(u,tmp,tmp_flag,1,work(:,:,:,1),work(:,:,:,2),&
                work(:,:,:,3),dir)

           call swpmom_stg(v,tmp,2,work(:,:,:,1),work(:,:,:,2), &
                work(:,:,:,3),momentum,dir,t)
           call swp_stg(v,tmp,tmp_flag,2,work(:,:,:,1),work(:,:,:,2),&
                work(:,:,:,3),dir)

           call swpmom_stg(w,tmp,3,work(:,:,:,1),work(:,:,:,2), &
                work(:,:,:,3),momentum,dir,t)
           call swp_stg(w,tmp,tmp_flag,3,work(:,:,:,1),work(:,:,:,2),&
                work(:,:,:,3),dir)

        endif

        if (dir.eq.1) then
           call get_velocity_from_momentum_staggered (momentum,u,du,1)
        elseif (dir.eq.2) then
           call get_velocity_from_momentum_staggered (momentum,v,dv,2)
        elseif (dir.eq.3) then
           call get_velocity_from_momentum_staggered (momentum,w,dw,3)
        endif

     enddo

  else
     do dir=1,3
        if (dir.eq.1) then
           call get_momentum_centered(u,1,momentum)
        elseif (dir.eq.2) then
           call get_momentum_centered(v,2,momentum)
        elseif (dir.eq.3) then
           call get_momentum_centered(w,3,momentum)
        endif

        tmp_flag = vof_flag

        if (VOF_advect=='Dick_Yue') call c_mask(tmp, work(:,:,:,2))
        if (MOD(tswap,3).eq.0) then  ! do z then x then y 
           call swpmom(w,tmp,3,work(:,:,:,1),work(:,:,:,2), &
                work(:,:,:,3),momentum,t)
           call swp(w,tmp,tmp_flag,3,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))

           call swpmom(u,tmp,1,work(:,:,:,1),work(:,:,:,2), &
                work(:,:,:,3),momentum,t)
           call swp(u,tmp,tmp_flag,1,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))

           call swpmom(v,tmp,2,work(:,:,:,1),work(:,:,:,2), &
                work(:,:,:,3),momentum,t)
           call swp(v,tmp,tmp_flag,2,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))

        elseif (MOD(tswap,3).eq.1) then ! do y z x

           call swpmom(v,tmp,2,work(:,:,:,1),work(:,:,:,2), &
                work(:,:,:,3),momentum,t)
           call swp(v,tmp,tmp_flag,2,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))

           call swpmom(w,tmp,3,work(:,:,:,1),work(:,:,:,2), &
                work(:,:,:,3),momentum,t)
           call swp(w,tmp,tmp_flag,3,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))

           call swpmom(u,tmp,1,work(:,:,:,1),work(:,:,:,2), &
                work(:,:,:,3),momentum,t)
           call swp(u,tmp,tmp_flag,1,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))

        else ! do x y z

           call swpmom(u,tmp,1,work(:,:,:,1),work(:,:,:,2),&
                work(:,:,:,3),momentum,t)
           call swp(u,tmp,tmp_flag,1,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))

           call swpmom(v,tmp,2,work(:,:,:,1),work(:,:,:,2), &
                work(:,:,:,3),momentum,t)
           call swp(v,tmp,tmp_flag,2,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))

           call swpmom(w,tmp,3,work(:,:,:,1),work(:,:,:,2), &
                work(:,:,:,3),momentum,t)
           call swp(w,tmp,tmp_flag,3,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))

        endif

        if (dir.eq.1) then
           call get_velocity_from_momentum_centered (momentum,u,du,1)
        elseif (dir.eq.2) then
           call get_velocity_from_momentum_centered (momentum,v,dv,2)
        elseif (dir.eq.3) then
           call get_velocity_from_momentum_centered (momentum,w,dw,3)
        endif


     enddo

  endif

end subroutine vofandmomsweepsstaggered
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
! subroutine SetVOFBC: Sets the VOF fraction boundary condition
!-------------------------------------------------------------------------------------------------
  subroutine SetVOFBC(cv,fl)
    use module_grid
    use module_BC
    use module_2phase
    use module_solid
    implicit none
    include 'mpif.h'
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: cv  ! cvof
    integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: fl  ! vof_flag
    integer :: fb(6),d,l,m,n,c(3),try(2:3),sign,orientation,flag,flhere
    real(8) :: cvhere

    do orientation=1,6
       cond = vofbdry_cond(orientation)
       if(cond=='wet') then
          fb(orientation)=1
       else if(cond=='dry') then
          fb(orientation)=0
       else if(cond=='periodic') then
          fb(orientation) = 3  !  will do nothing
       else if(cond=='outflow') then
          fb(orientation) = 4  ! will copy inflow
       else if(cond=='jet') then
          if(orientation /= 1) call pariserror("jet only at x-")
          fb(orientation) = 2
       else if(cond=='90deg') then 
          fb(orientation) = 5
       else
          call pariserror("this vofbc not implemented")
       endif
    enddo
! 
    do orientation=1,6
       ! orientation order:    
       ! x- y- z- x+ y+ z+
       d = orientation
       sign=-1
       if(orientation>3) then
          d = orientation-3
          sign=1
       endif
       flag=fb(orientation)
! sort directions so that try(1) = d and try(2),try(3) are any other two directions. 
       m=1
       n=2  
       do while (m.le.3)
          if(m.ne.d) then
             try(n) = m
             n=n+1
          endif
          m=m+1
       enddo
! end sort
       if(coords(d)==proclimit(d,sign)) then
          do l=coordstart(try(2))-Ng,coordend(try(2))+Ng
             do m=coordstart(try(3))-Ng,coordend(try(3))+Ng
                c(try(2)) = l; c(try(3)) = m
                c(d)=coordlimit(d,sign) + sign
                if(flag<2) then
                   cv(c(1),c(2),c(3))=dble(flag)
                   fl(c(1),c(2),c(3))=flag
                   c(d) = c(d) + sign
                   cv(c(1),c(2),c(3))=dble(flag)
                   fl(c(1),c(2),c(3))=flag
                elseif(flag==2) then  ! jet
                   cv(c(1),c(2),c(3))=dble(inject(c(2),c(3)))
                   fl(c(1),c(2),c(3))=inject(c(2),c(3))
                elseif(flag==4) then
                   c(d)=coordlimit(d,sign)
                   cvhere=cv(c(1),c(2),c(3))
                   flhere=fl(c(1),c(2),c(3))
                   c(d) = c(d) + sign
                   cv(c(1),c(2),c(3))=cvhere
                   fl(c(1),c(2),c(3))=flhere
                   c(d) = c(d) + sign
                   cv(c(1),c(2),c(3))=cvhere
                   fl(c(1),c(2),c(3))=flhere
                elseif(flag==5) then !90deg 
                   if (d==1 .and. sign==-1) then 
                      cv(c(1),c(2),c(3))=cv(c(1)+1,c(2),c(3))
                      fl(c(1),c(2),c(3))=fl(c(1)+1,c(2),c(3))
                   else if (d==1 .and. sign==1) then
                      cv(c(1),c(2),c(3))=cv(c(1)-1,c(2),c(3))
                      fl(c(1),c(2),c(3))=fl(c(1)-1,c(2),c(3))
                   else if (d==2 .and. sign==-1) then 
                      cv(c(1),c(2),c(3))=cv(c(1),c(2)+1,c(3))
                      fl(c(1),c(2),c(3))=fl(c(1),c(2)+1,c(3))
                   else if (d==2 .and. sign==1) then
                      cv(c(1),c(2),c(3))=cv(c(1),c(2)-1,c(3))
                      fl(c(1),c(2),c(3))=fl(c(1),c(2)-1,c(3))
                   else if (d==3 .and. sign==-1) then 
                      cv(c(1),c(2),c(3))=cv(c(1),c(2),c(3)+1)
                      fl(c(1),c(2),c(3))=fl(c(1),c(2),c(3)+1)
                   else if (d==3 .and. sign==1) then
                      cv(c(1),c(2),c(3))=cv(c(1),c(2),c(3)-1)
                      fl(c(1),c(2),c(3))=fl(c(1),c(2),c(3)-1)
                   end if !d
                endif
             enddo
          enddo
       endif
    enddo
  contains
    function inject(j,k)
      use module_grid
      implicit none
      integer :: j,k
      integer :: inject
      inject=0
      if ( inject_type == 2 .or. inject_type == 4) then 
         if ((y(j) - jetcenter_yc)**2 + (z(k) - jetcenter_zc)**2.lt.radius_liq_inject**2) inject=1
      else if ( inject_type == 3 ) then
         if ((y(j) - jetcenter_yc) <= radius_liq_inject ) inject = 1 
      else if ( inject_type == 5 ) then
         if ((y(j) - jetcenter_yc) > radius_gas_inject .and. &  
             (y(j) - jetcenter_yc) < radius_liq_inject ) inject = 1 
      end if ! inject_type
    end function inject
    !
    function xcoord(d,i)
      implicit none
      real(8) :: xcoord
      integer, intent(in) :: d,i
      if(d==1) then
         xcoord = x(i)
      else if(d==2) then
         xcoord = y(i)
      else
         xcoord = z(i)
      endif
    end function xcoord
  end subroutine SetVOFBC
!=================================================================================================
  subroutine test_cell_size()
    implicit none
    if(dabs(dyh(js)-dxh(is))*1d12/dxh(is)>1d0.or.dabs(dzh(ks)-dxh(is))*1d12/dxh(is)>1d0) then
       print *, "non-cubic cells", dxh(is),dyh(js),dzh(ks)
       stop
    endif
  end subroutine test_cell_size
!=================================================================================================

  subroutine get_staggered_fractions (cstag,d)
  use module_grid
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(out) :: cstag
  integer :: i,j,k,d
  integer :: i0,j0,k0

  call init_i0j0k0 (d,i0,j0,k0)

  work(:,:,:,1) = 0.d0
  work(:,:,:,2) = 0.d0
  do k=ks-1,ke+1
     do j=js-1,je+1
        do i=is-1,ie+1
            call get_half_fractions(cvof,d,i,j,k,work(i,j,k,1),work(i,j,k,2))
        enddo
     enddo
  enddo

  call do_all_ghost(work(:,:,:,1))
  call do_all_ghost(work(:,:,:,2))

  do k=ks-1,ke+1
     do j=js-1,je+1
        do i=is-1,ie+1
           cstag(i,j,k) = 0.5d0*(work(i,j,k,2)+work(i+i0,j+j0,k+k0,1))
        enddo
     enddo
  enddo

  call do_all_ghost(cstag)

  end subroutine get_staggered_fractions

end module module_vof
!=================================================================================================
!-------------------------------------------------------------------------------------------------
module module_output_vof
  use module_IO
  use module_flow
  use module_grid
  use module_solid
  use module_vof
  implicit none
  integer ,save :: vof_opened=0;
contains
  subroutine append_VOF_visit_file(rootname)
    implicit none
    character(*) :: rootname
    character(len=30), save :: file_name
    integer prank
    if(rank.ne.0) call pariserror('rank.ne.0 in append_VOF')
    if(vof_opened==0) then
       if(output_format==4) then
          file_name='fields.visit'
       else
          file_name='vof.visit'
       endif
       OPEN(UNIT=88,FILE=TRIM(file_name))
       write(88,10) nPdomain
10     format('!NBLOCKS ',I4)
       vof_opened=1
    else
       OPEN(UNIT=88,FILE=TRIM(file_name),position='append')
    endif
    do prank=0,NpDomain-1
       write(88,11) rootname//TRIM(int2text(prank,padding))//'.vtk'
11     format(A)
    enddo
    close(88)
  end subroutine  append_VOF_visit_file
!=================================================================================================
subroutine output_VOF(nf,i1,i2,j1,j2,k1,k2)
  implicit none
  integer :: nf,i1,i2,j1,j2,k1,k2
  if(output_format==2) call output_VOF_VTKSG(nf,i1,i2,j1,j2,k1,k2)
  if(output_format==3) call output_VOF_VTKSP(nf,i1,i2,j1,j2,k1,k2)
end subroutine output_VOF
!=================================================================================================
  subroutine output_VOF_VTKSP(nf,i1,i2,j1,j2,k1,k2)
    implicit none
    integer ::nf,i1,i2,j1,j2,k1,k2,i,j,k
    real(8) :: cfiltered
    character(len=30) :: rootname,filename
    rootname=trim(out_path)//'/VTK/VOF'//TRIM(int2text(nf,padding))//'-'
    if(rank==0) call append_VOF_visit_file(TRIM(rootname))

    OPEN(UNIT=8,FILE=TRIM(rootname)//TRIM(int2text(rank,padding))//'.vtk')
    write(8,10)
    write(8,11) time
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
    write(8,17)'VOF'
    write(8,18)
19  format('POINT_DATA ',I17)
17  format('SCALARS ',A20,' float 1')
18  format('LOOKUP_TABLE default')

    do k=k1,k2; do j=j1,j2; do i=i1,i2;
      if ( output_filtered_VOF ) then 
         cfiltered = b1*cvof(i,j,k) + & 
               b2*( cvof(i-1,j,k) + cvof(i,j-1,k) + cvof(i,j,k-1) + &
                    cvof(i+1,j,k) + cvof(i,j+1,k) + cvof(i,j,k+1) ) + &
               b3*( cvof(i+1,j+1,k) + cvof(i+1,j-1,k) + cvof(i-1,j+1,k) + cvof(i-1,j-1,k) + &
                    cvof(i+1,j,k+1) + cvof(i+1,j,k-1) + cvof(i-1,j,k+1) + cvof(i-1,j,k-1) + &
                    cvof(i,j+1,k+1) + cvof(i,j+1,k-1) + cvof(i,j-1,k+1) + cvof(i,j-1,k-1) ) + &
               b4*( cvof(i+1,j+1,k+1) + cvof(i+1,j+1,k-1) + cvof(i+1,j-1,k+1) + cvof(i+1,j-1,k-1) +  &
                    cvof(i-1,j+1,k+1) + cvof(i-1,j+1,k-1) + cvof(i-1,j-1,k+1) + cvof(i-1,j-1,k-1) )
         write(8,210) cfiltered
      else 
         write(8,210) cvof(i,j,k)
      end if ! output_filtered_VOF
    enddo; enddo; enddo
210 format(e14.5)
! 310 format(e14.5,e14.5,e14.5)
    close(8)
! TEMPORARY 
    if ( zip_data ) then 
      filename = TRIM(rootname)//TRIM(int2text(rank,padding))//'.vtk'
      call system('gzip '//trim(filename))
    end if ! zip_data
! END TEMPORARY 
end subroutine output_VOF_VTKSP
!=================================================================================================
  subroutine output_VOF_VTKSG(nf,i1,i2,j1,j2,k1,k2)
    implicit none
    integer ::nf,i1,i2,j1,j2,k1,k2,i,j,k
    real(8) :: cfiltered
    character(len=30) :: rootname,filename
    rootname=trim(out_path)//'/VTK/VOF'//TRIM(int2text(nf,padding))//'-'
    if(rank==0) call append_VOF_visit_file(TRIM(rootname))

    OPEN(UNIT=8,FILE=TRIM(rootname)//TRIM(int2text(rank,padding))//'.vtk')
    write(8,10)
    write(8,11) time
    write(8,12)
    write(8,13)
    write(8,14)i2-i1+1,j2-j1+1,k2-k1+1
    write(8,15)(i2-i1+1)*(j2-j1+1)*(k2-k1+1)
10  format('# vtk DataFile Version 2.0')
11  format('grid, time ',F16.8)
12  format('ASCII')
13  format('DATASET STRUCTURED_GRID')
14  format('DIMENSIONS ',I5,I5,I5)
15  format('POINTS ',I17,' double' )

    do k=k1,k2; do j=j1,j2; do i=i1,i2;
      write(8,320) x(i),y(j),z(k)
    enddo; enddo; enddo
320 format(e14.5,e14.5,e14.5)

    write(8,19)(i2-i1+1)*(j2-j1+1)*(k2-k1+1)
    write(8,17)'VOF'
    write(8,18)
19  format('POINT_DATA ',I17)
17  format('SCALARS ',A20,' float 1')
18  format('LOOKUP_TABLE default')

    do k=k1,k2; do j=j1,j2; do i=i1,i2
      if ( output_filtered_VOF ) then 
         cfiltered = b1*cvof(i,j,k) + & 
               b2*( cvof(i-1,j,k) + cvof(i,j-1,k) + cvof(i,j,k-1) + &
                    cvof(i+1,j,k) + cvof(i,j+1,k) + cvof(i,j,k+1) ) + &
               b3*( cvof(i+1,j+1,k) + cvof(i+1,j-1,k) + cvof(i-1,j+1,k) + cvof(i-1,j-1,k) + &
                    cvof(i+1,j,k+1) + cvof(i+1,j,k-1) + cvof(i-1,j,k+1) + cvof(i-1,j,k-1) + &
                    cvof(i,j+1,k+1) + cvof(i,j+1,k-1) + cvof(i,j-1,k+1) + cvof(i,j-1,k-1) ) + &
               b4*( cvof(i+1,j+1,k+1) + cvof(i+1,j+1,k-1) + cvof(i+1,j-1,k+1) + cvof(i+1,j-1,k-1) +  &
                    cvof(i-1,j+1,k+1) + cvof(i-1,j+1,k-1) + cvof(i-1,j-1,k+1) + cvof(i-1,j-1,k-1) )
         write(8,210) cfiltered
      else 
         write(8,210) cvof(i,j,k)
      end if ! output_filtered_VOF
    enddo; enddo; enddo
210 format(e14.5)
! 310 format(e14.5,e14.5,e14.5)
    close(8)
! TEMPORARY 
    if ( zip_data ) then 
      filename = TRIM(rootname)//TRIM(int2text(rank,padding))//'.vtk'
      call system('gzip '//trim(filename))
    end if ! zip_data
! END TEMPORARY 
end subroutine output_VOF_VTKSG
!=================================================================================================
!-------------------------------------------------------------------------------------------------
subroutine backup_VOF_write
  use module_freesurface
  implicit none
  integer ::i,j,k
  character(len=100) :: filename
  if (OrganizeOutFolder) then
    filename = trim(out_path)//'/BACKUP/backup_'//int2text(rank,padding)
  else 
    filename = trim(out_path)//'/backup_'//int2text(rank,padding)
  end if ! OrganizeOutFolder 
  call system('touch '//trim(filename)//'; mv '//trim(filename)//' '//trim(filename)//'.old')
  OPEN(UNIT=7,FILE=trim(filename),status='REPLACE',ACTION='write')
  !Note: p at ghost layers are needed for possion solver 
  write(7,1100)time,itimestep,imin,imax,jmin,jmax,kmin,kmax
  do k=kmin,kmax; do j=jmin,jmax; do i=imin,imax
    write(7,1200) u(i,j,k), v(i,j,k), w(i,j,k), p(i,j,k), cvof(i,j,k)
  enddo; enddo; enddo
  if (FreeSurface) then
     write(7,1201)R_RK,dR_RK,ddR_RK,V_0
  endif
  CLOSE(7)
  if(rank==0)print*,'Backup written at t=',time
  1100 FORMAT(es17.8e3,7I10)
  !Note: to guarantee identical results, 16 digits are needed for real8 
  1200 FORMAT(5es25.16e3)
  1201 FORMAT(4es25.16e3)
end subroutine backup_VOF_write
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
subroutine backup_VOF_read
  use module_freesurface
  implicit none
  integer ::i,j,k,i1,i2,j1,j2,k1,k2
  character(len=100) :: filename
  if (OrganizeOutFolder) then
    filename = trim(out_path)//'/BACKUP/backup_'//int2text(rank,padding)
  else 
    filename = trim(out_path)//'/backup_'//int2text(rank,padding)
  end if ! OrganizeOutFolder 
  OPEN(UNIT=7,FILE=trim(filename),status='old',action='read')
  read(7,*)time,itimestep,i1,i2,j1,j2,k1,k2
  if(i1/=imin .or. i2/=imax .or. j1/=jmin .or. j2/=jmax .or. k1/=kmin .or. k2/=kmax) &
    call pariserror("Error: backup_read")
  do k=kmin,kmax; do j=jmin,jmax; do i=imin,imax
    read(7,*) u(i,j,k), v(i,j,k), w(i,j,k), p(i,j,k), cvof(i,j,k)
  enddo; enddo; enddo
  if (FreeSurface) then
     read(7,*)R_RK,dR_RK,ddR_RK,V_0
  endif
  CLOSE(7)
end subroutine backup_VOF_read
!=================================================================================================
!-------------------------------------------------------------------------------------------------
  subroutine clean_debris()
    implicit none
    include 'mpif.h'
    integer, dimension(imin:imax,jmin:jmax,kmin:kmax) :: phase_id, clean_flag
    integer :: i,j,k,i0,j0,k0
    integer :: ni,nj,nk
    real(8) :: cvof_min,sum_cvof_min,sum1,sum_cvof
    integer :: m,n,l,inr,jnr,knr
    integer :: nposit
    integer :: cleaned_mix,cleaned_tot,itr,it_max,ierr
    real(8) :: d_vof, d_vof_tot
    logical :: all_clean_local, all_clean_global
   
    ! -1:gas, 1:liquid
    phase_id(:,:,:) = -1 
    ! -1:set to pure gas; 0:not changed; 1:set to pure liquid 
    clean_flag(:,:,:) = 0

    if ( clean_debris_method == 1 ) then 
       ! **********************************************************
       ! Stephane's approach
       ! **********************************************************
      
       ! Binary map of liquid and gas cells
       do i=imin,imax; do j=jmin,jmax; do k=kmin,kmax
         if ( cvof(i,j,k) > 0.5 ) phase_id(i,j,k) = 1
       end do; end do; end do

       ! flag debris
       ni=clean_debris_neighbours; nj=clean_debris_neighbours; nk=clean_debris_neighbours
       if (Nz == 2) nk = 0 !assume Nz=2 for 2d case 
       do i=is,ie; do j=js,je; do k=ks,ke
         clean_flag(i,j,k) = phase_id(i,j,k)
         do i0=-ni,ni; do j0=-nj,nj; do k0=-nk,nk
            if ( phase_id(i,j,k)*phase_id(i+i0,j+j0,k+k0) == -1 ) then
               clean_flag(i,j,k) = 0 
               exit
            end if ! phase_id
         enddo; enddo; enddo
       end do; end do; end do

    else if ( clean_debris_method == 2 ) then 
       ! **********************************************************
       ! Ruben's approach
       ! **********************************************************
       ! User defined parameter
       cvof_min      = 0.1d0 
       sum_cvof_min  = 0.3d0

       ni=clean_debris_neighbours;nj=clean_debris_neighbours;nk=clean_debris_neighbours
       if (Nz == 2) nk = 0 !assume Nz=2 for 2d case
       sum1 = dble((1+ni*2)*(1+nj*2)*(1+nk*2))
       do i=is,ie; do j=js,je; do k=ks,ke
         sum_cvof = 0.d0 
         do i0=-ni,ni; do j0=-nj,nj; do k0=-nk,nk
            sum_cvof = sum_cvof + cvof(i+i0,j+j0,k+k0) 
         enddo; enddo; enddo
         if ( cvof(i,j,k) <      cvof_min .and. & 
              sum_cvof    <  sum_cvof_min ) clean_flag(i,j,k) = -1 
         if ( cvof(i,j,k) > 1.d0-cvof_min .and. & 
              sum_cvof    > sum1-sum_cvof_min ) clean_flag(i,j,k) = 1 
       end do; end do; end do
    end if
    
! DEBUG
!!$    write(101,*) 'Cells to be cleaned at time',time,'for rank',rank
!!$    do i=is,ie; do j=js,je; do k=ks,ke
!!$      if ( clean_flag(i,j,k) == -1 .and. cvof(i,j,k) > 0.d0 ) then
!!$         write(101,*) '**********************************'
!!$         write(101,*) 'vof set to 0',i,j,k,time,cvof(i,j,k) 
!!$         write(101,'(9(E11.2))') cvof(i-1,j-1:j+1,k-1),cvof(i,j-1:j+1,k-1),cvof(i+1,j-1:j+1,k-1)
!!$         write(101,'(9(E11.2))') cvof(i-1,j-1:j+1,k  ),cvof(i,j-1:j+1,k  ),cvof(i+1,j-1:j+1,k  )
!!$         write(101,'(9(E11.2))') cvof(i-1,j-1:j+1,k+1),cvof(i,j-1:j+1,k+1),cvof(i+1,j-1:j+1,k+1)
!!$      else if ( clean_flag(i,j,k) ==  1 .and. cvof(i,j,k) < 1.d0 ) then 
!!$         write(101,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
!!$         write(101,*) 'vof set to 1',i,j,k,time,cvof(i,j,k) 
!!$         write(101,'(9(E11.2))') cvof(i-1,j-1:j+1,k-1),cvof(i,j-1:j+1,k-1),cvof(i+1,j-1:j+1,k-1)
!!$         write(101,'(9(E11.2))') cvof(i-1,j-1:j+1,k  ),cvof(i,j-1:j+1,k  ),cvof(i+1,j-1:j+1,k  )
!!$         write(101,'(9(E11.2))') cvof(i-1,j-1:j+1,k+1),cvof(i,j-1:j+1,k+1),cvof(i+1,j-1:j+1,k+1)
!!$      end if ! clean_flag 
!!$    end do; end do; end do
! END DEBUG

    if ( clean_debris_method == 1 .or. clean_debris_method == 2 ) then! remove debris
       do i=is,ie; do j=js,je; do k=ks,ke
          if ( clean_flag(i,j,k) == -1 .and. cvof(i,j,k) > 0.d0 ) then
             cvof(i,j,k) = 0.d0
          else if ( clean_flag(i,j,k) ==  1 .and. cvof(i,j,k) < 1.d0 ) then 
             cvof(i,j,k) = 1.d0
          end if ! clean_flag 
       end do; end do; end do
    endif
    
    if ( clean_debris_method == 3 ) then
       
       all_clean_local = .false.
       all_clean_global = .false.
       itr = 0
       cleaned_mix = 0

       do while ((.not.all_clean_global) .and. (itr<=10))
          itr = itr + 1
          do k=ks,ke; do j=js,je; do i=is,ie 
             if (vof_flag(i,j,k) == 2 ) then  ! mixed cell
                nposit=0
                do m=-1,1; do n=-1,1; do l=-1,1
                   inr=i+m
                   jnr=j+n
                   knr=k+l
                   if(vof_flag(inr,jnr,knr) == 2) then
                      nposit = nposit + 1
                   endif ! vof_flag
                enddo; enddo; enddo ! do m,n,l
                ! arrange coordinates so height direction is closest to normal
                ! try(:) array contains direction closest to normal first

                if(nposit<6) then
                   !clean cell
                   if (vof_phase(i,j,k)==1) then             
                      !calc mass diff
                      d_vof = d_vof + 1.0d0-cvof(i,j,k)
                      cvof(i,j,k) = 1.0d0
                      vof_flag(i,j,k) = 1
                   else
                      d_vof = d_vof - cvof(i,j,k)
                      cvof(i,j,k) = 0.0d0
                      vof_flag(i,j,k) = 0
                   endif
                   cleaned_mix = cleaned_mix + 1
                endif

             endif ! mixed cell loop
          enddo;enddo;enddo
          if (cleaned_mix==0) then
             all_clean_local=.true.
          endif
          call MPI_ALLREDUCE(all_clean_local, all_clean_global, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_Active, ierr)
          if (.not.all_clean_global) then
             call do_all_ghost(cvof)
             call get_flags_and_clip(cvof,vof_flag)
             call get_vof_phase(cvof)
          endif
       enddo ! all_clean
       call MPI_ALLREDUCE(cleaned_mix, cleaned_tot, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_Active, ierr)
       call MPI_ALLREDUCE(itr, it_max, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_Active, ierr)
       call MPI_ALLREDUCE(d_vof, d_vof_tot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_Active, ierr)
    endif
  end subroutine clean_debris

end module module_output_vof
