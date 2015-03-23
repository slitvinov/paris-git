!================================================================================================
!=================================================================================================
! Paris-0.1
! Extended from Code: FTC3D2011 (Front Tracking Code for 3D simulations)
! and Surfer. 
! 
! Authors: Sadegh Dabiri, Gretar Tryggvason
! Contact: sdabiri@gmail.com
! Author for Lagrangian particles extenstions: 
! Yue (Stanley) Ling (yueling@dalembert.upmc.fr), 
! Stephane Zaleski (zaleski@dalembert.upmc.fr) 
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
! module_lag_particle: Contains definition of variables for Lagrangian particle  
!   model and conversion between VOF and Lagrangian particle model
! Note: 1, Continuous phase is refered as fluid, and dispersed phase is refered
! as droplets or particles.
!       2, Droplets refered to dispersed phase resolved with Volume-of-Fluid
!       method (VOF); while Particles refered to dispersed phase represented by
!       Lagrangian Point-Particle Model (LPP).
!-------------------------------------------------------------------------------------------------
   module module_Lag_part
   use module_grid
   use module_BC
   use module_IO
   use module_tmpvar
   use module_VOF
   use module_2phase
   use module_flow
   use module_output_vof
   implicit none
   integer, dimension(:,:,:), allocatable :: tag_id
   integer, dimension(:,:,:), allocatable :: tag_flag
      ! 0 marked as untagged or reference fluid
      ! 1 marked as tagged droplet 
      ! 2 marked as S node
      ! 3 marked as C node
      ! 4 marked as reference fluid 
      ! 5 marked as ghost layer
   integer,parameter :: maxnum_diff_tag  = 50  ! ignore cases droplet spread over more than 1 block
   integer :: total_num_tag,totalnum_drop,totalnum_drop_indep,num_new_drop
   integer, dimension(:), allocatable :: num_drop
   integer, dimension(:), allocatable :: num_drop_merge
   integer, dimension(:), allocatable :: num_element
   integer, dimension(:), allocatable :: num_tag, tagmin, tagmax
   integer, dimension(:), allocatable :: tag_dropid
   integer, dimension(:), allocatable :: tag_rank
   integer, dimension(:), allocatable :: tag_mergeflag

   integer, dimension(:), allocatable :: new_drop_id
   real(8) :: tag_threshold=1d-9

   type element 
      real(8) :: xc,yc,zc,uc,vc,wc,duc,dvc,dwc,vol
      integer :: id 
   end type element
   type (element), dimension(:,:), allocatable :: element_stat

   type drop
      type(element) :: element
      integer :: num_cell_drop
      real(8) :: AspRatio
   end type drop
   type (drop), dimension(:), allocatable :: drops
   integer, dimension(:,:), allocatable :: drops_cell_list
   
   type drop_merge
      type(element) :: element
      integer :: num_cell_drop
      integer :: num_gcell
      integer :: num_diff_tag
      integer :: diff_tag_list(maxnum_diff_tag)
      integer :: flag_center_mass
   end type drop_merge
   type (drop_merge), dimension(:), allocatable :: drops_merge
   integer, dimension(:,:), allocatable :: drops_merge_cell_list
   integer, dimension(:,:,:), allocatable :: drops_merge_gcell_list

   type drop_merge_comm
      real(8) :: xc,yc,zc,uc,vc,wc,duc,dvc,dwc,vol
      integer :: id
      integer :: num_diff_tag 
      integer :: diff_tag_list(maxnum_diff_tag)
      integer :: flag_center_mass
   end type drop_merge_comm
   type (drop_merge_comm), dimension(:,:), allocatable :: drops_merge_comm

   logical :: LPP_initialized = .false.
   integer, dimension(:), allocatable :: num_part

   type particle
      type(element) :: element
      real(8) :: xcOld,ycOld,zcOld,ucOld,vcOld,wcOld,uf,vf,wf
      integer :: ic,jc,kc,tstepConvert,dummyint  
      ! Note: open_mpi sometimes failed to communicate the last varialbe in 
      !       MPI_TYPE_STRUCt correctly, dummyint is included to go around 
      !       this bug in mpi
   end type particle 
   type (particle), dimension(:,:), allocatable :: parts

   type particle_collect
      real(8) :: xc,yc,zc,hfx,hfy,hfz,vol
      integer :: ic,jc,kc,dummyint  
   end type particle_collect 
   type (particle_collect), dimension(:,:), allocatable :: parts_collect
   type (particle_collect), dimension(:), allocatable :: parts_collect_rank

   integer, dimension(:,:), allocatable :: parts_cross_id
   integer, dimension(:), allocatable :: parts_cross_id_rank
   integer, dimension(:,:), allocatable :: parts_cross_newrank
   integer, dimension(:), allocatable :: parts_cross_newrank_rank
   integer, dimension(:), allocatable :: num_part_cross

   ! substantial derivative of velocity
   real(8), dimension(:,:,:), allocatable :: sdu,sdv,sdw, & 
                                             sdu_work,sdv_work,sdw_work

   logical, dimension(:,:,:), allocatable :: RegAwayInterface

   integer, parameter :: CRAZY_INT = 3483129 
   real(8), parameter :: CRAZY_REAl = 123456789.123456789d-16
   real(8), parameter :: OneThird = 0.3333333333333333d0

   integer, parameter :: CriteriaRectangle = 1
   integer, parameter :: CriteriaCylinder  = 2
   integer, parameter :: CriteriaSphere    = 3
   integer, parameter :: CriteriaJet       = 4
   integer, parameter :: CriteriaInterface = 5
   integer, parameter :: CriteriaIsoDrop   = 6

   ! Stokes drag covers bubbles to particles limit. Finite Reynolds number
   ! extensions SN & CG are only for particles while MKL is only for bubbles
   integer, parameter :: dragmodel_Stokes = 1
   integer, parameter :: dragmodel_SN = 2     ! Schiller & Nauman
   integer, parameter :: dragmodel_CG = 3     ! Clift & Gauvin
   integer, parameter :: dragmodel_MKL = 4    ! Mei, Klausner,Lawrence 1994

   integer, parameter :: DropStatistics_PlotElement      = 1
   integer, parameter :: DropStatistics_ElementSizePDF   = 2
   integer, parameter :: DropStatistics_WriteElementData = 3
   integer, parameter :: DropStatistics_WriteElementDataEachRank = 4

   integer :: DropStatisticsMethod 
   logical :: DoConvertVOF2LPP 
   logical :: DoConvertLPP2VOF 
   integer :: dragmodel 
   integer :: ntimesteptag
   integer :: CriteriaConvertCase
   real(8) :: vol_cut, xlpp_min,ylpp_min,zlpp_min, & 
                       xlpp_max,ylpp_max,zlpp_max
   character(20) :: lppbdry_cond(6)

   integer :: max_num_drop, max_num_drop_merge 
   integer :: maxnum_cell_drop
   integer :: max_num_part
   integer :: max_num_part_cross

   integer :: outputlpp_format
   real(8) :: ConvertRegSizeToDiam
   
   integer :: nStepConverion

   real(8), parameter :: AspRatioSphere = 1.d0
   real(8) :: AspRatioTol

   logical :: WallEffectSettling
   logical :: output_lpp_evolution
   
   integer :: TwoWayCouplingFlag
   integer, parameter :: TwoWayCouplingIgnore      = 0
   integer, parameter :: TwoWayCouplingFilterForce = 1
   integer, parameter :: TwoWayCouplingFilterVel   = 2
   real(8) :: LengthLPP2dp, LengthLPP,mfLPP_ref

   logical :: UnsteadyPartForce

   integer :: SeedParticlesFlag
   logical :: SeedWithFluidVel
   logical :: AvoidSeedAtBlockBdry,AvoidSeedPartTooClose
   integer, parameter :: SeedParticlesNone      = 0
   integer, parameter :: SeedParticlesOneTime   = 1
   integer, parameter :: SeedParticlesRegular   = 2
   real(8) :: xmin_part_seed, ymin_part_seed, zmin_part_seed, dmin_part_seed, & 
              xmax_part_seed, ymax_part_seed, zmax_part_seed, dmax_part_seed, &
              umin_part_seed, vmin_part_seed, wmin_part_seed, & 
              umax_part_seed, vmax_part_seed, wmax_part_seed
   integer :: num_part_seed
   real(8) :: time_part_seed
   integer :: tag_opened=0
contains
!=================================================================================================
   subroutine initialize_LPP()

      implicit none
      call ReadLPPParameters

      allocate( parts(max_num_part,0:nPdomain-1) )
      allocate( num_part(0:nPdomain-1) )
      allocate( num_drop (0:nPdomain-1) )
      allocate( num_drop_merge (0:nPdomain-1) )
      allocate( num_element (0:nPdomain-1) )
      allocate( num_tag (0:nPdomain-1) )
      allocate( tagmin (0:nPdomain-1) )
      allocate( tagmax (0:nPdomain-1) )
      allocate( tag_flag(imin:imax,jmin:jmax,kmin:kmax) )
      allocate( tag_id  (imin:imax,jmin:jmax,kmin:kmax) )
      allocate( drops(max_num_drop) )
      allocate( drops_merge(max_num_drop_merge) )
      allocate( drops_cell_list(3,max_num_drop) )
      allocate( drops_merge_cell_list(3,max_num_drop_merge) )
      allocate( drops_merge_gcell_list(3,maxnum_cell_drop,max_num_drop_merge) )
      allocate( sdu(imin:imax,jmin:jmax,kmin:kmax), & 
                sdv(imin:imax,jmin:jmax,kmin:kmax), &
                sdw(imin:imax,jmin:jmax,kmin:kmax), & 
                sdu_work(imin:imax,jmin:jmax,kmin:kmax), & 
                sdv_work(imin:imax,jmin:jmax,kmin:kmax), &
                sdw_work(imin:imax,jmin:jmax,kmin:kmax) )
             
      if ( CriteriaConvertCase == CriteriaInterface ) & 
         allocate( RegAwayInterface(imin:imax,jmin:jmax,kmin:kmax) )

      if ( TwoWayCouplingFlag == TwoWayCouplingFilterForce) then  
         allocate( parts_collect(max_num_part,0:nPdomain-1) )
         allocate( parts_collect_rank(max_num_part) )
      end if 

      ! set default values
      num_tag  = 0
      num_part = 0
      num_drop = 0
      num_drop_merge = 0
      num_element = 0

      LPP_initialized = .true.

      sdu = 0.d0; sdv = 0.d0; sdw =0.d0
      sdu_work = 0.d0; sdv_work = 0.d0; sdw_work =0.d0

   end subroutine initialize_LPP

   subroutine ReadLPPParameters

      implicit none
      include 'mpif.h'

      integer ierr,in
      integer :: i
      logical file_is_there
      namelist /lppparameters/ DropStatisticsMethod, dragmodel, nTimeStepTag,   &
         DoConvertVOF2LPP,DoConvertLPP2VOF,CriteriaConvertCase, ConvertRegSizeToDiam, & 
         vol_cut,xlpp_min,xlpp_max,ylpp_min,ylpp_max,zlpp_min,zlpp_max,    &
         max_num_drop, max_num_drop_merge, maxnum_cell_drop, max_num_part, max_num_part_cross, &
         outputlpp_format,nStepConverion,AspRatioTol, &
         WallEffectSettling,output_lpp_evolution,lppbdry_cond,& 
         TwoWayCouplingFlag,LengthLPP2dp,mfLPP_ref,UnsteadyPartForce,&  
         SeedParticlesFlag, SeedWithFluidVel, & 
         AvoidSeedAtBlockBdry,AvoidSeedPartTooClose, & 
         num_part_seed, time_part_seed, & 
         xmin_part_seed, ymin_part_seed, zmin_part_seed, dmin_part_seed, & 
         xmax_part_seed, ymax_part_seed, zmax_part_seed, dmax_part_seed, &
         umin_part_seed, vmin_part_seed, wmin_part_seed, & 
         umax_part_seed, vmax_part_seed, wmax_part_seed

      in=32

      ! Set default values 
      DropStatisticsMethod = 0 
      dragmodel    = 1
      nTimeStepTag = 1
      CriteriaConvertCase = 1
      vol_cut  = 1.d-9
      xlpp_min = xh(   Ng)
      xlpp_max = xh(Nx+Ng)
      ylpp_min = yh(   Ng)
      ylpp_max = yh(Ny+Ng)
      zlpp_min = zh(   Ng)
      zlpp_max = zh(Nz+Ng)
      max_num_drop = 100
      max_num_drop_merge = 20
      maxnum_cell_drop = 10000
      max_num_part = 100
      max_num_part_cross = 10
      outputlpp_format = 3
      ConvertRegSizeToDiam = 2.d0 
      nStepConverion = 0
      AspRatioTol = 1.5d0
      WallEffectSettling = .false.
      output_lpp_evolution = .false.
      lppbdry_cond=['undefined','undefined','undefined','undefined','undefined','undefined']
      TwoWayCouplingFlag = 0 
      LengthLPP2dp = 8.d0
      mfLPP_ref = 0.5d0
      UnsteadyPartForce = .false.
      SeedParticlesFlag = 0
      SeedWithFluidVel =.false. 
      AvoidSeedAtBlockBdry =.false. 
      AvoidSeedPartTooClose =.false. 
      num_part_seed = 0; time_part_seed = 0.0;
      xmin_part_seed=0.0; ymin_part_seed=0.0; zmin_part_seed=0.0; dmin_part_seed=0.0 
      xmax_part_seed=0.0; ymax_part_seed=0.0; zmax_part_seed=0.0; dmax_part_seed=0.0
      umin_part_seed=0.0; vmin_part_seed=0.0; wmin_part_seed=0.0
      umax_part_seed=0.0; vmax_part_seed=0.0; wmax_part_seed=0.0; 

      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
      inquire(file='inputlpp',exist=file_is_there)
      open(unit=in, file='inputlpp', status='old', action='read', iostat=ierr)
      if (file_is_there) then
         if(ierr == 0) then
            read(UNIT=in,NML=lppparameters)
            if(rank==0) write(out,*)'Largranian point-particle parameters read successfully'
         else
            print *, 'rank=',rank,' has error ',ierr,' opening file inputlpp'
         endif
      else
         if (rank == 0) call pariserror("ReadLPPParameters: no 'inputlpp' file.")
      endif
      close(in)
      if (rank == 0) then
         write(UNIT=out,NML=lppparameters)
      end if ! rank

      ! Check consistence of inputs
      if ( (dragmodel == dragmodel_CG .or. dragmodel == dragmodel_SN) .and. mu1>mu2 ) & 
         call pariserror("Particle drag law is used for bubbles!")
      if ( (dragmodel == dragmodel_MKL) .and. mu1<mu2 ) & 
         call pariserror("Bubble drag law is used for particles!")
      if ( MOD(nstats,nTimeStepTag) /= 0 ) & 
         call pariserror("nstats must be integer times fo nTimeStepTag!")

      do i=1,6
         if(lppbdry_cond(i) == 'undefined') then 
            lppbdry_cond(i) = 'periodic'
            if ( rank == 0 ) write(*,*) 'undefined lpp condition in direction', i,' now is set to periodic by default!' 
         end if ! lppbdry_cond 
      end do !i

   end subroutine ReadLPPParameters


   subroutine lppsweeps(tswap,time,iStage)
      implicit none
    
      integer, intent(in) :: tswap,iStage
      real(8), intent(in) :: time

      call ComputePartForce(tswap,time)
      if ( TwoWayCouplingFlag /= TwoWayCouplingIgnore ) & 
         call TwoWayCouplingForce(TwoWayCouplingFlag) 
      call UpdatePartSol(iStage)

   end subroutine lppsweeps

   subroutine lppvofsweeps(tswap,time)
      implicit none
      include 'mpif.h'
    
      integer, intent(in) :: tswap
      real(8), intent(in) :: time

      ! Only do tagging and conversion in specific time steps
      if ( MOD(tswap,ntimestepTag) == 0 ) then 
         call tag_drop()
         if ( nPdomain > 1 ) call tag_drop_all
         call CreateTag2DropTable
         if ( nPdomain > 1 ) call merge_drop_pieces 

         if ( DropStatisticsMethod > 0 .and. & 
              MOD(tswap,nstats) == 0 ) call drop_statistics(tswap,time) 

         if ( (DoConvertVOF2LPP .or. DoConvertLPP2VOF) .and. & 
               CriteriaConvertCase == CriteriaInterface )  call MarkRegAwayInterface()


         if ( DoConvertVOF2LPP ) call ConvertVOF2LPP(tswap)
         if ( DoConvertLPP2VOF ) call ConvertLPP2VOF(tswap)   
         call ReleaseTag2DropTable
         if (DoConvertVOF2LPP .or. DoConvertLPP2VOF) then 
            if ( nPdomain > 1 ) & 
               call VOFCommGhost    ! Communicate vof field after conversion
            call linfunc(rho,rho1,rho2,DensMean)
            call linfunc(mu,mu1,mu2,ViscMean)
         end if ! DoConvertVOF2LPP,DoConvertLPP2VOF 
         
         if ( rank == 0 ) write(*,*) 'Total num of particles/elements:',time,sum(num_part),sum(num_element)
      end if ! tswap

!      if ( MOD(tswap,termout) == 0 ) then 
!         if ( rank == 0 ) write(*,*) 'Numbers of drops/particles:', time,sum(num_element)-sum(num_part),sum(num_part)
!      end if ! tswap

   end subroutine lppvofsweeps

  subroutine tag_drop()
    implicit none
    include 'mpif.h'
    integer :: i,j,k, i0,j0,k0
    integer :: isq,jsq,ksq
    integer :: current_id
    integer :: s_queue(Nx*Ny,3),c_queue(Nx*Ny,3) ! record i,j,k
    integer :: ns_queue,is_queue,nc_queue
    
    real(8) :: volcell,cvof_scaled

    integer :: num_cell_drop,cell_list(3)
    real(8) :: xc,yc,zc,uc,vc,wc,duc,dvc,dwc,vol 

    logical :: merge_drop

    integer :: imin_drop,imax_drop,jmin_drop,jmax_drop,kmin_drop,kmax_drop
    integer :: idif_drop,jdif_drop,kdif_drop

    if (.not. LPP_initialized) then 
      call initialize_LPP()
      LPP_initialized = .true.
    end if ! LPP_initialized

    tag_id  (:,:,:) = 0
    tag_flag(:,:,:) = 0
    current_id = 1

    drops_merge(:)%num_gcell = 0

    ns_queue = 0
    s_queue(:,:) = 0
    num_drop(:) = 0
    num_drop_merge(:) = 0
    do i=imin,imax; do j=jmin,jmax; do k=kmin,kmax
    !do i=is,ie; do j=js,je; do k=ks,ke
      if ( cvof(i,j,k) > tag_threshold .and. tag_flag(i,j,k) == 0 ) then 
        tag_id  (i,j,k) = current_id
        tag_flag(i,j,k) = 2 ! mark as S node
        num_drop(rank) = num_drop(rank) + 1
        if ( num_drop(rank) > max_num_drop ) & 
           call lpperror("Drop number exceeds maximum number!") 
        ! put the present node into S queue
        ns_queue = ns_queue + 1 
        s_queue(ns_queue,1) = i
        s_queue(ns_queue,2) = j
        s_queue(ns_queue,3) = k
  
        vol = 0.d0
        xc  = 0.d0
        yc  = 0.d0
        zc  = 0.d0
        uc  = 0.d0
        vc  = 0.d0
        wc  = 0.d0
        duc  = 0.d0
        dvc  = 0.d0
        dwc  = 0.d0
        num_cell_drop = 0
        cell_list = 0

        merge_drop = .false.
        imin_drop = imax
        imax_drop = imin
        jmin_drop = jmax
        jmax_drop = jmin
        kmin_drop = kmax
        kmax_drop = kmin
        do while ( ns_queue > 0 )
          nc_queue = 0 
          c_queue(:,:) = 0
          do is_queue=1,ns_queue
            isq = s_queue(is_queue,1)
            jsq = s_queue(is_queue,2)
            ksq = s_queue(is_queue,3)
       
            do i0=-1,1; do j0=-1,1; do k0=-1,1
               if ( isq+i0 >= imin .and. isq+i0 <= imax .and. & 
                    jsq+j0 >= jmin .and. jsq+j0 <= jmax .and. &
                    ksq+k0 >= kmin .and. ksq+k0 <= kmax ) then  
                  if ( cvof(isq+i0,jsq+j0,ksq+k0)      > tag_threshold .and. & 
                      tag_flag(isq+i0,jsq+j0,ksq+k0) == 0 ) then 
                     tag_id  (isq+i0,jsq+j0,ksq+k0) = current_id  ! tag node with id
                     tag_flag(isq+i0,jsq+j0,ksq+k0) = 3  ! mark as C node
                     ! put current node into C queue
                     nc_queue = nc_queue + 1
                     c_queue(nc_queue,1) = isq+i0
                     c_queue(nc_queue,2) = jsq+j0
                     c_queue(nc_queue,3) = ksq+k0
                  end if 
               end if ! isq
            enddo;enddo;enddo ! i0,j0,k0

            if (  isq >= is .and. isq <= ie .and. &     ! internal cells 
                  jsq >= js .and. jsq <= je .and. & 
                  ksq >= ks .and. ksq <= ke ) then
               tag_flag(isq,jsq,ksq) = 1 ! mark S node as tagged
               ! perform droplet calculation
               volcell = dx(isq)*dy(jsq)*dz(ksq)
               cvof_scaled = cvof(isq,jsq,ksq)*volcell
               vol = vol + cvof_scaled
               xc  = xc  + cvof_scaled*x(isq)
               yc  = yc  + cvof_scaled*y(jsq)
               zc  = zc  + cvof_scaled*z(ksq)
               uc  = uc  + cvof_scaled*u(isq,jsq,ksq)
               vc  = vc  + cvof_scaled*v(isq,jsq,ksq)
               wc  = wc  + cvof_scaled*w(isq,jsq,ksq)
               duc  = duc  + cvof_scaled*sdu(isq,jsq,ksq)
               dvc  = dvc  + cvof_scaled*sdv(isq,jsq,ksq)
               dwc  = dwc  + cvof_scaled*sdw(isq,jsq,ksq)
               if ( num_cell_drop == 0 ) then
                  num_cell_drop = num_cell_drop + 1
                  if ( num_cell_drop == 1) cell_list(1:3) = [isq,jsq,ksq]
               end if ! num_cell_drop
            else  ! ghost cells
               tag_flag(isq,jsq,ksq) = 5
               if ( merge_drop .eqv. .false.) then 
                  merge_drop = .true.
                  num_drop      (rank) = num_drop      (rank) - 1
                  num_drop_merge(rank) = num_drop_merge(rank) + 1
                  if ( num_drop_merge(rank) > max_num_drop_merge ) & 
                     call lpperror("Drop merge number exceeds maximum number!") 
               end if ! merge_drop
               if ( drops_merge(num_drop_merge(rank))%num_gcell < maxnum_cell_drop) then 
                  drops_merge(num_drop_merge(rank))%num_gcell = drops_merge(num_drop_merge(rank))%num_gcell + 1
                  drops_merge_gcell_list(1,drops_merge(num_drop_merge(rank))%num_gcell,num_drop_merge(rank)) = isq
                  drops_merge_gcell_list(2,drops_merge(num_drop_merge(rank))%num_gcell,num_drop_merge(rank)) = jsq
                  drops_merge_gcell_list(3,drops_merge(num_drop_merge(rank))%num_gcell,num_drop_merge(rank)) = ksq
               else
                  write(*,*) 'Warning: ghost cell number of tag',current_id,'at rank',rank,'reaches max value!'
               end if ! drops_merge(num_drop_merge(rank))%num_gcell
            end if ! isq, jsq, ksq

            imin_drop = MIN(isq,imin_drop)
            imax_drop = MAX(isq,imax_drop)
            jmin_drop = MIN(jsq,jmin_drop)
            jmax_drop = MAX(jsq,jmax_drop)
            kmin_drop = MIN(ksq,kmin_drop)
            kmax_drop = MAX(ksq,kmax_drop)
            idif_drop = imax_drop-imin_drop+1
            jdif_drop = jmax_drop-jmin_drop+1
            kdif_drop = kmax_drop-kmin_drop+1
          end do ! is_queue
          ! mark all C nodes as S nodes
          if ( nc_queue >= 0 ) then 
            s_queue(:,:) = c_queue(:,:)   ! mark all C nodes as S nodes
            ns_queue = nc_queue
          end if ! nc_queue
        end do ! ns_queue>0

        if ( vol > 0.d0 ) then  
          if ( merge_drop ) then
            drops_merge(num_drop_merge(rank))%element%id  = current_id
            drops_merge(num_drop_merge(rank))%element%vol = vol 
            drops_merge(num_drop_merge(rank))%element%xc  = xc/vol 
            drops_merge(num_drop_merge(rank))%element%yc  = yc/vol 
            drops_merge(num_drop_merge(rank))%element%zc  = zc/vol 
            drops_merge(num_drop_merge(rank))%element%uc  = uc/vol 
            drops_merge(num_drop_merge(rank))%element%vc  = vc/vol
            drops_merge(num_drop_merge(rank))%element%wc  = wc/vol
            drops_merge(num_drop_merge(rank))%element%duc  = duc/vol 
            drops_merge(num_drop_merge(rank))%element%dvc  = dvc/vol
            drops_merge(num_drop_merge(rank))%element%dwc  = dwc/vol
            drops_merge(num_drop_merge(rank))%num_cell_drop = num_cell_drop
            drops_merge_cell_list(:,num_drop_merge(rank)) = cell_list
          else 
            drops(num_drop(rank))%element%id  = current_id
            drops(num_drop(rank))%element%vol = vol 
            drops(num_drop(rank))%element%xc  = xc/vol
            drops(num_drop(rank))%element%yc  = yc/vol
            drops(num_drop(rank))%element%zc  = zc/vol
            drops(num_drop(rank))%element%uc  = uc/vol
            drops(num_drop(rank))%element%vc  = vc/vol
            drops(num_drop(rank))%element%wc  = wc/vol
            drops(num_drop(rank))%element%duc  = duc/vol
            drops(num_drop(rank))%element%dvc  = dvc/vol
            drops(num_drop(rank))%element%dwc  = dwc/vol
            drops(num_drop(rank))%num_cell_drop = num_cell_drop
            drops_cell_list(:,num_drop(rank)) = cell_list

            drops(num_drop(rank))%AspRatio = DBLE(MAX(idif_drop,jdif_drop,kdif_drop)) & 
                                                 /DBLE(MIN(idif_drop,jdif_drop,kdif_drop))
          end if ! merge_drop
          current_id = current_id+1
        else ! no need to merge droplet piece only contains ghost cells  
          if ( merge_drop ) then
             drops_merge(num_drop_merge(rank))%num_gcell = 0
             num_drop_merge(rank) = num_drop_merge(rank) - 1
          else 
             num_drop      (rank) = num_drop      (rank) - 1
          end if ! merge_drop
        end if ! vol
      else if ( cvof(i,j,k) == 0.d0 .and. tag_flag(i,j,k) == 0 ) then 
         tag_flag(i,j,k) = 4
      end if ! cvof(i,j,k)
    enddo; enddo; enddo

    num_tag(rank) = num_drop(rank) + num_drop_merge(rank)
    if (npDomain == 1 ) total_num_tag = num_tag(rank)

   end subroutine tag_drop

   subroutine tag_drop_all()
      implicit none
      include 'mpif.h'
      integer :: i,j,k
      integer :: ierr,irank,num_tag_accu
      integer :: req(4),sta(MPI_STATUS_SIZE,4)
      integer , parameter :: ngh=2
      integer , parameter :: root=0
      integer :: num_drop_rank, num_drop_merge_rank

      ! Broadcast num_drop(rank) to all processes
      num_drop_rank       = num_drop      (rank)
      num_drop_merge_rank = num_drop_merge(rank)
      call MPI_ALLGATHER(num_drop_rank       , 1, MPI_INTEGER, &
                         num_drop      (:)   , 1, MPI_INTEGER, MPI_Comm_World, ierr)
      call MPI_ALLGATHER(num_drop_merge_rank , 1, MPI_INTEGER, &
                         num_drop_merge(:)   , 1, MPI_INTEGER, MPI_Comm_World, ierr)
      num_tag(:) = num_drop(:) + num_drop_merge(:)
      total_num_tag =  sum(num_tag)
      tagmin(0) = min(1,num_tag(0))
      tagmax(0) = num_tag(0)
      do irank = 1,nPdomain-1
         tagmin(irank) = sum(num_tag(0:irank-1)) + min(1,num_tag(irank))
         tagmax(irank) = tagmin(irank) + max(num_tag(irank)-1,0)
      end do ! irank

      ! update tag_id from local to global id (Note: no need to change domain 0)
      if ( rank > 0 ) then 
         num_tag_accu = sum(num_tag(0:rank-1),1)
         do i=is,ie; do j=js,je; do k=ks,ke
            if ( tag_id(i,j,k) > 0 ) tag_id(i,j,k)=tag_id(i,j,k)+num_tag_accu
         end do; end do; end do
      end if ! rank

      call ighost_x(tag_id,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
      call ighost_y(tag_id,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
      call ighost_z(tag_id,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)

   end subroutine tag_drop_all

   subroutine merge_drop_pieces
      implicit none

      integer :: idrop,iCell,idrop1
      integer :: idiff_tag,tag,tag1,idiff_tag1,idiff_tag2,tag2
      real(8) :: vol_merge,xc_merge,yc_merge,zc_merge,uc_merge,vc_merge,wc_merge,vol1, &
                                                      duc_merge,dvc_merge,dwc_merge
      integer :: irank, irank1
      real(8) :: max_drop_merge_vol
      integer :: tag_max_drop_merge_vol
      logical :: TagAlreadyListed,FinishUpdateTag

      ! Check ghost cells of droplet pieces
      if ( num_drop_merge(rank) > 0 ) then 
      do idrop = 1, num_drop_merge(rank) 
         drops_merge(idrop)%num_diff_tag = 1 
         drops_merge(idrop)%diff_tag_list(drops_merge(idrop)%num_diff_tag) &
            = tag_id(drops_merge_gcell_list(1,1,idrop), &
                     drops_merge_gcell_list(2,1,idrop), &
                     drops_merge_gcell_list(3,1,idrop))
         if ( drops_merge(idrop)%num_gcell > 1 ) then 
            ! Go through the ghost cells and list all different tags
            do iCell = 2,drops_merge(idrop)%num_gcell
               ! Check if the tag of the current cell is already list 
               TagAlreadyListed = .false.
               do idiff_tag = 1, drops_merge(idrop)%num_diff_tag
                  if ( tag_id(drops_merge_gcell_list(1,iCell,idrop), & 
                              drops_merge_gcell_list(2,iCell,idrop), & 
                              drops_merge_gcell_list(3,iCell,idrop)) &
                    == drops_merge(idrop)%diff_tag_list(idiff_tag)) then 
                     TagAlreadyListed = .true.
                     exit
                  end if ! tag_id
               end do ! idiff_tag
               
               if ( TagAlreadyListed ) then 
                  cycle 
               else  ! a new differet tag is found
                  drops_merge(idrop)%num_diff_tag = &
                  drops_merge(idrop)%num_diff_tag + 1
                  if ( drops_merge(idrop)%num_diff_tag <= maxnum_diff_tag ) then  
!                     call lpperror("Number of different tags of a droplet pieces exceeds the max number!") 
                     drops_merge(idrop)%diff_tag_list(drops_merge(idrop)%num_diff_tag) &
                        = tag_id(drops_merge_gcell_list(1,iCell,idrop), &
                                 drops_merge_gcell_list(2,iCell,idrop), &
                                 drops_merge_gcell_list(3,iCell,idrop))
                  else
                     write(*,*) "Number of diff tags exceeds the max number!",time,rank,idrop 
                     call lpperror("Number of diff tags exceeds the max number!") 
                  end if ! maxnum_diff_tag
               end if ! idiff_tag
            end do ! iCell
         end if ! drops_merge(idrop)%num_gcell
      end do ! idrop
      end if ! num_drop_merge(rank)

      call CollectDropMerge

      ! merge droplets pieces & calculate droplet properties
      if ( rank == 0 ) then 
         ! Find complete list of different tags belonging to a drop
         FinishUpdateTag = .false. 
         do while (.not.FinishUpdateTag)
            FinishUpdateTag = .true. 
            do tag = 1,total_num_tag
               if ( tag_mergeflag(tag) == 0 ) cycle
               idrop = tag_dropid(tag)
               irank = tag_rank  (tag)
               do idiff_tag = 1,drops_merge_comm(idrop,irank)%num_diff_tag
                  tag1   = drops_merge_comm(idrop,irank)%diff_tag_list(idiff_tag)
                  idrop1 = tag_dropid(tag1)
                  irank1 = tag_rank  (tag1)
                  do idiff_tag1 = 1,drops_merge_comm(idrop1,irank1)%num_diff_tag
                     TagAlreadylisted = .false.
                     tag2   = drops_merge_comm(idrop1,irank1)%diff_tag_list(idiff_tag1)
                     if ( tag2 == tag ) TagAlreadyListed = .true.  
                     do idiff_tag2 = 1,drops_merge_comm(idrop,irank)%num_diff_tag
                        if ( tag2 == drops_merge_comm(idrop,irank)%diff_tag_list(idiff_tag2) ) & 
                           TagAlreadylisted = .true.
                     end do ! idiff_tag2
                     if (.not. TagAlreadylisted) then 
                        drops_merge_comm(idrop,irank)%num_diff_tag = & 
                           drops_merge_comm(idrop,irank)%num_diff_tag + 1
                        drops_merge_comm(idrop,irank)%diff_tag_list( & 
                           drops_merge_comm(idrop,irank)%num_diff_tag) = tag2
                        if ( FinishUpdateTag ) FinishUpdateTag = .false. 
                     end if ! TagAlreadylist
                  end do ! idiff_tag1
               end do ! idff_tag
            end do ! tag
         end do ! update_tag_finish

         ! Calculate merge drop properties by summing all pieces
         allocate( new_drop_id(1:total_num_tag) )
         new_drop_id(:) = 0
         num_new_drop = 0
         drops_merge_comm(:,:)%flag_center_mass = 0 
         do tag = 1,total_num_tag
            if ( new_drop_id(tag) == 0 .and. tag_mergeflag(tag) == 1 ) then 
               num_new_drop = num_new_drop + 1
               new_drop_id(tag) = num_new_drop
               idrop = tag_dropid(tag)
               irank = tag_rank  (tag)
               vol1   = drops_merge_comm(idrop,irank)%vol
               vol_merge = drops_merge_comm(idrop,irank)%vol
               xc_merge  = drops_merge_comm(idrop,irank)%xc*vol1
               yc_merge  = drops_merge_comm(idrop,irank)%yc*vol1
               zc_merge  = drops_merge_comm(idrop,irank)%zc*vol1
               uc_merge  = drops_merge_comm(idrop,irank)%uc*vol1
               vc_merge  = drops_merge_comm(idrop,irank)%vc*vol1
               wc_merge  = drops_merge_comm(idrop,irank)%wc*vol1
               duc_merge  = drops_merge_comm(idrop,irank)%duc*vol1
               dvc_merge  = drops_merge_comm(idrop,irank)%dvc*vol1
               dwc_merge  = drops_merge_comm(idrop,irank)%dwc*vol1

               max_drop_merge_vol = drops_merge_comm(idrop,irank)%vol
               tag_max_drop_merge_vol = tag
               do idiff_tag = 1,drops_merge_comm(idrop,irank)%num_diff_tag
                  tag1   = drops_merge_comm(idrop,irank)%diff_tag_list(idiff_tag)
                  new_drop_id(tag1) = num_new_drop 
                  idrop1 = tag_dropid(tag1)
                  irank1 = tag_rank  (tag1)
                  vol1   = drops_merge_comm(idrop1,irank1)%vol
                  vol_merge = vol_merge + vol1
                  xc_merge  = xc_merge  + drops_merge_comm(idrop1,irank1)%xc*vol1
                  yc_merge  = yc_merge  + drops_merge_comm(idrop1,irank1)%yc*vol1
                  zc_merge  = zc_merge  + drops_merge_comm(idrop1,irank1)%zc*vol1
                  uc_merge  = uc_merge  + drops_merge_comm(idrop1,irank1)%uc*vol1
                  vc_merge  = vc_merge  + drops_merge_comm(idrop1,irank1)%vc*vol1
                  wc_merge  = wc_merge  + drops_merge_comm(idrop1,irank1)%wc*vol1
                  duc_merge  = duc_merge  + drops_merge_comm(idrop1,irank1)%duc*vol1
                  dvc_merge  = dvc_merge  + drops_merge_comm(idrop1,irank1)%dvc*vol1
                  dwc_merge  = dwc_merge  + drops_merge_comm(idrop1,irank1)%dwc*vol1
                  if (drops_merge_comm(idrop1,irank1)%vol > max_drop_merge_vol) then
                     max_drop_merge_vol = drops_merge_comm(idrop1,irank1)%vol
                     tag_max_drop_merge_vol = tag1
                  end if ! max_drop_merge_vol
               end do ! idiff_tag
               xc_merge = xc_merge/vol_merge
               yc_merge = yc_merge/vol_merge
               zc_merge = zc_merge/vol_merge
               uc_merge = uc_merge/vol_merge
               vc_merge = vc_merge/vol_merge
               wc_merge = wc_merge/vol_merge
               duc_merge = duc_merge/vol_merge
               dvc_merge = dvc_merge/vol_merge
               dwc_merge = dwc_merge/vol_merge

               idrop1 = tag_dropid(tag_max_drop_merge_vol)
               irank1 = tag_rank  (tag_max_drop_merge_vol)
               drops_merge_comm(idrop1,irank1)%flag_center_mass = 1

               drops_merge_comm(idrop,irank)%vol = vol_merge
               drops_merge_comm(idrop,irank)%xc  = xc_merge
               drops_merge_comm(idrop,irank)%yc  = yc_merge
               drops_merge_comm(idrop,irank)%zc  = zc_merge
               drops_merge_comm(idrop,irank)%uc  = uc_merge
               drops_merge_comm(idrop,irank)%vc  = vc_merge
               drops_merge_comm(idrop,irank)%wc  = wc_merge
               drops_merge_comm(idrop,irank)%duc  = duc_merge
               drops_merge_comm(idrop,irank)%dvc  = dvc_merge
               drops_merge_comm(idrop,irank)%dwc  = dwc_merge
               do idiff_tag = 1,drops_merge_comm(idrop,irank)%num_diff_tag
                  tag1   = drops_merge_comm(idrop,irank)%diff_tag_list(idiff_tag)
                  idrop1 = tag_dropid(tag1)
                  irank1 = tag_rank  (tag1)
                  drops_merge_comm(idrop1,irank1)%vol = vol_merge
                  drops_merge_comm(idrop1,irank1)%xc  = xc_merge
                  drops_merge_comm(idrop1,irank1)%yc  = yc_merge
                  drops_merge_comm(idrop1,irank1)%zc  = zc_merge
                  drops_merge_comm(idrop1,irank1)%uc  = uc_merge
                  drops_merge_comm(idrop1,irank1)%vc  = vc_merge
                  drops_merge_comm(idrop1,irank1)%wc  = wc_merge
                  drops_merge_comm(idrop1,irank1)%duc  = duc_merge
                  drops_merge_comm(idrop1,irank1)%dvc  = dvc_merge
                  drops_merge_comm(idrop1,irank1)%dwc  = dwc_merge
               end do ! idiff_tag
            end if ! newdropid(tag)
         end do ! tag
         totalnum_drop_indep = sum(num_drop)
         new_drop_id = new_drop_id + totalnum_drop_indep
         totalnum_drop = totalnum_drop_indep + num_new_drop
      end if ! rank 

      call DistributeDropMerge
   
      ! finalize
      if ( rank == 0 ) deallocate( new_drop_id )

   end subroutine merge_drop_pieces

   subroutine drop_statistics(tswap,time)
      implicit none
      
      include 'mpif.h'
      integer, intent(in) :: tswap
      real(8), intent(in) :: time
      integer :: req(2),sta(MPI_STATUS_SIZE,2),MPI_Comm,ireq,ierr
      integer :: MPI_element_type, oldtypes(0:1), blockcounts(0:1), & 
                 offsets(0:1), extent,r8extent 
      integer :: maxnum_element, total_num_element
      integer :: irank, idrop, ielement, ielem_plot
      integer :: num_element_estimate(0:nPdomain-1)
      integer :: num_element_rank

      integer, parameter :: num_gaps = 1000
      real(8) :: gap,dmax,d
      integer :: igap,count_element(num_gaps)
      logical :: file_exist
      character(len=30) :: filename 
      type(element) :: element_NULL
      element_NULL%xc = 0.d0;element_NULL%yc = 0.d0;element_NULL%zc = 0.d0
      element_NULL%uc = 0.d0;element_NULL%vc = 0.d0;element_NULL%wc = 0.d0
      element_NULL%duc = 0.d0;element_NULL%dvc = 0.d0;element_NULL%dwc = 0.d0
      element_NULL%vol = 0.d0;element_NULL%id = CRAZY_INT

      num_element_estimate = num_drop+num_drop_merge+num_part
      maxnum_element = maxval(num_element_estimate) 
      allocate( element_stat(maxnum_element,0:nPdomain-1) )

      !  initialize element_stat
      num_element(:) = 0
      element_stat(:,:) = element_NULL
      if ( num_drop(rank) > 0 ) then 
         do idrop = 1,num_drop(rank) 
            num_element(rank) = num_element(rank) + 1
            element_stat(num_element(rank),rank)%vol = drops(idrop)%element%vol
            element_stat(num_element(rank),rank)%xc  = drops(idrop)%element%xc
            element_stat(num_element(rank),rank)%yc  = drops(idrop)%element%yc
            element_stat(num_element(rank),rank)%zc  = drops(idrop)%element%zc
            element_stat(num_element(rank),rank)%uc  = drops(idrop)%element%uc
            element_stat(num_element(rank),rank)%vc  = drops(idrop)%element%vc
            element_stat(num_element(rank),rank)%wc  = drops(idrop)%element%wc
            element_stat(num_element(rank),rank)%duc  = drops(idrop)%element%duc
            element_stat(num_element(rank),rank)%dvc  = drops(idrop)%element%dvc
            element_stat(num_element(rank),rank)%dwc  = drops(idrop)%element%dwc
            element_stat(num_element(rank),rank)%id  = drops(idrop)%element%id
         end do ! idrop
      end if ! num_drop(rank)
      
      if ( num_drop_merge(rank) > 0 ) then 
         do idrop = 1,num_drop_merge(rank)
            if ( drops_merge(idrop)%flag_center_mass == 1 ) then  
               num_element(rank) = num_element(rank) + 1
               element_stat(num_element(rank),rank)%vol = drops_merge(idrop)%element%vol
               element_stat(num_element(rank),rank)%xc  = drops_merge(idrop)%element%xc
               element_stat(num_element(rank),rank)%yc  = drops_merge(idrop)%element%yc
               element_stat(num_element(rank),rank)%zc  = drops_merge(idrop)%element%zc
               element_stat(num_element(rank),rank)%uc  = drops_merge(idrop)%element%uc
               element_stat(num_element(rank),rank)%vc  = drops_merge(idrop)%element%vc
               element_stat(num_element(rank),rank)%wc  = drops_merge(idrop)%element%wc
               element_stat(num_element(rank),rank)%duc  = drops_merge(idrop)%element%duc
               element_stat(num_element(rank),rank)%dvc  = drops_merge(idrop)%element%dvc
               element_stat(num_element(rank),rank)%dwc  = drops_merge(idrop)%element%dwc
               element_stat(num_element(rank),rank)%id  = drops_merge(idrop)%element%id
            end if !drops_merge(idrop)%flag_center_mass
         end do ! idrop
      end if ! num_drop(rank)

      if ( num_part(rank) > 0 ) then  
         do idrop = 1,num_part(rank) 
            num_element(rank) = num_element(rank) + 1
            element_stat(num_element(rank),rank)%vol = parts(idrop,rank)%element%vol
            element_stat(num_element(rank),rank)%xc  = parts(idrop,rank)%element%xc
            element_stat(num_element(rank),rank)%yc  = parts(idrop,rank)%element%yc
            element_stat(num_element(rank),rank)%zc  = parts(idrop,rank)%element%zc
            element_stat(num_element(rank),rank)%uc  = parts(idrop,rank)%element%uc
            element_stat(num_element(rank),rank)%vc  = parts(idrop,rank)%element%vc
            element_stat(num_element(rank),rank)%wc  = parts(idrop,rank)%element%wc
            element_stat(num_element(rank),rank)%duc  = parts(idrop,rank)%element%duc
            element_stat(num_element(rank),rank)%dvc  = parts(idrop,rank)%element%dvc
            element_stat(num_element(rank),rank)%dwc  = parts(idrop,rank)%element%dwc
            element_stat(num_element(rank),rank)%id  = parts(idrop,rank)%element%id
         end do ! idrop
      end if ! num_drop(rank)

      ! Output element stats without collective operation
      if (DropStatisticsMethod == DropStatistics_WriteElementDataEachRank ) then
         filename=TRIM(out_path)//'/element-stats_'//TRIM(int2text(rank,padding))//'.dat'
         inquire(FILE=filename,EXIST=file_exist)
         if ( file_exist ) then 
            OPEN(UNIT=102,FILE=filename, STATUS='old', POSITION='append', ACTION='write')
         else 
            OPEN(UNIT=102,FILE=filename, STATUS='new', ACTION='write')
         end if ! file_exist
         if ( num_element(rank) > 0 ) then
            do ielement = 1, num_element(rank)
               write(102,'(11(E15.6,1X))')   time, & 
                                             element_stat(ielement,rank)%xc, &
                                             element_stat(ielement,rank)%yc, &
                                             element_stat(ielement,rank)%zc, &
                                             element_stat(ielement,rank)%uc, &
                                             element_stat(ielement,rank)%vc, &
                                             element_stat(ielement,rank)%wc, &
                                             element_stat(ielement,rank)%duc, &
                                             element_stat(ielement,rank)%dvc, &
                                             element_stat(ielement,rank)%dwc, &
                                             element_stat(ielement,rank)%vol
            end do ! ielement
         end if ! num_element(rank) 
         close(102)
        
         deallocate( element_stat )
         return
      end if ! DropStatisticsMethod 

      ! Collect all discrete elements to rank 0
      
      !  Setup MPI derived type for element_type 
      offsets (0) = 0 
      oldtypes(0) = MPI_REAL8 
      blockcounts(0) = 10 
      call MPI_TYPE_EXTENT(MPI_REAL8, r8extent, ierr) 
      offsets    (1) = blockcounts(0)*r8extent 
      oldtypes   (1) = MPI_INTEGER  
      blockcounts(1) = 1  

      call MPI_TYPE_STRUCT(2, blockcounts, offsets, oldtypes, & 
                           MPI_element_type, ierr) 
      call MPI_TYPE_COMMIT(MPI_element_type, ierr)

      num_element_rank= num_element(rank)
      call MPI_ALLGATHER(num_element_rank, 1, MPI_INTEGER, &
                         num_element(:),   1, MPI_INTEGER, MPI_Comm_World, ierr)

      if ( rank > 0 ) then
         call MPI_ISEND(element_stat(1,rank),num_element(rank), & 
                        MPI_element_type, 0, 14, MPI_COMM_WORLD, req(1), ierr)
         call MPI_WAIT(req(1),sta(:,1),ierr)
      else
         do irank = 1,nPdomain-1
            call MPI_IRECV(element_stat(1,irank),num_element(irank), & 
                           MPI_element_type, irank, 14, MPI_COMM_WORLD, req(2), ierr)
            call MPI_WAIT(req(2),sta(:,2),ierr)
         end do ! irank
      end if ! rank

      ! Statisitcs or other operation at rank 0 
      if ( rank == 0 ) then
         select case(DropStatisticsMethod) 
            ! Plot each dispersed element
            case(DropStatistics_PlotElement)
               ielem_plot = 0 
               do irank = 0,nPdomain-1
                  if ( num_element(irank) > 0 ) then 
                     do ielement = 1, num_element(irank) 
                        ielem_plot = ielem_plot + 1 
                        OPEN(UNIT=200+ielem_plot,FILE=TRIM(out_path)//'/element-'//TRIM(int2text(ielem_plot,padding))//'.dat')
                        write(200+ielem_plot,'(11(E15.6,1X))') time,element_stat(ielement,irank)%xc, & 
                                                      element_stat(ielement,irank)%yc, &
                                                      element_stat(ielement,irank)%zc, &
                                                      element_stat(ielement,irank)%uc, &
                                                      element_stat(ielement,irank)%vc, &
                                                      element_stat(ielement,irank)%wc, &
                                                      element_stat(ielement,irank)%duc, &
                                                      element_stat(ielement,irank)%dvc, &
                                                      element_stat(ielement,irank)%dwc, &
                                                      element_stat(ielement,irank)%vol
                     end do ! ielement
                  end if ! num_element_irank) 
               end do ! irank
            ! compute element size pdf
            case(DropStatistics_ElementSizePDF)
               gap = yLength/dble(Ny) 
               dmax = dble(num_gaps)*gap
               count_element(:) = 0
               do irank = 0,nPdomain-1
                  if ( num_element(irank) > 0 ) then 
                     do ielement = 1, num_element(irank)
                        ! sum particles 
                        d =(element_stat(ielement,irank)%vol*6.d0/PI)**OneThird
                        if ( d < dmax ) then 
                           igap = int(d/gap)+1
                           count_element(igap) = count_element(igap) + 1
                        end if ! d
                     end do ! ielement
                  end if ! num_element(irank)
               end do ! irank
               open(unit=101,file=TRIM(out_path)//'/element-size-pdf_'//TRIM(int2text(tswap,padding))//'.dat')

               do igap = 1,num_gaps
                  write(101,*) (dble(igap)-0.5d0)*gap,count_element(igap)
               end do ! igap
               close(101)
            case(DropStatistics_WriteElementData)
               total_num_element = sum(num_element)
               open(unit=102,file=TRIM(out_path)//'/element-stats_'//TRIM(int2text(tswap,padding))//'.dat')
               write(102,'(E15.6,1X,I7)')   time, total_num_element 
               do irank = 0,nPdomain-1
                  if ( num_element(irank) > 0 ) then 
                     do ielement = 1, num_element(irank)
                        write(102,'(11(E15.6,1X))')   element_stat(ielement,irank)%xc, & 
                                                      element_stat(ielement,irank)%yc, &
                                                      element_stat(ielement,irank)%zc, &
                                                      element_stat(ielement,irank)%uc, &
                                                      element_stat(ielement,irank)%vc, &
                                                      element_stat(ielement,irank)%wc, &
                                                      element_stat(ielement,irank)%duc, &
                                                      element_stat(ielement,irank)%dvc, &
                                                      element_stat(ielement,irank)%dwc, &
                                                      element_stat(ielement,irank)%vol
                     end do ! ielement
                  end if ! num_element_irank) 
               end do ! irank
               close(102)
            case default
               call lpperror("unknown drop statistics method!")
         end select 
      end if ! rank

      ! finalize
      call MPI_TYPE_FREE(MPI_element_type, ierr)
      deallocate( element_stat )

   end subroutine drop_statistics

   subroutine CollectDropMerge
      implicit none

      include 'mpif.h'
      integer :: req(2),sta(MPI_STATUS_SIZE,2),MPI_Comm,ireq,ierr
      integer :: MPI_drop_merge_comm_type, oldtypes(0:1), blockcounts(0:1), & 
                 offsets(0:1), extent,r8extent 
      integer :: max_num_drop_merge_use 
      integer :: irank, idrop

      max_num_drop_merge_use = maxval(num_drop_merge)
      allocate( drops_merge_comm(max_num_drop_merge_use,0:nPdomain-1) )

      !  Setup MPI derived type for drop_merge_comm
      offsets (0) = 0 
      oldtypes(0) = MPI_REAL8 
      blockcounts(0) = 10 
      call MPI_TYPE_EXTENT(MPI_REAL8, r8extent, ierr) 
      offsets    (1) = blockcounts(0)*r8extent 
      oldtypes   (1) = MPI_INTEGER  
      blockcounts(1) = 1+1+maxnum_diff_tag+1  

      call MPI_TYPE_STRUCT(2, blockcounts, offsets, oldtypes, & 
                                  MPI_drop_merge_comm_type, ierr) 
      call MPI_TYPE_COMMIT(MPI_drop_merge_comm_type, ierr)

      !  initialize drops_merge_comm 
      if ( num_drop_merge(rank) > 0 ) then 
         do idrop = 1, num_drop_merge(rank)
            drops_merge_comm(idrop,rank)%id            = drops_merge(idrop)%element%id
            drops_merge_comm(idrop,rank)%num_diff_tag  = drops_merge(idrop)%num_diff_tag
            drops_merge_comm(idrop,rank)%diff_tag_list = drops_merge(idrop)%diff_tag_list
            drops_merge_comm(idrop,rank)%vol           = drops_merge(idrop)%element%vol
            drops_merge_comm(idrop,rank)%xc            = drops_merge(idrop)%element%xc
            drops_merge_comm(idrop,rank)%yc            = drops_merge(idrop)%element%yc
            drops_merge_comm(idrop,rank)%zc            = drops_merge(idrop)%element%zc
            drops_merge_comm(idrop,rank)%uc            = drops_merge(idrop)%element%uc
            drops_merge_comm(idrop,rank)%vc            = drops_merge(idrop)%element%vc
            drops_merge_comm(idrop,rank)%wc            = drops_merge(idrop)%element%wc
            drops_merge_comm(idrop,rank)%duc            = drops_merge(idrop)%element%duc
            drops_merge_comm(idrop,rank)%dvc            = drops_merge(idrop)%element%dvc
            drops_merge_comm(idrop,rank)%dwc            = drops_merge(idrop)%element%dwc
         end do ! idrop
      end if ! num_drop_merge(rank)

      ! Send all droplet pieces to rank 0
      if ( rank > 0 ) then
         if ( num_drop_merge(rank) > 0 ) then 
            call MPI_ISEND(drops_merge_comm(1,rank),num_drop_merge(rank), & 
                           MPI_drop_merge_comm_type, 0,    13, MPI_COMM_WORLD, req(1), ierr)
            call MPI_WAIT(req(1),sta(:,1),ierr)
         end if ! num_drop_merge(rank)
      else
         do irank = 1,nPdomain-1
            if ( num_drop_merge(irank) > 0 ) then 
               call MPI_IRECV(drops_merge_comm(1,irank),num_drop_merge(irank), & 
                              MPI_drop_merge_comm_type, irank, 13, MPI_COMM_WORLD, req(2), ierr)
               call MPI_WAIT(req(2),sta(:,2),ierr)
            end if ! num_drop_merge(irank)
         end do ! irank
      end if ! rank

      ! finalize
      call MPI_TYPE_FREE(MPI_drop_merge_comm_type, ierr)

   end subroutine CollectDropMerge
   
   subroutine DistributeDropMerge
      implicit none

      include 'mpif.h'
      integer :: irank,idrop
      integer :: req(2),sta(MPI_STATUS_SIZE,2),MPI_Comm,ireq,ierr
      integer :: MPI_drop_merge_comm_type, oldtypes(0:1), blockcounts(0:1), & 
                 offsets(0:1), extent,r8extent 

      !  Setup MPI derived type for drop_merge_comm
      offsets (0) = 0 
      oldtypes(0) = MPI_REAL8 
      blockcounts(0) = 10 
      call MPI_TYPE_EXTENT(MPI_REAL8, r8extent, ierr) 
      offsets    (1) = blockcounts(0)*r8extent 
      oldtypes   (1) = MPI_INTEGER  
      blockcounts(1) = 1+1+maxnum_diff_tag+1  

      call MPI_TYPE_STRUCT(2, blockcounts, offsets, oldtypes, & 
                                  MPI_drop_merge_comm_type, ierr) 
      call MPI_TYPE_COMMIT(MPI_drop_merge_comm_type, ierr)

      ! Send merged droplet from rank 0 to all ranks 
      if ( rank > 0 ) then
         if ( num_drop_merge(rank) > 0 ) then 
            call MPI_IRECV(drops_merge_comm(1,rank),num_drop_merge(rank), & 
                           MPI_drop_merge_comm_type, 0,    13, MPI_COMM_WORLD, req(1), ierr)
            call MPI_WAIT(req(1),sta(:,1),ierr)
         end if ! num_drop_merge(rank)
      else
         do irank = 1,nPdomain-1
               if ( num_drop_merge(irank) > 0 ) then 
               call MPI_ISEND(drops_merge_comm(1,irank),num_drop_merge(irank), & 
                              MPI_drop_merge_comm_type, irank, 13, MPI_COMM_WORLD, req(2), ierr)
               call MPI_WAIT(req(2),sta(:,2),ierr)
            end if ! num_drop_merge(irank)
         end do ! irank
      end if ! rank

      if ( num_drop_merge(rank) > 0 ) then
         do idrop = 1, num_drop_merge(rank)
            drops_merge(idrop)%element%vol = drops_merge_comm(idrop,rank)%vol
            drops_merge(idrop)%element%xc  = drops_merge_comm(idrop,rank)%xc
            drops_merge(idrop)%element%yc  = drops_merge_comm(idrop,rank)%yc
            drops_merge(idrop)%element%zc  = drops_merge_comm(idrop,rank)%zc
            drops_merge(idrop)%element%uc  = drops_merge_comm(idrop,rank)%uc
            drops_merge(idrop)%element%vc  = drops_merge_comm(idrop,rank)%vc
            drops_merge(idrop)%element%wc  = drops_merge_comm(idrop,rank)%wc
            drops_merge(idrop)%element%duc  = drops_merge_comm(idrop,rank)%duc
            drops_merge(idrop)%element%dvc  = drops_merge_comm(idrop,rank)%dvc
            drops_merge(idrop)%element%dwc  = drops_merge_comm(idrop,rank)%dwc
            drops_merge(idrop)%flag_center_mass  = drops_merge_comm(idrop,rank)%flag_center_mass
         end do ! idrop
      end if ! num_drop_merge(rank)

      ! finalize
      deallocate(drops_merge_comm)
      call MPI_TYPE_FREE(MPI_drop_merge_comm_type, ierr) 

   end subroutine DistributeDropMerge

   subroutine CreateTag2DropTable
      implicit none

      include 'mpif.h'
      integer :: req(6),sta(MPI_STATUS_SIZE,6),MPI_Comm,ireq,ierr
      integer :: idrop, irank

      allocate( tag_dropid   (1:total_num_tag) )
      allocate( tag_rank     (1:total_num_tag) )
      allocate( tag_mergeflag(1:total_num_tag) )

      if ( num_drop(rank) > 0 ) then 
      do idrop = 1,num_drop(rank)
         drops(idrop)%element%id = tag_id( drops_cell_list(1,idrop), &
                                                drops_cell_list(2,idrop), &
                                                drops_cell_list(3,idrop) )
         tag_dropid   (drops(idrop)%element%id) = idrop
         tag_rank     (drops(idrop)%element%id) = rank 
         tag_mergeflag(drops(idrop)%element%id) = 0 
      end do ! idrop
      end if ! num_drop(rank)
      
      if ( num_drop_merge(rank) > 0 ) then 
      do idrop = 1,num_drop_merge(rank)
         drops_merge(idrop)%element%id = tag_id( drops_merge_cell_list(1,idrop), &
                                                      drops_merge_cell_list(2,idrop), &
                                                      drops_merge_cell_list(3,idrop) )
         tag_dropid   (drops_merge(idrop)%element%id) = idrop
         tag_rank     (drops_merge(idrop)%element%id) = rank 
         tag_mergeflag(drops_merge(idrop)%element%id) = 1 
      end do ! idrop
      end if ! num_drop_merge(rank) 

      if ( nPdomain > 1 ) then 
         ! MPI communication for tag tables
         if ( rank > 0 ) then
            if ( num_tag(rank) > 0 ) then 
               call MPI_ISEND(tag_dropid   (tagmin(rank)),num_tag(rank), & 
                              MPI_INTEGER, 0,    10, MPI_COMM_WORLD, req(1), ierr)
               call MPI_ISEND(tag_rank     (tagmin(rank)),num_tag(rank), & 
                              MPI_INTEGER, 0,    11, MPI_COMM_WORLD, req(2), ierr)
               call MPI_ISEND(tag_mergeflag(tagmin(rank)),num_tag(rank), & 
                              MPI_INTEGER, 0,    12, MPI_COMM_WORLD, req(3), ierr)
               call MPI_WAITALL(3,req(1:3),sta(:,1:3),ierr)
            end if ! num_tag(rank)
         else
            ireq = 0 
            do irank = 1,nPdomain-1
               if ( num_tag(irank) > 0 ) then
                  call MPI_IRECV(tag_dropid   (tagmin(irank)),num_tag(irank), & 
                                 MPI_INTEGER, irank, 10, MPI_COMM_WORLD, req(4), ierr)
                  call MPI_IRECV(tag_rank     (tagmin(irank)),num_tag(irank), & 
                                 MPI_INTEGER, irank, 11, MPI_COMM_WORLD, req(5), ierr)
                  call MPI_IRECV(tag_mergeflag(tagmin(irank)),num_tag(irank), & 
                                 MPI_INTEGER, irank, 12, MPI_COMM_WORLD, req(6), ierr)
                  call MPI_WAITALL(3,req(4:6),sta(:,4:6),ierr)
               end if ! num_tag(irank)
            end do ! irank
         end if ! rank
      end if ! nPdomain 

   end subroutine CreateTag2DropTable

   subroutine ReleaseTag2DropTable
      implicit none

      deallocate(tag_dropid   )
      deallocate(tag_rank     )
      deallocate(tag_mergeflag)

   end subroutine ReleaseTag2DropTable

! ==============================================
! test tagging drops
! ==============================================
   subroutine test_tagDrop()
      implicit none
      integer :: i,j,k
      call output_VOF(0,is,ie+1,js,je+1,ks,ke+1)
      call tag_drop()
      if ( nPdomain > 1 ) call tag_drop_all
      call hello_coucou  ! 12
      call CreateTag2DropTable
      if ( nPdomain > 1 ) call merge_drop_pieces 
      call ReleaseTag2DropTable
      if ( DropStatisticsMethod > 0 ) call drop_statistics(0,0.d0)
   end subroutine test_tagDrop

  subroutine output_tag(nf,i1,i2,j1,j2,k1,k2)
    implicit none
    integer ::nf,i1,i2,j1,j2,k1,k2,i,j,k
    real(8) :: cfiltered
    character(len=30) :: rootname,filename
    rootname=trim(out_path)//'/VTK/tag'//TRIM(int2text(nf,padding))//'-'
    if(rank==0) call append_tag_visit_file(TRIM(rootname))

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
    write(8,17)'tag'
    write(8,18)
19  format('POINT_DATA ',I17)
17  format('SCALARS ',A20,' float 1')
18  format('LOOKUP_TABLE default')

    do k=k1,k2; do j=j1,j2; do i=i1,i2
      write(8,210) tag_id(i,j,k)
    enddo; enddo; enddo
210 format(i5)
    close(8)
  end subroutine output_tag
!-------------------------------------------------------------------------------------------------

  subroutine append_tag_visit_file(rootname)
    implicit none
    character(*) :: rootname
    integer prank
    if(rank.ne.0) call pariserror('rank.ne.0 in append_tag')
    if(tag_opened==0) then
       OPEN(UNIT=88,FILE='tag.visit')
       write(88,10) nPdomain
10     format('!NBLOCKS ',I4)
       tag_opened=1
    else
       OPEN(UNIT=88,FILE='tag.visit',position='append')
    endif
    do prank=0,NpDomain-1
       write(88,11) rootname//TRIM(int2text(prank,padding))//'.vtk'
11     format(A)
    enddo
    close(88)
  end subroutine  append_tag_visit_file

!   subroutine output_tag()
!      implicit none
!      integer :: i,j,k

!      OPEN(UNIT=90,FILE=TRIM(out_path)//'/tag-tecplot'//TRIM(int2text(rank,padding))//'.dat')

!      write(90,*) 'title= " 3d tag "'
!      write(90,*) 'variables = "x", "y", "z", "tag", "c" '
!      !write(90,*) 'zone i=,',Nx/nPx+Ng*2, 'j=',Ny/nPy+Ng*2, 'k=',Nz/nPz+Ng*2,'f=point'
!      write(90,*) 'zone i=',imax-imin+1, ',j=',jmax-jmin+1, ',k=',kmax-kmin+1,'f=point'
!      do k = kmin,kmax
!         do j=jmin,jmax
!            do i=imin,imax 
!               write(90,'(4(I5,1X),(E15.8))') i,j,k,tag_id(i,j,k),cvof(i,j,k)
!            end do ! i
!         end do ! j
!      end do ! k
!
!   end subroutine output_tag

! ==============================================
! test drop-particle conversion
! ==============================================
   subroutine test_DropPartConv(tswap)
      implicit none

      integer, intent(in) :: tswap
      integer :: i,j,k
      integer :: ib,ipart

      type(drop), dimension(NumBubble) :: drops_ex

      ! tag droplets and calculate drop properties
      call tag_drop()
      if ( nPdomain > 1 ) then 
         call tag_drop_all
         call merge_drop_pieces 
      end if ! nPdomain

      OPEN(UNIT=90,FILE=TRIM(out_path)//'/VOF_before_'//TRIM(int2text(rank,padding))//'.dat')
      write(90,*) 'title= " VOF drops before conversion "'
      write(90,*) 'variables = "x", "y", "z", "c", "tag" '
      write(90,*) 'zone i=,',Nx/nPx, 'j=',Ny/nPy, 'k=',Nz/nPz,'f=point'
      do k = ks,ke
         do j=js,je
            do i=is,ie 
               write(90,'(4(E15.8,1X),(I5))') x(i),y(j),z(k),cvof(i,j,k),tag_id(i,j,k)
            end do ! i
         end do ! j
      end do ! k
      CLOSE(90)

      if ( rank == 0 ) then 
      OPEN(UNIT=91,FILE=TRIM(out_path)//'/dropvol-'//TRIM(int2text(rank,padding))//'.dat')
      OPEN(UNIT=92,FILE=TRIM(out_path)//'/dropvol_ex-'//TRIM(int2text(rank,padding))//'.dat')
!      call QSort(drops,NumBubble)
      do ib = 1, NumBubble
         drops_ex(ib)%element%xc  = xc(ib)
         drops_ex(ib)%element%vol = 4.0d0*rad(ib)**3.d0*PI/3.d0
      end do ! i 
!      call QSort(drops_ex,NumBubble)
      do ib = 1, NumBubble 
         write(91,*) drops   (ib)%element%xc, drops   (ib)%element%vol
         write(92,*) drops_ex(ib)%element%xc, drops_ex(ib)%element%vol
      end do ! i 
      CLOSE(91)
      CLOSE(92)
      end if ! rank

      ! convert droplets to particles
      call ConvertVOF2LPP(tswap)

      ! output droplets & particles
      OPEN(UNIT=93,FILE=TRIM(out_path)//'/VOF_after_'//TRIM(int2text(rank,padding))//'.dat')
      write(93,*) 'title= " VOF drops after conversion "'
      write(93,*) 'variables = "x", "y", "z", "c", "tag" '
      write(93,*) 'zone i=,',Nx/nPx, 'j=',Ny/nPy, 'k=',Nz/nPz,'f=point'
      do k = ks,ke
         do j=js,je
            do i=is,ie
               write(93,'(4(E15.8,1X),(I5))') x(i),y(j),z(k),cvof(i,j,k),tag_id(i,j,k)
            end do ! i
         end do ! j
      end do ! k
      CLOSE(93)

!      call output_VOF(0,imin,imax,jmin,jmax,kmin,kmax)

      if ( num_part(rank) > 0 ) then 
      OPEN(UNIT=94,FILE=TRIM(out_path)//'/LPP_after_'//TRIM(int2text(rank,padding))//'.dat')
!      write(94,*) 'title= " Lagrangian particles after conversion "'
!      write(94,*) 'variables = "x", "y", "z"'
      do ipart = 1,num_part(rank)
         write(94,*)  & 
            parts(ipart,rank)%element%xc,parts(ipart,rank)%element%yc,&
            parts(ipart,rank)%element%zc,parts(ipart,rank)%element%vol
      end do ! ipart
      CLOSE(94)
      end if ! num_part(rank)

      if ( rank == 0 ) then 
         OPEN(UNIT=95,FILE=TRIM(out_path)//'/LPP_ex.dat')
         do ib = 1,NumBubble
            if ( 4.d0*PI*rad(ib)**3.d0/3.d0 < vol_cut ) &  
               write(95,*) xc(ib), yc(ib), zc(ib), 4.d0*PI*rad(ib)**3.d0/3.d0
         end do ! idrop
         CLOSE(95)
      end if ! rank
   end subroutine test_DropPartConv

! ===============================================
! Testing section
! ===============================================
   subroutine test_Lag_part(tswap)
      implicit none
      include 'mpif.h'
      integer, intent(in) :: tswap
!      integer :: ierr
                     
      if(test_tag) then
         call test_tagDrop()
      else if ( test_D2P ) then
         call test_DropPartConv(tswap)
      end if
                                                                             
   end subroutine test_Lag_part

! ===============================================
! Sort droplets according volume 
! with quicksort method
! 
! Modified based on the source code given in 
! http://rosettacode.org
! ===============================================
 
   recursive subroutine QSort(a,na)

   ! DUMMY ARGUMENTS
   integer, intent(in) :: nA
   type (drop), dimension(nA), intent(in out) :: A

   ! LOCAL VARIABLES
   integer :: left, right
   real(8) :: random
   real(8) :: pivot
   type (drop) :: temp
   integer :: marker

   if (nA > 1) then

      ! random pivor (not best performance, but avoids   worst-case)
      call random_number(random)
      pivot = A(int(random*real(nA-1))+1)%element%vol

      left = 0
      right = nA + 1

      do while (left < right)
         right = right - 1
         do while (A(right)%element%vol > pivot)
            right = right - 1
         end do
         left = left + 1
         do while (A(left)%element%vol < pivot)
            left = left + 1
         end do
         if (left < right) then
            temp = A(left)
            A(left) = A(right)
            A(right) = temp
         end if
      end do

      if (left == right) then
         marker = left + 1
      else
         marker = left
      end if

      call QSort(A(:marker-1),marker-1)
      call QSort(A(marker:),nA-marker+1)

    end if ! nA
   end subroutine QSort

   subroutine ConvertVOF2LPP(tswap)
      implicit none

      include 'mpif.h'

      integer, intent(in) :: tswap

      integer :: idrop,ilist
      logical :: ConvertDropFlag,convertDoneFlag
      real(8) :: MinDistPart2CellCenter, DistPart2CellCenter
      integer :: ierr
      real(8) :: uf,vf,wf,dp,ConvertRegSize,sduf,sdvf,sdwf
      integer :: i,j,k
      integer :: i1,ic,i2,j1,jc,j2,k1,kc,k2
      real(8) :: x1,x2,y1,y2,z1,z2
      real(8) :: ufp,vfp,wfp,dist2,wt
      logical :: ConvertMergeDrop=.false.
      real(8) :: volcell
      integer :: num_part_rank

      if ( num_drop(rank) > 0 ) then 
      do idrop = 1,num_drop(rank)

         ! Find cell index and check conversion criteria 
         call FindPartLocCell(drops(idrop)%element%xc,  & 
                              drops(idrop)%element%yc,  & 
                              drops(idrop)%element%zc,  &
                              ic,jc,kc) 
         call CheckConvertDropCriteria(drops(idrop)%element%vol, & 
                                       drops(idrop)%element%xc,  & 
                                       drops(idrop)%element%yc,  & 
                                       drops(idrop)%element%zc,  &
                                       ic,jc,kc,                      & 
                                       ConvertDropFlag,               & 
                                       CriteriaConvertCase,           &
                                       drops(idrop)%AspRatio)

         volcell = dx(ic)*dy(jc)*dz(kc)
         if ( ConvertDropFlag ) then
            write(*,*) 'Drop is converted to particle', idrop,rank,tswap, &
               drops(idrop)%element%vol/volcell
            
            ! transfer droplet properties to particle
            num_part(rank) = num_part(rank) + 1
            parts(num_part(rank),rank)%element = drops(idrop)%element

            parts(num_part(rank),rank)%ic = ic 
            parts(num_part(rank),rank)%jc = jc 
            parts(num_part(rank),rank)%kc = kc

            ! Record the time step the drop is converted from VOF to LPP
            parts(num_part(rank),rank)%tstepConvert = tswap

            ! Define conversion region
            dp = (6.d0*drops(idrop)%element%vol/PI)**OneThird
            ConvertRegSize = ConvertRegSizeToDiam*dp
            call FindCellIndexBdryConvertRegUnifMesh(ConvertRegSize, & 
                                             i1,ic,i2,j1,jc,j2,k1,kc,k2)
            ! reconstruct undisturbed flow field within the conversion region
            do k=k1,k2-1; do j=j1,j2-1; do i=i1,i2-1
               ! compute undisturbed flow field from flow property outside
               x1 = xh(i)-xh(i1-1)
               x2 = xh(i2)-xh(i)
               y1 = yh(j)-yh(j1-1)
               y2 = yh(j2)-yh(j)
               z1 = zh(k)-zh(k1-1)
               z2 = zh(k2)-zh(k)
               call Surface2VolumeIntrpl(u(i1-1,j,k),u(i2,j,k), &
                                         u(i,j1-1,k),u(i,j2,k), &
                                         u(i,j,k1-1),u(i,j,k2), &
                                         x1,x2,y1,y2,z1,z2,uf)
               call Surface2VolumeIntrpl(v(i1-1,j,k),v(i2,j,k), & 
                                         v(i,j1-1,k),v(i,j2,k), & 
                                         v(i,j,k1-1),v(i,j,k2), &
                                         x1,x2,y1,y2,z1,z2,vf)
               call Surface2VolumeIntrpl(w(i1-1,j,k),w(i2,j,k), & 
                                         w(i,j1-1,k),w(i,j2,k), & 
                                         w(i,j,k1-1),w(i,j,k2), &
                                         x1,x2,y1,y2,z1,z2,wf)
                                         
               ! compute undisturbed flow field from droplet acceleration
!               call ComputeUndisturbedVelEE(drops(idrop)%element%uc,  & 
!                                            drops(idrop)%element%duc, & 
!                                            drops(idrop)%element%vc,  & 
!                                            drops(idrop)%element%dvc, & 
!                                            drops(idrop)%element%wc,  & 
!                                            drops(idrop)%element%dwc, & 
!                                            dp,rho1,rho2,mu1,mu2,Gx,Gy,Gz,ufp,vfp,wfp)

!               dist2 = (x(i)-drops(idrop)%element%xc)**2.d0 & 
!                     + (y(j)-drops(idrop)%element%yc)**2.d0 & 
!                     + (z(k)-drops(idrop)%element%zc)**2.d0
!               wt = exp(-dist2*4.d0/dp/dp)
!               u(i,j,k) = uf*(1.d0-wt) + ufp*wt
!               v(i,j,k) = vf*(1.d0-wt) + vfp*wt
!               w(i,j,k) = wf*(1.d0-wt) + wfp*wt
               
               u(i,j,k) = uf 
               v(i,j,k) = vf
               w(i,j,k) = wf
            end do; end do; end do

            ! remove droplet vof structure
            cvof(i1:i2,j1:j2,k1:k2) = 0.d0 

         end if !ConvertDropFlag
      end do ! idrop
      end if ! num_drop(rank) 

      ! XXX Note: drop_merge_converge is not working yet
      if ( num_drop_merge(rank) > 0 .and. ConvertMergeDrop ) then
      do idrop = 1,num_drop_merge(rank)

         call CheckConvertDropCriteria(drops_merge(idrop)%element%vol, & 
                                       drops_merge(idrop)%element%xc,  & 
                                       drops_merge(idrop)%element%yc,  & 
                                       drops_merge(idrop)%element%zc,  &
                                       ic,jc,kc,                            & 
                                       ConvertDropFlag,CriteriaConvertCase, & 
                                       AspRatioSphere)

         if ( ConvertDropFlag ) then
            write(*,*) 'Drop_merge is converted to particle', idrop,rank,tswap

            ! compute average fluid quantities
            ! Note: XXX temporary, need to be improved later 
            uf = 0.d0 
            vf = 0.d0 
            wf = 0.d0 

            ! remove droplet vof structure
!            do ilist = 1,drops_merge(idrop)%num_cell_drop
!               cvof(drops_merge_cell_list(1,ilist,idrop), &
!                    drops_merge_cell_list(2,ilist,idrop), &
!                    drops_merge_cell_list(3,ilist,idrop)) = 0.0
!                  u(drops_merge_cell_list(1,ilist,idrop), &
!                    drops_merge_cell_list(2,ilist,idrop), &
!                    drops_merge_cell_list(3,ilist,idrop)) = uf
!                  v(drops_merge_cell_list(1,ilist,idrop), &
!                    drops_merge_cell_list(2,ilist,idrop), &
!                    drops_merge_cell_list(3,ilist,idrop)) = vf
!                  w(drops_merge_cell_list(1,ilist,idrop), &
!                    drops_merge_cell_list(2,ilist,idrop), &
!                    drops_merge_cell_list(3,ilist,idrop)) = wf
!            end do ! ilist

            ! remove droplet vof structure
!            do ilist = 1,drops_merge(idrop)%num_gcell
!               cvof(drops_merge_gcell_list(1,ilist,idrop), &
!                    drops_merge_gcell_list(2,ilist,idrop), &
!                    drops_merge_gcell_list(3,ilist,idrop)) = 0.0
!                  u(drops_merge_gcell_list(1,ilist,idrop), &
!                    drops_merge_gcell_list(2,ilist,idrop), &
!                    drops_merge_gcell_list(3,ilist,idrop)) = uf
!                  v(drops_merge_gcell_list(1,ilist,idrop), &
!                    drops_merge_gcell_list(2,ilist,idrop), &
!                    drops_merge_gcell_list(3,ilist,idrop)) = vf
!                  w(drops_merge_gcell_list(1,ilist,idrop), &
!                    drops_merge_gcell_list(2,ilist,idrop), &
!                    drops_merge_gcell_list(3,ilist,idrop)) = wf
!            end do ! ilist

            ! transfer droplet properties to particle if center of mass located
            ! in this droplet piece
            if ( drops_merge(idrop)%flag_center_mass == 1 ) then 
               num_part(rank) = num_part(rank) + 1
               parts(num_part(rank),rank)%element = drops_merge(idrop)%element
            
               ! Find particle location cell
!               MinDistPart2CellCenter = 1.0d10
!               do ilist = 1,drops_merge(idrop)%num_cell_drop
!                  DistPart2CellCenter = ( drops_merge(idrop)%element%xc  & 
!                              - x(drops_merge_cell_list(1,ilist,idrop)))**2.d0 & 
!                                      + ( drops_merge(idrop)%element%yc  & 
!                              - y(drops_merge_cell_list(2,ilist,idrop)))**2.d0 &
!                                      + ( drops_merge(idrop)%element%zc  & 
!                              - z(drops_merge_cell_list(3,ilist,idrop)))**2.d0
!                  if ( DistPart2CellCenter < MinDistPart2CellCenter ) then 
!                     MinDistPart2CellCenter = DistPart2CellCenter
!                     parts(num_part(rank),rank)%ic = drops_merge_cell_list(1,ilist,idrop)
!                     parts(num_part(rank),rank)%jc = drops_merge_cell_list(2,ilist,idrop)
!                     parts(num_part(rank),rank)%kc = drops_merge_cell_list(3,ilist,idrop)
!                  end if !DistPart2CellCenter
!               end do ! ilist
            end if ! flag_center_mass 

         end if !ConvertDropFlag
      end do ! idrop
      end if ! num_drop_merge(rank) 

      ! Update num_part to all ranks. Note: no need for num_drop &
      ! num_drop_merge since they will be regenerated next step
      num_part_rank = num_part(rank)
      call MPI_ALLGATHER(num_part_rank, 1, MPI_INTEGER, &
                         num_part(:)  , 1, MPI_INTEGER, MPI_Comm_World, ierr)

   end subroutine ConvertVOF2LPP

   subroutine CheckConvertDropCriteria(vol,xc,yc,zc,ic,jc,kc,ConvertDropFlag,CriteriaConvertCase,AspRatio)
      implicit none

      real(8), intent(in ) :: vol,xc,yc,zc,AspRatio
      integer, intent(in ) :: ic,jc,kc
      integer, intent(in ) :: CriteriaConvertCase
      logical, intent(out) :: ConvertDropFlag

      real(8) :: volcell,dp,ConvertRegSize
      integer :: i1,i,i2,j1,j,j2,k1,k,k2
      integer :: tag0
      logical :: ConvertDropRegion 
      logical :: ConvertDropShape 

      ConvertDropFlag = .false.
      ConvertDropShape = .false.
      ConvertDropRegion = .false.

      ! Shape criteria
      volcell = dx(ic)*dy(jc)*dz(kc) 
      if ( vol < volcell*3.375d0 ) then   
         ConvertDropShape = .true. 
      ! Check aspect ratio for drop diameter larger 1.5dx
      else if ( (vol < vol_cut) .and. (AspRatio < AspRatioTol) ) then 
         ConvertDropShape = .true.
      end if ! vol, AspRatio

      ! Location criteria, only checked if shape criteria is satisfied
      if ( ConvertDropShape ) then 
      select case (CriteriaConvertCase) 
         case (CriteriaRectangle)
            if ( (xc  > xlpp_min) .and. (xc  < xlpp_max) .and. & 
                 (yc  > ylpp_min) .and. (yc  < ylpp_max) .and. & 
                 (zc  > zlpp_min) .and. (zc  < zlpp_max) ) then 
               ConvertDropRegion = .true.
            end if !xlpp_min...
         case (CriteriaCylinder )   ! Note: assuming axis along x-direction
                                    ! radius indicated by ylpp_min & ylpp_max
            if ( ( (yc-0.5d0)**2.d0 + (zc-0.5d0)**2.d0 > ylpp_min**2.d0)  .and. & 
                 ( (yc-0.5d0)**2.d0 + (zc-0.5d0)**2.d0 < ylpp_max**2.d0) ) then 
               ConvertDropRegion = .true.
            end if !xlpp_min... 
         case (CriteriaSphere   )   ! radius indicated by ylpp_min & ylpp_max 
            if ( ( (xc**2.d0 + yc**2.d0 + zc**2.d0) > ylpp_min**2.d0)   .and. & 
                 ( (xc**2.d0 + yc**2.d0 + zc**2.d0) < ylpp_max**2.d0) ) then 
               ConvertDropRegion = .true.
            end if !xlpp_min... 
         case (CriteriaJet      )   ! Note: assuming axis along x-direction
                                    ! xlpp_max varies in time & given as 
                                    ! xlpp_max = xlpp_min + zlpp_max*t
                                    ! radius indicated by ylpp_min & ylpp_max
            xlpp_max = xlpp_min + zlpp_max*time
            if ( (  xc  > xlpp_min) .and. (xc  < xlpp_max)                .and. & 
                 ( (yc-0.5d0)**2.d0 + (zc-0.5d0)**2.d0 > ylpp_min**2.d0)  .and. & 
                 ( (yc-0.5d0)**2.d0 + (zc-0.5d0)**2.d0 < ylpp_max**2.d0) ) then 
               ConvertDropRegion = .true.
            end if !xlpp_min... 
         case (CriteriaInterface)
            if ( RegAwayInterface(ic,jc,kc) ) then 
               ConvertDropRegion = .true.
            end if !RegAwayInterface...
         case (CriteriaIsoDrop)
            dp = (6.d0*vol/PI)**OneThird
            ConvertRegSize = ConvertRegSizeToDiam*dp
            call FindCellIndexBdryConvertRegUnifMesh(ConvertRegSize,&
                                                     i1,ic,i2,j1,jc,j2,k1,kc,k2)
            tag0 = tag_id(ic,jc,kc)
            ConvertDropRegion = .true.
            do i=i1,i2; do j=j1,j2; do k=k1,k2
               if ( tag_id(i,j,k) /= 0 .and. tag_id(i,j,k) /= tag0 ) then 
                  ConvertDropRegion = .false.
                  exit 
               end if ! tag_id
            end do; end do; end do
      end select
      end if ! ConvertDropShape 

      ConvertDropFlag = ConvertDropRegion .and. ConvertDropShape

      ! Always keep drops with diameter smaller than dx as LPP (mostly generated by VOF
      ! error)
      ! if ( vol < volcell*PI/6.d0 ) ConvertDropFlag = .true. 
      
   end subroutine CheckConvertDropCriteria

   subroutine MarkRegAwayInterface

      implicit none
      include 'mpif.h'
      integer :: i,j,k
      integer :: idrop,droptag,droprank
      real(8) :: vol_drop,d_cut
      integer :: shift,i1,i2,j1,j2,k1,k2   

      RegAwayInterface = .true.
      d_cut = (6.d0*vol_cut/PI)**OneThird
      shift = INT(dble(ConvertRegSizeToDiam)*0.5d0*d_cut/(xh(is+1)-xh(is)))

      do k=ks,ke; do j=js,je; do i=is,ie
         if ( cvof(i,j,k) > 0.d0 .and. cvof(i,j,k) < 1.d0 ) then 
            ! Check volume of drop corresponding to the current cell
            droptag  = tag_id(i,j,k) 
            idrop    = tag_dropid   (droptag) 
!            droprank = tag_rank     (droptag)  
            if ( tag_mergeflag(droptag) == 0 ) then 
               vol_drop = drops(idrop)%element%vol
            else if ( tag_mergeflag(droptag) == 1 ) then 
               vol_drop = drops_merge(idrop)%element%vol
            else 
               call lpperror("Unknown tag_mergeflag!")
            end if ! tag_mergeflag

            ! if interface belongs to big liquid drop, mark the neighbor cells false
            if ( vol_drop > vol_cut ) then
               i1 = MAX(i-shift,imin)
               i2 = MIN(i+shift,imax)
               j1 = MAX(j-shift,jmin)
               j2 = MIN(j+shift,jmax)
               k1 = MAX(k-shift,kmin)
               k2 = MIN(k+shift,kmax)
               RegAwayInterface(i1:i2,j1:j2,k1:k2) = .false.
            end if ! vol_drop
         end if ! cvof(i,j,k)
      end do; end do; end do
   end subroutine MarkRegAwayInterface

   subroutine ComputeUndisturbedVelEE(uc,duc,vc,dvc,wc,dwc,dp,rhof,rhop,muf,mup,Gx,Gy,Gz,uf,vf,wf)

      implicit none
      real(8), intent (in) :: uc,duc,vc,dvc,wc,dwc,dp,rhof,rhop,muf,mup,Gx,Gy,Gz
      real(8), intent(out) :: uf,vf,wf

      real(8) :: u0,v0,w0,u1,v1,w1,Re1,u2,v2,w2,Re2,fu1,fv1,fw1,fu2,fv2,fw2
      real(8) :: taup, a,bx,by,bz, nuf,tol,velc

      integer :: iter
      integer, parameter :: iter_max = 10
      real(8), parameter :: toler = 1.e-9

      taup = rhop*dp*dp/18.0d0/muf & 
            *(3.d0 + 3.d0*muf/mup)/(3.d0 + 2.d0*muf/mup)
      a  = 1.d0 - rhof/rhop
      bx = duc - Gx
      by = dvc - Gy
      bz = dwc - Gz
      nuf = muf/rhof

      u2 = uc; v2 = vc; w2 = wc
      Re2 = Rep(u2,uc,v2,vc,w2,wc,nuf,dp)
      fu2 = fsol(u2,uc,taup,Re2,a,bx)
      fv2 = fsol(v2,vc,taup,Re2,a,by)
      fw2 = fsol(w2,wc,taup,Re2,a,bz)

      u1 = fu2 + u2; v1 = fv2 + v2; w1 = fw2 + w2
      Re1 = Rep(u1,uc,v1,vc,w1,wc,nuf,dp)
      fu1 = fsol(u1,uc,taup,Re1,a,bx)
      fv1 = fsol(v1,vc,taup,Re1,a,by)
      fw1 = fsol(w1,wc,taup,Re1,a,bz)

      do iter = 1, iter_max
         u0 = (u2*fu1-u1*fu2)/(fu1-fu2)  
         v0 = (v2*fv1-v1*fv2)/(fv1-fv2)  
         w0 = (w2*fw1-w1*fw2)/(fw1-fw2)  
         velc = sqrt(u1*u1+v1*v1+w1*w1) + 1.d-16
         tol = sqrt( (u0-u1)**2.d0 + (v0-v1)**2.d0 + (w0-w1)**2.d0 )/velc 
         if ( tol < toler ) then  
            exit 
         else
            u2=u1; v2=v1; w2=w1
            fu2=fu1;fv2=fv1;fw2=fw1
            u1=u0; v1=v0; w1=w0
            Re1 = Rep(u1,uc,v1,vc,w1,wc,nuf,dp)
            fu1 = fsol(u1,uc,taup,Re1,a,bx)
            fv1 = fsol(v1,vc,taup,Re1,a,by)
            fw1 = fsol(w1,wc,taup,Re1,a,bz)
         end if ! tol
      end do ! iter
      uf = u0; vf = v0; wf = w0

   end subroutine ComputeUndisturbedVelEE

   function fsol(uf,up,taup,Rep,a,b)
      implicit none
      real(8), intent(in) :: uf,up,taup,Rep,a,b
      real(8) :: fsol

      fsol = up + taup/phi(dragmodel,Rep)*a*b - uf
   end function fsol

   function Rep(uf,up,vf,vp,wf,wp,nuf,dp)
      implicit none
      real(8), intent(in) :: uf,up,vf,vp,wf,wp,nuf,dp 
      real(8) :: Rep

      Rep = sqrt((uf-up)**2.d0 + (vf-vp)**2.d0 + (wf-wp)**2.d0)*dp/nuf

   end function Rep

   subroutine Surface2VolumeIntrpl(ux1,ux2,uy1,uy2,uz1,uz2,x1,x2,y1,y2,z1,z2,uf)
      
      implicit none
      real(8), intent(in)  :: ux1,ux2,uy1,uy2,uz1,uz2,x1,x2,y1,y2,z1,z2
      real(8), intent(out) :: uf
   
      real(8) :: wtsum,x1p,x2p,y1p,y2p,z1p,z2p

      x1p = 1.d0/x1 
      x2p = 1.d0/x2  
      y1p = 1.d0/y1
      y2p = 1.d0/y2
      z1p = 1.d0/z1
      z2p = 1.d0/z2
      wtsum = x1p + x2p + y1p + y2p + z1p + z2p 
       
      uf =(x1p*ux1 + x2p*ux2 & 
         + y1p*uy1 + y2p*uy2 & 
         + z1p*uz1 + z2p*uz2)/wtsum

   end subroutine Surface2VolumeIntrpl

   subroutine FindPartLocCell(xp,yp,zp,ip,jp,kp)
      implicit none

      real(8), intent (in) :: xp,yp,zp
      integer, intent(out) :: ip,jp,kp

      if ( xp < xh(imin)-dx(imin) .or. xp > xh(imax) .or. & 
           yp < yh(jmin)-dy(jmin) .or. yp > yh(jmax) .or. & 
           zp < zh(kmin)-dz(kmin) .or. zp > zh(kmax) )    & 
           call lpperror("Fail to find cell index for particle outside of domain!") 

      do ip = imin,imax
         if (xp <= xh(ip)) exit
      end do ! ip

      do jp = jmin,jmax
         if (yp <= yh(jp)) exit  
      end do ! jp

      do kp = kmin,kmax
         if (zp <= zh(kp)) exit  
      end do ! kp

   end subroutine 

   subroutine FindCellIndexBdryConvertReg(xc,yc,zc,l,i1,ic,i2,j1,jc,j2,k1,kc,k2)
      implicit none

      real(8), intent (in) :: xc,yc,zc,l
      integer, intent (in) :: ic,jc,kc
      integer, intent(out) :: i1,i2,j1,j2,k1,k2

      real(8) :: l2,xl,xr,yl,yr,zl,zr

      l2 = 0.5d0*l

      xl = xc - l2
      do i1 = ic,is,-1
         if ( xl > xh(i1-1) ) exit 
      end do ! i
      
      xr = xc + l2
      do i2 = ic,ie
         if ( xr < xh(i2  ) ) exit 
      end do ! i

      yl = yc - l2
      do j1 = jc,js,-1
         if ( yl > yh(i1-1) ) exit 
      end do ! j

      yr = yc + l2
      do j2 = jc,je
         if ( yr < yh(j2  ) ) exit 
      end do ! j

      zl = zc - l2
      do k1 = kc,ks,-1
         if ( zl > zh(k1-1) ) exit 
      end do ! k

      zr = zc + l2
      do k2 = kc,ke
         if ( zr < zh(k2  ) ) exit 
      end do ! k

   end subroutine FindCellIndexBdryConvertReg

   subroutine FindCellIndexBdryConvertRegUnifMesh(l,i1,ic,i2,j1,jc,j2,k1,kc,k2)
      implicit none

      real(8), intent (in) :: l
      integer, intent (in) :: ic,jc,kc
      integer, intent(out) :: i1,i2,j1,j2,k1,k2

      real(8) :: l2
      integer :: ncl2

      l2 = 0.5d0*l
      ncl2 = NINT(l2/dx(is))
      i1 = ic - ncl2
      i2 = ic + ncl2
      j1 = jc - ncl2
      j2 = jc + ncl2
      k1 = kc - ncl2
      k2 = kc + ncl2

      i1 = max(is,i1)
      i2 = min(ie,i2)
      j1 = max(js,j1)
      j2 = min(je,j2)
      k1 = max(ks,k1)
      k2 = min(ke,k2)
   end subroutine FindCellIndexBdryConvertRegUnifMesh

   subroutine ConvertLPP2VOF(tswap) 
      implicit none
      include 'mpif.h'

      integer, intent (in) :: tswap

      real(8) :: rp,dp,uf,vf,wf,dummyreal,ConvertRegSize,mu2tomu1
      logical :: ConvertDropFlag
      integer :: i1,ic,i2,j1,jc,j2,k1,kc,k2
      integer :: ipart,ipart1,num_part_old
      integer :: i,j,k
      logical :: PartAtBlockEdge = .true.
      real(8) :: Re,vmf,vmp
      integer :: ierr
      real(8) :: volcell
      integer :: num_part_rank

      if ( num_part(rank) > 0 ) then
      num_part_old = num_part(rank)
      do ipart = 1,num_part_old
         ic = parts(ipart,rank)%ic
         jc = parts(ipart,rank)%jc
         kc = parts(ipart,rank)%kc

         ! Check if LPP locates at "LPP region"
         call CheckConvertDropCriteria(parts(ipart,rank)%element%vol, & 
                                       parts(ipart,rank)%element%xc,  & 
                                       parts(ipart,rank)%element%yc,  & 
                                       parts(ipart,rank)%element%zc,  &
                                       ic,jc,kc,          &
                                       ConvertDropFlag,               & 
                                       CriteriaConvertCase,AspRatioSphere)

         ! Check if LPP locates at block boundary 
         dp = (6.d0*parts(ipart,rank)%element%vol/PI)**(1.d0/3.d0)
         rp = dp*0.5d0
         if ( parts(ipart,rank)%element%xc - rp > xh(is)-dx(is) .and. & 
              parts(ipart,rank)%element%xc + rp < xh(ie)        .and. &
              parts(ipart,rank)%element%yc - rp > yh(js)-dy(js) .and. &
              parts(ipart,rank)%element%yc + rp < yh(je)        .and. &
              parts(ipart,rank)%element%zc - rp > zh(ks)-dz(ks) .and. &
              parts(ipart,rank)%element%zc + rp < zh(ke)       )     & 
              PartAtBlockEdge = .false.

         volcell = dx(ic)*dy(jc)*dz(kc)
         if ( .not.ConvertDropFlag .and. .not.PartAtBlockEdge ) then
            write(*,*) 'Particle is converted to drop', ipart,rank,tswap,parts(ipart,rank)%element%vol/volcell
            ! transfer droplet properties to particle
            num_drop(rank) = num_drop(rank) + 1
            drops(num_drop(rank))%element = parts(ipart,rank)%element

            ! Define conversion region
            ConvertRegSize = ConvertRegSizeToDiam*dp

            call FindCellIndexBdryConvertRegUnifMesh(ConvertRegSize, & 
                                             i1,ic,i2,j1,jc,j2,k1,kc,k2)

            call GetFluidProp(parts(ipart,rank)%ic, &
                              parts(ipart,rank)%jc, &
                              parts(ipart,rank)%kc, & 
                              parts(ipart,rank)%element%xc, & 
                              parts(ipart,rank)%element%yc, & 
                              parts(ipart,rank)%element%zc, & 
                              uf,vf,wf,dummyreal,dummyreal,dummyreal,& 
                              parts(ipart,rank)%element%vol,&
                              tswap-parts(ipart,rank)%tstepConvert,& 
                              TwoWayCouplingFlag)

            mu2tomu1 = mu2/mu1
            vmf = sqrt(uf*uf+vf*vf+wf*wf)
            vmp = sqrt(parts(ipart,rank)%element%uc**2.d0 & 
                     + parts(ipart,rank)%element%vc**2.d0 & 
                     + parts(ipart,rank)%element%wc**2.d0)
            Re = rho1*ABS(vmf-vmp)*dp/mu1
            
            call BuildFlowFieldVOF(i1,i2,j1,j2,k1,k2,dp,mu2tomu1, & 
                                   parts(ipart,rank)%element%xc,& 
                                   parts(ipart,rank)%element%yc,& 
                                   parts(ipart,rank)%element%zc,&
                                   parts(ipart,rank)%element%uc,& 
                                   parts(ipart,rank)%element%vc,& 
                                   parts(ipart,rank)%element%wc,&
                                   uf,vf,wf,Re)
            
            ! Remove the particle from the array  
            do ipart1 = ipart+1,num_part(rank)
               parts(ipart1-1,rank) = parts(ipart1,rank)
            end do ! ipart1
            num_part(rank) = num_part(rank) - 1
         end if ! ConvertDropFlag

! TEMPORARY - remove particle within liquid (failed to capture its impact on interface)
!         if ( cvof(ic,jc,kc) == 1.d0 ) then
!            write(*,*) 'Particle is removed since within liquid', ipart,rank,parts(ipart,rank)%element%vol/volcell
!            ! Remove the particle from the array  
!            do ipart1 = ipart+1,num_part(rank)
!               parts(ipart1-1,rank) = parts(ipart1,rank)
!            end do ! ipart1
!            num_part(rank) = num_part(rank) - 1
!         end if ! cvof(ic,jc,kc)
! END TEMPORARY

         ! Note: num_part(rank) has been updated
         if ( ipart > num_part(rank) ) exit
      end do ! ipart
      end if ! num_part(rank)

      ! Update num_part to all ranks
      num_part_rank = num_part(rank)
      call MPI_ALLGATHER(num_part_rank, 1, MPI_INTEGER, &
                         num_part(:)  , 1, MPI_INTEGER, MPI_Comm_World, ierr)
   end subroutine ConvertLPP2VOF

   subroutine BuildFlowFieldVOF(i1,i2,j1,j2,k1,k2,dp,mu,xp,yp,zp,up,vp,wp,uf,vf,wf,Re)
      implicit none

      integer, intent (in) :: i1,i2,j1,j2,k1,k2
      real(8), intent (in) :: dp,mu,xp,yp,zp,up,vp,wp,uf,vf,wf,Re

      integer :: i,j,k,i0,j0,k0
      real(8) :: radp,c,stencil3x3(-1:1,-1:1,-1:1)
      real(8) :: rm,rm2,rm3,rx,ry,rz,rx2,ry2,rz2
      real(8) :: um,ux,uy,uz
      real(8) :: tm,tx,ty,tz,theta
      real(8) :: a1,a2,a3,ur,ut,urp,utp
      integer :: nflag
      real(8) :: ufnew,vfnew,wfnew,wt
      logical :: UseCreepingFlow, UseAddedMass
! TEMPORARY
      real(8) :: expterm,term,cam
! END TEMPORARY 

! TEMPORARY  
         UseCreepingFlow = .false.
!         write(*,*) ' Creeping flow solution is used for LPP2VOF conversion!'
         UseAddedMass = .true.
!         write(*,*) ' Added mass compensation is used for LPP2VOF conversion!'
! END TEMPORARY

      radp = dp/2.d0
      a1 = radp*(2.d0+3.d0*mu)/2.d0/(1.d0+mu)
      a2 = radp**3.d0*mu      /2.d0/(1.d0+mu)
      a3 = 2.d0*(1.d0+mu)
      ux=uf-up; uy=vf-vp; uz=wf-wp
      um = sqrt(ux*ux+uy*uy+uz*uz)

      do k=k1+1,k2-1; do j=j1+1,j2-1; do i=i1+1,i2-1 
         
         ! Build VOF field of droplet
         do i0=-1,1; do j0=-1,1; do k0=-1,1
            stencil3x3(i0,j0,k0) = radp**2.d0 & 
               - ((x(i+i0)-xp)**2.d0+(y(j+j0)-yp)**2.d0+(z(k+k0)-zp)**2.d0)
         enddo; enddo; enddo
         call ls2vof_in_cell(stencil3x3,c,nflag)
         vof_flag(i,j,k) = nflag
         cvof(i,j,k)     = c
         if ( cvof(i,j,k) > 1.d0 ) cvof(i,j,k) = 1.d0
         if ( cvof(i,j,k) < 0.d0 ) cvof(i,j,k) = 0.d0
        
         ! Build flow field 
         if ( UseCreepingFlow ) then 
            ! Use creeping flow solution for outside of droplet
            ! XXX Note: the velocity is now computed at cell centers, need to be
            ! improved later

            ! transfer to spherical polar coordinate attached to droplet
            rx=x(i)-xp; ry=y(j)-yp; rz=z(k)-zp
            rm = sqrt(rx*rx+ry*ry+rz*rz)
            theta = acos( (-rx*ux-ry*uy-rz*uz)/rm/um )
            rm2 = rm*rm
            rm3 = rm2*rm
         
            ! Creeping flow solution (Bubbles,Drops, & Particles, Clift et al)
            if ( cvof(i,j,k) == 0.d0 ) then ! outside of droplet
               !Stokes solution 
               ur = -um*cos(theta)*(1.d0 - a1/rm      + a2/rm3)
               ut =  um*sin(theta)*(1.d0 - a1/rm/2.d0 - a2/rm3/2.d0)
               
               ! Oseen solution
               !term = rm*Re/4.d0/radp
               !expterm = exp(-term*(1.d0+cos(theta)))
               !ur = -um*cos(theta)*(1.d0 + a2/rm3) & 
               !    + um*a1/rm/2.d0*(term*(1.d0-expterm) - (1.d0-cos(theta))*expterm)
               !ut =  um*sin(theta)*(1.d0 - a2/rm3/2.d0 - a1/rm/2.d0*expterm )
            !else ! inside of droplet
            !   urp = um*cos(theta)/a3*(1.d0-     rm2/radp/radp)
            !   utp =-um*sin(theta)/a3*(1.d0-2.d0*rm2/radp/radp)
            !end if ! rm

               ! Transfer back to cartestian coordinate
               rx2 = rx*rx
               ry2 = ry*ry
               rz2 = rz*rz
               tx  = (rz2+ry2)*ux - rx*(rz*uz+ry*uy)
               ty  = (rx2+rz2)*uy - ry*(rx*ux+rz*uz)
               tz  = (ry2+rx2)*uz - rz*(ry*uy+rx*ux)
               tm  = sqrt(tx*tx + ty*ty + tz*tz)

               ufnew = ur*rx/rm + ut*tx/tm + up 
               vfnew = ur*ry/rm + ut*ty/tm + vp
               wfnew = ur*rz/rm + ut*tz/tm + wp

               u(i,j,k) = ufnew
               v(i,j,k) = vfnew
               w(i,j,k) = wfnew
               
               ! smoothening the transition from the original value
!               wt = exp(-(rm-radp)**2.d0/(0.5*radp)**2.d0)
!               u(i,j,k) = ufnew*wt + u(i,j,k)*(1.d0-wt)
!               v(i,j,k) = vfnew*wt + v(i,j,k)*(1.d0-wt)
!               w(i,j,k) = wfnew*wt + w(i,j,k)*(1.d0-wt)
            else 
               u(i,j,k) = up
               v(i,j,k) = vp
               w(i,j,k) = wp
            end if ! cvof
         else if ( UseAddedMass ) then
            cam = 1.d0 + 4.d0*rho1/rho2
            if((cvof(i,j,k) + cvof(i+1,j,k)) > 0.0d0) u(i,j,k) = cam*(up-uf)+uf
            if((cvof(i,j,k) + cvof(i,j+1,k)) > 0.0d0) v(i,j,k) = cam*(vp-vf)+vf
            if((cvof(i,j,k) + cvof(i,j,k+1)) > 0.0d0) w(i,j,k) = cam*(wp-wf)+wf
         else ! Simply replace inside of droplet with point velocity
            if((cvof(i,j,k) + cvof(i+1,j,k)) > 0.0d0) u(i,j,k) = up
            if((cvof(i,j,k) + cvof(i,j+1,k)) > 0.0d0) v(i,j,k) = vp
            if((cvof(i,j,k) + cvof(i,j,k+1)) > 0.0d0) w(i,j,k) = wp
         end if ! UseCreepSolution
      end do; end do; end do
      
   end subroutine BuildFlowFieldVOF

   subroutine VOFCommGhost

      implicit none
      include 'mpif.h'
      integer :: ierr, req(12),sta(MPI_STATUS_SIZE,12)
      integer , parameter :: ngh=2

      call ighost_x(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
      call ighost_y(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
      call ighost_z(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       
      call ghost_x(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
      call ghost_y(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
      call ghost_z(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
   end subroutine VOFCommGhost

   subroutine ComputePartForce(tswap,time)

      implicit none
      integer, intent(in) :: tswap
      real(8), intent(in) :: time

      real(8), parameter :: Cm = 0.5d0

      real(8) :: relvel(4), partforce(3)
      real(8) :: dp, mp,Rep, muf, mup, rhof, rhop, taup
      real(8) :: xp,yp,zp, up,vp,wp, uf,vf,wf, DufDt,DvfDt,DwfDt
      real(8) :: fhx,fhy,fhz
      
      real(8) :: dp2Lx

      integer :: ipart
      if ( num_part(rank) > 0 ) then
         do ipart = 1,num_part(rank)
            xp = parts(ipart,rank)%element%xc
            yp = parts(ipart,rank)%element%yc 
            zp = parts(ipart,rank)%element%zc 
            up = parts(ipart,rank)%element%uc
            vp = parts(ipart,rank)%element%vc
            wp = parts(ipart,rank)%element%wc

            rhof = rho1
            rhop = rho2
            muf  = mu1
            mup  = mu2

            dp = (parts(ipart,rank)%element%vol*6.d0/PI)**(1.d0/3.d0)
            taup = rhop *dp*dp/18.0d0/muf & 
                 *(3.d0 + 3.d0*muf/mup)/(3.d0 + 2.d0*muf/mup)

            call GetFluidProp(parts(ipart,rank)%ic, &
                              parts(ipart,rank)%jc, &
                              parts(ipart,rank)%kc, &
                              xp,yp,zp, & 
                              uf,vf,wf,DufDt,DvfDt,DwfDt,   & 
                              parts(ipart,rank)%element%vol,&
                              tswap-parts(ipart,rank)%tstepConvert,& 
                              TwoWayCouplingFlag)

            if ( taup < 5.d0*dt ) then 
               ! If taup is much smaller than dt, no need to compute force,
               ! particle velocity is set to fluid velocity directly
               partforce(1:3) = CRAZY_REAL
               parts(ipart,rank)%element%uc = uf 
               parts(ipart,rank)%element%vc = vf
               parts(ipart,rank)%element%wc = wf
            else 

            relvel(1) = uf - up
            relvel(2) = vf - vp
            relvel(3) = wf - wp 
            relvel(4) = sqrt(relvel(1)**2.0d0 + relvel(2)**2.0d0 + relvel(3)**2.0)
            Rep  = rhof*relvel(4)*dp/muf

            if ( WallEffectSettling ) then 
               dp2Lx = dp/xLength   ! assuming Lx=Lz, y is settling direction
               taup  = taup*((1.d0-dp2Lx)/(1.d0-0.33d0*dp2Lx))**2.7d0
               ! correction correlation by De Felice 1996
            end if ! WallEffectSettling

            ! Note: set history to be zero for now
            fhx=0.d0; fhy=0.d0; fhz=0.d0

            if ( UnsteadyPartForce ) then 
               partforce(1) =(relvel(1)/taup*phi(dragmodel,Rep) + (1.d0-rhof/rhop)*Gx  &   
                            + (1.d0+Cm)*rhof/rhop*DufDt                 &  
                            + fhx )/(1.d0+Cm*rhof/rhop)
               partforce(2) =(relvel(2)/taup*phi(dragmodel,Rep) + (1.d0-rhof/rhop)*Gy  &
                            + (1.d0+Cm)*rhof/rhop*DvfDt                 &
                            + fhy )/(1.d0+Cm*rhof/rhop)
               partforce(3) =(relvel(3)/taup*phi(dragmodel,Rep) + (1.d0-rhof/rhop)*Gz  &
                            + (1.d0+Cm)*rhof/rhop*DwfDt                 &
                            + fhz )/(1.d0+Cm*rhof/rhop)
            else 
               partforce(1) = relvel(1)/taup*phi(dragmodel,Rep) + (1.d0-rhof/rhop)*Gx     
               partforce(2) = relvel(2)/taup*phi(dragmodel,Rep) + (1.d0-rhof/rhop)*Gy  
               partforce(3) = relvel(3)/taup*phi(dragmodel,Rep) + (1.d0-rhof/rhop)*Gz 
            end if ! UnsteadyPartForce 

            if ( output_lpp_evolution .and. mod(tswap,nstats)==0 ) then
               call output_LPP_parameter(rank,ipart,xp,yp,zp,up,vp,wp,uf,vf,wf,dp,time)
            end if ! 

            end if ! taup 

            if ( partforce(1) /= partforce(1) ) then 
               write(*,*) ipart,rank,num_part(rank),UnsteadyPartForce,partforce(1:3),relvel(1:3),taup
               write(*,*) parts(ipart,rank)%ic,parts(ipart,rank)%jc,parts(ipart,rank)%kc,parts(ipart,rank)%element%vol
               write(*,*) xp,yp,zp,up,vp,wp
               call lpperror("particle force in x direction is NaN!")
            else if ( partforce(2) /= partforce(2) ) then 
               write(*,*) ipart,rank,num_part(rank),UnsteadyPartForce,partforce(1:3),relvel(1:3),taup
               write(*,*) parts(ipart,rank)%ic,parts(ipart,rank)%jc,parts(ipart,rank)%kc,parts(ipart,rank)%element%vol
               write(*,*) xp,yp,zp,up,vp,wp
               call lpperror("particle force in y direction is NaN!")
            else if ( partforce(3) /= partforce(3) ) then 
               write(*,*) ipart,rank,num_part(rank),UnsteadyPartForce,partforce(1:3),relvel(1:3),taup
               write(*,*) parts(ipart,rank)%ic,parts(ipart,rank)%jc,parts(ipart,rank)%kc,parts(ipart,rank)%element%vol
               write(*,*) xp,yp,zp,up,vp,wp
               call lpperror("particle force in z direction is NaN!")
            else 
               parts(ipart,rank)%element%duc = partforce(1) 
               parts(ipart,rank)%element%dvc = partforce(2)
               parts(ipart,rank)%element%dwc = partforce(3)
            end if ! partforce(1)

            parts(ipart,rank)%uf = uf
            parts(ipart,rank)%vf = vf
            parts(ipart,rank)%wf = wf

            if ( TwoWayCouplingFlag == TwoWayCouplingFilterForce) then 
               parts_collect(ipart,rank)%xc = xp
               parts_collect(ipart,rank)%yc = yp
               parts_collect(ipart,rank)%zc = zp
               parts_collect(ipart,rank)%ic = parts(ipart,rank)%ic
               parts_collect(ipart,rank)%jc = parts(ipart,rank)%jc
               parts_collect(ipart,rank)%kc = parts(ipart,rank)%kc
               mp =  rho2*parts(ipart,rank)%element%vol
               parts_collect(ipart,rank)%hfx = &
                  mp*(partforce(1) - (1.d0-rho1/rho2)*Gx )  
               parts_collect(ipart,rank)%hfy = &
                  mp*(partforce(2) - (1.d0-rho1/rho2)*Gy )  
               parts_collect(ipart,rank)%hfz = &
                  mp*(partforce(3) - (1.d0-rho1/rho2)*Gz ) 
               parts_collect(ipart,rank)%vol = parts(ipart,rank)%element%vol
            end if ! TwoWayCouplingFlag
         end do ! ipart
      end if ! num_part(rank) 
   end subroutine ComputePartForce

   subroutine output_LPP_parameter(rank,ipart,xp,yp,zp,up,vp,wp,uf,vf,wf,dp,time)

      integer, intent(in) :: rank, ipart
      real(8), intent(in) :: xp,yp,zp,up,vp,wp,uf,vf,wf,dp,time

      character(len=30) :: rootname
      real(8) :: urel2, Re,We

      urel2 = (up-uf)*(up-uf) + (vp-vf)*(vp-vf) + (wp-wf)*(wp-wf)
      Re = sqrt(urel2)*rho1*dp/(mu1+1.d-16)
      We = urel2*rho1*dp/(sigma+1.d-16)

      rootname=trim(out_path)//'/LPP-PARA'//TRIM(int2text(rank,padding))//'-'
      OPEN(UNIT=11,FILE=TRIM(rootname)//TRIM(int2text(ipart,padding))//'.dat',position='append')
      write(11,'(13(E15.8,1X))') time,xp,yp,zp,up,vp,wp,uf,vf,wf,Re,We,dp
      CLOSE(11)
   
   end subroutine output_LPP_parameter

   function phi(dragmodel,Rep)

      ! Note: Compute finite Reynolds number correction on quasi-steady drag
      implicit none
      integer, intent(in) :: dragmodel
      real(8), intent(in) :: Rep
      real(8) :: phi

      select case ( dragmodel ) 
         case ( dragmodel_Stokes ) 
            phi = 1.d0
         case ( dragmodel_SN ) 
            phi = 1.d0+0.15d0*Rep**0.687d0
         case ( dragmodel_CG ) 
            phi = 1.d0+0.15d0*Rep**0.687d0 & 
                  + 1.75d-2*Rep/(1.0d0 + 4.25d4/Rep**1.16d0)
         case ( dragmodel_MKL )
            phi = 1.d0 + 1.0d0/(8.d0/Rep & 
                              + 0.5d0 *(1.d0 + 3.315d0/Rep**0.5d0))
         case default
            call lpperror("wrong quasi-steady drag model!")
      end select ! dragmodel
   end function phi

   subroutine GetFluidProp(ip,jp,kp,xp,yp,zp,uf,vf,wf,DufDt,DvfDt,DwfDt,volp,iStepConversion,TwoWayCouplingFlag) 

      integer, intent(in)  :: ip,jp,kp,iStepConversion
      real(8), intent(in)  :: xp,yp,zp 
      real(8), intent(in)  :: volp 
      real(8), intent(out) :: uf,vf,wf,DufDt,DvfDt,DwfDt
      integer, intent(in)  :: TwoWayCouplingFlag
      
      real(8) :: dp, dxi,dyj,dzk,max_gridsize,min_gridsize
      real(8) :: Lx,Ly,Lz
      real(8), parameter :: small = 1.d0-40
      real(8), parameter :: large = 1.d0+40

      real(8) :: ConvertRegSize
      real(8) :: xh1,xh2,x1,x2,yh1,yh2,y1,y2,zh1,zh2,z1,z2
      integer :: i1,i2,j1,j2,k1,k2

      ! Check particle location
      if ( xp > xh(imax) .or. xp < (xh(imin)-dx(is)) ) then 
         write(*,*) ip,jp,kp,xp,yp,zp,(xh(imin)-dx(is)),xh(imax)
         write(*,*) rank,is,ie,js,je,ks,ke,volp
         call lpperror("Failed to get fluid properties as particle out x range")
      else if ( yp > yh(jmax) .or. yp < (yh(jmin)-dy(js)) ) then
         write(*,*) ip,jp,kp,xp,yp,zp,(yh(jmin)-dy(js)),yh(jmax)
         write(*,*) rank,is,ie,js,je,ks,ke,volp
         call lpperror("Failed to get fluid properties as particle out y range")
      else if ( zp > zh(kmax) .or. zp < (zh(kmin)-dz(ks)) ) then 
         write(*,*) ip,jp,kp,xp,yp,zp,(zh(kmin)-dz(ks)),zh(kmax)
         write(*,*) rank,is,ie,js,je,ks,ke,volp
         call lpperror("Failed to get fluid properties as particle out z range")
      end if ! ip,jp,kp

      dp = (6.d0*volp/PI)**OneThird

      if ( TwoWayCouplingFlag == TwoWayCouplingIgnore .or. & 
           TWoWayCouplingFlag == TwoWayCouplingFilterForce )  then 

         ConvertRegSize = ConvertRegSizeToDiam*dp
         if (  nStepConverion > 0 .and. & 
              iStepConversion < nStepConverion .and. &  
              ConvertRegSize > 2.d0*dx(ip) ) then
!            call ComputeAveFluidVel(xp,yp,zp,ip,jp,kp,uf,vf,wf,DufDt,DvfDt,DwfDt,ConvertRegSize)
            call FindCellIndexBdryConvertRegUnifMesh(ConvertRegSize,i1,ip,i2,j1,jp,j2,k1,kp,k2)
            xh1 = xp-xh(i1); xh2 = xh(i2)-xp; x1  = xp-x(i1); x2  = x(i2)-xp
            yh1 = yp-yh(j1); yh2 = yh(j2)-yp; y1  = yp-y(j1); y2  = y(j2)-yp            
            zh1 = zp-zh(k1); zh2 = zh(k2)-zp; z1  = zp-z(k1); z2  = z(k2)-zp
            call Surface2VolumeIntrpl(u(i1,jp,kp),u(i2,jp,kp), &
                                      u(ip,j1,kp),u(ip,j2,kp), &
                                      u(ip,jp,k1),u(ip,jp,k2), &
                                      xh1,xh2,y1,y2,z1,z2,uf)
            call Surface2VolumeIntrpl(v(i1,jp,kp),v(i2,jp,kp), &
                                      v(ip,j1,kp),v(ip,j2,kp), &
                                      v(ip,jp,k1),v(ip,jp,k2), &
                                      x1,x2,yh1,yh2,z1,z2,vf)
            call Surface2VolumeIntrpl(w(i1,jp,kp),w(i2,jp,kp), &
                                      w(ip,j1,kp),w(ip,j2,kp), &
                                      w(ip,jp,k1),w(ip,jp,k2), &
                                      x1,x2,y1,y2,zh1,zh2,wf)
            if ( UnsteadyPartForce ) then 
               call Surface2VolumeIntrpl(sdu(i1,jp,kp),sdu(i2,jp,kp), &
                                         sdu(ip,j1,kp),sdu(ip,j2,kp), &
                                         sdu(ip,jp,k1),sdu(ip,jp,k2), &
                                         xh1,xh2,y1,y2,z1,z2,DufDt)
               call Surface2VolumeIntrpl(sdv(i1,jp,kp),sdv(i2,jp,kp), &
                                         sdv(ip,j1,kp),sdv(ip,j2,kp), &
                                         sdv(ip,jp,k1),sdv(ip,jp,k2), &
                                         x1,x2,yh1,yh2,z1,z2,DvfDt)
               call Surface2VolumeIntrpl(sdw(i1,jp,kp),sdw(i2,jp,kp), &
                                         sdw(ip,j1,kp),sdw(ip,j2,kp), &
                                         sdw(ip,jp,k1),sdw(ip,jp,k2), &
                                         x1,x2,y1,y2,zh1,zh2,DwfDt)
            end if ! UnsteadyPartForce 
         else 
            call TrilinearIntrplFluidVel(xp,yp,zp,ip,jp,kp,uf,DufDt,1)
            call TrilinearIntrplFluidVel(xp,yp,zp,ip,jp,kp,vf,DvfDt,2)
            call TrilinearIntrplFluidVel(xp,yp,zp,ip,jp,kp,wf,DwfDt,3)
         end if
      else if ( TwoWayCouplingFlag == TwoWayCouplingFilterVel ) then 
         call lpperror("Two way coupling with filtered fluid veloicty is not ready yet!")
      end if ! dp 

   end subroutine GetFluidProp

   subroutine TwoWayCouplingForce(TwoWayCouplingFlag) 

      implicit none 
      include 'mpif.h'
      integer, intent(in) :: TwoWayCouplingFlag
     
      integer :: ierr,irank
      integer :: MPI_part_collect_type,MPI_part_collect_row, & 
                 oldtypes(0:3), blockcounts(0:3), & 
                 offsets(0:3), intextent,r8extent

      integer :: i,j,k,ipart,ip,jp,kp
      real(8) :: xp,yp,zp,fx,fy,fz,x2,y2,z2,s0,s1,dp
      real(8) :: massLPP, massVOF, mfLPP
      real(8), parameter :: mfLPPsmall = 1.0d-40
      real :: erf

      if (TwoWayCouplingFlag == TwoWayCouplingFilterForce) then
         
         !  Setup MPI derived type for particle 
         call MPI_TYPE_EXTENT(MPI_REAL8,   r8extent,  ierr) 
         call MPI_TYPE_EXTENT(MPI_INTEGER, intextent, ierr)

         offsets    (0) = 0 
         oldtypes   (0) = MPI_REAL8 
         blockcounts(0) = 7 
         offsets    (1) = offsets(0) + blockcounts(0)*r8extent 
         oldtypes   (1) = MPI_INTEGER  
         blockcounts(1) = 4  
         
         call MPI_TYPE_STRUCT(2, blockcounts, offsets, oldtypes, & 
                              MPI_part_collect_type, ierr) 
         call MPI_TYPE_COMMIT(MPI_part_collect_type, ierr)
         call MPI_TYPE_CONTIGUOUS (max_num_part, MPI_part_collect_type, MPI_part_collect_row, ierr)
         call MPI_TYPE_COMMIT(MPI_part_collect_row, ierr)

         parts_collect_rank(1:max_num_part) = parts_collect(1:max_num_part,rank)
         call MPI_ALLGATHER(parts_collect_rank(1:max_num_part), 1, MPI_part_collect_row, &
                            parts_collect(1:max_num_part,:),    1, MPI_part_collect_row, & 
                            MPI_Comm_World, ierr)

         ! Compute LPP mass fraction
         massLPP = 0.d0 
         if ( num_part(rank) > 0 ) then 
            do i = 1,num_part(rank)
               massLPP = massLPP + parts(ipart,rank)%element%vol 
            end do ! i
            massLPP = massLPP*rho2
         end if ! num_part(rank)
         massVOF = 0.d0 
         do k=ks,ke; do j=js,je; do i=is,ie
            massVOF = massVOF + rho(i,j,k)*dx(i)*dy(j)*dz(k) 
         end do; end do; end do

         do k=ks,ke; do j=js,je; do i=is,ie
            do irank = 0,nPdomain-1
               if ( num_part(irank) > 0 ) then
               do ipart = 1,num_part(irank) 
                  !ip = parts_collect(ipart,irank)%ic
                  !jp = parts_collect(ipart,irank)%jc
                  !kp = parts_collect(ipart,irank)%kc

                  mfLPP = (massLPP-parts_collect(ipart,irank)%vol*rho2)/massVOF
                  mfLPP = max(mfLPP,mfLPPsmall) 
                  dp = (6.d0*parts_collect(ipart,irank)%vol/PI)**OneThird 
                  LengthLPP = LengthLPP2dp*dp/erf(mfLPP/mfLPP_ref)
                  s1 = LengthLPP
                  
                  if ( i > Ng .and. i < Ng+Nx .and. & 
                       j > Ng .and. j < Ng+Ny .and. &
                       k > Ng .and. k < Ng+Nz ) then 
                     xp = parts_collect(ipart,irank)%xc
                     yp = parts_collect(ipart,irank)%yc
                     zp = parts_collect(ipart,irank)%zc
                     fx = parts_collect(ipart,irank)%hfx
                     fy = parts_collect(ipart,irank)%hfy
                     fz = parts_collect(ipart,irank)%hfz
   
                     x2 = (xh(i)-xp)*(xh(i)-xp) + (y(j)-yp)*(y(j)-yp) + (z(k)-zp)*(z(k)-zp)
                     y2 = (x(i)-xp)*(x(i)-xp) + (yh(j)-yp)*(yh(j)-yp) + (z(k)-zp)*(z(k)-zp)
                     z2 = (x(i)-xp)*(x(i)-xp) + (y(j)-yp)*(y(j)-yp) + (zh(k)-zp)*(zh(k)-zp)

                     s0 = (2.d0*PI*s1*s1)**(-1.5d0)
                     du(i,j,k) = du(i,j,k) - fx*s0*exp(-x2*0.5d0/s1/s1)
                     dv(i,j,k) = dv(i,j,k) - fy*s0*exp(-y2*0.5d0/s1/s1)
                     dw(i,j,k) = dw(i,j,k) - fz*s0*exp(-z2*0.5d0/s1/s1)
                  end if ! i,j,k
               end do ! ipart
               end if ! num_part(irank)
            end do ! irank
         end do; end do; end do
      else 
         call lpperror("Unknow two-way coupling flag!")
      end if ! TwoWayCouplingFlag

   end subroutine TwoWayCouplingForce

   subroutine UpdatePartSol(iStage)
      implicit none
      integer, intent (in) :: iStage
      integer :: ipart
      real(8) :: volcell

      if ( num_part(rank) > 0 ) then 
         do ipart = 1,num_part(rank)
! TEMPORARY - Not update particle position if volume too small   
            volcell = dx(parts(ipart,rank)%ic) &
                     *dy(parts(ipart,rank)%jc) &
                     *dz(parts(ipart,rank)%kc)
            if ( parts(ipart,rank)%element%vol > volcell/512.d0 ) then 
               parts(ipart,rank)%element%xc = parts(ipart,rank)%element%xc +& 
                                              parts(ipart,rank)%element%uc*dt 
               parts(ipart,rank)%element%yc = parts(ipart,rank)%element%yc +& 
                                              parts(ipart,rank)%element%vc*dt 
               parts(ipart,rank)%element%zc = parts(ipart,rank)%element%zc +& 
                                              parts(ipart,rank)%element%wc*dt
            end if ! parts(ipart,rank)%element%vol
! END TEMPORARY    

! TEMPORARY - Check particle location inside domain or not  
            if ( parts(ipart,rank)%element%xc < xh(1)-dx(1) .or. & 
                 parts(ipart,rank)%element%xc > xh(Nxt) ) then
               write(*,*) ipart,rank,iStage,parts(ipart,rank)%xcOld,parts(ipart,rank)%element%xc,& 
                           parts(ipart,rank)%element%uc,parts(ipart,rank)%uf, & 
                           parts(ipart,rank)%element%duc, parts(ipart,rank)%element%vol
               call lpperror("Particle location is out of bound in x direction!") 
            else if ( parts(ipart,rank)%element%yc < yh(1)-dy(1) .or. & 
                      parts(ipart,rank)%element%yc > yh(Nyt) ) then
               write(*,*) ipart,rank,iStage,parts(ipart,rank)%ycOld,parts(ipart,rank)%element%yc,& 
                           parts(ipart,rank)%element%vc,parts(ipart,rank)%element%dvc,&
                           parts(ipart,rank)%element%vol
               call lpperror("Particle location is out of bound in y direction!") 
            else if ( parts(ipart,rank)%element%zc < zh(1)-dz(1) .or. & 
                      parts(ipart,rank)%element%zc > zh(Nzt) ) then
               write(*,*) ipart,rank,iStage,parts(ipart,rank)%zcOld,parts(ipart,rank)%element%zc,& 
                           parts(ipart,rank)%element%wc,parts(ipart,rank)%element%dwc,&
                           parts(ipart,rank)%element%vol
               call lpperror("Particle location is out of bound in z direction!") 
            end if ! i1
! END TEMPORARY
            if ( parts(ipart,rank)%element%duc /= CRAZY_REAL ) then 
               parts(ipart,rank)%element%uc = parts(ipart,rank)%element%uc + &
                                              parts(ipart,rank)%element%duc*dt 
               parts(ipart,rank)%element%vc = parts(ipart,rank)%element%vc +&
                                              parts(ipart,rank)%element%dvc*dt 
               parts(ipart,rank)%element%wc = parts(ipart,rank)%element%wc +&
                                              parts(ipart,rank)%element%dwc*dt
            end if ! parts(ipart.rank)%element%duc
         end do ! ipart
        
         if ( iStage == 1 ) call UpdatePartLocCell   
      end if ! num_part(rank)
   end subroutine UPdatePartSol

   subroutine StoreOldPartSol()
      implicit none

      if ( num_part(rank) > 0 ) then 
         parts(1:num_part(rank),rank)%xcOld = parts(1:num_part(rank),rank)%element%xc 
         parts(1:num_part(rank),rank)%ycOld = parts(1:num_part(rank),rank)%element%yc 
         parts(1:num_part(rank),rank)%zcOld = parts(1:num_part(rank),rank)%element%zc 
   
         parts(1:num_part(rank),rank)%ucOld = parts(1:num_part(rank),rank)%element%uc 
         parts(1:num_part(rank),rank)%vcOld = parts(1:num_part(rank),rank)%element%vc 
         parts(1:num_part(rank),rank)%wcOld = parts(1:num_part(rank),rank)%element%wc 
      end if ! num_part(rank)

   end subroutine StoreOldPartSol

   subroutine AveragePartSol()
      implicit none
     
      ! Update particle solution for 2nd order time integration
      parts(1:num_part(rank),rank)%element%xc = & 
         0.5d0*( parts(1:num_part(rank),rank)%element%xc + parts(1:num_part(rank),rank)%xcOld ) 
      parts(1:num_part(rank),rank)%element%yc = & 
         0.5d0*( parts(1:num_part(rank),rank)%element%yc + parts(1:num_part(rank),rank)%ycOld ) 
      parts(1:num_part(rank),rank)%element%zc = & 
         0.5d0*( parts(1:num_part(rank),rank)%element%zc + parts(1:num_part(rank),rank)%zcOld )  
      
      parts(1:num_part(rank),rank)%element%uc = & 
         0.5d0*( parts(1:num_part(rank),rank)%element%uc + parts(1:num_part(rank),rank)%ucOld ) 
      parts(1:num_part(rank),rank)%element%vc = & 
         0.5d0*( parts(1:num_part(rank),rank)%element%vc + parts(1:num_part(rank),rank)%vcOld ) 
      parts(1:num_part(rank),rank)%element%wc = & 
         0.5d0*( parts(1:num_part(rank),rank)%element%wc + parts(1:num_part(rank),rank)%wcOld )  

      call UpdatePartLocCell   
      
   end subroutine AveragePartSol

   subroutine PartBCWrapper
      ! Transfer particles crossing blocks 
      if ( nPdomain > 1 ) then  
         call CollectPartCrossBlocks   
         call TransferPartCrossBlocks   
      end if ! nPdomain

      call SetPartBC
   end subroutine PartBCWrapper

   subroutine UpdatePartLocCell   
      implicit none
      
      integer :: i,j,k,ipart,i1,j1,k1
      real(8) :: xp,yp,zp

      ! A fast version for uniform mesh
      do ipart = 1,num_part(rank)
         ! x direction 
         i  = parts(ipart,rank)%ic
         xp = parts(ipart,rank)%element%xc
         i1 = INT((xp + dble(Ng)*dx(2))/dx(2)) + 1
         if ( i1 < 1 .or. i1 > Nxt ) call lpperror("Cell index is out of bound in x direction!") 
!         if      ( i <= Ng+1 ) then
!            i1 = INT((xp + dble(Ng)*dx(2))/dx(2)) + 1
!         else if ( i > Nx+Ng ) then 
!            i1 = Nx+2*Ng - INT((xh(Nx+2*Ng)-xp)/dx(Nx+Ng+1)) 
!         else if ( xp > xh(i-1) .and. xp <= xh(i  ) ) then
!            i1 = i
!         else if ( xp > xh(i  ) .and. xp <= xh(i+1) ) then 
!            i1 = i+1
!         else if ( xp > xh(i+1) .and. xp <= xh(i+2) ) then 
!            i1 = i+2
!         else if ( xp > xh(i-2) .and. xp <= xh(i-1) ) then 
!            i1 = i-1
!         else if ( xp > xh(i-3) .and. xp <= xh(i-2) ) then 
!            i1 = i-2
!         else 
!            call lpperror("Particle moves out of tracking range in x direction !")
!         end if !xp
        
         ! y direction 
         j  = parts(ipart,rank)%jc
         yp = parts(ipart,rank)%element%yc 
         j1 = INT((yp + dble(Ng)*dy(2))/dy(2)) + 1
         if ( j1 < 1 .or. j1 > Nyt ) call lpperror("Cell index is out of bound in y direction!") 
!         if      ( j <= Ng+1 ) then
!            j1 = INT((yp + dble(Ng)*dy(2))/dy(2)) + 1
!         else if ( j > Ny+Ng ) then 
!            j1 = Ny+2*Ng - INT((yh(Ny+2*Ng)-yp)/dy(Ny+Ng+1)) 
!         else if ( yp > yh(j-1) .and. yp <= yh(j  ) ) then
!            j1 = j
!         else if ( yp > yh(j  ) .and. yp <= yh(j+1) ) then 
!            j1 = j+1
!         else if ( yp > yh(j+1) .and. yp <= yh(j+2) ) then 
!            j1 = j+2
!         else if ( yp > yh(j-2) .and. yp <= yh(j-1) ) then 
!            j1 = j-1
!         else if ( yp > yh(j-3) .and. yp <= yh(j-2) ) then 
!            j1 = j-2
!         else
!            call lpperror("Particle moves out of tracking range in y direction !")
!         end if !yp

         ! z direction 
         k  = parts(ipart,rank)%kc
         zp = parts(ipart,rank)%element%zc 
         k1 = INT((zp + dble(Ng)*dz(2))/dz(2)) + 1
         if ( k1 < 1 .or. k1 > Nzt ) call lpperror("Cell index is out of bound in z direction!") 
!         if      ( k <= Ng+1 ) then
!            k1 = INT((zp + dble(Ng)*dz(2))/dz(2)) + 1
!         else if ( k > Nz+Ng ) then 
!            k1 = Nz+2*Ng - INT((zh(Nz+2*Ng)-zp)/dz(Nz+Ng+1)) 
!         else if ( zp > zh(k-1) .and. zp <= zh(k  ) ) then
!            k1 = k
!         else if ( zp > zh(k  ) .and. zp <= zh(k+1) ) then 
!            k1 = k+1
!         else if ( zp > zh(k+1) .and. zp <= zh(k+2) ) then 
!            k1 = k+2
!         else if ( zp > zh(k-2) .and. zp <= zh(k-1) ) then 
!            k1 = k-1
!         else if ( zp > zh(k-3) .and. zp <= zh(k-2) ) then 
!            k1 = k-2
!         else 
!            call lpperror("Particle moves out of tracking range in z direction !")
!         end if !zp

         parts(ipart,rank)%ic = i1 
         parts(ipart,rank)%jc = j1
         parts(ipart,rank)%kc = k1 

      end do ! ipart
   end subroutine UpdatePartLocCell   
   
   subroutine CollectPartCrossBlocks
      implicit none
       
      integer :: ipart,i,j,k
      integer :: ranknew
      integer :: c1,c2,c3
      logical :: PartNeedTransfer, PartExitDomain 

      allocate( num_part_cross(0:nPdomain-1) )
      allocate( parts_cross_id     (max_num_part_cross,0:nPdomain-1) )
      allocate( parts_cross_id_rank(max_num_part_cross) )
      allocate( parts_cross_newrank     (max_num_part_cross,0:nPdomain-1) )
      allocate( parts_cross_newrank_rank(max_num_part_cross) )
      num_part_cross(:) = 0
      parts_cross_id     (:,:) = CRAZY_INT 
      parts_cross_newrank(:,:) = CRAZY_INT 
      if ( num_part(rank) > 0 ) then 
         do ipart = 1,num_part(rank)
            i = parts(ipart,rank)%ic
            j = parts(ipart,rank)%jc
            k = parts(ipart,rank)%kc

            ! Check if the particle need to be transferred across blocks
            PartNeedTransfer = .false.
            PartExitDomain = .false.
            if ( i < is ) then 
               if ( i > Ng ) then 
                  PartNeedTransfer = .true.
               else 
                  if ( lppbdry_cond(1) == 'periodic' ) then 
                     i = i + Nx
                     PartNeedTransfer = .true.
                  else if ( lppbdry_cond(1) == 'exit' ) then 
                     PartExitDomain = .true. 
                  end if ! lppbdry
               end if ! i
            else if ( i > ie ) then 
               if ( i <= Ng+Nx ) then 
                  PartNeedTransfer = .true.
               else 
                  if ( lppbdry_cond(4) == 'periodic' ) then
                     i = i - Nx
                     PartNeedTransfer = .true.
                  else if ( lppbdry_cond(4) == 'exit' ) then 
                     PartExitDomain = .true. 
                  end if ! lppbdry
               end if ! i
            end if ! i

            if ( j < js ) then 
               if ( j > Ng ) then 
                  PartNeedTransfer = .true.
               else 
                  if ( lppbdry_cond(2) == 'periodic' ) then 
                     j = j + Ny
                     PartNeedTransfer = .true.
                  else if ( lppbdry_cond(2) == 'exit' ) then
                     PartExitDomain = .true. 
                  end if ! lppbdry
               end if ! j
            else if ( j > je ) then 
               if ( j <= Ng+Ny ) then 
                  PartNeedTransfer = .true.
               else 
                  if ( lppbdry_cond(5) == 'periodic' ) then 
                     j = j - Ny
                     PartNeedTransfer = .true.
                  else if ( lppbdry_cond(5) == 'exit' ) then
                     PartExitDomain = .true. 
                  end if ! lppbdry
               end if ! j
            end if ! j

            if ( k < ks ) then 
               if ( k > Ng ) then 
                  PartNeedTransfer = .true.
               else 
                  if ( lppbdry_cond(3) == 'periodic' ) then 
                     k = k + Nz
                     PartNeedTransfer = .true.
                  else if ( lppbdry_cond(3) == 'exit' ) then
                     PartExitDomain = .true. 
                  end if ! lppbdry
               end if ! k
            else if ( k > ke ) then 
               if ( k <= Ng+Nz ) then 
                  PartNeedTransfer = .true.
               else 
                  if ( lppbdry_cond(6) == 'periodic' ) then 
                     k = k - Nz
                     PartNeedTransfer = .true.
                  else if ( lppbdry_cond(6) == 'exit' ) then
                     PartExitDomain = .true. 
                  end if ! lppbdry
               end if ! k 
            end if ! k

            ! if particle locates at boundary cell in any one direction and the
            ! BC is not periodic, then no need to transfer particle 
            if ( PartExitDomain ) PartNeedTransfer = .false.

            ! Note: here only collect and transfer particles which cross blocks 
            !       due to periodic BC, the location and cell index stored in 
            !       parts(ipart,rank) will not be changed until SetPartBC is called

            ! Compute new rank for the particle if transfer is needed
            if ( PartNeedTransfer ) then  
               c1 = (i-Ng-1)/Mx  
               c2 = (j-Ng-1)/My  
               c3 = (k-Ng-1)/Mz  
               ranknew = c1*npy*npz + c2*npz + c3
               if ( ranknew > nPdomain-1 .or. ranknew < 0 ) then
                  call lpperror("new rank of particle out of range!")
               else if ( ranknew /= rank ) then 
                  num_part_cross(rank)  = num_part_cross(rank) + 1
                  parts_cross_id     (num_part_cross(rank),rank) = ipart 
                  parts_cross_newrank(num_part_cross(rank),rank) = ranknew 
               end if ! ranknew
            end if ! PartNeedTransfer

         end do ! ipart
      end if ! num_part(rank)
   end subroutine CollectPartCrossBlocks

   subroutine TransferPartCrossBlocks
      implicit none

      include 'mpif.h'

      integer :: ipart,ipart_cross,ipart1
      integer :: ierr,irank
      integer :: ranknew
      integer :: req(4),sta(MPI_STATUS_SIZE,4)
      integer :: MPI_particle_type, oldtypes(0:3), blockcounts(0:3), & 
                 offsets(0:3), intextent,r8extent
      integer :: maxnum_part_cross, MPI_int_row
      type(element) :: element_NULL
      integer :: num_part_cross_rank,num_part_rank
   
      element_NULL%xc = 0.d0;element_NULL%yc = 0.d0;element_NULL%zc = 0.d0
      element_NULL%uc = 0.d0;element_NULL%vc = 0.d0;element_NULL%wc = 0.d0
      element_NULL%duc = 0.d0;element_NULL%dvc = 0.d0;element_NULL%dwc = 0.d0
      element_NULL%vol = 0.d0;element_NULL%id = CRAZY_INT

      num_part_cross_rank = num_part_cross(rank)
      call MPI_ALLGATHER(num_part_cross_rank, 1, MPI_INTEGER, &
                         num_part_cross,    1, MPI_INTEGER, MPI_Comm_World, ierr)
      maxnum_part_cross = maxval(num_part_cross)
      if ( maxnum_part_cross  > 0 ) then 

         call MPI_TYPE_CONTIGUOUS (maxnum_part_cross, MPI_INTEGER, MPI_int_row, ierr)
         call MPI_TYPE_COMMIT(MPI_int_row, ierr)

         parts_cross_id_rank(1:maxnum_part_cross)=parts_cross_id(1:maxnum_part_cross,rank)
         call MPI_ALLGATHER(parts_cross_id_rank(1:maxnum_part_cross), 1, MPI_int_row, &
                            parts_cross_id(1:maxnum_part_cross,:),    1, MPI_int_row, & 
                            MPI_Comm_World, ierr)
         parts_cross_newrank_rank(1:maxnum_part_cross) = parts_cross_newrank(1:maxnum_part_cross,rank)
         call MPI_ALLGATHER(parts_cross_newrank_rank(1:maxnum_part_cross), 1, MPI_int_row, &
                            parts_cross_newrank(1:maxnum_part_cross,:),    1, MPI_int_row, & 
                            MPI_Comm_World, ierr)
      !  Setup MPI derived type for particle 
      call MPI_TYPE_EXTENT(MPI_REAL8,   r8extent,  ierr) 
      call MPI_TYPE_EXTENT(MPI_INTEGER, intextent, ierr) 
      offsets    (0) = 0 
      oldtypes   (0) = MPI_REAL8 
      blockcounts(0) = 10 
      offsets    (1) = offsets(0) + blockcounts(0)*r8extent 
      oldtypes   (1) = MPI_INTEGER  
      blockcounts(1) = 1  
      offsets    (2) = offsets(1) + blockcounts(1)*intextent 
      oldtypes   (2) = MPI_REAL8  
      blockcounts(2) = 9
      offsets    (3) = offsets(2) + blockcounts(2)*r8extent 
      oldtypes   (3) = MPI_INTEGER  
      blockcounts(3) = 5  

      call MPI_TYPE_STRUCT(4, blockcounts, offsets, oldtypes, & 
                           MPI_particle_type, ierr) 
      call MPI_TYPE_COMMIT(MPI_particle_type, ierr)

      ! Transfer particle to the destination block
      do irank = 0,nPdomain-1
         if ( num_part_cross(irank) > 0 ) then
            do ipart_cross = 1,num_part_cross(irank)
               ipart   = parts_cross_id     (ipart_cross,irank)
               ranknew = parts_cross_newrank(ipart_cross,irank)
               if ( rank == irank ) then 
                  call MPI_ISEND(parts(ipart,irank),1, MPI_particle_type, & 
                                 ranknew, 15, MPI_COMM_WORLD, req(1), ierr)
                  call MPI_WAIT(req(1),sta(:,1),ierr)
               else if ( rank == ranknew ) then 
                  call MPI_IRECV(parts(num_part(ranknew)+1,ranknew),1,MPI_particle_type, & 
                                 irank, 15, MPI_COMM_WORLD, req(2), ierr)
                  call MPI_WAIT(req(2),sta(:,2),ierr)
                  num_part(ranknew) = num_part(ranknew) + 1
               end if ! rank 
            end do ! ipart_cross 
         end if ! num_part_cross(irank)
      end do ! irank

      ! Remove particle in the origin block
      if ( num_part_cross(rank) > 0 ) then
         do ipart_cross = 1,num_part_cross(rank)
            ipart = parts_cross_id(ipart_cross,rank)
            do ipart1 = ipart,num_part(rank)-1
               parts(ipart1,rank) = parts(ipart1+1,rank)
            end do ! ipart1
            parts(num_part(rank),rank)%element = element_NULL
            num_part(rank) = num_part(rank) - 1
            if ( ipart_cross < num_part_cross(rank) ) & 
               parts_cross_id(ipart_cross+1:num_part_cross(rank),rank) = &
               parts_cross_id(ipart_cross+1:num_part_cross(rank),rank) - 1
         end do ! ipart_cross 
      end if ! num_part_cross(irank)

      ! Allgather updated number of partilces in every block
      num_part_rank = num_part(rank)
      call MPI_ALLGATHER(num_part_rank, 1, MPI_INTEGER, &
                         num_part(:)   , 1, MPI_INTEGER, MPI_Comm_World, ierr)
      end if ! maxnum_part_cross

      ! final
      deallocate(num_part_cross)
      deallocate(parts_cross_id)
      deallocate(parts_cross_id_rank)
      deallocate(parts_cross_newrank)
      deallocate(parts_cross_newrank_rank)

   end subroutine TransferPartCrossBlocks

   subroutine SetPartBC
      implicit none
      include 'mpif.h'
      integer :: ipart,ipart1,ierr
      type(element) :: element_NULL
      integer :: num_part_rank
   
      element_NULL%xc = 0.d0;element_NULL%yc = 0.d0;element_NULL%zc = 0.d0
      element_NULL%uc = 0.d0;element_NULL%vc = 0.d0;element_NULL%wc = 0.d0
      element_NULL%duc = 0.d0;element_NULL%dvc = 0.d0;element_NULL%dwc = 0.d0
      element_NULL%vol = 0.d0;element_NULL%id = CRAZY_INT

      if ( num_part(rank) > 0 ) then 
         do ipart = 1,num_part(rank)
            if ( parts(ipart,rank)%ic <= Ng ) then
               if ( lppbdry_cond(1) == 'exit') then 
                  call PartBC_exit(ipart,rank)
                  cycle 
               else 
                  call ImposePartBC(ipart,rank,1)
               end if ! lppbdry_cond
            else if ( parts(ipart,rank)%ic > Ng+Nx ) then 
               if ( lppbdry_cond(4) == 'exit') then 
                  call PartBC_exit(ipart,rank)
                  cycle 
               else 
                  call ImposePartBC(ipart,rank,4)
               end if ! lppbdry_cond
            end if ! ic

            if ( parts(ipart,rank)%jc <= Ng ) then 
               if ( lppbdry_cond(2) == 'exit') then 
                  call PartBC_exit(ipart,rank)
                  cycle 
               else 
                  call ImposePartBC(ipart,rank,2)
               end if ! lppbdry_cond
            else if ( parts(ipart,rank)%jc > Ng+Ny ) then 
               if ( lppbdry_cond(5) == 'exit') then 
                  call PartBC_exit(ipart,rank)
                  cycle 
               else 
                  call ImposePartBC(ipart,rank,5)
               end if ! lppbdry_cond
            end if ! jc 

            if ( parts(ipart,rank)%kc <= Ng ) then  
               if ( lppbdry_cond(3) == 'exit') then 
                  call PartBC_exit(ipart,rank)
                  cycle 
               else 
                  call ImposePartBC(ipart,rank,3)
               end if ! lppbdry_cond
            else if ( parts(ipart,rank)%kc > Ng+Nz ) then 
               if ( lppbdry_cond(6) == 'exit') then 
                  call PartBC_exit(ipart,rank)
                  cycle 
               else 
                  call ImposePartBC(ipart,rank,6)
               end if ! lppbdry_cond
            end if ! kc
         
         end do ! ipart

         ! Remove particles that exit the domain from the list (vol=-1)
         ipart = 0 
         do
            ipart = ipart + 1
            if ( parts(ipart,rank)%element%vol < 0.d0 ) then
               if ( ipart < num_part(rank) ) then 
                  do ipart1 = ipart+1,num_part(rank)
                     parts(ipart1-1,rank) = parts(ipart1,rank)
                  end do ! ipart1
               end if ! ipart
               parts(num_part(rank),rank)%element = element_NULL
               num_part(rank) = num_part(rank) - 1
               ipart          = ipart          - 1 ! Note: the index also needs to be shifted
            end if ! parts(ipart,rank)%element%vol
            if ( ipart >= num_part(rank) ) exit
         end do ! ipart

      end if ! num_part(rank)

      ! Allgather updated number of particles after imposing BC
      num_part_rank = num_part(rank)
      call MPI_ALLGATHER(num_part_rank, 1, MPI_INTEGER, &
                         num_part(:)   , 1, MPI_INTEGER, MPI_Comm_World, ierr)

   end subroutine SetPartBC

   subroutine ImposePartBC(ipart,rank,d)
      implicit none
      integer, intent (in) :: ipart,rank,d

      if ( lppbdry_cond(d) == 'periodic' ) then
         call PartBC_periodic(ipart,rank,d)
      else if ( lppbdry_cond(d) == 'reflect' ) then 
         call PartBC_reflect(ipart,rank,d)
      else 
         call lpperror("Unknown particle bondary condition!")
      end if ! lppbdry_cond
   end subroutine ImposePartBC

   subroutine PartBC_periodic(ipart,rank,d)
      implicit none
      integer, intent (in) :: ipart,rank,d
      
      if ( d == 1 ) then 
         parts(ipart,rank)%ic = parts(ipart,rank)%ic + Nx
         parts(ipart,rank)%element%xc = parts(ipart,rank)%element%xc + xLength
      else if ( d == 4 ) then 
         parts(ipart,rank)%ic = parts(ipart,rank)%ic - Nx
         parts(ipart,rank)%element%xc = parts(ipart,rank)%element%xc - xLength
      else if ( d == 2 ) then 
         parts(ipart,rank)%jc = parts(ipart,rank)%jc + Ny
         parts(ipart,rank)%element%yc = parts(ipart,rank)%element%yc + yLength
      else if ( d == 5 ) then 
         parts(ipart,rank)%jc = parts(ipart,rank)%jc - Ny
         parts(ipart,rank)%element%yc = parts(ipart,rank)%element%yc - yLength
      else if ( d == 3 ) then 
         parts(ipart,rank)%kc = parts(ipart,rank)%kc + Nz
         parts(ipart,rank)%element%zc = parts(ipart,rank)%element%zc + zLength
      else if ( d == 6 ) then 
         parts(ipart,rank)%kc = parts(ipart,rank)%kc - Nz
         parts(ipart,rank)%element%zc = parts(ipart,rank)%element%zc - zLength
      end if ! d
   end subroutine PartBC_periodic

   subroutine PartBC_exit(ipart,rank)
      implicit none
      integer, intent (in) :: ipart,rank
      
      ! For particle that will exit the domain, mark volume=-1 
      parts(ipart,rank)%element%vol = -1.d0
   end subroutine PartBC_exit

   subroutine PartBC_reflect(ipart,rank,d)
      implicit none
      integer, intent (in) :: ipart,rank,d
      
      if ( d == 1 ) then 
         parts(ipart,rank)%ic = 2* Ng    +1-parts(ipart,rank)%ic 
         parts(ipart,rank)%element%xc = -parts(ipart,rank)%element%xc 
         parts(ipart,rank)%element%uc = -parts(ipart,rank)%element%uc 
      else if ( d == 4 ) then 
         parts(ipart,rank)%ic = 2*(Nx+Ng)+1-parts(ipart,rank)%ic 
         parts(ipart,rank)%element%xc = 2.d0*xLength - parts(ipart,rank)%element%xc 
         parts(ipart,rank)%element%uc = -parts(ipart,rank)%element%uc 
      else if ( d == 2 ) then 
         parts(ipart,rank)%jc = 2* Ng    +1-parts(ipart,rank)%jc 
         parts(ipart,rank)%element%yc = -parts(ipart,rank)%element%yc 
         parts(ipart,rank)%element%vc = -parts(ipart,rank)%element%vc 
      else if ( d == 5 ) then 
         parts(ipart,rank)%jc = 2*(Ny+Ng)+1-parts(ipart,rank)%jc 
         parts(ipart,rank)%element%yc = 2.d0*yLength - parts(ipart,rank)%element%yc 
         parts(ipart,rank)%element%vc = -parts(ipart,rank)%element%vc 
      else if ( d == 3 ) then 
         parts(ipart,rank)%kc = 2* Ng    +1-parts(ipart,rank)%kc 
         parts(ipart,rank)%element%zc = -parts(ipart,rank)%element%zc 
         parts(ipart,rank)%element%wc = -parts(ipart,rank)%element%wc 
      else if ( d == 6 ) then 
         parts(ipart,rank)%kc = 2*(Nz+Ng)+1-parts(ipart,rank)%kc 
         parts(ipart,rank)%element%zc = 2.d0*zLength - parts(ipart,rank)%element%zc 
         parts(ipart,rank)%element%wc = -parts(ipart,rank)%element%wc 
      end if ! d
   end subroutine PartBC_reflect
   
   subroutine LinearIntrpl(x,x0,x1,f0,f1,f)
      implicit none
      real(8), intent (in) :: x,x0,x1,f0,f1
      real(8), intent(out) :: f      
      real(8) :: xl,xr

      xl = (x-x0)/(x1-x0)
      xr = 1.d0 - xl
      f  = f0*xr + f1*xl
   end subroutine LinearIntrpl

   subroutine BilinearIntrpl(x,y,x0,y0,x1,y1,f00,f01,f10,f11,f)
      implicit none
      real(8), intent (in) :: x,y,x0,y0,x1,y1,f00,f01,f10,f11
      real(8), intent(out) :: f      
      real(8) :: f0,f1

      call LinearIntrpl(x,x0,x1,f00,f10,f0)
      call LinearIntrpl(x,x0,x1,f01,f11,f1)
      call LinearIntrpl(y,y0,y1,f0 ,f1 ,f)
   end subroutine BilinearIntrpl

   subroutine TrilinearIntrpl(x,y,z,x0,y0,z0,x1,y1,z1,f000,f001,f010,f011,f100,f101,f110,f111,f)
      implicit none
      real(8), intent (in) :: x,y,z,x0,y0,z0,x1,y1,z1, & 
                              f000,f001,f010,f011,f100,f101,f110,f111
      real(8), intent(out) :: f      
      real(8) :: f0,f1,f00,f01,f10,f11

      call LinearIntrpl(x,x0,x1,f000,f100,f00)
      call LinearIntrpl(x,x0,x1,f010,f110,f10)
      call LinearIntrpl(x,x0,x1,f001,f101,f01)
      call LinearIntrpl(x,x0,x1,f011,f111,f11)
      call LinearIntrpl(y,y0,y1,f00,f10,f0)
      call LinearIntrpl(y,y0,y1,f01,f11,f1)
      call LinearIntrpl(z,z0,z1,f0 ,f1 ,f)
   end subroutine TrilinearIntrpl

   subroutine TrilinearIntrplFluidVel(xp,yp,zp,ip,jp,kp,vel,sdvel,dir)
      implicit none
      real(8), intent (in) :: xp,yp,zp
      integer, intent (in) :: ip,jp,kp,dir
      real(8), intent(out) :: vel,sdvel

      integer :: si,sj,sk

      if ( dir == 1 ) then ! Trilinear interpolation for u
         si = -1;sj = -1;sk = -1 
         if ( yp > y(jp) ) sj =  0 
         if ( zp > z(kp) ) sk =  0
         call TrilinearIntrpl(xp,yp,zp,xh(ip  +si),yh(jp  +sj),zh(kp  +sk),        & 
                                       xh(ip+1+si),yh(jp+1+sj),zh(kp+1+sk),        &
                                       u(ip  +si,jp  +sj,kp  +sk), & 
                                       u(ip  +si,jp  +sj,kp+1+sk), & 
                                       u(ip  +si,jp+1+sj,kp  +sk), & 
                                       u(ip  +si,jp+1+sj,kp+1+sk), & 
                                       u(ip+1+si,jp  +sj,kp  +sk), & 
                                       u(ip+1+si,jp  +sj,kp+1+sk), & 
                                       u(ip+1+si,jp+1+sj,kp  +sk), & 
                                       u(ip+1+si,jp+1+sj,kp+1+sk), vel)

         if ( UnsteadyPartForce ) & 
         call TrilinearIntrpl(xp,yp,zp,xh(ip  +si),yh(jp  +sj),zh(kp  +sk),        & 
                                       xh(ip+1+si),yh(jp+1+sj),zh(kp+1+sk),        &
                                       sdu(ip  +si,jp  +sj,kp  +sk), & 
                                       sdu(ip  +si,jp  +sj,kp+1+sk), & 
                                       sdu(ip  +si,jp+1+sj,kp  +sk), & 
                                       sdu(ip  +si,jp+1+sj,kp+1+sk), & 
                                       sdu(ip+1+si,jp  +sj,kp  +sk), & 
                                       sdu(ip+1+si,jp  +sj,kp+1+sk), & 
                                       sdu(ip+1+si,jp+1+sj,kp  +sk), & 
                                       sdu(ip+1+si,jp+1+sj,kp+1+sk), sdvel)
      else if ( dir == 2 ) then ! Trilinear interpolation for v
         si = -1;sj = -1;sk = -1 
         if ( xp > x(ip) ) si =  0 
         if ( zp > z(kp) ) sk =  0
         call TrilinearIntrpl(xp,yp,zp,xh(ip  +si),yh(jp  +sj),zh(kp  +sk),        & 
                                       xh(ip+1+si),yh(jp+1+sj),zh(kp+1+sk),        &
                                       v(ip  +si,jp  +sj,kp  +sk), & 
                                       v(ip  +si,jp  +sj,kp+1+sk), & 
                                       v(ip  +si,jp+1+sj,kp  +sk), & 
                                       v(ip  +si,jp+1+sj,kp+1+sk), & 
                                       v(ip+1+si,jp  +sj,kp  +sk), & 
                                       v(ip+1+si,jp  +sj,kp+1+sk), & 
                                       v(ip+1+si,jp+1+sj,kp  +sk), & 
                                       v(ip+1+si,jp+1+sj,kp+1+sk), vel)

         if ( UnsteadyPartForce ) & 
         call TrilinearIntrpl(xp,yp,zp,xh(ip  +si),yh(jp  +sj),zh(kp  +sk),        & 
                                       xh(ip+1+si),yh(jp+1+sj),zh(kp+1+sk),        &
                                       sdv(ip  +si,jp  +sj,kp  +sk), & 
                                       sdv(ip  +si,jp  +sj,kp+1+sk), & 
                                       sdv(ip  +si,jp+1+sj,kp  +sk), & 
                                       sdv(ip  +si,jp+1+sj,kp+1+sk), & 
                                       sdv(ip+1+si,jp  +sj,kp  +sk), & 
                                       sdv(ip+1+si,jp  +sj,kp+1+sk), & 
                                       sdv(ip+1+si,jp+1+sj,kp  +sk), & 
                                       sdv(ip+1+si,jp+1+sj,kp+1+sk), sdvel)
      else if ( dir == 3 ) then ! Trilinear interpolation for w
         si = -1;sj = -1;sk = -1 
         if ( xp > x(ip) ) si =  0 
         if ( yp > y(jp) ) sj =  0
         call TrilinearIntrpl(xp,yp,zp,xh(ip  +si),yh(jp  +sj),zh(kp  +sk),        & 
                                       xh(ip+1+si),yh(jp+1+sj),zh(kp+1+sk),        &
                                       w(ip  +si,jp  +sj,kp  +sk), & 
                                       w(ip  +si,jp  +sj,kp+1+sk), & 
                                       w(ip  +si,jp+1+sj,kp  +sk), & 
                                       w(ip  +si,jp+1+sj,kp+1+sk), & 
                                       w(ip+1+si,jp  +sj,kp  +sk), & 
                                       w(ip+1+si,jp  +sj,kp+1+sk), & 
                                       w(ip+1+si,jp+1+sj,kp  +sk), & 
                                       w(ip+1+si,jp+1+sj,kp+1+sk), vel)

         if ( UnsteadyPartForce ) & 
         call TrilinearIntrpl(xp,yp,zp,xh(ip  +si),yh(jp  +sj),zh(kp  +sk),        & 
                                       xh(ip+1+si),yh(jp+1+sj),zh(kp+1+sk),        &
                                       sdw(ip  +si,jp  +sj,kp  +sk), & 
                                       sdw(ip  +si,jp  +sj,kp+1+sk), & 
                                       sdw(ip  +si,jp+1+sj,kp  +sk), & 
                                       sdw(ip  +si,jp+1+sj,kp+1+sk), & 
                                       sdw(ip+1+si,jp  +sj,kp  +sk), & 
                                       sdw(ip+1+si,jp  +sj,kp+1+sk), & 
                                       sdw(ip+1+si,jp+1+sj,kp  +sk), & 
                                       sdw(ip+1+si,jp+1+sj,kp+1+sk), sdvel)
      else 
         call lpperror("Wrong direction in velocity interploation!")
      end if ! dir 
   end subroutine TrilinearIntrplFluidVel

   subroutine ComputeAveFluidVel(xp,yp,zp,ip,jp,kp,um,vm,wm,sdum,sdvm,sdwm,L)

      implicit none 
      real(8), intent(in) :: xp,yp,zp,L
      integer, intent(in) :: ip,jp,kp
      real(8), intent(out):: um,vm,wm,sdum,sdvm,sdwm

      integer :: i,j,k
      integer :: i1,i2,j1,j2,k1,k2
      real(8) :: numcell

      call FindCellIndexBdryConvertRegUnifMesh(L,i1,ip,i2,j1,jp,j2,k1,kp,k2)

      um = 0.d0
      vm = 0.d0 
      wm = 0.d0 
      sdum = 0.d0 
      sdvm = 0.d0 
      sdwm = 0.d0 
      do k=k1,k2; do j=j1,j2; do i=i1,i2
         um = um + u(i,j,k)
         vm = vm + v(i,j,k)
         wm = wm + w(i,j,k)
         sdum = sdum + sdu(i,j,k)
         sdvm = sdvm + sdv(i,j,k)
         sdwm = sdwm + sdw(i,j,k)
      end do; end do; end do
      numcell = dble(i2-i1+1)*dble(j2-j1+1)*dble(k2-k1+1)
      um = um/numcell
      vm = vm/numcell
      wm = wm/numcell
      sdum = sdum/numcell
      sdvm = sdvm/numcell
      sdwm = sdwm/numcell

   end subroutine ComputeAveFluidVel 

   subroutine StoreBeforeConvectionTerms()
      ! Store du, dv, dw before convection terms are computed and added
      implicit none
      if ( UnsteadyPartForce ) then 
         sdu_work = du
         sdv_work = dv
         sdw_work = dw
      end if ! UnsteadyPartForce
   end subroutine StoreBeforeConvectionTerms

   subroutine StoreAfterConvectionTerms()
      ! Store the portions of du, dv, dw due to convection terms
      implicit none
      if ( UnsteadyPartForce ) then 
         sdu_work = du - sdu_work
         sdv_work = dv - sdv_work
         sdw_work = dw - sdw_work
      end if ! UnsteadyPartForce
   end subroutine StoreAfterConvectionTerms

   subroutine ComputeSubDerivativeVel()
      implicit none

      integer :: i, j, k
     
      ! Note: XXX The current way to compute substantial derivatives of velocity
      ! only works when diffusion terms are solved EXPLICITLY

      ! Subtract the portions of du, dv, dw due to convection terms from the
      ! overall values
      if ( UnsteadyPartForce ) then 
         sdu = du - sdu_work
         sdv = dv - sdv_work
         sdw = dw - sdw_work

         ! Correct du, dv, dw with pressure gradients
         do k=ks,ke;  do j=js,je; do i=is,ieu     
            sdu(i,j,k)=sdu(i,j,k)-(2.0/dxh(i))*(p(i+1,j,k)-p(i,j,k))/(rho(i+1,j,k)+rho(i,j,k))
         enddo; enddo; enddo

         do k=ks,ke;  do j=js,jev; do i=is,ie    
            sdv(i,j,k)=sdv(i,j,k)-(2.0/dyh(j))*(p(i,j+1,k)-p(i,j,k))/(rho(i,j+1,k)+rho(i,j,k))
         enddo; enddo; enddo

         do k=ks,kew;  do j=js,je; do i=is,ie   
            sdw(i,j,k)=sdw(i,j,k)-(2.0/dzh(k))*(p(i,j,k+1)-p(i,j,k))/(rho(i,j,k+1)+rho(i,j,k))
         enddo; enddo; enddo
      end if ! UnsteadyPartForce

   end subroutine ComputeSubDerivativeVel

   subroutine SeedParticles()

      implicit none
      include 'mpif.h'
#ifdef __INTEL_COMPILER
!      use IFPORT
#endif

      logical :: SeedParticlesNow
      integer :: ipart_seed,itime,itime_part_seed
      real(8) :: xp,yp,zp,up,vp,wp,dp,volp
      real(8) :: rand
      integer :: ip,jp,kp
!      integer :: c1,c2,c3,rankpart
      integer :: ierr
      integer :: num_part_rank
      real(8) :: dist,dp1
      logical :: PartSeedTooClose
      integer :: ipart,ipart1,iter,niters
      real(8) :: xb1,xb2,yb1,yb2,zb1,zb2
      real(8) :: Mxt,Myt,Mzt
      integer :: num_part_seeded(0:nPdomain-1),rankseed,c1,c2,c3
      real(8) :: dummyreal,uf,vf,wf
      integer,parameter :: largeint = 1234567890

      SeedParticlesNow = .false.
      if ( SeedParticlesFlag == SeedParticlesNone ) then 
         return 
      else if ( SeedParticlesFlag == SeedParticlesOneTime ) then
         itime           = nint(time/dt)
         itime_part_seed = nint(time_part_seed/dt)
         if ( itime == itime_part_seed ) then 
            SeedParticlesNow = .true.
         end if ! t
      else if ( SeedParticlesFlag == SeedParticlesRegular ) then
         itime           = nint(time/dt)
         itime_part_seed = nint(time_part_seed/dt)
         if ( mod(itime,itime_part_seed) == 0 ) then 
            SeedParticlesNow = .true.
         end if ! itime_part_seed
      else 
         call lpperror("Unknown seeding particle flag!")
      end if ! SeedParticlesFlag

      num_part_seeded = 0 
      if ( SeedParticlesNow ) then
         call random_seed()
         if ( rank == 0 ) write(*,*) num_part_seed,'particles to seed at t=',time
         ipart_seed = 0
         do 
            ipart_seed = ipart_seed + 1
            call generate_random_particle(dmin_part_seed,dmax_part_seed,   & 
               xmin_part_seed,xmax_part_seed,ymin_part_seed,ymax_part_seed,&
               zmin_part_seed,zmax_part_seed,umin_part_seed,umax_part_seed,&
               vmin_part_seed,vmax_part_seed,wmin_part_seed,wmax_part_seed,& 
               dp,volp,xp,yp,zp,up,vp,wp)

            if ( test_injectdrop ) then
               Mxt = NINT((xmax_part_seed - xmin_part_seed)/dx(is))/npx
               Myt = NINT((ymax_part_seed - ymin_part_seed)/dy(js))/npy
               Mzt = NINT((zmax_part_seed - zmin_part_seed)/dz(ks))/npz
               c1 = INT((xp-xmin_part_seed)/dx(is))/Mxt  
               c2 = INT((yp-ymin_part_seed)/dy(js))/Myt 
               c3 = INT((zp-zmin_part_seed)/dz(ks))/Mzt 
               rankseed = c1*npy*npz + c2*npz + c3
               niters = 200
               do iter = 1,niters
                  if ( num_part_seeded(rankseed) >= NINT(dble(num_part_seed)/dble(npx*npy*npz)*1.1d0) ) then   
                     call generate_random_particle(dmin_part_seed,dmax_part_seed,   & 
                        xmin_part_seed,xmax_part_seed,ymin_part_seed,ymax_part_seed,&
                        zmin_part_seed,zmax_part_seed,umin_part_seed,umax_part_seed,&
                        vmin_part_seed,vmax_part_seed,wmin_part_seed,wmax_part_seed,& 
                        dp,volp,xp,yp,zp,up,vp,wp)
                     c1 = INT((xp-xmin_part_seed)/dx(is))/Mxt  
                     c2 = INT((yp-ymin_part_seed)/dy(js))/Myt 
                     c3 = INT((zp-zmin_part_seed)/dz(ks))/Mzt 
                     rankseed = c1*npy*npz + c2*npz + c3
                  else 
                     exit
                  end if ! num_part_seeded(rankseed)
               end do ! iter
               num_part_seeded(rankseed) = num_part_seeded(rankseed) + 1
            end if ! test_injectdrop

            if ( xp > xh(is-1) .and. xp <= xh(ie) .and. & 
                 yp > yh(js-1) .and. yp <= yh(je) .and. & 
                 zp > zh(ks-1) .and. zp <= zh(ke) ) then   
               num_part(rank) = num_part(rank) + 1
               parts(num_part(rank),rank)%element%xc = xp
               parts(num_part(rank),rank)%element%yc = yp
               parts(num_part(rank),rank)%element%zc = zp
               parts(num_part(rank),rank)%element%vol = volp
               call FindPartLocCell(xp,yp,zp,ip,jp,kp)
               parts(num_part(rank),rank)%ic = ip
               parts(num_part(rank),rank)%jc = jp
               parts(num_part(rank),rank)%kc = kp
               if ( SeedWithFluidVel ) then 
                  call GetFluidProp(ip,jp,kp,xp,yp,zp, & 
                              uf,vf,wf,dummyreal,dummyreal,dummyreal,& 
                              volp,largeint,TwoWayCouplingFlag)
                  up = uf; vp = vf; wp = wf
               end if ! SeedWithFluidVel
               parts(num_part(rank),rank)%element%uc = up
               parts(num_part(rank),rank)%element%vc = vp
               parts(num_part(rank),rank)%element%wc = wp
            end if ! xp,yp,zp 
            if ( ipart_seed >= num_part_seed ) exit
         end do ! ipart_seed
      
         ! Adjust particle positions to avoid locating at block bdry 
         if ( AvoidSeedAtBlockBdry .and. num_part(rank) > 0 ) then
            do ipart = 1,num_part(rank)
               dp = (parts(ipart,rank)%element%vol*6.d0/PI)**OneThird
               xp = parts(ipart,rank)%element%xc
               yp = parts(ipart,rank)%element%yc
               zp = parts(ipart,rank)%element%zc
               if ( xp - 0.65d0*dp < xh(is-1) ) xp = xh(is-1) + 0.65d0*dp
               if ( xp + 0.65d0*dp > xh(ie)   ) xp = xh(ie)   - 0.65d0*dp
               if ( yp - 0.65d0*dp < yh(js-1) ) yp = yh(js-1) + 0.65d0*dp
               if ( yp + 0.65d0*dp > yh(je)   ) yp = yh(je)   - 0.65d0*dp
               if ( zp - 0.65d0*dp < zh(ks-1) ) zp = zh(ks-1) + 0.65d0*dp
               if ( zp + 0.65d0*dp > zh(ke)   ) zp = zh(ke)   - 0.65d0*dp
               parts(ipart,rank)%element%xc = xp 
               parts(ipart,rank)%element%yc = yp
               parts(ipart,rank)%element%zc = zp
            end do !ipart
         end if !AvoidSeedAtBlockBdry
           
         ! Adjust particle positions to avoid particle too close to each other 
         if ( AvoidSeedPartTooClose .and. num_part(rank) > 1 ) then
            niters = 5000 
            do iter = 1,niters
               do ipart = 1,num_part(rank)
                  ! Check if a particle is too close with others
                  PartSeedTooClose = .false.
                  dp = (parts(ipart,rank)%element%vol*6.d0/PI)**OneThird
                  do ipart1 = 1,num_part(rank)
                     if ( ipart1 /= ipart ) then 
                        dist = (parts(ipart ,rank)%element%xc - & 
                                parts(ipart1,rank)%element%xc)**2.d0 & 
                             + (parts(ipart ,rank)%element%yc - &
                                parts(ipart1,rank)%element%yc)**2.d0 &
                             + (parts(ipart ,rank)%element%zc - &
                                parts(ipart1,rank)%element%zc)**2.d0
                        dp1  = (parts(ipart1,rank)%element%vol*6.d0/PI)**OneThird
                        if ( sqrt(dist) < (dp+dp1)*1.3d0 ) then 
                           PartSeedTooClose = .true.
                           exit
                        end if ! sqrt(dist)
                     end if ! ipart1 
                  end do ! ipart1
                  ! relocate the particle within the block
                  if ( PartSeedTooClose ) then
                     xb1 = max(xh(is-1)+0.65d0*dp,xmin_part_seed)
                     xb2 = min(xh(ie  )-0.65d0*dp,xmax_part_seed)
                     call random_number(rand) 
                     parts(ipart,rank)%element%xc = xb1 + rand*(xb2-xb1)
                                 
                     yb1 = max(yh(js-1)+0.65d0*dp,ymin_part_seed) 
                     yb2 = min(yh(je  )-0.65d0*dp,ymax_part_seed)
                     call random_number(rand) 
                     parts(ipart,rank)%element%yc = yb1 + rand*(yb2-yb1)

                     zb1 = max(zh(ks-1)+0.65d0*dp,zmin_part_seed) 
                     zb2 = min(zh(ke  )-0.65d0*dp,zmax_part_seed) 
                     call random_number(rand) 
                     parts(ipart,rank)%element%zc = zb1 + rand*(zb2-zb1)
                     exit
                  end if ! PartSeedTooClose
               end do ! ipart
               if ( .NOT.PartSeedTooClose ) exit
            end do ! iter 
            write(*,*) 'iteration on adjusting particle locations',rank,iter,num_part(rank)
         end if !AvoidSeedPartTooClose

         ! recompute particle velocity if seed with fluid vel and part positions
         ! have been adjusted
         if ( SeedWithFluidVel .and. (AvoidSeedAtBlockBdry .or. AvoidSeedPartTooClose) ) then
            do ipart = 1, num_part(rank)
               xp = parts(ipart,rank)%element%xc
               yp = parts(ipart,rank)%element%yc
               zp = parts(ipart,rank)%element%zc
               call FindPartLocCell(xp,yp,zp,ip,jp,kp)
               parts(ipart,rank)%ic = ip
               parts(ipart,rank)%jc = jp
               parts(ipart,rank)%kc = kp
               call GetFluidProp(ip,jp,kp,xp,yp,zp, &
                           uf,vf,wf,dummyreal,dummyreal,dummyreal,& 
                           volp,largeint,TwoWayCouplingFlag)
               up = uf; vp = vf; wp = wf
               parts(ipart,rank)%element%uc = up
               parts(ipart,rank)%element%vc = vp
               parts(ipart,rank)%element%wc = wp
            end do ! ipart 
         end if ! SeedWithFluidVel

         ! Allgather updated number of partilces after seeding
         num_part_rank = num_part(rank)
         call MPI_ALLGATHER(num_part_rank, 1, MPI_INTEGER, &
                            num_part(:)   , 1, MPI_INTEGER, MPI_Comm_World, ierr)
         if ( rank == 0 ) write(*,*) sum(num_part), "particles seeded!"
      end if ! SeedParticlesNow
   
      contains
         subroutine generate_random_particle(dmin_part_seed,dmax_part_seed,   & 
               xmin_part_seed,xmax_part_seed,ymin_part_seed,ymax_part_seed,&
               zmin_part_seed,zmax_part_seed,umin_part_seed,umax_part_seed,&
               vmin_part_seed,vmax_part_seed,wmin_part_seed,wmax_part_seed,& 
               dp,volp,xp,yp,zp,up,vp,wp)

            implicit none

            real(8),intent(in) :: dmin_part_seed,dmax_part_seed,   &
                                 xmin_part_seed,xmax_part_seed,ymin_part_seed,ymax_part_seed,&
                                 zmin_part_seed,zmax_part_seed,umin_part_seed,umax_part_seed,&
                                 vmin_part_seed,vmax_part_seed,wmin_part_seed,wmax_part_seed
            real(8),intent(out) :: dp,volp,xp,yp,zp,up,vp,wp

            call random_number(rand) 
            dp = dmin_part_seed + rand*(dmax_part_seed-dmin_part_seed)
            volp = PI*dp*dp*dp/6.d0
            call random_number(rand) 
            xp = xmin_part_seed + rand*(xmax_part_seed-xmin_part_seed)
            call random_number(rand) 
            yp = ymin_part_seed + rand*(ymax_part_seed-ymin_part_seed)
            call random_number(rand) 
            zp = zmin_part_seed + rand*(zmax_part_seed-zmin_part_seed)
            call random_number(rand) 
            up = umin_part_seed + rand*(umax_part_seed-umin_part_seed)
            call random_number(rand) 
            vp = vmin_part_seed + rand*(vmax_part_seed-vmin_part_seed)
            call random_number(rand) 
            wp = wmin_part_seed + rand*(wmax_part_seed-wmin_part_seed)
         end subroutine generate_random_particle
   end subroutine seedparticles

   subroutine lpperror(message) 
      use module_IO
      implicit none
      include 'mpif.h'
      integer ierr, MPI_errorcode
      character(*) :: message
      MPI_errorcode=1
      write(*,*) "LPP ERROR *** ",message, " *** STOP at rank: ", rank
      ! Exit MPI gracefully
      close(out)
      call MPI_ABORT(MPI_COMM_WORLD, MPI_errorcode, ierr)
      call MPI_finalize(ierr)
      stop 
   end subroutine lpperror

end module module_Lag_part

! ==================================================================================================
! module_output_lpp: I/O module for Lagrangian particle 
! ==================================================================================================
module module_output_lpp
   use module_IO
   use module_flow
   use module_Lag_part
   implicit none
   integer :: lpp_opened=0

   contains

   subroutine append_LPP_visit_file(rootname)
      implicit none
      character(*) :: rootname
      integer prank
      integer, parameter :: LPPformatPlot3D = 1
      integer, parameter :: LPPformatVOFVTK = 2
      integer, parameter :: LPPformatVTK = 3
      if(rank.ne.0) call pariserror('rank.ne.0 in append_LPP')
      if(lpp_opened==0) then
         OPEN(UNIT=88,FILE='lpp.visit')
         write(88,10) NpDomain
10       format('!NBLOCKS ',I4)
         lpp_opened=1
      else
         OPEN(UNIT=88,FILE='lpp.visit',position='append')
      endif
      do prank=0,NpDomain-1
         if ( outputlpp_format == LPPformatPlot3D ) then 
            write(88,11) rootname//TRIM(int2text(prank,padding))//'.3D'
         else if ( outputlpp_format == LPPformatVOFVTK .or. & 
                   outputlpp_format == LPPformatVTK ) then
            write(88,11) rootname//TRIM(int2text(prank,padding))//'.vtk'
         else 
            call pariserror("Unknow LPP output format!")
         end if ! outputlpp_format
11       format(A)
      enddo
      close(88)
   end subroutine  append_LPP_visit_file

   subroutine output_LPP(nf,i1,i2,j1,j2,k1,k2)

      implicit none
      integer,intent(in)  :: nf,i1,i2,j1,j2,k1,k2
      integer, parameter :: LPPformatPlot3D = 1
      integer, parameter :: LPPformatVOFVTK = 2
      integer, parameter :: LPPformatVTK = 3

      if ( maxval(num_part) < 1 .and. .not.DoConvertVOF2LPP ) return 

      if ( outputlpp_format == LPPformatPlot3D ) then 
         call output_LPP_Plot3D(nf,i1,i2,j1,j2,k1,k2)
      else if ( outputlpp_format == LPPformatVOFVTK ) then
         call output_LPP_VOFVTK(nf,i1,i2,j1,j2,k1,k2)
      else if ( outputlpp_format == LPPformatVTK ) then
         call output_LPP_VTK(nf,i1,i2,j1,j2,k1,k2)
      else 
         call pariserror("Unknow LPP output format!")
      end if ! outputlpp_format

   end subroutine output_LPP

   subroutine output_LPP_Plot3D(nf,i1,i2,j1,j2,k1,k2)
      implicit none
      integer,intent(in)  :: nf,i1,i2,j1,j2,k1,k2
      character(len=30) :: rootname
      integer :: ipart
      real(8), parameter :: volsmall=1.0d-60

      ! output lpp data in plot3d formate 
      rootname=trim(out_path)//'/VTK/LPP'//TRIM(int2text(nf,padding))//'-'
      if(rank==0) call append_LPP_visit_file(TRIM(rootname))

      OPEN(UNIT=8,FILE=TRIM(rootname)//TRIM(int2text(rank,padding))//'.3D')
!      write(8,10)
      write(8,11)
!10    format('# plot3D data file')
11    format('x y z vol')

      if ( num_part(rank) > 0 ) then 
         do ipart = 1,num_part(rank) 
            write(8,320) parts(ipart,rank)%element%xc,& 
            parts(ipart,rank)%element%yc, & 
            parts(ipart,rank)%element%zc, &  
            parts(ipart,rank)%element%vol
         enddo
      end if ! num_part(rank)
      ! add virtual particles at the 8 corners of the block, so that
      ! moledular plot in VisIt would work in the same scale of vof
      write(8,320) x(i1),y(j1),z(k1),volsmall 
      write(8,320) x(i1),y(j1),z(k2),volsmall 
      write(8,320) x(i1),y(j2),z(k1),volsmall 
      write(8,320) x(i1),y(j2),z(k2),volsmall 
      write(8,320) x(i2),y(j1),z(k1),volsmall 
      write(8,320) x(i2),y(j1),z(k2),volsmall 
      write(8,320) x(i2),y(j2),z(k1),volsmall 
      write(8,320) x(i2),y(j2),z(k2),volsmall 
320   format(e14.5,e14.5,e14.5,e14.5)
      close(8)

      ! write complete lpp data for other post-processing treatment
      rootname=trim(out_path)//'/VTK/LPPDATA'//TRIM(int2text(nf,padding))//'-'
      OPEN(UNIT=8,FILE=TRIM(rootname)//TRIM(int2text(rank,padding))//'.dat')
      write(8,*) time,num_part(rank) 
      if ( num_part(rank) > 0 ) then 
         do ipart = 1,num_part(rank) 
            write(8,'(10(E15.8,1X))') & 
            parts(ipart,rank)%element%xc,& 
            parts(ipart,rank)%element%yc, & 
            parts(ipart,rank)%element%zc, &  
            parts(ipart,rank)%element%uc, &  
            parts(ipart,rank)%element%vc, &  
            parts(ipart,rank)%element%wc, &  
            parts(ipart,rank)%uf, &  
            parts(ipart,rank)%vf, &  
            parts(ipart,rank)%wf, &  
            parts(ipart,rank)%element%vol
         enddo
      end if ! num_part(rank)
      close(8)

   end subroutine output_LPP_Plot3D

   subroutine output_LPP_VOFVTK(nf,il,ir,jl,jr,kl,kr)
      implicit none
      integer,intent(in)  :: nf,il,ir,jl,jr,kl,kr
      character(len=30) :: rootname,filename
      integer :: ipart
      real(8) :: lppvof(imin:imax,jmin:jmax,kmin:kmax)
      integer :: i1,ic,i2,j1,jc,j2,k1,kc,k2,i,j,k,i0,j0,k0
      real(8) :: dp,rp,ConvertRegSize,xp,yp,zp,c,stencil3x3(-1:1,-1:1,-1:1)
      integer :: nflag
      real(8) :: volcell


      lppvof = 0.d0

      if ( num_part(rank) > 0 ) then 
         do ipart = 1,num_part(rank)
            ic = parts(ipart,rank)%ic
            jc = parts(ipart,rank)%jc
            kc = parts(ipart,rank)%kc
            volcell = dx(ic)*dy(jc)*dz(kc)
            if ( parts(ipart,rank)%element%vol < volcell ) cycle  ! Note: not plot drops to small
            dp = (6.d0*parts(ipart,rank)%element%vol/PI)**(1.d0/3.d0)
            rp = 0.5d0*dp
            ConvertRegSize = ConvertRegSizeToDiam*dp
            xp = parts(ipart,rank)%element%xc
            yp = parts(ipart,rank)%element%yc
            zp = parts(ipart,rank)%element%zc
            call FindCellIndexBdryConvertRegUnifMesh(ConvertRegSize, & 
                                             i1,ic,i2,j1,jc,j2,k1,kc,k2)
         
            do i=i1+1,i2-1; do j=j1+1,j2-1; do k=k1+1,k2-1
               ! Build VOF field of droplet
               do i0=-1,1; do j0=-1,1; do k0=-1,1
                  stencil3x3(i0,j0,k0) = rp**2.d0 & 
                     - ((x(i+i0)-xp)**2.d0+(y(j+j0)-yp)**2.d0+(z(k+k0)-zp)**2.d0)
               enddo; enddo; enddo
               call ls2vof_in_cell(stencil3x3,c,nflag)
               lppvof(i,j,k) = c
               if ( lppvof(i,j,k) > 1.d0 ) lppvof(i,j,k) = 1.d0
               if ( lppvof(i,j,k) < 0.d0 ) lppvof(i,j,k) = 0.d0
            end do; end do; end do
         end do ! ipart
      end if ! num_part(rank)

      rootname=trim(out_path)//'/VTK/LPPVOF'//TRIM(int2text(nf,padding))//'-'
      if(rank==0) call append_LPP_visit_file(TRIM(rootname))

      OPEN(UNIT=8,FILE=TRIM(rootname)//TRIM(int2text(rank,padding))//'.vtk')
      write(8,10)
      write(8,11) time
      write(8,12)
      write(8,13)
      write(8,14)ir-il+1,jr-jl+1,kr-kl+1
      write(8,15)(ir-il+1)*(jr-jl+1)*(kr-kl+1)
10    format('# vtk DataFile Version 2.0')
11    format('grid, time ',F16.8)
12    format('ASCII')
13    format('DATASET STRUCTURED_GRID')
14    format('DIMENSIONS ',I5,I5,I5)
15    format('POINTS ',I17,' double' )

      do k=kl,kr; do j=jl,jr; do i=il,ir;
         write(8,320) x(i),y(j),z(k)
      enddo; enddo; enddo
320   format(e14.5,e14.5,e14.5)

      write(8,16)(ir-il+1)*(jr-jl+1)*(kr-kl+1)
      write(8,17)'LPPVOF'
      write(8,18)
16    format('POINT_DATA ',I17)
17    format('SCALARS ',A20,' float 1')
18    format('LOOKUP_TABLE default')

      do k=kl,kr; do j=jl,jr; do i=il,ir;
         write(8,210) lppvof(i,j,k)
      enddo; enddo; enddo
210   format(e14.5)
      close(8)
         
      if ( zip_data ) then 
         filename = TRIM(rootname)//TRIM(int2text(rank,padding))//'.vtk'
         call system('gzip '//trim(filename))
      end if ! zip_data
      
      ! write complete lpp data for other post-processing treatment
      rootname=trim(out_path)//'/VTK/LPPDATA'//TRIM(int2text(nf,padding))//'-'
      OPEN(UNIT=8,FILE=TRIM(rootname)//TRIM(int2text(rank,padding))//'.dat')
      write(8,*) time,num_part(rank) 
      if ( num_part(rank) > 0 ) then 
         do ipart = 1,num_part(rank) 
            write(8,'(10(E15.8,1X))') & 
            parts(ipart,rank)%element%xc,& 
            parts(ipart,rank)%element%yc, & 
            parts(ipart,rank)%element%zc, &  
            parts(ipart,rank)%element%uc, &  
            parts(ipart,rank)%element%vc, &  
            parts(ipart,rank)%element%wc, &  
            parts(ipart,rank)%uf, &  
            parts(ipart,rank)%vf, &  
            parts(ipart,rank)%wf, &  
            parts(ipart,rank)%element%vol
         enddo
      end if ! num_part(rank)
      close(8)

   end subroutine output_LPP_VOFVTK

   subroutine output_LPP_VTK(nf,i1,i2,j1,j2,k1,k2)
      implicit none
      integer,intent(in)  :: nf,i1,i2,j1,j2,k1,k2
      character(len=30) :: rootname
      integer :: ipart
      real(8), parameter :: smallradius=1.0d-60

      ! output lpp data in plot3d formate 
      rootname=trim(out_path)//'/VTK/LPP'//TRIM(int2text(nf,padding))//'-'
      if(rank==0) call append_LPP_visit_file(TRIM(rootname))

      OPEN(UNIT=8,FILE=TRIM(rootname)//TRIM(int2text(rank,padding))//'.vtk')
      write(8,10)
      write(8,11)
      write(8,12)
      write(8,13)
      write(8,14)
      write(8,15) num_part(rank)+2
10    format('# vtk DataFile Version 3.0')
11    format('vtk output')
12    format('ASCII')
13    format('DATASET POLYDATA')
14    format('')
15    format('POINTS ',I17,' float' )

      do ipart = 1,num_part(rank) 
         write(8,320) parts(ipart,rank)%element%xc,& 
         parts(ipart,rank)%element%yc, & 
         parts(ipart,rank)%element%zc  
      enddo
      ! add virtual particles at the 8 corners of the block, so that
      ! moledular plot in VisIt would work in the same scale of vof
      write(8,320) x(i1),y(j1),z(k1) 
      write(8,320) x(i2),y(j2),z(k2) 
320   format(e14.5,e14.5,e14.5)

      write(8,16) num_part(rank)+2
      write(8,17)'radius'
      write(8,18)
16    format('POINT_DATA ',I17)
17    format('SCALARS ',A20,' float')
18    format('LOOKUP_TABLE default')

      do ipart = 1,num_part(rank) 
         write(8,330) (parts(ipart,rank)%element%vol*0.75d0/PI)**0.3333333333333333d0 
      enddo
      write(8,330) smallradius 
      write(8,330) smallradius 
330   format(e14.5)

      write(8,17)'up'
      write(8,18)
      do ipart = 1,num_part(rank) 
         write(8,330) parts(ipart,rank)%element%uc
      enddo
      write(8,330) smallradius 
      write(8,330) smallradius

      write(8,17)'vp'
      write(8,18)
      do ipart = 1,num_part(rank) 
         write(8,330) parts(ipart,rank)%element%vc
      enddo
      write(8,330) smallradius 
      write(8,330) smallradius

      write(8,17)'wp'
      write(8,18)
      do ipart = 1,num_part(rank) 
         write(8,330) parts(ipart,rank)%element%wc
      enddo
      write(8,330) smallradius 
      write(8,330) smallradius

      write(8,17)'uf'
      write(8,18)
      do ipart = 1,num_part(rank) 
         write(8,330) parts(ipart,rank)%uf
      enddo
      write(8,330) smallradius 
      write(8,330) smallradius

      write(8,17)'vf'
      write(8,18)
      do ipart = 1,num_part(rank) 
         write(8,330) parts(ipart,rank)%vf
      enddo
      write(8,330) smallradius 
      write(8,330) smallradius

      write(8,17)'wf'
      write(8,18)
      do ipart = 1,num_part(rank) 
         write(8,330) parts(ipart,rank)%wf
      enddo
      write(8,330) smallradius 
      write(8,330) smallradius

!19    format('VECTORS ',A20,' float')
!      write(8,19)'velf'
!      write(8,18)
!      do ipart = 1,num_part(rank) 
!         write(8,320) parts(ipart,rank)%uf, &
!            parts(ipart,rank)%vf, &
!            parts(ipart,rank)%wf
!      enddo
!      write(8,320) smallradius,smallradius,smallradius 
!      write(8,320) smallradius,smallradius,smallradius

      close(8)
   end subroutine output_LPP_VTK

!-------------------------------------------------------------------------------------------------
   subroutine backup_LPP_write
      implicit none
      integer ::ipart
      character(len=100) :: filename
      filename = trim(out_path)//'/backuplpp_'//int2text(rank,padding)
      call system('mv '//trim(filename)//' '//trim(filename)//'.old')
      OPEN(UNIT=7,FILE=trim(filename),status='unknown',action='write')
      write(7,1100)time,itimestep,num_part(rank)
      if ( num_part(rank) > 0 ) then 
         do ipart=1,num_part(rank)
            write(7,1200) parts(ipart,rank)%element%xc, & 
                          parts(ipart,rank)%element%yc, & 
                          parts(ipart,rank)%element%zc, & 
                          parts(ipart,rank)%element%uc, & 
                          parts(ipart,rank)%element%vc, & 
                          parts(ipart,rank)%element%wc, & 
                          parts(ipart,rank)%element%duc, & 
                          parts(ipart,rank)%element%dvc, & 
                          parts(ipart,rank)%element%dwc, & 
                          parts(ipart,rank)%element%vol, &  
                          parts(ipart,rank)%ic, & 
                          parts(ipart,rank)%jc, & 
                          parts(ipart,rank)%kc 
         end do! ipart
      end if ! num_part(rank)
      if(rank==0)print*,'Backup LPP written at t=',time
      1100 FORMAT(es17.8e3,2I10)
      1200 FORMAT(10es17.8e3,3I5)
      CLOSE(7)
   end subroutine backup_LPP_write

!-------------------------------------------------------------------------------------------------
   subroutine backup_LPP_read
      implicit none
      integer ::ipart,ierr
      OPEN(UNIT=7,FILE=trim(out_path)//'/backuplpp_'//int2text(rank,padding),status='old',action='read')
      read(7,*)time,itimestep,num_part(rank)
      if ( num_part(rank) < 0 ) &
         call pariserror("Error: backuplpp_read")
      if ( num_part(rank) > 0 ) then 
         do ipart=1,num_part(rank)
            read(7,*    ) parts(ipart,rank)%element%xc, & 
                          parts(ipart,rank)%element%yc, & 
                          parts(ipart,rank)%element%zc, & 
                          parts(ipart,rank)%element%uc, & 
                          parts(ipart,rank)%element%vc, & 
                          parts(ipart,rank)%element%wc, & 
                          parts(ipart,rank)%element%duc, & 
                          parts(ipart,rank)%element%dvc, & 
                          parts(ipart,rank)%element%dwc, & 
                          parts(ipart,rank)%element%vol, &  
                          parts(ipart,rank)%ic, & 
                          parts(ipart,rank)%jc, & 
                          parts(ipart,rank)%kc 
         end do !ipart
      end if ! num_part(rank)
      CLOSE(7)
   end subroutine backup_LPP_read
!-------------------------------------------------------------------------------------------------

end module module_output_LPP
!=================================================================================================

