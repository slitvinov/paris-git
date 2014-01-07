!=================================================================================================
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
   integer,parameter :: maxnum_diff_tag  = 7   ! ignore cases droplet spread over more than 1 block
   integer :: total_num_tag,totalnum_drop,totalnum_drop_indep,num_new_drop
   integer, dimension(:), allocatable :: num_drop
   integer, dimension(:), allocatable :: num_drop_merge
   integer, dimension(:), allocatable :: num_element
   integer, dimension(:), allocatable :: num_tag, tagmin, tagmax
   integer, dimension(:), allocatable :: tag_dropid
   integer, dimension(:), allocatable :: tag_rank
   integer, dimension(:), allocatable :: tag_mergeflag

   integer, dimension(:), allocatable :: new_drop_id

   type element 
      real(8) :: xc,yc,zc,uc,vc,wc,vol
      integer :: id 
   end type element
   type (element), dimension(:,:), allocatable :: element_stat

   type drop
      type(element) :: element
      integer :: num_cell_drop
   end type drop
   type (drop), dimension(:,:), allocatable :: drops
   integer, dimension(:,:,:), allocatable :: drops_cell_list
   
   type drop_merge
      type(element) :: element
      integer :: num_cell_drop
      integer :: num_gcell
      integer :: num_diff_tag
      integer :: diff_tag_list(maxnum_diff_tag)
      integer :: flag_center_mass
   end type drop_merge
   type (drop_merge), dimension(:,:), allocatable :: drops_merge
   integer, dimension(:,:,:), allocatable :: drops_merge_cell_list
   integer, dimension(:,:,:), allocatable :: drops_merge_gcell_list

   type drop_merge_comm
      real(8) :: xc,yc,zc,uc,vc,wc,vol
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
      real(8) :: fx,fy,fz
      real(8) :: xcOld,ycOld,zcOld,ucOld,vcOld,wcOld
      integer :: ic,jc,kc,dummyint  
      ! Note: open_mpi sometimes failed to communicate the last varialbe in 
      !       MPI_TYPE_STRUCt correctly, dummyint is included to go around 
      !       this bug in mpi
   end type particle 
   type (particle), dimension(:,:), allocatable :: parts

   integer, dimension(:,:), allocatable :: parts_cross_id
   integer, dimension(:,:), allocatable :: parts_cross_newrank
   integer, dimension(:), allocatable :: num_part_cross

   ! substantial derivative of fluid velocity, used for unsteady force calulation
   real(8), dimension(:,:,:), allocatable :: sdu,sdv,sdw

   real(8), parameter :: PI = 3.14159265359d0
   integer, parameter :: CRAZY_INT = 3483129 

   integer, parameter :: CriteriaRectangle = 1
   integer, parameter :: CriteriaCylinder  = 2
   integer, parameter :: CriteriaSphere    = 3

   logical :: DoDropStatistics 
   logical :: DoConvertVOF2LPP 
   logical :: DoConvertLPP2VOF 
   integer :: dragmodel 
   integer :: ntimesteptag
   integer :: CriteriaConvertCase
   real(8) :: vol_cut, xlpp_min,ylpp_min,zlpp_min, & 
                       xlpp_max,ylpp_max,zlpp_max  

   integer :: max_num_drop 
   integer :: maxnum_cell_drop
   integer :: max_num_part
   integer :: max_num_part_cross
contains
!=================================================================================================
   subroutine initialize_LPP()

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
      allocate( drops(max_num_drop,0:nPdomain-1) )
      allocate( drops_merge(max_num_drop,0:nPdomain-1) )
      allocate( drops_cell_list(3,maxnum_cell_drop,max_num_drop) )
      allocate( drops_merge_cell_list(3,maxnum_cell_drop,max_num_drop) )
      allocate( drops_merge_gcell_list(3,maxnum_cell_drop,max_num_drop) )
      allocate( sdu(imin:imax,jmin:jmax,kmin:kmax), & 
                sdv(imin:imax,jmin:jmax,kmin:kmax), &
                sdw(imin:imax,jmin:jmax,kmin:kmax) )

      ! set default values
      num_tag  = 0
      num_part = 0
      num_drop = 0
      num_drop_merge = 0
      num_element = 0

      LPP_initialized = .true.

      sdu = 0.d0; sdv = 0.d0; sdw =0.d0

   end subroutine initialize_LPP

   subroutine ReadLPPParameters

      use module_flow
      use module_BC
      implicit none
      include 'mpif.h'

      integer ierr,in
      logical file_is_there
      namelist /lppparameters/ DoDropStatistics, dragmodel, nTimeStepTag,  &
         DoConvertVOF2LPP,DoConvertLPP2VOF,CriteriaConvertCase,            & 
         vol_cut,xlpp_min,xlpp_max,ylpp_min,ylpp_max,zlpp_min,zlpp_max,    &
         max_num_drop, maxnum_cell_drop, max_num_part, max_num_part_cross

      in=32

      ! Set default values 
      DoDropStatistics = .false. 
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
      max_num_drop = 10
      maxnum_cell_drop = 10000
      max_num_part = 100
      max_num_part_cross = 10

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
         if (rank == 0) STOP "ReadLPPParameters: no 'inputlpp' file."
      endif
      close(in)
   end subroutine ReadLPPParameters


   subroutine lppsweeps(tswap)
    
      integer, intent(in) :: tswap
      ! Only do tagging and conversion in specific time steps
      if ( MOD(tswap,ntimestepTag) == 0 ) then 
         call tag_drop()
         if ( nPdomain > 1 ) then 
            call tag_drop_all
            call merge_drop_pieces 
         end if ! nPdomain

         if ( DoDropStatistics ) call drop_statistics(tswap) 

         if ( DoConvertVOF2LPP ) call ConvertDrop2Part(tswap)
!         if ( DoConvertLPP2VOF ) call ConvertPart2Drop()   ! TBA
      end if ! tswap
      call ComputeFluidAccel()
      call ComputePartForce(tswap)
      call UpdatePartSol(tswap)
   end subroutine lppsweeps

  subroutine tag_drop()
    include 'mpif.h'
    integer :: i,j,k, i0,j0,k0
    integer :: isq,jsq,ksq
    integer :: current_id
    integer :: s_queue(Nx*Ny,3),c_queue(Nx*Ny,3) ! record i,j,k
    integer :: ns_queue,is_queue,nc_queue
    
    real(8) :: volcell,cvof_scaled

    integer :: num_cell_drop,cell_list(3,maxnum_cell_drop)
    real(8) :: xc,yc,zc,uc,vc,wc,vol 

    logical :: merge_drop

    if (.not. LPP_initialized) then 
      call initialize_LPP()
      LPP_initialized = .true.
    end if ! LPP_initialized

    tag_id  (:,:,:) = 0
    tag_flag(:,:,:) = 0
    current_id = 1

    drops_merge(:,rank)%num_gcell = 0

    ns_queue = 0
    s_queue(:,:) = 0
    num_drop(:) = 0
    num_drop_merge(:) = 0
    do i=imin,imax; do j=jmin,jmax; do k=kmin,kmax
    !do i=is,ie; do j=js,je; do k=ks,ke
      if ( cvof(i,j,k) > 0.d0 .and. tag_flag(i,j,k) == 0 ) then 
        tag_id  (i,j,k) = current_id
        tag_flag(i,j,k) = 2 ! mark as S node
        num_drop(rank) = num_drop(rank) + 1
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
        num_cell_drop = 0
        cell_list = 0

        merge_drop = .false.
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
                  if ( cvof(isq+i0,jsq+j0,ksq+k0)      > 0.d0 .and. & 
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
               if ( num_cell_drop < maxnum_cell_drop ) then
                  num_cell_drop = num_cell_drop + 1
                  cell_list(1:3,num_cell_drop) = [isq,jsq,ksq]
               else
                  call pariserror('Number of cells of droplet is larger than the maxinum value!')
               end if ! num_cell_drop
            else if ( isq >= Ng+1 .and. isq <= Ng+Nx .and. &     ! block ghost cells 
                      jsq >= Ng+1 .and. jsq <= Ng+Ny .and. & 
                      ksq >= Ng+1 .and. ksq <= Ng+Nz ) then
               tag_flag(isq,jsq,ksq) = 5
               if ( merge_drop .eqv. .false.) then 
                  merge_drop = .true.
                  num_drop      (rank) = num_drop      (rank) - 1
                  num_drop_merge(rank) = num_drop_merge(rank) + 1
               end if ! merge_drop
               if ( drops_merge(num_drop_merge(rank),rank)%num_gcell < maxnum_cell_drop) then 
                  drops_merge(num_drop_merge(rank),rank)%num_gcell = drops_merge(num_drop_merge(rank),rank)%num_gcell + 1
                  drops_merge_gcell_list(1,drops_merge(num_drop_merge(rank),rank)%num_gcell,num_drop_merge(rank)) = isq
                  drops_merge_gcell_list(2,drops_merge(num_drop_merge(rank),rank)%num_gcell,num_drop_merge(rank)) = jsq
                  drops_merge_gcell_list(3,drops_merge(num_drop_merge(rank),rank)%num_gcell,num_drop_merge(rank)) = ksq
               else
                  call pariserror('Number of ghost cells of droplet is larger than the maxinum value!')
               end if ! drops_merge(num_drop_merge(rank),rank)%num_gcell
            else                                                        ! domain ghost cells 
               ! Note: periodic bdry cond, to be added later
            end if ! isq, jsq, ksq
            
          end do ! is_queue
          ! mark all C nodes as S nodes
          if ( nc_queue >= 0 ) then 
            s_queue(:,:) = c_queue(:,:)   ! mark all C nodes as S nodes
            ns_queue = nc_queue
          end if ! nc_queue
        end do ! ns_queue>0

        if ( vol <= 0.d0 ) call pariserror("zero or negative droplet volume!")
        if ( merge_drop ) then 
            drops_merge(num_drop_merge(rank),rank)%element%id  = current_id
            drops_merge(num_drop_merge(rank),rank)%element%vol = vol 
            drops_merge(num_drop_merge(rank),rank)%element%xc  = xc/vol 
            drops_merge(num_drop_merge(rank),rank)%element%yc  = yc/vol 
            drops_merge(num_drop_merge(rank),rank)%element%zc  = zc/vol 
            drops_merge(num_drop_merge(rank),rank)%element%uc  = uc/vol 
            drops_merge(num_drop_merge(rank),rank)%element%vc  = vc/vol
            drops_merge(num_drop_merge(rank),rank)%element%wc  = wc/vol
            drops_merge(num_drop_merge(rank),rank)%num_cell_drop = num_cell_drop
            drops_merge_cell_list(:,:,num_drop_merge(rank)) = cell_list
         else 
            drops(num_drop(rank),rank)%element%id  = current_id
            drops(num_drop(rank),rank)%element%vol = vol 
            drops(num_drop(rank),rank)%element%xc  = xc/vol
            drops(num_drop(rank),rank)%element%yc  = yc/vol
            drops(num_drop(rank),rank)%element%zc  = zc/vol
            drops(num_drop(rank),rank)%element%uc  = uc/vol
            drops(num_drop(rank),rank)%element%vc  = vc/vol
            drops(num_drop(rank),rank)%element%wc  = wc/vol
            drops(num_drop(rank),rank)%num_cell_drop = num_cell_drop
            drops_cell_list(:,:,num_drop(rank)) = cell_list
         end if ! merge_drop
         current_id = current_id+1
      else if ( cvof(i,j,k) == 0.d0 .and. tag_flag(i,j,k) == 0 ) then 
         tag_flag(i,j,k) = 4
      end if ! cvof(i,j,k)
    enddo; enddo; enddo

    num_tag(rank) = num_drop(rank) + num_drop_merge(rank)
   end subroutine tag_drop

   subroutine tag_drop_all()
      include 'mpif.h'
      integer :: i,j,k
      integer :: ierr,irank,num_tag_accu
      integer :: req(4),sta(MPI_STATUS_SIZE,4),MPI_Comm,ireq
      integer , parameter :: ngh=2
      integer , parameter :: root=0

      ! Broadcast num_drop(rank) to all processes
      call MPI_ALLGATHER(num_drop      (rank), 1, MPI_INTEGER, &
                         num_drop      (:)   , 1, MPI_INTEGER, MPI_Comm_World, ierr)
      call MPI_ALLGATHER(num_drop_merge(rank), 1, MPI_INTEGER, &
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
      if ( nPdomain > 1 .and. rank > 0 ) then 
         num_tag_accu = sum(num_drop(0:rank-1),1) + sum(num_drop_merge(0:rank-1),1)
         do i=is,ie; do j=js,je; do k=ks,ke
            if ( tag_id(i,j,k) > 0 ) tag_id(i,j,k)=tag_id(i,j,k)+num_tag_accu
         end do; end do; end do
      end if ! rank

      call ighost_x(tag_id,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
      call ighost_y(tag_id,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
      call ighost_z(tag_id,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)

   end subroutine tag_drop_all

   subroutine merge_drop_pieces

      integer :: idrop,iCell,idrop1
      integer :: idiff_tag,tag,tag1
      real(8) :: vol_merge,xc_merge,yc_merge,zc_merge,uc_merge,vc_merge,wc_merge,vol1
      integer :: irank, irank1
      real(8) :: max_drop_merge_vol
      integer :: tag_max_drop_merge_vol

      call CreateTag2DropTable

      ! Check ghost cells of droplet pieces
      if ( num_drop_merge(rank) > 0 ) then 
      do idrop = 1, num_drop_merge(rank) 
         drops_merge(idrop,rank)%num_diff_tag = 1 
         drops_merge(idrop,rank)%diff_tag_list(drops_merge(idrop,rank)%num_diff_tag) &
            = tag_id(drops_merge_gcell_list(1,1,idrop), &
                     drops_merge_gcell_list(2,1,idrop), &
                     drops_merge_gcell_list(3,1,idrop))
         if ( drops_merge(idrop,rank)%num_gcell > 1 ) then  
            do iCell = 2,drops_merge(idrop,rank)%num_gcell
               do idiff_tag = 1, drops_merge(idrop,rank)%num_diff_tag
                  if ( tag_id(drops_merge_gcell_list(1,iCell,idrop), & 
                              drops_merge_gcell_list(2,iCell,idrop), & 
                              drops_merge_gcell_list(3,iCell,idrop)) &
                    == drops_merge(idrop,rank)%diff_tag_list(idiff_tag)) exit
               end do ! idiff_tag
               if ( idiff_tag == drops_merge(idrop,rank)%num_diff_tag + 1 ) then 
                  drops_merge(idrop,rank)%num_diff_tag = &
                  drops_merge(idrop,rank)%num_diff_tag + 1
                  drops_merge(idrop,rank)%diff_tag_list(drops_merge(idrop,rank)%num_diff_tag) &
                  = tag_id(drops_merge_gcell_list(1,iCell,idrop), &
                           drops_merge_gcell_list(2,iCell,idrop), &
                           drops_merge_gcell_list(3,iCell,idrop))
               end if ! idiff_tag
            end do ! iCell
         end if ! drops_merge(idrop,irank)%num_gcell
      end do ! idrop
      end if ! num_drop_merge(rank)

      call CollectDropMerge

      ! merge droplets pieces & calculate droplet properties
      if ( rank == 0 ) then 
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
               end do ! idiff_tag
            end if ! newdropid(tag)
         end do ! tag
         totalnum_drop_indep = sum(num_drop)
         new_drop_id = new_drop_id + totalnum_drop_indep
         totalnum_drop = totalnum_drop_indep + num_new_drop
      end if ! rank 

      call distributeDropMerge
   
      ! finalize
      call ReleaseTag2DropTable
      if ( rank == 0 ) deallocate( new_drop_id )

   end subroutine merge_drop_pieces

   subroutine drop_statistics(tswap)
      
      include 'mpif.h'
      integer, intent(in) :: tswap
      integer :: req(2),sta(MPI_STATUS_SIZE,2),MPI_Comm,ireq,ierr
      integer :: MPI_element_type, oldtypes(0:1), blockcounts(0:1), & 
                 offsets(0:1), extent,r8extent, MPI_element_row 
      integer :: maxnum_element 
      integer :: irank, idrop, ielement, ielem_plot
      integer :: num_element_estimate(0:nPdomain-1)
      type(element) :: element_NULL

      num_element_estimate = num_drop+num_drop_merge+num_part
      maxnum_element = maxval(num_element_estimate) 
      allocate( element_stat(maxnum_element,0:nPdomain-1) )

      element_NULL%xc = 0.d0;element_NULL%yc = 0.d0;element_NULL%zc = 0.d0
      element_NULL%uc = 0.d0;element_NULL%vc = 0.d0;element_NULL%wc = 0.d0
      element_NULL%vol = 0.d0;element_NULL%id = CRAZY_INT

      !  Setup MPI derived type for element_type 
      offsets (0) = 0 
      oldtypes(0) = MPI_REAL8 
      blockcounts(0) = 7 
      call MPI_TYPE_EXTENT(MPI_REAL8, r8extent, ierr) 
      offsets    (1) = blockcounts(0)*r8extent 
      oldtypes   (1) = MPI_INTEGER  
      blockcounts(1) = 1  

      call MPI_TYPE_STRUCT(2, blockcounts, offsets, oldtypes, & 
                           MPI_element_type, ierr) 
      call MPI_TYPE_COMMIT(MPI_element_type, ierr)

      !  initialize element_stat
      num_element(:) = 0
      element_stat(:,:) = element_NULL
      if ( num_drop(rank) > 0 ) then 
         do idrop = 1,num_drop(rank) 
            num_element(rank) = num_element(rank) + 1
            element_stat(num_element(rank),rank)%vol = drops(idrop,rank)%element%vol
            element_stat(num_element(rank),rank)%xc  = drops(idrop,rank)%element%xc
            element_stat(num_element(rank),rank)%yc  = drops(idrop,rank)%element%yc
            element_stat(num_element(rank),rank)%zc  = drops(idrop,rank)%element%zc
            element_stat(num_element(rank),rank)%uc  = drops(idrop,rank)%element%uc
            element_stat(num_element(rank),rank)%vc  = drops(idrop,rank)%element%vc
            element_stat(num_element(rank),rank)%wc  = drops(idrop,rank)%element%wc
            element_stat(num_element(rank),rank)%id  = drops(idrop,rank)%element%id
         end do ! idrop
      end if ! num_drop(rank)
      
      if ( num_drop_merge(rank) > 0 ) then 
         do idrop = 1,num_drop_merge(rank)
            if ( drops_merge(idrop,rank)%flag_center_mass == 1 ) then  
               num_element(rank) = num_element(rank) + 1
               element_stat(num_element(rank),rank)%vol = drops_merge(idrop,rank)%element%vol
               element_stat(num_element(rank),rank)%xc  = drops_merge(idrop,rank)%element%xc
               element_stat(num_element(rank),rank)%yc  = drops_merge(idrop,rank)%element%yc
               element_stat(num_element(rank),rank)%zc  = drops_merge(idrop,rank)%element%zc
               element_stat(num_element(rank),rank)%uc  = drops_merge(idrop,rank)%element%uc
               element_stat(num_element(rank),rank)%vc  = drops_merge(idrop,rank)%element%vc
               element_stat(num_element(rank),rank)%wc  = drops_merge(idrop,rank)%element%wc
               element_stat(num_element(rank),rank)%id  = drops_merge(idrop,rank)%element%id
            end if !drops_merge(idrop,rank)%flag_center_mass
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
            element_stat(num_element(rank),rank)%id  = parts(idrop,rank)%element%id
         end do ! idrop
      end if ! num_drop(rank)
 
      ! Collect all discrete elements to rank 0
      call MPI_ALLGATHER(num_element(rank), 1, MPI_INTEGER, &
                         num_element(:),    1, MPI_INTEGER, MPI_Comm_World, ierr)
      call MPI_TYPE_CONTIGUOUS (maxnum_element, MPI_element_type, MPI_element_row, ierr)
      call MPI_TYPE_COMMIT(MPI_element_row, ierr)

      if ( rank > 0 ) then
         call MPI_ISEND(element_stat(1:maxnum_element,rank),1, & 
                        MPI_element_row, 0, 14, MPI_COMM_WORLD, req(1), ierr)
         call MPI_WAIT(req(1),sta(:,1),ierr)
      else
         do irank = 1,nPdomain-1
            call MPI_IRECV(element_stat(1:maxnum_element,irank),1, & 
                           MPI_element_row, irank, 14, MPI_COMM_WORLD, req(2), ierr)
            call MPI_WAIT(req(2),sta(:,2),ierr)
         end do ! irank
      end if ! rank

      ! Statisitcs or other operation at rank 0 
      if ( rank == 0 ) then
         ielem_plot = 0 
         do irank = 0,nPdomain-1
            if ( num_element(irank) > 0 ) then 
               do ielement = 1, num_element(irank) 
                  ielem_plot = ielem_plot + 1 
                  OPEN(UNIT=200+ielem_plot,FILE=TRIM(out_path)//'/element-'//TRIM(int2text(ielem_plot,padding))//'.dat')
                  write(200+ielem_plot,*) tswap,element_stat(ielement,irank)%xc, & 
                                                element_stat(ielement,irank)%yc, &
                                                element_stat(ielement,irank)%yc, &
                                                element_stat(ielement,irank)%uc, &
                                                element_stat(ielement,irank)%vc, &
                                                element_stat(ielement,irank)%wc, &
                                                element_stat(ielement,irank)%vol
!                  CLOSE(UNIT=200+ielem_plot)
               end do ! ielement
            end if ! num_element_irank) 
         end do ! irank
      end if ! rank

      ! finalize
      call MPI_TYPE_FREE(MPI_element_type, ierr)
      deallocate( element_stat )

   end subroutine drop_statistics

   subroutine CollectDropMerge

      include 'mpif.h'
      integer :: req(2),sta(MPI_STATUS_SIZE,2),MPI_Comm,ireq,ierr
      integer :: MPI_drop_merge_comm_type, oldtypes(0:1), blockcounts(0:1), & 
                 offsets(0:1), extent,r8extent 
      integer :: maxnum_drop_merge 
      integer :: irank, idrop

      maxnum_drop_merge = maxval(num_drop_merge)
      allocate( drops_merge_comm(maxnum_drop_merge,0:nPdomain-1) )

      !  Setup MPI derived type for drop_merge_comm
      offsets (0) = 0 
      oldtypes(0) = MPI_REAL8 
      blockcounts(0) = 7 
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
            drops_merge_comm(idrop,rank)%id            = drops_merge(idrop,rank)%element%id
            drops_merge_comm(idrop,rank)%num_diff_tag  = drops_merge(idrop,rank)%num_diff_tag
            drops_merge_comm(idrop,rank)%diff_tag_list = drops_merge(idrop,rank)%diff_tag_list
            drops_merge_comm(idrop,rank)%vol           = drops_merge(idrop,rank)%element%vol
            drops_merge_comm(idrop,rank)%xc            = drops_merge(idrop,rank)%element%xc
            drops_merge_comm(idrop,rank)%yc            = drops_merge(idrop,rank)%element%yc
            drops_merge_comm(idrop,rank)%zc            = drops_merge(idrop,rank)%element%zc
            drops_merge_comm(idrop,rank)%uc            = drops_merge(idrop,rank)%element%uc
            drops_merge_comm(idrop,rank)%vc            = drops_merge(idrop,rank)%element%vc
            drops_merge_comm(idrop,rank)%wc            = drops_merge(idrop,rank)%element%wc
         end do ! idrop
      end if ! num_drop_merge(rank)

      ! Send all droplet pieces to rank 0
      if ( rank > 0 ) then
         if ( num_drop_merge(rank) > 0 ) then 
            call MPI_ISEND(drops_merge_comm(1:num_drop_merge(rank),rank),num_drop_merge(rank), & 
                           MPI_drop_merge_comm_type, 0,    13, MPI_COMM_WORLD, req(1), ierr)
            call MPI_WAIT(req(1),sta(:,1),ierr)
         end if ! num_drop_merge(rank)
      else
         do irank = 1,nPdomain-1
            if ( num_drop_merge(irank) > 0 ) then 
               call MPI_IRECV(drops_merge_comm(1:num_drop_merge(irank),irank),num_drop_merge(irank), & 
                              MPI_drop_merge_comm_type, irank, 13, MPI_COMM_WORLD, req(2), ierr)
               call MPI_WAIT(req(2),sta(:,2),ierr)
            end if ! num_drop_merge(irank)
         end do ! irank
      end if ! rank

      ! finalize
      call MPI_TYPE_FREE(MPI_drop_merge_comm_type, ierr)

   end subroutine CollectDropMerge
   
   subroutine DistributeDropMerge

      include 'mpif.h'
      integer :: irank,idrop
      integer :: req(2),sta(MPI_STATUS_SIZE,2),MPI_Comm,ireq,ierr
      integer :: MPI_drop_merge_comm_type, oldtypes(0:1), blockcounts(0:1), & 
                 offsets(0:1), extent,r8extent 

      !  Setup MPI derived type for drop_merge_comm
      offsets (0) = 0 
      oldtypes(0) = MPI_REAL8 
      blockcounts(0) = 7 
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
            call MPI_IRECV(drops_merge_comm(1:num_drop_merge(rank),rank),num_drop_merge(rank), & 
                           MPI_drop_merge_comm_type, 0,    13, MPI_COMM_WORLD, req(1), ierr)
            call MPI_WAIT(req(1),sta(:,1),ierr)
         end if ! num_drop_merge(rank)
      else
         do irank = 1,nPdomain-1
               if ( num_drop_merge(irank) > 0 ) then 
               call MPI_ISEND(drops_merge_comm(1:num_drop_merge(irank),irank),num_drop_merge(irank), & 
                              MPI_drop_merge_comm_type, irank, 13, MPI_COMM_WORLD, req(2), ierr)
               call MPI_WAIT(req(2),sta(:,2),ierr)
            end if ! num_drop_merge(irank)
         end do ! irank
      end if ! rank

      if ( num_drop_merge(rank) > 0 ) then
         do idrop = 1, num_drop_merge(rank)
            drops_merge(idrop,rank)%element%vol = drops_merge_comm(idrop,rank)%vol
            drops_merge(idrop,rank)%element%xc  = drops_merge_comm(idrop,rank)%xc
            drops_merge(idrop,rank)%element%yc  = drops_merge_comm(idrop,rank)%yc
            drops_merge(idrop,rank)%element%zc  = drops_merge_comm(idrop,rank)%zc
            drops_merge(idrop,rank)%element%uc  = drops_merge_comm(idrop,rank)%uc
            drops_merge(idrop,rank)%element%vc  = drops_merge_comm(idrop,rank)%vc
            drops_merge(idrop,rank)%element%wc  = drops_merge_comm(idrop,rank)%wc
            drops_merge(idrop,rank)%flag_center_mass  = drops_merge_comm(idrop,rank)%flag_center_mass
         end do ! idrop
      end if ! num_drop_merge(rank)

      ! finalize
      deallocate(drops_merge_comm)
      call MPI_TYPE_FREE(MPI_drop_merge_comm_type, ierr) 

   end subroutine DistributeDropMerge

   subroutine CreateTag2DropTable

      include 'mpif.h'
      integer :: req(6),sta(MPI_STATUS_SIZE,6),MPI_Comm,ireq,ierr
      integer :: idrop, irank

      allocate( tag_dropid   (1:total_num_tag) )
      allocate( tag_rank     (1:total_num_tag) )
      allocate( tag_mergeflag(1:total_num_tag) )

      if ( num_drop(rank) > 0 ) then 
      do idrop = 1,num_drop(rank)
         drops(idrop,rank)%element%id = tag_id( drops_cell_list(1,1,idrop), &
                                                drops_cell_list(2,1,idrop), &
                                                drops_cell_list(3,1,idrop) )
         tag_dropid   (drops(idrop,rank)%element%id) = idrop
         tag_rank     (drops(idrop,rank)%element%id) = rank 
         tag_mergeflag(drops(idrop,rank)%element%id) = 0 
      end do ! idrop
   end if ! num_drop(rank)
      
      if ( num_drop_merge(rank) > 0 ) then 
      do idrop = 1,num_drop_merge(rank)
         drops_merge(idrop,rank)%element%id = tag_id( drops_merge_cell_list(1,1,idrop), &
                                                      drops_merge_cell_list(2,1,idrop), &
                                                      drops_merge_cell_list(3,1,idrop) )
         tag_dropid   (drops_merge(idrop,rank)%element%id) = idrop
         tag_rank     (drops_merge(idrop,rank)%element%id) = rank 
         tag_mergeflag(drops_merge(idrop,rank)%element%id) = 1 
      end do ! idrop
   end if ! num_drop_merge(rank) 

      ! MPI communication for tag tables
      if ( rank > 0 ) then
         if ( num_tag(rank) > 0 ) then 
            call MPI_ISEND(tag_dropid   (tagmin(rank):tagmax(rank)),num_tag(rank), & 
                           MPI_INTEGER, 0,    10, MPI_COMM_WORLD, req(1), ierr)
            call MPI_ISEND(tag_rank     (tagmin(rank):tagmax(rank)),num_tag(rank), & 
                           MPI_INTEGER, 0,    11, MPI_COMM_WORLD, req(2), ierr)
            call MPI_ISEND(tag_mergeflag(tagmin(rank):tagmax(rank)),num_tag(rank), & 
                           MPI_INTEGER, 0,    12, MPI_COMM_WORLD, req(3), ierr)
            call MPI_WAITALL(3,req(1:3),sta(:,1:3),ierr)
         end if ! num_tag(rank)
      else
         ireq = 0 
         do irank = 1,nPdomain-1
            if ( num_tag(irank) > 0 ) then
               call MPI_IRECV(tag_dropid   (tagmin(irank):tagmax(irank)),num_tag(irank), & 
                              MPI_INTEGER, irank, 10, MPI_COMM_WORLD, req(4), ierr)
               call MPI_IRECV(tag_rank     (tagmin(irank):tagmax(irank)),num_tag(irank), & 
                              MPI_INTEGER, irank, 11, MPI_COMM_WORLD, req(5), ierr)
               call MPI_IRECV(tag_mergeflag(tagmin(irank):tagmax(irank)),num_tag(irank), & 
                              MPI_INTEGER, irank, 12, MPI_COMM_WORLD, req(6), ierr)
               call MPI_WAITALL(3,req(4:6),sta(:,4:6),ierr)
            end if ! num_tag(irank)
         end do ! irank
      end if ! rank

   end subroutine CreateTag2DropTable

   subroutine ReleaseTag2DropTable

      deallocate(tag_dropid   )
      deallocate(tag_rank     )
      deallocate(tag_mergeflag)

   end subroutine ReleaseTag2DropTable

! ==============================================
! output tag of droplets
! ==============================================
   subroutine output_tag()
      integer :: i,j,k

      call output_VOF(0,imin,imax,jmin,jmax,kmin,kmax)
      call tag_drop()
      if (nPdomain > 1 ) then 
         call tag_drop_all()
         call merge_drop_pieces 
      end if ! nPdomain

!      OPEN(UNIT=89,FILE=TRIM(out_path)//'/tag-'//TRIM(int2text(rank,padding))//'.txt')
      OPEN(UNIT=90,FILE=TRIM(out_path)//'/tag-tecplot'//TRIM(int2text(rank,padding))//'.dat')

!      do j=jmax,jmin,-1
!         write(89,'(36(I5,1X))') tag_id(imin:imax,j,Nz/4+Ng)
!      end do ! j
!      CLOSE(89)

      write(90,*) 'title= " 3d tag "'
      write(90,*) 'variables = "x", "y", "z", "tag", "c" '
      write(90,*) 'zone i=,',Nx/nPx+Ng*2, 'j=',Ny/nPy+Ng*2, 'k=',Nz/nPz+Ng*2,'f=point'
      do k = kmin,kmax
         do j=jmin,jmax
            do i=imin,imax 
               write(90,'(4(I5,1X),(E15.8))') i,j,k,tag_id(i,j,k),cvof(i,j,k)
            end do ! i
         end do ! j
      end do ! k

!      write(90,*) 'title= " 2d tag "'
!      write(90,*) 'variables = "x", "y", "tag", "c" '
!      write(90,*) 'zone i=,',Nx/nPx+Ng*2, 'j=',Ny/nPy+Ng*2, 'f=point'
!      k = Nz/4+Ng
!         do j=jmin,jmax
!            do i=imin,imax 
!               write(90,'(2(E15.8,1X),(I5,1X),(E15.8))') x(i),y(j),tag_id(i,j,k),cvof(i,j,k)
!            end do ! i
!        end do ! j
!      CLOSE(90)
!      if ( kmin < k .and. kmax > k ) then 
!      write(*,*) ' ************************ ', rank
!      do j = jmax,jmin,-1
!         write(*,'(20(I2,1X))') tag_id(imin:imax,j,k)
!      end do ! j
!      write(*,*) ' ************************ '
!      end if ! k 
   end subroutine output_tag

! ==============================================
! output droplets & particles 
! ==============================================
   subroutine output_DP(tswap)

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
         write(91,*) drops   (ib,rank)%element%xc, drops   (ib,rank)%element%vol
         write(92,*) drops_ex(ib)%element%xc, drops_ex(ib)%element%vol
      end do ! i 
      CLOSE(91)
      CLOSE(92)
      end if ! rank

      ! convert droplets to particles
      call ConvertDrop2Part(tswap)

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
   end subroutine output_DP

! ===============================================
! Testing section
! ===============================================
   subroutine test_Lag_part(tswap)
      include 'mpif.h'
      integer, intent(in) :: tswap
!      integer :: ierr
                     
      if(test_tag) then
         call output_tag()
      else if ( test_D2P ) then
         call output_DP(tswap)
      end if
                                                                             
!      ! Exit MPI gracefully
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      call MPI_finalize(ierr)
!      stop
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

   subroutine ConvertDrop2Part(tswap)

      include 'mpif.h'

      integer, intent(in) :: tswap

      integer :: idrop,ilist
      logical :: ConvertDropFlag,convertDoneFlag
      real(8) :: MinDistPart2CellCenter, DistPart2CellCenter
      integer :: ierr
      real(8) :: uf,vf,wf
      integer :: i,j,k

      if (.not. LPP_initialized) then 
         call initialize_LPP()
         LPP_initialized = .true.
      end if ! LPP_initialized

      if ( num_drop(rank) > 0 ) then 
      do idrop = 1,num_drop(rank)

         call CheckConvertDropCriteria(drops(idrop,rank)%element%vol, & 
                                       drops(idrop,rank)%element%xc,  & 
                                       drops(idrop,rank)%element%yc,  & 
                                       drops(idrop,rank)%element%zc,  &
                                       ConvertDropFlag,CriteriaConvertCase)

         if ( ConvertDropFlag ) then
! TEMPORARY
            write(*,*) 'Drop is converted to particle', idrop,rank,tswap
! END TEMPORARY
            ! transfer droplet properties to particle
            num_part(rank) = num_part(rank) + 1
            parts(num_part(rank),rank)%element = drops(idrop,rank)%element

            ! compute average fluid quantities
            ! Note: XXX temporary, need to be improved later 
            uf = drops(idrop,rank)%element%uc
            vf = drops(idrop,rank)%element%vc
            wf = drops(idrop,rank)%element%wc

            ! remove droplet vof structure
            MinDistPart2CellCenter = 1.0d10
            do ilist = 1,drops(idrop,rank)%num_cell_drop
               cvof(drops_cell_list(1,ilist,idrop), &
                    drops_cell_list(2,ilist,idrop), &
                    drops_cell_list(3,ilist,idrop)) = 0.0
                  u(drops_cell_list(1,ilist,idrop), &
                    drops_cell_list(2,ilist,idrop), &
                    drops_cell_list(3,ilist,idrop)) = uf 
                  v(drops_cell_list(1,ilist,idrop), &
                    drops_cell_list(2,ilist,idrop), &
                    drops_cell_list(3,ilist,idrop)) = vf 
                  w(drops_cell_list(1,ilist,idrop), &
                    drops_cell_list(2,ilist,idrop), &
                    drops_cell_list(3,ilist,idrop)) = wf 
               DistPart2CellCenter = ( drops(idrop,rank)%element%xc  & 
                           - x(drops_cell_list(1,ilist,idrop)))**2.d0 & 
                                   + ( drops(idrop,rank)%element%yc  & 
                           - y(drops_cell_list(2,ilist,idrop)))**2.d0 &
                                   + ( drops(idrop,rank)%element%zc  & 
                           - z(drops_cell_list(3,ilist,idrop)))**2.d0
               if ( DistPart2CellCenter < MinDistPart2CellCenter ) then 
                  MinDistPart2CellCenter = DistPart2CellCenter
                  parts(num_part(rank),rank)%ic = drops_cell_list(1,ilist,idrop)
                  parts(num_part(rank),rank)%jc = drops_cell_list(2,ilist,idrop)
                  parts(num_part(rank),rank)%kc = drops_cell_list(3,ilist,idrop)
               end if !DistPart2CellCenter
            end do ! ilist
         end if !ConvertDropFlag
      end do ! idrop
      end if ! num_drop(rank) 
      
      if ( num_drop_merge(rank) > 0 ) then
      do idrop = 1,num_drop_merge(rank)

         call CheckConvertDropCriteria(drops_merge(idrop,rank)%element%vol, & 
                                       drops_merge(idrop,rank)%element%xc,  & 
                                       drops_merge(idrop,rank)%element%yc,  & 
                                       drops_merge(idrop,rank)%element%zc,  &
                                       ConvertDropFlag,CriteriaConvertCase)

         if ( ConvertDropFlag ) then
! TEMPORARY
            write(*,*) 'Drop_merge is converted to particle', idrop,rank,tswap
! END TEMPORARY

            ! compute average fluid quantities
            ! Note: XXX temporary, need to be improved later 
            uf = drops_merge(idrop,rank)%element%uc
            vf = drops_merge(idrop,rank)%element%vc
            wf = drops_merge(idrop,rank)%element%wc

            ! remove droplet vof structure
            do ilist = 1,drops_merge(idrop,rank)%num_cell_drop
               cvof(drops_merge_cell_list(1,ilist,idrop), &
                    drops_merge_cell_list(2,ilist,idrop), &
                    drops_merge_cell_list(3,ilist,idrop)) = 0.0
                  u(drops_merge_cell_list(1,ilist,idrop), &
                    drops_merge_cell_list(2,ilist,idrop), &
                    drops_merge_cell_list(3,ilist,idrop)) = uf
                  v(drops_merge_cell_list(1,ilist,idrop), &
                    drops_merge_cell_list(2,ilist,idrop), &
                    drops_merge_cell_list(3,ilist,idrop)) = vf
                  w(drops_merge_cell_list(1,ilist,idrop), &
                    drops_merge_cell_list(2,ilist,idrop), &
                    drops_merge_cell_list(3,ilist,idrop)) = wf
            end do ! ilist
      
            ! remove droplet vof structure
            do ilist = 1,drops_merge(idrop,rank)%num_gcell
               cvof(drops_merge_gcell_list(1,ilist,idrop), &
                    drops_merge_gcell_list(2,ilist,idrop), &
                    drops_merge_gcell_list(3,ilist,idrop)) = 0.0
                  u(drops_merge_gcell_list(1,ilist,idrop), &
                    drops_merge_gcell_list(2,ilist,idrop), &
                    drops_merge_gcell_list(3,ilist,idrop)) = uf
                  v(drops_merge_gcell_list(1,ilist,idrop), &
                    drops_merge_gcell_list(2,ilist,idrop), &
                    drops_merge_gcell_list(3,ilist,idrop)) = vf
                  w(drops_merge_gcell_list(1,ilist,idrop), &
                    drops_merge_gcell_list(2,ilist,idrop), &
                    drops_merge_gcell_list(3,ilist,idrop)) = wf
            end do ! ilist

            ! transfer droplet properties to particle if center of mass located
            ! in this droplet piece
            if ( drops_merge(idrop,rank)%flag_center_mass == 1 ) then 
               num_part(rank) = num_part(rank) + 1
               parts(num_part(rank),rank)%element = drops_merge(idrop,rank)%element
            
               ! Find particle location cell
               MinDistPart2CellCenter = 1.0d10
               do ilist = 1,drops_merge(idrop,rank)%num_cell_drop
                  DistPart2CellCenter = ( drops_merge(idrop,rank)%element%xc  & 
                              - x(drops_merge_cell_list(1,ilist,idrop)))**2.d0 & 
                                      + ( drops_merge(idrop,rank)%element%yc  & 
                              - y(drops_merge_cell_list(2,ilist,idrop)))**2.d0 &
                                      + ( drops_merge(idrop,rank)%element%zc  & 
                              - z(drops_merge_cell_list(3,ilist,idrop)))**2.d0
                  if ( DistPart2CellCenter < MinDistPart2CellCenter ) then 
                     MinDistPart2CellCenter = DistPart2CellCenter
                     parts(num_part(rank),rank)%ic = drops_merge_cell_list(1,ilist,idrop)
                     parts(num_part(rank),rank)%jc = drops_merge_cell_list(2,ilist,idrop)
                     parts(num_part(rank),rank)%kc = drops_merge_cell_list(3,ilist,idrop)
                  end if !DistPart2CellCenter
               end do ! ilist
            end if ! flag_center_mass 

         end if !ConvertDropFlag
      end do ! idrop
      end if ! num_drop_merge(rank) 

      ! Update num_part to all ranks. Note: no need for num_drop &
      ! num_drop_merge since they will be regenerated next step  
      call MPI_ALLGATHER(num_part(rank), 1, MPI_INTEGER, &
                         num_part(:)   , 1, MPI_INTEGER, MPI_Comm_World, ierr)

   end subroutine ConvertDrop2Part

   subroutine CheckConvertDropCriteria(vol,xc,yc,zc,ConvertDropFlag,CriteriaConvertCase)

      real(8), intent(in ) :: vol,xc,yc,zc
      integer, intent(in ) :: CriteriaConvertCase
      logical, intent(out) :: ConvertDropFlag

      select case (CriteriaConvertCase) 
         case (CriteriaRectangle)
            if ( (vol < vol_cut)  .and. &
                 (xc  > xlpp_min) .and. (xc  < xlpp_max) .and. & 
                 (yc  > ylpp_min) .and. (yc  < ylpp_max) .and. & 
                 (zc  > zlpp_min) .and. (zc  < zlpp_max) ) then 
               ConvertDropFlag = .true.
            else 
               ConvertDropFlag = .false.
            end if !vol_cut, xlpp_min...
         case (CriteriaCylinder )   ! Note: assuming axis along x-direction
                                    ! radius indicated by ylpp_min & ylpp_max
            if ( (  vol < vol_cut)                          .and. &
                 ( (yc**2.d0 + zc**2.d0) > ylpp_min**2.d0)  .and. & 
                 ( (yc**2.d0 + zc**2.d0) < ylpp_max**2.d0) ) then 
               ConvertDropFlag = .true.
            else 
               ConvertDropFlag = .false.
            end if !vol_cut, xlpp_min... 
         case (CriteriaSphere   )   ! radius indicated by ylpp_min & ylpp_max 
            if ( (  vol < vol_cut)                                      .and. &
                 ( (xc**2.d0 + yc**2.d0 + zc**2.d0) > ylpp_min**2.d0)   .and. & 
                 ( (xc**2.d0 + yc**2.d0 + zc**2.d0) < ylpp_max**2.d0) ) then 
               ConvertDropFlag = .true.
            else 
               ConvertDropFlag = .false.
            end if !vol_cut, xlpp_min... 
      end select
      
   end subroutine CheckConvertDropCriteria

   subroutine ComputePartForce(tswap)

      integer, intent(in) :: tswap

      integer, parameter :: drag_model_Stokes = 1
      integer, parameter :: drag_model_SN = 2     ! Schiller & Nauman
      integer, parameter :: drag_model_CG = 3     ! Clift & Gauvin

      real(8), parameter :: Cm = 0.5d0

      real(8) :: relvel(4), partforce(3)
      real(8) :: dp, Rep, muf, phi, rhof, rhop, taup
      real(8) :: up,vp,wp, uf,vf,wf, DufDt,DvfDt,DwfDt
      real(8) :: fhx,fhy,fhz

      integer :: ipart
      if ( num_part(rank) > 0 ) then
         do ipart = 1,num_part(rank)
            up = parts(ipart,rank)%element%uc
            vp = parts(ipart,rank)%element%vc
            wp = parts(ipart,rank)%element%wc

            call GetFluidProp(parts(ipart,rank)%ic, &
                              parts(ipart,rank)%jc, &
                              parts(ipart,rank)%kc, & 
                              parts(ipart,rank)%element%xc, & 
                              parts(ipart,rank)%element%yc, & 
                              parts(ipart,rank)%element%zc, & 
                              uf,vf,wf,DufDt,DvfDt,DwfDt,   & 
                              parts(ipart,rank)%element%vol)

            relvel(1) = uf - up
            relvel(2) = vf - vp
            relvel(3) = wf - wp 
            relvel(4) = sqrt(relvel(1)**2.0d0 + relvel(2)**2.0d0 + relvel(3)**2.0)

            dp = (parts(ipart,rank)%element%vol*6.d0/PI)**(1.d0/3.d0)

            rhof = rho1
            rhop = rho2
            muf  = mu1
            Rep  = rhof*relvel(4)*dp/muf
            taup = rhop *dp*dp/18.0d0/muf

            select case ( dragmodel ) 
               case ( drag_model_Stokes ) 
                  phi = 1.d0
               case ( drag_model_SN ) 
                  phi = 1.d0+0.15d0*Rep**0.687d0
               case ( drag_model_CG ) 
                  phi = 1.d0+0.15d0*Rep**0.687d0 & 
                      + 1.75d-2*Rep/(1.0d0 + 4.25d4/Rep**1.16d0)
               case default
                  call pariserror("wrong quasi-steady drag model!")
            end select ! dragmodel
! TEMPORARY 
!            phi = phi * (1.d0+2.d0*0.03351)/(1.d0-0.03351)**2.d0
! END TEMPORARY 
            ! Note: set history to be zero for now
            fhx=0.d0; fhy=0.d0; fhz=0.d0

            partforce(1) =(relvel(1)/taup*phi + (1.d0-rhof/rhop)*Gx  &   
                         + (1.d0+Cm)*rhof/rhop*DufDt                 &  
                         + fhx )/(1.d0+Cm*rhof/rhop)
            partforce(2) =(relvel(2)/taup*phi + (1.d0-rhof/rhop)*Gy  &
                         + (1.d0+Cm)*rhof/rhop*DvfDt                 &
                         + fhy )/(1.d0+Cm*rhof/rhop)
            partforce(3) =(relvel(3)/taup*phi + (1.d0-rhof/rhop)*Gz  &
                         + (1.d0+Cm)*rhof/rhop*DwfDt                 &
                         + fhz )/(1.d0+Cm*rhof/rhop)

            parts(ipart,rank)%fx = partforce(1)
            parts(ipart,rank)%fy = partforce(2)
            parts(ipart,rank)%fz = partforce(3)
         end do ! ipart
      end if ! num_part(rank) 
   end subroutine ComputePartForce

   subroutine GetFluidProp(ip,jp,kp,xp,yp,zp,uf,vf,wf,DufDt,DvfDt,DwfDt,volp) 

      integer, intent(in)  :: ip,jp,kp
      real(8), intent(in)  :: xp,yp,zp 
      real(8), intent(in)  :: volp 
      real(8), intent(out) :: uf,vf,wf,DufDt,DvfDt,DwfDt
      
      real(8) :: dp, dx,dy,dz,max_gridsize,min_gridsize
      real(8) :: Lx,Ly,Lz

      dp = (6.d0*volp/PI)**(1.d0/3.d0)
      dx = xh(ip)-xh(ip-1)
      dy = yh(jp)-yh(jp-1)
      dz = zh(kp)-zh(kp-1)
      max_gridsize = max(dx,dy,dz)
      min_gridsize = min(dx,dy,dz)

      if ( dp < min_gridsize ) then ! Interploation
         call TrilinearIntrplFluidVel(xp,yp,zp,ip,jp,kp,uf,DufDt,1)
         call TrilinearIntrplFluidVel(xp,yp,zp,ip,jp,kp,vf,DvfDt,2)
         call TrilinearIntrplFluidVel(xp,yp,zp,ip,jp,kp,wf,DwfDt,3)
      else if ( dp > max_gridsize ) then  ! Check flow scale scale
         Lx = dx*(u(ip-1,jp,kp)+u(ip,jp,kp))/(u(ip,jp,kp)-u(ip-1,jp,kp))
         Ly = dy*(v(ip,jp-1,kp)+v(ip,jp,kp))/(v(ip,jp,kp)-v(ip,jp-1,kp))
         Lz = dz*(w(ip,jp,kp-1)+w(ip,jp,kp))/(w(ip,jp,kp)-w(ip,jp,kp-1))
         if ( Lx > dp ) then 
            call TrilinearIntrplFluidVel(xp,yp,zp,ip,jp,kp,uf,DufDt,1)
         else 
            !ComputeAveFluidVel
         end if ! Lx
         if ( Ly > dp ) then 
            call TrilinearIntrplFluidVel(xp,yp,zp,ip,jp,kp,vf,DvfDt,2)
         else 
            !ComputeAveFluidVel
         end if ! Lx
         if ( Lz > dp ) then 
            call TrilinearIntrplFluidVel(xp,yp,zp,ip,jp,kp,wf,DwfDt,3)
         else 
            !ComputeAveFluidVel
         end if ! Lz
      else ! interploation & averaging
      end if ! dp 
   end subroutine GetFluidProp

   subroutine UpdatePartSol(tswap)
  
      integer, intent(in) :: tswap

      if ( num_part(rank) > 0 ) then 
         parts(1:num_part(rank),rank)%element%xc = parts(1:num_part(rank),rank)%element%xc +& 
                                                   parts(1:num_part(rank),rank)%element%uc*dt 
         parts(1:num_part(rank),rank)%element%yc = parts(1:num_part(rank),rank)%element%yc +& 
                                                   parts(1:num_part(rank),rank)%element%vc*dt 
         parts(1:num_part(rank),rank)%element%zc = parts(1:num_part(rank),rank)%element%zc +& 
                                                   parts(1:num_part(rank),rank)%element%wc*dt 

         parts(1:num_part(rank),rank)%element%uc = parts(1:num_part(rank),rank)%element%uc +&
                                                   parts(1:num_part(rank),rank)%fx*dt 
         parts(1:num_part(rank),rank)%element%vc = parts(1:num_part(rank),rank)%element%vc +&
                                                   parts(1:num_part(rank),rank)%fy*dt 
         parts(1:num_part(rank),rank)%element%wc = parts(1:num_part(rank),rank)%element%wc +&
                                                   parts(1:num_part(rank),rank)%fz*dt
        
         call UpdatePartLocCell   
      end if ! num_part(rank)
      if ( nPdomain > 1 ) then  
         call CollectPartCrossBlocks   
         call TransferPartCrossBlocks   
      end if ! nPdomain
      call SetPartBC
   end subroutine UPdatePartSol

   subroutine StoreOldPartSol()

      parts(:,rank)%xcOld = parts(:,rank)%element%xc 
      parts(:,rank)%ycOld = parts(:,rank)%element%yc 
      parts(:,rank)%zcOld = parts(:,rank)%element%zc 

      parts(:,rank)%ucOld = parts(:,rank)%element%uc 
      parts(:,rank)%vcOld = parts(:,rank)%element%vc 
      parts(:,rank)%wcOld = parts(:,rank)%element%wc 

   end subroutine StoreOldPartSol

   subroutine AveragePartSol()
      
      parts(:,rank)%element%uc = 0.5d0*( parts(:,rank)%element%uc + parts(:,rank)%ucOld ) 
      parts(:,rank)%element%vc = 0.5d0*( parts(:,rank)%element%vc + parts(:,rank)%vcOld ) 
      parts(:,rank)%element%wc = 0.5d0*( parts(:,rank)%element%wc + parts(:,rank)%wcOld )  

   end subroutine AveragePartSol

   subroutine UpdatePartLocCell   
      
      integer :: i,j,k,ipart
      real(8) :: xp,yp,zp

      do ipart = 1,num_part(rank)
         ! x direction 
         i  = parts(ipart,rank)%ic
         xp = parts(ipart,rank)%element%xc 
         if ( ( parts(ipart,rank)%element%uc > 0.d0 .and. &
                ABS(x(i)-xp) > (0.5d0*(x(i+1)+x(i+2))-x(i)) ) .or. & 
              ( parts(ipart,rank)%element%uc < 0.d0 .and. &
                ABS(x(i)-xp) > (x(i)-0.5d0*(x(i-1)+x(i-2))) )) then
            call pariserror("Particles move more than dx in dt!")
         else if ( parts(ipart,rank)%element%uc > 0.d0 .and. &
                    ABS(x(i)-xp) > ABS(x(i+1)-xp) ) then 
            i = i+1
         else if ( parts(ipart,rank)%element%uc < 0.d0 .and. &
                    ABS(x(i)-xp) > ABS(x(i-1)-xp) ) then 
            i = i-1
         end if ! parts(ipart,rank)%element%uc
         
         ! y direction 
         j  = parts(ipart,rank)%jc
         yp = parts(ipart,rank)%element%yc 
         if ( ( parts(ipart,rank)%element%vc > 0.d0 .and. &
                ABS(y(j)-yp) > (0.5d0*(y(j+1)+y(j+2))-y(j)) ) .or. & 
              ( parts(ipart,rank)%element%vc < 0.d0 .and. &
                ABS(y(j)-yp) > (y(j)-0.5d0*(y(j-1)+y(j-2))) )) then
            call pariserror("Particles move more than dy in dt!")
         else if ( parts(ipart,rank)%element%vc > 0.d0 .and. &
                    ABS(y(j)-yp) > ABS(y(j+1)-yp) ) then 
            j = j+1
         else if ( parts(ipart,rank)%element%vc < 0.d0 .and. &
                    ABS(y(j)-yp) > ABS(y(j-1)-yp) ) then 
            j = j-1
         end if ! parts(ipart,rank)%element%vc

         ! z direction 
         k  = parts(ipart,rank)%kc
         zp = parts(ipart,rank)%element%zc 
         if ( ( parts(ipart,rank)%element%wc > 0.d0 .and. &
                ABS(z(k)-zp) > (0.5d0*(z(k+1)+z(k+2))-z(k)) ) .or. & 
              ( parts(ipart,rank)%element%wc < 0.d0 .and. &
                ABS(z(k)-zp) > (z(k)-0.5d0*(z(k-1)+z(k-2))) )) then
             call pariserror("Particles move more than dz in dt!")
         else if ( parts(ipart,rank)%element%wc > 0.d0 .and. &
                    ABS(z(k)-zp) > ABS(z(k+1)-zp) ) then 
            k = k+1
         else if ( parts(ipart,rank)%element%wc < 0.d0 .and. &
                    ABS(z(k)-zp) > ABS(z(k-1)-zp) ) then 
            k = k-1
         end if ! parts(ipart,rank)%element%wc

         parts(ipart,rank)%ic = i 
         parts(ipart,rank)%jc = j 
         parts(ipart,rank)%kc = k 

      end do ! ipart
   end subroutine UpdatePartLocCell   
   
   subroutine CollectPartCrossBlocks
       
      integer :: ipart,ipart_cross,ipart1,i,j,k
      integer :: ranknew
      integer :: c1,c2,c3

      allocate( num_part_cross(0:nPdomain-1) )
      allocate( parts_cross_id     (max_num_part_cross,0:nPdomain-1) )
      allocate( parts_cross_newrank(max_num_part_cross,0:nPdomain-1) )
      num_part_cross(:) = 0
      parts_cross_id     (:,:) = CRAZY_INT 
      parts_cross_newrank(:,:) = CRAZY_INT 
      if ( num_part(rank) > 0 ) then 
         do ipart = 1,num_part(rank)
            i = parts(ipart,rank)%ic
            j = parts(ipart,rank)%jc
            k = parts(ipart,rank)%kc

            if ( vofbdry_cond(1) == 'periodic' ) then
               if ( i < Ng ) then 
                  i = i + Nx
               else if ( i > Ng+Nx ) then
                  i = i + Nx
               end if !i
            end if ! vofbrdy_cond(1)

            if ( vofbdry_cond(2) == 'periodic' ) then
               if ( j < Ng ) then 
                  j = j + Ny
               else if ( j > Ng+Ny ) then
                  j = j + Ny
               end if !i
            end if ! vofbrdy_cond(2)
            
            if ( vofbdry_cond(3) == 'periodic' ) then
               if ( k < Ng ) then 
                  k = k + Nz
               else if ( k > Ng+Nz ) then
                  k = k + Nz
               end if !i
            end if ! vofbrdy_cond(3)
            ! Note: here only collect and transfer particles which cross blocks 
            !        due to periodic BC, the location and cell information will 
            !        not be changed until SetPartBC is called

            if ( i > ie .or. j > je .or. k > ke .or. &
                 i < is .or. j < js .or. k < ks ) then
               c1 = (i-Ng-1)/Mx  
               c2 = (j-Ng-1)/My  
               c3 = (k-Ng-1)/Mz  
               ranknew = c1*npy*npz + c2*npz + c3
               if ( ranknew > nPdomain-1 .or. ranknew < 0 ) then
                  call pariserror("new rank of particle out of range!")
               else if ( ranknew /= rank ) then 
                  num_part_cross(rank)  = num_part_cross(rank) + 1
                  parts_cross_id     (num_part_cross(rank),rank) = ipart 
                  parts_cross_newrank(num_part_cross(rank),rank) = ranknew 
               end if ! ranknew
            end if ! i,j,k
         end do ! ipart
      end if ! num_part(rank)
   end subroutine CollectPartCrossBlocks

   subroutine TransferPartCrossBlocks

      include 'mpif.h'

      integer :: ipart,ipart_cross,ipart1,i,j,k
      integer :: ierr,irank
      integer :: ranknew
      integer :: req(4),sta(MPI_STATUS_SIZE,4),MPI_Comm,ireq
      integer :: MPI_particle_type, oldtypes(0:3), blockcounts(0:3), & 
                 offsets(0:3), intextent,r8extent
      integer :: maxnum_part_cross, MPI_int_row

      call MPI_ALLGATHER(num_part_cross(rank), 1, MPI_INTEGER, &
                         num_part_cross,    1, MPI_INTEGER, MPI_Comm_World, ierr)
      maxnum_part_cross = maxval(num_part_cross)
      if ( maxnum_part_cross  > 0 ) then 

         call MPI_TYPE_CONTIGUOUS (maxnum_part_cross, MPI_INTEGER, MPI_int_row, ierr)
         call MPI_TYPE_COMMIT(MPI_int_row, ierr)

         call MPI_ALLGATHER(parts_cross_id(1:maxnum_part_cross,rank), 1, MPI_int_row, &
                            parts_cross_id(1:maxnum_part_cross,:),    1, MPI_int_row, & 
                            MPI_Comm_World, ierr)
         call MPI_ALLGATHER(parts_cross_newrank(1:maxnum_part_cross,rank), 1, MPI_int_row, &
                            parts_cross_newrank(1:maxnum_part_cross,:),    1, MPI_int_row, & 
                            MPI_Comm_World, ierr)
      !  Setup MPI derived type for drop_merge_comm
      call MPI_TYPE_EXTENT(MPI_REAL8,   r8extent,  ierr) 
      call MPI_TYPE_EXTENT(MPI_INTEGER, intextent, ierr) 
      offsets    (0) = 0 
      oldtypes   (0) = MPI_REAL8 
      blockcounts(0) = 7 
      offsets    (1) = offsets(0) + blockcounts(0)*r8extent 
      oldtypes   (1) = MPI_INTEGER  
      blockcounts(1) = 1  
      offsets    (2) = offsets(1) + blockcounts(1)*intextent 
      oldtypes   (2) = MPI_REAL8  
      blockcounts(2) = 9
      offsets    (3) = offsets(2) + blockcounts(2)*r8extent 
      oldtypes   (3) = MPI_INTEGER  
      blockcounts(3) = 4  

      call MPI_TYPE_STRUCT(4, blockcounts, offsets, oldtypes, & 
                           MPI_particle_type, ierr) 
      call MPI_TYPE_COMMIT(MPI_particle_type, ierr)

      do irank = 0,nPdomain-1
         if ( num_part_cross(irank) > 0 ) then
            do ipart_cross = 1,num_part_cross(irank)
               ipart   = parts_cross_id     (ipart_cross,irank)
               ranknew = parts_cross_newrank(ipart_cross,irank)
               if ( rank == irank ) then 
                  call MPI_ISEND(parts(ipart,irank),1, MPI_particle_type, & 
                                 ranknew, 15, MPI_COMM_WORLD, req(1), ierr)
                  call MPI_WAIT(req(1),sta(:,1),ierr)
                  do ipart1 = ipart,num_part(irank)-1
                     parts(ipart1,irank) = parts(ipart1+1,irank)
                  end do ! ipart1
                  num_part(irank) = num_part(irank) - 1
               else if ( rank == ranknew ) then 
                  call MPI_IRECV(parts(num_part(ranknew)+1,ranknew),1,MPI_particle_type, & 
                                 irank, 15, MPI_COMM_WORLD, req(2), ierr)
                  call MPI_WAIT(req(2),sta(:,2),ierr)
                  num_part(ranknew) = num_part(ranknew) + 1
               end if ! rank 
            end do ! ipart_cross 
         end if ! num_part_cross(irank)
      end do ! irank
      call MPI_ALLGATHER(num_part(rank), 1, MPI_INTEGER, &
                         num_part(:)   , 1, MPI_INTEGER, MPI_Comm_World, ierr)
      end if ! maxnum_part_cross

      ! final
      deallocate(num_part_cross)
      deallocate(parts_cross_id)
      deallocate(parts_cross_newrank)

   end subroutine TransferPartCrossBlocks

   subroutine SetPartBC

      integer :: ipart

      if ( num_part(rank) > 0 ) then 
         do ipart = 1,num_part(rank)
            if ( parts(ipart,rank)%ic < Ng .or. parts(ipart,rank)%ic > Ng+Nx ) then 
               call ImposePartBC(ipart,rank,1)
            end if ! parts(ipart,rank)%ic

            if ( parts(ipart,rank)%jc < Ng .or. parts(ipart,rank)%jc > Ng+Ny ) then 
               call ImposePartBC(ipart,rank,2)
            end if ! parts(ipart,rank)%jc

            if ( parts(ipart,rank)%kc < Ng .or. parts(ipart,rank)%kc > Ng+Nz ) then 
               call ImposePartBC(ipart,rank,3)
            end if ! parts(ipart,rank)%kc
         end do ! ipart
      end if ! num_part(rank)

   end subroutine SetPartBC

   subroutine ImposePartBC(ipart,rank,d)
      integer, intent (in) :: ipart,rank,d

      if ( vofbdry_cond(d) == 'periodic' ) then
         call PartBC_periodic(ipart,rank,d)
      else 
         call pariserror("unknown particle bondary condition!")
      end if ! vofbdry_cond
   end subroutine ImposePartBC

   subroutine PartBC_periodic(ipart,rank,d)
      integer, intent (in) :: ipart,rank,d
      
      if ( d == 1 ) then 
         if ( parts(ipart,rank)%ic < Ng ) then 
            parts(ipart,rank)%ic = parts(ipart,rank)%ic + Nx
            parts(ipart,rank)%element%xc = parts(ipart,rank)%element%xc + xLength
         else if ( parts(ipart,rank)%ic > Ng+Nx ) then 
            parts(ipart,rank)%ic = parts(ipart,rank)%ic - Nx
            parts(ipart,rank)%element%xc = parts(ipart,rank)%element%xc - xLength
         end if ! parts(ipart,rank)%ic
      else if ( d == 2 ) then 
         if ( parts(ipart,rank)%jc < Ng ) then 
            parts(ipart,rank)%jc = parts(ipart,rank)%jc + Ny
            parts(ipart,rank)%element%yc = parts(ipart,rank)%element%yc + yLength
         else if ( parts(ipart,rank)%jc > Ng+Ny ) then 
            parts(ipart,rank)%jc = parts(ipart,rank)%jc - Ny
            parts(ipart,rank)%element%yc = parts(ipart,rank)%element%yc - yLength
         end if ! parts(ipart,rank)%jc
      else if ( d == 3 ) then 
         if ( parts(ipart,rank)%kc < Ng ) then 
            parts(ipart,rank)%kc = parts(ipart,rank)%kc + Nz
            parts(ipart,rank)%element%zc = parts(ipart,rank)%element%zc + zLength
         else if ( parts(ipart,rank)%kc > Ng+Nz ) then 
            parts(ipart,rank)%kc = parts(ipart,rank)%kc - Nz
            parts(ipart,rank)%element%zc = parts(ipart,rank)%element%zc - zLength
         end if ! parts(ipart,rank)%kc
      end if ! d

   end subroutine PartBC_periodic

   subroutine LinearIntrpl(x,x0,x1,f0,f1,f)
      real(8), intent (in) :: x,x0,x1,f0,f1
      real(8), intent(out) :: f      
      real(8) :: xl,xr

      xl = (x-x0)/(x1-x0)
      xr = 1.d0 - xl
      f  = f0*xr + f1*xl
   end subroutine LinearIntrpl

   subroutine BilinearIntrpl(x,y,x0,y0,x1,y1,f00,f01,f10,f11,f)
      real(8), intent (in) :: x,y,x0,y0,x1,y1,f00,f01,f10,f11
      real(8), intent(out) :: f      
      real(8) :: f0,f1

      call LinearIntrpl(x,x0,x1,f00,f10,f0)
      call LinearIntrpl(x,x0,x1,f01,f11,f1)
      call LinearIntrpl(y,y0,y1,f0 ,f1 ,f)
   end subroutine BilinearIntrpl

   subroutine TrilinearIntrpl(x,y,z,x0,y0,z0,x1,y1,z1,f000,f001,f010,f011,f100,f101,f110,f111,f)
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
         call pariserror("Wrong direction in velocity interploation!")
      end if ! dir 
   end subroutine TrilinearIntrplFluidVel

   subroutine StoreDiffusionTerms()
      sdu = du
      sdv = dv
      sdw = dw
   end subroutine StoreDiffusionTerms

   subroutine ComputeFluidAccel()

      integer :: i,j,k
      real(8) :: dpdx,dpdy,dpdz

      ! Compute substantial derivatives of fluid velocity (fluid acceleration)
      ! Note: only needed in fluid phase
      do i=is,ie; do j=js,je; do k=ks,ke
         if ( cvof (i,j,k) == 0.d0 ) then 
            dpdx = (p(i,j,k) - p(i-1,j,k))/(x(i)-x(i-1))
            dpdy = (p(i,j,k) - p(i,j-1,k))/(y(j)-y(j-1))
            dpdz = (p(i,j,k) - p(i,j,k-1))/(z(k)-z(k-1))
            
            sdu(i,j,k) = sdu(i,j,k)-dpdx/rho1
            sdv(i,j,k) = sdv(i,j,k)-dpdy/rho1
            sdw(i,j,k) = sdw(i,j,k)-dpdz/rho1
         end if ! cvof 
      end do; end do; end do ! i,j,k

   end subroutine ComputeFluidAccel

end module module_Lag_part
