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
   integer,parameter :: max_num_drop = 10000
   integer,parameter :: maxnum_cell_drop = 2000
   integer,parameter :: maxnum_diff_tag  = 7   ! ignore cases droplet spread over more than 1 block
   integer,parameter :: ntimesteptag = 1
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
      integer :: id 
      real(8) :: xc,yc,zc,uc,vc,wc,vol
   end type element
   type (element), dimension(:,:), allocatable :: element_stat

   type drop
      type(element) :: element
      integer :: num_cell_drop
      integer :: cell_list(3,maxnum_cell_drop)
   end type drop
   type (drop), dimension(:,:), allocatable :: drops
   
   type drop_merge
      type(element) :: element
      integer :: num_cell_drop
      integer :: cell_list(3,maxnum_cell_drop)
      integer :: num_gcell
      integer :: gcell_list(3,maxnum_cell_drop)
      integer :: num_diff_tag
      integer :: diff_tag_list(maxnum_diff_tag)
      integer :: flag_center_mass
   end type drop_merge
   type (drop_merge), dimension(:,:), allocatable :: drops_merge

   type drop_merge_comm
      real(8) :: xc,yc,zc,uc,vc,wc,vol
      integer :: id
      integer :: num_diff_tag 
      integer :: diff_tag_list(maxnum_diff_tag)
      integer :: flag_center_mass
   end type drop_merge_comm
   type (drop_merge_comm), dimension(:,:), allocatable :: drops_merge_comm

   integer :: max_num_part = 10000
   logical :: LPP_initialized = .false.
   integer, dimension(:), allocatable :: num_part

   type particle
      type(element) :: element
      real(8) :: fx,fy,fz
      real(8) :: xcOld,ycOld,zcOld,ucOld,vcOld,wcOld 
   end type particle 
   type (particle), dimension(:,:), allocatable :: parts

   real(8) :: vol_cut, y_cut   ! Note: convert to input parameter later
   real(8), parameter :: PI = 3.14159265359d0

   logical :: DoDropStatistics = .false.
   integer :: dragmodel = 1

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

      num_tag  = 0
      num_part = 0
      num_drop = 0
      num_drop_merge = 0
      num_element = 0
   end subroutine initialize_LPP

   subroutine ReadLPPParameters

      include 'mpif.h'
      integer ierr,in
      logical file_is_there
      namelist /lppparameters/ DoDropStatistics, dragmodel 
      in=32

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

         call convertDrop2Part()
!         call convertPart2Drop()   ! TBA
      end if ! tswap
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

    logical :: merge_drop

    if (.not. LPP_initialized) then 
      call initialize_LPP()
      LPP_initialized = .true.
    end if ! LPP_initialized

    tag_id  (:,:,:) = 0
    tag_flag(:,:,:) = 0
    current_id = 1

    drops(:,rank)%element%vol = 0.d0
    drops(:,rank)%element%xc  = 0.d0
    drops(:,rank)%element%yc  = 0.d0
    drops(:,rank)%element%zc  = 0.d0
    drops(:,rank)%element%uc  = 0.d0
    drops(:,rank)%element%vc  = 0.d0
    drops(:,rank)%element%wc  = 0.d0
    drops(:,rank)%num_cell_drop = 0

    drops_merge(:,rank)%element%vol = 0.d0
    drops_merge(:,rank)%element%xc  = 0.d0
    drops_merge(:,rank)%element%yc  = 0.d0
    drops_merge(:,rank)%element%zc  = 0.d0
    drops_merge(:,rank)%element%uc  = 0.d0
    drops_merge(:,rank)%element%vc  = 0.d0
    drops_merge(:,rank)%element%wc  = 0.d0
    drops_merge(:,rank)%num_cell_drop = 0
    drops_merge(:,rank)%num_gcell = 0

    ns_queue = 0
    s_queue(:,:) = 0
    num_drop(:) = 0
    num_drop_merge(:) = 0
    do i=is,ie; do j=js,je; do k=ks,ke
      merge_drop = .false.
      if ( cvof(i,j,k) > 0.d0 .and. tag_flag(i,j,k) == 0 ) then 
        tag_id  (i,j,k) = current_id
        tag_flag(i,j,k) = 2 ! mark as S node
        num_drop(rank) = num_drop(rank) + 1
        ! put the present node into S queue
        ns_queue = ns_queue + 1 
        s_queue(ns_queue,1) = i
        s_queue(ns_queue,2) = j
        s_queue(ns_queue,3) = k
  
        do while ( ns_queue > 0 )
          nc_queue = 0 
          c_queue(:,:) = 0
          do is_queue=1,ns_queue
            isq = s_queue(is_queue,1)
            jsq = s_queue(is_queue,2)
            ksq = s_queue(is_queue,3)
       
            do i0=-1,1; do j0=-1,1; do k0=-1,1
              if ( cvof(isq+i0,jsq+j0,ksq+k0)      > 0.d0 .and. & 
                   tag_flag(isq+i0,jsq+j0,ksq+k0) == 0 ) then 
                  if (  isq+i0 >= is .and. isq+i0 <= ie .and. &     ! internal cells 
                        jsq+j0 >= js .and. jsq+j0 <= je .and. & 
                        ksq+k0 >= ks .and. ksq+k0 <= ke ) then
                     tag_id  (isq+i0,jsq+j0,ksq+k0) = current_id  ! tag node with id
                     tag_flag(isq+i0,jsq+j0,ksq+k0) = 3  ! mark as C node
                     ! put current node into C queue
                     nc_queue = nc_queue + 1
                     c_queue(nc_queue,1) = isq+i0
                     c_queue(nc_queue,2) = jsq+j0
                     c_queue(nc_queue,3) = ksq+k0
                  else                                               ! ghost cells
                     ! XXX Note: need to distinguish block ghost cells and domain
                     ! ghost cells!
                     tag_flag(isq+i0,jsq+j0,ksq+k0) = 5
                     if ( merge_drop .eqv. .false.) then 
                        merge_drop = .true.
                        num_drop      (rank) = num_drop      (rank) - 1
                        num_drop_merge(rank) = num_drop_merge(rank) + 1
                     end if ! merge_drop
                     drops_merge(num_drop_merge(rank),rank)%num_gcell = drops_merge(num_drop_merge(rank),rank)%num_gcell + 1
                     drops_merge(num_drop_merge(rank),rank)%gcell_list(1,drops_merge(num_drop_merge(rank),rank)%num_gcell) = isq+i0
                     drops_merge(num_drop_merge(rank),rank)%gcell_list(2,drops_merge(num_drop_merge(rank),rank)%num_gcell) = jsq+j0
                     drops_merge(num_drop_merge(rank),rank)%gcell_list(3,drops_merge(num_drop_merge(rank),rank)%num_gcell) = ksq+k0
                  end if ! isq+i0 > is
              end if 
            enddo;enddo;enddo ! i0,j0,k0
            tag_flag(isq,jsq,ksq) = 1 !unmark S node and marked as tagged
            
            ! perform droplet calculation
            volcell = dx(isq)*dy(jsq)*dz(ksq) 
            cvof_scaled = cvof(isq,jsq,ksq)*volcell
            if ( merge_drop ) then 
               drops_merge(num_drop_merge(rank),rank)%element%id = current_id
               drops_merge(num_drop_merge(rank),rank)%element%vol = &
                  drops_merge(num_drop_merge(rank),rank)%element%vol + cvof_scaled
               drops_merge(num_drop_merge(rank),rank)%element%xc  = & 
                  drops_merge(num_drop_merge(rank),rank)%element%xc  + cvof_scaled*x(isq)
               drops_merge(num_drop_merge(rank),rank)%element%yc  = & 
                  drops_merge(num_drop_merge(rank),rank)%element%yc  + cvof_scaled*y(jsq)
               drops_merge(num_drop_merge(rank),rank)%element%zc  = & 
                  drops_merge(num_drop_merge(rank),rank)%element%zc  + cvof_scaled*z(ksq)
               drops_merge(num_drop_merge(rank),rank)%element%uc  = &
                  drops_merge(num_drop_merge(rank),rank)%element%uc  + cvof_scaled*u(isq,jsq,ksq)
               drops_merge(num_drop_merge(rank),rank)%element%vc  = &
                  drops_merge(num_drop_merge(rank),rank)%element%vc  + cvof_scaled*v(isq,jsq,ksq)
               drops_merge(num_drop_merge(rank),rank)%element%wc  = &
                  drops_merge(num_drop_merge(rank),rank)%element%wc  + cvof_scaled*w(isq,jsq,ksq)
               drops_merge(num_drop_merge(rank),rank)%num_cell_drop = &
                  drops_merge(num_drop_merge(rank),rank)%num_cell_drop + 1
               drops_merge(num_drop_merge(rank),rank) &
                  %cell_list(1:3,drops_merge(num_drop_merge(rank),rank)%num_cell_drop) = [isq,jsq,ksq]
            else 
               drops(num_drop(rank),rank)%element%id = current_id
               drops(num_drop(rank),rank)%element%vol = drops(num_drop(rank),rank)%element%vol & 
                                                      + cvof_scaled
               drops(num_drop(rank),rank)%element%xc  = drops(num_drop(rank),rank)%element%xc  & 
                                                      + cvof_scaled*x(isq)
               drops(num_drop(rank),rank)%element%yc  = drops(num_drop(rank),rank)%element%yc  & 
                                                      + cvof_scaled*y(jsq)
               drops(num_drop(rank),rank)%element%zc  = drops(num_drop(rank),rank)%element%zc  & 
                                                      + cvof_scaled*z(ksq)
               drops(num_drop(rank),rank)%element%uc  = drops(num_drop(rank),rank)%element%uc  & 
                                                      + cvof_scaled*u(isq,jsq,ksq)
               drops(num_drop(rank),rank)%element%vc  = drops(num_drop(rank),rank)%element%vc  & 
                                                      + cvof_scaled*v(isq,jsq,ksq)
               drops(num_drop(rank),rank)%element%wc  = drops(num_drop(rank),rank)%element%wc  & 
                                                      + cvof_scaled*w(isq,jsq,ksq)
               drops(num_drop(rank),rank)%num_cell_drop = drops(num_drop(rank),rank)%num_cell_drop + 1
               drops(num_drop(rank),rank)%cell_list(1:3,drops(num_drop(rank),rank)%num_cell_drop) = [isq,jsq,ksq]
            end if ! merge_drop
          end do ! is_queue
          ! mark all C nodes as S nodes
          if ( nc_queue >= 0 ) then 
            s_queue(:,:) = c_queue(:,:)   ! mark all C nodes as S nodes
            ns_queue = nc_queue
          end if ! nc_queue
        end do ! ns_queue>0
        current_id = current_id+1
      else if ( cvof(i,j,k) == 0.d0 .and. tag_flag(i,j,k) == 0 ) then 
         tag_flag(i,j,k) = 4
      end if ! cvof(i,j,k)
    enddo; enddo; enddo

    if ( num_drop(rank) > 0 ) then 
      drops(:,rank)%element%xc = drops(:,rank)%element%xc/(drops(:,rank)%element%vol+1.0d-60)
      drops(:,rank)%element%yc = drops(:,rank)%element%yc/(drops(:,rank)%element%vol+1.0d-60)
      drops(:,rank)%element%zc = drops(:,rank)%element%zc/(drops(:,rank)%element%vol+1.0d-60)
      drops(:,rank)%element%uc = drops(:,rank)%element%uc/(drops(:,rank)%element%vol+1.0d-60)
      drops(:,rank)%element%vc = drops(:,rank)%element%vc/(drops(:,rank)%element%vol+1.0d-60)
      drops(:,rank)%element%wc = drops(:,rank)%element%wc/(drops(:,rank)%element%vol+1.0d-60)
    end if ! num_drop(rank)

    if ( num_drop_merge(rank) > 0 ) then 
      drops_merge(:,rank)%element%xc = drops_merge(:,rank)%element%xc/(drops_merge(:,rank)%element%vol+1.0d-60)
      drops_merge(:,rank)%element%yc = drops_merge(:,rank)%element%yc/(drops_merge(:,rank)%element%vol+1.0d-60)
      drops_merge(:,rank)%element%zc = drops_merge(:,rank)%element%zc/(drops_merge(:,rank)%element%vol+1.0d-60)
      drops_merge(:,rank)%element%uc = drops_merge(:,rank)%element%uc/(drops_merge(:,rank)%element%vol+1.0d-60)
      drops_merge(:,rank)%element%vc = drops_merge(:,rank)%element%vc/(drops_merge(:,rank)%element%vol+1.0d-60)
      drops_merge(:,rank)%element%wc = drops_merge(:,rank)%element%wc/(drops_merge(:,rank)%element%vol+1.0d-60)
    end if ! num_drop_merge(rank)

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
                         num_drop            , 1, MPI_INTEGER, MPI_Comm_World, ierr)
      call MPI_ALLGATHER(num_drop_merge(rank), 1, MPI_INTEGER, &
                         num_drop_merge      , 1, MPI_INTEGER, MPI_Comm_World, ierr)
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
            = tag_id(drops_merge(idrop,rank)%gcell_list(1,1), &
                     drops_merge(idrop,rank)%gcell_list(2,1), &
                     drops_merge(idrop,rank)%gcell_list(3,1))
         if ( drops_merge(idrop,rank)%num_gcell > 1 ) then  
            do iCell = 2,drops_merge(idrop,rank)%num_gcell
               do idiff_tag = 1, drops_merge(idrop,rank)%num_diff_tag
                  if ( tag_id(drops_merge(idrop,rank)%gcell_list(1,iCell), & 
                              drops_merge(idrop,rank)%gcell_list(2,iCell), & 
                              drops_merge(idrop,rank)%gcell_list(3,iCell)) &
                    == drops_merge(idrop,rank)%diff_tag_list(idiff_tag)) exit
               end do ! idiff_tag
               if ( idiff_tag == drops_merge(idrop,rank)%num_diff_tag + 1 ) then 
                  drops_merge(idrop,rank)%num_diff_tag = &
                  drops_merge(idrop,rank)%num_diff_tag + 1
                  drops_merge(idrop,rank)%diff_tag_list(drops_merge(idrop,rank)%num_diff_tag) &
                  = tag_id(drops_merge(idrop,rank)%gcell_list(1,iCell), &
                           drops_merge(idrop,rank)%gcell_list(2,iCell), &
                           drops_merge(idrop,rank)%gcell_list(3,iCell))
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
                 offsets(0:1), extent,r8extent 
      integer :: maxnum_tag 
      integer :: irank, idrop, ielement

! DEBUG
      call pariserror("drop-statisitics called!")
! END DEBUG

      maxnum_tag = maxval(num_tag)
      allocate( element_stat(maxnum_tag,0:nPdomain-1) )

      !  Setup MPI derived type for drop_merge_comm
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
      num_element(rank) = 0 
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
 
      ! Send all droplet pieces to rank 0
      if ( rank > 0 ) then
         if ( num_element(rank) > 0 ) then 
            call MPI_ISEND(element_stat(1:num_element(rank),rank),num_element(rank), & 
                           MPI_element_type, 0, 14, MPI_COMM_WORLD, req(1), ierr)
            call MPI_WAIT(req(1),sta(:,1),ierr)
         end if ! num_element(rank)
      else
         do irank = 1,nPdomain-1
            if ( num_element(irank) > 0 ) then 
               call MPI_IRECV(element_stat(1:num_element(irank),irank),num_element(irank), & 
                              MPI_element_type, irank, 14, MPI_COMM_WORLD, req(2), ierr)
               call MPI_WAIT(req(2),sta(:,2),ierr)
            end if ! num_drop_merge(irank)
         end do ! irank
      end if ! rank

      ! Statisitcs or other operation at rank 0 
      if ( rank == 0 ) then 
         do irank = 1,nPdomain-1
            do ielement = 1, num_element(irank) 
               write(200+ielement,*) tswap,  element_stat(ielement,irank)%xc, & 
                                             element_stat(ielement,irank)%yc, &
                                             element_stat(ielement,irank)%yc, &
                                             element_stat(ielement,irank)%uc, &
                                             element_stat(ielement,irank)%vc, &
                                             element_stat(ielement,irank)%wc, &
                                             element_stat(ielement,irank)%vol
            end do ! ielement
         end do ! irank
      end if ! rank

      ! finalize
      call MPI_TYPE_FREE(MPI_element_type, ierr)

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

      do idrop = 1,num_drop(rank)
         drops(idrop,rank)%element%id = tag_id( drops(idrop,rank)%cell_list(1,1), &
                                                drops(idrop,rank)%cell_list(2,1), &
                                                drops(idrop,rank)%cell_list(3,1) )
         tag_dropid   (drops(idrop,rank)%element%id) = idrop
         tag_rank     (drops(idrop,rank)%element%id) = rank 
         tag_mergeflag(drops(idrop,rank)%element%id) = 0 
      end do ! idrop
      do idrop = 1,num_drop_merge(rank)
         drops_merge(idrop,rank)%element%id = tag_id( drops_merge(idrop,rank)%cell_list(1,1), &
                                                      drops_merge(idrop,rank)%cell_list(2,1), &
                                                      drops_merge(idrop,rank)%cell_list(3,1) )
         tag_dropid   (drops_merge(idrop,rank)%element%id) = idrop
         tag_rank     (drops_merge(idrop,rank)%element%id) = rank 
         tag_mergeflag(drops_merge(idrop,rank)%element%id) = 1 
      end do ! idrop

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
      if (nPdomain > 1 ) call tag_drop_all()

!      OPEN(UNIT=89,FILE=TRIM(out_path)//'/tag-'//TRIM(int2text(rank,padding))//'.txt')
      OPEN(UNIT=90,FILE=TRIM(out_path)//'/tag-tecplot'//TRIM(int2text(rank,padding))//'.dat')

!      do j=jmax,jmin,-1
!         write(89,'(36(I5,1X))') tag_id(imin:imax,j,Nz/4+Ng)
!      end do ! j
!      CLOSE(89)

!      write(90,*) 'title= " 3d tag "'
!      write(90,*) 'variables = "x", "y", "z", "tag", "c" '
!      write(90,*) 'zone i=,',Nx/nPx+Ng*2, 'j=',Ny/nPy+Ng*2, 'k=',Nz/nPz+Ng*2,'f=point'
!      do k = kmin,kmax
!         do j=jmin,jmax
!            do i=imin,imax 
!               write(90,'(4(I5,1X),(E15.8))') i,j,k,tag_id(i,j,k),cvof(i,j,k)
!            end do ! i
!         end do ! j
!      end do ! k

!      write(90,*) 'title= " 2d tag "'
!      write(90,*) 'variables = "x", "y", "tag", "c" '
!      write(90,*) 'zone i=,',Nx/nPx+Ng*2, 'j=',Ny/nPy+Ng*2, 'f=point'
!      k = Nz/4+Ng
!         do j=jmin,jmax
!            do i=imin,imax 
!               write(90,'(2(E15.8,1X),(I5,1X),(E15.8))') x(i),y(j),tag_id(i,j,k),cvof(i,j,k)
!            end do ! i
!         end do ! j
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
   subroutine output_DP()
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
      call convertDrop2Part()

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
   subroutine test_Lag_part()
      implicit none
      include 'mpif.h'
!      integer :: ierr
                     
      if(test_tag) then
         call output_tag()
      else if ( test_D2P ) then
         call output_DP()
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

   subroutine convertDrop2Part()

      integer :: idrop,ilist
      logical :: convertDropFlag,convertDoneFlag

      if (.not. LPP_initialized) then 
         call initialize_LPP()
         LPP_initialized = .true.
      end if ! LPP_initialized

      ! Define criteria for conversion 
      if ( test_D2P ) then 
         vol_cut = 4.0d0*0.6d0**3.d0*PI/3.d0
      else if ( test_bubbles ) then 
         vol_cut = 4.0d0*0.05d-1**3.d0*PI/3.d0
         y_cut   = 0.3d-1
      end if ! test_type

      do idrop = 1,num_drop(rank)
         convertDropFlag = .false.
         if ( test_D2P ) then 
            if ( drops(idrop,rank)%element%vol < vol_cut ) then 
               convertDropFlag = .true.
            end if
         else if ( test_bubbles ) then 
            if ( drops(idrop,rank)%element%vol < vol_cut .and. &
                 drops(idrop,rank)%element%yc  > y_cut ) then 
               convertDropFlag = .true.
            end if
         end if ! test_D2P

         if ( convertDropFlag ) then
            ! transfer droplet properties to particle
            num_part(rank) = num_part(rank) + 1
            parts(num_part(rank),rank)%element = drops(idrop,rank)%element

            ! compute average fluid quantities

            ! remove droplet vof structure
            do ilist = 1,drops(idrop,rank)%num_cell_drop
               cvof(drops(idrop,rank)%cell_list(1,ilist), &
                    drops(idrop,rank)%cell_list(2,ilist), &
                    drops(idrop,rank)%cell_list(3,ilist)) = 0.0
            end do ! ilist
         end if !convertDropFlag
      end do ! idrop
          
      do idrop = 1,num_drop_merge(rank)
         if ( drops_merge(idrop,rank)%element%vol < vol_cut ) then 
            convertDropFlag = .true.
            ! other criteria can be added here
         else 
            convertDropFlag = .false.
         end if 

         if ( convertDropFlag ) then
            ! compute average fluid quantities

            ! remove droplet vof structure
            do ilist = 1,drops_merge(idrop,rank)%num_cell_drop
               cvof(drops_merge(idrop,rank)%cell_list(1,ilist), &
                    drops_merge(idrop,rank)%cell_list(2,ilist), &
                    drops_merge(idrop,rank)%cell_list(3,ilist)) = 0.0
            end do ! ilist
      
            ! transfer droplet properties to particle if center of mass located
            ! in this droplet piece
            if ( drops_merge(idrop,rank)%flag_center_mass == 1 ) then 
               num_part(rank) = num_part(rank) + 1
               parts(num_part(rank),rank)%element = drops_merge(idrop,rank)%element
            end if ! flag_center_mass 

         end if !convertDropFlag
      end do ! idrop

      ! transfer merged droplet to particles in rank 0 

   end subroutine convertDrop2Part

   subroutine ComputePartForce(tswap)

      integer, intent(in) :: tswap

      integer, parameter :: drag_model_Stokes = 1
      integer, parameter :: drag_model_SN = 2     ! Schiller & Nauman
      integer, parameter :: drag_model_CG = 3     ! Clift & Gauvin

      real(8) :: relvel(4), partforce(3)
      real(8) :: dp, Rep, muf, phi, rhof, rhop, taup, up,vp,wp, uf,vf,wf

      integer :: ipart
      if ( num_part(rank) > 0 ) then
         do ipart = 1,num_part(rank)
            up = parts(ipart,rank)%element%uc
            vp = parts(ipart,rank)%element%vc
            wp = parts(ipart,rank)%element%wc

            call GetFluidProp(parts(ipart,rank)%element%xc, &
                              parts(ipart,rank)%element%yc, &
                              parts(ipart,rank)%element%zc, & 
                              uf,vf,wf)

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

            partforce(1) = relvel(1)/taup*phi + (1.d0-rhof/rhop)*Gx   
            partforce(2) = relvel(2)/taup*phi + (1.d0-rhof/rhop)*Gy  
            partforce(3) = relvel(3)/taup*phi + (1.d0-rhof/rhop)*Gz  
            ! XXX: include other forces later

            parts(ipart,rank)%fx = partforce(1)
            parts(ipart,rank)%fy = partforce(2)
            parts(ipart,rank)%fz = partforce(3)
         end do ! ipart
      end if ! num_part(rank) 
   end subroutine ComputePartForce

   subroutine GetFluidProp(xp,yp,zp,uf,vf,wf) 
      real(8), intent(in) :: xp,yp,zp 
      real(8), intent(out) :: uf,vf,wf

! TEMPORARY
      uf = 0.d0 
      vf = 0.d0
      wf = 0.d0
! END TEMPORARY
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
! DEBUG
         write(101,*) tswap, & 
         parts(1,rank)%element%xc,parts(1,rank)%element%yc,parts(1,rank)%element%zc, &
         parts(1,rank)%element%uc,parts(1,rank)%element%vc,parts(1,rank)%element%wc, & 
         parts(1,rank)%fx,parts(1,rank)%fy,parts(1,rank)%fz
! END DEBUG
      end if ! num_part(rank)
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

end module module_Lag_part
