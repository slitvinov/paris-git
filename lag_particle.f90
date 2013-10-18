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
   implicit none
   integer, dimension(:,:,:), allocatable :: tag_id
   integer, dimension(:,:,:), allocatable :: tag_flag
      ! 0 marked as untagged or reference fluid
      ! 1 marked as tagged droplet 
      ! 2 marked as S node
      ! 3 marked as C node
      ! 4 marked as reference fluid 
   integer :: max_num_drop = 10000
   logical :: tag_drop_initialized = .false.
   integer :: totalnum_drop
   integer, dimension(:), allocatable :: num_drop

   !   real(8), dimension(:,:), allocatable :: xdrop,udrop
!      ! (:,:) = (dimension,particle_id)
!   real(8), dimension(:,:), allocatable :: prop_drop
!      ! (:,:) = (property-index,particle_id)
!      ! property index: 1, temperature 2, volume
!   integer, dimension(:),   allocatable :: droptag
!      ! (:) = (particle_id)
   type drop
      integer :: tag
      real(8) :: xc,yc,zc,uc,vc,wc,vol,tem
   end type drop
   type (drop), dimension(:), allocatable :: drops

   
   integer :: max_num_part = 10000
   logical :: LP_initialized = .false.
   integer :: num_part 
   real(8), dimension(:,:), allocatable :: xpart,upart
      ! (:,:) = (dimension,particle_id)
   real(8), dimension(:,:), allocatable :: prop_part
      ! (:,:) = (property-index,particle_id)
      ! property index: 1, temperature 2, volume

contains
!=================================================================================================
   subroutine initialize_Lag_part()
      allocate( xpart    (3,max_num_part) )
      allocate( upart    (3,max_num_part) )
      allocate( prop_part(2,max_num_part) )
   end subroutine initialize_Lag_part

   subroutine initialize_tag_drop()
      allocate( num_drop (0:nPdomain-1) )
      num_drop = 0
      allocate( tag_flag(imin:imax,jmin:jmax,kmin:kmax) )
      allocate( tag_id  (imin:imax,jmin:jmax,kmin:kmax) )
!      allocate( xdrop    (3,max_num_drop) )
!      allocate( udrop    (3,max_num_drop) )
!      allocate( prop_drop(2,max_num_drop) )
!      allocate( droptag  (  max_num_drop) )
      allocate( drops(max_num_drop) )
   end subroutine initialize_tag_drop

! 
! Purpose: tag isolated droplets and record properties
! 
  subroutine tag_drop()
    include 'mpif.h'
    integer :: i,j,k, i0,j0,k0
    integer :: isq,jsq,ksq
    integer :: current_id
    integer :: s_queue(Nx*Ny,3),c_queue(Nx*Ny,3) ! record i,j,k
    integer :: ns_queue,is_queue,nc_queue
    integer :: ierr,irank
    
    real :: volscale

    if (.not. tag_drop_initialized) then 
      call initialize_tag_drop()
      tag_drop_initialized = .true.
    end if ! tag_drop_initialized

    tag_id  (:,:,:) = 0
    tag_flag(:,:,:) = 0
    current_id = 1

    drops(:)%vol = 0.d0
    volscale = 1.d0/(real(Nx)*real(Ny)*real(Nz))

    ns_queue = 0
    s_queue(:,:) = 0
    num_drop(rank) = 0
    WRITE(*,*) ' start tagging ... ',rank
    do i=is,ie; do j=js,je; do k=ks,ke
      if ( cvof(i,j,k) > 0.d0 .and. tag_flag(i,j,k) == 0 ) then 
        tag_id  (i,j,k) = current_id
        tag_flag(i,j,k) = 2 ! mark as S node
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
              if ( cvof(isq+i0,jsq+j0,ksq+k0) > 0.d0 .and. tag_flag(isq+i0,jsq+j0,ksq+k0) == 0 ) then 
                tag_id  (isq+i0,jsq+j0,ksq+k0) = current_id  ! tag node with id
                tag_flag(isq+i0,jsq+j0,ksq+k0) = 3  ! mark as C node
                ! put current node into C queue
                nc_queue = nc_queue + 1
                c_queue(nc_queue,1) = isq+i0
                c_queue(nc_queue,2) = jsq+j0
                c_queue(nc_queue,3) = ksq+k0
              end if 
            enddo;enddo;enddo ! i0,j0,k0
            tag_flag(isq,jsq,ksq) = 1 !unmark S node and marked as tagged
            ! perform droplet calculation
            drops(current_id)%tag = current_id
            drops(current_id)%vol = drops(current_id)%vol & 
                                    + cvof(isq,jsq,ksq)*volscale
            drops(current_id)%xc  = drops(current_id)%xc  & 
                                    + cvof(isq,jsq,ksq)*volscale*x(isq)
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
    drops(:)%xc = drops(:)%xc/(drops(:)%vol+1.0d-60)

    ! Broadcast num_drop(rank) to all processes
    if ( nPdomain > 1 ) then 
      num_drop(rank) = current_id-1
      do irank = 0, nPdomain-1
         call MPI_BCAST(num_drop(irank), 1, MPI_INTEGER, &
                        irank, MPI_Comm_Cart, ierr)
      end do ! irank
      totalnum_drop = sum(num_drop,1)
    end if 

  end subroutine tag_drop

   subroutine tag_drop_all()
      include 'mpif.h'
      integer :: i,j,k
      integer :: ierr,numaccu

      ! update tag_id from local to global id
      !    Note: no need to change domain 0
      if ( nPdomain > 1 .and. rank > 0 ) then 
         numaccu = sum(num_drop(0:rank-1),1)
         do i=is,ie; do j=js,je; do k=ks,ke
            if ( tag_id(i,j,k) > 0 ) tag_id(i,j,k)=tag_id(i,j,k)+numaccu
         end do; end do; end do
      end if ! rank
      ! 3, merge droplet pieces crossing block boundary

! DEBUG
      WRITE(*,*) ' ********************************'
      WRITE(*,*) rank,numaccu
      write(*,*) num_drop
      write(*,*) totalnum_drop
! END DEBUG
   end subroutine tag_drop_all

! ==============================================
! output tag of droplets
! ==============================================
   subroutine output_tag()
      integer :: i,j,k

      call tag_drop()
      call tag_drop_all()

      OPEN(UNIT=89,FILE=TRIM(out_path)//'/tag-'//TRIM(int2text(rank,padding))//'.txt')
      OPEN(UNIT=90,FILE=TRIM(out_path)//'/tag-tecplot'//TRIM(int2text(rank,padding))//'.dat')

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
      do j=jmax,jmin,-1
         write(89,'(36(I5,1X))') tag_id(imin:imax,j,Nz/4+Ng)
      end do ! j
      CLOSE(89)
      CLOSE(90)
   end subroutine output_tag

! ==============================================
! output droplets & particles 
! ==============================================
   subroutine output_DP()
      integer :: i,j,k
      integer :: ib
      real(8), parameter :: PI = 3.14159265359d0

      type(drop), dimension(NumBubble) :: drops_ex

      call tag_drop()
      call tag_drop_all()

      OPEN(UNIT=90,FILE=TRIM(out_path)//'/dropvof_before_'//TRIM(int2text(rank,padding))//'.dat')
      write(90,*) 'title= " VOF drops before conversion "'
      write(90,*) 'variables = "x", "y", "z", "c", "tag" '
      write(90,*) 'zone i=,',Nx/nPx+Ng*2, 'j=',Ny/nPy+Ng*2, 'k=',Nz/nPz+Ng*2,'f=point'
      do k = kmin,kmax
         do j=jmin,jmax
            do i=imin,imax 
               write(90,'(3(I5,1X),(E15.8),(I5))') i,j,k,cvof(i,j,k),tag_id(i,j,k)
            end do ! i
         end do ! j
      end do ! k
      CLOSE(90)

      ! convert droplets to particles
      OPEN(UNIT=91,FILE=TRIM(out_path)//'/dropvol-'//TRIM(int2text(rank,padding))//'.dat')
      OPEN(UNIT=92,FILE=TRIM(out_path)//'/dropvol_ex-'//TRIM(int2text(rank,padding))//'.dat')
      call QSort(drops,NumBubble)
      do ib = 1, NumBubble
         drops_ex(ib)%xc  = xc(ib)
         drops_ex(ib)%vol = 4.0d0*rad(ib)**3.d0*PI/3.d0
      end do ! i 
      call QSort(drops_ex,NumBubble)
      do ib = 1, NumBubble 
         write(91,*) drops   (ib)%xc, drops   (ib)%vol
         write(92,*) drops_ex(ib)%xc, drops_ex(ib)%vol
      end do ! i 
      CLOSE(91)
      CLOSE(92)

      ! output droplets & particles

   end subroutine output_DP

! ===============================================
! Testing section
! ===============================================
   subroutine test_Lag_part()
      implicit none
      include 'mpif.h'
      integer :: ierr
                     
      if(test_tag) then
         call output_tag()
      else if ( test_D2P ) then
         call output_DP()
      end if
                                                                             
      ! Exit MPI gracefully
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_finalize(ierr)
      stop
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
      pivot = A(int(random*real(nA-1))+1)%vol

      left = 0
      right = nA + 1

      do while (left < right)
         right = right - 1
         do while (A(right)%vol > pivot)
            right = right - 1
         end do
         left = left + 1
         do while (A(left)%vol < pivot)
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
 
end module module_Lag_part
