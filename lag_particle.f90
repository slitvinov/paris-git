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
   integer,parameter :: max_num_drop = 10000
   integer,parameter :: maxnum_cell_drop = 500
   logical :: tag_drop_initialized = .false.
   integer :: totalnum_drop
   integer, dimension(:), allocatable :: num_drop

   type element 
      integer :: id 
      real(8) :: xc,yc,zc,uc,vc,wc,vol,tem
   end type element

   type drop
      type(element) :: element
      integer :: num_cell_drop
      integer :: cell_list(3,maxnum_cell_drop)
   end type drop
   type (drop), dimension(:), allocatable :: drops
   
   integer :: max_num_part = 10000
   logical :: LPP_initialized = .false.
   integer, dimension(:), allocatable :: num_part

   type particle
      type(element) :: element
   end type particle 
   type (particle), dimension(:), allocatable :: parts

   real(8) :: vol_cut
contains
!=================================================================================================
   subroutine initialize_LPP()
      allocate( parts(max_num_part) )
      allocate( num_part(0:nPdomain-1) )
      num_part = 0
   end subroutine initialize_LPP

   subroutine initialize_tag_drop()
      allocate( num_drop (0:nPdomain-1) )
      num_drop = 0
      allocate( tag_flag(imin:imax,jmin:jmax,kmin:kmax) )
      allocate( tag_id  (imin:imax,jmin:jmax,kmin:kmax) )
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
    
    real :: volscale,cvof_scaled

    if (.not. tag_drop_initialized) then 
      call initialize_tag_drop()
      tag_drop_initialized = .true.
    end if ! tag_drop_initialized

    tag_id  (:,:,:) = 0
    tag_flag(:,:,:) = 0
    current_id = 1

    drops(:)%element%vol = 0.d0
    drops(:)%element%xc = 0.d0
    drops(:)%element%yc = 0.d0
    drops(:)%element%zc = 0.d0
    drops(:)%element%uc = 0.d0
    drops(:)%element%vc = 0.d0
    drops(:)%element%wc = 0.d0
    drops(:)%num_cell_drop = 0
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
            drops(current_id)%element%id = current_id
            cvof_scaled = cvof(isq,jsq,ksq)*volscale
            drops(current_id)%element%vol = drops(current_id)%element%vol & 
                                          + cvof_scaled
            drops(current_id)%element%xc  = drops(current_id)%element%xc  & 
                                          + cvof_scaled*x(isq)
            drops(current_id)%element%yc  = drops(current_id)%element%yc  & 
                                          + cvof_scaled*y(jsq)
            drops(current_id)%element%zc  = drops(current_id)%element%zc  & 
                                          + cvof_scaled*z(ksq)
            drops(current_id)%element%uc  = drops(current_id)%element%uc  & 
                                          + cvof_scaled*u(isq,jsq,ksq)
            drops(current_id)%element%vc  = drops(current_id)%element%vc  & 
                                          + cvof_scaled*v(isq,jsq,ksq)
            drops(current_id)%element%wc  = drops(current_id)%element%wc  & 
                                          + cvof_scaled*w(isq,jsq,ksq)
            drops(current_id)%num_cell_drop = drops(current_id)%num_cell_drop + 1
            drops(current_id)%cell_list(1:3,drops(current_id)%num_cell_drop) = [isq,jsq,ksq]
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
    drops(:)%element%xc = drops(:)%element%xc/(drops(:)%element%vol+1.0d-60)
    drops(:)%element%yc = drops(:)%element%yc/(drops(:)%element%vol+1.0d-60)
    drops(:)%element%zc = drops(:)%element%zc/(drops(:)%element%vol+1.0d-60)
    drops(:)%element%uc = drops(:)%element%uc/(drops(:)%element%vol+1.0d-60)
    drops(:)%element%vc = drops(:)%element%vc/(drops(:)%element%vol+1.0d-60)
    drops(:)%element%wc = drops(:)%element%wc/(drops(:)%element%vol+1.0d-60)

    num_drop(rank) = current_id-1
    ! Broadcast num_drop(rank) to all processes
    if ( nPdomain > 1 ) then 
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
      integer :: ib,ipart
      real(8), parameter :: PI = 3.14159265359d0

      type(drop), dimension(NumBubble) :: drops_ex

      ! tag droplets and calculate drop properties
      call tag_drop()
      call tag_drop_all()

      OPEN(UNIT=90,FILE=TRIM(out_path)//'/VOF_before_'//TRIM(int2text(rank,padding))//'.dat')
      write(90,*) 'title= " VOF drops before conversion "'
      write(90,*) 'variables = "x", "y", "z", "c", "tag" '
      write(90,*) 'zone i=,',Nx/nPx+Ng*2, 'j=',Ny/nPy+Ng*2, 'k=',Nz/nPz+Ng*2,'f=point'
      do k = kmin,kmax
         do j=jmin,jmax
            do i=imin,imax 
               write(90,'(4(E15.8,1X),(I5))') x(i),y(j),z(k),cvof(i,j,k),tag_id(i,j,k)
            end do ! i
         end do ! j
      end do ! k
      CLOSE(90)

      OPEN(UNIT=91,FILE=TRIM(out_path)//'/dropvol-'//TRIM(int2text(rank,padding))//'.dat')
      OPEN(UNIT=92,FILE=TRIM(out_path)//'/dropvol_ex-'//TRIM(int2text(rank,padding))//'.dat')
      call QSort(drops,NumBubble)
      do ib = 1, NumBubble
         drops_ex(ib)%element%xc  = xc(ib)
         drops_ex(ib)%element%vol = 4.0d0*rad(ib)**3.d0*PI/3.d0
      end do ! i 
      call QSort(drops_ex,NumBubble)
      do ib = 1, NumBubble 
         write(91,*) drops   (ib)%element%xc, drops   (ib)%element%vol
         write(92,*) drops_ex(ib)%element%xc, drops_ex(ib)%element%vol
      end do ! i 
      CLOSE(91)
      CLOSE(92)

      ! convert droplets to particles
      vol_cut = 4.0d0*0.04d0**3.d0*PI/3.d0
      call convertDrop2Part()

      ! output droplets & particles
      OPEN(UNIT=93,FILE=TRIM(out_path)//'/VOF_after_'//TRIM(int2text(rank,padding))//'.dat')
      write(93,*) 'title= " VOF drops after conversion "'
      write(93,*) 'variables = "x", "y", "z", "c", "tag" '
      write(93,*) 'zone i=,',Nx/nPx+Ng*2, 'j=',Ny/nPy+Ng*2, 'k=',Nz/nPz+Ng*2,'f=point'
      do k = kmin,kmax
         do j=jmin,jmax
            do i=imin,imax
               write(93,'(4(E15.8,1X),(I5))') x(i),y(j),z(k),cvof(i,j,k),tag_id(i,j,k)
            end do ! i
         end do ! j
      end do ! k
      CLOSE(93)

!      call output_VOF(0,imin,imax,jmin,jmax,kmin,kmax)

      OPEN(UNIT=94,FILE=TRIM(out_path)//'/LPP_after_'//TRIM(int2text(rank,padding))//'.dat')
      write(94,*) 'title= " Lagrangian particles after conversion "'
      write(94,*) 'variables = "x", "y", "z"'
      do ipart = 1,num_part(rank)
         write(94,*) parts(ipart)%element%xc,parts(ipart)%element%yc,parts(ipart)%element%zc
      end do ! ipart
      CLOSE(94)

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

      convertDropFlag = .false.
      convertDoneFlag = .false.
      do idrop = 1,num_drop(rank)
         if ( drops(idrop)%element%vol < vol_cut ) then 
            convertDropFlag = .true.
            ! other criteria can be added here
         else 
            convertDropFlag = .false.
! TEMPORARY ! works only with the vol-cut criteria
            convertDoneflag = .true. 
! END TEMPORARY
         end if 

         if ( convertDoneFlag ) exit

         if ( convertDropFlag ) then
            ! transfer droplet properties to particle
            num_part(rank) = num_part(rank) + 1
            parts(num_part(rank))%element = drops(idrop)%element
            ! compute average fluid quantities

            ! remove droplet vof structure
            do ilist = 1,drops(idrop)%num_cell_drop
               cvof(drops(idrop)%cell_list(1,ilist), &
                    drops(idrop)%cell_list(2,ilist), &
                    drops(idrop)%cell_list(3,ilist)) = 0.0
            end do ! ilist
         end if !convertDropFlag
      end do ! idrop
           
   end subroutine convertDrop2Part

end module module_Lag_part
