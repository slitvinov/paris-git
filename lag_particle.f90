!=================================================================================================
!=================================================================================================
! Paris-0.1
! Extended from Code: FTC3D2011 (Front Tracking Code for 3D simulations)
! and Surfer. 
! 
! Authors:
! 
!   Sadegh Dabiri (sdabiri@gmail.com), Gretar Tryggvason
!   Stephane Zaleski (zaleski@dalembert.upmc.fr) and Yue Ling (ling.stanley@gmail.com)
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
   implicit none
   integer, dimension(:,:,:), allocatable :: tag_id
   integer, dimension(:,:,:), allocatable :: tag_flag
      ! 0 marked as untagged or reference fluid
      ! 1 marked as tagged droplet 
      ! 2 marked as S node
      ! 3 marked as C node
   integer :: max_num_drop = 10000
   logical :: tag_drop_initialized = .false.
   integer :: num_drop
   real(8), dimension(:,:), allocatable :: xdrop,udrop
      ! (:,:) = (dimension,particle_id)
   real(8), dimension(:,:), allocatable :: prop_drop
      ! (:,:) = (property-index,particle_id)
      ! property index: 1, temperature 2, volume
   integer, dimension(:),   allocatable :: droptag
      ! (:) = (particle_id)

   
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
      allocate( tag_flag(imin:imax,jmin:jmax,kmin:kmax) )
      allocate( tag_id  (imin:imax,jmin:jmax,kmin:kmax) )
      allocate( xdrop    (3,max_num_drop) )
      allocate( udrop    (3,max_num_drop) )
      allocate( prop_drop(2,max_num_drop) )
      allocate( droptag  (  max_num_drop) )
   end subroutine initialize_tag_drop

! 
! Purpose: tag isolated droplets and record properties
! 
  subroutine tag_drop()
    integer :: i,j,k, i0,j0,k0
    integer :: isq,jsq,ksq
    integer :: current_id
    integer :: s_queue(200,3),c_queue(200,3) ! max 8 neighbors, record i,j,k
    integer :: ns_queue,is_queue,ic_queue,nc_queue

    if (.not. tag_drop_initialized) call initialize_tag_drop()
    tag_drop_initialized = .true.

    tag_id  (:,:,:) = 0
    tag_flag(:,:,:) = 0
    current_id = 1

    ns_queue = 0
    s_queue(:,:) = 0
    WRITE(*,*) ' start tagging ... '
    do i=imin+1,imax-1; do j=jmin+1,jmax-1; do k=kmin+1,kmax-1
      if ( cvof(i,j,k) > 0.d0 .and. tag_flag(i,j,k) == 0 ) then 
! DEBUG
      WRITE(*,*) current_id, i,j,k
! END DEBUG
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
          end do ! is_queue
          ! mark all C nodes as S nodes
          if ( nc_queue >= 0 ) then 
            s_queue(:,:) = c_queue(:,:)   ! mark all C nodes as S nodes
            ns_queue = nc_queue
          end if ! nc_queue
        end do ! ns_queue>0
        current_id = current_id+1
      end if ! cvof(i,j,k)
    enddo; enddo; enddo

  end subroutine tag_drop

! ==============================================
! output tag of droplets
! ==============================================
   subroutine output_tag()
      implicit none
      integer :: i,j,k

      call tag_drop()

      OPEN(UNIT=89,FILE=TRIM(out_path)//'/tag-'//TRIM(int2text(rank,padding))//'.txt')
      OPEN(UNIT=90,FILE=TRIM(out_path)//'/tag-tecplot'//TRIM(int2text(rank,padding))//'.dat')

      write(*,*) 
      write(*,*) "Print tag_id" 
      write(90,*) 'title= " 3d tag "'
      write(90,*) 'variables = "x", "y", "z", "tag", "c" '
      write(90,*) 'zone i=,',Nx+Ng*2, 'j=',Ny+Ng*2, 'k=',Nz+Ng*2,'f=point'
      do k = kmin,kmax
         do j=jmin,jmax
            do i=imin,imax 
               write(90,'(4(I5,1X),(E15.8))') i,j,k,tag_id(i,j,k),cvof(i,j,k)
            end do ! i
         end do ! j
      end do ! k
      do j=jmax,jmin,-1
         write(89,'(20(I5,1X))') tag_id(imin:imax,j,(Nz+4)/2)
      end do ! j
      CLOSE(89)
      CLOSE(90)
   end subroutine output_tag

! ===============================================
! Testing section
! ===============================================
   subroutine test_Lag_part()
      implicit none
      include 'mpif.h'
      integer :: ierr
                     
      if(test_tag) then
         call output_tag()
      end if
                                                                             
      ! Exit MPI gracefully
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_finalize(ierr)
      stop
   end subroutine test_Lag_part

end module module_Lag_part
