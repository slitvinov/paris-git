!=================================================================================================
!=================================================================================================
! Paris-0.1
! Extended from Code: FTC3D2011 (Front Tracking Code for 3D simulations)
! and Surfer. 
! 
! Authors: Sadegh Dabiri, Gretar Tryggvason
! author for VOF extenstions Stephane Zaleski (zaleski@dalembert.upmc.fr) 
! Contact: sdabiri@gmail.com
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
! module_surface_tension: Contains definition of variables for surface tension from
!  Volume of Fluid interface tracking.
!-------------------------------------------------------------------------------------------------
module module_surface_tension
  use module_grid
  use module_IO
  use module_tmpvar
  use module_VOF
  implicit none
  real(8), dimension(:,:,:), allocatable :: n1,n2,n3 ! normals
  real(8), dimension(:,:,:,:), allocatable :: height ! normals
   ! 4th index: 1 for positive height in x, 2 for negative height in x
   !            3 for positive height in y, 4 for negative height in y, etc... 
 integer, dimension(:,:,:), allocatable :: vof_flag ! 
  !   0 empty
  !   1 full
  !   2 fractional
  integer, dimension(:,:,:,:), allocatable :: height_flag ! 
  !   0 undecided (not fully tested yet)
  !   1 height found
  !   2 no height found
  !   3 other cases (for instance empty cell)
  ! 4th index: 1 for positive height in x, 2 for negative height in x, etc... 
  !            3 for positive height in y, 4 for negative height in y, etc... 
  logical st_initialized = .false.
contains
!=================================================================================================
  subroutine initialize_surface_tension()
    allocate(  n1(imin:imax,jmin:jmax,kmin:kmax), n2(imin:imax,jmin:jmax,kmin:kmax),  &
               n3(imin:imax,jmin:jmax,kmin:kmax))
  end subroutine initialize_surface_tension

  subroutine get_normals()
    if(.not.st_initialized) call initialize_surface_tension()
    st_initialized=.true.
    if(ng.lt.2) stop "wrong ng"
    do k=ks-1,ke+1
       do j=js-1,je+1
          do i=is-1,ie+1
             do i0=-1,1; do j0=-1,1; do k0=-1,1
                 stencil3x3(i0,j0,k0) = c(i+i0,j+j0,k+k0)
              enddo;enddo;enddo
              call mycs(stencil3x3,mxyz)
              n1(i,j,k) = mxyz(1)
              n2(i,j,k) = mxyz(2)
              n3(i,j,k) = mxyz(3)
           enddo
        enddo
     enddo
   end subroutine get_normals

   subroutine get_all_heights
     integer :: direction
     integer :: req(24),sta(MPI_STATUS_SIZE,24)
     do direction=1,3
        call get_heights(direction)
     enddo
     call ghost_x(height(:,:,:,1),2,req( 1: 4));  
     call ghost_x(height(:,:,:,2),2,req( 5: 8)); 
     call ghost_x(height(:,:,:,3),2,req( 9:12));
     call ghost_x(height(:,:,:,4),2,req(13:16));
     call ghost_x(height(:,:,:,5),2,req(17:20));
     call ghost_x(height(:,:,:,6),2,req(21:24));
     call MPI_WAITALL(24,req(1:24),sta(:,1:24),ierr)

     call ghost_y(height(:,:,:,1),2,req( 1: 4));  
     call ghost_y(height(:,:,:,2),2,req( 5: 8)); 
     call ghost_y(height(:,:,:,3),2,req( 9:12));
     call ghost_y(height(:,:,:,4),2,req(13:16));
     call ghost_y(height(:,:,:,5),2,req(17:20));
     call ghost_y(height(:,:,:,6),2,req(21:24));
     call MPI_WAITALL(24,req(1:24),sta(:,1:24),ierr)

     call ghost_z(height(:,:,:,1),2,req( 1: 4));  
     call ghost_z(height(:,:,:,2),2,req( 5: 8)); 
     call ghost_z(height(:,:,:,3),2,req( 9:12));
     call ghost_z(height(:,:,:,4),2,req(13:16));
     call ghost_z(height(:,:,:,5),2,req(17:20));
     call ghost_z(height(:,:,:,6),2,req(21:24));
     call MPI_WAITALL(24,req(1:24),sta(:,1:24),ierr)

   end subroutine get_all_heights

   subroutine get_heights(direction)
     integer, intent=in :: direction
     integer :: index
     logical base_not_found
     real(8) height_p, height_n
     integer NDEPTH=3
     integer si,sj,sk
     ! NDEPTH is the depth of layers tested above or below the reference cell. 
     ! including the central layer
     ! NDEPTH*2 - 1 = 5 means a 5 x 3 stencil. 
     !
     !  height_p : height for a normal pointing up (reference phase under the other phase)
     !  height_n : height for a normal pointing down (reference phase above the other phase)
    
     if(direction.eq.1) then
        si=1; sj=0; sk=0;
     else if (direction.eq.2) then
        si=0; sj=1; sk=0;
     else if (direction.eq.3) then
        si=0; sj=0; sk=1
     endif
   
     do k=ks,ke
        do j=js,je
           do i=is,ie
!
!  first case: cell is either full or empty
!  odd index: normal pointing up
!
              index = 2*(direction-1) + 1
              if(vof_flag(i,j,k).eq.0.and.vof_flag(i-si,j-sj,k-sk).eq.1) then
                 height(i,j,k,index) =  - 0.5d0
                 height_flag(i,j,k,index) = 1
              else if(vof_flag(i,j,k).eq.1.and.vof_flag(i+si,j+sj,k+sk).eq.0) then
                 height(i,j,k,index) = 0.5d0
                 height_flag(i,j,k,index) = 1
              else
                 height_flag(i,j,k,index) = 2
              endif

!  even index: normal pointing down

              index = 2*(direction-1) + 2
              if(vof_flag(i,j,k).eq.1.and.vof_flag(i-si,j-sj,k-sk).eq.0) then
                 height(i,j,k,index) =  - 0.5d0
                 height_flag(i,j,k,index) = 1
               else if(vof_flag(i,j,k).eq.0.and.vof_flag(i+si,j+sj,k+sk).eq.1) then
                 height(i,j,k,index) = 0.5d0
                 height_flag(i,j,k,index) = 1
               endif
!
!  end empty/full case
!
!  second case: cell is fractional
!
              if(vof_flag(i,j,k).eq.2) then
                 i0 = i; j0 = j; k0 = k
                 height_p = c(i0,j0,k0) - 0.5d0
                 height_n = c(i0,j0,k0) - 0.5d0
                 i0 = i0 - si; j0 = j0 -sj; k0 = k0 - sk
                 base_not_found = .true.
                 nd = 1
                 do while ( base_not_found ) 
                    height_p = height_p + c(i0,j0,k0) - 1.d0
                    i0 = i0 - si; j0 = j0 -sj; k0 = k0 - sk
                    nd = nd + 1
                    if(vof_flag(i0,j0,k0).eq.1) then
                       bottom_p_found = .true.
                       base_not_found = .false.
                    endif
                    if(nd .eq. NDEPTH) base_not_found = .false.
                 end do

                 i0 = i - si; j0 = j -sj; k0 = k - sk
                 nd = 1
                 base_not_found = .true.
                 do while ( base_not_found ) 
                    height_n = height_n + c(i0,j0,k0) - 1.d0
                    i0 = i0 - si; j0 = j0 -sj; k0 = k0 - sk
                    nd = nd + 1
                    if(vof_flag(i0,j0,k0).eq.0) then
                       bottom_n_found = .true.
                       base_not_found = .false.
                    endif
                    if(nd .eq. NDEPTH) base_not_found = .false.
                 end do

                 ! same thing going up

                 i0 = i + si; j0 = j + sj; k0 = k + sk
                 nd = 1
                 base_not_found = .true.
                 do while ( base_not_found ) 
                    height_p = height_p + c(i0,j0,k0) - 1.d0
                    i0 = i0 + si; j0 = j0 + sj; k0 = k0 + sk
                    nd = nd + 1
                    if(vof_flag(i0,j0,k0).eq.0) then
                       top_p_found = .true.
                       base_not_found = .false.
                    endif
                    if(nd .eq. NDEPTH) base_not_found = .false.
                 end do

                 i0 = i + si; j0 = j + sj; k0 = k + sk
                 nd = 1
                 base_not_found = .true.
                 do while ( base_not_found ) 
                    height_n = height_n + c(i0,j0,k0) - 1.d0
                    i0 = i0 + si; j0 = j0 + sj; k0 = k0 + sk
                    nd = nd + 1
                    if(vof_flag(i0,j0,k0).eq.0) then
                       top_n_found = .true.
                       base_not_found = .false.
                    endif
                    if(nd .eq. NDEPTH) base_not_found = .false.
                 end do

! put everything in the height array. 
                 index = 2*(direction-1) + 1
                 if (bottom_p_found.and.top_p_found) then 
                    height_flag(i,j,k,index) = 1
                    height(i,j,k,index) = height_p
                 else
                    height_flag(i,j,k,index) = 2
                 endif

                 index = 2*(direction-1) + 2
                 if (bottom_n_found.and.top_n_found) then 
                    height_flag(i,j,k,index) = 1
                    height(i,j,k,index) = height_n
                 else
                    height_flag(i,j,k,index) = 2
                 endif
              endif
! end of fractional cell case
! note: we could check that height orientation and normal orientation agree. 
           enddo
        enddo
     enddo
   end subroutine get_heights    

   subroutine output_heights()
     implicit none
     integer i,j,k
     OPEN(UNIT=89,FILE=TRIM(out_path)//'/height-'//TRIM(int2text(rank,padding))//'.txt')
     k = nz/2
     j = ny/2
     do i=is,ie
        if (height_flag(i,j,k,1).eq.1) then
           write(89,100) x(i),height(i,j,k,5)
        else
           write(89,100) x(i),'-'
        endif
     enddo
   end subroutine output_heights
        
