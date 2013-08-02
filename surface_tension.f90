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
  use module_BC
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

  logical :: st_initialized = .false.
contains
!=================================================================================================
  subroutine initialize_surface_tension()
    allocate(  n1(imin:imax,jmin:jmax,kmin:kmax), n2(imin:imax,jmin:jmax,kmin:kmax),  &
               n3(imin:imax,jmin:jmax,kmin:kmax), vof_flag(imin:imax,jmin:jmax,kmin:kmax), &
               height(imin:imax,jmin:jmax,kmin:kmax,6))
      if(nx.ge.500000.or.ny.gt.500000.or.nz.gt.500000) then
         stop 'nx too large'
      endif
! TEMPORARY - Stanley
      height = 2.d6
! END TEMPORARY
   end subroutine initialize_surface_tension

   subroutine get_normals()
      real(8) :: stencil3x3(-1:1,-1:1,-1:1)
      integer :: i,j,k
      integer :: i0,j0,k0
      real(8) :: mxyz(3)
      if(.not.st_initialized) call initialize_surface_tension()
      st_initialized=.true.
      if(ng.lt.2) stop "wrong ng"
      do k=ks-1,ke+1
         do j=js-1,je+1
            do i=is-1,ie+1
               do i0=-1,1; do j0=-1,1; do k0=-1,1
                  stencil3x3(i0,j0,k0) = cvof(i+i0,j+j0,k+k0)
               enddo;enddo;enddo
               call mycs(stencil3x3,mxyz)
               n1(i,j,k) = mxyz(1)
               n2(i,j,k) = mxyz(2)
               n3(i,j,k) = mxyz(3)
            enddo
         enddo
      enddo
   end subroutine get_normals

  subroutine get_flags()
     integer :: i,j,k
     if(.not.st_initialized) call initialize_surface_tension()
     st_initialized=.true.
     if(ng.lt.2) stop "wrong ng"
     do k=kmin,kmax
        do j=jmin,jmax
           do i=imin,imax
              if(cvof(i,j,k).le.0.d0) then
                 vof_flag(i,j,k) = 0
              else if(cvof(i,j,k).ge.1.d0) then
                 vof_flag(i,j,k) = 1
              else
                 vof_flag(i,j,k) = 2
              endif
           enddo
        enddo
     enddo
   end subroutine get_flags

   subroutine get_all_heights
     include 'mpif.h'
     integer :: direction, ierr
     integer :: req(24),sta(MPI_STATUS_SIZE,24)
     if(.not.st_initialized) call initialize_surface_tension()
     st_initialized=.true.

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
     integer, intent(in) :: direction
     integer :: index
     logical :: base_not_found, bottom_n_found, bottom_p_found, top_n_found, top_p_found
     real(8) :: height_p, height_n
     integer :: NDEPTH
     parameter (ndepth=3)
     integer :: si,sj,sk
     integer :: i,j,k
     integer :: i0,j0,k0
     integer :: nd, ierr
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
     else
       stop "bad direction"
     endif
   
     do k=ks,ke
        do j=js,je
           do i=is,ie
!
!  first case: cell is either full or empty
!  odd index: reference fluid (C=1) is below. p-height.
!
              index = 2*(direction-1) + 1
              if(vof_flag(i,j,k).eq.0.and.vof_flag(i-si,j-sj,k-sk).eq.1) then
                 height(i,j,k,index) =  - 0.5d0
              else if(vof_flag(i,j,k).eq.1.and.vof_flag(i+si,j+sj,k+sk).eq.0) then
                 height(i,j,k,index) = 0.5d0
              else
                 height(i,j,k,index) = 2d6
              endif

!  even index: reference fluid (C=1) is above. n-height.

              index = 2*(direction-1) + 2
              if(vof_flag(i,j,k).eq.1.and.vof_flag(i-si,j-sj,k-sk).eq.0) then
                 height(i,j,k,index) =  - 0.5d0
               else if(vof_flag(i,j,k).eq.0.and.vof_flag(i+si,j+sj,k+sk).eq.1) then
                 height(i,j,k,index) = 0.5d0
               else
                 height(i,j,k,index) = 2d6
               endif
               !
               !  end empty/full case
               !
               !  second case: cell is fractional
               !
              if(vof_flag(i,j,k).eq.2) then
                 
                 height_p = cvof(i,j,k) - 0.5d0
                 height_n = - height_p

                 bottom_p_found=.false.
                 bottom_n_found=.false.
                 top_p_found=.false.
                 top_n_found=.false.
                 
                 ! going down
                 ! First pass of "p-heights"

                 i0 = i - si; j0 = j - sj; k0 = k - sk
                 base_not_found = .true.
                 nd = 1
                 do while ( base_not_found ) 
                    height_p = height_p + cvof(i0,j0,k0) - 1.d0
                    i0 = i0 - si; j0 = j0 - sj; k0 = k0 - sk
                    nd = nd + 1
                    if(vof_flag(i0,j0,k0).eq.1) then
                       bottom_p_found = .true.
                       base_not_found = .false.
                    endif
                    if(nd .eq. NDEPTH) base_not_found = .false.
                 end do

                 ! First pass of "n-heights"

                 i0 = i - si; j0 = j - sj; k0 = k - sk
                 base_not_found = .true.
                 nd = 1
                 do while ( base_not_found ) 
                    height_n = height_n - cvof(i0,j0,k0) 
                    i0 = i0 - si; j0 = j0 - sj; k0 = k0 - sk
                    nd = nd + 1
                    if(vof_flag(i0,j0,k0).eq.0) then
                       bottom_n_found = .true.
                       base_not_found = .false.
                    endif
                    if(nd .eq. NDEPTH) base_not_found = .false.
                 end do

                 ! same thing going up
                 ! Second pass of "p-heights"

                 i0 = i + si; j0 = j + sj; k0 = k + sk
                 nd = 1
                 base_not_found = .true.
                 do while ( base_not_found ) 
                    height_p = height_p + cvof(i0,j0,k0) 
                    i0 = i0 + si; j0 = j0 + sj; k0 = k0 + sk
                    nd = nd + 1
                    if(vof_flag(i0,j0,k0).eq.0) then
                       top_p_found = .true.
                       base_not_found = .false.
                    endif
                    if(nd .eq. NDEPTH) base_not_found = .false.
                 end do

                 ! Second pass of "n-heights"

                 i0 = i + si; j0 = j + sj; k0 = k + sk
                 nd = 1
                 base_not_found = .true.
                 do while ( base_not_found ) 
                    height_n = height_n - cvof(i0,j0,k0) + 1.d0
                    i0 = i0 + si; j0 = j0 + sj; k0 = k0 + sk
                    nd = nd + 1
                    if(vof_flag(i0,j0,k0).eq.1) then
                       top_n_found = .true.
                       base_not_found = .false.
                    endif
                    if(nd .eq. NDEPTH) base_not_found = .false.
                 end do

! put everything in the height array.

! "p-heights" have indexes 1, 3, 5
 
                 index = 2*(direction-1) + 1
                 if (bottom_p_found.and.top_p_found) then 
                    height(i,j,k,index) = height_p
                 else
                    height(i,j,k,index) = 2d6
                 endif

! "n-heights" have indexes 2, 4, 6

                 index = 2*(direction-1) + 2
                 if (bottom_n_found.and.top_n_found) then 
                    height(i,j,k,index) = height_n
                 else
                    height(i,j,k,index) = 2d6
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
     integer i,j,k,d,index
     real(8) h, th
     k = nz/2 + 2  ! +2 because of ghost layers
     j = ny/2

     if(k<ks.or.k>ke) return
     if(j<js.or.j>je) return

     OPEN(UNIT=89,FILE=TRIM(out_path)//'/height-'//TRIM(int2text(rank,padding))//'.txt')
     OPEN(UNIT=90,FILE=TRIM(out_path)//'/reference-'//TRIM(int2text(rank,padding))//'.txt')

!     write(89,*) " "
!     write(89,*) " p heights at  z = L/2 as a function of x"
!     write(89,*) " "
 
     index=5
     do i=is,ie
        th = (- z(k) + zlength/2.d0)/dx(i) + A_h*cos(2.*3.14159*x(i)/xlength)
        write(90,100) x(i),th
        if (height(i,j,k,index).lt.1d6) then
           h = height(i,j,k,index)
        else
           ! search for height
           d=0
           do while(d.lt.3)
              d = d + 1
              if (height(i,j,k+d,index).lt.1d6) then
                 h = height(i,j,k+d,index) + d
                 height(i,j,k,index) = h 
              else if (height(i,j,k-d,index).lt.1d6) then
                 h = height(i,j,k-d,index) - d
                 height(i,j,k,index) = h 
              endif
           enddo
        endif
        if(height(i,j,k,index).gt.1d6) then
           write(89,101) x(i),' -'
        else
           write(89,100) x(i), h 
        endif
     enddo

     !     write(89,*) " "
     !     write(89,*) " n heights at  z = L/2 as a function of x, not searched"
     !     write(89,*) " "

     if(1==0) then
        do i=is,ie
           write(89,103) x(i),height(i,j,k,1:6),height(i,j,k,1:6)
        enddo
        do k=ks,ke
           write(*,104) cvof(is:ie,ny/2,k)
        enddo
        print *, " "
        do k=ks,ke
           write(*,105) cvof(is:ie,ny/2,k)
        enddo

        write(89,*) " "
        write(89,*) " n heights in z"
        write(89,*) " "
        do k=ke,ks,-1
           write(89,1041) k, height(is:ie,ny/2,k,6)
        enddo

        write(89,*) " "
        do k=ke,ks,-1
           write(89,1051) k, height(is:ie,ny/2,k,6)
        enddo

        write(89,*) " "
        write(89,*) " p heights in z"
        write(89,*) " "
        do k=ke,ks,-1
           write(89,1041) k, height(is:ie,ny/2,k,5)
        enddo

        write(89,*) " "
        do k=ke,ks,-1
           write(89,1051) k, height(is:ie,ny/2,k,5)
        enddo
     endif  ! 1==0

100  format(2(f24.16))
101  format(f24.16,A2)
103  format(7(f16.8),6(I2))
104  format(16(f5.1,' '))
1041 format(I3,16(f5.1,' '))
105  format(16(I2,' '))
1051 format(I3,16(I2,' '))

     close(89)
     close(90)

   end subroutine output_heights
!=================================================================================================
!   Check if we find heights in the neighboring cells
!=================================================================================================
   subroutine get_local_heights(i0,j0,k0,nfound,indexfound,hloc)
      implicit none
      integer :: NDEPTH
      parameter (NDEPTH=3)
      integer, intent(in) :: i0,j0,k0
      integer, intent(out) :: nfound, indexfound
      real(8), intent(out) :: hloc(3,3)  ! last element is index 
      integer :: d,s
      integer :: i,j,k,m,n
      integer :: i1(3,3,3), j1(3,3,3), k1(3,3,3)
      integer :: index
      logical :: notfound
      integer :: si,sj,sk
!
! mapping
!
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
!
!  Loop over directions
! 
      d=0
      notfound=.true.
      do while (d.lt.3.and.notfound)
         d = d+1
         if(d.eq.1) then
            si=1; sj=0; sk=0;
         else if (d.eq.2) then
            si=0; sj=1; sk=0;
         else if (d.eq.3) then
            si=0; sj=0; sk=1
         else
            stop "bad direction"
         endif
         index = 2*(d-1)
         do while (index.lt.2*(d-1)+2.and.notfound)
            index = index + 1
            hloc = 2d6
            nfound = 0
            do m=-1,1 
               do n=-1,1
                  if(height(i1(m,n,d),j1(m,n,d),k1(m,n,d),index).lt.1d6) then
                     ! one height found
                     hloc(m,n) = height(i1(m,n,d),j1(m,n,d),k1(m,n,d),index)
                     nfound = nfound + 1
                  else
                     do while(s.lt.NDEPTH)
                        s = s + 1
                        if (height(i+si*s,j+sj*s,k+sk*s,index).lt.1d6) then
                           hloc(m,n) = height(i+si*s,j+sj*s,k+sk*s,index) + s
                           nfound = nfound + 1
                           s = NDEPTH  ! to exit loop
                        else if  (height(i-si*s,j-sj*s,k-sk*s,index).lt.1d6) then
                          hloc(m,n) = height(i-si*s,j-sj*s,k-sk*s,index) - s
                          nfound = nfound + 1
                          s = NDEPTH  ! to exit loop
                        endif
                     end do
                  end if ! found at same level
               end do ! n
            end do ! m 
            if(nfound.ge.6) then
               notfound = .false.
               indexfound = index
            endif
         end do ! index
      end do ! d
   end subroutine get_local_heights

   subroutine get_curvature(i0,j0,k0,kappa,indexCurv)
      implicit none
      integer, intent(in) :: i0,j0,k0
      real(8), intent(out) :: kappa  
      integer, intent(out) :: indexCurv

      integer :: nfound,d,indexfound
      real(8) :: h(-1:1,-1:1),hm,hn,hmm,hnn,hmn
      integer :: nCentroids
      integer :: si,sj,sk

      call get_local_heights(i0,j0,k0,nfound,indexfound,h)
      d=(indexfound-1)/2+1

      kappa = 0.d0
      if ( nfound == 9 ) then
         hm  = (h(1,0)-h(-1,0))/2.d0
         hn  = (h(0,1)-h(0,-1))/2.d0
         hmm =  h(1,0)-2.d0*h(0,0)+h(-1,0)
         hnn =  h(0,1)-2.d0*h(0,0)+h(0,-1)
         hmn = (h(1,1)-h(-1,1)-h(1,-1)+h(-1,-1))/4.d0
         kappa = (hmm + hnn + hmm*hn*hn + hnn*hm*hm - 2.d0*hmn*hm*hn) &
                /(1.d0+hm*hm+hn*hn)**(1.5d0)
         indexCurv = indexfound
      else if ( nfound >= 6 ) then 
         ! least square fit
      else if ( nfound <6 ) then 
         nCentroids = 6 - nfound
         ! find centroid
         ! least squre fit
      end if ! nfound

   end subroutine get_curvature

   subroutine output_curvature()
      implicit none
      
      integer :: i,j,k,indexCurv
      real(8) :: kappa, kappa_exact, radius

      do i = is,ie; do j=js,je; do k=ks,ke
         call get_curvature(i,j,k,kappa,indexCurv)
         kappa_exact = 1.0d0/radius
!         kappa_norm2 
      end do; end do; end do
   end subroutine output_curvature
 
end module module_surface_tension

