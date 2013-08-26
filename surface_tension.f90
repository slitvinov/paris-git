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
! module_surface_tension: Contains definition of variables for surface tension from
!  Volume of Fluid interface tracking.
!-------------------------------------------------------------------------------------------------
module module_surface_tension
  use module_grid
  use module_BC
  use module_IO
  use module_tmpvar
  use module_2phase
  use module_VOF
  use module_flow ! for curvature test only
  implicit none
  real(8), dimension(:,:,:), allocatable :: n1,n2,n3 ! normals
  real(8), dimension(:,:,:,:), allocatable :: height ! normals

   ! 4th index: 1 for positive height in x, 2 for negative height in x
   !            3 for positive height in y, 4 for negative height in y, etc... 
  integer, dimension(:,:,:,:), allocatable :: ixheight ! HF flags


!  integer, dimension(:,:,:,:), allocatable :: height_flag ! 
  !   0 undecided (not fully tested yet)
  !   1 height found
  !   2 no height found
  !   3 other cases (for instance empty cell)
  ! 4th index: 1 for positive height in x, 2 for negative height in x, etc... 
  !            3 for positive height in y, 4 for negative height in y, etc... 
  logical :: st_initialized = .false.
contains
!=================================================================================================
  subroutine initialize_surface_tension()
    allocate(  n1(imin:imax,jmin:jmax,kmin:kmax), n2(imin:imax,jmin:jmax,kmin:kmax),  &
               n3(imin:imax,jmin:jmax,kmin:kmax), &
               height(imin:imax,jmin:jmax,kmin:kmax,6))
      if(nx.ge.500000.or.ny.gt.500000.or.nz.gt.500000) then
         stop 'nx too large'
      endif
      height = 2.d6
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
           enddo ! k
        enddo ! j
     enddo ! i
   end subroutine get_heights    

   subroutine output_heights()
     implicit none
     integer i,j,k,d,index
     real(8) h, th
     k = nz/2 + 2  ! +2 because of ghost layers
     j = ny/2 + 2  ! +2 because of ghost layers

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
   subroutine get_local_heights(i0,j0,k0,nfoundmax,indexfound,hlocmax)
      implicit none
      integer :: NDEPTH
      parameter (NDEPTH=3)
      integer, intent(in) :: i0,j0,k0
      integer, intent(out) :: nfoundmax, indexfound
      real(8), intent(out) :: hlocmax(-1:1,-1:1)   
      real(8) :: hloc(-1:1,-1:1)   
      integer :: d,s
      integer :: i,j,k,m,n
      integer :: i1(-1:1,-1:1,3), j1(-1:1,-1:1,3), k1(-1:1,-1:1,3)
      integer :: index,nfound
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
      hloc = 2d6
      nfoundmax = 0 
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
                     s = 0 
                     do while(s.lt.NDEPTH)
                        s = s + 1
                        if (height(i1(m,n,d)+si*s,j1(m,n,d)+sj*s,k1(m,n,d)+sk*s,index).lt.1d6) then
                           hloc(m,n) = height(i1(m,n,d)+si*s,j1(m,n,d)+sj*s,k1(m,n,d)+sk*s,index) + s
                           nfound = nfound + 1
                           s = NDEPTH  ! to exit loop
                        else if  (height(i1(m,n,d)-si*s,j1(m,n,d)-sj*s,k1(m,n,d)-sk*s,index).lt.1d6) then
                           hloc(m,n) = height(i1(m,n,d)-si*s,j1(m,n,d)-sj*s,k1(m,n,d)-sk*s,index) - s
                           nfound = nfound + 1
                           s = NDEPTH  ! to exit loop
                        endif
                     end do
                  end if ! found at same level
               end do ! n
            end do ! m 
            if(nfound.ge.9) then
               notfound = .false.
               hlocmax = hloc
               nfoundmax = nfound
               indexfound = index
            else if (nfound > nfoundmax) then ! return maximum number of heights 
               hlocmax = hloc
               nfoundmax = nfound
               indexfound = index
            end if ! nfound
         end do ! index
      end do ! d
   end subroutine get_local_heights

   subroutine get_curvature(i0,j0,k0,kappa,indexCurv,nfound,nindepend)
      implicit none
      integer, intent(in) :: i0,j0,k0
      real(8), intent(out) :: kappa  
      integer, intent(out) :: indexCurv, nfound
      integer, intent(out) :: nindepend

      integer :: d,indexfound
      real(8) :: h(-1:1,-1:1),a(6)
      integer :: nCentroids

      integer :: m,n,ifit 
      real(8) :: xfit(9),yfit(9),hfit(9)
      logical :: fit_success = .false.

      logical :: independ_flag(9)
      integer :: i1(-1:1,-1:1,3), j1(-1:1,-1:1,3), k1(-1:1,-1:1,3)
      integer :: si,sj,sk,hsign,s
      
      integer :: NDEPTH
      parameter (ndepth=3)
      real(8) :: centroid(3)

      call get_local_heights(i0,j0,k0,nfound,indexfound,h)

      d=(indexfound-1)/2+1

      kappa = 0.d0
      nindepend = 0 
      ! if all nine heights found 
      if ( nfound == 9 ) then
         a(1) = (h(1,0)-2.d0*h(0,0)+h(-1,0))/2.d0
         a(2) = (h(0,1)-2.d0*h(0,0)+h(0,-1))/2.d0
         a(3) = (h(1,1)-h(-1,1)-h(1,-1)+h(-1,-1))/4.d0
         a(4) = (h(1,0)-h(-1,0))/2.d0
         a(5) = (h(0,1)-h(0,-1))/2.d0
         kappa = 2.d0*(a(1)*(1.d0+a(5)*a(5)) + a(2)*(1.d0+a(4)*a(4)) - a(3)*a(4)*a(5)) &
               /(1.d0+a(4)*a(4)+a(5)*a(5))**(1.5d0)
         indexCurv = indexfound
      ! if more than six heights found 
      else if ( nfound > 6 ) then 
         do ifit = 1,9
            m = (ifit-1)/3-1
            n = mod((ifit-1),3)-1
            xfit(ifit) = real(m)
            yfit(ifit) = real(n)
            hfit(ifit) = h(m,n)
            if( hfit(ifit) < 2.d6 ) independ_flag(ifit) = .true.
         end do ! ifit
         call parabola_fit(xfit,yfit,hfit,independ_flag,a,fit_success)
         kappa = 2.d0*(a(1)*(1.d0+a(5)*a(5)) + a(2)*(1.d0+a(4)*a(4)) - a(3)*a(4)*a(5)) &
               /(1.d0+a(4)*a(4)+a(5)*a(5))**(1.5d0)
         indexCurv = indexfound
      ! if less than six heights found
      else
         ! check independent interface points
         do ifit = 1,9
            m = (ifit-1)/3-1
            n = mod((ifit-1),3)-1
            xfit(ifit) = real(m)
            yfit(ifit) = real(n)
            hfit(ifit) = h(m,n)
            ! use height if exist
            if ( hfit(ifit) < 2.d6 ) then
               independ_flag(ifit) = .true.
               nindepend = nindepend + 1
            ! find centroid if no height
            else 
               d = (indexfound-1)/2 + 1
               hsign = mod((indexfound-1),2)*2-1
 
               if(d.eq.1) then
                  si=1; sj=0; sk=0;
               else if (d.eq.2) then
                  si=0; sj=1; sk=0;
               else if (d.eq.3) then
                  si=0; sj=0; sk=1
               else
                  stop "bad direction"
               end if

               ! create mapping for i1,j1,k1
               ! d=1
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

               ! search cut cell in stencil
               !  first search at the same layer
               if(vof_flag(i1(m,n,d),j1(m,n,d),k1(m,n,d)) == 2 ) then
                  call FindCutAreaCentroid(i1(m,n,d),j1(m,n,d),k1(m,n,d),centroid)
                  xfit(ifit) = centroid(1) + m
                  yfit(ifit) = centroid(2) + n
                  hfit(ifit) = centroid(3)
                  independ_flag(ifit) = .true.
                  nindepend = nindepend + 1
               !  search the next two layers along indexfound direction
               else 
                  s = 0 
                  do while(s.lt.NDEPTH)
                     s = s + 1
                     if ( vof_flag(i1(m,n,d)+si*hsign*s,j1(m,n,d)+sj*hsign*s,&
                                   k1(m,n,d)+sk*hsign*s) == 2 ) then
                        call FindCutAreaCentroid(i1(m,n,d)+si*hsign*s,  &
                                                 j1(m,n,d)+sj*hsign*s,  &
                                                 k1(m,n,d)+sk*hsign*s,centroid)
                        xfit(ifit) = centroid(1) + m
                        yfit(ifit) = centroid(2) + n
                        hfit(ifit) = centroid(3) + s*hsign
                        independ_flag(ifit) = .true.
                        nindepend = nindepend + 1
                        s = NDEPTH  ! to exit loop
                     endif
                  end do
               end if ! found at same level
            end if ! hfit(ifit)
         end do !ifit

         ! check if there are six independent positions
         if ( nindepend < 6 ) then 
            kappa = 0.d0
            indexCurv = 0
            return
         else 
            call parabola_fit(xfit,yfit,hfit,independ_flag,a,fit_success)
            kappa = 2.d0*(a(1)*(1.d0+a(5)*a(5)) + a(2)*(1.d0+a(4)*a(4)) - a(3)*a(4)*a(5)) &
                  /(1.d0+a(4)*a(4)+a(5)*a(5))**(1.5d0)
            indexCurv = indexfound
         end if ! nindepend
      end if ! nfound 
! TEMPORARY 
      if ( nfound < 9 ) then
            WRITE(*,*) ' WARNING: nfound < 9:',i0,j0,k0,nfound,nindepend
      end if ! nfound
! END TEMPORARY

   end subroutine get_curvature

   subroutine output_curvature()
      implicit none
      
      integer :: i,j,k,indexCurv
      integer :: ib
      real(8) :: kappa
      real(8) :: angle 
      real(8) :: kappamin
      real(8) :: kappamax
      real(8) :: kappa_exact
      real(8) :: hex,hpex,hmex,dhex,d2hex,kappa_hex
      real(8) :: hnum,hpnum,hmnum,dhnum,d2hnum,hloc(-1:1,-1:1)
      real(8) :: err_h,err_hp,err_hm,err_dh,err_d2h,err_K
      real(8) :: L2_err_h,L2_err_hp,L2_err_hm,L2_err_dh,L2_err_d2h,L2_err_K
      real(8) :: S2_err_h,S2_err_hp,S2_err_hm,S2_err_dh,S2_err_d2h,S2_err_K
      real(8) :: Lm_err_h,Lm_err_hp,Lm_err_hm,Lm_err_dh,Lm_err_d2h,Lm_err_K
      integer :: sumCount,nfound,indexfound,nindepend
      real(8), parameter :: PI= 3.14159265359d0

      OPEN(UNIT=89,FILE=TRIM(out_path)//'/curvature-'//TRIM(int2text(rank,padding))//'.txt')
      OPEN(UNIT=90,FILE=TRIM(out_path)//'/reference-'//TRIM(int2text(rank,padding))//'.txt')
      OPEN(UNIT=91,FILE=TRIM(out_path)//'/bigerror-'//TRIM(int2text(rank,padding))//'.txt')
      ib = 1
      kappamin = 1d20
      kappamax = -1d20
      if ( test_curvature ) then 
         kappa_exact = 2.d0/rad(ib)
         do i=is,ie; do j=js,je; do k=ks,ke
            ! find curvature only for cut cells
            if (vof_flag(i,j,k) == 2 ) then 
               call get_curvature(i,j,k,kappa,indexCurv,nfound,nindepend)
               kappa = kappa*dble(Nx)
               kappamax = max(ABS(kappa),kappamax)
               kappamin = min(ABS(kappa),kappamin)
               angle = atan2(y(j)-yc(ib),x(i)-xc(ib))/PI*180.d0
               write(89,*) angle,ABS(kappa)
               write(90,*) angle,kappa_exact
               err_K = ABS(ABS(kappa)-kappa_exact)/kappa_exact
               if ( err_K > 0.1d0 ) &
                  write(91,'(4(I3,1X),2(E15.8,1X),I2)') i,j,k,indexCurv,kappa,kappa_exact,nfound
            end if ! cvof(i,j,k)
         end do; end do; end do
      else if ( test_curvature_2D) then 
         kappa_exact = 1.d0/rad(ib)
         sumCount = 0
         S2_err_h=0.d0;S2_err_hp=0.d0;S2_err_hm=0.d0;S2_err_dh=0.d0;S2_err_d2h=0.d0;S2_err_K=0.d0
         Lm_err_h=0.d0;Lm_err_hp=0.d0;Lm_err_hm=0.d0;Lm_err_dh=0.d0;Lm_err_d2h=0.d0;Lm_err_K=0.d0
         do i=is,ie; do j=js,je; do k=ks,ke
            if (vof_flag(i,j,k) == 2 .and. k==(Nz+4)/2) then 
               call get_curvature(i,j,k,kappa,indexCurv,nfound,nindepend)
               kappa = kappa*dble(Nx)
               kappamax = max(ABS(kappa),kappamax)
               kappamin = min(ABS(kappa),kappamin)
               angle = atan2(y(j)-yc(ib),x(i)-xc(ib))/PI*180.d0
               write(89,*) angle,ABS(kappa)
               write(90,*) angle,kappa_exact
               call CalExactHeight_Circle(x(i),y(j),1.d0/dble(Nx),1.d0/dble(Ny),&
                  xc(ib),yc(ib),rad(ib),indexCurv,hex,hpex,hmex,dhex,d2hex)
               call get_local_heights(i,j,k,nfound,indexfound,hloc)
               hnum  = hloc( 0,0)/dble(Nx)
               hpnum = hloc( 1,0)/dble(Nx)
               hmnum = hloc(-1,0)/dble(Nx)
               dhnum = (hpnum-hmnum)/(2.d0/dble(Nx))
               d2hnum = (hpnum-2.d0*hnum+hmnum)/(1.d0/dble(Nx))**2.d0
               err_h    = ABS(hnum   - hex )
               err_hp   = ABS(hpnum  - hpex)
               err_hm   = ABS(hmnum  - hmex)
               err_dh   = ABS(dhnum  - dhex)
               err_d2h  = ABS(d2hnum - d2hex)/ABS(d2hex)
               err_K    = ABS(ABS(kappa)-kappa_exact)/kappa_exact
               write(91,'(7(E15.8,1X),6(I5,1X))') angle,err_h,err_hp,err_hm,& 
                                         err_dh,err_d2h,err_K,indexCurv,i,j,nfound,nindepend
               sumCount = sumCount + 1
               S2_err_h    = S2_err_h  + err_h  **2.d0
               S2_err_hp   = S2_err_hp + err_hp **2.d0
               S2_err_hm   = S2_err_hm + err_hm **2.d0
               S2_err_dh   = S2_err_dh + err_dh **2.d0
               S2_err_d2h  = S2_err_d2h+ err_d2h**2.d0
               S2_err_K    = S2_err_K  + err_K  **2.d0

               Lm_err_h    = MAX(Lm_err_h,   err_h ) 
               Lm_err_hp   = MAX(Lm_err_hp,  err_hp) 
               Lm_err_hm   = MAX(Lm_err_hm,  err_hm) 
               Lm_err_dh   = MAX(Lm_err_dh,  err_dh) 
               Lm_err_d2h  = MAX(Lm_err_d2h, err_d2h) 
               Lm_err_K    = MAX(Lm_err_K,   err_K) 
            end if ! cvof(i,j,k)
         end do; end do; end do
         L2_err_h    = sqrt(S2_err_h)  /dble(sumCount)
         L2_err_hp   = sqrt(S2_err_hp) /dble(sumCount)
         L2_err_hm   = sqrt(S2_err_hm) /dble(sumCount)
         L2_err_dh   = sqrt(S2_err_dh) /dble(sumCount)
         L2_err_d2h  = sqrt(S2_err_d2h)/dble(sumCount)
         L2_err_K    = sqrt(S2_err_K)  /dble(sumCount)
         write(*,*) 'L2 Norm:'
         write(*,'(I5,1X,6(E15.8,1X))') Nx,L2_err_h,L2_err_hp,L2_err_hm,L2_err_dh,L2_err_d2h,L2_err_K
         write(*,*) 'Linfty Norm:'
         write(*,'(I5,1X,6(E15.8,1X))') Nx,Lm_err_h,Lm_err_hp,Lm_err_hm,Lm_err_dh,Lm_err_d2h,Lm_err_K
      end if ! test_curvature
      write(*,*) 'max, min, and exact ABS(kappa)', kappamax, kappamin,kappa_exact
      write(*,*) 'max error', MAX(ABS(kappamax-kappa_exact), ABS(kappamin-kappa_exact))/kappa_exact
      CLOSE(89)
      CLOSE(90)
      CLOSE(91)

      contains
! TEMPORARY

   subroutine CalExactHeight_Circle(x1,y1,dx,dy,xc,yc,R,ih,h,hp,hm,dh,d2h)

      implicit none

      real(8), intent(in)  :: x1,y1,dx,dy,xc,yc,R
      integer, intent(in)  :: ih
      real(8), intent(out) :: h,hp,hm,dh,d2h

      integer :: s,d
      real(8) :: xm,ym,xp,yp

      s = -2*mod(ih,2)+1
      d = (ih+1)/2

      xm = x1-dx
      xp = x1+dx
      ym = y1-dy
      yp = y1+dy

      if ( d == 1 ) then 
         h  = s*(ABS(x1-xc)-sqrt(R**2.d0-(y1-yc)**2.d0))
         hp = s*(ABS(x1-xc)-sqrt(R**2.d0-(yp-yc)**2.d0))
         hm = s*(ABS(x1-xc)-sqrt(R**2.d0-(ym-yc)**2.d0))
         dh = s*(y1-yc)/sqrt(R**2.d0-(y1-yc)**2.d0)
         d2h= s*(  1.d0/sqrt(R**2.d0-(y1-yc)**2.d0) &
            + (y1-yc)**2.d0/(R**2.d0-(y1-yc)**2.d0)**1.5d0)
      else if ( d == 2 ) then 
         h  = s*(ABS(y1-yc)-sqrt(R**2.d0-(x1-xc)**2.d0))
         hp = s*(ABS(y1-yc)-sqrt(R**2.d0-(xp-xc)**2.d0))
         hm = s*(ABS(y1-yc)-sqrt(R**2.d0-(xm-xc)**2.d0))
         dh = s*(x1-xc)/sqrt(R**2.d0-(x1-xc)**2.d0)
         d2h= s*(  1.d0/sqrt(R**2.d0-(x1-xc)**2.d0) &
            + (x1-xc)**2.d0/(R**2.d0-(x1-xc)**2.d0)**1.5d0)
      end if ! d

   end subroutine CalExactHeight_Circle
! END TEMPORARY
   end subroutine output_curvature

   subroutine parabola_fit(xfit,yfit,hfit,independ_flag,a,fit_success)
      implicit none

      real(8), intent(in)  :: xfit(9),yfit(9),hfit(9)
      logical, intent(in)  :: independ_flag(9)  
      real(8), intent(out) :: a(6)
      logical, intent(out) :: fit_success

      real(8) :: m(6,6), invm(6,6)
      real(8) :: rhs(6)
      integer :: ifit, im,jm
      logical :: inv_success
      real(8) :: x1,x2,x3,x4,y1,y2,y3,y4

      a = 0.d0

      ! evaluate the linear system for least-square fit
      m   = 0.d0
      rhs = 0.d0
      do ifit = 1, 9
         ! sum over independent points with height
         if ( independ_flag(ifit) ) then
            x1 =    xfit(ifit)
            x2 = x1*xfit(ifit)
            x3 = x2*xfit(ifit)
            x4 = x3*xfit(ifit)
            y1 =    yfit(ifit)
            y2 = y1*yfit(ifit)
            y3 = y2*yfit(ifit)
            y4 = y3*yfit(ifit)
            
            m(1,1) = m(1,1) + x4
            m(2,2) = m(2,2) + y4
            m(3,3) = m(3,3) + x2*y2
            m(4,4) = m(4,4) + x2
            m(5,5) = m(5,5) + y2
            m(6,6) = m(6,6) + 1.d0

            m(1,3) = m(1,3) + x3*y1
            m(1,4) = m(1,4) + x3
            m(1,5) = m(1,5) + x2*y1
            m(2,3) = m(2,3) + x1*y3
            m(2,4) = m(2,4) + x1*y2
            m(2,5) = m(2,5) + y3
            m(3,6) = m(3,6) + x1*y1
            m(4,6) = m(4,6) + x1
            m(5,6) = m(5,6) + y1

            rhs(1) = rhs(1) + x2   *hfit(ifit)
            rhs(2) = rhs(2) + y2   *hfit(ifit)
            rhs(3) = rhs(3) + x1*y1*hfit(ifit)
            rhs(4) = rhs(4) + x1   *hfit(ifit)
            rhs(5) = rhs(5) + y1   *hfit(ifit)
            rhs(6) = rhs(6) +       hfit(ifit)
         end if ! hfit
      end do ! ifit
      m(1,2) = m(3,3)
      m(1,6) = m(4,4)
      m(2,6) = m(5,5)
      m(3,4) = m(1,5)
      m(3,5) = m(2,4)
      m(4,5) = m(3,6)

      do im = 1,6; do jm = 1,6
         if ( im > jm ) m(im,jm) = m(jm,im)
      end do; end do 

      ! Solve linear system
      call FindInverseMatrix(m,invm,6,inv_success)
      if ( inv_success ) then 
         do im=1,6
            do jm=1,6
               a(im) = a(im) + invm(im,jm)*rhs(jm)
            end do
         end do 
         fit_success = .true.
      else
         fit_success = .false.
      end if ! inv_success
   end subroutine parabola_fit
!
!  Subroutine to find the inverse of a square matrix
!  From Stanley's previous code
!
   SUBROUTINE FindInverseMatrix(matrix,inverse,n,inverse_success)
      implicit none

         !---Declarations
        INTEGER, INTENT(IN ) :: n
        real(8), INTENT(IN ), DIMENSION(n,n) :: matrix  !Input A matrix
        real(8), INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix
        logical, INTENT(OUT) :: inverse_success 

        integer :: i, j, k, l
        real(8) :: m
        real(8), DIMENSION(n,2*n) :: augmatrix !augmented matrix


        !Augment input matrix with an identity matrix
        DO i = 1,n
          DO j = 1,2*n
            IF (j <= n ) THEN
              augmatrix(i,j) = matrix(i,j)
            ELSE IF ((i+n) == j) THEN
              augmatrix(i,j) = 1.0d0
            Else
              augmatrix(i,j) = 0.0d0
            ENDIF
          END DO
        END DO
                
        !Ensure diagonal elements are non-zero
        DO k = 1,n-1
          DO j = k+1,n
            IF (augmatrix(k,k) == 0) THEN
               DO i = k+1, n
                 IF (augmatrix(i,k) /= 0) THEN
                   DO  l = 1, 2*n
                     augmatrix(k,l) = augmatrix(k,l)+augmatrix(i,l)
                   END DO
                 ENDIF
               END DO
            ENDIF
          END DO
        END DO
                
        !Reduce augmented matrix to upper traingular form
        DO k =1, n-1
          DO j = k+1, n
            m = augmatrix(j,k)/augmatrix(k,k)
            DO i = k, 2*n
              augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
            END DO
          END DO
        END DO

        !Test for invertibility
        DO i = 1, n
          IF (augmatrix(i,i) == 0) THEN
            write(*,*) "ERROR-Matrix is non-invertible"
            inverse = 0.d0
            inverse_success = .false.
            return
          ENDIF
        END DO
                
        !Make diagonal elements as 1
        DO i = 1 , n
          m = augmatrix(i,i)
          DO j = i , (2 * n)
            augmatrix(i,j) = (augmatrix(i,j) / m)
          END DO
        END DO
                
        !Reduced right side half of augmented matrix to identity matrix
        DO k = n-1, 1, -1
          DO i =1, k
            m = augmatrix(i,k+1)
            DO j = k, (2*n)
              augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
            END DO
          END DO
        END DO
                
        !store answer
        DO i =1, n
          DO j = 1, n
            inverse(i,j) = augmatrix(i,j+n)
          END DO
        END DO
        inverse_success = .true.
                
   END SUBROUTINE FindInverseMatrix

   subroutine FindCutAreaCentroid(i,j,k,centroid)
      implicit none

      integer, intent(in)  :: i,j,k
      real(8), intent(out) :: centroid(3)

      integer :: i0,j0,k0
      real(8) :: dmx,dmy,dmz, mxyz(3)
      real(8) :: invx,invy,invz
      real(8) :: alpha, alphamax,al3d
      real(8) :: stencil3x3(-1:1,-1:1,-1:1)

      ! find cut area centroid 
      !***
      !     (1) normal vector: dmx,dmy,dmz, and |dmx|+|dmy|+|dmz| = 1.
      !     (2) dmx,dmy,dmz>0 and record sign
      !     (3) get alpha;               
      !     (4) compute centroid with dmx,dmy,dmz and alpha;
      !     (5) transfer to local coordinate;
      !*(1)*
      do i0=-1,1; do j0=-1,1; do k0=-1,1
         stencil3x3(i0,j0,k0) = cvof(i+i0,j+j0,k+k0)
      enddo;enddo;enddo
      call mycs(stencil3x3,mxyz)
      dmx = mxyz(1)
      dmy = mxyz(2)
      dmz = mxyz(3)
      !*(2)*  
      invx = 1.d0
      invy = 1.d0
      invz = 1.d0
      if (dmx .lt. 0.0d0) then
         dmx = -dmx
         invx = -1.d0
      endif
      if (dmy .lt. 0.0d0) then
         dmy = -dmy
         invy = -1.d0
      endif
      if (dmz .lt. 0.0d0) then
         dmz = -dmz
         invz = -1.d0
      endif
      !*(3)*  
      alpha = al3d(dmx,dmy,dmz,cvof(i,j,k))
      !*(4)*  
      call ComputeCentroid(dmx,dmy,dmz,alpha,centroid)
      !*(5)*
      ! rotate
      centroid(1) = centroid(1)*invx
      centroid(2) = centroid(2)*invy
      centroid(3) = centroid(3)*invz
      ! shift
      centroid(1) = centroid(1) - invx*0.5d0
      centroid(2) = centroid(2) - invy*0.5d0
      centroid(3) = centroid(3) - invz*0.5d0

!      dmx = invx*dmx
!      dmy = invy*dmy
!      dmz = invz*dmz

   end subroutine FindCutAreaCentroid

   subroutine ComputeCentroid(m1,m2,m3,alpha,centroid)
      implicit none

      real(8), intent(in)  :: m1,m2,m3,alpha
      real(8), intent(out) :: centroid(3)
      
      real(8) :: m123,alphamax,Area
      real(8) :: x0(3),x1(3),x2(3),x3(3),x4(3),x5(3),x6(3)

      alphamax = m1+m2+m3

      m123 = sqrt(m1*m1+m2*m2+m3*m3)/2.d0/m1/m2/m3
      Area   = m123*( alpha*alpha &
                  - func(2.d0,(alpha-m1)) &
                  - func(2.d0,(alpha-m2)) &
                  - func(2.d0,(alpha-m3)) &
                  + func(2.d0,(alpha-alphamax+m1)) &
                  + func(2.d0,(alpha-alphamax+m2)) &
                  + func(2.d0,(alpha-alphamax+m3)) )

      x0(:) = [alpha/m1/3.d0,alpha/m2/3.d0,alpha/m3/3.d0]
      x1(:) = [(alpha/m1+2.d0)/3.d0,(alpha-m1)/m2/3.d0,(alpha-m1)/m3/3.d0]
      x2(:) = [(alpha-m2)/m1/3.d0,(alpha/m2+2.d0)/3.d0,(alpha-m2)/m3/3.d0]
      x3(:) = [(alpha-m3)/m1/3.d0,(alpha-m3)/m2/3.d0,(alpha/m3+2.d0)/3.d0]
      x4(:) = [(alpha-alphamax+m1)/m1/3.d0,((alpha-m3)/m2+2.d0)/3.d0,((alpha-m2)/m3+2.d0)/3.d0]
      x5(:) = [((alpha-m3)/m1+2.d0)/3.d0,(alpha-alphamax+m2)/m2/3.d0,((alpha-m1)/m3+2.d0)/3.d0]
      x6(:) = [((alpha-m2)/m1+2.d0)/3.d0,((alpha-m1)/m2+2.d0)/3.d0,(alpha-alphamax+m3)/m3/3.d0]

      centroid(:) = m123/Area*( alpha*alpha*x0(:) &
                              - func(2.d0,(alpha-m1))*x1(:) &
                              - func(2.d0,(alpha-m2))*x2(:) &
                              - func(2.d0,(alpha-m3))*x3(:) &
                              + func(2.d0,(alpha-alphamax+m1))*x4(:) &
                              + func(2.d0,(alpha-alphamax+m2))*x5(:) &
                              + func(2.d0,(alpha-alphamax+m3))*x6(:) )

! -------------------------------------------------------------------------------------------------------
      contains
         function func(n,z)
            implicit none
   
            real(8) :: n,z
            real(8) :: func
   
            if ( z > 0.d0 ) then
               func = z**n
            else
               func = 0.d0
            end if !z
         end function func

   end subroutine ComputeCentroid
!=========================================================================================================
!
!  Testing section
! 
!=========================================================================================================

   subroutine test_VOF_HF()
     implicit none
     integer :: calc_imax

     if(calc_imax(vof_flag)/=2) then
        write(*,*) calc_imax(vof_flag), "expecting maximum flag = 2"
        call pariserror("bad flags")
     endif
     call get_all_heights()
     if(test_heights) then
        call output_heights()
     else if(test_curvature .or. test_curvature_2D) then
        call output_curvature()
     end if
  end subroutine test_VOF_HF

!=========================================================================================================
!
!  Begin Ruben + Phil routines
! 
!=========================================================================================================

  subroutine get_flags_rph()
     integer :: i,j,k,q
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

              do q=1,3
                 ixheight(i,j,k,q)=0
                 height(i,j,k,q)=0.d0
              enddo

           enddo
        enddo
     enddo

!!$     call ghost_x(height(:,:,:,1),2,req( 1: 4));  
!!$     call ghost_x(height(:,:,:,2),2,req( 5: 8)); 
!!$     call ghost_x(height(:,:,:,3),2,req( 9:12));
!!$     call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
!!$
!!$     call ghost_y(height(:,:,:,1),2,req( 1: 4));  
!!$     call ghost_y(height(:,:,:,2),2,req( 5: 8)); 
!!$     call ghost_y(height(:,:,:,3),2,req( 9:12));
!!$     call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
!!$
!!$     call ghost_z(height(:,:,:,1),2,req( 1: 4));  
!!$     call ghost_z(height(:,:,:,2),2,req( 5: 8)); 
!!$     call ghost_z(height(:,:,:,3),2,req( 9:12));
!!$     call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
!!$
!!$     call ighost_x(ixheight(:,:,:,1),2,req( 1: 4));  
!!$     call ighost_x(ixheight(:,:,:,2),2,req( 5: 8)); 
!!$     call ighost_x(ixheight(:,:,:,3),2,req( 9:12));
!!$     call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
!!$
!!$     call ighost_y(ixheight(:,:,:,1),2,req( 1: 4));  
!!$     call ighost_y(ixheight(:,:,:,2),2,req( 5: 8)); 
!!$     call ighost_y(ixheight(:,:,:,3),2,req( 9:12));
!!$     call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
!!$
!!$     call ighost_z(ixheight(:,:,:,1),2,req( 1: 4));  
!!$     call ighost_z(ixheight(:,:,:,2),2,req( 5: 8)); 
!!$     call ighost_z(ixheight(:,:,:,3),2,req( 9:12));
!!$     call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)

   end subroutine get_flags_rph


   subroutine get_all_heights_rph()
     include 'mpif.h'
     integer :: direction, ierr
     integer :: req(48),sta(MPI_STATUS_SIZE,48)
     if(.not.st_initialized) call initialize_surface_tension()
     st_initialized=.true.

!!$     call ghost_x(height(:,:,:,1),2,req( 1: 4));  
!!$     call ghost_x(height(:,:,:,2),2,req( 5: 8)); 
!!$     call ghost_x(height(:,:,:,3),2,req( 9:12));
!!$     call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
!!$
!!$     call ghost_y(height(:,:,:,1),2,req( 1: 4));  
!!$     call ghost_y(height(:,:,:,2),2,req( 5: 8)); 
!!$     call ghost_y(height(:,:,:,3),2,req( 9:12));
!!$     call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
!!$
!!$     call ghost_z(height(:,:,:,1),2,req( 1: 4));  
!!$     call ghost_z(height(:,:,:,2),2,req( 5: 8)); 
!!$     call ghost_z(height(:,:,:,3),2,req( 9:12));
!!$     call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
!!$
!!$     call ighost_x(ixheight(:,:,:,1),2,req( 1: 4));  
!!$     call ighost_x(ixheight(:,:,:,2),2,req( 5: 8)); 
!!$     call ighost_x(ixheight(:,:,:,3),2,req( 9:12));
!!$     call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
!!$
!!$     call ighost_y(ixheight(:,:,:,1),2,req( 1: 4));  
!!$     call ighost_y(ixheight(:,:,:,2),2,req( 5: 8)); 
!!$     call ighost_y(ixheight(:,:,:,3),2,req( 9:12));
!!$     call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
!!$
!!$     call ighost_z(ixheight(:,:,:,1),2,req( 1: 4));  
!!$     call ighost_z(ixheight(:,:,:,2),2,req( 5: 8)); 
!!$     call ighost_z(ixheight(:,:,:,3),2,req( 9:12));
!!$     call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)

     do direction=1,3
        call get_heights_rph(direction)
     enddo

     call ghost_x(height(:,:,:,1),2,req(1: 4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr) 
     call ghost_x(height(:,:,:,2),2,req(1: 4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)  
     call ghost_x(height(:,:,:,3),2,req(1: 4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr) 


     call ghost_y(height(:,:,:,1),2,req(1: 4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)   
     call ghost_y(height(:,:,:,2),2,req(1: 4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)  
     call ghost_y(height(:,:,:,3),2,req(1: 4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr) 


     call ghost_z(height(:,:,:,1),2,req(1: 4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)  
     call ghost_z(height(:,:,:,2),2,req(1: 4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr) 
     call ghost_z(height(:,:,:,3),2,req(1: 4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)


     call ighost_x(ixheight(:,:,:,1),2,req(1: 4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)  
     call ighost_x(ixheight(:,:,:,2),2,req(1: 4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr) 
     call ighost_x(ixheight(:,:,:,3),2,req(1: 4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)


     call ighost_y(ixheight(:,:,:,1),2,req(1: 4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)  
     call ighost_y(ixheight(:,:,:,2),2,req(1: 4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr) 
     call ighost_y(ixheight(:,:,:,3),2,req(1: 4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)

     call ighost_z(ixheight(:,:,:,1),2,req(1: 4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)  
     call ighost_z(ixheight(:,:,:,2),2,req(1: 4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr) 
     call ighost_z(ixheight(:,:,:,3),2,req(1: 4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
 
   end subroutine get_all_heights_rph



   subroutine get_heights_rph(direction)
     integer, intent(in) :: direction
     integer :: index
     logical :: base_not_found, bottom_n_found, bottom_p_found, top_n_found, top_p_found
     real(8) :: height_p, height_n
     integer, parameter  :: NDEPTH=2
     integer :: si,sj,sk
     integer :: i,j,k
     integer :: i0,j0,k0
     integer :: nd, ierr
     ! NDEPTH is the depth of layers tested above or below the reference cell. 
     ! including the central layer
     ! NDEPTH*2 + 1 = 5 means a 5 x 3 stencil. 
     !
     !  height_p : height for a normal pointing up (reference phase under the other phase)
     !  height_n : height for a normal pointing down (reference phase above the other phase)
    
     real(8) :: cc, cp, cm, cpm, dc, res, sumcc
     !     real(8) :: 

     integer :: icell, np, nm, q, qs, npmax, nmmax

     logical :: inp, inm, inpm, icpm

   
     dc = 1.d-8 ! tolerance for a full or empty cell
     npmax = 2
     nmmax = 2

     if (direction .eq. 1) then
     do k=ks,ke
        do j=js,je
           do i=is,ie

                 if(vof_flag(i,j,k).eq.2) then
                    res = mod(ixheight(i,j,k,1),10)
                 endif
                 
                 if (res.eq.0) then

                       np = 0
                       sumcc = 0.d0
                       inp = .true.
                       cp = 0.d0 ; cm = 0.d0
!                       npmax = min(NDEPTH, ie-i)
!                       nmmax = min(NDEPTH, i-1)
                       
                       do while (inp.and.(np.lt.npmax) ) 
                          cc =  cvof(i+np+1,j,k)
                          sumcc = sumcc + cc
                          if ( (cc .le. dc).or.(cc.ge.(1.-dc) )) then
                             inp = .false.
                             cp = cc
                          end if
                          np = np+1
                       end do
                       
                       inpm = .not.inp
                       nm = 0
                       inm = .true.
                       
                       do while (inpm.and.(nm.lt.nmmax) ) 
                          cc =  cvof(i-nm-1,j,k)
                          sumcc = sumcc + cc
                          if ( (cc .le. dc).or.(cc.ge.(1.-dc) )) then
                             inpm = .false.
                             inm = .false.
                             cm = cc
                          endif
                          nm = nm+1
                       end do
                    
                       cpm = cp + cm
                    
                       icpm = .false.
                       if ( (cpm.gt. 0.9).and.(cpm.lt.1.1)) then
                          icpm = .true.
                       endif

                       inpm = .not.(inp.and.inm)
                    
                       if ( inpm.and.icpm ) then
                          sumcc = sumcc + cvof(i,j,k)
                          icell = int(sumcc)
                          
                          if ( cm.gt.cp ) then
                             ixheight(i-nm+icell,j,k,1) = 1
                             qs = icell - 1
                             do q=1,qs
                                if ( ixheight(i-nm+q,j,k,1).eq.0 ) then
                                   ixheight(i-nm+q,j,k,1)=-1
                                endif
                             enddo
                          
                             qs = np+nm-1-icell
                             
                             do q=1,qs
                                if ( ixheight(i-nm+icell+q,j,k,1).eq.0 ) then
                                   ixheight(i-nm+icell+q,j,k,1)=-1
                                endif
                             enddo
                             
                             height(i-nm+icell,j,k,1) = sumcc - icell - 0.5d0

!!$	      xx[1] = (i-nl-1.)*h;
!!$	      yy[1] = (j-0.5)*h;
!!$	      xx[2] = (i-nl-1.+sumcc)*h;
!!$	      yy[2] = (j-0.5)*h;
!!$	      segments(fp,xx,yy,1,2,1,4,1,0,0);
!!$	      }
                             
                          else
                          
                             ixheight(i+np-icell,j,k,1) = 11
                             qs = icell -1
                             
                             do q=1,qs
                                if ( ixheight(i+np-q,j,k,1).eq.0 ) then
                                   ixheight(i+np-q,j,k,1) = -1
                                endif
                             enddo
                             qs = np+nm-1-icell
                           
                             do q=1,qs
                                if (ixheight(i+np-icell-q,j,k,1).eq.0 ) then
                                   ixheight(i+np-icell-q,j,k,1) = -1
                                endif
                             enddo
                           
                             height(i+np-icell,j,k,1) = 0.5d0 - (sumcc-icell)

!!$ 
!!$	      {                       /* next 5 lines: red lines in the graph */
!!$	      xx[1] = (i+nr)*h;
!!$	      yy[1] = (j-0.5)*h;
!!$	      xx[2] = (i+nr-sumcc)*h;
!!$	      yy[2] = (j-0.5)*h;
!!$	      segments(fp,xx,yy,1,2,1,4,1,0,0);
!!$	      }

                          endif
                       endif
                    endif

                 enddo
              enddo
           enddo

        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! end of x - direction 1

     if (direction .eq. 2) then

     do k=ks,ke
        do j=js,je
           do i=is,ie

                 if(vof_flag(i,j,k).eq.2) then
                    res = mod(ixheight(i,j,k,2),10)
                 endif
                 
                 if (res.eq.0) then

                       np = 0
                       sumcc = 0.d0
                       inp = .true.
                       cp = 0.d0 ; cm = 0.d0
!                       npmax = min(NDEPTH, je-j)
!                       nmmax = min(NDEPTH, j-1)
                       
                       do while (inp.and.(np.lt.npmax) ) 
                          cc =  cvof(i,j+np+1,k)
                          sumcc = sumcc + cc
                          if ( (cc .le. dc).or.(cc.ge.(1.-dc) )) then
                             inp = .false.
                             cp = cc
                          end if
                          np = np+1
                       end do
                       
                       inpm = .not.inp
                       nm = 0
                       inm = .true.
                       
                       do while (inpm.and.(nm.lt.nmmax) ) 
                          cc =  cvof(i,j-nm-1,k)
                          sumcc = sumcc + cc
                          if ( (cc .le. dc).or.(cc.ge.(1.-dc) )) then
                             inpm = .false.
                             inm = .false.
                             cm = cc
                          endif
                          nm = nm+1
                       end do
                    
                       cpm = cp + cm
                    
                       icpm = .false.
                       if ( (cpm.gt. 0.9).and.(cpm.lt.1.1)) then
                          icpm = .true.
                       endif

                       inpm = .not.(inp.and.inm)
                    
                       if ( inpm.and.icpm ) then
                          sumcc = sumcc + cvof(i,j,k)
                          icell = int(sumcc)
                          
                          if ( cm.gt.cp ) then
                             ixheight(i,j-nm+icell,k,2) = 2
                             qs = icell - 1
                             do q=1,qs
                                if ( ixheight(i,j-nm+q,k,2).eq.0 ) then
                                   ixheight(i,j-nm+q,k,2)=-2
                                endif
                             enddo
                          
                             qs = np+nm-1-icell
                             
                             do q=1,qs
                                if ( ixheight(i,j-nm+icell+q,k,2).eq.0 ) then
                                   ixheight(i,j-nm+icell+q,k,2)=-2
                                endif
                             enddo
                             
                             height(i,j-nm+icell,k,2) = sumcc - icell - 0.5d0

!!$	      xx[1] = (i-nl-1.)*h;
!!$	      yy[1] = (j-0.5)*h;
!!$	      xx[2] = (i-nl-1.+sumcc)*h;
!!$	      yy[2] = (j-0.5)*h;
!!$	      segments(fp,xx,yy,1,2,1,4,1,0,0);
!!$	      }
                             
                          else
                          
                             ixheight(i,j+np-icell,k,2) = 22
                             qs = icell -1
                             
                             do q=1,qs
                                if ( ixheight(i,j+np-q,k,2).eq.0 ) then
                                   ixheight(i,j+np-q,k,2) = -2
                                endif
                             enddo
                             qs = np+nm-1-icell
                           
                             do q=1,qs
                                if (ixheight(i,j+np-icell-q,k,2).eq.0 ) then
                                   ixheight(i,j+np-icell-q,k,2) = -2
                                endif
                             enddo
                           
                             height(i,j+np-icell,k,2) = 0.5d0 - (sumcc-icell)

!!$ 
!!$	      {                       /* next 5 lines: red lines in the graph */
!!$	      xx[1] = (i+nr)*h;
!!$	      yy[1] = (j-0.5)*h;
!!$	      xx[2] = (i+nr-sumcc)*h;
!!$	      yy[2] = (j-0.5)*h;
!!$	      segments(fp,xx,yy,1,2,1,4,1,0,0);
!!$	      }

                          endif
                       endif
                    endif

                 enddo
              enddo
           enddo
           
           end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  end of y - direction 2

     if (direction .eq. 3) then 
       if (rank==2) then
           write(*,*) 'pippa:'
        endif
     do k=ks,ke
        do j=js,je
           do i=is,ie

                 if(vof_flag(i,j,k).eq.2) then
                    res = mod(ixheight(i,j,k,3),10)
                 endif
                 
                 if (res.eq.0) then

                       np = 0
                       sumcc = 0.d0
                       inp = .true.
                       cp = 0.d0 ; cm = 0.d0
!                       npmax = min(NDEPTH, ke-k)
!                       nmmax = min(NDEPTH, k-1)
                       
                       do while (inp.and.(np.lt.npmax) ) 
                          cc =  cvof(i,j,k+np+1)
                          sumcc = sumcc + cc
                          if ( (cc .le. dc).or.(cc.ge.(1.-dc) )) then
                             inp = .false.
                             cp = cc
                          end if
                          np = np+1
                       end do
                       
                       inpm = .not.inp
                       nm = 0
                       inm = .true.
                       
                       do while (inpm.and.(nm.lt.nmmax) ) 
                          cc =  cvof(i,j,k-nm-1)
                          sumcc = sumcc + cc
                          if ( (cc .le. dc).or.(cc.ge.(1.-dc) )) then
                             inpm = .false.
                             inm = .false.
                             cm = cc
                          endif
                          nm = nm+1
                       end do
                    
                       cpm = cp + cm
                    
                       icpm = .false.
                       if ( (cpm.gt. 0.9).and.(cpm.lt.1.1)) then
                          icpm = .true.
                       endif

                       inpm = .not.(inp.and.inm)
                    
                       if ( inpm.and.icpm ) then
                          sumcc = sumcc + cvof(i,j,k)
                          icell = int(sumcc)
                          
                          if ( cm.gt.cp ) then
                             ixheight(i,j,k-nm+icell,3) = 3
                             qs = icell - 1
                             do q=1,qs
                                if ( ixheight(i,j,k-nm+q,3).eq.0 ) then
                                   ixheight(i,j,k-nm+q,3)=-3
                                endif
                             enddo
                          
                             qs = np+nm-1-icell
                             
                             do q=1,qs
                                if ( ixheight(i,j,k-nm+icell+q,3).eq.0 ) then
                                   ixheight(i,j,k-nm+icell+q,3)=-3
                                endif
                             enddo
                             
                             height(i,j,k-nm+icell,3) = sumcc - icell - 0.5d0

!!$	      xx[1] = (i-nl-1.)*h;
!!$	      yy[1] = (j-0.5)*h;
!!$	      xx[2] = (i-nl-1.+sumcc)*h;
!!$	      yy[2] = (j-0.5)*h;
!!$	      segments(fp,xx,yy,1,2,1,4,1,0,0);
!!$	      }
                             
                          else
                          
                             ixheight(i,j,k+np-icell,3) = 33
                             qs = icell -1
                             
                             do q=1,qs
                                if ( ixheight(i,j,k+np-q,3).eq.0 ) then
                                   ixheight(i,j,k+np-q,3) = -3
                                endif
                             enddo
                             qs = np+nm-1-icell
                           
                             do q=1,qs
                                if (ixheight(i,j,k+np-icell-q,3).eq.0 ) then
                                   ixheight(i,j,k+np-icell-q,3) = -3
                                endif
                             enddo
                           
                             height(i,j,k+np-icell,3) = 0.5d0 - (sumcc-icell)

!!$ 
!!$	      {                       /* next 5 lines: red lines in the graph */
!!$	      xx[1] = (i+nr)*h;
!!$	      yy[1] = (j-0.5)*h;
!!$	      xx[2] = (i+nr-sumcc)*h;
!!$	      yy[2] = (j-0.5)*h;
!!$	      segments(fp,xx,yy,1,2,1,4,1,0,0);
!!$	      }

                          endif
                       endif
                    endif

                 enddo
              enddo
           enddo


        if (rank==0) then
           i=8
           k=6
           j = ny/2          
           write(*,*)  'pipc',i,k,cvof(i,j,k-2),cvof(i,j,k-1),cvof(i,j,k),cvof(i,j,k+1),cvof(i,j,k+2)
           write(*,*)  'pipf',i,k,vof_flag(i,j,k-2),vof_flag(i,j,k-1),vof_flag(i,j,k),vof_flag(i,j,k+1),vof_flag(i,j,k+2)
           write(*,*)  'piph',i,k,height(i,j,k,1),height(i,j,k,3)
!           write(*,*)  'pip',i,k,ixheight(i,j,k-2,3),ixheight(i,j,k-1,3),ixheight(i,j,k,3),ixheight(i,j,k+1,3),ixheight(i,j,k+2,3)
        endif
        endif

        index =3
        k = (ke-ks)/2+3
        j = (je-js)/2
        do i=is,ie
           print *, x(i), height(i,j,k-3,index),height(i,j,k-2,index),height(i,j,k-1,index) &
                ,height(i,j,k,index),height(i,j,k+1,index),height(i,j,k+2,index),height(i,j,k+3,index)
           print *, x(i), ixheight(i,j,k-3,index),ixheight(i,j,k-2,index),ixheight(i,j,k-1,index) &
                ,ixheight(i,j,k,index),ixheight(i,j,k+1,index),ixheight(i,j,k+2,index),ixheight(i,j,k+3,index)
        enddo
      end subroutine get_heights_rph

   ! at fixed j,k look at heights (index, p:5, n:6)
   subroutine output_heights_rph()
     implicit none
     integer i,j,k,d,index
     real(8) h, th

     j = ny/2

     OPEN(UNIT=89,FILE=TRIM(out_path)//'/height-'//TRIM(int2text(rank,padding))//'.txt')
     OPEN(UNIT=90,FILE=TRIM(out_path)//'/reference-'//TRIM(int2text(rank,padding))//'.txt')

     !     write(89,*) " "
     !     write(89,*) " p heights at  z = L/2 as a function of x"
     !     write(89,*) " "
 
       do k=kmin,kmax
         do i=imin,imax
           th = (- z(k) + zlength/2.d0)/dx(i) + A_h*cos(2.*3.14159*x(i)/xlength)
!           write(90,100) x(i),th
           write(89,100) i,k,height(i,j,k,1),height(i,j,k,3),ixheight(i,j,k,1),ixheight(i,j,k,3)
      enddo
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

100  format(2(I3),2(f17.12),2(I3))
101  format(f24.16,A2)
103  format(7(f16.8),6(I2))
104  format(16(f5.1,' '))
1041 format(I3,16(f5.1,' '))
105  format(16(I2,' '))
1051 format(I3,16(I2,' '))

     close(89)
     close(90)

   end subroutine output_heights_rph
 
end module module_surface_tension

