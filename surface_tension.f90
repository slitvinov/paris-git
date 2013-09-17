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
  real(8), dimension(:,:,:,:), allocatable :: height ! 

   ! 4th index: 1 for positive height in x, 2 for negative height in x
   !            3 for positive height in y, 4 for negative height in y, 
   !              etc... 
  integer, dimension(:,:,:,:), allocatable :: ixheight ! HF flags for rph (Ruben-Phil) routines
  logical :: st_initialized = .false.
  logical :: recomputenormals = .true.
contains
!=================================================================================================
  subroutine initialize_surface_tension()
    implicit none
    if(.not.recomputenormals) then
       allocate(n1(imin:imax,jmin:jmax,kmin:kmax), n2(imin:imax,jmin:jmax,kmin:kmax),  &
               n3(imin:imax,jmin:jmax,kmin:kmax))
    endif
    allocate(height(imin:imax,jmin:jmax,kmin:kmax,6))
      if(nx.ge.500000.or.ny.gt.500000.or.nz.gt.500000) then
         stop 'nx too large'
      endif
      height = 2.d6
      st_initialized=.true.
   end subroutine initialize_surface_tension
!=================================================================================================
! 
!  Put normals in a common array. Absolutely not sure this is efficient
!
!=================================================================================================
   subroutine get_normals()
     implicit none
     real(8) :: stencil3x3(-1:1,-1:1,-1:1)
     integer :: i,j,k
     integer :: i0,j0,k0
     real(8) :: mxyz(3)
     if(recomputenormals) call pariserror("normals not allocated")
     if(.not.st_initialized) call initialize_surface_tension()

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
!=================================================================================================
!
! the core of HF computation
!
!=================================================================================================
   subroutine get_all_heights
     implicit none
     include 'mpif.h'
     integer :: direction, ierr
     integer :: req(24),sta(MPI_STATUS_SIZE,24)
     if(.not.st_initialized) call initialize_surface_tension()

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
!=================================================================================================
! 
!   the actual HF
! 
!=================================================================================================
   subroutine get_heights(direction)
     implicit none
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
     integer :: direction=2
     k = nz/2 + 2  ! +2 because of ghost layers
     j = ny/2 + 2  ! +2 because of ghost layers

     if(k<ks.or.k>ke) return
     if(j<js.or.j>je) return

     OPEN(UNIT=89,FILE=TRIM(out_path)//'/height-'//TRIM(int2text(rank,padding))//'.txt')
     OPEN(UNIT=90,FILE=TRIM(out_path)//'/reference-'//TRIM(int2text(rank,padding))//'.txt')

     index=5
     do i=is,ie
        th = wave2ls(x(i),y(j),z(k),direction)/dx(nx/2+2)
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
100  format(2(f24.16))
101  format(f24.16,A2)
     close(89)
     close(90)
   end subroutine output_heights
!=======================================================================================================
!   Check if we find nine heights in the neighboring cells, if not collect all heights in all directions
!=======================================================================================================
   subroutine get_local_heights(i1,j1,k1,try,nfound,indexfound,hloc,points,nposit)
      implicit none
      integer :: NDEPTH,NPOS,NOR
      parameter  (NOR=6) ! number of orientations
      parameter (NDEPTH=3,NPOS=NOR*27)
      integer, intent(in) :: i1(-1:1,-1:1,3), j1(-1:1,-1:1,3), k1(-1:1,-1:1,3)
      integer, intent(out) :: nfound, indexfound
      real(8), intent(out) :: hloc(-1:1,-1:1)   
      real(8), intent(out) :: points(NPOS,3)
      integer, intent(out) :: nposit
      integer, intent(inout) :: try(3)
      integer :: i0,j0,k0
      real(8) :: deltax
      integer :: d,s
      integer :: i,j,k,m,n,l
      integer :: index
      logical :: dirnotfound
      integer :: si,sj,sk
!
!  Loop over directions until an orientation with 9 heights is found. 
! 
      hloc = 2d6
      indexfound = -1
      nposit = 0 
      l=0
      dirnotfound=.true.
      deltax=dx(nx/2)
      do while (l.lt.3.and.dirnotfound)
         l = l+1
         d = try(l)    ! on entry, try(l) sorts the directions , closest to normal first. 
         if(d.eq.1) then
            si=1; sj=0; sk=0
         else if (d.eq.2) then
            si=0; sj=1; sk=0;
         else if (d.eq.3) then
            si=0; sj=0; sk=1
         else
            stop "bad direction"
         endif
         index = 2*(d-1)
         do while (index.lt.2*(d-1)+2.and.dirnotfound)   ! try both orientations of direction d
            ! fixme: change to try only orientation given by normal. 
            index = index + 1
            hloc = 2d6
            nfound = 0
            do m=-1,1 
               do n=-1,1
                  if(height(i1(m,n,d),j1(m,n,d),k1(m,n,d),index).lt.1d6) then  ! search at same level
                     ! one height found
                     hloc(m,n) = height(i1(m,n,d),j1(m,n,d),k1(m,n,d),index)
                     nfound = nfound + 1
                     nposit = nposit + 1
                     points(nposit,1) = x(i1(m,n,d))/deltax + hloc(m,n)*si
                     points(nposit,2) = y(j1(m,n,d))/deltax + hloc(m,n)*sj
                     points(nposit,3) = z(k1(m,n,d))/deltax + hloc(m,n)*sk
                  else
                     s = 0 
                     do while(s.lt.NDEPTH) ! search at other levels
                        s = s + 1
                        if (height(i1(m,n,d)+si*s,j1(m,n,d)+sj*s,k1(m,n,d)+sk*s,index).lt.1d6) then
                           hloc(m,n) = height(i1(m,n,d)+si*s,j1(m,n,d)+sj*s,k1(m,n,d)+sk*s,index) + s
                           nfound = nfound + 1
                           nposit = nposit + 1
                           points(nposit,1) = x(i1(m,n,d)+si*s)/deltax + hloc(m,n)*si
                           points(nposit,2) = y(j1(m,n,d)+sj*s)/deltax + hloc(m,n)*sj
                           points(nposit,3) = z(k1(m,n,d)+sk*s)/deltax + hloc(m,n)*sk
                           s = NDEPTH  ! to exit loop
                        else if  (height(i1(m,n,d)-si*s,j1(m,n,d)-sj*s,k1(m,n,d)-sk*s,index).lt.1d6) then
                           hloc(m,n) = height(i1(m,n,d)-si*s,j1(m,n,d)-sj*s,k1(m,n,d)-sk*s,index) - s
                           nfound = nfound + 1
                           nposit = nposit + 1
                           points(nposit,1) = x(i1(m,n,d)+si*s)/deltax + hloc(m,n)*si
                           points(nposit,2) = y(j1(m,n,d)+sj*s)/deltax + hloc(m,n)*sj
                           points(nposit,3) = z(k1(m,n,d)+sk*s)/deltax + hloc(m,n)*sk
                           s = NDEPTH  ! to exit loop
                        endif
                     end do ! searcher at other levels
                  end if ! search at same level
               end do ! n
            end do ! m 
            if(nfound.eq.9) then
               dirnotfound = .false.
               indexfound = index
               ! on exit, let try(3) be the h direction found
               m=1
               n=1
               do while (m.le.3)
                  if(m.ne.d) then
                     try(n) = m
                     n=n+1
                  endif
                  m=m+1
               enddo
               try(3)=d  ! then exit
            end if ! nfound = 9
         end do ! index and dirnotfound
      end do ! d and dirnotfound
   end subroutine get_local_heights
!
   subroutine get_curvature(i0,j0,k0,kappa,indexCurv,nfound,nposit)
      implicit none
      integer, intent(in) :: i0,j0,k0
      real(8), intent(out) :: kappa  
      integer, intent(out) :: indexCurv, nfound
      integer, intent(out) :: nposit

      integer :: d,index,indexfound
      real(8) :: h(-1:1,-1:1),a(6)
      integer :: nCentroids

      integer :: m,n,l,i,j,k
      logical :: fit_success = .false.
      logical :: height_found
      integer :: i1(-1:1,-1:1,3), j1(-1:1,-1:1,3), k1(-1:1,-1:1,3),try(3)
      integer :: si,sj,sk,hsign,s,c(3)
      
      integer :: NDEPTH,NOR,NPOS
      parameter (NOR=6) ! number of orientations
      parameter (NDEPTH=3,NPOS=NOR*27)
      real(8) :: points(NPOS,3)
      real(8) :: xfit(NPOS),yfit(NPOS),hfit(NPOS),fit(NPOS,3)
      real(8) :: centroid(3),mx(3),stencil3x3(-1:1,-1:1,-1:1)

      call map3x3in2x2(i1,j1,k1,i0,j0,k0)
!   define in which order directions will be tried 
!   direction closest to normal first
      if(recomputenormals) then
         do m=-1,1; do n=-1,1; do l=-1,1
            stencil3x3(m,n,l) = cvof(i0+m,j0+n,k0+l)
         enddo;enddo;enddo
         call fd32(stencil3x3,mx)
      else
         mx(1) = n1(i0,j0,k0)      
         mx(2) = n2(i0,j0,k0)      
         mx(3) = n3(i0,j0,k0)
      endif
      call orientation(mx,try)
      call get_local_heights(i1,j1,k1,try,nfound,indexfound,h,points,nposit)
      kappa = 0.d0
      
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
         ! This stops the code in case kappa becomes NaN.
         if(kappa.ne.kappa) call pariserror("HF9: Invalid Curvature")               
         ! if more than six independent heights found in all directions 
         return
      else 
         nfound = - ind_pos(points,nposit) 
      endif ! nfound == 9
      if ( (-nfound) > 6 )  then  ! 6 points to avoid special 2D degeneracy. 
         xfit=points(:,try(2))
         yfit=points(:,try(3))
         hfit=points(:,try(1))
         ! fit over all positions, not only independent ones. 
         call parabola_fit(xfit,yfit,hfit,nposit,a,fit_success) 
         if(.not.fit_success) call pariserror("no fit success")
         kappa = 2.d0*(a(1)*(1.d0+a(5)*a(5)) + a(2)*(1.d0+a(4)*a(4)) - a(3)*a(4)*a(5)) &
               /(1.d0+a(4)*a(4)+a(5)*a(5))**(1.5d0)
         indexCurv = indexfound
         ! This stops the code in case kappa becomes NaN.
         if(kappa.ne.kappa) call pariserror("HF6: Invalid Curvature")
         nposit=0
         return
      else  ! ind_pos <= 5
         ! Find all centroids in 3**3
         ! use direction closest to normal
         nfound = -100 + nfound  ! encode number of independent positions into nfound
         nposit=1
         do m=-1,1; do n=-1,1; do l=-1,1
            i=i0+m
            j=j0+n
            k=k0+l
            c(1)=m
            c(2)=n
            c(3)=l
            if(vof_flag(i,j,k) == 2) then
               height_found=.false.
               do index=1,6
                  if(height(i,j,k,index) < 1d6) then   !  Use height if it exists
                     height_found=.true.
                     d = (index-1)/2 + 1
                     do s=1,3  ! positions relative to center if i0,j0,k0 in cell units
                        if(d==s) then
                           fit(nposit,s) = c(s)
                        else
                           fit(nposit,s) = c(s) + height(i,j,k,index)
                        endif
                     enddo ! do s
                     nposit = nposit + 1
                  endif ! height exists
               enddo ! do index
               if(.not.height_found) then
                  call FindCutAreaCentroid(i,j,k,centroid)
                  do s=1,3 
                     fit(nposit,s) = centroid(s) + c(s)
                  end do
                  nposit = nposit + 1
               endif 
            endif ! vof_flag
         enddo; enddo; enddo ! do m
         ! arrange coordinates so height direction is closest to normal
         ! try(:) array contains direction closest to normal first

         xfit=fit(:,try(2))
         yfit=fit(:,try(3))
         hfit=fit(:,try(1))
         call parabola_fit(xfit,yfit,hfit,nposit,a,fit_success)
         kappa = 2.d0*(a(1)*(1.d0+a(5)*a(5)) + a(2)*(1.d0+a(4)*a(4)) - a(3)*a(4)*a(5)) &
              /(1.d0+a(4)*a(4)+a(5)*a(5))**(1.5d0)
         ! This stops the code in case kappa becomes NaN.
         if(kappa.ne.kappa) then
            call pariserror("PF6: Invalid Curvature")  
         endif
      end if ! -nfound > 5
   end subroutine get_curvature

   subroutine output_curvature()
      implicit none      
      integer :: i,j,k,indexCurv,l,m,n
      integer :: ib,dir,try(3)
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
      integer :: i1(-1:1,-1:1,3), j1(-1:1,-1:1,3), k1(-1:1,-1:1,3)
      real(8) :: stencil3x3(-1:1,-1:1,-1:1),mx(3)
      integer :: NOR,NPOS,nposit
      parameter  (NOR=6) ! number of orientations
      parameter (NPOS=NOR*27)
      real(8) :: points(NPOS,3)
      real(8), parameter :: PI= 3.14159265359d0

      OPEN(UNIT=89,FILE=TRIM(out_path)//'/curvature-'//TRIM(int2text(rank,padding))//'.txt')
      OPEN(UNIT=90,FILE=TRIM(out_path)//'/reference-'//TRIM(int2text(rank,padding))//'.txt')
      OPEN(UNIT=91,FILE=TRIM(out_path)//'/bigerror-'//TRIM(int2text(rank,padding))//'.txt')
      OPEN(UNIT=92,FILE=TRIM(out_path)//'/debug-'//TRIM(int2text(rank,padding))//'.txt')
      ib = 1
      kappamin = 1d20
      kappamax = -1d20
      if ( test_curvature ) then 
         kappa_exact = 2.d0/rad(ib)
         do i=is,ie; do j=js,je; do k=ks,ke
            ! find curvature only for cut cells
            if (vof_flag(i,j,k) == 2 ) then 
               call get_curvature(i,j,k,kappa,indexCurv,nfound,nposit)
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
         S2_err_K=0.d0
         Lm_err_K=0.d0

         if(1==0) then ! begin remove h to h_exact comparison
            S2_err_h=0.d0;S2_err_hp=0.d0;S2_err_hm=0.d0;S2_err_dh=0.d0;S2_err_d2h=0.d0;
            Lm_err_h=0.d0;Lm_err_hp=0.d0;Lm_err_hm=0.d0;Lm_err_dh=0.d0;Lm_err_d2h=0.d0;
         endif
         
         k = (Nz+4)/2

         do j=js,je
            write(92,'(I2)',advance='no') j
            do i=is,ie
               if (vof_flag(i,j,k) == 2) then 
                  call get_curvature(i,j,k,kappa,indexCurv,nfound,nposit)
                  write(92,'(1X,E10.2)',advance='no') kappa
               else
                  write(92,'(2X,"*",I1,"*",2X)',advance='no') vof_flag(i,j,k)
               endif
            enddo
            write(92,*) "  "
         enddo
         write(92,*) "  "
         
         do j=js,je
            write(92,'(I2)',advance='no') j
            do i=is,ie
               if (vof_flag(i,j,k) == 2) then 
                  call get_curvature(i,j,k,kappa,indexCurv,nfound,nposit)
!                  write(92,'(1X,"<",I1,"-",I1,">")',advance='no') nfound,nindepend
                  write(92,'(I4,1X)',advance='no') nfound
               else
                  write(92,'(1X,"*",I1,"*",1X)',advance='no') vof_flag(i,j,k)
               endif
            enddo
            write(92,*) "  "
         enddo
         write(92,*) "  "

        do j=js,je
            write(92,'(I2)',advance='no') j
            do i=is,ie
               if (vof_flag(i,j,k) == 2) then 
                  call get_curvature(i,j,k,kappa,indexCurv,nfound,nposit)
!                  write(92,'(1X,"<",I1,"-",I1,">")',advance='no') nfound,nindepend
                  write(92,'(I4,1X)',advance='no') nposit
               else
                  write(92,'(1X,"*",I1,"*",1X)',advance='no') vof_flag(i,j,k)
               endif
            enddo
            write(92,*) "  "
         enddo
         write(92,*) "  "

         do n=1,6
            write(92,*) " "
            write(92,*) "orientation",n
            write(92,*) " "
            do i=is,ie
               write(92,'(3X,I1,3X)',advance='no') i
            enddo
            write(92,*) " "
            do j=je,js,-1
               write(92,'(I2)',advance='no') j
               do i=is,ie
                  if(height(i,j,k,n)>1d6) then
                     write(92,'(3X,".",3X)',advance='no')
                  else
                     write(92,'(1X,E10.2)',advance='no') height(i,j,k,n)
                  endif
               enddo
               write(92,*) " "
            enddo
         enddo
         write(92,*) " "
 
         do i=is,ie; do j=js,je
            if (vof_flag(i,j,k) == 2) then 
               call get_curvature(i,j,k,kappa,indexCurv,nfound,nposit)
               ! This stops the code in case kappa becomes NaN.
               if(kappa.ne.kappa) call pariserror("Invalid Curvature")  
               if(nfound==-1.or.abs(kappa)<1d-4) then
                  write(6,*) "i,j,k,nfound,nindepend,kappa ",i,j,k,nfound,nindepend,kappa
                  call pariserror("OC: curvature not found")
               else
                  kappa = kappa*dble(Nx)
                  kappamax = max(ABS(kappa),kappamax)
                  kappamin = min(ABS(kappa),kappamin)
                  angle = atan2(y(j)-yc(ib),x(i)-xc(ib))/PI*180.d0
                  write(89,*) angle,ABS(kappa)
                  write(90,*) angle,kappa_exact
                  
                  err_K    = ABS(ABS(kappa)-kappa_exact)/kappa_exact
                  S2_err_K    = S2_err_K  + err_K  **2.d0
                  Lm_err_K    = MAX(Lm_err_K,   err_K) 
                  sumCount = sumCount + 1

                  if(1==0) then ! remove h to h_exact comparison
                     call CalExactHeight_Circle(x(i),y(j),1.d0/dble(Nx),1.d0/dble(Ny),&
                          xc(ib),yc(ib),rad(ib),indexCurv,hex,hpex,hmex,dhex,d2hex)
                     call map3x3in2x2(i1,j1,k1,i,j,k)
                     if(recomputenormals) then
                        do m=-1,1; do n=-1,1; do l=-1,1
                           stencil3x3(m,n,l) = cvof(i+m,j+n,k+l)
                        enddo;enddo;enddo
                        call fd32(stencil3x3,mx)
                     else
                        mx(1) = n1(i,j,k)      
                        mx(2) = n2(i,j,k)      
                        mx(3) = n3(i,j,k)
                     endif
                     call orientation(mx,try)

                     call get_local_heights(i1,j1,k1,try,nfound,indexCurv,hloc,points,nposit)
                     hnum  = hloc( 0,0)/dble(Nx)
                     hpnum = hloc( 1,0)/dble(Nx)
                     hmnum = hloc(-1,0)/dble(Nx)
                     if(hnum+hpnum+hmnum<1d6) then 
                        dhnum = (hpnum-hmnum)/(2.d0/dble(Nx))
                        d2hnum = (hpnum-2.d0*hnum+hmnum)/(1.d0/dble(Nx))**2.d0
                        err_h    = ABS(hnum   - hex )
                        err_hp   = ABS(hpnum  - hpex)
                        err_hm   = ABS(hmnum  - hmex)
                        err_dh   = ABS(dhnum  - dhex)
                        err_d2h  = ABS(d2hnum - d2hex)/ABS(d2hex)
                        write(91,'(7(E15.8,1X),6(I5,1X))') angle,err_h,err_hp,err_hm,& 
                             err_dh,err_d2h,err_K,indexCurv,i,j,nfound,nindepend
                        !@@                     sumCount = sumCount + 1
                        S2_err_h    = S2_err_h  + err_h  **2.d0
                        S2_err_hp   = S2_err_hp + err_hp **2.d0
                        S2_err_hm   = S2_err_hm + err_hm **2.d0
                        S2_err_dh   = S2_err_dh + err_dh **2.d0
                        S2_err_d2h  = S2_err_d2h+ err_d2h**2.d0
                        
                        Lm_err_h    = MAX(Lm_err_h,   err_h ) 
                        Lm_err_hp   = MAX(Lm_err_hp,  err_hp) 
                        Lm_err_hm   = MAX(Lm_err_hm,  err_hm) 
                        Lm_err_dh   = MAX(Lm_err_dh,  err_dh) 
                        Lm_err_d2h  = MAX(Lm_err_d2h, err_d2h) 
                     endif ! h < 1d6
                  end if ! end remove h to hexact
               endif ! valid curvature
            end if ! cvof(i,j,k)
         end do; end do
         
         if(1==0) then
            if(SumCount>0) then
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
            endif
         endif
         L2_err_K    = sqrt(S2_err_K)  /dble(sumCount)
         write(*,*) 'L2 Norm:'
         write(*,'(I5,I5,1X,(E15.8,1X))') Nx,rank,L2_err_K
         write(*,*) 'Linfty Norm:'
         write(*,'(I5,I5,1X,(E15.8,1X))') Nx,rank,Lm_err_K
      end if ! test_curvature
      write(*,*) 'max, min, and exact ABS(kappa)', kappamax, kappamin,kappa_exact
      write(*,*) 'max relative error', MAX(ABS(kappamax-kappa_exact), ABS(kappamin-kappa_exact))/kappa_exact
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

   subroutine parabola_fit(xfit,yfit,hfit,nposit,a,fit_success)
      implicit none

      real(8), intent(in)  :: xfit(nposit),yfit(nposit),hfit(nposit)
      real(8), intent(out) :: a(6)
      logical, intent(out) :: fit_success

      real(8) :: m(6,6), invm(6,6)
      real(8) :: rhs(6)
      integer :: ifit, im,jm, nposit
      logical :: inv_success
      real(8) :: x1,x2,x3,x4,y1,y2,y3,y4

      a = 0.d0

      ! evaluate the linear system for least-square fit
      m   = 0.d0
      rhs = 0.d0
      do ifit = 1, nposit
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
      real(8) :: dmx,dmy,dmz, mxyz(3),px,py,pz
      real(8) :: invx,invy,invz
      real(8) :: alpha, alphamax,al3d, minnormcomp
      real(8) :: stencil3x3(-1:1,-1:1,-1:1)

      ! find cut area centroid 
      !***
      !     (1) normal vector: dmx,dmy,dmz, and |dmx|+|dmy|+|dmz| = 1.
      !     (2) dmx,dmy,dmz>0 and record sign
      !     (3) get alpha;               
      !     (4) compute centroid with dmx,dmy,dmz and alpha;
      !     (5) transfer to local coordinate;
      !*(1)*

      if(recomputenormals) then
         do i0=-1,1; do j0=-1,1; do k0=-1,1
            stencil3x3(i0,j0,k0) = cvof(i+i0,j+j0,k+k0)
         enddo;enddo;enddo
         call mycs(stencil3x3,mxyz)
         dmx = mxyz(1)
         dmy = mxyz(2)
         dmz = mxyz(3)      
      else
         dmx = n1(i0,j0,k0)      
         dmy = n2(i0,j0,k0)      
         dmz = n3(i0,j0,k0)
      endif
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
      call PlaneAreaCenter(dmx,dmy,dmz,alpha,px,py,pz)
      !*(5)*
      ! trap NaNs
      if(px.ne.px) call pariserror("FCAC:invalid px")
      if(py.ne.py) call pariserror("FCAC:invalid py")
      if(pz.ne.pz) call pariserror("FCAC:invalid pz")

      ! rotate
      centroid(1) = px*invx
      centroid(2) = py*invy
      centroid(3) = pz*invz
      ! shift
      centroid(1) = centroid(1) - invx*0.5d0
      centroid(2) = centroid(2) - invy*0.5d0
      centroid(3) = centroid(3) - invz*0.5d0

!      dmx = invx*dmx
!      dmy = invy*dmy
!      dmz = invz*dmz

   end subroutine FindCutAreaCentroid
! 
!   Computes the centroid as in gerris
! 
! * Fills p with the position of the center of mass of the polygon
! * obtained by interseectiing the plane  (m,alpha).
! * with the reference cell.
!
!  assumptions: mx,my,mz > 0 and |mx| + |my| + |mz| = 1
!
   subroutine PlaneAreaCenter (mx,my,mz, alpha, px,py,pz)
     implicit none
     real(8), intent(in) :: mx,my,mz,alpha
     real(8), intent(out) :: px,py,pz
     real(8) :: nx,ny,qx,qy
     real(8) :: area,b,amax

     if(mx<0.d0.or.my<0.d0.or.mz<0.d0) call pariserror("invalid mx-my-mz")
     if(abs(mx+my+mz-1d0)>1d-4) call pariserror("invalid mx+my+mz")

     if (mx < 1d-4) then
        nx = my
        ny = mz
        call LineCenter (nx,ny, alpha, qx,qy)
        px = 0.5d0
        py = qx
        pz = qy
        return
     endif
     if (my < 1d-4) then
        nx = mz
        ny = mx
        call LineCenter (nx,ny, alpha, qx,qy)
        px = qy
        py = 0.5d0
        pz = qx
        return
     endif
     if (mz < 1d-4) then
        call LineCenter (mx,my, alpha, px,py)
        pz = 0.5
        return
     endif

     if (alpha < 0.d0 .or. alpha > 1.d0) call pariserror("PAC: invalid alpha")

     area = alpha*alpha
     pz = area*alpha
     py = area*alpha
     pz = area*alpha
     b = alpha - mx
     if (b > 0.) then
        area = area - b*b
        px = px - b*b*(2.*mx + alpha)
        py = py - b*b*b
        pz = pz - b*b*b
     endif
     b = alpha - my
     if (b > 0.) then
        area = area - b*b
        py = py - b*b*(2.*my + alpha)
        px = px - b*b*b
        pz = pz - b*b*b
     endif
     b = alpha - mz
     if (b > 0.) then
        area = area - b*b
        pz = pz - b*b*(2.*mz + alpha)
        px = px - b*b*b
        py = py - b*b*b
     endif

     amax = alpha - 1.d0
     b = amax + mx
     if (b > 0.) then
        area = area + b*b
        py = py + b*b*(2.*my + alpha - mz)
        pz = pz + b*b*(2.*mz + alpha - my)
        px = px + b*b*b
     endif
     b = amax + my
     if (b > 0.) then
        area = area + b*b
        px = px + b*b*(2.*mx + alpha - mz)
        pz = pz + b*b*(2.*mz + alpha - mx)
        py = py + b*b*b
     endif
     b = amax + mz
     if (b > 0.) then
        area = area + b*b
        px = px + b*b*(2.*mx + alpha - my)
        py = py + b*b*(2.*my + alpha - mx)
        pz = pz + b*b*b
     endif

     area  = 3.d0*area
     px = px/(area*mx)
     py = py/(area*my)
     pz = pz/(area*mz)

     call THRESHOLD (px)
     call THRESHOLD (py)
     call THRESHOLD (pz)

   end subroutine PlaneAreaCenter

!-------------------------------------------------------------------------------------------------------
   subroutine LineCenter (mx,my, alpha, px,py)
     implicit none
     real(8), intent(in) :: mx,my,alpha
     real(8), intent(out) :: px,py
      
     if (alpha <= 0.d0 .or. alpha >= 1.d0) call pariserror("LC: invalid alpha")

     if (mx < 1d-4) then
        px = 0.5
        py = alpha;
        return
     endif

     if (my < 1d-4) then
        py = 0.5;
        px = alpha
        return
     endif

     px = 0.; py = 0.

     if (alpha >= mx) then
        px = px +  1.
        py = py +  (alpha - mx)/my
     else
        px = px +  alpha/mx
     endif

     if (alpha >= my) then
        py = py +  1.
        px = px +  (alpha - my)/mx
     else
        py = py +  alpha/my
     endif

     px = px/2.
     py = py/2.

     call THRESHOLD (px)
     call THRESHOLD (py)

     if(px.ne.px) call pariserror("LAC:invalid px")
     if(py.ne.py) call pariserror("LAC:invalid px")

   end subroutine

!   direction closest to normal first
   subroutine orientation (m,c)
     implicit none
     real(8), intent(in) :: m(3)
     integer, intent(inout) :: c(3)
     integer :: i,j,tmp
     do i = 1,3
        c(i) = i 
     enddo
     do i = 1,2
        do j=1,3-i
           if(abs(m(c(j+1))) > abs(m(c(j)))) then
              tmp = c(j)
              c(j) = c(j+1)
              c(j+1) = tmp
           endif
        enddo
     enddo
   end subroutine orientation

   function ind_pos (points, n)
     implicit none
     integer :: ind_pos
     integer, intent(in) :: n
     real(8), intent(in) :: points(n,3)
     integer :: i,j,ni,c
     real(8) :: d2
     logical :: depends

     if (n < 2) then
        ind_pos = n
        return
     endif

     ni=1
     
     do j=2,n
        depends = .false.
        do i=1,j-1
           if(.not.depends) then
              d2 = 0d0
              do c=1,3
                 d2 = d2 + (points(i,c) - points(j,c))**2
              enddo
              depends = (d2 < 0.5**2)
           endif
        enddo
        if(.not.depends) ni = ni + 1
     enddo
     ind_pos = ni
   end function ind_pos

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

!!$     call ghost_x(height(:,:,:,1),2,req( 1: 4))  
!!$     call ghost_x(height(:,:,:,2),2,req( 5: 8)) 
!!$     call ghost_x(height(:,:,:,3),2,req( 9:12))
!!$     call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
!!$
!!$     call ghost_y(height(:,:,:,1),2,req( 1: 4))  
!!$     call ghost_y(height(:,:,:,2),2,req( 5: 8)) 
!!$     call ghost_y(height(:,:,:,3),2,req( 9:12))
!!$     call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
!!$
!!$     call ghost_z(height(:,:,:,1),2,req( 1: 4))  
!!$     call ghost_z(height(:,:,:,2),2,req( 5: 8)) 
!!$     call ghost_z(height(:,:,:,3),2,req( 9:12))
!!$     call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
!!$
!!$     call ighost_x(ixheight(:,:,:,1),2,req( 1: 4))  
!!$     call ighost_x(ixheight(:,:,:,2),2,req( 5: 8)) 
!!$     call ighost_x(ixheight(:,:,:,3),2,req( 9:12))
!!$     call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
!!$
!!$     call ighost_y(ixheight(:,:,:,1),2,req( 1: 4))  
!!$     call ighost_y(ixheight(:,:,:,2),2,req( 5: 8)) 
!!$     call ighost_y(ixheight(:,:,:,3),2,req( 9:12))
!!$     call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
!!$
!!$     call ighost_z(ixheight(:,:,:,1),2,req( 1: 4))  
!!$     call ighost_z(ixheight(:,:,:,2),2,req( 5: 8)) 
!!$     call ighost_z(ixheight(:,:,:,3),2,req( 9:12))
!!$     call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)

   end subroutine get_flags_rph


   subroutine get_all_heights_rph()
     include 'mpif.h'
     integer :: direction, ierr
     integer :: req(48),sta(MPI_STATUS_SIZE,48)
     if(.not.st_initialized) call initialize_surface_tension()
     st_initialized=.true.

!!$     call ghost_x(height(:,:,:,1),2,req( 1: 4))  
!!$     call ghost_x(height(:,:,:,2),2,req( 5: 8)) 
!!$     call ghost_x(height(:,:,:,3),2,req( 9:12))
!!$     call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
!!$
!!$     call ghost_y(height(:,:,:,1),2,req( 1: 4))  
!!$     call ghost_y(height(:,:,:,2),2,req( 5: 8)) 
!!$     call ghost_y(height(:,:,:,3),2,req( 9:12))
!!$     call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
!!$
!!$     call ghost_z(height(:,:,:,1),2,req( 1: 4))  
!!$     call ghost_z(height(:,:,:,2),2,req( 5: 8)) 
!!$     call ghost_z(height(:,:,:,3),2,req( 9:12))
!!$     call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
!!$
!!$     call ighost_x(ixheight(:,:,:,1),2,req( 1: 4))  
!!$     call ighost_x(ixheight(:,:,:,2),2,req( 5: 8)) 
!!$     call ighost_x(ixheight(:,:,:,3),2,req( 9:12))
!!$     call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
!!$
!!$     call ighost_y(ixheight(:,:,:,1),2,req( 1: 4))  
!!$     call ighost_y(ixheight(:,:,:,2),2,req( 5: 8)) 
!!$     call ighost_y(ixheight(:,:,:,3),2,req( 9:12))
!!$     call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
!!$
!!$     call ighost_z(ixheight(:,:,:,1),2,req( 1: 4))  
!!$     call ighost_z(ixheight(:,:,:,2),2,req( 5: 8)) 
!!$     call ighost_z(ixheight(:,:,:,3),2,req( 9:12))
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
!!$	      endif}
                             
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

