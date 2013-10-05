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
  integer, parameter :: NDEPTH=2
  integer, parameter :: NOR=6 ! number of orientations
  integer, parameter :: NPOS=NOR*27
  real(8), parameter :: EPS_GEOM = 1d-4
  real(8), dimension(:,:,:), allocatable :: n1,n2,n3 ! normals
  real(8), dimension(:,:,:,:), allocatable :: height ! 

  ! 4th index: 1 for normal vector pointing towards positive x "positive height", 
  ! 2 for "negative" height in x
  ! 3 for positive height in y, 4 for negative height in y, 
  !  etc... 
  integer, dimension(:,:,:,:), allocatable :: ixheight ! HF flags for rph (Ruben-Phil) routines
  logical :: st_initialized = .false.
  logical :: recomputenormals = .true.
  logical :: debug_curvature = .false.
  logical :: debug_ij55 = .false.
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
     integer :: si,sj,sk
     integer :: i,j,k
     integer :: i0,j0,k0
     integer :: nd
     ! NDEPTH is the depth of layers tested above or below the reference cell. 
     ! including the central layer
     ! NDEPTH*2 - 1 = 5 means a 5 x 3 stencil. 
     !
     !  Note the normal is - grad C
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
   subroutine get_local_heights(i1,j1,k1,mx,try,nfound,indexfound,hloc,points,nposit)
      implicit none
      integer, intent(in) :: i1(-1:1,-1:1,3), j1(-1:1,-1:1,3), k1(-1:1,-1:1,3)  ! i1(:,:,d) 3x3 plane rotated in direction d
      integer, intent(out) :: nfound, indexfound
      real(8), intent(in)  :: mx(3)
      real(8), intent(out) :: hloc(-1:1,-1:1)   
      real(8), intent(out) :: points(NPOS,3)
      integer, intent(out) :: nposit
      integer, intent(inout) :: try(3)
!      integer :: i0,j0,k0
      real(8) :: deltax
      integer :: d,s
      integer :: i,j,k,m,n,l
      integer :: index
      logical :: dirnotfound,heightnotfound
      integer :: si,sj,sk

      i = i1(0,0,1)
      j = j1(0,0,1)
      k = k1(0,0,1)
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
         index =  2*(d-1)+2
         if(mx(d).gt.0) index = 2*(d-1)+1
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
                     s = 1 
                     heightnotfound=.true.
                     do while(s.le.NDEPTH.and.heightnotfound) ! search at other levels
                         if (height(i1(m,n,d)+si*s,j1(m,n,d)+sj*s,k1(m,n,d)+sk*s,index).lt.1d6) then
                           hloc(m,n) = height(i1(m,n,d)+si*s,j1(m,n,d)+sj*s,k1(m,n,d)+sk*s,index) + s
                           nfound = nfound + 1
                           nposit = nposit + 1
                           points(nposit,1) = x(i1(m,n,d)+si*s)/deltax + hloc(m,n)*si
                           points(nposit,2) = y(j1(m,n,d)+sj*s)/deltax + hloc(m,n)*sj
                           points(nposit,3) = z(k1(m,n,d)+sk*s)/deltax + hloc(m,n)*sk
                           heightnotfound=.false.  ! to exit loop
                        else if  (height(i1(m,n,d)-si*s,j1(m,n,d)-sj*s,k1(m,n,d)-sk*s,index).lt.1d6) then
                           hloc(m,n) = height(i1(m,n,d)-si*s,j1(m,n,d)-sj*s,k1(m,n,d)-sk*s,index) - s
                           nfound = nfound + 1
                           nposit = nposit + 1
                           points(nposit,1) = x(i1(m,n,d)+si*s)/deltax + hloc(m,n)*si
                           points(nposit,2) = y(j1(m,n,d)+sj*s)/deltax + hloc(m,n)*sj
                           points(nposit,3) = z(k1(m,n,d)+sk*s)/deltax + hloc(m,n)*sk
                           heightnotfound=.false.  ! to exit loop
                        endif
                        s = s + 1
                    end do ! while s lt ndepth 
                  end if ! search at same level
               end do ! n
            end do ! m 
!            write(*,*) "GLH: end m,n: nposit,d,index,nfound ",nposit,d,index,nfound
            if(nfound.eq.9) then
               dirnotfound = .false.
               indexfound = index
               ! on exit, let try(1) be the h direction found
               m=1
               n=2
               do while (m.le.3)
                  if(m.ne.d) then
                     try(n) = m
                     n=n+1
                  endif
                  m=m+1
               enddo
               try(1)=d  ! then exit
!               write(*,*) "GLH: 9F: i,j,k,try,nposit,nfound ",i,j,k,try,nposit,nfound
               return
            end if ! nfound = 9
      end do ! d and dirnotfound
      indexfound =  2*(try(1)-1)+2
      if(mx(try(1)).gt.0) indexfound = 2*(try(1)-1)+1
      if(nposit.gt.NPOS) call pariserror("GLH: nposit")
!      write(*,*) "GLH: not F:i,j,k, try,nposit,nfound ",i,j,k,try,nposit,nfound

   end subroutine get_local_heights
!
   subroutine get_curvature(i0,j0,k0,kappa,indexCurv,nfound,nposit,a)
      implicit none
      integer, intent(in) :: i0,j0,k0
      real(8), intent(out) :: kappa,a(6)  
      integer, intent(out) :: indexCurv, nfound
      integer, intent(out) :: nposit

      integer :: indexfound,index,d
      real(8) :: h(-1:1,-1:1)
!      integer :: nCentroids

      integer :: m,n,l,i,j,k
      logical :: fit_success = .false.
      logical :: height_found
      integer :: i1(-1:1,-1:1,3), j1(-1:1,-1:1,3), k1(-1:1,-1:1,3),try(3)
      integer :: s,c(3)
      
      real(8) :: points(NPOS,3),origin(3)
      real(8) :: xfit(NPOS),yfit(NPOS),hfit(NPOS),fit(NPOS,3)
      real(8) :: centroid(3),mx(3),stencil3x3(-1:1,-1:1,-1:1)

      ! debug
      integer :: k00
      real(8) :: debugfit(-1:1,-1:1,2)
      real(8) :: deltax
      integer :: iem,jem
      real(8) :: x1,y1,xvec,yvec,xpoint(0:2),ypoint(0:2),pc(12,12,2),diff(0:2)
      if(debug_curvature) then
         k00 = (Nz+4)/2
         jem = je - 2
         iem = ie - 2
         deltax=dx(nx/2)
         do i=js,jem
            write(79,'(4(E15.8,1X))') xh(is),yh(i),xh(iem)-xh(is),0.d0
         enddo
         do i=is,iem
            write(79,'(4(E15.8,1X))') xh(i),yh(js),0.,yh(jem)-yh(js)
         enddo
      endif
      ! end debug

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
      call get_local_heights(i1,j1,k1,mx,try,nfound,indexfound,h,points,nposit)
      kappa = 0.d0
      ! if all nine heights found 
      if ( nfound == 9 ) then
!
!  h = a6  + a4 x + a5 y + a3 xy + a1 x**2 + a2 y**2
!
         a(1) = (h(1,0)-2.d0*h(0,0)+h(-1,0))/2.d0
         a(2) = (h(0,1)-2.d0*h(0,0)+h(0,-1))/2.d0
         a(3) = (h(1,1)-h(-1,1)-h(1,-1)+h(-1,-1))/4.d0
         a(4) = (h(1,0)-h(-1,0))/2.d0
         a(5) = (h(0,1)-h(0,-1))/2.d0
         kappa = 2.d0*(a(1)*(1.d0+a(5)*a(5)) + a(2)*(1.d0+a(4)*a(4)) - a(3)*a(4)*a(5)) &
               /(1.d0+a(4)*a(4)+a(5)*a(5))**(1.5d0)
         kappa = sign(1.d0,mx(try(1)))*kappa
         indexCurv = indexfound
         ! This stops the code in case kappa becomes NaN.
!--debug         if(kappa.ne.kappa) call pariserror("HF9: Invalid Curvature")               
         ! if more than six independent heights found in all directions 
         return
      else 
         nfound = - ind_pos(points,nposit) 
      endif ! nfound == 9
! fixme: index not treated as it should (no reference to normal)
      height_found=.false.
      index=1
      do while(index<=6.and..not.height_found)
         if(height(i0,j0,k0,index) < 1d6) then   !  Use height if it exists
            height_found=.true.
            d = (index+1)/2
            origin=0.d0
            origin(d) = height(i0,j0,k0,index)
         endif
         index=index+1
      enddo
      if(.not.height_found) then
         call FindCutAreaCentroid(i0,j0,k0,centroid)
         do n=1,3
            origin(n) = centroid(n)
         enddo
      endif
      if ( (-nfound) > 6 )  then  ! more than 6 points to avoid special 2D degeneracy. 
         xfit=fit(:,try(2)) - origin(try(2))
         yfit=fit(:,try(3)) - origin(try(3))
         hfit=fit(:,try(1)) - origin(try(1))
         ! fit over all positions, not only independent ones. 
         call parabola_fit(xfit,yfit,hfit,nposit,a,fit_success) 
         if(.not.fit_success) call pariserror("no fit success")
         kappa = 2.d0*(a(1)*(1.d0+a(5)*a(5)) + a(2)*(1.d0+a(4)*a(4)) - a(3)*a(4)*a(5)) &
               /(1.d0+a(4)*a(4)+a(5)*a(5))**(1.5d0)
         indexCurv = indexfound
         kappa = sign(1.d0,mx(try(1)))*kappa
         ! This stops the code in case kappa becomes NaN.
!--debu         if(kappa.ne.kappa) call pariserror("HF6: Invalid Curvature")
         nposit=0
         return
      else  ! ind_pos <= 5
         ! Find all centroids in 3**3
         ! use direction closest to normal
         nfound = -100 + nfound  ! encode number of independent positions into nfound
         indexcurv=indexfound
         nposit=0
         debugfit=2e6
         do m=-1,1; do n=-1,1; do l=-1,1
            i=i0+m
            j=j0+n
            k=k0+l
            c(1)=m
            c(2)=n
            c(3)=l
            if(vof_flag(i,j,k) == 2) then
!               height_found=.false.
!
! comment the following as index not treated as it should (no reference to normal)
!
!               do index=1,6
!                   if(height(i,j,k,index) < 1d6) then   !  Use height if it exists
!                      height_found=.true.
!                      d = (index-1)/2 + 1
!                      nposit = nposit + 1
!                      do s=1,3  ! positions relative to center if i0,j0,k0 in cell units
!                         if(d==s) then
!                            fit(nposit,s) = c(s)
!                         else
!                            fit(nposit,s) = c(s) + height(i,j,k,index)
!                         endif
!                      enddo ! do s
!                   endif ! height exists
!               enddo ! do index
!               if(.not.height_found) then
                  nposit = nposit + 1
                  call FindCutAreaCentroid(i,j,k,centroid)
                  do s=1,3 
                     fit(nposit,s) = centroid(s) + c(s)
                     if(debug_curvature) then
                        if(k0==k00.and.l==0.and.s<3) then
                           debugfit(m,n,s) = centroid(s) + c(s)
                        endif
                     endif
                  end do
!                  write(92,'(1X,E10.2)',advance='no') fit(:,1)
!               endif 
            endif ! vof_flag
         enddo; enddo; enddo ! do m,n,l
         ! arrange coordinates so height direction is closest to normal
         ! try(:) array contains direction closest to normal first

         xfit=fit(:,try(2)) - origin(try(2))
         yfit=fit(:,try(3)) - origin(try(3))
         hfit=fit(:,try(1)) - origin(try(1))
         if(nposit.gt.NPOS) call pariserror("GLH: nposit")
         call parabola_fit(xfit,yfit,hfit,nposit,a,fit_success)
         if(debug_curvature) then
            if(rank==0) then
               if(is<=i0.and.i0<=ie.and.js<=j0.and.j0<=je) then
                  call PlotCutAreaCentroid(i0,j0,k0,centroid,x1,y1,xvec,yvec)
                  write(80,'(4(E15.8,1X))') x1,y1,xvec,yvec
                  ! rescale
                  do n=1,2
                     centroid(n) = deltax*debugfit(0,0,n)
                  enddo
                  centroid(1) = centroid(1) + x(i0) 
                  centroid(2) = centroid(2) + y(j0)
                  write(81,'(2(E15.8,1X))') centroid
               endif
               if(i0==5.and.j0==5.and.k0==k00) then
                  if(try(3).ne.3) then
                     write(*,*) "try: ",try
                     call pariserror("unexpected try(3)")
                  endif
                  if(try(1).ne.1) then
                     call pariserror("unexpected try(1)")
                  endif
                  xpoint(0)=debugfit(0,0,1)
                  ypoint(0)=debugfit(0,0,2)
                  xpoint(1)=debugfit(1,0,1)
                  ypoint(1)=debugfit(1,0,2)
                  xpoint(2)=debugfit(0,1,1)
                  ypoint(2)=debugfit(0,1,2)
                  !  h  = a4 t + a1 t**2 + a6   h=x, t=y
                  do n=0,2
                     diff(n) = abs(xpoint(n) - a(4)*ypoint(n) - a(1)*ypoint(n)**2 - a(6))
                     write(*,*) n, diff(n), a, nposit
                  enddo
                  WRITE(*,*) "debugfit ",debugfit(:,1,1)
               endif ! i0,j0=5,5
            endif ! rank=0
         endif ! debug
!         if(rank==0) write(*,*) "PF6: kappa,try,mx,nposit ",kappa,try,mx,nposit
!         if(rank==0) write(*,*) "PF6: a, fit_success ", a,fit_success
         kappa = 2.d0*(a(1)*(1.d0+a(5)*a(5)) + a(2)*(1.d0+a(4)*a(4)) - a(3)*a(4)*a(5)) &
             /sqrt(1.d0+a(4)*a(4)+a(5)*a(5))**3
         kappa = sign(1.d0,mx(try(1)))*kappa
         ! This stops the code in case kappa becomes NaN.
!--debu         if(kappa.ne.kappa) then
!           if(rank==0) write(*,*) "PF6: kappa,try,mx,nposit ",kappa,try,mx,nposit
!--debu            call pariserror("PF6: Invalid Curvature")  
!--debu         endif
      end if ! -nfound > 5
   end subroutine get_curvature

   subroutine output_curvature()
      implicit none      
      integer :: i,j,k,indexCurv ! ,l,m,n
      integer :: ib 
      real(8) :: kappa,a(6)
      real(8) :: angle 
      real(8) :: kappamin
      real(8) :: kappamax
      real(8) :: kappa_exact
      real(8) :: L2_err_K, err_K
      real(8) :: S2_err_K
      real(8) :: Lm_err_K
      integer :: sumCount,nfound,nindepend
      integer :: nposit
      real(8), parameter :: PI= 3.14159265359d0

      OPEN(UNIT=89,FILE=TRIM(out_path)//'/curvature-'//TRIM(int2text(rank,padding))//'.txt')
      OPEN(UNIT=90,FILE=TRIM(out_path)//'/reference-'//TRIM(int2text(rank,padding))//'.txt')
      OPEN(UNIT=91,FILE=TRIM(out_path)//'/bigerror-'//TRIM(int2text(rank,padding))//'.txt')
      OPEN(UNIT=92,FILE=TRIM(out_path)//'/debug-'//TRIM(int2text(rank,padding))//'.txt')
      if(rank==0) then
         OPEN(UNIT=79,FILE=TRIM(out_path)//'/grid.txt')
         OPEN(UNIT=80,FILE=TRIM(out_path)//'/segments.txt')
         OPEN(UNIT=81,FILE=TRIM(out_path)//'/points.txt')
         OPEN(UNIT=82,FILE=TRIM(out_path)//'/parabola.txt')
      endif
      ib = 1
      kappamin = 1d20
      kappamax = -1d20
      sumCount = 0
      S2_err_K=0.d0
      Lm_err_K=0.d0
      if ( test_curvature ) then 
         kappa_exact = 2.d0/rad(ib)
         do i=is,ie; do j=js,je; do k=ks,ke
            ! find curvature only for cut cells
            if (vof_flag(i,j,k) == 2 ) then 
               call get_curvature(i,j,k,kappa,indexCurv,nfound,nposit,a)
               kappa = kappa*dble(Nx) ! Nx = L/deltax
               kappamax = max(ABS(kappa),kappamax)
               kappamin = min(ABS(kappa),kappamin)
               angle = atan2(y(j)-yc(ib),x(i)-xc(ib))/PI*180.d0
               write(89,'(2(E15.8,1X))') angle,kappa
               write(92,'(2(E15.8,1X),I4)') angle,kappa,nfound
               write(90,*) angle,kappa_exact
               err_K = ABS(ABS(kappa)-kappa_exact)/kappa_exact
               if ( err_K > 0.1d0 ) &
                    write(91,'(4(I3,1X),2(E15.8,1X),I4)') i,j,k,indexCurv,kappa,kappa_exact,nfound
            end if ! cvof(i,j,k)
         end do; end do; end do
      else if ( test_curvature_2D) then 
         k = (Nz+4)/2
         kappa_exact = 1.d0/rad(ib)
          do i=is,ie; do j=js,je
            if (vof_flag(i,j,k) == 2) then 
               if(debug_ij55) then
                  debug_curvature = .true.
               endif
               call get_curvature(i,j,k,kappa,indexCurv,nfound,nposit,a)
               debug_curvature = .false.
               ! This stops the code in case kappa becomes NaN.
               if(kappa.ne.kappa) call pariserror("OC: Invalid Curvature")  
               if(nfound==-1.or.abs(kappa)<EPS_GEOM) then
                  write(6,*) "i,j,k,nfound,nindepend,kappa ",i,j,k,nfound,nindepend,kappa
                  call pariserror("OC: curvature not found")
               else
                  kappa = kappa*dble(Nx)  ! Nx = L/deltax
                  kappamax = max(ABS(kappa),kappamax)
                  kappamin = min(ABS(kappa),kappamin)
                  angle = atan2(y(j)-yc(ib),x(i)-xc(ib))/PI*180.d0
                  write(89,*) angle,ABS(kappa)
                  write(90,*) angle,kappa_exact
                  err_K    = ABS(ABS(kappa)-abs(kappa_exact))/kappa_exact
                  S2_err_K    = S2_err_K  + err_K**2
                  Lm_err_K    = MAX(Lm_err_K,   err_K) 
                  sumCount = sumCount + 1
               endif ! valid curvature
            end if ! cvof(i,j,k)
         end do; end do

         L2_err_K    = sqrt(S2_err_K/dble(sumCount))
         write(*,*) 'L2 Norm:'
         write(*,'(I5,I5,1X,(E15.8,1X))') Nx,rank,L2_err_K
         write(*,*) 'Linfty Norm:'
         write(*,'(I5,I5,1X,(E15.8,1X))') Nx,rank,Lm_err_K
      end if ! test_curvature
      write(*,*) 'max, min, and exact ABS(kappa)', kappamax, kappamin,kappa_exact
      write(*,*) 'max relative error', MAX(ABS(kappamax-kappa_exact), ABS(kappamin-kappa_exact))/kappa_exact
      if(rank==0) then
         CLOSE(79)
         CLOSE(80)
         CLOSE(81)
         CLOSE(82)
      endif
      CLOSE(89)
      CLOSE(90)
      CLOSE(91)
      CLOSE(92)
    end subroutine output_curvature

    subroutine parabola_fit(xfit,yfit,hfit,nposit,a,fit_success)
      implicit none

      real(8), intent(in)  :: xfit(nposit),yfit(nposit),hfit(nposit)
      real(8), intent(out) :: a(6)
      logical, intent(out) :: fit_success

      real(8) :: m(6,6), invm(6,6),error
      real(8) :: rhs(6)
      integer :: ifit, im,jm, nposit
      logical :: inv_success
      real(8) :: x1,x2,x3,x4,y1,y2,y3,y4

      fit_success=.false.
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
            
      ! The matrix is m_ij = sum_n alpha^n_i alpha^n_j
      ! and the "alpha_i" are the factors of the a_i coefficients:
      ! 
      !   x^2, y^2, xy, x^2, y^2, 1

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

      integer :: l,m,n
      real(8) :: dmx,dmy,dmz, mxyz(3),px,py,pz
      real(8) :: invx,invy,invz
      real(8) :: alpha, al3d
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
         do l=-1,1; do m=-1,1; do n=-1,1
            stencil3x3(l,m,n) = cvof(i+l,j+m,k+n)
         enddo;enddo;enddo
         call youngs(stencil3x3,mxyz)
         dmx = mxyz(1)
         dmy = mxyz(2)
         dmz = mxyz(3)      
      else
         dmx = n1(i,j,k)      
         dmy = n2(i,j,k)      
         dmz = n3(i,j,k)
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
!      if(px.ne.px) call pariserror("FCAC:invalid px")
!      if(py.ne.py) call pariserror("FCAC:invalid py")
!      if(pz.ne.pz) call pariserror("FCAC:invalid pz")

      ! rotate
      centroid(1) = px*invx
      centroid(2) = py*invy
      centroid(3) = pz*invz
      ! shift to cell-center coordinates
      centroid(1) = centroid(1) - invx*0.5d0
      centroid(2) = centroid(2) - invy*0.5d0
      centroid(3) = centroid(3) - invz*0.5d0
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

     if(mx<0.d0.or.my<0.d0.or.mz<0.d0) call pariserror("invalid mx my mz")
     if(abs(mx+my+mz-1d0)>EPS_GEOM) call pariserror("invalid mx+my+mz")

     if (mx < EPS_GEOM) then
        nx = my
        ny = mz
        call LineCenter (nx,ny, alpha, qx,qy)
        px = 0.5d0
        py = qx
        pz = qy
        return
     endif
     if (my < EPS_GEOM) then
        nx = mz
        ny = mx
        call LineCenter (nx,ny, alpha, qx,qy)
        px = qy
        py = 0.5d0
        pz = qx
        return
     endif
     if (mz < EPS_GEOM) then
        call LineCenter (mx,my, alpha, px,py)
        pz = 0.5
        return
     endif

     if (alpha < 0.d0 .or. alpha > 1.d0) call pariserror("PAC: invalid alpha")

     area = alpha*alpha
     px = area*alpha
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

     if (mx < EPS_GEOM) then
        px = 0.5
        py = alpha;
        return
     endif

     if (my < EPS_GEOM) then
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

 !    if(px.ne.px) call pariserror("LAC:invalid px")
 !    if(py.ne.py) call pariserror("LAC:invalid px")

   end subroutine

!   direction closest to normal first
   subroutine orientation (m,c)
     implicit none
     real(8), intent(in) :: m(3)
     integer, intent(out) :: c(3)
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
     if(.not.debug_ij55.and.test_curvature_2D.and.nx<=8.and.ny<=8.and.nz<=2) then
        call  plot_curvature()
     endif
  end subroutine test_VOF_HF
 
  subroutine plot_curvature()
    implicit none
    integer :: i,j,k,iem,jem,n,i0,j0,k0
    real(8) :: centroid(3),x1,y1,xvec,yvec,kappa,a(6),xpoint(0:2),ypoint(0:2),pc(12,12,2),diff(0:2)
    integer :: direction,indexcurv,nfound,nposit
    real(8) :: centroid_scaled(2), deltax
    k0 = (Nz+4)/2
    deltax=dx(nx/2)

    if(rank==0) then
      OPEN(UNIT=79,FILE=TRIM(out_path)//'/grid.txt')
      OPEN(UNIT=80,FILE=TRIM(out_path)//'/segments.txt')
      OPEN(UNIT=81,FILE=TRIM(out_path)//'/points.txt')
      OPEN(UNIT=82,FILE=TRIM(out_path)//'/parabola.txt')
      jem = je - 2
      iem = ie - 2
      do i=js,jem
         write(79,'(4(E15.8,1X))') xh(is),yh(i),xh(iem)-xh(is),0.d0
      enddo
      do i=is,iem
         write(79,'(4(E15.8,1X))') xh(i),yh(js),0.,yh(jem)-yh(js)
      enddo
      do i=is,ie; do j=js,je
         if(vof_flag(i,j,k).eq.2) then
            call PlotCutAreaCentroid(i,j,k,centroid,x1,y1,xvec,yvec)
            do n=1,2
               pc(i,j,n) = centroid(n)
            enddo
            write(80,'(4(E15.8,1X))') x1,y1,xvec,yvec
            do n=1,2 
               centroid_scaled(n) = deltax*centroid(n) 
            enddo
            centroid_scaled(1) = centroid_scaled(1) + x(i)
            centroid_scaled(2) = centroid_scaled(2) + y(j)
            write(81,'(2(E15.8,1X))') centroid_scaled(1),centroid_scaled(2) 
         endif
      enddo; enddo
      i0=5; j0=5
      call get_curvature(i0,j0,k0,kappa,indexCurv,nfound,nposit,a)
      i=i0
      j=j0
      xpoint(0)=0.d0
      ypoint(0)=0.d0
      i=i0+1; j=j0
      xpoint(1)=pc(i,j,1) + 1 - pc(i0,j0,1)
      ypoint(1)=pc(i,j,2) + 0 - pc(i0,j0,2)
      i=i0; j=j0+1
      xpoint(2)=pc(i,j,1) + 0 - pc(i0,j0,1)
      ypoint(2)=pc(i,j,2) + 1 - pc(i0,j0,2)
!  h  = a4 x + a1 x**2 + a6
      do n=0,2
          diff(n) = abs(xpoint(n) - a(4)*ypoint(n) - a(1)*ypoint(n)**2 - a(6))
         write(*,*) n, diff(n), a, nposit
!         diff(n) = abs(xpoint(n) - a(5)*ypoint(n) - a(2)*ypoint(n)**2 - a(6))
!         write(*,*) n, diff(n), "second orientation", nposit
      enddo
      CLOSE(79)
      CLOSE(80)
      CLOSE(81)
      CLOSE(82)
   endif
 end subroutine plot_curvature
 
subroutine PlotCutAreaCentroid(i,j,k,centroid,x1,y1,xvec,yvec)
      implicit none
      integer, intent(in)  :: i,j,k
      real(8), intent(out) :: centroid(3),x1,y1,xvec,yvec
      integer :: l,m,n
      real(8) :: dmx,dmy,dmz, mxyz(3),px,py,pz
      real(8) :: invx,invy
      real(8) :: alpha, al3d
      real(8) :: stencil3x3(-1:1,-1:1,-1:1)
      logical :: inv, swap
      real(8) :: deltax, tmp

      deltax=dx(nx/2)
      ! plot cut area centroid 
      !***
      !     (1) normal vector: dmx,dmy,dmz, and |dmx|+|dmy|+|dmz| = 1.
      !     (2) dmx,dmy,dmz>0 and record sign
      !     (3) get alpha;               
      !     (4) compute centroid with dmx,dmy,dmz and alpha;
      !     (5) transfer to local coordinate;
      !*(1)*

      if(recomputenormals) then
         do l=-1,1; do m=-1,1; do n=-1,1
            stencil3x3(l,m,n) = cvof(i+l,j+m,k+n)
         enddo;enddo;enddo
         call youngs(stencil3x3,mxyz)
         dmx = mxyz(1)
         dmy = mxyz(2)
         dmz = mxyz(3)      
      else
         dmx = n1(i,j,k)      
         dmy = n2(i,j,k)      
         dmz = n3(i,j,k)
      endif
      if(abs(dmz).gt.EPS_GEOM) call pariserror("PCAC: invalid dmz.")
      !*(2)*  
      invx = 1.d0
      invy = 1.d0
      if (dmx .lt. 0.0d0) then
         dmx = -dmx
         invx = -1.d0
      endif
      if (dmy .lt. 0.0d0) then
         dmy = -dmy
         invy = -1.d0
      endif
      alpha = al3d(dmx,dmy,dmz,cvof(i,j,k))
      !*(4)*  
      call PlaneAreaCenter(dmx,dmy,dmz,alpha,px,py,pz)
      !*(5)*
      ! rotate
      centroid(1) = px*invx
      centroid(2) = py*invy
      ! shift to cell-center coordinates
      centroid(1) = centroid(1) - invx*0.5d0
      centroid(2) = centroid(2) - invy*0.5d0
      ! rescale
!      do n=1,2; centroid(n) = deltax*centroid(n); enddo
      !*(6) 
      ! test alpha
      inv = .false.
      swap = .false.
      if(dmy > dmx ) then
         swap = .true.
         tmp = dmy
         dmy = dmx
         dmx = tmp
      endif

      if(alpha > 0.5d0) then
         inv = .true.
         alpha = 1.d0 - alpha
      endif
      x1 = alpha/dmx
      y1 = 0.d0
      if(alpha < dmy) then
         xvec = - x1
         yvec = alpha/dmy
      else
         xvec = - dmy/dmx
         yvec = 1.d0
      endif
      if(inv) then
         x1 = 1.d0 - x1; y1 = 1.d0 - y1
         xvec = -xvec; yvec = - yvec
      endif
      if(swap) then
         tmp = y1
         y1 = x1
         x1 = tmp
         tmp = yvec
         yvec = xvec
         xvec = tmp
      endif

      ! shift to cell center coordinates
      x1 = x1 - 0.5d0; y1 = y1 - 0.5d0
      ! reverse if needed and rescale
      x1 = x1*invx*deltax
      xvec = xvec*invx*deltax
      y1 = y1*invy*deltax
      yvec = yvec*invy*deltax
      ! shift
      x1 = x1 + x(i)
      centroid(1) = centroid(1)
      y1 = y1 + y(j)
      centroid(2) = centroid(2)
   end subroutine PlotCutAreaCentroid

end module module_surface_tension



