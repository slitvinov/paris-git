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
module module_st_testing
  use module_grid
  use module_BC
  use module_IO
  use module_tmpvar
  use module_2phase
  use module_VOF
  use module_surface_tension
  implicit none
  integer :: testint
contains
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
        method_count=0
        call output_heights()
        call output_curvature()
     end if
     if(test_curvature_2D.and.nx<=16.and.ny<=16.and.nz<=2) then
        call plot_curvature()
     end if
  end subroutine test_VOF_HF

   subroutine output_heights()
     implicit none
     integer i,j,k,d,index, direction
     real(8) h, th
     integer :: normalsign

     if(normal_up) then
        normalsign=1
     else
        normalsign=-1
     endif

     k = nz/2 + 2  ! +2 because of ghost layers
     j = ny/2 + 2  ! +2 because of ghost layers

     if(k<ks.or.k>ke) return
     if(j<js.or.j>je) return

     OPEN(UNIT=89,FILE=TRIM(out_path)//'/height-'//TRIM(int2text(rank,padding))//'.txt')
     OPEN(UNIT=90,FILE=TRIM(out_path)//'/reference-'//TRIM(int2text(rank,padding))//'.txt')

     if(cylinder_dir==2) then
        ! search in z direction
        direction = 3
        index = 2*(direction-1) + 1 + (-normalsign+1)/2
        do i=is,ie
           th = normalsign*wave2ls(x(i),y(j),z(k),cylinder_dir)/dx(nx/2+2)
           write(90,100) x(i),th
           if (height(i,j,k,index).lt.1d6) then
              h = height(i,j,k,index)
           else
              ! search for height
              d=0
              h = 2d6
              do while(h.gt.1d6.and.k+d<ke.and.k-d>ks)
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
     else if(cylinder_dir==3) then
        ! search in y direction
        direction=2
        index = 2*(direction-1) + 1 + (-normalsign+1)/2
        do i=is,ie
           th = normalsign*wave2ls(x(i),y(j),z(k),cylinder_dir)/dx(nx/2+2)
           write(90,100) x(i),th
           if (height(i,j,k,index).lt.1d6) then
              h = height(i,j,k,index)
           else
              ! search for height
              d=0
              h = 2d6
              do while(h.gt.1d6.and.j+d<je.and.j-d>js)
                 d = d + 1
                 if (height(i,j+d,k,index).lt.1d6) then
                    h = height(i,j+d,k,index) + d
                    height(i,j,k,index) = h 
                 else if (height(i,j-d,k,index).lt.1d6) then
                    h = height(i,j-d,k,index) - d
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
     endif

100  format(2(f24.16))
101  format(f24.16,A2)
     close(89)
     close(90)
   end subroutine output_heights
!==================================================================================================
   subroutine output_curvature()
      implicit none      
      integer :: i,j,k! ,l,m,n
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

      OPEN(UNIT=89,FILE=TRIM(out_path)//'/curvature-'//TRIM(int2text(rank,padding))//'.txt')
      OPEN(UNIT=90,FILE=TRIM(out_path)//'/reference-'//TRIM(int2text(rank,padding))//'.txt')
      OPEN(UNIT=91,FILE=TRIM(out_path)//'/bigerror-'//TRIM(int2text(rank,padding))//'.txt')
      OPEN(UNIT=92,FILE=TRIM(out_path)//'/debug-'//TRIM(int2text(rank,padding))//'.txt')
      ib = 1
      kappamin = 1d20
      kappamax = -1d20
      sumCount = 0
      S2_err_K=0.d0
      Lm_err_K=0.d0
      method_count=0
      if ( test_curvature ) then 
         kappa_exact = - 2.d0/rad(ib)
         do i=is,ie; do j=js,je; do k=ks,ke
            ! find curvature only for cut cells
            if (vof_flag(i,j,k) == 2 ) then 
               call get_curvature(i,j,k,kappa,nfound,nposit,a,.false.)
               if(nfound > 0) then
                  method_count(1) = method_count(1) + 1  ! nine heights
               else if( -nfound < 50) then 
                  method_count(2) = method_count(2) + 1  ! mixed heights
               else if (-nfound > 50) then
                  method_count(3) = method_count(3) + 1  ! centroids
               else
                  call pariserror("OC: unknown method_count") 
               endif
               kappa = kappa*dble(Nx) ! Nx = L/deltax
               kappamax = max(ABS(kappa),kappamax)
               kappamin = min(ABS(kappa),kappamin)
               angle = atan2(y(j)-yc(ib),x(i)-xc(ib))/PI*180.d0
               write(89,'(2(E15.8,1X))') angle,kappa
               write(92,'(2(E15.8,1X),I4)') angle,kappa,nfound
               write(90,*) angle,kappa_exact
               err_K = ABS(ABS(kappa)-kappa_exact)/kappa_exact
               if ( err_K > 0.1d0 ) &
                    write(91,'(3(I3,1X),2(E15.8,1X),I4)') i,j,k,kappa,kappa_exact,nfound
            end if ! cvof(i,j,k)
         end do; end do; end do
      else if ( test_curvature_2D) then 
         k = (Nz+4)/2
         kappa_exact = - 1.d0/rad(ib)
          do i=is,ie; do j=js,je
            if (vof_flag(i,j,k) == 2) then 
               call get_curvature(i,j,k,kappa,nfound,nposit,a,.false.)
               if(nfound > 0) then
                  method_count(1) = method_count(1) + 1  ! nine heights
               else if( -nfound < 50) then 
                  method_count(2) = method_count(2) + 1  ! mixed heights
               else if (-nfound > 50) then
                  method_count(3) = method_count(3) + 1  ! centroids
               endif
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
                  write(89,*) angle,kappa
                  write(90,*) angle,kappa_exact
                  err_K    = ABS(kappa-kappa_exact)/kappa_exact
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
      CLOSE(89)
      CLOSE(90)
      CLOSE(91)
      CLOSE(92)
      call print_method() 
contains
  subroutine print_method()
    integer :: total=0
    integer :: n
    real(8) :: fraction(3)
    do n=1,3
       total = method_count(n) + total
    enddo
    do n=1,3
       fraction(n) = float(method_count(n)) / float(total)   
    enddo

    OPEN(UNIT=89,FILE='mcount.tmp')
    write(89,*) fraction
    close(89)
end subroutine print_method
    end subroutine output_curvature
!=========================================================================================================
!
!  Testing section
! 
!=========================================================================================================
  subroutine plot_curvature()
    implicit none
    integer :: i,j,k,iem,jem,n
    real(8) :: centroid(3),x1,y1
    real(8), allocatable :: pc(:,:,:)
    real(8) :: centroid_scaled(2), deltax
    k = (Nz+4)/2
    deltax=dx(nx/2)

    allocate(pc(imin:imax,jmin:jmax,3))
    if(rank==0.and.cylinder_dir==3) then
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
            call PlotCutAreaCentroid(i,j,k,centroid,x1,y1)
            do n=1,2
               pc(i,j,n) = centroid(n)
            enddo
            write(80,'(2(E15.8,1X))') x1,y1
            do n=1,2 
               centroid_scaled(n) = deltax*centroid(n) 
            enddo
            centroid_scaled(1) = centroid_scaled(1) + x(i)
            centroid_scaled(2) = centroid_scaled(2) + y(j)
            write(81,'(2(E15.8,1X))') centroid_scaled(1),centroid_scaled(2) 
         endif
      enddo; enddo
      CLOSE(79)
      CLOSE(80)
      CLOSE(81)
      CLOSE(82)
   endif
 end subroutine plot_curvature
 
subroutine PlotCutAreaCentroid(i,j,k,centroid,x1,y1)
      implicit none
      integer, intent(in)  :: i,j,k
      real(8), intent(out) :: centroid(3),x1,y1
      integer :: l,m,n
      real(8) :: nr(3),dmx,dmy, al3dnew
      real(8) :: stencil3x3(-1:1,-1:1,-1:1)

      if(recomputenormals) then
         do l=-1,1; do m=-1,1; do n=-1,1
            stencil3x3(l,m,n) = cvof(i+l,j+m,k+n)
         enddo;enddo;enddo
         call youngs(stencil3x3,nr)
      else
         nr(1) = n1(i,j,k)      
         nr(2) = n2(i,j,k)      
         nr(3) = n3(i,j,k)
      endif
      dmx = nr(1)
      dmy = nr(2)
      if(abs(nr(3)).gt.EPS_GEOM) call pariserror("PCAC: invalid dmz.")
      call cent3D(nr,cvof(i,j,k),centroid)
      centroid = centroid - 0.5d0
      x1 = - al3dnew(nr,cvof(i,j,k))/dmx
      y1 = - al3dnew(nr,cvof(i,j,k))/dmy
      ! shift to cell center coordinates
      x1 = x1 - 0.5d0; y1 = y1 - 0.5d0
      ! shift
      x1 = x1 + x(i)
      y1 = y1 + y(j)
! some stuff is missing here
   end subroutine PlotCutAreaCentroid

  end module module_st_testing



