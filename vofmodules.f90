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
! module_VOF: Contains definition of variables for the Volume of Fluid interface tracking.
!-------------------------------------------------------------------------------------------------
module module_VOF
  use module_grid
  use module_IO
  use module_tmpvar
  implicit none
  real(8), dimension(:,:,:), allocatable, target :: cvof ! VOF tracer variable
  integer, dimension(:,:,:), allocatable :: vof_flag ! 
  !   0 empty
  !   1 full
  !   2 fractional
  !   3 unknown

  real(8), parameter  :: A_h = 2d0  ! For initialisation of height test
  real(8), parameter  :: TINY = 1d-50
  real(8), parameter  :: EPSC = 1.d-12  ! for clipping vof and setting flags
  real(8), parameter  :: EPSDP = 1d-16  ! assumed precision of double precision computations

  character(20) :: vofbdry_cond(6),test_type,vof_advect,cond
  integer :: parameters_read=0, refinement=-1 
  integer :: cylinder_dir = 0
  logical :: test_heights = .false.  
  logical :: normal_up = .true.    ! used for the h
  logical :: test_curvature = .false.  
  logical :: test_curvature_2D = .false.  
  logical :: test_HF = .false.
  logical :: test_LP = .false.
  logical :: test_tag = .false.
  logical :: test_D2P = .false.
  logical :: test_bubbles = .false.
  logical :: linfunc_initialized = .false.
  real(8) :: b1,b2,b3,b4
  integer :: nfilter=0

  logical :: DoLPP = .false.
contains
!=================================================================================================
!=================================================================================================
!__|__|__|__|
!__|__|_2.5_|    Maximum possible filter value
!__|__|_1.5_|
!__|__|_x|__|
!__|1.9__|__|
!2.9__|__|__|
!
  subroutine initialize_linfunc()
    implicit none
! Gerris method
! Each vertex averages 8 cells  --> 1/8
! Each Cell averages 8 vertices --> 1/8
! 8 vertex cells belong to only one vertex ->1/64
! 12 edge cells belong to 2 vertices -> 2/64
! 6 face cells belong to 4 vertices -> 4/64
! 8/4 + 6*4/64 + 12*2/64 + 8/64 = 1
    b1=8D0/64
    b2=4D0/64
    b3=2D0/64
    b4=1D0/64
    linfunc_initialized = .true.
  end subroutine initialize_linfunc
!------------------------------------------------------------------------
  subroutine filter(field)
    use module_tmpvar
    implicit none
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(out) :: field
    integer :: i,j,k,n
    field=cvof
    do n=1,nfilter
       do k=ks-1,ke+1; do j=js-1,je+1; do i=is-1,ie+1
          tmp(i,j,k) = b1*field(i,j,k) + & 
               b2*( field(i-1,j,k) + field(i,j-1,k) + field(i,j,k-1) + &
                    field(i+1,j,k) + field(i,j+1,k) + field(i,j,k+1) ) + &
               b3*( field(i+1,j+1,k) + field(i+1,j-1,k) + field(i-1,j+1,k) + field(i-1,j-1,k) + &
                    field(i+1,j,k+1) + field(i+1,j,k-1) + field(i-1,j,k+1) + field(i-1,j,k-1) + &
                    field(i,j+1,k+1) + field(i,j+1,k-1) + field(i,j-1,k+1) + field(i,j-1,k-1) ) + &
               b4*( field(i+1,j+1,k+1) + field(i+1,j+1,k-1) + field(i+1,j-1,k+1) + field(i+1,j-1,k-1) +  &
                    field(i-1,j+1,k+1) + field(i-1,j+1,k-1) + field(i-1,j-1,k+1) + field(i-1,j-1,k-1) )
! fixme: ghostx missing
       enddo; enddo; enddo
       field=tmp
    enddo
  end subroutine filter
!------------------------------------------------------------------------
  subroutine linfunc(field,a1,a2)
    implicit none
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(out) :: field
    real(8) :: cfiltered
    real(8), intent(in) :: a1,a2
    integer :: i,j,k
    if(.not.linfunc_initialized) call initialize_linfunc
    if(nfilter==0) then
       field = cvof*(a2-a1)+a1
    else if (nfilter==1) then
       do k=ks-1,ke+1; do j=js-1,je+1; do i=is-1,ie+1
          cfiltered = b1*cvof(i,j,k) + & 
               b2*( cvof(i-1,j,k) + cvof(i,j-1,k) + cvof(i,j,k-1) + &
                    cvof(i+1,j,k) + cvof(i,j+1,k) + cvof(i,j,k+1) ) + &
               b3*( cvof(i+1,j+1,k) + cvof(i+1,j-1,k) + cvof(i-1,j+1,k) + cvof(i-1,j-1,k) + &
                    cvof(i+1,j,k+1) + cvof(i+1,j,k-1) + cvof(i-1,j,k+1) + cvof(i-1,j,k-1) + &
                    cvof(i,j+1,k+1) + cvof(i,j+1,k-1) + cvof(i,j-1,k+1) + cvof(i,j-1,k-1) ) + &
               b4*( cvof(i+1,j+1,k+1) + cvof(i+1,j+1,k-1) + cvof(i+1,j-1,k+1) + cvof(i+1,j-1,k-1) +  &
                    cvof(i-1,j+1,k+1) + cvof(i-1,j+1,k-1) + cvof(i-1,j-1,k+1) + cvof(i-1,j-1,k-1) )
          field(i,j,k) = cfiltered*(a2-a1)+a1
       enddo; enddo; enddo
    endif
  end subroutine linfunc
!=================================================================================================

  subroutine ReadVOFParameters
    use module_flow
    use module_BC
    implicit none
    include 'mpif.h'
    integer ierr,in
    logical file_is_there
    namelist /vofparameters/ vofbdry_cond,test_type,VOF_advect,refinement, &
       cylinder_dir, normal_up, DoLPP
    in=31

    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    inquire(file='inputvof',exist=file_is_there)
    open(unit=in, file='inputvof', status='old', action='read', iostat=ierr)
    if (file_is_there) then
       if(ierr == 0) then
          read(UNIT=in,NML=vofparameters)
          if(rank==0) write(out,*)'VOF parameters read successfully'
       else
          print *, 'rank=',rank,' has error ',ierr,' opening file inputsolids'
       endif
    else
       if (rank == 0) STOP "ReadVOFParameters: no 'inputvof' file."
    endif
    close(in)
    if(refinement==-1) then
       refinement=8
       if(rank==0) write(*,*) "using default value for refinement"
    endif
  end subroutine ReadVOFParameters
!
  subroutine initialize_VOF
    use module_grid
    use module_BC
    implicit none
    include 'mpif.h'
    integer :: ierr,dir,orientation
    character :: dc(3) = (/ "x","y","z" /)
    call ReadVOFParameters
! Check grid
    if(read_x.or.read_y.or.read_z) then
       if((xform.ne.0.d0).or.(yform.ne.0.d0).or.(zform.ne.0.d0)) then
          if (rank == 0) print *, "VOF does not yet work with variable grids"
          call MPI_Finalize(ierr)
          stop
       endif
    endif
    allocate(cvof(imin:imax,jmin:jmax,kmin:kmax),vof_flag(imin:imax,jmin:jmax,kmin:kmax))
    cvof = 0.d0
    vof_flag = 3
    if(test_type=='droplet') then
       test_heights = .false.
    else if(test_type=='height_test') then
       test_heights = .true.
    else if(test_type=='Curvature_test') then
       test_curvature = .true.
    else if(test_type=='Curvature2D') then
       test_curvature_2D = .true.
    else if(test_type=='tag_test') then
       test_tag = .true.
    else if(test_type=='D2P_test') then
       test_D2P = .true.
    else if(test_type=='bubbles_test') then
       test_bubbles = .true.
    else
       write(*,*) test_type, rank
       stop 'unknown initialization'
    endif
    test_HF = test_heights .or. test_curvature .or. test_curvature_2D
    test_LP = test_tag .or. test_D2P 

! Post-read initialization of boundary conditions
    ! orientation order:    
    ! x- y- z- x+ y+ z+

    if(.not.bdry_read) stop "bdry not read"
    do orientation=1,6
       cond = vofbdry_cond(orientation)
       dir = orientation
       if(orientation>3) dir = orientation-3
       if(iachar(cond(1:1))==0) then
          if(orientation<=3) stop "missing vof bdry condition"
          ! mirror the color or periodicity except if outflow is set
          vofbdry_cond(orientation) = vofbdry_cond(dir)
          if(bdry_cond(orientation)==4) vofbdry_cond(orientation) = 'outflow'   ! outflow +
       endif
    enddo

    if(rank==0) then
       do dir=1,3
          write(6,'(A1,"-: ",(A)," / ",(A)T32," ",A1,"+: ",(A)," / ",(A))') &
               dc(dir),TRIM(expl(bdry_cond(dir))),TRIM(vofbdry_cond(dir)), &
               dc(dir),TRIM(expl(bdry_cond(dir+3))),TRIM(vofbdry_cond(dir+3)) ! ,rank
       enddo
    endif
  end subroutine initialize_VOF
!=================================================================================================
!   a hack to get the flags quickly (temporary)
!=================================================================================================
  subroutine get_flags_and_clip()
    integer :: i,j,k
    if(ng.lt.2) stop "wrong ng"
    do k=kmin,kmax
       do j=jmin,jmax
          do i=imin,imax
             if(cvof(i,j,k).le.EPSC) then
                vof_flag(i,j,k) = 0
                cvof(i,j,k)=0.d0
             else if(cvof(i,j,k).ge.1.d0-EPSC) then
                vof_flag(i,j,k) = 1
                cvof(i,j,k) = 1.d0
             else
                vof_flag(i,j,k) = 2
             endif
          enddo
       enddo
    enddo
  end subroutine get_flags_and_clip
  !=================================================================================================
  !  Initialize vof field and flags
  !=================================================================================================
  subroutine initconditions_VOF()
    use module_hello
    use module_flow
    use module_BC
    use module_2phase

    implicit none
    include 'mpif.h'
    integer :: ierr, req(12),sta(MPI_STATUS_SIZE,12)
    integer , parameter :: ngh=2
    integer :: ipar
    integer, parameter :: root_rank = 0
    
    if( test_D2P) then 
       if ( rank == root_rank ) call random_bubbles
       call MPI_BCAST(rad, NumBubble, MPI_REAL8, &
                      root_rank, MPI_Comm_Cart, ierr)
       call MPI_BCAST(xc , NumBubble, MPI_REAL8, &
                      root_rank, MPI_Comm_Cart, ierr)
       call MPI_BCAST(yc , NumBubble, MPI_REAL8, &
                      root_rank, MPI_Comm_Cart, ierr)
       call MPI_BCAST(zc , NumBubble, MPI_REAL8, &
                      root_rank, MPI_Comm_Cart, ierr)
    end if ! test_D2P

    if(test_heights) then 
       if(cylinder_dir==0) then
          write(*,*) "IVOF: Warning: cylinder_dir=0 set to 2"
          cylinder_dir=2
       endif
       ipar=cylinder_dir
       call levelset2vof(wave2ls,ipar)
    else if(NumBubble>0) then
       ! one cylinder in -ipar>0 direction otherwise spheres
       ipar=-cylinder_dir
      call levelset2vof(shapes2ls,ipar)
    else 
       if(bdry_cond(1) /= 3 ) write(*,*) &
            "IVOF: Warning: Nothing set. cylinder_dir=0 set to 2"
       cvof=0.d0
       vof_flag=0
    endif
    call ghost_x(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
    call ghost_y(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
    call ghost_z(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
    call ighost_x(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
    call ighost_y(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
    call ighost_z(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
    call setVOFBC(cvof,vof_flag)
    return
  end subroutine initconditions_VOF
  !=================================================================================================
  !   Generate random bubbles 
  !=================================================================================================
  subroutine random_bubbles()
    use module_2phase
    implicit none
    integer ib

    if(NumBubble>2) then 
      do ib=1,NumBubble
         rad(ib) = 0.02 + rand()*0.04
         xc(ib)  = 0.15 + rand()*0.7
         yc(ib)  = 0.15 + rand()*0.7
         zc(ib)  = 0.15 + rand()*0.7
      end do
       
    end if ! NumBubble
  end subroutine random_bubbles 
  !=================================================================================================
  !   Spheres and cylinders
  !=================================================================================================
  function shapes2ls(xx,yy,zz,ipar)
    use module_2phase
    implicit none
    real(8), intent(in) :: xx,zz,yy
    integer, intent(in) :: ipar
    real(8) :: a, cdir(0:3), shapes2ls
    integer ib

    if(.not.(-3<=ipar.and.ipar<=0)) call pariserror("S: invalid ipar")
    cdir(1:3) = 1.d0/(1.d0 + excentricity(1:3))
    cdir(-ipar) = 0.d0
    shapes2ls = -2.d6
    ! ipar < 0 cylinder case
    ! ipar = 0 spheres
    if(ipar < 0.and.NumBubble/=1) call pariserror("S: invalid NumBubbles")
    do ib=1,NumBubble
       a = rad(ib)**2 - (cdir(1)*(xx-xc(ib))**2+cdir(2)*(yy-yc(ib))**2+cdir(3)*(zz-zc(ib))**2) !fixme
       shapes2ls = MAX(shapes2ls,a)
    end do
  end function shapes2ls
  !=================================================================================================
  !  sine-wave interface
  !=================================================================================================
  function wave2ls(xx,yy,zz,ipar)
    use module_2phase
    implicit none
    real(8) wave2ls
    real(8), intent(in) :: xx,zz,yy
    integer, intent(in) :: ipar
    integer :: hpar,vpar,normalsign
    real(8) ::  vdir(13),hdir(3)

    if(normal_up) then
       normalsign=1
    else
       normalsign=-1
    endif

    if(.not.(1<=ipar.and.ipar<=3)) call pariserror("invalid ipar")
    if (ipar == 1) then
       hpar = 2
       vpar = 3
    else if (ipar == 2) then
       hpar = 1
       vpar = 3
    else
       hpar = 1
       vpar = 2
    endif
       
    vdir = 0.d0
    vdir(vpar) = 1.d0
    hdir = 0.d0
    hdir(hpar) = 1.d0

    wave2ls = vdir(3)*(- zz + zLength/2.d0) + vdir(2)*(-yy + yLength/2.d0) &
         + vdir(1)*(-xx + xLength/2.d0) &
         + A_h*dx(nx/2+2)*cos(2.*3.14159*(hdir(1)*xx/xLength &
         + hdir(2)*yy/yLength + hdir(3)*zz/zLength))
    wave2ls = wave2ls*dble(normalsign)
  end function wave2ls
  !=================================================================================================
  !   Converts a level-set field into a VOF field
  !=================================================================================================
  subroutine levelset2vof(lsfunction,ipar)
    use module_BC
    implicit none
    real(8), external :: lsfunction
    integer, intent(in) :: ipar
    include 'mpif.h'
    integer :: ierr, req(12),sta(MPI_STATUS_SIZE,12)
    integer , parameter :: ngh=2

    call ls2vof_refined(lsfunction,ipar,1)
    call ighost_x(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
    call ighost_y(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
    call ighost_z(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)

    call ghost_x(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
    call ghost_y(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
    call ghost_z(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
    call setVOFBC(cvof,vof_flag)
    call ls2vof_refined(lsfunction,ipar,refinement)
  end subroutine levelset2vof
  !=================================================================================================
  !   Finds cells surrounded by all empty or all full cells
  !=================================================================================================
  function notisolated(istencil3x3,is2D)
    implicit none
    logical :: notisolated
    integer, intent(in) :: is2D
    integer :: istencil3x3(-1:1,-1:1,-1:1)
    integer :: nfrac
    integer :: nfull
    integer :: i0,j0,k0
    integer :: isnot2D, ncells

    nfrac=0; nfull=0; ncells=0
    isnot2D = 1 - is2D
!     write(*,*) "isnot2D = ", isnot2D
!     if(is2D==1) write(*,'(9I1)',advance='no') istencil3x3(:,:,0)
    do i0=-1,1; 
       do j0=-1,1; 
          do k0=-isnot2D,isnot2D
             nfrac = nfrac + istencil3x3(i0,j0,k0)/2
             nfull = nfull + mod(istencil3x3(i0,j0,k0),2)
             ncells = ncells + 1
          enddo; enddo; enddo
    notisolated = .not.(nfrac==0.and.(nfull==0.or.nfull==ncells))
    ! write (*,*) " nfrac,nfull,ncells ",  nfrac,nfull,ncells
  end function notisolated
!=================================================================================================
!
! mapping
!
!=================================================================================================
  subroutine map3x3in2x2(i1,j1,k1,i0,j0,k0)
    implicit none
    integer, intent(in) :: i0,j0,k0
    integer, intent(out) :: i1(-1:1,-1:1,3), j1(-1:1,-1:1,3), k1(-1:1,-1:1,3)
    integer m,n
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
  end subroutine map3x3in2x2

  subroutine ls2vof_refined(lsfunction,ipar,n1)
    implicit none
    real(8), external :: lsfunction
    integer, intent(in) :: ipar,n1
    real(8) :: stencil3x3(-1:1,-1:1,-1:1),dx1,dy1,dz1,x0,y0,z0,x1,y1,z1,a,b
    integer :: i,j,k,i0,j0,k0,l,m,n,s
    integer :: nfrac,nflag,nfull
    integer :: istencil3x3(-1:1,-1:1,-1:1)
    integer :: i1(-1:1,-1:1,3), j1(-1:1,-1:1,3), k1(-1:1,-1:1,3)
    logical :: refinethis 
    real(8) :: count
    integer :: calc_imax
    integer :: dirselect(0:3), d, is2D,max_flag

! initialization
    count=0.d0
    refinethis = .false.

! Initialize 2D/3D switch
    dirselect = 1  ! spheres: all directions selected: default. 
    d = max(ipar,-ipar)
    dirselect(d)=0

!    Some error checking
    if(d>3) call pariserror("wrong ipar")
    if(min(min(nx,ny),nz)<2) call pariserror("minimum dimension nx ny nz too small")
    max_flag=calc_imax(vof_flag)
    if(n1>1.and.max_flag/=2.and.A_h>1d-16) then
       if(rank==0) then
          if(max_flag==0) then
             write(*,*) "ls2vof_refined: error: single phase flow ? Nothing to initialize !?"
          else
             write(*,*) "ls2vof_refined: maximum vof_flag = ", max_flag, "but expecting maximum flag = 2"
          endif
       endif
       call pariserror("bad flag")
    endif
 
! main loop
    do k=ks,ke; do j=js,je; do i=is,ie
       is2D=0  ! Default
       if(n1>1) then  ! refinement on second pass only
          if(d==0) then ! check for isolated cells
             do i0=-1,1; do j0=-1,1; do k0=-1,1
                istencil3x3(i0,j0,k0) = vof_flag(i+i0,j+j0,k+k0)
             enddo; enddo; enddo
          else if(d>0) then 
             call map3x3in2x2(i1,j1,k1,i,j,k)
             do m=-1,1; do n=-1,1
                istencil3x3(m,n,0) = vof_flag(i1(m,n,d),j1(m,n,d),k1(m,n,d))
             enddo; enddo
             is2D=1
          else
             call pariserror("bad d")
          endif
          refinethis = notisolated(istencil3x3,is2D)
          if(refinethis) count = count + 1.
       endif
! refine and initialize subcells
       if(n1==1.or.refinethis) then  ! if n1>1 and no refine leave cell as is. 
          dx1 = dx(i)/n1; dy1 = dy(j)/n1; dz1 = dz(k)/n1
          nfrac=0; nfull=0
          b=0.d0
          s=n1/2
          do l=0,n1-1; do m=0,n1-1; do n=0,n1-1
             x0 = x(i) - 0.5d0*dx(i) + 0.5d0*dx1 + dx1*l
             y0 = y(j) - 0.5d0*dy(j) + 0.5d0*dy1 + dy1*m
             z0 = z(k) - 0.5d0*dz(k) + 0.5d0*dz1 + dz1*n
             do i0=-1,1; do j0=-1,1; do k0=-1,1
                x1 = x0 + i0*dx1; y1 = y0 + j0*dy1; z1 = z0 + k0*dz1 
                stencil3x3(i0,j0,k0) = lsfunction(x1,y1,z1,ipar)
             enddo; enddo; enddo
             call ls2vof_in_cell(stencil3x3,a,nflag)
             if(nflag==2) then 
                nfrac = nfrac + 1
             else if(nflag==1) then
                nfull = nfull + 1
             endif
             b=b+a   ! *(1)*
          enddo; enddo; enddo
          cvof(i,j,k) = b/(n1**3)
          if(nfrac > 0) then
             vof_flag(i,j,k) = 2
             ! now either all full, all empty, or mix full/empty : 
          else if(nfull==n1**3) then ! all full
             vof_flag(i,j,k) = 1
             cvof(i,j,k) = 1.d0  ! because arithmetic at (1) may cause round off errors
          else if(nfull==0) then ! all empty
             vof_flag(i,j,k) = 0
             cvof(i,j,k) = 0.d0  ! paranoid programming.
          else ! mix of full and empty
             vof_flag(i,j,k) = 2
          end if
       endif
    enddo; enddo; enddo
!    IF(N1>1) write(*,*) "proportion refined ",100.*count/(nx*ny*nz),"%"
    return
  end subroutine ls2vof_refined
  !=================================================================================================
  subroutine c_mask(cbinary)
    implicit none
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(out) :: cbinary

    where (cvof > 0.5d0)
       cbinary = 1.d0
    elsewhere
       cbinary = 0.d0
    end where
  end subroutine c_mask
  !=================================================================================================
  subroutine vofsweeps(tswap)
    use module_BC
    use module_flow
    use module_tmpvar
    implicit none
    include 'mpif.h'
    integer, intent(in) :: tswap
    integer :: req(48),sta(MPI_STATUS_SIZE,48)
    integer, parameter :: ngh=2
    integer :: ierr

    if (VOF_advect=='Dick_Yue') call c_mask(work(:,:,:,2))
    if (MOD(tswap,3) .eq. 0) then
       call swp(w,cvof,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3),vof_flag,3)
       call ghost_x(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ghost_y(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ghost_z(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ighost_x(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ighost_y(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ighost_z(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)

       call swp(u,cvof,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3),vof_flag,1)
       call ghost_x(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ghost_y(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ghost_z(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ighost_x(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ighost_y(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ighost_z(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)

       call swp(v,cvof,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3),vof_flag,2)
       call ghost_x(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ghost_y(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ghost_z(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ighost_x(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ighost_y(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ighost_z(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)

    elseif (MOD(tswap,2) .eq. 0) then
       call swp(v,cvof,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3),vof_flag,2)
       call ghost_x(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ghost_y(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ghost_z(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ighost_x(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ighost_y(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ighost_z(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)

       call swp(w,cvof,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3),vof_flag,3)
       call ghost_x(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ghost_y(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ghost_z(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ighost_x(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ighost_y(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ighost_z(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)

       call swp(u,cvof,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3),vof_flag,1)
       call ghost_x(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ghost_y(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ghost_z(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ighost_x(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ighost_y(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ighost_z(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
    else 
       call swp(u,cvof,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3),vof_flag,1)
       call ghost_x(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ghost_y(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ghost_z(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ighost_x(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ighost_y(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ighost_z(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)

       call swp(v,cvof,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3),vof_flag,2)
       call ghost_x(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ghost_y(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ghost_z(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ighost_x(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ighost_y(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ighost_z(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)

       call swp(w,cvof,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3),vof_flag,3)
       call ghost_x(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ghost_y(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ghost_z(cvof,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ighost_x(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ighost_y(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ighost_z(vof_flag,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
    endif
  end subroutine vofsweeps
!=================================================================================================
! subroutine SetVOFBC: Sets the VOF fraction boundary condition
!-------------------------------------------------------------------------------------------------
  subroutine SetVOFBC(cv,fl)
    use module_grid
    use module_BC
    implicit none
    include 'mpif.h'
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: cv  ! cvof
    integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: fl  ! vof_flag
    integer :: fb(6),d,l,m,n,c(3),try(2:3),sign,orientation,flag,flhere
    real(8) :: xi,eta,cvhere

    do orientation=1,6
       cond = vofbdry_cond(orientation)
       if(cond=='wet') then
          fb(orientation)=1
       else if(cond=='dry') then
          fb(orientation)=0
       else if(cond=='periodic') then
          fb(orientation) = 3
       else if(cond=='outflow') then
          fb(orientation) = 4 
        else if(cond=='jet') then
           if(orientation /= 1) stop "jet only at x-"
           fb(orientation) = 2
       else
          call pariserror("this vofbc not implemented")
       endif
    enddo
! 
    do orientation=1,6
    ! orientation order:    
    ! x- y- z- x+ y+ z+
       d = orientation
       sign=-1
       if(orientation>3) then
          d = orientation-3
          sign=1
       endif
       flag=fb(orientation)
! sort directions so that try(1) = d and try(2),try(3) are any other two directions. 
       m=1
       n=2  
       do while (m.le.3)
          if(m.ne.d) then
             try(n) = m
             n=n+1
          endif
          m=m+1
       enddo
! end sort
       if(coords(d)==proclimit(d,sign)) then
          do l=coordstart(try(2))-Ng,coordend(try(2))+Ng
             do m=coordstart(try(3))-Ng,coordend(try(3))+Ng
                c(try(2)) = l; c(try(3)) = m
                c(d)=coordlimit(d,sign) + sign
                if(flag<2) then
                      cv(c(1),c(2),c(3))=dble(flag)
                      fl(c(1),c(2),c(3))=flag
                      c(d) = c(d) + sign
                      cv(c(1),c(2),c(3))=dble(flag)
                      fl(c(1),c(2),c(3))=flag
                   elseif(flag==2) then  ! jet
                      cv(c(1),c(2),c(3))=dble(inject(c(2),c(3)))
                      fl(c(1),c(2),c(3))=inject(c(2),c(3))
                   elseif(flag==4) then
                      c(d)=coordlimit(d,sign)
                      cvhere=cv(c(1),c(2),c(3))
                      flhere=fl(c(1),c(2),c(3))
                      c(d) = c(d) + sign
                      cv(c(1),c(2),c(3))=cvhere
                      fl(c(1),c(2),c(3))=flhere
                      c(d) = c(d) + sign
                      cv(c(1),c(2),c(3))=cvhere
                      fl(c(1),c(2),c(3))=flhere
                   endif
                enddo
             enddo
          endif
       enddo
  contains
  function inject(j,k)
    use module_grid
    implicit none
    integer :: j,k
    integer :: inject
    inject=0d0
    if((y(j) - 0.5d0)**2 + (z(k) - 0.5d0)**2.lt.0.25d0*0.25d0) inject=1
  end function inject
!
   function xcoord(d,i)
      implicit none
      real(8) :: xcoord
      integer, intent(in) :: d,i
      if(d==1) then
         xcoord = x(i)
      else if(d==2) then
         xcoord = y(i)
      else
         xcoord = z(i)
      endif
    end function xcoord
  end subroutine SetVOFBC
!=================================================================================================
  subroutine test_cell_size()
    implicit none
    if(dabs(dyh(js)-dxh(is))*1d14/dxh(is)>1d0.or.dabs(dzh(ks)-dxh(is))*1d14/dxh(is)>1d0) then
       print *, "non-cubic cells"
       stop
    endif
  end subroutine test_cell_size
end module module_vof
!=================================================================================================
!-------------------------------------------------------------------------------------------------
module module_output_vof
  use module_IO
  use module_flow
  use module_grid
  use module_solid
  use module_vof
  implicit none
  integer :: vof_opened=0;
contains
  subroutine append_VOF_visit_file(rootname)
    implicit none
    character(*) :: rootname
    integer prank
    if(rank.ne.0) call pariserror('rank.ne.0 in append_VOF')
    if(vof_opened==0) then
       OPEN(UNIT=88,FILE='vof.visit')
       write(88,10) nPdomain
10     format('!NBLOCKS ',I4)
       vof_opened=1
    else
       OPEN(UNIT=88,FILE='vof.visit',access='append')
    endif
    do prank=0,NpDomain-1
       write(88,11) rootname//TRIM(int2text(prank,padding))//'.vtk'
11     format(A)
    enddo
    close(88)
  end subroutine  append_VOF_visit_file
!=================================================================================================
  subroutine output_VOF(nf,i1,i2,j1,j2,k1,k2)
    implicit none
    integer ::nf,i1,i2,j1,j2,k1,k2,i,j,k
    character(len=30) :: rootname
    rootname=trim(out_path)//'/VTK/VOF'//TRIM(int2text(nf,padding))//'-'
    if(rank==0) call append_VOF_visit_file(TRIM(rootname))

    OPEN(UNIT=8,FILE=TRIM(rootname)//TRIM(int2text(rank,padding))//'.vtk')
    write(8,10)
    write(8,11) time
    write(8,12)
    write(8,13)
    write(8,14)i2-i1+1,j2-j1+1,k2-k1+1
    write(8,15)(i2-i1+1)*(j2-j1+1)*(k2-k1+1)
10  format('# vtk DataFile Version 2.0')
11  format('grid, time ',F16.8)
12  format('ASCII')
13  format('DATASET STRUCTURED_GRID')
14  format('DIMENSIONS ',I5,I5,I5)
15  format('POINTS ',I17,' float' )

    do k=k1,k2; do j=j1,j2; do i=i1,i2;
      write(8,320) x(i),y(j),z(k)
    enddo; enddo; enddo
320 format(e14.5,e14.5,e14.5)

    write(8,16)(i2-i1+1)*(j2-j1+1)*(k2-k1+1)
    write(8,17)'VOF'
    write(8,18)
16  format('POINT_DATA ',I17)
17  format('SCALARS ',A20,' float 1')
18  format('LOOKUP_TABLE default')

    do k=k1,k2; do j=j1,j2; do i=i1,i2;
      write(8,210) cvof(i,j,k)
    enddo; enddo; enddo
210 format(e14.5)
! 310 format(e14.5,e14.5,e14.5)
    close(8)
end subroutine output_VOF
!=================================================================================================
!-------------------------------------------------------------------------------------------------
end module module_output_vof
!-------------------------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------------------------
subroutine swp(us,c,vof1,vof2,vof3,f,d)
  use module_vof
  implicit none
  real (8)  , dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: us
  integer, intent(in) :: d
  real (8)  , dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: c,vof1,vof2,vof3
  integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: f
  if (VOF_advect=='Dick_Yue') then
     call swpr(us,c,vof1,vof2,vof3,f,d)
  elseif (VOF_advect=='CIAM') then
     call swpz(us,c,vof1,vof2,vof3,f,d)
  else
     call pariserror("*** unknown vof scheme")
  endif
  call get_flags_and_clip()
end subroutine swp
!
!  Implements the CIAM (Lagrangian Explicit, onto square)
!  advection method of Jie Li. 
! 
! ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
! split advection of the interface along the x (d=1), y (d=2) and z (d=3)
! directions
! ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
subroutine swpz(us,c,vof1,vof2,vof3,f,d)
  !***
  use module_grid
  use module_flow
  use module_vof
  implicit none
  integer i,j,k,invx,invy,invz
  real (8)  , dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: us
  integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: f
  integer, intent(in) :: d
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: c,vof1,vof2,vof3
  real(8) dmx,dmy,dmz,mm1,mm2
  real(8) a1,a2,alpha,al3d,fl3d
  real(8) mxyz(3),stencil3x3(-1:1,-1:1,-1:1)
  integer i0,j0,k0
  intrinsic dmax1,dmin1
  !***
  if(ng.lt.2) call pariserror("wrong ng")
  do k=ks-1,ke+1
     do j=js-1,je+1
        do i=is-1,ie+1
           if (d.eq.1) then
              a2 = us(i,j,k)*dt/dxh(i)
              a1 = us(i-1,j,k)*dt/dxh(i-1)
           elseif (d.eq.2) then
              a2 = us(i,j,k)*dt/dyh(j)
              a1 = us(i,j-1,k)*dt/dyh(j-1)
           elseif (d.eq.3) then
              a2 = us(i,j,k)*dt/dzh(k)
              a1 = us(i,j,k-1)*dt/dzh(k-1)
           endif
           !***
           !     3 cases: 1: default (c=0. and fluxes=0.); 2: c=1.; 3:c>0.
           !***
           vof1(i,j,k) = 0.0d0
           vof2(i,j,k) = 0.0d0
           vof3(i,j,k) = 0.0d0

           ! we need to introduce full/empty flags

           if (c(i,j,k) .eq. 1.0d0) then
              vof1(i,j,k) = dmax1(-a1,0.d0)
              vof2(i,j,k) = 1.d0 - dmax1(a1,0.d0) + dmin1(a2,0.d0)
              vof3(i,j,k) = dmax1(a2,0.d0)

           else if (c(i,j,k) .gt. 0.d0) then
              !***
              !     (1) normal vector: dmx,dmy,dmz, and |dmx|+|dmy|+|dmz| = 1.
              !     (2) dmx,dmy,dmz>0.
              !     (3) get alpha;               (4) back to original plane;
              !     (5) lagrangian advection;    (6) get fluxes
              !*(1)*

              do i0=-1,1; do j0=-1,1; do k0=-1,1
                 stencil3x3(i0,j0,k0) = c(i+i0,j+j0,k+k0)
              enddo;enddo;enddo
              call mycs(stencil3x3,mxyz)
              dmx = mxyz(1)
              dmy = mxyz(2)
              dmz = mxyz(3)
              !*(2)*  
              invx = 1
              invy = 1
              invz = 1
              if (dmx .lt. 0.0d0) then
                 dmx = -dmx
                 invx = -1
              endif
              if (dmy .lt. 0.0d0) then
                 dmy = -dmy
                 invy = -1
              endif
              if (dmz .lt. 0.0d0) then
                 dmz = -dmz
                 invz = -1
              endif
              !*(3)*  
              alpha = al3d(dmx,dmy,dmz,c(i,j,k))
              !*(4)*  
              dmx = invx*dmx
              dmy = invy*dmy
              dmz = invz*dmz
              alpha = alpha + dmin1(0.d0,dmx) + dmin1(0.d0,dmy) + &
                   dmin1(0.d0,dmz)
              !*(5)*  
              mm1 = dmax1(a1,0.0d0)
              mm2 = 1.d0 - mm1 + dmin1(0.d0,a2)
              if (d.eq.1) then
                 dmx = dmx/(1.0d0 - a1 + a2)
                 alpha = alpha + dmx*a1
                 if (a1 .lt. 0.d0) &
                      vof1(i,j,k) = fl3d(dmx,dmy,dmz,alpha,a1  ,-a1)
                 if (a2 .gt. 0.d0) &
                      vof3(i,j,k) = fl3d(dmx,dmy,dmz,alpha,1.d0,a2)
                 vof2(i,j,k) = fl3d(dmx,dmy,dmz,alpha,mm1,mm2)
              elseif (d.eq.2) then
                 dmy = dmy/(1.0d0 - a1 + a2)
                 alpha = alpha + dmy*a1
                 if (a1 .lt. 0.d0) &
                      vof1(i,j,k) = fl3d(dmy,dmz,dmx,alpha,a1  ,-a1)
                 if (a2 .gt. 0.d0) &
                      vof3(i,j,k) = fl3d(dmy,dmz,dmx,alpha,1.d0,a2)
                 vof2(i,j,k) = fl3d(dmy,dmz,dmx,alpha,mm1,mm2)
              elseif (d.eq.3) then
                 dmz = dmz/(1.0d0 - a1 + a2)
                 alpha = alpha + dmz*a1
                 if (a1 .lt. 0.d0) &
                      vof1(i,j,k) = fl3d(dmz,dmx,dmy,alpha,a1  ,-a1)
                 if (a2 .gt. 0.d0) &
                      vof3(i,j,k) = fl3d(dmz,dmx,dmy,alpha,1.d0,a2)
                 vof2(i,j,k) = fl3d(dmz,dmx,dmy,alpha,mm1,mm2)
              endif
              !           elseif (c(i,j,k).ne.0.d0) then
              !              call pariserror("case not allowed")
           endif
        enddo
     enddo
  enddo
  !
  ! assume that ghost layers take care of the boundary conditions. 
  ! so i-1, i+1 needs to be computed. 
  ! at least the ghost layers is-2, is-1, ie+1,ie+2  need to be there
  ! at the beginning of the subroutine, so that fluxes vof1,vof3 are computed
  ! for is-1, ie+1. 
  !    (1) new values of c and  clip it: 0. <= c <= 1.
  !    (2) apply proper boundary conditions to c
  !*(1)* 
  do k=ks,ke
     do j=js,je
        do i=is,ie
           if (d.eq.1) then
              c(i,j,k) = vof3(i-1,j,k) + vof2(i,j,k) + vof1(i+1,j,k)
           elseif (d.eq.2) then
              c(i,j,k) = vof3(i,j-1,k) + vof2(i,j,k) + vof1(i,j+1,k)
           elseif (d.eq.3) then
              c(i,j,k) = vof3(i,j,k-1) + vof2(i,j,k) + vof1(i,j,k+1)
           endif
           c(i,j,k) = dmax1(0.0d0,dmin1(1.0d0,c(i,j,k)))
        enddo
     enddo
  enddo
  !*(2)*
  call setvofbc(c,f)
  !***
end subroutine swpz
!
!=================================================================================================
! split 1D advection of the interface along the x,y,z (d=1,2,3) directions
!
! Following the advection method of Weymouth & Yue : Weymouth, G D, and Dick K P Yue,
! "Conservative Volume-of-Fluid Method for Free-Surface Simulations on Cartesian-Grids."
! Journal of Computational Physics 229, no. 8 (April 2010): 2853-2865. doi:10.1016/j.jcp.2009.12.018.
!=================================================================================================
!=================================================================================================
SUBROUTINE swpr(us,c,vof1,cg,vof3,f,dir)
!***
    USE module_grid
    USE module_flow
    USE module_vof
    use module_hello
    IMPLICIT NONE
    include 'mpif.h'
    INTEGER :: i,j,k
    INTEGER :: invx,invy,invz,ii,jj,kk,i0,j0,k0
    INTEGER, INTENT(IN) :: dir
    REAL (8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), INTENT(IN) :: us,cg
    REAL (8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), INTENT(INOUT) :: c,vof1,vof3
    integer, dimension(imin:imax,jmin:jmax,kmin:kmax),  intent(inout) :: f
    REAL(8), TARGET :: dmx,dmy,dmz,dxyz
    REAL(8), POINTER :: dm1,dm2,dm3
    REAL(8) :: a1,a2,alpha,AL3D,FL3D
    real(8) :: mxyz(3),stencil3x3(-1:1,-1:1,-1:1)
    INTRINSIC DMAX1,DMIN1
!
  if(ng < 2) call pariserror("wrong ng")
  ii=0; jj=0; kk=0
  if (dir == 1) then
     ii=1; dm1 => dmx;  dm2 => dmy;  dm3 => dmz 
  else if (dir == 2) then
     jj=1; dm1 => dmy;  dm2 => dmz;  dm3 => dmx 
  else if (dir == 3) then
     kk=1; dm1 => dmz;  dm2 => dmx;  dm3 => dmy 
  endif
  dxyz = dxh(is)
  call test_cell_size()

  do k=ks-1,ke+1
     do j=js-1,je+1
        do i=is-1,ie+1
           a2 = us(i,j,k)*dt/dxyz
           a1 = us(i-ii,j-jj,k-kk)*dt/dxyz
           !  default: fluxes=0. (good also for c=0.)
           vof1(i,j,k) = 0.d0
           vof3(i,j,k) = 0.d0
           !  c = 1.
           if (c(i,j,k) == 1.0d0) then
              vof1(i,j,k) = DMAX1(-a1,0.d0)
              vof3(i,j,k) = DMAX1(a2,0.d0)
           ! 0. < c < 1.
           else if (c(i,j,k) > 0.d0) then
              ! local stencil and normal vector: |dmx|+|dmy|+|dmz| = 1.
              do i0=-1,1; do j0=-1,1; do k0=-1,1
                 stencil3x3(i0,j0,k0) = c(i+i0,j+j0,k+k0)
              enddo;enddo;enddo
              call mycs(stencil3x3,mxyz)
              dmx = mxyz(1); dmy = mxyz(2); dmz = mxyz(3)
              ! positive dmx,dmy,dmz
              invx = 1; invy = 1; invz = 1
              if (dmx < 0.d0) then
                 dmx = -dmx; invx = -1
              endif
              if (dmy < 0.d0) then
                 dmy = -dmy; invy = -1
              endif
              if (dmz < 0.0d0) then
                 dmz = -dmz; invz = -1
              endif
              ! get alpha
              alpha = AL3D(dmx,dmy,dmz,c(i,j,k))
              ! back to the original plane
              dmx = invx*dmx
              dmy = invy*dmy
              dmz = invz*dmz
              alpha = alpha + DMIN1(0.d0,dmx) + DMIN1(0.d0,dmy) + DMIN1(0.d0,dmz)
              ! Eulerian advection
                 if (a1 < 0.d0) &
                      vof1(i,j,k) = FL3D(dm1,dm2,dm3,alpha,0.d0,-a1)
                 if (a2 > 0.d0) &
                      vof3(i,j,k) = FL3D(dm1,dm2,dm3,alpha,1.d0-a2,a2)
           endif
        enddo
     enddo
  enddo
  ! assume that ghost layers take care of the boundary conditions, then 
  ! fluxes vof1,vof3 must be computed for is-1, ie+1 
  ! new clipped values of c (0. <= c <= 1)
  do k=ks,ke
     do j=js,je
        do i=is,ie
           a2 = us(i,j,k)*dt/dxyz
           a1 = us(i-ii,j-jj,k-kk)*dt/dxyz
           c(i,j,k) = c(i,j,k) - (vof3(i,j,k) - vof1(i+ii,j+jj,k+kk)) + & 
                      (vof3(i-ii,j-jj,k-kk) - vof1(i,j,k)) + cg(i,j,k)*(a2-a1);
!!$           c(i,j,k) = DMAX1(0.d0,DMIN1(1.d0,c(i,j,k)))
           if (c(i,j,k) < EPSC) then
              c(i,j,k) = 0.d0
           elseif (c(i,j,k) >  (1.d0 - EPSC)) then
              c(i,j,k) = 1.d0
           endif
        enddo
     enddo
  enddo
  ! apply proper boundary conditions to c
  call setvofbc(c,f)
end subroutine swpr
!=================================================================================================
!=================================================================================================
! ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
! PROGRAM TO FIND alpha IN: m1 x1 + m2 x2 + m3 x3 = alpha,
! GIVEN m1+m2+m3=1 (all > 0) AND THE VOLUMETRIC FRACTION cc
! ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
function al3d(b1,b2,b3,cc)
  !***
  implicit none
  real(8) m1,m2,m3,cc,b1,b2,b3,tmp,pr,ch,mm,m12
  real(8) p,p12,q,teta,cs,al3d
  real(8) untier,v1,v2,v3
  parameter (untier=1.d0/3.d0)
  intrinsic dmax1,dmin1,dsqrt,dacos,dcos
  !***  
  !     (1) order coefficients: m1<m2<m3; (2) get ranges: v1<v2<v3;
  !     (3) limit ch (0.d0 < ch < 0.5d0); (4) calculate alpha
  !*(1)* 
  m1 = dmin1(b1,b2)
  m3 = dmax1(b1,b2)
  m2 = b3
  if (m2 .lt. m1) then
     tmp = m1
     m1 = m2
     m2 = tmp
  else if (m2 .gt. m3) then
     tmp = m3
     m3 = m2
     m2 = tmp
  endif
  !*(2)*
  m12 = m1 + m2 
  pr  = DMAX1(6.d0*m1*m2*m3,1.d-50)
  V1  = m1*m1*m1/pr
  V2  = V1 + 0.5d0*(m2-m1)/m3
  if (m3 .LT. m12) then
     mm = m3
     V3 = (m3*m3*(3.d0*m12-m3) + m1*m1*(m1-3.d0*m3) +&
          m2*m2*(m2-3.d0*m3))/pr
  else
     mm = m12
     V3 = 0.5d0*mm/m3
  endif
  !*(3)*
  ch = DMIN1(cc,1.d0-cc)
  !*(4)*      
  if (ch .LT. V1) then
     !***         AL3D = cbrt(pr*ch)
     AL3D = (pr*ch)**UNTIER
  else if (ch .LT. V2) then
     AL3D = 0.5d0*(m1 + DSQRT(m1*m1 + 8.d0*m2*m3*(ch-V1)))
  else if (ch .LT. V3) then
     p = 2.d0*m1*m2
     q = 1.5d0*m1*m2*(m12 - 2.d0*m3*ch)
     p12 = DSQRT(p)
     teta = DACOS(q/(p*p12))/3.d0
     cs = DCOS(teta)
     AL3D = p12*(DSQRT(3.d0*(1.d0-cs*cs)) - cs) + m12
  else if (m12 .LT. m3) then
     AL3D = m3*ch + 0.5d0*mm
  else 
     p = m1*(m2+m3) + m2*m3 - 0.25d0
     q = 1.5d0*m1*m2*m3*(0.5d0-ch)
     p12 = DSQRT(p)
     teta = DACOS(q/(p*p12))/3.0
     cs = DCOS(teta)
     AL3D = p12*(DSQRT(3.d0*(1.d0-cs*cs)) - cs) + 0.5d0
  endif

  if (cc .GT. 0.5d0)  AL3D = 1.d0 - AL3D
  !***
  return
end function al3d
! ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
! PROGRAM TO FIND THE "CUT VOLUME" V0 GIVEN r0, dr0 AND
! m1 x1 + m2 x2 + m3 x3 = alpha
! ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
function fl3d(m1,m2,m3,alpha,r0,dr0)
  !***
  implicit none
  real(8) m1,m2,m3,alpha,r0,dr0,fl3D
  real(8) al,al0,n1,n2,n3,b1,b2,b3,b12,bm,tmp,pr
  INTRINSIC DMAX1,DMIN1,DABS
  !***
  !     (1) move origin to r0 along r ;  (2) reflect parallelepiped;
  !     (3) limit alpha (0<= al0 <=0.5); (4) order coefficients: b1<b2<b3;
  !     (5) calculate volume (NOTE: it is assumed:s0=t0=0; ds0=dt0=1.)
  !*(1)*
  al = alpha - m1*r0
  !*(2)*
  al = al + DMAX1(0.d0,-m1*dr0)+DMAX1(0.d0,-m2)+DMAX1(0.d0,-m3)
  tmp = DABS(m1)*dr0 + DABS(m2) + DABS(m3)
  n1 = DABS(m1)/tmp
  n2 = DABS(m2)/tmp
  n3 = DABS(m3)/tmp
  al = DMAX1(0.d0,DMIN1(1.d0,al/tmp))
  !*(3)*
  al0 = DMIN1(al,1.d0-al)
  !*(4)*
  b1 = DMIN1(n1*dr0,n2)
  b3 = DMAX1(n1*dr0,n2)
  b2 = n3
  if (b2 .LT. b1) then
     tmp = b1
     b1 = b2
     b2 = tmp
  else if (b2 .GT. b3) then
     tmp = b3
     b3 = b2
     b2 = tmp
  endif
  b12 = b1 + b2
  bm = DMIN1(b12,b3)
  pr = DMAX1(6.d0*b1*b2*b3,1.0d-50)
  !*5*     
  if (al0 .LT. b1) then
     tmp = al0*al0*al0/pr
  else if (al0 .LT. b2) then
     tmp = 0.5d0*al0*(al0-b1)/(b2*b3) +  b1*b1*b1/pr
  else if (al0 .LT. bm) then
     tmp = (al0*al0*(3.d0*b12-al0) + b1*b1*(b1-3.d0*al0) +&
          b2*b2*(b2-3.d0*al0))/pr
  else if (b12 .LT. b3) then
     tmp = (al0 - 0.5d0*bm)/b3
  else
     tmp = (al0*al0*(3.d0-2.d0*al0) + b1*b1*(b1-3.d0*al0) +&
          b2*b2*(b2-3.d0*al0) + b3*b3*(b3-3.d0*al0))/pr
  endif

  if (al .LE. 0.5d0) then
     FL3D = tmp*dr0
  else
     FL3D = (1.d0-tmp)*dr0
  endif
  !***  
  return
end function fl3d

subroutine ls2vof_in_cell(stencil3x3,c,nflag)
  implicit none
  real(8), intent(out):: c
  integer, intent(out):: nflag
  real(8) :: zero, one, norml1
  real(8) :: mx,my,mz,alpha
  real(8) :: fl3d
  real(8) :: mxyz(3),stencil3x3(-1:1,-1:1,-1:1)

  zero=0.d0
  one=1.d0
  !***
  !     (1) gradient*32: mx,my,mz; (2) mx,my,mz>0. and mx+my+mz = 1.;
  !     (3) normalize alpha = level set at center. Cell units. 
  !     (4) shift alpha to origin=vertex;   (5) get volume from alpha.  
  !
  !     *(1)*  
  !***
  call fd32(stencil3x3,mxyz)
  !***
  !     *(2)*  
  !***
  mx = dabs(mxyz(1)); my = dabs(mxyz(2)); mz = dabs(mxyz(3))
  norml1 = mx+my+mz
  mx = mx/norml1;     my = my/norml1;     mz = mz/norml1
  !***
  !     *(3)*  
  !***
  ! the factor is 32 because grad ls=(1,0,0) gives mx=32.
  alpha = 32.d0*stencil3x3(0,0,0)/norml1   
  !***
  !     *(4)*  
  !***
  alpha = alpha + 0.5d0
  !***
  !     *(5)*  
  !***
  if(alpha.ge.1.d0) then 
     c = 1.d0
     nflag = 1
  else if (alpha.le.0.d0) then
     c = 0.d0
     nflag = 0 
  else 
     c = fl3d(mx,my,mz,alpha,zero,one)
     nflag = 2
  end if
  return
end subroutine ls2vof_in_cell
!
! *-----------------------------------------------------* 
! *  MYC - Mixed Youngs and Central Scheme              *
! * returns normal normalized so that |mx|+|my|+|mz| = 1* 
! *-----------------------------------------------------*
! 
!
!Known problems: the index (1,1,1), i.e. the central cell
!in the block, never occurs: neither in the central scheme
!nor in Youngs' method. Therefore an isolated droplet will have
!a normal with all components to zero. I took care of the
!division-by-zero issue, but not of this one.
!
!Ruben
!
!
! Translated into f90 by Stephane Z.
!
subroutine mycs(c,mxyz)
  !***
  implicit none
  real(8) c(0:2,0:2,0:2)
  real(8) mxyz(0:2)
  real(8) m1,m2,m(0:3,0:2),t0,t1,t2
  integer cn
  real(8), parameter  :: NOT_ZERO=1.e-30

  ! write the plane as: sgn(mx) X =  my Y +  mz Z + alpha 
  !                           m00 X = m01 Y + m02 Z + alpha 

  m1 = c(0,1,0) + c(0,1,2) + c(0,0,1) + c(0,2,1) + &
       c(0,1,1)
  m2 = c(2,1,0) + c(2,1,2) + c(2,0,1) + c(2,2,1) + &
       c(2,1,1)

  if(m1>m2) then
     m(0,0) = 1.
  else
     m(0,0) = -1.
  end if

  m1 = c(0,0,1)+ c(2,0,1)+ c(1,0,1)
  m2 = c(0,2,1)+ c(2,2,1)+ c(1,2,1)
  m(0,1) = 0.5*(m1-m2)

  m1 = c(0,1,0)+ c(2,1,0)+ c(1,1,0)
  m2 = c(0,1,2)+ c(2,1,2)+ c(1,1,2)
  m(0,2) = 0.5*(m1-m2)

  ! write the plane as: sgn(my) Y =  mx X +  mz Z + alpha, 
  !                          m11 Y = m10 X + m12 Z + alpha.

  m1 = c(0,0,1) + c(0,2,1) + c(0,1,1)
  m2 = c(2,0,1) + c(2,2,1) + c(2,1,1)
  m(1,0) = 0.5*(m1-m2)

  m1 = c(1,0,0) + c(1,0,2) + c(2,0,1) + c(0,0,1) +&
       c(1,0,1)
  m2 = c(1,2,0) + c(1,2,2) + c(2,2,1) + c(0,2,1) +&
       c(1,2,1)


  if(m1>m2) then
     m(1,1) = 1.
  else
     m(1,1) = -1.
  end if

  m1 = c(1,0,0)+ c(1,1,0)+ c(1,2,0)
  m2 = c(1,0,2)+ c(1,1,2)+ c(1,2,2)
  m(1,2) = 0.5*(m1-m2)

  ! write the plane as: sgn(mz) Z =  mx X +  my Y + alpha 
  !                          m22 Z = m20 X + m21 Y + alpha

  m1 = c(0,1,0)+ c(0,1,2)+ c(0,1,1)
  m2 = c(2,1,0)+ c(2,1,2)+ c(2,1,1)
  m(2,0) = 0.5*(m1-m2)

  m1 = c(1,0,0)+ c(1,0,2)+ c(1,0,1)
  m2 = c(1,2,0)+ c(1,2,2)+ c(1,2,1)
  m(2,1) = 0.5*(m1-m2)

  m1 = c(0,1,0) + c(2,1,0) + c(1,0,0) + c(1,2,0) +&
       c(1,1,0)
  m2 = c(0,1,2) + c(2,1,2) + c(1,0,2) + c(1,2,2) +&
       c(1,1,2)

  if(m1>m2) then
     m(2,2) = 1.
  else
     m(2,2) = -1.
  end if

  ! normalize each set (mx,my,mz): |mx|+|my|+|mz| = 1

  t0 = DABS(m(0,0)) + DABS(m(0,1)) + DABS(m(0,2))
  m(0,0) = m(0,0)/t0
  m(0,1) = m(0,1)/t0
  m(0,2) = m(0,2)/t0

  t0 = DABS(m(1,0)) + DABS(m(1,1)) + DABS(m(1,2))
  m(1,0) = m(1,0)/t0
  m(1,1) = m(1,1)/t0
  m(1,2) = m(1,2)/t0

  t0 = DABS(m(2,0)) + DABS(m(2,1)) + DABS(m(2,2))
  m(2,0) = m(2,0)/t0
  m(2,1) = m(2,1)/t0
  m(2,2) = m(2,2)/t0

  ! choose among the three central schemes */ 
  t0 = DABS(m(0,0))
  t1 = DABS(m(1,1))
  t2 = DABS(m(2,2))

  cn = 0
  if (t1 > t0) then
    t0 = t1
    cn = 1
  endif

  if (t2 > t0) cn = 2

  ! Youngs-CIAM scheme */  
  
  call fd32(c,m(3,0:2))

  ! normalize the set (mx,my,mz): |mx|+|my|+|mz| = 1 

  t0 = DABS(m(3,0)) + DABS(m(3,1)) + DABS(m(3,2)) + NOT_ZERO
  m(3,0) = m(3,0)/t0
  m(3,1) = m(3,1)/t0
  m(3,2) = m(3,2)/t0

  ! choose between the previous choice and Youngs-CIAM 
  t0 = DABS (m(3,0))
  t1 = DABS (m(3,1))
  t2 = DABS (m(3,2))
  if (t1 > t0)  t0 = t1
  if (t2 > t0)  t0 = t2

  if (DABS(m(cn,cn)) > t0)  cn = 3

  ! components of the normal vector */
  mxyz(0) = m(cn,0)
  mxyz(1) = m(cn,1)
  mxyz(2) = m(cn,2)

  return 
  end subroutine mycs

! *----------------------------------------------------------------* 
! *  FD32 - Youngs Finite Difference Gradient Scheme               *
! *  the gradient is computed with a multiplicative factor of -32: *
! *  mm = - 32 * grad (c)                                          *
! *----------------------------------------------------------------*
! 
!Known problems: the index (1,1,1), i.e. the central cell
!in the block, never occurs:
!Therefore an isolated droplet will have
!a normal with all components to zero. 
!
!Ruben
!
!
! Translated into f90 by Stephane Z.
!
subroutine fd32(c,mm)
  !***
  implicit none
  real(8), intent(in) :: c(0:2,0:2,0:2)
  real(8), intent(out) :: mm(0:2)
  real(8) :: m1,m2

  m1 = c(0,0,0) + c(0,2,0) + c(0,0,2) + c(0,2,2) +&
       2.d0*(c(0,0,1) + c(0,2,1) + c(0,1,0) + c(0,1,2)) +&
       4.d0*c(0,1,1)
  m2 = c(2,0,0) + c(2,2,0) + c(2,0,2) + c(2,2,2) +&
       2.d0*(c(2,0,1) + c(2,2,1) + c(2,1,0) + c(2,1,2)) +&
       4.d0*c(2,1,1)
  mm(0) = m1-m2

  m1 = c(0,0,0) + c(0,0,2) + c(2,0,0) + c(2,0,2) +&
       2.d0*(c(0,0,1) + c(2,0,1) + c(1,0,0) + c(1,0,2)) +&
       4.d0*c(1,0,1)
  m2 = c(0,2,0) + c(0,2,2) + c(2,2,0) + c(2,2,2) +&
       2.d0*(c(0,2,1) + c(2,2,1) + c(1,2,0) + c(1,2,2)) +&
       4.d0*c(1,2,1)
  mm(1) = m1-m2

  m1 = c(0,0,0) + c(0,2,0) + c(2,0,0) + c(2,2,0) +&
       2.d0*(c(0,1,0) + c(2,1,0) + c(1,0,0) + c(1,2,0)) +&
       4.d0*c(1,1,0)
  m2 = c(0,0,2) + c(0,2,2) + c(2,0,2) + c(2,2,2) +&
       2.d0*(c(0,1,2) + c(2,1,2) + c(1,0,2) + c(1,2,2)) +&
       4.d0*c(1,1,2)
  mm(2) = m1-m2

  return 
  end subroutine fd32

!
! *----------------------------------------------------------------* 
! *  youngs - Youngs Finite Difference Gradient Scheme             *
! *  the gradient is normed so that |mx|+|my|+|mz| = 1             *
! *----------------------------------------------------------------*
! 
! Known problems: the index (1,1,1), i.e. the central cell
! in the block, never occurs:
! Therefore an isolated droplet will have
! a normal with all components to zero. 
!
!
subroutine youngs(c,mm)
  !***
  implicit none
  real(8), intent(in) :: c(0:2,0:2,0:2)
  real(8), intent(out) :: mm(0:2)
  integer :: i
  real(8) :: norm

  call fd32(c,mm)
  norm = abs(mm(0))+abs(mm(1))+abs(mm(2)) ! + TINY
  do i=0,2
     mm(i) = mm(i)/norm
  enddo
  return 
  end subroutine youngs
!
