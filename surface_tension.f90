!=================================================================================================
!=================================================================================================
! PARIS  Parallel Robust Interface Simulator 
!=================================================================================================
! module_surface_tension: Contains definition of variables for surface tension from
!  Volume of Fluid interface tracking.
!
! Contact: Stephane Zaleski zaleski@dalembert.upmc.fr
! 
! Authors:
!         Daniel Fuster
! 	  Yue "Stanley" Ling 
!         Leon Malan
!         Ruben Scardovelli  
!         Phil Yecko         
!         Stephane Zaleski   
!
! GPL Licence
!
!     This file is part of PARIS.
!
!     PARIS is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
! 
!     PARIS is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with PARIS.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------------------------
module module_surface_tension
  use module_grid
  use module_BC
  use module_IO
  use module_2phase
  use module_freesurface
  use module_VOF
  implicit none

  ! choice of method  
  logical :: recomputenormals = .true.

  ! initial state
  logical :: st_initialized = .false.
  real(8), parameter :: kappamax = 4.d0
  integer, parameter :: nfound_min =  6 
  ! Caution the line below is read by test scripts, so do not move it below line 60. Do not 
  ! write ndepth= anywhere else in the first 60 lines. 
  ! Do not change the lettercase for ndepth. 
  integer, parameter :: NDEPTH=3 
  integer, parameter :: BIGINT=100
  real(8), parameter :: D_HALF_BIGINT = DBLE(BIGINT/2)
  integer, parameter :: MAX_EXT_H = 1
  integer, parameter :: NOR=6 ! number of orientations
  integer, parameter :: NPOS=27+1
  real(8), parameter :: EPS_GEOM = 1d-4
  real(8), dimension(:,:,:), allocatable :: n1,n2,n3 ! normals
  real(8), dimension(:,:,:,:), allocatable :: height ! 
  real(8), parameter :: UNCOMPUTED=2D6

  ! 4th index: 1 for normal vector pointing towards positive x "positive height", 
  ! 2 for "negative" height in x
  ! 3 for positive height in y, 4 for negative height in y, 
  !  etc... 
  integer, dimension(:,:,:,:), allocatable :: ixheight ! Height-Function flags for Ruben-Phil routines
  integer :: method_count(5)
  integer, parameter :: NGC=21
  integer :: geom_case_count(NGC)
  integer :: nkcomp
  integer, parameter :: NFOUND_PURE=NGC+1

contains
  !=================================================================================================
  subroutine initialize_surface_tension()
    implicit none
    integer :: minM
    allocate(height(imin:imax,jmin:jmax,kmin:kmax,6))
    if(nx.ge.500000.or.ny.gt.500000.or.nz.gt.500000) call pariserror("nx too large")
    if(NDEPTH.gt.20) call pariserror("ndepth too large")
    if(NDEPTH>BIGINT/2-2) call pariserror("BIGINT too small")
    minM=min(Mx,My,Mz)
    ! detect thin slice
    if(Ny<=2) minM=min(Mx,Mz)
    if(Nz<=2) minM=min(Mx,My)
    ! thin slice expected in y or z. 
    if(MinM<=2) call pariserror("unexpected thin slice")
    ! detect small subdomains in which the stencil will straddle two boundaries. 
    if(2*NDEPTH+1.gt.minM) then
       if(rank==0) then
          print *, 'Mx=',Mx,'My=',My,'Mz=',Mz,'minM=',minM
       endif
       call pariserror("ndepth too large for current box size nx/npx")
    endif
    if(MAX_EXT_H>BIGINT/2-2) call pariserror("MAX_EXT > BIGINT/2")
    if(MAX_EXT_H>nx/2) call pariserror("MAX_EXT > nx/2")
    !    allocate(geom_case_list(10))
    geom_case_count = 0
    height = 2.d6
    st_initialized=.true.
  end subroutine initialize_surface_tension
  !
  subroutine print_st_stats(iout)
    use module_BC
    use module_IO
    implicit none
    include 'mpif.h'
    integer :: ierr,i,iout
    integer :: glob_count(ngc)
    character(len=85) :: glob_desc(ngc)
    if(st_initialized) then
       call MPI_ALLREDUCE(geom_case_count, glob_count, ngc, MPI_INTEGER, MPI_SUM, MPI_COMM_Cart, ierr)  
       if(rank==0) then
          open(unit=101, file=trim(out_path)//'/st_stats-'//TRIM(int2text(iout,padding)), action='write', iostat=ierr)
          glob_desc(1)="mixed w/less than 4 mixed neighbors (quasi-isolated mixed, unfittable by sphere)"
          glob_desc(2)="mixed w/less than 6 mixed neighbors (quasi-isolated mixed, unfittable by paraboloid)"
          glob_desc(3)="pure cells w/more than 2 other color pure neighbors (grid-aligned interfaces)"
          glob_desc(4)="non-bulk pure cells w 0 valid neighbors"
          glob_desc(5)="                      1  " 
          glob_desc(6)="                      2  " 
          glob_desc(7)="                      3  " 
          glob_desc(8)="                      4  "
          glob_desc(9)="                      5 valid neighbors (unfittable by paraboloid)"
          glob_desc(10)="no fit success with 0 points"
          glob_desc(11)="                    1  " 
          glob_desc(12)="                    2  " 
          glob_desc(13)="                    3  " 
          glob_desc(14)="                    4  " 
          glob_desc(15)="                    5 points"
          glob_desc(16)="                    6 or more points"
          glob_desc(17)="large kappa"
          glob_desc(18)="no surface tension force in x direction"
          glob_desc(19)="no surface tension force in y direction"
          glob_desc(20)="no surface tension force in z direction"
          glob_desc(21)="mxyz vector is TINY in get_curvature (indicates debris cell)"
          do i=1,ngc
             write(101,'(I10," ",A85)') glob_count(i), glob_desc(i)
          enddo
          close(101)
       endif
       geom_case_count=0
    endif
  end subroutine print_st_stats

  !=================================================================================================
  ! 
  !  Put normals in a common array. Absolutely not sure this is efficient
  !
  !=================================================================================================
  subroutine get_normals()
    implicit none
    include 'mpif.h'
    real(8) :: stencil3x3(-1:1,-1:1,-1:1)
    integer :: i,j,k,ierr
    integer :: i0,j0,k0
    real(8) :: mxyz(3)
    integer :: req(36),sta(MPI_STATUS_SIZE,36)
    if(.not.st_initialized) call initialize_surface_tension()
    if(recomputenormals) call pariserror("recomputenormals is true, normals not allocated")

    if(ng.lt.2) call pariserror("wrong ng")
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
    call do_ghost_vector(n1,n2,n3)
  end subroutine get_normals
  !=================================================================================================
  !
  ! the core of HF computation
  !
  !=================================================================================================
  subroutine get_all_heights(iout)
    !use module_grid
    !use module_IO
    use module_timer
    implicit none
    include 'mpif.h'
    integer :: direction, ierr, i
    integer :: req(24),sta(MPI_STATUS_SIZE,24)
    logical :: debug
    integer :: iout,l,j,k
    integer :: prank
    if(.not.st_initialized) call initialize_surface_tension()

    !*** Initialize
    height=2d6

    do direction=1,3
       call get_heights_pass1(direction)
    enddo
    call my_timer(5)

    do i=1,6
       call ghost_x(height(:,:,:,i),2,req(4*(i-1)+1:4*i))
    enddo
    call MPI_WAITALL(24,req(1:24),sta(:,1:24),ierr)
    do i=1,6
       call ghost_y(height(:,:,:,i),2,req(4*(i-1)+1:4*i))
    enddo
    call MPI_WAITALL(24,req(1:24),sta(:,1:24),ierr)
    do i=1,6
       call ghost_z(height(:,:,:,i),2,req(4*(i-1)+1:4*i))
    enddo
    call MPI_WAITALL(24,req(1:24),sta(:,1:24),ierr)
    call my_timer(6)
    do direction=1,3
       call get_heights_pass2(direction)
       !        call get_heights_pass3(direction)
    enddo
    call my_timer(5)
    do i=1,6
       call ghost_x(height(:,:,:,i),2,req(4*(i-1)+1:4*i))
    enddo
    call MPI_WAITALL(24,req(1:24),sta(:,1:24),ierr)
    do i=1,6
       call ghost_y(height(:,:,:,i),2,req(4*(i-1)+1:4*i))
    enddo
    call MPI_WAITALL(24,req(1:24),sta(:,1:24),ierr)
    do i=1,6
       call ghost_z(height(:,:,:,i),2,req(4*(i-1)+1:4*i))
    enddo
    call MPI_WAITALL(24,req(1:24),sta(:,1:24),ierr)
    call my_timer(6)

    !! Parallel debugging output for pariscompare3D
    if (debug_par) then
       do l=1,6
          if (rank==0) then
             OPEN(UNIT=21,FILE='Height'//TRIM(int2text(l,1))//'holder-'//TRIM(int2text(iout,padding))//'.txt')
             write(21,12)NpDomain
             do prank=0,NpDomain-1
                write(21,13)TRIM(out_path)//'/Height'//TRIM(int2text(l,1))//'-'//&
                     TRIM(int2text(prank,padding))//'-'//TRIM(int2text(iout,padding))//'.txt'
             enddo
          endif
          OPEN(UNIT=20,FILE=TRIM(out_path)//'/Height'//TRIM(int2text(l,1))//'-'//&
               TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout,padding))//'.txt')
          do k=ks,ke; do j=js,je; do i=is,ie
             write(20,14)x(i),y(j),z(k),height(i,j,k,l)
          enddo;enddo;enddo
          close(20)
       enddo
    endif
12  format('NPROCS',I4)
13  format(A)  
14  format(4e14.5)
  end subroutine get_all_heights
  !=================================================================================================
  ! 
  !   the actual HF
  ! 
  !=================================================================================================
  subroutine get_heights_pass1(d)
    implicit none
    integer, intent(in) :: d
    integer :: index
    logical :: same_flag, abandon_search, height_found_at_start
    real(8) :: height_p     !  partial height 
    integer :: i,j,k,s,c0,c1,c(3),ctop
    integer :: sign, flag_other_end, climitp1, normalsign
    ! NDEPTH is the depth of layers tested above or below the reference cell. 
    ! including the central layer and the empty/full cells
    ! NDEPTH*2 + 1 = 7 means a 7 x 3^2 stencil including full/empty. 
    ! NDEPTH*2 + 1 = 5 (ndepth = 2) is the usual 3x3 stencil with checking of the cells
    ! above and below
    ! 2*NDEPTH + 1 = N_H + 2 in the notation of Mark Owkes & Olivier Desjardins
    ! JCP, Volume 281, 2015, Pages 285-300. 
    !  Note the normal is - grad C

    do k=ks,ke; do j=js,je; do i=is,ie
       if(vof_flag(i,j,k)/2==0) then ! flag is 0 or 1
          ! loop over search directions. sign = -1 search down. sign = +1 search up. 
          do sign=-1,1,2
             c(1)=i; c(2)=j; c(3)=k
             !  vof_flag=1 and sign = +  positive normal orientation 
             !  vof_flag=1 and sign = -  negative normal orientation 
             !  vof_flag=0 and sign = +  negative normal orientation 
             !  vof_flag=0 and sign = -  positive normal orientation 
             normalsign = (2*vof_flag(i,j,k)-1) * sign
             !  index: 2*(d-1) + 1 for normal pointing up (reference phase under the other phase)
             !  index: 2*(d-1) + 2 for normal pointing down
             index = 2*(d-1) + 1 + (-normalsign+1)/2
             flag_other_end = 1 - vof_flag(i,j,k)
             climitp1 = coordlimit(d,sign) + sign   ! first ghost layer
             height_p = 0.d0
             s = 0
             c0 = c(d) ! start of stack
             c1 = c0 + sign*ndepth ! middle of stack starting at c0 in direction sign and having maximum extent
             abandon_search=.false.
             !call verify_indices(c(1),c(2),c(3),index,0)
             height_found_at_start  = height(c(1),c(2),c(3),index)<D_HALF_BIGINT
             do while (.not.abandon_search) 
                !call verify_indices(c(1),c(2),c(3),index,1)
                !call verify_indices(i,j,k,1,2)
                same_flag = s>0.and.vof_flag(c(1),c(2),c(3))==vof_flag(i,j,k)
                height_p = height_p + (cvof(c(1),c(2),c(3)) - 0.5d0)*normalsign
                abandon_search =        &
                     same_flag          & ! case (1) no height, do nothing
                     .or.height_found_at_start                    & ! case (2) height already found, do nothing
                     .or.vof_flag(c(1),c(2),c(3))==flag_other_end & ! case (3) found the full height in previous pass
                     .or.c(d)==climitp1 & ! case (4) reached top but : not full height since checked above in (3)
                     .or.s==2*ndepth      ! case (5) stack too long : now c(d) = c0 + 2*ndepth
                if(.not.abandon_search) then
                   s = s + 1
                   c(d) = c(d) + sign ! go forward. Now c(d) = c0 + sign*s
                   ! abandon search, but first perform operations that record the height
                else if(vof_flag(c(1),c(2),c(3))==flag_other_end) then 
                   ! (3) *found the full height* !
                   ! there may be missing terms in the sum since the maximum extent
                   ! (s=2*ndepth) of the stack was not
                   ! necessarily reached. Add these terms. Here s = c(d) - c0
                   height_p = height_p + (2*ndepth-s)*(cvof(c(1),c(2),c(3))-0.5d0)*normalsign
                   ! height is now computed with respect to c1
                   do while (c(d)/=(c0-sign))
                      ! call verify_indices(c(1),c(2),c(3),index,3)
                      ! correct the height to give it with respect to current position c(d)
                      height(c(1),c(2),c(3),index) = height_p + c1 - c(d)
                      !                    call check_all(c(1),c(2),c(3),index)
                      c(d) = c(d) - sign ! go back down
                   enddo
                else if(c(d)==climitp1) then 
                   ! (4) reached top but : not full height since checked above in (3)
                   !     save partial height in the last normal layer (normal as opposite to ghost)
                   height_p = height_p + (- cvof(c(1),c(2),c(3)) + 0.5d0)*normalsign ! remove last addition
                   c(d) = c(d) - sign ! go back one step to climit
                   ! (**) here s is the algebraic distance to c0 thus
                   !      s = (c(d) - c0)*sign + 1 and s=1 for c(d)=c0=climit
                   !call verify_indices(c(1),c(2),c(3),index,4)
                   height(c(1),c(2),c(3),index) = height_p + BIGINT*s 
                   !                call check_all(c(1),c(2),c(3),index)
                endif ! abandon_search
             enddo ! while not abandon search
          enddo ! sign
       endif ! vof_flag
    enddo; enddo; enddo;  ! i,j,k
    !      contains
    !        subroutine verify_indices(i,j,k,index,pass)
    !          implicit none
    !          include 'mpif.h'
    !          integer, intent(in) :: i,j,k
    !          integer, intent(in) :: index,pass
    !          integer :: ierr, MPI_errorcode
    !          if(i.lt.imin.or.i.gt.imax.or.   &
    !               j.lt.jmin.or.j.gt.jmax.or. &
    !               k.lt.kmin.or.k.gt.kmax.or. &
    !               index.lt.1.or.index.gt.6) then 
    !             OPEN(UNIT=88,FILE=TRIM(out_path)//'/error-rank-'//TRIM(int2text(rank,padding))//'.txt')
    !             write(88,*) "imin,imax,jmin,jmax,kmin,kmax",imin,imax,jmin,jmax,kmin,kmax
    !             write(88,*) "i,j,k,index,pass",i,j,k,index,pass
    !             close(88)
    !             close(out)
    !             if(rank==0) print *, "index error in get_heights"
    !             call MPI_ABORT(MPI_COMM_WORLD, MPI_errorcode, ierr)
    !             call MPI_finalize(ierr)
    !             stop 
    !          end if
    !        end subroutine verify_indices
  end subroutine get_heights_pass1
  !
  !  Enable parallel computation: exchange information accross boundaries. 
  !
  subroutine get_heights_pass2(d)
    implicit none
    integer, intent(in) :: d
    integer :: index,i,j,k
    real(8) :: ha,hb,cmiddle
    integer :: l,m,n,c0,cb,c(3),try(3)
    integer :: sign, sabove, sbelow
    ! NDEPTH is the depth of layers tested above or below the reference cell. 
    try(1)=d 
    m=1
    n=2
    do while (m.le.3)
       if(m.ne.d) then
          try(n) = m
          n=n+1
       endif
       m=m+1
    enddo

    do l=coordstart(try(2)),coordend(try(2))
       do m=coordstart(try(3)),coordend(try(3))
          c(try(2)) = l; c(try(3)) = m
          do sign=-1,1,2  ! search in both directions
             do index = 2*(d-1) + 1, 2*(d-1) + 2  ! and for both indexes. 
                cb = coordlimit(d,sign)   ! coordinate "below" boundary
                c(d) = cb
                hb = height(c(1),c(2),c(3),index)
                if(hb>D_HALF_BIGINT.and.hb<1d6) then ! partial height in cell below
                   c(d) = cb + sign  
                   ha = height(c(1),c(2),c(3),index)
                   c(d) = cb
                   if(ha<D_HALF_BIGINT) then ! height already found above
                      height(c(1),c(2),c(3),index) = ha + sign   ! set height below accordingly  @@@ check no logical issues
                   else if(ha>D_HALF_BIGINT.and.ha<1d6) then ! try to match
                      sbelow = FLOOR(REAL(hb + D_HALF_BIGINT)/REAL(BIGINT)) ! "below" is assuming sign=1
                      hb = hb - BIGINT*sbelow  ! above, below in direction of sign
                      sabove = FLOOR(REAL(ha + D_HALF_BIGINT)/REAL(BIGINT))
                      ha = ha - BIGINT*sabove
                      ! c(d) = c0          index of bottom of stack
                      !        cmiddle     index of center of stack (for this stack, can be half integer)
                      !        ctop        index of top of this stack
                      !        ctop-c0+1 = length of stack
                      ! see (**) in pass 1. s is a positive quantity
                      !             cb-c0 =  sign*sbelow - sign
                      !             ca-ctop= - sign*sabove + sign
                      ! hence
                      !  |cb-c0| +  |ca-ctop| + 2 = ctop-c0 + 1 = sabove + sbelow
                      !  cmiddle = c0 + (ctop-c0)/2
                      if(sabove + sbelow - 1 <= 2*ndepth+1) then  ! 
                         ! bottom is at 
                         c0   = cb - (sbelow-1)*sign  
                         cmiddle   = dble(c0) + sign*(sabove+sbelow-1)*0.5d0
                         c(d) = cb + 2*sign 
                         do while (c(d)/=(c0-sign)) 
                            height(c(1),c(2),c(3),index) = ha + hb + cmiddle - c(d)
                            c(d) = c(d) - sign ! go back to c0 
                         enddo
                      endif ! not over stack height
                   endif ! partial height above
                endif ! partial height in cell below: if not, either full height or no-height, leave as is
                ! need to correct cell above accordingly if full height below and if correct height wanted in the 
                ! first ghost layer.
             enddo ! index
          enddo ! sign
       enddo ! l
    enddo ! m

    do index = 2*(d-1) + 1, 2*(d-1) + 2  ! for both indexes. 
       do k=kmin,kmax;do j=jmin,jmax;do i=imin,imax
          if(height(i,j,k,index)>D_HALF_BIGINT) height(i,j,k,index)=2d6
       enddo;enddo;enddo
    enddo

  end subroutine get_heights_pass2

  !    subroutine get_heights_pass3(d)   ! needs fixing
  !      implicit none
  !      integer, intent(in) :: d
  !      integer :: index
  !      logical :: limit_not_found
  !      integer :: i,j,k,c0,c(3)
  !      integer :: sign, climitp2, oppnormalsign
  !      ! need to extend heights
  !      ! start from full cells for which the height is defined
  !      ! and go the opposite way (towards the opposite interface); 
  !      do i=is-1,ie+1; do j=js-1,je+1; do k=ks-1,ke+1
  !         if(vof_flag(i,j,k)/2==0) then
  !            ! loop over search directions
  !            do sign=-1,1,2; 
  !               ! We want the opposite of search direction in pass 1 so 
  !               ! negative normal orientation if vof_flag=1 and sign = +, etc...
  !               ! oppnormalsign is the opposite of the sign of the normal
  !               oppnormalsign = - (2*vof_flag(i,j,k)-1) * sign
  !               index = 2*(d-1) + 1 + (-oppnormalsign+1)/2
  !               if(ABS(height(i,j,k,index))<MAX_EXT_H) then 
  !                  ! height has been properly computed in passes 1 and 2 and we are not more than max_ext_h from stack center. 
  !                  climitp2 = coordlimit(d,sign) + 2*sign
  !                  c(1) = i; c(2) = j; c(3) = k
  !                  c0 = c(d)
  !                  c(d) = c0 + sign ! start of region to be filled
  !                  limit_not_found=.not.(c0==climitp2) 
  !                  do while (limit_not_found) 
  !                     limit_not_found = .not.((vof_flag(c(1),c(2),c(3))==2) &
  !                          .or.(c(d)==climitp2).or.(abs(c(d)-c0)>=MAX_EXT_H))
  !                     height(c(1),c(2),c(3),index) = height(i,j,k,index) + c0 - c(d)
  !                     c(d) = c(d) + sign 
  !                  enddo
  !               endif
  !            enddo!; enddo
  !         endif
  !      enddo; enddo; enddo
  !    end subroutine get_heights_pass3
  !=======================================================================================================
  !   Check if we find nine heights in the neighboring cells, if not collect all heights in all directions
  !=======================================================================================================
  subroutine get_local_heights(i1,j1,k1,mxyz,try,nfound,hloc,points,nposit)
    implicit none
    integer, intent(in) :: i1(-1:1,-1:1,3), j1(-1:1,-1:1,3), k1(-1:1,-1:1,3)  
    ! i1(:,:,d) 3x3 plane rotated in direction d
    integer, intent(out) :: nfound
    real(8), intent(in)  :: mxyz(3)
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
    points = 0.d0
    hloc = 2d6
    nposit = 0 
    l=0
    dirnotfound=.true.
    deltax=dx(nx/2)
    do while (l.lt.3.and.dirnotfound)
       l = l+1
       d = try(l)    ! on entry, try(l) should sort the directions , closest to normal first. 
       if(d.eq.1) then
          si=1; sj=0; sk=0
       else if (d.eq.2) then
          si=0; sj=1; sk=0;
       else if (d.eq.3) then
          si=0; sj=0; sk=1
       else
          call pariserror("bad direction")
       endif
       index =  2*(d-1)+2
       if(mxyz(d).gt.0) index = 2*(d-1)+1
       hloc = 2d6
       nfound = 0
       do m=-1,1 
          do n=-1,1
             if(height(i1(m,n,d),j1(m,n,d),k1(m,n,d),index).lt.1d6) then  ! search at same level
                ! one height found
                hloc(m,n) = height(i1(m,n,d),j1(m,n,d),k1(m,n,d),index)
                nfound = nfound + 1
                nposit = nposit + 1
                points(nposit,1) = hloc(m,n)*si + i1(m,n,d)-i
                points(nposit,2) = hloc(m,n)*sj + j1(m,n,d)-j
                points(nposit,3) = hloc(m,n)*sk + k1(m,n,d)-k
             else
                s = 1 
                heightnotfound=.true.
                do while(s.le.Ng.and.heightnotfound) ! search at other levels
                   if (height(i1(m,n,d)+si*s,j1(m,n,d)+sj*s,k1(m,n,d)+sk*s,index).lt.1d6) then
                      hloc(m,n) = height(i1(m,n,d)+si*s,j1(m,n,d)+sj*s,k1(m,n,d)+sk*s,index) + s
                      nfound = nfound + 1
                      nposit = nposit + 1
                      points(nposit,1) = hloc(m,n)*si + i1(m,n,d)-i
                      points(nposit,2) = hloc(m,n)*sj + j1(m,n,d)-j
                      points(nposit,3) = hloc(m,n)*sk + k1(m,n,d)-k
                      heightnotfound=.false.  ! to exit loop
                   else if  (height(i1(m,n,d)-si*s,j1(m,n,d)-sj*s,k1(m,n,d)-sk*s,index).lt.1d6) then
                      hloc(m,n) = height(i1(m,n,d)-si*s,j1(m,n,d)-sj*s,k1(m,n,d)-sk*s,index) - s
                      nfound = nfound + 1
                      nposit = nposit + 1
                      points(nposit,1) = hloc(m,n)*si + i1(m,n,d)-i
                      points(nposit,2) = hloc(m,n)*sj + j1(m,n,d)-j
                      points(nposit,3) = hloc(m,n)*sk + k1(m,n,d)-k
                      heightnotfound=.false.  ! to exit loop
                   endif
                   s = s + 1
                end do ! while s lt ndepth 
             end if ! search at same level
          end do ! n
       end do ! m 
       if(nfound.eq.9) then
          dirnotfound = .false.
          ! on exit, redefine try() so that try(1) be the h direction found
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
          return
       end if ! nfound = 9
    end do ! d and dirnotfound
    if(nposit.gt.NPOS) call pariserror("GLH: nposit")
  end subroutine get_local_heights
  !
  !=================================================================================================
  !
  ! the core of HF computation
  !
  !=================================================================================================
  subroutine print_cvof_3x3x3(i0,j0,k0)
    implicit none
    integer, intent(in) :: i0,j0,k0
    integer :: l,m
    print *, "vof 3^3 cube at i j k ", i0,j0,k0
    print *, "                x y z ", x(i0), y(j0), z(k0)
    print *, "cvof()"
    do l=-1,1
       print *, " "
       do  m=-1,1
          print *, cvof(i0-1:i0+1,j0+m,k0+l)
       enddo
    enddo
    print *, " "
    print *, "flags"
    do l=-1,1
       print *, " "
       do  m=-1,1
          print *, vof_flag(i0-1:i0+1,j0+m,k0+l)
       enddo
    enddo
    print *, " "
  end subroutine print_cvof_3x3x3

  subroutine get_all_curvatures(kapparray,iout)
    implicit none
    real(8), intent(out) :: kapparray(imin:imax,jmin:jmax,kmin:kmax)
    integer, intent(in) :: iout
    call get_all_curvatures_pop(kapparray,iout)
  end subroutine get_all_curvatures

  !=================================================================================================
  !
  ! get_all_curvature_pop
  !
  !=================================================================================================

  ! This is a global function.
  ! It does a two-pass curvature computation. 
  ! In the first pass, curvature is computed from nine heights 
  ! or from mixed heights. 
  ! 
  ! In the second pass if curvature is not found it is 
  ! computed from averages if possible and otherwise from centroids. 
  !

  ! kappa_flag=0  9 heights
  ! kappa_flag=1  mixed heights
  ! kappa_flag=2  average from the above
  ! kappa_flag=3  centroids
  ! kappa_flag=4  no curvature possible  UNCOMPUTABLE
  ! kappa_flag=5  undecided yet   UNDECIDED

  !=================================================================================================

  subroutine get_all_curvatures_pop(kapparray,iout)
    use module_IO
    implicit none
    include 'mpif.h'
    real(8), intent(out) :: kapparray(imin:imax,jmin:jmax,kmin:kmax)
    real(8) :: afit(6), kappa
    integer :: ierr, i,j,k, nfound, nposit
    integer :: req(24),sta(MPI_STATUS_SIZE,24)
    integer, intent(in) :: iout
    integer, allocatable, dimension(:,:,:) :: kappaflag
    integer, parameter :: UNDECIDED = 5
    integer, parameter :: UNCOMPUTABLE = 4

    if(.not.st_initialized) call initialize_surface_tension()
    allocate(kappaflag(imin:imax,jmin:jmax,kmin:kmax))

    !*** Initialize
    kapparray=UNCOMPUTED
    kappaflag=UNDECIDED
    method_count=0
    do k=ks,ke; do j=js,je; do i=is,ie
       kappa = UNCOMPUTED
       if (vof_flag(i,j,k) == 2 ) then  ! mixed cell
          call curvature_pop_pass1(i,j,k,kappa,nfound)
          kapparray(i,j,k) = kappa
          if(nfound==0) then
             kappaflag(i,j,k) = 0 
             method_count(1) = method_count(1) + 1  ! nine heights
          endif
          if(nfound>-50.and.nfound<0) then
             kappaflag(i,j,k) = 1  !* could be faster
             method_count(2) = method_count(2) + 1  ! mixed heights
          endif
       endif
    enddo;enddo;enddo

    call do_all_ghost(kapparray)
    call do_all_ighost(kappaflag) 

    ! *** Pass 2
    do k=ks,ke; do j=js,je; do i=is,ie
       if (vof_flag(i,j,k) == 2 ) then
          if (kappaflag(i,j,k) > 1) then  ! failed curvature computation, 
             ! try averaging 
             call get_average(kapparray,kappaflag,i,j,k)
             if(kappaflag(i,j,k)==2)  method_count(3) = method_count(3) + 1  ! average
          endif
       endif
    enddo;enddo;enddo

    do k=ks,ke; do j=js,je; do i=is,ie
       kappa=UNCOMPUTED
       if (vof_flag(i,j,k) == 2 ) then
          if (kappaflag(i,j,k) > 2) then 
             call curvature_pop_pass2(i,j,k,kappa,nfound)
             if(nfound>0) then
                kappaflag(i,j,k) = UNCOMPUTABLE ! could not find centroid curvature
                method_count(5)  = method_count(5) + 1 
             else
                kapparray(i,j,k) = kappa
                kappaflag(i,j,k) = 3 
                method_count(4)  = method_count(4) + 1  ! centroids
             endif
          endif
       else ! all pure cells even if close to the phase boundary are discarded (may have a pure face)
          kappaflag(i,j,k) = UNCOMPUTABLE
       endif
    enddo;enddo;enddo
!     if(debug_curvature) print *, "mcount",  method_count

    if(iout==-1) then  ! test run
       call print_method()
       deallocate(kappaflag)
       return
    endif

    call do_all_ghost(kapparray)
    call do_all_ighost(kappaflag) 

    ! Clip and check curvature
    do k=kmin,kmax; do j=jmin,jmax; do i=imin,imax
       if(kappaflag(i,j,k) < UNCOMPUTABLE) then
          kappa=kapparray(i,j,k) 
          if(abs(kappa)>kappamax) then
             geom_case_count(17) = geom_case_count(17) + 1
             kappa = sign(1d0,kappa)*kappamax
          endif
          kapparray(i,j,k) = kappa  
          if(kappaflag(i,j,k) == UNDECIDED) then 
             call pariserror("Undecided curvature")
          endif
       endif
       if(debug_curvature) then !debugging
          if (kappa /= kappa) then 
             write(*,'("Kappa read into array is NaN, Kapparay, Kappa: ",2e14.5,3I8)') &
                  kapparray(i,j,k), kappa, i,j,k 
             call pariserror("Kappa read into array is NaN")
          endif
        endif                    ! end debugging
    enddo;enddo;enddo
    if (debug_par) then
       call write_par_var("Kappa     ",iout,kapparray)
    endif
    deallocate(kappaflag)
  end subroutine get_all_curvatures_pop

  ! Compute averages in a rather slow way (no vectorization possible)
  subroutine get_average(kapparray,kappaflag,i0,j0,k0)
    integer, intent(in) :: i0,j0,k0
    real(8), intent(inout) :: kapparray(imin:imax,jmin:jmax,kmin:kmax)
    integer, intent(inout) :: kappaflag(imin:imax,jmin:jmax,kmin:kmax)
    integer :: i,j,k,m,n,l,s
    real(8) :: kappa
    s=0; kappa=0d0
    do m=-1,1; do n=-1,1; do l=-1,1
       i=i0+m
       j=j0+n
       k=k0+l
       if(kappaflag(i,j,k) < 2) then
          s=s+1
          kappa=kappa+kapparray(i,j,k)
       endif
    enddo;enddo;enddo
    if(s>0) then
       kapparray(i0,j0,k0) = kappa/s
       kappaflag(i0,j0,k0) = 2
    endif
  end subroutine get_average
  !
  subroutine print_method()
    include "mpif.h"
    integer :: total
    integer :: n,ierr
    real(8) :: fraction(5)
    integer :: glob_m_count(5)
    call MPI_REDUCE(method_count, glob_m_count, 5, MPI_INTEGER , MPI_SUM, 0, MPI_COMM_Domain, ierr)
    if(rank==0) then
       total=0
       do n=1,5
          total = glob_m_count(n) + total
       enddo
       do n=1,5
          fraction(n) = float(glob_m_count(n)) / float(total)   
       enddo

       OPEN(UNIT=89,FILE='mcount.tmp')
       write(89,*) fraction
       close(89)
    endif
  end subroutine print_method
  !
  !=================================================================================================
  !
  ! curvature_pop_pass1(i0,j0,k0,kappa,nfound,nposit,a,pure_non_bulk)
  !
  !=================================================================================================
  
  ! This is a local function for a given cell. 
  ! It does the first pass of a two-pass curvature computation. 
  !
  ! A) The curvature is  returned by the variable kappa,
  ! except when the algorithm fails. If that is the 
  ! case, then kappa is unchanged from it initialisation
  ! value. 
  
  ! B) The condition of the calculation is returned by the variable nfound. 
  
  ! Success: nfound <= 0
  ! nfound=0         : nine heights found
  ! -50 < nfound < 0 : -nfound mixed heights were used. 
  
  ! Failure: nfound > 0
  ! 0 < nfound < nfound_min : not enough mixed heights
  ! nfound=16        : no fit success 
  ! nfound=21        : failure for cells with exactly grad C = 0 
  !=================================================================================================
  subroutine curvature_pop_pass1(i0,j0,k0,kappa,nfound)
    implicit none
    integer, intent(in) :: i0,j0,k0
    real(8), intent(out) :: kappa
    integer, intent(out) :: nfound

    real(8) :: h(-1:1,-1:1),a(6)  
    integer :: m,n,l,i,j,k
    logical :: fit_success 
    integer :: i1(-1:1,-1:1,3), j1(-1:1,-1:1,3), k1(-1:1,-1:1,3),try(3)
    integer :: s,c(3),d,neighbor,esign,nposit

    real(8) :: points(NPOS,3),bpoints(NPOS,3),origin(3)
    real(8) :: fit(NPOS,3),weights(NPOS)
    real(8) :: mxyz(3),mv(3),stencil3x3(-1:1,-1:1,-1:1)
    real(8) :: wg, kappasign, area3d

    call map3x3in2x2(i1,j1,k1,i0,j0,k0)
    !   define in which order directions will be tried 
    !   direction closest to normal first
    !   a) first determine normal
    if(recomputenormals) then
       do m=-1,1; do n=-1,1; do l=-1,1
          stencil3x3(m,n,l) = cvof(i0+m,j0+n,k0+l)
       enddo;enddo;enddo
       call mycs(stencil3x3,mxyz)
    else
       mxyz(1) = n1(i0,j0,k0)      
       mxyz(2) = n2(i0,j0,k0)      
       mxyz(3) = n3(i0,j0,k0)
    endif
    !   b) order the orientations
    call orientation(mxyz,try)
    call get_local_heights(i1,j1,k1,mxyz,try,nfound,h,points,nposit)
    ! Begin computation. If it does not succeed, return with kappa unchanged.
    ! Avoid finding curvature for cells with exactly grad C = 0 
    ! Note: a perfect plane will have 9 heights and grad C = 0
    if ( ABS(mxyz(1)) < TINY .and. ABS(mxyz(2)) < TINY .and. ABS(mxyz(3)) < TINY) then
       geom_case_count(21) = geom_case_count(21) + 1
       nfound = 21
       return 
    end if
    ! ***
    ! first prize: all nine heights found 
    ! ***
    if ( nfound == 9 .and.use_full_heights) then
       !
       !  h = a6  + a4 x + a5 y + a3 xy + a1 x**2 + a2 y**2
       !
       a(1) = h(1,0)-2.d0*h(0,0)+h(-1,0)
       a(2) = h(0,1)-2.d0*h(0,0)+h(0,-1)
       a(3) = (h(1,1)-h(-1,1)-h(1,-1)+h(-1,-1))/4.d0
       a(4) = (h(1,0)-h(-1,0))/2.d0
       a(5) = (h(0,1)-h(0,-1))/2.d0
       kappa = (a(1)*(1.d0+a(5)*a(5)) + a(2)*(1.d0+a(4)*a(4)) - 2d0*a(3)*a(4)*a(5)) &
            /(1.d0+a(4)*a(4)+a(5)*a(5))**(1.5d0)
       kappa = sign(1.d0,mxyz(try(1)))*kappa
       nfound=0
       return
    endif ! endif for first prize. 

    ! if not succesful, continue search
    ! ***
    ! Second prize: determine curvature from mixed heights and fits
    ! ***
    ! 1) determine the origin. 
    call FindCutAreaCentroid(i0,j0,k0,origin)
    ! 2) search for independent positions from mixed heights. Include origin.
    if(nposit>NPOS-1) call pariserror("nposit too large") ! paranoid programming
    nfound = ind_pos_sorted_add(origin,points,bpoints,nposit) 
    ! 3) Rotate coordinate sytems by permutation of x,y,z
    !  x_i' = x_k(i)
    !  m'_i = m_k(i)
    mv(1) = mxyz(try(2))
    mv(2) = mxyz(try(3))
    mv(3) = mxyz(try(1))
    if ( nfound > nfound_min + 1 )  then  ! At least 7 points + origin centroid
                                       !  to avoid special 2D degeneracy. 
       ! rotate and shift origin
       !  x_i' = x_k(i)
       !  m'_i = m_k(i)
       points(:,1) = bpoints(:,try(2))  - origin(try(2))
       points(:,2) = bpoints(:,try(3))  - origin(try(3))
       points(:,3) = bpoints(:,try(1))  - origin(try(1))   
       weights=1d0
       weights(1) = 1d2*area3d(mxyz,cvof(i0,j0,k0))
       ! 4) fit over all positions returned by ind_pos 
       ! call parabola_fit_with_rotation(points,fit,weights,mv,nposit,a,kappasign,fit_success) 
       ! 4) fit only over positions separated by a minimum distance
       call parabola_fit_with_rotation(points,fit,weights,mv,nfound,a,kappasign,fit_success) 
       ! call parabola_fit(points,weights,nposit,a,fit_success) 
       if(fit_success) then
          kappa = 2.d0*(a(1)*(1.d0+a(5)*a(5)) + a(2)*(1.d0+a(4)*a(4)) - a(3)*a(4)*a(5)) &
               /(1.d0+a(4)*a(4)+a(5)*a(5))**(1.5d0)
          kappa = kappasign*kappa
          nfound = - nfound  ! encode the fact that mixed-heights were used 
          return
       else
          geom_case_count(16) = geom_case_count(16) + 1
          nfound = 16  ! encode the fact that no curvature set and no fit success. kappa is unchanged from initialisation value.
          return
       endif ! fit_success
    endif !  nfound > nfound_min  
    ! if not succesful, return with positive nfound
    nfound = nfound + 1
  end subroutine curvature_pop_pass1
  !=================================================================================================
  !
  ! curvature_pop_pass2(kapparray1,i0,j0,k0,kappa,nfound,nposit,a)
  !
  !=================================================================================================
  
  ! This is a local function for a given cell. It attempts to compute the curvature from
  ! a paraboloid fit of VOF-PLIC facet centroid positions. 
  
  ! On exit:
  
  ! A) The curvature is  returned by the variable kappa,
  ! except when the algorithm fails. If that is the 
  ! case, then kappa is unchanged from it initialisation
  ! value. 
  
  ! B) The condition of the calculation is returned by the variable nfound. 
  
  ! Success: 
  ! nfound < -50     : -nfound - 50 centroids were used. 
  
  ! Failure: nfound > 0
  ! 11 < nfound < 15 : the paraboloid fit in 3D cannot be made due to insufficient number of centroids.
  ! nfound=16        : no fit success with centroids
  ! nfound=21        : failure for cells with exactly grad C = 0 
  !=================================================================================================
  subroutine curvature_pop_pass2(i0,j0,k0,kappa,nfound)
    implicit none
    integer, intent(in) :: i0,j0,k0
    real(8), intent(out) :: kappa
    integer, intent(out) :: nfound

    integer :: m,n,l,i,j,k
    logical :: fit_success 
    integer :: s,c(3),d,neighbor,esign, nposit,try(3)
    
    real(8) :: points(NPOS,3),bpoints(NPOS,3),origin(3),a(6)
    real(8) :: fit(NPOS,3),weights(NPOS)
    real(8) :: centroid(3),mxyz(3),mv(3),stencil3x3(-1:1,-1:1,-1:1)
    real(8) :: wg, kappasign
    
    if(recomputenormals) then
       do m=-1,1; do n=-1,1; do l=-1,1
          stencil3x3(m,n,l) = cvof(i0+m,j0+n,k0+l)
       enddo;enddo;enddo
       call fd32(stencil3x3,mxyz)
    else
       mxyz(1) = n1(i0,j0,k0)      
       mxyz(2) = n2(i0,j0,k0)      
       mxyz(3) = n3(i0,j0,k0)
    endif
    !   b) order the orientations
    call orientation(mxyz,try)
    ! Begin computation. If it does not succeed, return with kappa unchanged.
    ! Avoid finding curvature for cells with exactly grad C = 0 
    if ( ABS(mxyz(1)) < TINY .and. ABS(mxyz(2)) < TINY .and. ABS(mxyz(3)) < TINY) then
       geom_case_count(21) = geom_case_count(21) + 1
       nfound = 21
       return 
    end if
    ! 2) determine the origin. 
    call FindCutAreaCentroid(i0,j0,k0,origin)
    ! 3) Rotate coordinate sytems by permutation of x,y,z
    !  x_i' = x_k(i)
    !  m'_i = m_k(i)
    mv(1) = mxyz(try(2))
    mv(2) = mxyz(try(3))
    mv(3) = mxyz(try(1))
    ! ***
    ! Determine curvature from centroids
    ! ***
    ! 1) Find all centroids in 3**3
    ! filter out small cells if threshold set. 
    nposit=0
    do m=-1,1; do n=-1,1; do l=-1,1
       i=i0+m
       j=j0+n
       k=k0+l
       c(1)=m
       c(2)=n
       c(3)=l
       wg = cvof(i,j,k)*(1d0-cvof(i,j,k))
       if(vof_flag(i,j,k) == 2.and.wg>cwg_threshold) then
          nposit = nposit + 1
          call NewFindCutAreaCentroid(i,j,k,centroid)
          do s=1,3 
             fit(nposit,s) = centroid(s) + c(s)
          end do
          if(wg.lt.0d0) call pariserror("w<0")
          weights(nposit) = 1d0 ! wg  ! sqrt(wg)
       endif ! vof_flag
    enddo; enddo; enddo ! do m,n,l
    ! permute coordinates so z direction is closest to normal
    ! try(:) array contains direction closest to normal first
    points(:,1) = fit(:,try(2)) - origin(try(2))
    points(:,2) = fit(:,try(3)) - origin(try(3))
    points(:,3) = fit(:,try(1)) - origin(try(1))
    if(nposit.gt.NPOS) call pariserror("GLH: nposit")
    if(nposit < 6) then
       geom_case_count(nposit+10) = geom_case_count(nposit+10) + 1
       nfound = nposit+10  ! encode the fact that curvature was not set due to insufficient number of centroids. 
       return
    else
       call parabola_fit_with_rotation(points,fit,weights,mv,nposit,a,kappasign,fit_success) 
       if(.not.fit_success) then
          geom_case_count(16) = geom_case_count(16) + 1
          nfound = 16  ! encode the fact that no curvature set and no fit success. Curvature remains at
          ! initial (subroutine argument) value
          return
       else
          kappa = 2.d0*(a(1)*(1.d0+a(5)*a(5)) + a(2)*(1.d0+a(4)*a(4)) - a(3)*a(4)*a(5)) &
               /sqrt(1.d0+a(4)*a(4)+a(5)*a(5))**3
!           print *,"I,j,k,kappa,kappasign,nposit ", i0,j0,k0,kappa,kappasign,nposit
          kappa = kappasign*kappa
          nfound = - nposit - 50  ! encode the fact that nposit distinct centroids were used 
          return
       endif
    endif
    call pariserror("GLH: this statement should be unreachable")
  end subroutine curvature_pop_pass2

    ! Performs a rotation to align z axis with normal, then calls fit

    subroutine parabola_fit_with_rotation(bfit,fit,weights,mv,nposit,a,kappasign,fit_success)
      implicit none
      real(8), intent(in)  :: weights(NPOS)
      integer, intent(in)  :: nposit
      real(8), intent(inout) :: mv(3),bfit(NPOS,3)
      real(8), intent(out) :: a(6),fit(NPOS,3),kappasign
      logical, intent(out) :: fit_success
      logical :: inv_success
      real(8) :: invm(3,3), rhs(3), norm(3), mv1(3)
      real(8) :: ev(3,3)   ! New basis vector expressed in old (canonical) basis: 
      ! ev(coord i, basis vector j) = ev(i,j) = e_ij
      ! x_old_i = e_ij x_new_j
      ! x_new_i = e_ji x_old_j
      real(8) :: testm(3,3), id(3,3), error
      integer :: i,j

      if(do_rotation) then
         ! normal = direction z
         ! e_z' = m
         ev(:,3) = mv
         ! now the x'_3 = z' direction is aligned with the normal so 
         kappasign = 1d0 
         ! mv(3) is the vector of components of the normal vector
         ! let mv1 == m1  be the canonical basis vector furthest from normal 
         mv1 = 0d0    
         mv1(2) = 1d0
         ! direction 1 orthogonal to normal and mv1
         ! e_x' = m1 x m
         ! full expression:
         ! ev(1,1) =  mv1(2)*mv(3) - mv1(3)*mv(2)
         ! ev(2,1) = -mv1(1)*mv(3) + mv1(3)*mv(1)
         ! ev(3,1) =  mv1(1)*mv(2) - mv1(2)*mv(1)
         ev(1,1) =   mv(3)
         ev(2,1) =   0d0
         ev(3,1) =  -mv(1)
         ! e_y' = - e_x' x e_z' 
         ! full expression:
         ! ev(1,2) =  - ev(2,1)*ev(3,3) + ev(3,1)*ev(2,3) =   - m1 m2
         ! ev(2,2) =    ev(1,1)*ev(3,3) - ev(3,1)*ev(1,3) =  m3^2 + m1^2
         ! ev(3,2) =  - ev(1,1)*ev(2,3) + ev(2,1)*ev(1,3) =   - m3 m2 

         ev(1,2) =  - ev(2,1)*ev(3,3) + ev(3,1)*ev(2,3)
         ev(2,2) =    ev(1,1)*ev(3,3) - ev(3,1)*ev(1,3)
         ev(3,2) =  - ev(1,1)*ev(2,3) + ev(2,1)*ev(1,3)
         norm = sqrt(ev(1,:)**2 + ev(2,:)**2 + ev(3,:)**2)
         do i=1,3
            ev(i,:) = ev(i,:)/norm
         enddo
         if(debug_curvature) then
            if(min(norm(1),norm(2),norm(3)) < 1d-12) call pariserror("small rotation matrix in curvature")
            invm = transpose(ev)
            testm = matmul(invm,ev)
            id=0d0; do i=1,3; id(i,i) = 1d0; enddo  ! define identity matrix
            invm = testm - id
            error =  0d0
            do i=1,3
               do j=1,3
                  error = error + abs(invm(i,j))
               enddo
            enddo
            if(error>1d-14) then
               call pariserror("non orthogonal rotation matrix")
            end if
         endif
         fit = bfit
         do i=1,3
            bfit(:,i) = ev(1,i)*fit(:,1) + ev(2,i)*fit(:,2) + ev(3,i)*fit(:,3)
         enddo
      else
         !      if no rotation then usual sign calculation
         kappasign = sign(1.d0,mv(3))
      endif ! do_rotation
      call parabola_fit(bfit,weights,nposit,a,fit_success)
    end subroutine parabola_fit_with_rotation
    !-------------------------------------------------------------------------------------------------------
    !
    !  Parabolic fit
    !
    !-------------------------------------------------------------------------------------------------------

    subroutine parabola_fit(fit,weights,nposit,a,fit_success)
      implicit none
      real(8), intent(in)  :: fit(NPOS,3)
      real(8), intent(in)  :: weights(NPOS)
      real(8), intent(out) :: a(6)
      logical, intent(out) :: fit_success
      real(8) :: m(6,6), invm(6,6)
      real(8) :: rhs(6)
      integer :: ifit, im,jm, nposit
      logical :: inv_success
      real(8) :: x1,x2,x3,x4,y1,y2,y3,y4,wg

      fit_success=.false.
      a = 0.d0
      ! evaluate the linear system for least-square fit
      m   = 0.d0
      rhs = 0.d0

      do ifit = 1, nposit
         x1 =    fit(ifit,1)
         x2 = x1*fit(ifit,1)
         x3 = x2*fit(ifit,1)
         x4 = x3*fit(ifit,1)
         y1 =    fit(ifit,2)
         y2 = y1*fit(ifit,2)
         y3 = y2*fit(ifit,2)
         y4 = y3*fit(ifit,2)
         wg = weights(ifit)

         ! The matrix is m_ij = sum_n alpha^n_i alpha^n_j
         ! and the "alpha_i" are the factors of the a_i coefficients:
         ! 
         !   x^2, y^2, xy, x^2, y^2, 1

         m(1,1) = m(1,1) + x4*wg
         m(2,2) = m(2,2) + y4*wg
         m(3,3) = m(3,3) + x2*y2*wg
         m(4,4) = m(4,4) + x2*wg
         m(5,5) = m(5,5) + y2*wg
         m(6,6) = m(6,6) + wg

         m(1,3) = m(1,3) + x3*y1*wg
         m(1,4) = m(1,4) + x3*wg
         m(1,5) = m(1,5) + x2*y1*wg
         m(2,3) = m(2,3) + x1*y3*wg
         m(2,4) = m(2,4) + x1*y2*wg
         m(2,5) = m(2,5) + y3*wg
         m(3,6) = m(3,6) + x1*y1*wg
         m(4,6) = m(4,6) + x1*wg
         m(5,6) = m(5,6) + y1*wg

         rhs(1) = rhs(1) + x2   *fit(ifit,3)*wg
         rhs(2) = rhs(2) + y2   *fit(ifit,3)*wg
         rhs(3) = rhs(3) + x1*y1*fit(ifit,3)*wg
         rhs(4) = rhs(4) + x1   *fit(ifit,3)*wg
         rhs(5) = rhs(5) + y1   *fit(ifit,3)*wg
         rhs(6) = rhs(6) +       fit(ifit,3)*wg
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
         !         print *, "WARNING: no fit success"
      end if ! inv_success
    end subroutine parabola_fit
    !-------------------------------------------------------------------------------------------------------
    !
    ! Matrix inversion 
    !
    !  Subroutine to find the inverse of a square matrix with non-zero diagonal coeeficients. 
    !  From Stanley's previous code
    !
    !-------------------------------------------------------------------------------------------------------
    SUBROUTINE FindInverseMatrix(matrix,inverse,n,inverse_success)
      implicit none
      integer, intent(in ) :: n
      real(8), intent(in ), dimension(n,n) :: matrix 
      real(8), intent(out), dimension(n,n) :: inverse
      logical, intent(out) :: inverse_success 
      integer :: i, j, k, l
      real(8) :: m
      real(8), dimension(n,2*n) :: augmatrix ! augmented matrix
      real(8), dimension(:,:), allocatable :: test

      ! augment input matrix with an identity matrix
      do i = 1,n
         do j = 1,2*n
            if (j <= n ) then
               augmatrix(i,j) = matrix(i,j)
            else if ((i+n) == j) then
               augmatrix(i,j) = 1.0d0
            else
               augmatrix(i,j) = 0.0d0
            endif
         end do
      end do
      ! ensure initial diagonal elements are non-zero
      do k = 1,n-1
         if (abs(augmatrix(k,k)) <  EPSC) then
            do i = k+1, n
               if (abs(augmatrix(i,k)) > EPSC ) then
                  do  l = 1, 2*n
                     augmatrix(k,l) = augmatrix(k,l)+augmatrix(i,l)
                  end do
               endif
            end do
         endif
      end do

      ! reduce augmented matrix to upper triangular form
      do k =1, n-1
         if (abs(augmatrix(k,k)) <  EPSC ) then 
            inverse_success = .false.
            return
         end if
         do j = k+1, n
            m = augmatrix(j,k)/augmatrix(k,k)
            do i = k, 2*n
               augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
            end do
         end do
      end do

      ! test for invertibility
      if( abs(augmatrix(n,n)) <  EPSC) then
         ! write(*,*) "error-matrix has A_",n,n,"=",augmatrix(n,n)," and is non-invertible"
         inverse_success = .false.
         return
      endif

      ! rescale lines to make diagonal elements unity
      do i = 1 , n
         m = augmatrix(i,i)
         do j = i , 2*n
            augmatrix(i,j) = augmatrix(i,j)/m
         end do
      end do

      ! reduce left side of augmented matrix to identity matrix
      do k = n-1, 1, -1
         do i =1, k
            m = augmatrix(i,k+1)
            do j = k, (2*n)
               augmatrix(i,j) = augmatrix(i,j) - augmatrix(k+1,j)*m
            end do
         end do
      end do

      ! store answer
      inverse(:,1:n) = augmatrix(:,n+1:2*n)

      if(debug_curvature) then
         allocate(test(n,n))
         test = matmul(matrix,inverse) 
         do i=1,n
            test(i,i) = test(i,i) - 1d0
         enddo
         if (maxval(abs(test)).gt.1d-5) print *, "max error of matrix inversion > 1d-5:", maxval(abs(test))
         deallocate(test)
      endif
      inverse_success = .true.
    end subroutine FindInverseMatrix
    !-------------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------------
    subroutine NewFindCutAreaCentroid(i,j,k,centroid)
      implicit none

      integer, intent(in)  :: i,j,k
      real(8), intent(out) :: centroid(3)

      integer :: l,m,n
      real(8) :: nr(3)
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
      call area_centroid(nr,cvof(i,j,k),centroid)
      centroid = centroid - 0.5d0
    end subroutine NewFindCutAreaCentroid
    !-------------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------------
    subroutine FindCutAreaCentroid(i,j,k,centroid)
      implicit none

      integer, intent(in)  :: i,j,k
      real(8), intent(out) :: centroid(3)

      integer :: l,m,n
      real(8) :: dmx,dmy,dmz, mxyz(3),px,py,pz
      real(8) :: invx,invy,invz
      real(8) :: alpha, al3d, nr(3)
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
         nr = mxyz
      else
         nr(1) = n1(i,j,k)      
         nr(2) = n2(i,j,k)      
         nr(3) = n3(i,j,k)
      endif

      dmx = nr(1); dmy = nr(2); dmz = nr(3)
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
      alpha = al3d([dmx,dmy,dmz],cvof(i,j,k))
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
    !-------------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------------
    ! 
    !   Computes the centroid as in gerris
    ! 
    ! * Fills p with the position of the center of mass of the polygon
    ! * obtained by interseectiing the plane  (m,alpha).
    ! * with the reference cell.
    !
    !  assumptions: dmx,dmy,dmz > 0 and |dmx| + |dmy| + |dmz| = 1
    !
    subroutine PlaneAreaCenter (dmx,dmy,dmz, alpha, px,py,pz)
      implicit none
      real(8), intent(in) :: dmx,dmy,dmz,alpha
      real(8), intent(out) :: px,py,pz
      real(8) :: nx,ny,qx,qy
      real(8) :: area,b,amax

      if(dmx<0.d0.or.dmy<0.d0.or.dmz<0.d0) call pariserror("invalid dmx dmy dmz")
      if(abs(dmx+dmy+dmz-1d0)>EPS_GEOM) call pariserror("invalid dmx+dmy+dmz")
      if (dmx < EPS_GEOM) then
         nx = dmy
         ny = dmz
         call LineCenter (nx,ny, alpha, qx,qy)
         px = 0.5d0
         py = qx
         pz = qy
         return
      endif
      if (dmy < EPS_GEOM) then
         nx = dmz
         ny = dmx
         call LineCenter (nx,ny, alpha, qx,qy)
         px = qy
         py = 0.5d0
         pz = qx
         return
      endif
      if (dmz < EPS_GEOM) then
         call LineCenter (dmx,dmy, alpha, px,py)
         pz = 0.5
         return
      endif

      if (alpha < 0.d0 .or. alpha > 1.d0) then
         print *, "alpha =", alpha
         call pariserror("PAC: invalid alpha")
      endif

      ! Debris patch. Ideally, should not be looking for centroids of debris. 
      if(alpha < EPS_GEOM) then
         px=0d0; py=0d0; pz=0d0
         return
      endif

      if(alpha > 1d0 - EPS_GEOM) then
         px=1d0; py=1d0; pz=1d0
         return
      endif
      ! end debris patch

      area = alpha*alpha
      px = area*alpha
      py = area*alpha
      pz = area*alpha
      b = alpha - dmx
      if (b > 0.) then
         area = area - b*b
         px = px - b*b*(2.d0*dmx + alpha)
         py = py - b*b*b
         pz = pz - b*b*b
      endif
      b = alpha - dmy
      if (b > 0.d0) then
         area = area - b*b
         py = py - b*b*(2.d0*dmy + alpha)
         px = px - b*b*b
         pz = pz - b*b*b
      endif
      b = alpha - dmz
      if (b > 0.d0) then
         area = area - b*b
         pz = pz - b*b*(2.d0*dmz + alpha)
         px = px - b*b*b
         py = py - b*b*b
      endif

      amax = alpha - 1.d0
      b = amax + dmx
      if (b > 0.d0) then
         area = area + b*b
         py = py + b*b*(2.d0*dmy + alpha - dmz)
         pz = pz + b*b*(2.d0*dmz + alpha - dmy)
         px = px + b*b*b
      endif
      b = amax + dmy
      if (b > 0.d0) then
         area = area + b*b
         px = px + b*b*(2.d0*dmx + alpha - dmz)
         pz = pz + b*b*(2.d0*dmz + alpha - dmx)
         py = py + b*b*b
      endif
      b = amax + dmz
      if (b > 0.d0) then
         area = area + b*b
         px = px + b*b*(2.d0*dmx + alpha - dmy)
         py = py + b*b*(2.d0*dmy + alpha - dmx)
         pz = pz + b*b*b
      endif

      area  = 3.d0*area
      px = px/(area*dmx)
      py = py/(area*dmy)
      pz = pz/(area*dmz)

      call THRESHOLD (px)
      call THRESHOLD (py)
      call THRESHOLD (pz)

    end subroutine PlaneAreaCenter
    !-------------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------------
    subroutine LineCenter (dmx,dmy, alpha, px,py)
      implicit none
      real(8), intent(in) :: dmx,dmy,alpha
      real(8), intent(out) :: px,py

      !  if (alpha <= 0.d0 .or. alpha >= 1.d0) 
      ! call pariserror("LC: invalid alpha")

      if (alpha < 0.d0 .or. alpha > 1.d0) then
         print *, "alpha =", alpha
         call pariserror("LC: invalid alpha")
      endif

      if (dmx < EPS_GEOM) then
         px = 0.5d0
         py = alpha;
         return
      endif

      if (dmy < EPS_GEOM) then
         py = 0.5d0;
         px = alpha
         return
      endif

      px = 0.; py = 0.d0

      if (alpha >= dmx) then
         px = px +  1.d0
         py = py +  (alpha - dmx)/dmy
      else
         px = px +  alpha/dmx
      endif

      if (alpha >= dmy) then
         py = py +  1.d0
         px = px +  (alpha - dmy)/dmx
      else
         py = py +  alpha/dmy
      endif

      px = px/2.d0
      py = py/2.d0

      call THRESHOLD (px)
      call THRESHOLD (py)
      !    if(px.ne.px) call pariserror("LAC:invalid px")
      !    if(py.ne.py) call pariserror("LAC:invalid px")
    end subroutine LineCenter
    !------------------------------------------------------------------------

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
      real(8), intent(in) :: points(NPOS,3)
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
               depends = (d2 < 0.5d0**2)
            endif
         enddo
         if(.not.depends) ni = ni + 1
      enddo
      ind_pos = ni
    end function ind_pos

    function ind_pos_sorted (points, bpoints, n)
      implicit none
      integer :: ind_pos_sorted
      integer, intent(in) :: n
      real(8), intent(inout) :: points(NPOS,3), bpoints(NPOS,3)
      integer :: i,j,ni,c
      real(8) :: d2,d1
      logical :: reject
      if (n < 2) then
         ind_pos_sorted = n
         return
      endif
      ni=1
      bpoints(ni,:) = points(1,:)
      do j=2,n      ! j: new point
         d1 =  sum(points(j,:)**2)                  ! distance (O,x_j)
         ! reject = (d1>3d0)
         reject=.false.
         do i=1,j-1 ! i: old point
            if(.not.reject) then
               d2 = 0d0
               do c=1,3
                  d2 = d2 + (points(i,c) - points(j,c))**2  ! distance (x_i,x_j)
               enddo
               reject = (d2 < 0.5d0**2)
            endif
         enddo
         if(.not.reject)  then ! add new point to list
            ni = ni + 1
            bpoints(ni,:) = points(j,:)
         endif
      enddo
      ind_pos_sorted = ni
    end function ind_pos_sorted

    function ind_pos_sorted_add (centroid, points, bpoints, n)
      implicit none
      integer :: ind_pos_sorted_add
      integer, intent(in) :: n
      real(8), intent(in) :: points(NPOS,3), centroid(3)
      real(8), intent(out) :: bpoints(NPOS,3)
      integer :: i,j,ni,c
      real(8) :: d2,d1
      logical :: reject
      bpoints(1,:) = centroid
      if (n < 2) then
         ind_pos_sorted_add = n
         return
      endif
      ni=2
      bpoints(ni,:) = points(1,:)
      do j=2,n      ! j: new point
         d1 =  sum(points(j,:)**2)                  ! distance (O,x_j)
         ! reject = (d1>3d0)
         reject=.false.
         do i=1,j-1 ! i: old point
            if(.not.reject) then
               d2 = 0d0
               do c=1,3
                  d2 = d2 + (points(i,c) - points(j,c))**2  ! distance (x_i,x_j)
               enddo
               reject = (d2 < 0.5d0**2)
            endif
         enddo
         if(.not.reject)  then ! add new point to list
            ni = ni + 1
            bpoints(ni,:) = points(j,:)
         endif
      enddo
      ind_pos_sorted_add = ni
    end function ind_pos_sorted_add

    subroutine project_velocity_staggered (vel,velmask,dtp,pres,d)
      use module_grid
      use module_flow

      implicit none
      real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: vel
      real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: velmask,pres
      real(8), intent(in) :: dtp
      real(8) :: rhoavg,invdeltax
      integer :: i,j,k,d
      integer :: i0,j0,k0

      call init_i0j0k0 (d,i0,j0,k0)
      call get_staggered_fractions (tmp,d) 

      if (d.eq.1) then
         invdeltax = 1.d0/dx(nx)
      elseif (d.eq.2) then
         invdeltax = 1.d0/dy(ny)
      else
         invdeltax = 1.d0/dz(nz)
      endif

      if (STGhost) tmp = NINT(tmp)

      do k=ks,ke;  do j=js,je; do i=is,ieu  
         rhoavg = tmp(i,j,k)*rho2 + (1.d0-tmp(i,j,k))*rho1
         vel(i,j,k)=vel(i,j,k)-dtp*velmask(i,j,k)*invdeltax*(pres(i+i0,j+j0,k+k0)-pres(i,j,k))/rhoavg
      enddo; enddo; enddo

      if (STGhost) then 

         call get_all_curvatures(work(:,:,:,1),0)

         do k=ks,ke; do j=js,je; do i=is,ie

            rhoavg = tmp(i,j,k)*rho2 + (1.d0-tmp(i,j,k))*rho1
            vel(i,j,k) = vel(i,j,k) - dtp/rhoavg*sigma* &
                 (1.d0-tmp(i,j,k))*(work(i+i0,j+j0,k+k0,1)*cvof(i+i0,j+j0,k+k0)-&
                 work(i,j,k,1)*cvof(i,j,k))*invdeltax**2

            vel(i,j,k) = vel(i,j,k) + dtp/rhoavg*sigma* &
                 tmp(i,j,k)*(work(i+i0,j+j0,k+k0,1)*(1.d0-cvof(i+i0,j+j0,k+k0)) &
                 -work(i,j,k,1)*(1.d0-cvof(i,j,k)))*invdeltax**2

         enddo; enddo; enddo
      endif

    end subroutine project_velocity_staggered

    !=================================================================================================
    ! The Poisson equation for the pressure is setup with matrix A as
    ! A7*Pijk = A1*Pi-1jk + A2*Pi+1jk + A3*Pij-1k + 
    !           A4*Pij+1k + A5*Pijk-1 + A6*Pijk+1 + A8
    !-------------------------------------------------------------------------------------------------
    subroutine get_Poisson_matrix_x (A1,A2,A8, velmask,dtp,d)

      use module_grid
      use module_BC
      use module_2phase
      use module_flow

      implicit none
      real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: velmask
      real(8), dimension(is:ie,js:je,ks:ke), intent(out) :: A1,A2
      real(8), dimension(is:ie,js:je,ks:ke), intent(inout) :: A8
      real(8), intent(in) :: dtp
      real(8) :: rhoavg, invdeltax
      integer :: i,j,k,d
      integer :: i0,j0,k0

      call init_i0j0k0 (d,i0,j0,k0)
      call get_staggered_fractions (tmp,d) 
      if (d.eq.1) then
         invdeltax = 1.d0/dx(nx)
      elseif (d.eq.2) then
         invdeltax = 1.d0/dy(ny)
      else
         invdeltax = 1.d0/dz(nz)
      endif

      ! I use STGhost flag because it allows to call this function without
      ! the ghost fluid method if required in the future
      if (STGhost) tmp = NINT(tmp)

      do k=ks,ke; do j=js,je; do i=is,ie
         rhoavg = tmp(i-i0,j-j0,k-k0)*rho2 + (1.d0-tmp(i-i0,j-j0,k-k0))*rho1
         A1(i,j,k) = dtp*velmask(i-i0,j-j0,k-k0)*invdeltax**2/rhoavg
         rhoavg = tmp(i,j,k)*rho2 + (1.d0-tmp(i,j,k))*rho1
         A2(i,j,k) = dtp*velmask(i,j,k)*invdeltax**2/rhoavg
      enddo; enddo; enddo

      if (STGhost) then 

         call get_all_curvatures(work(:,:,:,1),0)

         do k=ks,ke; do j=js,je; do i=is,ie

            rhoavg    = tmp(i-i0,j-j0,k-k0)*rho2 + (1.d0-tmp(i-i0,j-j0,k-k0))*rho1
            A8(i,j,k) = A8(i,j,k) + dtp/rhoavg*sigma*invdeltax**3*( &
                 (1.d0-tmp(i-i0,j-j0,k-k0))*work(i-i0,j-j0,k-k0,1)*cvof(i-i0,j-j0,k-k0) &
                 -(1.d0-tmp(i-i0,j-j0,k-k0))*work(i,j,k,1)*cvof(i,j,k) &
                 -tmp(i-i0,j-j0,k-k0) *work(i-i0,j-j0,k-k0,1)*(1.d0-cvof(i-i0,j-j0,k-k0)) &
                 +tmp(i-i0,j-j0,k-k0) *work(i,j,k,1)         *(1.d0-cvof(i,j,k))  )

            rhoavg = tmp(i,j,k)*rho2 + (1.d0-tmp(i,j,k))*rho1
            A8(i,j,k) = A8(i,j,k) + dtp/rhoavg*sigma*invdeltax**3*( &
                 (1.d0-tmp(i,j,k))*work(i+i0,j+j0,k+k0,1)*cvof(i+i0,j+j0,k+k0) &
                 - (1.d0-tmp(i,j,k))*work(i,j,k,1)*cvof(i,j,k) &
                 - tmp(i,j,k)*work(i+i0,j+j0,k+k0,1)*(1.d0-cvof(i+i0,j+j0,k+k0)) &
                 + tmp(i,j,k)*work(i,j,k,1)         *(1.d0-cvof(i,j,k))    ) 

         enddo; enddo; enddo

      endif

    end subroutine get_Poisson_matrix_x

    subroutine SetupPoissonGhost(utmp,vtmp,wtmp,umask,vmask,wmask,dt,A,pmask,VolumeSource) 
      use module_grid
      use module_BC
      use module_2phase
      use module_IO

      implicit none
      real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: utmp,vtmp,wtmp
      real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: umask,vmask,wmask
      real(8), dimension(is:ie,js:je,ks:ke), intent(inout) :: pmask
      real(8), dimension(is:ie,js:je,ks:ke,8), intent(out) :: A
      real(8), intent(in) :: dt, VolumeSource
      integer :: i,j,k,l

      A(:,:,:,8) = 0.d0

      call get_Poisson_matrix_x(A(:,:,:,1),A(:,:,:,2),A(:,:,:,8),umask,dt,1)
      call get_Poisson_matrix_x(A(:,:,:,3),A(:,:,:,4),A(:,:,:,8),vmask,dt,2)
      call get_Poisson_matrix_x(A(:,:,:,5),A(:,:,:,6),A(:,:,:,8),wmask,dt,3)

      do k=ks,ke; do j=js,je; do i=is,ie
         A(i,j,k,7) = sum(A(i,j,k,1:6))
         A(i,j,k,8) =  A(i,j,k,8) - (  VolumeSource +(utmp(i,j,k)-utmp(i-1,j,k))/dx(i) &
              +  (vtmp(i,j,k)-vtmp(i,j-1,k))/dy(j) &
              +  (wtmp(i,j,k)-wtmp(i,j,k-1))/dz(k) )
      enddo; enddo; enddo

      call Poisson_BCs(A)
      call check_and_debug_Poisson(A,umask,vmask,wmask,1.,pmask,dt,VolumeSource)

    end subroutine SetupPoissonGhost

    subroutine output_centroids(nf)
      use module_grid
      implicit none
      real(8) :: cent(3)
      integer :: i,j,k,nf
      character(len=30) :: rootname

      rootname=TRIM(out_path)//'/centroid-'//TRIM(int2text(nf,padding))//'-'
      if (rank==0) call append_3d_file(TRIM(rootname))

      OPEN(UNIT=111,FILE=TRIM(rootname)//TRIM(int2text(rank,padding))//'.3D')
      write(111,'("X Y Z Centroid")')
      do i=is,ie; do j=js,je; do k=ks,ke
         if (vof_flag(i,j,k)==2) then
            call NewFindCutAreaCentroid(i,j,k,cent)
            write(111,'(4es17.8e3)')x(i)+cent(1)*dx(i),y(j)+cent(2)*dy(j),z(k)+cent(3)*dz(k),1.000
         endif
      enddo; enddo; enddo
      CLOSE(111)
    end subroutine output_centroids
    !=================================================================================================

  end module module_surface_tension



