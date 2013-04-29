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
! module_solid: Contains definition of variables for the solids.
!-------------------------------------------------------------------------------------------------
module module_solid
  use module_grid
  use module_flow
  use module_IO
  use module_BC
  implicit none
  real(8), dimension(:,:,:), allocatable :: solids
  character(20) :: solid_type
  integer il,ih,jl,jh,kl,kh
  integer :: solid_opened=0 
!***********************************************************************
contains
!***********************************************************************
!=================================================================================================
  SUBROUTINE append_solid_visit_file(rootname)
    implicit none
    character(*) :: rootname
    integer prank
    if(rank.ne.0) stop 'rank.ne.0 in append_solid'

    if(solid_opened==0) then
       OPEN(UNIT=89,FILE='solid.visit')
       write(89,10) nPdomain
10     format('!NBLOCKS ',I4)
       solid_opened=1
    endif

    do prank=0,NpDomain-1
       write(89,11) rootname//TRIM(int2text(prank,padding))//'.vtk'
 11 format(A)
    enddo
  end subroutine  append_solid_visit_file

  subroutine close_solid_visit_file()
    close(89)
  end subroutine close_solid_visit_file
!***********************************************************************
!=================================================================================================
  SUBROUTINE initialize_solids()
    implicit none
    include 'mpif.h'
    integer :: i,j,k
    real(8) :: s1
    integer :: ierr
    integer :: req(48),sta(MPI_STATUS_SIZE,48)
    
    call ReadSolidParameters
    if(dosolids) then
       allocate(solids(imin:imax,jmin:jmax,kmin:kmax))
       solids=0.
!       if(rank==0) print*, "***x ",x 
       do i=imin,imax; do j=jmin,jmax; do k=kmin,kmax; 
          if(solid_type == 'CFC') then
             s1 = solid_func_CFC(x(i),y(j),z(k),xlength) 
          else if (solid_type == 'SingleSphere') then
             s1 = solid_func_one_sphere(x(i),y(j),z(k),xlength) 
           else
             stop 'invalid type'
          endif
          if(s1 > 0.) then
             solids(i,j,k) = 1.d0
          endif
       enddo; enddo; enddo
!       if(rank==0) print*, "+++x ",x 

       call output_solids(0,is,ie,js,je,ks,ke)
       if(rank==0) write(out,*)'solids type ', solid_type, ' initialized'
       if(rank==0) write(6,*)  'solids type ', solid_type, ' initialized'
       call ghost_x(solids,2,req(1:4));  call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ghost_y(solids,2,req(1:4));  call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ghost_z(solids,2,req(1:4));  call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
    endif
  END SUBROUTINE initialize_solids
  ! sphere definition function 
  FUNCTION sphere_func(x1,y1,z1,x0,y0,z0,radius)
    !***
    implicit none
    real(8) :: sphere_func
    real(8) , intent(in) :: x1,y1,z1,x0,y0,z0,radius
    sphere_func = -((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0) + (z1 - z0)*(z1 - z0)) + radius*radius
    return
  end function sphere_func

  ! example implicit solid definition function 
  FUNCTION solid_func_one_sphere(x1, y1, z1,boxL)
    implicit none
    real(8) :: solid_func_one_sphere
    real(8) , intent(in) :: x1,y1,z1,boxL
    real(8) :: x0,y0,z0,x2,y2,z2
    real(8) :: radius
    x2 = x1/boxL
    y2 = y1/boxL
    z2 = z1/boxL
    x0=0.5;y0=0.5;z0=0.5
    radius = 0.45
    solid_func_one_sphere = sphere_func(x2,y2,z2,x0,y0,z0,radius)
    return 
  end function  solid_func_one_sphere

  ! One basic CFC cell: one vertex + three faces
  FUNCTION solid_func_CFC_scaled( x1,  y1,  z1)
    implicit none
    intrinsic dmax1
    real(8), intent(in) :: x1,  y1,  z1
    real(8) :: solid_func_CFC_scaled
    real(8) :: a,b,x2,y2,z2
    real(8) :: radius 
    radius = 0.25*sqrt(2.d0)
    x2=x1;y2=y1;z2=z1
    ! 1 lower vertex and x=0 face 
    a = DMAX1(sphere_func(x2,y2,z2,0.d0,0.d0,0.d0,radius),sphere_func(x2,y2,z2,0.d0,0.5d0,0.5d0,radius))
    b = DMAX1(sphere_func(x2,y2,z2,0.5d0,0.d0,0.5d0,radius),sphere_func(x2,y2,z2,0.5d0,0.5d0,0.d0,radius))
    solid_func_CFC_scaled =  DMAX1(a,b)
    return 
  end function solid_func_CFC_scaled

  ! One basic cube with CFC array and all vertices and faces 
  FUNCTION solid_func_CFC( x1,  y1,  z1,  boxL)
    implicit none
    intrinsic dmax1
    real(8), intent(in) :: x1,  y1,  z1,  boxL
    real(8) :: solid_func_CFC
    real(8) ::  a,b,c,x2,y2,z2
    !  rescale by boxL
    x2 = x1/boxL
    y2 = y1/boxL
    z2 = z1/boxL
    !  shift origin from lower left corner to center
    x2 = x2+0.5 
    y2 = y2+0.5 
    z2 = z2+0.5 
    ! cells at root vertex and three first neighbors */
    a = DMAX1(solid_func_CFC_scaled(x2,y2,z2),solid_func_CFC_scaled(x2-1.d0,y2,z2))
    b = DMAX1(solid_func_CFC_scaled(x2,y2-1.d0,z2),solid_func_CFC_scaled(x2,y2,z2-1.d0))
    a = DMAX1(a,b)
    ! three next lower vertices 
    c = DMAX1(solid_func_CFC_scaled(x2-1.d0,y2-1.d0,z2-1.d0),solid_func_CFC_scaled(x2-1.d0,y2-1.d0,z2))
    b = DMAX1(solid_func_CFC_scaled(x2,y2-1.d0,z2-1.d0),solid_func_CFC_scaled(x2-1.d0,y2,z2-1.d0))
    a = DMAX1(c,a)
    solid_func_CFC =  DMAX1(b,a)
    return 
  end function solid_func_CFC
!***********************************************************************
    SUBROUTINE outfarray(carray) 
      INTEGER CARRAY(2,2) 
      
      WRITE(6,*) 'OUTFARRAY'
      WRITE(6,*) '1: ', CARRAY(1,1),' , ' , CARRAY(1,2) 
      WRITE(6,*) '2: ', CARRAY(2,1),' , ' , CARRAY(2,2) 
      WRITE(6,*) ' '
      STOP 'OUTFARRAY'
    END SUBROUTINE OUTFARRAY
!***********************************************************************
    subroutine ReadSolidParameters
      use module_flow
      use module_BC
      use module_IO
      implicit none
      include 'mpif.h'
      integer ierr,in
      logical file_is_there
      namelist /solidparameters/ dosolids, solid_type
      in=1
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
      inquire(file='inputsolids',exist=file_is_there)
      open(unit=in, file='inputsolids', status='old', action='read', iostat=ierr)
      if (file_is_there) then
         if(ierr == 0) then
            read(UNIT=in,NML=solidparameters)
            if(rank==0) write(out,*)'Solid parameters read successfully'
         else
            print *, 'rank=',rank,' has error ',ierr,' opening file inputsolids'
         endif
      else
         if (rank == 0) then 
            print *, "ReadSolidParameters: no 'inputsolids' file."
            print *, 'solids will not be activated.'
         endif
      endif
      close(in)
!      if(rank==0) print  * , "dosolids =",dosolids
    end subroutine ReadSolidParameters
!
! Output the velocity profile
!
  function test_point_in(i,j,k)
    integer, intent(in) :: i,j,k
    logical :: test_point_in
    test_point_in = (imin < i).and.(i < imax).and.(jmin < j).and.(j < jmax).and.(kmin < k).and.(k < kmax)
  end function test_point_in
!
  subroutine output_at_location()
    integer :: j,jproc
    integer :: imidline,jmidline,kmidline
    jproc = (jmin+jmax)*npy/(2*Nyt)
    jl=jmin
    if(jproc.eq.0) jl=ng
    jh=jmax
    if(jproc.eq.(Npy-1)) jh=Nyt-ng+1
    imidline = Nxt/2
    jmidline = Nyt/2   ! possible bug: why test j, Nyt need not be in processor's subdomain. 
    kmidline = Nzt/2
    if(test_point_in(imidline,jmidline,kmidline)) then
       OPEN(UNIT=11,FILE=trim(out_path)//'/output_location'//TRIM(int2text(jproc,padding)),status='unknown',action='write')
       do j=jl,jh
          write(11,1100) y(j),u(imidline,j,kmidline)
       enddo
       close(11)
    endif
    if(rank==0) then
       OPEN(UNIT=111,FILE=trim(out_path)//'/Poiseuille_theory',status='unknown',action='write')
       do j=ng,Nyt-ng+1
          write(111,1100) y(j), 0.5*y(j)*(1-y(j))
       enddo
       close(111)
    endif
1100 FORMAT(es25.16e3,es25.16e3)
  end subroutine output_at_location
!=================================================================================================
!-------------------------------------------------------------------------------------------------
  subroutine output_solids(nf,i1,i2,j1,j2,k1,k2)
    integer ::nf,i1,i2,j1,j2,k1,k2,i,j,k
    character(len=30) :: rootname
    rootname=trim(out_path)//'/VTK/solid'//TRIM(int2text(nf,padding))//'-'
    if(rank==0) call append_solid_visit_file(TRIM(rootname))

    OPEN(UNIT=8,FILE=TRIM(rootname)//TRIM(int2text(rank,padding))//'.vtk')
    write(8,10)
    write(8,11)time
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
    write(8,17)'solids'
    write(8,18)
16  format('POINT_DATA ',I17)
17  format('SCALARS ',A20,' float 1')
20  format('VECTORS uv float')
18  format('LOOKUP_TABLE default')

    do k=k1,k2; do j=j1,j2; do i=i1,i2;
      write(8,210) solids(i,j,k)
    enddo; enddo; enddo
210 format(e14.5)
310 format(e14.5,e14.5,e14.5)

    close(8)
    if(rank==0) call close_solid_visit_file()
end subroutine output_solids
!***********************************************************************
end module module_solid
