
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
  use module_2phase
  implicit none
  real(8), dimension(:,:,:), allocatable :: solids
  real(8) :: solid_radius
  character(20) :: solid_type
  integer :: il,ih,jl,jh,kl,kh
  integer :: solid_opened=0 
  real(8) :: sxyzrad(4,10000)
  integer :: NumSpheres
  logical :: bitmap_opened=.false.
  integer :: remove_layers=0
  real(8) :: porosity
  logical :: DoOutputSolids
!***********************************************************************
contains
!***********************************************************************
  subroutine open_bitmap()
    use module_tmpvar
    implicit none
    include 'mpif.h'
    integer bit,i,j,k
    open(unit=89,file='bitmap-'//trim(int2text(rank,padding))//'.txt')
    do i=is,ie; do j=js,je; do k=ks,ke; 
       read(89,'(I1)') bit
       if(bit<0.or.bit>1) then
          call pariserror("wrong bit")
       endif
       bit = 1 - bit
       tmp(i,j,k) = dble(bit)
    enddo;enddo;enddo
    close(89)
  end subroutine open_bitmap
!=================================================================================================
  function read_bitmap(i,j,k)
    use module_grid
    use module_tmpvar
    implicit none
    integer, intent(in) :: i,j,k
    real(8) read_bitmap
    if(.not.bitmap_opened) then 
       call open_bitmap
       bitmap_opened=.true.
    endif
    read_bitmap=tmp(i,j,k)
    if(i<remove_layers+Ng.or.i>Nx+Ng-remove_layers) read_bitmap=0d0
  end function read_bitmap
!=================================================================================================
  SUBROUTINE append_solid_visit_file(rootname)
    implicit none
    character(*) :: rootname
    integer prank
    if(rank.ne.0) call pariserror("rank.ne.0 in append_solid")

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
    use module_timer
    use module_2phase
    implicit none
    include 'mpif.h'
    integer :: i,j,k,index
    real(8) :: s1
    integer :: ierr
    integer :: req(48),sta(MPI_STATUS_SIZE,48)
    real(8) :: x0,y0,z0,radius,x2,y2,z2
    ! for 2d nozzle
    real(8), parameter :: PI = 3.14159265359d0
    real(8) :: NozzleThickness,ryz

    call ReadSolidParameters
    if(dosolids) then
       allocate(solids(imin:imax,jmin:jmax,kmin:kmax))
       call add_2_my_sizer(1,8)
       solids=0d0
       !call calcsum(solids)
       !if(rank==0) print *, "all 1 : ", s1
       !call print_small_solid(solids)
       do i=is,ie; do j=js,je; do k=ks,ke; 
          if(solid_type == 'CFC') then
             s1 = solid_func_CFC(x(i),y(j),z(k),xlength) 
          else if (solid_type == 'SingleSphere') then
             s1 = solid_func_one_sphere(x(i),y(j),z(k),xlength) 
          else if (solid_type == 'SingleDisk') then
             s1 = solid_func_one_disk(x(i),y(j),xlength) 
          else if (solid_type == 'SingleSquare') then
             s1 = solid_func_one_square(x(i),y(j),xlength) 
          else if (solid_type == 'List_of_Spheres') then
             s1=-1d0
             do index=1,NumSpheres
                x0 = sxyzrad(1,index)/xlength
                y0 = sxyzrad(2,index)/xlength
                z0 = sxyzrad(3,index)/xlength
                x2 = x(i)/xlength
                y2 = y(j)/xlength
                z2 = z(k)/xlength
                radius = sxyzrad(4,index)/xlength
                s1 = MAX(s1, sphere_func(x2,y2,z2,x0,y0,z0,radius))
             enddo
          else if (solid_type == 'BitMap') then
             s1 = read_bitmap(i,j,k)
          else if (solid_type == '2D_nozzle') then
             NozzleThickness = NozzleThick2Cell*dx(is)
             if ( x(i) < NozzleLength .and. & 
                  y(j) < radius_liq_inject+NozzleThickness .and. & 
                  y(j) > radius_liq_inject ) then 
               s1 = 1.d0 
             else
               s1 =-1.d0
             end if ! x(i),y(j)
          else if (solid_type == '3D_nozzle') then
             ryz = sqrt( (y(j) - jetcenter_yc)**2.d0 + (z(k) - jetcenter_zc)**2.d0 )
             if(radius_gap_liqgas==0d0) then
                NozzleThickness = NozzleThick2Cell*dx(is)
                if ( x(i) < NozzleLength .and. & 
                     ryz < radius_liq_inject+NozzleThickness .and. & 
                     ryz > radius_liq_inject ) then 
                  s1 = 1.d0 
                else
                  s1 =-1.d0
                end if ! x(i),y(j)
             else
                if ( x(i) < NozzleLength .and. & 
                     ryz < radius_gap_liqgas .and. & 
                     ryz > radius_liq_inject ) then 
                  s1 = 1.d0 
                else
                  s1 =-1.d0
                end if ! x(i),y(j)
             end if ! radius_gap_liqgas 
          else if (solid_type == 'pipe') then
             radius = 0.5d0*xLength
             s1 = ((x(i) - radius)*(x(i) -radius) + (z(k) - radius)*(z(k) - radius)) - radius*radius
          else
             call pariserror("invalid type")
             s1 =0d0 ! to avoid -Wmaybe-uninitialized
          endif
          if(s1 > 0.) then
             solids(i,j,k) = 1d0
          endif
!          solids(i,j,k) = s1
       enddo; enddo; enddo

       if(rank==0) write(out,*)'solids type ', solid_type, ' initialized'
       if(rank==0) write(6,*)  'solids type ', solid_type, ' initialized'
       call ghost_x(solids,2,req(1:4));  call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ghost_y(solids,2,req(1:4));  call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       call ghost_z(solids,2,req(1:4));  call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
       if (DoOutputSolids) call output_solids(0,is,ie+1,js,je+1,ks,ke+1)

       ! For solid objects set mask according to placement of solids
       umask=1d0; vmask=1d0; wmask=1d0
       do i=imin,imax-1; do j=jmin,jmax-1; do k=kmin,kmax-1
          if((solids(i,j,k) + solids(i+1,j,k)) > 0.5d0) umask(i,j,k) = 0d0
          if((solids(i,j,k) + solids(i,j+1,k)) > 0.5d0) vmask(i,j,k) = 0d0
          if((solids(i,j,k) + solids(i,j,k+1)) > 0.5d0) wmask(i,j,k) = 0d0
       enddo; enddo; enddo
       if (DoOutputSolids) call printpor(solids)
       if (DoOutputSolids) call calcpor(solids,0,porosity)
   endif
  end subroutine initialize_solids
!============================================================================================================
  subroutine printpor(solids)
    implicit none
    include 'mpif.h'
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: solids
    integer i,j,k
    open(unit=105,file=trim(out_path)//'/bitmap-'//trim(int2text(rank,padding))//'.txt')!,&
!         status='new',action='write')
    do i=is,ie
       do j=js,je; 
          do k=ks,ke; 
             if(solids(i,j,k) <= 0.5d0) then   ! free space
                write(105,'("1")') ! 1 is free (fluid) space
             else
                write(105,'("0")') 
             endif
          enddo
       enddo
    enddo
    close(105)
  end subroutine printpor
  !============================================================================================================
!
! sphere definition function 
!
!============================================================================================================
  FUNCTION sphere_func(x1,y1,z1,x0,y0,z0,radius)
    !***
    implicit none
    real(8) :: sphere_func
    real(8) , intent(in) :: x1,y1,z1,x0,y0,z0,radius
    sphere_func = -((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0) + (z1 - z0)*(z1 - z0)) + radius*radius
    return
  end function sphere_func

  FUNCTION disk_func(x1,y1,x0,y0,radius)
    !***
    implicit none
    real(8) :: disk_func
    real(8) , intent(in) :: x1,y1,x0,y0,radius
    disk_func = -((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0)) + radius*radius
    return
  end function disk_func

  FUNCTION square_func(x1,y1,x0,y0,radius)
    !***
    implicit none
    real(8) :: square_func
    real(8) , intent(in) :: x1,y1,x0,y0,radius
    square_func = - max(abs(x1-x0),abs(y1-y0)) + radius
    return
  end function square_func


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
    radius = solid_radius
    solid_func_one_sphere = sphere_func(x2,y2,z2,x0,y0,z0,radius)
    return 
  end function  solid_func_one_sphere

 FUNCTION solid_func_one_disk(x1, y1, boxL)
    implicit none
    real(8) :: solid_func_one_disk
    real(8) , intent(in) :: x1,y1,boxL
    real(8) :: x0,y0,x2,y2
    real(8) :: radius
    x2 = x1/boxL
    y2 = y1/boxL
    x0=0.5;y0=0.5
    radius = solid_radius
    solid_func_one_disk = disk_func(x2,y2,x0,y0,radius)
    return 
  end function  solid_func_one_disk

 FUNCTION solid_func_one_square(x1, y1, boxL)
    implicit none
    real(8) :: solid_func_one_square
    real(8) , intent(in) :: x1,y1,boxL
    real(8) :: x0,y0,x2,y2
    real(8) :: radius
    x2 = x1/boxL
    y2 = y1/boxL
    x0=0.5;y0=0.5
    radius = solid_radius
    solid_func_one_square = square_func(x2,y2,x0,y0,radius)
    return 
  end function  solid_func_one_square

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
      CALL PARISERROR("OUTFARRAY")
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
      namelist /solidparameters/ dosolids, solid_type, solid_radius,NumSpheres,sxyzrad, &
          remove_layers,DoOutputSolids
      dosolids=.false.
      solid_type='SingleSphere'
      solid_radius=0.1d0
      NumSpheres=1;
      sxyzrad=0d0
      remove_layers=0
      DoOutputSolids=.true.

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
    end subroutine ReadSolidParameters
!
! Output the velocity profile
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
    if(test_point_in_wide(imidline,jmidline,kmidline)) then
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
    INCLUDE 'mpif.h'
    integer ::nf,i1,i2,j1,j2,k1,k2,i,j,k
    integer :: ierr
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
15  format('POINTS ',I17,' double' )

    do k=k1,k2; do j=j1,j2; do i=i1,i2;
      write(8,320) x(i),y(j),z(k)
    enddo; enddo; enddo
320 format(e14.5,e14.5,e14.5)

    write(8,16)(i2-i1+1)*(j2-j1+1)*(k2-k1+1)
    write(8,17)'solids'
    write(8,18)
16  format('POINT_DATA ',I17)
17  format('SCALARS ',A20,' float 1')
! 20  format('VECTORS uv float')
18  format('LOOKUP_TABLE default')

    do k=k1,k2; do j=j1,j2; do i=i1,i2;
      write(8,210) solids(i,j,k)
    enddo; enddo; enddo
210 format(e14.5)
! 310 format(e14.5,e14.5,e14.5)

    close(8)
    if(rank==0) call close_solid_visit_file()
       
    ! zip and tar data
    if ( zip_data ) then
       call system('gzip '//TRIM(rootname)//TRIM(int2text(rank,padding))//'.vtk')
       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       if ( rank == 0 ) then 
         call system('tar cvf solid.tar '//TRIM(out_path)//'/VTK/solid*.vtk.gz')
         call system('rm '//TRIM(out_path)//'/VTK/solid*.vtk.gz')
       end if ! rank 
    end if ! zip_data

end subroutine output_solids

  subroutine print_small_solid(solids)
    use module_grid
    implicit none
    INCLUDE 'mpif.h'
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: solids
    integer i,j,k,ins(4),jins, IERR

       do k=ks-1,ke+1
          if(rank==0) write(6,*) '-----' , k
          do i=is-1,ie+1
             ins = 0
             jins = 1
             do j=js-1,je+1; 
                if(solids(i,j,k) > 0.5)  ins(jins)=1
                jins = jins+1
             enddo
             if(rank==0) write(6,'(I1,"  ",4(I1," "))') i, ins
          enddo
          call mpi_barrier(MPI_COMM_CART, ierr)
          if(rank==1) write(6,*) '*****' , k
          do i=is-1,ie+1
             ins = 0
             jins = 1
             do j=js-1,je+1; 
                if(solids(i,j,k) > 0.5)  ins(jins)=1
                jins = jins+1
             enddo
             if(rank==1) write(6,'(I1,"  ",4(I1," "))') i, ins
          enddo
       enddo

    if(rank==0) write(6,*) '@@@@@@@@@@@@@@' , k
  end subroutine print_small_solid
!=================================================================================================
subroutine final_output(flowrate)
  use module_grid
  use module_IO
  implicit none
  real(8) flowrate
  open(unit=121,file=TRIM(out_path)//'/flowrate.txt',position='append')
  write(121,'(e16.8,e16.8)') dx(3),flowrate
  close(121)
 end subroutine final_output
!=================================================================================================
end module module_solid

subroutine calcpor(smask,type,porosity)
  use module_grid
  use module_BC
  use module_IO
  implicit none
  include 'mpif.h'
  integer, intent(in) :: type
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: smask
  real(8) :: sum, totalsum
  real(8), intent(out) :: porosity
  real(8) volume
  integer :: i,j,k, ierr
  sum=0.d0
  volume=Nx*Ny*Nz
  do k=ks,ke; do j=js,je; do i=is,ie;
     sum=sum+smask(i,j,k)
  enddo;enddo;enddo
  call MPI_REDUCE(sum, totalsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_Domain, ierr)
  porosity=type*totalsum/volume + (1-type)*(1d0 - totalsum/volume)
  if(rank==0) then 
     open(UNIT=89,file=trim(out_path)//'/porosity.txt')
     write(89,310)  porosity
     close(unit=89)
  endif
  310 format(F17.11)
end subroutine calcpor
