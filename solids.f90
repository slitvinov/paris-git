!=================================================================================================
! module_solids: Contains definition of variables for the solids.
!-------------------------------------------------------------------------------------------------
module module_solids
  use module_grid
  use module_hello
  implicit none
  real(8), dimension(:,:,:), allocatable :: solids
  logical :: dosolids = .false.
!***********************************************************************
  interface
    real(c_double) FUNCTION solid_func_CFC(xdd,ydd,zdd,ldd) bind(C, name="solid_func_CFC_")
      use iso_c_binding, only: c_double
      real(c_double) , VALUE :: xdd,ydd,zdd,ldd
    END FUNCTION solid_func_CFC
 end interface
contains
!***********************************************************************
    SUBROUTINE initsolids()
      integer :: i,j,k
      real(8) :: ttt
      real(8) :: xx,xy,xz,xl
!***
!   test the C and Fortran array indexing- for development only
!***
!      call comparef2c()
!     STOP
!*** 
!     end test section
!***
      allocate(solids(imin:imax,jmin:jmax,kmin:kmax))
      do i=imin,imax; do j=jmin,jmax; do k=kmin,kmax; 
         solids(i,j,k) = solid_func_CFC(x(i),y(j),z(k),xlength)
      enddo; enddo; enddo
    END SUBROUTINE initsolids
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
  integer ierr,in
  namelist /solidparameters/ dosolids
  in=1
  open(unit=in, file='input', status='old', action='read', iostat=ierr)
  if (ierr .ne. 0) stop 'ReadSolidParameters: error opening input file'
  read(UNIT=in,NML=solidparameters)
  close(in)
  print  * , "dosolids =",dosolids
end subroutine ReadSolidParameters
end module module_solids

module module_output_solids
    use module_IO
    use module_flow
    use module_grid
    use module_solids
  implicit none
  integer il,ih,jl,jh,kl,kh
  integer :: padding=3
contains
!=================================================================================================
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
    jl=jmin
    jh=jmax
    imidline = Nxt/2
    jmidline = Nyt/2
    kmidline = Nzt/2
    if(test_point_in(imidline,jmidline,kmidline)) then
       jproc = (jmin+jmax)*npy/(2*Nyt)
       OPEN(UNIT=11,FILE=trim(out_path)//'/output_location'//TRIM(int2text(jproc,padding)),status='unknown',action='write')
       do j=jl,jh
          write(11,1100) y(j),u(imidline,j,kmidline)
       enddo
       close(11)
    endif
1100 FORMAT(es25.16e3,es25.16e3)
  end subroutine output_at_location
!=================================================================================================
!-------------------------------------------------------------------------------------------------
  subroutine output_solids(nf,i1,i2,j1,j2,k1,k2)
    integer ::nf,i1,i2,j1,j2,k1,k2,i,j,k
    character(len=30) :: rootname
    rootname=trim(out_path)//'/solid'//TRIM(int2text(nf,padding))//'-'
    if(rank==0) call append_visit_file(TRIM(rootname),padding)

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
end subroutine output_solids
!***********************************************************************
end module module_output_solids
