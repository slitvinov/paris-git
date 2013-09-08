!=================================================================================================
!=================================================================================================
!=================================================================================================
! Paris-0.1
! Extended from Code: FTC3D2011 (Front Tracking Code for 3D simulations)
! 
! Authors: Sadegh Dabiri, Gretar Tryggvason
! author for minor extenstions Stephane Zaleski(zaleski@dalembert.upmc.fr) 
! Contact: sdabiri@gmail.com
! A three-dimensional Navier-Stokes flow solver for modeling of multiphase 
! flows. Flow can be driven by wall motion, density difference or pressure gradient.
! Boundary conditions supported: wall and periodic
! Version 1.0   1/21/2011   The 3D flow solver for variable density/viscosity is written. 
!                           The density is advected by an ENO scheme.
! Version 2.0   2/25/2011   Parallel implementation.
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
!=================================================================================================
!=================================================================================================
! module_grid: Contains definition of variables for the grid.
!-------------------------------------------------------------------------------------------------
module module_grid
  implicit none
  save
  integer :: Nx, Ny, Nz, Ng ! Ng is the number of ghost cells
  integer :: Nxt, Nyt, Nzt ! total number of cells
  integer :: is, ie, js, je, ks, ke
  integer :: ieu, jev, kew
  real(8), dimension(:), allocatable :: x, xh, dx, dxh
  real(8), dimension(:), allocatable :: y, yh, dy, dyh
  real(8), dimension(:), allocatable :: z, zh, dz, dzh
  real(8) :: xLength, yLength, zLength, xform, yform, zform !non-uniformity of grid
  logical :: TwoPhase
  integer :: nPx, nPy, nPz, Mx, My, Mz, rank, ndim=3, nPdomain, NumProcess
  integer, dimension(:), allocatable :: dims, coords, periodic, reorder
  integer :: MPI_Comm_Cart, MPI_Comm_Domain, MPI_Comm_Active
  integer :: imin, imax, jmin, jmax, kmin, kmax
end module module_grid

!=================================================================================================
! module_hello: Contains definition of variables and subroutines to say hello
! This is useful for debugging. 
!-------------------------------------------------------------------------------------------------
module module_hello
  implicit none
  integer :: hello_count = 1
  contains
subroutine hello_coucou
  use module_grid
  integer, parameter  :: debug=1
  if(debug == 1) then 
  if(rank==0) write(6,*) 'coucou ',hello_count, "Process0"
  if(rank==nPdomain) write(6,*) 'coucou ',hello_count, "Front"
  hello_count = hello_count + 1
  end if
end subroutine hello_coucou
end module module_hello
!=================================================================================================
!=================================================================================================
! module_flow: Contains definition of variables for the flow solver.
!-------------------------------------------------------------------------------------------------
module module_flow
  implicit none
  save
  real(8), dimension(:,:,:), allocatable :: u, v, w, uold, vold, wold, fx, fy, fz, color
  real(8), dimension(:,:,:), allocatable :: p, rho, rhoo, muold, mu, dIdx, dIdy, dIdz
  real(8), dimension(:,:,:), allocatable :: umask,vmask,wmask
  real(8), dimension(:,:,:), allocatable :: du,dv,dw,drho
  real(8), dimension(:,:), allocatable :: averages,oldaverages, allaverages
  logical, allocatable, dimension(:,:,:) :: mask

  real(8) :: gx, gy, gz, mu1, mu2, r_avg, dt, dtFlag, rho_ave, p_ave, vdt
  real(8) :: max_velocity, maxTime, Time, EndTime, MaxDt, CFL, mystats(16), stats(16)
  logical :: ZeroReynolds,DoVOF, DoFront, Implicit, hypre, GetPropertiesFromFront
  logical :: dosolids = .false.
  real(8) :: rho1, rho2, s
  real(8) :: U_init
  real(8) :: dpdx, dpdy, dpdz, W_ave  !pressure gradients in case of pressure driven channel flow
  real(8) :: dpdx_stat, dpdy_stat, dpdz_stat
  real(8) :: beta, MaxError
  integer :: maxit, it, itime_scheme, BuoyancyCase, drive
  integer :: sbx, sby, Nstep
  integer :: maxStep, itmax, iTimeStep
end module module_flow
!=================================================================================================
!=================================================================================================
! Temporary variables
!-------------------------------------------------------------------------------------------------
module module_tmpvar
  real(8), dimension(:,:,:,:), allocatable :: work, A
  real(8), dimension(:,:,:), allocatable :: tmp
  real(8) :: tcpu(100),t0
end module module_tmpvar

!=================================================================================================
!=================================================================================================
! module_2phase: Contains variables for two-phase flow
!-------------------------------------------------------------------------------------------------
module module_2phase
  real(8), dimension( : ), allocatable :: rad, xc, yc, zc, vol
  real(8) :: sigma
  integer :: NumBubble
end module module_2phase
!=================================================================================================
!=================================================================================================
! module_IO: Contains input/output variables and procedures
!-------------------------------------------------------------------------------------------------
module module_IO
  implicit none
  save
  integer :: padding=5
  integer :: opened=0;
  integer ::nout, out, output_format, nbackup, nstats, termout
  character(len=20) :: out_path, x_file, y_file, z_file
  logical :: read_x, read_y, read_z, restart, ICOut, restartFront, restartAverages
  contains
!=================================================================================================
subroutine write_vec_gnuplot(u,v,w,iout)
  use module_grid
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: u, v, w
  integer, intent(in) :: iout
  integer :: i,j,k
  real(8) :: norm=0.d0, coeff
  intrinsic dsqrt
  OPEN(UNIT=89,FILE=TRIM(out_path)//'/UV-'//TRIM(int2text(rank,padding))//'-'//TRIM(int2text(iout,padding))//'.txt')
  norm=0.d0
  k=nz/2+ng
  do i=imin,imax; do j=jmin,jmax
     norm = norm + u(i,j,k)**2 + v(i,j,k)**2
  enddo; enddo
  norm = dsqrt(norm/((jmax-jmin+1)*(imax-imin+1)))
  coeff = 0.8/norm
  ! print *," norm, rank ",norm,rank
  do i=is,ie; do j=js,je
     write(89,310) x(i),y(j),coeff*dx(i)*u(i,j,k),coeff*dx(i)*v(i,j,k)
     ! write(89,311) i,j,k, u(i,j,k),v(i,j,k)
  enddo; enddo
  close(unit=89)
310 format(e14.5,e14.5,e14.5,e14.5)
311 format(I2,I2,I2,e14.5,e14.5)
end subroutine  write_vec_gnuplot
!=================================================================================================
! append
    SUBROUTINE append_visit_file(rootname)
      use module_flow
      use module_grid
    implicit none
    character(*) :: rootname
    integer prank
    if(rank.ne.0) stop 'rank.ne.0 in append'

    if(opened==0) then
       OPEN(UNIT=90,FILE='parallel.visit')
       write(90,10) nPdomain
10     format('!NBLOCKS ',I4)
       opened=1
    endif

    do prank=0,NpDomain-1
       write(90,11) rootname//TRIM(int2text(prank,padding))//'.vtk'
 11 format(A)
    enddo
  end subroutine  append_visit_file

  subroutine close_visit_file()
    close(90)
  end subroutine close_visit_file
!=================================================================================================
! function int2text
!   Returns 'number' as a string with length of 'length'
!   called in:    function output
!-------------------------------------------------------------------------------------------------
function int2text(number,length)
  integer :: number, length, i
  character(len=length) :: int2text
  character, dimension(0:9) :: num = (/'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/)
  if(number>=10**length)print*, 'Warning: int2text: large input number. Increase "length".'
  do i=1,length
    int2text(length+1-i:length+1-i) = num(mod(number/(10**(i-1)),10))
  enddo
end function
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
subroutine backup_write
  use module_flow
  use module_grid
  implicit none
  integer ::i,j,k
  character(len=100) :: filename
  filename = trim(out_path)//'/backup_'//int2text(rank,3)
  call system('mv '//trim(filename)//' '//trim(filename)//'.old')
  OPEN(UNIT=7,FILE=trim(filename),status='unknown',action='write')
  write(7,1100)time,itimestep,is,ie,js,je,ks,ke
  do k=ks,ke; do j=js,je; do i=is,ie
    write(7,1200) u(i,j,k), v(i,j,k), w(i,j,k), p(i,j,k), color(i,j,k)
  enddo; enddo; enddo
  do i=1,10; do j=js,je
    write(7,1200) averages(i,j)
  enddo; enddo
  close(7)
  if(rank==0)print*,'Backup written at t=',time
  1100 FORMAT(es17.8e3,7I10)
  1200 FORMAT(5es17.8e3)
end subroutine backup_write
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
subroutine backup_read
  use module_flow
  use module_grid
  use module_hello
  implicit none
  integer ::i,j,k,i1,i2,j1,j2,k1,k2,ierr
  OPEN(UNIT=7,FILE=trim(out_path)//'/backup_'//int2text(rank,3),status='old',action='read')
  read(7,*)time,itimestep,i1,i2,j1,j2,k1,k2
  if(i1/=is .or. i2/=ie .or. j1/=js .or. j2/=je .or. k1/=ks .or. k2/=ke) &
    stop 'Error: backup_read'
  do k=ks,ke; do j=js,je; do i=is,ie
    read(7,*) u(i,j,k), v(i,j,k), w(i,j,k), p(i,j,k), color(i,j,k)
  enddo; enddo; enddo
  if(restartAverages)then
    do i=1,20; do j=js,je
      read(7,*,iostat=ierr) averages(i,j)
      if(ierr/=0)return !no average data, return
    enddo; enddo
    !oldaverages=averages
  endif
  close(7)
  1100 FORMAT(es25.16e3,7I10)
  1200 FORMAT(5es25.16e3)
end subroutine backup_read
!=================================================================================================
!=================================================================================================
! subroutine output
!   Writes the solution in the 'out_path' directory in the tecplot or vtk format
!   called in:    program main
!-------------------------------------------------------------------------------------------------
subroutine output(nf,i1,i2,j1,j2,k1,k2)
  implicit none
  integer :: nf,i1,i2,j1,j2,k1,k2
  if(output_format==1) call output1(nf,i1,i2,j1,j2,k1,k2)
  if(output_format==2) call output2(nf,i1,i2,j1,j2,k1,k2)
end subroutine output
!-------------------------------------------------------------------------------------------------
subroutine output1(nf,i1,i2,j1,j2,k1,k2)
  use module_flow
  use module_grid
  use module_hello
!  use module_tmpvar
  !use IO_mod
  implicit none
  integer :: nf,i1,i2,j1,j2,k1,k2,i,j,k
  logical, save :: first_time=.true.
  
  if(first_time)then
    first_time = .false.
    OPEN(UNIT=7,FILE=trim(out_path)//'/grid_'//int2text(rank,3)//'.dat')
    write(7,*) 'FILETYPE = GRID'
    write(7,*) 'variables = "x","y","z"'
    write(7,2100) i2-i1+1, j2-j1+1, k2-k1+1
    do k=k1,k2; do j=j1,j2; do i=i1,i2;
      write(7,1200) x(i),y(j),z(k)
    enddo; enddo; enddo
    close(7)
  endif

  OPEN(UNIT=7,FILE=trim(out_path)//'/plot'//int2text(nf,3)//'_'//int2text(rank,3)//'.dat')
  write(7,*) 'FILETYPE = SOLUTION'
  write(7,*) 'variables = "u","v","w","p","color"'
  write(7,1100) time, i2-i1+1, j2-j1+1, k2-k1+1
  do k=k1,k2; do j=j1,j2; do i=i1,i2;
    write(7,1200) 0.5d0*(u(i,j,k)+u(i-1,j,k)), &
                  0.5d0*(v(i,j,k)+v(i,j-1,k)), &
                  0.5d0*(w(i,j,k)+w(i,j,k-1)), &
                  p(i,j,k) , color(i,j,k)
  enddo; enddo; enddo
  close(7)
  if(rank==0)then
    OPEN(UNIT=7,FILE=trim(out_path)//'/averages.dat',access='append')
    write(7,1110) time
    do j=Ng,Ng+Ny+1
      write(7,1200) y(j),real(allaverages(:,j)-oldaverages(:,j))
    enddo
    close(7)
    oldaverages=allaverages
  endif
  1000 FORMAT('variables = "x","y","z","u","v","w","p","color"')
  1100 FORMAT('ZONE solutiontime=', 1PG15.7e2, ', I=',I4, 2X, ', J=',I4, 2X, ', K=',I4, 2X)
  1110 FORMAT('ZONE solutiontime=', 1PG15.7e2)
  1200 FORMAT(21es14.6e2)
  2100 FORMAT('ZONE I=',I4, 2X, ', J=',I4, 2X, ', K=',I4, 2X)
end subroutine output1
!-------------------------------------------------------------------------------------------------
subroutine output2(nf,i1,i2,j1,j2,k1,k2)
  use module_flow
  use module_grid
  use module_hello
  !use IO_mod
  implicit none
  integer ::nf,i1,i2,j1,j2,k1,k2,i,j,k, itype=5
!  logical, save :: first_time=.true.
  character(len=30) :: rootname
  rootname=TRIM(out_path)//'/VTK/plot'//TRIM(int2text(nf,padding))//'-'

  if(rank==0) call append_visit_file(TRIM(rootname))

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
    if (itype .le. 4)then
      write(8,17)'density'
      write(8,18)
    else
      write(8,20)
    endif
16  format('POINT_DATA ',I17)
17  format('SCALARS ',A20,' float 1')
20  format('VECTORS uv float')
18  format('LOOKUP_TABLE default')

    do k=k1,k2; do j=j1,j2; do i=i1,i2;
      if (itype .eq. 1)write(8,210)rho(i,j,k)
      if (itype .eq. 5)write(8,310)0.5*(u(i,j,k)+u(i-1,j,k)), &
         0.5*(v(i,j,k)+v(i,j-1,k)),0.5*(w(i,j,k)+w(i,j,k-1))
    enddo; enddo; enddo
210 format(e14.5)
310 format(e14.5,e14.5,e14.5)

    close(8)
end subroutine output2
!-------------------------------------------------------------------------------------------------
end module module_IO
!=================================================================================================
!=================================================================================================
! module module_BC: Sets the boundary conditions of the problem.
!-------------------------------------------------------------------------------------------------
module module_BC
  use module_grid
  use module_hello
  implicit none
  integer :: bdry_cond(6)
  ! bdry_cond(i) = is the type if boundary condition in i'th direction
  ! 0:wall;  1:periodic
  real(8) :: WallVel(6,3), WallShear(6,3)
  ! Tangential velocities on the surfaces of domain. First index represent the 
  ! side on which the velocity in the direction of the second index is specified.
  ! The sides are in this order: -x,+x,-y,+y,-z,+z.
  ! Example: WallVel(4,3) represent the W velocity on +y side of the domain.
  contains
!=================================================================================================
!=================================================================================================
! subroutine SetPressureBC: Sets the pressure boundary condition
!-------------------------------------------------------------------------------------------------
    subroutine SetPressureBC(umask,vmask,wmask)
    use module_grid
    implicit none
    include 'mpif.h'
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: umask,vmask,wmask

  ! for walls set the mask to zero
    if(bdry_cond(1)==0)then
      if(coords(1)==0    ) umask(is-1,js-1:je+1,ks-1:ke+1)=0d0
      if(coords(1)==nPx-1) umask(ie,js-1:je+1,ks-1:ke+1)=0d0
    endif

    if(bdry_cond(2)==0)then
      if(coords(2)==0    ) vmask(is-1:ie+1,js-1,ks-1:ke+1)=0d0
      if(coords(2)==nPy-1) vmask(is-1:ie+1,je,ks-1:ke+1)=0d0
    endif

    if(bdry_cond(3)==0)then
      if(coords(3)==0    ) wmask(is-1:ie+1,js-1:je+1,ks-1)=0d0
      if(coords(3)==nPz-1) wmask(is-1:ie+1,js-1:je+1,ke)=0d0
    endif
  end subroutine SetPressureBC
!=================================================================================================
!=================================================================================================
! subroutine SetVelocityBC: Sets the velocity boundary condition
!-------------------------------------------------------------------------------------------------
  subroutine SetVelocityBC(u,v,w,umask,vmask,wmask)
    use module_grid
    use module_hello
    implicit none
    include 'mpif.h'
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: u, v, w
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: umask,vmask,wmask
    ! solid obstacles
    u = u*umask
    v = v*vmask
    w = w*wmask

    ! wall boundary condition
    if(bdry_cond(1)==0 .and. coords(1)==0    ) then
        u(is-1,:,:)=0d0
        u(is-2,:,:)=-u(is,:,:)
        v(is-1,:,:)=2*WallVel(1,2)-v(is,:,:)
        w(is-1,:,:)=2*WallVel(1,3)-w(is,:,:)
    endif
    if(bdry_cond(4)==0 .and. coords(1)==nPx-1) then
        u(ie  ,:,:)=0d0
        u(ie+1,:,:)=-u(ie-1,:,:)
        v(ie+1,:,:)=2*WallVel(2,2)-v(ie,:,:)
        w(ie+1,:,:)=2*WallVel(2,3)-w(ie,:,:)
    endif
    if(bdry_cond(2)==0 .and. coords(2)==0    ) then
        v(:,js-1,:)=0d0
        v(:,js-2,:)=-v(:,js,:)
        u(:,js-1,:)=2*WallVel(3,1)-u(:,js,:)
        w(:,js-1,:)=2*WallVel(3,3)-w(:,js,:)
    endif
    if(bdry_cond(5)==0 .and. coords(2)==nPy-1) then
        v(:,je  ,:)=0d0
        v(:,je+1,:)=-v(:,je-1,:)
        u(:,je+1,:)=2*WallVel(4,1)-u(:,je,:)
        w(:,je+1,:)=2*WallVel(4,3)-w(:,je,:)
    endif
    if(bdry_cond(3)==0 .and. coords(3)==0    ) then
        w(:,:,ks-1)=0d0
        w(:,:,ks-2)=-w(:,:,ks)
        u(:,:,ks-1)=2*WallVel(5,1)-u(:,:,ks)
        v(:,:,ks-1)=2*WallVel(5,2)-v(:,:,ks)
    endif
    if(bdry_cond(6)==0 .and. coords(3)==nPz-1) then
        w(:,:,ke  )=0d0
        w(:,:,ke+1)=-w(:,:,ke-1)
        u(:,:,ke+1)=2*WallVel(6,1)-u(:,:,ke)
        v(:,:,ke+1)=2*WallVel(6,2)-v(:,:,ke)
    endif
    ! wall boundary condition: shear
    if(bdry_cond(1)==2 .and. coords(1)==0    ) then
        v(is-1,:,:) = -dxh(is-1)*WallShear(1,2)+v(is,:,:)
        w(is-1,:,:) = -dxh(is-1)*WallShear(1,3)+w(is,:,:)
    endif
    if(bdry_cond(4)==2 .and. coords(1)==nPx-1) then
        v(ie+1,:,:) = dxh(ie)*WallShear(2,2)+v(ie,:,:)
        w(ie+1,:,:) = dxh(ie)*WallShear(2,3)+w(ie,:,:)
    endif
    if(bdry_cond(2)==2 .and. coords(2)==0    ) then
        u(:,js-1,:) = -dyh(js-1)*WallShear(3,1)+u(:,js,:)
        w(:,js-1,:) = -dyh(js-1)*WallShear(3,3)+w(:,js,:)
    endif
    if(bdry_cond(5)==2 .and. coords(2)==nPy-1) then
        u(:,je+1,:) = dyh(je)*WallShear(4,1)+u(:,je,:)
        w(:,je+1,:) = dyh(je)*WallShear(4,3)+w(:,je,:)
    endif
    if(bdry_cond(3)==2 .and. coords(3)==0    ) then
        u(:,:,ks-1) = -dzh(ks-1)*WallShear(5,1)+u(:,:,ks)
        v(:,:,ks-1) = -dzh(ks-1)*WallShear(5,2)+v(:,:,ks)
    endif
    if(bdry_cond(6)==2 .and. coords(3)==nPz-1) then
        u(:,:,ke+1) = dzh(ke)*WallShear(6,1)+u(:,:,ke)
        v(:,:,ke+1) = dzh(ke)*WallShear(6,2)+v(:,:,ke)
    endif
  end subroutine SetVelocityBC
!=================================================================================================
!=================================================================================================
! subroutine SetVelocityBC: Sets the velocity boundary condition
!-------------------------------------------------------------------------------------------------
  subroutine SetVectorBC(fx,fy,fz)
    use module_grid
    use module_hello
    implicit none
    include 'mpif.h'
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: fx, fy, fz
    ! Add the vector outside the domain to the neighboring cell inside the domain
    ! This is used for color gradient vector and surface tension forces
    if(bdry_cond(1)==0 .and. coords(1)==0    ) then
        fx(is-1,:,:)=fx(is-1,:,:)+fx(is-2,:,:)
        fy(is  ,:,:)=fy(is  ,:,:)+fy(is-1,:,:)+fy(is-2,:,:)
        fz(is  ,:,:)=fz(is  ,:,:)+fz(is-1,:,:)+fz(is-2,:,:)
    endif
    if(bdry_cond(4)==0 .and. coords(1)==nPx-1) then
        fx(ie,:,:)=fx(ie,:,:)+fx(ie+1,:,:)
        fy(ie,:,:)=fy(ie,:,:)+fy(ie+1,:,:)+fy(ie+2,:,:)
        fz(ie,:,:)=fz(ie,:,:)+fz(ie+1,:,:)+fz(ie+2,:,:)
    endif
    if(bdry_cond(2)==0 .and. coords(2)==0    ) then
        fy(:,js-1,:)=fy(:,js-1,:)+fy(:,js-2,:)
        fx(:,js  ,:)=fx(:,js  ,:)+fx(:,js-1,:)+fx(:,js-2,:)
        fz(:,js  ,:)=fz(:,js  ,:)+fz(:,js-1,:)+fz(:,js-2,:)
    endif
    if(bdry_cond(5)==0 .and. coords(2)==nPy-1) then
        fy(:,je,:)=fy(:,je,:)+fy(:,je+1,:)
        fx(:,je,:)=fx(:,je,:)+fx(:,je+1,:)+fx(:,je+2,:)
        fz(:,je,:)=fz(:,je,:)+fz(:,je+1,:)+fz(:,je+2,:)
    endif
    if(bdry_cond(3)==0 .and. coords(3)==0    ) then
        fz(:,:,ks-1)=fz(:,:,ks-1)+fz(:,:,ks-2)
        fx(:,:,ks  )=fx(:,:,ks  )+fx(:,:,ks-1)+fx(:,:,ks-2)
        fy(:,:,ks  )=fy(:,:,ks  )+fy(:,:,ks-1)+fy(:,:,ks-2)
    endif
    if(bdry_cond(6)==0 .and. coords(3)==nPz-1) then
        fz(:,:,ke)=fz(:,:,ke)+fz(:,:,ke+1)
        fx(:,:,ke)=fx(:,:,ke)+fx(:,:,ke+1)+fx(:,:,ke+2)
        fy(:,:,ke)=fy(:,:,ke)+fy(:,:,ke+1)+fy(:,:,ke+2)
    endif
  end subroutine SetVectorBC
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
  subroutine ghost_x(F,ngh,req)
    use module_grid
    use module_hello
    implicit none
    include 'mpif.h'
    integer, intent(in) :: ngh ! number of ghost cell layers to fill
    integer, intent(out) :: req(4)
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: F
    integer :: jlen, klen, ierr !,sta(MPI_STATUS_SIZE,4)
    integer, save :: srcL, srcR, destL, destR, face(2)
    logical, save :: first_time=.true.

    if(ngh>Ng) stop 'ghost error: not enough ghost layers to fill'
    if(first_time) then
      first_time=.false.
      jlen=jmax-jmin+1; klen=kmax-kmin+1; !ilen=ngh
      call para_type_block3a(imin, imax, jmin, jmax, 1, jlen, klen, MPI_DOUBLE_PRECISION, face(1))
      call para_type_block3a(imin, imax, jmin, jmax, 2, jlen, klen, MPI_DOUBLE_PRECISION, face(2))
      call MPI_CART_SHIFT(MPI_COMM_CART, 0, 1, srcR, destR, ierr)
      call MPI_CART_SHIFT(MPI_COMM_CART, 0,-1, srcL, destL, ierr)
    endif

    call MPI_IRECV(F(is-ngh  ,jmin,kmin),1,face(ngh),srcR ,0,MPI_COMM_Cart,req(1),ierr)
    call MPI_ISEND(F(ie-ngh+1,jmin,kmin),1,face(ngh),destR,0,MPI_COMM_Cart,req(2),ierr)
    call MPI_IRECV(F(ie+1    ,jmin,kmin),1,face(ngh),srcL ,0,MPI_COMM_Cart,req(3),ierr)
    call MPI_ISEND(F(is      ,jmin,kmin),1,face(ngh),destL,0,MPI_COMM_Cart,req(4),ierr)
!    call MPI_WAITALL(4,req,sta,ierr)
  end subroutine ghost_x
!-------------------------------------------------------------------------------------------------
  subroutine ghost_y(F,ngh,req)
    use module_grid
    use module_hello
    implicit none
    include 'mpif.h'
    integer, intent(in) :: ngh
    integer, intent(out) :: req(4)
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: F
    integer :: ilen, klen, ierr !,sta(MPI_STATUS_SIZE,4)
    integer, save :: srcL, srcR, destL, destR, face(2)
    logical, save :: first_time=.true.

    if(ngh>Ng) stop 'ghost error: not enough ghost layers to fill'
    if(first_time)then
      first_time=.false.
      klen=kmax-kmin+1; ilen=imax-imin+1; !jlen=ngh
      call para_type_block3a(imin, imax, jmin, jmax, ilen, 1, klen, MPI_DOUBLE_PRECISION, face(1))
      call para_type_block3a(imin, imax, jmin, jmax, ilen, 2, klen, MPI_DOUBLE_PRECISION, face(2))
      call MPI_CART_SHIFT(MPI_COMM_CART, 1, 1, srcR, destR, ierr)
      call MPI_CART_SHIFT(MPI_COMM_CART, 1,-1, srcL, destL, ierr)
    endif

    call MPI_IRECV(F(imin,js-ngh  ,kmin),1,face(ngh),srcR ,0,MPI_COMM_Cart,req(1),ierr)
    call MPI_ISEND(F(imin,je-ngh+1,kmin),1,face(ngh),destR,0,MPI_COMM_Cart,req(2),ierr)
    call MPI_IRECV(F(imin,je+1    ,kmin),1,face(ngh),srcL ,0,MPI_COMM_Cart,req(3),ierr)
    call MPI_ISEND(F(imin,js      ,kmin),1,face(ngh),destL,0,MPI_COMM_Cart,req(4),ierr)
!    call MPI_WAITALL(4,req,sta,ierr)
  end subroutine ghost_y
!-------------------------------------------------------------------------------------------------
  subroutine ghost_z(F,ngh,req)
    use module_grid
    use module_hello
    implicit none
    include 'mpif.h'
    integer, intent(in) :: ngh
    integer, intent(out) :: req(4)
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: F
    integer :: ilen, jlen, ierr !,sta(MPI_STATUS_SIZE,4)
    integer, save :: srcL, srcR, destL, destR, face(2)
    logical, save :: first_time=.true.

    if(ngh>Ng) stop 'ghost error: not enough ghost layers to fill'
    if(first_time)then
      first_time=.false.
      ilen=imax-imin+1; jlen=jmax-jmin+1; !klen=ngh
      call para_type_block3a(imin, imax, jmin, jmax, ilen, jlen, 1, MPI_DOUBLE_PRECISION, face(1))
      call para_type_block3a(imin, imax, jmin, jmax, ilen, jlen, 2, MPI_DOUBLE_PRECISION, face(2))
      call MPI_CART_SHIFT(MPI_COMM_CART, 2, 1, srcR, destR, ierr)
      call MPI_CART_SHIFT(MPI_COMM_CART, 2,-1, srcL, destL, ierr)
    endif

    call MPI_IRECV(F(imin,jmin,ks-ngh  ),1,face(ngh),srcR ,0,MPI_COMM_Cart,req(1),ierr)
    call MPI_ISEND(F(imin,jmin,ke-ngh+1),1,face(ngh),destR,0,MPI_COMM_Cart,req(2),ierr)
    call MPI_IRECV(F(imin,jmin,ke+1    ),1,face(ngh),srcL ,0,MPI_COMM_Cart,req(3),ierr)
    call MPI_ISEND(F(imin,jmin,ks      ),1,face(ngh),destL,0,MPI_COMM_Cart,req(4),ierr)
!    call MPI_WAITALL(4,req,sta,ierr)
  end subroutine ghost_z
!=================================================================================================
!=================================================================================================
  subroutine ighost_x(F,ngh,req)
    use module_grid
    use module_hello
    implicit none
    include 'mpif.h'
    integer, intent(in) :: ngh ! number of ghost cell layers to fill
    integer, intent(out) :: req(4)
    integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: F
    integer :: jlen, klen, ierr !,sta(MPI_STATUS_SIZE,4)
    integer, save :: srcL, srcR, destL, destR, face(2)
    logical, save :: first_time=.true.

    if(ngh>Ng) stop 'ghost error: not enough ghost layers to fill'
    if(first_time) then
      first_time=.false.
      jlen=jmax-jmin+1; klen=kmax-kmin+1; !ilen=ngh
      call para_type_block3a(imin, imax, jmin, jmax, 1, jlen, klen, MPI_INTEGER, face(1))
      call para_type_block3a(imin, imax, jmin, jmax, 2, jlen, klen, MPI_INTEGER, face(2))
      call MPI_CART_SHIFT(MPI_COMM_CART, 0, 1, srcR, destR, ierr)
      call MPI_CART_SHIFT(MPI_COMM_CART, 0,-1, srcL, destL, ierr)
    endif

    call MPI_IRECV(F(is-ngh  ,jmin,kmin),1,face(ngh),srcR ,0,MPI_COMM_Cart,req(1),ierr)
    call MPI_ISEND(F(ie-ngh+1,jmin,kmin),1,face(ngh),destR,0,MPI_COMM_Cart,req(2),ierr)
    call MPI_IRECV(F(ie+1    ,jmin,kmin),1,face(ngh),srcL ,0,MPI_COMM_Cart,req(3),ierr)
    call MPI_ISEND(F(is      ,jmin,kmin),1,face(ngh),destL,0,MPI_COMM_Cart,req(4),ierr)
!    call MPI_WAITALL(4,req,sta,ierr)
  end subroutine ighost_x
!-------------------------------------------------------------------------------------------------
  subroutine ighost_y(F,ngh,req)
    use module_grid
    use module_hello
    implicit none
    include 'mpif.h'
    integer, intent(in) :: ngh
    integer, intent(out) :: req(4)
    integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: F
    integer :: ilen, klen, ierr !,sta(MPI_STATUS_SIZE,4)
    integer, save :: srcL, srcR, destL, destR, face(2)
    logical, save :: first_time=.true.

    if(ngh>Ng) stop 'ghost error: not enough ghost layers to fill'
    if(first_time)then
      first_time=.false.
      klen=kmax-kmin+1; ilen=imax-imin+1; !jlen=ngh
      call para_type_block3a(imin, imax, jmin, jmax, ilen, 1, klen, MPI_INTEGER, face(1))
      call para_type_block3a(imin, imax, jmin, jmax, ilen, 2, klen, MPI_INTEGER, face(2))
      call MPI_CART_SHIFT(MPI_COMM_CART, 1, 1, srcR, destR, ierr)
      call MPI_CART_SHIFT(MPI_COMM_CART, 1,-1, srcL, destL, ierr)
    endif

    call MPI_IRECV(F(imin,js-ngh  ,kmin),1,face(ngh),srcR ,0,MPI_COMM_Cart,req(1),ierr)
    call MPI_ISEND(F(imin,je-ngh+1,kmin),1,face(ngh),destR,0,MPI_COMM_Cart,req(2),ierr)
    call MPI_IRECV(F(imin,je+1    ,kmin),1,face(ngh),srcL ,0,MPI_COMM_Cart,req(3),ierr)
    call MPI_ISEND(F(imin,js      ,kmin),1,face(ngh),destL,0,MPI_COMM_Cart,req(4),ierr)
!    call MPI_WAITALL(4,req,sta,ierr)
  end subroutine ighost_y
!-------------------------------------------------------------------------------------------------
  subroutine ighost_z(F,ngh,req)
    use module_grid
    use module_hello
    implicit none
    include 'mpif.h'
    integer, intent(in) :: ngh
    integer, intent(out) :: req(4)
    integer, dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: F
    integer :: ilen, jlen, ierr !,sta(MPI_STATUS_SIZE,4)
    integer, save :: srcL, srcR, destL, destR, face(2)
    logical, save :: first_time=.true.

    if(ngh>Ng) stop 'ghost error: not enough ghost layers to fill'
    if(first_time)then
      first_time=.false.
      ilen=imax-imin+1; jlen=jmax-jmin+1; !klen=ngh
      call para_type_block3a(imin, imax, jmin, jmax, ilen, jlen, 1, MPI_INTEGER, face(1))
      call para_type_block3a(imin, imax, jmin, jmax, ilen, jlen, 2, MPI_INTEGER, face(2))
      call MPI_CART_SHIFT(MPI_COMM_CART, 2, 1, srcR, destR, ierr)
      call MPI_CART_SHIFT(MPI_COMM_CART, 2,-1, srcL, destL, ierr)
    endif

    call MPI_IRECV(F(imin,jmin,ks-ngh  ),1,face(ngh),srcR ,0,MPI_COMM_Cart,req(1),ierr)
    call MPI_ISEND(F(imin,jmin,ke-ngh+1),1,face(ngh),destR,0,MPI_COMM_Cart,req(2),ierr)
    call MPI_IRECV(F(imin,jmin,ke+1    ),1,face(ngh),srcL ,0,MPI_COMM_Cart,req(3),ierr)
    call MPI_ISEND(F(imin,jmin,ks      ),1,face(ngh),destL,0,MPI_COMM_Cart,req(4),ierr)
!    call MPI_WAITALL(4,req,sta,ierr)
  end subroutine ighost_z
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
  subroutine ghost_xAdd(F,ir1,is1,iwork,req)
    use module_grid
    use module_hello
    use module_tmpvar
    implicit none
    include 'mpif.h'
    integer, intent(in) :: ir1, is1, iwork
    integer, intent(out) :: req(4)
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: F
    integer :: jlen, klen, ierr !,sta(MPI_STATUS_SIZE,4)
    integer, save :: srcL, srcR, destL, destR, face
    logical, save :: first_time=.true.

    if(first_time)then
      first_time=.false.
      jlen=jmax-jmin+1; klen=kmax-kmin+1; !ilen=ngh
      call para_type_block3a(imin, imax, jmin, jmax, 2, jlen, klen, MPI_DOUBLE_PRECISION, face)
      call MPI_CART_SHIFT(MPI_COMM_CART, 0, 1, srcR, destR, ierr)
      call MPI_CART_SHIFT(MPI_COMM_CART, 0,-1, srcL, destL, ierr)
    endif

    work(ir1:ir1+1,jmin:jmax,kmin:kmax,iwork) = 0d0
    work(ie-1:ie  ,jmin:jmax,kmin:kmax,iwork) = 0d0
    call MPI_IRECV(work(ir1 ,jmin,kmin,iwork),1,face,srcR ,0,MPI_Comm_Cart,req(1),ierr)
    call MPI_ISEND(F   (is1 ,jmin,kmin      ),1,face,destR,0,MPI_Comm_Cart,req(2),ierr)
    call MPI_IRECV(work(ie-1,jmin,kmin,iwork),1,face,srcL ,0,MPI_Comm_Cart,req(3),ierr)
    call MPI_ISEND(F   (is-2,jmin,kmin      ),1,face,destL,0,MPI_Comm_Cart,req(4),ierr)
!    call MPI_WAITALL(4,req,sta,ierr)

  end subroutine ghost_xAdd
!-------------------------------------------------------------------------------------------------
  subroutine ghost_yAdd(F,jr1,js1,iwork,req)
    use module_grid
    use module_hello
    use module_tmpvar
    implicit none
    include 'mpif.h'
    integer, intent(in) :: jr1, js1, iwork
    integer, intent(out) :: req(4)
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: F
    integer :: ilen, klen, ierr !,sta(MPI_STATUS_SIZE,4)
    integer, save :: srcL, srcR, destL, destR, face
    logical, save :: first_time=.true.

    if(first_time)then
      first_time=.false.
      ilen=imax-imin+1; klen=kmax-kmin+1; !ilen=ngh
      call para_type_block3a(imin, imax, jmin, jmax, ilen, 2, klen, MPI_DOUBLE_PRECISION, face)
      call MPI_CART_SHIFT(MPI_COMM_CART, 1, 1, srcR, destR, ierr)
      call MPI_CART_SHIFT(MPI_COMM_CART, 1,-1, srcL, destL, ierr)
    endif

    work(imin:imax,jr1:jr1+1,kmin:kmax,iwork) = 0d0
    work(imin:imax,je-1:je  ,kmin:kmax,iwork) = 0d0
    call MPI_IRECV(work(imin,jr1 ,kmin,iwork),1,face,srcR ,0,MPI_Comm_Cart,req(1),ierr)
    call MPI_ISEND(F   (imin,js1 ,kmin      ),1,face,destR,0,MPI_Comm_Cart,req(2),ierr)
    call MPI_IRECV(work(imin,je-1,kmin,iwork),1,face,srcL ,0,MPI_Comm_Cart,req(3),ierr)
    call MPI_ISEND(F   (imin,js-2,kmin      ),1,face,destL,0,MPI_Comm_Cart,req(4),ierr)
!    call MPI_WAITALL(4,req,sta,ierr)
  end subroutine ghost_yAdd
!-------------------------------------------------------------------------------------------------
  subroutine ghost_zAdd(F,kr1,ks1,iwork,req)
    use module_grid

      use module_hello
    use module_tmpvar
    implicit none
    include 'mpif.h'
    integer, intent(in) :: kr1, ks1, iwork
    integer, intent(out) :: req(4)
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: F
    integer :: ilen, jlen, ierr !,sta(MPI_STATUS_SIZE,4)
    integer, save :: srcL, srcR, destL, destR, face
    logical, save :: first_time=.true.

    if(first_time)then
      first_time=.false.
      ilen=imax-imin+1; jlen=jmax-jmin+1; !ilen=ngh
      call para_type_block3a(imin, imax, jmin, jmax, ilen, jlen, 2, MPI_DOUBLE_PRECISION, face)
      call MPI_CART_SHIFT(MPI_COMM_CART, 2, 1, srcR, destR, ierr)
      call MPI_CART_SHIFT(MPI_COMM_CART, 2,-1, srcL, destL, ierr)
    endif

    work(imin:imax,jmin:jmax,kr1:kr1+1,iwork) = 0d0
    work(imin:imax,jmin:jmax,ke-1:ke  ,iwork) = 0d0
    call MPI_IRECV(work(imin,jmin,kr1 ,iwork),1,face,srcR ,0,MPI_Comm_Cart,req(1),ierr)
    call MPI_ISEND(F   (imin,jmin,ks1       ),1,face,destR,0,MPI_Comm_Cart,req(2),ierr)
    call MPI_IRECV(work(imin,jmin,ke-1,iwork),1,face,srcL ,0,MPI_Comm_Cart,req(3),ierr)
    call MPI_ISEND(F   (imin,jmin,ks-2      ),1,face,destL,0,MPI_Comm_Cart,req(4),ierr)
!    call MPI_WAITALL(4,req,sta,ierr)
  end subroutine ghost_zAdd
!-------------------------------------------------------------------------------------------------
end module module_BC
!=================================================================================================
!=================================================================================================
! module_poisson:
! Using hypre library, solves linear equations with matrix A and solution P as
! A7*Pijk = A1*Pi-1jk + A2*Pi+1jk + A3*Pij-1k + 
!           A4*Pij+1k + A5*Pijk-1 + A6*Pijk+1 + A8
! 
! Syntax: call poi_initialize(mpi_comm_in,is,ie,js,je,ks,ke,Nx,Ny,Nz,bdry_cond)
!
!   Input:  integer :: mpi_comm_in               MPI communicator
!                      is,ie,js,je,ks,ke         bounds of each subdomain
!                      Nx,Ny,Nz                  total size of computational domain
!                      bdry_cond(3)              boundary conditions: use 1 for periodic
!
! Syntax: call poi_solve(A,p,maxError,maxIteration,numIteration)
!
!   Input:  real(8) :: A(is:ie,js:je,ks:ke,1:8)  coefficients and source term
!                      maxError                  criterion of convergence
!                      p(is:ie,js:je,ks:ke)      initial guess
!           integer :: maxIteration              maximum number of iterations
!
!   Output: real(8) :: p(is:ie,js:je,ks:ke)      solution
!           integer :: numIteration              number of iterations performed
!
! Syntax: call poi_finalize
!
! written by Sadegh Dabiri sdabiri@nd.edu 10/2011
!-------------------------------------------------------------------------------------------------
module module_poisson
  implicit none
  save
  private
  public :: poi_initialize, poi_solve, poi_finalize

  integer :: nstencil
  integer*8 :: grid_obj, stencil, Amat, Bvec, Xvec, solver
  integer, dimension(:), allocatable :: stencil_indices, ilower, iupper
  integer :: mpi_comm_poi,is,ie,js,je,ks,ke,Mx,My,Mz
  integer, parameter :: ndim=3
  contains
!=================================================================================================
!=================================================================================================
  subroutine poi_initialize(mpi_comm_in,iis,iie,jjs,jje,kks,kke,Nx,Ny,Nz,bdry_cond)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: mpi_comm_in,iis,iie,jjs,jje,kks,kke,Nx,Ny,Nz,bdry_cond(3)
    integer :: ierr, periodic_array(3), i
    integer, dimension(:,:), allocatable :: offsets
    
    mpi_comm_poi = mpi_comm_in
    is = iis;   ie = iie;   Mx = ie-is+1
    js = jjs;   je = jje;   My = je-js+1
    ks = kks;   ke = kke;   Mz = ke-ks+1

    nstencil = 2 * ndim + 1
    allocate(ilower(ndim), iupper(ndim), stencil_indices(nstencil) )
    ilower = (/is, js, ks/);  iupper = (/ie, je, ke/)
    periodic_array = 0
    if( bdry_cond(1)==1) periodic_array(1) = Nx
    if( bdry_cond(2)==1) periodic_array(2) = Ny
    if( bdry_cond(3)==1) periodic_array(3) = Nz

    call HYPRE_StructGridCreate(mpi_comm_poi, ndim, grid_obj, ierr)    ! create a 3d grid object
    call HYPRE_StructGridSetExtents(grid_obj, ilower, iupper, ierr)    ! add a box to the grid
    call HYPRE_StructGridSetPeriodic(grid_obj, periodic_array, ierr)   ! set periodic
    call HYPRE_StructGridAssemble(grid_obj, ierr)                      ! assemble the grid
    call HYPRE_StructStencilCreate(ndim, nstencil, stencil, ierr)      ! create a 3d 7-pt stencil

    allocate(offsets(ndim,nstencil), stat=ierr) ! define the geometry of the stencil
    offsets(:,1) = (/ 0, 0, 0/)
    offsets(:,2) = (/-1, 0, 0/)
    offsets(:,3) = (/ 1, 0, 0/)
    offsets(:,4) = (/ 0,-1, 0/)
    offsets(:,5) = (/ 0, 1, 0/)
    offsets(:,6) = (/ 0, 0,-1/)
    offsets(:,7) = (/ 0, 0, 1/)
    do i = 1, nstencil
      stencil_indices(i) = i-1
      call HYPRE_StructStencilSetElement(stencil, stencil_indices(i), offsets(:,i), ierr)
    enddo

    call HYPRE_StructMatrixCreate(mpi_comm_poi, grid_obj, stencil, Amat, ierr)! set up matrix Amat
    call HYPRE_StructMatrixInitialize(Amat, ierr)

    call HYPRE_StructVectorCreate(mpi_comm_poi, grid_obj, Bvec, ierr)! create vector object Bvec
    call HYPRE_StructVectorInitialize(Bvec, ierr)

    call HYPRE_StructVectorCreate(mpi_comm_poi, grid_obj, Xvec, ierr)! create vector object Xvec
    call HYPRE_StructVectorInitialize(Xvec, ierr)

  end subroutine poi_initialize
!=================================================================================================
! Solve Poisson equation. p is initial guess. 
!=================================================================================================
  subroutine poi_solve(A,p,maxError,maxit,num_iterations)
    implicit none
    include 'mpif.h'
    integer :: ierr, nvalues, ijk, i,j,k
    real(8), dimension(:), allocatable :: values
    real(8), dimension(is:ie,js:je,ks:ke,8), intent(in) :: A
    real(8), dimension(is:ie,js:je,ks:ke), intent(inout) :: p
!    real(8), dimension(is-2:ie+2,js-2:je+2,ks-2:ke+2), intent(inout) :: p
    real(8), intent(in) :: maxError
    integer, intent(in) :: maxit
    integer, intent(out):: num_iterations
    real(8) :: final_res_norm
!----------------------------------------Fill in matrix Amat--------------------------------------
    num_iterations = 0
    nvalues = mx * my * mz * nstencil
    allocate(values(nvalues), stat=ierr)
    if(ierr/=0)stop '**** poi_solve: allocation error ****'

    ijk = 1
    do k=ks,ke;  do j=js,je;  do i=is,ie
      values(ijk+1) = -A(i,j,k,1)
      values(ijk+2) = -A(i,j,k,2)
      values(ijk+3) = -A(i,j,k,3)
      values(ijk+4) = -A(i,j,k,4)
      values(ijk+5) = -A(i,j,k,5)
      values(ijk+6) = -A(i,j,k,6)
      values(ijk  ) =  A(i,j,k,7)
      ijk = ijk + 7
    enddo;  enddo;  enddo
    call HYPRE_StructMatrixSetBoxValues(Amat, ilower, iupper, nstencil, stencil_indices, &
                                        values, ierr)
    deallocate(values, stat=ierr)
!------------------------------Fill in source term and initial guess------------------------------
    nvalues = mx * my * mz
    allocate(values(nvalues), stat=ierr)

    ijk = 1
    do k=ks,ke;  do j=js,je;  do i=is,ie
      values(ijk) = A(i,j,k,8)
      ijk = ijk + 1
    enddo;  enddo;  enddo
    call HYPRE_StructVectorSetBoxValues(Bvec, ilower, iupper, values, ierr)

    ijk = 1
    do k=ks,ke;  do j=js,je;  do i=is,ie
      values(ijk) = p(i,j,k)
      ijk = ijk + 1
    enddo;  enddo;  enddo
    call HYPRE_StructVectorSetBoxValues(Xvec, ilower, iupper, values, ierr)
    deallocate(values, stat=ierr)
!---------------------------Assemble matrix Amat and vectors Bvec,Xvec----------------------------
    call HYPRE_StructMatrixAssemble(Amat, ierr)
    call HYPRE_StructVectorAssemble(Bvec, ierr)
    call HYPRE_StructVectorAssemble(Xvec, ierr)
!---------------------------------------Solve the equations---------------------------------------
    call HYPRE_StructPFMGCreate(mpi_comm_poi, solver, ierr)
    call HYPRE_StructPFMGSetMaxIter(solver, maxit, ierr)
    call HYPRE_StructPFMGSetTol(solver, MaxError, ierr)

    call HYPRE_StructPFMGSetRelaxType(solver, 3, ierr)

    call hypre_structPFMGsetLogging(solver, 1, ierr)
    call HYPRE_StructPFMGSetup(solver, Amat, Bvec, Xvec, ierr)
    call HYPRE_StructPFMGSolve(solver, Amat, Bvec, Xvec, ierr)
    !call hypre_structPFMGgetnumiterations(solver, num_iterations, ierr)
    !call hypre_structPFMGgetfinalrelative(solver, final_res_norm, ierr)
!--------------------------------------Retrieve the solution--------------------------------------
    nvalues = mx * my * mz
    allocate(values(nvalues), stat=ierr)
    call HYPRE_StructVectorGetBoxValues(Xvec, ilower, iupper, values, ierr)
    call HYPRE_StructPFMGDestroy(solver, ierr)

    ijk = 1
    do k=ks,ke;  do j=js,je;  do i=is,ie
      p(i,j,k) = values(ijk)
      ijk = ijk + 1
    enddo;  enddo;  enddo
    deallocate(values, stat=ierr)
  end subroutine poi_solve
!=================================================================================================
!=================================================================================================
  subroutine poi_finalize
    implicit none
    include 'mpif.h'
    integer :: ierr
    call HYPRE_StructGridDestroy(grid_obj, ierr)
    call HYPRE_StructStencilDestroy(stencil, ierr)
    call HYPRE_StructMatrixDestroy(Amat, ierr)
    call HYPRE_StructVectorDestroy(Bvec, ierr)
    call HYPRE_StructVectorDestroy(Xvec, ierr)
    deallocate(ilower, iupper, stencil_indices, stat=ierr)
  end subroutine poi_finalize
!-------------------------------------------------------------------------------------------------
end module module_poisson
!=================================================================================================
!=================================================================================================
! creates new data type for passing non-contiguous 
!-------------------------------------------------------------------------------------------------
SUBROUTINE para_type_block3a(imin, imax, jmin, jmax, ilen, jlen, klen, ioldtype,inewtype)
  implicit none
  INCLUDE 'mpif.h'
  integer :: imin, imax, jmin, jmax, ilen, jlen, klen, ioldtype,inewtype,isize, ierr, itemp, idist
  CALL MPI_TYPE_EXTENT(ioldtype, isize, ierr)
  CALL MPI_TYPE_VECTOR(jlen, ilen, imax - imin + 1, ioldtype, itemp, ierr)
  idist = (imax - imin + 1) * (jmax - jmin + 1) * isize
  CALL MPI_TYPE_HVECTOR(klen, 1, idist, itemp, inewtype, ierr)
  CALL MPI_TYPE_COMMIT(inewtype, ierr)
END
!=================================================================================================
!=================================================================================================
! The Poisson equation for the density is setup with matrix A as
! A7*Pijk = A1*Pi-1jk + A2*Pi+1jk + A3*Pij-1k + 
!           A4*Pij+1k + A5*Pijk-1 + A6*Pijk+1 + A8
!-------------------------------------------------------------------------------------------------
subroutine SetupDensity(dIdx,dIdy,dIdz,A,color) !,mask)
  use module_grid
  use module_BC
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: dIdx,dIdy,dIdz, color
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(out) :: A
  integer :: i,j,k
  logical, save :: first=.true.
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax) :: dI
  do k=ks,ke; do j=js,je; do i=is,ie
      A(i,j,k,1) = 1d0/dx(i)/dxh(i-1)
      A(i,j,k,2) = 1d0/dx(i)/dxh(i  )
      A(i,j,k,3) = 1d0/dy(j)/dyh(j-1)
      A(i,j,k,4) = 1d0/dy(j)/dyh(j  )
      A(i,j,k,5) = 1d0/dz(k)/dzh(k-1)
      A(i,j,k,6) = 1d0/dz(k)/dzh(k  )
      A(i,j,k,7) = sum(A(i,j,k,1:6))
      A(i,j,k,8) = -(dIdx(i,j,k)-dIdx(i-1,j,k))/dx(i) &
                   -(dIdy(i,j,k)-dIdy(i,j-1,k))/dy(j) &
                   -(dIdz(i,j,k)-dIdz(i,j,k-1))/dz(k)
  enddo; enddo; enddo

  if(bdry_cond(1)==0)then
    if(coords(1)==0    ) A(is,:,:,8)=A(is,:,:,8)+A(is,:,:,1)
    if(coords(1)==nPx-1) A(ie,:,:,8)=A(ie,:,:,8)+A(ie,:,:,2)
    if(coords(1)==0    ) A(is,:,:,1) = 0d0
    if(coords(1)==nPx-1) A(ie,:,:,2) = 0d0
  endif
  if(bdry_cond(2)==0)then
    if(coords(2)==0    ) A(:,js,:,8)=A(:,js,:,8)+A(:,js,:,3)
    if(coords(2)==nPy-1) A(:,je,:,8)=A(:,je,:,8)+A(:,je,:,4)
    if(coords(2)==0    ) A(:,js,:,3) = 0d0
    if(coords(2)==nPy-1) A(:,je,:,4) = 0d0
  endif
  if(bdry_cond(3)==0)then
    if(coords(3)==0    ) A(:,:,ks,8)=A(:,:,ks,8)+A(:,:,ks,5)
    if(coords(3)==nPz-1) A(:,:,ke,8)=A(:,:,ke,8)+A(:,:,ke,6)
    if(coords(3)==0    ) A(:,:,ks,5) = 0d0
    if(coords(3)==nPz-1) A(:,:,ke,6) = 0d0
  endif

  if(.not. first)then
  	first=.false.
    do k=ks,ke; do j=js,je; do i=is,ie
      dI(i,j,k) = ( (dIdx(i,j,k)+dIdx(i-1,j,k))**2 + &
                    (dIdx(i,j,k)+dIdx(i-1,j,k))**2 + &
                    (dIdx(i,j,k)+dIdx(i-1,j,k))**2 )
      if( dI(i,j,k)<1e-6 )then
        A(i,j,k,1:6) = 0d0
        A(i,j,k,7) = 1d0
        A(i,j,k,8) = float(floor(color(i,j,k)+0.5))
      endif
    enddo; enddo; enddo
 endif

!  do k=ks,ke; do j=js,je; do i=is,ie
!      A(i,j,k,7) = sum(A(i,j,k,1:6))
!  enddo; enddo; enddo
  ! Anchor a point to 1
!  if(coords(1)==0 .and. coords(2)==0 .and. coords(3)==0) then
!    A(is,js,ks,:)=0d0
!    A(is,js,ks,7:8)=1d0
!  endif
end subroutine SetupDensity
!=================================================================================================
!=================================================================================================
! The Poisson equation for the pressure is setup with matrix A as
! A7*Pijk = A1*Pi-1jk + A2*Pi+1jk + A3*Pij-1k + 
!           A4*Pij+1k + A5*Pijk-1 + A6*Pijk+1 + A8
!-------------------------------------------------------------------------------------------------
subroutine SetupPoisson(utmp,vtmp,wtmp,umask,vmask,wmask,rhot,dt,A) !,mask)
  use module_grid
  use module_hello
  use module_BC
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: utmp,vtmp,wtmp,rhot
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: umask,vmask,wmask
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(out) :: A
  real(8) :: dt
  integer :: i,j,k
  do k=ks,ke; do j=js,je; do i=is,ie;
!    if(mask(i,j,k))then
      A(i,j,k,1) = 2d0*dt*umask(i-1,j,k)/(dx(i)*dxh(i-1)*(rhot(i-1,j,k)+rhot(i,j,k)))
      A(i,j,k,2) = 2d0*dt*umask(i,j,k)/(dx(i)*dxh(i  )*(rhot(i+1,j,k)+rhot(i,j,k)))
      A(i,j,k,3) = 2d0*dt*vmask(i,j-1,k)/(dy(j)*dyh(j-1)*(rhot(i,j-1,k)+rhot(i,j,k)))
      A(i,j,k,4) = 2d0*dt*vmask(i,j,k)/(dy(j)*dyh(j  )*(rhot(i,j+1,k)+rhot(i,j,k)))
      A(i,j,k,5) = 2d0*dt*wmask(i,j,k-1)/(dz(k)*dzh(k-1)*(rhot(i,j,k-1)+rhot(i,j,k)))
      A(i,j,k,6) = 2d0*dt*wmask(i,j,k)/(dz(k)*dzh(k  )*(rhot(i,j,k+1)+rhot(i,j,k)))
      A(i,j,k,7) = sum(A(i,j,k,1:6)) + 1d-49
      A(i,j,k,8) = -( (utmp(i,j,k)-utmp(i-1,j,k))/dx(i) &
                     +(vtmp(i,j,k)-vtmp(i,j-1,k))/dy(j) &
                     +(wtmp(i,j,k)-wtmp(i,j,k-1))/dz(k) )
!    endif
  enddo; enddo; enddo
end subroutine SetupPoisson
!=================================================================================================
!=================================================================================================
! The equation for the U velocity is setup with matrix A as
! A7*Uijk = A1*Ui-1jk + A2*Ui+1jk + A3*Uij-1k + 
!           A4*Uij+1k + A5*Uijk-1 + A6*Uijk+1 + A8
!-------------------------------------------------------------------------------------------------
subroutine SetupUvel(u,du,rho,mu,rho1,mu1,dt,A)
  use module_grid
  use module_hello
  use module_BC

  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: u,du,rho,mu
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(out) :: A
  real(8), intent(in) :: dt,rho1,mu1
  real(8) :: rhom
  integer :: i,j,k
  if(TwoPhase) then
     do k=ks,ke; do j=js,je; do i=is,ie;
        rhom = 0.5d0*(rho(i+1,j,k)+rho(i,j,k))
        A(i,j,k,1) = dt/(dx(i  )*dxh(i)*rhom)*2d0*mu(i  ,j,k)
        A(i,j,k,2) = dt/(dx(i+1)*dxh(i)*rhom)*2d0*mu(i+1,j,k)
        A(i,j,k,3) = dt/(dy(j)*dyh(j-1)*rhom)*0.25d0*(mu(i,j,k)+mu(i+1,j,k)+mu(i,j-1,k)+mu(i+1,j-1,k))
        A(i,j,k,4) = dt/(dy(j)*dyh(j  )*rhom)*0.25d0*(mu(i,j,k)+mu(i+1,j,k)+mu(i,j+1,k)+mu(i+1,j+1,k))
        A(i,j,k,5) = dt/(dz(k)*dzh(k-1)*rhom)*0.25d0*(mu(i,j,k)+mu(i+1,j,k)+mu(i,j,k-1)+mu(i+1,j,k-1))
        A(i,j,k,6) = dt/(dz(k)*dzh(k  )*rhom)*0.25d0*(mu(i,j,k)+mu(i+1,j,k)+mu(i,j,k+1)+mu(i+1,j,k+1))
        A(i,j,k,7) = 1d0+sum(A(i,j,k,1:6))
        A(i,j,k,8) = u(i,j,k) + dt*du(i,j,k)
     enddo; enddo; enddo
  else
     do k=ks,ke; do j=js,je; do i=is,ie;
        rhom = rho1
        A(i,j,k,1) = dt/(dx(i  )*dxh(i)*rhom)*mu1
        A(i,j,k,2) = dt/(dx(i+1)*dxh(i)*rhom)*mu1
        A(i,j,k,3) = dt/(dy(j)*dyh(j-1)*rhom)*mu1
        A(i,j,k,4) = dt/(dy(j)*dyh(j  )*rhom)*mu1
        A(i,j,k,5) = dt/(dz(k)*dzh(k-1)*rhom)*mu1
        A(i,j,k,6) = dt/(dz(k)*dzh(k  )*rhom)*mu1
        A(i,j,k,7) = 1d0+sum(A(i,j,k,1:6))
        A(i,j,k,8) = u(i,j,k) + dt*du(i,j,k)
     enddo; enddo; enddo
  endif

!-------------------------------------------------------------------------------------------------
  !wall boundary conditions
  if(bdry_cond(2)==0 .and. coords(2)==0    ) then
    do k=ks,ke; do i=is,ie;
      A(i,js,k,7) = A(i,js,k,7) + A(i,js,k,3)
      A(i,js,k,3) = 0d0
    enddo; enddo
  endif
  if(bdry_cond(5)==0 .and. coords(2)==nPy-1) then
    do k=ks,ke; do i=is,ie;
      A(i,je,k,7) = A(i,je,k,7) + A(i,je,k,4)
      A(i,je,k,4) = 0d0
    enddo; enddo
  endif

  if(bdry_cond(3)==0 .and. coords(3)==0    ) then
    do j=js,je; do i=is,ie;
      A(i,j,ks,7) = A(i,j,ks,7) + A(i,j,ks,5)
      A(i,j,ks,5) = 0d0
    enddo; enddo
  endif
  if(bdry_cond(6)==0 .and. coords(3)==nPz-1) then
    do j=js,je; do i=is,ie;
      A(i,j,ke,7) = A(i,j,ke,7) + A(i,j,ke,6)
      A(i,j,ke,6) = 0d0
    enddo; enddo
  endif
end subroutine SetupUvel
!=================================================================================================
!=================================================================================================
subroutine SetupVvel(v,dv,rho,mu,rho1,mu1,dt,A) 
  use module_grid
  use module_hello
  use module_BC
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: v,dv,rho,mu
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(out) :: A
  real(8), intent(in) :: dt,rho1,mu1
  real(8) :: rhom
  integer :: i,j,k
  if(TwoPhase) then
     do k=ks,ke; do j=js,je; do i=is,ie;
        rhom = 0.5d0*(rho(i,j+1,k)+rho(i,j,k))
        if((dy(j  )*dyh(j)*rhom)==0d0)then
           write(*,'(10I4,3f7.3)'),rank,i,j,k,is,ie,js,je,ks,ke,dy(j),dyh(j),rhom
        endif
        A(i,j,k,3) = dt/(dy(j  )*dyh(j)*rhom)*2d0*mu(i,j  ,k)
        A(i,j,k,4) = dt/(dy(j+1)*dyh(j)*rhom)*2d0*mu(i,j+1,k)
        A(i,j,k,5) = dt/(dz(k)*dzh(k-1)*rhom)*0.25d0*(mu(i,j,k)+mu(i,j+1,k)+mu(i,j,k-1)+mu(i,j+1,k-1))
        A(i,j,k,6) = dt/(dz(k)*dzh(k  )*rhom)*0.25d0*(mu(i,j,k)+mu(i,j+1,k)+mu(i,j,k+1)+mu(i,j+1,k+1))
        A(i,j,k,1) = dt/(dx(i)*dxh(i-1)*rhom)*0.25d0*(mu(i,j,k)+mu(i,j+1,k)+mu(i-1,j,k)+mu(i-1,j+1,k))
        A(i,j,k,2) = dt/(dx(i)*dxh(i  )*rhom)*0.25d0*(mu(i,j,k)+mu(i,j+1,k)+mu(i+1,j,k)+mu(i+1,j+1,k))
        A(i,j,k,7) = 1d0+sum(A(i,j,k,1:6))
        A(i,j,k,8) = v(i,j,k) + dt*dv(i,j,k)
     enddo; enddo; enddo
  else
     do k=ks,ke; do j=js,je; do i=is,ie;
        rhom = rho1
        A(i,j,k,3) = dt/(dy(j  )*dyh(j)*rhom)*mu1
        A(i,j,k,4) = dt/(dy(j+1)*dyh(j)*rhom)*mu1
        A(i,j,k,5) = dt/(dz(k)*dzh(k-1)*rhom)*mu1
        A(i,j,k,6) = dt/(dz(k)*dzh(k  )*rhom)*mu1
        A(i,j,k,1) = dt/(dx(i-1)*dxh(i)*rhom)*mu1
        A(i,j,k,2) = dt/(dx(i)*dxh(i)*rhom)*mu1
        A(i,j,k,7) = 1d0+sum(A(i,j,k,1:6))
        A(i,j,k,8) = v(i,j,k) + dt*dv(i,j,k)
     enddo; enddo; enddo
  endif
  !-------------------------------------------------------------------------------------------------
  !wall boundary conditions
  if(bdry_cond(1)==0 .and. coords(1)==0    ) then
     do k=ks,ke; do j=js,je;
        A(is,j,k,7) = A(is,j,k,7) + A(is,j,k,1)
      A(is,j,k,1) = 0d0
    enddo; enddo
  endif
  if(bdry_cond(4)==0 .and. coords(1)==nPx-1) then
    do k=ks,ke; do j=js,je;
      A(ie,j,k,7) = A(ie,j,k,7) + A(ie,j,k,2)
      A(ie,j,k,2) = 0d0
    enddo; enddo
  endif

  if(bdry_cond(3)==0 .and. coords(3)==0    ) then
    do j=js,je; do i=is,ie;
      A(i,j,ks,7) = A(i,j,ks,7) + A(i,j,ks,5)
      A(i,j,ks,5) = 0d0
    enddo; enddo
  endif
  if(bdry_cond(6)==0 .and. coords(3)==nPz-1) then
    do j=js,je; do i=is,ie;
      A(i,j,ke,7) = A(i,j,ke,7) + A(i,j,ke,6)
      A(i,j,ke,6) = 0d0
    enddo; enddo
  endif
end subroutine SetupVvel
!=================================================================================================
!=================================================================================================
subroutine SetupWvel(w,dw,rho,mu,rho1,mu1,dt,A)
  use module_grid
  use module_hello
  use module_BC
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: w,dw,rho,mu
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(out) :: A
  real(8), intent(in) :: dt,rho1,mu1
  real(8) :: rhom
  integer :: i,j,k
  if(TwoPhase) then
     do k=ks,ke; do j=js,je; do i=is,ie;
        rhom = 0.5d0*(rho(i,j,k+1)+rho(i,j,k))
        A(i,j,k,5) = dt/(dz(k  )*dzh(k)*rhom)*2d0*mu(i,j,k  )
        A(i,j,k,6) = dt/(dz(k+1)*dzh(k)*rhom)*2d0*mu(i,j,k+1)
        A(i,j,k,1) = dt/(dx(i)*dxh(i-1)*rhom)*0.25d0*(mu(i,j,k)+mu(i,j,k+1)+mu(i-1,j,k)+mu(i-1,j,k+1))
        A(i,j,k,2) = dt/(dx(i)*dxh(i  )*rhom)*0.25d0*(mu(i,j,k)+mu(i,j,k+1)+mu(i+1,j,k)+mu(i+1,j,k+1))
        A(i,j,k,3) = dt/(dy(j)*dyh(j-1)*rhom)*0.25d0*(mu(i,j,k)+mu(i,j,k+1)+mu(i,j-1,k)+mu(i,j-1,k+1))
        A(i,j,k,4) = dt/(dy(j)*dyh(j  )*rhom)*0.25d0*(mu(i,j,k)+mu(i,j,k+1)+mu(i,j+1,k)+mu(i,j+1,k+1))
        A(i,j,k,7) = 1d0+sum(A(i,j,k,1:6))
        A(i,j,k,8) = w(i,j,k) + dt*dw(i,j,k)
     enddo; enddo; enddo
  else
      do k=ks,ke; do j=js,je; do i=is,ie;
        rhom = rho1
        A(i,j,k,5) = dt/(dz(k  )*dzh(k)*rhom)*mu1
        A(i,j,k,6) = dt/(dz(k+1)*dzh(k)*rhom)*mu1
        A(i,j,k,1) = dt/(dx(i)*dxh(i-1)*rhom)*mu1
        A(i,j,k,2) = dt/(dx(i)*dxh(i  )*rhom)*mu1
        A(i,j,k,3) = dt/(dy(j)*dyh(j-1)*rhom)*mu1
        A(i,j,k,4) = dt/(dy(j)*dyh(j  )*rhom)*mu1
        A(i,j,k,7) = 1d0+sum(A(i,j,k,1:6))
        A(i,j,k,8) = w(i,j,k) + dt*dw(i,j,k)
     enddo; enddo; enddo
  endif
  !-------------------------------------------------------------------------------------------------
  !wall boundary conditions
  if(bdry_cond(1)==0 .and. coords(1)==0    ) then
     do k=ks,ke; do j=js,je;
        A(is,j,k,7) = A(is,j,k,7) + A(is,j,k,1)
        A(is,j,k,1) = 0d0
     enddo; enddo
  endif
  if(bdry_cond(4)==0 .and. coords(1)==nPx-1) then
    do k=ks,ke; do j=js,je;
      A(ie,j,k,7) = A(ie,j,k,7) + A(ie,j,k,2)
      A(ie,j,k,2) = 0d0
    enddo; enddo
  endif

  if(bdry_cond(2)==0 .and. coords(2)==0    ) then
    do k=ks,ke; do i=is,ie;
      A(i,js,k,7) = A(i,js,k,7) + A(i,js,k,3)
      A(i,js,k,3) = 0d0
    enddo; enddo
  endif
  if(bdry_cond(5)==0 .and. coords(2)==nPy-1) then
    do k=ks,ke; do i=is,ie;
      A(i,je,k,7) = A(i,je,k,7) + A(i,je,k,4)
      A(i,je,k,4) = 0d0
    enddo; enddo
  endif
end subroutine SetupWvel
!=================================================================================================
!=================================================================================================
! Solves the following linear equiation:
! A7*Pijk = A1*Pi-1jk + A2*Pi+1jk + A3*Pij-1k + 
!           A4*Pij+1k + A5*Pijk-1 + A6*Pijk+1 + A8
!-------------------------------------------------------------------------------------------------
subroutine LinearSolver(A,p,maxError,beta,maxit,it,ierr)
  use module_grid
  use module_hello
  use module_BC
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: p
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(in) :: A
  real(8), intent(in) :: beta, maxError
  integer, intent(in) :: maxit
  integer, intent(out) :: it, ierr
  real(8) :: res, totalres
  integer :: i,j,k
  integer :: req(12),sta(MPI_STATUS_SIZE,12)
  logical :: mask(imin:imax,jmin:jmax,kmin:kmax)
!--------------------------------------ITERATION LOOP--------------------------------------------  
  do it=1,maxit
    do k=ks,ke; do j=js,je; do i=is,ie
      p(i,j,k)=(1d0-beta)*p(i,j,k)+beta* 1d0/A(i,j,k,7)*(              &
        A(i,j,k,1) * p(i-1,j,k) + A(i,j,k,2) * p(i+1,j,k) +            &
        A(i,j,k,3) * p(i,j-1,k) + A(i,j,k,4) * p(i,j+1,k) +            &
        A(i,j,k,5) * p(i,j,k-1) + A(i,j,k,6) * p(i,j,k+1) + A(i,j,k,8))
    enddo; enddo; enddo
!---------------------------------CHECK FOR CONVERGENCE-------------------------------------------
    res = 0d0
    call ghost_x(p,1,req( 1: 4)); call ghost_y(p,1,req( 5: 8)); call ghost_z(p,1,req( 9:12))
    do k=ks+1,ke-1; do j=js+1,je-1; do i=is+1,ie-1
      res=res+abs(-p(i,j,k) * A(i,j,k,7) +                             &
        A(i,j,k,1) * p(i-1,j,k) + A(i,j,k,2) * p(i+1,j,k) +            &
        A(i,j,k,3) * p(i,j-1,k) + A(i,j,k,4) * p(i,j+1,k) +            &
        A(i,j,k,5) * p(i,j,k-1) + A(i,j,k,6) * p(i,j,k+1) + A(i,j,k,8) )
    enddo; enddo; enddo
    call MPI_WAITALL(12,req,sta,ierr)
    mask=.true.
    mask(is+1:ie-1,js+1:je-1,ks+1:ke-1)=.false.
    do k=ks,ke; do j=js,je; do i=is,ie
      if(mask(i,j,k))res=res+abs(-p(i,j,k) * A(i,j,k,7) +              &
        A(i,j,k,1) * p(i-1,j,k) + A(i,j,k,2) * p(i+1,j,k) +            &
        A(i,j,k,3) * p(i,j-1,k) + A(i,j,k,4) * p(i,j+1,k) +            &
        A(i,j,k,5) * p(i,j,k-1) + A(i,j,k,6) * p(i,j,k+1) + A(i,j,k,8) )
    enddo; enddo; enddo
    res = res/float(Nx*Ny*Nz)
    call MPI_ALLREDUCE(res, totalres, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr)
    if (.not.(totalres<1e10)) then
      ierr=1 !stop '***** solution has diverged *****'
      if(rank==0) print*,'Solver 0 diverged after',it,'iterations.'
      return
    endif
    if (totalres<maxError) exit
  enddo
  if(it==maxit+1 .and. rank==0) write(*,*) 'Warning: LinearSolver reached maxit.'
end subroutine LinearSolver
!=================================================================================================
!=================================================================================================
! Solves the following linear equiation:
! A7*Uijk = umask*(A1*Ui-1jk + A2*Ui+1jk + A3*Uij-1k + 
!           A4*Uij+1k + A5*Uijk-1 + A6*Uijk+1 + A8)
!-------------------------------------------------------------------------------------------------
subroutine LinearSolver1(A,u,umask,maxError,beta,maxit,it,ierr)
  use module_grid
  use module_hello
  use module_BC
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: u
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: umask
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(in) :: A
  real(8), intent(in) :: beta, maxError
  integer, intent(in) :: maxit
  integer, intent(out) :: it, ierr
  real(8) :: res, totalres
  integer :: i,j,k
  integer :: req(12),sta(MPI_STATUS_SIZE,12)
  logical :: mask(imin:imax,jmin:jmax,kmin:kmax)
!--------------------------------------ITERATION LOOP--------------------------------------------  
  do it=1,maxit
    do k=ks,ke; do j=js,je; do i=is,ie
      u(i,j,k)=umask(i,j,k)*((1d0-beta)*u(i,j,k)+beta* 1d0/A(i,j,k,7)*(              &
        A(i,j,k,1) * u(i-1,j,k) + A(i,j,k,2) * u(i+1,j,k) +            &
        A(i,j,k,3) * u(i,j-1,k) + A(i,j,k,4) * u(i,j+1,k) +            &
        A(i,j,k,5) * u(i,j,k-1) + A(i,j,k,6) * u(i,j,k+1) + A(i,j,k,8)))
    enddo; enddo; enddo
!---------------------------------CHECK FOR CONVERGENCE-------------------------------------------
    res = 0d0
    call ghost_x(u,1,req( 1: 4)); call ghost_y(u,1,req( 5: 8)); call ghost_z(u,1,req( 9:12))
    do k=ks+1,ke-1; do j=js+1,je-1; do i=is+1,ie-1
      res=res+ umask(i,j,k)*abs(-u(i,j,k) * A(i,j,k,7) +                             &
        A(i,j,k,1) * u(i-1,j,k) + A(i,j,k,2) * u(i+1,j,k) +            &
        A(i,j,k,3) * u(i,j-1,k) + A(i,j,k,4) * u(i,j+1,k) +            &
        A(i,j,k,5) * u(i,j,k-1) + A(i,j,k,6) * u(i,j,k+1) + A(i,j,k,8) )
    enddo; enddo; enddo
    call MPI_WAITALL(12,req,sta,ierr)
    mask=.true.
    mask(is+1:ie-1,js+1:je-1,ks+1:ke-1)=.false.
    do k=ks,ke; do j=js,je; do i=is,ie
      if(mask(i,j,k))res=res+umask(i,j,k)*abs(-u(i,j,k) * A(i,j,k,7) +              &
        A(i,j,k,1) * u(i-1,j,k) + A(i,j,k,2) * u(i+1,j,k) +            &
        A(i,j,k,3) * u(i,j-1,k) + A(i,j,k,4) * u(i,j+1,k) +            &
        A(i,j,k,5) * u(i,j,k-1) + A(i,j,k,6) * u(i,j,k+1) + A(i,j,k,8) )
    enddo; enddo; enddo
    res = res/float(Nx*Ny*Nz)
    call MPI_ALLREDUCE(res, totalres, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr)
    if (.not.(totalres<1e10)) then
      ierr=1 !stop '***** solution has diverged *****'
      if(rank==0) print*,'Solver 1 diverged after',it,'iterations.'
      return
    endif
    if (totalres<maxError) exit
  enddo
  if(it==maxit+1 .and. rank==0) write(*,*) 'Warning: LinearSolver reached maxit.'
end subroutine LinearSolver1
!=================================================================================================
!=================================================================================================
! Returns the residual
!-------------------------------------------------------------------------------------------------
subroutine calcResidual(A,p, Residual)
  use module_grid
  use module_BC
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: p
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(in) :: A
  real(8) :: res, totalres, Residual
  integer :: i,j,k, ierr
  res = 0d0
  do k=ks,ke; do j=js,je; do i=is,ie
    res=res+abs(-p(i,j,k) * A(i,j,k,7) +                             &
      A(i,j,k,1) * p(i-1,j,k) + A(i,j,k,2) * p(i+1,j,k) +            &
      A(i,j,k,3) * p(i,j-1,k) + A(i,j,k,4) * p(i,j+1,k) +            &
      A(i,j,k,5) * p(i,j,k-1) + A(i,j,k,6) * p(i,j,k+1) + A(i,j,k,8) )**2
  enddo; enddo; enddo
  res = res/float(Nx*Ny*Nz)
  call MPI_ALLREDUCE(res, totalres, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr)
  Residual = sqrt(totalres)
end subroutine calcResidual
!=================================================================================================
!=================================================================================================
! Returns the residual
!-------------------------------------------------------------------------------------------------
subroutine calcResidual1(A,p,pmask,Residual)
  use module_grid
  use module_BC
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: p
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: pmask
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(in) :: A
  real(8) :: res, totalres, Residual
  integer :: i,j,k, ierr
  res = 0d0
  do k=ks,ke; do j=js,je; do i=is,ie
    res=res+pmask(i,j,k)*abs(-p(i,j,k) * A(i,j,k,7) +                             &
      A(i,j,k,1) * p(i-1,j,k) + A(i,j,k,2) * p(i+1,j,k) +            &
      A(i,j,k,3) * p(i,j-1,k) + A(i,j,k,4) * p(i,j+1,k) +            &
      A(i,j,k,5) * p(i,j,k-1) + A(i,j,k,6) * p(i,j,k+1) + A(i,j,k,8) )**2
  enddo; enddo; enddo
  res = res/float(Nx*Ny*Nz)
  call MPI_ALLREDUCE(res, totalres, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr)
  Residual = sqrt(totalres)
end subroutine calcResidual1
!=================================================================================================
!=================================================================================================
! Returns the flow rate
!-------------------------------------------------------------------------------------------------
subroutine calcsum(p)
  use module_grid
  use module_BC
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: p
  real(8) :: flow
  call calcsum_shift(p,flow,0,0,0)
end subroutine calcsum

subroutine calcsum2(p,flowrate)
  use module_grid
  use module_BC
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: p
  real(8) :: flow, totalflow, flowrate
  integer :: i,j,k, ierr
  call calcsum_shift(p,flow,1,0,0)
  call calcsum_shift(p,flow,0,1,0)
  call calcsum_shift(p,flow,0,0,1)
  call calcsum_shift(p,flow,-1,0,0)
  call calcsum_shift(p,flow,0,-1,0)
  call calcsum_shift(p,flow,0,0,-1)
  call calcsum_shift(p,flow,0,0,0)
  flowrate = -1.
end subroutine calcsum2

subroutine calcsum_shift(p,porosity,sx,sy,sz)
  use module_grid
  use module_BC
  use module_IO
  implicit none
  include 'mpif.h'
  integer, intent(in) :: sx,sy,sz
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: p
  real(8) :: flow, totalflow
  real(8), intent(out) :: porosity
  real(8) volume
  integer :: i,j,k, ierr
  flow=0.d0
  volume=Nx*Ny*Nz
  do k=ks,ke; do j=js,je; do i=is,ie;
     flow=flow+p(i+sx,j+sy,k+sz)
  enddo;enddo;enddo
  call MPI_REDUCE(flow, totalflow, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_Domain, ierr)
  porosity=totalflow/volume
  if(rank==0) then 
     open(UNIT=89,file=trim(out_path)//'/porosity.txt')
     write(89,310)  porosity
     close(unit=89)
  endif
  310 format(F17.11)
end subroutine calcsum_shift
!-------------------------------------------------------------------------------------------------

function calc_imax(f)
  use module_grid
  integer calc_imax
  integer, dimension(imin:imax,jmin:jmax,kmin:kmax),  intent(in) :: f
  integer :: i,j,k
  calc_imax=-2147483647
  do k=kmin,kmax
     do j=jmin,jmax
        do i=imin,imax
           calc_imax = max(calc_imax,f(i,j,k))
        enddo
     enddo
  enddo
  end function calc_imax
 
  subroutine THRESHOLD(z)
    implicit none
    real(8),intent(inout) :: z
    if ( z < 0.d0 ) then
       z = 0.d0
    else if ( z > 1.d0 ) then
       z = 1.d0
    end if
  end subroutine THRESHOLD
