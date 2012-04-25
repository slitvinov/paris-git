!=================================================================================================
!=================================================================================================
! Code: FTC3D2011 (Front Tracking Code for 3D simulations)
! Authors: Sadegh Dabiri, Gretar Tryggvason
! Contact: sdabiri@gmail.com
! A three dimensional Navier-Stokes flow solver for modeling of multiphase 
! flows. Flow can be driven by wall motion, density difference or pressure gradient.
! Boundary conditions supported: wall and periodic
! Version 1.0   1/21/2011   The 3D flow solver for variable density/viscosity is written. 
!                           The desity is advected by ENO scheme.
! Version 2.0   2/25/2011   Parallel implementation.
!=================================================================================================
!=================================================================================================
! module_grid: Contains definition of variables for the grid.
!-------------------------------------------------------------------------------------------------
module module_grid
  implicit none
  integer :: Nx, Ny, Nz, Ng ! Ng isnumber of ghost cells
  integer :: Nxt, Nyt, Nzt ! total number of cells
  integer :: is, ie, js, je, ks, ke
  integer :: ieu, jev, kew
  real(8), dimension(:), allocatable :: x, xh, dx, dxh
  real(8), dimension(:), allocatable :: y, yh, dy, dyh
  real(8), dimension(:), allocatable :: z, zh, dz, dzh
  real(8) :: xLength, yLength, zLength, xform, yform, zform

  integer :: nPx, nPy, nPz, Mx, My, Mz, rank, ndim=3, numProcess
  integer, dimension(:), allocatable :: dims, coords, periodic, reorder
  integer :: MPI_Comm_Cart
  integer :: imin, imax, jmin, jmax, kmin, kmax
  logical :: hypre
end module module_grid
!=================================================================================================
!=================================================================================================
! module_flow: Contains definition of variables for the flow solver.
!-------------------------------------------------------------------------------------------------
module module_flow
  implicit none
  real(8), dimension(:,:,:), allocatable :: u, v, w, uold, vold, wold, fx, fy, fz
  real(8), dimension(:,:,:), allocatable :: p, rho, rhoo, muold, mu
  real(8), dimension(:,:,:), allocatable :: du,dv,dw,drho
  real(8) gx, gy, gz, mu1, mu2, r_avg, dt, dtFlag
  real(8) max_velocity, maxTime, Time, EndTime, MaxDt, CFL, myfr, flowrate
  logical :: TwoPhase
  real(8) :: rho1, rho2, sigma, s
  real(8), dimension(:), allocatable :: rad, xc, yc, zc
  real(8) :: dpdx, dpdy, dpdz  !pressure gradients in case of pressure driven channel flow
  real(8) :: beta, MaxError
  integer :: maxit, it, itime_scheme, ii, BuoyancyCase, numBubble
  integer :: sbx, sby, Nstep
  integer :: maxStep, itmax, iTimeStep
end module module_flow
!=================================================================================================
!=================================================================================================
! Temporary variables
!-------------------------------------------------------------------------------------------------
module module_tmpvar
  real(8), dimension(:,:,:,:), allocatable :: work
  real(8), dimension(:,:,:), allocatable :: tmp
  real(8) :: tcpu(100),t0
end module module_tmpvar
!=================================================================================================
!=================================================================================================
! module_IO: Contains input/output variables and procedures
!-------------------------------------------------------------------------------------------------
module module_IO
  implicit none
  save
  integer :: nout, out, output_format, nbackup
  character(len=20) :: out_path, x_file, y_file, z_file
  logical :: read_x, read_y, read_z, restart, ICOut
  contains
!=================================================================================================
!=================================================================================================
! function int2text
!   Returns 'number' as a string with length of 'length'
!   called in:    function output
!-------------------------------------------------------------------------------------------------
function int2text(number,length)
  integer :: number, length, i
  character(len=length) :: int2text
  character, dimension(0:9) :: num = (/'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/)
  if(number>=10**length)stop 'int2text error: not enough large string'
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
  OPEN(UNIT=7,FILE=trim(out_path)//'/backup_'//int2text(rank,3),status='unknown',action='write')
  write(7,1100)time,itimestep,is,ie,js,je,ks,ke
  do k=ks,ke; do j=js,je; do i=is,ie
    write(7,1200) u(i,j,k), v(i,j,k), w(i,j,k), p(i,j,k), rho(i,j,k)
  enddo; enddo; enddo
  close(7)
  1100 FORMAT(es25.16e3,7I5)
  1200 FORMAT(5es25.16e3)
end subroutine backup_write
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
subroutine backup_read
  use module_flow
  use module_grid
  implicit none
  integer ::i,j,k,i1,i2,j1,j2,k1,k2
  OPEN(UNIT=7,FILE=trim(out_path)//'/backup_'//int2text(rank,3),status='old',action='read')
  read(7,1100)time,itimestep,i1,i2,j1,j2,k1,k2
  if(i1/=is .or. i2/=ie .or. j1/=js .or. j2/=je .or. k1/=ks .or. k2/=ke) &
    stop 'Error: backup_read'
  do k=ks,ke; do j=js,je; do i=is,ie
    read(7,1200) u(i,j,k), v(i,j,k), w(i,j,k), p(i,j,k), rho(i,j,k)
  enddo; enddo; enddo
  close(7)
  1100 FORMAT(es25.16e3,7I5)
  1200 FORMAT(5es25.16e3)
end subroutine backup_read
!=================================================================================================
!=================================================================================================
! subroutine output
!   Writes the solution in the 'out_path' directory in the tecplot or vtk format
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
!  use module_tmpvar
  !use IO_mod
  implicit none
  integer ::nf,i1,i2,j1,j2,k1,k2,i,j,k
!  logical, save :: first_time=.true.

  OPEN(UNIT=7,FILE=trim(out_path)//'/plot'//int2text(nf,3)//'_'//int2text(rank,3)//'.dat')
  write(7,1000)
  write(7,1100) time, i2-i1+1, j2-j1+1, k2-k1+1
  do k=k1,k2; do j=j1,j2; do i=i1,i2;
  write(7,1200) x(i),y(j),z(k),0.5d0*(u(i,j,k)+u(i-1,j,k)), &
      0.5d0*(v(i,j,k)+v(i,j-1,k)),0.5d0*(w(i,j,k)+w(i,j,k-1)), &
!  write(7,1200) x(i),y(j),z(k),0.5d0*(dIdx(i,j,k)+dIdx(i,j,k)), &
!      0.5d0*(dIdy(i,j,k)+dIdy(i,j,k)),0.5d0*(dIdz(i,j,k)+dIdz(i,j,k)), &
      p(i,j,k) !, rho(i,j,k)
  enddo; enddo; enddo
  close(7)
  1000 FORMAT('variables = "x","y","z","u","v","w","p"') !,"rho"')
  1100 FORMAT('ZONE solutiontime=', 1PG15.7e2, ', I=',I4, 2X, ', J=',I4, 2X, ', K=',I4, 2X)
  1200 FORMAT(20es14.6e2)
end subroutine output1
!-------------------------------------------------------------------------------------------------
subroutine output2(nf,i1,i2,j1,j2,k1,k2)
  use module_flow
  use module_grid
  !use IO_mod
  implicit none
  integer ::nf,i1,i2,j1,j2,k1,k2,i,j,k, itype=1
!  logical, save :: first_time=.true.
  
    OPEN(UNIT=8,FILE=trim(out_path)//'/plot'//int2text(nf,3)//'_'//int2text(rank,3)//'.vtk')
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

    close(7)
end subroutine output2
!-------------------------------------------------------------------------------------------------
end module module_IO
!=================================================================================================
!=================================================================================================
! module module_BC: Sets the boundary conditions of the problem.
!-------------------------------------------------------------------------------------------------
module module_BC
  use module_grid
  implicit none
  character(20) :: bdry_cond(3)
  ! bdry_cond(i) = is the type if boundary condition in i'th direction
  real(8) :: WallVel(6,3)
  ! Tangential velocities on the surfaces of domain. First index represent the 
  ! side on which the velocity in the direction of the second index is specified.
  ! The sides are in this order: -x,+x,-y,+y,-z,+z.
  ! Example: WallVel(4,3) represent the W velocity on +y side of the domain.
  contains
!=================================================================================================
!=================================================================================================
! subroutine SetPressureBC: Sets the pressure boundary condition
!-------------------------------------------------------------------------------------------------
  subroutine SetPressureBC(density)
    implicit none
    include 'mpif.h'
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: density
    real(8), save :: Large=1.0e20
  ! for walls set the density to a large value in the ghost cells
    if(bdry_cond(1)=='wall')then
      if(coords(1)==0    ) density(is-1,js-1:je+1,ks-1:ke+1)=Large
      if(coords(1)==nPx-1) density(ie+1,js-1:je+1,ks-1:ke+1)=Large
    endif

    if(bdry_cond(2)=='wall')then
      if(coords(2)==0    ) density(is-1:ie+1,js-1,ks-1:ke+1)=Large
      if(coords(2)==nPy-1) density(is-1:ie+1,je+1,ks-1:ke+1)=Large
    endif

    if(bdry_cond(3)=='wall')then
      if(coords(3)==0    ) density(is-1:ie+1,js-1:je+1,ks-1)=Large
      if(coords(3)==nPz-1) density(is-1:ie+1,js-1:je+1,ke+1)=Large
    endif
  end subroutine SetPressureBC
!=================================================================================================
!=================================================================================================
! subroutine SetVelocityBC: Sets the velocity boundary condition
!-------------------------------------------------------------------------------------------------
  subroutine SetVelocityBC(u,v,w)
    use module_grid
    implicit none
    include 'mpif.h'
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: u, v, w
    ! tangential velocity at boundaries
    if(bdry_cond(1)=='wall')then
      if(coords(1)==0    ) then
        v(is-1,:,:)=2*WallVel(1,2)-v(is,:,:)
        w(is-1,:,:)=2*WallVel(1,3)-w(is,:,:)
      endif
      if(coords(1)==nPx-1) then
        v(ie+1,:,:)=2*WallVel(2,2)-v(ie,:,:)
        w(ie+1,:,:)=2*WallVel(2,3)-w(ie,:,:)
      endif
    endif
    if(bdry_cond(2)=='wall')then
      if(coords(2)==0    ) then
        u(:,js-1,:)=2*WallVel(3,1)-u(:,js,:)
        w(:,js-1,:)=2*WallVel(3,3)-w(:,js,:)
      endif
      if(coords(2)==nPy-1) then
        u(:,je+1,:)=2*WallVel(4,1)-u(:,je,:)
        w(:,je+1,:)=2*WallVel(4,3)-w(:,je,:)
      endif
    endif
    if(bdry_cond(3)=='wall')then
      if(coords(3)==0    ) then
        u(:,:,ks-1)=2*WallVel(5,1)-u(:,:,ks)
        v(:,:,ks-1)=2*WallVel(5,2)-v(:,:,ks)
      endif
      if(coords(3)==nPz-1) then
        u(:,:,ke+1)=2*WallVel(6,1)-u(:,:,ke)
        v(:,:,ke+1)=2*WallVel(6,2)-v(:,:,ke)
      endif
    endif
  end subroutine SetVelocityBC
!=================================================================================================
!=================================================================================================
! subroutine ghost_x: fills the ghost cells in x-direction from the adjacent core
!-------------------------------------------------------------------------------------------------
  subroutine ghost_x(F,ngh,req)
    use module_grid
    implicit none
    include 'mpif.h'
    integer, intent(in) :: ngh ! number of ghost cell layers to fill
    integer, intent(out) :: req(4)
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: F
    integer :: jlen, klen, ierr !,sta(MPI_STATUS_SIZE,4)
    integer, save :: srcL, srcR, destL, destR, face(2)
    logical, save :: first_time=.true.

    if(ngh>Ng) stop 'ghost error: not enought ghost layers to fill'
    if(first_time)then
      first_time=.false.
      jlen=jmax-jmin+1; klen=kmax-kmin+1; !ilen=ngh
      call para_type_block3a(imin, imax, jmin, jmax, 1, jlen, klen, MPI_DOUBLE_PRECISION, face(1))
      call para_type_block3a(imin, imax, jmin, jmax, 2, jlen, klen, MPI_DOUBLE_PRECISION, face(2))
      call MPI_CART_SHIFT(MPI_COMM_CART, 0, 1, srcR, destR, ierr)
      call MPI_CART_SHIFT(MPI_COMM_CART, 0,-1, srcL, destL, ierr)
    endif

    call MPI_IRECV(F(is-ngh  ,jmin,kmin),1,face(ngh),srcR ,0,MPI_COMM_WORLD,req(1),ierr)
    call MPI_ISEND(F(ie-ngh+1,jmin,kmin),1,face(ngh),destR,0,MPI_COMM_WORLD,req(2),ierr)
    call MPI_IRECV(F(ie+1    ,jmin,kmin),1,face(ngh),srcL ,0,MPI_COMM_WORLD,req(3),ierr)
    call MPI_ISEND(F(is      ,jmin,kmin),1,face(ngh),destL,0,MPI_COMM_WORLD,req(4),ierr)
!    call MPI_WAITALL(4,req,sta,ierr)
  end subroutine ghost_x
!-------------------------------------------------------------------------------------------------
  subroutine ghost_y(F,ngh,req)
    use module_grid
    implicit none
    include 'mpif.h'
    integer, intent(in) :: ngh
    integer, intent(out) :: req(4)
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: F
    integer :: ilen, klen, ierr !,sta(MPI_STATUS_SIZE,4)
    integer, save :: srcL, srcR, destL, destR, face(2)
    logical, save :: first_time=.true.

    if(ngh>Ng) stop 'ghost error: not enought ghost layers to fill'
    if(first_time)then
      first_time=.false.
      klen=kmax-kmin+1; ilen=imax-imin+1; !jlen=ngh
      call para_type_block3a(imin, imax, jmin, jmax, ilen, 1, klen, MPI_DOUBLE_PRECISION, face(1))
      call para_type_block3a(imin, imax, jmin, jmax, ilen, 2, klen, MPI_DOUBLE_PRECISION, face(2))
      call MPI_CART_SHIFT(MPI_COMM_CART, 1, 1, srcR, destR, ierr)
      call MPI_CART_SHIFT(MPI_COMM_CART, 1,-1, srcL, destL, ierr)
    endif

    call MPI_IRECV(F(imin,js-ngh  ,kmin),1,face(ngh),srcR ,0,MPI_COMM_WORLD,req(1),ierr)
    call MPI_ISEND(F(imin,je-ngh+1,kmin),1,face(ngh),destR,0,MPI_COMM_WORLD,req(2),ierr)
    call MPI_IRECV(F(imin,je+1    ,kmin),1,face(ngh),srcL ,0,MPI_COMM_WORLD,req(3),ierr)
    call MPI_ISEND(F(imin,js      ,kmin),1,face(ngh),destL,0,MPI_COMM_WORLD,req(4),ierr)
!    call MPI_WAITALL(4,req,sta,ierr)
  end subroutine ghost_y
!-------------------------------------------------------------------------------------------------
  subroutine ghost_z(F,ngh,req)
    use module_grid
    implicit none
    include 'mpif.h'
    integer, intent(in) :: ngh
    integer, intent(out) :: req(4)
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: F
    integer :: ilen, jlen, ierr !,sta(MPI_STATUS_SIZE,4)
    integer, save :: srcL, srcR, destL, destR, face(2)
    logical, save :: first_time=.true.

    if(ngh>Ng) stop 'ghost error: not enought ghost layers to fill'
    if(first_time)then
      first_time=.false.
      ilen=imax-imin+1; jlen=jmax-jmin+1; !klen=ngh
      call para_type_block3a(imin, imax, jmin, jmax, ilen, jlen, 1, MPI_DOUBLE_PRECISION, face(1))
      call para_type_block3a(imin, imax, jmin, jmax, ilen, jlen, 2, MPI_DOUBLE_PRECISION, face(2))
      call MPI_CART_SHIFT(MPI_COMM_CART, 2, 1, srcR, destR, ierr)
      call MPI_CART_SHIFT(MPI_COMM_CART, 2,-1, srcL, destL, ierr)
    endif

    call MPI_IRECV(F(imin,jmin,ks-ngh  ),1,face(ngh),srcR ,0,MPI_COMM_WORLD,req(1),ierr)
    call MPI_ISEND(F(imin,jmin,ke-ngh+1),1,face(ngh),destR,0,MPI_COMM_WORLD,req(2),ierr)
    call MPI_IRECV(F(imin,jmin,ke+1    ),1,face(ngh),srcL ,0,MPI_COMM_WORLD,req(3),ierr)
    call MPI_ISEND(F(imin,jmin,ks      ),1,face(ngh),destL,0,MPI_COMM_WORLD,req(4),ierr)
!    call MPI_WAITALL(4,req,sta,ierr)
  end subroutine ghost_z
!-------------------------------------------------------------------------------------------------
end module module_BC
!=================================================================================================
!=================================================================================================
! creates new data type for passing non-contiguous data in ghost_x
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
! module_poisson: solves the pressure poisson equation using hypre
!-------------------------------------------------------------------------------------------------
module module_poisson
  use module_grid
  use module_BC
  use module_tmpvar
  implicit none
  integer :: nstencil
  integer*8 :: grid_obj, stencil, Amat, Bvec, Xvec, solver, precond
  integer, dimension(:), allocatable :: stencil_indices, ilower, iupper
  contains
!-------------------------------------------------------------------------------------------------
  subroutine poi_initialize
    implicit none
    include 'mpif.h'
    integer :: ierr, periodic_array(3), i
    integer, dimension(:,:), allocatable :: offsets

    nstencil = 2 * ndim + 1 !7
    allocate(ilower(ndim), iupper(ndim), stencil_indices(nstencil) )
    ilower = (/is, js, ks/);  iupper = (/ie, je, ke/)
    periodic_array = 0
    if( bdry_cond(1)=='periodic') periodic_array(1) = Nx
    if( bdry_cond(2)=='periodic') periodic_array(2) = Ny
    if( bdry_cond(3)=='periodic') periodic_array(3) = Nz

    allocate(offsets(ndim,nstencil), stat=ierr)
    offsets(:,1) = (/ 0, 0, 0/)
    offsets(:,2) = (/-1, 0, 0/)
    offsets(:,3) = (/ 1, 0, 0/)
    offsets(:,4) = (/ 0,-1, 0/)
    offsets(:,5) = (/ 0, 1, 0/)
    offsets(:,6) = (/ 0, 0,-1/)
    offsets(:,7) = (/ 0, 0, 1/)

    call HYPRE_StructGridCreate(mpi_comm_world, ndim, grid_obj, ierr)  ! create a 3d grid object
    call HYPRE_StructGridSetExtents(grid_obj, ilower, iupper, ierr)    ! add a box to the grid
    call HYPRE_StructGridSetPeriodic(grid_obj, periodic_array, ierr)   ! set periodic
    call HYPRE_StructGridAssemble(grid_obj, ierr)                      ! assemble the grid
    call HYPRE_StructStencilCreate(ndim, nstencil, stencil, ierr)      ! create a 3d 7-pt stencil

    do i = 1, nstencil          ! define the geometry of the stencil
      stencil_indices(i) = i-1
      call HYPRE_StructStencilSetElement(stencil, stencil_indices(i), offsets(:,i), ierr)
    enddo

    call HYPRE_StructMatrixCreate(mpi_comm_world, grid_obj, stencil, Amat, ierr) ! set up matrix A
    call HYPRE_StructMatrixInitialize(Amat, ierr)

    call HYPRE_StructVectorCreate(mpi_comm_world, grid_obj, Bvec, ierr)! create vector object B
    call HYPRE_StructVectorInitialize(Bvec, ierr)

    call HYPRE_StructVectorCreate(mpi_comm_world, grid_obj, Xvec, ierr)! create vector object X
    call HYPRE_StructVectorInitialize(Xvec, ierr)

  end subroutine poi_initialize
!-------------------------------------------------------------------------------------------------
  subroutine poi_solve(utmp,vtmp,wtmp,rhot,dt,p,maxError,maxit)
    use module_grid
    implicit none
    include 'mpif.h'
    integer :: ierr, nvalues, ijk, i,j,k
    real(8), dimension(:), allocatable :: values
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: utmp,vtmp,wtmp,rhot
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: p
    real(8), intent(in) :: dt,maxError
    integer, intent(in) :: maxit
    !real(8) :: final_res_norm
!---------------------------------------------Setup matrix A--------------------------------------
    nvalues = mx * my * mz * nstencil
    allocate(values(nvalues), stat=ierr)
    if(ierr/=0)stop '**** poi_solve: allocation error ****'

    ijk = 1
    do k=ks,ke;  do j=js,je;  do i=is,ie
!      ijk = 1 + i-is + mx*( j-js + my*( k-ks ) )
      values(ijk+1) = -(1./dx(i))*( 1./(dxh(i-1)*(rhot(i-1,j,k)+rhot(i,j,k))) )
      values(ijk+2) = -(1./dx(i))*( 1./(dxh(i  )*(rhot(i+1,j,k)+rhot(i,j,k))) )
      values(ijk+3) = -(1./dy(j))*( 1./(dyh(j-1)*(rhot(i,j-1,k)+rhot(i,j,k))) )
      values(ijk+4) = -(1./dy(j))*( 1./(dyh(j  )*(rhot(i,j+1,k)+rhot(i,j,k))) )
      values(ijk+5) = -(1./dz(k))*( 1./(dzh(k-1)*(rhot(i,j,k-1)+rhot(i,j,k))) )
      values(ijk+6) = -(1./dz(k))*( 1./(dzh(k  )*(rhot(i,j,k+1)+rhot(i,j,k))) )
      values(ijk  ) = -sum(values(ijk+1:ijk+6))
      ijk = ijk + 7
    enddo;  enddo;  enddo
    call HYPRE_StructMatrixSetBoxValues(Amat, ilower, iupper, nstencil, stencil_indices, &
                                        values, ierr)
    deallocate(values, stat=ierr)
    !if(rank==0)print*,'Matix A done.'
!------------------------------------------SET UP THE SOURCE TERM and initial guess---------------
    nvalues = mx * my * mz
    allocate(values(nvalues), stat=ierr)

    ijk = 1
    do k=ks,ke;  do j=js,je;  do i=is,ie
      values(ijk) = -0.5/dt *( (utmp(i,j,k)-utmp(i-1,j,k))/dx(i) &
                              +(vtmp(i,j,k)-vtmp(i,j-1,k))/dy(j) &
                              +(wtmp(i,j,k)-wtmp(i,j,k-1))/dz(k)   )
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
    !if(rank==0)print*,'Vector B done.'
!-----------------------------------------------Anchor point p(1,1,1)=0---------------------------
    if(coords(1)==0 .and. coords(2)==0 .and. coords(3) == 0)then
      nvalues = nstencil
      allocate(values(nvalues), stat=ierr)
      values = 0d0
      values(1) = 1d0
      call HYPRE_StructMatrixAddToBoxValues(Amat, ilower, ilower, nstencil, stencil_indices, &
                                            values, ierr)
      call HYPRE_StructVectorSetBoxValues(Bvec, ilower, ilower, 0d0, ierr)
      deallocate(values, stat=ierr)
    endif
    !if(rank==0)print*,'Anchor done.'
!-----------------------------------------------Assemble matrix A and vectors B,X-----------------
    call HYPRE_StructMatrixAssemble(Amat, ierr)
    call HYPRE_StructVectorAssemble(Bvec, ierr)
    call HYPRE_StructVectorAssemble(Xvec, ierr)
    !if(rank==0)print*,'Assemble done.'
!-----------------------------------------SOLVE THE EQUATIONS-------------------------------------
    call HYPRE_StructPFMGCreate(mpi_comm_world, solver, ierr)
    call HYPRE_StructPFMGSetMaxIter(solver, maxit, ierr)
    call HYPRE_StructPFMGSetTol(solver, MaxError, ierr)

!    call HYPRE_StructPFMGSetMaxLevels(solver, 4, ierr)    
    call HYPRE_StructPFMGSetRelaxType(solver, 3, ierr)    
!    call HYPRE_StructPFMGSetNumPreRelax(solver, 3, ierr)
!    call HYPRE_StructPFMGSetNumPostRelax(solver, 0, ierr)

    call hypre_structPFMGsetLogging(solver, 1, ierr)
    call HYPRE_StructPFMGSetup(solver, Amat, Bvec, Xvec, ierr)
    call HYPRE_StructPFMGSolve(solver, Amat, Bvec, Xvec, ierr)
!    call hypre_structPFMGgetnumiterations(solver, num_iterations, ierr)
!    call hypre_structPFMGgetfinalrelative(solver, final_res_norm, ierr)

    nvalues = mx * my * mz
    allocate(values(nvalues), stat=ierr)
    call HYPRE_StructVectorGetBoxValues(Xvec, ilower, iupper, values, ierr)
    call HYPRE_StructPFMGDestroy(solver, ierr)

    ! set the values of b vector
    ijk = 1
    do k=ks,ke;  do j=js,je;  do i=is,ie
      p(i,j,k) = values(ijk)
      ijk = ijk + 1
    enddo;  enddo;  enddo
    deallocate(values, stat=ierr)
  end subroutine poi_solve
!-------------------------------------------------------------------------------------------------
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
! Main program to solve the 3D NS equations for multiphase flows
!-------------------------------------------------------------------------------------------------
Program ftc3d2011
  use module_flow
  use module_grid
  use module_BC
  use module_tmpvar
  use module_poisson
  use module_IO
  implicit none
  include 'mpif.h'
  integer :: ierr, i,j,k
  ! Locals for marching and timing
  real(8) :: start_time, end_time !, IC_time, wall_time
  integer :: req(48),sta(MPI_STATUS_SIZE,48)
  logical, allocatable, dimension(:,:,:) :: mask
!---------------------------------------INITIALIZATION--------------------------------------------
  ! Initialize MPI
  call MPI_INIT(ierr)
  If (ierr /= 0) STOP '*** Main: unsuccessful MPI-initialization'
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numProcess, ierr)
  start_time = MPI_Wtime(ierr)

  call ReadParameters
  if(rank==0) write(out,*)'Parameters read successfully'
  ! check number of processors
  If (numProcess .NE. nPx*nPy*nPz) STOP '*** Main: Problem with nPx!'

  call initialize
  call InitCondition
  if(rank==0) write(out,*)'initialized'

  if(HYPRE) call poi_initialize
  if(HYPRE .and. rank==0) write(out,*)'hypre initialized'

  allocate(mask(imin:imax,jmin:jmax,kmin:kmax))
  ! output initial condition
  if(ICOut)call output(0,is,ie,js,je,ks,ke)
  if(rank==0)start_time = MPI_WTIME()
!---------------------------------------MAIN TIME LOOP--------------------------------------------
  do while(time<EndTime .and. itimestep<nstep)
    if(dtFlag==2)call TimeStepSize(dt)
    time=time+dt
    itimestep=itimestep+1
    if(rank==0)write(out,'("Step: ",I6," dt=",es16.5e2," time=",es16.5e2)')itimestep,dt,time
    if(itime_scheme==2) then
      uold = u
      vold = v
      wold = w
      rhoo = rho
      muold  = mu
    endif
!------------------------------------ADVECTION & DIFFUSION----------------------------------------
    do ii=1, itime_scheme
      call momentumDiffusion(u,v,w,rho,mu,du,dv,dw)
      call momentumConvection(u,v,w,du,dv,dw)
      call volumeForce(rho,rho1,rho2,dpdx,dpdy,dpdz,BuoyancyCase,fx,fy,fz,gx,gy,gz,du,dv,dw)
      if(TwoPhase) call density(rho,u,v,w,drho)
      u = u + dt * du
      v = v + dt * dv
      w = w + dt * dw
      if(TwoPhase) rho = rho + dt * drho

      call SetVelocityBC(u,v,w)
      
      call ghost_x(u  ,2,req( 1: 4));  call ghost_x(v,2,req( 5: 8)); call ghost_x(w,2,req( 9:12)); 
      call ghost_x(rho,2,req(13:16));  call MPI_WAITALL(16,req(1:16),sta(:,1:16),ierr)
      call ghost_y(u  ,2,req( 1: 4));  call ghost_y(v,2,req( 5: 8)); call ghost_y(w,2,req( 9:12)); 
      call ghost_y(rho,2,req(13:16));  call MPI_WAITALL(16,req(1:16),sta(:,1:16),ierr)
      call ghost_z(u  ,2,req( 1: 4));  call ghost_z(v,2,req( 5: 8)); call ghost_z(w,2,req( 9:12)); 
      call ghost_z(rho,2,req(13:16));  call MPI_WAITALL(16,req(1:16),sta(:,1:16),ierr)

      if(TwoPhase) mu  = (rho-rho1)/(rho2-rho1)*(mu2-mu1)+mu1
    enddo

    if(itime_scheme==2) then
      u = 0.5*(u+uold)
      v = 0.5*(v+vold)
      w = 0.5*(w+wold)
      if(TwoPhase) then
        rho = 0.5*(rho+rhoo)
        mu  = 0.5*(mu +muold )
      endif
    endif
!--------------------------------------PROJECTION STEP--------------------------------------------
    ! Compute source term and the coefficient do p(i,j,k)
!    tmp(is:ie,js:je,ks:ke)=rho(is:ie,js:je,ks:ke)
!    tmp = rho
!    call SetPressureBC(density=tmp,pressure=p)
!    call SetupPoisson(u,v,w,tmp,dt,work)

    tmp = rho
    call SetPressureBC(tmp)
    if(HYPRE)then
      call poi_solve(u,v,w,tmp,dt,p,maxError,maxit)
      call ghost_x(p,1,req(1:4 ))
      call ghost_y(p,1,req(5:8 ))
      call ghost_z(p,1,req(9:12))
      call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
    else
      mask=.true.
      call SetupPoisson(u,v,w,tmp,dt,work,mask)
      call SOR_Solver(work,p,maxError,beta,maxit,it)
    endif
    do k=ks,ke;  do j=js,je; do i=is,ieu;    ! CORRECT THE u-velocity 
      u(i,j,k)=u(i,j,k)-dt*(2.0/dxh(i))*(p(i+1,j,k)-p(i,j,k))/(rho(i+1,j,k)+rho(i,j,k))
    enddo; enddo; enddo
      
    do k=ks,ke;  do j=js,jev; do i=is,ie;    ! CORRECT THE v-velocity
      v(i,j,k)=v(i,j,k)-dt*(2.0/dyh(j))*(p(i,j+1,k)-p(i,j,k))/(rho(i,j+1,k)+rho(i,j,k))
    enddo; enddo; enddo
    
    do k=ks,kew;  do j=js,je; do i=is,ie;   ! CORRECT THE v-velocity
      w(i,j,k)=w(i,j,k)-dt*(2.0/dzh(k))*(p(i,j,k+1)-p(i,j,k))/(rho(i,j,k+1)+rho(i,j,k))
    enddo; enddo; enddo

    call SetVelocityBC(u,v,w)

    call ghost_x(u,2,req( 1: 4)); call ghost_x(v,2,req( 5: 8)); call ghost_x(w,2,req( 9:12))
    call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
    call ghost_y(u,2,req( 1: 4)); call ghost_y(v,2,req( 5: 8)); call ghost_y(w,2,req( 9:12))
    call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)
    call ghost_z(u,2,req( 1: 4)); call ghost_z(v,2,req( 5: 8)); call ghost_z(w,2,req( 9:12))
    call MPI_WAITALL(12,req(1:12),sta(:,1:12),ierr)

!--------------------------------------------OUTPUT-----------------------------------------------
    ! calculate the volumetric flow rate in x-direction
    i=is
    myfr=0d0
    do k=ks,ke;  do j=js,je
      myfr=myfr+u(i,j,k)*dy(j)*dz(k)
    enddo; enddo
    call MPI_ALLREDUCE(myfr, flowrate, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    flowrate = flowrate/dfloat(nPx)
    if(mod(itimestep,nbackup)==0)call backup_write
    if(mod(itimestep,nout)==0)call output(ITIMESTEP/nout,is,ie,js,je,ks,ke)
    if(rank==0)then
      end_time =  MPI_WTIME()
      write(out,'("Step:",I6," Iterations:",I5," cpu(s):",f10.2)')itimestep,it,end_time-start_time
      open(unit=121,file='track')
      write(121,'("Step:",I6," dt=",es16.5e2," time=",es16.5e2)')itimestep,dt,time
      write(121,'("            Iterations:",I5," cpu(s):",f10.2)')it,end_time-start_time
      close(121)
      write(*,'("Step:",I6," dt=",es16.5e2," time=",es16.5e2)')itimestep,dt,time
      write(*,'("            Iterations:",I5," cpu(s):",f10.2)')it,end_time-start_time
      open(unit=121,file='stats',access='append')
      write(121,'(5es16.5e2)')time,flowrate,end_time-start_time
      close(121)
    endif
  enddo
  if(HYPRE) call poi_finalize
  call MPI_FINALIZE(ierr)
  stop
end program ftc3d2011
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
subroutine TimeStepSize(deltaT)
  use module_grid
  use module_flow
  implicit none
  include "mpif.h"
  real(8) :: deltaT, h, vmax, dtadv, mydt !, vmaxglob
  real(8), save :: dtvisc, nu_min, nu_max !, hglob
!  logical, save :: first_time=.true.
  integer :: ierr
  h  = min(minval(dx),minval(dy),minval(dz))
  nu_min = mu1/rho1;   if(TwoPhase) nu_min = min(nu_min,mu2/rho2)
  nu_max = mu1/rho1;   if(TwoPhase) nu_max = max(nu_max,mu2/rho2)
  dtvisc=h**2/(6.0d0*nu_max)
  vmax = maxval(sqrt(u(is:ie,js:je,ks:ke)**2+v(is:ie,js:je,ks:ke)**2+w(is:ie,js:je,ks:ke)**2))
  vmax = max(vmax,1e-3)
!  dtadv  = min(h/vmax,2d0*nu_min/vmax**2)
  dtadv  = h/vmax
  mydt = CFL*min(dtvisc,dtadv)
  mydt = min(mydt,MaxDt)

  call MPI_ALLREDUCE(mydt, deltat, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
end subroutine TimeStepSize
!=================================================================================================
!=================================================================================================
! subroutine momentumConvection
! calculates the convection terms in mumentum equation using ENO scheme
! and returns them in du, dv, dw
!-------------------------------------------------------------------------------------------------
subroutine momentumConvection(u,v,w,du,dv,dw)
  use module_grid
  use module_tmpvar
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: u, v, w
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: du, dv, dw
  real(8), external :: minabs
  integer :: i,j,k
!-------------------------------------ENO interpolation u-velocity--------------------------------
  do k=ks,ke; do j=js,je; do i=is,ieu+1
    if (u(i-1,j,k)+u(i,j,k)>0.0) then
      work(i,j,k,1)=u(i-1,j,k)+0.5*minabs((u(i,j,k)-u(i-1,j,k)),(u(i-1,j,k)-u(i-2,j,k)))
    else
      work(i,j,k,1)=u(i,j,k)-0.5*minabs((u(i+1,j,k)-u(i,j,k)),(u(i,j,k)-u(i-1,j,k)))
    endif
  enddo; enddo; enddo
  do k=ks,ke; do j=js-1,je; do i=is-1,ie
    if(v(i,j,k)+v(i+1,j,k)>0.0) then
      work(i,j,k,2)=u(i,j,k)+0.5*minabs((u(i,j+1,k)-u(i,j,k)),(u(i,j,k)-u(i,j-1,k)))
    else
      work(i,j,k,2)=u(i,j+1,k)-0.5*minabs((u(i,j+2,k)-u(i,j+1,k)),(u(i,j+1,k)-u(i,j,k)))
    endif
  enddo; enddo; enddo
  do k=ks-1,ke; do j=js,je; do i=is-1,ie
    if(w(i,j,k)+w(i+1,j,k)>0.0) then
      work(i,j,k,3)=u(i,j,k)+0.5*minabs((u(i,j,k+1)-u(i,j,k)),(u(i,j,k)-u(i,j,k-1)))
    else
      work(i,j,k,3)=u(i,j,k+1)-0.5*minabs((u(i,j,k+2)-u(i,j,k+1)),(u(i,j,k+1)-u(i,j,k)))
    endif
  enddo; enddo; enddo
  do k=ks,ke;  do j=js,je; do i=is,ieu
    du(i,j,k)=du(i,j,k)-0.5*((u(i,j  ,k  )+u(i+1,j  ,k  ))*work(i+1,j  ,k  ,1)-&
                    (u(i,j  ,k  )+u(i-1,j  ,k  ))*work(i  ,j  ,k  ,1))/dxh(i) &
              -0.5*((v(i,j  ,k  )+v(i+1,j  ,k  ))*work(i  ,j  ,k  ,2)-&
                    (v(i,j-1,k  )+v(i+1,j-1,k  ))*work(i  ,j-1,k  ,2))/dy(j)  &
              -0.5*((w(i,j  ,k  )+w(i+1,j  ,k  ))*work(i  ,j  ,k  ,3)-&
                    (w(i,j  ,k-1)+w(i+1,j  ,k-1))*work(i  ,j  ,k-1,3))/dz(k)
  enddo; enddo; enddo

!-------------------------------------ENO interpolation v-velocity--------------------------------
  do k=ks,ke; do j=js,jev+1; do i=is,ie
    if (v(i,j-1,k)+v(i,j,k)>0.0) then
      work(i,j,k,2)=v(i,j-1,k)+0.5*minabs((v(i,j,k)-v(i,j-1,k)),(v(i,j-1,k)-v(i,j-2,k)))
    else
      work(i,j,k,2)=v(i,j,k)-0.5*minabs((v(i,j+1,k)-v(i,j,k)),(v(i,j,k)-v(i,j-1,k)))
    endif
  enddo; enddo; enddo
  do k=ks,ke; do j=js-1,je; do i=is-1,ie
    if(u(i,j,k)+u(i,j+1,k)>0.0) then
      work(i,j,k,1)=v(i,j,k)+0.5*minabs((v(i+1,j,k)-v(i,j,k)),(v(i,j,k)-v(i-1,j,k)))
    else
      work(i,j,k,1)=v(i+1,j,k)-0.5*minabs((v(i+2,j,k)-v(i+1,j,k)),(v(i+1,j,k)-v(i,j,k)))
    endif
  enddo; enddo; enddo
  do k=ks-1,ke; do j=js-1,je; do i=is,ie
    if(w(i,j,k)+w(i,j+1,k)>0.0) then
      work(i,j,k,3)=v(i,j,k)+0.5*minabs((v(i,j,k+1)-v(i,j,k)),(v(i,j,k)-v(i,j,k-1)))
    else
      work(i,j,k,3)=v(i,j,k+1)-0.5*minabs((v(i,j,k+2)-v(i,j,k+1)),(v(i,j,k+1)-v(i,j,k)))
    endif
  enddo; enddo; enddo
  do k=ks,ke;  do j=js,jev; do i=is,ie
    dv(i,j,k)=dv(i,j,k)-0.5*((u(i  ,j,k  )+u(i  ,j+1,k  ))*work(i  ,j  ,k  ,1)-&
                    (u(i-1,j,k  )+u(i-1,j+1,k  ))*work(i-1,j  ,k  ,1))/dx(i)  &
              -0.5*((v(i  ,j,k  )+v(i  ,j+1,k  ))*work(i  ,j+1,k  ,2)-&
                    (v(i  ,j,k  )+v(i  ,j-1,k  ))*work(i  ,j  ,k  ,2))/dyh(j) &
              -0.5*((w(i  ,j,k  )+w(i  ,j+1,k  ))*work(i  ,j  ,k  ,3)-&
                    (w(i  ,j,k-1)+w(i  ,j+1,k-1))*work(i  ,j  ,k-1,3))/dz(k)
  enddo; enddo; enddo
!-------------------------------------ENO interpolation w-velocity--------------------------------
  do k=ks,kew+1; do j=js,je; do i=is,ie
    if (w(i,j,k-1)+w(i,j,k)>0.0) then
      work(i,j,k,3)=w(i,j,k-1)+0.5*minabs((w(i,j,k)-w(i,j,k-1)),(w(i,j,k-1)-w(i,j,k-2)))
    else
      work(i,j,k,3)=w(i,j,k)-0.5*minabs((w(i,j,k+1)-w(i,j,k)),(w(i,j,k)-w(i,j,k-1)))
    endif
  enddo; enddo; enddo
  do k=ks-1,ke; do j=js,je; do i=is-1,ie
    if(u(i,j,k)+u(i,j,k+1)>0.0) then
      work(i,j,k,1)=w(i,j,k)+0.5*minabs((w(i+1,j,k)-w(i,j,k)),(w(i,j,k)-w(i-1,j,k)))
    else
      work(i,j,k,1)=w(i+1,j,k)-0.5*minabs((w(i+2,j,k)-w(i+1,j,k)),(w(i+1,j,k)-w(i,j,k)))
    endif
  enddo; enddo; enddo
  do k=ks-1,ke; do j=js-1,je; do i=is,ie
    if(v(i,j,k)+v(i,j,k+1)>0.0) then
      work(i,j,k,2)=w(i,j,k)+0.5*minabs((w(i,j+1,k)-w(i,j,k)),(w(i,j,k)-w(i,j-1,k)))
    else
      work(i,j,k,2)=w(i,j+1,k)-0.5*minabs((w(i,j+2,k)-w(i,j+1,k)),(w(i,j+1,k)-w(i,j,k)))
    endif
  enddo; enddo; enddo
  do k=ks,kew;  do j=js,je; do i=is,ie
    dw(i,j,k)=dw(i,j,k)-0.5*((u(i  ,j  ,k)+u(i  ,j  ,k+1))*work(i  ,j  ,k  ,1)-&
                    (u(i-1,j  ,k)+u(i-1,j  ,k+1))*work(i-1,j  ,k  ,1))/dx(i)  &
              -0.5*((v(i  ,j  ,k)+v(i  ,j  ,k+1))*work(i  ,j  ,k  ,2)-&
                    (v(i  ,j-1,k)+v(i  ,j-1,k+1))*work(i  ,j-1,k  ,2))/dy(j)  &
              -0.5*((w(i  ,j  ,k)+w(i  ,j  ,k+1))*work(i  ,j  ,k+1,3)-&
                    (w(i  ,j  ,k)+w(i  ,j  ,k-1))*work(i  ,j  ,k  ,3))/dzh(k)
  enddo; enddo; enddo

end subroutine momentumConvection
!=================================================================================================
!=================================================================================================
! Calculates the diffusion terms in the momentum equation and adds them to du,dv,dw
!-------------------------------------------------------------------------------------------------
subroutine momentumDiffusion(u,v,w,rho,mu,du,dv,dw)
  use module_grid
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: u, v, w, rho, mu
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(out) :: du, dv, dw
  integer :: i,j,k
!-------------------------------------PREDICTED u-velocity-DIFFUSION------------------------------
  do k=ks,ke; do j=js,je; do i=is,ieu
    du(i,j,k)= ( & ! du(i,j,k)+(   &
      (1./dxh(i))*2.*(   mu(i+1,j,k)*(1./dx(i+1))*(u(i+1,j,k)-u(i,j,k)) -                &
                         mu(i,j,k)  *(1./dx(i))  *(u(i,j,k)-u(i-1,j,k)) )                &
     +(1./dy(j))*( 0.25*(mu(i,j,k)+mu(i+1,j,k)+mu(i+1,j+1,k)+mu(i,j+1,k))*               &
         ((1./dyh(j))*  (u(i,j+1,k)-u(i,j,k)) + (1./dxh(i))*(v(i+1,j,k)-v(i,j,k)) ) -    &
                   0.25*(mu(i,j,k)+mu(i+1,j,k)+mu(i+1,j-1,k)+mu(i,j-1,k))*               &
         ((1./dyh(j-1))*(u(i,j,k)-u(i,j-1,k)) + (1./dxh(i))*(v(i+1,j-1,k)-v(i,j-1,k))))  &
     +(1./dz(k))*( 0.25*(mu(i,j,k)+mu(i+1,j,k)+mu(i+1,j,k+1)+mu(i,j,k+1))*               &
         ((1./dzh(k))*  (u(i,j,k+1)-u(i,j,k)) + (1./dxh(i))*(w(i+1,j,k)-w(i,j,k)) ) -    &
                   0.25*(mu(i,j,k)+mu(i+1,j,k)+mu(i+1,j,k-1)+mu(i,j,k-1))*               &
         ((1./dzh(k-1))*(u(i,j,k)-u(i,j,k-1)) + (1./dxh(i))*(w(i+1,j,k-1)-w(i,j,k-1))))  &
                          )/(0.5*(rho(i+1,j,k)+rho(i,j,k))   )
  enddo; enddo; enddo
!-------------------------------------PREDICTED v-velocity-DIFFUSION------------------------------
  do k=ks,ke; do j=js,jev; do i=is,ie
    dv(i,j,k)= ( & ! dv(i,j,k)+(   &
      (1./dx(i))*( 0.25*(mu(i,j,k)+mu(i+1,j,k)+mu(i+1,j+1,k)+mu(i,j+1,k))*               &
         ((1./dyh(j))*(u(i,j+1,k)-u(i,j,k)) + (1./dxh(i))*(v(i+1,j,k)-v(i,j,k)) ) -      &
                   0.25*(mu(i,j,k)+mu(i,j+1,k)+mu(i-1,j+1,k)+mu(i-1,j,k))*               &
         ((1./dyh(j))*(u(i-1,j+1,k)-u(i-1,j,k))+ (1./dxh(i-1))*(v(i,j,k)- v(i-1,j,k))))  &
     +(1./dyh(j))*2.*(   mu(i,j+1,k)*(1./dy(j+1))*(v(i,j+1,k)-v(i,j,k)) -                &
                         mu(i,j,k)  *(1./dy(j))*  (v(i,j,k)-v(i,j-1,k)) )                &
     +(1./dz(k))*( 0.25*(mu(i,j,k)+mu(i,j+1,k)+mu(i,j+1,k+1)+mu(i,j,k+1))*               &
         ((1./dzh(k))*  (v(i,j,k+1)-v(i,j,k)) + (1./dyh(j))*(w(i,j+1,k)-w(i,j,k)) ) -    &
                   0.25*(mu(i,j,k)+mu(i,j+1,k)+mu(i,j+1,k-1)+mu(i,j,k-1))*               &
         ((1./dzh(k-1))*(v(i,j,k)-v(i,j,k-1)) + (1./dyh(j))*(w(i,j+1,k-1)-w(i,j,k-1))))  &
                           )/(0.5*(rho(i,j+1,k)+rho(i,j,k))   )
  enddo; enddo; enddo
!-------------------------------------PREDICTED w-velocity-DIFFUSION------------------------------
  do k=ks,kew; do j=js,je; do i=is,ie
    dw(i,j,k)= ( & ! dw(i,j,k)+(   &
      (1./dx(i))*( 0.25*(mu(i,j,k)+mu(i+1,j,k)+mu(i+1,j,k+1)+mu(i,j,k+1))*               &
         ((1./dzh(k))*(u(i,j,k+1)-u(i,j,k)) + (1./dxh(i))*(w(i+1,j,k)-w(i,j,k)) ) -      &
                   0.25*(mu(i,j,k)+mu(i-1,j,k)+mu(i-1,j,k+1)+mu(i,j,k+1))*               &
         ((1./dzh(k))*(u(i-1,j,k+1)-u(i-1,j,k))+ (1./dxh(i-1))*(w(i,j,k)-w(i-1,j,k))) )  &
     +(1./dy(j))*( 0.25*(mu(i,j,k)+mu(i,j+1,k)+mu(i,j+1,k+1)+mu(i,j,k+1))*               &
         ((1./dzh(k))*(v(i,j,k+1)-v(i,j,k)) + (1./dyh(j))*(w(i,j+1,k)-w(i,j,k)) ) -      &
                   0.25*(mu(i,j,k)+mu(i,j-1,k)+mu(i,j-1,k+1)+mu(i,j,k+1))*               &
         ((1./dzh(k))*(v(i,j-1,k+1)-v(i,j-1,k))+ (1./dyh(j-1))*(w(i,j,k)-w(i,j-1,k))) )  &
     +(1./dzh(k))*2.*(   mu(i,j,k+1)*(1./dz(k+1))*(w(i,j,k+1)-w(i,j,k)) -                &
                         mu(i,j,k)  *(1./dz(k))*  (w(i,j,k)-w(i,j,k-1)) )                &
                          )/(0.5*(rho(i,j,k+1)+rho(i,j,k))   )
  enddo; enddo; enddo
end subroutine momentumDiffusion
!=================================================================================================
!=================================================================================================
! Calculates the volume force in the momentum equations and adds them to du,dv,dw
!-------------------------------------------------------------------------------------------------
subroutine volumeForce(rho,rho1,rho2,dpdx,dpdy,dpdz,BuoyancyCase,fx,fy,fz,gx,gy,gz,du,dv,dw)
  use module_grid
  implicit none
  include "mpif.h"
  integer, intent(in) :: BuoyancyCase
  real(8), intent(in) :: gx,gy,gz,rho1,rho2, dpdx, dpdy, dpdz
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: rho, fx,fy,fz
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: du, dv, dw
  real(8) :: vol, rro, ro
  integer :: i,j,k, ierr
  if(BuoyancyCase==0) then
    rro = 0.0
  elseif(BuoyancyCase==1) then
    rro=rho1
  elseif(BuoyancyCase==2) then
    rro=rho2
  elseif(BuoyancyCase==3) then
    ! calculate average density
    vol=sum(dx(is:ie))*sum(dy(js:je))*sum(dz(ks:ke))
    ro=0.0
    do k=ks,ke; do j=js,je; do i=is,ie
      ro = ro+dx(i)*dy(j)*dz(k)*rho(i,j,k)
    enddo; enddo; enddo
    ro=ro/(xLength*yLength*zLength)
    call MPI_ALLREDUCE(ro, rro, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_World, ierr)
  else
    stop 'volumeForce: invalid buoyancy option'
  endif
   
  do k=ks,ke;  do j=js,je; do i=is,ieu
    du(i,j,k)=du(i,j,k)+(fx(i,j,k)-dpdx)/(0.5*(rho(i+1,j,k)+rho(i,j,k))) + &
        (1.0 -rro/(0.5*(rho(i+1,j,k)+rho(i,j,k))) )*gx
  enddo; enddo; enddo

  do k=ks,ke;  do j=js,jev; do i=is,ie
    dv(i,j,k)=dv(i,j,k)+(fy(i,j,k)-dpdy)/(0.5*(rho(i,j+1,k)+rho(i,j,k))) + &
        (1.0 -rro/(0.5*(rho(i,j+1,k)+rho(i,j,k))) )*gy
  enddo; enddo; enddo

  do k=ks,kew;  do j=js,je; do i=is,ie
    dw(i,j,k)=dw(i,j,k)+(fz(i,j,k)-dpdz)/(0.5*(rho(i,j,k+1)+rho(i,j,k))) + &
        (1.0 -rro/(0.5*(rho(i,j,k+1)+rho(i,j,k))) )*gz
  enddo; enddo; enddo
end subroutine volumeForce
!=================================================================================================
!=================================================================================================
! The Poisson equation for the pressure is setup with matrix A as
! A7*Pijk = A1*Pi-1jk + A2*Pi+1jk + A3*Pij-1k + 
!           A4*Pij+1k + A5*Pijk-1 + A6*Pijk+1 + A8
!-------------------------------------------------------------------------------------------------
subroutine SetupPoisson(utmp,vtmp,wtmp,rhot,dt,A,mask)
  use module_grid
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: utmp,vtmp,wtmp,rhot
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax,8), intent(out) :: A
  logical, intent(in) :: mask(imin:imax,jmin:jmax,kmin:kmax)
  real(8) :: dt
  integer :: i,j,k
  do k=ks,ke; do j=js,je; do i=is,ie;
    if(mask(i,j,k))then
      A(i,j,k,1) = (1./dx(i))*( 1./(dxh(i-1)*(rhot(i-1,j,k)+rhot(i,j,k))) )
      A(i,j,k,2) = (1./dx(i))*( 1./(dxh(i  )*(rhot(i+1,j,k)+rhot(i,j,k))) )
      A(i,j,k,3) = (1./dy(j))*( 1./(dyh(j-1)*(rhot(i,j-1,k)+rhot(i,j,k))) )
      A(i,j,k,4) = (1./dy(j))*( 1./(dyh(j  )*(rhot(i,j+1,k)+rhot(i,j,k))) )
      A(i,j,k,5) = (1./dz(k))*( 1./(dzh(k-1)*(rhot(i,j,k-1)+rhot(i,j,k))) )
      A(i,j,k,6) = (1./dz(k))*( 1./(dzh(k  )*(rhot(i,j,k+1)+rhot(i,j,k))) )
      A(i,j,k,7) = sum(A(i,j,k,1:6))
      A(i,j,k,8) = -0.5/dt *( (utmp(i,j,k)-utmp(i-1,j,k))/dx(i) &
                             +(vtmp(i,j,k)-vtmp(i,j-1,k))/dy(j) &
                             +(wtmp(i,j,k)-wtmp(i,j,k-1))/dz(k)   )
    endif
  enddo; enddo; enddo
end subroutine SetupPoisson
!=================================================================================================
!=================================================================================================
! SOR sover to solve the following linear equiation:
! A7*Pijk = A1*Pi-1jk + A2*Pi+1jk + A3*Pij-1k + 
!           A4*Pij+1k + A5*Pijk-1 + A6*Pijk+1 + A8
!-------------------------------------------------------------------------------------------------
subroutine SOR_Solver(A,p,maxError,beta,maxit,it)
  use module_grid
  use module_BC
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: p
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax,8), intent(in) :: A
  real(8), intent(in) :: beta, maxError
  integer, intent(in) :: maxit
  integer, intent(out) :: it
  real(8) :: res, totalres
  integer :: i,j,k, ierr
  integer :: req(12),sta(MPI_STATUS_SIZE,12)
  logical :: mask(imin:imax,jmin:jmax,kmin:kmax)
!--------------------------------------ITTERATION LOOP--------------------------------------------  
  do it=1,maxit
    do k=ks,ke; do j=js,je; do i=is,ie
      p(i,j,k)=(1.0-beta)*p(i,j,k)+beta* 1.0/A(i,j,k,7)*(              &
        A(i,j,k,1) * p(i-1,j,k) + A(i,j,k,2) * p(i+1,j,k) +            &
        A(i,j,k,3) * p(i,j-1,k) + A(i,j,k,4) * p(i,j+1,k) +            &
        A(i,j,k,5) * p(i,j,k-1) + A(i,j,k,6) * p(i,j,k+1) + A(i,j,k,8))
    enddo; enddo; enddo
!---------------------------------CHECK FOR CONVERGENCE-------------------------------------------
    res = 0.0
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
        if(mask(i,j,k)) res=res+abs(-p(i,j,k) *  A(i,j,k,7) +            &
          A(i,j,k,1) * p(i-1,j,k) + A(i,j,k,2) * p(i+1,j,k) +            &
          A(i,j,k,3) * p(i,j-1,k) + A(i,j,k,4) * p(i,j+1,k) +            &
          A(i,j,k,5) * p(i,j,k-1) + A(i,j,k,6) * p(i,j,k+1) + A(i,j,k,8) )
      enddo; enddo; enddo
    call MPI_ALLREDUCE(res, totalres, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    totalres = totalres/float(Nx*Ny*Nz)
    if (.not.(totalres<1e10)) stop '***** solution has diverged *****'
    if (totalres<maxError) exit
  enddo
  if(it==maxit+1 .and. rank==0) write(*,*) 'Warning: SOR_Solver reached maxit.'
end subroutine SOR_Solver
!=================================================================================================
!=================================================================================================
! subroutine initialize
!   Initializes the flow solver
!   called in:    program ftc3d2011
!-------------------------------------------------------------------------------------------------
subroutine initialize
  use module_grid
  use module_flow
  use module_BC
  use module_IO
  use module_tmpvar
  implicit none
  include 'mpif.h'
  integer :: ierr, i,j,k

  allocate(dims(ndim),periodic(ndim),reorder(ndim),coords(ndim),STAT=ierr)
  dims(1) = nPx; dims(2) = nPy; dims(3) = nPz
  periodic = (bdry_cond=='periodic')
  reorder = .true.

  call MPI_Cart_Create(MPI_Comm_World,ndim,dims,periodic,reorder,MPI_Comm_Cart,ierr)
  if (ierr /= 0) STOP '*** Grid: unsuccessful Cartesian MPI-initialization'

  call MPI_Cart_Coords(MPI_Comm_Cart,rank,ndim,coords,ierr)
  if (ierr /= 0) STOP '*** Grid: unsuccessful in getting topology coordinates'

  is=coords(1)*Mx+1+Ng; ie=coords(1)*Mx+Mx+Ng; imin=is-Ng; imax=ie+Ng
  js=coords(2)*My+1+Ng; je=coords(2)*My+My+Ng; jmin=js-Ng; jmax=je+Ng
  ks=coords(3)*Mz+1+Ng; ke=coords(3)*Mz+Mz+Ng; kmin=ks-Ng; kmax=ke+Ng
  ieu=ie; if(bdry_cond(1)/='periodic' .and. coords(1)==nPx-1) ieu=ie-1
  jev=je; if(bdry_cond(2)/='periodic' .and. coords(2)==nPy-1) jev=je-1
  kew=ke; if(bdry_cond(3)/='periodic' .and. coords(3)==nPz-1) kew=ke-1

  Nxt=Nx+2*Ng; Nyt=Ny+2*Ng; Nzt=Nz+2*Ng ! total number of cells

!  allocate( x(imin:imax),dx(imin:imax),dxh(imin:imax),&
!            y(jmin:jmax),dy(jmin:jmax),dyh(jmin:jmax),&
!            z(kmin:kmax),dz(kmin:kmax),dzh(kmin:kmax) )
  allocate(x(Nxt), xh(Nxt), dx(Nxt), dxh(Nxt), &
           y(Nyt), yh(Nyt), dy(Nyt), dyh(Nyt), &
           z(Nzt), zh(Nzt), dz(Nzt), dzh(Nzt)  )

  allocate(  u(imin:imax,jmin:jmax,kmin:kmax), uold(imin:imax,jmin:jmax,kmin:kmax), &
            du(imin:imax,jmin:jmax,kmin:kmax),   fx(imin:imax,jmin:jmax,kmin:kmax), &
             v(imin:imax,jmin:jmax,kmin:kmax), vold(imin:imax,jmin:jmax,kmin:kmax), &
            dv(imin:imax,jmin:jmax,kmin:kmax),   fy(imin:imax,jmin:jmax,kmin:kmax), &
             w(imin:imax,jmin:jmax,kmin:kmax), wold(imin:imax,jmin:jmax,kmin:kmax), &
            dw(imin:imax,jmin:jmax,kmin:kmax),   fz(imin:imax,jmin:jmax,kmin:kmax), &
           rho(imin:imax,jmin:jmax,kmin:kmax), rhoo(imin:imax,jmin:jmax,kmin:kmax), &
          drho(imin:imax,jmin:jmax,kmin:kmax),    p(imin:imax,jmin:jmax,kmin:kmax), &
            mu(imin:imax,jmin:jmax,kmin:kmax),muold(imin:imax,jmin:jmax,kmin:kmax)  )

  allocate(tmp(imin:imax,jmin:jmax,kmin:kmax), work(imin:imax,jmin:jmax,kmin:kmax,8))

  x=0.0;y=0.0;z=0.0;dx=0.0;dy=0.0;dz=0.0;dxh=0.0;dyh=0.0;dzh=0.0;du=0.0;dv=0.0;dw=0.0
  u=0.0;v=0.0;w=0.0;p=0.0;tmp=0.0;fx=0.0;fy=0.0;fz=0.0;drho=0.0;rho=0.0;mu=0.0;work=00.

! Set the stretched grid
  do i=1,Nxt; s=dfloat(i-Ng)/dfloat(Nx); xh(i)=xLength*(xform*s*(0.5d0-s)*(1d0-s)+s); enddo
  do j=1,Nyt; s=dfloat(j-Ng)/dfloat(Ny); yh(j)=yLength*(yform*s*(0.5d0-s)*(1d0-s)+s); enddo
  do k=1,Nzt; s=dfloat(k-Ng)/dfloat(Nz); zh(k)=zLength*(zform*s*(0.5d0-s)*(1d0-s)+s); enddo

  if(read_x)then
    open(unit=12,file=trim(x_file),status='old',action='read')
  	read(12,*) (xh(i), i=Ng,Nx+Ng)
    close(12)
    xLength = xh(Nx+Ng)-xh(Ng)
    i=Nx+Ng+1
    xh(i) = 2d0*xh(i-1)-xh(i-2)
    if(Ng==2)xh(1) = 2d0*xh(2)-xh(3)
    if(Ng==2)xh(i+1) = 2d0*xh(i)-xh(i-1)
  endif
  if(read_y)then
    open(unit=12,file=trim(y_file),status='old',action='read')
  	read(12,*) (yh(i), i=Ng,Ny+Ng)
    close(12)
    yLength = yh(Ny+Ng)-yh(Ng)
    i=Ny+Ng+1
    yh(i) = 2d0*yh(i-1)-yh(i-2)
    if(Ng==2)yh(1) = 2d0*yh(2)-yh(3)
    if(Ng==2)yh(i+1) = 2d0*yh(i)-yh(i-1)
  endif
  if(read_z)then
    open(unit=12,file=trim(z_file),status='old',action='read')
  	read(12,*) (zh(i), i=Ng,Nz+Ng)
    close(12)
    zLength = zh(Nz+Ng)-zh(Ng)
    i=Nz+Ng+1
    zh(i) = 2d0*zh(i-1)-zh(i-2)
    if(Ng==2)zh(1) = 2d0*zh(2)-zh(3)
    if(Ng==2)zh(i+1) = 2d0*zh(i)-zh(i-1)
  endif

  do i=2,Nxt; x(i)=0.5d0*(xh(i)+xh(i-1)); enddo; x(1)=2d0*xh(1)-x(2)
  do j=2,Nyt; y(j)=0.5d0*(yh(j)+yh(j-1)); enddo; y(1)=2d0*yh(1)-y(2)
  do k=2,Nzt; z(k)=0.5d0*(zh(k)+zh(k-1)); enddo; z(1)=2d0*zh(1)-z(2)

  do i=1,Nxt-1; dxh(i)=x(i+1)-x(i); enddo; dxh(Nxt)=dxh(Nxt-1)
  do j=1,Nyt-1; dyh(j)=y(j+1)-y(j); enddo; dyh(Nyt)=dyh(Nyt-1)
  do k=1,Nzt-1; dzh(k)=z(k+1)-z(k); enddo; dzh(Nzt)=dzh(Nzt-1)
  
  do i=2,Nxt; dx(i)=xh(i)-xh(i-1); enddo; dx(1)=dx(2);
  do j=2,Nyt; dy(j)=yh(j)-yh(j-1); enddo; dy(1)=dy(2);
  do k=2,Nzt; dz(k)=zh(k)-zh(k-1); enddo; dz(1)=dz(2);

end subroutine initialize

!=================================================================================================
!=================================================================================================
! subroutine InitCondition
!   Sets the initial conditions
!   called in:    program ftc3d2011
!-------------------------------------------------------------------------------------------------
subroutine InitCondition
  use module_grid
  use module_flow
  use module_BC
  use module_IO
  use module_tmpvar
  implicit none
  include 'mpif.h'
  integer :: ierr, i,j,k, ib
  integer :: req(48),sta(MPI_STATUS_SIZE,48)

  ! Set the wall velocities
  WallVel=0.0

  if(restart)then
    call backup_read
    call ghost_x(u  ,2,req( 1: 4));  call ghost_x(v,2,req( 5: 8)); call ghost_x(w,2,req( 9:12)); 
    call ghost_x(rho,2,req(13:16));  call MPI_WAITALL(16,req(1:16),sta(:,1:16),ierr)
    call ghost_y(u  ,2,req( 1: 4));  call ghost_y(v,2,req( 5: 8)); call ghost_y(w,2,req( 9:12)); 
    call ghost_y(rho,2,req(13:16));  call MPI_WAITALL(16,req(1:16),sta(:,1:16),ierr)
    call ghost_z(u  ,2,req( 1: 4));  call ghost_z(v,2,req( 5: 8)); call ghost_z(w,2,req( 9:12)); 
    call ghost_z(rho,2,req(13:16));  call MPI_WAITALL(16,req(1:16),sta(:,1:16),ierr)
    mu = mu1
    if(TwoPhase) mu  = (rho-rho1)/(rho2-rho1)*(mu2-mu1)+mu1
  else
    time = 0d0
    itimestep = 0
    rho=rho1; mu=mu1
  ! Set density and viscosity in the domain and the drop
    if(TwoPhase)then
      do ib=1,numBubble
        do i=imin,imax; do j=jmin,jmax; do k=kmin,kmax; 
          if ( (x(i)-xc(ib))**2+(y(j)-yc(ib))**2+(z(k)-zc(ib))**2 < rad(ib)**2) then 
             rho(i,j,k)=rho2; mu(i,j,k)=mu2; 
          endif 
        enddo; enddo; enddo
      enddo
    endif
  endif
end subroutine InitCondition
!=================================================================================================
!=================================================================================================
! subroutine density
!   calculates drho/dt
!   called in:    program ftc3d2011
!-------------------------------------------------------------------------------------------------
subroutine density(rho,u,v,w,drho)
  use module_grid
  use module_tmpvar
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: u, v, w,rho !,D
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(out) :: drho
  real(8), external :: minabs
  integer :: i,j,k

  do k=ks,ke; do j=js,je; do i=is-1,ie
    if (u(i,j,k)>0.0) then
      work(i,j,k,1)=rho(i  ,j,k)+0.5*minabs((rho(i+1,j,k)-rho(i  ,j,k)),(rho(i  ,j,k)-rho(i-1,j,k)))
    else
      work(i,j,k,1)=rho(i+1,j,k)-0.5*minabs((rho(i+2,j,k)-rho(i+1,j,k)),(rho(i+1,j,k)-rho(i  ,j,k)))
    endif
  enddo; enddo; enddo
  do k=ks,ke; do j=js-1,je; do i=is,ie
    if(v(i,j,k)>0.0) then
      work(i,j,k,2)=rho(i,j  ,k)+0.5*minabs((rho(i,j+1,k)-rho(i,j  ,k)),(rho(i,j  ,k)-rho(i,j-1,k)))
    else
      work(i,j,k,2)=rho(i,j+1,k)-0.5*minabs((rho(i,j+2,k)-rho(i,j+1,k)),(rho(i,j+1,k)-rho(i,j  ,k)))
    endif
  enddo; enddo; enddo
  do k=ks-1,ke; do j=js,je; do i=is,ie
    if(w(i,j,k)>0.0) then
      work(i,j,k,3)=rho(i,j,k  )+0.5*minabs((rho(i,j,k+1)-rho(i,j,k  )),(rho(i,j,k  )-rho(i,j,k-1)))
    else
      work(i,j,k,3)=rho(i,j,k+1)-0.5*minabs((rho(i,j,k+2)-rho(i,j,k+1)),(rho(i,j,k+1)-rho(i,j,k  )))
    endif
  enddo; enddo; enddo
  do k=ks,ke;  do j=js,je; do i=is,ie
    drho(i,j,k) = -( u(i  ,j  ,k  )*work(i  ,j  ,k  ,1)          &
                    -u(i-1,j  ,k  )*work(i-1,j  ,k  ,1) )/dx(i)  &
                  -( v(i  ,j  ,k  )*work(i  ,j  ,k  ,2)          &
                    -v(i  ,j-1,k  )*work(i  ,j-1,k  ,2) )/dy(j)  &
                  -( w(i  ,j  ,k  )*work(i  ,j  ,k  ,3)          &
                    -w(i  ,j  ,k-1)*work(i  ,j  ,k-1,3) )/dz(k)  !&
!          + 0.5*( ( (D(i,j,k)+D(i+1,j,k))*(rho(i+1,j,k)-rho(i,j,k))/dxh(i  )         &
!                  - (D(i,j,k)+D(i-1,j,k))*(rho(i,j,k)-rho(i-1,j,k))/dxh(i-1) )/dx(i) &
!                 +( (D(i,j,k)+D(i,j+1,k))*(rho(i,j+1,k)-rho(i,j,k))/dyh(j  )         &
!                  - (D(i,j,k)+D(i,j-1,k))*(rho(i,j,k)-rho(i,j-1,k))/dyh(j-1) )/dy(j) &
!                 +( (D(i,j,k)+D(i,j,k+1))*(rho(i,j,k+1)-rho(i,j,k))/dzh(k  )         &
!                  - (D(i,j,k)+D(i,j,k-1))*(rho(i,j,k)-rho(i,j,k-1))/dzh(k-1) )/dz(k) )
  enddo; enddo; enddo

end subroutine density
!=================================================================================================
!=================================================================================================
! function minabs
!   used for ENO interpolations
!   called in:    subroutine momentumConvection
!                 subroutine density
!-------------------------------------------------------------------------------------------------
function minabs(a,b)
implicit none
real(8) :: minabs, a, b
  if(abs(a)<abs(b)) then
    minabs=a
  else
    minabs=b
  endif
!  minabs = 0.5*(sign(1.0,abs(b)-abs(a))*(a-b)+a+b)
end function minabs
!=================================================================================================
!=================================================================================================
! subroutine ReadParameters
!   Reads the parameters in "input" file for the calculation
!   called in:    program ftc3d2011
!-------------------------------------------------------------------------------------------------
subroutine ReadParameters
  use module_grid
  use module_flow
  use module_BC
  use module_IO
  implicit none
  include 'mpif.h'
  integer :: in, ierr
  integer, parameter :: MaxFront=1000
  real(8) :: xyzrad(4,MaxFront)
  namelist /parameters/ out_path , Nx, Ny, Nz, Ng, xLength, yLength, zLength, gx, gy, gz, hypre, &
                        bdry_cond, dpdx, dpdy, dpdz, itime_scheme, nstep, maxit, maxError,  &
                        beta, nout, TwoPhase, rho1, mu1, rho2, mu2, sigma,BuoyancyCase, &
                        nPx,nPy,nPz, EndTime, MaxDt, CFL, xform, yform, zform, output_format, &
                        nbackup, x_file, y_file, z_file, read_x, read_y, read_z, restart, dt, &
                        dtFlag, ICOut, NumBubble, xyzrad

  in=1
  out=2

  open(unit=in, file='input', status='old', action='read', iostat=ierr)
  if (ierr .ne. 0) stop 'ReadParameters: error openning input file'
  read(UNIT=in,NML=parameters)
  close(in)

  if(numBubble>MaxFront) stop 'Error: ReadParameters: increase size of xyzrad array (MaxFront)'
  allocate(rad(MaxFront), xc(MaxFront), yc(MaxFront), zc(MaxFront))
  xc (1:NumBubble) = xyzrad(1,1:NumBubble)
  yc (1:NumBubble) = xyzrad(2,1:NumBubble)
  zc (1:NumBubble) = xyzrad(3,1:NumBubble)
  rad(1:NumBubble) = xyzrad(4,1:NumBubble)

  if(rank==0)then
    call system('mkdir            '//trim(out_path))
    call system('cp input         '//trim(out_path))
    !call system('cp ftc3d2011.f90 '//trim(out_path))
    open(unit=out, file=trim(out_path)//'/output', action='write', iostat=ierr)
    if (ierr .ne. 0) stop 'ReadParameters: error openning output file'
    write(UNIT=out,NML=parameters)
  endif

  ! Number of grid points in streamwise and transverse directions must
  ! be integral multiples of total number of processors
  if(mod(Nx,nPx) /= 0) Stop 'ReadParameters: Nx not divisible by nPx!'
  if(mod(Ny,nPy) /= 0) Stop 'ReadParameters: Ny not divisible by nPy!'
  if(mod(Nz,nPz) /= 0) Stop 'ReadParameters: Nz not divisible by nPz!'
  Mx = Nx/nPx; My = Ny/nPy; Mz = Nz/nPz
  if(numProcess /= nPx*nPy*nPz)then
    if(rank==0) write(out,*) 'Error: incorrect number of processors.'
    call MPI_Finalize(ierr)
    stop
  endif
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
end subroutine ReadParameters
!=================================================================================================
!=================================================================================================
