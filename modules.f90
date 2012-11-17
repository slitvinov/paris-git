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
!=================================================================================================
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
  real(8) :: xLength, yLength, zLength, xform, yform, zform !non-uniformity of grid

  integer :: nPx, nPy, nPz, Mx, My, Mz, rank, ndim=3, numProcess
  integer, dimension(:), allocatable :: dims, coords, periodic, reorder
  integer :: MPI_Comm_Cart
  integer :: imin, imax, jmin, jmax, kmin, kmax
  logical :: hypre
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
  if(rank==0) write(6,*) 'coucou ',hello_count
  hello_count = hello_count + 1
end subroutine hello_coucou
end module module_hello
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
  interface
     SUBROUTINE append_visit_file(rootname,padding)
       character :: rootname(*)
       integer :: padding
     END SUBROUTINE append_visit_file
  end interface
  contains
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
  integer :: padding=5
  OPEN(UNIT=7,FILE=trim(out_path)//'/backup_'//int2text(rank,padding),status='unknown',action='write')
  write(7,1100)time,itimestep,is,ie,js,je,ks,ke
  do k=ks,ke; do j=js,je; do i=is,ie
    write(7,1200) u(i,j,k), v(i,j,k), w(i,j,k), p(i,j,k), rho(i,j,k)
  enddo; enddo; enddo
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
  implicit none
  integer ::i,j,k,i1,i2,j1,j2,k1,k2
  integer :: padding=5
  OPEN(UNIT=7,FILE=trim(out_path)//'/backup_'//int2text(rank,padding),status='old',action='read')
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
!  use module_tmpvar
  !use IO_mod
  implicit none
  integer ::nf,i1,i2,j1,j2,k1,k2,i,j,k
  integer :: padding=5
!  logical, save :: first_time=.true.

  OPEN(UNIT=7,FILE=trim(out_path)//'/plot'//int2text(nf,padding)//'_'//int2text(rank,3)//'.dat')
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
  1110 FORMAT('ZONE solutiontime=', 1PG15.7e2)
  1200 FORMAT(21es14.6e2)
end subroutine output1
!-------------------------------------------------------------------------------------------------
subroutine output2(nf,i1,i2,j1,j2,k1,k2)
  use module_flow
  use module_grid
  !use IO_mod
  implicit none
  integer ::nf,i1,i2,j1,j2,k1,k2,i,j,k, itype=5
!  logical, save :: first_time=.true.
  character(len=30) :: rootname
  integer :: padding=5
  rootname=TRIM(out_path)//'/plot'//TRIM(int2text(nf,padding))//'-'

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

    if(ngh>Ng) stop 'ghost error: not enough ghost layers to fill'
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

    if(ngh>Ng) stop 'ghost error: not enough ghost layers to fill'
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

    if(ngh>Ng) stop 'ghost error: not enough ghost layers to fill'
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
  use module_hello
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
