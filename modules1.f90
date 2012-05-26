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
