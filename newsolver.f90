!=================================================================================================
!=================================================================================================
!=================================================================================================
! Paris-0.1
! Extended from Code: FTC3D2011 (Front Tracking Code for 3D simulations)
! 
! Authors: Sadegh Dabiri, Gretar Tryggvason
! author for this file: Stephane Zaleski(zaleski@dalembert.upmc.fr) 
! 
! A three-dimensional Navier-Stokes flow solver for modeling of multiphase 
! flows. Flow can be driven by wall motion, density difference or pressure gradient.
! Boundary conditions supported: wall and periodic
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
! Solves the following linear equiation:
! A7*Pijk = A1*Pi-1jk + A2*Pi+1jk + A3*Pij-1k + 
!           A4*Pij+1k + A5*Pijk-1 + A6*Pijk+1 + A8
!-------------------------------------------------------------------------------------------------
subroutine NewSolver(A,p,maxError,beta,maxit,it,ierr)
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
  real(8) :: res1,res2,resinf,intvol,intsource,maxerrsrc
  real(8) :: tres1,tres2,tresinf,tintvol,tintsource
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
    res1 = 0d0; res2=0.d0; resinf=0.d0; intvol=0.d0
    call ghost_x(p,1,req( 1: 4)); call ghost_y(p,1,req( 5: 8)); call ghost_z(p,1,req( 9:12))
    do k=ks+1,ke-1; do j=js+1,je-1; do i=is+1,ie-1
      res2=res2+abs(-p(i,j,k) * A(i,j,k,7) +                           &
        A(i,j,k,1) * p(i-1,j,k) + A(i,j,k,2) * p(i+1,j,k) +            &
        A(i,j,k,3) * p(i,j-1,k) + A(i,j,k,4) * p(i,j+1,k) +            &
        A(i,j,k,5) * p(i,j,k-1) + A(i,j,k,6) * p(i,j,k+1) + A(i,j,k,8) )**2
    enddo; enddo; enddo
    !write(*,'(e14.5)')res2
    !write(*,'(2i4)')ks,ke
    call MPI_WAITALL(12,req,sta,ierr)
    mask=.true.
    mask(is+1:ie-1,js+1:je-1,ks+1:ke-1)=.false.
    do k=ks,ke; do j=js,je; do i=is,ie
      if(mask(i,j,k))res2=res2+abs(-p(i,j,k) * A(i,j,k,7) +            &
        A(i,j,k,1) * p(i-1,j,k) + A(i,j,k,2) * p(i+1,j,k) +            &
        A(i,j,k,3) * p(i,j-1,k) + A(i,j,k,4) * p(i,j+1,k) +            &
        A(i,j,k,5) * p(i,j,k-1) + A(i,j,k,6) * p(i,j,k+1) + A(i,j,k,8) )**2
    enddo; enddo; enddo
    res2 = res2/dble(Nx*Ny*Nz)
    call catch_divergence(res2,ierr)
    call MPI_ALLREDUCE(res2, tres2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr)
    tres2=sqrt(tres2)
    if (tres2<maxError) exit
  enddo
  if(it==maxit+1 .and. rank==0) write(*,*) 'Warning: LinearSolver reached maxit: ||res||_2: ',tres2
contains
  subroutine catch_divergence(res2,ierr)
    real(8), intent(in) :: res2
    integer, intent(out) :: ierr
    if ((res2*npx*npy*npz)>1.d16 ) then
       if(rank==0) print*,'Pressure solver diverged after',it,'iterations at rank ',rank
       call pariserror("newsolver error")
    else if (res2 .ne. res2) then 
       if(rank==0) print*, 'it:',it,'Pressure residual value is invalid at rank', rank
       call pariserror("newsolver error")
    else
       ierr=0
    endif
  end subroutine catch_divergence
end subroutine NewSolver
!=================================================================================================
!=================================================================================================
!---------------------------------CHECK FOR WELL POSEDNESS --------------------------------------  
!=================================================================================================
subroutine check_poisson_setup(A,tmp)
  use module_grid
  use module_BC
  implicit none
  include 'mpif.h'
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(in) :: A
  real(8), dimension(is:ie,js:je,ks:ke), intent(out) :: tmp  ! dims: 
  ! array used only for temp computations sp redimensioning to smaller size is allowed
  integer :: ierr
  real(8) :: intsource,maxerr
  real(8) :: tintsource,maxtmp,mintmp
  integer :: i,j,k
! verify that the rhs is orthogonal to the kernel of the adjoint

  maxerr=1d-10
  intsource = sum(A(:,:,:,8))
  call MPI_ALLREDUCE(intsource, tintsource, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr)
  tintsource=tintsource/dble(Nx*Ny*Nz)
!   if(rank==0) then 
!      print *, "tintsource=",tintsource,"ierrmpi=",ierr
!  endif
  if (tintsource>maxerr) call pariserror("large volume source")

! Verify that 1 is in the kernel of the adjoint. 

  tmp=0d0
  maxtmp=-1d50
  mintmp=1d50
  do k=ks,ke; do j=js,je; do i=is,ie
     tmp(i,j,k) = sum(A(i,j,k,1:6)) - A(i,j,k,7)
     maxtmp=max(maxtmp, tmp(i,j,k))
     mintmp=min(mintmp, tmp(i,j,k))
  enddo; enddo; enddo
  if(max(maxtmp,-mintmp)>maxerr) then
     print *, "maxmin ",maxtmp,mintmp,rank
     call pariserror("constant p not in kernel")
  endif
end subroutine check_poisson_setup
! 
!=================================================================================================
!-------------------------------------------------------------------------------------------------
!
!  write a file with the volume flux in each x-plane
!
!-------------------------------------------------------------------------------------------------
!=================================================================================================
!
subroutine check_corrected_vel(u,v,w,umask,vmask,wmask,iout)
  use module_grid
  use module_hello
  use module_BC
  use module_IO
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: u,v,w
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: umask,vmask,wmask
  integer, intent(in) :: iout
  integer :: i,j,k,ierr,i1,i2
  real(8), dimension(:), allocatable, save :: flux,fullflux,fluxm,fullfluxm
  logical :: initialized = .false.
  if(.not. initialized) then 
     allocate(flux(nx+2*Ng),fullflux(nx+2*Ng))
     allocate(fluxm(nx+2*Ng),fullfluxm(nx+2*Ng))
     initialized = .true.
  endif

!  print*,"xh(is-1),xh(ie)",xh(Ng),xh(Nx+Ng)


  flux = 0d0; fluxm = 0d0
  i1=is
  i2=ieu
  if(coords(1)==0) i1=1
  if(coords(1)==Npx-1) i2=Nx+2*Ng
  do i=i1,i2
     do k=ks,ke
        do j=js,je
           flux(i) = flux(i) + u(i,j,k)*dz(k)*dy(j)
           fluxm(i) = fluxm(i) + umask(i,j,k)*u(i,j,k)*dz(k)*dy(j)
        enddo
     enddo
  enddo
  call MPI_REDUCE(flux(1), fullflux(1),nx+2*Ng, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_Domain, ierr)
  call MPI_REDUCE(fluxm(1), fullfluxm(1),nx+2*Ng, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_Domain, ierr)
  if(rank==0) then
     OPEN(UNIT=89,FILE=TRIM(out_path)//'/volume_flux-'//TRIM(int2text(iout,padding))//'.txt')
     do i=1,Nx+2*Ng
        write(89,310) i, fullflux(i)/(Ylength*Zlength), xh(i) , fullfluxm(i)/(Ylength*Zlength)
     enddo
     close(89)
  endif
310 format(I3,'  ',4(e14.5))
end subroutine check_corrected_vel
