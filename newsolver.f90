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
subroutine NewSolver(A,p,maxError,beta,maxit,it,ierr,norm)
  use module_mgsolver
  use module_grid
  use module_BC
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: p
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(in) :: A
  real(8), intent(in) :: beta, maxError
  real(8) :: tres2
  integer, intent(in) :: maxit
  integer, intent(out) :: it, ierr
  integer, intent(in) :: norm

  if (MultiGrid) then
      call NewSolverMG(A,p,maxError,beta,maxit,it,ierr,norm,tres2)
  else
      call NewSolver_std(A,p,maxError,beta,maxit,it,ierr,norm,tres2)
  endif

  if (test_MG.AND.rank==0) print *, 'residual', tres2
  if(it==maxit+1 .and. rank==0) write(*,*) 'Warning: LinearSolver reached maxit: ||res||: ',tres2

end subroutine NewSolver

subroutine relax_step(A,p,beta,L)
  use module_grid
  use module_BC
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: p
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(in) :: A
  real(8), intent(in) :: beta
  integer :: L
  integer, parameter :: relaxtype=1

  if(relaxtype==2) then 
      call LineRelax(A,p,beta,L)
  elseif(relaxtype==1) then
      call RedBlackRelax(A,p,beta,L)
  endif

end subroutine relax_step

subroutine NewSolver_std(A,p,maxError,beta,maxit,it,ierr,norm,tres2)
  use module_mgsolver
  use module_grid
  use module_BC
  use module_IO
  use module_freesurface
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: p
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(in) :: A
  real(8), intent(in)  :: beta, maxError
  real(8), intent(out) :: tres2
  integer, intent(in) :: maxit
  integer, intent(out) :: it, ierr
  real(8) :: res1,res2,resinf,intvol
  integer :: i,j,k
  integer :: req(12),sta(MPI_STATUS_SIZE,12)
  logical :: mask(imin:imax,jmin:jmax,kmin:kmax)
  integer, intent(in) :: norm
  integer, save :: itime=0
! Open file for convergence history
  if(rank==0.and.recordconvergence) then
     OPEN(UNIT=89,FILE=TRIM(out_path)//'/convergence_history-'//TRIM(int2text(itime,padding))//'.txt')
  endif
  itime=itime+1
  !--------------------------------------ITERATION LOOP--------------------------------------------  
  do it=1,maxit

      call relax_step(A,p,beta,-1)
!---------------------------------CHECK FOR CONVERGENCE-------------------------------------------
     res1 = 0d0; res2=0.d0; resinf=0.d0; intvol=0.d0
     call ghost_x(p,1,req( 1: 4)); call ghost_y(p,1,req( 5: 8)); call ghost_z(p,1,req( 9:12))
     do k=ks+1,ke-1; do j=js+1,je-1; do i=is+1,ie-1
        res2=res2+abs(-p(i,j,k) * A(i,j,k,7) +                           &
             A(i,j,k,1) * p(i-1,j,k) + A(i,j,k,2) * p(i+1,j,k) +            &
             A(i,j,k,3) * p(i,j-1,k) + A(i,j,k,4) * p(i,j+1,k) +            &
             A(i,j,k,5) * p(i,j,k-1) + A(i,j,k,6) * p(i,j,k+1) + A(i,j,k,8) )**norm 
     enddo; enddo; enddo
     call MPI_WAITALL(12,req,sta,ierr)
    mask=.true.
    mask(is+1:ie-1,js+1:je-1,ks+1:ke-1)=.false.
    do k=ks,ke; do j=js,je; do i=is,ie
       if(mask(i,j,k)) res2=res2+abs(-p(i,j,k) * A(i,j,k,7) +            &
            A(i,j,k,1) * p(i-1,j,k) + A(i,j,k,2) * p(i+1,j,k) +            &
            A(i,j,k,3) * p(i,j-1,k) + A(i,j,k,4) * p(i,j+1,k) +            &
            A(i,j,k,5) * p(i,j,k-1) + A(i,j,k,6) * p(i,j,k+1) + A(i,j,k,8) )**norm
    enddo; enddo; enddo
    res2 = res2/dble(Nx*Ny*Nz)
    call catch_divergence(res2,ierr)
    call MPI_ALLREDUCE(res2, tres2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr)
    if(norm==2) tres2=sqrt(tres2)
    if(rank==0.and.mod(it,10) == 0.and.recordconvergence) write(89,310) it, tres2
310 format(I6,'  ',(e14.5))
    if (tres2<maxError) then 
       if(rank==0.and.recordconvergence) close(89)
       exit
    endif
  enddo
  if(rank==0.and.recordconvergence) close(89)
contains
  subroutine catch_divergence(res2,ierr)
    real(8), intent(in) :: res2
    integer, intent(out) :: ierr
    logical :: diverged=.false.
    logical :: extended=.true.
    integer :: l
    if(extended) then
       do k=ks,ke; do j=js,je; do i=is,ie
          do l=1,8
             if(A(i,j,k,l)/=A(i,j,k,l).or.p(i,j,k)/=p(i,j,k)) then
                diverged=.true.
             endif
          enddo
          if(diverged) then
             OPEN(UNIT=88,FILE=TRIM(out_path)//'/message-rank-'//TRIM(int2text(rank,padding))//'.txt')
             write(88,*) "ijk rank",i,j,k,rank
             write(88,*) "A",  A(i,j,k,:)
             write(88,*) "p",  p(i,j,k)
             write(88,*) 'A or p is NaN after',it,'iterations at rank ',rank
             close(88)
             if(rank<=30) print*,'A or p is NaN after',it,'iterations at rank ',rank
             call pariserror("A or p is NaN")
             exit
          endif
       end do; end do; end do
    endif
    if ((res2*npx*npy*npz)>1.d16 ) then
       if(rank<=30) print*,'Pressure solver diverged after',it,'iterations at rank ',rank
       call pariserror("newsolver error")
    else if (res2 .ne. res2) then 
       if(rank<=30) print*, 'it:',it,'Pressure residual value is invalid at rank', rank
       call pariserror("newsolver error")
    else
       ierr=0
    endif
  end subroutine catch_divergence
end subroutine NewSolver_std

subroutine apply_BC(p)
  use module_grid
  use module_BC
  implicit none
  include 'mpif.h'
  real(8) :: p(imin:imax,jmin:jmax,kmin:kmax)
  integer :: L,ierr
  integer :: req(12),sta(MPI_STATUS_SIZE,12)

  call ghost_x(p,1,req( 1: 4))
  call ghost_y(p,1,req( 5: 8))
  call ghost_z(p,1,req( 9:12))
  call MPI_WAITALL(12,req,sta,ierr)

end subroutine

subroutine apply_BC_MG(p,L)
  use module_grid
  use module_BC
  implicit none
  include 'mpif.h'
  real(8) :: p(imin:imax,jmin:jmax,kmin:kmax)
  integer :: L,ierr
  integer :: req(12),sta(MPI_STATUS_SIZE,12)

  call ghost_MG_x(p,1,req( 1: 4),L)
  call ghost_MG_y(p,1,req( 5: 8),L)
  call ghost_MG_z(p,1,req( 9:12),L)
  call MPI_WAITALL(12,req,sta,ierr)

end subroutine

!--------------------------------------ONE RELAXATION ITERATION (SMOOTHER)----------------------------------
subroutine RedBlackRelax(A,p,beta,L)
  use module_grid
  use module_BC
  use module_IO
  use module_freesurface
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: p
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(in) :: A
  real(8), intent(in) :: beta
  integer :: L
  integer :: i,j,k,isw,jsw,ksw,ipass

  ksw=1
  do ipass=1,2
     jsw=ksw
     do k=ks,ke
        isw=jsw
        do j=js,je
           do i=isw+is-1,ie,2
              p(i,j,k)=(1d0-beta)*p(i,j,k) + beta/A(i,j,k,7)*(             &
                   A(i,j,k,1) * p(i-1,j,k) + A(i,j,k,2) * p(i+1,j,k) +       &
                   A(i,j,k,3) * p(i,j-1,k) + A(i,j,k,4) * p(i,j+1,k) +       &
                   A(i,j,k,5) * p(i,j,k-1) + A(i,j,k,6) * p(i,j,k+1) + A(i,j,k,8))
           enddo
           isw=3-isw
        enddo
        jsw=3-jsw
     enddo
     ksw=3-ksw
     if (ipass==1) then
       if (L.GT.0) then
           call apply_BC_MG(p,L)
       else
           call apply_BC(p)
       endif
     endif
  enddo
end subroutine RedBlackRelax
!--------------------------------------ONE RELAXATION ITERATION (SMOOTHER)----------------------------------
subroutine LineRelax(A,p,beta,L)
  use module_grid
  use module_freesurface
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: p
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(in) :: A
  real(8), intent(in) :: beta
  integer :: i,j,k
  integer :: L
!--------------------------------------ITERATION LOOP--------------------------------------------  
  do k=ks,ke; do j=js,je; do i=is,ie
     p(i,j,k)=(1d0-beta)*p(i,j,k) + (beta/A(i,j,k,7))*(              &
          A(i,j,k,1) * p(i-1,j,k) + A(i,j,k,2) * p(i+1,j,k) +        &
          A(i,j,k,3) * p(i,j-1,k) + A(i,j,k,4) * p(i,j+1,k) +        &
          A(i,j,k,5) * p(i,j,k-1) + A(i,j,k,6) * p(i,j,k+1) + A(i,j,k,8))
  enddo; enddo; enddo

  if (L.GT.0) then
      call apply_BC_MG(p,L)
  else
      call apply_BC(p)
  endif
end subroutine LineRelax
!=================================================================================================
!---------------------------------CHECK FOR WELL POSEDNESS --------------------------------------  
!=================================================================================================
subroutine check_poisson_setup(A,tmp,umask,vmask,wmask)
  use module_grid
  use module_BC
  use module_IO
  use module_freesurface
  implicit none
  include 'mpif.h'
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(in) :: A
  real(8), dimension(is:ie,js:je,ks:ke), intent(out) :: tmp  ! dims: 
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: umask,vmask,wmask

  ! array used only for temp computations sp redimensioning to smaller size is allowed
  integer :: ierr
  real(8) :: intsource,maxerr,A7
  real(8) :: tintsource,maxtmp,mintmp
  integer :: i,j,k,l
  logical :: diverged=.false.

! verify that we shall not divide by zero

  do k=ks,ke; do j=js,je; do i=is,ie
     do l=1,8
        if(A(i,j,k,l)/=A(i,j,k,l).or.A(i,j,k,7).lt.1d-55) then
           diverged=.true.
        endif
     enddo
     if(diverged) then
        OPEN(UNIT=88,FILE=TRIM(out_path)//'/message-rank-'//TRIM(int2text(rank,padding))//'.txt')
        write(88,*) "ijk rank",i,j,k,rank
        write(88,*) "A",  A(i, j, k,:)
        write(88,*) 'A trouble after setup at rank ',rank
        close(88)
        call pariserror("A is NaN or A7 is too close to 0")
        exit
     endif
  end do; end do; end do

! verify that the rhs is orthogonal to the kernel of the adjoint

  maxerr=1d-10
  intsource = sum(A(:,:,:,8))
  call MPI_ALLREDUCE(intsource, tintsource, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr)
  tintsource=tintsource/dble(Nx*Ny*Nz)
!   if(rank==0) then 
!      print *, "tintsource=",tintsource,"ierrmpi=",ierr
!  endif
  if (tintsource>maxerr) call pariserror("large volume source")

! Verify that 1 is in the kernel of the adjoint in the non solid regions. 

  tmp=0d0
  maxtmp=-1d50
  mintmp=1d50
  do k=ks,ke; do j=js,je; do i=is,ie
     A7 = A(i,j,k,7)
     if(umask(i-1,j,k).lt.0.5d0.and.umask(i,j,k).lt.0.5d0.and.     &
          vmask(i,j-1,k).lt.0.5d0.and.vmask(i,j,k).lt.0.5d0.and.   &
          wmask(i,j,k-1).lt.0.5d0.and.wmask(i,j,k).lt.0.5d0 ) then 
        ! we are in the solid or in an isolated fluid cell. 
        A7 = 0d0
     endif
     tmp(i,j,k) = sum(A(i,j,k,1:6)) - A7
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
subroutine check_corrected_vel(u,umask,iout)
  use module_grid
  
  use module_BC
  use module_IO
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: u
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: umask
  integer, intent(in) :: iout
  integer :: i,j,k,ierr,i1,i2
  real(8), dimension(:), allocatable, save :: flux,fullflux,fluxm,fullfluxm
  logical :: initialized = .false.
  if(.not. initialized) then 
     allocate(flux(nx+2*Ng),fullflux(nx+2*Ng))
     allocate(fluxm(nx+2*Ng),fullfluxm(nx+2*Ng))
     initialized = .true.
  endif
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

