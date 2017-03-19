!================================================================================================
!=================================================================================================
! Paris-0.1
!
! Multigrid solver
! written by Daniel Fuster
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

MODULE module_mgsolver

  TYPE equation
    real(8), ALLOCATABLE :: K(:,:,:,:)
  END TYPE equation
  TYPE(equation), ALLOCATABLE :: DataMG(:), pMG(:), EMG(:)
  integer :: nd(3), LMG

CONTAINS

subroutine init_MG(nvar)
  use module_grid
  use module_BC
  implicit none
  include 'mpif.h'
  integer, intent(in) :: nvar
  integer :: Ld(3), nL(3)
  integer :: i,ierr

  nd(1) = ie-is+1; nd(2) = je-js+1; nd(3) = ke-ks+1
  IF (IAND (nd(1), nd(1)-1).NE.0) STOP 'multigrid requires 2^n+1 nodes per direction'
  IF (IAND (nd(2), nd(2)-1).NE.0) STOP 'multigrid requires 2^n+1 nodes per direction'
  IF (IAND (nd(3), nd(3)-1).NE.0) STOP 'multigrid requires 2^n+1 nodes per direction'

  Ld(:) = int(log(float(nd(:)))/log(2.d0))
  LMG = minval(Ld(:))
  call MPI_ALLREDUCE(LMG, LMG, 1, MPI_INTEGER, MPI_MIN, MPI_Comm_Cart, ierr) 

  ALLOCATE (DataMG(LMG)); ALLOCATE (pMG(LMG)); ALLOCATE (EMG(LMG))

  DO i=LMG,1,-1
      nL(:) = nd(:)/2**(LMG-i)
      ALLOCATE(DataMG(i)%K(nL(1),nL(2),nL(3),8))
      nL(:) = nd(:)/2**(LMG-i)+2*Ng
      ALLOCATE(pMG(i)%K(nL(1),nL(2),nL(3),nvar)) 
      ALLOCATE(EMG(i)%K(nL(1),nL(2),nL(3),nvar)) 
      DataMG(i)%K = 0.d0
      pMG(i)%K    = 0.d0
      EMG(i)%K    = 0.d0
  enddo

end subroutine init_MG

subroutine finalize_MG
  implicit none
  include 'mpif.h'
  integer :: i

  !free memory
  DO i=1,LMG
      DEALLOCATE(DataMG(i)%K) 
      DEALLOCATE(pMG(i)%K) 
      DEALLOCATE(EMG(i)%K) 
  ENDDO
  DEALLOCATE(DataMG) 
  DEALLOCATE(pMG) 
  DEALLOCATE(EMG) 

end subroutine finalize_MG
 
subroutine get_residual(A,p,norm,error_norm)
  use module_grid
  use module_BC
  implicit none
  include 'mpif.h'
  real(8) :: A(is:ie,js:je,ks:ke,8), p(imin:imax,jmin:jmax,kmin:kmax)
  real(8), intent(out) :: error_norm
  integer, intent(in)  :: norm
  integer :: i,j,k

  error_norm = 0.d0

  do k=ks,ke; do j=js,je; do i=is,ie
      A(i,j,k,8) = A(i,j,k,8)     - A(i,j,k,7)*p(i,j,k)    + &
        A(i,j,k,1)*p(i-1,j,k) + A(i,j,k,2)*p(i+1,j,k)  + &
        A(i,j,k,3)*p(i,j-1,k) + A(i,j,k,4)*p(i,j+1,k)  + &
        A(i,j,k,5)*p(i,j,k-1) + A(i,j,k,6)*p(i,j,k+1) 
    error_norm = error_norm + abs(A(i,j,k,8))**norm
  enddo; enddo; enddo

end subroutine get_residual

subroutine compute_residual(A,p,norm,error_norm)
  use module_grid
  use module_BC
  implicit none
  include 'mpif.h'
  real(8) :: A(is:ie,js:je,ks:ke,8), p(imin:imax,jmin:jmax,kmin:kmax)
  real(8), intent(out) :: error_norm
  integer, intent(in)  :: norm
  integer :: i,j,k

  error_norm = 0.d0

  do k=ks,ke; do j=js,je; do i=is,ie
      error_norm = error_norm + abs(A(i,j,k,8)     - &
        A(i,j,k,7)*p(i,j,k)    + &
        A(i,j,k,1)*p(i-1,j,k) + A(i,j,k,2)*p(i+1,j,k)  + &
        A(i,j,k,3)*p(i,j-1,k) + A(i,j,k,4)*p(i,j+1,k)  + &
        A(i,j,k,5)*p(i,j,k-1) + A(i,j,k,6)*p(i,j,k+1) )**norm
  enddo; enddo; enddo

end subroutine compute_residual

subroutine coarse_fine_interp(pf,pc,imin1,jmin1,kmin1)
  use module_grid
  use module_BC
    implicit none
    integer :: imin1,jmin1,kmin1
    real(8) :: pf(imin1:,jmin1:,kmin1:,:), pc(imin:,jmin:,kmin:,:)
    real(8) :: residtot
    integer :: i,j,k,is1,js1,ks1
    
    is1=imin1+Ng; js1=jmin1+Ng; ks1=kmin1+Ng
    !interpolation from coarse level

  IF (sym_MG) THEN
    do i=is,ie; do j=js,je; do k=ks,ke
    !symmetric operator
      residtot = abs(pf(2*(i-is)+is1,  2*(j-js)+js1,  2*(k-ks)+ks1,1))   + & 
                 abs(pf(2*(i-is)+is1+1,2*(j-js)+js1,  2*(k-ks)+ks1,1))   + & 
                 abs(pf(2*(i-is)+is1,  2*(j-js)+js1+1,2*(k-ks)+ks1,1))   + & 
                 abs(pf(2*(i-is)+is1+1,2*(j-js)+js1+1,2*(k-ks)+ks1,1))  + & 
                 abs(pf(2*(i-is)+is1,  2*(j-js)+js1,  2*(k-ks)+ks1+1,1)) + & 
                 abs(pf(2*(i-is)+is1+1,2*(j-js)+js1,  2*(k-ks)+ks1+1,1)) + & 
                 abs(pf(2*(i-is)+is1,  2*(j-js)+js1+1,2*(k-ks)+ks1+1,1)) + & 
                 abs(pf(2*(i-is)+is1+1,2*(j-js)+js1+1,2*(k-ks)+ks1+1,1)) + 1.d-15

      pf(2*(i-is)+is1,  2*(j-js)+js1,  2*(k-ks)+ks1,1)   = & 
          abs(pf(2*(i-is)+is1,  2*(j-js)+js1,  2*(k-ks)+ks1,1))/residtot * pc(i,j,k,1)
      pf(2*(i-is)+is1+1,2*(j-js)+js1,  2*(k-ks)+ks1,1)   = & 
          abs(pf(2*(i-is)+is1+1,2*(j-js)+js1,  2*(k-ks)+ks1,1))/residtot * pc(i,j,k,1)
      pf(2*(i-is)+is1,  2*(j-js)+js1+1,2*(k-ks)+ks1,1)   = & 
          abs(pf(2*(i-is)+is1,  2*(j-js)+js1+1,2*(k-ks)+ks1,1))/residtot * pc(i,j,k,1)
      pf(2*(i-is)+is1+1,2*(j-js)+js1+1,2*(k-ks)+ks1,1)   = & 
          abs(pf(2*(i-is)+is1+1,2*(j-js)+js1+1,2*(k-ks)+ks1,1))/residtot * pc(i,j,k,1)
      pf(2*(i-is)+is1,  2*(j-js)+js1,  2*(k-ks)+ks1+1,1) = & 
          abs(pf(2*(i-is)+is1,  2*(j-js)+js1,  2*(k-ks)+ks1+1,1))/residtot * pc(i,j,k,1)
      pf(2*(i-is)+is1+1,2*(j-js)+js1,  2*(k-ks)+ks1+1,1) = & 
          abs(pf(2*(i-is)+is1+1,2*(j-js)+js1,  2*(k-ks)+ks1+1,1))/residtot * pc(i,j,k,1)
      pf(2*(i-is)+is1,  2*(j-js)+js1+1,2*(k-ks)+ks1+1,1) = & 
          abs(pf(2*(i-is)+is1,  2*(j-js)+js1+1,2*(k-ks)+ks1+1,1))/residtot * pc(i,j,k,1)
      pf(2*(i-is)+is1+1,2*(j-js)+js1+1,2*(k-ks)+ks1+1,1) = & 
          abs(pf(2*(i-is)+is1+1,2*(j-js)+js1+1,2*(k-ks)+ks1+1,1))/residtot * pc(i,j,k,1)
    enddo; enddo; enddo
  ELSE
    do i=is,ie; do j=js,je; do k=ks,ke
    !original non-symmetric operator
      pf(2*(i-is)+is1,  2*(j-js)+js1,  2*(k-ks)+ks1,1)   = pc(i,j,k,1)
      pf(2*(i-is)+is1+1,2*(j-js)+js1,  2*(k-ks)+ks1,1)   = pc(i,j,k,1)
      pf(2*(i-is)+is1,  2*(j-js)+js1+1,2*(k-ks)+ks1,1)   = pc(i,j,k,1)
      pf(2*(i-is)+is1+1,2*(j-js)+js1+1,2*(k-ks)+ks1,1)   = pc(i,j,k,1)
      pf(2*(i-is)+is1,  2*(j-js)+js1,  2*(k-ks)+ks1+1,1) = pc(i,j,k,1)
      pf(2*(i-is)+is1+1,2*(j-js)+js1,  2*(k-ks)+ks1+1,1) = pc(i,j,k,1)
      pf(2*(i-is)+is1,  2*(j-js)+js1+1,2*(k-ks)+ks1+1,1) = pc(i,j,k,1)
      pf(2*(i-is)+is1+1,2*(j-js)+js1+1,2*(k-ks)+ks1+1,1) = pc(i,j,k,1)
    enddo; enddo; enddo
  ENDIF

end subroutine coarse_fine_interp

subroutine coarse_from_fine(Af,Ac,is1,js1,ks1)
  use module_grid
  use module_BC
    implicit none
    integer :: is1,js1,ks1
    integer :: i,j,k,i1,j1,k1
    real(8) :: Af(is1:,js1:,ks1:,:), Ac(is:,js:,ks:,:)

    Ac(:,:,:,8) = 0.d0
    DO i=is,ie; DO j=js,je; DO k=ks,ke
    DO i1=0,1; DO j1=0,1; DO k1=0,1
      Ac(i,j,k,8) = Ac(i,j,k,8) + &
      Af(2*(i-is)+is1+i1,  2*(j-js)+js1+j1,  2*(k-ks)+ks1+k1,8)*0.125d0
    ENDDO; ENDDO; ENDDO
    ENDDO; ENDDO; ENDDO

end subroutine coarse_from_fine

subroutine fill_coefficients(Af, Ac, is1, js1, ks1)
  use module_grid
  use module_BC
  implicit none
  integer :: i,j,k,is1,js1,ks1
  real(8) :: Ac(is:,js:,ks:,:), Af(is1:,js1:,ks1:,:)

    DO k=ks,ke; DO j=js,je; DO i=is,ie
        Ac(i,j,k,:) =  0.125d0*(                       &
                Af(2*(i-is)+is1,  2*(j-js)+js1,  2*(k-ks)+ks1,:) + &
                Af(2*(i-is)+is1,  2*(j-js)+js1+1,2*(k-ks)+ks1,:) + &
                Af(2*(i-is)+is1,  2*(j-js)+js1,  2*(k-ks)+ks1+1,:) + &
                Af(2*(i-is)+is1,  2*(j-js)+js1+1,2*(k-ks)+ks1+1,:) + &
                Af(2*(i-is)+is1+1,2*(j-js)+js1,  2*(k-ks)+ks1,:) + &
                Af(2*(i-is)+is1+1,2*(j-js)+js1+1,2*(k-ks)+ks1,:) + &
                Af(2*(i-is)+is1+1,2*(j-js)+js1,  2*(k-ks)+ks1+1,:) + &
                Af(2*(i-is)+is1+1,2*(j-js)+js1+1,2*(k-ks)+ks1+1,:))
     ENDDO; ENDDO; ENDDO

!    DO k=ks,ke; DO j=js,je; DO i=is,ie
!        Ac(i,j,k,1:2) =  0.25d0*(                       &
!                Af(2*(i-is)+is1,  2*(j-js)+js1,  2*(k-ks)+ks1,1:2) + &
!                Af(2*(i-is)+is1,  2*(j-js)+js1+1,2*(k-ks)+ks1,1:2) + &
!                Af(2*(i-is)+is1,  2*(j-js)+js1,  2*(k-ks)+ks1+1,1:2) + &
!                Af(2*(i-is)+is1,  2*(j-js)+js1+1,2*(k-ks)+ks1+1,1:2))
!     ENDDO; ENDDO; ENDDO
!
!    DO k=ks,ke; DO j=js,je; DO i=is,ie
!        Ac(i,j,k,3:4) = 0.25d0*(                       &
!                Af(2*(i-is)+is1,  2*(j-js)+js1,  2*(k-ks)+ks1,3:4) + &
!                Af(2*(i-is)+is1+1,2*(j-js)+js1,  2*(k-ks)+ks1,3:4) + &
!                Af(2*(i-is)+is1,  2*(j-js)+js1,  2*(k-ks)+ks1+1,3:4) + &
!                Af(2*(i-is)+is1+1,2*(j-js)+js1,  2*(k-ks)+ks1+1,3:4) )
!     ENDDO; ENDDO; ENDDO
!
!     DO k=ks,ke; DO j=js,je; DO i=is,ie
!        Ac(i,j,k,5:6) =  0.25d0*(                         &
!                Af(2*(i-is)+is1,  2*(j-js)+js1,  2*(k-ks)+ks1,5:6) + &
!                Af(2*(i-is)+is1+1,2*(j-js)+js1,  2*(k-ks)+ks1,5:6) + &
!                Af(2*(i-is)+is1,  2*(j-js)+js1+1,2*(k-ks)+ks1,5:6) + &
!                Af(2*(i-is)+is1+1,2*(j-js)+js1+1,2*(k-ks)+ks1,5:6) )
!     ENDDO; ENDDO; ENDDO
!
!     DO k=ks,ke; DO j=js,je; DO i=is,ie
!        Ac(i,j,k,7) = sum(Ac(i,j,k,1:6))
!     ENDDO; ENDDO; ENDDO


end subroutine fill_coefficients

subroutine MG_iteration(beta, norm)
  use module_grid
  use module_BC
  implicit none
  include 'mpif.h'
  real(8), intent(in) :: beta
  integer, intent(in) :: norm
  real(8) :: resMax, start_time, end_time
  integer :: i,j,k,n, Level
  integer :: is1, js1, ks1, imin1, jmin1, kmin1
  integer :: ncall=1
  integer :: nL(3)

  !-------------
  DO Level=LMG,2,-1

    pMG(Level)%K = 0.d0  !now pMG is the error
    nL(:) = nd(:)/2**(LMG-Level)
    call update_bounds(nL)
    
    DO i=1,nrelax
        call relax_step(DataMG(Level)%K,pMG(Level)%K(:,:,:,1),beta,Level)
        call apply_BC_MG(pMG(Level)%K(:,:,:,1),Level)
    ENDDO
   
    call get_residual(DataMG(Level)%K,pMG(Level)%K(:,:,:,1),norm,resMax)

    nL(:) = nd(:)/2**(LMG-Level+1)
    call update_bounds(nL)
    is1=coords(1)*nd(1)/2**(LMG-Level)+1+Ng
    js1=coords(2)*nd(2)/2**(LMG-Level)+1+Ng
    ks1=coords(3)*nd(3)/2**(LMG-Level)+1+Ng
    call coarse_from_fine(DataMG(Level)%K,DataMG(Level-1)%K,is1,js1,ks1) !projection at the next level

  ENDDO

  nL(:) = nd(:)/2**(LMG-2)
  call update_bounds(nL)
  call relax_step(DataMG(2)%K,pMG(2)%K(:,:,:,1),beta,2)
  call apply_BC_MG(pMG(2)%K(:,:,:,1),2)

  DO Level=3,LMG,1

    EMG(Level)%K = 0.d0  !initial guess of the error
    nL(:) = nd(:)/2**(LMG-Level+1)
    call update_bounds(nL)
    imin1=coords(1)*nd(1)/2**(LMG-Level)+1
    jmin1=coords(2)*nd(2)/2**(LMG-Level)+1
    kmin1=coords(3)*nd(3)/2**(LMG-Level)+1
    call coarse_fine_interp(EMG(Level)%K,pMG(Level-1)%K,imin1,jmin1,kmin1) !interpolation from coarse level
    nL(:) = nd(:)/2**(LMG-Level)
    call update_bounds(nL)
    call apply_BC_MG(EMG(Level)%K(:,:,:,1),Level)

    call get_residual(DataMG(Level)%K,EMG(Level)%K(:,:,:,1),norm,resMax) !new residual
    pMG(Level)%K = pMG(Level)%K + EMG(Level)%K

    EMG(Level)%K = 0.d0  !initial guess of the error
    DO i=1,nrelax
      call relax_step(DataMG(Level)%K,EMG(Level)%K(:,:,:,1),beta,Level)
      call apply_BC_MG(EMG(Level)%K(:,:,:,1),Level)
    ENDDO
    pMG(Level)%K = pMG(Level)%K + EMG(Level)%K

    call get_residual(DataMG(Level)%K,EMG(Level)%K(:,:,:,1),norm,resMax)
  ENDDO

end subroutine MG_iteration

subroutine get_MGmatrix_coef(A)
  use module_grid
  use module_BC
  implicit none
  include 'mpif.h'
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(in) :: A
  integer :: is1, js1, ks1, i, j, k, Level
  integer :: nL(3)

  do k=ks,ke; do j=js,je; do i=is,ie
      DataMG(LMG)%K(i-is+1,j-js+1,k-ks+1,:) = A(i,j,k,:) !fine level
  enddo; enddo; enddo

  DO Level=LMG-1,1,-1
    nL(:) = nd(:)/2**(LMG-Level)
    call update_bounds(nL)
    is1=coords(1)*nd(1)/2**(LMG-Level-1)+1+Ng
    js1=coords(2)*nd(2)/2**(LMG-Level-1)+1+Ng
    ks1=coords(3)*nd(3)/2**(LMG-Level-1)+1+Ng
    call fill_coefficients(DataMG(Level+1)%K,DataMG(Level)%K,is1,js1,ks1)
  ENDDO

end subroutine get_MGmatrix_coef

subroutine NewSolverMG(A,p,maxError,beta,maxit,it,ierr,norm,tres2)
  use module_grid
  use module_BC
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: p
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(in) :: A
  real(8), intent(in) :: beta, maxError
  real(8), intent(out) :: tres2
  integer, intent(in) :: maxit
  integer, intent(out) :: it, ierr
  real(8) :: resMax, end_time, start_time
  integer :: ncall=1
  integer :: nL(3)
  integer, intent(in) :: norm

  IF (rank==0.and.recordconvergence) THEN
    OPEN(UNIT=89,FILE='convergenceMG_history.txt',position='append')
    start_time =  MPI_WTIME()
  ENDIF

  call get_MGmatrix_coef(A)

  !compute residual at the finest level
  call update_bounds(nd)
  pMG(LMG)%K(:,:,:,1)    = p
  call get_residual(DataMG(LMG)%K,pMG(LMG)%K(:,:,:,1),norm,resMax)

  resMax = resMax/dble(Nx*Ny*Nz)
  call MPI_ALLREDUCE(resMax, tres2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr) 

  pMG(LMG)%K(:,:,:,1)    = 0.d0

  DO it=1,maxit,1

    call MG_iteration(beta, norm)

    p = p + pMG(LMG)%K(:,:,:,1)
    call update_bounds(nd)
    call compute_residual(A,p,norm,resMax)

    resMax = resMax/dble(Nx*Ny*Nz)
    call MPI_ALLREDUCE(resMax, tres2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr) 
    if(norm==2) tres2=sqrt(tres2)

    if(rank==0.and.recordconvergence) THEN
      end_time =  MPI_WTIME()
      WRITE(89,*) ncall, it, end_time-start_time, tres2
    ENDIF

    if (tres2<maxError) exit

  ENDDO

  call update_bounds(nd) !restore default values (MG finished)

  if(it==maxit+1 .and. rank==0) write(*,*) 'Warning: LinearSolver reached maxit: ||res||: ',tres2
  ncall = ncall + 1

end subroutine NewSolverMG

subroutine get_MGtest_err(p)
  use module_grid
  use module_BC
  implicit none
  include 'mpif.h'
  integer :: i,j,k,ierr
  real(8) :: res, pi=3.14159265359, ref
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: p

  res = 0.d0; ref = 0.d0
  do k=ks,ke; do j=js,je; do i=is,ie
      res = res + &
            abs(p(i,j,k)-sin(2.*pi*x(i))*sin(2.*pi*y(j))*sin(2.*pi*z(k)))
      ref = ref + &
            abs(sin(2.*pi*x(i))*sin(2.*pi*y(j))*sin(2.*pi*z(k)))
  enddo; enddo; enddo
  call MPI_ALLREDUCE(res, res, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Cart, ierr) 
  if (rank==0) print *, 'error', res/ref

end subroutine get_MGtest_err

END MODULE module_mgsolver
