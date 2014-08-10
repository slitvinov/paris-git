!================================================================================================
!=================================================================================================
! Paris-0.1
!
! Free surface extensions
! written by Leon Malan and Stephane Zaleski
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
!-------------------------------------------------------------------------------------------------
subroutine setuppoisson_fs(utmp,vtmp,wtmp,umask,vmask,wmask,rhot,dt,A,pmask,cvof,n1,n2,n3)
  use module_grid
  use module_BC
  use module_2phase
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: utmp,vtmp,wtmp,rhot
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: umask,vmask,wmask
  real(8), dimension(is:ie,js:je,ks:ke), intent(inout) :: pmask
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: cvof,n1,n2,n3
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(out) :: A
  real(8), dimension(is:ie,js:je,ks:ke) :: x_int, y_int, z_int
  real(8) :: alpha, x_test, y_test, z_test, n_x, n_y, n_z
  real(8) :: x_l, x_r, y_b, y_t, z_r, z_f
  real(8) :: nr(3),al3dnew
  real(8) :: dt, limit
  integer :: i,j,k,l

  x_int = 0d0; y_int = 0d0; z_int = 0d0
  pmask = 1d0

  do k=ks,ke; do j=js,je; do i=is,ie
     limit = 0.10*min(dx(i),dy(j),dz(k))
     if(cvof(i,j,k) >= 0.5d0) then ! pressure 0 in the cvof=1 phase. 
        pmask(i,j,k) = 0d0

        nr(1) = n1(i,j,k);         nr(2) = n2(i,j,k);         nr(3) = n3(i,j,k)
        alpha = al3dnew(nr,cvof(i,j,k))

        if (n_x .ne. 0.d0) then
           x_test = (alpha - (n_y+n_z)/2d0)/n_x
           if (x_test<1.d0 .and. x_test>0d0) then
              if (n1(i,j,k) > 0d0) then 
                 x_int(i,j,k) = xh(i-1) + x_test*dx(i)
                 x_l = x(i+1) - x_int(i,j,k)
                 x_mod(i,j,k) = dxh(i) - x_l
                 A(i+1,j,k,1) = 2d0*dt*umask(i,j,k)/(dx(i+1)*x_l*(rhot(i,j,k)+rhot(i+1,j,k)))
              else 
                 x_int(i,j,k) = xh(i) - x_test*dx(i)
                 x_r = x_int(i,j,k) - x(i-1) 
                 x_mod(i-1,j,k) = dxh(i-1) - x_r
                 A(i-1,j,k,2) = 2d0*dt*umask(i-1,j,k)/(dx(i-1)*x_r*(rhot(i,j,k)+rhot(i-1,j,k)))
              endif
           endif
        endif

        if (n_y .ne. 0.d0) then
           y_test = (alpha - (n_x+n_z)/2d0)/n_y
           if (y_test<1.d0 .and. y_test>0d0) then
              if (n2(i,j,k) > 0d0) then 
                 y_int(i,j,k) = yh(j-1) + y_test*dy(j)
                 y_b = y(j+1) - y_int(i,j,k)
                 y_mod(i,j,k) = dyh(j)-y_b
                 A(i,j+1,k,3) = 2d0*dt*vmask(i,j,k)/(dy(j+1)*y_b*(rhot(i,j,k)+rhot(i,j+1,k)))
              else 
                 y_int(i,j,k) = yh(j) - y_test*dy(j)
                 y_t = y_int(i,j,k) - y(j-1)
                 y_mod(i,j-1,k) = dyh(j-1)-y_t
                 A(i,j-1,k,4) = 2d0*dt*vmask(i,j-1,k)/(dy(j-1)*y_t*(rhot(i,j,k)+rhot(i,j-1,k)))
              endif
           endif
        endif

        if (n_z .ne. 0.d0) then
           z_test = (alpha - (n_x+n_y)/2d0)/n_z
           if (z_test<1.d0 .and. z_test>0d0) then
              if (n3(i,j,k) > 0d0) then 
                 z_int(i,j,k) = zh(k-1) + z_test*dz(k)
                 z_r = z(k+1) - z_int(i,j,k)
                 z_mod(i,j,k) = dzh(k)-z_r
                 A(i,j,k+1,5) = 2d0*dt*wmask(i,j,k)/(dz(k+1)*z_r*(rhot(i,j,k)+rhot(i,j,k+1)))
              else 
                 z_int(i,j,k) = zh(k) - z_test*dz(k)
                 z_f = z_int(i,j,k) - z(k-1)
                 z_mod(i,j,k-1) = dzh(k-1)-z_f
                 A(i,j,k-1,6) = 2d0*dt*wmask(i,j,k-1)/(dz(k-1)*z_f*(rhot(i,j,k)+rhot(i,j,k-1)))
              endif
           endif
        endif
     endif

     if (cvof(i,j,k)>0d0 .and. cvof(i,j,k)<0.5d0) then

        nr(1) = n1(i,j,k);         nr(2) = n2(i,j,k);         nr(3) = n3(i,j,k)
        alpha = al3dnew(nr,cvof(i,j,k))

        if (n_x .ne. 0.d0) then
           x_test = (alpha - (n_y+n_z)/2d0)/n_x
           if (x_test<1.d0 .and. x_test>0d0) then
              if (n1(i,j,k) > 0d0) then 
                 x_int(i,j,k) = xh(i-1) + x_test*dx(i)
                 x_l = x(i) - x_int(i,j,k)
                 if (x_l < limit) x_l = limit !arbitrary small limit, to be evaluated
                 A(i,j,k,1) = 2d0*dt*umask(i-1,j,k)/((dx(i)+x_l)/2d0*x_l*(rhot(i-1,j,k)+rhot(i,j,k)))
                 A(i,j,k,2) = 2d0*dt*umask(i,j,k)/((dx(i)+x_l)/2d0*dxh(i+1)*(rhot(i,j,k)+rhot(i+1,j,k)))
              else 
                 x_int(i,j,k) = xh(i) - x_test*dx(i)
                 x_r = x_int(i,j,k) - x(i)
                 if (x_r < limit) x_r = limit !arbitrary small limit, to be evaluated
                 A(i,j,k,1) = 2d0*dt*umask(i-1,j,k)/((dx(i)+x_r)/2d0*dxh(i-1)*(rhot(i,j,k)+rhot(i-1,j,k)))
                 A(i,j,k,2) = 2d0*dt*umask(i,j,k)/((dx(i)+x_r)/2d0*x_r*(rhot(i+1,j,k)+rhot(i,j,k)))
              endif
           endif
        endif

        if (n_y .ne. 0.d0) then
           y_test = (alpha - (n_x+n_z)/2d0)/n_y
           if (y_test<1.d0 .and. y_test>0d0) then
              if (n2(i,j,k) > 0d0) then 
                 y_int(i,j,k) = yh(j-1) + y_test*dy(j)
                 y_b = y(j) - y_int(i,j,k)
                 if (y_b < limit) y_b = limit !arbitrary small limit, to be evaluated
                 A(i,j,k,3) = 2d0*dt*vmask(i,j-1,k)/((dy(j)+y_b)/2d0*y_b*(rhot(i,j-1,k)+rhot(i,j,k)))
                 A(i,j,k,4) = 2d0*dt*vmask(i,j,k)/((dy(j)+y_b)/2d0*dyh(j+1)*(rhot(i,j,k)+rhot(i,j+1,k)))
              else 
                 y_int(i,j,k) = yh(j) - y_test*dy(j)
                 y_t = y_int(i,j,k) - y(j) 
                 if (y_t < limit) y_t = limit !arbitrary small limit, to be evaluated
                 A(i,j,k,3) = 2d0*dt*vmask(i,j-1,k)/((dy(j)+y_t)/2d0*dyh(j-1)*(rhot(i,j,k)+rhot(i,j-1,k)))
                 A(i,j,k,4) = 2d0*dt*vmask(i,j,k)/((dy(j)+y_t)/2d0*y_t*(rhot(i,j+1,k)+rhot(i,j,k)))
              endif
           endif
        endif

        if (n_z .ne. 0.d0) then
           z_test = (alpha - (n_x+n_y)/2d0)/n_z
           if (z_test<1.d0 .and. z_test>0d0) then
              if (n3(i,j,k) > 0d0) then 
                 z_int(i,j,k) = zh(k-1) + z_test*dz(k)
                 z_r = z(k) - z_int(i,j,k)
                 if (z_r < limit) z_r = limit !arbitrary small limit, to be evaluated
                 A(i,j,k,5) = 2d0*dt*wmask(i,j,k-1)/((dz(k)+z_r)/2d0*z_r*(rhot(i,j,k-1)+rhot(i,j,k)))
                 A(i,j,k,6) = 2d0*dt*wmask(i,j,k)/((dz(k)+z_r)/2d0*dzh(k+1)*(rhot(i,j,k)+rhot(i,j,k+1)))
              else 
                 z_int(i,j,k) = zh(k) - z_test*dz(k)
                 z_f = z_int(i,j,k) - z(k)
                 if (z_f < limit) z_f = limit !arbitrary small limit, to be evaluated
                 A(i,j,k,5) = 2d0*dt*wmask(i,j,k-1)/((dz(k)+z_f)/2d0*dzh(k-1)*(rhot(i,j,k)+rhot(i,j,k-1)))
                 A(i,j,k,6) = 2d0*dt*wmask(i,j,k)/((dz(k)+z_f)/2d0*z_f*(rhot(i,j,k+1)+rhot(i,j,k)))
              endif
           endif
        endif
     endif
  enddo;enddo;enddo

  do k=ks,ke; do j=js,je; do i=is,ie
     do l=1,6
        A(i,j,k,l) = pmask(i,j,k)*A(i,j,k,l)
     enddo
     A(i,j,k,7) = sum(A(i,j,k,1:6)) + (1d0-pmask(i,j,k)) + 1d-49
     A(i,j,k,8) = pmask(i,j,k)*A(i,j,k,8)
  enddo;enddo;enddo
end subroutine setuppoisson_fs
