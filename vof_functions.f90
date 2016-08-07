!===============================================================================
!   Copyright (C) 2014  
!
!   Author: Ruben Scardovelli (ruben.scardovelli@unibo.it)
!
!   This file contains a few routines in FORTRAN 90 for standard
!   calculations with the Volume-of-Fluid (VOF) method
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <http://www.gnu.org/licenses/>.
!===============================================================================
!                              ***  NOTE  *** 
!
! The normal coefficients MUST satisfy the following relation 
! |nr1|+|nr2|+|nr3|=1
!
! HOWEVER, in function FL3D and in subroutine FLUX_CENTROID the normalization 
! is done internally as
! |nr1*dx1| + |nr2*dx2| + |nr3*dx3| = 1
! since a right hexahedron of sides (dx1,dx2,dx3) is considered
!
! FURTHERMORE, it would be better to avoid situations where 
! |nr_i| ~ 1, |nr_j| and/or |nr_k| < 1.d-14
! the user should take care of this issue by setting these tiny values 
! to zero before calling the routines 
!===============================================================================
! DESCRIPTION OF FUNCTION AL3D:
! compute in the unit cube the plane constant alpha satisfying
! [nr]*[x] = nr1*x1 + nr2*x2 + nr3*x3 = alpha
! where the normal is pointing outward from the reference phase
! INPUT: normal coefficients in nr(3) and volume fraction cc 
!-------------------------------------------------------------------------------
FUNCTION AL3D(nr,cc)

  IMPLICIT NONE
  REAL(8), INTENT(IN):: nr(3),cc
  REAL(8) :: AL3D
  REAL(8) :: cch,c01,c02,c03,np1,np2,np3
  REAL(8) :: m1,m2,m3,m12,numer,denom,p,pst,q,arc,csarc
  REAL(8), PARAMETER :: athird=1.d0/3.d0,eps0=1.d-50 
  INTRINSIC DABS,DMIN1,DMAX1,DSQRT,DACOS,DCOS

  np1 = DABS(nr(1))                                 ! need positive coefficients
  np2 = DABS(nr(2))
  np3 = DABS(nr(3))
  m1 = DMIN1(np1,np2)                              ! order positive coefficients
  m3 = DMAX1(np1,np2)
  if (np3 < m1) then
     m2 = m1
     m1 = np3
  else if (np3 >= m3) then
     m2 = m3
     m3 = np3
   else
     m2 = np3
  endif

  denom = DMAX1(6.d0*m1*m2*m3,eps0)                           
  cch = DMIN1(cc,1.d0-cc)                              ! limit to: 0 < cch < 1/2
  c01 = m1*m1*m1/denom                                          ! get cch ranges
  c02  = c01 + 0.5d0*(m2-m1)/m3
  m12 = m1 + m2
  if (m12 <= m3) then
     c03 = 0.5d0*m12/m3
  else
     numer = m3*m3*(3.d0*m12-m3) + m1*m1*(m1-3.d0*m3) + m2*m2*(m2-3.d0*m3)
     c03 = numer/denom
  endif

! 1: C<=C1; 2: C1<=C<=C2; 3: C2<=C<=C3; 4: C3<=C<=1/2 (a: m12<=m3; b: m3<m12)) 
  if (cch <= c01) then
     AL3D = (denom*cch)**athird                                       ! case (1)
  else if (cch <= c02) then 
     AL3D = 0.5d0*(m1 + DSQRT(m1*m1 + 8.d0*m2*m3*(cch-c01)))          ! case (2)
  else if (cch <= c03) then
     p = 2.d0*m1*m2
     q = 1.5d0*m1*m2*(m12 - 2.d0*m3*cch)
     pst = DSQRT(p)
     arc = athird*DACOS(q/(p*pst))
     csarc = DCOS(arc)
     AL3D = pst*(DSQRT(3.d0*(1.d0-csarc*csarc)) - csarc) + m12        ! case (3)
  else if (m12 <= m3) then
     AL3D = m3*cch + 0.5d0*m12                                       ! case (4a)
  else                                  
     p = m12*m3 + m1*m2 - 0.25d0                                     
     q = 1.5d0*m1*m2*m3*(0.5d0-cch)
     pst = DSQRT(p)
     arc = athird*DACOS(q/(p*pst))
     csarc = DCOS(arc)
     AL3D = pst*(DSQRT(3.d0*(1.d0-csarc*csarc)) - csarc) + 0.5d0     ! case (4b)
  endif

  if (cc > 0.5d0)  AL3D = 1.d0 - AL3D
  
! compute alpha for the given coefficients  
  AL3D = AL3D + DMIN1(0.d0,nr(1)) + DMIN1(0.d0,nr(2)) + DMIN1(0.d0,nr(3))
  
END FUNCTION AL3D

!===============================================================================
! DESCRIPTION OF FUNCTION CC3D:
! compute in the unit cube the volume cut by the plane
! [nr]*[x] = nr1*x1 + nr2*x2 + nr3*x3 = alpha
! where the normal is pointing outward from the reference phase
! INPUT: normal coefficients in nr(3) and plane constant alpha 
!-------------------------------------------------------------------------------
FUNCTION CC3D(nr,alpha)

  IMPLICIT NONE
  REAL(8), INTENT(IN):: nr(3),alpha
  REAL(8) :: CC3D
  REAL(8) :: al,alh,np1,np2,np3,m1,m2,m3,m12,mm,denom,top
  REAL(8), PARAMETER :: eps0=1.d-50 
  INTRINSIC DMAX1,DMIN1,DABS

  np1 = DABS(nr(1))                                 ! need positive coefficients
  np2 = DABS(nr(2))
  np3 = DABS(nr(3))
  m1 = DMIN1(np1,np2)                              ! order positive coefficients
  m3 = DMAX1(np1,np2)
  if (np3 < m1) then
     m2 = m1
     m1 = np3
  else if (np3 >= m3) then
     m2 = m3
     m3 = np3
   else
     m2 = np3
  endif

  al = alpha + DMAX1(0.d0,-nr(1))+DMAX1(0.d0,-nr(2))+DMAX1(0.d0,-nr(3))
  al = DMAX1(0.d0,DMIN1(1.d0,al))                                  ! safe limits
  alh = DMIN1(al,1.d0-al)                              ! limit to: 0 < alh < 1/2
  m12 = m1 + m2
  mm = DMIN1(m12,m3)
  denom = DMAX1(6.d0*m1*m2*m3,eps0)

! 1: al<=m1; 2: m1<=al<=m2; 3: m2<=al<=mm; 4: mm<=al<=1/2 (a:m12<=m3; b:m3<m12)) 
  if (alh <= m1) then
    CC3D = alh*alh*alh/denom                                          ! case (1)
  else if (alh <= m2) then
    CC3D = 0.5d0*alh*(alh-m1)/(m2*m3) +  m1*m1*m1/denom               ! case (2)
  else if (alh <= mm) then
    top = alh*alh*(3.d0*m12-alh) + m1*m1*(m1-3.d0*alh)              
    CC3D = (top + m2*m2*(m2-3.d0*alh))/denom                          ! case (3)
  else if (m12 <= m3) then
    CC3D = (alh - 0.5d0*m12)/m3                                      ! case (4a)
  else
    top = alh*alh*(3.d0-2.d0*alh) + m1*m1*(m1-3.d0*alh)             
    CC3D = (top + m2*m2*(m2-3.d0*alh) + m3*m3*(m3-3.d0*alh))/denom   ! case (4b)
  endif
  
  if (al > 0.5d0) CC3D = 1.d0 - CC3D

END FUNCTION CC3D

!===============================================================================
! DESCRIPTION OF FUNCTION FL3D:
! compute in the right hexahedron starting at (x01,x02,x03) and of sides 
! (dx1,dx2,dx3) the volume cut by the plane
! [nr]*[x] = nr1*x1 + nr2*x2 + nr3*x3 = alpha
! INPUT: normal coefficients in nr(3), plane constant alpha, starting
!        point x0(3), sides dx(3)
!-------------------------------------------------------------------------------
FUNCTION FL3D(nr,alpha,x0,dx)

  IMPLICIT NONE
  REAL(8), INTENT(IN):: nr(3),x0(3),dx(3),alpha
  REAL(8) :: FL3D
  REAL(8) :: al,almax,alh,np1,np2,np3,m1,m2,m3,m12,mm,denom,frac,top
  REAL(8), PARAMETER :: eps0=1.d-50 
  INTRINSIC DMAX1,DMIN1,DABS

! move origin to x0 
  al = alpha - nr(1)*x0(1) - nr(2)*x0(2) - nr(3)*x0(3)
! reflect the figure when negative coefficients
  al = al + DMAX1(0.d0,-nr(1)*dx(1)) + DMAX1(0.d0,-nr(2)*dx(2)) &
          + DMAX1(0.d0,-nr(3)*dx(3)) 
  np1 = DABS(nr(1))                                 ! need positive coefficients
  np2 = DABS(nr(2))
  np3 = DABS(nr(3))
  almax = np1*dx(1) + np2*dx(2) + np3*dx(3)                       
  al = DMAX1(0.d0,DMIN1(1.d0,al/almax))           !get new al within safe limits
  alh = DMIN1(al,1.d0-al)                              ! limit to: 0 < alh < 1/2


! normalized equation: m1*y1 + m2*y2 + m3*y3 = alh, with 0 <= m1 <= m2 <= m3
! the problem is then solved again in the unit cube
  np1 = np1/almax;
  np2 = np2/almax;
  np3 = np3/almax;
  m1 = DMIN1(np1*dx(1),np2*dx(2))                           ! order coefficients
  m3 = DMAX1(np1*dx(1),np2*dx(2))
  top = np3*dx(3)
  if (top < m1) then
     m2 = m1
     m1 = top
  else if (top >= m3) then
     m2 = m3
     m3 = top
   else
     m2 = top
  endif

  m12 = m1 + m2
  mm = DMIN1(m12,m3)
  denom = DMAX1(6.d0*m1*m2*m3,eps0)

! 1: al<=m1; 2: m1<=al<=m2; 3: m2<=al<=mm; 4: mm<=al<=1/2 (a:m12<=m3; b:m3<m12)) 
  if (alh <= m1) then
     frac = alh*alh*alh/denom                                         ! case (1)
  else if (alh <= m2) then
     frac = 0.5d0*alh*(alh-m1)/(m2*m3) +  m1*m1*m1/denom              ! case (2)
  else if (alh <= mm) then
     top = alh*alh*(3.d0*m12-alh) + m1*m1*(m1-3.d0*alh)              
     frac = (top + m2*m2*(m2-3.d0*alh))/denom                         ! case (3)
  else if (m12 <= m3) then
     frac = (alh - 0.5d0*m12)/m3                                     ! case (4a)
  else
     top = alh*alh*(3.d0-2.d0*alh) + m1*m1*(m1-3.d0*alh)             
     frac = (top + m2*m2*(m2-3.d0*alh) + m3*m3*(m3-3.d0*alh))/denom  ! case (4b)
  endif

  top = dx(1)*dx(2)*dx(3)
  if (al <= 0.5d0) then
     FL3D = frac*top
  else
     FL3D = (1.d0-frac)*top
  endif

END FUNCTION FL3D

!===============================================================================
! DESCRIPTION OF FUNCTION AREA3D:
! compute in the unit cube the area cut by the plane
! [nr]*[x] = nr1*x1 + nr2*x2 + nr3*x3 = alpha
! INPUT: normal coefficients in nr(3) and volume fraction cc 
! NOTE : the cut area A is invariant with respects to reflections, ordering  
!        of the coefficients and midpoint (i.e., A(cc) = A(1-cc))
!-------------------------------------------------------------------------------
FUNCTION AREA3D(nr,cc)

  IMPLICIT NONE
  REAL(8), INTENT(IN):: nr(3),cc
  REAL(8) :: AREA3D
  REAL(8) :: cch,ccr,al,c00,c01,c02,c03,np1,np2,np3
  REAL(8) :: m1,m2,m3,m12,numer,denom,p,pst,q,arc,csarc
  REAL(8), PARAMETER :: athird=1.d0/3.d0,eps0=1.d-50
  INTRINSIC DABS,DMIN1,DMAX1,DSQRT,DACOS,DCOS

  np1 = DABS(nr(1))                                 ! need positive coefficients
  np2 = DABS(nr(2))
  np3 = DABS(nr(3))
  m1 = DMIN1(np1,np2)                              ! order positive coefficients
  m3 = DMAX1(np1,np2)
  if (np3 < m1) then
     m2 = m1
     m1 = np3
  else if (np3 >= m3) then
     m2 = m3
     m3 = np3
   else
     m2 = np3
  endif

  denom = DMAX1(2.d0*m1*m2*m3,eps0)                           
  cch = DMIN1(cc,1.d0-cc)                              ! limit to: 0 < cch < 1/2
  ccr = DMAX1(cch,eps0)                                     ! avoid cch = m1 = 0  
  c01 = m1*m1*m1/(3.d0*denom)                                   ! get cch ranges
  c02  = c01 + 0.5d0*(m2-m1)/m3
  m12 = m1 + m2
  if (m12 <= m3) then
     c03 = 0.5d0*m12/m3
  else
     numer = m3*m3*(3.d0*m12-m3) + m1*m1*(m1-3.d0*m3) + m2*m2*(m2-3.d0*m3)
     c03 = numer/(3.d0*denom)
  endif
  c00 = DSQRT(m1*m1 + m2*m2 + m3*m3)

! 1: C<=C1; 2: C1<=C<=C2; 3: C2<=C<=C3; 4: C3<=C<=1/2 (a: m12<=m3; b: m3<m12)) 
  if (ccr <= c01) then
     al = (3.d0*denom*cch)**athird                                         
     AREA3D = c00*al*al/denom                                         ! case (1)
  else if (ccr <= c02) then 
    al = 0.5d0*(m1 + DSQRT(m1*m1 + 8.d0*m2*m3*(cch - c01)))           
    AREA3D = c00*(2.d0*al-m1)/(2.d0*m2*m3)                            ! case (2)
  else if (ccr <= c03) then
     p = 2.d0*m1*m2                                                   
     q = 1.5d0*m1*m2*(m12 - 2.d0*m3*cch)
     pst = DSQRT(p)
     arc = athird*DACOS(q/(p*pst))
     csarc = DCOS(arc)
     al = pst*(DSQRT(3.d0*(1.d0-csarc*csarc)) - csarc) + m12
     AREA3D = c00*(-al*al + m1*(2.d0*al-m1) + m2*(2.d0*al-m2))/denom  ! case (3)
  else if (m12 <= m3) then                                  
     al = m3*cch + 0.5d0*m12                                         
     AREA3D = c00/m3                                                 ! case (4a)
  else                                  
     p = m12*m3 + m1*m2 - 0.25d0                                     
     q = 1.5d0*m1*m2*m3*(0.5d0-cch)
     pst = DSQRT(p)
     arc = athird*DACOS(q/(p*pst))
     csarc = DCOS(arc)
     al = pst*(DSQRT(3.d0*(1.d0-csarc*csarc)) - csarc) + 0.5d0 
     AREA3D = c00*(2.d0*al*(1.d0-al) - c00*c00)/denom                ! case (4b)
  endif

END FUNCTION AREA3D

!===============================================================================
! DESCRIPTION OF SUBROUTINE CENT3D:
! NOTE: simply a wrapper to the new subroutine AREA_CENTROID
!-------------------------------------------------------------------------------
SUBROUTINE CENT3D(nr,cc,xc0)

  IMPLICIT NONE
  REAL(8), INTENT(IN) :: nr(3),cc
  REAL(8), INTENT(OUT) :: xc0(3)

  call AREA_CENTROID(nr,cc,xc0)

END SUBROUTINE CENT3D

!===============================================================================
! DESCRIPTION OF SUBROUTINE AREA_CENTROID:
! compute in the unit cube the coordinates of the centroid of the area 
! cut by the plane
! [nr]*[x] = nr1*x1 + nr2*x2 + nr3*x3 = alpha
! INPUT : normal coefficients in nr(3) and volume fraction cc 
! OUTPUT: centroid coordinates in xc0(3)
! NOTE  : the centroid coordinates change with respect to reflections, ordering 
!         of the coefficients and central point (i.e., xc0(cc) = 1 - xc0(1-cc))
!-------------------------------------------------------------------------------
SUBROUTINE AREA_CENTROID(nr,cc,xc0)

  IMPLICIT NONE
  REAL(8), INTENT(IN) :: nr(3),cc
  REAL(8), INTENT(OUT) :: xc0(3)
  REAL(8) :: ctd0(3)
  REAL(8) :: cch,ccr,al,c01,c02,c03,np1,np2,np3
  REAL(8) :: m1,m2,m3,m12,numer,denom,p,pst,q,arc,csarc,top,bot
  INTEGER :: ind(3)
  REAL(8), PARAMETER :: athird=1.d0/3.d0,eps0=1.d-50
  INTRINSIC DMAX1,DMIN1,DSQRT,DACOS,DCOS

  np1 = DABS(nr(1))                                 ! need positive coefficients
  np2 = DABS(nr(2))
  np3 = DABS(nr(3))
! coefficients in ascending order
  if (np1 <= np2) then                             
    m1 = np1
    m3 = np2
    ind(1) = 1                                            
    ind(3) = 2
  else
    m1 = np2
    m3 = np1
    ind(1) = 2
    ind(3) = 1
  endif
  
  if (np3 < m1) then
    m2 = m1
    m1 = np3
    ind(2) = ind(1)
    ind(1) = 3
  else if (np3 >= m3) then
    m2 = m3
    m3 = np3
    ind(2) = ind(3)
    ind(3) = 3
  else
    m2 = np3
    ind(2) = 3
  endif

  denom = DMAX1(6.d0*m1*m2*m3,eps0)                           
  cch = DMIN1(cc,1.d0-cc)                              ! limit to: 0 < cch < 1/2
  ccr = DMAX1(cch,eps0)                                     ! avoid cch = m1 = 0  
  c01 = m1*m1*m1/denom                                          ! get cch ranges
  c02  = c01 + 0.5d0*(m2-m1)/m3
  m12 = m1 + m2
  if (m12 <= m3) then
    c03 = 0.5d0*m12/m3
  else
    numer = m3*m3*(3.d0*m12-m3) + m1*m1*(m1-3.d0*m3) + m2*m2*(m2-3.d0*m3)
    c03 = numer/denom
  endif
  
! 1: C<=C1; 2: C1<=C<=C2; 3: C2<=C<=C3; 4: C3<=C<=1/2 (a: m12<=m3; b: m3<m12)) 
  if (ccr <= c01) then
    al = (denom*cch)**athird 
    ctd0(1) = athird*al/m1
    ctd0(2) = athird*al/m2
    ctd0(3) = athird*al/m3                                            ! case (1)
  else if (ccr <= c02) then 
    al = 0.5d0*(m1 + DSQRT(m1*m1 + 8.d0*m2*m3*(cch - c01)))    
    top = m1*m1 + 3.*al*(al-m1)
    bot = DMAX1(3.*(2.*al-m1),eps0)
    ctd0(1) = (3.*al-2.*m1)/bot
    ctd0(2) = top/(m2*bot)
    ctd0(3) = top/(m3*bot)                                            ! case (2)
  else if (ccr <= c03) then
    p = 2.d0*m1*m2                                                   
    q = 1.5d0*m1*m2*(m12 - 2.d0*m3*cch)
    pst = DSQRT(p)
    arc = athird*DACOS(q/(p*pst))
    csarc = DCOS(arc)
    al = pst*(DSQRT(3.d0*(1.d0-csarc*csarc)) - csarc) + m12
    bot = 3.*((al-m2)*(al-m2) + (al-m1)*(al-m1) - al*al)
    top = (al-m2)*(al-m2)*(al-m2) + m1*m1*(2.*m1-3.*al)
    ctd0(1) = top/(m1*bot)
    top = (al-m1)*(al-m1)*(al-m1) + m2*m2*(2.*m2-3.*al)
    ctd0(2) = top/(m2*bot)
!   top = (al-m1)*(al-m1)*(al-m1) + (al-m2)*(al-m2)*(al-m2) - al*al*al
!   ctd0(3) = top/(m3*bot)                                           
    ctd0(3) = (al - m1*ctd0(1) - m2*ctd0(2))/m3                       ! case (3)
  else if (m12 <= m3) then                                  
    al = m3*cch + 0.5d0*m12                                         
    ctd0(1) = 0.5d0
    ctd0(2) = 0.5d0
    ctd0(3) = (al-0.5d0*m12)/m3                                      ! case (4a) 
  else                                  
    p = m12*m3 + m1*m2 - 0.25d0                                     
    q = 1.5d0*m1*m2*m3*(0.5d0-cch)
    pst = DSQRT(p)
    arc = athird*DACOS(q/(p*pst))
    csarc = DCOS(arc)
    al = pst*(DSQRT(3.d0*(1.d0-csarc*csarc)) - csarc) + 0.5d0 
    bot = 2.d0*al*(1.d0-al) - (m1*m1 + m2*m2 + m3*m3)
    top = (m1*m1*m1 + m2*m2*m2 + m3*m3*m3) -3.*al*(m1*m1 + m2*m2 + m3*m3) 
    top = athird*(top + 3.*al*al - 2.*al*al*al)
    ctd0(1) = (top/m1 - (al-m1)*(al-m1))/bot
    ctd0(2) = (top/m2 - (al-m2)*(al-m2))/bot
!   ctd0(3) = (top/m3 - (al-m3)*(al-m3))/bot                         
    ctd0(3) = (al - m1*ctd0(1) - m2*ctd0(2))/m3                      ! case (4b)
  endif

! back to actual cc value
  if (cc > 0.5d0) then
    ctd0(1) = 1.d0 - ctd0(1)
    ctd0(2) = 1.d0 - ctd0(2)
    ctd0(3) = 1.d0 - ctd0(3)
  endif

! get correct indexing
  xc0(ind(1)) = ctd0(1)                              
  xc0(ind(2)) = ctd0(2)
  xc0(ind(3)) = ctd0(3)

! take care of negative coefficients
  if (nr(1) < 0.d0) xc0(1) = 1.d0 - xc0(1) 
  if (nr(2) < 0.d0) xc0(2) = 1.d0 - xc0(2)
  if (nr(3) < 0.d0) xc0(3) = 1.d0 - xc0(3)

END SUBROUTINE AREA_CENTROID

!===============================================================================
! DESCRIPTION OF SUBROUTINE VOLUME_CENTROID:
! compute in the unit cube the coordinates of the centroid of the volume 
! cut by the plane
! [nr]*[x] = nr1*x1 + nr2*x2 + nr3*x3 = alpha
! INPUT : normal coefficients in nr(3) and volume fraction cc
! OUTPUT: centroid coordinates in xc0(3)
! NOTE  : the centroid coordinates change with respect to reflections, ordering 
!         of the coefficients and central point (i.e., xc0(cc) = 1 - xc0(1-cc))
!-------------------------------------------------------------------------------
SUBROUTINE VOLUME_CENTROID(nr,cc,xc0)

  IMPLICIT NONE
  REAL(8), INTENT(IN) :: nr(3),cc
  REAL(8), INTENT(OUT) :: xc0(3)
  REAL(8) :: ctd0(3)
  REAL(8) :: cch,ccr,alh,c01,c02,c03,np1,np2,np3
  REAL(8) :: m1,m2,m3,m12,numer,denom,p,pst,q,arc,csarc
  REAL(8) :: tmp1,tmp2,tmp3,tmp4,a2,bot1
  INTEGER :: ind(3)
  REAL(8), PARAMETER :: athird=1.d0/3.d0,eps0=1.d-50 
  INTRINSIC DMAX1,DMIN1,DSQRT,DACOS,DCOS

  np1 = DABS(nr(1))                                 ! need positive coefficients
  np2 = DABS(nr(2))
  np3 = DABS(nr(3))
! coefficients in ascending order
  if (np1 <= np2) then                             
    m1 = np1
    m3 = np2
    ind(1) = 1                                            
    ind(3) = 2
  else
    m1 = np2
    m3 = np1
    ind(1) = 2
    ind(3) = 1
  endif

  if (np3 < m1) then
    m2 = m1
    m1 = np3
    ind(2) = ind(1)
    ind(1) = 3
  else if (np3 >= m3) then
    m2 = m3
    m3 = np3
    ind(2) = ind(3)
    ind(3) = 3
  else
    m2 = np3
    ind(2) = 3
  endif

  denom = DMAX1(6.d0*m1*m2*m3,eps0)                           
  cch = DMIN1(cc,1.d0-cc)                              ! limit to: 0 < cch < 1/2
  ccr = DMAX1(cch,eps0)                                     ! avoid cch = m1 = 0
  c01 = m1*m1*m1/denom                                          ! get cch ranges 
  c02  = c01 + 0.5d0*(m2-m1)/m3
  m12 = m1 + m2
  if (m12 <= m3) then
    c03 = 0.5d0*m12/m3
  else
    numer = m3*m3*(3.d0*m12-m3) + m1*m1*(m1-3.d0*m3) + m2*m2*(m2-3.d0*m3)
    c03 = numer/denom
  endif

! 1: C<=C1; 2: C1<=C<=C2; 3: C2<=C<=C3; 4: C3<=C<=1/2 (a: m12<=m3; b: m3<m12)) 
  if (ccr <= c01) then
    alh = 0.25d0*(denom*cch)**athird 
    ctd0(1) = alh/m1 
    ctd0(2) = alh/m2
    ctd0(3) = alh/m3                                                ! case (1)
  else if (ccr <= c02) then 
    alh = 0.5d0*(m1 + DSQRT(m1*m1 + 8.d0*m2*m3*(cch - c01)))    
    tmp1 = (2.d0*alh-m1)
    tmp2 = (2.d0*alh*alh - 2.d0*alh*m1 + m1*m1)
    bot1 = DMAX1(4.d0*(3.d0*alh*alh - 3.d0*alh*m1 + m1*m1),eps0)
    tmp2 = tmp1*tmp2/bot1
    ctd0(1) = 0.5d0 - m1*tmp1/bot1 
    ctd0(2) = tmp2/m2
    ctd0(3) = tmp2/m3                                                ! case (2)
  else if (ccr <= c03) then
    p = 2.d0*m1*m2                                                   
    q = 1.5d0*m1*m2*(m12 - 2.d0*m3*cch)
    pst = DSQRT(p)
    arc = athird*DACOS(q/(p*pst))
    csarc = DCOS(arc)
    alh = pst*(DSQRT(3.d0*(1.d0-csarc*csarc)) - csarc) + m12
    a2 = alh*alh
    tmp1 = (alh-m1)*(alh-m1)*(alh-m1)
    tmp2 = (alh-m2)*(alh-m2)*(alh-m2)
    bot1 = 4.d0*(a2*alh - tmp1 - tmp2)
    ctd0(1) = (a2*a2 - tmp1*(alh+3.d0*m1) - tmp2*(alh-m2))/(m1*bot1) 
    ctd0(2) = (a2*a2 - tmp1*(alh-m1) - tmp2*(alh+3.d0*m2))/(m2*bot1) 
    ctd0(3) = (a2*a2 - tmp1*(alh-m1) - tmp2*(alh-m2))/(m3*bot1)      ! case (3)
  else if (m12 <= m3) then
    alh = m3*cch + 0.5d0*m12                                         
    bot1 = DMAX1((2.d0*alh - m12),eps0)
    ctd0(1) = 0.5d0 - m1/(6.d0*bot1)
    ctd0(2) = 0.5d0 - m2/(6.d0*bot1)                                ! case (4a) 
    ctd0(3) = ((3.d0*alh - 2.d0*m12)*bot1 + alh*m12 - m1*m2)/(6.d0*m3*bot1)
  else
    p = m12*m3 + m1*m2 - 0.25d0                                     
    q = 1.5d0*m1*m2*m3*(0.5d0-cch)
    pst = DSQRT(p)
    arc = athird*DACOS(q/(p*pst))
    csarc = DCOS(arc)
    alh = pst*(DSQRT(3.d0*(1.d0-csarc*csarc)) - csarc) + 0.5d0 
    tmp1 = m1 + m2 + m3
    tmp2 = m1*m1 + m2*m2 + m3*m3
    tmp3 = m1*m1*m1 + m2*m2*m2 + m3*m3*m3
    tmp4 = m1*m1*m1*m1 + m2*m2*m2*m2 + m3*m3*m3*m3
    a2 = alh*alh
    bot1 = 4.d0*(2.d0*a2*alh - 3.d0*a2*tmp1 + 3.d0*alh*tmp2 - tmp3)
    tmp1 = 2.d0*a2*a2 - 4.d0*a2*alh*tmp1 + 6.d0*a2*tmp2 - 4.d0*alh*tmp3 + tmp4
    ctd0(1) = (tmp1 + 4.d0*m1*(alh-m1)*(alh-m1)*(alh-m1))/(m1*bot1)
    ctd0(2) = (tmp1 + 4.d0*m2*(alh-m2)*(alh-m2)*(alh-m2))/(m2*bot1)
    ctd0(3) = (tmp1 + 4.d0*m3*(alh-m3)*(alh-m3)*(alh-m3))/(m3*bot1) ! case (4b)
  endif

! back to actual c value, but the centroid of a full cell is the cell center
  if (cc > 0.5d0) then
    ctd0(1) = (0.5d0 - (1.d0 - cc)*(1.d0-ctd0(1)))/cc
    ctd0(2) = (0.5d0 - (1.d0 - cc)*(1.d0-ctd0(2)))/cc
    ctd0(3) = (0.5d0 - (1.d0 - cc)*(1.d0-ctd0(3)))/cc
  endif

! get correct indexing
  xc0(ind(1)) = ctd0(1)                              
  xc0(ind(2)) = ctd0(2)
  xc0(ind(3)) = ctd0(3)

! take care of negative coefficients
  if (nr(1) < 0.d0) xc0(1) = 1.d0 - xc0(1) 
  if (nr(2) < 0.d0) xc0(2) = 1.d0 - xc0(2)
  if (nr(3) < 0.d0) xc0(3) = 1.d0 - xc0(3)

END SUBROUTINE VOLUME_CENTROID

!===============================================================================
! DESCRIPTION OF SUBROUTINE FLUX_CENTROID:
! compute in the right hexahedron starting at (x01,x02,x03) and of sides 
! (dx1,dx2,dx3) the centroid of the volume cut by the plane
! [nr]*[x] = nr1*x1 + nr2*x2 + nr3*x3 = alpha
! INPUT : normal coefficients in nr(3), plane constant alpha, starting
!        point x0(3), sides dx(3)
! OUTPUT: centroid coordinates in xc0(3)
! NOTE  : the centroid coordinates change with respect to reflections, ordering 
!         of the coefficients and central point (i.e., xc0(cc) = 1 - xc0(1-cc))
!-------------------------------------------------------------------------------
SUBROUTINE FLUX_CENTROID(nr,alpha,x0,dx,xc0)

  IMPLICIT NONE
  REAL(8), INTENT(IN) :: nr(3),x0(3),dx(3),alpha
  REAL(8), INTENT(OUT) :: xc0(3)
  REAL(8) :: ctd0(3)
  REAL(8) :: al,almax,alh,alr,np1,np2,np3,m1,m2,m3,m12,mm,denom
  REAL(8) :: tmp1,tmp2,tmp3,tmp4,a2,bot1,frac
  REAL(8), PARAMETER :: eps0=1.d-50
  INTEGER :: ind(3)
  INTRINSIC DMAX1,DMIN1,DABS

! move origin to x0 
  al = alpha - nr(1)*x0(1) - nr(2)*x0(2) - nr(3)*x0(3)
! reflect the figure when negative coefficients
  al = al + DMAX1(0.d0,-nr(1)*dx(1)) + DMAX1(0.d0,-nr(2)*dx(2)) &
          + DMAX1(0.d0,-nr(3)*dx(3)) 
  np1 = DABS(nr(1))                                 ! need positive coefficients
  np2 = DABS(nr(2))
  np3 = DABS(nr(3))
  almax = np1*dx(1) + np2*dx(2) + np3*dx(3)                       
  al = DMAX1(0.d0,DMIN1(1.d0,al/almax))           !get new al within safe limits
  alh = DMIN1(al,1.d0-al)                              ! limit to: 0 < alh < 1/2
  alr = DMAX1(alh,eps0)                                     ! avoid alh = m1 = 0

! normalized equation: m1*y1 + m2*y2 + m3*y3 = alh, with 0 <= m1 <= m2 <= m3
! the problem is then solved again in the unit cube
  np1 = np1*dx(1)/almax;
  np2 = np2*dx(2)/almax;
  np3 = np3*dx(3)/almax;
! coefficients in ascending order
  if (np1 <= np2) then                             
    m1 = np1
    m3 = np2
    ind(1) = 1                                            
    ind(3) = 2
  else
    m1 = np2
    m3 = np1
    ind(1) = 2
    ind(3) = 1
  endif

  if (np3 < m1) then
    m2 = m1
    m1 = np3
    ind(2) = ind(1)
    ind(1) = 3
  else if (np3 >= m3) then
    m2 = m3
    m3 = np3
    ind(2) = ind(3)
    ind(3) = 3
  else
    m2 = np3
    ind(2) = 3
  endif

  m12 = m1 + m2
  mm = DMIN1(m12,m3)
  denom = DMAX1(6.d0*m1*m2*m3,eps0)

! 1: al<=m1; 2: m1<=al<=m2; 3: m2<=al<=mm; 4: mm<=al<=1/2 (a:m12<=m3; b:m3<m12)) 
  if (alr <= m1) then
    tmp1 = 0.25d0*alh
    ctd0(1) = tmp1/m1 
    ctd0(2) = tmp1/m2
    ctd0(3) = tmp1/m3                                                
    frac = alh*alh*alh/denom                                          ! case (1)
  else if (alr <= m2) then
    tmp1 = (2.d0*alh-m1)
    tmp2 = (2.d0*alh*alh - 2.d0*alh*m1 + m1*m1)
    bot1 = 4.d0*(3.d0*alr*alr - 3.d0*alh*m1 + m1*m1)
    tmp2 = tmp1*tmp2/bot1
    ctd0(1) = 0.5d0 - m1*tmp1/bot1 
    ctd0(2) = tmp2/m2
    ctd0(3) = tmp2/m3                                                
    frac = 0.5d0*alh*(alh-m1)/(m2*m3) +  m1*m1*m1/denom               ! case (2)
  else if (alr <= mm) then
    a2 = alh*alh
    tmp1 = (alh-m1)*(alh-m1)*(alh-m1)
    tmp2 = (alh-m2)*(alh-m2)*(alh-m2)
    bot1 = 4.d0*(a2*alh - tmp1 - tmp2)
    ctd0(1) = (a2*a2 - tmp1*(alh+3.d0*m1) - tmp2*(alh-m2))/(m1*bot1) 
    ctd0(2) = (a2*a2 - tmp1*(alh-m1) - tmp2*(alh+3.d0*m2))/(m2*bot1) 
    ctd0(3) = (a2*a2 - tmp1*(alh-m1) - tmp2*(alh-m2))/(m3*bot1)      
    tmp3 = alh*alh*(3.d0*m12-alh) + m1*m1*(m1-3.d0*alh)              
    frac = (tmp3 + m2*m2*(m2-3.d0*alh))/denom                         ! case (3)
  else if (m12 <= m3) then
    bot1 = DMAX1((2.d0*alh - m12),eps0)
    ctd0(1) = 0.5d0 - m1/(6.d0*bot1)
    ctd0(2) = 0.5d0 - m2/(6.d0*bot1)                                
    ctd0(3) = ((3.d0*alh - 2.d0*m12)*bot1 + alh*m12 - m1*m2)/(6.d0*m3*bot1)
    frac = (alh - 0.5d0*m12)/m3                                      ! case (4a)
  else
    tmp1 = m1 + m2 + m3
    tmp2 = m1*m1 + m2*m2 + m3*m3
    tmp3 = m1*m1*m1 + m2*m2*m2 + m3*m3*m3
    tmp4 = m1*m1*m1*m1 + m2*m2*m2*m2 + m3*m3*m3*m3
    a2 = alh*alh
    bot1 = 4.d0*(2.d0*a2*alh - 3.d0*a2*tmp1 + 3.d0*alh*tmp2 - tmp3)
    tmp1 = 2.d0*a2*a2 - 4.d0*a2*alh*tmp1 + 6.d0*a2*tmp2 - 4.d0*alh*tmp3 + tmp4
    ctd0(1) = (tmp1 + 4.d0*m1*(alh-m1)*(alh-m1)*(alh-m1))/(m1*bot1)
    ctd0(2) = (tmp1 + 4.d0*m2*(alh-m2)*(alh-m2)*(alh-m2))/(m2*bot1)
    ctd0(3) = (tmp1 + 4.d0*m3*(alh-m3)*(alh-m3)*(alh-m3))/(m3*bot1) 
    tmp1 = alh*alh*(3.d0-2.d0*alh) + m1*m1*(m1-3.d0*alh)             
    frac = (tmp1 + m2*m2*(m2-3.d0*alh) + m3*m3*(m3-3.d0*alh))/denom  ! case (4b)
  endif

! back to actual al value, but the centroid of a full cell is the cell center
! must be careful that frac = cc when al < 0.5, otherwise frac = 1 - cc
  if (al > 0.5d0) then
    ctd0(1) = (0.5d0 - frac*(1.d0-ctd0(1)))/(1.d0-frac)
    ctd0(2) = (0.5d0 - frac*(1.d0-ctd0(2)))/(1.d0-frac)
    ctd0(3) = (0.5d0 - frac*(1.d0-ctd0(3)))/(1.d0-frac)
  endif

! get correct indexing
  xc0(ind(1)) = ctd0(1)                              
  xc0(ind(2)) = ctd0(2)
  xc0(ind(3)) = ctd0(3)

! take care of negative coefficients
  if (nr(1) < 0.d0) xc0(1) = 1.d0 - xc0(1) 
  if (nr(2) < 0.d0) xc0(2) = 1.d0 - xc0(2)
  if (nr(3) < 0.d0) xc0(3) = 1.d0 - xc0(3)

! final position of the centroid with respect to the cell origin
! and by considering the actual sides of the hexahedron 
  xc0(1) = x0(1) + xc0(1)*dx(1) 
  xc0(2) = x0(2) + xc0(2)*dx(2) 
  xc0(3) = x0(3) + xc0(3)*dx(3) 

END SUBROUTINE FLUX_CENTROID

!===============================================================================
! to make it a module add the following lines at the top

!!$MODULE VOF_FUNCTIONS
!!$
!!$  IMPLICIT NONE
!!$
!!$CONTAINS

! and the following one at the bottom

!!$END MODULE VOF_FUNCTIONS
!===============================================================================
