!=================================================================================================
!=================================================================================================
! PARIS  Parallel Robust Interface Simulator 
! NEW VOF FUNCTIONS
!
! author:  Ruben Scardovelli ruben.scardovelli@unibo.it
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
!===============================================================================
! DESCRIPTION OF FUNCTION AL3DNEW:
! compute in the unit cube the coefficient alpha that satisfies
! [nr]*[x] = nr1*x1 + nr2*x2 + nr3*x3 = alpha
! INPUT: normal coefficients (nr1,nr2,nr3) and volume fraction cc 
! NOTE : the normal coefficients MUST satisfy the relation |nr1|+|nr2|+|nr3|=1
!-------------------------------------------------------------------------------
FUNCTION AL3DNEW(nr,cc)

  IMPLICIT NONE
  REAL(8), INTENT(IN):: nr(3),cc
  REAL(8) :: AL3DNEW
  REAL(8) :: cch,c01,c02,c03,np1,np2,np3
  REAL(8) :: m1,m2,m3,m12,numer,denom,p,pst,q,arc,csarc
  REAL(8), PARAMETER :: athird=1.d0/3.d0 
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
  cch = DMIN1(cc,1.d0-cc)                              ! symmetry: 0 < cch < 1/2
  denom = DMAX1(6.d0*m1*m2*m3,1.d-50)                           ! get cch ranges
  c01 = m1*m1*m1/denom
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
     AL3DNEW = (denom*cch)**athird                                       ! case (1)
  else if (cch <= c02) then 
     AL3DNEW = 0.5d0*(m1 + DSQRT(m1*m1 + 8.d0*m2*m3*(cch-c01)))          ! case (2)
  else if (cch <= c03) then
     p = 2.d0*m1*m2
     q = 1.5d0*m1*m2*(m12 - 2.d0*m3*cch)
     pst = DSQRT(p)
     arc = athird*DACOS(q/(p*pst))
     csarc = DCOS(arc)
     AL3DNEW = pst*(DSQRT(3.d0*(1.d0-csarc*csarc)) - csarc) + m12        ! case (3)
  else if (m12 <= m3) then
     AL3DNEW = m3*cch + 0.5d0*m12                                       ! case (4a)
  else                                  
     p = m12*m3 + m1*m2 - 0.25d0                                     
     q = 1.5d0*m1*m2*m3*(0.5d0-cch)
     pst = DSQRT(p)
     arc = athird*DACOS(q/(p*pst))
     csarc = DCOS(arc)
     AL3DNEW = pst*(DSQRT(3.d0*(1.d0-csarc*csarc)) - csarc) + 0.5d0     ! case (4b)
  endif

  if (cc > 0.5d0)  AL3DNEW = 1.d0 - AL3DNEW
  
! compute alpha for the real coefficients  
  AL3DNEW = AL3DNEW + DMIN1(0.d0,nr(1)) + DMIN1(0.d0,nr(2)) + DMIN1(0.d0,nr(3))
  
END FUNCTION AL3DNEW

!===============================================================================
! DESCRIPTION OF FUNCTION CC3DNEW:
! compute in the unit cube the volume cut by the plane
! [nr]*[x] = nr1*x1 + nr2*x2 + nr3*x3 = alpha
! where the normal is pointing outward from the reference phase
! INPUT: normal coefficients (nr1,nr2,nr3) and plane constant alpha 
! NOTE : the normal coefficients MUST satisfy the relation |nr1|+|nr2|+|nr3|=1
!-------------------------------------------------------------------------------
FUNCTION CC3DNEW(nr,alpha)

  IMPLICIT NONE
  REAL(8), INTENT(IN):: nr(3),alpha
  REAL(8) :: CC3DNEW
  REAL(8) :: al,alh,np1,np2,np3,m1,m2,m3,m12,mm,denom,top
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
  alh = DMIN1(al,1.d0-al)                              ! symmetry: 0 < alh < 1/2
  m12 = m1 + m2
  mm = DMIN1(m12,m3)
  denom = DMAX1(6.d0*m1*m2*m3,1.0d-50)

! 1: al<=m1; 2: m1<=al<=m2; 3: m2<=al<=mm; 4: mm<=al<=1/2 (a:m12<=m3; b:m3<m12)) 
  if (alh <= m1) then
    CC3DNEW = alh*alh*alh/denom                                          ! case (1)
  else if (alh <= m2) then
    CC3DNEW = 0.5d0*alh*(alh-m1)/(m2*m3) +  m1*m1*m1/denom               ! case (2)
  else if (alh <= mm) then
    top = alh*alh*(3.d0*m12-alh) + m1*m1*(m1-3.d0*alh)              
    CC3DNEW = (top + m2*m2*(m2-3.d0*alh))/denom                          ! case (3)
  else if (m12 <= m3) then
    CC3DNEW = (alh - 0.5d0*m12)/m3                                      ! case (4a)
  else
    top = alh*alh*(3.d0-2.d0*alh) + m1*m1*(m1-3.d0*alh)             
    CC3DNEW = (top + m2*m2*(m2-3.d0*alh) + m3*m3*(m3-3.d0*alh))/denom   ! case (4b)
  endif
  
  if (al > 0.5d0) CC3DNEW = 1.d0 - CC3DNEW

END FUNCTION CC3DNEW

!===============================================================================
! DESCRIPTION OF FUNCTION FL3DNEW:
! compute in the right hexahedron starting at (x01,x02,x03) and of sides 
! (dx1,dx2,dx3) the volume cut by the plane
! [nr]*[x] = nr1*x1 + nr2*x2 + nr3*x3 = alpha
! INPUT: normal coefficients (nr1,nr2,nr3), plane constant alpha, starting
!        point (x01,x02,x03), sides (dx1,dx2,dx3)
! NOTE : the normal coefficients MUST satisfy the relation |nr1|+|nr2|+|nr3|=1
!-------------------------------------------------------------------------------
FUNCTION FL3DNEW(nr,alpha,x0,dx)

  IMPLICIT NONE
  REAL(8), INTENT(IN):: nr(3),x0(3),dx(3),alpha
  REAL(8) :: FL3DNEW
  REAL(8) :: al,almax,alh,np1,np2,np3,m1,m2,m3,m12,mm,denom,frac,top
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
  alh = DMIN1(al,1.d0-al)                              ! symmetry: 0 < alh < 1/2


! normalized equation: m1*y1 + m2*y2 + m3*y3 = alh, with m1 <= m2 <= m3
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
  denom = DMAX1(6.d0*m1*m2*m3,1.0d-50)

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
     FL3DNEW = frac*top
  else
     FL3DNEW = (1.d0-frac)*top
  endif

END FUNCTION FL3DNEW

!===============================================================================
! DESCRIPTION OF FUNCTION AREA3D:
! compute in the unit cube the area cut by the plane
! [nr]*[x] = nr1*x1 + nr2*x2 + nr3*x3 = alpha
! INPUT: normal coefficients (nr1,nr2,nr3) and volume fraction cc 
! NOTE : the normal coefficients MUST satisfy the relation |nr1|+|nr2|+|nr3|=1;
!        the cut area A is invariant with respects to reflections, ordering  
!        of the coefficients and midpoint (i.e. A(alpha) = A(1-alpha))
!-------------------------------------------------------------------------------
FUNCTION AREA3D(nr,cc)

  IMPLICIT NONE
  REAL(8), INTENT(IN):: nr(3),cc
  REAL(8) :: AREA3D
  REAL(8) :: cch,al,c00,c01,c02,c03,np1,np2,np3
  REAL(8) :: m1,m2,m3,m12,numer,denom,p,pst,q,arc,csarc
  REAL(8), PARAMETER :: athird=1.d0/3.d0 
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
  cch = DMIN1(cc,1.d0-cc)                              ! symmetry: 0 < cch < 1/2
  denom = DMAX1(6.d0*m1*m2*m3,1.d-50)                           ! get cch ranges
  c01 = m1*m1*m1/denom
  c02  = c01 + 0.5d0*(m2-m1)/m3
  m12 = m1 + m2
  if (m12 <= m3) then
     c03 = 0.5d0*m12/m3
  else
     numer = m3*m3*(3.d0*m12-m3) + m1*m1*(m1-3.d0*m3) + m2*m2*(m2-3.d0*m3)
     c03 = numer/denom
  endif
  c00 = 3.d0*DSQRT(m1*m1 + m2*m2 + m3*m3)/denom

! 1: C<=C1; 2: C1<=C<=C2; 3: C2<=C<=C3; 4: C3<=C<=1/2 (a: m12<=m3; b: m3<m12)) 
  if (cch <= c01) then
     al = (denom*cch)**athird                                         
     AREA3D = c00*al*al                                               ! case (1)
  else if (cch <= c02) then 
    al = 0.5d0*(m1 + DSQRT(m1*m1 + 8.d0*m2*m3*(cch - c01)))           
    AREA3D = c00*m1*(2.d0*al-m1)                                      ! case (2)
  else if (cch <= c03) then
     p = 2.d0*m1*m2                                                   
     q = 1.5d0*m1*m2*(m12 - 2.d0*m3*cch)
     pst = DSQRT(p)
     arc = athird*DACOS(q/(p*pst))
     csarc = DCOS(arc)
     al = pst*(DSQRT(3.d0*(1.d0-csarc*csarc)) - csarc) + m12
     AREA3D = c00*(-al*al + m1*(2.d0*al-m1) + m2*(2.d0*al-m2))        ! case (3)
  else if (m12 <= m3) then                                  
     al = m3*cch + 0.5d0*m12                                         
     AREA3D = c00*2.d0*m1*m2                                         ! case (4a) 
  else                                  
     p = m12*m3 + m1*m2 - 0.25d0                                     
     q = 1.5d0*m1*m2*m3*(0.5d0-cch)
     pst = DSQRT(p)
     arc = athird*DACOS(q/(p*pst))
     csarc = DCOS(arc)
     al = pst*(DSQRT(3.d0*(1.d0-csarc*csarc)) - csarc) + 0.5d0 
     AREA3D = c00*(2.d0*al*(1.d0-al) - (m1*m1 + m2*m2 + m3*m3))      ! case (4b)
  endif

END FUNCTION AREA3D

!===============================================================================
! DESCRIPTION OF FUNCTION CENT3D:
! compute in the unit cube the coordinates of the centroid of the area 
! cut by the plane
! [nr]*[x] = nr1*x1 + nr2*x2 + nr3*x3 = alpha
! INPUT: normal coefficients nr1, nr2, nr3 and volume fraction cc 
! NOTE : the normal coefficients MUST satisfy the relation |nr1|+|nr2|+|nr3|=1;
!        the centroid coordinates change with the reflections required for 
!        negative normal coefficients, the midpoint symmetry and the ordering 
!        of the coefficients    
!-------------------------------------------------------------------------------
subroutine CENT3D(nr,cc,centroid)

  IMPLICIT NONE
  REAL(8), INTENT(IN):: nr(3),cc
  REAL(8), DIMENSION(3) :: centroid(3),ctd0
  REAL(8) :: cch,al,c01,c02,c03,np1,np2,np3
  REAL(8) :: m1,m2,m3,m12,numer,denom,p,pst,q,arc,csarc,top,bot
  REAL(8), PARAMETER :: athird=1.d0/3.d0 
  INTEGER, DIMENSION(3) :: index
  INTRINSIC DMAX1,DMIN1,DSQRT,DACOS,DCOS

  np1 = DABS(nr(1))                                 ! need positive coefficients
  np2 = DABS(nr(2))
  np3 = DABS(nr(3))
! order positive coefficients and get index ordering
  if (np1 <= np2) then                             
    m1 = np1
    m3 = np2
    index(1) = 1                                            
    index(3) = 2
  else
    m1 = np2
    m3 = np1
    index(1) = 2
    index(3) = 1
  endif
  
  if (np3 < m1) then
    m2 = m1
    m1 = np3
    index(2) = index(1)
    index(1) = 3
  else if (np3 >= m3) then
    m2 = m3
    m3 = np3
    index(2) = index(3)
    index(3) = 3
  else
    m2 = np3
    index(2) = 3
  endif
  cch = DMIN1(cc,1.d0-cc)                              ! symmetry: 0 < cch < 1/2
  denom = DMAX1(6.d0*m1*m2*m3,1.d-50)                           ! get cch ranges
  c01 = m1*m1*m1/denom
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
    al = (denom*cch)**athird 
    ctd0(1) = athird*al/m1
    ctd0(2) = athird*al/m2
    ctd0(3) = athird*al/m3                                            ! case (1)
  else if (cch <= c02) then 
    al = 0.5d0*(m1 + DSQRT(m1*m1 + 8.d0*m2*m3*(cch - c01)))    
    top = m1*m1 + 3.*al*(al-m1)
    bot = 3.*(2.*al-m1)
    ctd0(1) = (3.*al-2.*m1)/bot
    ctd0(2) = top/(m2*bot)
    ctd0(3) = top/(m3*bot)                                            ! case (2)
  else if (cch <= c03) then
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

! back from symmetry
  if (cc > 0.5d0) then
    ctd0(1) = 1.d0 - ctd0(1)
    ctd0(2) = 1.d0 - ctd0(2)
    ctd0(3) = 1.d0 - ctd0(3)
  endif

! get correct indexing
  CENTROID(index(1)) = ctd0(1)                              
  CENTROID(index(2)) = ctd0(2)
  CENTROID(index(3)) = ctd0(3)

! take care of negative coefficients
  if (nr(1) < 0.d0) CENTROID(1) = 1.d0 - CENTROID(1) 
  if (nr(2) < 0.d0) CENTROID(2) = 1.d0 - CENTROID(2)
  if (nr(3) < 0.d0) CENTROID(3) = 1.d0 - CENTROID(3)

END subroutine CENT3D
!=================================================================================================
! ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
! PROGRAM TO FIND alpha IN: m1 x1 + m2 x2 + m3 x3 = alpha,
! GIVEN m1+m2+m3=1 (all > 0) AND THE VOLUMETRIC FRACTION cc
! ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
function al3d(b1,b2,b3,cc)
  !***
  implicit none
  real(8) m1,m2,m3,cc,b1,b2,b3,tmp,pr,ch,mm,m12
  real(8) p,p12,q,teta,cs,al3d
  real(8) untier,v1,v2,v3
  parameter (untier=1.d0/3.d0)
  intrinsic dmax1,dmin1,dsqrt,dacos,dcos
  !***  
  !     (1) order coefficients: m1<m2<m3; (2) get ranges: v1<v2<v3;
  !     (3) limit ch (0.d0 < ch < 0.5d0); (4) calculate alpha
  !*(1)* 
  m1 = dmin1(b1,b2)
  m3 = dmax1(b1,b2)
  m2 = b3
  if (m2 .lt. m1) then
     tmp = m1
     m1 = m2
     m2 = tmp
  else if (m2 .gt. m3) then
     tmp = m3
     m3 = m2
     m2 = tmp
  endif
  !*(2)*
  m12 = m1 + m2 
  pr  = DMAX1(6.d0*m1*m2*m3,1.d-50)
  V1  = m1*m1*m1/pr
  V2  = V1 + 0.5d0*(m2-m1)/m3
  if (m3 .LT. m12) then
     mm = m3
     V3 = (m3*m3*(3.d0*m12-m3) + m1*m1*(m1-3.d0*m3) +&
          m2*m2*(m2-3.d0*m3))/pr
  else
     mm = m12
     V3 = 0.5d0*mm/m3
  endif
  !*(3)*
  ch = DMIN1(cc,1.d0-cc)
  !*(4)*      
  if (ch .LT. V1) then
     !***         AL3D = cbrt(pr*ch)
     AL3D = (pr*ch)**UNTIER
  else if (ch .LT. V2) then
     AL3D = 0.5d0*(m1 + DSQRT(m1*m1 + 8.d0*m2*m3*(ch-V1)))
  else if (ch .LT. V3) then
     p = 2.d0*m1*m2
     q = 1.5d0*m1*m2*(m12 - 2.d0*m3*ch)
     p12 = DSQRT(p)
     teta = DACOS(q/(p*p12))/3.d0
     cs = DCOS(teta)
     AL3D = p12*(DSQRT(3.d0*(1.d0-cs*cs)) - cs) + m12
  else if (m12 .LT. m3) then
     AL3D = m3*ch + 0.5d0*mm
  else 
     p = m1*(m2+m3) + m2*m3 - 0.25d0
     q = 1.5d0*m1*m2*m3*(0.5d0-ch)
     p12 = DSQRT(p)
     teta = DACOS(q/(p*p12))/3.0
     cs = DCOS(teta)
     AL3D = p12*(DSQRT(3.d0*(1.d0-cs*cs)) - cs) + 0.5d0
  endif

  if (cc .GT. 0.5d0)  AL3D = 1.d0 - AL3D
  !***
  return
end function al3d
! ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
! PROGRAM TO FIND THE "CUT VOLUME" V0 GIVEN r0, dr0 AND
! m1 x1 + m2 x2 + m3 x3 = alpha
! ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
function fl3d(m1,m2,m3,alpha,r0,dr0)
  !***
  implicit none
  real(8) m1,m2,m3,alpha,r0,dr0,fl3D
  real(8) al,al0,n1,n2,n3,b1,b2,b3,b12,bm,tmp,pr
  INTRINSIC DMAX1,DMIN1,DABS
  !***
  !     (1) move origin to r0 along r ;  (2) reflect parallelepiped;
  !     (3) limit alpha (0<= al0 <=0.5); (4) order coefficients: b1<b2<b3;
  !     (5) calculate volume (NOTE: it is assumed:s0=t0=0; ds0=dt0=1.)
  !*(1)*
  al = alpha - m1*r0
  !*(2)*
  al = al + DMAX1(0.d0,-m1*dr0)+DMAX1(0.d0,-m2)+DMAX1(0.d0,-m3)
  tmp = DABS(m1)*dr0 + DABS(m2) + DABS(m3)
  n1 = DABS(m1)/tmp
  n2 = DABS(m2)/tmp
  n3 = DABS(m3)/tmp
  al = DMAX1(0.d0,DMIN1(1.d0,al/tmp))
  !*(3)*
  al0 = DMIN1(al,1.d0-al)
  !*(4)*
  b1 = DMIN1(n1*dr0,n2)
  b3 = DMAX1(n1*dr0,n2)
  b2 = n3
  if (b2 .LT. b1) then
     tmp = b1
     b1 = b2
     b2 = tmp
  else if (b2 .GT. b3) then
     tmp = b3
     b3 = b2
     b2 = tmp
  endif
  b12 = b1 + b2
  bm = DMIN1(b12,b3)
  pr = DMAX1(6.d0*b1*b2*b3,1.0d-50)
  !*5*     
  if (al0 .LT. b1) then
     tmp = al0*al0*al0/pr
  else if (al0 .LT. b2) then
     tmp = 0.5d0*al0*(al0-b1)/(b2*b3) +  b1*b1*b1/pr
  else if (al0 .LT. bm) then
     tmp = (al0*al0*(3.d0*b12-al0) + b1*b1*(b1-3.d0*al0) +&
          b2*b2*(b2-3.d0*al0))/pr
  else if (b12 .LT. b3) then
     tmp = (al0 - 0.5d0*bm)/b3
  else
     tmp = (al0*al0*(3.d0-2.d0*al0) + b1*b1*(b1-3.d0*al0) +&
          b2*b2*(b2-3.d0*al0) + b3*b3*(b3-3.d0*al0))/pr
  endif

  if (al .LE. 0.5d0) then
     FL3D = tmp*dr0
  else
     FL3D = (1.d0-tmp)*dr0
  endif
  !***  
  return
end function fl3d
