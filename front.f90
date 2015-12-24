!=================================================================================================
!=================================================================================================
! Module module_front:
! Includes the data for the points and the linked list definitions and procedures for adding and
! deleting from the list.
!-------------------------------------------------------------------------------------------------
module module_front
  use module_2phase
  implicit none
  save
  private
  public :: InitFront, backup_front_write, backup_front_read, print_fronts, CalcVolume, &
            CalcSurfaceTension, Front2GridVector, AdvanceFront2, GetFront, DistributeFront, &
            InitConditionFront, StoreOldFront, AverageFront, RegridFront, SmoothFront, &
            CorrectVolume, &
            amin, amax, aspmax, sigma, smooth, nsmooth, MaxPoint, MaxElem, MaxFront, &
            maxErrorVol, FrontProps, nregrid, &
            test_frdroplet

  integer :: MaxPoint, MaxElem, MaxFront, FirstEmptyPoint, FirstEmptyElem, FirstEmptyFront, &
             FrontLength, FrontFirst, NumEmptyPoint, NumEmptyElem, NumEmptyFront, NumLocPoint

  integer, dimension( : ), allocatable :: PointLength,  PointFirst,     ElemLength, ElemFirst, &
                                          LocalPointIndex, LocalElemIndex
  integer, dimension(:,:), allocatable :: FrontConnect, ElemConnect,    ElemNgbr, &
                                          PointConnect, ElemCorner
  real(8), dimension(:,:), allocatable :: PointCoords,  PointCoordsOld, FrontForce, gradIFront, &
                                          PointCoordsBuff,  PointCoordsBuffo, FrontProps
  real(8), dimension( : ), allocatable :: SurfaceTension ! , rad, xc, yc, zc, vol
! Two-phase properties now defined in Two-Phase module

  integer, dimension( : ), allocatable :: GlobPointIndex,   GlobElemIndex
  integer, dimension(:,:), allocatable :: ElemNgbrBuff,     ElemCornerBuff
  real(8) :: amin, amax, aspmax, maxErrorVol
  real(8), dimension(:,:), allocatable :: smin, smax
  real(8), dimension( 3 ) :: mysmin, mysmax
  logical :: smooth
  integer :: nsmooth, nregrid
  integer, allocatable, dimension(:,:) :: request
  integer, allocatable, dimension(:,:,:) :: status
  logical :: test_frdroplet
  contains
!=================================================================================================
!=================================================================================================
! subroutine to add object NewObject to a linked list. The object comes from the unused part of 
! the array.
!-------------------------------------------------------------------------------------------------
  subroutine add_obj(NewObject, connect, NumObjects, FirstObject, NumEmptyObjects, &
                     FirstEmptyObject, MaxObjects)
    implicit none
    integer :: NewObject, NumObjects, FirstObject, NumEmptyObjects, FirstEmptyObject, MaxObjects
    integer :: connect(MaxObjects,2),m
    if(NumEmptyObjects.eq.0)then 
      write(*,*)'NO MORE OBJECTS IN LIST 1' 
      return 
    endif 
    if(NumObjects==0)then
      NewObject   = FirstEmptyObject
      FirstObject = NewObject
      FirstEmptyObject = connect(FirstEmptyObject,1)
      connect(NewObject,1) = NewObject
      connect(NewObject,2) = NewObject
      NumObjects = NumObjects+1
      NumEmptyObjects = NumEmptyObjects-1
    else
      NewObject = FirstEmptyObject
      m = connect(FirstObject,2)
      FirstEmptyObject = connect(FirstEmptyObject,1)
      connect(NewObject,1) = FirstObject
      connect(m,1) = NewObject
      connect(NewObject,2) = m
      connect(FirstObject,2) = NewObject
      NumObjects = NumObjects+1
      NumEmptyObjects = NumEmptyObjects-1
    endif
    return 
  end subroutine add_obj
!=================================================================================================
!=================================================================================================
! subroutine to delete object OldObject from a linked list. The object remains in the unused part
! of the array.
!-------------------------------------------------------------------------------------------------
  subroutine delete_obj(OldObject, connect, NumObjects, FirstObject, NumEmptyObjects, &
                        FirstEmptyObject, MaxObjects)
    implicit none
    integer OldObject, NumObjects, FirstObject, NumEmptyObjects, FirstEmptyObject, MaxObjects
    integer connect(MaxObjects,2)
    if(OldObject.eq.firstObject) firstObject = connect(firstObject,1)
    if(NumObjects.eq.1) firstObject = 0
    NumObjects = NumObjects-1
    connect(connect(OldObject,2),1)=connect(OldObject,1)
    connect(connect(OldObject,1),2)=connect(OldObject,2)
!    connect(FirstEmptyObject,2) = OldObject
    connect(OldObject,1) = FirstEmptyObject
    connect(OldObject,2) = 0
    FirstEmptyObject = OldObject
    NumEmptyObjects = NumEmptyObjects+1
    return
  end subroutine delete_obj
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
  subroutine InitFront
    use module_grid
    implicit none
    include 'mpif.h'
    integer :: i, ierr, root !, numProcess

    test_frdroplet=.true.
    FirstEmptyPoint = 1;  NumEmptyPoint = MaxPoint
    FirstEmptyElem  = 1;  NumEmptyElem  = MaxElem
    FirstEmptyFront = 1;  NumemptyFront = MaxFront
    allocate( FrontConnect(MaxFront,2), ElemConnect   (MaxElem,2), PointConnect   (MaxPoint,2), &
              PointLength (MaxFront  ), ElemCorner    (MaxElem,3), PointCoords    (MaxPoint,3), &
              PointFirst  (MaxFront  ), ElemNgbr      (MaxElem,3), LocalPointIndex(MaxPoint  ), &
              ElemLength  (MaxFront  ), LocalElemIndex(MaxElem  ), FrontForce     (3,MaxPoint), &
              ElemFirst   (MaxFront  ), SurfaceTension(MaxPoint)  )
    allocate( GlobPointIndex(MaxPoint), GlobElemIndex(MaxElem), ElemNgbrBuff(MaxElem,3), &
              PointCoordsBuff(9,MaxPoint),PointCoordsBuffo(9,MaxPoint),ElemCornerBuff(MaxElem,3),&
              PointCoordsOld(MaxPoint,3), gradIFront(3,MaxPoint) )
    FrontLength = 0;  FrontFirst  = 0
    PointLength = 0;  PointFirst  = 0
    ElemLength  = 0;  ElemFirst   = 0
    ElemCorner  = 0;  ElemNgbr    = 0
!   Set initial connectivity
    PointConnect = 0
    do i=1,MaxPoint-1
      PointConnect(i,1)=i+1
    enddo
    ElemConnect = 0
    do i=1,MaxElem-1
      ElemConnect(i,1)=i+1
    enddo
    FrontConnect = 0
    do i=1,MaxFront-1
      FrontConnect(i,1)=i+1
    enddo

    ! smin ans smax are used to decide which points on the front belongs to which process
    if(rank<nPdomain)then
      mysmin = (/ dfloat(is)-0.5d0 , dfloat(js)-0.5d0 , dfloat(ks)-0.5d0 /)
      mysmax = (/ dfloat(ie)+0.5d0 , dfloat(je)+0.5d0 , dfloat(ke)+0.5d0 /)
    else
      mysmin = 0d0
      mysmax = 0d0
    endif
    allocate(smin(1:3,0:numProcess-1),smax(1:3,0:numProcess-1))

    root = nPdomain
    call MPI_Gather(mysmin, 3, MPI_DOUBLE_PRECISION, smin, 3, MPI_DOUBLE_PRECISION, &
                    root , MPI_Comm_Active, ierr)
    call MPI_Gather(mysmax, 3, MPI_DOUBLE_PRECISION, smax, 3, MPI_DOUBLE_PRECISION, &
                    root , MPI_Comm_Active, ierr)

    allocate(request(4,0:nPdomain-1),status(MPI_STATUS_SIZE,4,0:nPdomain-1))
  end subroutine InitFront
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
  subroutine InitConditionFront
    use module_grid
    implicit none
    integer :: i, npoints, front
    do i=1, NumBubble
      npoints=int(rad(i)*3.1415/2.*float(Nx)/xLength/(sqrt(amin*amax)/1d0)/2.)*2 !must be even
      print*,'npoints=',npoints
      call add_obj(front, FrontConnect, FrontLength, FrontFirst, NumEmptyFront, &
                          FirstEmptyFront, MaxFront)
      call add_sphere(front,FrontProps(5,i),FrontProps(6,i),FrontProps(7,i),rad(i),npoints)
      print*,'sphere added'
    enddo
  end subroutine InitConditionFront
!=================================================================================================
!=================================================================================================
! subroutine add_sphere: Adds a sphere to the front 
!-------------------------------------------------------------------------------------------------
  subroutine add_sphere(front,xc,yc,zc,radin,nps)
    use module_grid
    implicit none
    integer :: front,nps
    real(8) :: xc,yc,zc,radin, ee, radx,rady,radz
    real(8) :: pi, dph, phi, theta, pt(4*nps*nps+2,3), xp(3)
    integer :: iq,i2,i1,iip,ist,iie,ia,ib,ic,id,icp(4*2*nps*nps,3),ine(4*2*nps*nps,3),iqq,ne,np, &
               indpt(4*nps*nps+2), indel(4*2*nps*nps), point, elem, i
    pi=4d0*ATAN(1d0)
    dph=0.5d0*pi/dfloat(nps)
!    ee is nonzero to create deformed bubbles
    ee=0.0d0
!   set north and south pole----CHECK
    radx=radin*sqrt(1.+excentricity(1))
    rady=radin*sqrt(1.+excentricity(2))
    radz=radin*sqrt(1.+excentricity(3))
    pt(1,1)=xc
    pt(1,2)=yc
    pt(1,3)=zc+radz
    pt(2,1)=xc
    pt(2,2)=yc
    pt(2,3)=zc-radz
    do iq = 1,4
      do i2 = 1,nps
        do i1 = 1,nps
! set points
          iip = 2+(iq-1)*nps*nps+(i2-1)*nps+i1
          phi = dph*dfloat(i1-i2)
!          phi = phi+(phi-pi/2d0)*(phi+pi/2d0)*0.1
          ist=i2-1
          if((i1-i2) .lt. 0)ist=i1-1    
          theta = 0.5d0*pi*(dfloat(iq-1)+dfloat(ist)/dfloat(nps-abs(i1-i2)) )
          pt(iip,1)=xc+radx*cos(phi)*cos(theta)
          pt(iip,2)=yc+rady*cos(phi)*sin(theta)
          pt(iip,3)=zc+radz*sin(phi)
! set elements
          iie=2*i1+2*nps*(i2-1)+2*nps*nps*(iq-1)
          ia=iip
          ib=iip+nps
          ic=ib+1
          id=ia+1
          if(i1 .eq. nps)then
            iqq=iq
            if(iqq .eq. 4)iqq=0
            ic=2+iqq*nps*nps+nps+1-i2
            id=ic+1
          end if
          if(i2 .eq. nps) then
            iqq=iq
            if(iqq .eq. 4)iqq=0
            ib=2+iqq*nps*nps+(nps+1-i1)*nps+1
            ic=ib-nps
          end if
          if((i1 .eq. nps) .and. (i2 .eq. 1)) id=1
          if((i2 .eq. nps) .and. (i1 .eq. 1)) ib=2
          icp(iie-1,1)=ia
          icp(iie-1,2)=ib
          icp(iie-1,3)=ic
          icp(iie  ,1)=ia
          icp(iie  ,2)=ic
          icp(iie  ,3)=id
          ine(iie-1,1)=iie-2
          ine(iie-1,2)=iie+2*nps
          ine(iie-1,3)=iie
          ine(iie  ,1)=iie-1
          ine(iie  ,2)=iie+1
          ine(iie  ,3)=iie-2*nps-1
          if(i1 .eq. 1  ) then
            iqq=iq-1
            if(iqq .eq. 0)iqq=4
            ine(iie-1,1)=iqq*2*nps*nps-2*i2+1
          end if
          if(i1 .eq. nps) then
            iqq=iq
            if(iqq .eq. 4)iqq=0
            ine(iie  ,2)=iqq*2*nps*nps+2*(nps+1-i2)
          end if
          if(i2 .eq. 1  ) then
            iqq=iq-1
            if(iqq .eq. 0)iqq=4
            ine(iie  ,3)=iqq*2*nps*nps-2*nps*(i1-1)
          end if 
          if(i2 .eq. nps) then
            iqq=iq
            if(iqq .eq. 4)iqq=0
            ine(iie-1,2)=iqq*2*nps*nps+2*nps*(nps-i1)+1
          end if
!          elprop(iie-1,1)=prop1
!          elprop(iie-1,2)=prop2
!          elprop(iie,1)=prop1
!          elprop(iie,2)=prop2
!          iprop(iie-1)=is
!          iprop(iie)=is
        enddo
      enddo
    enddo
    np=4*nps*nps+2
    ne=4*2*nps*nps

    if(PointLength(front)/=0 .or. ElemLength(front)/=0) &
      call pariserror("Error: Only one sphere per front is allowed.")

    do i = 1,np
      call add_obj(point, PointConnect, PointLength(front), PointFirst(front), &
                   NumEmptyPoint, FirstEmptyPoint, MaxPoint)
      indpt(i) = point

      xp = pt(i,:)
      call MapCoords(xp,point)

      SurfaceTension(point) = sigma
    enddo

    do i = 1,ne
      call add_obj(elem, ElemConnect, ElemLength(front), ElemFirst(front), &
                   NumEmptyElem, FirstEmptyElem, MaxElem)
      indel(i) = elem
      ElemCorner(elem,:) = indpt(icp(i,:))
    enddo
    elem = ElemFirst(front)
    do i = 1,ne
      ElemNgbr  (elem,:) = indel(ine(i,:))
      elem = ElemConnect(elem,1)
    enddo
  end subroutine add_sphere
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
  subroutine StoreOldFront
    implicit none
    integer :: front, ifr, i, point
    front = FrontFirst
    do ifr = 1, FrontLength
      point = PointFirst(front)
      do i = 1, PointLength(front)
        PointCoordsOld(point,:)=PointCoords(point,:)
        point = PointConnect(point,1)
      enddo
      front = FrontConnect(front,1)
    enddo
  end subroutine StoreOldFront
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
  subroutine AverageFront
    implicit none
    integer :: front, ifr, i, point
    front = FrontFirst
    do ifr = 1, FrontLength
      point = PointFirst(front)
      do i = 1, PointLength(front)
        PointCoords(point,:)=0.5d0*(PointCoords(point,:)+PointCoordsOld(point,:))
        point = PointConnect(point,1)
      enddo
      front = FrontConnect(front,1)
    enddo
  end subroutine AverageFront
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
  subroutine backup_front_write(time,iTimeStep)
    use module_grid
    use module_IO
    implicit none
    integer :: iunit, front, elem, point, ifr, i, j, iTimeStep
    real(8), intent(in) :: time
    character(len=100) :: filename
    filename = trim(out_path)//'/backup_front'
    call system('mv '//trim(filename)//' '//trim(filename)//'.old')
    iunit=7
    OPEN(UNIT=iunit,FILE=trim(filename))
    write(iunit,"('time=',es25.16e3,' iTime= ',I8,' Fronts= ',I8)")time, itimestep, FrontLength
    front = FrontFirst
    do ifr = 1, FrontLength
      write(iunit,"('NODES=',I8,' ,ELEMENTS=', I8, ' ,Volume=', e25.16e3)") &
                    PointLength(front), ElemLength(front), FrontProps(2,front)
      point = PointFirst(front)
      do i = 1, PointLength(front)
        write(iunit,'(I6,3e17.8e3)') point, (PointCoords(point,j),j=1,3) !, SurfaceTension(point)
        point = PointConnect(point,1)
      enddo
      elem = ElemFirst(front)
      do i = 1, ElemLength(front)
        call write_integer(iunit,elem,ElemCorner(elem,1) &
                                     ,ElemCorner(elem,2) &
                                     ,ElemCorner(elem,3) )
        elem = ElemConnect(elem,1)
      enddo
      elem = ElemFirst(front)
      do i = 1, ElemLength(front)
        call write_integer(iunit,ElemNgbr(elem,1) &
                                ,ElemNgbr(elem,2) &
                                ,ElemNgbr(elem,3) )
        elem = ElemConnect(elem,1)
      enddo
      front = FrontConnect(front,1)
    enddo
    close(iunit)
  end subroutine backup_front_write
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
  subroutine backup_front_read(time,iTimeStep)
    use module_grid
    use module_IO
    implicit none
    integer :: iunit, front, elem, point, ifr, i, j, FrontLength1, ElemLength1, PointLength1, &
               elem1, point1, indpt(MaxPoint), indel(MaxPoint), ind(3), iTimeStep
    real(8), intent(out) :: time
    real(8) :: volume
    iunit=7
    OPEN(UNIT=iunit,FILE=trim(out_path)//'/backup_front',status='old',action='read')
    read(iunit,"('time=',es25.16e3,' iTime= ',I8,' Fronts= ',I8)")time, itimestep, FrontLength1
    do ifr = 1, FrontLength1
      call add_obj(front, FrontConnect, FrontLength, FrontFirst, NumEmptyFront, &
                         FirstEmptyFront, MaxFront)
      !read(iunit,"('NODES=',I6,' ,ELEMENTS=', I6)")PointLength1, ElemLength1
      read(iunit,"('NODES=',I8,' ,ELEMENTS=', I8, ' ,Volume=', e25.16e3)") &
                    PointLength1, ElemLength1, volume
      FrontProps(2,front) = volume
      do i = 1, PointLength1
        call add_obj(point, PointConnect, PointLength(front), PointFirst(front), &
                     NumEmptyPoint, FirstEmptyPoint, MaxPoint)
        read(iunit,*) point1, (PointCoords(point,j),j=1,3)
        indpt(point1) = point
        SurfaceTension(point) = sigma
      enddo
      do i = 1, ElemLength1
        call add_obj(elem, ElemConnect, ElemLength(front), ElemFirst(front), &
                     NumEmptyElem, FirstEmptyElem, MaxElem)
        read(iunit,*) elem1, ind(1:3)
        ElemCorner(elem,1:3) = indpt(ind(1:3))
        indel(elem1) = elem
      enddo
      elem = ElemFirst(front)
      do i = 1, ElemLength(front)
        read(iunit,*) ind(1:3)
        ElemNgbr(elem,1:3) = indel(ind(1:3))
        elem = ElemConnect(elem,1)
      enddo
      front = FrontConnect(front,1)
    enddo
    close(iunit)
  end subroutine backup_front_read
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
  subroutine print_fronts(nf,time)
    use module_IO
    implicit none
    integer :: nf
    real(8) :: time
    if(output_format==1) call print_fronts1(nf,time) !tecplot format
    if(output_format==2) call print_fronts2(nf,time) !vtk format
    if(output_format==3) call print_fronts3(nf) !silo format
  end subroutine print_fronts
!-------------------------------------------------------------------------------------------------
  subroutine print_fronts1(nf,time)
    use module_grid
    use module_IO
    implicit none
    integer :: nf, iunit, front, elem, point, ifr, i
    real(8) :: time, xp(3)
    iunit=7
    OPEN(UNIT=iunit,FILE=trim(out_path)//'/front'//int2text(nf,3)//'.dat')
    write(iunit,*)'variables = "x","y","z"' !,"u","v","w"'
    front = FrontFirst
    do ifr = 1, FrontLength
      write(iunit,"('ZONE solutiontime=',1PG15.7E2,' DATAPACKING=POINT, ZONETYPE=FETRIANGLE')")time
      write(iunit,"('NODES=',I6,' ,ELEMENTS=', I6)")PointLength(front), ElemLength(front)
      point = PointFirst(front)
      do i = 1, PointLength(front)
        call GetCoords(xp,point)
        write(iunit,'(6e16.8e2)') xp !, gradIFront(:,point)
        LocalPointIndex(point)=i
        point = PointConnect(point,1)
      enddo
      elem = ElemFirst(front)
      do i = 1, ElemLength(front)
        call write_integer(iunit,LocalPointIndex(ElemCorner(elem,1)) &
                                ,LocalPointIndex(ElemCorner(elem,2)) &
                                ,LocalPointIndex(ElemCorner(elem,3)) )
        elem = ElemConnect(elem,1)
      enddo
      front = FrontConnect(front,1)
    enddo
    close(iunit)
  end subroutine print_fronts1
!-------------------------------------------------------------------------------------------------
  subroutine print_fronts2(nf,time)
    use module_grid
    use module_IO
    implicit none
    integer :: nf,iunit, front, elem, point, ifr, i
    real(8) :: time, xp(3)
    iunit=7
    OPEN(UNIT=iunit,FILE=trim(out_path)//'/front'//int2text(nf,3)//'.vtk')
    front = FrontFirst
    do ifr = 1, FrontLength

      write(iunit,10)
      write(iunit,11)time
      write(iunit,12)
      write(iunit,13)
      write(iunit,15)PointLength(front)
10    format('# vtk DataFile Version 2.0')
11    format('grid, time ',F16.8)
12    format('ASCII')
13    format('DATASET POLYDATA')
15    format('POINTS ',I17,' float' )

      point = PointFirst(front)
      do i = 1, PointLength(front)
        call GetCoords(xp,point)
        WRITE (iunit,255) xp
255     FORMAT (' ',3E13.5)
        LocalPointIndex(point)=i-1
        point = PointConnect(point,1)
      enddo

      write(iunit,16) ElemLength(front),4*ElemLength(front)
16    format('POLYGONS ',2I17)
      elem = ElemFirst(front)
      do i = 1, ElemLength(front)
        call write_integer(iunit,3,LocalPointIndex(ElemCorner(elem,1)) &
                                  ,LocalPointIndex(ElemCorner(elem,2)) &
                                  ,LocalPointIndex(ElemCorner(elem,3)) )
! 266     FORMAT (' 3',3I6)
        elem = ElemConnect(elem,1)
      enddo
      front = FrontConnect(front,1)
    enddo
    close(iunit)
  end subroutine print_fronts2

!-------------------------------------------------------------------------------------------------
  subroutine print_fronts3(nf)
    use module_grid
    use module_IO
    implicit none
#ifdef HAVE_SILO
    include 'silo_f9x.inc'
#endif   

    integer :: nf, front, elem, point, ifr, i, j
    real(8) :: xp0(3)

#ifdef HAVE_SILO
    integer :: dbfile, err, ierr, lname
    integer :: ndims, nzones, nnodes, Lnodelist
    integer :: shapesize, shapecount
    integer, parameter :: TotalPoint=4000000, TotalElem=8000000
    real   , dimension(TotalPoint) :: xp, yp, zp
    integer, dimension(3*TotalElem)  :: nodelist
    character(len=30) :: outfile

    nnodes = 0
    nzones = 0
    Lnodelist=0

    front = FrontFirst
    do ifr = 1, FrontLength

       point = PointFirst(front)
       do i = 1, PointLength(front)
          call GetCoords(xp0,point)
          xp(i +nnodes) = REAL(xp0(1))
          yp(i +nnodes) = REAL(xp0(2))
          zp(i +nnodes) = REAL(xp0(3))
          LocalPointIndex(point) = i +nnodes
          point = PointConnect(point,1)
       enddo

       elem = ElemFirst(front)
       do i = 1, ElemLength(front)
          do j = 1, 3
             Lnodelist=Lnodelist+1
             nodelist(Lnodelist) = LocalPointIndex(ElemCorner(elem,j))
          enddo
          elem = ElemConnect(elem,1)
       enddo

       nnodes = nnodes + PointLength(front) 
       nzones = nzones + ElemLength(front)
       front = FrontConnect(front,1)
    enddo

    ndims = 3
    shapesize = 3
    shapecount=nzones

    outfile=trim(out_path)//'/front'//int2text(nf,3)//'.silo'
    outfile=trim(adjustl(outfile))
    lname=len_trim(outfile)

    ierr= dbcreate(outfile, lname, 0, DB_LOCAL, "Front3D", 7, DB_PDB, dbfile)
    if(dbfile== -1) then
       write(*,*) 'Could not create Silo file!\n'
       stop
    endif

    err = dbputzl2(dbfile, "zonelist", 8, nzones, ndims, nodelist, Lnodelist, 1, 0, 0, &
                   DB_ZONETYPE_TRIANGLE, shapesize, shapecount, 1, DB_F77NULL, ierr)
    err = dbputum(dbfile, "front", 5, ndims, xp, yp, zp, "X", 1, "Y", 1, "Z", 1, DB_FLOAT, &
                  nnodes, nzones, "zonelist", 8, DB_F77NULL, 0, DB_F77NULL, ierr)
    ierr = dbclose(dbfile)
#endif   

  end subroutine print_fronts3
!=================================================================================================
!=================================================================================================
! writes 3 or 4 integer number to an output file without extra blanks
subroutine write_integer(iunit,n1,n2,n3,n4)
  implicit none
  integer :: iunit,n1,n2,n3,n4
  optional :: n4
  character(len=20) :: c1,c2,c3,c4
  
  write(c1,*)n1; c1=adjustl(c1)
  write(c2,*)n2; c2=adjustl(c2)
  write(c3,*)n3; c3=adjustl(c3)

  if(present(n4))then
    write(c4,*)n4; c4=adjustl(c4)
    write(iunit,'(A)') trim(c1)//' '//trim(c2)//' '//trim(c3)//' '//trim(c4)
  else
    write(iunit,'(A)') trim(c1)//' '//trim(c2)//' '//trim(c3)
  endif

  return
end subroutine write_integer
!=================================================================================================
!=================================================================================================
! 1 -> volume=V=int(dV)
! 2 -> initial volume
! 3 -> volume error
! 4 -> area=int(dA)
! First moments of inertia:
! 5 -> xc=int(x*dV)/V
! 6 -> yc=int(y*dV)/V
! 7 -> zc=int(z*dV)/V
! Second moments of inertia:
! 8 -> x3c=int(x*x*dV)-xc^2*V
! 9 -> y3c=int(y*y*dV)-yc^2*V
! 10-> z3c=int(z*z*dV)-zc^2*V
! 11-> xy3c=int(x*y*dV)-xc*yc*V
! 12-> xz3c=int(x*z*dV)-xc*zc*V
! 13-> yz3c=int(y*z*dV)-yc*zc*V
!-------------------------------------------------------------------------------------------------
  subroutine CalcVolume
    implicit none
    real(8) :: xx1(3), xx2(3), xx3(3), p1x,p1y,p1z, p2x,p2y,p2z, p3x,p3y,p3z, xc,yc,zc, x1,y1,z1,&
               x2,y2,z2, rx,ry,rz, ar, rnx,rny,rnz, xxr,yyr,zzr, M(3,3), a,b,c,d, x0(3)
    complex(8) :: roots(3)
    real(8), parameter :: r3=1.0d0/3.0d0
    integer :: elem, i, p1, p2, p3, front, ifr

    front = FrontFirst
    do ifr = 1,FrontLength
      p1 = PointFirst(front)
      call GetCoords(x0,p1)
      FrontProps(1,front) = 0d0
      FrontProps(3:14,front) = 0d0
      elem = ElemFirst(front)
      do i = 1,ElemLength(front)
        p1 = ElemCorner(elem,1)
        p2 = ElemCorner(elem,2)
        p3 = ElemCorner(elem,3)
        call GetCoords(xx1,p1)
        call GetCoords(xx2,p2)
        call GetCoords(xx3,p3)
        p1x = xx1(1)-x0(1)
        p1y = xx1(2)-x0(2)
        p1z = xx1(3)-x0(3)
        p2x = xx2(1)-x0(1)
        p2y = xx2(2)-x0(2)
        p2z = xx2(3)-x0(3)
        p3x = xx3(1)-x0(1)
        p3y = xx3(2)-x0(2)
        p3z = xx3(3)-x0(3)
        xc = r3*(p1x+p2x+p3x)
        yc = r3*(p1y+p2y+p3y)
        zc = r3*(p1z+p2z+p3z)
        x1 = p2x-p1x
        y1 = p2y-p1y
        z1 = p2z-p1z
        x2 = p3x-p1x
        y2 = p3y-p1y
        z2 = p3z-p1z
!
! components of n.dA and n
!
        rx = 0.5d0*(y1*z2-y2*z1)
        ry = 0.5d0*(z1*x2-z2*x1)
        rz = 0.5d0*(x1*y2-x2*y1)
        ar = sqrt(rx*rx+ry*ry+rz*rz)
        rnx = rx/ar
        rny = ry/ar
        rnz = rz/ar
!
! increments due to element elem
!
        FrontProps( 1,front) = FrontProps( 1,front)+xc*rx+yc*ry+zc*rz
        FrontProps( 4,front) = FrontProps( 4,front)+ar
        xxr = xc*xc*rx
        yyr = yc*yc*ry
        zzr = zc*zc*rz
        FrontProps( 5,front) = FrontProps( 5,front)+xxr
        FrontProps( 6,front) = FrontProps( 6,front)+yyr
        FrontProps( 7,front) = FrontProps( 7,front)+zzr
        FrontProps( 8,front) = FrontProps( 8,front)+xc*xxr
        FrontProps( 9,front) = FrontProps( 9,front)+yc*yyr
        FrontProps(10,front) = FrontProps(10,front)+zc*zzr
        FrontProps(11,front) = FrontProps(11,front)+yc*xxr+xc*yyr
        FrontProps(12,front) = FrontProps(12,front)+zc*xxr+xc*zzr
        FrontProps(13,front) = FrontProps(13,front)+zc*yyr+yc*zzr

        elem  =  ElemConnect(elem,1)
        !ke = elcon(ke,kf)
      end do
      FrontProps( 1,front) = FrontProps( 1,front) * r3
      FrontProps( 3,front) = FrontProps( 3,front) - FrontProps(2,front)
      FrontProps( 5,front) = FrontProps( 5,front) * 0.5d0  / FrontProps(1,front) + x0(1)
      FrontProps( 6,front) = FrontProps( 6,front) * 0.5d0  / FrontProps(1,front) + x0(2)
      FrontProps( 7,front) = FrontProps( 7,front) * 0.5d0  / FrontProps(1,front) + x0(3)
      FrontProps( 8,front) = FrontProps( 8,front) * r3
      FrontProps( 9,front) = FrontProps( 9,front) * r3
      FrontProps(10,front) = FrontProps(10,front) * r3
      FrontProps(11,front) = FrontProps(11,front) * 0.25d0
      FrontProps(12,front) = FrontProps(12,front) * 0.25d0
      FrontProps(13,front) = FrontProps(13,front) * 0.25d0

! calculate centered versions of second moments of inertia
      x1 = FrontProps(5,front) - x0(1)
      y1 = FrontProps(6,front) - x0(2)
      z1 = FrontProps(7,front) - x0(3)
      FrontProps( 8,front) = FrontProps( 8,front)-x1*x1*FrontProps(1,front)
      FrontProps( 9,front) = FrontProps( 9,front)-y1*y1*FrontProps(1,front)
      FrontProps(10,front) = FrontProps(10,front)-z1*z1*FrontProps(1,front)
      FrontProps(11,front) = FrontProps(11,front)-x1*y1*FrontProps(1,front)
      FrontProps(12,front) = FrontProps(12,front)-x1*z1*FrontProps(1,front)
      FrontProps(13,front) = FrontProps(13,front)-y1*z1*FrontProps(1,front)

      ! calculate bubble deformation
      M(1,1)=FrontProps(8 ,front)
      M(2,2)=FrontProps(9 ,front)
      M(3,3)=FrontProps(10,front)
      M(1,2)=FrontProps(11,front)
      M(2,1)=FrontProps(11,front)
      M(1,3)=FrontProps(12,front)
      M(3,1)=FrontProps(12,front)
      M(2,3)=FrontProps(13,front)
      M(3,2)=FrontProps(13,front)
      a = -1
      b = tr(M)
      c = 0.5*(tr(matmul(M,M))-(tr(M))**2)
      d = det(M)
      
      call cubic_solve(a,b,c,d,roots)
      a=real(roots(1))
      b=real(roots(2))
      c=real(roots(3))

      d=min(a,b,c)
      a=max(a,b,c)

      FrontProps(14,front) = sqrt(a/d)
  	  
      front = FrontConnect(front,1)
    enddo

  end subroutine CalcVolume
!=================================================================================================
!=================================================================================================
! subroutine CorrectVolume
! Corrects the volume of each bubble by putting a source at the centroid of each bubble
! proportional to its volume loss.
!-------------------------------------------------------------------------------------------------
  subroutine CorrectVolume
    use module_grid
    implicit none
    real(8) :: frontCenter(3), x1(3), dr(3), factor, r3, shift(3)
    integer :: point, i, front, ifr, ip, jfr, front2
    real(8) :: pi=3.141592653

    front = FrontFirst
    do ifr = 1, FrontLength
      if(abs(FrontProps(1,front)/FrontProps(2,front)-1d0)>maxErrorVol)then
        factor = (FrontProps(2,front)-FrontProps(1,front))/4d0/pi
        frontCenter = (/FrontProps(5,front), FrontProps(6,front), FrontProps(7,front)/)

        front2 = FrontFirst
        do jfr = 1, FrontLength
        shift(1) = xLength*floor((FrontProps(5,front2)-FrontProps(5,front))/xLength+0.5)
        shift(2) = yLength*floor((FrontProps(6,front2)-FrontProps(6,front))/yLength+0.5)
        shift(3) = zLength*floor((FrontProps(7,front2)-FrontProps(7,front))/zLength+0.5)

          point = PointFirst(front2)
          do ip = 1, PointLength(front2)
            call GetCoords(x1,point)
            dr = x1-frontCenter-shift
            r3 = (dr(1)**2+dr(2)**2+dr(3)**2)**1.5
            dr = factor * dr / r3

            i = floor(PointCoords(point,1));   i = Ng+modulo(i-Ng,Nx)
            PointCoords(point,1) = PointCoords(point,1) + dr(1)/dx(i+1)

            i = floor(PointCoords(point,2));   i = Ng+modulo(i-Ng,Ny)
            PointCoords(point,2) = PointCoords(point,2) + dr(2)/dy(i+1)

            i = floor(PointCoords(point,3));   i = Ng+modulo(i-Ng,Nz)
            PointCoords(point,3) = PointCoords(point,3) + dr(3)/dz(i+1)

            point = PointConnect(point,1)
          enddo

          front2 = FrontConnect(front2,1)
        enddo

      endif
      front = FrontConnect(front,1)
    enddo
  end subroutine CorrectVolume
!=================================================================================================
!=================================================================================================
! subroutine MapCoords
! Maps the point xp(:) in physical domain to PointCoords(point,:) in the front domain.
!-------------------------------------------------------------------------------------------------
  subroutine MapCoords(xp,point)
    use module_grid
    implicit none
    real(8) :: xp(3)
    integer :: point,i1, Ndomain
    Ndomain = floor(xp(1)/xLength)
    xp(1) = xp(1)-xLength*Ndomain
    do i1=Ng,Nx+Ng-1
      if(xp(1)>xh(i1+1))cycle
      PointCoords(point,1)=0.5+( dfloat(i1)*(xh(i1+1)-xp(1)) + dfloat(i1+1)*(xp(1)-xh(i1)) )/dx(i1+1)
      exit
    enddo
    PointCoords(point,1) = PointCoords(point,1)+dfloat(Nx*Ndomain)
    Ndomain = floor(xp(2)/yLength)
    xp(2) = xp(2)-yLength*Ndomain
    do i1=Ng,Ny+Ng-1
      if(xp(2)>yh(i1+1))cycle
      PointCoords(point,2)=0.5+( dfloat(i1)*(yh(i1+1)-xp(2)) + dfloat(i1+1)*(xp(2)-yh(i1)) )/dy(i1+1)
      exit
    enddo
    PointCoords(point,2) = PointCoords(point,2)+dfloat(Ny*Ndomain)
    Ndomain = floor(xp(3)/zLength)
    xp(3) = xp(3)-zLength*Ndomain
    do i1=Ng,Nz+Ng-1
      if(xp(3)>zh(i1+1))cycle
      PointCoords(point,3)=0.5+( dfloat(i1)*(zh(i1+1)-xp(3)) + dfloat(i1+1)*(xp(3)-zh(i1)) )/dz(i1+1)
      exit
    enddo
    PointCoords(point,3) = PointCoords(point,3)+dfloat(Nz*Ndomain)
  end subroutine MapCoords
!=================================================================================================
!=================================================================================================
! subroutine GetCoords
! Converts the mapped location of front points to the physical domain
!-------------------------------------------------------------------------------------------------
  subroutine GetCoords(xp,point)
    use module_grid
    implicit none
    real(8) :: xp(3)
    integer :: point,i,i1
    i = floor(PointCoords(point,1)-0.5);   i1 = Ng+modulo(i-Ng,Nx)
    !if(i>Nx .or. i<1)call pariserror("Error: GetCoords"); 
    xp(1)=xh(i1)+dx(i1+1)*(PointCoords(point,1)-dfloat(i)-0.5)+xLength*floor(dfloat(i-Ng)/dfloat(Nx))
    i = floor(PointCoords(point,2)-0.5);   i1 = Ng+modulo(i-Ng,Ny)
    !if(i>Ny .or. i<1)call pariserror("Error: GetCoords"); 
    xp(2)=yh(i1)+dy(i1+1)*(PointCoords(point,2)-dfloat(i)-0.5)+yLength*floor(dfloat(i-Ng)/dfloat(Ny))
    i = floor(PointCoords(point,3)-0.5);   i1 = Ng+modulo(i-Ng,Nz)
    !if(i>Nz .or. i<1)call pariserror("Error: GetCoords"); 
    xp(3)=zh(i1)+dz(i1+1)*(PointCoords(point,3)-dfloat(i)-0.5)+zLength*floor(dfloat(i-Ng)/dfloat(Nz))
  end subroutine GetCoords
!=================================================================================================
!=================================================================================================
! subroutine CalcSurfaceTension
! Calculates the surface tension and the color function gradient on the front grid
!-------------------------------------------------------------------------------------------------
  subroutine CalcSurfaceTension
    use module_flow
    implicit none
    real(8) :: l
    real(8), dimension(3) :: xcen, xn, xn1, sm, x1, x2, x3, x21, x32 !, area
    real(8), dimension(3,3) :: xm, xs
    integer :: front, elem, ifr, point, i, p1, p2, p3

    ! reset the surface tension force
    front = FrontFirst
    do ifr = 1, FrontLength
      point = PointFirst(front)
      do i = 1, PointLength(front)
        FrontForce(:,point)=0d0
        gradIFront(:,point)=0d0
        point = PointConnect(point,1)
      enddo
      front = FrontConnect(front,1)
    enddo

    ! add up the surface tension forces
    front = FrontFirst
    do ifr = 1, FrontLength
      elem = ElemFirst(front)
      do i = 1, ElemLength(front)
        p1=ElemCorner(elem,1)      
        p2=ElemCorner(elem,2)
        p3=ElemCorner(elem,3)

        call GetCoords(x1,p1)
        call GetCoords(x2,p2)
        call GetCoords(x3,p3)
        xm(:,1) = 0.5d0*(x2+x3)
        xm(:,2) = 0.5d0*(x3+x1)
        xm(:,3) = 0.5d0*(x1+x2)
        xcen = (x1+x2+x3)/3d0

        !normal vector to the element
        x21 = x2-x1
        x32 = x3-x2
        xn  = cross(x21,x32)
        l   = dsqrt(xn(1)**2+xn(2)**2+xn(3)**2)
        xn1 = xn/l !unit length

        x21 = xcen-xm(:,1);   xs(:,1) = cross(x21,xn1)
        x21 = xcen-xm(:,2);   xs(:,2) = cross(x21,xn1)
        x21 = xcen-xm(:,3);   xs(:,3) = cross(x21,xn1)

        !average surface tension
        sm(1) = (SurfaceTension(2)+SurfaceTension(3))/2.4d0+SurfaceTension(1)/6d0
        sm(2) = (SurfaceTension(3)+SurfaceTension(1))/2.4d0+SurfaceTension(2)/6d0
        sm(3) = (SurfaceTension(1)+SurfaceTension(2))/2.4d0+SurfaceTension(3)/6d0
        FrontForce(:,p1) = FrontForce(:,p1) + (xs(:,3)*sm(3) - xs(:,2)*sm(2))
        FrontForce(:,p2) = FrontForce(:,p2) + (xs(:,1)*sm(1) - xs(:,3)*sm(3))
        FrontForce(:,p3) = FrontForce(:,p3) + (xs(:,2)*sm(2) - xs(:,1)*sm(1))

        ! density gradient terms
        gradIFront(:,p1) = gradIFront(:,p1) + xn/6d0
        gradIFront(:,p2) = gradIFront(:,p2) + xn/6d0
        gradIFront(:,p3) = gradIFront(:,p3) + xn/6d0

        elem = ElemConnect(elem,1)
      enddo
      front = FrontConnect(front,1)
    enddo
  end subroutine CalcSurfaceTension
!=================================================================================================
!=================================================================================================
! subroutine Front2GridVector
! Maps the surface tension and color function gradient from the front to the flow grid
!-------------------------------------------------------------------------------------------------
  subroutine Front2GridVector(fx1, fy1, fz1, fx2, fy2, fz2)
    use module_grid
    use module_BC 
    use module_tmpvar
    implicit none
    include 'mpif.h'
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax) :: fx1, fy1, fz1, fx2, fy2, fz2
    real(8) :: xp,yp,zp, xt,yt,zt, xth,yth,zth, ffx1,ffy1,ffz1, ffx2,ffy2,ffz2, &
               drx,dry,drz,drxh,dryh,drzh, pd2
    integer :: point,ir,jr,kr,irh,jrh,krh,ii,jj,kk,iih,jjh,kkh,i1,j1,k1 !,i,front,ifr
    integer :: req(12),sta(MPI_STATUS_SIZE,12), ierr
    pd2 = 3.14159265358979323846/2d0

    ! add up the surface tension forces onto the fixed grid
!    front = FrontFirst
!    do ifr = 1, FrontLength
!      point = PointFirst(front)
!      do i = 1, PointLength(front)
    do point = 1, NumLocPoint

          xp   = PointCoordsBuff(1,point)
          yp   = PointCoordsBuff(2,point)
          zp   = PointCoordsBuff(3,point)
          ffx1 = PointCoordsBuff(4,point)
          ffy1 = PointCoordsBuff(5,point)
          ffz1 = PointCoordsBuff(6,point)
          ffx2 = PointCoordsBuff(7,point)
          ffy2 = PointCoordsBuff(8,point)
          ffz2 = PointCoordsBuff(9,point)

          ir = floor(xp);   irh = floor(xp-0.5d0);   xt = xp;   xth = xp-0.5d0
          jr = floor(yp);   jrh = floor(yp-0.5d0);   yt = yp;   yth = yp-0.5d0
          kr = floor(zp);   krh = floor(zp-0.5d0);   zt = zp;   zth = zp-0.5d0

          if(bdry_cond(1)==1) then
            ir  = is-1+modulo(ir-is+1,Nx)
            irh = is-1+modulo(irh-is+1,Nx)
            xt  = dfloat(is-1)+modulo(xp-dfloat(is-1),dfloat(Nx))
            xth = dfloat(is-1)+modulo(xp-dfloat(is-1)-0.5d0,dfloat(Nx))
          endif
          if(bdry_cond(2)==1) then
            jr  = js-1+modulo(jr-js+1,Ny)
            jrh = js-1+modulo(jrh-js+1,Ny)
            yt  = dfloat(js-1)+modulo(yp-dfloat(js-1),dfloat(Ny))
            yth = dfloat(js-1)+modulo(yp-dfloat(js-1)-0.5d0,dfloat(Ny))
          endif
          if(bdry_cond(3)==1) then
            kr  = ks-1+modulo(kr-ks+1,Nz)
            krh = ks-1+modulo(krh-ks+1,Nz)
            zt  = dfloat(ks-1)+modulo(zp-dfloat(ks-1),dfloat(Nz))
            zth = dfloat(ks-1)+modulo(zp-dfloat(ks-1)-0.5d0,dfloat(Nz))
          endif

          do i1=1,4;   do j1=1,4;   do k1=1,4
            ii  = ir -2+i1;   iih = irh-2+i1
            jj  = jr -2+j1;   jjh = jrh-2+j1
            kk  = kr -2+k1;   kkh = krh-2+k1

            drx  = 1d0 + cos((xt -float(ii ))*pd2)
            dry  = 1d0 + cos((yt -float(jj ))*pd2)
            drz  = 1d0 + cos((zt -float(kk ))*pd2)
            drxh = 1d0 + cos((xth-float(iih))*pd2)
            dryh = 1d0 + cos((yth-float(jjh))*pd2)
            drzh = 1d0 + cos((zth-float(kkh))*pd2)

            fx1(iih,jj,kk) = fx1(iih,jj,kk) + drxh*dry*drz*ffx1/64d0 /dxh(iih)/dy(jj)/dz(kk)
            fy1(ii,jjh,kk) = fy1(ii,jjh,kk) + drx*dryh*drz*ffy1/64d0 /dx(ii)/dyh(jjh)/dz(kk)
            fz1(ii,jj,kkh) = fz1(ii,jj,kkh) + drx*dry*drzh*ffz1/64d0 /dx(ii)/dy(jj)/dzh(kkh)

            fx2(iih,jj,kk) = fx2(iih,jj,kk) + drxh*dry*drz*ffx2/64d0 /dxh(iih)/dy(jj)/dz(kk)
            fy2(ii,jjh,kk) = fy2(ii,jjh,kk) + drx*dryh*drz*ffy2/64d0 /dx(ii)/dyh(jjh)/dz(kk)
            fz2(ii,jj,kkh) = fz2(ii,jj,kkh) + drx*dry*drzh*ffz2/64d0 /dx(ii)/dy(jj)/dzh(kkh)

          enddo;   enddo;   enddo
!        point = PointConnect(point,1)
!      enddo
!      front = FrontConnect(front,1)
    enddo
!-------------------------------------------------------------------------------------------------
    call ghost_xAdd(fx1,is-1,ie,1,req(1:4 ))
    call ghost_xAdd(fy1,is,ie+1,2,req(5:8 ))
    call ghost_xAdd(fz1,is,ie+1,3,req(9:12))
    call MPI_WAITALL(12,req,sta,ierr)
    fx1(is-1:is,:,:) = fx1(is-1:is,:,:) + work(is-1:is,:,:,1)
    fx1(ie-1:ie,:,:) = fx1(ie-1:ie,:,:) + work(ie-1:ie,:,:,1)
    fy1(is:is+1,:,:) = fy1(is:is+1,:,:) + work(is:is+1,:,:,2)
    fy1(ie-1:ie,:,:) = fy1(ie-1:ie,:,:) + work(ie-1:ie,:,:,2)
    fz1(is:is+1,:,:) = fz1(is:is+1,:,:) + work(is:is+1,:,:,3)
    fz1(ie-1:ie,:,:) = fz1(ie-1:ie,:,:) + work(ie-1:ie,:,:,3)

    call ghost_yAdd(fx1,js,je+1,1,req(1:4 ))
    call ghost_yAdd(fy1,js-1,je,2,req(5:8 ))
    call ghost_yAdd(fz1,js,je+1,3,req(9:12))
    call MPI_WAITALL(12,req,sta,ierr)
    fx1(:,js:js+1,:) = fx1(:,js:js+1,:) + work(:,js:js+1,:,1)
    fx1(:,je-1:je,:) = fx1(:,je-1:je,:) + work(:,je-1:je,:,1)
    fy1(:,js-1:js,:) = fy1(:,js-1:js,:) + work(:,js-1:js,:,2)
    fy1(:,je-1:je,:) = fy1(:,je-1:je,:) + work(:,je-1:je,:,2)
    fz1(:,js:js+1,:) = fz1(:,js:js+1,:) + work(:,js:js+1,:,3)
    fz1(:,je-1:je,:) = fz1(:,je-1:je,:) + work(:,je-1:je,:,3)

    call ghost_zAdd(fx1,ks,ke+1,1,req(1:4 ))
    call ghost_zAdd(fy1,ks,ke+1,2,req(5:8 ))
    call ghost_zAdd(fz1,ks-1,ke,3,req(9:12))
    call MPI_WAITALL(12,req,sta,ierr)
    fx1(:,:,ks:ks+1) = fx1(:,:,ks:ks+1) + work(:,:,ks:ks+1,1)
    fx1(:,:,ke-1:ke) = fx1(:,:,ke-1:ke) + work(:,:,ke-1:ke,1)
    fy1(:,:,ks:ks+1) = fy1(:,:,ks:ks+1) + work(:,:,ks:ks+1,2)
    fy1(:,:,ke-1:ke) = fy1(:,:,ke-1:ke) + work(:,:,ke-1:ke,2)
    fz1(:,:,ks-1:ks) = fz1(:,:,ks-1:ks) + work(:,:,ks-1:ks,3)
    fz1(:,:,ke-1:ke) = fz1(:,:,ke-1:ke) + work(:,:,ke-1:ke,3)
!-------------------------------------------------------------------------------------------------
    call ghost_xAdd(fx2,is-1,ie,1,req(1:4 ))
    call ghost_xAdd(fy2,is,ie+1,2,req(5:8 ))
    call ghost_xAdd(fz2,is,ie+1,3,req(9:12))
    call MPI_WAITALL(12,req,sta,ierr)
    fx2(is-1:is,:,:) = fx2(is-1:is,:,:) + work(is-1:is,:,:,1)
    fx2(ie-1:ie,:,:) = fx2(ie-1:ie,:,:) + work(ie-1:ie,:,:,1)
    fy2(is:is+1,:,:) = fy2(is:is+1,:,:) + work(is:is+1,:,:,2)
    fy2(ie-1:ie,:,:) = fy2(ie-1:ie,:,:) + work(ie-1:ie,:,:,2)
    fz2(is:is+1,:,:) = fz2(is:is+1,:,:) + work(is:is+1,:,:,3)
    fz2(ie-1:ie,:,:) = fz2(ie-1:ie,:,:) + work(ie-1:ie,:,:,3)

    call ghost_yAdd(fx2,js,je+1,1,req(1:4 ))
    call ghost_yAdd(fy2,js-1,je,2,req(5:8 ))
    call ghost_yAdd(fz2,js,je+1,3,req(9:12))
    call MPI_WAITALL(12,req,sta,ierr)
    fx2(:,js:js+1,:) = fx2(:,js:js+1,:) + work(:,js:js+1,:,1)
    fx2(:,je-1:je,:) = fx2(:,je-1:je,:) + work(:,je-1:je,:,1)
    fy2(:,js-1:js,:) = fy2(:,js-1:js,:) + work(:,js-1:js,:,2)
    fy2(:,je-1:je,:) = fy2(:,je-1:je,:) + work(:,je-1:je,:,2)
    fz2(:,js:js+1,:) = fz2(:,js:js+1,:) + work(:,js:js+1,:,3)
    fz2(:,je-1:je,:) = fz2(:,je-1:je,:) + work(:,je-1:je,:,3)

    call ghost_zAdd(fx2,ks,ke+1,1,req(1:4 ))
    call ghost_zAdd(fy2,ks,ke+1,2,req(5:8 ))
    call ghost_zAdd(fz2,ks-1,ke,3,req(9:12))
    call MPI_WAITALL(12,req,sta,ierr)
    fx2(:,:,ks:ks+1) = fx2(:,:,ks:ks+1) + work(:,:,ks:ks+1,1)
    fx2(:,:,ke-1:ke) = fx2(:,:,ke-1:ke) + work(:,:,ke-1:ke,1)
    fy2(:,:,ks:ks+1) = fy2(:,:,ks:ks+1) + work(:,:,ks:ks+1,2)
    fy2(:,:,ke-1:ke) = fy2(:,:,ke-1:ke) + work(:,:,ke-1:ke,2)
    fz2(:,:,ks-1:ks) = fz2(:,:,ks-1:ks) + work(:,:,ks-1:ks,3)
    fz2(:,:,ke-1:ke) = fz2(:,:,ke-1:ke) + work(:,:,ke-1:ke,3)

    !move the vectors that fall outside of a wall into the domain
    call SetVectorBC(fx1,fy1,fz1)
    call SetVectorBC(fx2,fy2,fz2)

  end subroutine Front2GridVector
!=================================================================================================
!=================================================================================================
! subroutine AdvanceFront
! Moves the points on the front by interpolating the velocity from flow to front grid
!-------------------------------------------------------------------------------------------------
  subroutine AdvanceFront2(u, v, w,color, dt)
    use module_grid
    use module_BC
    implicit none
!    real(8) :: PointCoordsBuff(9,Maxpoint)
    real(8), dimension(imin:imax,jmin:jmax,kmin:kmax) :: u, v, w, color
    real(8) :: xp,yp,zp, xt,yt,zt, xth,yth,zth, drx,dry,drz,drxh,dryh,drzh, pd2, dt, up, vp, wp, &
               colorp, xn(3),rn
    integer :: point,ir,jr,kr,irh,jrh,krh,ii,jj,kk,iih,jjh,kkh,i1,j1,k1
    pd2 = 3.14159265358979323846/2d0

    ! add up the surface tension forces onto the fixed grid
    do point = 1, NumLocPoint

          xp  = PointCoordsBuff(1,point)
          yp  = PointCoordsBuff(2,point)
          zp  = PointCoordsBuff(3,point)
          
          up=0d0;  vp=0d0;  wp=0d0;  colorp=0d0

          ir = floor(xp);   irh = floor(xp-0.5d0);   xt = xp;   xth = xp-0.5d0
          jr = floor(yp);   jrh = floor(yp-0.5d0);   yt = yp;   yth = yp-0.5d0
          kr = floor(zp);   krh = floor(zp-0.5d0);   zt = zp;   zth = zp-0.5d0

          if(bdry_cond(1)==1) then
            ir  = is-1+modulo(ir-is+1,Nx)
            irh = is-1+modulo(irh-is+1,Nx)
            xt  = dfloat(is-1)+modulo(xp-dfloat(is-1),dfloat(Nx))
            xth = dfloat(is-1)+modulo(xp-dfloat(is-1)-0.5d0,dfloat(Nx))
          endif
          if(bdry_cond(2)==1) then
            jr  = js-1+modulo(jr-js+1,Ny)
            jrh = js-1+modulo(jrh-js+1,Ny)
            yt  = dfloat(js-1)+modulo(yp-dfloat(js-1),dfloat(Ny))
            yth = dfloat(js-1)+modulo(yp-dfloat(js-1)-0.5d0,dfloat(Ny))
          endif
          if(bdry_cond(3)==1) then
            kr  = ks-1+modulo(kr-ks+1,Nz)
            krh = ks-1+modulo(krh-ks+1,Nz)
            zt  = dfloat(ks-1)+modulo(zp-dfloat(ks-1),dfloat(Nz))
            zth = dfloat(ks-1)+modulo(zp-dfloat(ks-1)-0.5d0,dfloat(Nz))
          endif

          do i1=1,4;   do j1=1,4;   do k1=1,4
            ii  = ir -2+i1;   iih = irh-2+i1
            jj  = jr -2+j1;   jjh = jrh-2+j1
            kk  = kr -2+k1;   kkh = krh-2+k1

            drx  = 1d0 + cos((xt -float(ii ))*pd2)
            dry  = 1d0 + cos((yt -float(jj ))*pd2)
            drz  = 1d0 + cos((zt -float(kk ))*pd2)
            drxh = 1d0 + cos((xth-float(iih))*pd2)
            dryh = 1d0 + cos((yth-float(jjh))*pd2)
            drzh = 1d0 + cos((zth-float(kkh))*pd2)

            up = up + drxh*dry*drz/64d0*u(iih,jj,kk) /dxh(iih)
            vp = vp + drx*dryh*drz/64d0*v(ii,jjh,kk) /dyh(jjh)
            wp = wp + drx*dry*drzh/64d0*w(ii,jj,kkh) /dzh(kkh)
            
            colorp = colorp + drx*dry*drz/64d0*color(ii,jj,kk)

          enddo;   enddo;   enddo

        ! avoid bubbles to overlap
        ! if color function is less than 0.1 move the points inward 0.01dx
        if(colorp<0.02)then
          !find the unit normal outward vector
          xn(1:3) = PointCoordsBuff(7:9,point)
          rn = sqrt(sum(xn**2))
          xn = xn / rn
          up = up - 0.01*xn(1)/dt
          vp = vp - 0.01*xn(2)/dt
          wp = wp - 0.01*xn(3)/dt
        endif

        PointCoordsBuff(1,point) = PointCoordsBuff(1,point) + up*dt
        PointCoordsBuff(2,point) = PointCoordsBuff(2,point) + vp*dt
        PointCoordsBuff(3,point) = PointCoordsBuff(3,point) + wp*dt

    enddo

    ! wall boundary condition
    if(bdry_cond(1)/=1 .and. coords(1)==0    ) then
      do point = 1, NumLocPoint
          PointCoordsBuff(1,point) = max(PointCoordsBuff(1,point), 2.51d0)
      enddo
    endif
    if(bdry_cond(4)/=1 .and. coords(1)==nPx-1) then
      do point = 1, NumLocPoint
          PointCoordsBuff(1,point) = min(PointCoordsBuff(1,point), dble(ie)+0.49d0)
      enddo
    endif
    if(bdry_cond(2)/=1 .and. coords(2)==0    ) then
      do point = 1, NumLocPoint
          PointCoordsBuff(2,point) = max(PointCoordsBuff(2,point), 2.51d0)
      enddo
    endif
    if(bdry_cond(5)/=1 .and. coords(2)==nPy-1) then
      do point = 1, NumLocPoint
          PointCoordsBuff(2,point) = min(PointCoordsBuff(2,point), dble(je)+0.49d0)
      enddo
    endif
    if(bdry_cond(3)/=1 .and. coords(3)==0    ) then
      do point = 1, NumLocPoint
          PointCoordsBuff(3,point) = max(PointCoordsBuff(3,point), 2.51d0)
      enddo
    endif
    if(bdry_cond(6)/=1 .and. coords(3)==nPz-1) then
      do point = 1, NumLocPoint
          PointCoordsBuff(3,point) = min(PointCoordsBuff(3,point), dble(ke)+0.49d0)
      enddo
    endif

  end subroutine AdvanceFront2
!=================================================================================================
!=================================================================================================
! subroutine DistributeFront
! Called by the front master, distributes the front nodes between flow processes and gets it back.
!-------------------------------------------------------------------------------------------------
  subroutine DistributeFront(dest,send_recv) !,request)
    use module_grid
    use module_BC
    
    implicit none
    include 'mpif.h'
    character(len=4) :: send_recv
    integer :: dest, request(1:2), front, point, ifr, i !, iunit, elem
    integer :: sta(MPI_STATUS_SIZE,1:2), ierr
    integer, allocatable, dimension(:), save :: ptcount !, send_type
    logical, save :: first_time=.true.
    real(8) :: xp, yp, zp
    if(first_time) then
      first_time = .false.
      allocate(ptcount(0:nPdomain) ) !,send_type(0:nPdomain) )
    endif

    if(send_recv=='send') then
      ptcount(dest)=0
      front = FrontFirst
      do ifr = 1, FrontLength
        point = PointFirst(front)
        do i = 1, PointLength(front)
          xp  = PointCoords(point,1)
          yp  = PointCoords(point,2)
          zp  = PointCoords(point,3)
          if(bdry_cond(1)==1)xp = dfloat(is)-0.5d0+modulo(xp-dfloat(is)+0.5d0,dfloat(Nx))
          if(bdry_cond(2)==1)yp = dfloat(js)-0.5d0+modulo(yp-dfloat(js)+0.5d0,dfloat(Ny))
          if(bdry_cond(3)==1)zp = dfloat(ks)-0.5d0+modulo(zp-dfloat(ks)+0.5d0,dfloat(Nz))

          if( xp>=smin(1,dest) .and. xp<smax(1,dest) .and. &
              yp>=smin(2,dest) .and. yp<smax(2,dest) .and. &
              zp>=smin(3,dest) .and. zp<smax(3,dest)     )then
            ptcount(dest)=ptcount(dest)+1
            PointCoordsBuff(1:3,ptcount(dest)) = PointCoords(point, :)
            PointCoordsBuff(4:6,ptcount(dest)) = FrontForce(:, point)
            PointCoordsBuff(7:9,ptcount(dest)) = gradIfront(:, point)
            LocalPointIndex(point)   = ptcount(dest)
            GlobPointIndex(ptcount(dest)) = point
          endif
          point = PointConnect(point,1)
        enddo
        front = FrontConnect(front,1)
      enddo
       call MPI_ISEND(PointCoordsBuff(1,1),9*ptcount(dest),MPI_DOUBLE_PRECISION, dest,0, &
                     MPI_Comm_Active,request(1),ierr)
      call MPI_ISEND(GlobPointIndex(1),ptcount(dest),MPI_INTEGER, dest,1, MPI_Comm_Active, &
                     request(2),ierr)
      call MPI_WAITALL(2,request(1:2),sta(:,1:2),ierr)

    elseif(send_recv=='recv')then

      call MPI_IRECV(PointCoordsBuff(1,1),9*ptcount(dest),MPI_DOUBLE_PRECISION, dest,0, &
                     MPI_Comm_Active,request(1),ierr)
      call MPI_IRECV(GlobPointIndex(1),ptcount(dest),MPI_INTEGER, dest,1, MPI_Comm_Active, &
                     request(2),ierr)
      call MPI_WAITALL(2,request(1:2),sta(:,1:2),ierr)
      do i=1,ptcount(dest)
        PointCoords(GlobPointIndex(i),1:3)=PointCoordsBuff(1:3,i)
      enddo

    else
      call pariserror("Error: Incorrect input for DistributeFront")
    endif
   end subroutine DistributeFront
!=================================================================================================
!=================================================================================================
! subroutine GetFront
! Called by thr flow processes, gets the front from front master and sends it back.
!-------------------------------------------------------------------------------------------------
  subroutine GetFront(send_recv)
    use module_grid
    implicit none
    save
    include 'mpif.h'
    character(len=4), intent(in) :: send_recv
    integer :: ierr, req(2), sta(MPI_STATUS_SIZE,2)

    if(send_recv=='recv')then
      ! get the front from the front master process
      call MPI_IRECV(PointCoordsBuff(1,1),9*MaxPoint   ,MPI_DOUBLE_PRECISION, nPdomain,0, &
                     MPI_Comm_Active,req(1),ierr)
      call MPI_IRECV(GlobPointIndex(1),MaxPoint   ,MPI_INTEGER, nPdomain,1, MPI_Comm_Active, &
                     req(2), ierr)
    elseif(send_recv=='send') then
      ! send the updated front back to the front master process
      call MPI_ISEND(PointCoordsBuff(1,1),9*NumLocPoint,MPI_DOUBLE_PRECISION, nPdomain,0, &
                     MPI_Comm_Active,req(1),ierr)
      call MPI_ISEND(GlobPointIndex(1),NumLocPoint,MPI_INTEGER, nPdomain,1, MPI_Comm_Active, &
                     req(2), ierr)
    elseif(send_recv=='wait') then
      ! wait for the communication to finish
      call MPI_WAITALL(2,req,sta,ierr)
      call MPI_GET_COUNT(sta(:,2), MPI_INTEGER, NumLocPoint, ierr)
    else
      call pariserror("Error: Incorrect input for DistributeFront")
    endif
  end subroutine GetFront
!=================================================================================================
!-------------------------------------------------------------------------------------------------
  subroutine SmoothFront
    implicit none
    integer :: front, elem, ifr, i, ielem, nn, im, i1, ip, ip0, i2, iter, kp, j, ll
    real(8) :: lambda, alfa
!-----------------------------------------------------
! local arrays
!    integer, parameter :: mxpt=100000
    Integer :: ntp(MaxPoint),ntpts(MaxPoint,25)
    real(8) :: delta(3),delta0(3,MaxPoint)
!-----------------------------------------------------

!      if(mxpt.lt.MaxPoint)then
!        write(*,*)'[fsmooth] warning: fix mxpt'
!        stop
!      endif

      ntp(1:MaxPoint) = 0
      ntpts(1:MaxPoint,1:25) = 0
      
      front = FrontFirst
      do ifr = 1, FrontLength
        elem = ElemFirst(front)
        do ielem = 1, ElemLength(front)

          do i=1,3
            ip=ElemCorner(elem,i)
            if(ntp(ip).eq.0) then
              do i1=1,2
                nn=ntp(ip)+1
                ntp(ip)=nn
                im=i-i1
                if(im.lt.1) im=im+3
                ip0=ElemCorner(elem,im)
!                ip0=icp(im,ke,kf)
                ntpts(ip,nn)=ip0
              enddo
            else
              nn=ntp(ip)
              do i1=1,2
                im=i-i1
                if(im.lt.1) im=im+3
                ip0=ElemCorner(elem,im)
!                ip0=icp(im,ke,kf)
                do i2=1,nn
                  if(ip0.eq.ntpts(ip,i2)) cycle !goto 10
                enddo
                nn=nn+1
                ntp(ip)=nn
                ntpts(ip,nn)=ip0
              enddo
            endif
          enddo
          elem = ElemConnect(elem,1)
        enddo

        lambda=0.1d0

        do iter=1,5 !20

          kp = PointFirst(front)
          do ll = 1, PointLength(front)
!          kp=ffp(kf)
!          do ll=1,np(kf)
            do i=1,3
              delta0(i,kp)=0.0d0
              do j=1,ntp(kp)
                im=ntpts(kp,j)
                delta0(i,kp)=delta0(i,kp)+(PointCoords(im,i)-PointCoords(kp,i))/dfloat(ntp(kp))
!                (pt(i,im,kf)-pt(i,kp,kf))/dfloat(ntp(kp))
              enddo
            enddo
            kp = PointConnect(kp,1)
!            kp=ptcon(kp,kf)
          enddo

          kp = PointFirst(front)
          do ll = 1, PointLength(front)
!          kp=ffp(kf)
!          do ll=1,np(kf)
            do i=1,3
              delta(i)=0.0d0
              alfa=0.0d0
              do j=1,ntp(kp)
                im=ntpts(kp,j)
                delta(i)=delta(i)+(delta0(i,im)-delta0(i,kp))/dfloat(ntp(kp))
                alfa=alfa+1.0d0/dble(ntp(im))
              enddo
              alfa=1.0d0+alfa/dfloat(ntp(kp))
              PointCoords(kp,i) = PointCoords(kp,i)-lambda*delta(i)/alfa
!              pt(i,kp,kf)=pt(i,kp,kf)-lambda*delta(i)/alfa
            enddo
            kp = PointConnect(kp,1)
!            kp=ptcon(kp,kf)
          enddo

        enddo

      front = FrontConnect(front,1)
    enddo


  end subroutine SmoothFront
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
  subroutine RegridFront
    implicit none
    integer :: ipass, front, elem, ifr, i, p1, p2, p3, n1, n2, n3, nn, iflg, irough
    real(8) :: s1, s2, s3, ss
    do ipass=1,6
!-------------------------------------------Add elements------------------------------------------
      front = FrontFirst
      ifr = 1
      Loop1: do while(ifr.le.FrontLength)
        elem = ElemFirst(front)
        do i = 1, ElemLength(front)
          p1=ElemCorner(elem,1)      
          p2=ElemCorner(elem,2)
          p3=ElemCorner(elem,3)
          s1=distance(p1,p2)
          s2=distance(p2,p3)
          s3=distance(p3,p1)
          n1=1
          n2=2
          n3=3
          if(s1.gt.s2)then                ! Order the lengths, s1 shortest, s3 longest 
            ss=s1;      s1=s2;      s2=ss
            nn=n1;      n1=n2;      n2=nn
          endif
          if(s2.gt.s3)then
            ss=s2;      s2=s3;      s3=ss
            nn=n2;      n2=n3;      n3=nn
          endif
          if(s1.gt.s2)then
            ss=s1;      s1=s2;      s2=ss
            nn=n1;      n1=n2;      n2=nn
          endif
          if((s3.gt.amax).or.((s1.gt.amin).and.(AspectRatio(elem).gt.aspmax)))then
            call AddElement(elem,front,n3)
            cycle Loop1                   ! Check current front from the beginning
          endif
          elem = ElemConnect(elem,1)
        enddo
        ifr = ifr+1
        front = FrontConnect(front,1)
      enddo Loop1
!-------------------------------------------Delete elements---------------------------------------
      front = FrontFirst
      ifr = 1
      Loop2: do while(ifr.le.FrontLength)
        elem = ElemFirst(front)
        do i = 1, ElemLength(front)
          p1=ElemCorner(elem,1)      
          p2=ElemCorner(elem,2)
          p3=ElemCorner(elem,3)
          s1=distance(p1,p2)
          s2=distance(p2,p3)
          s3=distance(p3,p1)
          n1=1
          n2=2
          n3=3
          if(s1.gt.s2)then                ! Order the lengths, s1 shortest, s3 longest 
            ss=s1;      s1=s2;      s2=ss
            nn=n1;      n1=n2;      n2=nn
          endif
          if(s2.gt.s3)then
            ss=s2;      s2=s3;      s3=ss
            nn=n2;      n2=n3;      n3=nn
          endif
          if(s1.gt.s2)then
            ss=s1;      s1=s2;      s2=ss
            nn=n1;      n1=n2;      n2=nn
          endif
          if(s1 .lt. amin)then
            iflg=0
            call DeleteElement(elem,front,n1,iflg)
            if(iflg.eq.0) cycle Loop2     ! Check current front from the beginning
          end if
          elem = ElemConnect(elem,1)
        enddo
        ifr = ifr+1
        front = FrontConnect(front,1)
      enddo Loop2
!-------------------------------------------Check for roughness-----------------------------------
      irough=0
      call FrontQuality(irough)
      if(irough.eq.0)return
    enddo
  end subroutine RegridFront
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
  subroutine FrontQuality(irough) 
    implicit none
    integer :: irough, front, elem, ifr, i
    real(8) :: aamin, aamax
    aamin=amin*amin*1.7/4.0
    aamax=amax*amax*1.7/4.0
    irough=0                              ! Check for rough spots
    front = FrontFirst
    do ifr = 1, FrontLength
      elem = ElemFirst(front)
      do i = 1, ElemLength(front)
        !c find normal and centroid of current element
        if(ElementArea(elem) .lt. aamin)then
!          write(*,*)' small element:', elem
          irough=irough+1
        end if
        if(ElementArea(elem) .gt. aamax)then
!          write(*,*)' large element:', elem
          irough=irough+1
        end if
        if(AspectRatio(elem) .gt. aspmax)then
!          write(*,*)' large element aspect rto:', elem
          irough=irough+1
        end if
        elem = ElemConnect(elem,1)
      enddo
      front = FrontConnect(front,1)
    enddo
!    write(*,*) ' QUALITY CHECK DONE'
  end subroutine FrontQuality
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
  subroutine DeleteElement(m,front,n1,iflg)
    implicit none
    integer :: m,front,n1,n2,n3,mc,nc1,nc2,nc3,p(8), iflg, i, j, m1, m2, m3, mc1, mc2, mo
    call nflocal(m,n1,n2,n3,mc,nc1,nc2,nc3,p)
    do i=1,3
      if( ( ElemCorner(ElemNgbr(mc,nc1),i).ne.p(1) ).and. &
          ( ElemCorner(ElemNgbr(mc,nc1),i).ne.p(4) ) ) p(7)=ElemCorner(ElemNgbr(mc,nc1),i)
      if( ( ElemCorner(ElemNgbr(mc,nc2),i).ne.p(4) ).and. &
          ( ElemCorner(ElemNgbr(mc,nc2),i).ne.p(2) ) ) p(8)=ElemCorner(ElemNgbr(mc,nc2),i)
    enddo
    ! abandon delete if the result is an element whose neighbors have a common side.
    m2=ElemNgbr(m,n2)
    m3=ElemNgbr(m,n3)
    mc1=ElemNgbr(mc,nc1)
    mc2=ElemNgbr(mc,nc2)
    do i=1,3;    do j=1,3
      if( ((ElemNgbr(mc1,i).eq.ElemNgbr(mc2,j)).and.(ElemNgbr(mc1,i).ne.mc)).or.&
          ((ElemNgbr(m2 ,i).eq.ElemNgbr(m3 ,j)).and.(ElemNgbr(m2 ,i).ne.m )) ) then
        !write(*,*)'Abandoned: Delete will result in a poor grid'
        iflg=1
        return
      endif
    enddo;    enddo
    ! set element properties, this is somewhat approximate since all element
    ! connected to p1 and p2 change their area.
    !       elprop(ine(mc,nc1),1)=elprop(ine(mc,nc1),1)+0.5*elprop(mc,1)
    !       elprop(ine(mc,nc2),1)=elprop(ine(mc,nc2),1)+0.5*elprop(mc,1)
    !       elprop(ine(m,n3),1)=elprop(ine(m,n3),1)+0.5*elprop(m,1)
    !       elprop(ine(m,n2),1)=elprop(ine(m,n2),1)+0.5*elprop(m,1)
    ! set the coordinates of the new point, incl. curvature effects.
    PointCoords(p(1),1) = 0.5d0   *( PointCoords(p(1),1)+PointCoords(p(2),1) ) + &
                          0.125d0 *( PointCoords(p(3),1)+PointCoords(p(4),1) ) - &
                          0.0625d0*( PointCoords(p(5),1)+PointCoords(p(6),1) +   &
                                     PointCoords(p(7),1)+PointCoords(p(8),1) )
    PointCoords(p(1),2) = 0.5d0   *( PointCoords(p(1),2)+PointCoords(p(2),2) ) + &
                          0.125d0 *( PointCoords(p(3),2)+PointCoords(p(4),2) ) - &
                          0.0625d0*( PointCoords(p(5),2)+PointCoords(p(6),2) +   &
                                     PointCoords(p(7),2)+PointCoords(p(8),2) )
    PointCoords(p(1),3) = 0.5d0   *( PointCoords(p(1),3)+PointCoords(p(2),3) ) + &
                          0.125d0 *( PointCoords(p(3),3)+PointCoords(p(4),3) ) - &
                          0.0625d0*( PointCoords(p(5),3)+PointCoords(p(6),3) +   &
                                     PointCoords(p(7),3)+PointCoords(p(8),3) )
    ! eliminate p2:
    call delete_obj(p(2), PointConnect, PointLength(front), PointFirst(front), &
                    NumEmptyPoint, FirstEmptyPoint, MaxPoint)
    ! eliminate m and mc:
    call delete_obj(m, ElemConnect, ElemLength(front), ElemFirst(front), &
                 NumEmptyElem, FirstEmptyElem, MaxElem)
    call delete_obj(mc, ElemConnect, ElemLength(front), ElemFirst(front), &
                 NumEmptyElem, FirstEmptyElem, MaxElem)
    ! set element pointers:
    do i=1,3
      if(ElemNgbr(ElemNgbr(m ,n2 ),i) .eq. m )ElemNgbr(ElemNgbr(m ,n2 ),i)=ElemNgbr(m ,n3 )
      if(ElemNgbr(ElemNgbr(m ,n3 ),i) .eq. m )ElemNgbr(ElemNgbr(m ,n3 ),i)=ElemNgbr(m ,n2 )
      if(ElemNgbr(ElemNgbr(mc,nc2),i) .eq. mc)ElemNgbr(ElemNgbr(mc,nc2),i)=ElemNgbr(mc,nc1)
      if(ElemNgbr(ElemNgbr(mc,nc1),i) .eq. mc)ElemNgbr(ElemNgbr(mc,nc1),i)=ElemNgbr(mc,nc2)
    enddo
    ! reset pointers to p2 to p1:
    m1=ElemNgbr(m,n2)
    do
       do i=1,3
          if(ElemCorner(m1,i).eq.p(2))then
             ElemCorner(m1,i)=p(1)
             mo=i
          end if
       enddo
       if(m1.eq.ElemNgbr(mc,nc2))exit
       m1=ElemNgbr(m1,mo)
    enddo
    m=m2
    return
  end subroutine DeleteElement
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
  subroutine AddElement(m,front,n1)
    implicit none
    integer :: m,front,n1,n2,n3,mc,nc1,nc2,nc3,p(8), i, newpt, newm, newmc, mm3, mmc1
    call nflocal(m,n1,n2,n3,mc,nc1,nc2,nc3,p)
    do i=1,3
      if( ( ElemCorner(ElemNgbr(mc,nc1),i).ne.p(1) ).and. &
          ( ElemCorner(ElemNgbr(mc,nc1),i).ne.p(4) ) ) p(7)=ElemCorner(ElemNgbr(mc,nc1),i)
      if( ( ElemCorner(ElemNgbr(mc,nc2),i).ne.p(4) ).and. &
          ( ElemCorner(ElemNgbr(mc,nc2),i).ne.p(2) ) ) p(8)=ElemCorner(ElemNgbr(mc,nc2),i)
    enddo
    ! set new point
    call add_obj(newpt, PointConnect, PointLength(front), PointFirst(front), &
                 NumEmptyPoint, FirstEmptyPoint, MaxPoint)
    ! set the coordinates of the new point, incl. curvature effects.
    PointCoords(newpt,1) = 0.5d0   *( PointCoords(p(1),1)+PointCoords(p(2),1) ) + &
                           0.125d0 *( PointCoords(p(3),1)+PointCoords(p(4),1) ) - &
                           0.0625d0*( PointCoords(p(5),1)+PointCoords(p(6),1) +   &
                                      PointCoords(p(7),1)+PointCoords(p(8),1) )
    PointCoords(newpt,2) = 0.5d0   *( PointCoords(p(1),2)+PointCoords(p(2),2) ) + &
                           0.125d0 *( PointCoords(p(3),2)+PointCoords(p(4),2) ) - &
                           0.0625d0*( PointCoords(p(5),2)+PointCoords(p(6),2) +   &
                                      PointCoords(p(7),2)+PointCoords(p(8),2) )
    PointCoords(newpt,3) = 0.5d0   *( PointCoords(p(1),3)+PointCoords(p(2),3) ) + &
                           0.125d0 *( PointCoords(p(3),3)+PointCoords(p(4),3) ) - &
                           0.0625d0*( PointCoords(p(5),3)+PointCoords(p(6),3) +   &
                                      PointCoords(p(7),3)+PointCoords(p(8),3) )
    SurfaceTension(newpt) = sigma
    ! set new elements
    call add_obj(newm, ElemConnect, ElemLength(front), ElemFirst(front), &
                 NumEmptyElem, FirstEmptyElem, MaxElem)
    call add_obj(newmc, ElemConnect, ElemLength(front), ElemFirst(front), &
                 NumEmptyElem, FirstEmptyElem, MaxElem)
    ElemCorner(m,n1)=newpt
    ElemCorner(mc,nc1)=newpt
    ElemCorner(newm,1)=newpt
    ElemCorner(newm,2)=p(3)
    ElemCorner(newm,3)=p(1)
    ElemCorner(newmc,1)=newpt
    ElemCorner(newmc,2)=p(1)
    ElemCorner(newmc,3)=p(4)
    mm3=ElemNgbr(m,n3)
    mmc1=ElemNgbr(mc,nc1)
    ElemNgbr(m,n3)=newm
    ElemNgbr(mc,nc1)=newmc
    ElemNgbr(newm,1)=m
    ElemNgbr(newm,2)=mm3
    ElemNgbr(newm,3)=newmc
    ElemNgbr(newmc,1)=newm
    ElemNgbr(newmc,2)=mmc1
    ElemNgbr(newmc,3)=mc
    do i=1,3
      if(ElemNgbr(mm3,i) .eq. m)ElemNgbr(mm3,i)=newm
      if(ElemNgbr(mmc1,i) .eq.mc)ElemNgbr(mmc1,i)=newmc
    enddo
!!c set element properties
!        iprop(newm)=iprop(m)
!        iprop(newmc)=iprop(mc)
!!c       elprop(newm,1)=0.5*elprop(m,1)
!!c       elprop(m,1)=0.5*elprop(m,1)
!!c       elprop(newmc,1)=0.5*elprop(mc,1)
!!c       elprop(mc,1)=0.5*elprop(mc,1)
    return
  end subroutine AddElement
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
  subroutine nflocal(m,n1,n2,n3,mc,nc1,nc2,nc3,p)
    implicit none
    integer :: m,n1,n2,n3,mc,nc1,nc2,nc3,p(8), i
    mc=ElemNgbr(m,n1)
    n2=n1+1;    if(n2.gt.3)n2=n2-3
    n3=n2+1;    if(n3.gt.3)n3=n3-3
    !c find the points
    p(1)=ElemCorner(m,n1)
    p(2)=ElemCorner(m,n2)
    p(3)=ElemCorner(m,n3)
    do i=1,3
      if(ElemCorner(mc,i) .eq. p(1))nc1=i
    enddo
    nc2=nc1+1;    if(nc2.gt.3)nc2=nc2-3
    nc3=nc2+1;    if(nc3.gt.3)nc3=nc3-3
    p(4)=ElemCorner(mc,nc2)
    do i=1,3
      if( ( ElemCorner(ElemNgbr(m,n2),i).ne.p(2) ).and. &
          ( ElemCorner(ElemNgbr(m,n2),i).ne.p(3) ) ) p(5)=ElemCorner(ElemNgbr(m,n2),i)
      if( ( ElemCorner(ElemNgbr(m,n3),i).ne.p(3) ).and. &
          ( ElemCorner(ElemNgbr(m,n3),i).ne.p(1) ) ) p(6)=ElemCorner(ElemNgbr(m,n3),i)
    enddo
  end subroutine nflocal
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
  function AspectRatio(elem)
    implicit none
    integer :: elem, p1,p2,p3
    real(8) :: AspectRatio, s1,s2,s3,s, area
    P1=ElemCorner(elem, 1)
    P2=ElemCorner(elem, 2)
    P3=ElemCorner(elem, 3)
    S1=distance(p1,p2)
    S2=distance(p2,p3)
    S3=distance(p3,p1)
    S=(S1+S2+S3)/3d0
    area = ElementArea(elem)
    AspectRatio =0.25d0*SQRT(3d0)*S*S/area
    return
  end function AspectRatio
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
  function ElementArea(elem)
    implicit none
    integer :: elem,p1,p2,p3
    real(8) :: ElementArea, x1, x2, y1, y2, z1, z2
    P1 = ElemCorner(elem, 1)
    P2 = ElemCorner(elem, 2)
    P3 = ElemCorner(elem, 3)
    x1 = PointCoords(p2,1)-PointCoords(p1,1)
    y1 = PointCoords(p2,2)-PointCoords(p1,2)
    z1 = PointCoords(p2,3)-PointCoords(p1,3)
    x2 = PointCoords(p3,1)-PointCoords(p1,1)
    y2 = PointCoords(p3,2)-PointCoords(p1,2)
    z2 = PointCoords(p3,3)-PointCoords(p1,3)
    ElementArea = 0.5d0*sqrt((y1*z2-y2*z1)**2+(x2*z1-x1*z2)**2+(x1*y2-x2*y1)**2)
    return
  end function ElementArea
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
  function distance(p1,p2)
    implicit none
    integer :: p1,p2
    real(8) :: distance
    distance = sqrt((PointCoords(p1,1)-PointCoords(p2,1))**2+ &
                    (PointCoords(p1,2)-PointCoords(p2,2))**2+ &
                    (PointCoords(p1,3)-PointCoords(p2,3))**2  )
  end function distance
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
  function cross(a,b)
    real(8), dimension(3) :: cross, a, b
    cross(1) = a(2)*b(3)-a(3)*b(2)
    cross(2) = a(3)*b(1)-a(1)*b(3)
    cross(3) = a(1)*b(2)-a(2)*b(1)
  end function cross
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
  function det(B)
    real(8)::B(3,3),det
    det=B(1,1)*B(2,2)*B(3,3)+B(1,2)*B(2,3)*B(3,1)+B(1,3)*B(2,1)*B(3,2) &
       -B(1,3)*B(2,2)*B(3,1)-B(1,2)*B(2,1)*B(3,3)-B(1,1)*B(2,3)*B(3,2)
  end function
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
  function tr(B)
    real(8)::B(3,3),tr
    tr=B(1,1)+B(2,2)+B(3,3)
  end function
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------
subroutine cubic_solve(a,b,c,d,roots)
implicit none
complex(8),dimension(3) :: roots
real(8) :: a,b,c,d
integer :: nroot

! ----------------------------------------------------------------------
! Solve a cubic equation where a, b, c, and d are real.
!   a*x**3 + b*x**2 + c*x + d = 0
!
! Variables used:
!   a, b, c, d  ... coefficients (input)
!   y1, y2, y3  ... three transformed solutions
!   y2r, y2i    ... real and imaginary parts of a pair of complex roots
!   x(i)        ... three (generally) complex solutions (outputdir)
!   nroot       ... number of roots
!
! Formula used are given in Tuma, "Engineering Mathematics Handbook", p7
!   (McGraw Hill, 1978).
!   Step 0: If a is 0. use the quadrati! formula to avoid dividing by 0.
!   Step 1: Calculate p and q
!           p = ( 3*c/a - (b/a)**2 ) / 3
!           q = ( 2*(b/a)**3 - 9*b*c/a/a + 27*d/a ) / 27
!   Step 2: Calculate discriminant D
!           D = (p/3)**3 + (q/2)**2
!   Step 3: Depending on the sign of D, we follow different strategy.
!           If D<0, three distinct real roots.
!           If D=0, three real roots of which at least two are equal.
!           If D>0, one real and two complex roots.
!   Step 3a: For D>0 and D=0,
!           Calculate u and v
!           u = cubic_root(-q/2 + sqrt(D))
!           v = cubic_root(-q/2 - sqrt(D))
!           Find the three transformed roots
!           y1 = u + v
!           y2 = -(u+v)/2 + i (u-v)*sqrt(3)/2
!           y3 = -(u+v)/2 - i (u-v)*sqrt(3)/2
!   Step 3b Alternately, for D<0, a trigonometri! formulation is more convenient
!           y1 =  2 * sqrt(|p|/3) * cos(phi/3)
!           y2 = -2 * sqrt(|p|/3) * cos((phi+pi)/3)
!           y3 = -2 * sqrt(|p|/3) * cos((phi-pi)/3)
!           where phi = acos(-q/2/sqrt(|p|**3/27))
!                 pi  = 3.141592654...
!   Step 4  Finally, find the three roots
!           x = y - b/a/3
!
! ----------------------------------------------------------------------
! Instructor: Nam Sun Wang
! ----------------------------------------------------------------------

! Declare variables
      complex(8) :: x(3)
      real(8) :: pi=3.141592654
      real(8) :: DD, p, q, phi, temp1, temp2, y1,y2,y3, u, v, y2r, y2i

DD=0d0;  p=0d0;  q=0d0;  phi=0d0;  temp1=0d0;  temp2=0d0;  y1=0d0; y2=0d0; y3=0d0;  u=0d0;  v=0d0;  y2r=0d0;  y2i=0d0

! Step 0: If a is 0 use the quadratic formula. -------------------------
      IF(a .eq. 0.)THEN
      if(b .eq. 0.)then
        if(c .eq. 0.)then
!         We have a non-equation; therefore, we have no valid solution
          nroot = 0
        else
!         We have a linear equation with 1 root.
          nroot = 1
          x(1) = cmplx(-d/c, 0.,kind=8)
        endif
      else
!     We have a true quadratic equation.  Apply the quadratic formula to find two roots.
      nroot = 2
        DD = c*c-4.*b*d
        if(DD .ge. 0.)then
          x(1) = cmplx((-c+sqrt(DD))/2./b, 0.,kind=8)
          x(2) = cmplx((-c-sqrt(DD))/2./b, 0.,kind=8)
        else
          x(1) = cmplx(-c/2./b, +sqrt(-DD)/2./b,kind=8)
          x(2) = cmplx(-c/2./b, -sqrt(-DD)/2./b,kind=8)
        endif
      endif

      ELSE

! Cubic equation with 3 roots
      nroot = 3

! Step 1: Calculate p and q --------------------------------------------
      p  = c/a - b*b/a/a/3.
      q  = (2.*b*b*b/a/a/a - 9.*b*c/a/a + 27.*d/a) / 27.

! Step 2: Calculate DD (discriminant) ----------------------------------
      DD = p*p*p/27. + q*q/4.

! Step 3: Branch to different algorithms based on DD -------------------
      if(DD .lt. 0.)then
!       Step 3b:
!       3 real unequal roots -- use the trigonometric formulation
        phi = acos(-q/2./sqrt(abs(p*p*p)/27.))
        temp1=2.*sqrt(abs(p)/3.)
        y1 =  temp1*cos(phi/3.)
        y2 = -temp1*cos((phi+pi)/3.)
        y3 = -temp1*cos((phi-pi)/3.)
      else
!       Step 3a:
!       1 real root & 2 conjugate complex roots OR 3 real roots (some are equal)
        temp1 = -q/2. + sqrt(DD)
        temp2 = -q/2. - sqrt(DD)
        u = abs(temp1)**(1./3.)
        v = abs(temp2)**(1./3.)
        if(temp1 .lt. 0.) u=-u
        if(temp2 .lt. 0.) v=-v
        y1  = u + v
        y2r = -(u+v)/2.
        y2i =  (u-v)*sqrt(3.)/2.
      endif

! Step 4: Final transformation -----------------------------------------
      temp1 = b/a/3.
      y1 = y1-temp1
      y2 = y2-temp1
      y3 = y3-temp1
      y2r=y2r-temp1

! Assign answers -------------------------------------------------------
      if(DD .lt. 0.)then
        x(1) = cmplx( y1,  0.,kind=8)
        x(2) = cmplx( y2,  0.,kind=8)
        x(3) = cmplx( y3,  0.,kind=8)
      elseif(DD .eq. 0.)then
        x(1) = cmplx( y1,  0.,kind=8)
        x(2) = cmplx(y2r,  0.,kind=8)
        x(3) = cmplx(y2r,  0.,kind=8)
      else
        x(1) = cmplx( y1,  0.,kind=8)
        x(2) = cmplx(y2r, y2i,kind=8)
        x(3) = cmplx(y2r,-y2i,kind=8)
      endif

    ENDIF

  roots = x

end subroutine cubic_solve
!-------------------------------------------------------------------------------------------------
end module module_front
!=================================================================================================
!=================================================================================================
!-------------------------------------------------------------------------------------------------


