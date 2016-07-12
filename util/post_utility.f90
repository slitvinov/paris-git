PROGRAM post_utility

  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INCLUDE 'silo_f9x.inc'

  INTEGER, DIMENSION(MPI_STATUS_SIZE)           :: status
  INTEGER                                       :: rank, code, fh, nx, ny, nz
  INTEGER :: npx, npy, npz
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :, :)    :: u_read
  INTEGER, PARAMETER :: array_rank=3
  INTEGER, DIMENSION(array_rank) :: shape_view_array, shape_sub_view_array, start_view_coord
  INTEGER :: is, ie, js, je, ks, ke
  INTEGER :: type_sub_array, icolor
  INTEGER, DIMENSION (3) :: coords, dims, comm_dom, comm_cart
  LOGICAL, DIMENSION (3) :: periodic
  LOGICAL :: reorder
  INTEGER(kind = MPI_OFFSET_KIND) :: initial_displacement

  LOGICAL,PARAMETER :: vof_only = .false.
  ! Silo variables************************************************************************
  REAL(8) :: deltaX
  INTEGER :: levarnames, lemeshnames, padd
  REAL(8) :: xLength, yLength, zLength, time, dt, time_start
  INTEGER :: startCounter, endCounter, counter
  INTEGER, DIMENSION(3) :: dims_mesh, dims_vof, ghoststart, ghostend
  INTEGER :: i, j, k, ierr2, dbfile, optlist, lfile_name
  !***************************************************************************************

  !
  ! Read input parameters of the program
  !
  if ( rank == 0 .and. vof_only ) write(*,*) "Warning: only plot VOF!"
  CALL readInput()
  padd = 5
  !
  ! Several parameters for MPI
  !
  dims(1) = npx
  dims(2) = npy
  dims(3) = npz
  reorder = .true.
  periodic = .false.
  icolor = 0

  !
  ! Initialise MPI
  !
  CALL MPI_INIT(code)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, code)
  CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, icolor, 0, comm_dom, code)
  CALL MPI_Cart_Create(comm_dom,array_rank,dims, &
       periodic,reorder,comm_cart,code)
  CALL MPI_Cart_Coords(comm_cart,rank,array_rank,coords,code)


  !
  ! Setting bounds, is,ie,js,je,ks,ke as a function of the core subdomain
  !
  CALL setBounds()

  !
  ! Allocate variables to read
  !
  ALLOCATE(u_read((ie-is+1), (je-js+1), (ke-ks+1)))
  u_read(:, :, :) = 0.d0

  time = time_start
  !write(*,*) SHAPE(u_read)
  DO counter=startCounter,endCounter

     !write(*, *) 'out/VTK/cvof-'//i2t(counter,padd)//'.data'

     !open (unit = 2, file ='out/VTK/test-'//i2t(rank,5)//'.out')  
     !  WRITE(2, 101)  u_read
     !close(2)    

     ! 
     ! Output the results in silo format**************************************************
     !

     CALL openSilo(counter, time)
     time = time + dt
     !CALL appendVariable()
     !CALL closeSilo()
     if (rank == 0) call progress(counter,endCounter)

  ENDDO

101 FORMAT  (E12.5)


  !
  ! Finish MPI process
  !
  CALL MPI_FINALIZE(code)

CONTAINS
  !
  ! Subroutine to read the input post processing parameters
  !
  SUBROUTINE readInput()

    OPEN(11, FILE='inputpost', STATUS='OLD')
    READ(11, *) nx
    READ(11, *) ny
    READ(11, *) nz
    READ(11, *) npx
    READ(11, *) npy
    READ(11, *) npz
    READ(11, *) startCounter
    READ(11, *) endCounter
    READ(11, *) time_start
    READ(11, *) dt
    READ(11, *) xLength
    READ(11, *) yLength
    READ(11, *) zLength
    CLOSE(11)
  END SUBROUTINE readInput

  !
  ! Subroutine to calculate the bounds of the local arrays
  !
  SUBROUTINE setBounds()

    !WRITE(*,*) 'Coords', coords(1), coords(2), coords(3)

    !
    ! Ghost vector initialisation
    !
    ghoststart = 2
    ghostend = 2
    !
    ! Setting limits in the x direction
    !
    is = (nx/npx)*coords(1) - 1
    IF(coords(1)==0 ) THEN
       is = 1
       ghoststart(1) = 0
    ENDIF
    ie = (nx/npx)*(coords(1) + 1) + 2
    IF(coords(1)== npx - 1) THEN 
       ie = nx
       ghostend(1) = 0
    ENDIF
    !
    ! Setting limits in the y direction
    !
    js = (ny/npy)*coords(2) - 1
    IF(coords(2)==0 ) THEN
       js = 1
       ghoststart(2) = 0
    ENDIF
    je = (ny/npy)*(coords(2) + 1) + 2
    IF(coords(2)== npy - 1) THEN
       je = ny
       ghostend(2) = 0
    ENDIF
    !
    ! Setting limits in the z direction
    !
    ks = (nz/npz)*coords(3) - 1
    IF(coords(3)==0 ) THEN
       ks = 1
       ghoststart(3) = 0
    ENDIF
    ke = (nz/npz)*(coords(3) + 1) + 2
    IF(coords(3)== npz - 1) THEN 
       ke = nz
       ghostend(3) = 0
    ENDIF
  END SUBROUTINE setBounds

  !
  ! Read mpi
  !
  FUNCTION read_MPI(v_name, index)
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    REAL(KIND=8), DIMENSION(ie-is+1,je-js+1,ke-ks+1) :: read_MPI
    CHARACTER(len=4), INTENT(in) :: v_name
    INTEGER, INTENT(IN) :: index 


    read_MPI = 0d0
    !
    ! Open files
    ! 
    CALL MPI_FILE_OPEN(MPI_COMM_WORLD, 'out/VTK/'//TRIM(v_name)//'-'//i2t(index,padd)//'.data', &
         MPI_MODE_RDONLY, &
         MPI_INFO_NULL, fh, code)

    !Error checking
    IF (code /= MPI_SUCCESS) THEN
       PRINT *, 'Error opening the file'
       CALL MPI_ABORT(MPI_COMM_WORLD, 2, code)
       CALL MPI_FINALIZE(code)
    END IF

    !
    ! Creation of the derived datatype type_sub_array corresponding to the matrix u without ghost
    !
    !Shape  of the array 
    shape_view_array(:)= (/nx, ny, nz/)

    !Shape of the subarray

    shape_sub_view_array(:) = (/ie - is + 1, je - js + 1, ke - ks + 1/)

    !Starting coordinates of the subarray
    start_view_coord(:) = (/ is - 1, js - 1, ks - 1 /)
    !WRITE(*, *) 'Start coords: ', start_view_coord

    !Creation of the derived datatype type_sub_array
    CALL MPI_TYPE_CREATE_SUBARRAY(array_rank, shape_view_array, shape_sub_view_array, &
         start_view_coord, MPI_ORDER_FORTRAN, MPI_REAL8, type_sub_array, code)

    !Commit type_sub_array
    CALL MPI_TYPE_COMMIT(type_sub_array, code)

    !
    !Change the file view
    !
    initial_displacement = 0
    CALL MPI_FILE_SET_VIEW(fh, initial_displacement, MPI_REAL8, &
         type_sub_array, "native", MPI_INFO_NULL, code)


    !
    ! Read the data in parallel from the binary file
    !
    CALL MPI_FILE_READ_ALL(fh, read_MPI, SIZE(read_MPI), &
         MPI_REAL8, status, code)
    !
    ! Close the binary file
    !     

    CALL MPI_FILE_CLOSE(fh, code)

  END FUNCTION read_MPI

  !
  ! Subroutine to output the loaded data in .silo format
  !
  SUBROUTINE openSilo(index, time)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: index
    INTEGER :: timestep
    INTEGER::narg,cptArg !#of arg & counter of arg
    CHARACTER(LEN=70)::name !Arg name
    CHARACTER(LEN=70)::path, file_name
    REAL(4), DIMENSION(:,:,:), ALLOCATABLE :: matrix_small
    REAL(8), DIMENSION(:), ALLOCATABLE :: x_axis, y_axis, z_axis
    INTEGER :: iee, ise, jee, jse, kee, kse
    REAL(KIND=8), INTENT(IN) :: time

    path = 'out/VTK'

    !Debugging messages
    !Write(*,*) 'Starting job in rank: ', rank

    ! Setting string lengths
    levarnames = 25 + 2*padd
    lemeshnames = 24 + 2*padd

    ! Voxel size of structured regular mesh
    deltaX = xLength/REAL(Nx)

    ! Setting limits of the domain to be output
    ise = is; jse = js; kse = ks
    iee = ie; jee = je; kee = ke 
    !ise=imin; iee=imax; if(coords(1)==0) ise=is; if(coords(1)==nPx-1) iee=ie;
    !jse=jmin; jee=jmax; if(coords(2)==0) jse=js; if(coords(2)==nPy-1) jee=je;
    !kse=kmin; kee=kmax; if(coords(3)==0) kse=ks; if(coords(3)==nPz-1) kee=ke;

    ! Debugging messages
    !WRITE(*,*) 'is In rank ', rank, is, ie, js, je, ks, ke
    !WRITE(*,*) 'Maximim In rank ', rank, imin, imax, jmin, jmax, kmin, kmax
    !WRITE(*,*) 'coords in rank ', rank, coords(1), coords(2), coords(3)
    !WRITE(*,*) 'ise In rank ', rank, ise, iee, jse, jee, kse, kee

    ! Allocating arrays

    !allocate(matrix_small(ise:iee,jse:jee,kse:kee))
    ALLOCATE(matrix_small(iee-ise+1,jee-jse+1,kee-kse+1))
    ALLOCATE(x_axis(ise:iee+1),y_axis(jse:jee+1),z_axis(kse:kee+1))

    ! Defining mesh axys
    DO i=ise,iee+1
       !x_axis(i)= REAL(i-3)*deltaX
       x_axis(i)= REAL(i-1)*deltaX
    ENDDO
    !Write(*,*) 'In rank', rank, ' x limits ', x_axis
    DO j=jse,jee+1
       !y_axis(j)= REAL(j-3)*deltaX
       y_axis(j)= REAL(j-1)*deltaX
    ENDDO
    DO k=kse,kee+1
       !z_axis(k)= REAL(k-3)*deltaX
       z_axis(k)= REAL(k-1)*deltaX
    enddo

    !WRITE(*,*) 'X axys' , x_axis, rank 
    !WRITE(*,*) 'Y axys' , y_axis, rank 
    !WRITE(*,*) 'Z axys' , z_axis, rank
    ! Defining dimensions of the mesh
    dims_mesh(1) = iee - ise + 2
    dims_mesh(2) = jee - jse + 2
    dims_mesh(3) = kee - kse + 2

    ! Defining dimensions of the cvof data
    dims_vof(1) = iee - ise + 1
    dims_vof(2) = jee - jse + 1
    dims_vof(3) = kee - kse + 1

    ! Debugging messages
    !WRITE(*,*) 'Limits In rank ', rank, iee, ise, jee, jse, kee, kse

    ! Defining ghost zones
    !ghostlow(1) = is-ise
    !ghostlow(2) = js-jse
    !ghostlow(3) = ks-kse
    !ghosttop(1) = iee-ie
    !ghosttop(2) = jee-je
    !ghosttop(3) = kee-ke

    ! Writing multi mesh file
    IF (rank == 0) CALL write_master(TRIM(path)//'/fbasic',index, time)

    ! Setting *.silo file path
    file_name = TRIM(path)//'/fbasic'//i2t(index,padd)//'-'//i2t(rank,padd)//".silo"
    ! Setting *.silo file path length
    lfile_name = 20 + 2*padd

    !Debugging message
    !write(*,*) 'Path for silo is ', file_name

    ! Generating .silo file
    ierr2 = dbcreate(TRIM(file_name), lfile_name, DB_CLOBBER, DB_LOCAL, &
         'Comment about the data', 22, DB_PDB, dbfile)

    ! Setting ghost layers
    ierr2 = dbmkoptlist(2, optlist)
    ierr2 = dbaddiopt(optlist, DBOPT_HI_OFFSET, ghostend)
    ierr2 = dbaddiopt(optlist, DBOPT_LO_OFFSET, ghoststart)


    ! Appending mesh to *.silo file
    ierr2 = dbputqm (dbfile, 'srm', 18, "x", 1, &
         "y", 1, "z", 1, x_axis, y_axis, z_axis, dims_mesh, 3, &
         DB_DOUBLE, DB_COLLINEAR, optlist, ierr2)
    !END SUBROUTINE openSilo

    !SUBROUTINE appendVariable()
    !IMPLICIT NONE
    !REAL(KIND=8), DIMENSION(ie-is+1,je-js+1,ke-ks+1), INTENT(IN) :: silo_var
    ! Appending cvof variable to *.silo file  

    u_read = read_MPI('cvof',counter)

    ierr2 = dbputqv1 (dbfile, 'cvof', 4, 'srm', 3, &
         REAL(u_read,4), dims_vof, &
         3, DB_F77NULL, 0, DB_FLOAT, DB_ZONECENT, DB_F77NULL, ierr2) 

    if ( .not. vof_only ) then 
    u_read = read_MPI('pres',counter) 

    ! Appending Pressure variable to *.silo file  
    ierr2 = dbputqv1 (dbfile, 'pres', 4, 'srm', 3, &
         REAL(u_read,4), dims_vof, &
         3, DB_F77NULL, 0, DB_FLOAT, DB_ZONECENT, DB_F77NULL, ierr2)

    u_read = read_MPI('uvel',counter)

    ! Appending u_component variable to *.silo file  
    ierr2 = dbputqv1 (dbfile, 'uvel', 4, 'srm', 3, &
         REAL(u_read,4), dims_vof, &
         3, DB_F77NULL, 0, DB_FLOAT, DB_ZONECENT, DB_F77NULL, ierr2) 

    u_read = read_MPI('vvel',counter)

    ! Appending v_component variable to *.silo file  
    ierr2 = dbputqv1 (dbfile, 'vvel', 4, 'srm', 3, &
         REAL(u_read,4), dims_vof, &
         3, DB_F77NULL, 0, DB_FLOAT, DB_ZONECENT, DB_F77NULL, ierr2)

    u_read = read_MPI('wvel',counter)

    ! Appending w_component variable to *.silo file  
    ierr2 = dbputqv1 (dbfile, 'wvel', 4, 'srm', 3, &
         REAL(u_read,4), dims_vof, &
         3, DB_F77NULL, 0, DB_FLOAT, DB_ZONECENT, DB_F77NULL, ierr2) 
    end if ! vof_only

    !END SUBROUTINE appendVariable

    !SUBROUTINE closeSilo()
    !IMPLICIT NONE
    ! Closing *.silo file		
    ierr2 = dbclose(dbfile)
  END SUBROUTINE openSilo

  !
  ! Function to transform an integer into a string of a given length
  !
  FUNCTION i2t(number,length)
    INTEGER :: number, length, i
    character(len=length) :: i2t
    character, dimension(0:9) :: num = (/'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/)
    IF(number>=10**length) WRITE(*,*) "i2t error: string is not large enough"
    DO i=1,length
       i2t(length+1-i:length+1-i) = num(mod(number/(10**(i-1)),10))
    ENDDO
  END FUNCTION i2t

  !
  ! Subroutine to write the multi root data file for visit
  !
  SUBROUTINE write_master(rootname, step, time)
    IMPLICIT NONE
    INCLUDE 'silo_f9x.inc'
    INTEGER Err, ierr, dbfile, nmesh, oldlen, m, step
    CHARACTER(*) :: rootname
    CHARACTER(LEN=50) :: file_n
    CHARACTER(LEN=40) :: fullname
    CHARACTER(LEN=levarnames), DIMENSION(npx*npy*npz) :: varnames, vec1n, vec2n
    CHARACTER(LEN=levarnames), DIMENSION(npx*npy*npz) :: vec3n, vecn, presn
    CHARACTER(LEN=lemeshnames),DIMENSION(npx*npy*npz) :: meshnames
    INTEGER, DIMENSION(npx*npy*npz) :: lmeshnames, lvarnames, meshtypes, vartypes
    INTEGER :: lfile_n, lsilon, optlist
    REAL(8) :: time

    ! Setting mesh types and variable types
    meshtypes = DB_QUAD_RECT
    vartypes = DB_QUADVAR

    lmeshnames = lemeshnames
    lvarnames = levarnames
    !WRITE(*, *) 'Processors', npx*npy*npz
    ! Setting multi mesh and variables paths
    file_n = TRIM(rootname)//i2t(step,padd)//'-'
    DO m=0,npx*npy*npz-1
       !WRITE(*,*) 'Loop n: ', m
       fullname = TRIM(file_n)//TRIM(i2t(m,padd))//'.silo:srm'
       !Debugging message
       !write(*,*) 'Paths1 ', fullname
       meshnames(m+1) = TRIM(fullname)
       varnames(m+1) = TRIM(file_n)//TRIM(i2t(m,padd))//'.silo:cvof'
       if ( .not. vof_only ) then 
       presn(m+1) = TRIM(file_n)//TRIM(i2t(m,padd))//'.silo:pres'
       vec1n(m+1) = TRIM(file_n)//TRIM(i2t(m,padd))//'.silo:uvel'
       vec2n(m+1) = TRIM(file_n)//TRIM(i2t(m,padd))//'.silo:vvel'
       vec3n(m+1) = TRIM(file_n)//TRIM(i2t(m,padd))//'.silo:wvel'
       end if ! vof_only  
    ENDDO

    ! Setting length of multi mesh file pash
    lsilon = 10 + padd

    ! Creating root file
    err = dbcreate('multi'//i2t(step,padd)//'.root', lsilon, DB_CLOBBER, DB_LOCAL, &
         "multimesh root", 14, DB_PDB, dbfile)

    IF(dbfile.eq.-1) WRITE(6,*) 'Could not create Silo file!'

    !Set the maximum string length
    oldlen = dbget2dstrlen()
    err = dbset2dstrlen(lemeshnames)

    ! Setting ghost layers, time step and physical time
    err = dbmkoptlist(2, optlist)
    err = dbaddiopt(optlist, DBOPT_CYCLE, step)
    err = dbaddiopt(optlist, DBOPT_DTIME, time)
    ! Append the multimesh object.
    err = dbputmmesh(dbfile, "srucmesh", 8, npx*npy*npz, meshnames, &
         lmeshnames, meshtypes, optlist, ierr)

    !Restore the previous value for maximum string length
    err = dbset2dstrlen(oldlen)

    ! Set the maximum string length
    oldlen = dbget2dstrlen()
    err = dbset2dstrlen(levarnames)

    ! Append the multivariable object.
    err = dbputmvar(dbfile, "cvof", 4, npx*npy*npz, varnames, lvarnames, &
         vartypes, DB_F77NULL, ierr)

    if ( .not. vof_only ) then 
    err = dbputmvar(dbfile, "p", 1, npx*npy*npz, presn, lvarnames, &
         vartypes, DB_F77NULL, ierr)

    err = dbputmvar(dbfile, "uvel", 4, npx*npy*npz, vec1n, lvarnames, &
         vartypes, DB_F77NULL, ierr)

    err = dbputmvar(dbfile, "vvel", 4, npx*npy*npz, vec2n, lvarnames, &
         vartypes, DB_F77NULL, ierr)

    err = dbputmvar(dbfile, "wvel", 4, npx*npy*npz, vec3n, lvarnames, &
         vartypes, DB_F77NULL, ierr)

    err = dbputdefvars(dbfile, 'defvars',7, 1, 'velocity',8, &
         DB_VARTYPE_VECTOR, '{uvel,vvel,wvel}', 16, DB_F77NULL, ierr)
    end if ! vof_only  

    ! Set maximum string length
    err = dbset2dstrlen(oldlen)

    ! Close file
    err = dbclose(dbfile)

    !End subroutine
  END SUBROUTINE write_master

  SUBROUTINE progress(j,max)
    IMPLICIT NONE
    INTEGER(KIND=4)::j,k,max
    CHARACTER(LEN=19)::bar="?????% |          |"
    WRITE(unit=bar(1:5),fmt="(f5.1)") (100/REAL(max))*REAL(j)
    DO k=1, INT((10/REAL(max))*REAL(j))
       bar(8+k:8+k)="*"
    ENDDO
    ! print the progress bar.
    WRITE(unit=6,fmt="(a1,a19)",advance="no") char(13), bar
    IF(j/=max) THEN
       FLUSH(unit=6)
    ELSE
       WRITE(unit=6,fmt=*)
    ENDIF
    RETURN
  END SUBROUTINE progress

END PROGRAM post_utility


