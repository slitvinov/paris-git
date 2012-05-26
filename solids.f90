!=================================================================================================
! module_solids: Contains definition of variables for the solids.
!-------------------------------------------------------------------------------------------------
module module_solids
  use module_grid
  implicit none
  real(8), dimension(:,:,:), allocatable :: solids
  contains
!***********************************************************************
    SUBROUTINE initsolids()
!***
!   test the C and Fortran array indexing- for development only
!***
!      call comparef2c()
!     STOP
!*** 
!     end test section
!***
      allocate(solids(imin:imax,jmin:jmax,kmin:kmax))
      solids = 0.d0
      call inisolids(solids,nx,ny,nz)
    END SUBROUTINE initsolids
!***********************************************************************
    SUBROUTINE outfarray(carray) 
      INTEGER CARRAY(2,2) 
      
      WRITE(6,*) 'OUTFARRAY'
      WRITE(6,*) '1: ', CARRAY(1,1),' , ' , CARRAY(1,2) 
      WRITE(6,*) '2: ', CARRAY(2,1),' , ' , CARRAY(2,2) 
      WRITE(6,*) ' '
      STOP 'OUTFARRAY'
    END SUBROUTINE OUTFARRAY
end module module_solids
