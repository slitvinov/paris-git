!=================================================================================================
! Paris-0.1
! Extended from Code: FTC3D2011 (Front Tracking Code for 3D simulations)
! and Surfer. 
! 
! Authors: Sadegh Dabiri, Gretar Tryggvason
! author for VOF extensions Stephane Zaleski (zaleski@dalembert.upmc.fr) 
! !
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
! module_averages: Contains subroutines to compute several averages of the flow
! author: Tomas Arrufat
!-------------------------------------------------------------------------------------------------
module module_averages
  implicit none 
  integer, parameter :: SIZEARR = 7
  save
  character(20) :: averages_to_do(SIZEARR), av_on_nodes(SIZEARR)
  logical :: averages_on_nodes(SIZEARR,3), do_averages= .true., output_vtk(SIZEARR)
  integer :: av_div(3), arr_div(3)
  integer :: selected_phase(SIZEARR)
  logical :: debug = .false.
!**************************************************************************************************
  contains
!**************************************************************************************************
  subroutine ComputeAverages(timeout)
    use module_solid
    use module_grid
    use module_flow
    use module_vof
    use module_io
    implicit none
    include 'mpif.h'
    integer :: i, timeout, ierr
    integer :: l, m, n, bounds(3)
    integer :: req(2),sta(MPI_STATUS_SIZE,2)
    real(8), dimension(:,:,:), allocatable :: loc_av_field, loc_av_field2, av_field1, av_field2, unitary_field
    
 
    if(do_averages) then
      if (rank==0 .and. debug) write(*,*) 'Starting averages'

      ALLOCATE(av_field1(1:arr_div(1),1:arr_div(2),1:arr_div(3)), av_field2(1:arr_div(1),1:arr_div(2),1:arr_div(3)), &
               loc_av_field(1:arr_div(1),1:arr_div(2),1:arr_div(3)), loc_av_field2(1:arr_div(1),1:arr_div(2),1:arr_div(3)), &
                unitary_field(is:ie,js:je,ks:ke))
      
      unitary_field(is:ie,js:je,ks:ke) = 1d0

      do i=1,SIZEARR
      
      av_field1(1:arr_div(1),1:arr_div(2),1:arr_div(3)) = 0d0      
      av_field2(1:arr_div(1),1:arr_div(2),1:arr_div(3)) = 0d0  
      loc_av_field(1:arr_div(1),1:arr_div(2),1:arr_div(3)) = 0d0 
      loc_av_field2(1:arr_div(1),1:arr_div(2),1:arr_div(3)) = 0d0

        if(averages_on_nodes(i,1)) then
          bounds(1) = av_div(1)
        else 
          bounds(1) = arr_div(1)
        endif
        if(averages_on_nodes(i,2)) then
          bounds(2) = av_div(2)
        else 
          bounds(2) = arr_div(2)
        endif
        if(averages_on_nodes(i,3)) then
          bounds(3) = av_div(3)
        else 
          bounds(3) = arr_div(3)
        endif
        ! Computing the averages of Darcy velocity
        if(averages_to_do(i)=='UDarcy') then
          if (rank==0 .and. debug) write(*,*) 'Attempting to compute UDarcy'
          if(dosolids) then
             call AverageField(u,loc_av_field,.true.,.false.,selected_phase(i),bounds)
          else
             call AverageField(u,loc_av_field,.false.,.false.,selected_phase(i),bounds)
          endif     
          call MPI_ALLREDUCE( loc_av_field,  av_field1, arr_div(1)*arr_div(2)*arr_div(3), &
                                       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_Domain,req(1), ierr)
          if(rank==0) then

            do n=1,bounds(3)
              do m=1,bounds(2)
                do l=1,bounds(1)
                  av_field1(l,m,n) = av_field1(l,m,n)*bounds(1)*bounds(2)*bounds(3)/(XLength*YLength*ZLength)
                enddo
              enddo
            enddo
            ! Printng results
            call print_results('UDarcyPhase'//int2text(selected_phase(i),2), av_field1, timeout, bounds)
            if(output_vtk(i))then
              call printVTK('UDarcyPhase'//int2text(selected_phase(i),2),timeout,bounds,av_field1)
            endif
          endif
        ! Computing averages of pressure
        elseif(averages_to_do(i)=='P') then
        if (rank==0 .and. debug) write(*,*) 'Attempting to compute Pressure'
          !Computing pressure when there is solid
          if(dosolids) then 
          if (debug) call outputmessage('Attempting to compute Pressure with solids')
             call AverageField(p,loc_av_field,.true.,.false.,selected_phase(i), bounds)
          if (rank==0 .and. debug) write(*,*) 'Step 1'
             call AverageField(solids,loc_av_field2,.false.,.true.,selected_phase(i),bounds)
          if (debug) call outputmessage('Step 2')
             call MPI_ALLREDUCE( loc_av_field,  av_field1, arr_div(1)*arr_div(2)*arr_div(3), &
                                          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_Domain, ierr)
          if (debug) call outputmessage('Step 3')
             call MPI_ALLREDUCE( loc_av_field2,  av_field2, arr_div(1)*arr_div(2)*arr_div(3), &
                                          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_Domain, ierr)
          if (debug) call outputmessage('Step 4')
             if(rank==0) then
               do n=1,bounds(3)
                 do m=1,bounds(2)
                   do l=1,bounds(1)
                     if(av_field2(l,m,n) > 1d-10) then
                       av_field1(l,m,n) = av_field1(l,m,n) / av_field2(l,m,n)
                     endif
                   enddo
                 enddo
               enddo           
             endif
          if (debug) call outputmessage('Step 5')
          !Computing pressure when there is no solid
          else
             if (rank==0 .and. debug) write(*,*) 'Attempting to compute Pressure with out solids'
             call AverageField(p,loc_av_field,.false.,.false.,selected_phase(i), bounds)
             if (rank==0 .and. debug) write(*,*) '1'
             call MPI_ALLREDUCE( loc_av_field,  av_field1, arr_div(1)*arr_div(2)*arr_div(3), &
                                          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_Domain, ierr)
             if (rank==0 .and. debug) write(*,*) '2'
             if(rank==0) then
               do n=1,bounds(3)
                 do m=1,bounds(2)
                   do l=1,bounds(1)
                      av_field1(l,m,n) = av_field1(l,m,n)*bounds(1)*bounds(2)*bounds(3)/(XLength*YLength*ZLength)
                   enddo
                 enddo
               enddo           
             endif       
          endif
          ! Printing results
          if(rank==0) then
             if (debug) call outputmessage('Printing results')
             call print_results('PressuPhase'//int2text(selected_phase(i),2), av_field1, timeout, bounds)
             if (debug) call outputmessage('Printing results in vtk')
             if(output_vtk(i))then
                call printVTK('PressuPhase'//int2text(selected_phase(i),2),timeout,bounds,av_field1)
             endif
          endif   
        ! Computing saturation   
        elseif(averages_to_do(i)=='Saturation') then
          if (rank==0 .and. debug) write(*,*) 'Attempting to compute Saturation'
          ! Computing saturation when there is solid.
          if(dosolids) then
             call AverageField(solids,loc_av_field,.false.,.true.,selected_phase(i), bounds)
             call AverageField(solids,loc_av_field2,.false.,.true.,0, bounds)
             call MPI_ALLREDUCE( loc_av_field,  av_field1, arr_div(1)*arr_div(2)*arr_div(3), &
                                          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_Domain,req(1), ierr)
             call MPI_ALLREDUCE( loc_av_field2,  av_field2, arr_div(1)*arr_div(2)*arr_div(3), &
                                          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_Domain,req(2), ierr)
             if(rank==0) then
               do n=1,bounds(3)
                 do m=1,bounds(2)
                   do l=1,bounds(1)
                     av_field1(l,m,n) = av_field1(l,m,n)/av_field2(l,m,n)
                   enddo
                 enddo
               enddo
             endif
          ! Computing saturation when there is no solid (To be corrected!)
          else
             call AverageField(unitary_field,loc_av_field,.false.,.false.,selected_phase(i), bounds)
             call MPI_ALLREDUCE( loc_av_field,  av_field1, arr_div(1)*arr_div(2)*arr_div(3), &
                                          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_Domain,req(1), ierr)
             if(rank==0) then
               do n=1,bounds(3)
                 do m=1,bounds(2)
                   do l=1,bounds(1)
                     av_field1(l,m,n) = av_field1(l,m,n)*bounds(1)*bounds(2)*bounds(3)/(XLength*YLength*ZLength)
                   enddo
                 enddo
               enddo
             endif
          endif
          ! Printing results
          if(rank==0) then
             call print_results('SaturaPhase'//int2text(selected_phase(i),2), av_field1, timeout, bounds)
             if(output_vtk(i))then
                call printVTK('SaturaPhase'//int2text(selected_phase(i),2),timeout,bounds,av_field1)
             endif
          endif
        endif  
          
      enddo
    endif
  end subroutine ComputeAverages
!======================================================================================================
  subroutine AverageField(original_field,localaveraged_field,solidFilter,complementary,colourFilter,bounds)
    use module_grid
    use module_solid
    use module_VOF
    use module_flow
    implicit none
    include 'mpif.h'
    real(8), dimension(:,:,:), allocatable, intent(in) :: original_field
    real(8), dimension(:,:,:), allocatable :: localaveraged_field
    real(8) :: vol, field, solidwall, colour

    integer :: colourFilter
    integer :: l, m, n, i, j, k, ierr, gap(3), bounds(3)
    integer :: cubedim(3)
    integer :: coreis, corejs, coreks, coreie, coreje, coreke
    integer :: isf, ief, jsf, jef, ksf, kef
    logical :: out_of_range, solidFilter, complementary
    !if (rank==0 .and. debug) write(*,*) 'Averages 1'
    if(colourFilter /= 0 .and. colourFilter /= 1 .and. colourFilter /= 2) then
      call parisError('Colour filter option not valid')
    endif
    !if (rank==0 .and. debug) write(*,*) 'Averages 2'
    !if((.not.dosolids) .and. solidFilter) then
    !   solidFilter = .false.
    !endif
    !if (rank==0 .and. debug) write(*,*) 'Averages 3'
    cubedim(1) = Nx / av_div(1)
    cubedim(2) = Ny / av_div(2)
    cubedim(3) = Nz / av_div(3)
    !if (rank==0 .and. debug) write(*,*) 'Averages 4'
    if(bounds(1) == av_div(1)) then
      gap(1) = 0    
    else 
      gap(1) = cubedim(1) / 2 - cubedim(1)
    endif
    if(bounds(2) == av_div(2)) then
      gap(2) = 0    
    else 
      gap(2) = cubedim(2) / 2 - cubedim(2)
    endif
    if(bounds(3) == av_div(3)) then
      gap(3) = 0    
    else 
      gap(3) = cubedim(3) / 2 - cubedim(3)
    endif
  
    !if (rank==0 .and. debug) write(*,*) 'Averages 5'

    coreis = coords(1)*Nx/nPx + 1
    coreie = (coords(1) + 1)*Nx/nPx

    corejs = coords(2)*Ny/nPy + 1
    coreje = (coords(2) + 1)*Ny/nPy

    coreks = coords(3)*Nz/nPz + 1
    coreke = (coords(3) + 1)*Nz/nPz
    !if (rank==0 .and. debug) write(*,*) 'Averages 6'
    localaveraged_field(1:arr_div(1),1:arr_div(2),1:arr_div(3)) = 0d0
    !if (rank==0 .and. debug) write(*,*) 'Averages 7'
    do n=1,bounds(3)

      do m=1,bounds(2)
  
        do l=1,bounds(1)
          if(l==1 .and. (bounds(1) /= av_div(1))) then
            isf =1; ief=gap(1) + cubedim(1)
          else
            isf = (l-1)*cubedim(1)+gap(1)+1
            ief = l*cubedim(1)+gap(1)
          endif
          if(m==1 .and. (bounds(2) /= av_div(2))) then
            jsf =1; jef=gap(2) + cubedim(2)
          else
            jsf = (m-1)*cubedim(2)+gap(2)+1
            jef = m*cubedim(2)+gap(2)
          endif
          if(n==1 .and. (bounds(3) /= av_div(3))) then
            ksf =1; kef=gap(3) + cubedim(3)
          else
            ksf = (n-1)*cubedim(3)+gap(3)+1
            kef = n*cubedim(3)+gap(3)
          endif

          !if (rank==0) write(*,*) 'First attempt of bounds xyz',isf, ief, jsf, jef, ksf, kef 
          out_of_range=.false.

          if( isf < coreis) isf = coreis
          if( ief <= coreis) out_of_range = .true.
          if( isf >= coreie) out_of_range = .true.
          if( ief > coreie) ief = coreie

          if( jsf < corejs) jsf = corejs
          if( jef <= corejs) out_of_range = .true.
          if( jsf >= coreje) out_of_range = .true.
          if( jef > coreje) jef = coreje

          if( ksf < coreks) ksf = coreks
          if( kef <= coreks) out_of_range = .true.
          if( ksf >= coreke) out_of_range = .true.
          if( kef > coreke) kef = coreke
 
          !if (rank==0) write(*,*) 'Correction of bounds xyz',isf, ief, jsf, jef, ksf, kef 
          isf = isf - coreis + is
          ief = ief - coreis + is
          jsf = jsf - corejs + js
          jef = jef - corejs + js
          ksf = ksf - coreks + ks 
          kef = kef - coreks + ks

          !if (rank==0) write(*,*) 'Adaptation of bounds xyz',isf, ief, jsf, jef, ksf, kef 
          if(out_of_range .eqv. .false.) then
            !if (rank==0) write(*,*) 'Bounds acceptied'
            !write(*,*) 'Rank ', rank,' operating in limits x,y,z ',isf, ief, jsf, jef, ksf, kef 
            do i=isf,ief
              do j=jsf,jef
                do k=ksf,kef 
                !if (rank==0 .and. debug) write(*,*) 'Averages 8'
                  vol = dx(i)*dy(j)*dz(k)
                  if(colourFilter == 0) then
                    colour = 1d0
                  elseif(colourFilter == 1) then
                    colour = cvof(i,j,k)
                  elseif(colourFilter == 2) then  
                    colour = 1d0 - cvof(i,j,k)
                  endif
                  !if (rank==0 .and. debug) write(*,*) 'Averages 9'
                  if(solidFilter) then
                    solidwall = 1d0 - solids(i,j,k)
                  else
                    solidwall = 1d0
                  endif
                  !if (rank==0 .and. debug) write(*,*) 'Averages 10'
                  if(complementary) then
                    field = 1d0 - original_field(i,j,k)
                  else 
                    field = original_field(i,j,k)
                  endif
                  !if (rank==0 .and. debug) write(*,*) 'Averages 11'
                  localaveraged_field(l,m,n) = localaveraged_field(l,m,n) + field*solidwall*vol*colour

                enddo
              enddo
            enddo

          endif

        enddo
      enddo
    enddo
    
  end subroutine AverageField
!======================================================================================================
  subroutine ReadAveParameters
    !use module_solid
    use module_grid
    use module_flow
    !use module_vof
    implicit none
    include 'mpif.h'
    integer :: in, ierr, i
    logical :: file_is_there


    namelist /aveparameters/ averages_to_do, av_on_nodes, av_div, selected_phase, output_vtk
    
    averages_to_do=['undefined','undefined','undefined','undefined','undefined','undefined','undefined']
    av_on_nodes=['none', 'none', 'none', 'none', 'none', 'none', 'none']
    av_div(1:3)=1
    selected_phase(1:SIZEARR)=0
    output_vtk(1:SIZEARR)=.false.
    
    in=31

    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    inquire(file='inputaverages',exist=file_is_there)
    open(unit=in, file='inputaverages', status='old', action='read', iostat=ierr)
    if (file_is_there) then
       if(ierr == 0) then
          read(UNIT=in,NML=aveparameters)
          if(rank==0) write(*,*)'Average parameters read successfully'
       else
          print *, 'rank=',rank,' has error ',ierr,' opening file inputaverages'
       endif
    else
       if (rank == 0) then
         write(*,*) "ReadAveParameters: no 'inputaverages' file."
         do_averages= .false.
       endif
    endif
    close(in)

    !if(.not.dosolids) then 
    !   if (rank == 0) then
    !      write(*,*) "ReadAveParameters: cannot work without solids."
    !      do_averages= .false.
    !   endif
    !endif

    if(mod(Nx,av_div(1)) /= 0) call pariserror("Number of nodes in X not divisible by requested average units in X")
    if(mod(Ny,av_div(2)) /= 0) call pariserror("Number of nodes in Y not divisible by requested average units in Y")
    if(mod(Nz,av_div(3)) /= 0) call pariserror("Number of nodes in Z not divisible by requested average units in Z")
 
    arr_div(1) = av_div(1) + 1
    arr_div(2) = av_div(2) + 1
    arr_div(3) = av_div(3) + 1
    
    do i=1,SIZEARR
      if((DoVOF .eqv. .false.)  .and. selected_phase(i) /= 0) call pariserror("Selected Phase not available")
      averages_on_nodes(i,1)=.false.;averages_on_nodes(i,2)=.false.;averages_on_nodes(i,3)=.false.
      if(av_on_nodes(i)=='x') then
        averages_on_nodes(i,1)=.true.
      elseif(av_on_nodes(i)=='y') then
        averages_on_nodes(i,2)=.true.
      elseif(av_on_nodes(i)=='z') then
        averages_on_nodes(i,3)=.true.
      elseif(av_on_nodes(i)=='xy') then
        averages_on_nodes(i,1)=.true.;averages_on_nodes(i,2)=.true.
      elseif(av_on_nodes(i)=='xz') then
        averages_on_nodes(i,1)=.true.;averages_on_nodes(i,3)=.true.
      elseif(av_on_nodes(i)=='yz') then
        averages_on_nodes(i,2)=.true.;averages_on_nodes(i,3)=.true.
      elseif(av_on_nodes(i)=='xyz') then
        averages_on_nodes(i,1)=.true.;averages_on_nodes(i,2)=.true.;averages_on_nodes(i,3)=.true.
      endif
    enddo

    ! Test
    if( rank==0 .and. debug )then
      write(*,*) 'Written parameters in averages'
      write(*,*) averages_to_do
      write(*,*) averages_on_nodes
      write(*,*) av_div
      write(*,*) selected_phase
      write(*,*) output_vtk
      write(*,*) Nx, Ny, Nz, nPx, nPy, nPz
    endif
  end subroutine ReadAveParameters
!=====================================================================================================
  SUBROUTINE print_results(file_name, vector, iout, bounds)
  use module_IO
  implicit none

  character(LEN=13) :: file_name
  integer :: x_div, iout, bounds(3)
  integer :: l, m, n
  real(8), dimension(:,:,:), allocatable :: vector
  integer, save :: unit_num = 89
  logical :: nodes
  

    OPEN(UNIT=unit_num,FILE=TRIM(out_path)//'/'//trim(file_name)//'-'//TRIM(int2text(iout,5))//'.dat') 

      do n=1,bounds(3); do m=1,bounds(2); do l=1,bounds(1);  
        write(unit_num,11) vector(l, m, n)
      enddo; enddo; enddo
   
    10  format(I5.5, E14.6)
    11  format(E14.6)
    close(unit=unit_num)

  END SUBROUTINE print_results

subroutine printVTK(filename,iout,bounds,vector)
  use module_flow
  use module_grid
  !use module_hello
  use module_IO
  !use IO_mod
  implicit none
  integer :: iout,bounds(3), l, m, n
  real(8) :: origin(3)
  real(8), dimension(:,:,:), allocatable :: vector
  character(len=13) :: rootname,filename
  rootname=TRIM(out_path)//'/VTK/'//TRIM(filename)//'-'//TRIM(int2text(iout,5))//'.vtk'

  !write(*,*) rootname
  if(bounds(1)==av_div(1)) then
    origin(1) = XLength / (2 * REAL(av_div(1)))
  else
    origin(1) = 0d0
  endif

  if(bounds(2)==av_div(2)) then
    origin(2) = YLength / (2 * REAL(av_div(2)))
  else
    origin(2) = 0d0
  endif

  if(bounds(3)==av_div(3)) then
    origin(3) = ZLength / (2 * REAL(av_div(3)))
  else
    origin(3) = 0d0
  endif

  OPEN(UNIT=89,FILE=TRIM(out_path)//'/VTK/'//TRIM(filename)//'-'//TRIM(int2text(iout,5))//'.vtk')
    write(89,10)
    write(89,11) time
    write(89,12)
    write(89,13)
    write(89,14) bounds(1),bounds(2),bounds(3)
    write(89,15) origin(1),origin(2),origin(3)
    write(89,16) XLength / REAL(av_div(1)), YLength / REAL(av_div(2)), ZLength / REAL(av_div(3))
10  format('# vtk DataFile Version 2.0')
11  format('grid, time ',F16.8)
12  format('ASCII')
13  format('DATASET STRUCTURED_POINTS')
14  format('DIMENSIONS ',I5,I5,I5)
15  format('ORIGIN ',F16.8,F16.8,F16.8)
16  format('SPACING ',F16.8,F16.8,F16.8)

    write(89,19) bounds(1)*bounds(2)*bounds(3)
    write(89,17) TRIM(filename)
    write(89,18)

19  format('POINT_DATA ',I17)
17  format('SCALARS ',A20,' double 1')
18  format('LOOKUP_TABLE default')

    do n=1,bounds(3); do m=1,bounds(2); do l=1,bounds(1); 
      write(89,210) vector(l,m,n)
    enddo; enddo; enddo
210 format(e14.5)


    close(89)

end subroutine printVTK
!=========================================================================
subroutine outputmessage(message) 
  use module_IO
  use module_grid
  implicit none
  include 'mpif.h'
  integer ierr
  character(*) :: message

  OPEN(UNIT=89,FILE=TRIM(out_path)//'/message-rank-'//TRIM(int2text(rank,padding))//'.txt')
  write(89,*) message
  close(89)
  
end subroutine outputmessage

end module module_averages
