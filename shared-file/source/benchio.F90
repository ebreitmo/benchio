! **************************************************************************
!                                                                          !
!    EPCC I/O benchmarking application.                                    !
!                                                                          !
!    Copyright (c) 2018, The University of Edinburgh                       !
!                            All rights reserved.                          !
!                     List of contributors in AUTHORS                      !
!                                                                          !
!    This program is free software: you can redistribute it and/or modify  !
!    it under the terms of the GNU General Public License as published by  !
!    the Free Software Foundation, either version 3 of the License, or     !
!    (at your option) any later version.                                   !
!                                                                          !
!    This program is distributed in the hope that it will be useful,       !
!    but WITHOUT ANY WARRANTY; without even the implied warranty of        !
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         !
!    GNU General Public License for more details.                          !
!                                                                          !
!    You should have received a copy of the GNU General Public License     !
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.!
!                                                                          !
! **************************************************************************

program benchio

  use benchclock
  use mpi
  use mpiio
  use iohdf5
  use ionetcdf

  implicit none

  integer, parameter :: numiolayer = 4
  integer, parameter :: maxlen = 64
  integer :: numrep
! To check preprocessor flag
  integer :: withserial,withmpiio,withhdf5,withnetcdf

  character*(maxlen), dimension(numiolayer)  :: iostring, iolayername

  character*(maxlen) :: filedir, filename

  integer :: iolayer, irep

! Set local array size - global sizes l1, l2 and l3 are scaled
! by number of processes in each dimension

  integer :: n1
  integer :: n2
  integer :: n3
  integer, parameter :: ndim = 3

  integer :: i1, i2, i3, j1, j2, j3, l1, l2, l3, p1, p2, p3

  double precision, allocatable :: iodata(:,:,:)

  integer :: rank, size, ierr, comm, cartcomm, dblesize
  integer, dimension(ndim) :: dims, coords

  integer, parameter :: iounit = 12
  integer, parameter :: mib = 1024*1024

  logical :: reorder = .false.
  logical, dimension(ndim) :: periods = [.false., .false., .false.]

  double precision :: t0, t1, time, iorate, mibdata
  double precision :: mintime, maxiorate, avgtime, avgiorate

! Default values, can be overwritten later by input file 'BenchIO-Input.txt'
 n1 = 256
 n2 = 256
 n3 = 256
 numrep=10
 filedir = 'benchio_files'
 ! Set default to false
  withserial = 0
  withmpiio = 0
  withhdf5 = 0
  withnetcdf = 0
 
  iostring(1) = 'Serial'
  iostring(2) = 'MPI-IO'
  iostring(3) = ' HDF5 '
  iostring(4) = 'NetCDF'

  iolayername(1) = 'serial.dat'
  iolayername(2) = 'mpiio.dat'
  iolayername(3) = 'hdf5.dat'
  iolayername(4) = 'netcdf.dat'

  call MPI_Init(ierr)

  comm = MPI_COMM_WORLD

  call MPI_Comm_size(comm, size, ierr)
  call MPI_Comm_rank(comm, rank, ierr)

! Read inpit file 'BenchIO-Input.txt' which contains all input parameters
 if (rank == 0) then
   open(50, FILE='BenchIO-Input.txt',FORM='FORMATTED',STATUS='old')
   read(50,*)n1          ! 3D grid
   read(50,*)n2
   read(50,*)n3
   read(50,*)numrep      ! Number of iterations
   read(50,*)filedir     ! Where to put all the files to be created
   read(50,*)withserial  ! Is serial benchmark to be executed: 0=no, 1=yes 
   read(50,*)withmpiio   ! Is mpiio benchmark to be executed
   read(50,*)withhdf5    ! IS hdf5 benchmark to be executed
   read(50,*)withnetcdf  ! IS netcdf benchmark to be executed

#ifndef WITH_SERIAL      ! In case the flag has not been set in the Makefile
#define WITH_SERIAL 0
#endif

#ifndef WITH_MPIIO
#define WITH_MPIIO 0
#endif

#ifndef WITH_HDf5
#define WITH_HDF5 0
#endif

#ifndef WITH_NETCDF
#define WITH_NETCDF 0
#endif

! Check if choice for executing the serial/mpiio/hdf5/netcdf benchmarks agrees in the Makefile (with compiler flags)
! and the input file 'BenchIO-Input.txt'
! Prints a warning if they don't agree
! The setting int he input file 'BenchIO-Input.txt' overwrites the setting from the Makefile
   if (withserial.eq.0.and.WITH_SERIAL.eq.1) then
    write(*,*)'Warning: Default variable WITH_SERIAL set differently in Makefile and input file!'
   endif
   if (withserial.eq.1.and.WITH_SERIAL.eq.0) then
    write(*,*)'Warning:	Default	variable WITH_SERIAL set differently in Makefile and input file!'
   endif
   if (withmpiio.eq.0.and.WITH_MPIIO.eq.1) then
    write(*,*)'Warning: Default variable WITH_MPIIO set differently in Makefile and input file!'
   endif
   if (withmpiio.eq.1.and.WITH_MPIIO.eq.0) then
    write(*,*)'Warning: Default variable WITH_MPIIO set differently in Makefile and input file!'
   endif
   if (withhdf5.eq.0.and.WITH_HDF5.eq.1) then
    write(*,*)'Warning: Default variable WITH_HDF5 set differently in Makefile and input file!'
   endif
   if (withhdf5.eq.1.and.WITH_HDF5.eq.0) then
    write(*,*)'Warning: Default variable WITH_HDF5 set differently in Makefile and input file!'
   endif
   if (withnetcdf.eq.0.and.WITH_NETCDF.eq.1) then
    write(*,*)'Warning: Default variable WITH_NETCDF set differently in Makefile and input file!'
   endif
   if (withnetcdf.eq.1.and.WITH_NETCDF.eq.0) then
    write(*,*)'Warning: Default variable WITH_NETCDF set differently in Makefile and input file!'
   endif
   close(50)
 endif
 
 !Broadcast the variables read in from 'BenchIO-Input.txt' to all processors
 call MPI_Bcast(n1,1,MPI_INTEGER,0,comm,ierr)
 call MPI_Bcast(n2,1,MPI_INTEGER,0,comm,ierr)
 call MPI_Bcast(n3,1,MPI_INTEGER,0,comm,ierr)
 call MPI_Bcast(numrep,1,MPI_INTEGER,0,comm,ierr)
 call MPI_Bcast(filedir,maxlen,MPI_CHARACTER,0,comm,ierr)
 call MPI_Bcast(withserial,1,MPI_INTEGER,0,comm,ierr)
 call MPI_Bcast(withmpiio,1,MPI_INTEGER,0,comm,ierr)
 call MPI_Bcast(withhdf5,1,MPI_INTEGER,0,comm,ierr)
 call MPI_Bcast(withnetcdf,1,MPI_INTEGER,0,comm,ierr)

   allocate(iodata(0:n1+1, 0:n2+1, 0:n3+1))

  dims(:) = 0

! Set 3D processor grid

  call MPI_Dims_create(size, ndim, dims, ierr)

! Reverse dimensions as MPI assumes C ordering (this is not essential)

  p1 = dims(3)
  p2 = dims(2)
  p3 = dims(1)

! Compute global sizes

  l1 = p1*n1
  l2 = p2*n2
  l3 = p3*n3

  call MPI_Type_size(MPI_DOUBLE_PRECISION, dblesize, ierr)

  mibdata = float(dblesize*n1*n2*n3)*float(p1*p2*p3)/float(mib)

  if (rank == 0) then
     write(*,*)
     write(*,*) 'Simple Parallel IO benchmark'
     write(*,*) '----------------------------'
     write(*,*)
     write(*,*) 'Running on ', size, ' process(es)'
     write(*,*) 'Process grid is (', p1, ', ', p2, ', ', p3, ')'
     write(*,*) 'Array size is   (', n1, ', ', n2, ', ', n3, ')'
     write(*,*) 'Global size is  (', l1, ', ', l2, ', ', l3, ')'
     write(*,*)
     write(*,*) 'Total amount of data = ', mibdata, ' MiB'
     write(*,*)
     write(*,*) 'Clock resolution is ', benchtick()*1.0e6, ', usecs'
  end if
  
  dims(1) = p1
  dims(2) = p2
  dims(3) = p3

  call MPI_Cart_create(comm, ndim, dims, periods, reorder, cartcomm, ierr)

! Set halos to illegal values

  iodata(:,:,:) = -1
  
! Set iodata core to have unique values 1, 2, ..., p1*n1*p2*n2*p3*n3

  call MPI_Cart_coords(cartcomm, rank, ndim, coords, ierr)
  
  do i3 = 1, n3
     do i2 = 1, n2
        do i1 = 1, n1

           j1 = coords(1)*n1 + i1
           j2 = coords(2)*n2 + i2
           j3 = coords(3)*n3 + i3

           iodata(i1,i2,i3) = (j3-1)*l1*l2 + (j2-1)*l1 + j1

        end do
     end do
  end do

  do iolayer = 1, numiolayer

!  Skip layer if support is not compiled in
!  Expects iolayers in order: serial, MPI-IO, HDF5, NetCDF
if (iolayer == 1) then
 if (withserial.ne.1) then
       cycle
     endif
endif

if (iolayer == 2) then
 if (withmpiio.ne.1) then
       cycle
     endif
endif

if (iolayer == 3) then
 if (withhdf5.ne.1) then
       cycle
     endif
endif

if (iolayer == 4) then
 if (withnetcdf.ne.1) then
        cycle
     endif
endif

     if (rank == 0) then
        write(*,*)
        write(*,*) '------'
        write(*,*) iostring(iolayer)
        write(*,*) '------'
        write(*,*)
     end if

     filename = trim(filedir)//'/'//trim(iolayername(iolayer))

     if (rank == 0) then
        write(*,*) 'Writing to ', filename
        mintime = 0
        maxiorate = 0
        avgtime = 0
        avgiorate = 0
     end if

     do irep = 1, numrep
       call MPI_Barrier(comm, ierr)
       if (rank == 0) then
         t0 = benchtime()
       end if

       select case (iolayer)

       case(1)
          call serialwrite(filename, iodata, n1, n2, n3, cartcomm)

       case(2)
          call mpiiowrite(filename, iodata, n1, n2, n3, cartcomm)

       case(3)
          call hdf5write(filename, iodata, n1, n2, n3, cartcomm)

       case(4)
          call netcdfwrite(filename, iodata, n1, n2, n3, cartcomm)

       case default
          write(*,*) 'Illegal value of iolayer = ', iolayer
          stop

       end select

       call MPI_Barrier(comm, ierr)
       if (rank == 0) then
          t1 = benchtime()

          time = t1 - t0
          iorate = mibdata/time
          avgtime = avgtime + time/numrep
          avgiorate = avgiorate + iorate/numrep

          if (maxiorate < iorate) then
            maxiorate = iorate
            mintime = time
          end if

          write(*,*) 'time = ', time, ', rate = ', iorate, ' MiB/s'
          call fdelete(filename)
       end if
     end do
     if (rank == 0) then
       write(*,*) 'mintime = ', mintime, ', maxrate = ', maxiorate, ' MiB/s'
       write(*,*) 'avgtime = ', avgtime, ', avgrate = ', avgiorate, ' MiB/s'
       write(*,*) 'Deleting: ', filename
       write(*,*)
     end if
  end do

  if (rank == 0) then
     write(*,*)
     write(*,*) '--------'
     write(*,*) 'Finished'
     write(*,*) '--------'
     write(*,*)
  end if

! Clean up allocated storage
  deallocate(iodata)

  call MPI_Finalize(ierr)
  
end program benchio

subroutine fdelete(filename)

  implicit none

  character *(*) :: filename
  integer, parameter :: iounit = 15
  integer :: stat

  open(unit=iounit, iostat=stat, file=filename, status='old')
  if (stat.eq.0) close(unit=iounit, status='delete')

end subroutine fdelete
