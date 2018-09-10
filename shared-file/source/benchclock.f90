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

module benchclock

  implicit none

  logical,          save, private :: firstcall = .true.
  double precision, save, private :: ticktime = 0.0

  integer, parameter :: int32kind = selected_int_kind( 9)
  integer, parameter :: int64kind = selected_int_kind(18)

!
!  Select high resolution clock
!

  integer, parameter :: intkind = int64kind
  integer(kind = intkind) :: count,rate

contains

double precision function benchtime()

  double precision :: dummy

! Ensure clock is initialised  

  if (firstcall) dummy = benchtick()

  call system_clock(count)

  benchtime  = dble(count)*ticktime

end function benchtime


double precision function benchtick()

  if (firstcall) then

     firstcall = .false.
     call system_clock(count, rate)
     ticktime = 1.0d0/dble(rate)

  end if

  benchtick = ticktime

end function benchtick

end module benchclock



