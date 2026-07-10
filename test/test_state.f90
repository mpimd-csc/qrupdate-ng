! Copyright (C) 2026 Martin Köhler <koehlerm(AT)mpi-magdeburg.mpg.de>
!
! This file is part of qrupdate-ng.
!
! qrupdate is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this software; see the file COPYING.  If not, see
! <http://www.gnu.org/licenses/>.
!
module test_state
    implicit none
    logical :: xerbla_called = .false.
    logical :: custom_handler_called = .false.
    character(len=64) :: last_srname = ''
    integer :: last_info = 0
    integer :: call_count = 0
contains
    subroutine reset()
        xerbla_called = .false.
        custom_handler_called = .false.
        last_srname = ''
        last_info = 0
        call_count = 0
    end subroutine reset
end module test_state
