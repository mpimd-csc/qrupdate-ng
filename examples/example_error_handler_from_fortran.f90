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
program example_error_handler_from_fortran
    use qrupdate_error
    use qrupdate
    use iso_fortran_env
    implicit none

    complex(real32) :: A(1,1), u(1)
    real(real32) :: rw(1)
    integer :: info
    integer :: lda
    integer :: n
    procedure(error_handler_if), pointer  :: current_handler


    n = -1
    lda = 1
    info = 0

    call qrupdate_set_error(my_error_handler)

    call cch1dn(n, A, lda, u, rw, info)

    current_handler => qrupdate_get_error()

    write(*,*) "reset the handler"

    call qrupdate_set_error()
    call cch1dn(n, A, lda, u, rw, info)

    write(*,*) "restore the saved handler"
    call qrupdate_set_error(current_handler)
    call cch1dn(n, A, lda, u, rw, info)

    stop

contains
    subroutine my_error_handler(srname, linfo, aux)
        character(len=*), intent(in) :: srname
        integer, intent(in) :: linfo
        class(*), optional, intent(in) :: aux

        write(*,*) "My error handler"
        write(*,*) "srname = ", srname
        write(*,*) "info = ", linfo
    end subroutine
end program
