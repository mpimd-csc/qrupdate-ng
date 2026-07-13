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
!> \brief Module for custom error handling.
!>
!> This module provides a mechanism to override the default BLAS/LAPACK
!> error handler (xerbla). Users can specify a custom subroutine to
!> handle errors, along with optional auxiliary data.
module qrupdate_error
    use qrupdate_blas
    implicit none

    !> \brief Abstract interface for custom error handlers.
    !> \ingroup error
    abstract interface
    subroutine error_handler_if(srname, info, aux)
        character(len=*), intent(in) :: srname
        integer, intent(in) :: info
    class(*), optional, intent(in) :: aux
end subroutine error_handler_if
    end interface

    procedure(error_handler_if), pointer :: global_error_handler => null()
class(*), pointer :: global_error_aux => null()

private :: global_error_handler, global_error_aux
contains

    !> \brief Sets the custom error handler.
    !>
    !> \param[in] p_handler
    !> \verbatim
    !>          p_handler is a procedure pointer conforming to the
    !>          error_handler_if interface.
    !> \endverbatim
    !> \ingroup error
    subroutine qrupdate_set_error(p_handler)
        procedure(error_handler_if), pointer, intent(in) :: p_handler
        global_error_handler => p_handler
    end subroutine qrupdate_set_error

    !> \brief Sets auxiliary data for the error handler.
    !>
    !> \param[in] p_aux
    !> \verbatim
    !>          p_aux is a pointer to a polymorphic object (class(*))
    !>          that will be passed to the error handler.
    !> \endverbatim
    !> \ingroup error
    subroutine qrupdate_set_error_data(p_aux)
    class(*), pointer :: p_aux
    global_error_aux => p_aux
end subroutine qrupdate_set_error_data

!> \brief Dispatches error reporting to the handler.
!>
!> This subroutine checks if a custom error handler has been set via
!> qrupdate_set_error. If so, it calls that handler with the provided
!> routine name, error code, and any set auxiliary data. Otherwise, it
!> falls back to the standard LAPACK xerbla routine.
!>
!> \param[in] srname
!> \verbatim
!>          srname is CHARACTER(LEN=*)
!>          The name of the routine that encountered the error.
!> \endverbatim
!> \param[in] info
!> \verbatim
!>          info is INTEGER
!>          The error code.
!> \endverbatim
!> \ingroup error
subroutine qrupdate_xerror(srname, info)
    character(len=*), intent(in) :: srname
    integer, intent(in) :: info
    if (associated(global_error_handler)) then
        call global_error_handler(srname, info, global_error_aux)
    else
        call xerbla(srname, info)
    end if
end subroutine qrupdate_xerror

end module qrupdate_error
