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
!> Since qrupdate 1.2.0, the way of handling errors changed a bit. Instead of
!> calling xerbla from LAPACK/BLAS, a user defined handler can registered. If this
!> is not done, xerbla is called as before. See \ref qrupdate_error::qrupdate_set_error for
!> details.
!>
!> ## Fortran error handling
!>
!> To set a custom error handler in Fortran, create a subroutine that matches
!> the \c error_handler_if interface and pass it to \c qrupdate_error::qrupdate_set_error:
!>
!> \code{.f90}
!> use qrupdate_error
!> implicit none
!>
!> subroutine my_error_handler(srname, info, aux)
!>     character(len=*), intent(in) :: srname
!>     integer, intent(in) :: info
!>     class(*), optional, intent(in) :: aux
!>
!>     print *, 'Error in routine:', srname
!>     print *, 'Error code:', info
!>     if (present(aux)) then
!>         print *, 'Auxiliary data available'
!>     end if
!> end subroutine my_error_handler
!>
!> call qrupdate_set_error(my_error_handler)
!> \endcode
!>
!> ## C error handling
!>
!> To use a custom error handler from C, implement a function with the following
!> signature and set it using \c qrupdate_set_error:
!>
!> \code{.c}
!> #include "qrupdate.h"
!>
!> void my_error_handler(const char *srname, int info, void *aux)
!> {
!>     printf("Error in routine %s (code: %d)\n", srname, info);
!>     if (aux != NULL) {
!>         printf("Auxiliary data: %p\n", aux);
!>     }
!> }
!>
!> int main() {
!>     qrupdate_set_error(my_error_handler);
!>     // ... call qrupdate functions
!>     return 0;
!> }
!> \endcode
!>
!> The C error handler receives:
!> - \p srname: The routine name as a character array (without null terminator)
!> - \p info: Pointer to the error code (integer)
!> - \p aux: Pointer to auxiliary data (NULL if not set)
!> - \p srname_len: Length of the routine name string
!>
!> To pass auxiliary data to the C error handler, use \c qrupdate_set_error_data_c:
!>
!> \code{.c}
!> int my_data = 42;
!> qrupdate_set_error_data_c(&my_data);
!> \endcode
module qrupdate_error
    use iso_c_binding
    implicit none

    !> \brief Abstract interface for Fortran custom error handlers.
    !> \ingroup error
    !>
    !> This interface defines the required signature for Fortran error handler
    !> subroutines passed to qrupdate_set_error. The handler is called when
    !> qrupdate_xerror is invoked with the routine name, error code, and
    !> optional auxiliary data.
    !>
    !> The handler must be a subroutine with the following signature:
    !>
    !> \code{.f90}
    !> subroutine my_handler(srname, info, aux)
    !>     character(len=*), intent(in) :: srname
    !>     integer, intent(in) :: info
    !>     class(*), optional, intent(in) :: aux
    !>     ! Implementation
    !> end subroutine
    !> \endcode
    !>
    !> \param[in] srname
    !> \verbatim
    !>          The name of the routine that encountered the error.
    !> \endverbatim
    !> \param[in] info
    !> \verbatim
    !>          The error code. Positive values indicate errors.
    !> \endverbatim
    !> \param[in] aux
    !> \verbatim
    !>          Optional pointer to polymorphic auxiliary data.
    !>          Use present(aux) to check if data was provided.
    !> \endverbatim
    abstract interface
        subroutine error_handler_if(srname, info, aux)
            character(len=*), intent(in) :: srname
            integer, intent(in) :: info
            class(*), optional, intent(in) :: aux
        end subroutine error_handler_if
    end interface

    !> \brief Abstract interface for C custom error handlers.
    !> \ingroup error
    !>
    !> This interface defines the required signature for C error handler
    !> functions when used through the C-binding interface qrupdate_set_error_c.
    !> The handler receives the routine name as a character array, error code
    !> as a pointer, and optional auxiliary data as a C pointer.
    !>
    !> The C function must have the following signature:
    !>
    !> \code{.c}
    !> void my_handler(const char *srname, int info, void *aux )
    !> {
    !>     // Implementation
    !> }
    !> \endcode
    !>
    !> \param[in] srname
    !> \verbatim
    !>          The name of the routine that encountered the error.
    !>          Passed as a Fortran-style character array (no null terminator).
    !>          The length is provided in \p srname_len.
    !> \endverbatim
    !> \param[in] info
    !> \verbatim
    !>          The error code (default C integer).
    !> \endverbatim
    !> \param[in] aux
    !> \verbatim
    !>          Optional C pointer to auxiliary data.
    !>          NULL if no data was set via qrupdate_set_error_data_c.
    !> \endverbatim
    abstract interface
        subroutine error_handler_if_c(srname, info, aux)
            import c_ptr, c_char, c_int
            character(kind=c_char), intent(in) :: srname(*)
            integer(kind = c_int), value, intent(in) :: info
            type(c_ptr), value, intent(in) :: aux
        end subroutine error_handler_if_c
    end interface

    !> \brief Set or reset the custom error handler.
    !> \ingroup error
    !>
    !> This generic interface provides two ways to manage the error handler:
    !>
    !> - \ref qrupdate_error::qrupdate_set_error_internal - Register a custom Fortran error handler
    !> - \ref qrupdate_error::qrupdate_reset_error - Reset to the default LAPACK xerbla behavior
    !>
    !> When a custom error handler is set via \ref qrupdate_set_error_internal,
    !> it will be called instead of xerbla when errors occur. Use
    !> \ref qrupdate_reset_error to restore the default behavior.
    !>
    !> \note The C interface uses \ref qrupdate_set_error_c instead.
    !>
    !> \code{.f90}
    !> use qrupdate_error
    !> implicit none
    !>
    !> call qrupdate_set_error(my_error_handler)     ! Set custom handler
    !> call qrupdate_set_error()                      ! Reset to default
    !> \endcode
    !>
    !> \sa qrupdate_error::qrupdate_set_error_internal
    !> \sa qrupdate_error::qrupdate_reset_error
    interface qrupdate_set_error
        procedure :: qrupdate_set_error_internal, qrupdate_reset_error
    end interface


    procedure(error_handler_if), pointer :: global_error_handler => null()
    procedure(error_handler_if_c), pointer :: global_error_handler_c => null()
    class(*), pointer :: global_error_aux => null()
    type(c_ptr) :: global_error_aux_c = c_null_ptr
    logical :: set_from_c = .false.

    private :: global_error_handler, global_error_handler_c , global_error_aux, global_error_aux_c, set_from_c
contains

    !> \brief Sets the custom error handler (Fortran interface).
    !>
    !> This subroutine registers a custom error handler that will be called
    !> instead of the default LAPACK xerbla routine when an error occurs.
    !>
    !> The error handler must conform to the \c error_handler_if interface,
    !> which requires the following signature:
    !>
    !> \code{.f90}
    !> subroutine my_error_handler(srname, info, aux)
    !>     use qrupdate_error
    !>     implicit none
    !>     character(len=*), intent(in) :: srname
    !>     integer, intent(in) :: info
    !>     class(*), optional, intent(in) :: aux
    !>     ! Handler implementation
    !> end subroutine my_error_handler
    !> \endcode
    !>
    !> where:
    !> - \p srname is the name of the routine that encountered the error
    !> - \p info is the error code (positive value indicates an error)
    !> - \p aux is optional auxiliary data passed via qrupdate_set_error_data
    !>
    !> To pass auxiliary data to the error handler:
    !>
    !> \code{.f90}
    !> type(my_data_type), target :: my_data
    !> call qrupdate_set_error_data(my_data)
    !> call qrupdate_set_error(my_error_handler)
    !> \endcode
    !>
    !> \param[in] p_handler
    !> \verbatim
    !>          p_handler is a procedure pointer conforming to the
    !>          error_handler_if interface.
    !> \endverbatim
    !> \ingroup error
    !> \sa qrupdate_error::qrupdate_set_error
    subroutine qrupdate_set_error_internal(p_handler)
        procedure(error_handler_if) :: p_handler
        set_from_c = .false.
        global_error_handler => p_handler
    end subroutine qrupdate_set_error_internal

    !> \brief Reset the custom error handler to the default
    !>
    !> The subroutine resets the custom error handler to the default
    !> behavior calling XERBLA. If the error handler was set from C
    !> inbetween, the reset is done as well.
    !>
    !> A more flexible call to this subroutine is done using the
    !> \ref qrupdate_error::qrupdate_set_error interface.
    !>
    !> \ingroup error
    subroutine qrupdate_reset_error()
        set_from_c = .false.
        global_error_handler => null()
    end subroutine

    !> \brief Return the current error handler
    !> \ingroup error
    !>
    !> Retrieves the currently registered Fortran error handler.
    !>
    !> \returns p_handler
    !> \verbatim
    !>          A procedure pointer to the current error handler subroutine,
    !>          or NULL if no custom handler is set (default behavior uses xerbla).
    !> \endverbatim
    !>
    !> Since the function returns a pointer, it must be used in the following way:
    !>
    !> \code{.f90}
    !> procedure(error_handler_if), pointer :: current_handler
    !>
    !> current_hanlder => qrupdate_get_error()
    !> \endcode
    !>
    function qrupdate_get_error() result(p_handler)
        procedure(error_handler_if), pointer  :: p_handler

        p_handler => global_error_handler
        return
    end function

    !> \brief Sets the custom error handler (C interface).
    !>
    !> This subroutine is the C-binding wrapper for setting an error handler.
    !> It should be called with a function pointer to a C function with the
    !> following signature:
    !>
    !> \code{.c}
    !> void c_error_handler(const char *srname, int info, void *aux)
    !> \endcode
    !>
    !> The C error handler receives:
    !> - \p srname: The name of the routine that encountered the error (char array),
    !>              in Fortran style without trailing '\0'. The caller must extract
    !>              the string length from \p srname_len.
    !> - \p info: Pointer to the error code (integer pointer)
    !> - \p aux: Optional auxiliary data passed via qrupdate_set_error_data_c
    !>           (void pointer, NULL if not set)
    !>
    !> The error handler is called when qrupdate_xerror is invoked, typically
    !> when a LAPACK/BLAS routine detects an error. If no custom handler is
    !> set, the standard xerbla routine is called. If the function is called
    !> with a NULL argument, the error handler is set back to XERBLA again.
    !>
    !> \code{.c}
    !> #include "qrupdate.h"
    !>
    !> void my_error_handler(const char *srname, int info, void *aux)
    !> {
    !>     printf("Error in %s: code %d\n", srname, info);
    !> }
    !>
    !> int main() {
    !>     qrupdate_set_error(my_error_handler);
    !>     // ... call qrupdate functions
    !>     qrupdate_set_error(NULL);
    !>     // Now xerbla is the error handler.
    !>     return 0;
    !> }
    !> \endcode
    !>
    !> \param[in] p_handler
    !> \verbatim
    !>          p_handler is a C function pointer to the error handler. If the pointer is NULL, the
    !>          behavior is reseted to the default.
    !> \endverbatim
    !> \ingroup error
    subroutine qrupdate_set_error_c(p_handler) bind(C, name = "qrupdate_set_error")
        use iso_c_binding
        type(c_funptr), intent(in), value :: p_handler
        if ( .not. c_associated(p_handler)) then
            set_from_c = .false.
            global_error_handler_c => null()
        else
            set_from_c = .true.
            call c_f_procpointer(p_handler, global_error_handler_c)
        end if
    end subroutine

    !> \brief Sets auxiliary data for the error handler (Fortran interface).
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

    !> \brief Sets auxiliary data for the C error handler.
    !>
    !> This subroutine is the C-binding wrapper for setting auxiliary data
    !> that will be passed to the C error handler.
    !>
    !> The auxiliary data is passed as a C pointer (void*) to the error handler.
    !>
    !> \param[in] p_aux
    !> \verbatim
    !>          p_aux is a C pointer (type(c_ptr)) to the auxiliary data.
    !> \endverbatim
    !> \ingroup error
    subroutine qrupdate_set_error_data_c(p_aux) bind(C, name="qrupdate_set_error_data")
        type(c_ptr), intent(in), value :: p_aux

        global_error_aux_c = p_aux
    end subroutine

    !> \brief Return the current C error handler
    !> \ingroup error
    !>
    !> Retrieves the currently registered C error handler.
    !>
    !> \returns p_handler
    !> \verbatim
    !>          A C function pointer to the current error handler function,
    !>          or NULL if no custom handler is set (default behavior uses xerbla).
    !> \endverbatim
    function qrupdate_get_error_c() bind(C, name = "qrupdate_get_error") result(p_handler)
        type(c_funptr) :: p_handler
        p_handler = c_funloc(global_error_handler_c)
    end function


    !> \brief Dispatches error reporting to the handler.
    !>
    !> This subroutine checks if a custom error handler has been set via
    !> qrupdate_set_error or qrupdate_set_error_c. If so, it calls that handler
    !> with the provided routine name, error code, and any set auxiliary data.
    !> Otherwise, it falls back to the standard LAPACK xerbla routine.
    !>
    !> The error handler is typically invoked automatically when a qrupdate
    !> routine encounters an error. The dispatch logic is:
    !>
    !> 1. If a Fortran error handler is set (via qrupdate_set_error), it is called.
    !> 2. If a C error handler is set (via qrupdate_set_error_c), it is called.
    !> 3. Otherwise, the standard xerbla routine is called.
    !>
    !> \param[in] srname
    !> \verbatim
    !>          srname is CHARACTER(LEN=*)
    !>          The name of the routine that encountered the error.
    !> \endverbatim
    !> \param[in] info
    !> \verbatim
    !>          info is INTEGER
    !>          The error code (positive value indicates an error).
    !> \endverbatim
    !> \ingroup error
    subroutine qrupdate_xerror(srname, info)
        character(len=*), intent(in) :: srname
        integer, intent(in) :: info
        integer :: c_info

        interface
            subroutine xerbla(name, code)
                character(len = *) :: name
                integer :: code
            end subroutine
        end interface

        if (.not. set_from_c .and.associated(global_error_handler)) then
            call global_error_handler(trim(srname), info, global_error_aux)
            return
        end if
        if (set_from_c .and.associated(global_error_handler_c)) then
            c_info = int(info, c_int)
            call global_error_handler_c(trim(srname) // c_null_char, c_info, global_error_aux_c)
            return
        end if

        call xerbla(trim(srname), info)
    end subroutine qrupdate_xerror

end module qrupdate_error
