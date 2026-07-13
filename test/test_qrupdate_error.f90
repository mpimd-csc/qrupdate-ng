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

module test_aux
    implicit none
    type :: aux_data
        character(len=32) :: message
        integer :: val
    end type aux_data
end module test_aux

program test_qrupdate_error
    use iso_fortran_env
    use qrupdate_error
    use test_state
    use test_aux
    implicit none

    type(aux_data), target :: my_data
    class(*), pointer :: p_aux_data
    procedure(error_handler_if), pointer :: p_handler

    print *, 'Running qrupdate_error tests...'

    ! Test 1: Default behavior (no custom handler)
    call reset()
    call qrupdate_xerror('TEST_ROUTINE', 42)
    if (.not. xerbla_called) then
        print *, 'Test 1 Failed: xerbla was not called'
        stop
    end if
    print *, 'Test 1 Passed: Default behavior calls xerbla'

    ! Test 2: Custom handler without auxiliary data
    call reset()
    p_handler => my_custom_handler
    call qrupdate_set_error(p_handler)
    call qrupdate_xerror('CUSTOM_ROUTINE', 123)
    if (.not. custom_handler_called) then
        print *, 'Test 2 Failed: Custom handler was not called'
        stop
    end if
    if (last_srname /= 'CUSTOM_ROUTINE' .or. last_info /= 123) then
        print *, 'Test 2 Failed: Incorrect arguments passed to handler'
        stop
    end if
    print *, 'Test 2 Passed: Custom handler called correctly'

    ! Test 3: Custom handler with auxiliary data
    call reset()
    my_data%message = 'Hello World'
    my_data%val = 999
    p_aux_data => my_data
    call qrupdate_set_error_data(p_aux_data)
    call qrupdate_xerror('DATA_ROUTINE', 7)
    if (.not. custom_handler_called) then
        print *, 'Test 3 Failed: Custom handler was not called'
        stop
    end if
    print *, 'Test 3 Passed: Custom handler called with auxiliary data'

    ! Test 4: Overwriting handler
    call reset()
    p_handler => my_other_handler
    call qrupdate_set_error(p_handler)
    call qrupdate_xerror('OTHER_ROUTINE', 1)
    if (call_count /= 1) then
        print *, 'Test 4 Failed: Incorrect call count for other handler'
        stop
    end if
    print *, 'Test 4 Passed: Handler overwritten successfully'

    print *, 'All qrupdate_error tests passed!'

contains

    subroutine my_custom_handler(srname, info, aux)
        use test_state
        character(len=*), intent(in) :: srname
        integer, intent(in) :: info
        class(*), optional, intent(in) :: aux

 
        if (present(aux)) then 
            call_count = call_count
        endif        
        print *, 'my_custom_handler called from ', srname, ' with ', info 

        custom_handler_called = .true.
        last_srname = srname
        last_info = info
        call_count = call_count + 1
    end subroutine my_custom_handler

    subroutine my_other_handler(srname, info, aux)
        use test_state
        character(len=*), intent(in) :: srname
        integer, intent(in) :: info
        class(*), optional, intent(in) :: aux

        if (present(aux)) then 
            call_count = call_count
        endif 
        print *, 'my_other_handler called from ', srname, ' with ', info 
        custom_handler_called = .true.
        call_count = call_count + 1
    end subroutine my_other_handler

end program test_qrupdate_error


subroutine xerbla(srname, info)
    use test_state
    implicit none
    character(len=*), intent(in) :: srname
    integer, intent(in) :: info
    xerbla_called = .true.
    print *, 'Mock xerbla called from: ', srname, ' with info: ', info
end subroutine xerbla
