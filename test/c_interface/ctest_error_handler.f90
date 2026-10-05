module ctest_error_handler
    use iso_c_binding
    use qrupdate_error
    use test_state
    implicit none
contains
    subroutine ctest_set_handler() bind(C, name='ctest_set_handler')
        procedure(error_handler_if), pointer :: p_handler
        p_handler => pxerbla
        call qrupdate_set_error(p_handler)
    end subroutine ctest_set_handler

    subroutine ctest_reset_state() bind(C, name='ctest_reset_state')
        call reset()
    end subroutine ctest_reset_state

    function ctest_get_xerbla_called() result(called) bind(C, name='ctest_get_xerbla_called')
        logical(C_BOOL) :: called
        called = xerbla_called
    end function ctest_get_xerbla_called

    function ctest_get_last_info() result(info) bind(C, name='ctest_get_last_info')
        integer(C_INT) :: info
        info = last_info
    end function ctest_get_last_info

    function ctest_srname_eq(expected) result(match) bind(C, name='ctest_srname_eq')
        character(C_CHAR), intent(in) :: expected(*)
        logical(C_BOOL) :: match
        character(len=64) :: fstr
        integer :: i, n
        n = 0
        do while (expected(n+1) /= C_NULL_CHAR)
            n = n + 1
        end do
        fstr = ' '
        fstr(1:n) = ''
        do i = 1, n
            fstr(i:i) = expected(i)
        end do
        match = (trim(last_srname) == trim(fstr(1:n)))
    end function ctest_srname_eq

    subroutine pxerbla(srname, info, aux)
        character(len=*), intent(in) :: srname
        integer, intent(in) :: info
        class(*), optional, intent(in) :: aux
        xerbla_called = .true.
        last_srname = trim(srname)
        last_info = info
    end subroutine pxerbla
end module ctest_error_handler
