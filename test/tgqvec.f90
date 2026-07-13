! Copyright (C) 2026 Martin K&ouml;hler <koehlerm(AT)mpi-magdeburg.mpg.de>
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
program tgqvec
    use iso_fortran_env
    implicit none
    integer :: m, n, passed, failed
    integer :: i
    real(real32) :: slamch
    real(real64) :: dlamch
    external :: slamch, dlamch

    common /stats/ passed, failed
    passed = 0
    failed = 0

    write(*,*)
    write(*,*) 'Testing gqvec routines.'
    write(*,*)
    write(*,1000)

    ! Test cases: (m, n) pairs
    ! - n = 0: trivial case, u should be e_1
    ! - n = 1: orthogonal complement is m-1 dimensional
    ! - n = m-1: orthogonal complement is 1 dimensional
    ! - n = m div 2: moderate complement dimension
    ! - n = m div 2 - 1: slightly smaller complement
    ! - large m: stress test

    do i = 1, 6
        select case (i)
        case (1)
            m = 10
            n = 0
        case (2)
            m = 10
            n = 1
        case (3)
            m = 10
            n = m - 1
        case (4)
            m = 10
            n = m / 2
        case (5)
            m = 10
            n = m / 2 - 1
        case (6)
            m = 100
            n = 50
        end select

        write(*,*) 'Test case (m,n) = (', m, ',', n, '):'
        call stest(m, n)
        call dtest(m, n)
        call ctest(m, n)
        call ztest(m, n)
    end do

    write(*,1000)
    write(*,1001) passed, failed
    write(*,*)
    if (failed .ne. 0) stop 1

    1000 format(70('-'))
    1001 format('total:', 1x, 'PASSED', 1x, I6, 1x, 'FAILED', 1x, I6)
end program tgqvec

! Convert integer to string
subroutine int2str(n, s)
    integer, intent(in) :: n
    character(len=*), intent(out) :: s
    character(len=12) :: tmp
    write(tmp, '(I12)') n
    s = adjustl(tmp)
end subroutine int2str

subroutine stest(m, n)
    use iso_fortran_env
    integer, intent(in) :: m, n
    real(real32), allocatable :: A(:,:), Q(:,:), R(:,:), u(:), Qcheck(:,:)
    real(real32) :: slamch, tol
    integer :: j
    character(len=40) :: lbl
    integer :: passed, failed
    common /stats/ passed, failed

    external :: sgqvec, srandg, sqrgen, slamch
    real(real32) :: sdot, snrm2

    allocate(A(m, m), Q(m, m), R(m, m), u(m), Qcheck(m+3, m))

    call srandg(m, n, A, m)
    call sqrgen(m, n, A, m, Q, m, R, m)

    ! Test with ldq = m (tight packing)
    call sgqvec(m, n, Q, m, u)
    tol = 5e2 * slamch('p')
    call int2str(m, lbl)
    if (n == 0) then
        call check_unit_vector_s(m, u, tol, 's (n=0, ldq=' // lbl // ')')
    else
        call check_orthogonal_s(m, n, Q, m, u, tol, 's (n>0, ldq=' // lbl // ')')
        call check_unit_norm_s(m, u, tol, 's (n>0, ldq=' // lbl // ')')
    end if

    ! Test with ldq > m (leading dimension padding)
    do j = 1, 3
        call srandg(m, n, A, m)
        call sqrgen(m, n, A, m, Q, m, R, m)
        if (n == 0) then
            call sgqvec(m, n, Q, m+j, u)
            tol = 5e2 * slamch('p')
            call int2str(m + j, lbl)
            call check_unit_vector_s(m, u, tol, 's (n=0, ldq=' // lbl // ')')
        else
            call sgqvec_ldq(m, n, Q, m+j, u, Qcheck)
            tol = 5e2 * slamch('p')
            call int2str(m + j, lbl)
            call check_orthogonal_s(m, n, Qcheck, m+j, u, tol, 's (n>0, ldq=' // lbl // ')')
            call check_unit_norm_s(m, u, tol, 's (n>0, ldq=' // lbl // ')')
        end if
    end do

    deallocate(A, Q, R, u)

end subroutine stest

subroutine dtest(m, n)
    use iso_fortran_env
    integer, intent(in) :: m, n
    real(real64), allocatable :: A(:,:), Q(:,:), R(:,:), u(:), Qcheck(:,:)
    real(real64) :: dlamch, tol
    integer :: j
    character(len=40) :: lbl
    integer :: passed, failed
    common /stats/ passed, failed

    external :: dgqvec, drandg, dqrgen, dlamch
    real(real64) :: ddot, dnrm2

    allocate(A(m, m), Q(m, m), R(m, m), u(m), Qcheck(m+3, m))

    call drandg(m, n, A, m)
    call dqrgen(m, n, A, m, Q, m, R, m)

    call dgqvec(m, n, Q, m, u)
    tol = 5d2 * dlamch('p')
    call int2str(m, lbl)
    if (n == 0) then
        call check_unit_vector_d(m, u, tol, 'd (n=0, ldq=' // lbl // ')')
    else
        call check_orthogonal_d(m, n, Q, m, u, tol, 'd (n>0, ldq=' // lbl // ')')
        call check_unit_norm_d(m, u, tol, 'd (n>0, ldq=' // lbl // ')')
    end if

    ! Test with ldq > m
    do j = 1, 3
        call drandg(m, n, A, m)
        call dqrgen(m, n, A, m, Q, m, R, m)
        if (n == 0) then
            call dgqvec(m, n, Q, m+j, u)
            tol = 5d2 * dlamch('p')
            call int2str(m + j, lbl)
            call check_unit_vector_d(m, u, tol, 'd (n=0, ldq=' // lbl // ')')
        else
            call dgqvec_ldq(m, n, Q, m+j, u, Qcheck)
            tol = 5d2 * dlamch('p')
            call int2str(m + j, lbl)
            call check_orthogonal_d(m, n, Qcheck, m+j, u, tol, 'd (n>0, ldq=' // lbl // ')')
            call check_unit_norm_d(m, u, tol, 'd (n>0, ldq=' // lbl // ')')
        end if
    end do

    deallocate(A, Q, R, u)

end subroutine dtest

subroutine ctest(m, n)
    use iso_fortran_env
    use qrupdate_blas
    integer, intent(in) :: m, n
    complex(real32), allocatable :: A(:,:), Q(:,:), R(:,:), u(:), Qcheck(:,:)
    real(real32) :: slamch, tol
    integer :: j
    character(len=40) :: lbl
    integer :: passed, failed
    common /stats/ passed, failed

    external :: cgqvec, crandg, cqrgen

    allocate(A(m, m), Q(m, m), R(m, m), u(m), Qcheck(m+3, m))

    call crandg(m, n, A, m)
    call cqrgen(m, n, A, m, Q, m, R, m)

    call cgqvec(m, n, Q, m, u)
    tol = 5e2 * slamch('p')
    call int2str(m, lbl)
    if (n == 0) then
        call check_unit_vector_c(m, u, tol, 'c (n=0, ldq=' // lbl // ')')
    else
        call check_orthogonal_c(m, n, Q, m, u, tol, 'c (n>0, ldq=' // lbl // ')')
        call check_unit_norm_c(m, u, tol, 'c (n>0, ldq=' // lbl // ')')
    end if

    ! Test with ldq > m
    do j = 1, 3
        call crandg(m, n, A, m)
        call cqrgen(m, n, A, m, Q, m, R, m)
        if (n == 0) then
            call cgqvec(m, n, Q, m+j, u)
            tol = 5e2 * slamch('p')
            call int2str(m + j, lbl)
            call check_unit_vector_c(m, u, tol, 'c (n=0, ldq=' // lbl // ')')
        else
            call cgqvec_ldq(m, n, Q, m+j, u, Qcheck)
            tol = 5e2 * slamch('p')
            call int2str(m + j, lbl)
            call check_orthogonal_c(m, n, Qcheck, m+j, u, tol, 'c (n>0, ldq=' // lbl // ')')
            call check_unit_norm_c(m, u, tol, 'c (n>0, ldq=' // lbl // ')')
        end if
    end do

    deallocate(A, Q, R, u)

end subroutine ctest

subroutine ztest(m, n)
    use iso_fortran_env
    use qrupdate_blas
    integer, intent(in) :: m, n
    complex(real64), allocatable :: A(:,:), Q(:,:), R(:,:), u(:), Qcheck(:,:)
    real(real64) :: dlamch, tol
    integer :: j
    character(len=40) :: lbl
    integer :: passed, failed
    common /stats/ passed, failed

    external :: zgqvec, zrandg, zqrgen

    allocate(A(m, m), Q(m, m), R(m, m), u(m), Qcheck(m+3, m))

    call zrandg(m, n, A, m)
    call zqrgen(m, n, A, m, Q, m, R, m)

    call zgqvec(m, n, Q, m, u)
    tol = 5d2 * dlamch('p')
    call int2str(m, lbl)
    if (n == 0) then
        call check_unit_vector_z(m, u, tol, 'z (n=0, ldq=' // lbl // ')')
    else
        call check_orthogonal_z(m, n, Q, m, u, tol, 'z (n>0, ldq=' // lbl // ')')
        call check_unit_norm_z(m, u, tol, 'z (n>0, ldq=' // lbl // ')')
    end if

    ! Test with ldq > m
    do j = 1, 3
        call zrandg(m, n, A, m)
        call zqrgen(m, n, A, m, Q, m, R, m)
        if (n == 0) then
            call zgqvec(m, n, Q, m+j, u)
            tol = 5d2 * dlamch('p')
            call int2str(m + j, lbl)
            call check_unit_vector_z(m, u, tol, 'z (n=0, ldq=' // lbl // ')')
        else
            call zgqvec_ldq(m, n, Q, m+j, u, Qcheck)
            tol = 5d2 * dlamch('p')
            call int2str(m + j, lbl)
            call check_orthogonal_z(m, n, Qcheck, m+j, u, tol, 'z (n>0, ldq=' // lbl // ')')
            call check_unit_norm_z(m, u, tol, 'z (n>0, ldq=' // lbl // ')')
        end if
    end do

    deallocate(A, Q, R, u)

end subroutine ztest

! Helper: call gqvec with a padded Q buffer for ldq > m.
! Creates Qpad(ldq, n) with Q in the first m rows, calls gqvec, returns Qpad via Qout.
subroutine sgqvec_ldq(m, n, Q, ldq, u, Qout)
    use iso_fortran_env
    integer, intent(in) :: m, n, ldq
    real(real32), intent(in) :: Q(m, *)
    real(real32), intent(out) :: u(m)
    real(real32), intent(out) :: Qout(ldq, n)
    external sgqvec
    if (n == 0) then
        call sgqvec(m, n, Q, ldq, u)
        return
    end if
    Qout = 0e0
    Qout(1:m, 1:n) = Q(1:m, 1:n)
    call sgqvec(m, n, Qout, ldq, u)
end subroutine sgqvec_ldq

subroutine dgqvec_ldq(m, n, Q, ldq, u, Qout)
    use iso_fortran_env
    integer, intent(in) :: m, n, ldq
    real(real64), intent(in) :: Q(m, *)
    real(real64), intent(out) :: u(m)
    real(real64), intent(out) :: Qout(ldq, n)
    external dgqvec
    if (n == 0) then
        call dgqvec(m, n, Q, ldq, u)
        return
    end if
    Qout = 0d0
    Qout(1:m, 1:n) = Q(1:m, 1:n)
    call dgqvec(m, n, Qout, ldq, u)
end subroutine dgqvec_ldq

subroutine cgqvec_ldq(m, n, Q, ldq, u, Qout)
    use iso_fortran_env
    integer, intent(in) :: m, n, ldq
    complex(real32), intent(in) :: Q(m, *)
    complex(real32), intent(out) :: u(m)
    complex(real32), intent(out) :: Qout(ldq, n)
    external cgqvec
    if (n == 0) then
        call cgqvec(m, n, Q, ldq, u)
        return
    end if
    Qout = (0e0, 0e0)
    Qout(1:m, 1:n) = Q(1:m, 1:n)
    call cgqvec(m, n, Qout, ldq, u)
end subroutine cgqvec_ldq

subroutine zgqvec_ldq(m, n, Q, ldq, u, Qout)
    use iso_fortran_env
    integer, intent(in) :: m, n, ldq
    complex(real64), intent(in) :: Q(m, *)
    complex(real64), intent(out) :: u(m)
    complex(real64), intent(out) :: Qout(ldq, n)
    external zgqvec
    if (n == 0) then
        call zgqvec(m, n, Q, ldq, u)
        return
    end if
    Qout = (0d0, 0d0)
    Qout(1:m, 1:n) = Q(1:m, 1:n)
    call zgqvec(m, n, Qout, ldq, u)
end subroutine zgqvec_ldq

! Check that u is approximately the first canonical unit vector e_1
subroutine check_unit_vector_s(m, u, tol, label)
    use iso_fortran_env
    integer, intent(in) :: m
    real(real32), intent(in) :: u(m)
    real(real32), intent(in) :: tol
    character(*), intent(in) :: label
    real(real32) :: err
    integer :: passed, failed
    common /stats/ passed, failed
    external :: snrm2
    real(real32) :: snrm2
    character(6) :: result

    err = abs(u(1) - 1e0)
    if (err > tol) then
        result = 'FAIL'
        failed = failed + 1
    else
        result = 'PASS'
        passed = passed + 1
    end if
    write(*, 1001) trim(label), 'u(1)=1  :', err, result

    if (m > 1) then
        err = snrm2(m - 1, u(2), 1)
    else
        err = 0e0
    end if
    if (err > tol) then
        result = 'FAIL'
        failed = failed + 1
    else
        result = 'PASS'
        passed = passed + 1
    end if
    write(*, 1002) trim(label), 'u(2:m)=0:', err, result

    return
    1001 format('   ', A, 1x, A10, E12.4, 1x, A4)
    1002 format('   ', A, 1x, A10, E12.4, 1x, A4)
end subroutine check_unit_vector_s

subroutine check_unit_vector_d(m, u, tol, label)
    use iso_fortran_env
    integer, intent(in) :: m
    real(real64), intent(in) :: u(m)
    real(real64), intent(in) :: tol
    character(*), intent(in) :: label
    real(real64) :: err
    integer :: passed, failed
    common /stats/ passed, failed
    external :: dnrm2
    real(real64) :: dnrm2
    character(6) :: result

    err = abs(u(1) - 1d0)
    if (err > tol) then
        result = 'FAIL'
        failed = failed + 1
    else
        result = 'PASS'
        passed = passed + 1
    end if
    write(*, 1001) trim(label), 'u(1)=1  :', err, result

    if (m > 1) then
        err = dnrm2(m - 1, u(2), 1)
    else
        err = 0d0
    end if
    if (err > tol) then
        result = 'FAIL'
        failed = failed + 1
    else
        result = 'PASS'
        passed = passed + 1
    end if
    write(*, 1002) trim(label), 'u(2:m)=0:', err, result

    return
    1001 format('   ', A, 1x, A10, E12.4, 1x, A4)
    1002 format('   ', A, 1x, A10, E12.4, 1x, A4)
end subroutine check_unit_vector_d

subroutine check_unit_vector_c(m, u, tol, label)
    use iso_fortran_env
    use qrupdate_blas
    integer, intent(in) :: m
    complex(real32), intent(in) :: u(m)
    real(real32), intent(in) :: tol
    character(*), intent(in) :: label
    real(real32) :: err
    integer :: passed, failed
    common /stats/ passed, failed
    external :: scnrm2
    real(real32) :: scnrm2
    character(6) :: result

    err = abs(abs(u(1)) - 1e0)
    if (err > tol) then
        result = 'FAIL'
        failed = failed + 1
    else
        result = 'PASS'
        passed = passed + 1
    end if
    write(*, 1001) trim(label), '|u(1)|=1:', err, result

    if (m > 1) then
        err = scnrm2(m - 1, u(2), 1)
    else
        err = 0e0
    end if
    if (err > tol) then
        result = 'FAIL'
        failed = failed + 1
    else
        result = 'PASS'
        passed = passed + 1
    end if
    write(*, 1002) trim(label), 'u(2:m)=0:', err, result

    return
    1001 format('   ', A, 1x, A10, E12.4, 1x, A4)
    1002 format('   ', A, 1x, A10, E12.4, 1x, A4)
end subroutine check_unit_vector_c

subroutine check_unit_vector_z(m, u, tol, label)
    use iso_fortran_env
    use qrupdate_blas
    integer, intent(in) :: m
    complex(real64), intent(in) :: u(m)
    real(real64), intent(in) :: tol
    character(*), intent(in) :: label
    real(real64) :: err
    integer :: passed, failed
    common /stats/ passed, failed
    external :: dznrm2
    real(real64) :: dznrm2
    character(6) :: result

    err = abs(abs(u(1)) - 1d0)
    if (err > tol) then
        result = 'FAIL'
        failed = failed + 1
    else
        result = 'PASS'
        passed = passed + 1
    end if
    write(*, 1001) trim(label), '|u(1)|=1:', err, result

    if (m > 1) then
        err = dznrm2(m - 1, u(2), 1)
    else
        err = 0d0
    end if
    if (err > tol) then
        result = 'FAIL'
        failed = failed + 1
    else
        result = 'PASS'
        passed = passed + 1
    end if
    write(*, 1002) trim(label), 'u(2:m)=0:', err, result

    return
    1001 format('   ', A, 1x, A10, E12.4, 1x, A4)
    1002 format('   ', A, 1x, A10, E12.4, 1x, A4)
end subroutine check_unit_vector_z

subroutine check_orthogonal_s(m, n, Q, ldq, u, tol, label)
    use iso_fortran_env
    integer, intent(in) :: m, n, ldq
    real(real32), intent(in) :: Q(ldq, n), u(m)
    real(real32), intent(in) :: tol
    character(*), intent(in) :: label
    real(real32) :: prod, err
    integer :: passed, failed
    common /stats/ passed, failed
    external :: sdot
    real(real32) :: sdot
    integer :: i
    character(6) :: result

    err = 0e0
    do i = 1, n
        prod = sdot(m, Q(1, i), 1, u, 1)
        if (abs(prod) > err) err = abs(prod)
    end do
    if (err > tol) then
        result = 'FAIL'
        failed = failed + 1
    else
        result = 'PASS'
        passed = passed + 1
    end if
    write(*, 1001) trim(label), 'max|Q''u|:', err, result

    return
    1001 format('   ', A, 1x, A10, E12.4, 1x, A4)
end subroutine check_orthogonal_s

subroutine check_orthogonal_d(m, n, Q, ldq, u, tol, label)
    use iso_fortran_env
    integer, intent(in) :: m, n, ldq
    real(real64), intent(in) :: Q(ldq, n), u(m)
    real(real64), intent(in) :: tol
    character(*), intent(in) :: label
    real(real64) :: prod, err
    integer :: passed, failed
    common /stats/ passed, failed
    external :: ddot
    real(real64) :: ddot
    integer :: i
    character(6) :: result

    err = 0d0
    do i = 1, n
        prod = ddot(m, Q(1, i), 1, u, 1)
        if (abs(prod) > err) err = abs(prod)
    end do
    if (err > tol) then
        result = 'FAIL'
        failed = failed + 1
    else
        result = 'PASS'
        passed = passed + 1
    end if
    write(*, 1001) trim(label), 'max|Q''u|:', err, result

    return
    1001 format('   ', A, 1x, A10, E12.4, 1x, A4)
end subroutine check_orthogonal_d

subroutine check_orthogonal_c(m, n, Q, ldq, u, tol, label)
    use iso_fortran_env
    use qrupdate_blas
    integer, intent(in) :: m, n, ldq
    complex(real32), intent(in) :: Q(ldq, n), u(m)
    real(real32), intent(in) :: tol
    character(*), intent(in) :: label
    real(real32) :: err, rr
    complex(real32) :: rc
    integer :: i
    integer :: passed, failed
    common /stats/ passed, failed
    character(6) :: result

    err = 0e0
    do i = 1, n
        call qrupdate_cdotc(rc, m, Q(1, i), 1, u, 1)
        rr = sqrt(real(rc, real32)**2 + aimag(rc)**2)
        if (rr > err) err = rr
    end do
    if (err > tol) then
        result = 'FAIL'
        failed = failed + 1
    else
        result = 'PASS'
        passed = passed + 1
    end if
    write(*, 1001) trim(label), 'max|Q''u|:', err, result

    return
    1001 format('   ', A, 1x, A10, E12.4, 1x, A4)
end subroutine check_orthogonal_c

subroutine check_orthogonal_z(m, n, Q, ldq, u, tol, label)
    use iso_fortran_env
    use qrupdate_blas
    integer, intent(in) :: m, n, ldq
    complex(real64), intent(in) :: Q(ldq, n), u(m)
    real(real64), intent(in) :: tol
    character(*), intent(in) :: label
    real(real64) :: err, rr
    complex(real64) :: rc
    integer :: i
    integer :: passed, failed
    common /stats/ passed, failed
    character(6) :: result

    err = 0d0
    do i = 1, n
        call qrupdate_zdotc(rc, m, Q(1, i), 1, u, 1)
        rr = dsqrt(real(rc, real64)**2 + aimag(rc)**2)
        if (rr > err) err = rr
    end do
    if (err > tol) then
        result = 'FAIL'
        failed = failed + 1
    else
        result = 'PASS'
        passed = passed + 1
    end if
    write(*, 1001) trim(label), 'max|Q''u|:', err, result

    return
    1001 format('   ', A, 1x, A10, E12.4, 1x, A4)
end subroutine check_orthogonal_z

subroutine check_unit_norm_s(m, u, tol, label)
    use iso_fortran_env
    integer, intent(in) :: m
    real(real32), intent(in) :: u(m)
    real(real32), intent(in) :: tol
    character(*), intent(in) :: label
    real(real32) :: unrm, err
    integer :: passed, failed
    common /stats/ passed, failed
    external :: snrm2
    real(real32) :: snrm2
    character(6) :: result

    unrm = snrm2(m, u, 1)
    err = abs(unrm - 1e0)
    if (err > tol) then
        result = 'FAIL'
        failed = failed + 1
    else
        result = 'PASS'
        passed = passed + 1
    end if
    write(*, 1001) trim(label), 'norm(u) :', unrm, result

    return
    1001 format('   ', A, 1x, A10, E12.4, 1x, A4)
end subroutine check_unit_norm_s

subroutine check_unit_norm_d(m, u, tol, label)
    use iso_fortran_env
    integer, intent(in) :: m
    real(real64), intent(in) :: u(m)
    real(real64), intent(in) :: tol
    character(*), intent(in) :: label
    real(real64) :: unrm, err
    integer :: passed, failed
    common /stats/ passed, failed
    external :: dnrm2
    real(real64) :: dnrm2
    character(6) :: result

    unrm = dnrm2(m, u, 1)
    err = abs(unrm - 1d0)
    if (err > tol) then
        result = 'FAIL'
        failed = failed + 1
    else
        result = 'PASS'
        passed = passed + 1
    end if
    write(*, 1001) trim(label), 'norm(u) :', unrm, result

    return
    1001 format('   ', A, 1x, A10, E12.4, 1x, A4)
end subroutine check_unit_norm_d

subroutine check_unit_norm_c(m, u, tol, label)
    use iso_fortran_env
    use qrupdate_blas
    integer, intent(in) :: m
    complex(real32), intent(in) :: u(m)
    real(real32), intent(in) :: tol
    character(*), intent(in) :: label
    real(real32) :: unrm, err
    integer :: passed, failed
    common /stats/ passed, failed
    external :: scnrm2
    real(real32) :: scnrm2
    character(6) :: result

    unrm = scnrm2(m, u, 1)
    err = abs(unrm - 1e0)
    if (err > tol) then
        result = 'FAIL'
        failed = failed + 1
    else
        result = 'PASS'
        passed = passed + 1
    end if
    write(*, 1001) trim(label), 'norm(u) :', unrm, result

    return
    1001 format('   ', A, 1x, A10, E12.4, 1x, A4)
end subroutine check_unit_norm_c

subroutine check_unit_norm_z(m, u, tol, label)
    use iso_fortran_env
    use qrupdate_blas
    integer, intent(in) :: m
    complex(real64), intent(in) :: u(m)
    real(real64), intent(in) :: tol
    character(*), intent(in) :: label
    real(real64) :: unrm, err
    integer :: passed, failed
    common /stats/ passed, failed
    external :: dznrm2
    real(real64) :: dznrm2
    character(6) :: result

    unrm = dznrm2(m, u, 1)
    err = abs(unrm - 1d0)
    if (err > tol) then
        result = 'FAIL'
        failed = failed + 1
    else
        result = 'PASS'
        passed = passed + 1
    end if
    write(*, 1001) trim(label), 'norm(u) :', unrm, result

    return
    1001 format('   ', A, 1x, A10, E12.4, 1x, A4)
end subroutine check_unit_norm_z
