! Copyright (C) 2008, 2009  VZLU Prague, a.s., Czech Republic, Jaroslav Hajek <highegg@gmail.com>
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
!> \brief Generates a unit vector orthogonal to the column space of a unitary matrix.
!>
!> \par Definition:
! =============
!> \verbatim
!>       subroutine sgqvec(m,n,Q,ldq,u)
!>
!>       .. Scalar Arguments ..
!>       integer            m, n, ldq
!>       ..
!>       .. Array Arguments ..
!>       real               Q(ldq,*), u(*)
!>       ..
!> \endverbatim
!>
!> \par Purpose:
! =============
!> \verbatim
!>
!> SGQVEC generates a vector u in the orthogonal complement of the
!> column space of an orthogonal matrix Q.  Given an m-by-n
!> orthogonal matrix Q with n < m, SGQVEC generates a
!> vector u of length m such that Q.'*u = 0 and norm(u) = 1, where
!> Q.' denotes the transpose of Q.
!>
!> The algorithm projects canonical unit vectors onto the orthogonal
!> complement of Q's column space until a nonzero result is found.
!> If n = 0, the first canonical unit vector is returned.
!> \endverbatim
!>
!> \param[in] m
!> \verbatim
!>          m is INTEGER
!>          The number of rows of the matrix Q.  m >= 0.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is INTEGER
!>          The number of columns of the matrix Q.  n >= 0 and
!>          n < m.
!> \endverbatim
!>
!> \param[in] Q
!> \verbatim
!>          Q is REAL array, dimension (ldq,n)
!>          The orthogonal m-by-n matrix Q.
!> \endverbatim
!>
!> \param[in] ldq
!> \verbatim
!>          ldq is INTEGER
!>          The leading dimension of the array Q.  ldq >= m.
!> \endverbatim
!>
!> \param[out] u
!> \verbatim
!>          u is REAL array, dimension (m)
!>          The generated vector such that Q.'*u = 0 and norm(u) = 1.
!> \endverbatim
!>
!> \ingroup qrdecomp
subroutine sgqvec(m,n,Q,ldq,u)
    use iso_fortran_env
    use qrupdate_error
    integer, intent(in) :: m, n, ldq
    real(real32), intent(in) :: Q(ldq,*)
    real(real32), intent(out) :: u(*)
    external sdot,saxpy,snrm2,sscal
    real(real32) sdot,snrm2,r
    integer info,i,j
    ! quick return if possible.
    if (m == 0) return
    if (n == 0) then
        u(1) = 1e0
        do i = 2,m
            u(i) = 0e0
        end do
        return
    end if
    ! check arguments.
    info = 0
    if (m < 0) then
        info = 1
    else if (n < 0) then
        info = 2
    else if (ldq < m) then
        info = 4
    end if
    if (info /= 0) then
        call qrupdate_xerror('SGQVEC',info)
        return
    end if
    j = 1
    r = 0e0
    do while ( r .eq. 0e0 )
        ! probe j-th canonical unit vector.
        do i = 1,m
            u(i) = 0e0
        end do
        u(j) = 1e0
        ! form u - Q*Q'*u
        do i = 1,n
            r = sdot(m,Q(1,i),1,u,1)
            call saxpy(m,-r,Q(1,i),1,u,1)
        end do
        r = snrm2(m,u,1)
        if (r == 0e0) then
            j = j + 1
            if (j > n) then
                ! this is fatal, and in theory, it can't happen.
                stop 'fatal: impossible condition in SGQVEC'
            end if
        end if
    end do
    call sscal(m,1e0/r,u,1)
end subroutine
