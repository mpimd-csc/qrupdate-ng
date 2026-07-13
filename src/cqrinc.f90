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
!> \brief Updates a QR factorization after inserting a new column.
!>
!> \par Definition:
! =============
!> \verbatim
!>       subroutine cqrinc(m,n,k,Q,ldq,R,ldr,j,x,rw)
!>
!>       .. Scalar Arguments ..
!>       integer             m, n, k, ldq, ldr, j
!>       ..
!>       .. Array Arguments ..
!>       complex             Q(ldq,*)
!>       complex             R(ldr,*)
!>       complex             x(*)
!>       real                rw(*)
!>       ..
!> \endverbatim
!>
!> \par Purpose:
! =============
!> \verbatim
!>
!> CQRINC updates a QR factorization after inserting a new column. i.e.,
!> given an m-by-k unitary matrix Q, an m-by-n upper trapezoidal matrix
!> R and index j in the range 1:n+1, CQRINC updates the matrix
!> Q -> Q1 and R -> R1 so that Q1 is again unitary, R1 upper
!> trapezoidal, and Q1*R1 = [A(:,1:j-1); x; A(:,j:n)], where A = Q*R.
!> (complex version)
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
!>          The number of columns of the matrix R.  n >= 0.
!> \endverbatim
!>
!> \param[in] k
!> \verbatim
!>          k is INTEGER
!>          The number of columns of Q, and rows of R.  Must be
!>          either k = m (full Q) or k = n <= m (economical form,
!>          basis dimension will increase).
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is COMPLEX array, dimension (ldq,*)
!>          On entry, the unitary m-by-k matrix Q.  On exit,
!>          the updated matrix Q1.
!> \endverbatim
!>
!> \param[in] ldq
!> \verbatim
!>          ldq is INTEGER
!>          The leading dimension of Q.  ldq >= m.
!> \endverbatim
!>
!> \param[in,out] R
!> \verbatim
!>          R is COMPLEX array, dimension (ldr,*)
!>          On entry, the original matrix R.  On exit, the
!>          updated matrix R1.
!> \endverbatim
!>
!> \param[in] ldr
!> \verbatim
!>          ldr is INTEGER
!>          The leading dimension of R.  ldr >= min(m,n+1).
!> \endverbatim
!>
!> \param[in] j
!> \verbatim
!>          j is INTEGER
!>          The position of the new column in R1.  1 <= j <= n+1.
!> \endverbatim
!>
!> \param[in] x
!> \verbatim
!>          x is COMPLEX array, dimension (*)
!>          The column being inserted.
!> \endverbatim
!>
!> \param[out] rw
!> \verbatim
!>          rw is REAL array, dimension (*)
!>          A real workspace vector of size k.
!> \endverbatim
!>
!> \ingroup qrdecomp
subroutine cqrinc(m,n,k,Q,ldq,R,ldr,j,x,rw)
    use iso_fortran_env
    use qrupdate_blas
    use qrupdate_error
    integer, intent(in) :: m, n, k, ldq, ldr, j
    complex(real32), intent(inout) :: Q(ldq,*), R(ldr,*)
    complex(real32), intent(in) :: x(*)
    real(real32), intent(out) :: rw(*)
    external cgqvec, cqrtv1,cqrqh,cqrot
    external ccopy,caxpy,csscal,scnrm2
    real(real32) scnrm2,rx
    integer info,i,k1
    logical full
    ! quick return if possible.
    if (m == 0) return
    ! check arguments.
    info = 0
    if (m < 0) then
        info = 1
    else if (n < 0) then
        info = 2
    else if (k /= m .and. (k /= n .or. n >= m)) then
        info = 3
    else if (ldq < m) then
        info = 5
    else if (ldr < min(m,k+1)) then
        info = 7
    else if (j < 1 .or. j > n+1) then
        info = 8
    end if
    if (info /= 0) then
        call qrupdate_xerror('CQRINC',info)
        return
    end if
    full = k == m
    ! insert empty column at j-th position
    do i = n,j,-1
        call ccopy(k,R(1,i),1,R(1,i+1),1)
    end do
    ! insert Q'*u into R. In the nonfull case, form also u-Q*Q'*u.
    if (full) then
        k1 = k
        do i = 1,k
            call qrupdate_cdotc(R(i,j), m,Q(1,i),1,x,1)
        end do
    else
        k1 = k + 1
        ! zero last row of R
        do i = 1,n+1
            R(k1,i) = 0e0
        end do
        call ccopy(m,x,1,Q(1,k1),1)
        do i = 1,k
            call qrupdate_cdotc(R(i,j), m,Q(1,i),1,Q(1,k1),1)
            call caxpy(m,-R(i,j),Q(1,i),1,Q(1,k1),1)
        end do
        ! get norm of the inserted column
        rx = scnrm2(m,Q(1,k1),1)
        R(k1,j) = rx
        if (rx == 0e0) then
            ! in the rare case when rx is exact zero, we still need to provide
            ! a valid orthogonal unit vector. The details are boring, so handle
            ! that elsewhere.
            call cgqvec(m,k,Q,ldq,Q(1,k1))
        else
            ! otherwise, just normalize the added column.
            call csscal(m,1e0/rx,Q(1,k1),1)
        end if
    end if
    ! maybe we're finished.
    if (j > k) return
    ! eliminate the spike.
    call cqrtv1(k1+1-j,R(j,j),rw)
    ! apply rotations to R(j:k,j:n).
    if (j <= n) call cqrqh(k1+1-j,n+1-j,R(j,j+1),ldr,rw,R(j+1,j))
    ! apply rotations to Q(:,j:k).
    call cqrot('B',m,k1+1-j,Q(1,j),ldq,rw,R(j+1,j))
    ! zero spike.
    do i = j+1,k1
        R(i,j) = 0e0
    end do
end subroutine
