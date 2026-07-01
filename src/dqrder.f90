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
!> \brief Updates a QR factorization after deleting a row.
!>
!> \par Definition:
! =============
!> \verbatim
!>       subroutine dqrder(m,n,Q,ldq,R,ldr,j,w)
!>
!>       .. Scalar Arguments ..
!>       integer             m, n, ldq, ldr, j
!>       ..
!>       .. Array Arguments ..
!>       double precision    Q(ldq,*)
!>       double precision    R(ldr,*)
!>       double precision    w(*)
!>       ..
!> \endverbatim
!>
!> \par Purpose:
! =============
!> \verbatim
!>
!> DQRDER updates a QR factorization after deleting a row. i.e., given
!> an m-by-m orthogonal matrix Q, an m-by-n upper trapezoidal matrix R
!> and index j in the range 1:m, DQRDER updates Q ->Q1 and an R
!> -> R1 so that Q1 is again orthogonal, R1 upper trapezoidal, and Q1*R1
!> = [A(1:j-1,:); A(j+1:m,:)], where A = Q*R. (real version)
!> \endverbatim
!>
!> \param[in] m
!> \verbatim
!>          m is INTEGER
!>          The number of rows of the matrix Q.  m >= 1.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is INTEGER
!>          The number of columns of the matrix R.  n >= 0.
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (ldq,*)
!>          On entry, the orthogonal matrix Q.  On exit, the
!>          updated matrix Q1.
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
!>          R is DOUBLE PRECISION array, dimension (ldr,*)
!>          On entry, the original matrix R.  On exit, the
!>          updated matrix R1.
!> \endverbatim
!>
!> \param[in] ldr
!> \verbatim
!>          ldr is INTEGER
!>          The leading dimension of R.  ldr >= m.
!> \endverbatim
!>
!> \param[in] j
!> \verbatim
!>          j is INTEGER
!>          The position of the deleted row.  1 <= j <= m.
!> \endverbatim
!>
!> \param[out] w
!> \verbatim
!>          w is DOUBLE PRECISION array, dimension (*)
!>          A workspace vector of size 2*m.
!> \endverbatim
!>
!> \ingroup qrdecomp
subroutine dqrder(m,n,Q,ldq,R,ldr,j,w)
    integer, intent(in) :: m, n, j, ldq, ldr
    double precision, intent(inout) :: Q(ldq,*)
    double precision, intent(inout) :: R(ldr,*)
    double precision, intent(out) :: w(*)
    external xerbla,dcopy,dqrtv1,dqrot,dqrqh
    integer info,i,k
    ! quick return if possible.
    if (m == 1) return
    ! check arguments
    info = 0
    if (m < 1) then
        info = 1
    else if (j < 1 .or. j > m) then
        info = 7
    end if
    if (info /= 0) then
        call xerbla('DQRDER',info)
        return
    end if
    ! eliminate Q(j,2:m).
    call dcopy(m,Q(j,1),ldq,w,1)
    call dqrtv1(m,w,w(m+1))
    ! apply rotations to Q.
    call dqrot('B',m,m,Q,ldq,w(m+1),w(2))
    ! form Q1.
    do k = 1,m-1
        if (j > 1) call dcopy(j-1,Q(1,k+1),1,Q(1,k),1)
        if (j < m) call dcopy(m-j,Q(j+1,k+1),1,Q(j,k),1)
    end do
    ! apply rotations to R.
    call dqrqh(m,n,R,ldr,w(m+1),w(2))
    ! form R1.
    do k = 1,n
        do i = 1,m-1
            R(i,k) = R(i+1,k)
        end do
    end do
end subroutine
