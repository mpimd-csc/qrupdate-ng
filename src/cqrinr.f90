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
!> \brief Updates a QR factorization after inserting a new row.
!>
!> \par Definition:
! =============
!> \verbatim
!>       subroutine cqrinr(m,n,Q,ldq,R,ldr,j,x,rw)
!>
!>       .. Scalar Arguments ..
!>       integer             m, n, ldq, ldr, j
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
!> CQRINR updates a QR factorization after inserting a new row. i.e.,
!> given an m-by-m unitary matrix Q, an m-by-n upper trapezoidal matrix
!> R and index j in the range 1:m+1, CQRINR updates Q -> Q1 and
!> R -> R1 so that Q1 is again unitary, R1 upper trapezoidal, and Q1*R1
!> = [A(1:j-1,:); x; A(j:m,:)], where A = Q*R. (complex version)
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
!> \param[in,out] Q
!> \verbatim
!>          Q is COMPLEX array, dimension (ldq,*)
!>          On entry, the unitary matrix Q.  On exit, the
!>          updated matrix Q1.
!> \endverbatim
!>
!> \param[in] ldq
!> \verbatim
!>          ldq is INTEGER
!>          The leading dimension of Q.  ldq >= m+1.
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
!>          The leading dimension of R.  ldr >= m+1.
!> \endverbatim
!>
!> \param[in] j
!> \verbatim
!>          j is INTEGER
!>          The position of the new row in R1.  1 <= j <= m+1.
!> \endverbatim
!>
!> \param[in,out] x
!> \verbatim
!>          x is COMPLEX array, dimension (*)
!>          On entry, the row being added.  On exit, x is
!>          destroyed.
!> \endverbatim
!>
!> \param[out] rw
!> \verbatim
!>          rw is REAL array, dimension (*)
!>          A real workspace vector of size min(m,n).
!> \endverbatim
!>
!> \ingroup qrdecomp
subroutine cqrinr(m,n,Q,ldq,R,ldr,j,x,rw)
    integer, intent(in) :: m, n, j, ldq, ldr
    complex, intent(inout) :: Q(ldq,*)
    complex, intent(inout) :: R(ldr,*)
    complex, intent(inout) :: x(*)
    real, intent(out) :: rw(*)
    external xerbla,ccopy,cqhqr,cqrot
    integer info,i,k
    ! check arguments
    info = 0
    if (n < 0) then
        info = 2
    else if (j < 1 .or. j > m+1) then
        info = 7
    end if
    if (info /= 0) then
        call xerbla('CQRINR',info)
        return
    end if
    ! permute the columns of Q1 and rows of R1 so that c the new row ends
    ! up being the topmost row of R1.
    do i = m,1,-1
        if (j > 1) then
            call ccopy(j-1,Q(1,i),1,Q(1,i+1),1)
        end if
        Q(j,i+1) = 0e0
        if (j <= m) then
            call ccopy(m+1-j,Q(j,i),1,Q(j+1,i+1),1)
        end if
    end do
    ! set up the 1st column
    do i = 1,j-1
        Q(i,1) = 0e0
    end do
    Q(j,1) = 1e0
    do i = j+1,m+1
        Q(i,1) = 0e0
    end do
    ! set up the new matrix R1
    do k = 1,n
        if (k < m) R(m+1,k) = 0e0
        do i = min(m,k),1,-1
            R(i+1,k) = R(i,k)
        end do
        R(1,k) = x(k)
    end do
    ! retriangularize R
    call cqhqr(m+1,n,R,ldr,rw,x)
    ! apply rotations to Q
    call cqrot('F',m+1,min(m,n)+1,Q,ldq,rw,x)
end subroutine
