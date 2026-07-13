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
!> \brief Updates a QR factorization after deleting a column.
!>
!> \par Definition:
! =============
!> \verbatim
!>       subroutine sqrdec(m,n,k,Q,ldq,R,ldr,j,w)
!>
!>       .. Scalar Arguments ..
!>       integer             m, n, k, ldq, ldr, j
!>       ..
!>       .. Array Arguments ..
!>       real                Q(ldq,*)
!>       real                R(ldr,*)
!>       real                w(*)
!>       ..
!> \endverbatim
!>
!> \par Purpose:
! =============
!> \verbatim
!>
!> SQRDEC updates a QR factorization after deleting a column. i.e.,
!> given an m-by-k orthogonal matrix Q, an k-by-n upper trapezoidal
!> matrix R and index j in the range 1:n+1, SQRDEC updates the
!> matrix Q -> Q1 and R -> R1 so that Q1 remains orthogonal, R1 is upper
!> trapezoidal, and Q1*R1 = [A(:,1:j-1) A(:,j+1:n)], where A = Q*R.
!> (real version)
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
!>          either k = m (full Q) or k = n < m (economical form,
!>          basis dimension will decrease).
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is REAL array, dimension (ldq,*)
!>          On entry, the orthogonal m-by-k matrix Q.  On exit,
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
!>          R is REAL array, dimension (ldr,*)
!>          On entry, the original matrix R.  On exit, the
!>          updated matrix R1.
!> \endverbatim
!>
!> \param[in] ldr
!> \verbatim
!>          ldr is INTEGER
!>          The leading dimension of R.  ldr >= k.
!> \endverbatim
!>
!> \param[in] j
!> \verbatim
!>          j is INTEGER
!>          The position of the deleted column in R.  1 <= j <= n.
!> \endverbatim
!>
!> \param[out] w
!> \verbatim
!>          w is REAL array, dimension (*)
!>          A workspace vector of size k-j.
!> \endverbatim
!>
!> \ingroup qrdecomp
subroutine sqrdec(m,n,k,Q,ldq,R,ldr,j,w)
    use iso_fortran_env
    use qrupdate_error
    integer, intent(in) :: m, n, k, ldq, ldr, j
    real(real32), intent(inout) :: Q(ldq,*), R(ldr,*)
    real(real32), intent(out) :: w(*)
    external scopy,sqhqr,sqrot
    integer info,i
    ! quick return if possible.
    if (m == 0 .or. n == 0 .or. j == n) return
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
    else if (ldr < k) then
        info = 7
    else if (j < 1 .or. j > n+1) then
        info = 8
    end if
    if (info /= 0) then
        call qrupdate_xerror('SQRDEC',info)
        return
    end if
    ! delete the j-th column.
    do i = j,n-1
        call scopy(k,R(1,i+1),1,R(1,i),1)
    end do
    ! retriangularize.
    if (j < k) then
        call sqhqr(k+1-j,n-j,R(j,j),ldr,w,R(1,n))
        ! apply rotations to Q.
        call sqrot('F',m,min(k,n)+1-j,Q(1,j),ldq,w,R(1,n))
    end if
end subroutine
