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
!> \brief Updates a Cholesky factorization after a symmetric shift of rows and columns.
!>
!> \par Definition:
! =============
!> \verbatim
!>       subroutine dchshx(n,R,ldr,i,j,w)
!>
!>       .. Scalar Arguments ..
!>       integer            n, ldr, i, j
!>       ..
!>       .. Array Arguments ..
!>       double precision   R(ldr,*), w(*)
!>       ..
!> \endverbatim
!>
!> \par Purpose:
! =============
!> \verbatim
!>
!> DCHSHX updates the Cholesky factorization of a symmetric
!> positive definite matrix A after a symmetric shift of rows and
!> columns.  Given an upper triangular matrix R that is a Cholesky
!> factor of A, i.e., A = R.'*R, where R.' denotes the transpose
!> of R, DCHSHX updates R -> R1 so that
!> R1.'*R1 = A(p,p), where p is the permutation
!> [1:i-1, shift(i:j,-1), j+1:n] if i < j, or
!> [1:j-1, shift(j:i,+1), i+1:n] if j < i.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is INTEGER
!>          The order of matrix R.  n >= 0.
!> \endverbatim
!>
!> \param[in,out] R
!> \verbatim
!>          R is DOUBLE PRECISION array, dimension (ldr,n)
!>          On entry, the upper triangular matrix R, the Cholesky
!>          factor of A.  On exit, the updated upper triangular
!>          matrix R1, the Cholesky factor of A(p,p).
!> \endverbatim
!>
!> \param[in] ldr
!> \verbatim
!>          ldr is INTEGER
!>          The leading dimension of the array R.  ldr >= n.
!> \endverbatim
!>
!> \param[in] i
!> \verbatim
!>          i is INTEGER
!>          The first index determining the range of the shift.
!>          1 <= i <= n.
!> \endverbatim
!>
!> \param[in] j
!> \verbatim
!>          j is INTEGER
!>          The second index determining the range of the shift.
!>          1 <= j <= n.
!> \endverbatim
!>
!> \param[out] w
!> \verbatim
!>          w is DOUBLE PRECISION array, dimension (2*n)
!>          Workspace vector used during the retriangularization.
!> \endverbatim
!>
!> \ingroup choldecomp
subroutine dchshx(n,R,ldr,i,j,w)
    use iso_fortran_env
    use qrupdate_error
    integer, intent(in) :: n, ldr, i, j
    real(real64), intent(inout) :: R(ldr,*)
    real(real64), intent(out) :: w(*)
    external dcopy,dqrtv1,dqrqh,dqhqr
    integer info,l
    ! quick return if possible.
    if (n == 0 .or. n == 1) return
    info = 0
    ! check arguments.
    if (n < 0) then
        info = 1
    else if (i < 1 .or. i > n) then
        info = 4
    else if (j < 1 .or. j > n) then
        info = 5
    end if
    if (info /= 0) then
        call qrupdate_xerror('DCHSHX',info)
        return
    end if
    if (i < j) then
        ! shift columns
        call dcopy(n,R(1,i),1,w,1)
        do l = i,j-1
            call dcopy(n,R(1,l+1),1,R(1,l),1)
        end do
        call dcopy(n,w,1,R(1,j),1)
        ! retriangularize
        call dqhqr(n+1-i,n+1-i,R(i,i),ldr,w(n+1),w)
    else if (j < i) then
        ! shift columns
        call dcopy(n,R(1,i),1,w,1)
        do l = i,j+1,-1
            call dcopy(n,R(1,l-1),1,R(1,l),1)
        end do
        call dcopy(n,w,1,R(1,j),1)
        ! eliminate the introduced spike.
        call dqrtv1(n+1-j,R(j,j),w(n+1))
        ! apply rotations to R
        call dqrqh(n+1-j,n-j,R(j,j+1),ldr,w(n+1),R(j+1,j))
        ! zero spike.
        do l = j+1,n
            R(l,j) = 0d0
        end do
    end if
end subroutine
