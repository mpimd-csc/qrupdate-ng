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
!> \brief Updates a Cholesky factorization after deleting a row and column.
!>
!> \par Definition:
! =============
!> \verbatim
!>       subroutine cchdex(n,R,ldr,j,rw)
!>
!>       .. Scalar Arguments ..
!>       integer            n, ldr, j
!>       ..
!>       .. Array Arguments ..
!>       complex            R(ldr,*)
!>       real               rw(*)
!>       ..
!> \endverbatim
!>
!> \par Purpose:
! =============
!> \verbatim
!>
!> CCHDEX updates the Cholesky factorization of a hermitian
!> positive definite matrix A after deleting a row/column.
!> Given an upper triangular matrix R that is a Cholesky
!> factor of A, i.e., A = R'*R, where R' denotes the conjugate
!> transpose of R, CCHDEX updates R -> R1 so that
!> R1'*R1 = A(jj,jj), where jj = [1:j-1, j+1:n+1].
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
!>          R is COMPLEX array, dimension (ldr,n)
!>          On entry, the upper triangular matrix R, the Cholesky
!>          factor of A.  On exit, the updated upper triangular
!>          matrix R1, the Cholesky factor of A(jj,jj).
!> \endverbatim
!>
!> \param[in] ldr
!> \verbatim
!>          ldr is INTEGER
!>          The leading dimension of the array R.  ldr >= n.
!> \endverbatim
!>
!> \param[in] j
!> \verbatim
!>          j is INTEGER
!>          The position of the deleted row/column.
!> \endverbatim
!>
!> \param[out] rw
!> \verbatim
!>          rw is REAL array, dimension (n)
!>          A workspace vector.
!> \endverbatim
!>
!> \ingroup choldecomp
subroutine cchdex(n,R,ldr,j,rw)
  use iso_fortran_env
  use qrupdate_error
    integer, intent(in) :: n, ldr, j
    complex(real32), intent(inout) :: R(ldr,*)
    real(real32), intent(out) :: rw(*)
    integer info,i
    external ccopy,cqhqr
    ! quick return if possible.
    if (n == 1) return
    ! check arguments
    info = 0
    if (n < 0) then
        info = 1
    else if (j < 1 .or. j > n) then
        info = 4
    end if
    if (info /= 0) then
        call qrupdate_xerror('CCHDEX',info)
        return
    end if
    ! delete the j-th column.
    do i = j,n-1
        call ccopy(n,R(1,i+1),1,R(1,i),1)
    end do
    ! retriangularize.
    if (j < n) then
        call cqhqr(n+1-j,n-j,R(j,j),ldr,rw,R(1,n))
    end if
end subroutine
