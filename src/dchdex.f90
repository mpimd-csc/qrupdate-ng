! Copyright (C) 2008, 2009  VZLU Prague, a.s., Czech Republic
!
! Author: Jaroslav Hajek <highegg@gmail.com>
!
! This file is part of qrupdate.
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
subroutine dchdex(n,R,ldr,j,w)
    ! purpose:      given an upper triangular matrix R that is a Cholesky
    !               factor of a symmetric positive definite matrix A, i.e.
    !               A = R'*R, this subroutine updates R -> R1 so that
    !               R1'*R1 = A(jj,jj), where jj = [1:j-1,j+1:n+1].
    !               (real version)
    ! arguments:
    ! n (in)        the order of matrix R.
    ! R (io)        on entry, the original upper trapezoidal matrix R.
    !               on exit, the updated matrix R1.
    ! ldr (in)      leading dimension of R. ldr >= n.
    ! j (in)        the position of the deleted row/column.
    ! w (out)       a workspace vector of size n.
    !
    integer n,ldr,j
    double precision R(ldr,*),w(*)
    integer info,i
    external xerbla,dcopy,dqhqr

    ! quick return if possible
    if (n == 1) return

    ! check arguments
    info = 0
    if (n < 0) then
        info = 1
    else if (j < 1 .or. j > n) then
        info = 4
    end if
    if (info /= 0) then
        call xerbla('DCHDEX',info)
        return
    end if

    ! delete the j-th column.
    do i = j,n-1
        call dcopy(n,R(1,i+1),1,R(1,i),1)
    end do
    ! retriangularize.
    if (j < n) then
        call dqhqr(n+1-j,n-j,R(j,j),ldr,w,R(1,n))
    end if
end subroutine
