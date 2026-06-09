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
subroutine dqrder(m,n,Q,ldq,R,ldr,j,w)
    ! purpose:      updates a QR factorization after deleting a row.
    !               i.e., given an m-by-m orthogonal matrix Q, an m-by-n
    !               upper trapezoidal matrix R and index j in the range
    !               1:m, this subroutine updates Q ->Q1 and an R -> R1
    !               so that Q1 is again orthogonal, R1 upper trapezoidal,
    !               and Q1*R1 = [A(1:j-1,:); A(j+1:m,:)], where A = Q*R.
    !               (real version)
    !
    ! arguments:
    ! m (in)        number of rows of the matrix Q.
    ! n (in)        number of columns of the matrix R.
    ! Q (io)        on entry, the orthogonal matrix Q.
    !               on exit, the updated matrix Q1.
    ! ldq (in)      leading dimension of Q. ldq >= m.
    ! R (io)        on entry, the original matrix R.
    !               on exit, the updated matrix R1.
    ! ldr (in)      leading dimension of R. ldr >= m.
    ! j (in)        the position of the deleted row.
    ! w (out)       a workspace vector of size 2*m.
    !
    integer m,n,j,ldq,ldr
    double precision Q(ldq,*),R(ldr,*),w(*)
    external xerbla,dcopy,dqrtv1,dqrot,dqrqh
    integer info,i,k
    ! quick return if possible
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
