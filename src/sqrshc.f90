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
subroutine sqrshc(m,n,k,Q,ldq,R,ldr,i,j,w)
    ! purpose:      updates a QR factorization after circular shift of
    !               columns.
    !               i.e., given an m-by-k orthogonal matrix Q, an k-by-n
    !               upper trapezoidal matrix R and index j in the range
    !               1:n+1, this subroutine updates the matrix Q -> Q1 and
    !               R -> R1 so that Q1 is again orthogonal, R1 upper
    !               trapezoidal, and
    !               Q1*R1 = A(:,p), where A = Q*R and p is the permutation
    !               [1:i-1,shift(i:j,-1),j+1:n] if i < j  or
    !               [1:j-1,shift(j:i,+1),i+1:n] if j < i.
    !               (real version)
    ! arguments:
    ! m (in)        number of rows of the matrix Q.
    ! n (in)        number of columns of the matrix R.
    ! k (in)        number of columns of Q1, and rows of R1. Must be
    !               either k = m (full Q) or k = n <= m (economical form).
    ! Q (io)        on entry, the unitary m-by-k matrix Q.
    !               on exit, the updated matrix Q1.
    ! ldq (in)      leading dimension of Q. ldq >= m.
    ! R (io)        on entry, the original matrix R.
    !               on exit, the updated matrix R1.
    ! ldr (in)      leading dimension of R. ldr >= k.
    ! i (in)        the first index determining the range (see above)
    ! j (in)        the second index determining the range (see above)
    ! w (o)         a workspace vector of size 2*k.
    !
    integer m,n,k,ldq,ldr,i,j
    real Q(ldq,*),R(ldr,*),w(*)
    external xerbla,scopy,sqrtv1,sqrqh,sqhqr,sqrot
    integer info,jj,kk,l
    ! quick return if possible.
    if (m == 0 .or. n == 1) return
    info = 0
    ! check arguments.
    if (m < 0) then
        info = 1
    else if (n < 0) then
        info = 2
    else if (k /= m .and. (k /= n .or. n > m)) then
        info = 3
    else if (i < 1 .or. i > n) then
        info = 6
    else if (j < 1 .or. j > n) then
        info = 7
    end if
    if (info /= 0) then
        call xerbla('SQRSHC',info)
        return
    end if

    if (i < j) then
        ! shift columns
        call scopy(k,R(1,i),1,w,1)
        do l = i,j-1
            call scopy(k,R(1,l+1),1,R(1,l),1)
        end do
        call scopy(k,w,1,R(1,j),1)
        ! retriangularize
        if (i < k) then
            kk = min(k,j)
            call sqhqr(kk+1-i,n+1-i,R(i,i),ldr,w(k+1),w)
            ! apply rotations to Q.
            call sqrot('F',m,kk+1-i,Q(1,i),ldq,w(k+1),w)
        end if
    else if (j < i) then
        ! shift columns
        call scopy(k,R(1,i),1,w,1)
        do l = i,j+1,-1
            call scopy(k,R(1,l-1),1,R(1,l),1)
        end do
        call scopy(k,w,1,R(1,j),1)
        ! retriangularize
        if (j < k) then
            jj = min(j+1,n)
            kk = min(k,i)
            ! eliminate the introduced spike.
            call sqrtv1(kk+1-j,R(j,j),w(k+1))
            ! apply rotations to R
            call sqrqh(kk+1-j,n-j,R(j,jj),ldr,w(k+1),R(j+1,j))
            ! apply rotations to Q
            call sqrot('B',m,kk+1-j,Q(1,j),ldq,w(k+1),R(j+1,j))
            ! zero spike.
            do l = j+1,kk
                R(l,j) = 0e0
            end do
        end if
    end if
end subroutine
