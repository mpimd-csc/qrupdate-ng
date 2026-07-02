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
!> \brief Updates a QR factorization after a circular shift of columns.
!>
!> \par Definition:
! =============
!> \verbatim
!>       subroutine cqrshc(m,n,k,Q,ldq,R,ldr,i,j,w,rw)
!>
!>       .. Scalar Arguments ..
!>       integer             m, n, k, ldq, ldr, i, j
!>       ..
!>       .. Array Arguments ..
!>       complex             Q(ldq,*)
!>       complex             R(ldr,*)
!>       complex             w(*)
!>       real                rw(*)
!>       ..
!> \endverbatim
!>
!> \par Purpose:
! =============
!> \verbatim
!>
!> CQRSHC updates a QR factorization after circular shift of columns.
!> i.e., given an m-by-k unitary matrix Q, an k-by-n upper trapezoidal
!> matrix R and index j in the range 1:n+1, CQRSHC updates the
!> matrix Q -> Q1 and R -> R1 so that Q1 is again unitary, R1 upper
!> trapezoidal, and Q1*R1 = A(:,p), where A = Q*R and p is the
!> permutation [1:i-1,shift(i:j,-1),j+1:n] if i < j or
!> [1:j-1,shift(j:i,+1),i+1:n] if j < i. (complex version)
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
!>          The number of columns of Q1, and rows of R1.  Must be
!>          either k = m (full Q) or k = n <= m (economical form).
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
!>          The leading dimension of R.  ldr >= k.
!> \endverbatim
!>
!> \param[in] i
!> \verbatim
!>          i is INTEGER
!>          The first index determining the range (see above).
!> \endverbatim
!>
!> \param[in] j
!> \verbatim
!>          j is INTEGER
!>          The second index determining the range (see above).
!> \endverbatim
!>
!> \param[out] w
!> \verbatim
!>          w is COMPLEX array, dimension (*)
!>          A workspace vector of size k.
!> \endverbatim
!>
!> \param[out] rw
!> \verbatim
!>          rw is REAL array, dimension (*)
!>          A real workspace vector of size k.
!> \endverbatim
!>
!> \ingroup qrdecomp
subroutine cqrshc(m,n,k,Q,ldq,R,ldr,i,j,w,rw)
  use iso_fortran_env
    integer, intent(in) :: m, n, k, ldq, ldr, i, j
    complex(real32), intent(inout) :: Q(ldq,*)
    complex(real32), intent(inout) :: R(ldr,*)
    complex(real32), intent(out) :: w(*)
    real(real32), intent(out) :: rw(*)
    external xerbla,ccopy,cqrtv1,cqrqh,cqhqr,cqrot
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
        call xerbla('CQRSHC',info)
        return
    end if

    if (i < j) then
        ! shift columns
        call ccopy(k,R(1,i),1,w,1)
        do l = i,j-1
            call ccopy(k,R(1,l+1),1,R(1,l),1)
        end do
        call ccopy(k,w,1,R(1,j),1)
        ! retriangularize
        if (i < k) then
            kk = min(k,j)
            call cqhqr(kk+1-i,n+1-i,R(i,i),ldr,rw,w)
            ! apply rotations to Q.
            call cqrot('F',m,kk+1-i,Q(1,i),ldq,rw,w)
        end if
    else if (j < i) then
        ! shift columns
        call ccopy(k,R(1,i),1,w,1)
        do l = i,j+1,-1
            call ccopy(k,R(1,l-1),1,R(1,l),1)
        end do
        call ccopy(k,w,1,R(1,j),1)
        ! retriangularize
        if (j < k) then
            jj = min(j+1,n)
            kk = min(k,i)
            ! eliminate the introduced spike.
            call cqrtv1(kk+1-j,R(j,j),rw)
            ! apply rotations to R
            call cqrqh(kk+1-j,n-j,R(j,jj),ldr,rw,R(j+1,j))
            ! apply rotations to Q
            call cqrot('B',m,kk+1-j,Q(1,j),ldq,rw,R(j+1,j))
            ! zero spike.
            do l = j+1,kk
                R(l,j) = 0e0
            end do
        end if
    end if
end subroutine
