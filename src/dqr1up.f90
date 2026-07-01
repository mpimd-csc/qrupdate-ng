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
!> \brief Updates a QR factorization after a rank-1 modification.
!>
!> \par Definition:
! =============
!> \verbatim
!>       subroutine dqr1up(m,n,k,Q,ldq,R,ldr,u,v,w)
!>
!>       .. Scalar Arguments ..
!>       integer             m, n, k, ldq, ldr
!>       ..
!>       .. Array Arguments ..
!>       double precision    Q(ldq,*)
!>       double precision    R(ldr,*)
!>       double precision    u(*)
!>       double precision    v(*)
!>       double precision    w(*)
!>       ..
!> \endverbatim
!>
!> \par Purpose:
! =============
!> \verbatim
!>
!> DQR1UP updates a QR factorization after rank-1 modification i.e.,
!> given a m-by-k orthogonal Q and m-by-n upper trapezoidal R, an
!> m-vector u and n-vector v, DQR1UP updates Q -> Q1 and R ->
!> R1 so that Q1*R1 = Q*R + u*v', and Q1 is again orthonormal and R1
!> upper trapezoidal. (real version)
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
!>          either k = m (full Q) or k = n < m (economical form).
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (ldq,*)
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
!>          R is DOUBLE PRECISION array, dimension (ldr,*)
!>          On entry, the upper trapezoidal m-by-n matrix R.  On
!>          exit, the updated matrix R1.
!> \endverbatim
!>
!> \param[in] ldr
!> \verbatim
!>          ldr is INTEGER
!>          The leading dimension of R.  ldr >= k.
!> \endverbatim
!>
!> \param[in,out] u
!> \verbatim
!>          u is DOUBLE PRECISION array, dimension (*)
!>          On entry, the left m-vector.  On exit, if k < m,
!>          u is destroyed.
!> \endverbatim
!>
!> \param[in,out] v
!> \verbatim
!>          v is DOUBLE PRECISION array, dimension (*)
!>          On entry, the right n-vector.  On exit, v is
!>          destroyed.
!> \endverbatim
!>
!> \param[out] w
!> \verbatim
!>          w is DOUBLE PRECISION array, dimension (*)
!>          A workspace vector of size 2*k.
!> \endverbatim
!>
!> \ingroup qrdecomp
subroutine dqr1up(m,n,k,Q,ldq,R,ldr,u,v,w)
    integer, intent(in) :: m, n, k, ldq, ldr
    double precision, intent(inout) :: Q(ldq,*)
    double precision, intent(inout) :: R(ldr,*)
    double precision, intent(inout) :: u(*)
    double precision, intent(inout) :: v(*)
    double precision, intent(out) :: w(*)
    external xerbla, dch1up, dqrqh,dqhqr,dqrot,dqrtv1
    external daxpy,ddot,dnrm2,dlamch,dscal,drot
    double precision ddot,dnrm2,dlamch,ru,ruu
    integer info,i
    logical full
    ! quick return if possible.
    if (k == 0 .or. n == 0) return
    ! check arguments.
    info = 0
    if (m < 0) then
        info = 1
    else if (n < 0) then
        info = 2
    else if (k /= m .and. (k /= n .or. n > m)) then
        info = 3
    else if (ldq < m) then
        info = 5
    else if (ldr < k) then
        info = 7
    endif
    if (info /= 0) then
        call xerbla('DQR1UP',info)
        return
    end if

    full = k == m
    ! in the non-full case, we shall need the norm of u.
    if (.not.full) ru = dnrm2(m,u,1)
    ! form Q'*u. In the non-full case, form also u - Q*Q'u.
    do i = 1,k
        w(i) = ddot(m,Q(1,i),1,u,1)
        if (.not.full) call daxpy(m,-w(i),Q(1,i),1,u,1)
    end do
    ! generate rotations to eliminate Q'*u.
    call dqrtv1(k,w,w(k+1))
    ! apply rotations to R.
    call dqrqh(k,n,R,ldr,w(k+1),w(2))
    ! apply rotations to Q.
    call dqrot('B',m,k,Q,ldq,w(k+1),w(2))
    ! update the first row of R.
    call daxpy(n,w(1),v,1,R(1,1),ldr)
    ! retriangularize R.
    call dqhqr(k,n,R,ldr,w(k+1),w)
    ! apply rotations to Q.
    call dqrot('F',m,min(k,n+1),Q,ldq,w(k+1),w)
    ! in the full case, we're finished
    if (full) return
    ! compute relative residual norm
    ruu = dnrm2(m,u,1)
    ru = ru * dlamch('e')
    if (ruu <= ru) return
    ! update the orthogonal basis.
    call dscal(n,ruu,v,1)
    call dscal(m,1d0/ruu,u,1)
    call dch1up(n,R,ldr,v,w(k+1))
    do i = 1,n
        call drot(m,Q(1,i),1,u,1,w(k+i),v(i))
    end do
end subroutine
