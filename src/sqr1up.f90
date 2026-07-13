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
!>       subroutine sqr1up(m,n,k,Q,ldq,R,ldr,u,v,w)
!>
!>       .. Scalar Arguments ..
!>       integer             m, n, k, ldq, ldr
!>       ..
!>       .. Array Arguments ..
!>       real                Q(ldq,*)
!>       real                R(ldr,*)
!>       real                u(*)
!>       real                v(*)
!>       real                w(*)
!>       ..
!> \endverbatim
!>
!> \par Purpose:
! =============
!> \verbatim
!>
!> SQR1UP updates a QR factorization after rank-1 modification i.e.,
!> given a m-by-k orthogonal Q and m-by-n upper trapezoidal R, an
!> m-vector u and n-vector v, SQR1UP updates Q -> Q1 and R ->
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
!>          u is REAL array, dimension (*)
!>          On entry, the left m-vector.  On exit, if k < m,
!>          u is destroyed.
!> \endverbatim
!>
!> \param[in,out] v
!> \verbatim
!>          v is REAL array, dimension (*)
!>          On entry, the right n-vector.  On exit, v is
!>          destroyed.
!> \endverbatim
!>
!> \param[out] w
!> \verbatim
!>          w is REAL array, dimension (*)
!>          A workspace vector of size 2*k.
!> \endverbatim
!>
!> \ingroup qrdecomp
subroutine sqr1up(m,n,k,Q,ldq,R,ldr,u,v,w)
  use iso_fortran_env
  use qrupdate_error
    integer, intent(in) :: m, n, k, ldq, ldr
    real(real32), intent(inout) :: Q(ldq,*), R(ldr,*)
    real(real32), intent(inout) :: u(*), v(*)
    real(real32), intent(out) :: w(*)
    external sqrqh,sqhqr,sqrot,sqrtv1,sch1up
    external saxpy,sdot,snrm2,slamch,sscal,srot
    real(real32) sdot,snrm2,slamch,ru,ruu
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
        call qrupdate_xerror('SQR1UP',info)
        return
    end if
    full = k == m
    ! in the non-full case, we shall need the norm of u.
    ru = 1.0
    if (.not.full) ru = snrm2(m,u,1)
    ! form Q'*u. In the non-full case, form also u - Q*Q'u.
    do i = 1,k
        w(i) = sdot(m,Q(1,i),1,u,1)
        if (.not.full) call saxpy(m,-w(i),Q(1,i),1,u,1)
    end do
    ! generate rotations to eliminate Q'*u.
    call sqrtv1(k,w,w(k+1))
    ! apply rotations to R.
    call sqrqh(k,n,R,ldr,w(k+1),w(2))
    ! apply rotations to Q.
    call sqrot('B',m,k,Q,ldq,w(k+1),w(2))
    ! update the first row of R.
    call saxpy(n,w(1),v,1,R(1,1),ldr)
    ! retriangularize R.
    call sqhqr(k,n,R,ldr,w(k+1),w)
    ! apply rotations to Q.
    call sqrot('F',m,min(k,n+1),Q,ldq,w(k+1),w)
    ! in the full case, we're finished
    if (full) return
    ! compute relative residual norm
    ruu = snrm2(m,u,1)
    ru = ru * slamch('e')
    if (ruu <= ru) return
    ! update the orthogonal basis.
    call sscal(n,ruu,v,1)
    call sscal(m,1e0/ruu,u,1)
    call sch1up(n,R,ldr,v,w(k+1))
    do i = 1,n
        call srot(m,Q(1,i),1,u,1,w(k+i),v(i))
    end do
end subroutine
