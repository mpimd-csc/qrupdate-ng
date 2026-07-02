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
!>       subroutine zqr1up(m,n,k,Q,ldq,R,ldr,u,v,w,rw)
!>
!>       .. Scalar Arguments ..
!>       integer             m, n, k, ldq, ldr
!>       ..
!>       .. Array Arguments ..
!>       double complex      Q(ldq,*)
!>       double complex      R(ldr,*)
!>       double complex      u(*)
!>       double complex      v(*)
!>       double complex      w(*)
!>       double precision    rw(*)
!>       ..
!> \endverbatim
!>
!> \par Purpose:
! =============
!> \verbatim
!>
!> ZQR1UP updates a QR factorization after rank-1 modification i.e.,
!> given a m-by-k unitary Q and m-by-n upper trapezoidal R, an m-vector
!> u and n-vector v, ZQR1UP updates Q -> Q1 and R -> R1 so that
!> Q1*R1 = Q*R + u*v', and Q1 is again unitary and R1 upper trapezoidal.
!> (complex version)
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
!>          Q is COMPLEX*16 array, dimension (ldq,*)
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
!>          R is COMPLEX*16 array, dimension (ldr,*)
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
!>          u is COMPLEX*16 array, dimension (*)
!>          On entry, the left m-vector.  On exit, if k < m,
!>          u is destroyed.
!> \endverbatim
!>
!> \param[in,out] v
!> \verbatim
!>          v is COMPLEX*16 array, dimension (*)
!>          On entry, the right n-vector.  On exit, v is
!>          destroyed.
!> \endverbatim
!>
!> \param[out] w
!> \verbatim
!>          w is COMPLEX*16 array, dimension (*)
!>          A workspace vector of size k.
!> \endverbatim
!>
!> \param[out] rw
!> \verbatim
!>          rw is DOUBLE PRECISION array, dimension (*)
!>          A real workspace vector of size k.
!> \endverbatim
!>
!> \ingroup qrdecomp
subroutine zqr1up(m,n,k,Q,ldq,R,ldr,u,v,w,rw)
  use iso_fortran_env
    integer, intent(in) :: m, n, k, ldq, ldr
    complex(real64), intent(inout) :: Q(ldq,*), R(ldr,*), u(*), v(*)
    complex(real64), intent(out) :: w(*)
    real(real64), intent(out) :: rw(*)
    external zqrqh,zqhqr,zqrot,zqrtv1,zaxpy,zaxcpy,xerbla,zch1up
    external zdotc,dznrm2,dlamch,zdscal,zrot
    complex(real64) zdotc
    real(real64) dznrm2,dlamch,ru,ruu
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
        call xerbla('ZQR1UP',info)
        return
    end if

    full = k == m
    ! in the non-full case, we shall need the norm of u.
    if (.not.full) ru = dznrm2(m,u,1)
    ! form Q'*u. In the non-full case, form also u - Q*Q'u.
    do i = 1,k
        w(i) = zdotc(m,Q(1,i),1,u,1)
        if (.not.full) call zaxpy(m,-w(i),Q(1,i),1,u,1)
    end do
    ! generate rotations to eliminate Q'*u.
    call zqrtv1(k,w,rw)
    ! apply rotations to R.
    call zqrqh(k,n,R,ldr,rw,w(2))
    ! apply rotations to Q.
    call zqrot('B',m,k,Q,ldq,rw,w(2))
    ! update the first row of R.
    call zaxcpy(n,w(1),v,1,R(1,1),ldr)
    ! retriangularize R.
    call zqhqr(k,n,R,ldr,rw,w)
    ! apply rotations to Q.
    call zqrot('F',m,min(k,n+1),Q,ldq,rw,w)
    ! in the full case, we're finished
    if (full) return
    ! compute relative residual norm
    ruu = dznrm2(m,u,1)
    ru = ru * dlamch('e')
    if (ruu <= ru) return
    ! update the orthogonal basis.
    call zdscal(n,ruu,v,1)
    call zdscal(m,1d0/ruu,u,1)
    call zch1up(n,R,ldr,v,rw)
    do i = 1,n
        call zrot(m,Q(1,i),1,u,1,rw(i),conjg(v(i)))
    end do
end subroutine
