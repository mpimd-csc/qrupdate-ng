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
!> \brief Updates a row-pivoted LU factorization after rank-1 modification.
!>
!> \par Definition:
! =============
!> \verbatim
!>       subroutine zlup1up(m,n,L,ldl,R,ldr,p,u,v,w)
!>
!>       .. Scalar Arguments ..
!>       integer            m, n, ldl, ldr
!>       ..
!>       .. Array Arguments ..
!>       integer            p(*)
!>       double complex     L(ldl,*), R(ldr,*), u(*), v(*), w(*)
!>       ..
!> \endverbatim
!>
!> \par Purpose:
! =============
!> \verbatim
!>
!> ZLUP1UP updates a row-pivoted LU factorization after rank-1
!> modification.  Given an m-by-k lower-triangular matrix L with
!> unit diagonal, a k-by-n upper-trapezoidal matrix R, and a
!> permutation vector p, where k = min(m,n), ZLUP1UP
!> updates L -> L1, R -> R1 and p -> p1 so that L1 is again
!> lower unit triangular, R1 upper trapezoidal, p1 a permutation,
!> and P1'*L1*R1 = P'*L*R + u*v', where v' denotes the conjugate
!> transpose of v and P is the permutation matrix corresponding
!> to p.
!> \endverbatim
!>
!> \param[in] m
!> \verbatim
!>          m is INTEGER
!>          The number of rows of the matrix L.  m >= 0.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is INTEGER
!>          The number of columns of the matrix R.  n >= 0.
!> \endverbatim
!>
!> \param[in,out] L
!> \verbatim
!>          L is COMPLEX*16 array, dimension (ldl,k)
!>          On entry, the unit lower triangular matrix L.  On exit,
!>          the updated unit lower triangular matrix L1.
!> \endverbatim
!>
!> \param[in] ldl
!> \verbatim
!>          ldl is INTEGER
!>          The leading dimension of the array L.  ldl >= m.
!> \endverbatim
!>
!> \param[in,out] R
!> \verbatim
!>          R is COMPLEX*16 array, dimension (ldr,n)
!>          On entry, the upper trapezoidal m-by-n matrix R.
!>          On exit, the updated upper trapezoidal matrix R1.
!> \endverbatim
!>
!> \param[in] ldr
!> \verbatim
!>          ldr is INTEGER
!>          The leading dimension of the array R.  ldr >= k,
!>          where k = min(m,n).
!> \endverbatim
!>
!> \param[in] p
!> \verbatim
!>          p is INTEGER array, dimension (m)
!>          The permutation vector representing the row pivoting.
!>          On exit, p is updated to reflect the new pivoting.
!> \endverbatim
!>
!> \param[in] u
!> \verbatim
!>          u is COMPLEX*16 array, dimension (m)
!>          The left m-vector defining the rank-1 modification.
!> \endverbatim
!>
!> \param[in] v
!> \verbatim
!>          v is COMPLEX*16 array, dimension (n)
!>          The right n-vector defining the rank-1 modification.
!> \endverbatim
!>
!> \param[out] w
!> \verbatim
!>          w is COMPLEX*16 array, dimension (m)
!>          Workspace vector used during the update computation.
!> \endverbatim
!>
!> \ingroup ludecomp
subroutine zlup1up(m,n,L,ldl,R,ldr,p,u,v,w)
    use iso_fortran_env
    use qrupdate_error
    integer, intent(in) :: m, n, ldl, ldr
    integer, intent(inout) :: p(*)
    complex(real64), intent(inout) :: L(ldl,*), R(ldr,*)
    complex(real64), intent(in) :: u(*), v(*)
    complex(real64), intent(out) :: w(*)
    complex(real64) one,tmp
    real(real64) tau
    parameter (one = 1d0, tau = 1d-1)
    integer k,info,i,j,itmp
    external zcopy,zaxpy,ztrsv,zgeru,zgemv,zswap
    ! quick return if possible.
    k = min(m,n)
    if (k == 0) return
    ! check arguments.
    info = 0
    if (m < 0) then
        info = 1
    else if (n < 0) then
        info = 2
    else if (ldl < m) then
        info = 4
    else if (ldr < k) then
        info = 6
    endif
    if (info /= 0) then
        call qrupdate_xerror('ZLUP1UP',info)
        return
    end if
    ! form L \ P*u.
    do i = 1,m
        w(i) = u(p(i))
    end do
    call ztrsv('L','N','U',k,L,ldl,w,1)
    ! if m > k = n, subtract the trailing part.
    if (m > k) then
        call zgemv('N',m-k,k,-one,L(k+1,1),ldl,w,1,one,w(k+1),1)
    end if
    ! work from bottom to top
    do j = k-1,1,-1
        if (abs(w(j)) < tau * abs(L(j+1,j)*w(j) + w(j+1))) then
            ! need pivoting. swap j and j+1
            tmp = w(j)
            w(j) = w(j+1)
            w(j+1) = tmp
            ! update p
            itmp = p(j)
            p(j) = p(j+1)
            p(j+1) = itmp
            ! update L
            call zswap(m-j+1,L(j,j),1,L(j,j+1),1)
            call zswap(j+1,L(j,1),ldl,L(j+1,1),ldl)
            ! update R
            call zswap(n-j+1,R(j,j),ldr,R(j+1,j),ldr)
            ! make L lower triangular again
            tmp = -L(j,j+1)
            call zaxpy(m-j+1,tmp,L(j,j),1,L(j,j+1),1)
            ! update R
            call zaxpy(n-j+1,-tmp,R(j+1,j),ldr,R(j,j),ldr)
            ! update w
            w(j) = w(j) - tmp*w(j+1)
        end if
        ! eliminate w(j+1)
        tmp = w(j+1)/w(j)
        w(j+1) = 0
        ! update R.
        call zaxpy(n-j+1,-tmp,R(j,j),ldr,R(j+1,j),ldr)
        ! update L.
        call zaxpy(m-j,tmp,L(j+1,j+1),1,L(j+1,j),1)
    end do
    ! add a multiple of v to R
    call zaxpy(n,w(1),v,1,R(1,1),ldr)
    ! forward sweep
    do j = 1,k-1
        if (abs(R(j,j)) < tau * abs(L(j+1,j)*R(j,j) + R(j+1,j))) then
            ! need pivoting. swap j and j+1
            ! update p
            itmp = p(j)
            p(j) = p(j+1)
            p(j+1) = itmp
            ! update L
            call zswap(m-j+1,L(j,j),1,L(j,j+1),1)
            call zswap(j+1,L(j,1),ldl,L(j+1,1),ldl)
            ! update R
            call zswap(n-j+1,R(j,j),ldr,R(j+1,j),ldr)
            ! make L lower triangular again
            tmp = -L(j,j+1)
            call zaxpy(m-j+1,tmp,L(j,j),1,L(j,j+1),1)
            ! update R
            call zaxpy(n-j+1,-tmp,R(j+1,j),ldr,R(j,j),ldr)
        end if
        ! eliminate R(j+1,j)
        tmp = R(j+1,j)/R(j,j)
        ! update R.
        R(j+1,j) = 0d0
        call zaxpy(n-j,-tmp,R(j,j+1),ldr,R(j+1,j+1),ldr)
        ! update L.
        call zaxpy(m-j,tmp,L(j+1,j+1),1,L(j+1,j),1)
    end do
    ! if m > k = n, complete the update by updating the lower part of L.
    if (m > k) then
        call zcopy(k,v,1,w,1)
        call ztrsv('U','T','N',k,R,ldr,w,1)
        call zgeru(m-k,k,one,w(k+1),1,w,1,L(k+1,1),ldl)
    endif
end subroutine
