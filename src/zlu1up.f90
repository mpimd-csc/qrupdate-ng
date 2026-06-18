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
!> \brief Updates an LU factorization after a rank-1 modification.
!>
!> \par Definition:
! =============
!> \verbatim
!>       subroutine zlu1up(m,n,L,ldl,R,ldr,u,v)
!>
!>       .. Scalar Arguments ..
!>       integer            m, n, ldl, ldr
!>       ..
!>       .. Array Arguments ..
!>       double complex     L(ldl,*), R(ldr,*), u(*), v(*)
!>       ..
!> \endverbatim
!>
!> \par Purpose:
! =============
!> \verbatim
!>
!> ZLU1UP updates an LU factorization after rank-1 modification.
!> Given an m-by-k lower-triangular matrix L with unit diagonal and
!> a k-by-n upper-trapezoidal matrix R, where k = min(m,n), this
!> ZLU1UP updates L -> L1 and R -> R1 so that L1 is again
!> lower unit triangular, R1 upper trapezoidal, and
!> L1*R1 = L*R + u*v', where v' denotes the conjugate transpose
!> of v.
!>
!> The update is performed using the Bennett algorithm with
!> column-major access, which processes the leading k-by-k block
!> first and then finishes the trailing part of R if needed.
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
!> \param[in,out] u
!> \verbatim
!>          u is COMPLEX*16 array, dimension (m)
!>          On entry, the left m-vector defining the rank-1
!>          modification.  On exit, if k < m, u is destroyed;
!>          otherwise, u contains the updated vector.
!> \endverbatim
!>
!> \param[in,out] v
!> \verbatim
!>          v is COMPLEX*16 array, dimension (n)
!>          On entry, the right n-vector defining the rank-1
!>          modification.  On exit, v is destroyed.
!> \endverbatim
!>
!> \ingroup ludecomp
subroutine zlu1up(m,n,L,ldl,R,ldr,u,v)
    integer m,n,ldl,ldr
    double complex L(ldl,*),R(ldr,*),u(*),v(*)
    double complex ui,vi
    integer k,info,i,j
    external xerbla
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
        call xerbla('ZLU1UP',info)
        return
    end if

    ! The Bennett algorithm, modified for column-major access.
    ! The leading part.
    do i = 1,k
        ! prefetch
        ui = u(i)
        vi = v(i)
        ! delayed R update
        do j = 1,i-1
            R(j,i) = R(j,i) + u(j)*vi
            vi = vi - v(j)*R(j,i)
        end do
        ! diagonal update
        R(i,i) = R(i,i) + ui*vi
        vi = vi/R(i,i)
        ! L update
        do j = i+1,m
            u(j) = u(j) - ui*L(j,i)
            L(j,i) = L(j,i) + u(j)*vi
        end do
        u(i) = ui
        v(i) = vi
    end do

    ! Finish the trailing part of R if needed.
    do i = k+1,n
        vi = v(i)
        do j = 1,k
            R(j,i) = R(j,i) + u(j)*vi
            vi = vi - v(j)*R(j,i)
        end do
        v(i) = vi
    end do
end subroutine