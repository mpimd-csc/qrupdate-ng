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
!> \brief Converts an upper trapezoidal matrix to upper Hessenberg form.
!>
!> \par Definition:
! =============
!> \verbatim
!>       subroutine zqrqh(m,n,R,ldr,c,s)
!>
!>       .. Scalar Arguments ..
!>       integer             m, n, ldr
!>       ..
!>       .. Array Arguments ..
!>       double complex      R(ldr,*)
!>       double precision    c(*)
!>       double complex      s(*)
!>       ..
!> \endverbatim
!>
!> \par Purpose:
! =============
!> \verbatim
!>
!> ZQRQH brings an upper trapezoidal matrix R into upper Hessenberg form
!> using min(m-1,n) Givens rotations. (complex version)
!> \endverbatim
!>
!> \param[in] m
!> \verbatim
!>          m is INTEGER
!>          The number of rows of the matrix R.  m >= 0.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is INTEGER
!>          The number of columns of the matrix R.  n >= 0.
!> \endverbatim
!>
!> \param[in,out] R
!> \verbatim
!>          R is COMPLEX*16 array, dimension (ldr,*)
!>          On entry, the upper Hessenberg matrix R.  On exit,
!>          the updated upper trapezoidal matrix.
!> \endverbatim
!>
!> \param[in] ldr
!> \verbatim
!>          ldr is INTEGER
!>          The leading dimension of R.  ldr >= m.
!> \endverbatim
!>
!> \param[in] c
!> \verbatim
!>          c is DOUBLE PRECISION array, dimension (*)
!>          The rotation cosines.  Must contain at least
!>          min(m-1,n) elements.
!> \endverbatim
!>
!> \param[in] s
!> \verbatim
!>          s is COMPLEX*16 array, dimension (*)
!>          The rotation sines.  Must contain at least
!>          min(m-1,n) elements.
!> \endverbatim
!>
!> \ingroup qrdecomp
subroutine zqrqh(m,n,R,ldr,c,s)
    integer, intent(in) :: m, n, ldr
    double complex, intent(inout) :: R(ldr,*)
    double precision, intent(in) :: c(*)
    double complex, intent(in) :: s(*)
    external xerbla
    double complex t
    integer info,i,ii,j
    ! quick return if possible.
    if (m == 0 .or. m == 1 .or. n == 0) return
    ! check arguments.
    info = 0
    if (m < 0) then
        info = 1
    else if (n < 0) then
        info = 2
    else if (ldr < m) then
        info = 4
    end if
    if (info /= 0) then
        call xerbla('ZQRQH',info)
        return
    end if
    do i = 1,n
        ii = min(m-1,i)
        ! apply stored rotations, column-wise
        t = R(ii+1,i)
        do j = ii,1,-1
            R(j+1,i) = c(j)*t - conjg(s(j))*R(j,i)
            t = c(j)*R(j,i) + s(j)*t
        end do
        R(1,i) = t
    end do
end subroutine
