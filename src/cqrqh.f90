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
!>       subroutine cqrqh(m,n,R,ldr,c,s)
!>
!>       .. Scalar Arguments ..
!>       integer            m, n, ldr
!>       ..
!>       .. Array Arguments ..
!>       complex            R(ldr,*), s(*)
!>       real               c(*)
!>       ..
!> \endverbatim
!>
!> \par Purpose:
! =============
!> \verbatim
!>
!> CQRQH brings an m-by-n upper trapezoidal matrix R into upper
!> Hessenberg form.  Given an m-by-n upper trapezoidal matrix R,
!> CQRQH applies min(m-1,n) inverse Givens rotations
!> from the right to introduce subdiagonal elements, producing an
!> upper Hessenberg matrix.
!>
!> On exit, c contains the cosine parts and s contains the sine
!> parts of the Givens rotations used in the transformation.
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
!>          R is COMPLEX array, dimension (ldr,n)
!>          On entry, the upper trapezoidal matrix R.  On exit, the
!>          upper Hessenberg matrix.
!> \endverbatim
!>
!> \param[in] ldr
!> \verbatim
!>          ldr is INTEGER
!>          The leading dimension of the array R.  ldr >= m.
!> \endverbatim
!>
!> \param[in] c
!> \verbatim
!>          c is REAL array, dimension (min(m-1,n))
!>          The cosine parts of the Givens rotations.
!> \endverbatim
!>
!> \param[in] s
!> \verbatim
!>          s is COMPLEX array, dimension (min(m-1,n))
!>          The sine parts of the Givens rotations.
!> \endverbatim
!>
!> \ingroup qrdecomp
subroutine cqrqh(m,n,R,ldr,c,s)
  use iso_fortran_env
    integer, intent(in) :: m, n, ldr
    complex(real32), intent(inout) :: R(ldr,*)
    complex(real32), intent(in) :: s(*)
    real(real32), intent(in) :: c(*)
    external xerbla
    complex(real32) t
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
        call xerbla('CQRQH',info)
        return
    end if
    do i = 1,n
        ! apply stored rotations, column-wise
        ii = min(m-1,i)
        t = R(ii+1,i)
        do j = ii,1,-1
            R(j+1,i) = c(j)*t - conjg(s(j))*R(j,i)
            t = c(j)*R(j,i) + s(j)*t
        end do
        R(1,i) = t
    end do
end subroutine