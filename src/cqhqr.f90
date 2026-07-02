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
!> \brief Reduces an upper Hessenberg matrix to upper trapezoidal form.
!>
!> \par Definition:
! =============
!> \verbatim
!>       subroutine cqhqr(m,n,R,ldr,c,s)
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
!> CQHQR reduces an m-by-n upper Hessenberg matrix R to upper
!> trapezoidal form.  Given an m-by-n upper Hessenberg matrix R,
!> CQHQR applies min(m-1,n) Givens rotations from the
!> left to eliminate the subdiagonal elements, producing an upper
!> trapezoidal matrix.
!>
!> On exit, c contains the cosine parts and s contains the sine
!> parts of the Givens rotations used in the reduction.
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
!>          On entry, the upper Hessenberg matrix R.  On exit, the
!>          updated upper trapezoidal matrix.
!> \endverbatim
!>
!> \param[in] ldr
!> \verbatim
!>          ldr is INTEGER
!>          The leading dimension of the array R.  ldr >= m.
!> \endverbatim
!>
!> \param[out] c
!> \verbatim
!>          c is REAL array, dimension (min(m-1,n))
!>          On exit, the cosine parts of the Givens rotations used
!>          to reduce R to upper trapezoidal form.
!> \endverbatim
!>
!> \param[out] s
!> \verbatim
!>          s is COMPLEX array, dimension (min(m-1,n))
!>          On exit, the sine parts of the Givens rotations used
!>          to reduce R to upper trapezoidal form.
!> \endverbatim
!>
!> \ingroup qrdecomp
subroutine cqhqr(m,n,R,ldr,c,s)
  use iso_fortran_env
    integer, intent(in) :: m, n, ldr
    complex(real32), intent(inout) :: R(ldr,*)
    real(real32), intent(out) :: c(*)
    complex(real32), intent(out) :: s(*)
    external xerbla,clartg
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
        call xerbla('CQHQR',info)
        return
    end if
    do i = 1,n
        ! apply stored rotations, column-wise
        t = R(1,i)
        ii = min(m,i)
        do j = 1,ii-1
            R(j,i) = c(j)*t + s(j)*R(j+1,i)
            t = c(j)*R(j+1,i) - conjg(s(j))*t
        end do
        if (ii < m) then
            ! generate next rotation
            call clartg(t,R(ii+1,i),c(i),s(i),R(ii,i))
            R(ii+1,i) = 0e0
        else
            R(ii,i) = t
        end if
    end do
end subroutine