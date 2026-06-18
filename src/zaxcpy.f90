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
!> \brief Performs scaled conjugate vector addition.
!>
!> \par Definition:
! =============
!> \verbatim
!>       subroutine zaxcpy(n,a,x,incx,y,incy)
!>
!>       .. Scalar Arguments ..
!>       integer            n, incx, incy
!>       ..
!>       .. Array Arguments ..
!>       double complex     a, x(*), y(*)
!>       ..
!> \endverbatim
!>
!> \par Purpose:
! =============
!> \verbatim
!>
!> ZAXCPY performs the operation y := y + a * conjg(x), where a is a
!> double complex scalar, x is a vector of length n, conjg(x) denotes
!> the element-wise complex conjugate of x, and y is a vector of the
!> same length.  On entry, y contains the existing values; on exit, y
!> is overwritten with the result.  This is the double precision
!> complex analogue of the BLAS zaxpy, with the x argument conjugated
!> before scaling.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is INTEGER
!>          The number of elements in vectors x and y.  If n <= 0,
!>          the subroutine returns immediately without modification.
!> \endverbatim
!>
!> \param[in] a
!> \verbatim
!>          a is COMPLEX*16
!>          The double complex scalar used to scale the conjugated
!>          vector conjg(x) before accumulation into y.
!> \endverbatim
!>
!> \param[in] x
!> \verbatim
!>          x is COMPLEX*16 array, dimension (*)
!>          The vector of length n whose complex conjugate is scaled
!>          by a and added to y.  x is not modified.
!> \endverbatim
!>
!> \param[in] incx
!> \verbatim
!>          incx is INTEGER
!>          The stride (increment) for elements of x.  If incx > 0,
!>          elements are accessed starting from x(1); if incx < 0,
!>          elements are accessed starting from
!>          x(1 + (-n+1)*incx).  A value of 1 accesses
!>          contiguous elements.
!> \endverbatim
!>
!> \param[in,out] y
!> \verbatim
!>          y is COMPLEX*16 array, dimension (*)
!>          On entry, the vector y of length n.  On exit, y is
!>          overwritten with y + a * conjg(x).
!> \endverbatim
!>
!> \param[in] incy
!> \verbatim
!>          incy is INTEGER
!>          The stride (increment) for elements of y.  If incy > 0,
!>          elements are accessed starting from y(1); if incy < 0,
!>          elements are accessed starting from
!>          y(1 + (-n+1)*incy).  A value of 1 accesses
!>          contiguous elements.
!> \endverbatim
!> \ingroup aux
subroutine zaxcpy(n,a,x,incx,y,incy)
    integer n,incx,incy
    double complex a,x(*),y(*)
    integer i,ix,iy
    ! quick return if possible.
    if (n <= 0) return
    if (incx /= 1 .or. incy /= 1) then
        ! code for unequal increments or equal increments not equal to 1
        ix = 1
        iy = 1
        if (incx.lt.0) ix = (-n+1)*incx + 1
        if (incy.lt.0) iy = (-n+1)*incy + 1
        do i = 1,n
            y(iy) = y(iy) + a*conjg(x(ix))
            ix = ix + incx
            iy = iy + incy
        end do
    else
        ! code for both increments equal to 1
        do i = 1,n
            y(i) = y(i) + a*conjg(x(i))
        end do
    end if
end subroutine
