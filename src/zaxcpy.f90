! Copyright (C) 2008, 2009  VZLU Prague, a.s., Czech Republic
!
! Author: Jaroslav Hajek <highegg@gmail.com>
!
! This file is part of qrupdate.
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
subroutine zaxcpy(n,a,x,incx,y,incy)
    ! purpose:      constant times a conjugated vector plus a vector.
    ! arguments:
    ! n (in)        vector length
    ! a (in)        complex factor
    ! x (in)        added vector
    ! incx (in)     x increments
    ! y (io)        accumulator vector
    ! incy (in)     y increments
    !
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
