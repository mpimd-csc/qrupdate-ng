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
!> \brief ZDOTC Complex dot product (x**H y version)
!
!> \par  Definition:
!  ===========
!> \verbatim
!>       SUBROUTINE QRUPDATE_ZDOTC(RET, N,CX,INCX,CY,INCY)
!>
!>       .. Scalar Arguments ..
!>       INTEGER INCX,INCY,N
!>       COMPLEX(real64) ret
!>       ..
!>       .. Array Arguments ..
!>       COMPLEX(real64) CX(*),CY(*)
!>       ..
!> \endverbatim
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZDOTC forms the dot product of two complex vectors
!>      ZDOTC = X^H * Y
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[out] RET
!> \verbatim
!>         RET is COMPLEX(real64)
!>         RET contains the dot product.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>         N is INTEGER
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] CX
!> \verbatim
!>          CX is COMPLEX(real64) array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>         INCX is INTEGER
!>         storage spacing between elements of CX
!> \endverbatim
!>
!> \param[in] CY
!> \verbatim
!>          CY is COMPLEX(real64) array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>         INCY is INTEGER
!>         storage spacing between elements of CY
!> \endverbatim
!
!> \remark This routine is equivalent to zdotc from BLAS but invariant with the compiler's Fortran ABI
!
!> \ingroup aux
!  =====================================================================
subroutine qrupdate_zdotc(ret,n,cx,incx,cy,incy)
    use iso_fortran_env
    implicit none

    !     .. scalar arguments ..
    integer, intent(in) :: incx,incy,n
    !     ..
    !     .. array arguments ..
    complex(real64), intent(in) :: cx(*),cy(*)
    complex(real64), intent(out) :: ret
    !     ..
    !
    !     .. local scalars ..
    complex(real64) ctemp
    integer i,ix,iy
    !     ..
    ctemp = (0.0,0.0)
    ret = (0.0,0.0)
    if (n.le.0) return
    if (incx.eq.1 .and. incy.eq.1) then
        !
        !        code for both increments equal to 1
        !
        do i = 1,n
            ctemp = ctemp + conjg(cx(i))*cy(i)
        end do
    else
        !
        !        code for unequal increments or equal increments
        !          not equal to 1
        !
        ix = 1
        iy = 1
        if (incx.lt.0) ix = (-n+1)*incx + 1
        if (incy.lt.0) iy = (-n+1)*incy + 1
        do i = 1,n
            ctemp = ctemp + conjg(cx(ix))*cy(iy)
            ix = ix + incx
            iy = iy + incy
        end do
    end if
    ret = ctemp
    return
    !
    !     end of zdotc
    !
end subroutine


