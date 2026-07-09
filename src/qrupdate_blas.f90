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

!
! This module contains the (c/z)dot(u/c) replacements to obtain a
! Fortran (gfortran/flang) ABI invariant implementation.
!
module qrupdate_blas
    use iso_fortran_env
    implicit none

    interface
        subroutine qrupdate_cdotc(ret,n,cx,incx,cy,incy)
            use iso_fortran_env
            integer, intent(in) :: incx, incy, n
            complex(real32), intent(in) :: cx(*),cy(*)
            complex(real32), intent(out) :: ret
        end subroutine qrupdate_cdotc
    end interface

    interface
        subroutine qrupdate_cdotu(ret,n,cx,incx,cy,incy)
            use iso_fortran_env
            integer, intent(in) :: incx, incy, n
            complex(real32), intent(in) :: cx(*),cy(*)
            complex(real32), intent(out) :: ret
        end subroutine qrupdate_cdotu
    end interface

    interface
        subroutine qrupdate_zdotc(ret,n,cx,incx,cy,incy)
            use iso_fortran_env
            integer, intent(in) :: incx, incy, n
            complex(real64), intent(in) :: cx(*),cy(*)
            complex(real64), intent(out) :: ret
        end subroutine qrupdate_zdotc
    end interface

    interface
        subroutine qrupdate_zdotu(ret,n,cx,incx,cy,incy)
            use iso_fortran_env
            integer, intent(in) :: incx, incy, n
            complex(real64), intent(in) :: cx(*),cy(*)
            complex(real64), intent(out) :: ret
        end subroutine qrupdate_zdotu
    end interface

end module qrupdate_blas
