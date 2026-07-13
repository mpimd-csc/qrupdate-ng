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
! Furthermore, I provides the interface to LAPACK's xerbla.
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

    interface
        subroutine xerbla( srname, info )
            character*(*), intent(in) :: srname
            integer, intent(in) ::  info
        end subroutine
    end interface

contains

    function lsame( ca, cb )
        character, intent(in) :: ca, cb
        logical :: lsame
        integer :: inta, intb, zcode

        lsame = ca == cb
        if ( lsame ) return

        zcode = ichar( 'Z' )
        inta = ichar( ca )
        intb = ichar( cb )

        if ( zcode == 90 .or. zcode == 122 ) then
            ! ASCII
            if ( inta >= 97 .and. inta <= 122 ) inta = inta - 32
            if ( intb >= 97 .and. intb <= 122 ) intb = intb - 32
        else if ( zcode == 233 .or. zcode == 169 ) then
            ! EBCDIC
            if ( ( inta >= 129 .and. inta <= 137 ) .or. &
                 ( inta >= 145 .and. inta <= 153 ) .or. &
                 ( inta >= 162 .and. inta <= 169 ) ) inta = inta + 64
            if ( ( intb >= 129 .and. intb <= 137 ) .or. &
                 ( intb >= 145 .and. intb <= 153 ) .or. &
                 ( intb >= 162 .and. intb <= 169 ) ) intb = intb + 64
        else if ( zcode == 218 .or. zcode == 250 ) then
            ! ASCII on Prime machines
            if ( inta >= 225 .and. inta <= 250 ) inta = inta - 32
            if ( intb >= 225 .and. intb <= 250 ) intb = intb - 32
        end if

        lsame = inta == intb
    end function lsame

end module qrupdate_blas
