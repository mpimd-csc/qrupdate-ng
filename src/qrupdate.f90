!
! SPDX-License-Identifier: GPL-3.0-or-later
!
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
module qrupdate
    use iso_fortran_env
    implicit none

    interface
        subroutine caxcpy(n, a, x, incx, y, incy)
            use iso_fortran_env
            integer, intent(in)                     :: n
            complex(kind = real32), intent(in)      :: a
            complex(kind = real32), intent(in)      :: x(*)
            integer, intent(in)                     :: incx
            complex(kind = real32), intent(inout)   :: y(*)
            integer, intent(in)                     :: incy
        end subroutine caxcpy
    end interface
    public :: caxcpy

    interface
        subroutine cch1dn(n, r, ldr, u, rw, info)
            use iso_fortran_env
            integer, intent(in)                     :: n
            complex(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            complex(kind = real32), intent(inout)   :: u(*)
            real(kind = real32), intent(out)        :: rw(*)
            integer, intent(out)                    :: info
        end subroutine cch1dn
    end interface
    public :: cch1dn

    interface
        subroutine cch1up(n, r, ldr, u, w)
            use iso_fortran_env
            integer, intent(in)                     :: n
            complex(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            complex(kind = real32), intent(inout)   :: u(*)
            real(kind = real32), intent(out)        :: w(*)
        end subroutine cch1up
    end interface
    public :: cch1up

    interface
        subroutine cchdex(n, r, ldr, j, rw)
            use iso_fortran_env
            integer, intent(in)                     :: n
            complex(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            integer, intent(in)                     :: j
            real(kind = real32), intent(out)        :: rw(*)
        end subroutine cchdex
    end interface
    public :: cchdex

    interface
        subroutine cchinx(n, r, ldr, j, u, rw, info)
            use iso_fortran_env
            integer, intent(in)                     :: n
            complex(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            integer, intent(in)                     :: j
            complex(kind = real32), intent(inout)   :: u(*)
            real(kind = real32), intent(out)        :: rw(*)
            integer, intent(out)                    :: info
        end subroutine cchinx
    end interface
    public :: cchinx

    interface
        subroutine cchshx(n, r, ldr, i, j, w, rw)
            use iso_fortran_env
            integer, intent(in)                     :: n
            complex(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            integer, intent(in)                     :: i
            integer, intent(in)                     :: j
            complex(kind = real32), intent(out)     :: w(*)
            real(kind = real32), intent(out)        :: rw(*)
        end subroutine cchshx
    end interface
    public :: cchshx

    interface
        subroutine cgqvec(m, n, q, ldq, u)
            use iso_fortran_env
            integer, intent(in)                   :: m
            integer, intent(in)                   :: n
            complex(kind = real32), intent(in)    :: q(ldq, *)
            integer, intent(in)                   :: ldq
            complex(kind = real32), intent(out)   :: u(*)
        end subroutine cgqvec
    end interface
    public :: cgqvec

    interface
        subroutine clu1up(m, n, l, ldl, r, ldr, u, v)
            use iso_fortran_env
            integer, intent(in)                     :: m
            integer, intent(in)                     :: n
            complex(kind = real32), intent(inout)   :: l(ldl, *)
            integer, intent(in)                     :: ldl
            complex(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            complex(kind = real32), intent(inout)   :: u(*)
            complex(kind = real32), intent(inout)   :: v(*)
        end subroutine clu1up
    end interface
    public :: clu1up

    interface
        subroutine clup1up(m, n, l, ldl, r, ldr, p, u, &
                       v, w)
            use iso_fortran_env
            integer, intent(in)                     :: m
            integer, intent(in)                     :: n
            complex(kind = real32), intent(inout)   :: l(ldl, *)
            integer, intent(in)                     :: ldl
            complex(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            integer, intent(inout)                  :: p(*)
            complex(kind = real32), intent(in)      :: u(*)
            complex(kind = real32), intent(in)      :: v(*)
            complex(kind = real32), intent(out)     :: w(*)
        end subroutine clup1up
    end interface
    public :: clup1up

    interface
        subroutine cqhqr(m, n, r, ldr, c, s)
            use iso_fortran_env
            integer, intent(in)                     :: m
            integer, intent(in)                     :: n
            complex(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            real(kind = real32), intent(out)        :: c(*)
            complex(kind = real32), intent(out)     :: s(*)
        end subroutine cqhqr
    end interface
    public :: cqhqr

    interface
        subroutine cqr1up(m, n, k, q, ldq, r, ldr, u, &
                       v, w, rw)
            use iso_fortran_env
            integer, intent(in)                     :: m
            integer, intent(in)                     :: n
            integer, intent(in)                     :: k
            complex(kind = real32), intent(inout)   :: q(ldq, *)
            integer, intent(in)                     :: ldq
            complex(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            complex(kind = real32), intent(inout)   :: u(*)
            complex(kind = real32), intent(inout)   :: v(*)
            complex(kind = real32), intent(out)     :: w(*)
            real(kind = real32), intent(out)        :: rw(*)
        end subroutine cqr1up
    end interface
    public :: cqr1up

    interface
        subroutine cqrdec(m, n, k, q, ldq, r, ldr, j, &
                       rw)
            use iso_fortran_env
            integer, intent(in)                     :: m
            integer, intent(in)                     :: n
            integer, intent(in)                     :: k
            complex(kind = real32), intent(inout)   :: q(ldq, *)
            integer, intent(in)                     :: ldq
            complex(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            integer, intent(in)                     :: j
            real(kind = real32), intent(out)        :: rw(*)
        end subroutine cqrdec
    end interface
    public :: cqrdec

    interface
        subroutine cqrder(m, n, q, ldq, r, ldr, j, w, &
                       rw)
            use iso_fortran_env
            integer, intent(in)                     :: m
            integer, intent(in)                     :: n
            complex(kind = real32), intent(inout)   :: q(ldq, *)
            integer, intent(in)                     :: ldq
            complex(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            integer, intent(in)                     :: j
            complex(kind = real32), intent(out)     :: w(*)
            real(kind = real32), intent(out)        :: rw(*)
        end subroutine cqrder
    end interface
    public :: cqrder

    interface
        subroutine cqrinc(m, n, k, q, ldq, r, ldr, j, &
                       x, rw)
            use iso_fortran_env
            integer, intent(in)                     :: m
            integer, intent(in)                     :: n
            integer, intent(in)                     :: k
            complex(kind = real32), intent(inout)   :: q(ldq, *)
            integer, intent(in)                     :: ldq
            complex(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            integer, intent(in)                     :: j
            complex(kind = real32), intent(in)      :: x(*)
            real(kind = real32), intent(out)        :: rw(*)
        end subroutine cqrinc
    end interface
    public :: cqrinc

    interface
        subroutine cqrinr(m, n, q, ldq, r, ldr, j, x, &
                       rw)
            use iso_fortran_env
            integer, intent(in)                     :: m
            integer, intent(in)                     :: n
            complex(kind = real32), intent(inout)   :: q(ldq, *)
            integer, intent(in)                     :: ldq
            complex(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            integer, intent(in)                     :: j
            complex(kind = real32), intent(inout)   :: x(*)
            real(kind = real32), intent(out)        :: rw(*)
        end subroutine cqrinr
    end interface
    public :: cqrinr

    interface
        subroutine cqrot(dir, m, n, q, ldq, c, s)
            use iso_fortran_env
            character(len=*), intent(in)                   :: dir
            integer, intent(in)                     :: m
            integer, intent(in)                     :: n
            complex(kind = real32), intent(inout)   :: q(ldq, *)
            integer, intent(in)                     :: ldq
            real(kind = real32), intent(in)         :: c(*)
            complex(kind = real32), intent(in)      :: s(*)
        end subroutine cqrot
    end interface
    public :: cqrot

    interface
        subroutine cqrqh(m, n, r, ldr, c, s)
            use iso_fortran_env
            integer, intent(in)                     :: m
            integer, intent(in)                     :: n
            complex(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            real(kind = real32), intent(in)         :: c(*)
            complex(kind = real32), intent(in)      :: s(*)
        end subroutine cqrqh
    end interface
    public :: cqrqh

    interface
        subroutine cqrshc(m, n, k, q, ldq, r, ldr, i, &
                       j, w, rw)
            use iso_fortran_env
            integer, intent(in)                     :: m
            integer, intent(in)                     :: n
            integer, intent(in)                     :: k
            complex(kind = real32), intent(inout)   :: q(ldq, *)
            integer, intent(in)                     :: ldq
            complex(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            integer, intent(in)                     :: i
            integer, intent(in)                     :: j
            complex(kind = real32), intent(out)     :: w(*)
            real(kind = real32), intent(out)        :: rw(*)
        end subroutine cqrshc
    end interface
    public :: cqrshc

    interface
        subroutine cqrtv1(n, u, w)
            use iso_fortran_env
            integer, intent(in)                     :: n
            complex(kind = real32), intent(inout)   :: u(*)
            real(kind = real32), intent(out)        :: w(*)
        end subroutine cqrtv1
    end interface
    public :: cqrtv1

    interface
        subroutine dch1dn(n, r, ldr, u, w, info)
            use iso_fortran_env
            integer, intent(in)                  :: n
            real(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            real(kind = real64), intent(inout)   :: u(*)
            real(kind = real64), intent(out)     :: w(*)
            integer, intent(out)                 :: info
        end subroutine dch1dn
    end interface
    public :: dch1dn

    interface
        subroutine dch1up(n, r, ldr, u, w)
            use iso_fortran_env
            integer, intent(in)                  :: n
            real(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            real(kind = real64), intent(inout)   :: u(*)
            real(kind = real64), intent(out)     :: w(*)
        end subroutine dch1up
    end interface
    public :: dch1up

    interface
        subroutine dchdex(n, r, ldr, j, w)
            use iso_fortran_env
            integer, intent(in)                  :: n
            real(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            integer, intent(in)                  :: j
            real(kind = real64), intent(out)     :: w(*)
        end subroutine dchdex
    end interface
    public :: dchdex

    interface
        subroutine dchinx(n, r, ldr, j, u, w, info)
            use iso_fortran_env
            integer, intent(in)                  :: n
            real(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            integer, intent(in)                  :: j
            real(kind = real64), intent(inout)   :: u(*)
            real(kind = real64), intent(out)     :: w(*)
            integer, intent(out)                 :: info
        end subroutine dchinx
    end interface
    public :: dchinx

    interface
        subroutine dchshx(n, r, ldr, i, j, w)
            use iso_fortran_env
            integer, intent(in)                  :: n
            real(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            integer, intent(in)                  :: i
            integer, intent(in)                  :: j
            real(kind = real64), intent(out)     :: w(*)
        end subroutine dchshx
    end interface
    public :: dchshx

    interface
        subroutine dgqvec(m, n, q, ldq, u)
            use iso_fortran_env
            integer, intent(in)                :: m
            integer, intent(in)                :: n
            real(kind = real64), intent(in)    :: q(ldq, *)
            integer, intent(in)                :: ldq
            real(kind = real64), intent(out)   :: u(*)
        end subroutine dgqvec
    end interface
    public :: dgqvec

    interface
        subroutine dlu1up(m, n, l, ldl, r, ldr, u, v)
            use iso_fortran_env
            integer, intent(in)                  :: m
            integer, intent(in)                  :: n
            real(kind = real64), intent(inout)   :: l(ldl, *)
            integer, intent(in)                  :: ldl
            real(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            real(kind = real64), intent(inout)   :: u(*)
            real(kind = real64), intent(inout)   :: v(*)
        end subroutine dlu1up
    end interface
    public :: dlu1up

    interface
        subroutine dlup1up(m, n, l, ldl, r, ldr, p, u, &
                       v, w)
            use iso_fortran_env
            integer, intent(in)                  :: m
            integer, intent(in)                  :: n
            real(kind = real64), intent(inout)   :: l(ldl, *)
            integer, intent(in)                  :: ldl
            real(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            integer, intent(inout)               :: p(*)
            real(kind = real64), intent(in)      :: u(*)
            real(kind = real64), intent(in)      :: v(*)
            real(kind = real64), intent(out)     :: w(*)
        end subroutine dlup1up
    end interface
    public :: dlup1up

    interface
        subroutine dqhqr(m, n, r, ldr, c, s)
            use iso_fortran_env
            integer, intent(in)                  :: m
            integer, intent(in)                  :: n
            real(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            real(kind = real64), intent(out)     :: c(*)
            real(kind = real64), intent(out)     :: s(*)
        end subroutine dqhqr
    end interface
    public :: dqhqr

    interface
        subroutine dqr1up(m, n, k, q, ldq, r, ldr, u, &
                       v, w)
            use iso_fortran_env
            integer, intent(in)                  :: m
            integer, intent(in)                  :: n
            integer, intent(in)                  :: k
            real(kind = real64), intent(inout)   :: q(ldq, *)
            integer, intent(in)                  :: ldq
            real(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            real(kind = real64), intent(inout)   :: u(*)
            real(kind = real64), intent(inout)   :: v(*)
            real(kind = real64), intent(out)     :: w(*)
        end subroutine dqr1up
    end interface
    public :: dqr1up

    interface
        subroutine dqrdec(m, n, k, q, ldq, r, ldr, j, &
                       w)
            use iso_fortran_env
            integer, intent(in)                  :: m
            integer, intent(in)                  :: n
            integer, intent(in)                  :: k
            real(kind = real64), intent(inout)   :: q(ldq, *)
            integer, intent(in)                  :: ldq
            real(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            integer, intent(in)                  :: j
            real(kind = real64), intent(out)     :: w(*)
        end subroutine dqrdec
    end interface
    public :: dqrdec

    interface
        subroutine dqrder(m, n, q, ldq, r, ldr, j, w)
            use iso_fortran_env
            integer, intent(in)                  :: m
            integer, intent(in)                  :: n
            real(kind = real64), intent(inout)   :: q(ldq, *)
            integer, intent(in)                  :: ldq
            real(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            integer, intent(in)                  :: j
            real(kind = real64), intent(out)     :: w(*)
        end subroutine dqrder
    end interface
    public :: dqrder

    interface
        subroutine dqrinc(m, n, k, q, ldq, r, ldr, j, &
                       x, w)
            use iso_fortran_env
            integer, intent(in)                  :: m
            integer, intent(in)                  :: n
            integer, intent(in)                  :: k
            real(kind = real64), intent(inout)   :: q(ldq, *)
            integer, intent(in)                  :: ldq
            real(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            integer, intent(in)                  :: j
            real(kind = real64), intent(in)      :: x(*)
            real(kind = real64), intent(out)     :: w(*)
        end subroutine dqrinc
    end interface
    public :: dqrinc

    interface
        subroutine dqrinr(m, n, q, ldq, r, ldr, j, x, &
                       w)
            use iso_fortran_env
            integer, intent(in)                  :: m
            integer, intent(in)                  :: n
            real(kind = real64), intent(inout)   :: q(ldq, *)
            integer, intent(in)                  :: ldq
            real(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            integer, intent(in)                  :: j
            real(kind = real64), intent(inout)   :: x(*)
            real(kind = real64), intent(out)     :: w(*)
        end subroutine dqrinr
    end interface
    public :: dqrinr

    interface
        subroutine dqrot(dir, m, n, q, ldq, c, s)
            use iso_fortran_env
            character(len=*), intent(in)                :: dir
            integer, intent(in)                  :: m
            integer, intent(in)                  :: n
            real(kind = real64), intent(inout)   :: q(ldq, *)
            integer, intent(in)                  :: ldq
            real(kind = real64), intent(in)      :: c(*)
            real(kind = real64), intent(in)      :: s(*)
        end subroutine dqrot
    end interface
    public :: dqrot

    interface
        subroutine dqrqh(m, n, r, ldr, c, s)
            use iso_fortran_env
            integer, intent(in)                  :: m
            integer, intent(in)                  :: n
            real(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            real(kind = real64), intent(in)      :: c(*)
            real(kind = real64), intent(in)      :: s(*)
        end subroutine dqrqh
    end interface
    public :: dqrqh

    interface
        subroutine dqrshc(m, n, k, q, ldq, r, ldr, i, &
                       j, w)
            use iso_fortran_env
            integer, intent(in)                  :: m
            integer, intent(in)                  :: n
            integer, intent(in)                  :: k
            real(kind = real64), intent(inout)   :: q(ldq, *)
            integer, intent(in)                  :: ldq
            real(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            integer, intent(in)                  :: i
            integer, intent(in)                  :: j
            real(kind = real64), intent(out)     :: w(*)
        end subroutine dqrshc
    end interface
    public :: dqrshc

    interface
        subroutine dqrtv1(n, u, w)
            use iso_fortran_env
            integer, intent(in)                  :: n
            real(kind = real64), intent(inout)   :: u(*)
            real(kind = real64), intent(out)     :: w(*)
        end subroutine dqrtv1
    end interface
    public :: dqrtv1

    interface
        subroutine qrupdate_cdotc(ret, n, cx, incx, cy, incy)
            use iso_fortran_env
            complex(kind = real32), intent(out)   :: ret
            integer, intent(in)                   :: n
            complex(kind = real32), intent(in)    :: cx(*)
            integer, intent(in)                   :: incx
            complex(kind = real32), intent(in)    :: cy(*)
            integer, intent(in)                   :: incy
        end subroutine qrupdate_cdotc
    end interface
    public :: qrupdate_cdotc

    interface
        subroutine qrupdate_cdotu(ret, n, cx, incx, cy, incy)
            use iso_fortran_env
            complex(kind = real32), intent(out)   :: ret
            integer, intent(in)                   :: n
            complex(kind = real32), intent(in)    :: cx(*)
            integer, intent(in)                   :: incx
            complex(kind = real32), intent(in)    :: cy(*)
            integer, intent(in)                   :: incy
        end subroutine qrupdate_cdotu
    end interface
    public :: qrupdate_cdotu

    interface
        subroutine qrupdate_zdotc(ret, n, cx, incx, cy, incy)
            use iso_fortran_env
            complex(kind = real64), intent(out)   :: ret
            integer, intent(in)                   :: n
            complex(kind = real64), intent(in)    :: cx(*)
            integer, intent(in)                   :: incx
            complex(kind = real64), intent(in)    :: cy(*)
            integer, intent(in)                   :: incy
        end subroutine qrupdate_zdotc
    end interface
    public :: qrupdate_zdotc

    interface
        subroutine qrupdate_zdotu(ret, n, cx, incx, cy, incy)
            use iso_fortran_env
            complex(kind = real32), intent(out)   :: ret
            integer, intent(in)                   :: n
            complex(kind = real32), intent(in)    :: cx(*)
            integer, intent(in)                   :: incx
            complex(kind = real32), intent(in)    :: cy(*)
            integer, intent(in)                   :: incy
        end subroutine qrupdate_zdotu
    end interface
    public :: qrupdate_zdotu

    interface
        subroutine sch1dn(n, r, ldr, u, w, info)
            use iso_fortran_env
            integer, intent(in)                  :: n
            real(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            real(kind = real32), intent(inout)   :: u(*)
            real(kind = real32), intent(out)     :: w(*)
            integer, intent(out)                 :: info
        end subroutine sch1dn
    end interface
    public :: sch1dn

    interface
        subroutine sch1up(n, r, ldr, u, w)
            use iso_fortran_env
            integer, intent(in)                  :: n
            real(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            real(kind = real32), intent(inout)   :: u(*)
            real(kind = real32), intent(out)     :: w(*)
        end subroutine sch1up
    end interface
    public :: sch1up

    interface
        subroutine schdex(n, r, ldr, j, w)
            use iso_fortran_env
            integer, intent(in)                  :: n
            real(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            integer, intent(in)                  :: j
            real(kind = real32), intent(out)     :: w(*)
        end subroutine schdex
    end interface
    public :: schdex

    interface
        subroutine schinx(n, r, ldr, j, u, w, info)
            use iso_fortran_env
            integer, intent(in)                  :: n
            real(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            integer, intent(in)                  :: j
            real(kind = real32), intent(inout)   :: u(*)
            real(kind = real32), intent(out)     :: w(*)
            integer, intent(out)                 :: info
        end subroutine schinx
    end interface
    public :: schinx

    interface
        subroutine schshx(n, r, ldr, i, j, w)
            use iso_fortran_env
            integer, intent(in)                  :: n
            real(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            integer, intent(in)                  :: i
            integer, intent(in)                  :: j
            real(kind = real32), intent(out)     :: w(*)
        end subroutine schshx
    end interface
    public :: schshx

    interface
        subroutine sgqvec(m, n, q, ldq, u)
            use iso_fortran_env
            integer, intent(in)                :: m
            integer, intent(in)                :: n
            real(kind = real32), intent(in)    :: q(ldq, *)
            integer, intent(in)                :: ldq
            real(kind = real32), intent(out)   :: u(*)
        end subroutine sgqvec
    end interface
    public :: sgqvec

    interface
        subroutine slu1up(m, n, l, ldl, r, ldr, u, v)
            use iso_fortran_env
            integer, intent(in)                  :: m
            integer, intent(in)                  :: n
            real(kind = real32), intent(inout)   :: l(ldl, *)
            integer, intent(in)                  :: ldl
            real(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            real(kind = real32), intent(inout)   :: u(*)
            real(kind = real32), intent(inout)   :: v(*)
        end subroutine slu1up
    end interface
    public :: slu1up

    interface
        subroutine slup1up(m, n, l, ldl, r, ldr, p, u, &
                       v, w)
            use iso_fortran_env
            integer, intent(in)                  :: m
            integer, intent(in)                  :: n
            real(kind = real32), intent(inout)   :: l(ldl, *)
            integer, intent(in)                  :: ldl
            real(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            integer, intent(inout)               :: p(*)
            real(kind = real32), intent(in)      :: u(*)
            real(kind = real32), intent(in)      :: v(*)
            real(kind = real32), intent(out)     :: w(*)
        end subroutine slup1up
    end interface
    public :: slup1up

    interface
        subroutine sqhqr(m, n, r, ldr, c, s)
            use iso_fortran_env
            integer, intent(in)                  :: m
            integer, intent(in)                  :: n
            real(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            real(kind = real32), intent(in)      :: c(*)
            real(kind = real32), intent(in)      :: s(*)
        end subroutine sqhqr
    end interface
    public :: sqhqr

    interface
        subroutine sqr1up(m, n, k, q, ldq, r, ldr, u, &
                       v, w)
            use iso_fortran_env
            integer, intent(in)                  :: m
            integer, intent(in)                  :: n
            integer, intent(in)                  :: k
            real(kind = real32), intent(inout)   :: q(ldq, *)
            integer, intent(in)                  :: ldq
            real(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            real(kind = real32), intent(inout)   :: u(*)
            real(kind = real32), intent(inout)   :: v(*)
            real(kind = real32), intent(out)     :: w(*)
        end subroutine sqr1up
    end interface
    public :: sqr1up

    interface
        subroutine sqrdec(m, n, k, q, ldq, r, ldr, j, &
                       w)
            use iso_fortran_env
            integer, intent(in)                  :: m
            integer, intent(in)                  :: n
            integer, intent(in)                  :: k
            real(kind = real32), intent(inout)   :: q(ldq, *)
            integer, intent(in)                  :: ldq
            real(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            integer, intent(in)                  :: j
            real(kind = real32), intent(out)     :: w(*)
        end subroutine sqrdec
    end interface
    public :: sqrdec

    interface
        subroutine sqrder(m, n, q, ldq, r, ldr, j, w)
            use iso_fortran_env
            integer, intent(in)                  :: m
            integer, intent(in)                  :: n
            real(kind = real32), intent(inout)   :: q(ldq, *)
            integer, intent(in)                  :: ldq
            real(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            integer, intent(in)                  :: j
            real(kind = real32), intent(out)     :: w(*)
        end subroutine sqrder
    end interface
    public :: sqrder

    interface
        subroutine sqrinc(m, n, k, q, ldq, r, ldr, j, &
                       x, w)
            use iso_fortran_env
            integer, intent(in)                  :: m
            integer, intent(in)                  :: n
            integer, intent(in)                  :: k
            real(kind = real32), intent(inout)   :: q(ldq, *)
            integer, intent(in)                  :: ldq
            real(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            integer, intent(in)                  :: j
            real(kind = real32), intent(in)      :: x(*)
            real(kind = real32), intent(out)     :: w(*)
        end subroutine sqrinc
    end interface
    public :: sqrinc

    interface
        subroutine sqrinr(m, n, q, ldq, r, ldr, j, x, &
                       w)
            use iso_fortran_env
            integer, intent(in)                  :: m
            integer, intent(in)                  :: n
            real(kind = real32), intent(inout)   :: q(ldq, *)
            integer, intent(in)                  :: ldq
            real(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            integer, intent(in)                  :: j
            real(kind = real32), intent(inout)   :: x(*)
            real(kind = real32), intent(out)     :: w(*)
        end subroutine sqrinr
    end interface
    public :: sqrinr

    interface
        subroutine sqrot(dir, m, n, q, ldq, c, s)
            use iso_fortran_env
            character(len=*), intent(in)                :: dir
            integer, intent(in)                  :: m
            integer, intent(in)                  :: n
            real(kind = real32), intent(inout)   :: q(ldq, *)
            integer, intent(in)                  :: ldq
            real(kind = real32), intent(in)      :: c(*)
            real(kind = real32), intent(in)      :: s(*)
        end subroutine sqrot
    end interface
    public :: sqrot

    interface
        subroutine sqrqh(m, n, r, ldr, c, s)
            use iso_fortran_env
            integer, intent(in)                  :: m
            integer, intent(in)                  :: n
            real(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            real(kind = real32), intent(in)      :: c(*)
            real(kind = real32), intent(in)      :: s(*)
        end subroutine sqrqh
    end interface
    public :: sqrqh

    interface
        subroutine sqrshc(m, n, k, q, ldq, r, ldr, i, &
                       j, w)
            use iso_fortran_env
            integer, intent(in)                  :: m
            integer, intent(in)                  :: n
            integer, intent(in)                  :: k
            real(kind = real32), intent(inout)   :: q(ldq, *)
            integer, intent(in)                  :: ldq
            real(kind = real32), intent(inout)   :: r(ldr, *)
            integer, intent(in)                  :: ldr
            integer, intent(in)                  :: i
            integer, intent(in)                  :: j
            real(kind = real32), intent(out)     :: w(*)
        end subroutine sqrshc
    end interface
    public :: sqrshc

    interface
        subroutine sqrtv1(n, u, w)
            use iso_fortran_env
            integer, intent(in)                  :: n
            real(kind = real32), intent(inout)   :: u(*)
            real(kind = real32), intent(out)     :: w(*)
        end subroutine sqrtv1
    end interface
    public :: sqrtv1

    interface
        subroutine zaxcpy(n, a, x, incx, y, incy)
            use iso_fortran_env
            integer, intent(in)                     :: n
            complex(kind = real64), intent(in)      :: a
            complex(kind = real64), intent(in)      :: x(*)
            integer, intent(in)                     :: incx
            complex(kind = real64), intent(inout)   :: y(*)
            integer, intent(in)                     :: incy
        end subroutine zaxcpy
    end interface
    public :: zaxcpy

    interface
        subroutine zch1dn(n, r, ldr, u, rw, info)
            use iso_fortran_env
            integer, intent(in)                     :: n
            complex(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            complex(kind = real64), intent(inout)   :: u(*)
            real(kind = real64), intent(out)        :: rw(*)
            integer, intent(out)                    :: info
        end subroutine zch1dn
    end interface
    public :: zch1dn

    interface
        subroutine zch1up(n, r, ldr, u, w)
            use iso_fortran_env
            integer, intent(in)                     :: n
            complex(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            complex(kind = real64), intent(inout)   :: u(*)
            real(kind = real64), intent(out)        :: w(*)
        end subroutine zch1up
    end interface
    public :: zch1up

    interface
        subroutine zchdex(n, r, ldr, j, rw)
            use iso_fortran_env
            integer, intent(in)                     :: n
            complex(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            integer, intent(in)                     :: j
            real(kind = real64), intent(out)        :: rw(*)
        end subroutine zchdex
    end interface
    public :: zchdex

    interface
        subroutine zchinx(n, r, ldr, j, u, rw, info)
            use iso_fortran_env
            integer, intent(in)                     :: n
            complex(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            integer, intent(in)                     :: j
            complex(kind = real64), intent(inout)   :: u(*)
            real(kind = real64), intent(out)        :: rw(*)
            integer, intent(out)                    :: info
        end subroutine zchinx
    end interface
    public :: zchinx

    interface
        subroutine zchshx(n, r, ldr, i, j, w, rw)
            use iso_fortran_env
            integer, intent(in)                     :: n
            complex(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            integer, intent(in)                     :: i
            integer, intent(in)                     :: j
            complex(kind = real64), intent(out)     :: w(*)
            real(kind = real64), intent(out)        :: rw(*)
        end subroutine zchshx
    end interface
    public :: zchshx

    interface
        subroutine zgqvec(m, n, q, ldq, u)
            use iso_fortran_env
            integer, intent(in)                   :: m
            integer, intent(in)                   :: n
            complex(kind = real64), intent(in)    :: q(ldq, *)
            integer, intent(in)                   :: ldq
            complex(kind = real64), intent(out)   :: u(*)
        end subroutine zgqvec
    end interface
    public :: zgqvec

    interface
        subroutine zlu1up(m, n, l, ldl, r, ldr, u, v)
            use iso_fortran_env
            integer, intent(in)                     :: m
            integer, intent(in)                     :: n
            complex(kind = real64), intent(inout)   :: l(ldl, *)
            integer, intent(in)                     :: ldl
            complex(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            complex(kind = real64), intent(inout)   :: u(*)
            complex(kind = real64), intent(inout)   :: v(*)
        end subroutine zlu1up
    end interface
    public :: zlu1up

    interface
        subroutine zlup1up(m, n, l, ldl, r, ldr, p, u, &
                       v, w)
            use iso_fortran_env
            integer, intent(in)                     :: m
            integer, intent(in)                     :: n
            complex(kind = real64), intent(inout)   :: l(ldl, *)
            integer, intent(in)                     :: ldl
            complex(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            integer, intent(inout)                  :: p(*)
            complex(kind = real64), intent(in)      :: u(*)
            complex(kind = real64), intent(in)      :: v(*)
            complex(kind = real64), intent(out)     :: w(*)
        end subroutine zlup1up
    end interface
    public :: zlup1up

    interface
        subroutine zqhqr(m, n, r, ldr, c, s)
            use iso_fortran_env
            integer, intent(in)                     :: m
            integer, intent(in)                     :: n
            complex(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            real(kind = real64), intent(out)        :: c(*)
            complex(kind = real64), intent(out)     :: s(*)
        end subroutine zqhqr
    end interface
    public :: zqhqr

    interface
        subroutine zqr1up(m, n, k, q, ldq, r, ldr, u, &
                       v, w, rw)
            use iso_fortran_env
            integer, intent(in)                     :: m
            integer, intent(in)                     :: n
            integer, intent(in)                     :: k
            complex(kind = real64), intent(inout)   :: q(ldq, *)
            integer, intent(in)                     :: ldq
            complex(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            complex(kind = real64), intent(inout)   :: u(*)
            complex(kind = real64), intent(inout)   :: v(*)
            complex(kind = real64), intent(out)     :: w(*)
            real(kind = real64), intent(out)        :: rw(*)
        end subroutine zqr1up
    end interface
    public :: zqr1up

    interface
        subroutine zqrdec(m, n, k, q, ldq, r, ldr, j, &
                       rw)
            use iso_fortran_env
            integer, intent(in)                     :: m
            integer, intent(in)                     :: n
            integer, intent(in)                     :: k
            complex(kind = real64), intent(inout)   :: q(ldq, *)
            integer, intent(in)                     :: ldq
            complex(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            integer, intent(in)                     :: j
            real(kind = real64), intent(out)        :: rw(*)
        end subroutine zqrdec
    end interface
    public :: zqrdec

    interface
        subroutine zqrder(m, n, q, ldq, r, ldr, j, w, &
                       rw)
            use iso_fortran_env
            integer, intent(in)                     :: m
            integer, intent(in)                     :: n
            complex(kind = real64), intent(inout)   :: q(ldq, *)
            integer, intent(in)                     :: ldq
            complex(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            integer, intent(in)                     :: j
            complex(kind = real64), intent(out)     :: w(*)
            real(kind = real64), intent(out)        :: rw(*)
        end subroutine zqrder
    end interface
    public :: zqrder

    interface
        subroutine zqrinc(m, n, k, q, ldq, r, ldr, j, &
                       x, rw)
            use iso_fortran_env
            integer, intent(in)                     :: m
            integer, intent(in)                     :: n
            integer, intent(in)                     :: k
            complex(kind = real64), intent(inout)   :: q(ldq, *)
            integer, intent(in)                     :: ldq
            complex(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            integer, intent(in)                     :: j
            complex(kind = real64), intent(in)      :: x(*)
            real(kind = real64), intent(out)        :: rw(*)
        end subroutine zqrinc
    end interface
    public :: zqrinc

    interface
        subroutine zqrinr(m, n, q, ldq, r, ldr, j, x, &
                       rw)
            use iso_fortran_env
            integer, intent(in)                     :: m
            integer, intent(in)                     :: n
            complex(kind = real64), intent(inout)   :: q(ldq, *)
            integer, intent(in)                     :: ldq
            complex(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            integer, intent(in)                     :: j
            complex(kind = real64), intent(inout)   :: x(*)
            real(kind = real64), intent(out)        :: rw(*)
        end subroutine zqrinr
    end interface
    public :: zqrinr

    interface
        subroutine zqrot(dir, m, n, q, ldq, c, s)
            use iso_fortran_env
            character(len=*), intent(in)                   :: dir
            integer, intent(in)                     :: m
            integer, intent(in)                     :: n
            complex(kind = real64), intent(inout)   :: q(ldq, *)
            integer, intent(in)                     :: ldq
            real(kind = real64), intent(in)         :: c(*)
            complex(kind = real64), intent(in)      :: s(*)
        end subroutine zqrot
    end interface
    public :: zqrot

    interface
        subroutine zqrqh(m, n, r, ldr, c, s)
            use iso_fortran_env
            integer, intent(in)                     :: m
            integer, intent(in)                     :: n
            complex(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            real(kind = real64), intent(in)         :: c(*)
            complex(kind = real64), intent(in)      :: s(*)
        end subroutine zqrqh
    end interface
    public :: zqrqh

    interface
        subroutine zqrshc(m, n, k, q, ldq, r, ldr, i, &
                       j, w, rw)
            use iso_fortran_env
            integer, intent(in)                     :: m
            integer, intent(in)                     :: n
            integer, intent(in)                     :: k
            complex(kind = real64), intent(inout)   :: q(ldq, *)
            integer, intent(in)                     :: ldq
            complex(kind = real64), intent(inout)   :: r(ldr, *)
            integer, intent(in)                     :: ldr
            integer, intent(in)                     :: i
            integer, intent(in)                     :: j
            complex(kind = real64), intent(out)     :: w(*)
            real(kind = real64), intent(out)        :: rw(*)
        end subroutine zqrshc
    end interface
    public :: zqrshc

    interface
        subroutine zqrtv1(n, u, w)
            use iso_fortran_env
            integer, intent(in)                     :: n
            complex(kind = real64), intent(inout)   :: u(*)
            real(kind = real64), intent(out)        :: w(*)
        end subroutine zqrtv1
    end interface
    public :: zqrtv1

end module qrupdate
