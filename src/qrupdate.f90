module qrupdate
    use iso_fortran_env
    implicit none

    interface
        subroutine caxcpy(n, a, x, incx, y, incy)
            use iso_fortran_env
            integer, intent(in) :: n
            complex(real32), intent(in) :: a
            complex(real32), intent(in) :: x
            integer, intent(in) :: incx
            complex(real32), intent(inout) :: y
            integer, intent(in) :: incy
        end subroutine caxcpy
    end interface

    interface
        subroutine cch1dn(n, R, ldr, u, rw, info)
            use iso_fortran_env
            integer, intent(in) :: n
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            complex(real32), intent(inout) :: u
            real(real32), intent(out) :: rw
            integer, intent(out) :: info
        end subroutine cch1dn
    end interface

    interface
        subroutine cch1up(n, R, ldr, u, w)
            use iso_fortran_env
            integer, intent(in) :: n
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            complex(real32), intent(inout) :: u
            real(real32), intent(out) :: w
        end subroutine cch1up
    end interface

    interface
        subroutine cchdex(n, R, ldr, j, rw)
            use iso_fortran_env
            integer, intent(in) :: n
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: j
            real(real32), intent(out) :: rw
        end subroutine cchdex
    end interface

    interface
        subroutine cchinx(n, R, ldr, j, u, rw, info)
            use iso_fortran_env
            integer, intent(in) :: n
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: j
            complex(real32), intent(inout) :: u
            real(real32), intent(out) :: rw
            integer, intent(out) :: info
        end subroutine cchinx
    end interface

    interface
        subroutine cchshx(n, R, ldr, i, j, w, rw)
            use iso_fortran_env
            integer, intent(in) :: n
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: i
            integer, intent(in) :: j
            complex(real32), intent(out) :: w
            real(real32), intent(out) :: rw
        end subroutine cchshx
    end interface

    interface
        subroutine cgqvec(m, n, Q, ldq, u)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            complex(real32), intent(in) :: Q
            integer, intent(in) :: ldq
            complex(real32), intent(out) :: u
        end subroutine cgqvec
    end interface

    interface
        subroutine clu1up(m, n, L, ldl, R, ldr, u, v)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            complex(real32), intent(inout) :: L
            integer, intent(in) :: ldl
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            complex(real32), intent(inout) :: u
            complex(real32), intent(inout) :: v
        end subroutine clu1up
    end interface

    interface
        subroutine clup1up(m, n, L, ldl, R, ldr, p, u, v, w)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            complex(real32), intent(inout) :: L
            integer, intent(in) :: ldl
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: p
            complex(real32), intent(in) :: u
            complex(real32), intent(in) :: v
            complex(real32), intent(out) :: w
        end subroutine clup1up
    end interface

    interface
        subroutine cqhqr(m, n, R, ldr, c, s)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            real(real32), intent(out) :: c
            complex(real32), intent(out) :: s
        end subroutine cqhqr
    end interface

    interface
        subroutine cqr1up(m, n, k, Q, ldq, R, ldr, u, v, w, rw)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            integer, intent(in) :: k
            complex(real32), intent(inout) :: Q
            integer, intent(in) :: ldq
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            complex(real32), intent(inout) :: u
            complex(real32), intent(inout) :: v
            complex(real32), intent(out) :: w
            real(real32), intent(out) :: rw
        end subroutine cqr1up
    end interface

    interface
        subroutine cqrdec(m, n, k, Q, ldq, R, ldr, j, rw)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            integer, intent(in) :: k
            complex(real32), intent(inout) :: Q
            integer, intent(in) :: ldq
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: j
            real(real32), intent(out) :: rw
        end subroutine cqrdec
    end interface

    interface
        subroutine cqrder(m, n, Q, ldq, R, ldr, j, w, rw)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            complex(real32), intent(inout) :: Q
            integer, intent(in) :: ldq
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: j
            complex(real32), intent(out) :: w
            real(real32), intent(out) :: rw
        end subroutine cqrder
    end interface

    interface
        subroutine cqrinc(m, n, k, Q, ldq, R, ldr, j, x, rw)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            integer, intent(in) :: k
            complex(real32), intent(inout) :: Q
            integer, intent(in) :: ldq
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: j
            complex(real32), intent(in) :: x
            real(real32), intent(out) :: rw
        end subroutine cqrinc
    end interface

    interface
        subroutine cqrinr(m, n, Q, ldq, R, ldr, j, x, rw)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            complex(real32), intent(inout) :: Q
            integer, intent(in) :: ldq
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: j
            complex(real32), intent(inout) :: x
            real(real32), intent(out) :: rw
        end subroutine cqrinr
    end interface

    interface
        subroutine cqrot(dir, m, n, Q, ldq, c, s)
            use iso_fortran_env
            character, intent(in) :: dir
            integer, intent(in) :: m
            integer, intent(in) :: n
            complex(real32), intent(inout) :: Q
            integer, intent(in) :: ldq
            real(real32), intent(in) :: c
            complex(real32), intent(in) :: s
        end subroutine cqrot
    end interface

    interface
        subroutine cqrqh(m, n, R, ldr, c, s)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            real(real32), intent(in) :: c
            complex(real32), intent(in) :: s
        end subroutine cqrqh
    end interface

    interface
        subroutine cqrshc(m, n, k, Q, ldq, R, ldr, i, j, w, rw)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            integer, intent(in) :: k
            complex(real32), intent(inout) :: Q
            integer, intent(in) :: ldq
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: i
            integer, intent(in) :: j
            complex(real32), intent(out) :: w
            real(real32), intent(out) :: rw
        end subroutine cqrshc
    end interface

    interface
        subroutine cqrtv1(n, u, w)
            use iso_fortran_env
            integer, intent(in) :: n
            complex(real32), intent(inout) :: u
            real(real32), intent(out) :: w
        end subroutine cqrtv1
    end interface

    interface
        subroutine dch1dn(n, R, ldr, u, w, info)
            use iso_fortran_env
            integer, intent(in) :: n
            real(real64), intent(inout) :: R
            integer, intent(in) :: ldr
            real(real64), intent(inout) :: u
            real(real64), intent(out) :: w
            integer, intent(out) :: info
        end subroutine dch1dn
    end interface

    interface
        subroutine dch1up(n, R, ldr, u, w)
            use iso_fortran_env
            integer, intent(in) :: n
            real(real64), intent(inout) :: R
            integer, intent(in) :: ldr
            real(real64), intent(inout) :: u
            real(real64), intent(out) :: w
        end subroutine dch1up
    end interface

    interface
        subroutine dchdex(n, R, ldr, j, w)
            use iso_fortran_env
            integer, intent(in) :: n
            real(real64), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: j
            real(real64), intent(out) :: w
        end subroutine dchdex
    end interface

    interface
        subroutine dchinx(n, R, ldr, j, u, w, info)
            use iso_fortran_env
            integer, intent(in) :: n
            real(real64), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: j
            real(real64), intent(inout) :: u
            real(real64), intent(out) :: w
            integer, intent(out) :: info
        end subroutine dchinx
    end interface

    interface
        subroutine dchshx(n, R, ldr, i, j, w)
            use iso_fortran_env
            integer, intent(in) :: n
            real(real64), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: i
            integer, intent(in) :: j
            real(real64), intent(out) :: w
        end subroutine dchshx
    end interface

    interface
        subroutine dgqvec(m, n, Q, ldq, u)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            real(real64), intent(in) :: Q
            integer, intent(in) :: ldq
            real(real64), intent(out) :: u
        end subroutine dgqvec
    end interface

    interface
        subroutine dlu1up(m, n, L, ldl, R, ldr, u, v)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            real(real64), intent(inout) :: L
            integer, intent(in) :: ldl
            real(real64), intent(inout) :: R
            integer, intent(in) :: ldr
            real(real64), intent(inout) :: u
            real(real64), intent(inout) :: v
        end subroutine dlu1up
    end interface

    interface
        subroutine dlup1up(m, n, L, ldl, R, ldr, p, u, v, w)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            real(real64), intent(inout) :: L
            integer, intent(in) :: ldl
            real(real64), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: p
            real(real64), intent(in) :: u
            real(real64), intent(in) :: v
            real(real64), intent(out) :: w
        end subroutine dlup1up
    end interface

    interface
        subroutine dqhqr(m, n, R, ldr, c, s)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            real(real64), intent(inout) :: R
            integer, intent(in) :: ldr
            real(real64), intent(out) :: c
            real(real64), intent(out) :: s
        end subroutine dqhqr
    end interface

    interface
        subroutine dqr1up(m, n, k, Q, ldq, R, ldr, u, v, w)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            integer, intent(in) :: k
            real(real64), intent(inout) :: Q
            integer, intent(in) :: ldq
            real(real64), intent(inout) :: R
            integer, intent(in) :: ldr
            real(real64), intent(inout) :: u
            real(real64), intent(inout) :: v
            real(real64), intent(out) :: w
        end subroutine dqr1up
    end interface

    interface
        subroutine dqrdec(m, n, k, Q, ldq, R, ldr, j, w)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            integer, intent(in) :: k
            real(real64), intent(inout) :: Q
            integer, intent(in) :: ldq
            real(real64), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: j
            real(real64), intent(out) :: w
        end subroutine dqrdec
    end interface

    interface
        subroutine dqrder(m, n, Q, ldq, R, ldr, j, w)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            real(real64), intent(inout) :: Q
            integer, intent(in) :: ldq
            real(real64), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: j
            real(real64), intent(out) :: w
        end subroutine dqrder
    end interface

    interface
        subroutine dqrinc(m, n, k, Q, ldq, R, ldr, j, x, w)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            integer, intent(in) :: k
            real(real64), intent(inout) :: Q
            integer, intent(in) :: ldq
            real(real64), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: j
            real(real64), intent(in) :: x
            real(real64), intent(out) :: w
        end subroutine dqrinc
    end interface

    interface
        subroutine dqrinr(m, n, Q, ldq, R, ldr, j, x, w)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            real(real64), intent(inout) :: Q
            integer, intent(in) :: ldq
            real(real64), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: j
            real(real64), intent(inout) :: x
            real(real64), intent(out) :: w
        end subroutine dqrinr
    end interface

    interface
        subroutine dqrot(dir, m, n, Q, ldq, c, s)
            use iso_fortran_env
            character, intent(in) :: dir
            integer, intent(in) :: m
            integer, intent(in) :: n
            real(real64), intent(inout) :: Q
            integer, intent(in) :: ldq
            real(real64), intent(in) :: c
            real(real64), intent(in) :: s
        end subroutine dqrot
    end interface

    interface
        subroutine dqrqh(m, n, R, ldr, c, s)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            real(real64), intent(inout) :: R
            integer, intent(in) :: ldr
            real(real64), intent(in) :: c
            real(real64), intent(in) :: s
        end subroutine dqrqh
    end interface

    interface
        subroutine dqrshc(m, n, k, Q, ldq, R, ldr, i, j, w)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            integer, intent(in) :: k
            real(real64), intent(inout) :: Q
            integer, intent(in) :: ldq
            real(real64), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: i
            integer, intent(in) :: j
            real(real64), intent(out) :: w
        end subroutine dqrshc
    end interface

    interface
        subroutine dqrtv1(n, u, w)
            use iso_fortran_env
            integer, intent(in) :: n
            real(real64), intent(inout) :: u
            real(real64), intent(out) :: w
        end subroutine dqrtv1
    end interface

    interface
        subroutine sch1dn(n, R, ldr, u, w, info)
            use iso_fortran_env
            integer, intent(in) :: n
            real(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            real(real32), intent(inout) :: u
            real(real32), intent(out) :: w
            integer, intent(out) :: info
        end subroutine sch1dn
    end interface

    interface
        subroutine sch1up(n, R, ldr, u, w)
            use iso_fortran_env
            integer, intent(in) :: n
            real(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            real(real32), intent(inout) :: u
            real(real32), intent(out) :: w
        end subroutine sch1up
    end interface

    interface
        subroutine schdex(n, R, ldr, j, w)
            use iso_fortran_env
            integer, intent(in) :: n
            real(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: j
            real(real32), intent(out) :: w
        end subroutine schdex
    end interface

    interface
        subroutine schinx(n, R, ldr, j, u, w, info)
            use iso_fortran_env
            integer, intent(in) :: n
            real(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: j
            real(real32), intent(inout) :: u
            real(real32), intent(out) :: w
            integer, intent(out) :: info
        end subroutine schinx
    end interface

    interface
        subroutine schshx(n, R, ldr, i, j, w)
            use iso_fortran_env
            integer, intent(in) :: n
            real(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: i
            integer, intent(in) :: j
            real(real32), intent(out) :: w
        end subroutine schshx
    end interface

    interface
        subroutine sgqvec(m, n, Q, ldq, u)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            real(real32), intent(in) :: Q
            integer, intent(in) :: ldq
            real(real32), intent(out) :: u
        end subroutine sgqvec
    end interface

    interface
        subroutine slu1up(m, n, L, ldl, R, ldr, u, v)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            real(real32), intent(inout) :: L
            integer, intent(in) :: ldl
            real(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            real(real32), intent(inout) :: u
            real(real32), intent(inout) :: v
        end subroutine slu1up
    end interface

    interface
        subroutine slup1up(m, n, L, ldl, R, ldr, p, u, v, w)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            real(real32), intent(inout) :: L
            integer, intent(in) :: ldl
            real(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: p
            real(real32), intent(in) :: u
            real(real32), intent(in) :: v
            real(real32), intent(out) :: w
        end subroutine slup1up
    end interface

    interface
        subroutine sqhqr(m, n, R, ldr, c, s)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            real(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            real(real32), intent(out) :: c
            real(real32), intent(out) :: s
        end subroutine sqhqr
    end interface

    interface
        subroutine sqr1up(m, n, k, Q, ldq, R, ldr, u, v, w)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            integer, intent(in) :: k
            real(real32), intent(inout) :: Q
            integer, intent(in) :: ldq
            real(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            real(real32), intent(inout) :: u
            real(real32), intent(inout) :: v
            real(real32), intent(out) :: w
        end subroutine sqr1up
    end interface

    interface
        subroutine sqrdec(m, n, k, Q, ldq, R, ldr, j, w)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            integer, intent(in) :: k
            real(real32), intent(inout) :: Q
            integer, intent(in) :: ldq
            real(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: j
            real(real32), intent(out) :: w
        end subroutine sqrdec
    end interface

    interface
        subroutine sqrder(m, n, Q, ldq, R, ldr, j, w)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            real(real32), intent(inout) :: Q
            integer, intent(in) :: ldq
            real(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: j
            real(real32), intent(out) :: w
        end subroutine sqrder
    end interface

    interface
        subroutine sqrinc(m, n, k, Q, ldq, R, ldr, j, x, w)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            integer, intent(in) :: k
            real(real32), intent(inout) :: Q
            integer, intent(in) :: ldq
            real(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: j
            real(real32), intent(in) :: x
            real(real32), intent(out) :: w
        end subroutine sqrinc
    end interface

    interface
        subroutine sqrinr(m, n, Q, ldq, R, ldr, j, x, w)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            real(real32), intent(inout) :: Q
            integer, intent(in) :: ldq
            real(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: j
            real(real32), intent(inout) :: x
            real(real32), intent(out) :: w
        end subroutine sqrinr
    end interface

    interface
        subroutine sqrot(dir, m, n, Q, ldq, c, s)
            use iso_fortran_env
            character, intent(in) :: dir
            integer, intent(in) :: m
            integer, intent(in) :: n
            real(real32), intent(inout) :: Q
            integer, intent(in) :: ldq
            real(real32), intent(in) :: c
            real(real32), intent(in) :: s
        end subroutine sqrot
    end interface

    interface
        subroutine sqrqh(m, n, R, ldr, c, s)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            real(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            real(real32), intent(in) :: c
            real(real32), intent(in) :: s
        end subroutine sqrqh
    end interface

    interface
        subroutine sqrshc(m, n, k, Q, ldq, R, ldr, i, j, w)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            integer, intent(in) :: k
            real(real32), intent(inout) :: Q
            integer, intent(in) :: ldq
            real(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: i
            integer, intent(in) :: j
            real(real32), intent(out) :: w
        end subroutine sqrshc
    end interface

    interface
        subroutine sqrtv1(n, u, w)
            use iso_fortran_env
            integer, intent(in) :: n
            real(real32), intent(inout) :: u
            real(real32), intent(out) :: w
        end subroutine sqrtv1
    end interface

    interface
        subroutine zaxcpy(n, a, x, incx, y, incy)
            use iso_fortran_env
            integer, intent(in) :: n
            complex(real32), intent(in) :: a
            complex(real32), intent(in) :: x
            integer, intent(in) :: incx
            complex(real32), intent(inout) :: y
            integer, intent(in) :: incy
        end subroutine zaxcpy
    end interface

    interface
        subroutine zch1dn(n, R, ldr, u, rw, info)
            use iso_fortran_env
            integer, intent(in) :: n
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            complex(real32), intent(inout) :: u
            real(real64), intent(out) :: rw
            integer, intent(out) :: info
        end subroutine zch1dn
    end interface

    interface
        subroutine zch1up(n, R, ldr, u, w)
            use iso_fortran_env
            integer, intent(in) :: n
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            complex(real32), intent(inout) :: u
            real(real64), intent(out) :: w
        end subroutine zch1up
    end interface

    interface
        subroutine zchdex(n, R, ldr, j, rw)
            use iso_fortran_env
            integer, intent(in) :: n
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: j
            real(real64), intent(out) :: rw
        end subroutine zchdex
    end interface

    interface
        subroutine zchinx(n, R, ldr, j, u, rw, info)
            use iso_fortran_env
            integer, intent(in) :: n
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: j
            complex(real32), intent(inout) :: u
            real(real64), intent(out) :: rw
            integer, intent(out) :: info
        end subroutine zchinx
    end interface

    interface
        subroutine zchshx(n, R, ldr, i, j, w, rw)
            use iso_fortran_env
            integer, intent(in) :: n
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: i
            integer, intent(in) :: j
            complex(real32), intent(out) :: w
            real(real64), intent(out) :: rw
        end subroutine zchshx
    end interface

    interface
        subroutine zgqvec(m, n, Q, ldq, u)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            complex(real32), intent(in) :: Q
            integer, intent(in) :: ldq
            complex(real32), intent(out) :: u
        end subroutine zgqvec
    end interface

    interface
        subroutine zlu1up(m, n, L, ldl, R, ldr, u, v)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            complex(real32), intent(inout) :: L
            integer, intent(in) :: ldl
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            complex(real32), intent(inout) :: u
            complex(real32), intent(inout) :: v
        end subroutine zlu1up
    end interface

    interface
        subroutine zlup1up(m, n, L, ldl, R, ldr, p, u, v, w)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            complex(real32), intent(inout) :: L
            integer, intent(in) :: ldl
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: p
            complex(real32), intent(in) :: u
            complex(real32), intent(in) :: v
            complex(real32), intent(out) :: w
        end subroutine zlup1up
    end interface

    interface
        subroutine zqhqr(m, n, R, ldr, c, s)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            real(real64), intent(out) :: c
            complex(real32), intent(out) :: s
        end subroutine zqhqr
    end interface

    interface
        subroutine zqr1up(m, n, k, Q, ldq, R, ldr, u, v, w, rw)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            integer, intent(in) :: k
            complex(real32), intent(inout) :: Q
            integer, intent(in) :: ldq
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            complex(real32), intent(inout) :: u
            complex(real32), intent(inout) :: v
            complex(real32), intent(out) :: w
            real(real64), intent(out) :: rw
        end subroutine zqr1up
    end interface

    interface
        subroutine zqrdec(m, n, k, Q, ldq, R, ldr, j, rw)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            integer, intent(in) :: k
            complex(real32), intent(inout) :: Q
            integer, intent(in) :: ldq
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: j
            real(real64), intent(out) :: rw
        end subroutine zqrdec
    end interface

    interface
        subroutine zqrder(m, n, Q, ldq, R, ldr, j, w, rw)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            complex(real32), intent(inout) :: Q
            integer, intent(in) :: ldq
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: j
            complex(real32), intent(out) :: w
            real(real64), intent(out) :: rw
        end subroutine zqrder
    end interface

    interface
        subroutine zqrinc(m, n, k, Q, ldq, R, ldr, j, x, rw)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            integer, intent(in) :: k
            complex(real32), intent(inout) :: Q
            integer, intent(in) :: ldq
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: j
            complex(real32), intent(in) :: x
            real(real64), intent(out) :: rw
        end subroutine zqrinc
    end interface

    interface
        subroutine zqrinr(m, n, Q, ldq, R, ldr, j, x, rw)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            complex(real32), intent(inout) :: Q
            integer, intent(in) :: ldq
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: j
            complex(real32), intent(inout) :: x
            real(real64), intent(out) :: rw
        end subroutine zqrinr
    end interface

    interface
        subroutine zqrot(dir, m, n, Q, ldq, c, s)
            use iso_fortran_env
            character, intent(in) :: dir
            integer, intent(in) :: m
            integer, intent(in) :: n
            complex(real32), intent(inout) :: Q
            integer, intent(in) :: ldq
            real(real64), intent(in) :: c
            complex(real32), intent(in) :: s
        end subroutine zqrot
    end interface

    interface
        subroutine zqrqh(m, n, R, ldr, c, s)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            real(real64), intent(in) :: c
            complex(real32), intent(in) :: s
        end subroutine zqrqh
    end interface

    interface
        subroutine zqrshc(m, n, k, Q, ldq, R, ldr, i, j, w, rw)
            use iso_fortran_env
            integer, intent(in) :: m
            integer, intent(in) :: n
            integer, intent(in) :: k
            complex(real32), intent(inout) :: Q
            integer, intent(in) :: ldq
            complex(real32), intent(inout) :: R
            integer, intent(in) :: ldr
            integer, intent(in) :: i
            integer, intent(in) :: j
            complex(real32), intent(out) :: w
            real(real64), intent(out) :: rw
        end subroutine zqrshc
    end interface

    interface
        subroutine zqrtv1(n, u, w)
            use iso_fortran_env
            integer, intent(in) :: n
            complex(real32), intent(inout) :: u
            real(real64), intent(out) :: w
        end subroutine zqrtv1
    end interface

end module qrupdate
