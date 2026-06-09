! Copyright (C) 2009  VZLU Prague, a.s., Czech Republic
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
subroutine slup1up(m,n,L,ldl,R,ldr,p,u,v,w)
    ! purpose:      updates a row-pivoted LU factorization after rank-1 modi
    !               i.e., given an m-by-k lower-triangular matrix L with uni
    !               diagonal, a k-by-n upper-trapezoidal matrix R, and a
    !               permutation matrix P, where k = min(m,n),
    !               this subroutine updates L -> L1, R -> R1 and P -> P1 so
    !               L is again lower unit triangular, R upper trapezoidal,
    !               P permutation and P1'*L1*R1 = P'*L*R + u*v.'.
    !               (real version)
    ! arguments:
    ! m (in)        order of the matrix L.
    ! n (in)        number of columns of the matrix U.
    ! L (io)        on entry, the unit lower triangular matrix L.
    !               on exit, the updated matrix L1.
    ! ldl (in)      the leading dimension of L. ldl >= m.
    ! R (io)        on entry, the upper trapezoidal m-by-n matrix R.
    !               on exit, the updated matrix R1.
    ! ldr (in)      the leading dimension of R. ldr >= min(m,n).
    ! p (in)        the permutation vector representing P
    ! u (in)        the left m-vector.
    ! v (in)        the right n-vector.
    ! w (work)      a workspace vector of size m.
    !
    ! REMARK:       Algorithm is due to
    !               A. Kielbasinski, H. Schwetlick, Numerische Lineare
    !               Algebra, Verlag Harri Deutsch, 1988
    !
    integer m,n,ldl,ldr,p(*)
    real L(ldl,*),R(ldr,*),u(*),v(*),w(*)
    real one,tau,tmp
    parameter (one = 1e0, tau = 1e-1)
    integer k,info,i,j,itmp
    external xerbla,scopy,saxpy,strsv,sger,sgemv,sswap

    ! quick return if possible.
    k = min(m,n)
    if (k == 0) return
    ! check arguments.
    info = 0
    if (m < 0) then
        info = 1
    else if (n < 0) then
        info = 2
    else if (ldl < m) then
        info = 4
    else if (ldr < k) then
        info = 6
    endif
    if (info /= 0) then
        call xerbla('SLU1UP',info)
        return
    end if

    ! form L \ P*u.
    do i = 1,m
        w(i) = u(p(i))
    end do
    call strsv('L','N','U',k,L,ldl,w,1)
    ! if m > k = n, subtract the trailing part.
    if (m > k) then
        call sgemv('N',m-k,k,-one,L(k+1,1),ldl,w,1,one,w(k+1),1)
    end if

    ! work from bottom to top
    do j = k-1,1,-1
        if (abs(w(j)) < tau * abs(L(j+1,j)*w(j) + w(j+1))) then
            ! need pivoting. swap j and j+1
            tmp = w(j)
            w(j) = w(j+1)
            w(j+1) = tmp
            ! update p
            itmp = p(j)
            p(j) = p(j+1)
            p(j+1) = itmp
            ! update L
            call sswap(m-j+1,L(j,j),1,L(j,j+1),1)
            call sswap(j+1,L(j,1),ldl,L(j+1,1),ldl)
            ! update R
            call sswap(n-j+1,R(j,j),ldr,R(j+1,j),ldr)
            ! make L lower triangular again
            tmp = -L(j,j+1)
            call saxpy(m-j+1,tmp,L(j,j),1,L(j,j+1),1)
            ! update R
            call saxpy(n-j+1,-tmp,R(j+1,j),ldr,R(j,j),ldr)
            ! update w
            w(j) = w(j) - tmp*w(j+1)
        end if
        ! eliminate w(j+1)
        tmp = w(j+1)/w(j)
        w(j+1) = 0
        ! update R.
        call saxpy(n-j+1,-tmp,R(j,j),ldr,R(j+1,j),ldr)
        ! update L.
        call saxpy(m-j,tmp,L(j+1,j+1),1,L(j+1,j),1)
    end do

    ! add a multiple of v to R
    call saxpy(n,w(1),v,1,R(1,1),ldr)

    ! forward sweep
    do j = 1,k-1
        if (abs(R(j,j)) < tau * abs(L(j+1,j)*R(j,j) + R(j+1,j))) then
            ! need pivoting. swap j and j+1
            ! update p
            itmp = p(j)
            p(j) = p(j+1)
            p(j+1) = itmp
            ! update L
            call sswap(m-j+1,L(j,j),1,L(j,j+1),1)
            call sswap(j+1,L(j,1),ldl,L(j+1,1),ldl)
            ! update R
            call sswap(n-j+1,R(j,j),ldr,R(j+1,j),ldr)
            ! make L lower triangular again
            tmp = -L(j,j+1)
            call saxpy(m-j+1,tmp,L(j,j),1,L(j,j+1),1)
            ! update R
            call saxpy(n-j+1,-tmp,R(j+1,j),ldr,R(j,j),ldr)
        end if
        ! eliminate R(j+1,j)
        tmp = R(j+1,j)/R(j,j)
        ! update R.
        R(j+1,j) = 0e0
        call saxpy(n-j,-tmp,R(j,j+1),ldr,R(j+1,j+1),ldr)
        ! update L.
        call saxpy(m-j,tmp,L(j+1,j+1),1,L(j+1,j),1)
    end do

    ! if m > k = n, complete the update by updating the lower part of L.
    if (m > k) then
        call scopy(k,v,1,w,1)
        call strsv('U','T','N',k,R,ldr,w,1)
        call sger(m-k,k,one,w(k+1),1,w,1,L(k+1,1),ldl)
    endif
end subroutine
