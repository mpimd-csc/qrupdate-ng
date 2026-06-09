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
subroutine clu1up(m,n,L,ldl,R,ldr,u,v)
    ! purpose:      updates an LU factorization after rank-1 modification
    !               i.e., given an m-by-k lower-triangular matrix L with uni
    !               diagonal and a k-by-n upper-trapezoidal matrix R,
    !               where k = min(m,n),
    !               this subroutine updates L -> L1 and R -> R1 so that
    !               L is again lower unit triangular, R upper trapezoidal,
    !               and L1*R1 = L*R + u*v.'.
    !               (complex version)
    ! arguments:
    ! m (in)        order of the matrix L.
    ! n (in)        number of columns of the matrix U.
    ! L (io)        on entry, the unit lower triangular matrix L.
    !               on exit, the updated matrix L1.
    ! ldl (in)      the leading dimension of L. ldl >= m.
    ! R (io)        on entry, the upper trapezoidal m-by-n matrix R.
    !               on exit, the updated matrix R1.
    ! ldr (in)      the leading dimension of R. ldr >= min(m,n).
    ! u (io)        the left m-vector. On exit, if k < m, u is destroyed.
    ! v (io)        the right n-vector. On exit, v is destroyed.
    !
    ! REMARK:       Algorithm is due to
    !               J. Bennett: Triangular factors of modified matrices,
    !                           Numerische Mathematik, 7 (1965)
    !
    integer m,n,ldl,ldr
    complex L(ldl,*),R(ldr,*),u(*),v(*)
    complex ui,vi
    integer k,info,i,j
    external xerbla
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
        call xerbla('CLU1UP',info)
        return
    end if

    ! The Bennett algorithm, modified for column-major access.
    ! The leading part.
    do i = 1,k
        ! prefetch
        ui = u(i)
        vi = v(i)
        ! delayed R update
        do j = 1,i-1
            R(j,i) = R(j,i) + u(j)*vi
            vi = vi - v(j)*R(j,i)
        end do
        ! diagonal update
        R(i,i) = R(i,i) + ui*vi
        vi = vi/R(i,i)
        ! L update
        do j = i+1,m
            u(j) = u(j) - ui*L(j,i)
            L(j,i) = L(j,i) + u(j)*vi
        end do
        u(i) = ui
        v(i) = vi
    end do

    ! Finish the trailing part of R if needed.
    do i = k+1,n
        vi = v(i)
        do j = 1,k
            R(j,i) = R(j,i) + u(j)*vi
            vi = vi - v(j)*R(j,i)
        end do
        v(i) = vi
    end do
end subroutine
