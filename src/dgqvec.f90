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
subroutine dgqvec(m,n,Q,ldq,u)
    ! purpose:      given an orthogonal m-by-n matrix Q, n < m, generates
    !               a vector u such that Q'*u = 0 and norm(u) = 1.
    ! arguments:
    ! m (in)        number of rows of matrix Q.
    ! n (in)        number of columns of matrix Q.
    ! Q (in)        the orthogonal matrix Q.
    ! ldq (in)      leading dimension of Q.
    ! u (out)       the generated vector.
    !
    integer m,n,ldq
    double precision Q(ldq,*),u(*)
    external ddot,daxpy,dnrm2,dscal,xerbla
    double precision ddot,dnrm2,r
    integer info,i,j
    ! quick return if possible.
    if (m == 0) return
    if (n == 0) then
        u(1) = 1d0
        do i = 2,m
            u(i) = 0d0
        end do
        return
    end if
    ! check arguments.
    info = 0
    if (m < 0) then
        info = 1
    else if (n < 0) then
        info = 2
    else if (ldq < m) then
        info = 4
    end if
    if (info /= 0) then
        call xerbla('DGQVEC',info)
        return
    end if

    j = 1
    r = 0d0
    do while ( r .eq. 0d0 )
        ! probe j-th canonical unit vector.
        do i = 1,m
            u(i) = 0d0
        end do
        u(j) = 1d0
        ! form u - Q*Q'*u
        do i = 1,n
            r = ddot(m,Q(1,i),1,u,1)
            call daxpy(m,-r,Q(1,i),1,u,1)
        end do
        r = dnrm2(m,u,1)
        if (r == 0d0) then
            j = j + 1
            if (j > n) then
                ! this is fatal, and in theory, it can't happen.
                stop 'fatal: impossible condition in DGQVEC'
            end if
        end if
    end do
    call dscal(m,1d0/r,u,1)
end subroutine
