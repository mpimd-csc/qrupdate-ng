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
subroutine sqhqr(m,n,R,ldr,c,s)
    ! purpose:      given an m-by-n upper Hessenberg matrix R, this
    !               subroutine updates R to upper trapezoidal form
    !               using min(m-1,n) Givens rotations.
    !               (real version)
    ! arguments:
    ! m (in)        number of rows of the matrix R
    ! n (in)        number of columns of the matrix R
    ! R (io)        on entry, the upper Hessenberg matrix R
    !               on exit, the updated upper trapezoidal matrix
    ! ldr (in)      leading dimension of R, >= m
    ! c(out)        rotation cosines, size at least min(m-1,n)
    ! s(out)        rotation sines, size at least min(m-1,n)
    !
    integer m,n,ldr
    real R(ldr,*),c(*),s(*)
    external xerbla,slartg
    real t
    integer info,i,ii,j
    ! quick return if possible.
    if (m == 0 .or. m == 1 .or. n == 0) return
    ! check arguments.
    info = 0
    if (m < 0) then
        info = 1
    else if (n < 0) then
        info = 2
    else if (ldr < m) then
        info = 4
    end if
    if (info /= 0) then
        call xerbla('SQHQR',info)
        return
    end if
    do i = 1,n
        ! apply stored rotations, column-wise
        t = R(1,i)
        ii = min(m,i)
        do j = 1,ii-1
            R(j,i) = c(j)*t + s(j)*R(j+1,i)
            t = c(j)*R(j+1,i) - s(j)*t
        end do
        if (ii < m) then
            ! generate next rotation
            call slartg(t,R(ii+1,i),c(i),s(i),R(ii,i))
            R(ii+1,i) = 0e0
        else
            R(ii,i) = t
        end if
    end do
end subroutine
