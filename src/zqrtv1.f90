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
subroutine zqrtv1(n,u,w)
    ! purpose:      generates a sequence of n-1 Givens rotations that
    !               eliminate all but the first element of a vector u.
    ! arguments:
    ! n (in)        the length of the vector u
    ! u (io)        on entry, the vector u.
    !               on exit, u(2:n) contains the rotation sines, u(1)
    !               contains the remaining element.
    ! w (o)         on exit, w contains the rotation cosines.
    !
    integer n
    double complex u(*)
    double precision w(*)
    external zlartg
    double complex rr,t
    integer i
    ! quick return if possible.
    if (n <= 0) return
    rr = u(n)
    do i = n-1,1,-1
        call zlartg(u(i),rr,w(i),u(i+1),t)
        rr = t
    end do
    u(1) = rr
end subroutine
