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
subroutine sch1up(n,R,ldr,u,w)
    ! purpose:      given an upper triangular matrix R that is a Cholesky
    !               factor of a symmetric positive definite matrix A, i.e.
    !               A = R'*R, this subroutine updates R -> R1 so that
    !               R1'*R1 = A + u*u'
    !               (real version)
    ! arguments:
    ! n (in)        the order of matrix R
    ! R (io)        on entry, the upper triangular matrix R
    !               on exit, the updated matrix R1
    ! ldr (in)      leading dimension of R. ldr >= n.
    ! u (io)        the vector determining the rank-1 update
    !               on exit, u contains the rotation sines
    !               used to transform R to R1.
    ! w (out)       cosine parts of rotations.
    !
    integer n,ldr
    real R(ldr,*),u(*)
    real w(*)
    external slartg
    real rr,ui,t
    integer i,j

    do i = 1,n
        ! apply stored rotations, column-wise
        ui = u(i)
        do j = 1,i-1
            t = w(j)*R(j,i) + u(j)*ui
            ui = w(j)*ui - u(j)*R(j,i)
            R(j,i) = t
        end do
        ! generate next rotation
        call slartg(R(i,i),ui,w(i),u(i),rr)
        R(i,i) = rr
    end do
end subroutine
