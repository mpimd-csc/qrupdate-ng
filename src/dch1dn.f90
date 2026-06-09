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
subroutine dch1dn(n,R,ldr,u,w,info)
    ! purpose:      given an upper triangular matrix R that is a Cholesky
    !               factor of a hermitian positive definite matrix A, i.e.
    !               A = R'*R, this subroutine downdates R -> R1 so that
    !               R1'*R1 = A - u*u'
    !               (real version)
    ! arguments:
    ! n (in)        the order of matrix R
    ! R (io)        on entry, the upper triangular matrix R
    !               on exit, the updated matrix R1
    ! ldr (in)      leading dimension of R. ldr >= n.
    ! u (io)        the vector determining the rank-1 update
    !               on exit, u contains the reflector sines
    !               used to transform R to R1.
    ! w (out)       cosine parts of reflectors.
    !
    ! info (out)    on exit, error code:
    !                info = 0: success.
    !                info = 1: update violates positive-definiteness.
    !                info = 2: R is singular.
    !
    integer n,ldr
    double precision R(ldr,*),u(*),w(*)
    integer info
    external xerbla,dtrsv,dlartg,dnrm2
    double precision dnrm2,rho,rr,ui,t
    integer i,j

    ! quick return if possible.
    if (n == 0) return
    ! check arguments.
    info = 0
    if (n < 0) then
        info = -1
    else if (ldr < n) then
        info = -3
    end if
    if (info /= 0) then
        call xerbla('DCH1DN',-info)
        return
    end if

    ! check for singularity of R.
    do i = 1,n
        if (R(i,i) == 0d0) goto 20
    end do
    ! form R' \ u
    call dtrsv('U','T','N',n,R,ldr,u,1)
    rho = dnrm2(n,u,1)
    ! check positive definiteness
    rho = 1 - rho**2
    if (rho <= 0d0) goto 10
    rho = sqrt(rho)
    ! eliminate R' \ u
    do i = n,1,-1
        ui = u(i)
        ! generate next rotation
        call dlartg(rho,ui,w(i),u(i),rr)
        rho = rr
    end do
    ! apply rotations
    do i = n,1,-1
        ui = 0d0
        do j = i,1,-1
            t = w(j)*ui + u(j)*R(j,i)
            R(j,i) = w(j)*R(j,i) - u(j)*ui
            ui = t
        end do
    end do

    ! normal return
    return
    ! error returns
    10 info = 1
    return
    20 info = 2
    return
end subroutine
