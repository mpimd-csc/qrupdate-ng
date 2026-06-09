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
subroutine cch1dn(n,R,ldr,u,rw,info)
    ! purpose:      given an upper triangular matrix R that is a Cholesky
    !               factor of a hermitian positive definite matrix A, i.e.
    !               A = R'*R, this subroutine downdates R -> R1 so that
    !               R1'*R1 = A - u*u'
    !               (complex version)
    ! arguments:
    ! n (in)        the order of matrix R
    ! R (io)        on entry, the upper triangular matrix R
    !               on exit, the updated matrix R1
    ! ldr (in)      leading dimension of R. ldr >= n.
    ! u (io)        the vector determining the rank-1 update
    !               on exit, u contains the reflector sines
    !               used to transform R to R1.
    ! rw (out)      cosine parts of reflectors.
    !
    ! info (out)    on exit, error code:
    !                info = 0: success.
    !                info = 1: update violates positive-definiteness.
    !                info = 2: R is singular.
    !
    integer n,ldr
    complex R(ldr,*),u(*)
    real rw(*)
    integer info
    external ctrsv,clartg,scnrm2,xerbla
    complex crho,rr,ui,t
    real scnrm2,rho
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
        call xerbla('CCH1DN',-info)
        return
    end if

    ! check for singularity of R.
    do i = 1,n
        if (R(i,i) == 0e0) goto 20
    end do
    ! form R' \ u
    call ctrsv('U','C','N',n,R,ldr,u,1)
    rho = scnrm2(n,u,1)
    ! check positive definiteness
    rho = 1 - rho**2
    if (rho <= 0e0) goto 10
    crho = sqrt(rho)
    ! eliminate R' \ u
    do i = n,1,-1
        ui = u(i)
        ! generate next rotation
        call clartg(crho,ui,rw(i),u(i),rr)
        crho = rr
    end do
    ! apply rotations
    do i = n,1,-1
        ui = 0e0
        do j = i,1,-1
            t = rw(j)*ui + u(j)*R(j,i)
            R(j,i) = rw(j)*R(j,i) - conjg(u(j))*ui
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
