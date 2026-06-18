! Copyright (C) 2008, 2009  VZLU Prague, a.s., Czech Republic, Jaroslav Hajek <highegg@gmail.com>
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
!> \brief Downdates a Cholesky factorization after a rank-1 modification.
!>
!> \par Definition:
! =============
!> \verbatim
!>       subroutine zch1dn(n,R,ldr,u,rw,info)
!>
!>       .. Scalar Arguments ..
!>       integer            n, ldr, info
!>       ..
!>       .. Array Arguments ..
!>       double complex     R(ldr,*), u(*)
!>       double precision   rw(*)
!>       ..
!> \endverbatim
!>
!> \par Purpose:
! =============
!> \verbatim
!>
!> ZCH1DN downdates the Cholesky factorization of a hermitian
!> positive definite matrix A after a rank-1 modification.  Given an
!> upper triangular matrix R that is a Cholesky factor of A, i.e.,
!> A = R'*R, where R' denotes the conjugate transpose of R, this
!> ZCH1DN downdates R -> R1 so that R1'*R1 = A - u*u', where u
!> is a given vector.
!>
!> The downdate is performed by applying a sequence of hyperbolic
!> rotations to restore the upper triangular structure of R.  On
!> exit, u contains the rotation sines and rw contains the rotation
!> cosines used in the transformation.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is INTEGER
!>          The order of matrix R.  n >= 0.
!> \endverbatim
!>
!> \param[in,out] R
!> \verbatim
!>          R is COMPLEX*16 array, dimension (ldr,n)
!>          On entry, the upper triangular matrix R, the Cholesky
!>          factor of A.  On exit, the updated upper triangular
!>          matrix R1, the Cholesky factor of A - u*u'.
!> \endverbatim
!>
!> \param[in] ldr
!> \verbatim
!>          ldr is INTEGER
!>          The leading dimension of the array R.  ldr >= n.
!> \endverbatim
!>
!> \param[in,out] u
!> \verbatim
!>          u is COMPLEX*16 array, dimension (n)
!>          On entry, the vector determining the rank-1 downdate.
!>          On exit, u contains the rotation sines used to
!>          transform R to R1.
!> \endverbatim
!>
!> \param[out] rw
!> \verbatim
!>          rw is DOUBLE PRECISION array, dimension (n)
!>          On exit, rw contains the cosine parts of the
!>          rotations used to transform R to R1.
!> \endverbatim
!>
!> \param[out] info
!> \verbatim
!>          info is INTEGER
!>          = 0:  successful exit
!>          = 1:  the update would violate positive-definiteness
!>          = 2:  R is singular
!> \endverbatim
!>
!> \ingroup choldecomp
subroutine zch1dn(n,R,ldr,u,rw,info)
    integer n,ldr
    double complex R(ldr,*),u(*)
    double precision rw(*)
    integer info
    external xerbla,ztrsv,zlartg,dznrm2
    double complex crho,rr,ui,t
    double precision dznrm2,rho
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
        call xerbla('ZCH1DN',-info)
        return
    end if

    ! check for singularity of R.
    do i = 1,n
        if (R(i,i) == 0d0) goto 20
    end do
    ! form R' \ u
    call ztrsv('U','C','N',n,R,ldr,u,1)
    rho = dznrm2(n,u,1)
    ! check positive definiteness
    rho = 1 - rho**2
    if (rho <= 0d0) goto 10
    crho = sqrt(rho)
    ! eliminate R' \ u
    do i = n,1,-1
        ui = u(i)
        ! generate next rotation
        call zlartg(crho,ui,rw(i),u(i),rr)
        crho = rr
    end do
    ! apply rotations
    do i = n,1,-1
        ui = 0d0
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