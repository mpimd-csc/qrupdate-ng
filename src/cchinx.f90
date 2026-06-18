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
!> \brief Updates a Cholesky factorization after inserting a row and column.
!>
!> \par Definition:
! =============
!> \verbatim
!>       subroutine cchinx(n,R,ldr,j,u,rw,info)
!>
!>       .. Scalar Arguments ..
!>       integer            n, j, ldr, info
!>       ..
!>       .. Array Arguments ..
!>       complex            R(ldr,*), u(*), rw(*)
!>       ..
!> \endverbatim
!>
!> \par Purpose:
! =============
!> \verbatim
!>
!> CCHINX updates the Cholesky factorization of a hermitian
!> positive definite matrix A after inserting a row and column.
!> Given an upper triangular matrix R that is a Cholesky factor of
!> A, i.e., A = R'*R, where R' denotes the conjugate transpose of
!> R, CCHINX updates R -> R1 so that R1'*R1 = A1, where
!> A1(jj,jj) = A, A1(j,:) = u', A1(:,j) = u, and
!> jj = [1:j-1, j+1:n+1].
!>
!> On exit, u is destroyed and R is extended by one row and column.
!> The insertion is performed by first solving R'*u = v, checking
!> positive definiteness, and then retriangularizing.
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
!>          R is COMPLEX array, dimension (ldr,n+1)
!>          On entry, the upper triangular matrix R, the Cholesky
!>          factor of A.  On exit, the updated upper triangular
!>          matrix R1, the Cholesky factor of A1.
!> \endverbatim
!>
!> \param[in] ldr
!> \verbatim
!>          ldr is INTEGER
!>          The leading dimension of the array R.  ldr >= n+1.
!> \endverbatim
!>
!> \param[in] j
!> \verbatim
!>          j is INTEGER
!>          The position of the inserted row and column.
!>          1 <= j <= n+1.
!> \endverbatim
!>
!> \param[in,out] u
!> \verbatim
!>          u is COMPLEX array, dimension (n+1)
!>          On entry, the vector defining the inserted row/column.
!>          On exit, u is destroyed.
!> \endverbatim
!>
!> \param[out] rw
!> \verbatim
!>          rw is REAL array, dimension (n+1)
!>          Workspace vector used to store rotation cosines during
!>          the retriangularization.
!> \endverbatim
!>
!> \param[out] info
!> \verbatim
!>          info is INTEGER
!>          = 0:  successful exit
!>          = 1:  the update would violate positive-definiteness
!>          = 2:  R is singular
!>          = 3:  the diagonal element of u is not real
!> \endverbatim
!>
!> \ingroup choldecomp
subroutine cchinx(n,R,ldr,j,u,rw,info)
    integer n,j,ldr,info
    complex R(ldr,*),u(*),rw(*)
    external xerbla,ccopy,scnrm2,ctrsv,cqrtv1,cqrqh
    complex t
    real scnrm2,rho
    integer i

    ! check arguments
    info = 0
    if (n < 0) then
        info = -1
    else if (j < 1 .or. j > n+1) then
        info = -4
    end if
    if (info /= 0) then
        call xerbla('CCHINX',info)
        return
    end if

    ! shift vector.
    t = u(j)
    do i = j,n
        u(i) = u(i+1)
    end do
    ! the diagonal element must be real.
    if (imag(t) /= 0e0) goto 30

    ! check for singularity of R.
    do i = 1,n
        if (R(i,i) == 0e0) goto 20
    end do
    ! form R' \ u
    call ctrsv('U','C','N',n,R,ldr,u,1)
    rho = scnrm2(n,u,1)
    ! check positive definiteness.
    rho = real(t) - rho**2
    if (rho <= 0e0) goto 10
    ! shift columns
    do i = n,j,-1
        call ccopy(i,R(1,i),1,R(1,i+1),1)
        R(i+1,i+1) = 0e0
    end do
    call ccopy(n,u,1,R(1,j),1)
    R(n+1,j) = sqrt(rho)
    ! retriangularize
    if (j < n+1) then
        ! eliminate the introduced spike.
        call cqrtv1(n+2-j,R(j,j),rw)
        ! apply rotations to R
        call cqrqh(n+2-j,n+1-j,R(j,j+1),ldr,rw,R(j+1,j))
        ! zero spike.
        do i = j+1,n+1
            R(i,j) = 0e0
        end do
    end if
    ! normal return.
    return
    ! error returns.
    10 info = 1
    return
    20 info = 2
    return
    30 info = 3
    return
end subroutine