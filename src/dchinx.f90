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
!>       subroutine dchinx(n,R,ldr,j,u,w,info)
!>
!>       .. Scalar Arguments ..
!>       integer            n, j, ldr, info
!>       ..
!>       .. Array Arguments ..
!>       double precision   R(ldr,*), u(*), w(*)
!>       ..
!> \endverbatim
!>
!> \par Purpose:
! =============
!> \verbatim
!>
!> DCHINX updates the Cholesky factorization of a symmetric
!> positive definite matrix A after inserting a row and column.
!> Given an upper triangular matrix R that is a Cholesky factor of
!> A, i.e., A = R.'*R, where R.' denotes the transpose of R, this
!> DCHINX updates R -> R1 so that R1.'*R1 = A1, where
!> A1(jj,jj) = A, A1(j,:) = u.', A1(:,j) = u, and
!> jj = [1:j-1, j+1:n+1].
!>
!> On exit, u is destroyed and R is extended by one row and column.
!> The insertion is performed by first solving R.'*u = v, checking
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
!>          R is DOUBLE PRECISION array, dimension (ldr,n+1)
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
!>          u is DOUBLE PRECISION array, dimension (n+1)
!>          On entry, the vector defining the inserted row/column.
!>          On exit, u is destroyed.
!> \endverbatim
!>
!> \param[out] w
!> \verbatim
!>          w is DOUBLE PRECISION array, dimension (n+1)
!>          Workspace vector used during the retriangularization.
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
subroutine dchinx(n,R,ldr,j,u,w,info)
    integer n,j,ldr,info
    double precision R(ldr,*),u(*),w(*)
    external xerbla,dcopy,dnrm2,dtrsv,dqrtv1,dqrqh
    double precision dnrm2,rho,t,rr
    integer i

    ! check arguments
    info = 0
    if (n < 0) then
        info = -1
    else if (j < 1 .or. j > n+1) then
        info = -4
    end if
    if (info /= 0) then
        call xerbla('DCHINX',-info)
        return
    end if

    ! shift vector.
    t = u(j)
    do i = j,n
        u(i) = u(i+1)
    end do

    ! check for singularity of R.
    do i = 1,n
        if (R(i,i) == 0d0) goto 20
    end do
    ! form R' \ u
    call dtrsv('U','T','N',n,R,ldr,u,1)
    rho = dnrm2(n,u,1)
    ! check positive definiteness.
    rho = t - rho**2
    if (rho <= 0d0) goto 10
    ! shift columns
    do i = n,j,-1
        call dcopy(i,R(1,i),1,R(1,i+1),1)
        R(i+1,i+1) = 0d0
    end do
    call dcopy(n,u,1,R(1,j),1)
    R(n+1,j) = sqrt(rho)
    ! retriangularize
    if (j < n+1) then
        ! eliminate the introduced spike.
        call dqrtv1(n+2-j,R(j,j),w)
        ! apply rotations to R
        call dqrqh(n+2-j,n+1-j,R(j,j+1),ldr,w,R(j+1,j))
        ! zero spike.
        do i = j+1,n+1
            R(i,j) = 0d0
        end do
    end if
    ! normal return.
    return
    ! error returns.
    10 info = 1
    return
    20 info = 2
    return
end subroutine
