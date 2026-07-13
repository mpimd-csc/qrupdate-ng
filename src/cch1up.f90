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
!> \brief Updates a Cholesky factorization after a rank-1 modification.
!>
!> \par Definition:
! =============
!> \verbatim
!>       subroutine cch1up(n,R,ldr,u,w)
!>
!>       .. Scalar Arguments ..
!>       integer            n, ldr
!>       ..
!>       .. Array Arguments ..
!>       complex            R(ldr,*), u(*)
!>       real               w(*)
!>       ..
!> \endverbatim
!>
!> \par Purpose:
! =============
!> \verbatim
!>
!> CCH1UP updates the Cholesky factorization of a hermitian
!> positive definite matrix A after a rank-1 modification.  Given an
!> upper triangular matrix R that is a Cholesky factor of A, i.e.,
!> A = R'*R, where R' denotes the conjugate transpose of R, this
!> CCH1UP updates R -> R1 so that R1'*R1 = A + u*u', where u is
!> a given vector.
!>
!> The update is performed by applying a sequence of Givens rotations
!> to restore the upper triangular structure of R.  On exit, u
!> contains the rotation sines and w contains the rotation cosines
!> used in the transformation.
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
!>          R is COMPLEX array, dimension (ldr,n)
!>          On entry, the upper triangular matrix R, the Cholesky
!>          factor of A.  On exit, the updated upper triangular
!>          matrix R1, the Cholesky factor of A + u*u'.
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
!>          u is COMPLEX array, dimension (n)
!>          On entry, the vector determining the rank-1 update.
!>          On exit, u contains the rotation sines used to
!>          transform R to R1.
!> \endverbatim
!>
!> \param[out] w
!> \verbatim
!>          w is REAL array, dimension (n)
!>          On exit, w contains the cosine parts of the Givens
!>          rotations used to transform R to R1.
!> \endverbatim
!>
!> \ingroup choldecomp
subroutine cch1up(n,R,ldr,u,w)
    use iso_fortran_env
    integer, intent(in) :: n, ldr
    complex(real32), intent(inout) :: R(ldr,*)
    complex(real32), intent(inout) :: u(*)
    real(real32), intent(out) :: w(*)
    external clartg
    complex(real32) rr,ui,t
    integer i,j

    do i = 1,n
        ! apply stored rotations, column-wise
        ui = conjg(u(i))
        do j = 1,i-1
            t = w(j)*R(j,i) + u(j)*ui
            ui = w(j)*ui - conjg(u(j))*R(j,i)
            R(j,i) = t
        end do
        ! generate next rotation
        call clartg(R(i,i),ui,w(i),u(i),rr)
        R(i,i) = rr
    end do
end subroutine
