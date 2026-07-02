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
!> \brief Generates Givens rotations to eliminate all but the first element of a vector.
!>
!> \par Definition:
! =============
!> \verbatim
!>       subroutine sqrtv1(n,u,w)
!>
!>       .. Scalar Arguments ..
!>       integer            n
!>       ..
!>       .. Array Arguments ..
!>       real               u(*), w(*)
!>       ..
!> \endverbatim
!>
!> \par Purpose:
! =============
!> \verbatim
!>
!> SQRTV1 generates a sequence of n-1 Givens rotations that
!> eliminate all but the first element of a real vector u.
!> On entry, u contains the vector to be reduced.  On exit,
!> u(1) contains the resulting element, u(2:n) contains the
!> rotation sines, and w contains the rotation cosines.
!>
!> The rotations are generated from the bottom up, so that the
!> first rotation eliminates u(n), the second eliminates u(n-1),
!> and so on.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is INTEGER
!>          The length of the vector u.  If n <= 0, the subroutine
!>          returns immediately without modification.
!> \endverbatim
!>
!> \param[in,out] u
!> \verbatim
!>          u is REAL array, dimension (n)
!>          On entry, the vector to be reduced.  On exit, u(1)
!>          contains the remaining element, and u(2:n) contains
!>          the sine parts of the Givens rotations.
!> \endverbatim
!>
!> \param[out] w
!> \verbatim
!>          w is REAL array, dimension (n)
!>          On exit, w contains the cosine parts of the Givens
!>          rotations.
!> \endverbatim
!>
!> \ingroup givens
subroutine sqrtv1(n,u,w)
  use iso_fortran_env
    integer, intent(in) :: n
    real(real32), intent(inout) :: u(*)
    real(real32), intent(out) :: w(*)
    external slartg
    real(real32) rr,t
    integer i
    ! quick return if possible.
    if (n <= 0) return
    rr = u(n)
    do i = n-1,1,-1
        call slartg(u(i),rr,w(i),u(i+1),t)
        rr = t
    end do
    u(1) = rr
end subroutine
