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
!> \brief Applies a sequence of Givens rotations from the right to a matrix.
!>
!> \par Definition:
! =============
!> \verbatim
!>       subroutine sqrot(dir,m,n,Q,ldq,c,s)
!>
!>       .. Scalar Arguments ..
!>       character           dir
!>       integer             m, n, ldq
!>       ..
!>       .. Array Arguments ..
!>       real                Q(ldq,*)
!>       real                c(*)
!>       real                s(*)
!>       ..
!> \endverbatim
!>
!> \par Purpose:
! =============
!> \verbatim
!>
!> SQROT applies a sequence of Givens rotations from the right
!> side to an m-by-n matrix Q.  Given a direction indicator
!> dir, the rotation cosine and sine vectors c and s, SQROT
!> applies the rotations to Q, updating it in place.  If dir
!> is 'F' (forward), rotations are applied from the first to
!> the last; if dir is 'B' (backward), from the last to the
!> first.
!> \endverbatim
!>
!> \param[in] dir
!> \verbatim
!>          dir is CHARACTER
!>          If 'B' or 'b', rotations are applied backwards
!>          (from the last to the first).  If 'F' or 'f',
!>          rotations are applied forwards (from the first to
!>          the last).
!> \endverbatim
!>
!> \param[in] m
!> \verbatim
!>          m is INTEGER
!>          The number of rows of matrix Q.  m >= 0.
!> \endverbatim
!>
!> \param[in] n
!> \verbatim
!>          n is INTEGER
!>          The number of columns of the matrix Q.  n >= 0.
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is REAL array, dimension (ldq,*)
!>          On entry, the matrix Q.  On exit, the updated
!>          matrix Q1.
!> \endverbatim
!>
!> \param[in] ldq
!> \verbatim
!>          ldq is INTEGER
!>          The leading dimension of Q.  ldq >= m.
!> \endverbatim
!>
!> \param[in] c
!> \verbatim
!>          c is REAL array, dimension (*)
!>          The rotation cosines.  Must contain at least
!>          n-1 elements.
!> \endverbatim
!>
!> \param[in] s
!> \verbatim
!>          s is REAL array, dimension (*)
!>          The rotation sines.  Must contain at least
!>          n-1 elements.
!> \endverbatim
!>
!> \ingroup givens
subroutine sqrot(dir,m,n,Q,ldq,c,s)
  use iso_fortran_env
    character, intent(in) :: dir
    integer, intent(in) :: m, n, ldq
    real(real32), intent(inout) :: Q(ldq,*)
    real(real32), intent(in) :: c(*)
    real(real32), intent(in) :: s(*)
    external xerbla, srot,lsame
    logical lsame,fwd
    integer info,i
    ! quick return if possible.
    if (m == 0 .or. n == 0 .or. n == 1) return
    ! check arguments.
    info = 0
    fwd = lsame(dir,'F')
    if (.not.(fwd .or. lsame(dir,'B'))) then
        info = 1
    else if (m < 0) then
        info = 2
    else if (n < 0) then
        info = 3
    else if (ldq < m) then
        info = 5
    end if
    if (info /= 0) then
        call xerbla('SQROT',info)
        return
    end if

    if (fwd) then
        do i = 1,n-1
            call srot(m,Q(1,i),1,Q(1,i+1),1,c(i),s(i))
        end do
    else
        do i = n-1,1,-1
            call srot(m,Q(1,i),1,Q(1,i+1),1,c(i),s(i))
        end do
    end if
end subroutine
