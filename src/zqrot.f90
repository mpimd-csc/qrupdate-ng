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
subroutine zqrot(dir,m,n,Q,ldq,c,s)
    ! purpose:      Apply a sequence of inv. rotations from right
    !
    ! arguments:
    ! dir (in)      if 'B' or 'b', rotations are applied from backwards
    !               if 'F' or 'f', from forwards.
    ! m (in)        number of rows of matrix Q
    ! n (in)        number of columns of the matrix Q
    ! Q (io)        on entry, the matrix Q
    !               on exit, the updated matrix Q1
    ! ldq (in)      the leading dimension of Q
    ! c (in)        n-1 rotation cosines
    ! s (in)        n-1 rotation sines
    !
    character dir
    integer m,n,ldq
    double complex Q(ldq,*),s(*)
    double precision c(*)
    external zrot,lsame,xerbla
    logical lsame,fwd
    integer info,i
    ! quick return if possible
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
        call xerbla('ZQROT',info)
        return
    end if

    if (fwd) then
        do i = 1,n-1
            call zrot(m,Q(1,i),1,Q(1,i+1),1,c(i),conjg(s(i)))
        end do
    else
        do i = n-1,1,-1
            call zrot(m,Q(1,i),1,Q(1,i+1),1,c(i),conjg(s(i)))
        end do
    end if
end subroutine
