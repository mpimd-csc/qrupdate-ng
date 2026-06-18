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
program tchdex
    integer n,j
    external stest, dtest, ctest, ztest, pstats


    write (*,*)
    write (*,*) 'testing Cholesky symmetric delete routines.'
    write (*,*) 'All residual errors are expected to be small.'
    write (*,*)

    n = 50
    j = 15
    write (*,*) 'schdex test:'
    call stest(n,j)
    write (*,*) 'dchdex test:'
    call dtest(n,j)
    write (*,*) 'cchdex test:'
    call ctest(n,j)
    write (*,*) 'zchdex test:'
    call ztest(n,j)

    call pstats
end program

subroutine stest(n,j)
    integer n,j
    real A(n,n),R(n,n),wrk(n)
    external srandg,schgen,schdex,schchk
    integer i,k
    ! set up random matrix & vectors
    call srandg(n,n,A,n)
    ! generate A'*A and its Cholesky decomposition
    call schgen(n,A,n,R,n)
    ! update the matrix A
    do k = j,n-1
        do i = 1,j-1
            A(i,k) = A(i,k+1)
        end do
        do i = j,k
            A(i,k) = A(i+1,k+1)
        end do
    end do
    ! update the Cholesky decomposition
    call schdex(n,R,n,j,wrk)
    ! check result
    call schchk(n-1,A,n,R,n)

end subroutine

subroutine dtest(n,j)
    integer n,j
    double precision A(n,n),R(n,n),wrk(n)
    external drandg,dchgen,dchdex,dchchk
    integer i,k
    ! set up random matrix & vectors
    call drandg(n,n,A,n)
    ! generate A'*A and its Cholesky decomposition
    call dchgen(n,A,n,R,n)
    ! update the matrix A
    do k = j,n-1
        do i = 1,j-1
            A(i,k) = A(i,k+1)
        end do
        do i = j,k
            A(i,k) = A(i+1,k+1)
        end do
    end do
    ! update the Cholesky decomposition
    call dchdex(n,R,n,j,wrk)
    ! check result
    call dchchk(n-1,A,n,R,n)

end subroutine

subroutine ctest(n,j)
    integer n,j
    complex A(n,n),R(n,n)
    real rwrk(n)
    external crandg,cchgen,cchdex,cchchk
    integer i,k
    ! set up random matrix & vectors
    call crandg(n,n,A,n)
    ! generate A'*A and its Cholesky decomposition
    call cchgen(n,A,n,R,n)
    ! update the matrix A
    do k = j,n-1
        do i = 1,j-1
            A(i,k) = A(i,k+1)
        end do
        do i = j,k
            A(i,k) = A(i+1,k+1)
        end do
    end do
    ! update the Cholesky decomposition
    call cchdex(n,R,n,j,rwrk)
    ! check result
    call cchchk(n-1,A,n,R,n)

end subroutine

subroutine ztest(n,j)
    integer n,j
    double complex A(n,n),R(n,n)
    double precision rwrk(n)
    external zrandg,zchgen,zchdex,zchchk
    integer i,k
    ! set up random matrix & vectors
    call zrandg(n,n,A,n)
    ! generate A'*A and its Cholesky decomposition
    call zchgen(n,A,n,R,n)
    ! update the matrix A
    do k = j,n-1
        do i = 1,j-1
            A(i,k) = A(i,k+1)
        end do
        do i = j,k
            A(i,k) = A(i+1,k+1)
        end do
    end do
    ! update the Cholesky decomposition
    call zchdex(n,R,n,j,rwrk)
    ! check result
    call zchchk(n-1,A,n,R,n)

    END
