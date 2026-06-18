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
program tch1up
    integer n
    external stest, dtest, ctest, ztest, pstats

    write (*,*)
    write (*,*) 'testing Cholesky rank-1 update routines.'
    write (*,*) 'All residual errors are expected to be small.'
    write (*,*)

    n = 50
    write (*,*) 'sch1up test:'
    call stest(n)
    write (*,*) 'dch1up test:'
    call dtest(n)
    write (*,*) 'cch1up test:'
    call ctest(n)
    write (*,*) 'zch1up test:'
    call ztest(n)

    call pstats
end program

subroutine stest(n)
    integer n
    real A(n,n),R(n,n),u(n),wrk(n)
    external srandg,schgen,ssyr,sch1up,schchk
    ! set up random matrix & vectors
    call srandg(n,n,A,n)
    call srandg(n,1,u,n)
    ! generate A'*A and its Cholesky decomposition
    call schgen(n,A,n,R,n)
    ! update the matrix A
    call ssyr('U',n,1e0,u,1,A,n)
    ! update the Cholesky decomposition
    call sch1up(n,R,n,u,wrk)
    ! check result
    call schchk(n,A,n,R,n)

end subroutine

subroutine dtest(n)
    integer n
    double precision A(n,n),R(n,n),u(n),wrk(n)
    external drandg,dchgen,dsyr,dch1up,dchchk
    ! set up random matrix & vectors
    call drandg(n,n,A,n)
    call drandg(n,1,u,n)
    ! generate A'*A and its Cholesky decomposition
    call dchgen(n,A,n,R,n)
    ! update the matrix A
    call dsyr('U',n,1d0,u,1,A,n)
    ! update the Cholesky decomposition
    call dch1up(n,R,n,u,wrk)
    ! check result
    call dchchk(n,A,n,R,n)

end subroutine

subroutine ctest(n)
    integer n
    complex A(n,n),R(n,n),u(n)
    real rwrk(n)
    external crandg,cchgen,cher,cch1up,cchchk
    ! set up random matrix & vectors
    call crandg(n,n,A,n)
    call crandg(n,1,u,n)
    ! generate A'*A and its Cholesky decomposition
    call cchgen(n,A,n,R,n)
    ! update the matrix A
    call cher('U',n,1e0,u,1,A,n)
    ! update the Cholesky decomposition
    call cch1up(n,R,n,u,rwrk)
    ! check result
    call cchchk(n,A,n,R,n)

end subroutine

subroutine ztest(n)
    integer n
    double complex A(n,n),R(n,n),u(n)
    double precision rwrk(n)
    external zrandg,zchgen,zher,zch1up,zchchk
    ! set up random matrix & vectors
    call zrandg(n,n,A,n)
    call zrandg(n,1,u,n)
    ! generate A'*A and its Cholesky decomposition
    call zchgen(n,A,n,R,n)
    ! update the matrix A
    call zher('U',n,1d0,u,1,A,n)
    ! update the Cholesky decomposition
    call zch1up(n,R,n,u,rwrk)
    ! check result
    call zchchk(n,A,n,R,n)

    END
