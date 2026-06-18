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
program tqrinr
    integer m,n,j
    external stest, dtest, ctest, ztest, pstats


    write (*,*)
    write (*,*) 'testing QR row insert routines.'
    write (*,*) 'All residual errors are expected to be small.'
    write (*,*)

    m = 60
    n = 40
    j = 30
    write (*,*) 'sqrinr test (full factorization):'
    call stest(m,n,j)
    write (*,*) 'dqrinr test (full factorization):'
    call dtest(m,n,j)
    write (*,*) 'cqrinr test (full factorization):'
    call ctest(m,n,j)
    write (*,*) 'zqrinr test (full factorization):'
    call ztest(m,n,j)

    call pstats
end program

subroutine stest(m,n,j)
    integer m,n,j
    real A(m+1,max(m+1,n)),Q(m+1,m+1),R(m+1,n),u(n),wrk(n)
    external srandg,sqrgen,scopy,sqrinr,sqrchk
    integer i
    ! set up random matrix & vector
    call srandg(m,n,A,m+1)
    call srandg(n,1,u,n)
    ! generate QR decomposition
    call sqrgen(m,n,A,m+1,Q,m+1,R,m+1)
    ! update A
    do i = m,j,-1
        call scopy(n,A(i,1),m+1,A(i+1,1),m+1)
    end do
    call scopy(n,u,1,A(j,1),m+1)
    ! update the QR decomposition
    call sqrinr(m,n,Q,m+1,R,m+1,j,u,wrk)
    ! check result
    call sqrchk(m+1,n,m+1,A,m+1,Q,m+1,R,m+1)

end subroutine

subroutine dtest(m,n,j)
    integer m,n,j
    double precision A(m+1,max(m+1,n)),Q(m+1,m+1),R(m+1,n),u(n),wrk(n)
    external drandg,dqrgen,dcopy,dqrinr,dqrchk
    integer i
    ! set up random matrix & vector
    call drandg(m,n,A,m+1)
    call drandg(n,1,u,n)
    ! generate QR decomposition
    call dqrgen(m,n,A,m+1,Q,m+1,R,m+1)
    ! update A
    do i = m,j,-1
        call dcopy(n,A(i,1),m+1,A(i+1,1),m+1)
    end do
    call dcopy(n,u,1,A(j,1),m+1)
    ! update the QR decomposition
    call dqrinr(m,n,Q,m+1,R,m+1,j,u,wrk)
    ! check result
    call dqrchk(m+1,n,m+1,A,m+1,Q,m+1,R,m+1)

end subroutine

subroutine ctest(m,n,j)
    integer m,n,j
    complex A(m+1,max(m+1,n)),Q(m+1,m+1),R(m+1,n),u(n),wrk(n)
    external crandg,cqrgen,ccopy,cqrinr,cqrchk
    integer i
    ! set up random matrix & vector
    call crandg(m,n,A,m+1)
    call crandg(n,1,u,n)
    ! generate QR decomposition
    call cqrgen(m,n,A,m+1,Q,m+1,R,m+1)
    ! update A
    do i = m,j,-1
        call ccopy(n,A(i,1),m+1,A(i+1,1),m+1)
    end do
    call ccopy(n,u,1,A(j,1),m+1)
    ! update the QR decomposition
    call cqrinr(m,n,Q,m+1,R,m+1,j,u,wrk)
    ! check result
    call cqrchk(m+1,n,m+1,A,m+1,Q,m+1,R,m+1)

end subroutine

subroutine ztest(m,n,j)
    integer m,n,j
    double complex A(m+1,max(m+1,n)),Q(m+1,m+1),R(m+1,n),u(n),wrk(n)
    external zrandg,zqrgen,zcopy,zqrinr,zqrchk
    integer i
    ! set up random matrix & vector
    call zrandg(m,n,A,m+1)
    call zrandg(n,1,u,n)
    ! generate QR decomposition
    call zqrgen(m,n,A,m+1,Q,m+1,R,m+1)
    ! update A
    do i = m,j,-1
        call zcopy(n,A(i,1),m+1,A(i+1,1),m+1)
    end do
    call zcopy(n,u,1,A(j,1),m+1)
    ! update the QR decomposition
    call zqrinr(m,n,Q,m+1,R,m+1,j,u,wrk)
    ! check result
    call zqrchk(m+1,n,m+1,A,m+1,Q,m+1,R,m+1)

    END
