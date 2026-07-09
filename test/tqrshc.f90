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
program tqrshc
    integer m,n,i,j
    external stest, dtest, ctest, ztest, pstats


    write (*,*)
    write (*,*) 'testing QR column shift routines.'
    write (*,*) 'All residual errors are expected to be small.'
    write (*,*)

    m = 60
    n = 50
    i = 20
    j = 40
    write (*,*) 'sqrshc test (left shift, full factorization):'
    call stest(m,n,i,j,0)
    write (*,*) 'dqrshc test (left shift, full factorization):'
    call dtest(m,n,i,j,0)
    write (*,*) 'cqrshc test (left shift, full factorization):'
    call ctest(m,n,i,j,0)
    write (*,*) 'zqrshc test (left shift, full factorization):'
    call ztest(m,n,i,j,0)

    i = 40
    j = 20
    write (*,*) 'sqrshc test (right shift, economized factorization):'
    call stest(m,n,i,j,1)
    write (*,*) 'dqrshc test (right shift, economized factorization):'
    call dtest(m,n,i,j,1)
    write (*,*) 'cqrshc test (right shift, economized factorization):'
    call ctest(m,n,i,j,1)
    write (*,*) 'zqrshc test (right shift, economized factorization):'
    call ztest(m,n,i,j,1)

    call pstats
end program

subroutine stest(m,n,i,j,ec)
    use iso_fortran_env
    integer m,n,i,j,ec
    real(real32) A(m,max(m,n)),Q(m,m),R(m,n),wrk(2*m)
    external srandg,sqrgen,scopy,sqrshc,sqrchk
    integer k
    ! set up random matrix & vector
    call srandg(m,n,A,m)
    ! generate QR decomposition
    call sqrgen(m,n,A,m,Q,m,R,m)
    ! update A
    if (i < j) then
        call scopy(m,A(1,i),1,wrk,1)
        do k = i,j-1
            call scopy(m,A(1,k+1),1,A(1,k),1)
        end do
        call scopy(m,wrk,1,A(1,j),1)
    else
        call scopy(m,A(1,i),1,wrk,1)
        do k = i,j+1,-1
            call scopy(m,A(1,k-1),1,A(1,k),1)
        end do
        call scopy(m,wrk,1,A(1,j),1)
    end if
    ! update the QR decomposition
    k = m
    if (ec == 1) k = n
    call sqrshc(m,n,k,Q,m,R,m,i,j,wrk)
    ! check result
    call sqrchk(m,n,k,A,m,Q,m,R,m)

end subroutine

subroutine dtest(m,n,i,j,ec)
    use iso_fortran_env
    integer m,n,i,j,ec
    real(real64) A(m,max(m,n)),Q(m,m),R(m,n),wrk(2*m)
    external drandg,dqrgen,dcopy,dqrshc,dqrchk
    integer k
    ! set up random matrix & vector
    call drandg(m,n,A,m)
    ! generate QR decomposition
    call dqrgen(m,n,A,m,Q,m,R,m)
    ! update A
    if (i < j) then
        call dcopy(m,A(1,i),1,wrk,1)
        do k = i,j-1
            call dcopy(m,A(1,k+1),1,A(1,k),1)
        end do
        call dcopy(m,wrk,1,A(1,j),1)
    else
        call dcopy(m,A(1,i),1,wrk,1)
        do k = i,j+1,-1
            call dcopy(m,A(1,k-1),1,A(1,k),1)
        end do
        call dcopy(m,wrk,1,A(1,j),1)
    end if
    ! update the QR decomposition
    k = m
    if (ec == 1) k = n
    call dqrshc(m,n,k,Q,m,R,m,i,j,wrk)
    ! check result
    call dqrchk(m,n,k,A,m,Q,m,R,m)

end subroutine

subroutine ctest(m,n,i,j,ec)
    use iso_fortran_env
    integer m,n,i,j,ec
    complex(real32) A(m,max(m,n)),Q(m,m),R(m,n),wrk(m)
    real(real32) rwrk(m)
    external crandg,cqrgen,ccopy,cqrshc,cqrchk
    integer k
    ! set up random matrix & vector
    call crandg(m,n,A,m)
    ! generate QR decomposition
    call cqrgen(m,n,A,m,Q,m,R,m)
    ! update A
    if (i < j) then
        call ccopy(m,A(1,i),1,wrk,1)
        do k = i,j-1
            call ccopy(m,A(1,k+1),1,A(1,k),1)
        end do
        call ccopy(m,wrk,1,A(1,j),1)
    else
        call ccopy(m,A(1,i),1,wrk,1)
        do k = i,j+1,-1
            call ccopy(m,A(1,k-1),1,A(1,k),1)
        end do
        call ccopy(m,wrk,1,A(1,j),1)
    end if
    ! update the QR decomposition
    k = m
    if (ec == 1) k = n
    call cqrshc(m,n,k,Q,m,R,m,i,j,wrk,rwrk)
    ! check result
    call cqrchk(m,n,k,A,m,Q,m,R,m)

end subroutine

subroutine ztest(m,n,i,j,ec)
    use iso_fortran_env
    integer m,n,i,j,ec
    complex(real64) A(m,max(m,n)),Q(m,m),R(m,n),wrk(m)
    real(real64) rwrk(m)
    external zrandg,zqrgen,zcopy,zqrshc,zqrchk
    integer k
    ! set up random matrix & vector
    call zrandg(m,n,A,m)
    ! generate QR decomposition
    call zqrgen(m,n,A,m,Q,m,R,m)
    ! update A
    if (i < j) then
        call zcopy(m,A(1,i),1,wrk,1)
        do k = i,j-1
            call zcopy(m,A(1,k+1),1,A(1,k),1)
        end do
        call zcopy(m,wrk,1,A(1,j),1)
    else
        call zcopy(m,A(1,i),1,wrk,1)
        do k = i,j+1,-1
            call zcopy(m,A(1,k-1),1,A(1,k),1)
        end do
        call zcopy(m,wrk,1,A(1,j),1)
    end if
    ! update the QR decomposition
    k = m
    if (ec == 1) k = n
    call zqrshc(m,n,k,Q,m,R,m,i,j,wrk,rwrk)
    ! check result
    call zqrchk(m,n,k,A,m,Q,m,R,m)

    END
