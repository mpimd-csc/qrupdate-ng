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
program tlu1up
    integer m,n
    external stest, dtest, ctest, ztest, pstats


    write (*,*)
    write (*,*) 'testing LU rank-1 update routines.'
    write (*,*) 'All residual errors are expected to be small.'
    write (*,*)

    m = 60
    n = 40
    write (*,*) 'slu1up test (rows > columns):'
    call stest(m,n)
    write (*,*) 'dlu1up test (rows > columns):'
    call dtest(m,n)
    write (*,*) 'clu1up test (rows > columns):'
    call ctest(m,n)
    write (*,*) 'zlu1up test (rows > columns):'
    call ztest(m,n)

    m = 40
    n = 60
    write (*,*) 'slu1up test (rows < columns):'
    call stest(m,n)
    write (*,*) 'dlu1up test (rows < columns):'
    call dtest(m,n)
    write (*,*) 'clu1up test (rows < columns):'
    call ctest(m,n)
    write (*,*) 'zlu1up test (rows < columns):'
    call ztest(m,n)

    call pstats
end program

subroutine stest(m,n)
    integer m,n
    real A(m,n),L(m,min(m,n)),R(min(m,n),n),u(m),v(n)
    external srandg,slugen,sger,slu1up,sluchk
    integer k
    ! set up random matrix & vectors
    call srandg(m,n,A,m)
    call srandg(m,1,u,m)
    call srandg(n,1,v,n)
    k = min(m,n)
    ! generate LU decomposition
    call slugen(m,n,A,m,L,m,R,k)
    ! update A
    call sger(m,n,1e0,u,1,v,1,A,m)
    ! update the LU decomposition
    call slu1up(m,n,L,m,R,k,u,v)
    ! check result
    call sluchk(m,n,A,m,L,m,R,k)

end subroutine

subroutine dtest(m,n)
    integer m,n
    double precision A(m,n),L(m,min(m,n)),R(min(m,n),n),u(m),v(n)
    external drandg,dlugen,dger,dlu1up,dluchk
    integer k
    ! set up random matrix & vectors
    call drandg(m,n,A,m)
    call drandg(m,1,u,m)
    call drandg(n,1,v,n)
    k = min(m,n)
    ! generate LU decomposition
    call dlugen(m,n,A,m,L,m,R,k)
    ! update A
    call dger(m,n,1d0,u,1,v,1,A,m)
    ! update the LU decomposition
    call dlu1up(m,n,L,m,R,k,u,v)
    ! check result
    call dluchk(m,n,A,m,L,m,R,k)

end subroutine

subroutine ctest(m,n)
    integer m,n
    complex A(m,n),L(m,min(m,n)),R(min(m,n),n),u(m),v(n)
    external crandg,clugen,cgeru,clu1up,cluchk
    integer k
    ! set up random matrix & vectors
    call crandg(m,n,A,m)
    call crandg(m,1,u,m)
    call crandg(n,1,v,n)
    k = min(m,n)
    ! generate LU decomposition
    call clugen(m,n,A,m,L,m,R,k)
    ! update A
    call cgeru(m,n,(1e0,0e0),u,1,v,1,A,m)
    ! update the LU decomposition
    call clu1up(m,n,L,m,R,k,u,v)
    ! check result
    call cluchk(m,n,A,m,L,m,R,k)

end subroutine

subroutine ztest(m,n)
    integer m,n
    double complex A(m,n),L(m,min(m,n)),R(min(m,n),n),u(m),v(n)
    external zrandg,zlugen,zgeru,zlu1up,zluchk
    integer k
    ! set up random matrix & vectors
    call zrandg(m,n,A,m)
    call zrandg(m,1,u,m)
    call zrandg(n,1,v,n)
    k = min(m,n)
    ! generate LU decomposition
    call zlugen(m,n,A,m,L,m,R,k)
    ! update A
    call zgeru(m,n,(1d0,0d0),u,1,v,1,A,m)
    ! update the LU decomposition
    call zlu1up(m,n,L,m,R,k,u,v)
    ! check result
    call zluchk(m,n,A,m,L,m,R,k)

    END
