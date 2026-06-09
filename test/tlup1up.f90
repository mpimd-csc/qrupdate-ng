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
program tlup1up
    integer m,n
    external stest, dtest, ctest, ztest, pstats



    write (*,*)
    write (*,*) 'testing pivoted LU rank-1 update routines.'
    write (*,*) 'All residual errors are expected to be small.'
    write (*,*)

    m = 60
    n = 40
    write (*,*) 'slup1up test (rows > columns):'
    call stest(m,n)
    write (*,*) 'dlup1up test (rows > columns):'
    call dtest(m,n)
    write (*,*) 'clup1up test (rows > columns):'
    call ctest(m,n)
    write (*,*) 'zlup1up test (rows > columns):'
    call ztest(m,n)

    m = 40
    n = 60
    write (*,*) 'slup1up test (rows < columns):'
    call stest(m,n)
    write (*,*) 'dlup1up test (rows < columns):'
    call dtest(m,n)
    write (*,*) 'clup1up test (rows < columns):'
    call ctest(m,n)
    write (*,*) 'zlup1up test (rows < columns):'
    call ztest(m,n)

    call pstats
end program

subroutine stest(m,n)
    integer m,n
    real A(m,n),L(m,min(m,n)),R(min(m,n),n),u(m),v(n),w(m)
    external srandg,slupgen,sger,slup1up,slupchk
    integer k,p(m)
    ! set up random matrix & vectors
    call srandg(m,n,A,m)
    call srandg(m,1,u,m)
    call srandg(n,1,v,n)
    k = min(m,n)
    ! generate LU decomposition
    call slupgen(m,n,A,m,L,m,R,k,p)
    ! update A
    call sger(m,n,1e0,u,1,v,1,A,m)
    ! update the pivoted LU decomposition
    call slup1up(m,n,L,m,R,k,p,u,v,w)
    ! check result
    call slupchk(m,n,A,m,L,m,R,k,p)
end subroutine

subroutine dtest(m,n)
    integer m,n
    double precision A(m,n),L(m,min(m,n)),R(min(m,n),n),u(m),v(n),w(m)
    external drandg,dlupgen,dger,dlup1up,dlupchk
    integer k,p(m)
    ! set up random matrix & vectors
    call drandg(m,n,A,m)
    call drandg(m,1,u,m)
    call drandg(n,1,v,n)
    k = min(m,n)
    ! generate LU decomposition
    call dlupgen(m,n,A,m,L,m,R,k,p)
    ! update A
    call dger(m,n,1d0,u,1,v,1,A,m)
    ! update the pivoted LU decomposition
    call dlup1up(m,n,L,m,R,k,p,u,v,w)
    ! check result
    call dlupchk(m,n,A,m,L,m,R,k,p)
end subroutine

subroutine ctest(m,n)
    integer m,n
    complex A(m,n),L(m,min(m,n)),R(min(m,n),n),u(m),v(n),w(m)
    external crandg,clupgen,cgeru,clup1up,clupchk
    integer k,p(m)
    ! set up random matrix & vectors
    call crandg(m,n,A,m)
    call crandg(m,1,u,m)
    call crandg(n,1,v,n)
    k = min(m,n)
    ! generate LU decomposition
    call clupgen(m,n,A,m,L,m,R,k,p)
    ! update A
    call cgeru(m,n,(1e0,0e0),u,1,v,1,A,m)
    ! update the pivoted LU decomposition
    call clup1up(m,n,L,m,R,k,p,u,v,w)
    ! check result
    call clupchk(m,n,A,m,L,m,R,k,p)
end subroutine

subroutine ztest(m,n)
    integer m,n
    double complex A(m,n),L(m,min(m,n)),R(min(m,n),n),u(m),v(n),w(m)
    external zrandg,zlupgen,zgeru,zlup1up,zlupchk
    integer k,p(m)
    ! set up random matrix & vectors
    call zrandg(m,n,A,m)
    call zrandg(m,1,u,m)
    call zrandg(n,1,v,n)
    k = min(m,n)
    ! generate LU decomposition
    call zlupgen(m,n,A,m,L,m,R,k,p)
    ! update A
    call zgeru(m,n,(1d0,0d0),u,1,v,1,A,m)
    ! update the pivoted LU decomposition
    call zlup1up(m,n,L,m,R,k,p,u,v,w)
    ! check result
    call zlupchk(m,n,A,m,L,m,R,k,p)
    END
