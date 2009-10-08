c Copyright (C) 2008, 2009  VZLU Prague, a.s., Czech Republic
c
c Author: Jaroslav Hajek <highegg@gmail.com>
c
c This file is part of qrupdate.
c
c qrupdate is free software; you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation; either version 3 of the License, or
c (at your option) any later version.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c
c You should have received a copy of the GNU General Public License
c along with this software; see the file COPYING.  If not, see
c <http://www.gnu.org/licenses/>.
c
      program tlup1up
      integer m,n

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
      external srandg,slugen,sger,slup1up,sluchk
      integer k,p(m)
c set up random matrix & vectors
      call srandg(m,n,A,m)
      call srandg(m,1,u,m)
      call srandg(n,1,v,n)
      k = min(m,n)
c generate LU decomposition
      call slupgen(m,n,A,m,L,m,R,k,p)
c update A
      call sger(m,n,1e0,u,1,v,1,A,m)
c update the pivoted LU decomposition
      call slup1up(m,n,L,m,R,k,p,u,v,w)
c check result
      call slupchk(m,n,A,m,L,m,R,k,p)
      end subroutine

      subroutine dtest(m,n)
      integer m,n
      double precision A(m,n),L(m,min(m,n)),R(min(m,n),n),u(m),v(n),w(m)
      external drandg,dlugen,dger,dlup1up,dluchk
      integer k,p(m)
c set up random matrix & vectors
      call drandg(m,n,A,m)
      call drandg(m,1,u,m)
      call drandg(n,1,v,n)
      k = min(m,n)
c generate LU decomposition
      call dlupgen(m,n,A,m,L,m,R,k,p)
c update A
      call dger(m,n,1d0,u,1,v,1,A,m)
c update the pivoted LU decomposition
      call dlup1up(m,n,L,m,R,k,p,u,v,w)
c check result
      call dlupchk(m,n,A,m,L,m,R,k,p)
      end subroutine

      subroutine ctest(m,n)
      integer m,n
      complex A(m,n),L(m,min(m,n)),R(min(m,n),n),u(m),v(n),w(m)
      external crandg,clugen,cgeru,clup1up,cluchk
      integer k,p(m)
c set up random matrix & vectors
      call crandg(m,n,A,m)
      call crandg(m,1,u,m)
      call crandg(n,1,v,n)
      k = min(m,n)
c generate LU decomposition
      call clupgen(m,n,A,m,L,m,R,k,p)
c update A
      call cgeru(m,n,(1e0,0e0),u,1,v,1,A,m)
c update the pivoted LU decomposition
      call clup1up(m,n,L,m,R,k,p,u,v,w)
c check result
      call clupchk(m,n,A,m,L,m,R,k,p)
      end subroutine

      subroutine ztest(m,n)
      integer m,n
      double complex A(m,n),L(m,min(m,n)),R(min(m,n),n),u(m),v(n),w(m)
      external zrandg,zlugen,zgeru,zlup1up,zluchk
      integer k,p(m)
c set up random matrix & vectors
      call zrandg(m,n,A,m)
      call zrandg(m,1,u,m)
      call zrandg(n,1,v,n)
      k = min(m,n)
c generate LU decomposition
      call zlupgen(m,n,A,m,L,m,R,k,p)
c update A
      call zgeru(m,n,(1d0,0d0),u,1,v,1,A,m)
c update the pivoted LU decomposition
      call zlup1up(m,n,L,m,R,k,p,u,v,w)
c check result
      call zlupchk(m,n,A,m,L,m,R,k,p)
      end subroutine
