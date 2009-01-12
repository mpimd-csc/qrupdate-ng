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
      program tqr1up
      integer m,n

      write (*,*)
      write (*,*) 'testing QR rank-1 update routines.'
      write (*,*) 'All residual errors are expected to be small.'
      write (*,*)

      m = 60
      n = 40
      write (*,*) 'sqr1up test (full factorization):'
      call stest(m,n,0)
      write (*,*) 'dqr1up test (full factorization):'
      call dtest(m,n,0)
      write (*,*) 'cqr1up test (full factorization):'
      call ctest(m,n,0)
      write (*,*) 'zqr1up test (full factorization):'
      call ztest(m,n,0)

      write (*,*) 'sqr1up test (economized factorization):'
      call stest(m,n,1)
      write (*,*) 'dqr1up test (economized factorization):'
      call dtest(m,n,1)
      write (*,*) 'cqr1up test (economized factorization):'
      call ctest(m,n,1)
      write (*,*) 'zqr1up test (economized factorization):'
      call ztest(m,n,1)

      m = 40
      n = 60
      write (*,*) 'sqr1up test (rows < columns):'
      call stest(m,n,0)
      write (*,*) 'dqr1up test (rows < columns):'
      call dtest(m,n,0)
      write (*,*) 'cqr1up test (rows < columns):'
      call ctest(m,n,0)
      write (*,*) 'zqr1up test (rows < columns):'
      call ztest(m,n,0)

      call pstats
      end program

      subroutine stest(m,n,ec)
      integer m,n,ec
      real A(m,max(m,n)),Q(m,m),R(m,n),u(m),v(n),wrk(2*m)
      external srandg,sqrgen,sger,sqr1up,sqrchk
      integer k
c set up random matrix & vectors
      call srandg(m,n,A,m)
      call srandg(m,1,u,m)
      call srandg(n,1,v,n)
c generate QR decomposition
      call sqrgen(m,n,A,m,Q,m,R,m)
c update A
      call sger(m,n,1e0,u,1,v,1,A,m)
c update the QR decomposition
      k = m
      if (ec == 1) k = n
      call sqr1up(m,n,k,Q,m,R,m,u,v,wrk)
c check result
      call sqrchk(m,n,k,A,m,Q,m,R,m)

      end subroutine

      subroutine dtest(m,n,ec)
      integer m,n,ec
      double precision A(m,max(m,n)),Q(m,m),R(m,n),u(m),v(n),wrk(2*m)
      external drandg,dqrgen,dger,dqr1up,dqrchk
      integer k
c set up random matrix & vectors
      call drandg(m,n,A,m)
      call drandg(m,1,u,m)
      call drandg(n,1,v,n)
c generate QR decomposition
      call dqrgen(m,n,A,m,Q,m,R,m)
c update A
      call dger(m,n,1d0,u,1,v,1,A,m)
c update the QR decomposition
      k = m
      if (ec == 1) k = n
      call dqr1up(m,n,k,Q,m,R,m,u,v,wrk)
c check result
      call dqrchk(m,n,k,A,m,Q,m,R,m)

      end subroutine

      subroutine ctest(m,n,ec)
      integer m,n,ec
      complex A(m,max(m,n)),Q(m,m),R(m,n),u(m),v(n),wrk(m)
      real rwrk(m)
      external crandg,cqrgen,cgerc,cqr1up,cqrchk
      integer k
c set up random matrix & vectors
      call crandg(m,n,A,m)
      call crandg(m,1,u,m)
      call crandg(n,1,v,n)
c generate QR decomposition
      call cqrgen(m,n,A,m,Q,m,R,m)
c update A
      call cgerc(m,n,(1e0,0e0),u,1,v,1,A,m)
c update the QR decomposition
      k = m
      if (ec == 1) k = n
      call cqr1up(m,n,k,Q,m,R,m,u,v,wrk,rwrk)
c check result
      call cqrchk(m,n,k,A,m,Q,m,R,m)

      end subroutine

      subroutine ztest(m,n,ec)
      integer m,n,ec
      double complex A(m,max(m,n)),Q(m,m),R(m,n),u(m),v(n),wrk(m)
      double precision rwrk(m)
      external zrandg,zqrgen,zgerc,zqr1up,zqrchk
      integer k
c set up random matrix & vectors
      call zrandg(m,n,A,m)
      call zrandg(m,1,u,m)
      call zrandg(n,1,v,n)
c generate QR decomposition
      call zqrgen(m,n,A,m,Q,m,R,m)
c update A
      call zgerc(m,n,(1d0,0d0),u,1,v,1,A,m)
c update the QR decomposition
      k = m
      if (ec == 1) k = n
      call zqr1up(m,n,k,Q,m,R,m,u,v,wrk,rwrk)
c check result
      call zqrchk(m,n,k,A,m,Q,m,R,m)

      end subroutine
