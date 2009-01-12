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
      program tqrder
      integer m,n,j

      write (*,*)
      write (*,*) 'testing QR row delete routines.'
      write (*,*) 'All residual errors are expected to be small.'
      write (*,*)

      m = 60
      n = 40
      j = 30
      write (*,*) 'sqrder test (full factorization):'
      call stest(m,n,j)
      write (*,*) 'dqrder test (full factorization):'
      call dtest(m,n,j)
      write (*,*) 'cqrder test (full factorization):'
      call ctest(m,n,j)
      write (*,*) 'zqrder test (full factorization):'
      call ztest(m,n,j)

      call pstats
      end program

      subroutine stest(m,n,j)
      integer m,n,j
      real A(m,max(m,n)),Q(m,m),R(m,n),wrk(2*m)
      external srandg,sqrgen,scopy,sqrder,sqrchk
      integer i
c set up random matrix & vector
      call srandg(m,n,A,m)
c generate QR decomposition
      call sqrgen(m,n,A,m,Q,m,R,m)
c update A
      do i = j,m-1
        call scopy(n,A(i+1,1),m,A(i,1),m)
      end do
c update the QR decomposition
      call sqrder(m,n,Q,m,R,m,j,wrk)
c check result
      call sqrchk(m-1,n,m-1,A,m,Q,m,R,m)

      end subroutine

      subroutine dtest(m,n,j)
      integer m,n,j
      double precision A(m,max(m,n)),Q(m,m),R(m,n),wrk(2*m)
      external drandg,dqrgen,dcopy,dqrder,dqrchk
      integer i
c set up random matrix & vector
      call drandg(m,n,A,m)
c generate QR decomposition
      call dqrgen(m,n,A,m,Q,m,R,m)
c update A
      do i = j,m-1
        call dcopy(n,A(i+1,1),m,A(i,1),m)
      end do
c update the QR decomposition
      call dqrder(m,n,Q,m,R,m,j,wrk)
c check result
      call dqrchk(m-1,n,m-1,A,m,Q,m,R,m)

      end subroutine

      subroutine ctest(m,n,j)
      integer m,n,j
      complex A(m,max(m,n)),Q(m,m),R(m,n),wrk(m)
      real rwrk(m)
      external crandg,cqrgen,ccopy,cqrder,cqrchk
      integer i
c set up random matrix & vector
      call crandg(m,n,A,m)
c generate QR decomposition
      call cqrgen(m,n,A,m,Q,m,R,m)
c update A
      do i = j,m-1
        call ccopy(n,A(i+1,1),m,A(i,1),m)
      end do
c update the QR decomposition
      call cqrder(m,n,Q,m,R,m,j,wrk,rwrk)
c check result
      call cqrchk(m-1,n,m-1,A,m,Q,m,R,m)

      end subroutine

      subroutine ztest(m,n,j)
      integer m,n,j
      double complex A(m,max(m,n)),Q(m,m),R(m,n),wrk(m)
      double precision rwrk(m)
      external zrandg,zqrgen,zcopy,zqrder,zqrchk
      integer i
c set up random matrix & vector
      call zrandg(m,n,A,m)
c generate QR decomposition
      call zqrgen(m,n,A,m,Q,m,R,m)
c update A
      do i = j,m-1
        call zcopy(n,A(i+1,1),m,A(i,1),m)
      end do
c update the QR decomposition
      call zqrder(m,n,Q,m,R,m,j,wrk,rwrk)
c check result
      call zqrchk(m-1,n,m-1,A,m,Q,m,R,m)

      end subroutine
