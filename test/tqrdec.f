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
      program tqrdec
      integer m,n,j

      write (*,*)
      write (*,*) 'testing QR column delete routines.'
      write (*,*) 'All residual errors are expected to be small.'
      write (*,*)

      m = 60
      n = 40
      j = 12
      write (*,*) 'sqrdec test (full factorization):'
      call stest(m,n,j,0)
      write (*,*) 'dqrdec test (full factorization):'
      call dtest(m,n,j,0)
      write (*,*) 'cqrdec test (full factorization):'
      call ctest(m,n,j,0)
      write (*,*) 'zqrdec test (full factorization):'
      call ztest(m,n,j,0)

      write (*,*) 'sqrdec test (economized factorization):'
      call stest(m,n,j,1)
      write (*,*) 'dqrdec test (economized factorization):'
      call dtest(m,n,j,1)
      write (*,*) 'cqrdec test (economized factorization):'
      call ctest(m,n,j,1)
      write (*,*) 'zqrdec test (economized factorization):'
      call ztest(m,n,j,1)

      call pstats
      end program

      subroutine stest(m,n,j,ec)
      integer m,n,j,ec
      real A(m,max(m,n)),Q(m,m),R(m,n),wrk(m)
      external srandg,sqrgen,scopy,sqrdec,sqrchk
      integer k,i
c set up random matrix & vector
      call srandg(m,n,A,m)
c generate QR decomposition
      call sqrgen(m,n,A,m,Q,m,R,m)
c update A
      do i = j,n-1
        call scopy(m,A(1,i+1),1,A(1,i),1)
      end do
c update the QR decomposition
      k = m
      if (ec == 1) k = n
      call sqrdec(m,n,k,Q,m,R,m,j,wrk)
c check result
      if (ec == 1) k = n+1
      call sqrchk(m,n-1,k,A,m,Q,m,R,m)

      end subroutine

      subroutine dtest(m,n,j,ec)
      integer m,n,j,ec
      double precision A(m,max(m,n)),Q(m,m),R(m,n),wrk(m)
      external drandg,dqrgen,dcopy,dqrdec,dqrchk
      integer k,i
c set up random matrix & vector
      call drandg(m,n,A,m)
c generate QR decomposition
      call dqrgen(m,n,A,m,Q,m,R,m)
c update A
      do i = j,n-1
        call dcopy(m,A(1,i+1),1,A(1,i),1)
      end do
c update the QR decomposition
      k = m
      if (ec == 1) k = n
      call dqrdec(m,n,k,Q,m,R,m,j,wrk)
c check result
      if (ec == 1) k = n+1
      call dqrchk(m,n-1,k,A,m,Q,m,R,m)

      end subroutine

      subroutine ctest(m,n,j,ec)
      integer m,n,j,ec
      complex A(m,max(m,n)),Q(m,m),R(m,n),wrk(m)
      external crandg,cqrgen,ccopy,cqrdec,cqrchk
      integer k,i
c set up random matrix & vector
      call crandg(m,n,A,m)
c generate QR decomposition
      call cqrgen(m,n,A,m,Q,m,R,m)
c update A
      do i = j,n-1
        call ccopy(m,A(1,i+1),1,A(1,i),1)
      end do
c update the QR decomposition
      k = m
      if (ec == 1) k = n
      call cqrdec(m,n,k,Q,m,R,m,j,wrk)
c check result
      if (ec == 1) k = n+1
      call cqrchk(m,n-1,k,A,m,Q,m,R,m)

      end subroutine

      subroutine ztest(m,n,j,ec)
      integer m,n,j,ec
      double complex A(m,max(m,n)),Q(m,m),R(m,n),wrk(m)
      external zrandg,zqrgen,zcopy,zqrdec,zqrchk
      integer k,i
c set up random matrix & vector
      call zrandg(m,n,A,m)
c generate QR decomposition
      call zqrgen(m,n,A,m,Q,m,R,m)
c update A
      do i = j,n-1
        call zcopy(m,A(1,i+1),1,A(1,i),1)
      end do
c update the QR decomposition
      k = m
      if (ec == 1) k = n
      call zqrdec(m,n,k,Q,m,R,m,j,wrk)
c check result
      if (ec == 1) k = n+1
      call zqrchk(m,n-1,k,A,m,Q,m,R,m)

      end subroutine
