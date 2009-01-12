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
      program tqrinc
      integer m,n,j

      write (*,*)
      write (*,*) 'testing QR column insert routines.'
      write (*,*) 'All residual errors are expected to be small.'
      write (*,*)

      m = 60
      n = 40
      j = 28
      write (*,*) 'sqrinc test (full factorization):'
      call stest(m,n,j,0)
      write (*,*) 'dqrinc test (full factorization):'
      call dtest(m,n,j,0)
      write (*,*) 'cqrinc test (full factorization):'
      call ctest(m,n,j,0)
      write (*,*) 'zqrinc test (full factorization):'
      call ztest(m,n,j,0)

      write (*,*) 'sqrinc test (economized factorization):'
      call stest(m,n,j,1)
      write (*,*) 'dqrinc test (economized factorization):'
      call dtest(m,n,j,1)
      write (*,*) 'cqrinc test (economized factorization):'
      call ctest(m,n,j,1)
      write (*,*) 'zqrinc test (economized factorization):'
      call ztest(m,n,j,1)

      call pstats
      end program

      subroutine stest(m,n,j,ec)
      integer m,n,j,ec
      real A(m,max(m,n+1)),Q(m,m),R(m,n+1),u(m),wrk(m)
      external srandg,sqrgen,scopy,sqrinc,sqrchk
      integer k,i
c set up random matrix & vector
      call srandg(m,n,A,m)
      call srandg(m,1,u,m)
c generate QR decomposition
      call sqrgen(m,n,A,m,Q,m,R,m)
c update A
      do i = n,j,-1
        call scopy(m,A(1,i),1,A(1,i+1),1)
      end do
      call scopy(m,u,1,A(1,j),1)
c update the QR decomposition
      k = m
      if (ec == 1) k = n
      call sqrinc(m,n,k,Q,m,R,m,j,u,wrk)
c check result
      if (ec == 1) k = n+1
      call sqrchk(m,n+1,k,A,m,Q,m,R,m)

      end subroutine

      subroutine dtest(m,n,j,ec)
      integer m,n,j,ec
      double precision A(m,max(m,n+1)),Q(m,m),R(m,n+1),u(m),wrk(m)
      external drandg,dqrgen,dcopy,dqrinc,dqrchk
      integer k,i
c set up random matrix & vector
      call drandg(m,n,A,m)
      call drandg(m,1,u,m)
c generate QR decomposition
      call dqrgen(m,n,A,m,Q,m,R,m)
c update A
      do i = n,j,-1
        call dcopy(m,A(1,i),1,A(1,i+1),1)
      end do
      call dcopy(m,u,1,A(1,j),1)
c update the QR decomposition
      k = m
      if (ec == 1) k = n
      call dqrinc(m,n,k,Q,m,R,m,j,u,wrk)
c check result
      if (ec == 1) k = n+1
      call dqrchk(m,n+1,k,A,m,Q,m,R,m)

      end subroutine

      subroutine ctest(m,n,j,ec)
      integer m,n,j,ec
      complex A(m,max(m,n+1)),Q(m,m),R(m,n+1),u(m),wrk(m)
      external crandg,cqrgen,ccopy,cqrinc,cqrchk
      integer k,i
c set up random matrix & vector
      call crandg(m,n,A,m)
      call crandg(m,1,u,m)
c generate QR decomposition
      call cqrgen(m,n,A,m,Q,m,R,m)
c update A
      do i = n,j,-1
        call ccopy(m,A(1,i),1,A(1,i+1),1)
      end do
      call ccopy(m,u,1,A(1,j),1)
c update the QR decomposition
      k = m
      if (ec == 1) k = n
      call cqrinc(m,n,k,Q,m,R,m,j,u,wrk)
c check result
      if (ec == 1) k = n+1
      call cqrchk(m,n+1,k,A,m,Q,m,R,m)

      end subroutine

      subroutine ztest(m,n,j,ec)
      integer m,n,j,ec
      double complex A(m,max(m,n+1)),Q(m,m),R(m,n+1),u(m),wrk(m)
      external zrandg,zqrgen,zcopy,zqrinc,zqrchk
      integer k,i
c set up random matrix & vector
      call zrandg(m,n,A,m)
      call zrandg(m,1,u,m)
c generate QR decomposition
      call zqrgen(m,n,A,m,Q,m,R,m)
c update A
      do i = n,j,-1
        call zcopy(m,A(1,i),1,A(1,i+1),1)
      end do
      call zcopy(m,u,1,A(1,j),1)
c update the QR decomposition
      k = m
      if (ec == 1) k = n
      call zqrinc(m,n,k,Q,m,R,m,j,u,wrk)
c check result
      if (ec == 1) k = n+1
      call zqrchk(m,n+1,k,A,m,Q,m,R,m)

      end subroutine
