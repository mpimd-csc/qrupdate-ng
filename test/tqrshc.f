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
      program tqrshc
      integer m,n,i,j

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
      integer m,n,i,j,ec
      real A(m,max(m,n)),Q(m,m),R(m,n),wrk(2*m)
      external srandg,sqrgen,scopy,sqrshc,sqrchk
      integer k
c set up random matrix & vector
      call srandg(m,n,A,m)
c generate QR decomposition
      call sqrgen(m,n,A,m,Q,m,R,m)
c update A
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
c update the QR decomposition
      k = m
      if (ec == 1) k = n
      call sqrshc(m,n,k,Q,m,R,m,i,j,wrk)
c check result
      call sqrchk(m,n,k,A,m,Q,m,R,m)

      end subroutine

      subroutine dtest(m,n,i,j,ec)
      integer m,n,i,j,ec
      double precision A(m,max(m,n)),Q(m,m),R(m,n),wrk(2*m)
      external drandg,dqrgen,dcopy,dqrshc,dqrchk
      integer k
c set up random matrix & vector
      call drandg(m,n,A,m)
c generate QR decomposition
      call dqrgen(m,n,A,m,Q,m,R,m)
c update A
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
c update the QR decomposition
      k = m
      if (ec == 1) k = n
      call dqrshc(m,n,k,Q,m,R,m,i,j,wrk)
c check result
      call dqrchk(m,n,k,A,m,Q,m,R,m)

      end subroutine

      subroutine ctest(m,n,i,j,ec)
      integer m,n,i,j,ec
      complex A(m,max(m,n)),Q(m,m),R(m,n),wrk(m)
      real rwrk(m)
      external crandg,cqrgen,ccopy,cqrshc,cqrchk
      integer k
c set up random matrix & vector
      call crandg(m,n,A,m)
c generate QR decomposition
      call cqrgen(m,n,A,m,Q,m,R,m)
c update A
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
c update the QR decomposition
      k = m
      if (ec == 1) k = n
      call cqrshc(m,n,k,Q,m,R,m,i,j,wrk,rwrk)
c check result
      call cqrchk(m,n,k,A,m,Q,m,R,m)

      end subroutine

      subroutine ztest(m,n,i,j,ec)
      integer m,n,i,j,ec
      double complex A(m,max(m,n)),Q(m,m),R(m,n),wrk(m)
      double precision rwrk(m)
      external zrandg,zqrgen,zcopy,zqrshc,zqrchk
      integer k
c set up random matrix & vector
      call zrandg(m,n,A,m)
c generate QR decomposition
      call zqrgen(m,n,A,m,Q,m,R,m)
c update A
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
c update the QR decomposition
      k = m
      if (ec == 1) k = n
      call zqrshc(m,n,k,Q,m,R,m,i,j,wrk,rwrk)
c check result
      call zqrchk(m,n,k,A,m,Q,m,R,m)

      end subroutine
