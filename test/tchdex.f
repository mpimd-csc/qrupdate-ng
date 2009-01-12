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
      program tchdex
      integer n,j

      write (*,*)
      write (*,*) 'testing Cholesky symmetric delete routines.'
      write (*,*) 'All residual errors are expected to be small.'
      write (*,*)

      n = 50
      j = 15
      write (*,*) 'schdex test:'
      call stest(n,j)
      write (*,*) 'dchdex test:'
      call dtest(n,j)
      write (*,*) 'cchdex test:'
      call ctest(n,j)
      write (*,*) 'zchdex test:'
      call ztest(n,j)

      call pstats
      end program

      subroutine stest(n,j)
      integer n,j
      real A(n,n),R(n,n),wrk(n)
      external srandg,schgen,schdex,schchk
      integer i,k
c set up random matrix & vectors
      call srandg(n,n,A,n)
c generate A'*A and its Cholesky decomposition
      call schgen(n,A,n,R,n)
c update the matrix A
      do k = j,n-1
        do i = 1,j-1
          A(i,k) = A(i,k+1)
        end do
        do i = j,k
          A(i,k) = A(i+1,k+1)
        end do
      end do
c update the Cholesky decomposition
      call schdex(n,R,n,j,wrk)
c check result
      call schchk(n-1,A,n,R,n)

      end subroutine

      subroutine dtest(n,j)
      integer n,j
      double precision A(n,n),R(n,n),wrk(n)
      external drandg,dchgen,dchdex,dchchk
      integer i,k
c set up random matrix & vectors
      call drandg(n,n,A,n)
c generate A'*A and its Cholesky decomposition
      call dchgen(n,A,n,R,n)
c update the matrix A
      do k = j,n-1
        do i = 1,j-1
          A(i,k) = A(i,k+1)
        end do
        do i = j,k
          A(i,k) = A(i+1,k+1)
        end do
      end do
c update the Cholesky decomposition
      call dchdex(n,R,n,j,wrk)
c check result
      call dchchk(n-1,A,n,R,n)

      end subroutine

      subroutine ctest(n,j)
      integer n,j
      complex A(n,n),R(n,n)
      real rwrk(n)
      external crandg,cchgen,cchdex,cchchk
      integer i,k
c set up random matrix & vectors
      call crandg(n,n,A,n)
c generate A'*A and its Cholesky decomposition
      call cchgen(n,A,n,R,n)
c update the matrix A
      do k = j,n-1
        do i = 1,j-1
          A(i,k) = A(i,k+1)
        end do
        do i = j,k
          A(i,k) = A(i+1,k+1)
        end do
      end do
c update the Cholesky decomposition
      call cchdex(n,R,n,j,rwrk)
c check result
      call cchchk(n-1,A,n,R,n)

      end subroutine

      subroutine ztest(n,j)
      integer n,j
      double complex A(n,n),R(n,n)
      double precision rwrk(n)
      external zrandg,zchgen,zchdex,zchchk
      integer i,k
c set up random matrix & vectors
      call zrandg(n,n,A,n)
c generate A'*A and its Cholesky decomposition
      call zchgen(n,A,n,R,n)
c update the matrix A
      do k = j,n-1
        do i = 1,j-1
          A(i,k) = A(i,k+1)
        end do
        do i = j,k
          A(i,k) = A(i+1,k+1)
        end do
      end do
c update the Cholesky decomposition
      call zchdex(n,R,n,j,rwrk)
c check result
      call zchchk(n-1,A,n,R,n)

      end subroutine
