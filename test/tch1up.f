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
      program tch1up
      integer n

      write (*,*)
      write (*,*) 'testing Cholesky rank-1 update routines.'
      write (*,*) 'All residual errors are expected to be small.'
      write (*,*)

      n = 50
      write (*,*) 'sch1up test:'
      call stest(n)
      write (*,*) 'dch1up test:'
      call dtest(n)
      write (*,*) 'cch1up test:'
      call ctest(n)
      write (*,*) 'zch1up test:'
      call ztest(n)

      call pstats
      end program

      subroutine stest(n)
      integer n
      real A(n,n),R(n,n),u(n),wrk(n)
      external srandg,schgen,ssyr,sch1up,schchk
c set up random matrix & vectors
      call srandg(n,n,A,n)
      call srandg(n,1,u,n)
c generate A'*A and its Cholesky decomposition
      call schgen(n,A,n,R,n)
c update the matrix A
      call ssyr('U',n,1e0,u,1,A,n)
c update the Cholesky decomposition
      call sch1up(n,R,n,u,wrk)
c check result
      call schchk(n,A,n,R,n)

      end subroutine

      subroutine dtest(n)
      integer n
      double precision A(n,n),R(n,n),u(n),wrk(n)
      external drandg,dchgen,dsyr,dch1up,dchchk
c set up random matrix & vectors
      call drandg(n,n,A,n)
      call drandg(n,1,u,n)
c generate A'*A and its Cholesky decomposition
      call dchgen(n,A,n,R,n)
c update the matrix A
      call dsyr('U',n,1d0,u,1,A,n)
c update the Cholesky decomposition
      call dch1up(n,R,n,u,wrk)
c check result
      call dchchk(n,A,n,R,n)

      end subroutine

      subroutine ctest(n)
      integer n
      complex A(n,n),R(n,n),u(n)
      real rwrk(n)
      external crandg,cchgen,cher,cch1up,cchchk
c set up random matrix & vectors
      call crandg(n,n,A,n)
      call crandg(n,1,u,n)
c generate A'*A and its Cholesky decomposition
      call cchgen(n,A,n,R,n)
c update the matrix A
      call cher('U',n,1e0,u,1,A,n)
c update the Cholesky decomposition
      call cch1up(n,R,n,u,rwrk)
c check result
      call cchchk(n,A,n,R,n)

      end subroutine

      subroutine ztest(n)
      integer n
      double complex A(n,n),R(n,n),u(n)
      double precision rwrk(n)
      external zrandg,zchgen,zher,zch1up,zchchk
c set up random matrix & vectors
      call zrandg(n,n,A,n)
      call zrandg(n,1,u,n)
c generate A'*A and its Cholesky decomposition
      call zchgen(n,A,n,R,n)
c update the matrix A
      call zher('U',n,1d0,u,1,A,n)
c update the Cholesky decomposition
      call zch1up(n,R,n,u,rwrk)
c check result
      call zchchk(n,A,n,R,n)

      end subroutine
