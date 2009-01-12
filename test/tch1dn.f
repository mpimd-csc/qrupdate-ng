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
      program tch1dn
      integer n

      write (*,*)
      write (*,*) 'testing Cholesky rank-1 downdate routines.'
      write (*,*) 'All residual errors are expected to be small.'
      write (*,*)

      n = 50
      write (*,*) 'sch1dn test:'
      call stest(n)
      write (*,*) 'dch1dn test:'
      call dtest(n)
      write (*,*) 'cch1dn test:'
      call ctest(n)
      write (*,*) 'zch1dn test:'
      call ztest(n)

      call pstats
      end program

      subroutine stest(n)
      integer n
      real A(n,n),R(n,n),u(n),wrk(2*n)
      external srandg,scopy,schgen,sch1up,schchk
      integer info
c set up random matrix & vectors
      call srandg(n,n,A,n)
      call srandg(n,1,u,n)
      call scopy(n,u,1,wrk,1)
c generate A'*A and its Cholesky decomposition
      call schgen(n,A,n,R,n)
c update the Cholesky decomposition
      call sch1up(n,R,n,u,wrk(1+n))
c downdate it back
      call sch1dn(n,R,n,wrk,wrk(1+n),info)
c check result
      call schchk(n,A,n,R,n)

      end subroutine

      subroutine dtest(n)
      integer n
      double precision A(n,n),R(n,n),u(n),wrk(2*n)
      external drandg,dcopy,dchgen,dch1up,dchchk
      integer info
c set up random matrix & vectors
      call drandg(n,n,A,n)
      call drandg(n,1,u,n)
      call dcopy(n,u,1,wrk,1)
c generate A'*A and its Cholesky decomposition
      call dchgen(n,A,n,R,n)
c update the Cholesky decomposition
      call dch1up(n,R,n,u,wrk(1+n))
c downdate it back
      call dch1dn(n,R,n,wrk,wrk(1+n),info)
c check result
      call dchchk(n,A,n,R,n)

      end subroutine

      subroutine ctest(n)
      integer n
      complex A(n,n),R(n,n),u(n),wrk(n)
      real rwrk(n)
      external crandg,ccopy,cchgen,cch1up,cchchk
      integer info
c set up random matrix & vectors
      call crandg(n,n,A,n)
      call crandg(n,1,u,n)
      call ccopy(n,u,1,wrk,1)
c generate A'*A and its Cholesky decomposition
      call cchgen(n,A,n,R,n)
c update the Cholesky decomposition
      call cch1up(n,R,n,u,rwrk)
c downdate it back
      call cch1dn(n,R,n,wrk,rwrk,info)
c check result
      call cchchk(n,A,n,R,n)

      end subroutine

      subroutine ztest(n)
      integer n
      double complex A(n,n),R(n,n),u(n),wrk(n)
      double precision rwrk(n)
      external zrandg,zcopy,zchgen,zch1up,zchchk
      integer info
c set up random matrix & vectors
      call zrandg(n,n,A,n)
      call zrandg(n,1,u,n)
      call zcopy(n,u,1,wrk,1)
c generate A'*A and its Cholesky decomposition
      call zchgen(n,A,n,R,n)
c update the Cholesky decomposition
      call zch1up(n,R,n,u,rwrk)
c downdate it back
      call zch1dn(n,R,n,wrk,rwrk,info)
c check result
      call zchchk(n,A,n,R,n)

      end subroutine
