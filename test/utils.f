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
      subroutine srandg(m,n,x,ldx)
      integer m,n,ldx
      real x(ldx,*)
      external slaruv
      integer seed(4),j
      common /xrand/ seed
      do j = 1,n
        call slaruv(seed,m,x(1,j))
      end do
      end subroutine

      subroutine drandg(m,n,x,ldx)
      integer m,n,ldx
      double precision x(ldx,*)
      external dlaruv
      integer seed(4),j
      common /xrand/ seed
      do j = 1,n
        call dlaruv(seed,m,x(1,j))
      end do
      end subroutine

      subroutine crandg(m,n,x,ldx)
      integer m,n,ldx
      complex x(ldx,*)
      external slaruv
      integer seed(4),j
      common /xrand/ seed
      do j = 1,n
        call slaruv(seed,2*m,x(1,j))
      end do
      end subroutine

      subroutine zrandg(m,n,x,ldx)
      integer m,n,ldx
      double complex x(ldx,*)
      external dlaruv
      integer seed(4),j
      common /xrand/ seed
      do j = 1,n
        call dlaruv(seed,2*m,x(1,j))
      end do
      end subroutine

      block data xrandi
      integer seed(4)
      common /xrand/ seed
      data seed /4*3/
      end block data

      subroutine sqrgen(m,n,A,lda,Q,ldq,R,ldr)
      integer m,n,lda,ldq,ldr
      real A(lda,n),Q(ldq,m),R(ldr,n)
      real work(max(m,n)),tau(min(m,n))
      integer info,i,j
      external slacpy,sgeqrf,sorgqr
      if (m == 0) return
      call slacpy('0',m,n,A,lda,R,ldr)
      call sgeqrf(m,n,R,ldr,tau,work,max(m,n),info)
      do i = 1,n
        do j = i+1,m
          Q(j,i) = R(j,i)
          R(j,i) = 0e0
        end do
      end do
      call sorgqr(m,m,min(m,n),Q,ldq,tau,work,max(m,n),info)
      end subroutine

      subroutine dqrgen(m,n,A,lda,Q,ldq,R,ldr)
      integer m,n,lda,ldq,ldr
      double precision A(lda,n),Q(ldq,m),R(ldr,n)
      double precision work(max(m,n)),tau(min(m,n))
      integer info,i,j
      external dlacpy,dgeqrf,dorgqr
      if (m == 0) return
      call dlacpy('0',m,n,A,lda,R,ldr)
      call dgeqrf(m,n,R,ldr,tau,work,max(m,n),info)
      do i = 1,n
        do j = i+1,m
          Q(j,i) = R(j,i)
          R(j,i) = 0d0
        end do
      end do
      call dorgqr(m,m,min(m,n),Q,ldq,tau,work,max(m,n),info)
      end subroutine

      subroutine cqrgen(m,n,A,lda,Q,ldq,R,ldr)
      integer m,n,lda,ldq,ldr
      complex A(lda,n),Q(ldq,m),R(ldr,n)
      complex work(max(m,n)),tau(min(m,n))
      integer info,i,j
      external clacpy,cgeqrf,cungqr
      if (m == 0) return
      call clacpy('0',m,n,A,lda,R,ldr)
      call cgeqrf(m,n,R,ldr,tau,work,max(m,n),info)
      do i = 1,n
        do j = i+1,m
          Q(j,i) = R(j,i)
          R(j,i) = 0e0
        end do
      end do
      call cungqr(m,m,min(m,n),Q,ldq,tau,work,max(m,n),info)
      end subroutine

      subroutine zqrgen(m,n,A,lda,Q,ldq,R,ldr)
      integer m,n,lda,ldq,ldr
      double complex A(lda,n),Q(ldq,m),R(ldr,n)
      double complex work(max(m,n)),tau(min(m,n))
      integer info,i,j
      external zlacpy,zgeqrf,zungqr
      if (m == 0) return
      call zlacpy('0',m,n,A,lda,R,ldr)
      call zgeqrf(m,n,R,ldr,tau,work,max(m,n),info)
      do i = 1,n
        do j = i+1,m
          Q(j,i) = R(j,i)
          R(j,i) = 0d0
        end do
      end do
      call zungqr(m,m,min(m,n),Q,ldq,tau,work,max(m,n),info)
      end subroutine

      subroutine smdump(name,m,n,A,lda)
      character(*) name
      integer m,n,lda
      real A(lda,n)
      integer i,j
      write (*,1001) name
      do i = 1,m
        do j = 1,n
          write(*,1002) A(i,j)
        end do
        write(*,*)
      end do

 1001 format (A,' = ')
 1002 format (1x,F6.3,$)
      end subroutine

      subroutine dmdump(name,m,n,A,lda)
      character(*) name
      integer m,n,lda
      double precision A(lda,n)
      integer i,j
      write (*,1001) name
      do i = 1,m
        do j = 1,n
          write(*,1002) A(i,j)
        end do
        write(*,*)
      end do

 1001 format (A,' = ')
 1002 format (1x,F6.3,$)
      end subroutine

      subroutine cmdump(name,m,n,A,lda)
      character(*) name
      integer m,n,lda
      complex A(lda,n)
      integer i,j
      write (*,1001) name
      do i = 1,m
        do j = 1,n
          write(*,1002) A(i,j)
        end do
        write(*,*)
      end do

 1001 format (A,' = ')
 1002 format (1x,F6.3,SP,F6.3,'i',$)
      end subroutine

      subroutine zmdump(name,m,n,A,lda)
      character(*) name
      integer m,n,lda
      double complex A(lda,n)
      integer i,j
      write (*,1001) name
      do i = 1,m
        do j = 1,n
          write(*,1002) A(i,j)
        end do
        write(*,*)
      end do

 1001 format (A,' = ')
 1002 format (1x,F6.3,SP,F6.3,'i',$)
      end subroutine

      subroutine sqrchk(m,n,k,A,lda,Q,ldq,R,ldr)
      integer m,n,k,lda,ldq,ldr
      real A(lda,max(n,k)),Q(ldq,k),R(ldr,n)
      real rnrm,slange,slansy
      external sgemm,ssyrk,slange,slansy
      real wrk(m)
      integer i

c get residual
      call sgemm('N','N',m,n,k,-1e0,Q,ldq,R,ldr,1e0,A,lda)
c get frobenius norm
      rnrm = slange('M',m,n,A,lda)
      write(*,1001) rnrm
c form Q'*Q - I
      call ssyrk('U','T',k,m,1e0,Q,ldq,0e0,A,lda)
      do i = 1,k
        A(i,i) = A(i,i) - 1e0
      end do
c get frobenius norm
      rnrm = slansy('M','U',k,A,lda,wrk)
      write(*,1002) rnrm
      return

 1001 format(6x,'residual error = ',10x,E21.12)
 1002 format('orth. residual error = ',10x,E21.12)
      end subroutine

      subroutine dqrchk(m,n,k,A,lda,Q,ldq,R,ldr)
      integer m,n,k,lda,ldq,ldr
      double precision A(lda,max(n,k)),Q(ldq,k),R(ldr,n)
      double precision rnrm,dlange,dlansy
      external dgemm,dsyrk,dlange,dlansy
      double precision wrk(m)
      integer i

c get residual
      call dgemm('N','N',m,n,k,-1d0,Q,ldq,R,ldr,1d0,A,lda)
c get frobenius norm
      rnrm = dlange('M',m,n,A,lda)
      write(*,1001) rnrm
c form Q'*Q - I
      call dsyrk('U','T',k,m,1d0,Q,ldq,0d0,A,lda)
      do i = 1,k
        A(i,i) = A(i,i) - 1d0
      end do
c get frobenius norm
      rnrm = dlansy('M','U',k,A,lda,wrk)
      write(*,1002) rnrm
      return

 1001 format(6x,'residual error = ',10x,E21.12)
 1002 format('orth. residual error = ',10x,E21.12)
      end subroutine

      subroutine cqrchk(m,n,k,A,lda,Q,ldq,R,ldr)
      integer m,n,k,lda,ldq,ldr
      complex A(lda,max(n,k)),Q(ldq,k),R(ldr,n)
      real rnrm,clange,clanhe
      external cgemm,csyrk,clange,clanhe
      real wrk(m)
      integer i

c get residual
      call cgemm('N','N',m,n,k,-(1e0,0e0),Q,ldq,R,ldr,(1e0,0e0),A,lda)
c get frobenius norm
      rnrm = clange('M',m,n,A,lda)
      write(*,1001) rnrm
c form Q'*Q - I
      call cherk('U','C',k,m,(1e0,0e0),Q,ldq,(0e0,0e0),A,lda)
      do i = 1,k
        A(i,i) = A(i,i) - 1e0
      end do
c get frobenius norm
      rnrm = clanhe('M','U',k,A,lda,wrk)
      write(*,1002) rnrm
      return

 1001 format(6x,'residual error = ',10x,E21.12)
 1002 format('orth. residual error = ',10x,E21.12)
      end subroutine

      subroutine zqrchk(m,n,k,A,lda,Q,ldq,R,ldr)
      integer m,n,k,lda,ldq,ldr
      double complex A(lda,max(n,k)),Q(ldq,k),R(ldr,n)
      double precision rnrm,zlange,zlanhe
      external zgemm,zsyrk,zlange,zlanhe
      double precision wrk(m)
      integer i

c get residual
      call zgemm('N','N',m,n,k,-(1d0,0d0),Q,ldq,R,ldr,(1d0,0d0),A,lda)
c get frobenius norm
      rnrm = zlange('M',m,n,A,lda)
      write(*,1001) rnrm
c form Q'*Q - I
      call zherk('U','C',k,m,(1d0,0d0),Q,ldq,(0d0,0d0),A,lda)
      do i = 1,k
        A(i,i) = A(i,i) - 1d0
      end do
c get frobenius norm
      rnrm = zlanhe('M','U',k,A,lda,wrk)
      write(*,1002) rnrm
      return

 1001 format(6x,'residual error = ',10x,E21.12)
 1002 format('orth. residual error = ',10x,E21.12)
      end subroutine

      subroutine schgen(n,A,lda,R,ldr)
      integer n,lda,ldr
      real A(lda,n),R(ldr,n)
      external ssyrk,slacpy,spotrf
      integer i,j,info
      call ssyrk('U','T',n,n,1e0,A,lda,0e0,R,ldr)
c augment to ensure strict positivity, zero below diag
      do i = 1,n
        R(i,i) = R(i,i) + 1e-3
c zero below diagonal
        do j = i+1,n
          R(j,i) = 0e0
        end do
      end do
      call slacpy('U',n,n,R,ldr,A,lda)
c symmetrize A
      do i = 1,n-1
        do j = i+1,n
          A(j,i) = A(i,j)
        end do
      end do
      call spotrf('U',n,R,ldr,info)
      if (info /= 0) stop 'fatal:error generating positive matrix'
      end subroutine

      subroutine dchgen(n,A,lda,R,ldr)
      integer n,lda,ldr
      double precision A(lda,n),R(ldr,n)
      external dsyrk,dlacpy,dpotrf
      integer i,j,info
      call dsyrk('U','T',n,n,1d0,A,lda,0d0,R,ldr)
c augment to ensure strict positivity
      do i = 1,n
        R(i,i) = R(i,i) + 1d-3
c zero below diagonal
        do j = i+1,n
          R(j,i) = 0d0
        end do
      end do
      call dlacpy('U',n,n,R,ldr,A,lda)
c symmetrize A
      do i = 1,n-1
        do j = i+1,n
          A(j,i) = A(i,j)
        end do
      end do
      call dpotrf('U',n,R,ldr,info)
      if (info /= 0) stop 'fatal:error generating positive matrix'
      end subroutine

      subroutine cchgen(n,A,lda,R,ldr)
      integer n,lda,ldr
      complex A(lda,n),R(ldr,n)
      external cherk,clacpy,cpotrf
      integer i,j,info
      call cherk('U','C',n,n,(1e0,0e0),A,lda,(0e0,0e0),R,ldr)
c augment to ensure strict positivity
      do i = 1,n
        R(i,i) = R(i,i) + 1e-3
c zero below diagonal
        do j = i+1,n
          R(j,i) = 0e0
        end do
      end do
      call clacpy('U',n,n,R,ldr,A,lda)
c symmetrize A
      do i = 1,n-1
        do j = i+1,n
          A(j,i) = conjg(A(i,j))
        end do
      end do
      call cpotrf('U',n,R,ldr,info)
      if (info /= 0) stop 'fatal:error generating positive matrix'
      end subroutine

      subroutine zchgen(n,A,lda,R,ldr)
      integer n,lda,ldr
      double complex A(lda,n),R(ldr,n)
      external zherk,zlacpy,zpotrf
      integer i,j,info
      call zherk('U','C',n,n,(1d0,0d0),A,lda,(0d0,0d0),R,ldr)
c augment to ensure strict positivity
      do i = 1,n
        R(i,i) = R(i,i) + 1d-3
c zero below diagonal
        do j = i+1,n
          R(j,i) = 0d0
        end do
      end do
      call zlacpy('U',n,n,R,ldr,A,lda)
c symmetrize A
      do i = 1,n-1
        do j = i+1,n
          A(j,i) = conjg(A(i,j))
        end do
      end do
      call zpotrf('U',n,R,ldr,info)
      if (info /= 0) stop 'fatal:error generating positive matrix'
      end subroutine

      subroutine schchk(n,A,lda,R,ldr)
      integer n,lda,ldr
      real A(lda,n),R(ldr,n)
      real rnrm,slansy
      external ssyrk,slansy
      real wrk(n)
      integer i,j

c zero lower triangle of R
      do j = 1,n-1
        do i = j+1,n
          R(i,j) = 0e0
        end do
      end do
c form A - R'*R
      call ssyrk('U','T',n,n,1e0,R,ldr,-1e0,A,lda)
c get frobenius norm
      rnrm = slansy('M','U',n,A,lda,wrk)
      write(*,1001) rnrm
      return

 1001 format(6x,'residual error = ',10x,E21.12)
      end subroutine

      subroutine dchchk(n,A,lda,R,ldr)
      integer n,lda,ldr
      double precision A(lda,n),R(ldr,n)
      double precision rnrm,dlansy
      external dsyrk,dlansy
      double precision wrk(n)
      integer i,j

c zero lower triangle of R
      do j = 1,n-1
        do i = j+1,n
          R(i,j) = 0d0
        end do
      end do
c form A - R'*R
      call dsyrk('U','T',n,n,1d0,R,ldr,-1d0,A,lda)
c get frobenius norm
      rnrm = dlansy('M','U',n,A,lda,wrk)
      write(*,1001) rnrm
      return

 1001 format(6x,'residual error = ',10x,E21.12)
      end subroutine

      subroutine cchchk(n,A,lda,R,ldr)
      integer n,lda,ldr
      complex A(lda,n),R(ldr,n)
      real rnrm,clanhe
      external cherk,clanhe
      real wrk(n)
      integer i,j

c zero lower triangle of R
      do j = 1,n-1
        do i = j+1,n
          R(i,j) = 0e0
        end do
      end do
c form A - R'*R
      call cherk('U','C',n,n,(1e0,0e0),R,ldr,(-1e0,0e0),A,lda)
c get frobenius norm
      rnrm = clanhe('M','U',n,A,lda,wrk)
      write(*,1001) rnrm
      return

 1001 format(6x,'residual error = ',10x,E21.12)
      end subroutine

      subroutine zchchk(n,A,lda,R,ldr)
      integer n,lda,ldr
      double complex A(lda,n),R(ldr,n)
      double precision rnrm,zlanhe
      external zherk,zlanhe
      double precision wrk(n)
      integer i,j

c zero lower triangle of R
      do j = 1,n-1
        do i = j+1,n
          R(i,j) = 0d0
        end do
      end do
c form A - R'*R
      call zherk('U','C',n,n,(1d0,0d0),R,ldr,(-1d0,0d0),A,lda)
c get frobenius norm
      rnrm = zlanhe('M','U',n,A,lda,wrk)
      write(*,1001) rnrm
      return

 1001 format(6x,'residual error = ',10x,E21.12)
      end subroutine
