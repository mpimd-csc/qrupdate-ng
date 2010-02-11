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
      integer seed(4),j,k
      common /xrand/ seed
      do j = 1,n
        do k = 1,m,128
          call slaruv(seed,min(m-k+1,128),x(k,j))
        end do
      end do
      do j = 1,n
        do k = 1,m
          if (.not. (x(k,j) > 0e0 .and. x(k,j) < 1e0)) then
            stop 'slaruv produced invalid number'
          end if
        end do
      end do
      end subroutine

      subroutine drandg(m,n,x,ldx)
      integer m,n,ldx
      double precision x(ldx,*)
      external dlaruv
      integer seed(4),j,k
      common /xrand/ seed
      do j = 1,n
        do k = 1,m,128
          call dlaruv(seed,min(m-k+1,128),x(k,j))
        end do
      end do
      do j = 1,n
        do k = 1,m
          if (.not. (x(k,j) > 0d0 .and. x(k,j) < 1d0)) then
            stop 'dlaruv produced invalid number'
          end if
        end do
      end do
      end subroutine

      subroutine crandg(m,n,x,ldx)
      integer m,n,ldx
      complex x(ldx,*)
      external srandg
      call srandg(2*m,n,x,2*ldx)
      end subroutine

      subroutine zrandg(m,n,x,ldx)
      integer m,n,ldx
      double complex x(ldx,*)
      external srandg
      call drandg(2*m,n,x,2*ldx)
      end subroutine

      block data xrandi
      integer seed(4)
      common /xrand/ seed
      data seed /4*1/
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

      character*4 function spftol(rnrm)
      real rnrm,slamch
      external slamch
      common /stats/ passed,failed
      integer passed,failed
      if (rnrm < 2e2*slamch('p')) then
        spftol = 'PASS'
        passed = passed + 1
      else
        spftol = 'FAIL'
        failed = failed + 1
      end if
      end function

      character*4 function dpftol(rnrm)
      double precision rnrm,dlamch
      external dlamch
      common /stats/ passed,failed
      integer passed,failed
      if (rnrm < 2d2*dlamch('p')) then
        dpftol = 'PASS'
        passed = passed + 1
      else
        dpftol = 'FAIL'
        failed = failed + 1
      end if
      end function

      subroutine pstats
      common /stats/ passed,failed
      integer passed,failed

      write(*,1001) 
      write(*,1002) passed,failed
      write(*,*)
 1001 format(70('-'))
 1002 format(1x,'total:',5x,'PASSED',1x,I3,5x,'FAILED',1x,I3)
      end subroutine

      subroutine sqrchk(m,n,k,A,lda,Q,ldq,R,ldr)
      integer m,n,k,lda,ldq,ldr
      real A(lda,max(n,k)),Q(ldq,k),R(ldr,n)
      real rnrm,slange,slansy
      external sgemm,ssyrk,slange,slansy,spftol
      character*4 spftol
      real wrk(m)
      integer i

c get residual
      call sgemm('N','N',m,n,k,-1e0,Q,ldq,R,ldr,1e0,A,lda)
c get frobenius norm
      rnrm = slange('M',m,n,A,lda,wrk)
      write(*,1001) rnrm,spftol(rnrm)
c form Q'*Q - I
      call ssyrk('U','T',k,m,1e0,Q,ldq,0e0,A,lda)
      do i = 1,k
        A(i,i) = A(i,i) - 1e0
      end do
c get frobenius norm
      rnrm = slansy('M','U',k,A,lda,wrk)
      write(*,1002) rnrm,spftol(rnrm)
      return

 1001 format(6x,'residual error = ',10x,E21.12,5x,A6)
 1002 format('orth. residual error = ',10x,E21.12,5x,A6)
      end subroutine

      subroutine dqrchk(m,n,k,A,lda,Q,ldq,R,ldr)
      integer m,n,k,lda,ldq,ldr
      double precision A(lda,max(n,k)),Q(ldq,k),R(ldr,n)
      double precision rnrm,dlange,dlansy
      external dgemm,dsyrk,dlange,dlansy,dpftol
      character*4 dpftol
      double precision wrk(m)
      integer i

c get residual
      call dgemm('N','N',m,n,k,-1d0,Q,ldq,R,ldr,1d0,A,lda)
c get frobenius norm
      rnrm = dlange('M',m,n,A,lda,wrk)
      write(*,1001) rnrm,dpftol(rnrm)
c form Q'*Q - I
      call dsyrk('U','T',k,m,1d0,Q,ldq,0d0,A,lda)
      do i = 1,k
        A(i,i) = A(i,i) - 1d0
      end do
c get frobenius norm
      rnrm = dlansy('M','U',k,A,lda,wrk)
      write(*,1002) rnrm,dpftol(rnrm)
      return

 1001 format(6x,'residual error = ',10x,E21.12,5x,A6)
 1002 format('orth. residual error = ',10x,E21.12,5x,A6)
      end subroutine

      subroutine cqrchk(m,n,k,A,lda,Q,ldq,R,ldr)
      integer m,n,k,lda,ldq,ldr
      complex A(lda,max(n,k)),Q(ldq,k),R(ldr,n)
      real rnrm,clange,clanhe
      external cgemm,cherk,clange,clanhe,spftol
      character*4 spftol
      real wrk(m)
      integer i

c get residual
      call cgemm('N','N',m,n,k,-(1e0,0e0),Q,ldq,R,ldr,(1e0,0e0),A,lda)
c get frobenius norm
      rnrm = clange('M',m,n,A,lda,wrk)
      write(*,1001) rnrm,spftol(rnrm)
c form Q'*Q - I
      call cherk('U','C',k,m,1e0,Q,ldq,0e0,A,lda)
      do i = 1,k
        A(i,i) = A(i,i) - 1e0
      end do
c get frobenius norm
      rnrm = clanhe('M','U',k,A,lda,wrk)
      write(*,1002) rnrm,spftol(rnrm)
      return

 1001 format(6x,'residual error = ',10x,E21.12,5x,A6)
 1002 format('orth. residual error = ',10x,E21.12,5x,A6)
      end subroutine

      subroutine zqrchk(m,n,k,A,lda,Q,ldq,R,ldr)
      integer m,n,k,lda,ldq,ldr
      double complex A(lda,max(n,k)),Q(ldq,k),R(ldr,n)
      double precision rnrm,zlange,zlanhe
      external zgemm,zherk,zlange,zlanhe,dpftol
      character*4 dpftol
      double precision wrk(m)
      integer i

c get residual
      call zgemm('N','N',m,n,k,-(1d0,0d0),Q,ldq,R,ldr,(1d0,0d0),A,lda)
c get frobenius norm
      rnrm = zlange('M',m,n,A,lda,wrk)
      write(*,1001) rnrm,dpftol(rnrm)
c form Q'*Q - I
      call zherk('U','C',k,m,1d0,Q,ldq,0d0,A,lda)
      do i = 1,k
        A(i,i) = A(i,i) - 1d0
      end do
c get frobenius norm
      rnrm = zlanhe('M','U',k,A,lda,wrk)
      write(*,1002) rnrm,dpftol(rnrm)
      return

 1001 format(6x,'residual error = ',10x,E21.12,5x,A6)
 1002 format('orth. residual error = ',10x,E21.12,5x,A6)
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
      external cherk,clacpy,cpotrf,cdotc
      complex cdotc,Rii
      integer i,j,info
      call cherk('U','C',n,n,1e0,A,lda,0e0,R,ldr)
c augment to ensure strict positivity
      do i = 1,n
c CHERK is often buggy. We'll recompute the diagonal elements and
c possibly warn about the bug.
        Rii = cdotc (n, A(1,i), 1, A(1,i), 1)
        if (.not. abs (Rii - R(i,i)) < 1e-5 * abs(Rii)) then
          write (*,1001)
          write (*,1002) R(i,i), Rii
        endif
        R(i,i) = Rii + 1e-3
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
 1001 format ('WARNING: Possible bug in BLAS CHERK:')
 1002 format ('WARNING: CHERK computed diagonal element (',
     +F6.3,F6.3,'), direct computation gives (',F6.3,F6.3,').')
      end subroutine

      subroutine zchgen(n,A,lda,R,ldr)
      integer n,lda,ldr
      double complex A(lda,n),R(ldr,n)
      external zherk,zlacpy,zpotrf,zdotc
      double complex zdotc,Rii
      integer i,j,info
      call zherk('U','C',n,n,1d0,A,lda,0d0,R,ldr)
c augment to ensure strict positivity
      do i = 1,n
c ZHERK is often buggy. We'll recompute the diagonal elements and
c possibly warn about the bug.
        Rii = zdotc (n, A(1,i), 1, A(1,i), 1)
        if (.not. abs (Rii - R(i,i)) < 1d-10 * abs(Rii)) then
          write (*,1001)
          write (*,1002) R(i,i), Rii 
        endif
        R(i,i) = Rii + 1d-3
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
 1001 format ('WARNING: Possible bug in BLAS CHERK:')
 1002 format ('WARNING: ZHERK computed diagonal element (',
     +F6.3,F6.3,'), direct computation gives (',F6.3,F6.3,').')
      end subroutine

      subroutine schchk(n,A,lda,R,ldr)
      integer n,lda,ldr
      real A(lda,n),R(ldr,n)
      real rnrm,slansy
      external ssyrk,slansy,spftol
      character*4 spftol
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
      write(*,1001) rnrm,spftol(rnrm)
      return

 1001 format(6x,'residual error = ',10x,E21.12,5x,A6)
      end subroutine

      subroutine dchchk(n,A,lda,R,ldr)
      integer n,lda,ldr
      double precision A(lda,n),R(ldr,n)
      double precision rnrm,dlansy
      external dsyrk,dlansy,dpftol
      character*4 dpftol
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
      write(*,1001) rnrm,dpftol(rnrm)
      return

 1001 format(6x,'residual error = ',10x,E21.12,5x,A6)
      end subroutine

      subroutine cchchk(n,A,lda,R,ldr)
      integer n,lda,ldr
      complex A(lda,n),R(ldr,n)
      real rnrm,clanhe
      external cherk,clanhe,spftol
      character*4 spftol
      real wrk(n)
      integer i,j

c zero lower triangle of R
      do j = 1,n-1
        do i = j+1,n
          R(i,j) = 0e0
        end do
      end do
c form A - R'*R
      call cherk('U','C',n,n,1e0,R,ldr,-1e0,A,lda)
c get frobenius norm
      rnrm = clanhe('M','U',n,A,lda,wrk)
      write(*,1001) rnrm,spftol(rnrm)
      return

 1001 format(6x,'residual error = ',10x,E21.12,5x,A6)
      end subroutine

      subroutine zchchk(n,A,lda,R,ldr)
      integer n,lda,ldr
      double complex A(lda,n),R(ldr,n)
      double precision rnrm,zlanhe
      external zherk,zlanhe,dpftol
      character*4 dpftol
      double precision wrk(n)
      integer i,j

c zero lower triangle of R
      do j = 1,n-1
        do i = j+1,n
          R(i,j) = 0d0
        end do
      end do
c form A - R'*R
      call zherk('U','C',n,n,1d0,R,ldr,-1d0,A,lda)
c get frobenius norm
      rnrm = zlanhe('M','U',n,A,lda,wrk)
      write(*,1001) rnrm,dpftol(rnrm)
      return

 1001 format(6x,'residual error = ',10x,E21.12,5x,A6)
      end subroutine

      subroutine slugen(m,n,A,lda,L,ldl,R,ldr)
      integer m,n,lda,ldl,ldr
      real A(lda,n),L(ldl,min(m,n)),R(ldr,n)
      integer ipiv(min(m,n)),info,i,j
      external sswap,slacpy,sgetrf
      if (m >= n) then
        call slacpy('0',m,n,A,lda,L,ldl)
        call sgetrf(m,n,L,ldl,ipiv,info)
        call slacpy('U',m,n,L,ldl,R,ldr)
      else
        call slacpy('0',m,n,A,lda,R,ldr)
        call sgetrf(m,n,R,ldr,ipiv,info)
        call slacpy('L',m,n,R,ldr,L,ldl)
      end if
      do i = 1,min(m,n)
        do j = 1,i-1
          L(j,i) = 0e0
        end do
        L(i,i) = 1e0
      end do
c permute the orig matrix      
      do i = 1,min(m,n)
        j = ipiv(i)
        if (i /= j) then
          call sswap(n,A(i,1),lda,A(j,1),lda)
        end if
      end do
      end subroutine

      subroutine dlugen(m,n,A,lda,L,ldl,R,ldr)
      integer m,n,lda,ldl,ldr
      double precision A(lda,n),L(ldl,min(m,n)),R(ldr,n)
      integer ipiv(min(m,n)),info,i,j
      external dswap,dlacpy,dgetrf
      if (m >= n) then
        call dlacpy('0',m,n,A,lda,L,ldl)
        call dgetrf(m,n,L,ldl,ipiv,info)
        call dlacpy('U',m,n,L,ldl,R,ldr)
      else
        call dlacpy('0',m,n,A,lda,R,ldr)
        call dgetrf(m,n,R,ldr,ipiv,info)
        call dlacpy('L',m,n,R,ldr,L,ldl)
      end if
      do i = 1,min(m,n)
        do j = 1,i-1
          L(j,i) = 0d0
        end do
        L(i,i) = 1d0
      end do
c permute the orig matrix      
      do i = 1,min(m,n)
        j = ipiv(i)
        if (i /= j) then
          call dswap(n,A(i,1),lda,A(j,1),lda)
        end if
      end do
      end subroutine

      subroutine clugen(m,n,A,lda,L,ldl,R,ldr)
      integer m,n,lda,ldl,ldr
      complex A(lda,n),L(ldl,min(m,n)),R(ldr,n)
      integer ipiv(min(m,n)),info,i,j
      external cswap,clacpy,cgetrf
      if (m >= n) then
        call clacpy('0',m,n,A,lda,L,ldl)
        call cgetrf(m,n,L,ldl,ipiv,info)
        call clacpy('U',m,n,L,ldl,R,ldr)
      else
        call clacpy('0',m,n,A,lda,R,ldr)
        call cgetrf(m,n,R,ldr,ipiv,info)
        call clacpy('L',m,n,R,ldr,L,ldl)
      end if
      do i = 1,min(m,n)
        do j = 1,i-1
          L(j,i) = 0d0
        end do
        L(i,i) = 1e0
      end do
c permute the orig matrix      
      do i = 1,min(m,n)
        j = ipiv(i)
        if (i /= j) then
          call cswap(n,A(i,1),lda,A(j,1),lda)
        end if
      end do
      end subroutine

      subroutine zlugen(m,n,A,lda,L,ldl,R,ldr)
      integer m,n,lda,ldl,ldr
      double complex A(lda,n),L(ldl,min(m,n)),R(ldr,n)
      integer ipiv(min(m,n)),info,i,j
      external zswap,zlacpy,zgetrf
      if (m >= n) then
        call zlacpy('0',m,n,A,lda,L,ldl)
        call zgetrf(m,n,L,ldl,ipiv,info)
        call zlacpy('U',m,n,L,ldl,R,ldr)
      else
        call zlacpy('0',m,n,A,lda,R,ldr)
        call zgetrf(m,n,R,ldr,ipiv,info)
        call zlacpy('L',m,n,R,ldr,L,ldl)
      end if
      do i = 1,min(m,n)
        do j = 1,i-1
          L(j,i) = 0d0
        end do
        L(i,i) = 1d0
      end do
c permute the orig matrix      
      do i = 1,min(m,n)
        j = ipiv(i)
        if (i /= j) then
          call zswap(n,A(i,1),lda,A(j,1),lda)
        end if
      end do
      end subroutine

      subroutine sluchk(m,n,A,lda,L,ldl,R,ldr)
      integer m,n,lda,ldl,ldr
      real A(lda,n),L(ldl,min(m,n)),R(ldr,n)
      real rnrm,slange
      external sgemm,slange,spftol
      character*4 spftol
      real wrk(1)
      integer i,j

c zero lower triangle of R
      do j = 1,n-1
        do i = j+1,min(m,n)
          R(i,j) = 0e0
        end do
      end do
c form A - L*R
      call sgemm('N','N',m,n,min(m,n),1e0,L,ldl,R,ldr,-1e0,A,lda)
c get frobenius norm
      rnrm = slange('M',m,n,A,lda,wrk)
      write(*,1001) rnrm,spftol(rnrm)
      return

 1001 format(6x,'residual error = ',10x,E21.12,5x,A6)
      end subroutine

      subroutine dluchk(m,n,A,lda,L,ldl,R,ldr)
      integer m,n,lda,ldl,ldr
      double precision A(lda,n),L(ldl,min(m,n)),R(ldr,n)
      double precision rnrm,dlange
      external dgemm,dlange,dpftol
      character*4 dpftol
      double precision wrk(1)
      integer i,j

c zero lower triangle of R
      do j = 1,n-1
        do i = j+1,min(m,n)
          R(i,j) = 0e0
        end do
      end do
c form A - L*R
      call dgemm('N','N',m,n,min(m,n),1d0,L,ldl,R,ldr,-1d0,A,lda)
c get frobenius norm
      rnrm = dlange('M',m,n,A,lda,wrk)
      write(*,1001) rnrm,dpftol(rnrm)
      return

 1001 format(6x,'residual error = ',10x,E21.12,5x,A6)
      end subroutine

      subroutine cluchk(m,n,A,lda,L,ldl,R,ldr)
      integer m,n,lda,ldl,ldr
      complex A(lda,n),L(ldl,min(m,n)),R(ldr,n)
      real rnrm,clange
      external cgemm,clange,spftol
      character*4 spftol
      real wrk(1)
      integer i,j

c zero lower triangle of R
      do j = 1,n-1
        do i = j+1,min(m,n)
          R(i,j) = 0e0
        end do
      end do
c form A - L*R
      call cgemm('N','N',m,n,min(m,n),(1e0,0e0),L,ldl,R,ldr,(-1e0,0e0),
     +A,lda)
c get frobenius norm
      rnrm = clange('M',m,n,A,lda,wrk)
      write(*,1001) rnrm,spftol(rnrm)
      return

 1001 format(6x,'residual error = ',10x,E21.12,5x,A6)
      end subroutine

      subroutine zluchk(m,n,A,lda,L,ldl,R,ldr)
      integer m,n,lda,ldl,ldr
      double complex A(lda,n),L(ldl,min(m,n)),R(ldr,n)
      double precision rnrm,zlange
      external zgemm,zlange,dpftol
      character*4 dpftol
      double precision wrk(1)
      integer i,j

c zero lower triangle of R
      do j = 1,n-1
        do i = j+1,min(m,n)
          R(i,j) = 0e0
        end do
      end do
c form A - L*R
      call zgemm('N','N',m,n,min(m,n),(1d0,0d0),L,ldl,R,ldr,(-1d0,0d0),
     +A,lda)
c get frobenius norm
      rnrm = zlange('M',m,n,A,lda,wrk)
      write(*,1001) rnrm,dpftol(rnrm)
      return

 1001 format(6x,'residual error = ',10x,E21.12,5x,A6)
      end subroutine

      subroutine slupgen(m,n,A,lda,L,ldl,R,ldr,p)
      integer m,n,lda,ldl,ldr,p(m)
      real A(lda,n),L(ldl,min(m,n)),R(ldr,n)
      integer ipiv(min(m,n)),info,i,j,tmp
      external sswap,slacpy,sgetrf
      if (m >= n) then
        call slacpy('0',m,n,A,lda,L,ldl)
        call sgetrf(m,n,L,ldl,ipiv,info)
        call slacpy('U',m,n,L,ldl,R,ldr)
      else
        call slacpy('0',m,n,A,lda,R,ldr)
        call sgetrf(m,n,R,ldr,ipiv,info)
        call slacpy('L',m,n,R,ldr,L,ldl)
      end if
      do i = 1,min(m,n)
        do j = 1,i-1
          L(j,i) = 0e0
        end do
        L(i,i) = 1e0
      end do
c generate permutation
      do i = 1,m
        p(i) = i
      end do
      do i = 1,min(m,n)
        j = ipiv(i)
        if (i /= j) then
          tmp = p(i)
          p(i) = p(j)
          p(j) = tmp
        end if
      end do
c zero lower triangle of R
      do j = 1,n-1
        do i = j+1,min(m,n)
          R(i,j) = 0e0
        end do
      end do
      end subroutine

      subroutine dlupgen(m,n,A,lda,L,ldl,R,ldr,p)
      integer m,n,lda,ldl,ldr,p(m)
      double precision A(lda,n),L(ldl,min(m,n)),R(ldr,n)
      integer ipiv(min(m,n)),info,i,j,tmp
      external dswap,dlacpy,dgetrf
      if (m >= n) then
        call dlacpy('0',m,n,A,lda,L,ldl)
        call dgetrf(m,n,L,ldl,ipiv,info)
        call dlacpy('U',m,n,L,ldl,R,ldr)
      else
        call dlacpy('0',m,n,A,lda,R,ldr)
        call dgetrf(m,n,R,ldr,ipiv,info)
        call dlacpy('L',m,n,R,ldr,L,ldl)
      end if
      do i = 1,min(m,n)
        do j = 1,i-1
          L(j,i) = 0d0
        end do
        L(i,i) = 1d0
      end do
c generate permutation
      do i = 1,m
        p(i) = i
      end do
      do i = 1,min(m,n)
        j = ipiv(i)
        if (i /= j) then
          tmp = p(i)
          p(i) = p(j)
          p(j) = tmp
        end if
      end do
c zero lower triangle of R
      do j = 1,n-1
        do i = j+1,min(m,n)
          R(i,j) = 0d0
        end do
      end do
      end subroutine

      subroutine clupgen(m,n,A,lda,L,ldl,R,ldr,p)
      integer m,n,lda,ldl,ldr,p(m)
      complex A(lda,n),L(ldl,min(m,n)),R(ldr,n)
      integer ipiv(min(m,n)),info,i,j,tmp
      external cswap,clacpy,cgetrf
      if (m >= n) then
        call clacpy('0',m,n,A,lda,L,ldl)
        call cgetrf(m,n,L,ldl,ipiv,info)
        call clacpy('U',m,n,L,ldl,R,ldr)
      else
        call clacpy('0',m,n,A,lda,R,ldr)
        call cgetrf(m,n,R,ldr,ipiv,info)
        call clacpy('L',m,n,R,ldr,L,ldl)
      end if
      do i = 1,min(m,n)
        do j = 1,i-1
          L(j,i) = 0e0
        end do
        L(i,i) = 1e0
      end do
c generate permutation
      do i = 1,m
        p(i) = i
      end do
      do i = 1,min(m,n)
        j = ipiv(i)
        if (i /= j) then
          tmp = p(i)
          p(i) = p(j)
          p(j) = tmp
        end if
      end do
c zero lower triangle of R
      do j = 1,n-1
        do i = j+1,min(m,n)
          R(i,j) = 0e0
        end do
      end do
      end subroutine

      subroutine zlupgen(m,n,A,lda,L,ldl,R,ldr,p)
      integer m,n,lda,ldl,ldr,p(m)
      double complex A(lda,n),L(ldl,min(m,n)),R(ldr,n)
      integer ipiv(min(m,n)),info,i,j,tmp
      external zswap,zlacpy,zgetrf
      if (m >= n) then
        call zlacpy('0',m,n,A,lda,L,ldl)
        call zgetrf(m,n,L,ldl,ipiv,info)
        call zlacpy('U',m,n,L,ldl,R,ldr)
      else
        call zlacpy('0',m,n,A,lda,R,ldr)
        call zgetrf(m,n,R,ldr,ipiv,info)
        call zlacpy('L',m,n,R,ldr,L,ldl)
      end if
      do i = 1,min(m,n)
        do j = 1,i-1
          L(j,i) = 0d0
        end do
        L(i,i) = 1d0
      end do
c generate permutation
      do i = 1,m
        p(i) = i
      end do
      do i = 1,min(m,n)
        j = ipiv(i)
        if (i /= j) then
          tmp = p(i)
          p(i) = p(j)
          p(j) = tmp
        end if
      end do
c zero lower triangle of R
      do j = 1,n-1
        do i = j+1,min(m,n)
          R(i,j) = 0d0
        end do
      end do
      end subroutine

c converts a linear permutation into LAPACK ipiv-style form
      subroutine p2ipiv(n,p)
      integer n,p(n)
      integer q(n),i,j,k
      do i = 1,n
        q(p(i)) = i
      end do
      do i = 1,n
        j = p(i)
        k = q(i)
        if (j /= i) then
          p(k) = j
          q(j) = k
          p(i) = j
        end if
      end do
      end subroutine

      subroutine slupchk(m,n,A,lda,L,ldl,R,ldr,p)
      integer m,n,lda,ldl,ldr,p(m)
      real A(lda,n),L(ldl,min(m,n)),R(ldr,n)
      real rnrm,slange
      external sswap,sgemm,slange,spftol
      character*4 spftol
      real wrk(1)
      integer i,j

c convert p into successive swaps      
      call p2ipiv(m,p)
c form A - L*R
      do i = 1,m
        j = p(i)
        if (i /= j) then
          call sswap(n,A(i,1),lda,A(j,1),lda)
        end if
      end do
      call sgemm('N','N',m,n,min(m,n),1e0,L,ldl,R,ldr,-1e0,A,lda)
c get frobenius norm
      rnrm = slange('M',m,n,A,lda,wrk)
      write(*,1001) rnrm,spftol(rnrm)
      return

 1001 format(6x,'residual error = ',10x,E21.12,5x,A6)
      end subroutine

      subroutine dlupchk(m,n,A,lda,L,ldl,R,ldr,p)
      integer m,n,lda,ldl,ldr,p(m)
      double precision A(lda,n),L(ldl,min(m,n)),R(ldr,n)
      double precision rnrm,dlange
      external dswap,dgemm,dlange,dpftol
      character*4 dpftol
      double precision wrk(1)
      integer i,j

c convert p into successive swaps      
      call p2ipiv(m,p)
c form A - L*R
      do i = 1,m
        j = p(i)
        if (i /= j) then
          call dswap(n,A(i,1),lda,A(j,1),lda)
        end if
      end do
      call dgemm('N','N',m,n,min(m,n),1d0,L,ldl,R,ldr,-1d0,A,lda)
c get frobenius norm
      rnrm = dlange('M',m,n,A,lda,wrk)
      write(*,1001) rnrm,dpftol(rnrm)
      return

 1001 format(6x,'residual error = ',10x,E21.12,5x,A6)
      end subroutine

      subroutine clupchk(m,n,A,lda,L,ldl,R,ldr,p)
      integer m,n,lda,ldl,ldr,p(m)
      complex A(lda,n),L(ldl,min(m,n)),R(ldr,n)
      real rnrm,clange
      external cswap,cgemm,clange,spftol
      character*4 spftol
      real wrk(1)
      integer i,j

c convert p into successive swaps      
      call p2ipiv(m,p)
c form A - L*R
      do i = 1,m
        j = p(i)
        if (i /= j) then
          call cswap(n,A(i,1),lda,A(j,1),lda)
        end if
      end do
      call cgemm('N','N',m,n,min(m,n),(1e0,0e0),L,ldl,R,ldr,(-1e0,0e0),
     +A,lda)
c get frobenius norm
      rnrm = clange('M',m,n,A,lda,wrk)
      write(*,1001) rnrm,spftol(rnrm)
      return

 1001 format(6x,'residual error = ',10x,E21.12,5x,A6)
      end subroutine

      subroutine zlupchk(m,n,A,lda,L,ldl,R,ldr,p)
      integer m,n,lda,ldl,ldr,p(m)
      double complex A(lda,n),L(ldl,min(m,n)),R(ldr,n)
      double precision rnrm,zlange
      external zswap,zgemm,zlange,dpftol
      character*4 dpftol
      double precision wrk(1)
      integer i,j

c convert p into successive swaps      
      call p2ipiv(m,p)
c form A - L*R
      do i = 1,m
        j = p(i)
        if (i /= j) then
          call zswap(n,A(i,1),lda,A(j,1),lda)
        end if
      end do
      call zgemm('N','N',m,n,min(m,n),(1d0,0d0),L,ldl,R,ldr,(-1d0,0d0),
     +A,lda)
c get frobenius norm
      rnrm = zlange('M',m,n,A,lda,wrk)
      write(*,1001) rnrm,dpftol(rnrm)
      return

 1001 format(6x,'residual error = ',10x,E21.12,5x,A6)
      end subroutine
