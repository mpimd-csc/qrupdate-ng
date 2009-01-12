      program tchshx
      integer n,i,j

      write (*,*)
      write (*,*) 'testing QR column shift routines.'
      write (*,*) 'All residual errors are expected to be small.'
      write (*,*)

      n = 50
      i = 20
      j = 40
      write (*,*) 'schshx test (left shift):'
      call stest(n,i,j)
      write (*,*) 'dchshx test (left shift):'
      call dtest(n,i,j)
      write (*,*) 'cchshx test (left shift):'
      call ctest(n,i,j)
      write (*,*) 'zchshx test (left shift):'
      call ztest(n,i,j)

      i = 40
      j = 20
      write (*,*) 'schshx test (right shift):'
      call stest(n,i,j)
      write (*,*) 'dchshx test (right shift):'
      call dtest(n,i,j)
      write (*,*) 'cchshx test (right shift):'
      call ctest(n,i,j)
      write (*,*) 'zchshx test (right shift):'
      call ztest(n,i,j)

      call pstats
      end program

      subroutine stest(n,i,j)
      integer n,i,j
      real A(n,n),R(n,n),wrk(2*n)
      external srandg,schgen,sswap,schshx,schchk
      integer k
c set up random matrix
      call srandg(n,n,A,n)
c generate A'*A and its Cholesky decomposition
      call schgen(n,A,n,R,n)
c update matrix
      if (i < j) then
        do k = i,j-1
          call sswap(n,A(1,k),1,A(1,k+1),1)
          call sswap(n,A(k,1),n,A(k+1,1),n)
        end do
      else if (i > j) then
        do k = i,j+1,-1
          call sswap(n,A(1,k),1,A(1,k-1),1)
          call sswap(n,A(k,1),n,A(k-1,1),n)
        end do
      end if
c update factorization
      call schshx(n,R,n,i,j,wrk)
c check result
      call schchk(n,A,n,R,n)

      end subroutine

      subroutine dtest(n,i,j)
      integer n,i,j
      double precision A(n,n),R(n,n),wrk(2*n)
      external drandg,dchgen,dswap,dchshx,dchchk
      integer k
c set up random matrix
      call drandg(n,n,A,n)
c generate A'*A and its Cholesky decomposition
      call dchgen(n,A,n,R,n)
c update matrix
      if (i < j) then
        do k = i,j-1
          call dswap(n,A(1,k),1,A(1,k+1),1)
          call dswap(n,A(k,1),n,A(k+1,1),n)
        end do
      else if (i > j) then
        do k = i,j+1,-1
          call dswap(n,A(1,k),1,A(1,k-1),1)
          call dswap(n,A(k,1),n,A(k-1,1),n)
        end do
      end if
c update factorization
      call dchshx(n,R,n,i,j,wrk)
c check result
      call dchchk(n,A,n,R,n)

      end subroutine

      subroutine ctest(n,i,j)
      integer n,i,j
      complex A(n,n),R(n,n),wrk(n)
      real rwrk(n)
      external crandg,cchgen,cswap,cchshx,cchchk
      integer k
c set up random matrix
      call crandg(n,n,A,n)
c generate A'*A and its Cholesky decomposition
      call cchgen(n,A,n,R,n)
c update matrix
      if (i < j) then
        do k = i,j-1
          call cswap(n,A(1,k),1,A(1,k+1),1)
          call cswap(n,A(k,1),n,A(k+1,1),n)
        end do
      else if (i > j) then
        do k = i,j+1,-1
          call cswap(n,A(1,k),1,A(1,k-1),1)
          call cswap(n,A(k,1),n,A(k-1,1),n)
        end do
      end if
c update factorization
      call cchshx(n,R,n,i,j,wrk,rwrk)
c check result
      call cchchk(n,A,n,R,n)

      end subroutine

      subroutine ztest(n,i,j)
      integer n,i,j
      double complex A(n,n),R(n,n),wrk(n)
      double precision rwrk(n)
      external zrandg,zchgen,zswap,zchshx,zchchk
      integer k
c set up random matrix
      call zrandg(n,n,A,n)
c generate A'*A and its Cholesky decomposition
      call zchgen(n,A,n,R,n)
c update matrix
      if (i < j) then
        do k = i,j-1
          call zswap(n,A(1,k),1,A(1,k+1),1)
          call zswap(n,A(k,1),n,A(k+1,1),n)
        end do
      else if (i > j) then
        do k = i,j+1,-1
          call zswap(n,A(1,k),1,A(1,k-1),1)
          call zswap(n,A(k,1),n,A(k-1,1),n)
        end do
      end if
c update factorization
      call zchshx(n,R,n,i,j,wrk,rwrk)
c check result
      call zchchk(n,A,n,R,n)

      end subroutine
