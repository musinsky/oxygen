      Function Random(n)         ! produces  0 <= random number < 1
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
C*
C***   The whole family of  Random number Generators
C*     x = Random(n)    where n=1,2,... - what kind of Generator
C*     For a while we have :
C*     n=1  : uniformly distributed numbers in [0,1) by G.Marsaglia
C*     n=2  :    ---"---       ---"---      in (0,1) by L'Ecuyer
C*
C========================================================================
C      n = 1
C This random number generator originally appeared in "Toward a Universal
C Random Number Generator" by George Marsaglia and Arif Zaman.
C Florida State University Report: FSU-SCRI-87-50 (1987)
C
C It passes ALL of the tests for random number generators and has a period
C of 2^144, is completely portable (gives bit identical results on all
C machines with at least 48-bit mantissas in the floating point representation).
C
C The algorithm is a combination of a Fibonacci sequence (with lags of 97
C and 33, and operation "subtraction plus one, modulo one") and an
C "arithmetic sequence" (using subtraction).
C For each random number only 5 float-point and 3 fixed-point additions
C are performed.
C
C Subroutine RandomInitiate(i,k)  initiates this Generator with 2 integer
C numbers  0 <= i <= 31328  and  0 <= k <= 30081.  It permits to get about
C 900000000 different random series.
C Default initialisation is :  i=1802,  k=9373
C RandomInitiate defines only serie, its members always start from the
C beginning.  Full state of this Generator contains in
C      common /RanGen_Parameters/ istat(4),rstat(100)
C      real*8 rstat
C Random(-1) prints istat(1) and istat(2) : current serie idents i and k.
C========================================================================
C                      n = 2
C         Portable random number generator proposed by L'Ecuyer
C              in Comm. ACM 31:743 (1988)
C Subroutine RandomInitiate(i,k)  initializes this Generator with 2 integer
C numbers i>0  and  k>0.   Default initialisation is :  i=1802,  k=9373.
C These 2 numbers completely define the state of this Generator.
C Random(-1) prints them.
C========================================================================
C
      parameter(ngener=3)         ! increase it for new member
      common /RanGen_Parameters/ ij,kl, i97,j97, c,cd,cm, u(97)
      data ij/-1/, kl/-1/, i97/-1/
C
      if(i97.le.0) call RandomInitiate(1802,9373)
      if(n.gt.0 .and. n.le.ngener) goto(1,2,3) n   ! add here for new member
      if(n.eq.-1) Write(*,*) ' Random(-1): i,k=',ij,kl
  1   continue
      r = U(I97) - U(J97)
      if( r .lt. 0.0d0 ) r = r + 1.0d0
      U(I97) = r
      I97 = I97 - 1
      if(I97 .eq. 0) I97 = 97
      J97 = J97 - 1
      if(J97 .eq. 0) J97 = 97
      C = C - CD
      if( C .lt. 0.0d0 ) C = C + CM
      r = r - C
      if( r .lt. 0.0d0 ) r = r + 1.0d0           !   [ 0,1 )
      Random = r
      return
C
  2   continue
      K = ij/53668
      ij = 40014*(ij - K*53668) - K*12211
      IF (ij .LT. 0) ij=ij+2147483563
      K = kl/52774
      kl = 40692*(kl - K*52774) - K* 3791
      IF (kl .LT. 0) kl=kl+2147483399
      IZ = ij - kl
      IF (IZ .LT. 1) IZ = IZ + 2147483562
      Random = IZ * 4.656613d-10                    ! ( 0,1 )
      return

  3   ij=314159261*ij+90763307
      if(ij.lt.0) ij=ij+2147483647
      r=ij/2.147483648d9
      if(r.gt.1) write(*,*) r,'  > 1 !!!'
      random=r
      end

      subroutine RandomInitiate(i0,k0)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      common /RanGen_Parameters/ ij,kl, i97,j97, c,cd,cm, u(97)
c---
      ij = i0
      kl = k0
        i = mod(IJ/177, 177) + 2
        j = mod(IJ    , 177) + 2
        k = mod(KL/169, 178) + 1
        l = mod(KL,     169)
        do ii = 1, 97
          s = 0.0d0
          t = 0.5d0
          do jj = 1, 48
            m = mod(mod(i*j, 179)*k, 179)
            i = j
            j = k
            k = m
            l = mod(53*l+1, 169)
            if (mod(l*m, 64) .ge. 32)  s = s + t
            t = 0.5d0 * t
          enddo
          U(ii) = s
        enddo
        C = 362436.0d0 / 16777216.0d0
        CD = 7654321.0d0 / 16777216.0d0
        CM = 16777213.0d0 /16777216.0d0
        I97 = 97
        J97 = 33
      return
      end
      
       FUNCTION RNDM(K)
       REAL *8 RANDOM
       RNDM=RANDOM(1)
       return
       END