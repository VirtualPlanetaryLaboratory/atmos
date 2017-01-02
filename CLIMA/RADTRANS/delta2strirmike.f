      SUBROUTINE DELTA2STRIR(SRFALB,ASY,TAULAM,OMG0,
     &  FUP,FDN,BPLANCK)
     
C
      INCLUDE 'CLIMA/INCLUDE/header.inc'
C
C   This is Jim Kasting's version of the Toon et al. 2-stream code.
C   (Ref.: JGR 94, 16287, 1989).  It vectorizes over height, rather
C   than wavelength.
C
C   The direct solar flux at the top (pi*Fs) is treated as 1.
C
C   For now, at least, it is hardwired as the delta-quadrature
C   approximation. The Eddington coefficients are in the code but are
C   commented out.
C
      PARAMETER(NZ=ND-1, NZ2=2*NZ)
      REAL NORM,B1
      DIMENSION TAU(NZ), G(NZ), GAM1(NZ), GAM2(NZ),
     1  ALAM(NZ), CGAM(NZ), E1(NZ), E2(NZ), E3(NZ), E4(NZ),
     2  CP0(NZ), CPB(NZ), CM0(NZ), CMB(NZ), Y1(NZ), Y2(NZ),
     3  W0(NZ)
      DIMENSION A(NZ2), B(NZ2), D(NZ2), E(NZ2), Y(NZ2)
      DIMENSION ASY(NZ),TAULAM(NZ), OMG0(NZ),
     &     FUP(ND),FDN(ND),BPLANCK(ND)
C
      DATA C,HP,BK,SIGMA,PI,SM/3.E10, 6.63E-27, 1.38E-16, 5.67E-5,
     2  3.14159, 1.67E-24/
     
      
      ALB = SRFALB
      PI = 3.14159
      emissivity = 1.
      BCON = 2.*HP/C/C
      HK = HP/BK
C      SQ3 = SQRT(3.)
C      U1 = 1./SQ3
C      U0M = 1./U0
C      U0M2 = U0M*U0M
      U1 = 0.5
      U1M = 1./U1
      NZM1 = NZ - 1
      NZP1 = NZ + 1
      MZ2 = NZ2
C
C   Delta Approximation
      DO 1 N=1,NZ
      TAU(N) = TAULAM(N)
      W0(N) = OMG0(N)
 1    G(N) = ASY(N)
C
C   Calculate the gamma's, lambda's, and e's
      DO 2 N=1,NZ
C      GAM1(N) = SQ3*(2. - W0(N)*(1.+G(N)))/2.
C      GAM2(N) = SQ3*W0(N)*(1.-G(N))/2.
      GAM1(N) = 2 - W0(N)*(1.+G(N)) 
      GAM2(N) = W0(N)*(1.-G(N))
C
      ALAM(N) = SQRT(GAM1(N)*GAM1(N) - GAM2(N)*GAM2(N))
      CGAM(N) = (GAM1(N) - ALAM(N))/GAM2(N)
      EMLT = EXP(-ALAM(N)*TAU(N))
C      PRINT*, ALAM(N),TAU(N),N
C
      E1(N) = 1. + CGAM(N)*EMLT
      E2(N) = 1. - CGAM(N)*EMLT
      E3(N) = CGAM(N) + EMLT
   2  E4(N) = CGAM(N) - EMLT
C
C   Calculate A, B, and D, i.e. the coefficients of the tridiagonal
C      matrix
C   Top of atmosphere
      A(1) = 0.
      B(1) = E1(1)
      D(1) = -E2(1)
C   Odd coefficients
      DO 3 N=1,NZM1
      L = 2*N + 1
      A(L) = E2(N)*E3(N) - E4(N)*E1(N)
      B(L) = E1(N)*E1(N+1) - E3(N)*E3(N+1)
   3  D(L) = E3(N)*E4(N+1) - E1(N)*E2(N+1)
C
C   Even coefficients
      DO 4 N=1,NZM1
      L = 2*N
      A(L) = E2(N+1)*E1(N) - E3(N)*E4(N+1)
      B(L) = E2(N)*E2(N+1) - E4(N)*E4(N+1)
   4  D(L) = E1(N+1)*E4(N+1) - E2(N+1)*E3(N+1)
C
C   Bottom of atmosphere
      A(NZ2) = E1(NZ) - ALB*E3(NZ)
      B(NZ2) = E2(NZ) - ALB*E4(NZ)
      D(NZ2) = 0.
      B(1) = E1(1)
C
C   Now, set up the RHS of the equation:
C
      NORM = 2*U1*PI
      DO 6 N=1,NZ
C
      B0n = BPLANCK(N)
      B1 = BPLANCK(N+1)
      B1n = (B1-B0n)/TAU(N)
      CP0(N) = NORM*(B0n + B1n*(1/(GAM1(N)+GAM2(N))))
      CPB(N) = NORM*(B0n + B1n*(TAU(N)+1/(GAM1(N)+GAM2(N))))
      CM0(N) = NORM*(B0n + B1n*(-1/(GAM1(N)+GAM2(N))))
   6  CMB(N) = NORM*(B0n + B1n*(TAU(N)-1/(GAM1(N)+GAM2(N))))
      SSFC = emissivity*PI*BPLANCK(ND)
C
C   Odd coefficients
      E(1) = - CM0(1)
      DO 7 N=1,NZM1
      L = 2*N + 1
   7  E(L) = (CP0(N+1)-CPB(N))*E3(N) + (CMB(N)-CM0(N+1))*E1(N)
C     
C   Even coefficients
      DO 8 N=1,NZM1
      L = 2*N
   8  E(L) = (CP0(N+1)-CPB(N))*E2(N+1) - (CM0(N+1)-CMB(N))*E4(N+1)
      E(NZ2) = SSFC - CPB(NZ) + ALB*CMB(NZ)
C     Call the tridiagonal solver.  Use Numerical Recipes for now.
       CALL TRIDAG(A,B,D,E,Y,MZ2)
C      WRITE(1,210)
C      WRITE(1,111) Y
C
      DO 9 N=1,NZ
      L = 2*N
      L1 = L-1
      Y1(N) = Y(L1)
   9  Y2(N) = Y(L)
C  Calculate upward and downward fluxes.
C
       FUP(1) = ((Y1(1)*E3(1) - Y2(1)*E4(1)) + CP0(1))
       FDN(1) = 0
      DO 10 N=1,NZ
         FUP(N+1) = (Y1(N)*E1(N) + Y2(N)*E2(N)
     &     + CPB(N))
         FDN(N+1) = (Y1(N)*E3(N) + Y2(N)*E4(N)
     &     + CMB(N))
 10   CONTINUE
      RETURN
      END
