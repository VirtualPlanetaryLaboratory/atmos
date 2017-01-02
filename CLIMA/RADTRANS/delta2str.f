      SUBROUTINE DELTA2STR(SRFALB,AMU0,ASY,TAULAM,OMG0,
     &  FUP,FDN)
C
      INCLUDE 'CLIMA/INCLUDE/header.inc'
C   (Jim Kasting's version)
C   This is Jim Kasting's version of the Toon et al. 2-stream code.
C   (R
C   The direct solar flux at the top (pi*Fs) is treated as 1.
C
C   For now, at least, it is hardwired as the delta-quadrature
C   approximation. The Eddington coefficients are in the code but are
C   commented out.
C
      PARAMETER(NZ=ND-1, NZ2=2*NZ)
      DIMENSION TAU(NZ), TAUC(NZ), G(NZ), GAM1(NZ), GAM2(NZ), GAM3(NZ),
     1  GAM4(NZ), ALAM(NZ), CGAM(NZ), E1(NZ), E2(NZ), E3(NZ), E4(NZ),
     2  CP0(NZ), CPB(NZ), CM0(NZ), CMB(NZ), Y1(NZ), Y2(NZ), ! EWS - AMEAN(ND) variable not used
     3  W0(NZ),DIRECT(ND)
      DIMENSION A(NZ2), B(NZ2), D(NZ2), E(NZ2), Y(NZ2)
      DIMENSION ASY(NZ),TAULAM(NZ),FMT(NZ),OMG0(NZ),
     &     FUP(ND),FDN(ND)
C
      ALB = SRFALB
      U0 = AMU0
C      U1 = 0.5  (Eddington value)
      SQ3 = SQRT(3.)
      U1 = 1./SQ3
      U0M = 1./U0
      U0M2 = U0M*U0M
      U1M = 1./U1
      NZM1 = NZ - 1
      NZP1 = NZ + 1
      MZ2 = NZ2
C
C   Delta Approximation
C-AP  I used approximation form Joseph et al. 1976
!GNA - is that the delta approx that Dave Crisp mentioned?
      DO 1 N=1,NZ
      FMT(N) = ASY(N)*ASY(N)
      TAU(N) = TAULAM(N)*(1-OMG0(N)*FMT(N))
      W0(N) = OMG0(N)*(1-FMT(N))/(1-OMG0(N)*FMT(N))
 1    G(N) = ASY(N)/(1+ASY(N))
C
C   Calculate the gamma's, lambda's, and e's
      DO 2 N=1,NZ
C      GAM1(N) = (7. - W0(N)*(4.+3.*G(N)))/4.
C      GAM2(N) = - (1. - W0(N)*(4.-3.*G(N)))/4.
C      GAM3(N) = (2. - 3.*G(N)*U0)/4.
C   (Eddington values above -- not currently used)
C
      GAM1(N) = SQ3*(2. - W0(N)*(1.+G(N)))/2.
      GAM2(N) = SQ3*W0(N)*(1.-G(N))/2.
      GAM3(N) = (1. - SQ3*G(N)*U0)/2.
      GAM4(N) = 1. - GAM3(N)
C
      ALAM(N) = SQRT(GAM1(N)*GAM1(N) - GAM2(N)*GAM2(N))
      CGAM(N) = (GAM1(N) - ALAM(N))/GAM2(N)
      EMLT = EXP(-ALAM(N)*TAU(N))
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
C
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
C
C   Now, set up the RHS of the equation:
C   TAUC(N) is the optical depth above layer N
      TAUC(1) = 0.
      DO 5 N=2,NZ
   5  TAUC(N) = TAUC(N-1) + TAU(N-1)
C
      DIRECT(1) = U0
      DO 6 N=1,NZ
      ET0 = EXP(-TAUC(N)/U0)
      ETB = ET0 * EXP(-TAU(N)/U0)
      DIRECT(N+1) = ETB*U0
      DENOM = ALAM(N)*ALAM(N) - U0M2
      FACP = W0(N) * ((GAM1(N)-U0M)*GAM3(N) + GAM4(N)*GAM2(N))
      FACM = W0(N) * ((GAM1(N)+U0M)*GAM4(N) + GAM2(N)*GAM3(N))
C
      CP0(N) = ET0*FACP/DENOM
      CPB(N) = ETB*FACP/DENOM
      CM0(N) = ET0*FACM/DENOM
   6  CMB(N) = ETB*FACM/DENOM
      SSFC = ALB*DIRECT(ND)
C   DIRECT(N) is the direct solar flux at the top of layer N
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
C
C      WRITE(1,104)
C 104  FORMAT(/5X,'CP0',7X,'CPB',7X,'CM0',7X,'CMB',7X,'TAU',7X,'TAUC')
C      WRITE(1,105) (CP0(N),CPB(N), CM0(N),CMB(N),TAU(N),TAUC(N),N=1,NZ)
C 105  FORMAT(1X,1P6E10.3)
C      WRITE(1,110)
C 110  FORMAT(/1X,'E')
C      WRITE(1,111) E
c  111  FORMAT(1X,1P5E10.3) !EWS - not used
C
C   Call the tridiagonal solver.  Use Numerical Recipes for now.
      CALL TRIDAG(A,B,D,E,Y,MZ2)
C      WRITE(1,210)
c 210  FORMAT(/1X,'Y (after TRIDAG)') !EWS - not used
C      WRITE(1,111) Y
C
      DO 9 N=1,NZ
      L = 2*N
      L1 = L-1
      Y1(N) = Y(L1)
   9  Y2(N) = Y(L)
C
C   Calculate the mean intensities, AMEAN(N), at the boundaries betweenC   
C   the layers.  AMEAN(N) is the intensity at the top of layer N.
C   (Top of atmosphere is different}
C      AMEAN(1) = U1M * (Y1(1)*E3(1) - Y2(1)*E4(1) + CP0(1)) + 1.
C      DO 10 N=1,NZ
C  10  AMEAN(N+1) = U1M * (Y1(N)*(E1(N)+E3(N)) + Y2(N)*(E2(N)+E4(N))
C     1  + CPB(N) + CMB(N)) + DIRECT(N+1)
C
C  Calculate upward and downward fluxes.
C
      FUP(1) = ((Y1(1)*E3(1) - Y2(1)*E4(1)) + CP0(1))
      FDN(1) = DIRECT(1)
      DO 10 N=1,NZ
         FUP(N+1) = (Y1(N)*E1(N) + Y2(N)*E2(N)
     &     + CPB(N))
         FDN(N+1) = (Y1(N)*E3(N) + Y2(N)*E4(N)
     &     + CMB(N)) + DIRECT(N+1)
 10   CONTINUE
C   Print out the results
C      WRITE(1,108)
C 108  FORMAT(/5X,'Y1',8X,'Y2',8X,'DIRECT',4X,'AMEAN')
C      WRITE(1,109) (Y1(N),Y2(N),DIRECT(N),AMEAN(N),N=1,NZ)
C 109  FORMAT(1X,1P4E10.3)
C      WRITE(1,112) DIRECT(ND),AMEAN(ND)
C 112  FORMAT(21X,1P2E10.3)
C
      RETURN
      END
