      SUBROUTINE MSCAT(SIGR,U0,WAV,sq,wl,columndepth,TTOT,S,USETD)

! WAV here is a single value at the center of a wavlength interval, wl is the index

C
C          THIS SUBROUTINE COMPUTES NORMALIZED SOURCE FUNCTIONS DUE
C     TO RAYLEIGH SCATTERING USING YUK YUNGS METHOD.
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      real*8 mass
      REAL*8 sq(kj,nz,kw),TTOT(NZ),S(NZ)
      integer wl
      real*8 columndepth(kj,NZ)
      DIMENSION B(ML2,ML2),IPVT(ML2),F(ML2),TAU(ML1),W0(ML1),
     2  D3(ML1),D5(ML1),D7(ML1),TA(NZ),TAB(NZ),W(NZ)
     2  ,SO3(NZ),SO2(NZ),W0P(NP),W0PS(NP),QEXT(NP),SIGS(NZ,NP),
     3  SIGA(NZ,NP),SIGAT(NZ),SIGST(NZ),TMS(NZ),TAB2(NZ)
      DIMENSION XI(ML1,2),EX1(ML1,2),E2(ML1,2),E3(ML1,2),E4(ML1,2),
     2  E5(ML1,2),E1S(ML1),XIS(ML1)
      REAL*8 L2(ML),L4(ML),G1(ML),G2(ML),G3(ML),G4(ML),M1(ML,ML),
     2  M3(ML,ML),M5(ML,ML)
      character*8 PLANET
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PHOTABLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/CBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/JBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/NBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/AERBLK.inc'
C
C   THIS SUBROUTINE DOES RAYLEIGH SCATTERING USING YUK YUNG'S
C   METHOD.  MIE SCATTERING BY PARTICLES IS COMPUTED BY USING VAN
C   DE HULST'S SIMILARITY RELATIONS TO SCALE FROM ANISOTROPIC TO
C   ISOTROPIC SCATTERING.  EQUATIONS ARE SOLVED ON A GRID THAT IS
C   EVENLY SPACED IN EXTINCTION OPTICAL DEPTH (TAU).
C
      PI = 3.14159
C
C   COMPUTE FACTORS FOR SCALING FORWARD TO ISOTROPIC SCATTERING
      G0 = 0.8

      do i=1,np

         LL=NQ-NP+i  !particles are the last elements in main loop
         if(USETD.EQ.1) LL=i+NQ !or are in tri-diag

         if (LL.EQ.LSO4AER .OR. LL.EQ.LSXO4AER) W0P(i) = 1.

         if (LL.EQ.LS8AER .OR. LL.EQ.LSXS7AER) then
           if (wav.gt.3500.) then
              W0P(i)=1.0
           else
            W0P(i) = 0.5
           endif 
         endif
      enddo


      DO 25 J=1,NP   
      W0PS(J) = (1.-G0)*W0P(J)/(1.-G0*W0P(J))
  25  QEXT(J) = 2.*(1.-G0*W0P(J))

C
C   COMPUTE EXTINCTION, ABSORPTION, AND SCATTERING CROSS SECTIONS FOR
C   PARTICLES (NOTE THAT THEY ARE SCALED TO COMPARE WITH MOLECULAR
C   ABSORPTION AND SCATTERING)
      DO 26 J=1,NP         
      DO 26 I=1,NZ
      SIGE = QEXT(J)*PI*RPAR(I,J)*RPAR(I,J)*AERSOL(I,J)/DEN(I)
      SIGS(I,J) = W0PS(J)*SIGE
  26  SIGA(I,J) = (1.-W0PS(J))*SIGE


C
C   COMPUTE SINGLE SCATTERING ALBEDO W(I) FOR GASES PLUS PARTICLES (cm^2)
      do I=1,NZ
       SIGAT(I)=0.0
       SIGST(I)=0.0
       do j=1,np
          SIGAT(I)=SIGAT(I)+SIGA(I,j)
          SIGST(I)=SIGST(I)+SIGS(I,j)
       enddo   

        sigabs=0.0
         do j=1,kj
          sigabs=sigabs + sq(j,i,wl)*absorbers(j,i)
         enddo
        sigabs=sigabs+sigat(I)
       
      W(I) = (SIGR + SIGST(I))/(SIGR + SIGST(I) + SIGABS)
      enddo

C


C   COMPUTE ABSORPTION (TAB) AND EXTINCTION (TA) OPTICAL DEPTHS
      do i=1,nz
         TAB(I)=0.0
      enddo

        do j=1,kj
           TAB(nz)=TAB(NZ) + sq(j,nz,wl)*columndepth(j,nz)
        enddo
       TAB(nz)=TAB(nz)+SIGAT(NZ)*DEN(NZ)*HSCALE(NZ)


      TMS(NZ) = SIGST(NZ)*DEN(NZ)*HSCALE(NZ)
      TA(NZ) = TAB(NZ) + SIGR*TTOT(NZ) + TMS(NZ)
!total extinction optical depth= absorption + rayleigh scattering + particle scattering

C
      do M=1,NZ1
       I = NZ - M   !run count from top down...
        do j=1,kj
           sfp=0.5*(sq(j,i,wl)+sq(j,i+1,wl)) * 
     $                    (columndepth(j,i) - columndepth(j,i+1)) 
          TAB(I)=TAB(I)+ sfp
        enddo
        TAB(I)=TAB(I)+TAB(I+1)+SIGAT(I)*DEN(I)

       TMS(I) = TMS(I+1) + SIGST(I)*DEN(I)*DZ(I)   !ACK check when grid goes variable
       TA(I) = TAB(I) + SIGR*TTOT(I) + TMS(I)
      enddo


C
C ***** DEFINE BOTTOM OF SCATTERING ATMOSPHERE WHERE THE PURE
C       ABSORPTION OPTICAL DEPTH IS EQUAL TO 5 *****
      DO 21 L=1,NZ
      I = NZ - L + 1
      IF(TAB(I).GT.5.) GO TO 22
  21  CONTINUE
      MZ = 1
      GO TO 23
  22  MZ = I + 1
  23  LZ = NZ - MZ + 1

C
      SDEPTH = TA(MZ)
      NL = SDEPTH/0.5
      NL = min(NL,ML)
      NL = max(NL,3)
      NL1 = NL + 1
      NLM1 = NL - 1
      NL2 = 2*NL
C
      TAU(1) = 0.
      TAU(NL1) = TA(MZ)
      DT = TAU(NL1)/FLOAT(NL)
      DO 2 J=2,NL
   2  TAU(J) = (J-1)*DT
C
C ***** FIND SINGLE SCATTERING ALBEDOS AT THE CENTERS OF THE EQUALLY-
C       SPACED TAU LEVELS *****
      W0(1) = W(NZ)
      W0(NL1) = W(MZ)
      IS = 1
      DO 3 J=2,NL
      DO 4 I=IS,LZ
      IN = NZ - I + 1
      IF(TA(IN).GT.TAU(J)) GO TO 5
   4  CONTINUE
   5  IS = NZ - IN - 1
C
C     TAU(J) LIES BETWEEN TA(IN) AND TA(IN+1)
C
      FR = (TAU(J) - TA(IN+1))/(TA(IN) - TA(IN+1))
      W0(J) = FR*W(IN) + (1.-FR)*W(IN+1)
   3  CONTINUE
C
C     CONVERT TO CENTERED VALUES
      DO 6 J=1,NL
   6  W0(J) = SQRT(W0(J)*W0(J+1))
C
C ***** COMPUTE THE VALUES OF LI AND MIJ *****
      DO 7 I=1,NL1
      DTAU = (I-1)*DT
      XIS(I) = DTAU
      E1S(I) = E1(DTAU)
   7  CONTINUE
C
      DO 30 I=1,NL1
      X = XIS(I)
      X2 = X*X
      X3 = X*X2
      X4 = X*X3
      X5 = X*X4
      X6 = X*X5
      E1X = E1S(I)
      PEX = EXP(-X)
      D3(I) = (X2*E1X + PEX*(1. - X))/2.
      D5(I) = (X4*E1X + PEX*(6. - 2.*X + X2 - X3))/24.
  30  D7(I) = (X6*E1X + PEX*(120. - 24.*X + 6.*X2 - 2.*X3 + X4
     2  - X5))/720.
C
      SM1 = 2.*(DT - 0.5 + D3(2))
      SM3 = 2.*(DT/3. - 0.25 + D5(2))
      SM5 = 2.*(DT/5. - 1./6. + D7(2))
      DO 8 I=1,NL
      L2(I) = D3(NL-I+1) - D3(NL-I+2)
      L4(I) = D5(NL-I+1) - D5(NL-I+2)
      M1(I,I) = SM1
      M3(I,I) = SM3
      M5(I,I) = SM5
   8  CONTINUE
C
      DO 16 I=2,NL
      M1(I,1) = D3(I-1) + D3(I+1) - 2.*D3(I)
      M3(I,1) = D5(I-1) + D5(I+1) - 2.*D5(I)
      M5(I,1) = D7(I-1) + D7(I+1) - 2.*D7(I)
  16  CONTINUE
C
      DO 17 I=2,NLM1
      I1 = I + 1
      DO 17 J=I1,NL
      K = J - I + 1
      M1(J,K) = M1(I,1)
      M3(J,K) = M3(I,1)
      M5(J,K) = M5(I,1)
  17  CONTINUE
C
      DO 18 I=2,NL
      I1 = I - 1
      DO 18 J=1,I1
      M1(J,I) = M1(I,J)
      M3(J,I) = M3(I,J)
      M5(J,I) = M5(I,J)
  18  CONTINUE
C
C ***** COMPUTE THE ELEMENTS OF THE MATRIX EQUATION *****
      U2 = U0*U0
      TG = TAU(NL1)
      C1 = 3.*U0/(32.*DT)
      C2 = 2.*ALB*EXP(-TG/U0)
C
      DO 9 I=1,NL
      D1 = EXP(-TAU(I)/U0)
      D2 = EXP(-TAU(I+1)/U0)
      F(I) = C1*W0(I)*((1.-U2)*(D1 - D2) + C2*(L2(I) - L4(I)))
      F(I+NL) = C1*W0(I)*((1.+U2)*(D1 - D2) + C2*(L2(I) + L4(I)))
   9  CONTINUE
C
      C1 = 0.75/DT
      C2 = 0.5*C1
      DO 10 I=1,NL
      DO 10 J=1,NL
      B(I,J) = -C1*W0(I)*(M1(I,J) - 2.*M3(I,J) + M5(I,J)
     2  + ALB*(L2(I) - L4(I))*(L2(J) - L4(J)))
      B(I,J+NL) = -C2*W0(I)*(M3(I,J) - M5(I,J)
     2  + ALB*(L2(I) - L4(I))*(L2(J) + L4(J)))
      B(I+NL,J) = - C1*W0(I)*(M3(I,J) - M5(I,J)
     2  + ALB*(L2(I) + L4(I))*(L2(J) - L4(J)))
      B(I+NL,J+NL) = - C2*W0(I)*(M1(I,J) + M5(I,J)
     2  + ALB*(L2(I) + L4(I))*(L2(J) + L4(J)))
  10  CONTINUE
C
      DO 11 I=1,NL2
  11  B(I,I) = B(I,I) + 1.
C
C ***** SOLVE THE MATRIX EQUATION *****
      CALL SGEFA(B,ML2,NL2,IPVT,INFO)
      IF(INFO.NE.0) PRINT 100,INFO,WAV
 100  FORMAT(1X,'INFO=',I3,5X,'WAV=',1PE11.4)
      CALL SGESL(B,ML2,NL2,IPVT,F,0)
C
C ***** COMPUTE THE NORMALIZED SOURCE FUNCTIONS AT THE ORIGINAL
C       GRID LEVELS *****
      C1 = 2.*ALB*U0*EXP(-TG/U0)
C
      DO 20 L=1,LZ
      K = NZ - L + 1
      TAK = TA(K)
      TU = TAK/U0
      IF(K.EQ.NZ.OR.TU.GT.0.02) GO TO 12
      S(K) = S(NZ)
      GO TO 20
  12  CONTINUE
C
      DO 13 I=1,NL
      IF(TAK.LT.TAU(I+1).AND.TAK.GT.TAU(I)) GO TO 14
      T1 = ABS(TAU(I) - TAK)
      T2 = ABS(TAU(I+1) - TAK)
      GO TO 31
  14  CONTINUE
      T1 = TAK - TAU(I)
      T2 = TAU(I+1) - TAK
  31  XI(I,1) = T1
      XI(I,2) = T2
      EX1(I,1) = E1(T1)
      EX1(I,2) = E1(T2)
  13  CONTINUE
C
      DO 32 J=1,2
      DO 32 I=1,NL
      X = XI(I,J)
      X2 = X*X
      X3 = X*X2
      X4 = X*X3
      PEX = EXP(-X)
      E1X = EX1(I,J)
      E2(I,J) = PEX - X*E1X
      E3(I,J) = (X2*E1X + PEX*(1. - X))/2.
      E4(I,J) = (-X3*E1X + PEX*(2. - X + X2))/6.
      E5(I,J) = (X4*E1X + PEX*(6. - 2.*X + X2 - X3))/24.
  32  CONTINUE
C
      DO 33 I=1,NL
      IF(TAK.LT.TAU(I+1) .AND. TAK.GT.TAU(I)) GO TO 34
      G1(I) = ABS(E2(I,1) - E2(I,2))
      G2(I) = ABS(E3(I,1) - E3(I,2))
      G3(I) = ABS(E4(I,1) - E4(I,2))
      G4(I) = ABS(E5(I,1) - E5(I,2))
      GO TO 33
  34  CONTINUE
      G1(I) = 2. - E2(I,1) - E2(I,2)
      G2(I) = 1. - E3(I,1) - E3(I,2)
      G3(I) = 2./3. - E4(I,1) - E4(I,2)
      G4(I) = 0.5 - E5(I,1) - E5(I,2)
  33  CONTINUE
C
      SUM1 = 0.
      SUM2 = 0.
      SUM3 = 0.
      SUM4 = 0.
      DO 15 I=1,NL
      J = I + NL
      SUM1 = SUM1 + (2.*F(I) + F(J))*G1(I)
      SUM2 = SUM2 + (2.*F(I) + F(J))*G2(I)
      SUM3 = SUM3 + (2.*F(I) - F(J))*G3(I)
      SUM4 = SUM4 + (2.*F(I) - F(J))*G4(I)
  15  CONTINUE
C
      DTG = TG - TAK
      D2 = EXP(-DTG) - DTG*E1(DTG)
      S(K) = EXP(-TAK/U0) + C1*D2 + 2.*(SUM1 - SUM3)
     2  + 4.*ALB*D2*(SUM2 - SUM4)
  20  CONTINUE
C
C ***** COMPUTE APPROXIMATE INTENSITIES AT LARGE OPTICAL DEPTHS *****
      IF(LZ.EQ.NZ) RETURN
      LZ1 = LZ + 1
      DO 24 L=LZ1,NZ
      K = NZ - L + 1
      DTAB = TAB(K) - TAB(MZ)
      S(K) = S(MZ)*EXP(-DTAB/U0)
  24  CONTINUE
C
      RETURN
      END

      real*8 FUNCTION SIGRAY(W)
      implicit real*8(A-H,O-Z)
      W1 = 1.E-4 * W
      W2 = W1 * W1
      W4 = W2 * W2
      SIGRAY = 4.006E-28*(1. + .0113/W2 + .00013/W4)/W4
      RETURN
      END

      real*8 FUNCTION E1(X)
      implicit real*8(A-H,O-Z)
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/EBLOK.inc'
      E1 = 0.
      IF(X.EQ.0.) RETURN
C
      X2 = X*X
      X3 = X2*X
      X4 = X3*X
      IF(X.GT.1) GO TO 1
      X5 = X4*X
      E1 = -LOG(X) + AI(1)*X + AI(2)*X2 + AI(3)*X3 + AI(4)*X4 + AI(5)*X5
     2  + AI(6)
      RETURN
C
   1  SUM1 = X4 + BI(1)*X3 + BI(2)*X2 + BI(3)*X + BI(4)
      SUM2 = X4 + CI(1)*X3 + CI(2)*X2 + CI(3)*X + CI(4)
      E1 = EXP(-X)*SUM1/SUM2/X
      RETURN
      END
