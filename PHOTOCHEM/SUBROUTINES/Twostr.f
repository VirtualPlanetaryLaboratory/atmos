      SUBROUTINE TWOSTR(SIGR,U0,sq,WAV,LL,S,NN,IKN)

!note i haven't yet done the work to fully abstract particles in this subroutine.
!hydrocarbons are using optical properties files, while sulfate and sulfur are hardcoded.
      
! in call, from Photo, send WAV(L) in the WAV position, and L in the LL position
!IKN is a printing option, NN is just N as sent to Photo
!S is the returned source function

C   This is my version of the Toon et al. 2-stream code.  (Ref.: JGR
C   94, 16287, 1989).  It vectorizes over height, rather than wavelength,
C   and is designed to work with PRIMS3 and its companion photochemical 
C   models.
C
C   For now, at least, it is hardwired as the quadrature approximation.
C   NP is the number of different types of particles

      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      PARAMETER(NZP1=NZ+1,NZ2=2*NZ)

      implicit real*8(A-H,O-Z)
      real*8 mass
      character*8 PLANET
      REAL*8 sq(kj,nz,kw),S(NZ)
 
      DIMENSION TAU(NZ),TAUCTSTR(NZP1),GT(NZ),GAM1(NZ),GAM2(NZ),
     1  GAM3(NZ),GAM4(NZ),ALAM(NZ),CGAM(NZ),E1(NZ),E2(NZ),E3(NZ),E4(NZ),
     2  CP0(NZ),CPB(NZ),CM0(NZ),CMB(NZ),Y1(NZ),Y2(NZ),W0(NZ),
     3  TAUSG(NZ),TAUSP(NZ),DIRECT(NZP1),AMEAN(NZP1),TAUG(NZ),FMT(NZ)
      DIMENSION A(NZ2),B(NZ2),D(NZ2),E(NZ2),FUP(NZP1),FDN(NZP1)
      DIMENSION W0P(NP),QEXT(NP),TAUSCAT_PART(NP,NZ),TAUP(NZ)
      dimension TAU_PART(NP,NZ)
      dimension SIGR(NZ)
!this FUP is different from the lower boundary flux terms in Output.f.  Should be fine here.

      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PHOTABLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/CBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/AERBLK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/MBLOK.inc'


C     U1 = 0.5  (Eddington value)
      SQ3 = SQRT(3.)
      PI = 3.14159
      U1 = 1./SQ3
      U0M = 1./U0
      U0M2 = U0M*U0M
      U1M = 1./U1
      NZM1 = NZ - 1
      MZ2 = NZ2
      GP = 0.8  !old

C
C   Calculate the optical depths of the different layers.  TAUA is absorption,
C   TAUSG is scattering by gases, TAUSP is scattering by particles, TAUG is
C   extinction due to gases, TAU is total extinction due to gases and
C   particles.  
C   Note that the grid for this subroutine is numbered from top to bottom,
c   whereas the main program is numbered from bottom to top.
C   First do gases 

      do I=1,NZ
      N = NZP1 - I
      TAUSP(N) = 0.
      TAUP(N) = 0.
      TAUSG(N) = SIGR(I)*DEN(I)*DZ(I)

! Can below be removed?
corig      TAUA = ( SO3(I)*O3(I) + SO2(I)*O2(I) + SCO2*CO2(I) + SH2O*H2O(I)
corig     2   + SSO2*FSO2(I) + SS2*S2(I) + SH2S*H2S(I) + SNH3*FNH3(I) ) 
corig     3   * DEN(I)*DZ
corig   TAUG(N) = TAUA + TAUSG(N)

!Sspecies(I) is the species cross section at the given height 
!Fspecies are absolute values of species mixing ratios
!so, in my scheme where everything absorbs, this would be:
!a loop over kj (to get all species that absorb (at each height)) - 
!of sq*absorbers
  
      taua=0.0
       do j=1,kj
          IF(PLANET.EQ.'WASP12B') THEN
           if(j.lt.9)taua = taua + sq(j,i,LL)*absorbers(j,i)
          ELSE
           taua = taua + sq(j,i,LL)*absorbers(j,i)
          ENDIF   
       enddo   

      TAUG(N) = TAUA*DEN(I)*DZ(I) + TAUSG(N)

      enddo

C   Now do particles.  Must combine their scattering optical depths into a
C   single array in order to calculate W0 and G.  (TAUSP(N))
C   Don't need any arrays for pure absorption.
C   Scale optical depth, W0, and G for the particles using the Delta-
C   Eddington approximation of Joseph et al. (Ref: J. Atmos. Sci. 33,
C   2452, 1976)

!gna - this is the correction for forward scattering: delta eddington approximation

      IF (NP.GT.O) THEN !REMOVE THE PARTICLE STUFF, YET AGAIN....
C   Particle 1 is sulfate, 2 is S8, 3 is HCAER
       DO K=1,NP
        W0P(K) = 0.0
        IF (K.LT.3) W0P(K) = 1.
        IF (K.EQ.2.AND.WAV.LE.3500.) W0P(K) = 0.5 !ACK hardcoded grid
 !      print*, 'K,W0P=', K,W0P
       ENDDO    
!here, Kevin uses a formula to compute WOPS and QEXT (fom G0=0.8 and WOP)

      		DO  J=1,NP
      		DO  I=1,NZ
      			N = NZP1 - I
C-MC Here, we are including correct value of Qext and W0P for hydrocarbons
      			IF (J .GE. 3) THEN 
       !hardcoded for HCAER and HCAER2  
       				QEXT(J) = QEXTT(LL,I,J)     
       				W0P(J) = W0T(LL,I,J)                          
      			ELSE
       !for sulfate and S8 : valid for large particles only, but this is status quo for now...
       				QEXT(J)=2.   
       !W0P set above
       !for sulfate and S8, set assymetry factor to GP=0.8 (alt and wl independent)
       				GFT(LL,I,J)= 0.8     
      			ENDIF

      		 TAUP1 = QEXT(J)*PI*RPAR(I,J)*RPAR(I,J)*AERSOL(I,J)*DZ(I)
      !particle extinction for each particle at each height
      !not used, but might as well keep  
             TAU_PART(J,I)=TAUP1  
      !particle scattering for each particle at each height
      		 TAUSCAT_PART(J,I)=W0P(J)*TAUP1  
      !TAUSP contains total particle scattering at each height
      		 TAUSP(N) = TAUSP(N) + TAUSCAT_PART(J,I)   
      !TAUP contains total particle extinction at each height 
      		 TAUP(N) = TAUP(N) + TAUP1  
      		ENDDO
      		ENDDO
      ENDIF

C   Calculate W0 and G by averaging over Rayleigh and Mie scatterers.  
C   (scattering due to gases vs. particles)
C   Avoid letting W0 equal exactly 1.
      DO  N=1,NZ
      !TAU is total extinction (Gas + Rayleigh + Particle)
      TAU(N)=TAUG(N) + TAUP(N) 
      !w0 = total scattering/total extinction
      W0(N) = (TAUSG(N) + TAUSP(N))/TAU(N)  
      W0(N) = AMIN1(W0(N), 0.99999) 
c-mab: Was "0.999" before - extra 99s needed for WASP12B convergence

      !GFT is still bottom to top, so needs a switch
      I= NZP1 - N  

      GPnew=0.0
      IF (NP.GT.O) THEN !REMOVE THE PARTICLE STUFF, YET AGAIN....
      do k=1,np
      !asymmetry factors weighted by particle scattering/total scattering
      GPnew=GPnew + GFT(LL,I,K)*TAUSCAT_PART(K,I)/(TAUSP(N) + TAUSG(N)) 
      enddo  
 
      !dont let assymetry factor get larger than 1
      	GT(N) = amin1(Gpnew,0.99999)
      ELSE
      	GT(N) = 0.0   
      ENDIF 

      enddo


C   Delta-Eddington scaling
C-AP I used approximation from Joseph et al. 1976
      DO 24 N=1,NZ
      FMT(N) = GT(N)*GT(N) 
      TAU(N) = TAU(N)*(1. - W0(N)*FMT(N))
      W0(N) = W0(N)*(1. - FMT(N))/(1. - W0(N)*FMT(N))
  24  GT(N) = GT(N)/(1. + GT(N))
C-AP**************************************************
C
C   Calculate the gamma's, lambda's, and e's
      DO 2 N=1,NZ
C     GAM1(N) = (7. - W0(N)*(4.+3.*GT(N)))/4.
C     GAM2(N) = - (1. - W0(N)*(4.-3.*GT(N)))/4.
C     GAM3(N) = (2. - 3.*GT(N)*U0)/4.
C   (Eddington values above; quadrature values below)
C
      GAM1(N) = SQ3*(2. - W0(N)*(1.+GT(N)))/2.
      GAM2(N) = SQ3*W0(N)*(1.-GT(N))/2.
      GAM3(N) = (1. - SQ3*GT(N)*U0)/2.
      IF(NP.EQ.0)GAM3(N) = (1. - SQ3*U0)/2.
      GAM4(N) = 1. - GAM3(N)

      ALAM(N) = SQRT(GAM1(N)*GAM1(N) - GAM2(N)*GAM2(N))
      CGAM(N) = (GAM1(N) - ALAM(N))/GAM2(N)
      EMLT = EXP(-ALAM(N)*TAU(N))

      E1(N) = 1. + CGAM(N)*EMLT
      E2(N) = 1. - CGAM(N)*EMLT
      E3(N) = CGAM(N) + EMLT
   2  E4(N) = CGAM(N) - EMLT

C   Calculate A, B, and D, i.e. the coefficients of the tridiagonal matrix 
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

C   Even coefficients
      DO 4 N=1,NZM1
      L = 2*N
      A(L) = E2(N+1)*E1(N) - E3(N)*E4(N+1)
      B(L) = E2(N)*E2(N+1) - E4(N)*E4(N+1)
   4  D(L) = E1(N+1)*E4(N+1) - E2(N+1)*E3(N+1)

C   Bottom of atmosphere
      A(NZ2) = E1(NZ) - ALB*E3(NZ)
      B(NZ2) = E2(NZ) - ALB*E4(NZ)
      D(NZ2) = 0.

C   Now, set up the RHS of the equation:
C   TAUCTSTR(N) is the optical depth above layer N
      TAUCTSTR(1) = 0.
      DO 5 N=2,NZP1
   5  TAUCTSTR(N) = TAUCTSTR(N-1) + TAU(N-1)
C On last call,
C Print out TAUCTSTR(N), TAU(N), W0(N), GT(N), TAUSG(N), TAUSP(N).
C Also print TAUC at the ground for all wavelengths.

      IF(NN.EQ.1 .AND. IKN.EQ.1) THEN
      !ACK - check LUN
      WRITE(22,114) WAV,TAUCTSTR(NZP1)   
 114  FORMAT(1X,F6.1,2X,1PE10.3)
      !ACK - what are these hardcoded wavelengths???
      IF (WAV.EQ.1860.5 .OR. WAV.EQ.2010. .OR. WAV.EQ.2116.5 .OR.  
     2    WAV.EQ.2211. .OR. WAV.EQ.2312.5 .OR. WAV.EQ.2516. .OR.
     3    WAV.EQ.3007.5 .OR. WAV.EQ.3900. .OR. WAV.EQ.4500.) THEN
      !ACK - check LUN
      WRITE(20,106)WAV,U0,ALB    
 106  FORMAT('# WAV = ',F6.1,2X,'U0 = ',F6.4,2X,'Rsfc = ',F5.3,
     2   /'# TWOSTR: TAUCTSTR(N) is the optical depth above layer N',
     3   /'#  Z',4X,'TAUCTSTR(N)',4X,'TAU(N)',5X,'W0(N)',6X,'G(N)',6X,
     4   'TAUSG(N)',3X,'TAUSP(N)',3X,'N')
      DO 113 N=1,NZ
      I = NZP1 - N
      !ACK - check LUN
      WRITE(20,107) I,TAUCTSTR(N),TAU(N),W0(N),GT(N),TAUSG(N),TAUSP(N),N
 107  FORMAT(1X,I3,1P6E11.3,1X,I3)
 113  CONTINUE
      !ACK - check LUN
      WRITE(20,105)TAUCTSTR(NZP1)   
 105  FORMAT('#  0',1PE11.3)
      ENDIF
      ENDIF

C   DIRECT(N) is the direct solar flux at the top of layer N.  Values
C   are normalized to unity.  DIRECT(NZP1) is the direct flux at the ground.
      DIRECT(1) = U0
      DO 6 N=1,NZ
      ET0 = EXP(-TAUCTSTR(N)/U0)
      ETB = ET0 * EXP(-TAU(N)/U0)
      DIRECT(N+1) = ETB*U0
      DENOM = ALAM(N)*ALAM(N) - U0M2
      FACP = W0(N) * ((GAM1(N)-U0M)*GAM3(N) + GAM4(N)*GAM2(N))
      FACM = W0(N) * ((GAM1(N)+U0M)*GAM4(N) + GAM2(N)*GAM3(N))

      CP0(N) = ET0*FACP/DENOM
      CPB(N) = ETB*FACP/DENOM
      CM0(N) = ET0*FACM/DENOM
   6  CMB(N) = ETB*FACM/DENOM
      SSFC = ALB*DIRECT(NZP1)

C   Odd coefficients
      E(1) = - CM0(1)
      DO 7 N=1,NZM1
      L = 2*N + 1
   7  E(L) = (CP0(N+1)-CPB(N))*E3(N) + (CMB(N)-CM0(N+1))*E1(N)

C   Even coefficients
      DO 8 N=1,NZM1
      L = 2*N
   8  E(L) = (CP0(N+1)-CPB(N))*E2(N+1) - (CM0(N+1)-CMB(N))*E4(N+1)
      E(NZ2) = SSFC - CPB(NZ) + ALB*CMB(NZ)

C   Call the tridiagonal solver (from LINPACK).  E is the RHS of the matrix
C   equation on input and is the solution vector Y on output
      CALL SGTSL(MZ2,A,B,D,E,NFLAG)
      IF (NFLAG .NE. 0) PRINT 100, NFLAG
 100  FORMAT(/1X,'Tridiagonal solver failed in TWOSTR, NFLAG =',I4)

      DO 9 N=1,NZ
      L = 2*N
      L1 = L-1
      Y1(N) = E(L1)
   9  Y2(N) = E(L)

C   Calculate the mean intensities, AMEAN(N), at the boundaries between
C   the layers.  AMEAN(N) is the intensity at the top of layer N.
      AMEAN(1) = U1M * (Y1(1)*E3(1) - Y2(1)*E4(1) + CP0(1)) + 1.
      DO 10 N=1,NZ
  10  AMEAN(N+1) = U1M * (Y1(N)*(E1(N)+E3(N)) + Y2(N)*(E2(N)+E4(N)) 
     1  + CPB(N) + CMB(N)) + DIRECT(N+1)/U0

C   Reset any AMEAN values that may go negative.  Check error file
C   to be sure this only happens near the ground where AMEAN ~ 0.
      DO 12 N=1,NZP1 
      IF(AMEAN(N).LT.0.0)THEN
         !ACK - check LUN
         WRITE(13,103) WAV,N,AMEAN(N)   
 103     FORMAT('WAVE =',F6.1,' AMEAN(',I3,')=',1PE11.3)
         AMEAN(N) = ABS(AMEAN(N))
      ENDIF
  12  CONTINUE

C  Calculate upward and downward fluxes.
 
      FUP(1) = ((Y1(1)*E3(1) - Y2(1)*E4(1)) + CP0(1))
      FDN(1) = DIRECT(1)
      DO 110 N=1,NZ
         FUP(N+1) = (Y1(N)*E1(N) + Y2(N)*E2(N)
     &     + CPB(N))
         FDN(N+1) = (Y1(N)*E3(N) + Y2(N)*E4(N)
     &     + CMB(N)) + DIRECT(N+1)
 110  CONTINUE

C   Convert back to main program grid.  S(I) is the mean intensity at the
C   midpoint of layer I.
      DO 11 I=1,NZ
      N = NZP1 - I
  11  S(I) = SQRT(AMEAN(N)*AMEAN(N+1))

C  Print out the results at a few wavelengths

      IF(NN.EQ.1 .AND. IKN.EQ.1) THEN
      IF (WAV.EQ.1860.5 .OR. WAV.EQ.2010. .OR. WAV.EQ.2116.5 .OR.
     2    WAV.EQ.2211. .OR. WAV.EQ.2312.5 .OR. WAV.EQ.2516. .OR.
     3    WAV.EQ.3007.5 .OR. WAV.EQ.3900. .OR. WAV.EQ.4500.) THEN
      !ACK - check LUN
      WRITE(21,108)WAV,U0,ALB  
 108  FORMAT(/'# WAV = ',F6.1,2X,'U0 = ',F6.4,2X,'Rsfc = ',F5.3,
     2  /'#  Z',4X,'S(Z)',7X,'Y1',9X,'Y2',
     3  9X,'DIRECT',5X,'AMEAN',6X,'FUP',8X,'FDN',7X,'N')
      DO 111 N=1,NZ
      I = NZP1 - N
      !ACK - check LUN
      WRITE(21,109) I,S(I),Y1(N),Y2(N),DIRECT(N),AMEAN(N), 
     2  FUP(N),FDN(N),N
 109  FORMAT(1X,I3,1P7E11.3,1X,I3)
 111  CONTINUE
      !ACK - check LUN
      WRITE(21,112) DIRECT(NZP1),AMEAN(NZP1),FUP(NZP1),FDN(NZP1) 
 112  FORMAT('#',36X,1P4E11.3)
      ENDIF
      ENDIF

      RETURN
      END
C-PK ***************************************
