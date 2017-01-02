

      SUBROUTINE IR(T,PF,P,CGAS,FI)
c
c jfk 6/27/08 The call sequence now contains both PF and P. This was
c     being done incorrectly before, because PF was being converted
c     into P during the call.
c
c  This subroutine calculates the infrared flux

c  This subroutine contains CH4

      INCLUDE 'CLIMA/INCLUDE/header.inc'
      PARAMETER (NF=55,NGS=5)

C new common block, von Paris, 21/04/2006
      COMMON/IRDATA/WEIGHT(8,3),xkappa(8,12,NF,8,3), 
     & CIA(7,NF), CPRW(ND,NF)  
c-rr !3/23/11 put CIA matrix in IRDATA
      COMMON/IRBLK/FUPIR(ND),FDNIR(ND),SRFALBIR,OMG0AIR(NF,ND-1),
     & ASYAIR(NF,ND-1),IO3,QEXTIR(NF,ND-1)
c     COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4,
      COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4,FC2H6,FNO2
      COMMON/ALTBLOK/DALT(ND-1),RADIUS(ND-1),PARTICLES(ND),RAER(ND),
     2 ALT(ND)
      COMMON/CONSS/C,BK,G,PI,SM,DM, DM2
      COMMON/WAVE/AV(NF),LAM(NF),W(NF)
c
c jfk 6/25/08 Temporary (PF should have been passed in the call
c     statement. Actually, this statement does remain once the code is
c     corrected.)
      DIMENSION PF(ND), MS1(ND),MS(ND), FX(ND)

c     DIMENSION AV(NF)
c      DIMENSION P(ND),BPLANCK(ND),CPR(NF),TPR(NF),TPRIND(ND), 
c     2 TAUGCO2(ND),TAUGH2O(ND),TAUGCH4(ND),TAUGIR(ND),FUP(ND),FDN(ND),
c     3 FUPA(ND),FDNA(ND),T(ND),CGAS(ND,NGS),TAUCONTIN(ND),WEIGHT(8,3),
c     & xkappa(8,12,55,8,3)
      DIMENSION P(ND),BPLANCK(ND),CPR(NF),TPR(NF),TPRIND(ND), 
     2 TAUGCO2(ND),TAUGH2O(ND),TAUGCH4(ND),TAUGIR(ND),FUP(ND),FDN(ND),
     3 FUPA(ND),FDNA(ND),T(ND),CGAS(ND,NGS),TAUCONTIN(ND)
      REAL KAPPALAYER(NF,8,3,ND),kappa(NF,8,3),KAPPALAYEROZ(8,ND),
     2 WEIGHTOZC(8),KAPPAOZC(8),TAUGOZ(ND),TAUCONTINT(NF),CNUT,FI(4,ND),
     3 TAUTOTAL(NF),TRANSLAYER(ND),TWGHTT(8),OMG0IR(ND-1),ASYIR(ND-1),
     4 TAULAMIR(ND-1),TAUAEXTIR(ND-1),TAUASIR(ND-1),TAUSIR(ND-1)
      REAL KAPPALAYERC2H6(NF,6)
      DATA HP,SIGMA/6.63E-27, 5.67E-5/
      DIMENSION WEIGHTC2H6(6),CGASC2H6(ND)
      DIMENSION TAUGC2H6(ND)
      REAL KAPPAC2H6(NF)
      REAL np
c-rr Added CIAMS1L, CIAMSL, and CPRL and CIA Common block items 3/24/2011
      REAL TAUSUM(6), CIAMS1L, CIAMSL, CPRL, CIA, CPRW
      

C   PRESSURE-INDUCED CO2 ABSORPTION (FROM JIM POLLACK)
      DATA CPR/4.3E-5, 3.8E-5, 1.2E-5, 2.8E-6, 7.6E-7, 4.5E-7, 2.3E-7,
     2  5.4E-7, 1.6E-6, 10*0., 4.E-7, 4.E-6, 1.4E-5, 1.0E-5,  ! 4e-7 value updates 7.5e-7 in Kasting et al. (1984)
     3  1.2E-6, 2.0E-7, 5.0E-8, 3.0E-8, 28*0./
c-jdh No 20-um band
c      DATA CPR/0, 0, 0, 0, 0, 0, 0,
c     2  0, 0, 10*0., 4.E-7, 4.E-6, 1.4E-5, 1.0E-5,
c     3  1.2E-6, 2.0E-7, 5.0E-8, 3.0E-8, 28*0./
c-jdh No 7.7-um band
c      DATA CPR/4.3E-5, 3.8E-5, 1.2E-5, 2.8E-6, 7.6E-7, 4.5E-7, 2.3E-7,
c     2  5.4E-7, 1.6E-6, 10*0., 0, 0, 0, 0,
c     3  0, 0, 0, 0, 28*0./

      DATA TPR/3.4, 2.2, 1.9, 52*1.7/

C   FREQUENCIES AT ENDS OF SPECTRAL INTERVALS (1/CM)
c-jdh commented out, declared common in Clima.f, converted to units of (1/S)
c      DATA AV/40., 100., 160., 220., 280., 330., 380., 440., 495.,
c     2  545., 617., 667., 720., 800., 875., 940., 1000., 1065.,
c     3  1108., 1200., 1275., 1350., 1450., 1550., 1650., 1750., 1850.,
c     4  1950., 2050., 2200., 2397., 2494., 2796., 3087., 3425., 3760.,
c     5  4030., 4540., 4950., 5370., 5925., 6390., 6990., 7650., 8315.,
c     6  8850., 9350., 9650., 10400., 11220., 11870., 12790., 13300.,
c     7  14470., 15000./
C
c-jdh Ethane data (for polynomial fit approximation)
c      DATA KAPPAC2H6/0,0,0,0,0,0,0,0,0,0,
c     &  0,0,0,8.627161E-21,1.036469E-20,1.098662E-21,0,0,0,0,
c     &  0,2.377869E-22,9.640977E-21,1.726593E-20,1.426406E-21,0,
c     &  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
c     &  0,0,0,0,0,0,0,0,0,0/

c-jdh Ethane k-coefficient weights
      DATA WEIGHTC2H6/0.08566,0.18038,0.23396,0.23396,0.18038,0.08566/

       INTEGER COUNTERIR

       
       COUNTERIR = 0
       SRFALBIR = 0.
       NLAYERS = ND - 1    
c      BCON = 2.*HP/C/C
       HK = HP/BK
C       np = 1.E+1
 
        
c-jdh Set C2H6 absorption coefficients
      DO NI=1,NF
        DO NK=1,6
          KAPPALAYERC2H6(NI,NK) = 0.0       
         END DO
      END DO
      KAPPALAYERC2H6(14,1) = 1.0140E-23
      KAPPALAYERC2H6(14,2) = 6.6404E-23
      KAPPALAYERC2H6(14,3) = 4.4811E-22
      KAPPALAYERC2H6(14,4) = 2.2630E-21
      KAPPALAYERC2H6(14,5) = 5.1880E-21
      KAPPALAYERC2H6(14,6) = 7.8162E-21
      KAPPALAYERC2H6(15,1) = 4.1940E-21
      KAPPALAYERC2H6(15,2) = 5.3313E-21
      KAPPALAYERC2H6(15,3) = 7.0534E-21
      KAPPALAYERC2H6(15,4) = 9.2389E-21
      KAPPALAYERC2H6(15,5) = 1.3508E-20
      KAPPALAYERC2H6(15,6) = 2.7746E-20
      KAPPALAYERC2H6(16,1) = 2.3159E-23
      KAPPALAYERC2H6(16,2) = 8.6746E-23
      KAPPALAYERC2H6(16,3) = 2.6336E-22
      KAPPALAYERC2H6(16,4) = 9.1975E-22
      KAPPALAYERC2H6(16,5) = 2.3556E-21
      KAPPALAYERC2H6(16,6) = 3.9631E-21
      KAPPALAYERC2H6(22,1) = 3.0909E-24
      KAPPALAYERC2H6(22,2) = 1.5773E-23
      KAPPALAYERC2H6(22,3) = 4.9806E-23
      KAPPALAYERC2H6(22,4) = 1.3673E-22
      KAPPALAYERC2H6(22,5) = 3.3965E-22
      KAPPALAYERC2H6(22,6) = 8.2701E-22
      KAPPALAYERC2H6(23,1) = 1.0923E-21
      KAPPALAYERC2H6(23,2) = 2.7999E-21
      KAPPALAYERC2H6(23,3) = 7.0277E-21
      KAPPALAYERC2H6(23,4) = 1.1359E-20
      KAPPALAYERC2H6(23,5) = 1.5960E-20
      KAPPALAYERC2H6(23,6) = 2.1332E-20
      KAPPALAYERC2H6(24,1) = 6.0680E-21
      KAPPALAYERC2H6(24,2) = 9.5224E-21
      KAPPALAYERC2H6(24,3) = 1.3251E-20
      KAPPALAYERC2H6(24,4) = 1.6683E-20
      KAPPALAYERC2H6(24,5) = 2.1862E-20
      KAPPALAYERC2H6(24,6) = 4.6329E-20
      KAPPALAYERC2H6(25,1) = 3.1607E-23
      KAPPALAYERC2H6(25,2) = 7.4713E-23
      KAPPALAYERC2H6(25,3) = 2.2452E-22
      KAPPALAYERC2H6(25,4) = 9.2931E-22
      KAPPALAYERC2H6(25,5) = 3.0558E-21
      KAPPALAYERC2H6(25,6) = 5.7388E-21      
      DO I=1,NLAYERS
c        CGASC2H6(I) = 1.3386e-3*BK*273.16/(SM*G)
c     &                * (P(I+1) - P(I))*FC2H6/DM
c        CGASC2H6(I) = 1.0e-5*BK*273.16/(SM*G)
c     &                * (P(I+1) - P(I))*FC2H6/DM
c jfk 6/27/08 P was changed to PF in the 2 lines below.
        CGASC2H6(I) = 1.0e-5*BK*273.16/(SM*G)
     &                * (PF(I+1) - PF(I))*FC2H6/DM
      END DO


C Read the IR exponential sums
c-jdh Moved to Clima.f 
c      CALL IREXPSUMS(WEIGHT,xkappa)

      DO 7 IL = 1,NLAYERS
c      CALL INTERP(T(IL),P(IL),xkappa,kappa)
c-rr   
       TIL = AMIN1(T(IL), 350.)
       TIL = AMAX1(TIL, 50.)
       CALL INTERP(TIL,P(IL),kappa)
c       if (IL.eq.95)print *, 'TIL=', TIL, 'kappa=', kappa

c       TIL = AMIN1(T(IL), 400.)
c       TIL = AMAX1(TIL, 100.)


c          iil= 95
c       if (IL.eq.iil)then
c      print *, 'T(IL)= ',T(IL)
c       endif               

       CALL INTERPCO2CIA(T(IL),IL,MS,MS1,FX)
c       CALL INTERPCO2CIA(TIL,IL,MS,MS1,FX)
       
       DO 8 I=1, 55
       DO 9 J=1, 8
       DO 10 K=1, 3
       KAPPALAYER(I,J,K,IL) = kappa(I,J,K)
  10  CONTINUE
  9   CONTINUE
  8   CONTINUE
	
  7   CONTINUE
c      print *, 'MS=', MS, 'MS1=', MS1
      DO IL = 1,NLAYERS
       CALL INTERPOZONE(P(IL),WEIGHTOZC,KAPPAOZC)
        DO J = 1, 8
        KAPPALAYEROZ(J,IL) = KAPPAOZC(J)
        ENDDO
      ENDDO
c       PRINT*, T,P
       DO 99 I=1, ND
        FUPIR(I) = 0.
        FDNIR(I) = 0.
 99   CONTINUE
       DO I=1,55
         TAUTOTAL(I) = 0.0
       ENDDO

****** Loop over frequency      
       DO 1 I=1, NF
       DO 20 J=1,ND 
C-jdh  VAC = C*AV(I)
       VAC = AV(I)
       BPLANCK(J) = PLANCK(VAC,T(J),HP,C,HK)
  20   CONTINUE
       DO 15 J=1,ND
         FUPA(J) = 0.0
         FDNA(J) = 0.0
         FUP(J) = 0.0
         FDN(J) = 0.0
         TRANSLAYER(J) = 0.0 
  15   CONTINUE
****** AEROSOLS

       DO IL=1,NLAYERS
        r = RAER(IL)
        np = PARTICLES(IL)
        TAUAEXTIR(IL) = QEXTIR(I,IL)*PI*r*r*DALT(IL)*np*1.E+5
        TAUASIR(IL) = OMG0AIR(I,IL)*TAUAEXTIR(IL)
       ENDDO
       IF (I.EQ.15) THEN
       TAUEXTIRTOTAL = 0.
       TAUASIRTOTAL = 0.
       DO IL=1,NLAYERS
       TAUEXTIRTOTAL = TAUEXTIRTOTAL + TAUAEXTIR(IL)
       TAUASIRTOTAL = TAUASIRTOTAL + TAUASIR(IL)
       ENDDO
c       PRINT *, '******************************'
c       PRINT *, 'TAUEXTIRTOTAL'
c      PRINT *, TAUEXTIRTOTAL
c       PRINT *, 'TAUSIRTOTAL'
c       PRINT *, TAUASIRTOTAL
       TAUAABSIRTOTAL = TAUEXTIRTOTAL - TAUASIRTOTAL
c       PRINT *, 'TAUAABSIRTOTAL'
c       PRINT *, TAUAABSIRTOTAL
c       PRINT *, '*******************************'
       ENDIF 
****  8 - 12 UM CONTINUUM
       IF ((I.GT.14).and.(I.LT.21)) THEN
C       PRINT *, 'TAUCONTIN'

       DO IL = 1,NLAYERS
c       DP = ABS(P(IL+1)-P(IL))*1.E6
c jfk 6/27/08 P was changed to PF in the line below.
       DP = ABS(PF(IL+1)-PF(IL))*1.E6
       TV = 1.25E-22 + 1.67E-19 * EXP(-2.62E-13*VAC)
       TF = EXP(1800. * (1./T(IL) - 1./296.))
       CNUT = TV*TF
c jfk 6/26/08 P should be OK in the next line because this is in the
c     middle of the layer
       TAUCONTIN(IL) = CNUT*FI(1,IL)*FI(1,IL)*P(IL)*DP/(DM*SM*G)
C       PRINT 19, TAUCONTIN(IL) 
       ENDDO
       ELSE
       DO IL = 1,NLAYERS
       TAUCONTIN(IL) = 0.
       ENDDO
       ENDIF
C******
C
c-jdh  If(I.lt.55) VAC1 = C*(AV(I+1)-AV(I))
c-jdh  If(I.eq.55) VAC1 = C*(AV(I) - AV(I-1))
c       IF(I.lt.NF) THEN
c         VAC1 = (AV(I+1)-AV(I))
c       ELSE
c         VAC1 = (AV(I)-AV(I-1))
c       ENDIF

       DO 2 K1=1, 8 
       DO 3 K2=1, 8 
       
       	
c	DO 4 K3=1, 6 ! This is methane loop
c       DO 5 K0=1, 6
c         TWGHT = WEIGHT(K1,1)*WEIGHT(K2,2)*WEIGHT(K3,3)
       
c        TWGHT = WEIGHT(K1,1)*WEIGHT(K2,2)*WEIGHT(K3,3)*WEIGHTC2H6(K0)
       	
        TWGHT = WEIGHT(K1,1)*WEIGHT(K2,2)

        SUM_TPRIND=0   ! summing TPRIND experiment. Sets SUMTPRIND when start a new frequency (I)
        DO 11 IL = 1,NLAYERS
         PPE = (1. + 0.3*FI(2,IL))*P(IL)     !CO2
         TPE = (300./T(IL))**TPR(I)
         TAUGCO2(IL) = KAPPALAYER(I,K1,1,IL)*CGAS(IL,3)/(44.*SM)
         TAUGH2O(IL) = KAPPALAYER(I,K2,2,IL)*CGAS(IL,2)/(18.*SM)
c-jdh    TAUGCH4(IL) = KAPPALAYER(I,K3,3,IL)*CGAS(IL,5)
c    & *(4.46*DM)/(16.*SM)
         TAUGCH4(IL) = KAPPALAYER(I,K3,3,IL)*CGAS(IL,5)*2.687E24
c         TAUGC2H6(IL) = KAPPALAYERC2H6(I,K0)*CGASC2H6(IL)*2.687E24
         TAUGC2H6(IL) = 0.
c-jdh    TAUGC2H6(IL) = KAPPAC2H6(I)*CGAS(IL,5)*2.687E24*(FC2H6/FCH4)

c-rr	This is the CO2 CIA loop that takes the saved values for FX, MS, and MS1 at a given height and calculates
c-rr	the correct interpolated CIA values. These are from the model of Wordsworth et al. (2010) that uses the GBB
c-rr	parametrization scheme. 3/24/11
c        print *, 'is this working?'
        
	CIAMS1L = alog(CIA(MS1(IL),I))
	CIAMSL = alog(CIA(MS(IL),I))
	CPRL = CIAMS1L + FX(IL)*(CIAMSL-CIAMS1L)
	CPRW(IL,I) = exp(CPRL)
       IF (CPRW(IL,I).lt.1.E-45)CPRW(IL,I)= 0. !c-rr Ensures bands where CIA is supposed to be zero, are zero. 5/28/2011
       
c        if ((IL.eq.iil).and.(I.eq.4).and.(k1.eq.2).and.(k2.eq.2)) 
c     & print *, 'CIA=',CPRW(IL,I), 'LAMB=',I

       IF(I .EQ. 38)THEN
          CPRW(IL,I) = 3.5E-8  !Rewriting CIA at intervals 37 and 40 with appropriate CIA at 2.3 and 1.73 microns, respectively(Tsang et al. 2008).
       ELSEIF (I.EQ.41)THEN
          CPRW(IL,I) = .57*6.0E-9  !Because the 1.73 micron window is only 57% of this bin
       ELSEIF ((I.EQ.45).OR.(I.EQ.46).OR.(I.EQ.47))THEN  !Rewriting CIA at intervals 45-47 with 1.2 micron complex
	      CPRW(IL,I) = 1.5E-9
       ENDIF
       

c	print *, ' I hope so'
****** Pressure induced absorption by CO2 
         CGAS1 = CGAS(IL,3)/1.963E-3  
c-rr 3/25/11 commenting out Kasting et al. (1984) parametrization
!        TPRIND(IL) = CPR(I)*PPE*TPE*CGAS1
!         TPRIND(IL) = 0.

c	 print *, TPRIND(IL)
c-rr 3/24/11 	 PUT NEW TPRIND parametrization here! 
          PCGS = P(IL)* 1.e6
!       TPRIND(IL) = CPRW(IL,I)*(PCGS/(BK*T(IL)*2.687E19))*CGAS1
!    &   *((1/1.3)+(FI(2,IL)/1.3))

!  c-rr Corrected bug in CO2 CIA. Removed extraneous pressure broadening factor of (1/1.3). 2/7/2012

        TPRIND(IL) = CPRW(IL,I)*(PCGS/(BK*T(IL)*2.687E19))*CGAS1
     &   *(1+FI(2,IL))

!	  if ((K1.eq.2).and.(K2.eq.2).and.(I.eq.4))then
!             print *, TPRIND(IL)
!          endif

c          if ((IL.eq.iil).and.(K1.eq.2).and.(K2.eq.2).and.(I.eq.4))then
c	print *, 'PCGS=', PCGS, 'T(IL)=', T(IL), 'FI=', FI(2,IL), 
c     &   'TPRIND(IL)=', TPRIND(IL), 'BK=', BK, 'CGAS1=', CGAS1, 
c     &   'CPRW=', CPRW(IL,I)
c	endif 

c	SUM_TPRIND = SUM_TPRIND + TPRIND(IL) !summing TPRIND experiment
        
         
c          print *, TPRIND(IL)
c        TPRIND(IL) = 0
*******
         TAUGIR(IL) = TAUGCH4(IL)+TAUGCO2(IL)+TAUGH2O(IL)+TPRIND(IL)
     & +TAUCONTIN(IL)+TAUGC2H6(IL)
  11    CONTINUE
	
c	print *, AV(I), SUM_TPRIND
c        read(*,*)

******* Ozone absorption
        IF (I.EQ.18) THEN
      DO K4=1,8
        TWGHTT(K4) = TWGHT*WEIGHTOZC(K4)
      DO IL=1,NLAYERS
        TAUGOZ(IL) =KAPPALAYEROZ(K4,IL)*CGAS(IL,4)*2.687E19
c       TAUGOZ(IL) =KAPPALAYEROZ(K4,IL)*CGAS(IL,4)*(4.46E-5*DM)/(48.*SM)
c -PJK There was a bug in the line of code below. The IR optical depth was
c      being compounded improperly in the ozone band. It's done right below.
c      TAUGIR(IL) =TAUGIR(IL) + TAUGOZ(IL)
      TAUGIR(IL) = TAUGCH4(IL)+TAUGCO2(IL)+TAUGH2O(IL)+TPRIND(IL)
     & +TAUCONTIN(IL) + TAUGOZ(IL) + TAUGC2H6(IL)
      ENDDO
      DO IL=1, NLAYERS
          TAULAMIR(IL) = TAUAEXTIR(IL) + TAUGIR(IL)
          TAUSIR(IL) = TAUASIR(IL) 
          ASYIR(IL) = ASYAIR(I,IL) 
          OMG0IR(IL) = TAUSIR(IL)/TAULAMIR(IL) 
          OMG0IR(IL) = AMIN1(OMG0IR(IL),0.99999)
          OMG0IR(IL) = AMAX1(OMG0IR(IL),1.E-5)
      ENDDO
c
c jfk  6/25/08  Include TAUTOP in the call sequence to do the upper BC.
c      This requires the use of PF, not P.
      TAUTOP = TAULAMIR(1)*PF(1)/(PF(2)-PF(1))
      TAUTOP = AMIN1(TAUTOP,1.)
      CALL DELTA2STRIR(SRFALBIR,ASYIR,TAULAMIR,OMG0IR,
     & FUP,FDN,BPLANCK,TAUTOP)
C
               J = 1
  13          IF (J .LE. ND) THEN
                  FUPA(J)=FUPA(J)+TWGHTT(K4)*FUP(J)
                  FDNA(J)=FDNA(J)+TWGHTT(K4)*FDN(J)
                  J=J+1
                  GOTO 13
               END IF
        

        ENDDO
        GOTO 5
        ENDIF
***************
      DO IL=1, NLAYERS
          TAULAMIR(IL) = TAUAEXTIR(IL) + TAUGIR(IL)
          TAUSIR(IL) = TAUASIR(IL) 
          ASYIR(IL) = ASYAIR(I,IL) 
          OMG0IR(IL) = TAUSIR(IL)/TAULAMIR(IL) 
          OMG0IR(IL) = AMIN1(OMG0IR(IL),0.99999)
          OMG0IR(IL) = AMAX1(OMG0IR(IL),1.E-5)
      ENDDO
c
c jfk  6/25/08  Include TAUTOP in the call sequence to do the upper BC
      TAUTOP = TAULAMIR(1)*PF(1)/(PF(2)-PF(1))
      TAUTOP = AMIN1(TAUTOP,1.)
      CALL DELTA2STRIR(SRFALBIR,ASYIR,TAULAMIR,OMG0IR,
     & FUP,FDN,BPLANCK,TAUTOP)
C
C
               J = 1
  16          IF (J .LE. ND) THEN
                  FUPA(J)=FUPA(J)+TWGHT*FUP(J)
                  FDNA(J)=FDNA(J)+TWGHT*FDN(J)
                  J=J+1
                  GOTO 16
               END IF
C

  5     CONTINUE
  4     CONTINUE  
                  
  3     CONTINUE            
  2     CONTINUE

c         write(89,*) AV(I)/C, SUM_TPRIND   !summing tprind experiment
	

         DO 14 J=1,ND
             FDNIR(J)=FDNIR(J)+W(I)*FDNA(J)
             FUPIR(J)=FUPIR(J)+W(I)*FUPA(J)
             IF (J.eq.1) THEN
             write(90, *)AV(I)/C,ABS(FDNIR(J) - FUPIR(J)) ! Outgoing IR at top of atmosphere   
             ENDIF
c-jdh uncomment to get band-by-band fluxes            
c             IF (J.EQ.1) THEN
c  17            FORMAT(1PE13.5,3x,1PE22.10)
c               PRINT 17, AV(I), W(I)*FUPA(J)
c             END IF

c-jdh        FDNIR(J)=FDNIR(J)+VAC1*FDNA(J)
c-jdh        FUPIR(J)=FUPIR(J)+VAC1*FUPA(J)
             BPLANCK(J) = 0.    
  14    CONTINUE
C        PRINT*, 'TRANSLAYER'
C        DO J=1,NLAYERS
C        PRINT 19, TRANSLAYER(J)
C         TAUTOTAL(I) = TAUTOTAL(I) - alog(TRANSLAYER(J))
C        ENDDO
   1     CONTINUE                   !**END LOOP over frequency**

19     FORMAT(/1X,1PE10.4)
100    FORMAT(1X,1P10E12.5) 
C      PRINT*,'TAUTOTAL'
C      PRINT 100, TAUTOTAL
C
      RETURN
      END
