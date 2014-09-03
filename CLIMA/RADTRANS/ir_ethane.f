C
C
      SUBROUTINE IR(T,PF,P,CGAS,FI)
c
c jfk 6/27/08 The call sequence now contains both PF and P. This was
c     being done incorrectly before, because PF was being converted
c     into P during the call.
c
c  This subroutine calculates the infrared flux

c  This subroutine contains CH4

      INCLUDE '../header.inc'
      PARAMETER (NF=55,NGS=5, NZ=100)

C new common block, von Paris, 21/04/2006
      COMMON/IRDATA/WEIGHT(8,3),xkappa(8,12,55,8,3)
      COMMON/IRBLK/FUPIR(ND),FDNIR(ND),SRFALBIR,OMG0AIR(NF,ND-1),
     & ASYAIR(NF,ND-1),IO3,QEXTIR(NF,ND-1)
c     COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4,
      COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4,FC2H6,FNO2
      COMMON/ALTBLOK/DALT(ND-1),RADIUS(ND-1),PARTICLES(ND),RAER(ND),
     2 ALT(ND)
      COMMON/CONSS/C,BK,G,PI,SM,DM 
      COMMON/WAVE/AV(NF),LAM(NF),W(NF)
c
      COMMON/TFBLOK11/FNO2_PHOTO(NZ), FSO2_PHOTO(NZ)
     | , FNO2_NEW(ND), FSO2_NEW(ND), FCH4_PHOTO(NZ)
     | , FCH4_NEW(ND), FC2H6_PHOTO(NZ), FC2H6_NEW(ND)
C
      REAL KAPPALAYERC2H6(NF,6), KAPPALAYERSO2(NF,6)
C      COMMON/TFBLOK12/KAPPALAYERC2H6(NF,6), KAPPALAYERSO2(NF,6)
      COMMON/TFBLOK12/KAPPALAYERC2H6, KAPPALAYERSO2
C
c jfk 6/25/08 Temporary (PF should have been passed in the call
c     statement. Actually, this statement does remain once the code is
c     corrected.)
      DIMENSION PF(ND)

c      DIMENSION P(ND),BPLANCK(ND),CPR(NF),TPR(NF),TPRIND(ND), 
c     2 TAUGCO2(ND),TAUGH2O(ND),TAUGCH4(ND),TAUGIR(ND),FUP(ND),FDN(ND),
c     3 FUPA(ND),FDNA(ND),T(ND),CGAS(ND,NGS),TAUCONTIN(ND),WEIGHT(8,3),
c     & xkappa(8,12,55,8,3)
      DIMENSION P(ND),BPLANCK(ND),CPR(NF),TPR(NF),TPRIND(ND), 
     2 TAUGCO2(ND),TAUGH2O(ND),TAUGCH4(ND),TAUGIR(ND),FUP(ND),FDN(ND),
     3 FUPA(ND),FDNA(ND),T(ND),CGAS(ND,NGS),TAUCONTIN(ND)
      REAL KAPPALAYER(55,8,3,ND),kappa(55,8,3),KAPPALAYEROZ(8,ND),
     2 WEIGHTOZC(8),KAPPAOZC(8),TAUGOZ(ND),TAUCONTINT(NF),CNUT,FI(4,ND),
     3 TAUTOTAL(NF),TRANSLAYER(ND),TWGHTT(8),OMG0IR(ND-1),ASYIR(ND-1),
     4 TAULAMIR(ND-1),TAUAEXTIR(ND-1),TAUASIR(ND-1),TAUSIR(ND-1)
C
      DATA HP,SIGMA/6.63E-27, 5.67E-5/
      DIMENSION WEIGHTC2H6(6),WEIGHTSO2(6),CGASC2H6(ND), CGASSO2(ND)
      DIMENSION TAUGC2H6(ND), TAUGSO2(ND)
      REAL np

      REAL TAUSUM(6)


C   PRESSURE-INDUCED CO2 ABSORPTION (FROM JIM POLLACK)
      DATA CPR/4.3E-5, 3.8E-5, 1.2E-5, 2.8E-6, 7.6E-7, 4.5E-7, 2.3E-7,
     2  5.4E-7, 1.6E-6, 10*0., 4.E-7, 4.E-6, 1.4E-5, 1.0E-5,
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
c-jdh Ethane k-coefficient weights
      DATA WEIGHTC2H6/0.08566,0.18038,0.23396,0.23396,0.18038,0.08566/
      DATA WEIGHTSO2/0.08566,0.18038,0.23396,0.23396,0.18038,0.08566/
C
       INTEGER COUNTERIR, INCLUDE_CH4, INCLUDE_C2H6, INCLUDE_SO2
       REAL TAU_TF=0.
C
C-TF   THESE FLAGS ARE FOR SELF-CONSISTENCY CHECK.
       INCLUDE_CH4 = 1
       INCLUDE_C2H6 =0
       INCLUDE_SO2 = 0
       K1=1
       K2=1
       K3=1
       K4=0
       K5=0
C
C
       COUNTERIR = 0
       SRFALBIR = 0.
       NLAYERS = ND - 1    
c      BCON = 2.*HP/C/C
       HK = HP/BK
C       np = 1.E+1
C
C      DO I=1,6
C      PRINT 99931, KAPPALAYERSO2(9:11,I),
C     | KAPPALAYERSO2(17:23,I), KAPPALAYERSO2(29:32,I)
C      ENDDO
99931 FORMAT(1P14E10.2)
C      DO I=1,6
C      PRINT 99938, KAPPALAYERC2H6(14:16,I),
C     | KAPPALAYERC2H6(22:25,I)
C      ENDDO
99938 FORMAT(1P7E10.2)
C
C-TF  ALTHOUGH CGASC2H6 IS COMPUTED HERE. IT IS NOT USED
C-TF  TO CALCULATE THE TAU IN JACOB'S OROGINAL CODE
      DO I=1,NLAYERS
c        CGASC2H6(I) = 1.3386e-3*BK*273.16/(SM*G)
c     &                * (P(I+1) - P(I))*FC2H6/DM
c        CGASC2H6(I) = 1.0e-5*BK*273.16/(SM*G)
c     &                * (P(I+1) - P(I))*FC2H6/DM
C
c jfk 6/27/08 P was changed to PF in the 2 lines below.
        CGASC2H6(I) = 1.0e-5*BK*273.16/(SM*G)
     &                * (PF(I+1) - PF(I))*FC2H6_NEW(I)/DM
C
C-TF   COMPUTE SO2 COLUMN DEPTH
C       CGASSO2(I) = 1.0e-5*BK*273.16/(SM*G)
C     &                * (PF(I+1) - PF(I))*FSO2_NEW(I)/DM
C
C-TF   CGASSO2 UNIT IS CM-2 BECAUSE THE CROSS SECTION DATA 
C-TF   ARE IN UNIT OF CM2
       CGASSO2(I) = (PF(I+1) - PF(I))/
     |  (SM*G*DM)*FSO2_NEW(I)*1e6
C
      TAUGCH4(I) = 0.
      TAUGC2H6(I) = 0.
      TAUGSO2(I) = 0.
C
      END DO
C
C Read the IR exponential sums, WEIGHT AND XKAPPA
c-jdh Moved to Clima.f 
c      CALL IREXPSUMS(WEIGHT,xkappa)

      DO 7 IL = 1,NLAYERS
c      CALL INTERP(T(IL),P(IL),xkappa,kappa)
       CALL INTERP(T(IL),P(IL),kappa)
       DO 8 I=1, 55
       DO 9 J=1, 8
       DO 10 K=1, 3
       KAPPALAYER(I,J,K,IL) = kappa(I,J,K)
  10  CONTINUE
  9   CONTINUE
  8   CONTINUE
  7   CONTINUE
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
       DO 100 I=1, NF
C
       DO 20 J=1,ND 
C-jdh  VAC = C*AV(I)
       VAC = AV(I)
       BPLANCK(J) = PLANCK(VAC,T(J),HP,C,HK)
  20   CONTINUE
C
       DO 15 J=1,ND
         FUPA(J) = 0.0
         FDNA(J) = 0.0
         FUP(J) = 0.0
         FDN(J) = 0.0
         TRANSLAYER(J) = 0.0 
  15   CONTINUE
C
****** AEROSOLS
C
       DO IL=1,NLAYERS
        r = RAER(IL)
        np = PARTICLES(IL)
        TAUAEXTIR(IL) = QEXTIR(I,IL)*PI*r*r*DALT(IL)*np*1.E+5
        TAUASIR(IL) = OMG0AIR(I,IL)*TAUAEXTIR(IL)
       ENDDO
C
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
C
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
C
       ELSE
       DO IL = 1,NLAYERS
       TAUCONTIN(IL) = 0.
       ENDDO
C
       ENDIF
C
C******
C
c-jdh  If(I.lt.55) VAC1 = C*(AV(I+1)-AV(I))
c-jdh  If(I.eq.55) VAC1 = C*(AV(I) - AV(I-1))
c       IF(I.lt.NF) THEN
c         VAC1 = (AV(I+1)-AV(I))
c       ELSE
c         VAC1 = (AV(I)-AV(I-1))
c       ENDIF
C
C
C
C-TF   START K-COEFFICIENT LOOPS
C-TF   KAPPALAYER(I,K,1:3,IL) ARE SET UP IN IREXPSUMS.F
C-TF   FOR H2O, CO2, AND CH4.
C
       TWGHT = 1.
       DO 1 K1=1, 8    ! H2O
       DO 2 K2=1, 8    ! CO2
       DO 3 K3=1, 6    ! CH4
c       DO 4 K4=1, 6    ! C2H6
C       DO 5 K5=1, 6    ! SO2
C
C-TF   THERE IS THE SELF-CONSISTENCY CHECK.
       IF (   ( (INCLUDE_CH4 .EQ. 0) .AND. (K3 .NE. 0) ) 
     |   .OR. ( (INCLUDE_C2H6 .EQ. 0).AND. (K4 .NE. 0) )
     |   .OR. ( (INCLUDE_SO2 .EQ. 0) .AND. (K5 .NE. 0) )
     |   .OR. ( (INCLUDE_SO2 .EQ. 1) .AND. (K5 .EQ. 0) )
     |   .OR. ( (INCLUDE_C2H6 .EQ. 1).AND. (K4 .EQ. 0) )
     |   .OR. ( (INCLUDE_CH4 .EQ. 1) .AND. (K3 .EQ. 0) )
     |    )  THEN 
       PRINT*, "in ir.f: one gas is included", 
     |   " in k-eff loop but its flag is off. Please check!!!"
       STOP
       ENDIF
C
        TWGHT = WEIGHT(K1,1)*WEIGHT(K2,2)
C        
        IF (INCLUDE_CH4 .EQ. 1)  TWGHT = TWGHT * WEIGHT(K3,3)
        IF (INCLUDE_C2H6 .EQ. 1) TWGHT = TWGHT * WEIGHTC2H6(K4)
        IF (INCLUDE_SO2 .EQ. 1)  TWGHT = TWGHT * WEIGHTSO2(K5)
C
        DO 11 IL = 1,NLAYERS
         PPE = (1. + 0.3*FI(2,IL))*P(IL)     !CO2
         TPE = (300./T(IL))**TPR(I)
         TAUGCO2(IL) = KAPPALAYER(I,K1,1,IL)*CGAS(IL,3)/(44.*SM)
         TAUGH2O(IL) = KAPPALAYER(I,K2,2,IL)*CGAS(IL,2)/(18.*SM)
C
C-TF     THE FOLLOWING WAS USED IN MIKE&WILL's MODEL
        IF (INCLUDE_CH4 .EQ. 1) THEN
         TAUGCH4(IL) = KAPPALAYER(I,K3,3,IL)*CGAS(IL,5)*2.687E24
        ENDIF 
C
C-TF     THE FOLLOWING EXPRESSION PROBABLY USES WRONG UNITS
        IF (INCLUDE_C2H6 .EQ. 1) THEN
         TAUGC2H6(IL) = KAPPALAYERC2H6(I,K4)*CGASC2H6(IL)*2.687E24
        ENDIF
C
C-TF    SO2
        IF (INCLUDE_SO2 .EQ. 1) THEN
         TAUGSO2(IL) = KAPPALAYERSO2(I,K5)*CGASSO2(IL)
        ENDIF 
C
****** Pressure induced absorption by CO2 
         CGAS1 = CGAS(IL,3)/1.963E-3  
         TPRIND(IL) = CPR(I)*PPE*TPE*CGAS1
c        TPRIND(IL) = 0
*******
         TAUGIR(IL) = TAUGCO2(IL)+TAUGH2O(IL)+TPRIND(IL)
     |               +TAUCONTIN(IL)
C
C
       IF (INCLUDE_CH4 .EQ. 1) TAUGIR(IL) = TAUGIR(IL) + TAUGCH4(IL)
       IF (INCLUDE_C2H6 .EQ. 1) TAUGIR(IL) = TAUGIR(IL) + TAUGC2H6(IL)
       IF (INCLUDE_SO2 .EQ. 1)  TAUGIR(IL) = TAUGIR(IL) + TAUGSO2(IL)
C
  11    CONTINUE
C
C        IF ( (K1 .EQ. 1) .AND. (K2 .EQ. 1) .AND. (K5 .EQ. 1) .AND.
C     |       (I .EQ. 9) ) THEN
C        PRINT*, KAPPALAYERSO2(9,1), TAU_TF
C        ENDIF
C
******* Ozone absorption
C-TF  THIS SPECIAL TREATMENT FOR OZONE IS DONE TO SAVE COMPUTER
C-TF  TIME BECAUSE ONLY ONE WAVELENGTH IS AFFECTED BY OZONE.
C-TF  HERE A NEW WEIGHT FUNCTION TWGHTT(K00) IS CALCULATED
C-TF  AND IS USED TO COMPUTE FUPA AND FDNA.
C
        IF (I.EQ.18) THEN
      DO K00=1,8
        TWGHTT(K00) = TWGHT*WEIGHTOZC(K00)
C
      DO IL=1,NLAYERS
        TAUGOZ(IL) =KAPPALAYEROZ(K00,IL)*CGAS(IL,4)*2.687E19
        TAUGIR(IL) =TAUGIR(IL) + TAUGOZ(IL)
      ENDDO
C
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
                  FUPA(J)=FUPA(J)+TWGHTT(K00)*FUP(J)
                  FDNA(J)=FDNA(J)+TWGHTT(K00)*FDN(J)
                  J=J+1
                  GOTO 13
               END IF
        ENDDO  ! END K00 LOOP
C
C-TF  JUMP OVER THE FOLLOWING CALCULATIONS BECAUSE THEY ARE ALREADY
C-TF  DONE HERE FOR OZONE.
        GOTO 1001
C        
        ENDIF  ! END IF (I.EQ.18)
C
*************** 
C-TF  IF THE WAVELENGTH IS NOT WHERE OZONE IS IMPORTANT
C-TF  SIMILAR CALCULATIONS NEED TO BE DONE BUT USING THE
C-TF  WEIGHT CALCULATIED RIGHT AFTER THE K[23456] LOOPS
C
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
1001    CONTINUE     ! JUMP POINT FOR THE OZONE CALCULATIONS
C
c  5     CONTINUE     ! END SO2 LOOP
c  4     CONTINUE     ! END C2H6 LOOP
  3     CONTINUE     ! END CH4 LOOP  
  2     CONTINUE     ! END CO2 LOOP       
  1     CONTINUE     ! END H2O LOOP
C
         DO 14 J=1,ND
             FDNIR(J)=FDNIR(J)+W(I)*FDNA(J)
             FUPIR(J)=FUPIR(J)+W(I)*FUPA(J)
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
100     CONTINUE                   !**END LOOP over frequency**
C
19     FORMAT(/1X,1PE10.4)
1005   FORMAT(1X,1P10E12.5) 
C      PRINT*,'TAUTOTAL'
C      PRINT 1005, TAUTOTAL
C
      RETURN
      END
