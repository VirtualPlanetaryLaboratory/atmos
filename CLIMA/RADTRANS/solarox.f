      SUBROUTINE SOLAROX(TF,LAST,FNC,NST)
C  Add NO2 to the mix of gases
C
C  FORMERLY TOM ACKERMANS SOLAR ROUTINE 'RADRIV'. ALL SUBROUTINES UNDER
C  'RADRIV' HAVE BEEN ELIMINATED EXCEPT FOR 'RAYLEY' (NOW ALSO CALLED
C  'RAYLEY'), 'SOLD2S' (NOW, THIS SUBROUTINE CALLED 'SOLAR'), AND
C  'SOLAR' (NOW PART OF KASTING'S MAIN CODE BEFORE THE IR-LOOP).
C  THE TWO-STREAM 'SOLVER' SUBROUTINES HAVE BEEN REPLACED BY A HARDIER
C  ROUTINE; ONE THAT DOESN'T BOMB FOR SMALL ZENITH ANGLES, CALLED 'DELTA2STR.'
C  ASSUME UNITS TO BE CGS UNLESS STATED OTHERWISE.   DMW (7-31-95)
C PF1


c  This subroutine contains CH4

      INCLUDE 'CLIMA/INCLUDE/header.inc'
      PARAMETER (NSOL=38,NGS=8,NF=55, IK=8)  !c-rr Adding IK=8, number of sums for CO2 and H2O absorption coefficients 8/27/2012
      PARAMETER(NS=3, NS1=NS+2, NS4=NS+5) ! Adding parameter statement needed for FI(NS1,ND) 5/23/2011
      REAL kmatrix_solh2o, kmatrix_solco2, weights, KAPPALAYERSOL_CO2,
     &   KAPPALAYERSOL_H2O ! EWS - BETA not used
C

       COMMON/SOLARBLK/AMU0,SRFALB,OMG0A(NSOL,ND-1),
     &  ASYA(NSOL,ND-1),TAUAER(NSOL),SIGERT(NSOL),FMA(NSOL),PF(ND),
     &  ALAMBDA(NSOL),CGAS(ND,NGS),FUPSOL(ND),FDNSOL(ND),
     &  NGAS(2,NSOL),WGHT(4,2,NSOL),NPR(2,NSOL),SOLINT(NSOL),
     &  TAULAM(ND-1),ASY(ND-1),OMG0(ND-1),FMT(ND-1),QEXT(NSOL,ND-1)

        COMMON/DATASOLAR/weightco2_h2oSOL(8), weights(3,NSOL,IK)   ! new common block for weights and interpolated coefficients for CO2, H2O, and methane. 8/26/2012



      COMMON/ALTBLOK/DALT(ND-1),RADIUS(ND-1),PARTICLES(ND),RAER(ND),
     2 ALT(ND)
       COMMON/PRESS/BETIR1(4,5,NSOL),BETIR2(4,5,NSOL),
     &  kappa_solh2o(NSOL,8,8,IK), kappa_solco2(NSOL,8,8,IK) ! Added new kappa matricies for each of CO2 and H2O coefficients. 8/26/2012

      COMMON/AOZONE/BETAO3(nsol), BETAO2(2),
     &   WGHTO2(NSOL,2)
!     Added O3 absorption coefficient vector(plus c and  d vectors), first 14 terms. 3/19/2012 ! Added new kappa matrix for mixed CO2-H2O coefficients
c     COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4
      COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4,FC2H6,FNO2,FI(NS1,ND),
     &   FH22
      COMMON/PRESSURE/P(ND),PLOG(ND)
      COMMON/CONSS/C,BK,G,GNEW(ND),PI,SM,DM,DM2
      COMMON/CH4BLOCK/ALPHACH4T188(4,17),BETACH4T188(4,17),
     & ALPHACH4T295(4,17),BETACH4T295(4,17),ALPHACH4Kark(4,21),
     & BETACH4Kark(4,21),GAMMAEXP188(17),GAMMAEXP295(17),
     & ALPHACH4NEW(6),BETACH4NEW(17,3,5,6),ALCH4(6,38)
      COMMON/WAVE/AV(NF),LAM(NF)
      COMMON/NO2BLOK/SIGNO2(NSOL)





!       COMMON/SOLARDATA/weightco2_h2oSOL(16), weights(2,NSOL,IK)
!     & ,kmatrix_sol(NSOL,IK) ! c-rr Added weightco2_h2O data array and weights array that combines co2_h2o and Ch4 weight arrays 3/19/2012! c-rr Added 16 weight coefficient array for CO2-H2O kspectrum mixed absorption coefficiets, 3/15/2012



      COMMON/CPHEAT/CPO2(ND),CPCO2(ND), CPN2(ND), CPH2O(ND),
     & CPN(ND), CPNT(ND)

      DIMENSION  TAUAS(ND-1),TAUS(ND-1), ! EWS - TAUA(ND-1) not used
     &  TAUR(ND-1),TAUG(ND-1),FRAC(ND-1),INDEX(ND-1),
     &  CRAY(ND-1),DPLYR(ND-1),PLAYR(ND),PMID(ND-1),
     &  FUP(ND),FDN(ND),FUPA(ND),FDNA(ND),ALPHACH4(4,38),
     &  BETACH4(4,38),TAUEXT(ND-1),TAUEXTID(ND-1),GAMMACH4(38)
      DIMENSION TF(ND),Fdiff(ND),SolHeat(NSOL,ND), PF1(ND)
      DIMENSION FNC(ND), KMATRIX_SOLH2O(NSOL,IK),
     &  KMATRIX_SOLCO2(NSOL,IK),self_absol(NSOL,ND),
     &   forn_absol(NSOL,ND),  ! Added self and foreign water continuum solar coefficients 8/30/2012 c-rr
     &  SOLAV(NSOL+1) ! Added solar intervals
      DIMENSION ALAM(NF),TAUCONTSOL(ND),  ! Added water continuum tau for solar 8/30/2012 ! W(NF), AVOLD(NF) not used - EWS 9/4/2015
     &  KAPPALAYERSOL_CO2(NSOL,ND,IK),KAPPALAYERSOL_H2O(NSOL,ND,IK) ! Turning SIGR into a scalar quantity not an array 5/23/2011
!     Added KAPPALAYERSOL which is KMATRIX_SOL looped over all heights
c     DIMENSION W(NF), ALAM(NF), AVOLD(NF),SIGR(NSOL)
C     DIMENSION ALAM(NF), AVOLD(NF),SIGR(NSOL)
C

       COMMON/BPS_SOL/s_absol(NSOL), f_absol(NSOL), TDsol(NSOL),  ! Added COMMON BLOCK FOR BPS CONTINUUM FOR SOLAR  8/30/2012 c-rr
     &  Bssol(NSOL), Bfsol(NSOL)
      REAL np, SIGR,ALCH4  ! removed extraneous kappa and kmatrix 3/23/2012

       DATA SOLAV/42087 ,36363 ,35087 ,32562 ,30376 ,29308 ,25641 ,
     & 22222 ,18518 ,18198 ,17649 ,16528 ,16000 ,15000 ,14470 ,13300 ,
     & 12790 ,11870 ,11220 ,10400 , 9650 , 9350 , 8850 , 8315 , 7650 ,
     & 6990 , 6390 , 5925 , 5370 , 4950 , 4540 , 4030 , 3760 , 3425 ,
     & 3087 , 2796 , 2494 , 2397 , 2200/ ! Solar intervals 9 - 38 have BPS water continuum. Converted to wavenumbers 8/31/2012

       C2 = 1.4388 ! constant for BPS continuum (in cmK) 8/31/2012 , c-rr
C


C  Define the NO2 mixing ratio
      COUNTERS = 0
      NLAYERS = ND - 1
      DO 1140 JS=1,ND
         PLAYR(JS) = PF(JS)*1.E6      !PF IN BARS; PLAYR IN DYNES/CM^2
 1140 CONTINUE
C



      DO 1145 IL = 1,NLAYERS
         PMID(IL) = 0.5*(PLAYR(IL+1)+PLAYR(IL))

 1145 CONTINUE

C


!      CALL RAYLEY(SIGR)   Commenting out this RAYLEY CALL AND MOVING IT DOWN INSIDE ALTITUDE LOOP. 5/23/2011

C
C      WEIGHT = 0.
      DO 1150 IL=1,NLAYERS
         DPLYR(IL)=PLAYR(IL+1)-PLAYR(IL)
C        AM = 44.*FCO2(IL) + DM*(1.-FCO2(IL))
c         AM = 44.*FCO2 + DM*FN2
          AM = DM

        CONS0=6.0255E23/GNEW(IL)
         CONS=CONS0/AM
       CRAY(IL)=CONS*DPLYR(IL)
         PL=LOG10(PMID(IL)) - 1.
         PL=AMAX1(PL,2.001)    !10^-3 BARS < P < 10 BARS =>
         PL=AMIN1(PL,5.999)    !10^3 DYNES/CM^2 < P < 10^7 DYNES/CM^2
         LP=PL
         INDEX(IL)=LP-1
         FRAC(IL)=PL-LP
 1150 CONTINUE
C

      DO 1152 I=1,ND
         FUPSOL(I) = 0.0
         FDNSOL(I) = 0.0
         FUP(I) = 0.0
         FDN(I) = 0.0
 1152 CONTINUE
      T =188
C SELECTION OF CH4 DATA USING THE CURRENT VALUE OF TEMPERATURE
      TCOMPARISON = (188+295)/2
      IF (T.LE.TCOMPARISON) THEN
C3340 FORMAT('Using data for 188 temperature')
C-rr Absorption coefficient units are in inverse ATM-KM for CH4 (according to gascon), so converting to cm^2/molecule by dividing by 2.687E24. 3/19/2012
       DO 3334 J=1,17
        GAMMACH4(39-J) = GAMMAEXP188(J)
       DO 3333 K3=1,4
        ALPHACH4(K3,39-J) = ALPHACH4T188(K3,J)
        BETACH4(K3,39-J) = BETACH4T188(K3,J)/2.687E24
 3333   CONTINUE
 3334  CONTINUE
       ELSE
C3341 FORMAT('Using data for 255 temperature')
       DO 3336 J=1,17
        GAMMACH4(39-J) = GAMMAEXP295(J)
       DO 3335 K3=1,4
        ALPHACH4(K3,39-J) = ALPHACH4T295(K3,J)
        BETACH4(K3,39-J) = BETACH4T295(K3,J)/2.687E24
 3335   CONTINUE
 3336  CONTINUE
      END IF
C PLUGGING IN CH4 DATA FROM KARKOSHKAs file
       DO 3338 J=1,21
        GAMMACH4(J) = 0.
       DO 3337 K3=1,4
       ALPHACH4(K3,J) = ALPHACH4Kark(K3,J)
       BETACH4(K3,J) = BETACH4Kark(K3,J)/2.687E24
 3337  CONTINUE
 3338 CONTINUE
        DO J = 1,ND
           PF1(J) = PF(J)*1.E6
        END DO
C

! c-rr modified interp routine for mixed kspectrum CO2-H2O absorption coefficients. Also, this loop feeds in altitude (j) for the interpolated kmatrix, 3/15/2012

!           print *,'GOING IN INTERPSOLAR....'

       DO j = 1,ND
             CALL interpsolar(TF(j),P(j),j,
     &      kmatrix_solco2, kmatrix_solh2o)
        DO K =1,IK

          DO I = 1, NSOL
          kappalayersol_co2(I,J,K) = kmatrix_solco2(I,K)  ! Interpolated CO2 coefficients stored in kappalayersolco2

!          if(I <=8)then

!             kappalayersol_h2o(I,J,K) = kappalayersol_h2o(1,J,K)  ! Interpolated H2O coefficients stored in kappalayersolh2o
!          else
          kappalayersol_h2o(I,J,K) = kmatrix_solh2o(I,K)  ! Interpolated H2O coefficients stored in kappalayersolh2o
!          endif
!          print *,kmatrix_solh2o(I,K),kappalayersol_h2o(I,J,K),I,K,j
   !       IF(I .EQ. 10)then
   !       print 3456,KAPPALAYERSOL(I,J,K),TF(J),P(j),FI(1,J),I,J,K
   !       pause
   !       endif
          ENDDO

       ENDDO
         ENDDO

c 3456     format(1p4e14.5,3(2x,i3)) !EWS - label not used

!         print *,'OUT OF INTERPSOLAR....'


!         print 2233, (ALAMBDA(I), i=1,NSOL)
      DO 1155 I=1,NSOL     !**BEGIN WAVELENGTH LOOP**

         DO 1160 J=1,ND
            FUPA(J) = 0.0
            FDNA(J) = 0.0
 1160    CONTINUE

             AL2 = ALAMBDA(I)**2  ! Took AL2 out of Rayley in order to do altitude loop and get AL2 for each individual mixing ratio.

!2233              format(6(f9.5))

! 5/23/2011
       DO 1165 IL=1,NLAYERS
c-rr Recalculation of mixing rations for the noncondensibles
c FI(1,J)= water
c FI(2,J)= carbon dioxide
c FI(3,J) = methane


c The condensibles are water and carbon dioxide. Water convects for planets closer in and CO2 condenses for planets further out
c FNC = noncondensible mixing ratio = FN2 + FO2 + FCH4 + FAR = 1 - FH20 - FCO2  5/3/2011
           Fwater = FI(1,IL) ! Needed for rayley
           FCO2 = FI(2,IL)  ! Needed for rayley
           FNCR = FNC(IL)  ! This is the FNC that is passed into rayley, one layer at a time  c-rr 6/7/2012

c-----------------------------------------------------------------------------------------
           CALL RAYLEY(SIGR, AL2,Fwater,FNCR) ! Call rayley to output SIGR at a given altitude and wavelength, inputting AL2. 5/23/2011
           ! 8/31/2011 removed FNC from rayley argument

c           TAUR(IL)= SIGR(I)*CRAY(IL)
        TAUR(IL)= SIGR*CRAY(IL)  ! SIGR is a scalar now 5/23/2011
        r = RAER(IL)
        np = PARTICLES(IL)
        TAUEXT(IL)=QEXT(I,IL)*PI*r*r*DALT(IL)*np*1.E+5
        TAUEXTID(IL) = 2*PI*r*r*DALT(IL)*np*1.E+5
        TAUAS(IL)=OMG0A(I,IL)*TAUEXT(IL)
        TAUS(IL)=TAUAS(IL)+TAUR(IL)
 1165  CONTINUE
C

       IF (I.EQ.9) THEN
       TAUEXTTOTAL = 0.
       TAUASTOTAL = 0.
       DO IL=1,NLAYERS
       TAUEXTTOTAL = TAUEXTTOTAL + TAUEXT(IL)
       TAUASTOTAL = TAUASTOTAL + TAUAS(IL)
       ENDDO
c       PRINT *, '***************************'
c       PRINT *, 'TAUEXTTOTALVIS'
c       PRINT *, TAUEXTTOTAL
c       PRINT *, 'TAUASTOTALVIS'
c       PRINT *, TAUASTOTAL
       TAUAABSTOTAL = TAUEXTTOTAL - TAUASTOTAL
c       PRINT *, 'TAUAABSTOTALVIS'
c       PRINT *, TAUAABSTOTAL
C       PRINT *, 'QEXTVIS      OMG0A      ALT'
C       DO IL=1,NLAYERS
C       PRINT *, QEXT(I,IL), OMG0A(I,IL), ALT(IL)
C       ENDDO
C      PRINT *, '***************************'
       ENDIF
C

       IF (I.EQ.1) THEN
       TAUEXTTOTAL = 0.
       TAUASTOTAL = 0.
       DO IL=1,NLAYERS
       TAUEXTTOTAL = TAUEXTTOTAL + TAUEXT(IL)
       TAUASTOTAL = TAUASTOTAL + TAUAS(IL)
       ENDDO
       ENDIF
c   19   FORMAT(/1X,1PE10.4) ! EWS - label not used
         NPR1=NPR(1,I)
         NPR2=NPR(2,I)
         NPR3=4  ! 4 terms for CH4 for intervals 1-21
         IF(I>21) THEN
         NPR3=6  ! 6 terms for Ch4 for intervals 22-38
         END IF



!----------------BPS WATER CONTINUUM 8/30/2012 c-rr
                 IF(I.le.8)THEN
                 DO IL = 1,NLAYERS
                 TAUCONTSOL(IL) = 0.0d0
                 ENDDO
                 ENDIF

                IF (I.ge.9) THEN  ! water continuum only in solar intervals 9 - 38
!                IF (((I.ge.30).and.(I.le.31)).or.((I.ge.34).and.(I.le.38)))THEN
!                IF ((I.eq.30).or.(I.eq.31))THEN

                DO IL = 1, NLAYERS

                PH2O = P(IL)*FI(1,IL)  ! water partial pressure (in bars)
                PDRY = P(IL)*(1.-FI(1,IL)) !partial pressure dry air (in bars)
                RHOW = (PH2O/1.)*(296./TF(IL))  ! RHO_water. Different from MT_CKD? Why does BPS use Ph2o instead of P?
                RHOF = (PDRY/1.)*(296./TF(IL))  ! RHO_foreign. Different from MT_CKD? Why does BPS use Pdry instead of P?

!                   Bssol(i) = 0.0d0
!                   Bfsol(i) = 0.0d0
                RADFLD = 1.  ! radiation field term is always negligible in solar because wavenumber always greater than 500cm-1.
                self_absol(I,IL) = RHOW*(s_absol(i)*
     &                            exp(TDsol(i)*(296.-TF(IL)))
     &         + RADFLD*Bssol(i)) ! Self broadening coefficient with radiation field and temperature dependence included (cm^2/molecule)
                forn_absol(I,IL) = RHOF*(f_absol(i) + Bfsol(i))*RADFLD
                  ! foreign broadening coefficient with radiation field and temperature dependence included (cm^2/molecule)

                ABSCONT = self_absol(I,IL) + forn_absol(I,IL)
                TAUCONTSOL(IL) = ABSCONT*CGAS(IL,6) !optical depth
!                TAUCONTSOL(IL) = 0.0d0
              ENDDO  ! ENDS LAYER LOOP IN CONTINUUM

              ENDIF! water continuum only in solar intervals 9 - 38

!------------------------------------------------------------------




C-rr  This section of code will do intervals 1-16 using Tom Ackerman's CO2 coefficients (WGHT arrays)
C     and 17-38 using Richard Freedman's (weight array). ALPHACH4(K3,I) for Intervals 1-21 and
C     ALPHACH4NEW(K3) for Intervals 22-38 have been combined into ALCH4(K3,I). 10/23/2010



C     ALPHACH4(K3,I) for Intervals 1-21 and
C     ALPHACH4NEW(K3) for Intervals 22-38 have been combined into ALCH4(K3,I). 10/23/2010


!        DO 1180 K3=1,NPR3              !**BEGIN LOOP by Ch4 exp sums' coef
!            IF (I.ge.17) THEN
!               NPR2=8
!            ENDIF
!         DO 1170 K1=1,NPR1              !**BEGIN NPR1 LOOP** K1 is gas 1 from solar38
!            DO 1175 K2=1,NPR2           !**BEGIN NPR2 LOOP** K2 is gas 2 from solar38
!               IF(I.lt.17) THEN
!               AP=WGHT(K1,1,I)*WGHT(K2,2,I)*ALCH4(K3,I)
!               ELSE
!               AP = WGHT(K1,1,I)*weight(K2)*ALCH4(K3,I)
!               END IF

               SUMAP = 0.
       DO 1190 KO2 = 1,2 ! For 2 oxygen weighting terms
!       DO 1180 K3 = 1, NPR3   !**BEGIN NPR3 LOOP** K3 is gas 2 (CH4) from nearIR_expsums.pdat
                DO 1175 K2 = 1, IK  ! K2 is kspectrum CO2 coefficients from solar38_h2o.dat
                     DO 1170 K1 = 1, IK  ! K1 is kspectrum H2O coefficients from solar38_h2o.dat

                   AP = weights(1, I, K1)*weights(2,I,K2)*
     &             WGHTO2(I,KO2)  ! Only convolve with O2



!               IG1=NGAS(1,I)  ! IG1 is mixed CO2/H2O and IG2 is CH4.
!               IG2=NGAS(2,I)
C
C     CONTRIBUTION OF CH4 absorption coefficients to the gas TAU
         DO 1181 IL= 1,NLAYERS
c            TAUG(IL)=0.
c  Add the NO2 absorption in here. Scale everything to CH4.


             CGAS(IL,2) = AMAX1(CGAS(IL,2),1.E-20)

             TAUG(IL) = SIGNO2(I)*CGAS(IL,2)*FNO2/FCH4 * 2.687E24  ! c-rr Change CGAS CH4 number from 5 to 2. 3/19/2012

c         IF (I .eq. 7) THEN
c         PRINT 99923, IL, FNO2, TAUG(IL), SIGNO2(I)
c         ENDIF
c 99923 FORMAT("in solar.f:",I5, 1P10E10.2) !EWS - label not used

            IF(I<22) THEN
            TAUG(IL)= TAUG(IL) + BETACH4(K3,I)*CGAS(IL,2) ! c-rr change gas CH4 number from 5 to 2. 3/19/2012
            END IF

         IF(I>21) THEN
         II = I - 21
C*********************ADDING STUFF HERE*****************************
C        Interpolation scheme to select the correct K-coefficients for BETACH4
C        Pressures are -4,-3,-2,-1,0
C        Temps are 112,188,295
C        Assumes Temp will never be lower than 112, log10(pressure) never lower than -4

      AFP = 0.0
      AFT = 0.0
      TP = AMAX1(TF(IL),112.)
      TP = AMIN1(TP,295.)
      PPLOG = AMIN1(PLOG(IL),0.)
      PPLOG = AMAX1(PPLOG,-4.)
      PP = AMAX1(P(IL),1.E-04)
      PP = AMIN1(PP,1.)
      ANEWBETA = 0.0


      IF(TP-112. < 76.) THEN
        AFT = (TP-112.)/76.
        LOGP = PPLOG-1
        SELECT CASE(LOGP)
        CASE(:-4)
        AFP = (PP-0.0001)/0.0009
        ANEWBETA= AFT*AFP*BETACH4NEW(II,2,2,K3)+(1-AFP)*(1-AFT)*
     1    BETACH4NEW(II,1,1,K3)+AFT*(1-AFP)*BETACH4NEW(II,2,1,K3)+
     2    AFP*(1-AFT)*BETACH4NEW(II,1,2,K3)

        CASE(-3)
        AFP = (PP-0.001)/0.009
        ANEWBETA= AFT*AFP*BETACH4NEW(II,2,3,K3)+(1-AFP)*(1-AFT)*
     1    BETACH4NEW(II,1,2,K3)+AFT*(1-AFP)*BETACH4NEW(II,2,2,K3)+
     2    AFP*(1-AFT)*BETACH4NEW(II,1,3,K3)

        CASE(-2)
        AFP = (PP-0.01)/0.09
        ANEWBETA= AFT*AFP*BETACH4NEW(II,2,4,K3)+(1-AFP)*(1-AFT)*
     1    BETACH4NEW(II,1,3,K3)+AFT*(1-AFP)*BETACH4NEW(II,2,3,K3)+
     2    AFP*(1-AFT)*BETACH4NEW(II,1,4,K3)

        CASE(-1:)
        AFP = (PP-0.1)/0.9
        ANEWBETA= AFT*AFP*BETACH4NEW(II,2,5,K3)+(1-AFP)*(1-AFT)*
     1    BETACH4NEW(II,1,4,K3)+AFT*(1-AFP)*BETACH4NEW(II,2,4,K3)+
     2    AFP*(1-AFT)*BETACH4NEW(II,1,5,K3)
        END SELECT

      ELSE IF(TP-188. <= 107.) THEN
        AFT = (TP-188.)/107.
        LOGP = PPLOG-1
        SELECT CASE(LOGP)
        CASE(:-4)
        AFP = (PP-0.0001)/0.0009
        ANEWBETA= AFT*AFP*BETACH4NEW(II,3,2,K3)+(1-AFP)*(1-AFT)*
     1    BETACH4NEW(II,2,1,K3)+AFT*(1-AFP)*BETACH4NEW(II,3,1,K3)+
     2    AFP*(1-AFT)*BETACH4NEW(II,2,2,K3)

        CASE(-3)
        AFP = (PP-0.001)/0.009
        ANEWBETA= AFT*AFP*BETACH4NEW(II,3,3,K3)+(1-AFP)*(1-AFT)*
     1    BETACH4NEW(II,2,2,K3)+AFT*(1-AFP)*BETACH4NEW(II,3,2,K3)+
     2    AFP*(1-AFT)*BETACH4NEW(II,2,3,K3)

        CASE(-2)
        AFP = (PP-0.01)/0.09
        ANEWBETA= AFT*AFP*BETACH4NEW(II,3,4,K3)+(1-AFP)*(1-AFT)*
     1    BETACH4NEW(II,2,3,K3)+AFT*(1-AFP)*BETACH4NEW(II,3,3,K3)+
     2    AFP*(1-AFT)*BETACH4NEW(II,2,4,K3)

        CASE(-1:)
        AFP = (PP-0.1)/0.9
        ANEWBETA= AFT*AFP*BETACH4NEW(II,3,5,K3)+(1-AFP)*(1-AFT)*
     1    BETACH4NEW(II,2,4,K3)+AFT*(1-AFP)*BETACH4NEW(II,3,4,K3)+
     2    AFP*(1-AFT)*BETACH4NEW(II,2,5,K3)
        END SELECT
      END IF
            TAUG(IL)= TAUG(IL) + ANEWBETA*CGAS(IL,2)  ! Change CH4 from 5 to 2
c            print 9999,TAUG(IL),ANEWBETA
c            print 9998,TP,LOGP
c           print 9999,AFT,pp
c 9998        format(1PE12.5,2x,I2) !EWS - label not used
      END IF

 1181    CONTINUE

!               IF (IG1.LE.4) THEN
C-rr           10/19/2010. Gas "3" is CO2 added
C               from Richard Freedmans eb.txt.
C              For intervals 1-16 use Tom Ackerman's kcoefficients to get BETAs
C              otherwise use Richard Freedman's.
!                  DO 1185 IL=1,NLAYERS
!                     L=INDEX(IL)
!                     FR=FRAC(IL)
!                     BETA=FR*BETIR1(K1,L+1,I)+(1.-FR)*BETIR1(K1,L,I) ! BETIR1 ozone interpolation
!                     TAUG(IL)=TAUG(IL)+BETA*CGAS(IL,IG1)
! 1185             CONTINUE
!               END IF
C


!  TAU expression simplified so that all species are done within a single loop 3/19/2012
                  DO 1182 IL = 1, NLAYERS
                       BETACO2 = KAPPALAYERSOL_CO2(I,IL,K2) ! kappas for CO2 coefficients
                       BETAH2O = KAPPALAYERSOL_H2O(I,IL,K1) ! kappas for H2O coefficients
!                      if(I .eq. 9)then
!                      print 9999, KAPPALAYERSOL_CO2(I,IL,K1)
!     &                         ,KAPPALAYERSOL_H2O(I,IL,K1),
!     &                          I,IL,K1, P(IL),TF(IL)
!                       pause
!                       endif

!                       IF (I.ge.9)THEN
!                       print *, BETACO2, BETAH2O, I, IL
!                       pause
!                       ENDIF

                       IF (I.le.14) THEN


                       CGAS(IL,3) = AMAX1(CGAS(IL,3),1.E-20)
                       CGAS(IL,4) = AMAX1(CGAS(IL,4),1.E-20)
                       TAUG(IL)= TAUG(IL) + BETAH2O*CGAS(IL,6)
     &                 + BETAO3(I)*CGAS(IL,4) ! FOR CGAS: Gas 3 is O2, Gas 4 is O3, Gas 5 is CO2, and Gas 6 is H2O. Gas 2 (CH4) is handled above)
!                         if (IL.eq.1) then
!                         print *, BETAO3(I)*2.687E19, I, CGAS(IL,1)
!                         pause
!                         endif
                       ELSEIF(I.eq.15) THEN
                       TAUG(IL)= TAUG(IL) + BETACO2*CGAS(IL,5) +
     &                 BETAH2O*CGAS(IL,6)
                       ELSEIF(I.eq.16)THEN
                       TAUG(IL) = TAUG(IL) + BETACO2*CGAS(IL,5) +
     &                 BETAH2O*CGAS(IL,6) +
     &                 BETAO2(KO2)*CGAS(IL,3)
!                          print *, BETAO2, IL, I
!                          pause
                       ELSE ! For intervals 17-38
                       TAUG(IL)= TAUG(IL) + BETACO2*CGAS(IL,5) +
     &                 BETAH2O*CGAS(IL,6)
                       ENDIF

!                       if((I.gt.17).and.(IL.eq.50))then
!                        print *, 'CGASH2O=', CGAS(IL,6),
!     &                  'CGASCO2=', CGAS(IL,5), IL
!                        pause
!                       print *, 'TAUG=', TAUG(IL)
!                       pause
!                       endif
 1182             CONTINUE

!               IF (IG2.LE.4) THEN
!                  DO 1190 IL=1,NLAYERS

!                    IF (I<17) THEN
!                      L = INDEX(IL)
!                      FR = FRAC(IL)
!                      BETA=FR*BETIR2(K2,L+1,I)+(1.-FR)*BETIR2(K2,L,I)  ! This interpolation for the first 16 intervals are done in solar not in a separate interpolation routine. KMATRIX, in contrast, is already interpolated. We probably don't need this because BETIR2 is water here and we have the new coefficients that already have these values for intervals 1-38.
!                      TAUG(IL)=TAUG(IL)+BETA*CGAS(IL,IG2)  ! CGAS is in units of what units for gas(from gascon) and BETA is inverse of that. We need to change them to molecules/CM2 and CM2/molcule, respectively.
!                     ELSE
!                       BETA = KMATRIX(I,IL,K1)
!                       TAUG(IL)= TAUG(IL) + BETA*CGAS(IL,3)


!                          ENDIF

c-rr 6/7/2011 near IR CO2 CIA section------------------------------------------------------------------------

!                   IF(IG2.EQ.3) THEN
                     IF  ((I .EQ. 22).OR.(I.EQ.23).OR.(I.EQ.24))THEN  ! Put CO2 CIA in 1.2 micron complex
                     DO IL =1,NLAYERS
                     TAUG(IL) = TAUG(IL) + 1.5E-9*(CGAS(IL,5)
     &                             /2.687E19)*FI(2,IL)*P(IL) ! Lochsmidt's value converts CGAS for CO2 from mol/cm^2 into atm-cm. The CIA values (i.e. 1.5E-9 are
!                                                                      in inverse amagat^2cm. So tau is unitless.
                     ENDDO


                     ELSEIF(I.EQ.28) THEN ! Put CO2 CIA in 1.73 micron band
                      DO IL = 1, NLAYERS
                        TAUG(IL) = TAUG(IL) + 0.57*6.0E-9*(CGAS(IL,5)
     &                            /2.687E19)*FI(2,IL)*P(IL)
                       ENDDO

                     ELSEIF(I.EQ.31) THEN ! Put CO2 CIA in 2.3 micron band
                     DO IL = 1, NLAYERS
                      TAUG(IL) = TAUG(IL) + 3.5E-8*(CGAS(IL,5)
     &                            /2.687E19)*FI(2,IL)*P(IL)
                     ENDDO
                     ENDIF
!                     ENDIF
C------------------------------------------------------------------------------------------


               DO 1200 IL=1,NLAYERS
!                TAUG(IL) = AMAX1(TAUG(IL), 1.e-17)

               TAULAM(IL) = TAUEXT(IL)+TAUR(IL)+TAUG(IL)
     &              + TAUCONTSOL(IL) ! Added tau for BPS continuum 8/28/2012

!              if(TF(ND) > 2200.0d0)then
!              print 9999, TAULAM(IL), TAUG(IL),CGAS(IL,2),CGAS(IL,3),
!     &                    CGAS(IL,4), CGAS(IL,5),IL,K1,K2,I
!              endif




c 9999   format(1P6E12.5,3(2x,I3), 2x,1p3e12.5) !EWS - label not used
c               TAULAM(IL) = AMIN1(TAULAM(IL),1000.)
               OMG0(IL) = TAUS(IL)/TAULAM(IL)
C
C  Do not let scattering albedo be larger than 0.99999.
                  OMG0(IL) = AMIN1(OMG0(IL),0.99999)
C

                  TSRAT = TAUAS(IL)/TAUS(IL)
                  ASY(IL) = TSRAT*ASYA(I,IL)
C-AP                  FMT(IL) = TSRAT*FMA(I)
                  OMG0(IL) = AMAX1(OMG0(IL),1.E-5)! changed lower limit of single-scattering albedo to 1e-12 (default 1e-5) c-rr 4/30/2012
 1200          CONTINUE
               COUNTERS = COUNTERS + 1
               CALL DELTA2STR(SRFALB,AMU0,ASY,TAULAM,OMG0,FUP,FDN)
C

               SUMAP = SUMAP + AP

               J = 1
 1205          IF (J .LE. ND) THEN
                  FUPA(J)=FUPA(J)+AP*FUP(J)
                  FDNA(J)=FDNA(J)+AP*FDN(J)

!                  print *, 'AP=',AP, 'FDN=',FDN(J),'FUP=', FUP(J)
!                  pause

                  J=J+1
                  GOTO 1205
               END IF
C         WEIGHT = WEIGHT + AP
C

 1170    CONTINUE             !**END NPR1 LOOP ! H2O
 1175    CONTINUE                !**END NPR2 LOOP ! CO2
! 1180    CONTINUE                !**END NPR3 LOOP ! CH4
 1190    CONTINUE                ! ENDS KO2 loop
!            print *, 'SUMAP=', SUMAP

         DO 1210 J=1,ND
             FDNSOL(J)=FDNSOL(J)+SOLINT(I)*FDNA(J)
             FUPSOL(J)=FUPSOL(J)+SOLINT(I)*FUPA(J)

             if (J.eq.1) then ! Print out fupsol and fdnsol to calculate planetary albedo per wavelength bin
C                          write(93,*) FUPSOL(J), FDNSOL(J)
                         elseif(J.eq.ND)then  ! Print out fluxes per wavelength
                         write(91,*) FDNSOL(J)-FUPSOL(J),SOLINT(I)
c                    print *, SOLINT(I)*FDNA(J) - SOLINT(I)*FUPA(J)
                         ENDIF
C-rr         PLAYING AROUND WITH FLUXES HERE
c             IF ((I==26).and.(J==95)) THEN
c             print *, 'FDNSOL=',FDNSOL(J),'FUPSOL=',FUPSOL(J)
c             ENDIF
             Fdiff(J) = SOLINT(I)*(FDNA(J)-FUPA(J))
c              CPCO2 = 7.7 + 5.3E-3*TF(J) - 8.3E-7*TF(J)*TF(J)
c-rr              Putting new CPCO2 parametrization here. Why is there no FCH4?
c              CPCO2 = 5.89 + 6.06E-3*TF(J) + 2.39E-5*TF(J)*TF(J)
c     &        -3.44E-8*TF(J)*TF(J)*TF(J)
c                if (j.eq.1)print *, 'CPCO2=', CPCO2
c                 CPN2 = 6.76 + 6.06E-4*TF(J) + 1.3E-7*TF(J)*TF(J)
c-rr              CPO2 = 8.27 + 2.58E-4*TF(J) - 1.877E5/TF(J)/TF(J)
c              Putting new CPO2 parametrization 3/28/11
c              CPO2 = 7.47 -4.84E-3*TF(J) + 1.38E-5*TF(J)*TF(J)
c     &        -8.73E-9*TF(J)*TF(J)*TF(J) - 1.76E-9/TF(J)/TF(J)
c                 CPO2 = AMAX1(CPO2,CPN2)
c              CPN = FCO2*CPCO2 + FN2*CPN2 + FO2*CPO2 + FAR*4.97
c              CPNT(J) = CPN*4.18*1.E7/DM
c            IF(J.eq.ND) CPNT(J)= 50.* 4.18*1.E7


!             print *,FDNSOL(J),FUPSOL(J),J


 1210    CONTINUE

        IF (LAST .EQ. 1) THEN
         DO J=1,NLAYERS
           DFdiff=Fdiff(J+1)-Fdiff(J)
           SolHeat(I,J)=-(DFdiff*GNEW(J)/CPNT(J)
     &                   /(PF1(J+1)-PF1(J))*86400)
         END DO
        END IF
c      print 99999, i,fupa(1)
c 99999 format(1x,'i =',i2,2x,1pe10.3) !EWS - label not used


c        IF (NST.eq.100)then
c        print *, TAUG(IL),ALAMBDA(I),I
c        endif

 1155 CONTINUE                   !**END WAVELENGTH LOOP**

        IF (LAST .EQ. 1) THEN
        do k=1,5
         istart = (k-1)*9 + 1
         istop = istart + 8
         istop = min0(istop,nsol)
C            WRITE(96,456)(i,i=istart,istop)
C         do j=1,nlayers
C           WRITE(96,457)PF(j),(SolHeat(i,j),i=istart,istop)
C         enddo
        enddo
        CLOSE(96)
        END IF
 456        FORMAT(9x,9I10)
 457    FORMAT(1P10E10.2)



      RETURN
      END
