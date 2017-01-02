      SUBROUTINE RAINOUT(JTROP,NRAIN,USETD)
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      real*8 mass
      CHARACTER*8 ISPEC
      DIMENSION HEFF(NQ),IPVT(NAQ),DJAC(NAQ,NAQ),
     2  F(NAQ),FP(NAQ),TAQ(NZ),X(NAQ)
      character*8 PLANET
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PHOTABLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/BBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/GBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/NBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/WBLOK.inc'
      DATA LSO2g,LH2COg,LSO2aq,LH2COaq,LHCO3_,LCO3_2,LHSO3_,
     1  LSO3_2,LH2COSO3,LOH_/1,2,3,4,5,6,7,8,9,10/




c  I am wanting to adjust rainout downwards for another planet
c  or for lower temperatures generally
c  Giorgi and Chameides use WH2O for the current rainout rate.
c   they quote 3.3e-6 g/cm2/s  as the global average when WH2O is
c  integrated over altitude 
c
C      THIS SUBROUTINE CALCULATES RAINOUT RATES USING GIORGI AND
C   CHAMEIDES (1985) MODEL.  FIRST, IT CALCULATES THE NORMAL HENRY'S
C   LAW COEFFICIENTS.  THEN IT SOLVES A SYSTEM OF NAQ AQUEOUS PHASE
C   REACTIONS TO FIND EFFECTIVE HENRY'S LAW COEFFICIENTS.  THE
C   REACTIONS ARE:
C
C     1)  (SO2)g + ALPHARAIN*[(SO2)aq  +  HSO3-  +  SO3=  +  CH2OHSO3-]
C                  =  (SO2)go
C     2)  (SO2)g   =  (SO2)aq
C     3)  (H2CO)g  =  CH2(OH)2
C     4)  (CO2)aq  =  HCO3-  +  H+
C     5)  (SO2)aq  =  HSO3-  +  H+
C     6)   HCO3-   =  CO3=   +  H+
C     7)   HSO3-   =  SO3=   +  H+
C     8)  CH2(OH)2  +  HSO3-  =  H2O  +  CH2OHSO3-
C     9)   H2O     =  H+  +  OH-
C    10)  (H2CO)g + ALPHARAIN*[CH2(OH)2  +  CH2OHSO3-]  =  (H2CO)go

C   ALONG WITH
C         (CO2)g   =  (CO2)aq
C         (H2SO4)g =  2H+  +  SO4=
C         H+  =  OH-  +  HCO3-  +  HSO3-  +  CH2OHSO3-  +
C                2*[CO3=  +  SO3=]
C
C     THE VARIABLES IN THE NEWTON STEP ARE:
C     1)  X(1)  =  (SO2)g
C     2)  X(2)  =  (H2CO)g
C     3)  X(3)  =  (SO2)aq
C     4)  X(4)  =  CH2(OH)2
C     5)  X(5)  =  HCO3-
C     6)  X(6)  =  CO3=
C     7)  X(7)  =  HSO3-
C     8)  X(8)  =  SO3=
C     9)  X(9)  =  CH2OHSO3-
C    10)  X(10) =  OH-
C
C   FIRST DEFINE RELEVANT CONSTANTS\



      EPS = 1.E-7         !Shawn has this as 1e-4
      INEWT = 20
      GAM15 = 8.64E+05/2.0
      GAM8 = 7.0E+06/2.0
      AV = 6.02E+23
      WL = 1.0
      R = 1.36E-22
      NH = JTROP  !height index of troposphere is called NH in this subroutine...
    
      NH1 = NH + 1



C   MODIFY TEMPERATURE PROFILE TO DO AQUEOUS CHEMISTRY
      T_triple = 273.15
      DO 7 I=1,NH
   7  TAQ(I) = max(T(I),T_triple)

 
C
C   CALCULATE NORMAL HENRY'S LAW COEFFICIENTS (PHYSICAL DISSOLUTION ONLY)
c   lets play guess the units.  looks like mol/liter/atm 
      IF (NRAIN.GT.0) GO TO 4

      do I=1,NH
       tfac=(1./TAQ(I) - 1./298.)
       HCO2(I) = 3.5E-2 * EXP(2400.*tfac) ! updated
      
       
      do J=1,NQ
       if(ISPEC(J).EQ.'O')   H(J,I) = 0.
       if(ISPEC(J).EQ.'H2O') H(J,I) = 0.
       if(ISPEC(J).EQ.'H')   H(J,I) = 0.
       if(ISPEC(J).EQ.'HCO') H(J,I) = 0. ! I think jim made this up
       if(ISPEC(J).EQ.'CH3') H(J,I) = 0.
       if(ISPEC(J).EQ.'S')          H(J,I) = 0.
       if(ISPEC(J).EQ.'S2')        H(J,I) = 0.
       if(ISPEC(J).EQ.'S4')       H(J,I) = 0.
       if(ISPEC(J).EQ.'S8')       H(J,I) = 0.
      !infinite (CHECK THIS against PAVLOV/KEVIN)
       if(ISPEC(J).EQ.'S8AER') H(J,I) = 7.E11
       if(ISPEC(J).EQ.'HNO') H(J,I) = 7.E11
       ! I think jim made the below up
       if(ISPEC(J).EQ.'HS')        H(J,I) = 1.E5
       ! I think jim made the below up
       if(ISPEC(J).EQ.'SO')        H(J,I) = 1.9E-3
       ! updated, infinite
       if(ISPEC(J).EQ.'SO3')      H(J,I) = 7.E11
       if(ISPEC(J).EQ.'H2SO4')  H(J,I) = 7.E11
       ! I think jim made the below up
       if(ISPEC(J).EQ.'HSO')      H(J,I) = 9.E3
       !infinite for particle species 
       if(ISPEC(J).EQ.'SO4AER')H(J,I) = 7.E11
       if(ISPEC(J).EQ.'HCAER') H(J,I) = 7.E11   !infinite for particle species 
       if(ISPEC(J).EQ.'HCAER2') H(J,I) = 7.E11   !infinite for particle species 
       if(ISPEC(J).EQ.'O2') H(J,I) = 1.3E-3* EXP(1500.*tfac) ! updated
       if(ISPEC(J).EQ.'OH') H(J,I) = 30. * EXP(4500.*tfac)  ! updated
       if(ISPEC(J).EQ.'HO2') H(J,I) = 4.E3 * EXP(5900.*tfac)  ! updated
       if(ISPEC(J).EQ.'H2O2') H(J,I) = 8.3E4 * EXP(7400.*tfac) ! updated
       if(ISPEC(J).EQ.'H2') H(J,I) = 7.8E-4 * EXP(500.*tfac)  ! updated
       if(ISPEC(J).EQ.'CO') H(J,I) = 1.0E-3 * EXP(1300.*tfac)  ! updated
       if(ISPEC(J).EQ.'H2CO') H(J,I) = 3.2E3 * EXP(6800.*tfac)  ! updated 
       if(ISPEC(J).EQ.'CH4') H(J,I) = 1.4E-3 * EXP(1600.*tfac)  ! updated
       if(ISPEC(J).EQ.'C2H6') H(J,I) = 1.9E-3 * EXP(2300.*tfac)  ! updated
       if(ISPEC(J).EQ.'NO') H(J,I) = 1.9E-3 * EXP(1500.*tfac)  ! updated
       if(ISPEC(J).EQ.'NO2') H(J,I) = 1.2E-2 * EXP(2500.*tfac)  ! updated
       if(ISPEC(J).EQ.'HNO2') H(J,I) = 50. * EXP(4900.*tfac)  ! added, updated
       if(ISPEC(J).EQ.'HNO3') H(J,I) = 2.1e5 * EXP(8700.*tfac)  ! added, updated
       if(ISPEC(J).EQ.'HO2NO2') H(J,I) = 1.2e4* EXP(6900.*tfac)  ! added, updated
       if(ISPEC(J).EQ.'H2S') H(J,I) = 
     $   0.1 * EXP(2000.*tfac) ! updated
       if(ISPEC(J).EQ.'SO2') H(J,I) = 
     $   1.4 * EXP(2900.*tfac)  ! updated
       if(ISPEC(J).EQ.'OCS') H(J,I) = 
     $   0.022 * EXP(2100.*tfac)  ! updated
       if(ISPEC(J).EQ.'CH3SH') H(J,I) =
     $    0.2 * EXP(2800.*tfac)  
       if(ISPEC(J).EQ.'C2H6S') H(J,I) = 
     $    0.48 * EXP(3100.*tfac)  
       if(ISPEC(J).EQ.'C2H6S2') H(J,I) = 
     $    0.96 * EXP(4000.*tfac)  
       if(ISPEC(J).EQ.'CS2') H(J,I) = 
     $    0.055 * EXP(2800.*tfac)  
c no info found for CS, CH3S, HCS. CO solubility is small, as is HCO, so 0's are probably OK. no info on CH3O either.

!none of the carbon species have any rainout terms in Shawn's model - how valid is this?
       if(ISPEC(J).EQ.'CH')   H(J,I) = 0.     
       if(ISPEC(J).EQ.'C2H2')   H(J,I) = 0.   

       if(ISPEC(J).EQ.'CH3O2')   H(J,I) = 0.  
       if(ISPEC(J).EQ.'CH3O')   H(J,I) = 0.   
       if(ISPEC(J).EQ.'CH2CO')   H(J,I) = 0.  
       if(ISPEC(J).EQ.'CH3CO')   H(J,I) = 0.  
       if(ISPEC(J).EQ.'CH3CHO')   H(J,I) = 0. 
       if(ISPEC(J).EQ.'C2H3')   H(J,I) = 0.   
       if(ISPEC(J).EQ.'C2H4')   H(J,I) = 0.   
       if(ISPEC(J).EQ.'C2H2OH')   H(J,I) = 0. 
       if(ISPEC(J).EQ.'C2H4OH')   H(J,I) = 0. 
       if(ISPEC(J).EQ.'C3H8')   H(J,I) = 0.   
       if(ISPEC(J).EQ.'C3H7')   H(J,I) = 0.   
       if(ISPEC(J).EQ.'C3H6')   H(J,I) = 0.   
       if(ISPEC(J).EQ.'C2H5HCO')   H(J,I) = 0.
       if(ISPEC(J).EQ.'C3H5')   H(J,I) = 0.   
       if(ISPEC(J).EQ.'CH2CCH2')   H(J,I) = 0.


c-mc should go through the henry table and look for chlorine as well as hazy species

      enddo
      enddo

C       print *, 'rain 3'

C
C   CALCULATE EQUILIBRIUM CONSTANTS FOR AQUEOUS PHASE REACTIONS
      DO 3 I=1,NH
C   4)  (CO2)aq  =  HCO3-  +  H+
      R4(I) = 4.3E-7 * EXP(-913.*(1./TAQ(I) - 1./298.))
C   5)  (SO2)aq  =  HSO3-  +  H+
      R5(I) = 1.7E-2 * EXP(2090.*(1./TAQ(I) - 1./298.))
C   6)   HCO3-   =  CO3=   +  H+
      R6(I) = 5.6E-11
C   7)   HSO3-   =  SO3=   +  H+
      R7(I) = 6.E-8 * EXP(1120.*(1./TAQ(I) - 1./298.))
C   8)  CH2(OH)2  +  HSO3-  =  H2O  +  CH2OHSO3-
      R8(I) = 1.E5
C   9)   H2O     =  H+  +  OH-
      R9(I) = 1.E-14 * EXP(-6716.*(1./TAQ(I) - 1./298.))
   3  CONTINUE
      
C
      DO 21 J=1,NQ
      DO 21 I=1,NZ
  21  ENHAN(J,I) = 1.


C   NOW ESTIMATE INITIAL CONCENTRATIONS AT GRID STEP 1 ON THE
C     FIRST CALL

        ALPHARAIN = WL * 1.E-9 * 6.02E23/DEN(1)  
        CO2aq = FCO2*DEN(1)*HCO2(1)*R*TAQ(1)      
        HPLUS = SQRT(R4(1)*CO2aq)
        X(LHCO3_) = R4(1)*CO2aq/HPLUS      !5
        X(LCO3_2) = R6(1)*X(LHCO3_)/HPLUS  !6
        X(LOH_) = R9(1)/HPLUS              !10
C      print *, 'rain 4'
        
        SO2g = USOL(LSO2,1)
        SO2aq = H(LSO2,1)*SO2g

        HSO3_ = R5(1)*SO2aq/HPLUS
        SO3_2 = R7(1)*HSO3_/HPLUS
        FAC = SO2g/(SO2g + ALPHARAIN*(SO2AQ + HSO3_ + SO3_2))
        X(LSO2g) = SO2g*FAC         !1
        X(LSO2aq) = SO2aq*FAC       !3
        X(LHSO3_) = HSO3_*FAC       !7
        X(LSO3_2) = SO3_2*FAC       !8
C

        H2COg = USOL(LH2CO,1)
        CH2OH2 = H(LH2CO,1)*H2COg

        FHSO3_ = R8(1)*CH2OH2*X(LHSO3_)
        FAC = H2COg/(H2COg + ALPHARAIN*(CH2OH2 + FHSO3_))
        X(LH2COg) = H2COg*FAC    !2
        X(LH2COaq) = CH2OH2*FAC  !4
        X(LH2COSO3) = FHSO3_*FAC !9


   4  CONTINUE     !end skip on first call section


C
C ***** LOOP OVER ALTITUDE *****
      DO 10 I=1,NH            !this is a big loop (NH is the tropopause height index JTROP elsewhere)

      IF (NRAIN .NE. 0) THEN 
        DO 9 K=1,NAQ
 9      X(K) = XSAVE(K,I)
      ENDIF 

      ALPHARAIN = WL * 1.E-9 * 6.02E23/DEN(I)

       if(USETD.EQ.1) then
        SO4_2 = (USOL(LH2SO4,I) + PARTICLES(I,LSO4AER-NQ))/ALPHARAIN  
       else
        SO4_2 = (USOL(LH2SO4,I) + USOL(LSO4AER,I))/ALPHARAIN  
       endif   

      SO4SAV(I) = SO4_2
      CO2aq = FCO2*DEN(I)*HCO2(I)*R*TAQ(I)
      
      SO2g0 = USOL(LSO2,I)
      H2COg0 = USOL(LH2CO,I)
      HSO2 = H(LSO2,I)
      HH2CO = H(LH2CO,I)





C
C   START NEWTON ITERATION
      DO 12 IN=1,INEWT
      CALL AQUEOUS(X,F,I)
C
      DO 13 J=1,NAQ
      XS = X(J)
      DX = EPS*X(J)
      X(J) = X(J) + DX
      CALL AQUEOUS(X,FP,I)
C
      DO 14 K=1,NAQ
  14  DJAC(K,J) = (FP(K) - F(K))/DX
  13  X(J) = XS
C     

      CALL SGEFA(DJAC,NAQ,NAQ,IPVT,INFO)
      IF (INFO .NE. 0) THEN 
        PRINT 100, INFO,I,NRAIN
 100  FORMAT(//1X,'NEWTON SOLVER FAILED IN AQUEOUS'/,5X,'INFO =',I3,
     2  '  GRID STEP =',I3,2X,'NRAIN =',I3)
        STOP
      ELSE 
        CALL SGESL(DJAC,NAQ,NAQ,IPVT,F,0)
      ENDIF
C

      LTEST = 0
      DO 22 J=1,NAQ
      X(J) = X(J) - F(J)
      TEST = ABS(F(J)/X(J))
      IF (TEST.GT.1.E-2) LTEST = 1
  22  CONTINUE
      IF (LTEST.EQ.0) GO TO 15    ! loop has converged; jump out of loop 
  12  CONTINUE
C
      PRINT 200,I,NRAIN
 200  FORMAT(//1X,'NEWTON SOLVER FAILED TO CONVERGE IN AQUEOUS'/,5X,
     2  'GRID STEP =',I3,2X,'NRAIN =',I3)
      STOP

  15  CONTINUE
C
C   CALCULATE EFFECTIVE HENRY'S LAW COEFFICIENTS (INCLUDING AQUEOUS
C      PHASE REACTIONS)
      DO 18 J=1,NQ
  18  HEFF(J) = H(J,I) + 1.E-30
C
      HEFF(LH2CO) =(X(LH2COaq) + X(LH2COSO3))/X(LH2COg)

      L1=LSO2
      L2=LH2S

      HEFF(L1) = (X(LSO2aq) + X(LHSO3_) + X(LSO3_2)
     1  + X(LH2COSO3)) / X(LSO2g)
      HEFF(L2) = H(L2,I)*(1. + 1.1E-7/HPLUS)
      PH(I) = - LOG10(HPLUS)

C
C   SAVE DENSITIES AND CALCULATE ENHANCEMENTS


      DO 16 J=1,NAQ
  16  XSAVE(J,I) = X(J)





      ENHAN(LH2CO,I) = HEFF(LH2CO)/H(LH2CO,I)
      ENHAN(L1,I) = HEFF(L1)/H(L1,I)
      ENHAN(L2,I) = HEFF(L2)/H(L2,I)
c  17  CONTINUE !EWS - not used
C
c what follows is incomprehensible to me
c  they state that the vertical integral of WH2O is the rainout rate
c   hence I should scale rainout by scaling WH2O
c  
C   NOW BEGIN GIORGI AND CHAMEIDES FORMULATION FOR RAINOUT RATES
      ZKM = Z(I)/1.E5   !convert vertical grid into kilometers

c this mod for high CO2
c      eleven = 11.0
c      ZKM = min(ZKM,eleven)
      TEMP = T(I)
C
C  Find appropriate GAMMA
      IF (ZKM.LE.1.51) THEN
         GAMMA = GAM15
      ELSE IF (ZKM.LT.8.) THEN
         GAMMA = GAM15 + (GAM8-GAM15)*((ZKM-1.5)/6.5)
      ELSE
         GAMMA = GAM8
      END IF
C
C  Find WH2O
      IF (ZKM.LE.1.) THEN
        Y = 11.35 + 0.1*ZKM
      ELSE
        Y = 11.5444 - 0.085333*ZKM - 9.1111E-03*ZKM*ZKM
      ENDIF
      WH2O = 10.**Y
 
      if (PLANET .EQ. 'EARTH') then 
       wh2o = 1.0 * wh2o  ! nominal earthly rain
c       wh2o=1e-9*wh2o  !ATACAMA
      else if (PLANET .EQ. 'MARS') then
       wh2o = 1e-9*wh2o   !turn off the rain  
      else if (PLANET .EQ. 'DRY') then
       wh2o = 1e-9*wh2o   !turn off rain for dry planet - EWS 9/14/2015
      endif   
 
C  Find F(Z)
      IF (ZKM.LE.1.51) THEN
        FZ = 0.1
      ELSE
        FZ = 0.16615 - 0.04916*ZKM + 3.37451E-3*ZKM*ZKM
      ENDIF
C
C  Loop over species
      DO 10 J=1,NQ
      RKJ = WH2O/55. /(AV*WL*1.E-9 + 1./(HEFF(J)*R*TEMP))
      QJ = 1. - FZ + FZ/(GAMMA*RKJ) * (1.0 - EXP(-RKJ*GAMMA))
  10  RAINGC(J,I) = (1. - EXP(-RKJ*GAMMA))/(GAMMA*QJ)
C ***** END ALTITUDE LOOP *****

C
c-mc set rainout rates to zero above the tropopause (nh=jtrop)

      DO 11 I=NH1,NZ
      DO 11 J=1,NQ
  11  RAINGC(J,I) = 0.


C
C ***** OLD (FISHMAN AND CRUTZEN) RAINOUT RATE *****
C   (USED FOR SCALING THE VERTICAL DISTRIBUTION OF LIGHTNING)
c   what is the 2.4e-6 constant? 
      DO 2 I=1,NH
        ZKM = Z(I)/1.E5
        RAIN(I) = 2.4E-6*EXP((6.-ZKM)/2.42)
        IF(ZKM.LT.6.) RAIN(I) = 2.4E-6
   2  CONTINUE
 
      DO 6 I=NH1,NZ
   6  RAIN(I) = 0.

      NRAIN = NRAIN + 1

      RETURN
      END
