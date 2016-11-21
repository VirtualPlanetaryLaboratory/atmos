      SUBROUTINE DENSTY
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      real*8  mass
      character*8 PLANET
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PHOTABLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/comPRESS1.inc'

C
C   THIS SUBROUTINE CALCULATES ATMOSPHERIC NUMBER DENSITIES, ASSUM-
C   ING HYDROSTATIC EQUILIBRIUM
      G0 = G
      RGAS = 8.3143E7
      BK = 1.38054E-16

c-mab: Using radius to distinguish terrestrial from giants here
      IF (R0.LT.(5E9)) THEN 
c-mab: giving a sub-saturn radius value as the marker here...
        FT = FO2 + FCO2 + FAR
c-mab: For error checking
        if(FO2.eq.0.)print*,"Warning: Check FO2 value, as FO2 =",FO2
        if(FCO2.eq.0.)print*,"Warning: Check FCO2 value, as FCO2 =",FCO2
        if(FAR.eq.0.)print*,"Warning: Check FAR value, as FAR =",FAR

        if(FO2.lt.0.OR.FCO2.lt.0.OR.FAR.lt.0) then
         print*,"Mixing ratios should not be negative!"
         print*,"FT,FO2,FCO2,FAR = ",FT,FO2,FCO2,FAR
         print*,"Aborting..."
         stop
        endif


        WT = FO2*32. + FCO2*44. + (1.-FT)*28. + FAR*40.   !assuming O2,CO2,N2 and Ar are the main players and pN2=1-pO2+pCO2+pAr

        print*, "Calculated WT for O2-N2 dominated atmosphere..."
C this is bad because CO and O can be important up high.
C (but maybe this doesn't really matter)
      ELSE
c-mab: For error checking
        FH2  = (1.0-FHE-FH2O-FH-FOH-FCO-FCO2-FCH4)
        if(FAR.ne.0.)then
         print*,"FAR =",FAR
         print*,"Warning: FAR should NOT exist in this template."
        endif
        if(FHE.eq.0.)print*,"Warning: Check FHE value, as FHE =",FHE
        if(FH2O.eq.0.)print*,"Warning: Check FH2O value, as FH2O =",FH2O
        if(FCO.eq.0.)print*,"Warning: Check FCO value, as FCO =",FCO
        if(FCO2.eq.0.)print*,"Warning: Check FCO2 value, as FCO2 =",FCO2
        if(FCH4.eq.0.)print*,"Warning: Check FCH4 value, as FCH4 =",FCH4

        if(FHE.lt.0.OR.FH2O.lt.0.OR.FCO.lt.0
     .             .OR.FCO2.lt.0.OR.FCH4.lt.0) then
         print*,"Mixing ratios should not be negative!"
         print*,"FH2,FHE,FH2O,FCO,FCO2,FCH4 =",FH2,FHE,FH2O,FCO,FCO2,FCH4
         print*,"Aborting..."
         stop
        endif
C-mab: Uncomment below for debugging
C        print*,"FH2,FHE,FCO2,FCO,FH2O,FH,FOH,FCH4",
C     .             FH2,FHE,FCO2,FCO,FH2O,FH,FOH,FCH4
 
        WT = FH2*2.0 + FHE*4 +  FCO*28.0 + FH2O*18.0 + FH + FCH4*16.0 
     .             + FCO2*44.0 
c=mab: Review above expression prior to release

        print*, "Calculated WT for He/H2 dominated atmosphere..."
      ENDIF
      
      	print*,"Molecular weight of atmosphere = ",WT

      ROVERM = RGAS/WT

!      PG = 1.013E6    !primary specification of pressure in the model
!      P0 = PG  (P0 now specified in PLANET.dat)
C     P0 = PG GIVES YOU A ONE BAR ATMOSPHERE

c-mc      DZ = Z(2) - Z(1)   !ACK

      T0 = T(1) + (T(1)-T(2))/2.
      HA = ROVERM*0.5*(T0 + T(1))/G0
      P1 = P0 *1e6 * EXP(-0.5*DZ(1)/HA)
      DEN(1) = P1/(BK*T(1))
C
C ***** FIND DENSITY FROM HYDROSTATIC EQUILIBRIUM *****
      DO 1 I=2,NZ
c-mc      DZ = Z(I) - Z(I-1)
      R = R0 + Z(I)
      GZ = G0 * (R0/R)*(R0/R)                              
      TAV = 0.5*(T(I) + T(I-1))
      HA = ROVERM*TAV/GZ
   1  DEN(I) = DEN(I-1)*EXP(-DZ(I)/HA)*T(I-1)/T(I) 
C


C ***** FIND PRESSURE FROM THIS DENSITY *********
      DO 334 I=1,NZ
        PRESS(I)=DEN(I)*BK*T(I)
 334  CONTINUE
C
      RETURN
      END
