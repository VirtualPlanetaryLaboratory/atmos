      SUBROUTINE DENSTY
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      real*8  mass
      character*8 PLANET
      dimension sumLL(NZ), ROVERM(NZ),HA(NZ)
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PHOTABLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/comPRESS1.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/BBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/NBLOK.inc'

C
C   THIS SUBROUTINE CALCULATES ATMOSPHERIC NUMBER DENSITIES, ASSUM-
C   ING HYDROSTATIC EQUILIBRIUM
      G0 = G
      RGAS = 8.3143E7
      BK = 1.380649E-16

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

        Do I=1,NZ  !STB
        sumLL(I)=0 !sum over molecular weights of LL species weighted by their mixing ratio
        sumUSOL(I)=0 !sum over mixing ratios for the LL species
        DO J=1,(NQ-NP)
        sumLL(I)= sumLL(I)+USOL(J,I)*mass(J)
        sumUSOL(I)=sumUSOL(I)+USOL(J,I)
        ENDDO
        If (LCO2.le.NQ) then !to ensure that CO2 is counted, e.g. for Mars
        WTa(I) = sumLL(I) + (1.-sumUSOL(I)-FAR)*28.+ FAR*40.
c        print*, WTa(I), 'This is WTa(i)'
        else
        WTa(I) = sumLL(I) + (1.-sumUSOL(I)-FAR-FCO2)*28.
     &         + FCO2*44 + FAR*40.
        ENDIF
        Enddo
        print*, "Calculated WTa for O2-N2 dominated atmosphere..."
        print*, "Surface WTa = ",WTa(1)
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

       Do I=1,NZ

        WTa(I) = FH2*2.0 + FHE*4 +  FCO*28.0 + FH2O*18.0 + FH + FCH4*16.0
     .             + FCO2*44.0
c=mab: Review above expression prior to release
      Enddo
        print*, "Calculated WT for He/H2 dominated atmosphere..."


    	print*,"Molecular weight of atmosphere = ",WTa(1)
      ENDIF
      WT=WTa(1) !insurance policy for when WT is used in a part of the code I'm unaware of

       Do I=1,NZ
       ROVERM(I) = RGAS/WTa(I)
       Enddo

!      PG = 1.013E6    !primary specification of pressure in the model
!      P0 = PG  (P0 now specified in PLANET.dat)
C     P0 = PG GIVES YOU A ONE BAR ATMOSPHERE

c-mc      DZ = Z(2) - Z(1)   !ACK

      T0 = T(1) + (T(1)-T(2))/2.
      HA(1) = ROVERM(1)*0.5*(T0 + T(1))/G0
      P1 = P0 *1e6 * EXP(-0.5*DZ(1)/HA(1))
      DEN(1) = P1/(BK*T(1))
      GZ(1)=G0
C
C ***** FIND DENSITY FROM HYDROSTATIC EQUILIBRIUM *****
      DO I=2,NZ
c-mc      DZ = Z(I) - Z(I-1)
      R = R0 + Z(I)
      GZ(I) = G0 * (R0/R)*(R0/R)
      TAV = 0.5*(T(I) + T(I-1))
      HA(I) = ROVERM(I)*TAV/GZ(I)
      DEN(I) = DEN(I-1)*EXP(-DZ(I)/HA(I))*T(I-1)/T(I)
      END DO
C


C ***** FIND PRESSURE FROM THIS DENSITY *********
      DO 334 I=1,NZ
        PRESS(I)=DEN(I)*BK*T(I)
 334  CONTINUE
C
      RETURN
      END
