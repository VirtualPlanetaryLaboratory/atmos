      SUBROUTINE DENSTY(FO2,poop3)
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
!      R0 = 6.371E8  (now set in PLANET.dat
      FT = FO2 + FCO2 + FAR
      WT = FO2*32. + FCO2*44. + (1.-FT)*28. + FAR*40.   !assuming O2,CO2,N2 and Ar are the main players and pN2=1-pO2+pCO2+pAr
      !this is bad because CO and O can be important up high.  (but maybe this doesn't really matter)

      ROVERM = RGAS/WT

!      PG = 1.013E6    !primary specification of pressure in the model
!      P0 = PG  (P0 now specified in PLANET.dat)
C     P0 = PG GIVES YOU A ONE BAR ATMOSPHERE

c-mc      DZ = Z(2) - Z(1)   !ACK

      T0 = T(1) + (T(1)-T(2))/2.
      HA = ROVERM*0.5*(T0 + T(1))/G0
      P1 = P0 *1e6 * poop3 * EXP(-0.5*DZ(1)/HA)
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
