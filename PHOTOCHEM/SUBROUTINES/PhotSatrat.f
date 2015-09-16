      SUBROUTINE PHOTSATRAT(JTROP,H2O)
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      real*8 mass
      DIMENSION HL(NZ),H2O(NZ),A(NZ)
      character*8 PLANET
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PHOTABLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/SATBLK.inc'

C
C   THIS SUBROUTINE CALCULATES THE SATURATION VAPOR PRESSURE OF WATER
C   FROM MAGNUS' EQUATION (SEE RUNAWAY GREENHOUSE PAPER IN ICARUS -
C   APPENDIX A).  IT THEN FIXES TROPOSPHERIC H2O USING A MANABE/
C   WETHERALD RELATIVE HUMIDITY DISTRIBUTION.
C
      T0 = 273.15
      Pzero = 6.103E-3          !vapor pressure at T0 I assume
      AMV = 18.
      VAPL = 597.3
      SUBL = 677.
      R = 1.9872
      A0 = 0.553
      BK = 1.38E-16
      PS = 1.E-6 * DEN(1)*BK*T(1)    !this is actually a good measure of pressure
C
      DO 1 J=1,NZ
      HL(J) = SUBL
      A(J) = 0.
      IF (T(J).LT.T0) GO TO 1
      HL(J) = VAPL + A0*T0
      A(J) = A0
   1  CONTINUE
C
C   FIND SATURATION VAPOR PRESSURE
      DO 2 J=1,NZ
      P1 = Pzero * (T0/T(J))**(AMV*A(J)/R)
      P2 = EXP(AMV*HL(J)/R * (1./T0 - 1./T(J)))
      PV = P1 * P2
      P(J) = 1.E-6 * DEN(J)*BK*T(J)    ! bars
   2  H2OSAT(J) = PV/P(J)
C
C   CALCULATE TROPOSPHERIC H2O CONCENTRATIONS
      DO 3 J=1,JTROP
      if (PLANET .EQ. 'EARTH') then 
       REL = 0.77 * (P(J)/PS - 0.02)/0.98   !manabe formula

       
c       REL = 1.0 ! Saturated troposphere - eventual should abstract this and ATACAMA

c       REL = 0.17   !test ATACAMA
      else if (PLANET .EQ. 'MARS') then
c      REL = 0.12                         ! 7 microns
       REL = 0.17                         ! this is good for 9.5 micrometers  
! this makes the model very dry  - warning nonstandard!!!      

      else if (PLANET .EQ. 'DRY') then !Dry planets. Implemented originally for post-runaway, O2-rich atmospheres
       REL= 1.e-8 ! For H2O-free atmospheres. Set above zero to prevent floating point errors. - EWS - 9/14/2015
      endif   

      RELH(J) = REL 
   3  H2O(J) = REL * H2OSAT(J)
C


      RETURN
      END
