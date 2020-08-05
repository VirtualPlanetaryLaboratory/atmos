      SUBROUTINE DIFCO
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      real*8 mass,BKMG(NZ)
      character*8 PLANET
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PHOTABLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/AERBLK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/CBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/NBLOK.inc'
      BK= 1.380649E-16! Boltzmann
      NZ1 = NZ - 1
      IF (FH2.GT.0.50) THEN
c-mab: Recomputing FH2 with the most major species VMR = (1-everything major) method:
        FH2  = (1.0-FHE-FH2O-FH-FOH-FCO-FCO2-FCH4)
c-mab: (Note: DENSTY.f is missing the x4 in HE in that version.)
	     WT = FH2*2.0 + FHE*4.0 !From Ravi's version DIFCO...
        print*,'DIFCO: WT for H2/HE dominated atmosphere....'
      ELSE
c-mab: assuming the O2,CO2,N2,and AR are the major species (identical computation in DENSTY)
      	WT = FO2*32. + FCO2*44. + FAR*40.+ (1.-FO2-FCO2-FAR)*28.
        print*,'DIFCO: WT for N2/O2 dominated atmosphere....'
      ENDIF
      DO  I=1,NZ
      BKMG(I) = BK/(1.660539E-24*WTa(I)*GZ(I))    !a good pressure Boltzmann divided by mean molecular weight, g, and inverse Avogadro constant
      ENDDO
C
C ***** DK(I) = K*N AT GRID STEP I+1/2 *****
C
      DO 1 I=1,NZ1
      EDDAV = SQRT(EDD(I)*EDD(I+1))  !average eddy diffusion in grid center
      DENAV = SQRT(DEN(I)*DEN(I+1))  !average density at grid center
   1  DK(I) = EDDAV*DENAV
C
C   COMPUTE DIFFUSION LIFETIME AT EACH HEIGHT (H*H/K)
      DO 2 I=1,NZ
      H = BKMG(I) * T(I)
      HSCALE(I) = H
   2  TAUEDD(I) = H*H/EDD(I)

      do i=1,NZ1
        TAV = SQRT(T(I)*T(I+1))  !average temperature at grid center
        H_atm(i) = BKMG(I) * TAV
c  compute scale heights of all the species in the atmosphere at all
c   heights
        do j=1,NQ
           scale_H(j,i) = BKMG(I) * TAV*WTa(I)/mass(j)
        enddo


c-mab: Proceed with this loop for giant planets instead...
c-mab: transferring the "b" portion of the "DI" expression from Ravi's version since DI = b/n
c-mab: where b is the binary molecular diffusion parameter and DI is the diffusion coefficient
c-mab: got relation from Jim's new book Chapter 5.
! (b has a constant*T^(some power) form for terrestrial gases in Chapter 5. Not established for giant planets?_?)
        do j=1,NQ
      	 bX1X2(j,i) = 1.52E18*(1./mass(j)+1./WTa(I))**0.5*(TAV**0.5)
        if (j.eq.LH)then
c-mab: One of the below options to be used for terrestrial planets.
              if(PLANET .EQ. 'EARTH') then
c                bX1X2(j,i) = 2.7E19*(TAV/200.)**0.75   ! correct for N2
                bX1X2(j,i)=4.8E17*TAV**0.70  !from catling book
c         bX1X2(j,i) = 4.0D18*(TAV/200.)**0.75
              else if (PLANET .EQ. 'MARS') then
c                bX1X2(j,i) = 0.8*1.8*1.4E19*(TAV/200.)**0.75  ! correct for CO2
              bX1X2(j,i)=1.8*31.4/BK*(TAV)**0.75*EXP(-11.7/TAV)   !from Zahnle 2008,using the factor 1.8 to estimate binary diffusion based on a 1.8 factor difference between H and H2 diffusion in He
c         bX1X2(j,i) = 4.0D18*(TAV/200.)**0.75
              else if (PLANET .EQ. 'DRY') then ! assume O2-dominated - this is the same as N2 for now
c                bX1X2(j,i) = 2.7E19*(TAV/200.)**0.75   ! correct for N2
c         bX1X2(j,i) = 4.0D18*(TAV/200.)**0.75
              endif
        endif
        if (j.eq.LH2)then
c-mab: One of the below options to be used for terrestrial planets.
              if(PLANET .EQ. 'EARTH') then
c                bX1X2(j,i) = 1.4E19*(TAV/200.)**0.75  ! correct for N2
                bX1X2(j,i)=2.7E17*TAV**0.75  !from catling book, for air
c        bX1X2(j,i) = 4.0D18*(TAV/200.)**0.75
              else if (PLANET .EQ. 'MARS') then
c               bX1X2(j,i) = 0.8*1.4E19*(TAV/200.)**0.75  ! correct for CO2
               bX1X2(j,i)=31.4/BK*(TAV)**0.75*EXP(-11.7/TAV) !from Zahnle 2008
c        bX1X2(j,i) = 4.0D18*(TAV/200.)**0.75
              else if (PLANET .EQ. 'DRY') then ! assume O2-dominated - this is the same as N2 for now
c               bX1X2(j,i) = 1.4E19*(TAV/200.)**0.75  ! correct for N2
c        bX1X2(j,i) = 4.0D18*(TAV/200.)**0.75
              endif
        endif
c-mab: Should we put an error check for cases that don't fall above or below?

        enddo
! Note: there will be a separate expression in output.f for this scenario as well
      enddo

      H_atm(nz) = BKMG(NZ) * TAV
      do j=1,NQ
        scale_H(j,NZ) = BKMG(NZ) * T(NZ)*WTa(NZ)/mass(j)
      enddo


        do j=1,NQ
      	 bX1X2(j,NZ) = 1.52E18*(1./mass(j)+1./WTa(NZ))**0.5*(T(NZ)**0.5)
        if (j.eq.LH)then
c-mab: One of the below options to be used for terrestrial planets.
              if(PLANET .EQ. 'EARTH') then
c                bX1X2(j,NZ) = 2.7E19*(T(NZ)/200.)**0.75   ! correct for N2
                bX1X2(j,NZ)=4.8E17*T(NZ)**0.70  !from catling book
c         bX1X2(j,NZ) = 4.0D18*(T(NZ)/200.)**0.75
              else if (PLANET .EQ. 'MARS') then
c                bX1X2(j,NZ) = 0.8*1.8*1.4E19*(T(nz)/200.)**0.75  ! correct for CO2
              bX1X2(j,NZ)=1.8*31.4/BK*(T(NZ))**0.75*EXP(-11.7/T(NZ))   !from Zahnle 2008,using the factor 1.8 to estimate binary diffusion based on a 1.8 factor difference between H and H2 diffusion in He
c         bX1X2(j,NZ) = 4.0D18*(T(NZ)/200.)**0.75
              else if (PLANET .EQ. 'DRY') then ! assume O2-dominated - this is the same as N2 for now
c                bX1X2(j,NZ) = 2.7E19*(T(NZ)/200.)**0.75   ! correct for N2
c         bX1X2(j,NZ) = 4.0D18*(T(NZ)/200.)**0.75
              endif
        endif
        if (j.eq.LH2)then
c-mab: One of the below options to be used for terrestrial planets.
              if(PLANET .EQ. 'EARTH') then
c                bX1X2(j,NZ) = 1.4E19*(T(NZ)200.)**0.75  ! correct for N2
                bX1X2(j,NZ)=2.7E17*T(NZ)**0.75  !from catling book, for air
c        bX1X2(j,NZ) = 4.0D18*(T(nz)/200.)**0.75
              else if (PLANET .EQ. 'MARS') then
c               bX1X2(j,NZ) = 0.8*1.4E19*(T(nz)/200.)**0.75  ! correct for CO2
               bX1X2(j,NZ)=31.4/BK*(T(NZ))**0.75*EXP(-11.7/T(NZ)) !from Zahnle 2008
c        bX1X2(j,NZ) = 4.0D18*(T(NZ)/200.)**0.75
              else if (PLANET .EQ. 'DRY') then ! assume O2-dominated - this is the same as N2 for now
C               bX1X2(j,NZ) = 1.4E19*(T(NZ)/200.)**0.75  ! correct for N2
c        bX1X2(j,NZ) = 4.0D18*(T(NZ)/200.)**0.75
              endif
        endif
        enddo
  !output.f impacted



      RETURN
      END
