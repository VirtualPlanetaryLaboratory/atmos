      SUBROUTINE DIFCO(FO2)
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      real*8 mass
      character*8 PLANET
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PHOTABLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/AERBLK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/CBLOK.inc'
      NZ1 = NZ - 1
      WT = FO2*32. + FCO2*44. + FAR*40.+ (1.-FO2-FCO2-FAR)*28.     !assuming the O2,CO2,N2,and AR are the major species (identical computation in DENSTY)
      BKMG = 1.38E-16/(1.67E-24*WT*G)    !a good pressure
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
      H = BKMG * T(I)
      HSCALE(I) = H
   2  TAUEDD(I) = H*H/EDD(I)
   
      do i=1,NZ1
        TAV = SQRT(T(I)*T(I+1))  !average temperature at grid center
        H_atm(i) = BKMG * TAV
c  compute scale heights of all the species in the atmosphere at all
c   heights
        do j=1,NQ
           scale_H(j,i) = BKMG * TAV*WT/mass(j)
        enddo

      if (PLANET .EQ. 'EARTH') then 
        bHN2(i) = 2.7E19*(TAV/200.)**0.75   ! correct for N2
        bH2N2(i) = 1.4E19*(TAV/200.)**0.75  ! correct for N2
c       bXN2(i) = 4.0D18*(TAV/200.)**0.75
      else if (PLANET .EQ. 'MARS') then
        bHN2(i) = 0.8*1.8*1.4E19*(TAV/200.)**0.75  ! correct for CO2
        bH2N2(i) = 0.8*1.4E19*(TAV/200.)**0.75  ! correct for CO2
c       bXN2(i) = 4.0D18*(TAV/200.)**0.75
      else if (PLANET .EQ. 'DRY') then ! assume O2-dominated - this is the same as N2 for now
        bHN2(i) = 2.7E19*(TAV/200.)**0.75   ! correct for N2
        bH2N2(i) = 1.4E19*(TAV/200.)**0.75  ! correct for N2
c       bXN2(i) = 4.0D18*(TAV/200.)**0.75
      endif   



      enddo
      
      H_atm(nz) = BKMG * TAV
      do j=1,NQ
        scale_H(j,NZ) = BKMG * T(NZ)*WT/mass(j)
      enddo

      if (PLANET .EQ. 'EARTH') then 
      bHN2(nz) = 2.7E19*(T(nz)/200.)**0.75   ! correct for N2
      bH2N2(nz) = 1.4E19*(T(nz)/200.)**0.75  ! correct for N2
c     bXN2(nz) = 4.0D18*(T(nz)/200.)**0.75
      else if (PLANET .EQ. 'MARS') then
      bHN2(nz) = 0.8*1.8*1.4E19*(T(nz)/200.)**0.75  ! correct for CO2
      bH2N2(nz) = 0.8*1.4E19*(T(nz)/200.)**0.75     ! correct for CO2
c     bXN2(nz) = 4.0D18*(T(nz)/200.)**0.75
      else if (PLANET .EQ. 'DRY') then ! assume O2-dominated - this is the same as N2 for now
        bHN2(i) = 2.7E19*(TAV/200.)**0.75   ! correct for N2
        bH2N2(i) = 1.4E19*(TAV/200.)**0.75  ! correct for N2
c       bXN2(i) = 4.0D18*(TAV/200.)**0.75
      endif    



      RETURN
      END
