      SUBROUTINE DIFCO
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      real*8 mass
      character*8 PLANET
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PHOTABLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/AERBLK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/CBLOK.inc'
      NZ1 = NZ - 1
      IF (FH2.GT.0.50) THEN
c-mab: Recomputing FH2 with the most major species VMR = (1-everything major) method:
        FH2  = (1.0-FHE-FH2O-FCO-FCO2-FCH4)
c-mab: (Note: DENSTY.f is missing the x4 in HE in that version.)
	WT = FH2*2.0 + FHE*4.0 !From Ravi's version DIFCO...
        print*,'DIFCO: WT for H2/HE dominated atmosphere....'      
      ELSE
c-mab: assuming the O2,CO2,N2,and AR are the major species (identical computation in DENSTY)
      	WT = FO2*32. + FCO2*44. + FAR*40.+ (1.-FO2-FCO2-FAR)*28.
        print*,'DIFCO: WT for N2/O2 dominated atmosphere....'      
      ENDIF
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
      else if (PLANET .EQ. 'WASP12B') then 
! assuming the "air" composition here is whatever went to the WT (H2 + HE)
! The general form for b here is from DI = b/n expression in Kopparapu et al. 2012 work
! (b has a A*T^(n) form for terrestrial gases in Chapter 5. 
! A and n are not presently established for giant planets...
       do j=1,NQ
      	bX1X2(j,i) = 1.52E18*(1./mass(j) + 1./WT)**0.5*(TAV**0.5)
       enddo 
! Note: Only works properly if diffusion expressions are updated in 
! TOTCtester and output.f for this scenario as well      
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
      else if (PLANET .EQ. 'WASP12B') then !see notes from previous loop on this
       do j=1,NQ
      	bX1X2(j,nz) = 1.52E18*(1./mass(j) + 1./WT)**0.5*(T(nz)**0.5)
       enddo       
      endif    



      RETURN
      END
