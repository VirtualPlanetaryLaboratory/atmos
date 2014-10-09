      SUBROUTINE PHOTO(ZY,AGL,LTIMES,ISEASON,IZYO2,IO2,INO,N,timega,
     $                 frak,msun)
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      real*8 mass
      Integer nw     !nw - number of wavelengths (nw<kw) set in gridw.f
      real*8 SQ(kj,nz,kw),TTOT(NZ),S(NZ),SALL(kw,NZ), STAU(kw),smax(kw) 
      real*8 temp(kj)
      character*8 REACTYPE,PLANET,CHEMJ,ISPEC
      CHARACTER*20 fmtstr
      DIMENSION SO2(NZ),SO3(NZ),SIGL(NZ),SIGNOL(NZ),CL(NZ),
     2  CNO(NZ),D0(2),CL2(NZ)
      dimension LLNO(35),ANO(9,2),BNO(5,2)

      dimension SIGR(NZ)

      dimension volmix(10,nz),icomp(10,nz),SIGR2(NZ),ncomp(nz)
!above new for Rayleigh Scatering

      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PHOTABLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/AERBLK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/BBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/CBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/JBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/QBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/RBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/SBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/NBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/ISOBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/MBLOK.inc'


      dimension columndepth(KJ,NZ)
      dimension PLOG(NZ)

c testing
      dimension ANEW(NR,NZ),SNEW(NZ)
      dimension PRATESNO(NZ)  !used in High resolution model



C   NO PREDISSOCIATION WAVELENGTHS (ALLEN AND FREDERICK, 1982)
c- this used only if INO=0
      DATA LLNO/3*0, 2*2, 3*0, 2*1, 25*0/

      PM = 1.67E-24               !used in S8 (???) - what is this?
      BK = 1.38E-16               !Boltzmann constant - in erg/K
      RGAS = 8.3143E7             !erg/mol/K
      
!compute mean molar mass (this could easily be abstracted so that N2 isn't hardcoded...
!updating for isotopes, but keeping rest as status quo for now...
      if(ISOTOPE.EQ.1) then
       pO2=UINERT(LO2-Loff,1)
      else
       pO2=USOL(LO2,1)
      endif 
      pCO2=FCO2


      WT = pO2*32. + pCO2*44. + (1.-pO2-pCO2)*28. + 0.4   !(g) mean molar mass
      RMG = RGAS/(WT*G)           !gm cm^2/s^2/mol/K  / g *s^2/cm ->  cm/mol/K

c-mc allows for high CO2 and O2, but forces N2 to decrease if CO2 increases.  Wrong as CO2 gets large
c-mc also makes one realize that fixed mixing ratio isn't the best way to approach high CO2 levels.
c-mc can never get more than 1 bar. - taking CO2 to 0.5 increase mean molar mass to 36.4 (from 28.56 at 0.01)


      PI = 3.14159
      ZYR = ZY*PI/180.            ! note ZY is passed in subroutine call - solar angle in radians
      U0 = COS(ZYR)
      AM = 1./U0
c -mc now in PHOTABLOK      ALB = 0.25                  ! albedo of surface 

      NPHOT = NPHOT + 1        ! counts calls to this subroutine
      HC = 6.625E-27 * 3.00E10     !planck constant*speed of light, used in EFLUX
c-mc erg s * cm/s -> erg cm

!below are used in old school NO calculations
      D0(1) = 1.2E-6              ! used for NO
      D0(2) = 2.6E-6              ! used for NO
      ALNO = 5.1E7                ! used for old NO
      QKNO = 1.5E-9               ! used for old NO
      DIJ = 1.65E9                ! used for old NO


      IF (N.EQ.1) write(14, 119)    !final step printout
 119  FORMAT(/1X,'ENERGY FLUXES IN W/M/M (NOT DIURNALLY AVERAGED)',//
     2  2X,'L',3X,'WAV',6X,'TAUR',6X,'EFLUX',5X,'GFLUX',5X,'S(1)')


c zero the photolysis rate for each species (j) at each height (i).

      do j=1,kj
         do i=1,nz
            prates(j,i)=0.
!temp
c       if(i.eq.1)print *,ISPEC(photoreac(j)),absorbers(j,i)
         enddo
      enddo
c      stop
C
C ***** CALCULATE COLUMN DEPTHS ABOVE EACH COLLOCATION POINT FOR
C ***** EACH SPECIES THAT ABSORBS PHOTONS

c Kevin's likes nothing above the top
      HA = 0.0 * RMG*T(NZ)       !cm/mol/K * K -> cm/mol
      HAD = 0.0 * HA*DEN(NZ)
      TTOT(NZ) = HAD 

      do k=1,kj   !kj=number of photolysis reactions
         columndepth(k,NZ)=absorbers(k,NZ)*HAD
      enddo

      do k=1,kj
       DO  M=1,NZ1
        I = NZ - M      !run through heights from the top down.

        HA = RMG*0.5*(T(I) + T(I+1))  !scale height RT/MG

c-mc        DZ = Z(I+1) - Z(I)  !ACK - this is good, but should already exist as a vector
c-mc in our new scheme DZ(I)=Z(I)-Z(I-1) so DZ(I+1)=Z(I+1)-Z(I)
C ACK - may have to return when I take this to a variable grid

        EFAC = (1. - EXP(-DZ(I+1)/HA))*DEN(I)*HA     !column depth of each layer
        TTOT(I) = TTOT(I+1) + EFAC     !total column depth above height level I   
        columndepth(k,I)= columndepth(k,I+1)    !species column depth above level I
     $                 + EFAC*SQRT(absorbers(k,i)*absorbers(k,i+1))
       enddo
      enddo

c - note that columndepth is indexed by photoreac so contains duplicate information 
c (i.e. columndepth(1,*) and (2,*)  are both the O2 column depth (assuming O2 is the first photoreaction)

      IF(LTIMES .EQ. 0) then
        CALL INITPHOTO(sq,columndepth,ZY,nw,timega,IO2,msun) 
 !this subroutine returns sq(nj,nz,nw) = cross section * quantum yield for each
 ! photolysis reaction at each height for each wavelength
! it also specifies the wavelength grid and the flux

c GNA
       !if (frak.eq.0) then 
         CALL INITMIE(nw,wavl,frak)
        ! else
        !CALL INITMIEFRAC(nw,wavl,frak)
      !endif


      endif  




C   REPEAT THIS SECTION ONLY IF PRESSURE AND TEMPERATURE VARY WITH TIME      
c-mc the following two sections require debugging before use...
      if (ISEASON .ge. 3) then  

        JO2_O1D=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'O2     ')
        JO3_O1D=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'O3     ')
        JNO=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'NO     ')


      CALL XS('O2      ',nw,wavl,wav,T,DEN,JO2_O1D,sq,columndepth,zy,
     $         IO2)         

      CALL XS('O3      ',nw,wavl,wav,T,DEN,JO3_O1D,sq,columndepth,zy)         

      CALL XS('NO      ',nw,wavl,wav,T,DEN,JNO,sq,columndepth,zy)         

       !put in calls to other P/T dependent cross sections here...
       !also note that above/below calls to NO and O2 should only be called if IO2=0,INO=0
       !if IO2=1, the O2 cross section is computed by Jim's undocumented exponential sums methodolgy.
       !hopefully this works, at least is computed each round, so should reflect changing O2 levels...

!no3?
      endif 
      

C    REPEAT THIS SECTION ONLY IF SOLAR ZENITH ANGLE OR O2 VARIES WITH TIME
      if (IZYO2 .GE. 1) then

       JO2_O1D=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'O2     ')
       JNO=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'NO     ')

      CALL XS('O2      ',nw,wavl,wav,T,DEN,JO2_O1D,sq,columndepth,zy,
     $        IO2)         

      CALL XS('NO      ',nw,wavl,wav,T,DEN,JNO,sq,columndepth,zy)         

      endif 


c - cross section 'J numbers' that are used in the loop below
c - for now, I am assuming that these species exist.  If they don't
c - the program will crash because the Jnumber will be returned as a '0'
c which will cause the sq(jnumber,..) or prate(jnumber,..) call to crash
c - consider some IF's here, but this would also entail changing the output files, etc.


      if (IO2.EQ.1) then
       JO2=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'O2     ')+1
       !note this is for the O2 + Hv -> O + O reaction,which is the second O2 reaction
      endif 

      if (INO.LE.2) then  !used if INO=0 or INO=1 - on JPL grid only... !actually for now using in high res too...
       JNO=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'NO     ')
      endif 
      
!return to this sulfur stuff in a bit....

!OK - i want to keep this as it is used to track wavelength dependence of photolysis
!in fact, will probably need to expand this down the road.
!for now, I just want to make sure it works for both the isotope and fine grid options...

       JSO2=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'SO2    ')
       if (JSO2.EQ.0) then
        JSO2=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'SXO2     ')
       endif

!so the above works if either SO2 or SXO2 has a photolysis reaction, but not if both (or neither) do
!for now, I don't anticipate this occuring, so leaving it for now as one or the other will probably always be true

!same goes for below

! a numbering scheme that works even if S8 isn't a species.  the key is to just not use JS8L,JS8R,JS8       
       JS8L=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'S8     ')
       if (JS8L.EQ.0) then
        JS8L=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'SXS7   ')
       endif
       JS8R=JS8L+1
       JS8=JS8L+2

C interpolate the optical properties from the Mie code onto the particle radii computed by photochemical code

      !this was in the time-stepping loop of the original code, but that seems like a mistake
C-AP Since all model is in cm we should convert RSTAND
c-mc test test test
cgna - uncommented why isn't it doing this still??
      !DO k=1,34
       !RSTAND(k) = RSTAND(k)/10000.
      !ENDDO


C  Calculate Qext , W0, G, for the current hydrocarbon aerosol distribution
C qext and W0 are hardcoded for sulfate and sulfur in Twostrm.f

c       do j=1,nz
c          print *, RPAR(j,3)
c       enddo
      if (NP.GE.3) then 
      do L=3,NP  !loop over hydrocarbon particles ONLY RIGHT now
      DO I=1,nw   
      DO J=1,NZ                                       !ISOHACK - i will need to interpolate the MIE parameters to the wavelength grid
      DO k=1,33  !ACK - hardcoded num particles (probably OK - this is how the HC grid was computed)  

         IF ((RPAR(J,L).GE.RSTAND(k)).and.(RPAR(J,L).LT.RSTAND(k+1)))   
     2 THEN

      drs = RSTAND(k+1) - RSTAND(k)
      dr  = RPAR(J,L) - RSTAND(k) 

      QEXTT(I,J,L) = QEXTHC(I,k) + ((QEXTHC(I,k+1) - 
     2 QEXTHC(I,k))/drs)*dr
 
      GFT(I,J,L) = GHC(I,k) + ((GHC(I,k+1) - 
     2 GHC(I,k))/drs)*dr

      W0T(I,J,L) = W0HC(I,k) + ((W0HC(I,k+1) - 
     2 W0HC(I,k))/drs)*dr

c      if (J.EQ.1)print *,QEXTT(I,J),W0T(I,J) 
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      enddo
      endif


C ***** ***** ***** START WAVELENGTH LOOP   ***** ***** *****
      do 19 L=1,nw   
         
      Lold=L-10   !hardcoded to JPL wavelength grid - only access if IO2<1, options which require this grid

      KN = 1      !exponentional sum index - reset below if IO2=1
      ALP = 1.    !exponentional sum coefficient - reset below if IO2=1

      IF (IO2.EQ.1 .AND. wavl(L).LE.2041. .AND. wavl(L).GE.1754.) then
         KN = NK(Lold) ! NK(L) are the number of exponential sum coefficients needed for O2
      endif  ! the coefficients are read in as ALPHAP(L,K) (where 1<K<4) and BETA(L,K)

c-mc zero out Rayleigh scattering vectors:
      do i=1,nz
         ncomp(i)=0  !number of "major" species at each height
       do j=1,10     !10 is the (arbitrary) number of major absorbers
        volmix(j,i)=0.0  !mixing ratio of major species
        icomp(j,i)=0     !hardcoded "index" number of major species
       enddo 
      enddo


                                
C   LOOP OVER K'S AT LOW O2 LEVELS (this is a long loop)
c     note that 19 is also target for loop over L 
c     so this is repeated once for L<1754 and L>2041A and NK(L) times for 1754<L<2041

      DO 19 K=1,KN    

        if (k.eq.1) then !compute Rayleigh scattering cross section (as a function of height?)
         do i=1,nz 

          SIGR(i) = SIGRAY(WAV(L)) * (1. + 1.5*pCO2)  !Old rayleigh cross section 
!note that SIGRAY is calculated in MSCAT. So it needs to stick around for a while. Or be moved into Twostr.f


       !set up new Rayleigh scattering vectors
         do j=1,NSP
           if (SL(j,i)/DEN(i).GE. 0.01) then  !if more than 0.1% of atmosphere, consider Rayleigh contribution
              ncomp(i)=ncomp(i)+1

c              if (Z(i)/1e5.eq.107.5) print *, ispec(j)

              volmix(ncomp(i),i)=SL(j,i)/DEN(i)
              
              if (ISPEC(j).eq.'CO2') then
                 icomp(ncomp(i),i)=2
              else if (ISPEC(j).eq.'N2') then
                 icomp(ncomp(i),i)=3
              else if (ISPEC(j).eq.'O2') then
                 icomp(ncomp(i),i)=4
              else if (ISPEC(j).eq.'H2') then
                 icomp(ncomp(i),i)=5
              else if (ISPEC(j).eq.'HE') then
                 icomp(ncomp(i),i)=6
              else               
                 icomp(ncomp(i),i)=1  !use Earth 'air' - better than nothing? hard to know...
               if (wavl(L).eq.2273) then
                  if (tempcount.eq.0) then
                    print *, ISPEC(j),'at ', Z(i)/1e5, 'km is major
     $  species without Rayleigh data - using AIR', SL(j,i)/DEN(i)
                    tempcount=1
                  endif

               endif 
              endif 

c              if (ncomp(i).eq.3) then
c                 print *, WAVL(L),Z(I)/1e5,ISPEC(J),ncomp(i)
c                 print *,ncomp
c                 print *, icomp
c              endif
c some test printing stuff...
c             if (wavl(L).eq.2273) then
c              print *, wavl(L),Z(i),ISPEC(j),ncomp(i),volmix(ncomp(i),i)
c     $                 ,icomp(ncomp(i),i)       
c             endif 

           endif !end loop over major species
         enddo  !end loop over all species in new Rayleigh setup loop

        enddo  !end loop over height in Rayleigh loop
                tempcount=0                 

        call RAYLEIGH(wavl(l)*1e-4,ncomp,icomp,volmix,SIGR)

       endif   !end case for Rayleigh loop



      if (IO2 .EQ. 1) then   !re-compute O2 cross section via exponential sums

        if (wavl(L) .LE. 2041. .AND. wavl(L).GE.1754.) then 
           ALP = ALPHAP(Lold,K)  !ALPHAP(17,4) are coeficients where 1<K<4
            do I=1,NZ
             sq(JO2,I,L)= SO2HZ(Lold) + BETA(Lold,K)
            enddo   
        endif 

        endif !end IO2=1 loop  
 


!mc - there used to be a mechanism for a Beer's law calculation when the optical depth was high or if IO2=0.
!removing this in favor of a permanent Twostr.f. Check the subversion archives if this ever needs to come back

       IKN=1
       IF (K.NE.KN) IKN=0   !output flag

       CALL TWOSTR(SIGR,U0,sq,WAV(L),L,S,N,IKN)  !two stream radiative tranfer
! this returns the Source function S to this code

       FLX = FLUX(L)*AGL*ALP*FSCALE  
      !AGL is diurnal averaging factor, ALP is 1 if out of the SR band or IO2.NE.1
      ! or is the exponential sum coeffiecent if in the SR band and IO2.EQ.1
      ! FLUX is already corrected based on solar age (timeGa set in INPUTFILES/PLANET.dat) 
      !FSACLE adjusts for position in the solar system (set in INPUTFILES/PLANET.dat) 
          !FSCALE=1 is Earth, FSCALE=0.43 is Mars

c compute photlysis rates for each reaction at each height (summed over wavelength)      

      do j=1,kj   
       do i=1,nz  
          prates(j,i) = prates(j,i) + FLX*sq(j,i,L)*S(i)
       enddo
      enddo
    
c save wavelength dependence of SO2 photolysis and optical depth
      do I=1,NZ
        PSO2MC(L,I) = FLX*sq(JSO2,I,L)*S(I)
        SALL(L,I) = S(I)
      enddo

      if (INO.LE.1) then   
C   NO PREDISSOCIATION IN THE D00 (1910 A) AND D10 (1830 A) BANDS
      if (wavl(L).LE.2500.0 .AND. wavl(L).GE.1754.) then
       NOL = LLNO(Lold)         !      DATA LLNO/3*0, 2*2, 3*0, 2*1, 25*0/
       IF (NOL .NE. 0) THEN         ! else bail out of loop over K (GOTO 19)

         IF (INO .EQ. 1) THEN 
C             old (cieslik and nicolet) method with intensities updated to
C              frederick and hudson (1979)
           do I=1,NZ
             RN2 = DIJ/(ALNO + DIJ + QKNO*DEN(I))
             prates(JNO,I)=prates(JNO,I)+ 0.5*D0(NOL)*S(I)*RN2*AGL*ALP  
           enddo
         ELSE
C            frederick and allen method (effective cross sections)
           do I=1,NZ
            prates(JNO,I)=prates(JNO,I) + FLX*SIGNO(I,NOL)*S(I)  
           enddo
        ENDIF

       ENDIF     !if NOL=0, then do nothing...

      endif !end NO wavlength loop
      
      else  !end if loop which restricts this behavior to INO=0 or INO=1
!so ww get here if INO=2
!for now just use the NO photo rate generated from the band model, even at high res (Jim's suggestion)
!the below is dumb - it should be removed from the wavelength loop to make it more clear...
!leaving in place for now as I try to get results by 5PM...
       if(LTIMES.EQ.0) then   
        if(L.eq.1) then
           print *, 'using hardcoded NO photorates'
           print *, 'be sure out.NOprates is valid for this atmosphere'
         open(60, file='out.NOprates',status='OLD')         ! formatted input
         read (60,*) (pratesNO(I),I=1,nz)
         close(60)
        endif   
       endif
       if(L.eq.1) then
         do i=1,nz
          prates(JNO,I)=pratesNO(I)
         enddo 
       endif  
      endif !end INO=2 loop

!print out on last timestep
      if (N .NE. 0 .AND. K.EQ.KN) then !N.NE.0 only on last timestep,  KN=1 or 1-4 for L<17
c        TAUR = SIGR*TTOT(1)
        TAUR = SIGR(1)*TTOT(1)  !ack temporary way of indicating total rayleigh optical depth. should be sum
        DELWAV = WAVU(L) - WAVL(L)
        EFLUX = 1.E6*HC*FLX/(WAV(L)*DELWAV*AGL)  !convert to W/m^2
        GFLUX = EFLUX*S(1)  !ground flux = TOA flux*optical path
        write(14, 120) L,WAV(L),TAUR,EFLUX,GFLUX,S(1)
 120    FORMAT(1X,I3,1X,F6.1,1X,1P4E10.3)
      endif


C ***** ***** ***** END WAVELENGTH LOOP ***** ***** *****  
  19  continue

c      print *, 'stopping after wavelength loop in PHOTO'
c      stop



      if (JS8L.GT.0) then  !if gaseous S8 is in the model, compute the photolysis rate by black magic
C
C ***** CALCULATE S8 PHOTORATE USING ANDY YOUNG'S METHOD *****
C     (ANC IS THE NUMBER OF COLLISIONS REQUIRED TO CLOSE THE RING,
C      QCOL IS THE COLLISION CROSS SECTION OF A MOLECULE, SCOL IS
C      THE COLLISION FREQUENCY, PS8R AND PS8L ARE THE PHOTOLYSIS
C      RATES OF THE RING AND LINEAR S8 MOLECULES, RESPECTIVELY.)
c   there is a bug in computing PS8 - the photolysis cross sections
c   go through "crises" that seem to be unrelated to column depths
      ANC = 1.
      QCOL = 3.E-15
      DO I=1,NZ
       VMEAN = SQRT(8.*BK*T(I)/(WT*PM*PI))
       SCOL = DEN(I)*QCOL*VMEAN
       prates(JS8,I) = prates(JS8R,I) * prates(JS8L,I)/   
     $                (prates(JS8L,I) + SCOL/ANC)       

       prates(JS8L,I)=0.0   !this keeps these predissociation reactions from
       prates(JS8R,I)=0.0   !factoring into the photoylsis/radiative transfer schemes
      
      enddo

      endif !end gaseous S8 loop


      if (N.NE.0) then  !on last timestep only...

!ACK - hardcoded wavelength grid

      DO 301 J=1,NZ
 301     write(27,399) (PSO2MC(L,J),L=11,30)   !ACK  11-30 is 1762-2116.5A
 399  format(20(1PE9.2,2X))

      DO 302 J=1,NZ
 302     write(27,499) (PSO2MC(L,J),L=31,45)   !ACK 31-45 is 2139.5-2516A
 499  format(15(1PE9.2,2X))


!now print out so2 photorates on the high resolution grid
      if (LGRID.eq.1) then
!L-13 is 1786.25 and L-934 is 2246.75
       L1=minloc(wavl,1,wavl.ge.1786.25)  !ACK hardcoded to high resolution grid 
       L2=minloc(wavl,1,wavl.ge.2246.75)  !ACK hardcoded to high resolution grid 
       write(61,*),L1,L2
        fmtstr='(    (1PE9.2,2x))'
        write(fmtstr(2:5),'(I4)')NZ
       do j=L1,L2
        write(61,fmtstr) (PSO2MC(j,i),i=1,nz)
       enddo   

      endif





c - td printout - here we are going to write out where tau=1
c - at each wavelength, need to find maximum height at which tau=1

      opticaldepth=1.0
      slev=EXP(-1.0*AM*opticaldepth)   !optical path length where tau=1 given zenith angle

      smax=MAXLOC(SALL,2, SALL .le. slev)

c orig      STAU=Z(MAXLOC(SALL,2, SALL .le. slev))/1e5   !original code not in loop      
      do L=1,nw
         if (smax(L) .eq. 0) smax(L)=1.      !for some reason MAXLOC started returning 0 for the lower bound
         STAU=Z(INT(smax(L)))/1e5
      enddo

      if(ISOTOPE.EQ.0) write(48,322) (STAU(L), L=1,nw) !write out tau=1 heights for the time-dependent codes

 322  format(118(F10.3))  !ACK a hardcoded wavelength grid


!print cross sections - should disable for production runs
!or better yet, just don't save them in evolve script...

       do L=1,nw  
        do I=1,nz
!         write(29,*) (sq(J,I,L),J=1,kj)   !not needed for now
        enddo
       enddo


c print wavlength and height grids...
c so analysis programs can pick up on nz and nw

       do i=1,nz 
        write(30,*) Z(i)
       enddo 
       

       do L=1,nw 
        if (ISOTOPE.EQ.0) write(31,*) wavl(L),wavu(L),wav(L)
         write(41,*), flux(L),relflux(L)
       enddo
       write(41,*) AM    !write out mu

      do L=1,nw  
       write(41,*) (SALL(L,I),I=1,NZ)
      enddo
      endif




C ***** FILL UP RATE MATRIX *****

      do j=1,kj
       do i=1,nz
          A(INT(photonums(j)),i)=prates(j,i)
       enddo   
      enddo




c      print *,'stopping in PHOTO'
c      stop

      LTIMES = LTIMES + 1
      RETURN
      END
