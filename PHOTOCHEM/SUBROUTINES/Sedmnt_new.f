      SUBROUTINE SEDMNT(FSULF,USETD,frak,HCDENS, monsize)
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      real*8 mass
      DIMENSION FSULF(NZ),TAUTRN(NZ),RHOP(NZ)
      DIMENSION TAURAN(NZ,NP),ALAM(NZ),TAUCPK(NZ,NP),ETA(NZ)
      character*8 PLANET,ISPEC
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PHOTABLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/CBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/GBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/NBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/AERBLK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/ISOBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc'
      DIMENSION CUNING(NZ,NP),amass(NZ,NP),THERMSP(NZ,NP),
     2  TAURELAXC(NZ,NP),TAURELAX(NZ,NP),AFPL(NZ,NP),delta(NZ,NP),
     3  BETAF(NZ,NP)
 
      REAL NMON,RMON,DF  !frachack 


C
C   THIS SUBROUTINE CALCULATES FALL VELOCITIES AND ESTIMATES PARTICLE SIZE
C   BASED ON THEIR COAGULATION LIFETIMES
C
C   CONSTANTS from Pruppacher and Klett page 450
      A = 1.257
      B = 0.4
      C = 1.1

      BK = 1.38*1.E-16
      PI = 3.14159
      NZ1 = NZ - 1

      IF (frak.eq.1) THEN
C **********************************************
C-EW  IMPLEMENTATION OF FRACTAL MICROPHYSICS
C-EW  ONLY FOR HYDROCARBONS, K=3 and 4
C
C-EW  IMPLEMENTATION USES A SIZE BIN DEPENDENT FRACTAL
C-EW  DIMENSION TO PARAMETERIZE AGGREGATE RESTRUCTURING 
C-EW  PERMEABILITY EFFECTS ARE NOT TREATED HERE
C
C-EW RPAR = equal mass spherical radii
C-EW RFRAC = fractal aggregate radii
C
C-EW  THIS IS THE MONOMER RADIUS [cm] 
C-GA aak! didn't realize this was hardcoded here! Changing now
C-   that we have the monsize keyword in input_photochem.dat
         IF (monsize .eq. 0)  RMON = 50.E-7 
         IF (monsize .eq. 1)  RMON = 10.E-7 
         IF (monsize .eq. 2)  RMON = 20.E-7 
         IF (monsize .eq. 3)  RMON = 70.E-7 
         IF (monsize .eq. 4)  RMON = 10.E-6 

         DO K=3,4               !hcaer1 and hcaer2
           DO J=1,NZ
             NMON = (RPAR(J,K)/RMON)**3.
             IF (NMON .LE. 1.) THEN
             !i.e. of number of monomers is <= 1, it's a 
             !spherical particle. By definition, DF for
             !spheres is 3.
               DF = 3.
            ELSE
             !calculate DF (fractal param) for nonspheres 
               DF = 2.4 - 0.9*EXP(-NMON/500.)
            ENDIF
           !calculates fractal radius based on 
           !equivalent mass spherial particles
           RFRAC(J,K) = RPAR(J,K)**(3./DF)*RMON**(1.-3./DF)
          ENDDO
        ENDDO


C-EW *******************************************   
C
      ENDIF

      DO J=1,NZ
        ALAM(J) = 1.63E14/DEN(J) !Eddy diffusion stuff?
        ETA(J)=ABS((1.718 + 0.0049*(T(J)-273.) - 1.2
     &  *(1.E-5)*(T(J)-273.)*(T(J)-273.))*1.E-4)
      ENDDO 


      DO 10 K=1,NP
C   (1 = SULFATE, 2 = S8, 3 = HCAER1, 4=HCAER2)

        if (ISOTOPE.EQ.0) L = LSO4AER   !so using the sulfate aersol for all particles?
        if (ISOTOPE.EQ.1) L = LSXO4AER  !so using the sulfate aersol for all particles?
       
        if(usetd.eq.1) then !Jim's code with particles in tri-diag uses H2SO4 rather than the aerosol
          if (ISOTOPE.EQ.0) L = LH2SO4  !using h2so4 for all particles? (sure why not - it's infinite rainout...)
          if (ISOTOPE.EQ.1) L = LH2SXO4  !using h2so4 for all particles?
        endif 

        if (ISOTOPE.EQ.0) then   !ISOHACK 


C-AP ESTIMATION OF THE AEROSOL FREE PATH LENGTH
        DO J = 1,NZ
          IF(K.GE.3 .AND. frak.eq.1) THEN  !ACK hardcoded particle number (getting HCAER and HCAER2)
C-EW  HYDROCARBONS: USE FRACTAL MICROPHYSICS
            ALPH = A + B*EXP(-C*RFRAC(J,K)/ALAM(J)) !giada - related to particle diffusion 
            CUNING(J,K) = 1 + ALPH*ALAM(J)/RFRAC(J,K)

C-GA: density of hc aerosols should be 0.63 g/cm3
C        see Trainer et al (2006)
C        --> density is now read in from input_photochem.dat
C-AP Notation is similar Fusch 1964
C THERMSP = thermal velocity of molecule

            amass(J,K) = (4./3.)*PI*RPAR(J,K)**3*HCDENS !particle mass
            THERMSP(J,K) = SQRT((8*BK*T(J))/(pi*amass(J,K))) !thermal velocity
            TAURELAXC(J,K)=2*RPAR(J,K)**3/(9*ETA(J)*RFRAC(J,K)) !relaxation time
            TAURELAX(J,K) = TAURELAXC(J,K)*CUNING(J,K) 
            AFPL(J,K) = THERMSP(J,K)*TAURELAX(J,K)           !related to particle mean free path
            delta(J,K) = (((2*RFRAC(J,K)+AFPL(J,K))**3 - (4* !delta is mean distance from ctr of sphere
     &          RFRAC(J,K)*RFRAC(J,K)+AFPL(J,K)*AFPL(J,K))**1.5)/ !reached by particle leaving the surface and traveling distance 
     &         (6*RFRAC(J,K)*AFPL(J,K)) - 2*RFRAC(J,K))*SQRT(2.) !equal to the mean free path (AFPL)
           ELSE
C-EW  S8,SO4,HC if frak=0: USE SPHERICAL MICROPHYSICS
             ALPH = A + B*EXP(-C*RPAR(J,K)/ALAM(J))
             CUNING(J,K) = 1 + ALPH*ALAM(J)/RPAR(J,K)

C-AP Notation is similar Fusch 1964
             adensity = HCDENS
             amass(J,K) = (4./3.)*PI*RPAR(J,K)**3*adensity !aerosol mass
             THERMSP(J,K) = SQRT((8*BK*T(J))/(pi*amass(J,K))) !thermal velocity
             TAURELAXC(J,K)=2*RPAR(J,K)*RPAR(J,K)/(9*ETA(J))  !relaxation timescale
             TAURELAX(J,K) = TAURELAXC(J,K)*CUNING(J,K) 
             AFPL(J,K) = THERMSP(J,K)*TAURELAX(J,K)            !related to MFP
             delta(J,K) = (((2*RPAR(J,K)+AFPL(J,K))**3 - (4*   !mean distance from center of sphere
     &            RPAR(J,K)*RPAR(J,K)+AFPL(J,K)*AFPL(J,K))**1.5)/
     &            (6*RPAR(J,K)*AFPL(J,K)) - 2*RPAR(J,K))*SQRT(2.) 
          ENDIF  !end frak loop
       ENDDO

C-AP Calculation of the correction to the coagulation kernel 
       DO J = 1,NZ
          IF(K.GE.3 .AND. frak.eq.1) THEN !get HCAER and HCAER2
C-EW  HYDROCARBONS: USE FRACTAL MICROPHYSICS 
            BETAF(J,K) = 1/(RFRAC(J,K)/(RFRAC(J,K)+delta(J,K)/2)
     &      +PI*AFPL(J,K)/(2*SQRT(2.)*RFRAC(J,K)))      
         ELSE
C-EW  S8,SO4,HC if frak=0: USE SPHERICAL MICROPHYSICS 
            BETAF(J,K) = 1/(RPAR(J,K)/(RPAR(J,K)+delta(J,K)/2)
     &      +PI*AFPL(J,K)/(2*SQRT(2.)*RPAR(J,K)))
         ENDIF
      ENDDO
C-AP ******************************************************  

C   ESTIMATE COAGULATION AND SEDIMENTATION LIFETIMES (TOON AND FARLOW, 1981)
      DO I=1,NZ
c-GA  this commented one was from Toon and Farlow 1981
c-GA  the new methodology below from Pavlov 2001
c     TAUC(I,K) = 1.E6/(AERSOL(I,K)*SQRT(RPAR(I,K)))      !T&F (81) p.41 - e-folding lifetime against coagulation
                                                          ! 1/tauc = 1/N DN/DT
c recompute tauc using new methodology 
         TAUCPK(I,K) = 3*ETA(I)/(4*AERSOL(I,K)*BK*
     &        T(I)*CUNING(I,K))
         TAUC(I,K)=TAUCPK(I,K)/BETAF(I,K) !coagulation lifetime
         TAURAN(I,K) = 1./(RAINGC(L,I) + 1.E-20) !rainout lifetime
C-GA Gotta love comments like the ones below (saving it for amusement and also b/c it may be important)...
! this is really wrong for S8 !!  !where does this come from?
! it looks like I never completly fixed the rainout thing.  TAURAN is using LSO4AER for both.  but isn't RAINGC really small for S8?
                                     
         TAUSED(I,K) = HSCALE(I)/WFALL(I,K) ! sedimentation lifetime
      ENDDO

C   FIND MINIMUM OF DIFFUSION AND SEDIMENTATION LIFETIMES, THEN SCALE PRTICLE SIZES
      e_minus_5 = 1.E-5
      e_hc=1.3E-7
      DO  I=1,NZ
         TAUTRN(I) = min(TAUSED(I,K),TAUEDD(I)) !find the minimum of the three destruction timescales 
         TAUTRN(I) = min(TAUTRN(I),TAURAN(I,K)) !lifetime against eddy diffusion is H*H/K 
                                                !where K is eddy diffusion coefficient, H is scale height
         RPAR(I,K) = RPAR(I,K) * (TAUTRN(I)/TAUC(I,K))**0.25 !grow the particle
         if (K.GE.3) then
            RPAR(I,K) = max(RPAR(I,K),e_hc)      !largest HC particles are smaller?
         else
            RPAR(I,K) = max(RPAR(I,K),e_minus_5) !largest particles are 1 micron
         endif
      ENDDO

C   DON'T ALLOW PARTICLES TO DECREASE IN SIZE AT LOW ALTITUDES
      DO 3 I=1,NZ1
         J = NZ - I
 3       RPAR(J,K) = max(RPAR(J,K),RPAR(J+1,K))

      endif  !END OF ISOTOPE SKIP LOOP
C
C   COMPUTE PARTICLE-TO-GAS CONVERSION FACTORS AND DENSITIES
      do I=1,nz
        R = RPAR(I,K)
        LL=NQ-NP+K  !ACK - assuming particles are last LL elements
C-GA: the above hardcoding deserves a comment in the inputs
C-GA: readme file, which I've added.
        if (USETD.EQ.1) LL=K+NQ !particles in tri-diag

        if (LL.EQ.LS8AER.OR.LL.EQ.LSXS7AER) then
           RHOP(I) = 2.07!density
           factor=2.03E7 !molecules per 0.1 um sphere 
        endif

        if (LL.EQ.LSO4AER.OR.LL.EQ.LSXO4AER) then
           RHOP(I) = 1. + 0.8*FSULF(I)!density
           factor = 4.6E7*FSULF(I)    !molecules per 0.1 um sphere 
        endif

       if (LL.EQ.LHCAER .OR. LL.EQ.LHCAER2) then
          RHOP(I) = HCDENS                     !hcaer density
          gpermolec = 2*1.66e-24+12*4*1.66e-24 !grams per molecule for HCAER
          factor=(4./3.)*PI*(1.E-5)**3*HCDENS/gpermolec !molecules per 0.1 um sphere 
       endif
          CONVER(I,K) = factor * (R/1.E-5)**3 !number of molecules/particle
                                              !used in main code to calculate aersol (num dens of aerosols)
C - GA: keeping this old comment because it's hilarious and caused hours of difficulties trying
C       to figure out what 'factor' meant
!"factor is related to density somehow (jim's words)"
!wow...thanks for that helpfulness :-P

       enddo

       if(ISOTOPE.EQ.0) then    !fall velocities should stay the same in isotope code
C   NOW COMPUTE FALL VELOCITIES
          DO  J=1,NZ
      
             IF(K.GE.3 .AND. frak.eq.1) THEN  
C-EW  HYDROCARBONS: USE FRACTAL MICROPHYSICS 
                R = RPAR(J,K)
                RF = RFRAC(J,K) 
                F1 = 2./9. * RHOP(J)*R*R*R*G/ETA(J)/RF !from stokes law F1 - settling velocity
C - GA it's computing terminal velocity  when frictional and buoyant forces
C- are equal to gravitational force
                ALPH = A + B*EXP(-C*RF/ALAM(J)) ! related to particle resistance to motion 
                WFALL(J,K) = F1*(1. + ALAM(J)*ALPH/RF) !wfall = fall velocity
                             !this term (alam*alph/rf) is particle diffusion?  I think -gna
             ELSE      
C-EW  S8, SO4,HC if frak=0: USE SPHERICAL MICROPHYSICS
                R = RPAR(J,K) 
                F1 = 2./9. * RHOP(J)*R*R*G/ETA(J)
                ALPH = A + B*EXP(-C*R/ALAM(J))
                WFALL(J,K) = F1*(1. + ALAM(J)*ALPH/R)
             ENDIF
          enddo

       endif                    !end ISOTOPE skip loop 

 10   CONTINUE


      RETURN
      END
