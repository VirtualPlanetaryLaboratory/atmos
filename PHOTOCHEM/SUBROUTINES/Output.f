      SUBROUTINE OUTPUT(N,NSTEPS,TIME,jtrop,vdep,USOLORIG,USETD,frak)
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      real*8 mass
      CHARACTER*8 ISPEC,REACTYPE,PLANET,CHEMJ
      CHARACTER*20 fmtstr,fmtdatastr
      CHARACTER*60 fmtheadstr

      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PHOTABLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/BBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/CBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/FBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/GBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/JBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/NBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/QBLOK.inc'  !this has flux in it - compare against below
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/RBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/SBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/WBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/ZBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/SATBLK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/SULBLK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/AERBLK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/RRATS.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/MBLOK.inc'

      COMMON/LifeTime/TAUO2,TAUCH4,TAUSO2

      DIMENSION FUP(NQ1),FLOW(NQ1),CON(NQ1),FLUXO(NQ1,NZ)
     2  ,ZF(NZ),TAUAER(NP+1)
      dimension TAUHCABS(kw),TAUHCEXT(kw), TAUHCEFF(kw)
      dimension vdep(nq), flux_H(nz),USOLORIG(NQ,NZ)
      dimension PSO2MIF(nz),PROSO2MIF(nz),PROSO(nz)


C   THIS SUBROUTINE PRINTS OUT ALL THE DATA.  THE VARIABLE ISKIP SAYS
C   HOW MANY GRID POINTS YOU WANT TO LOOK AT.
C
      ISKIP = 4
      JSKIP = ISKIP
      IF(N.EQ.NSTEPS) ISKIP = 2
      TIMEY = TIME/3600./24./365.
      write(14, 100) TIME,TIMEY

      write(14, 101) NPHOT
 !100  format(/1X,'TIME =',E11.4,5X,'TIMEY =',F10.5,1X,'YEARS')
 !101  format(/1X,'NPHOT =',I3)
c-mab above are old format, using the format below to avoid **** in out.out
 100  FORMAT(/1X,"TIME =",1PE9.2,5X,"TIMEY =",E9.2,1X,"YEARS")
 101  format(/1X,'NPHOT =',I5)

C
      write(14, 105)
 105  format(/1X,'MIXING RATIOS OF LONG-LIVED SPECIES'/)
      IROW = 18
      LR = NQ/IROW + 1
      RL = FLOAT(NQ)/IROW + 1
      DIF = RL - LR
      IF (DIF.LT.0.001) LR = LR - 1
C
      DO 8 L=1,LR
      K1 = 1 + (L-1)*IROW
      K2 = K1 + IROW - 1
      IF (L.EQ.LR) K2 = NQ
      write(14, 110) (ISPEC(K),K=K1,K2)
 110  format(/5X,'Z',8X,18(A8,1X))
      DO I=1,3
      write(14, 120) Z(I),(USOL(K,I),K=K1,K2)
      END DO
      DO I=4,NZ,ISKIP
      write(14, 120) Z(I),(USOL(K,I),K=K1,K2)
      END DO
 120  format(1X,1P19E9.2)
      IF (N.EQ.0) GO TO 8 !on first call don't print TP, TL
      write(14, 140)
 140  format(/1X,'TP, TL')
      write(14, 145) (TP(K),K=K1,K2)
      write(14, 145) (TL(K),K=K1,K2)
 145  format(10X,1P18E9.2)
   8  CONTINUE
C
c      write(14, 106)
c106  format(//1X,'S8AER MIXING RATIO'/)
c      write(14, 180) S8AER

      IF (N.EQ.0) RETURN
C     EVERYTHING BELOW PRINTED AFTER FINAL TIMESTEP

      write(14, 107) TP(LS8AER),TL(LS8AER)
 107  format(/5X,'TP =',1PE10.3,2X,'TL =',E10.3)


      write(14, 150) O3COL
 150  format(//1X,'OZONE COLUMN DEPTH = ',1PE11.4)
      write(19, 150) O3COL   !terse output
      CO_column = 0.0
      CH4_column = 0.0
      O2_column = 0.0
      do i=1,nz
        CO_column = dz(i)*usol(LCO,i)*den(i) + CO_column
        CH4_column = dz(i)*usol(LCH4,i)*den(i) + CH4_column
        O2_column = dz(i)*usol(LO2,i)*den(i) + O2_column
      enddo
      write(14, 159) O2_column, CO_column, CH4_column   !terse output
      write(19, 159) O2_column, CO_column, CH4_column   !terse output
 159  format(//1X,'O2, CO, CH4 COLUMN DEPTHs = ',1P3E14.3)
      write(14, 152) H2SCOL,SO2COL,S2COL,S4COL,S8COL
 152  format(/1X,'SULFUR COLUMN DEPTHS:  H2S =',1PE10.3,2X,'SO2 =',
     2  E10.3,2X,'S2 =',E10.3,2X,'S4 =',E10.3,2X,'S8 =',E10.3)


c-mab FOR THE FIRST TWO PARTICLES (ASSUMED SULPHATE AND S8 IN HARDCODING ORDER)
      IF (NP.GT.0) THEN
       DO J=1,NP
!ack hardcoding to retain Kevin's original scheme
!where tauaer 2 and 3 were both for elemental sulfur
!this should be integrated/updated below, mab: just did that below (after defining TAUAER with NP+1 above)
        IF (J.NE.3) TAUAER(J) = 0.
c-mab: hardcoding, designed to keep tauaer 3 value from previous J = 2 step
corig      TAUAER(NP+1) = 0.
        DO I=1,NZ
         R = RPAR(I,J)
!this should be the extinction efficiency at some wavelength - presumably visible? why 0.6?
!c-mab reconfigured this loop cause there was a bug previously...
         IF (J.LT.3) TAUAER(J) = TAUAER(J) +
     &                           0.6*3.14159*R*R*AERSOL(I,J)*DZ(I) !S8 VIS (0.6)
         IF (J.EQ.2) TAUAER(J+1) = TAUAER(J+1) +  !c-mab Hold S8(UV) in tauaer3
     2     1.2*3.14159*R*R*AERSOL(I,J)*DZ(I)
!so 1.2 is somehow Qext of S8 in the UV? - don't we have these all as 2?
         IF (J.EQ.3) TAUAER(J) = TAUAER(J)
!c-mab Double check to make sure tauaer3 doesn't get re-written by 0.0 or something
        ENDDO
       ENDDO
        write(14, 153) (TAUAER(J),J=1,3)
!ACK hardcoded to 2 particles + S8UV (mab note: S8UV held in next, i.e. 3rd level for now)
c-mab: IF TAUAER is used for addition particles, increment J by 1, i.e. particle 3 would go to J=4
 153    format(/1X,'SCALED AEROSOL EXTINCTION OPTICAL DEPTHS',/5X,
     2  'SULFATE =',1PE10.3,5X,'S8(VIS) =',E10.3,5X,'S8(UV) =',E10.3)
C
      ENDIF

!mc perhaps need to abstract the above with better wavelength behavior.
!am porting over some HCAER stuff from Shawn's code:

      do i=1,kw
      TAUHCABS(i) = 0.
      TAUHCEXT(i) = 0.
      TAUHCEFF(i) = 0.
c      print *, W0T(j,1)
      enddo
c      stop

      if (NP.gt.2) then
         L3 = 3 !EWS - prevents gfortran compilation warnings
      do L=L3,NP  !ACK - hardcoded for HCAER/HCAER2 to be the final particle
      do j=1,kw
      do I=1,NZ
      R = RPAR(I,L)  !hcaer  (ACK - should be an if here)
      extinction=QEXTT(j,I,L)*3.14159*R*R*AERSOL(I,L)*DZ(I)
      TAUHCABS(j) = TAUHCABS(j) + extinction*(1-W0T(j,I,L))
      TAUHCEXT(j) = TAUHCEXT(j) + extinction
      TAUHCEFF(j) = TAUHCEFF(j) + extinction *
     $              (1-W0T(j,I,L)/2.-GFT(j,I,L)*W0T(j,I,L)/2.)

c      if (I.EQ.1.and.j.lt.118) print *, L,j,GFT(J,I,L),W0T(J,I,L),
c     $   1-W0T(j,I,L)/2.,-GFT(j,I,L)*W0T(j,I,L)/2.,
c     $              (1-W0T(j,I,L)/2.-GFT(j,I,L)*W0T(j,I,L)/2.)
      enddo
      enddo
      enddo
c      stop

      i190=minloc(wavl,1,wavl.ge.1900)
c      i190=minloc(wavl,1,wavl.ge.1750)  !test
      i200=minloc(wavl,1,wavl.ge.2000)
      i210=minloc(wavl,1,wavl.ge.2100)
      i220=minloc(wavl,1,wavl.ge.2200)
      i258=minloc(wavl,1,wavl.ge.2580)
      i550=minloc(wavl,1,wavl.ge.5500)

      PRINT 353, TAUHCABS(i190),TAUHCABS(i200),TAUHCABS(i210),
     $  TAUHCABS(i220),TAUHCABS(i258),TAUHCABS(i550)

      write(14,353) TAUHCABS(i190),TAUHCABS(i200),TAUHCABS(i210),
     $  TAUHCABS(i220),TAUHCABS(i258),TAUHCABS(i550)

      write(63,*) TAUHCABS(i190),TAUHCABS(i200),TAUHCABS(i210),
     $  TAUHCABS(i220),TAUHCABS(i258),TAUHCABS(i550)

 353  FORMAT(/1X,"SCALED HC AEROSOL ABSORPTION OPTICAL DEPTHS",/5X,
     2  "190 nm =",1PE10.3,5X,"200 nm =",E10.3,5X,"210 nm =",
     3  E10.3,5X,"220 nm =",E10.3,5X,"258 nm =",E10.3,5X,"550 nm ="
     4  ,E10.3)

      PRINT 253, TAUHCEXT(i190),TAUHCEXT(i200),TAUHCEXT(i210),
     $  TAUHCEXT(i220),TAUHCEXT(i258),TAUHCEXT(i550)

      write (14,253) TAUHCEXT(i190),TAUHCEXT(i200),TAUHCEXT(i210),
     $  TAUHCEXT(i220),TAUHCEXT(i258),TAUHCEXT(i550)

      write (63,*) TAUHCEXT(i190),TAUHCEXT(i200),TAUHCEXT(i210),
     $  TAUHCEXT(i220),TAUHCEXT(i258),TAUHCEXT(i550)

 253  FORMAT(/1X,"SCALED HC AEROSOL EXTINCTION OPTICAL DEPTHS",/5X,
     2  "190 nm =",1PE10.3,5X,"200 nm =",E10.3,5X,"210 nm =",
     3  E10.3,5X,"220 nm =",E10.3,5X,"258 nm =",E10.3,5X,"550 nm ="
     4  ,E10.3)


      PRINT 453, TAUHCEFF(i190),TAUHCEFF(i200),TAUHCEFF(i210),
     $  TAUHCEFF(i220),TAUHCEFF(i258),TAUHCEFF(i550)

      write (14,453) TAUHCEFF(i190),TAUHCEFF(i200),TAUHCEFF(i210),
     $  TAUHCEFF(i220),TAUHCEFF(i258),TAUHCEFF(i550)

      write (63,*) TAUHCEFF(i190),TAUHCEFF(i200),TAUHCEFF(i210),
     $  TAUHCEFF(i220),TAUHCEFF(i258),TAUHCEFF(i550)


 453  FORMAT(/1X,"SCALED HC AEROSOL EFFECTIVE OPTICAL DEPTHS",/5X,
     2  "190 nm =",1PE10.3,5X,"200 nm =",E10.3,5X,"210 nm =",
     3  E10.3,5X,"220 nm =",E10.3,5X,"258 nm =",E10.3,5X,"550 nm ="
     4  ,E10.3)



!shawn created an outputfile for the extinction optical depths at 190,220,550
!this could be useful, but check for existing out.tau - although it probably doesn't do particles right...


c      ISKIP=10  !don't know if we really care about this, but porting for now

c      PRINT *, 'QEXTUV      OMG0A       ALT'
c      DO I=1,NZ,ISKIP
c       PRINT *, QEXTT(i258,I), W0T(i258,I), Z(I)
c      ENDDO
c      PRINT *, 'QEXTVIS      OMG0A       ALT'
c      DO I=1,NZ,ISKIP
c       PRINT *, QEXTT(i550,I), W0T(i550,I), Z(I)
c      ENDDO


      endif  !loop for particles 3 and 4


      write(14, 151) USOL(LH2O,NH)
 151  format(/1X,'FH2O AT COLD TRAP =',1PE10.3)


c in order to print fluxes I need to do this next bit
c     IF(N.LT.NSTEPS) RETURN
C
C ***** write on LAST ITERATION ONLY *****
c-mc      DO 1 I=1,NZ
c-mc   1  ZF(I) = Z(I) + 0.5*DZ    ! ZF(1) = 1.0 km  !this is just used for printing
c ACK this might need a little tire kicking as I go to a variable grid
      do i=1,nz
       ZF(i) = Z(i) + DZ(i)/2    !this is just used for printing
      enddo


C
c  what follows gives the fluxes
      DO 3 K=1,NQ
      DO I=1,NZ
      SL(K,I) = USOL(K,I)*DEN(I)
      END DO
      DO I=1,NZ-1
      FLUXO(K,I) = - DK(I)*(USOL(K,I+1) - USOL(K,I))/DZ(I)
      END DO
   3  CONTINUE
!diffusion added below

      vturb=0.01
c      vturb=0.0

c fluxes for particles
      if(NP.GT.0) THEN
      do k=1,np
       do I=1,NZ1

       if (USETD.EQ.1) then
       if (I.EQ.1) then
       FLUXO(K+NQ,I) = - DK(I)*(PARTICLES(I+1,K) - PARTICLES(I,K))/DZ(I)
     1  - 0.5*((WFALL(I,K)+vturb)*DEN(I)*PARTICLES(I,K)
     2          + WFALL(I+1,K)*DEN(I+1)*PARTICLES(I+1,K))
!turbulent diffusion velocity of 0.01 at the lower boundary (Pavlov 02 appendix)
!hmm. maybe the wfall(I+1) shouldn't be here... - or vturb should - but the fluxes for LL don't seem to know about the boundary conditions. maybe this should just be transport...
       else

       FLUXO(K+NQ,I) = - DK(I)*(PARTICLES(I+1,K) - PARTICLES(I,K))/DZ(I)
     1  - 0.5*(WFALL(I,K)*DEN(I)*PARTICLES(I,K)
     2          + WFALL(I+1,K)*DEN(I+1)*PARTICLES(I+1,K))

      endif
      else  !particles are in the main matrix
         !for now, assume they are in the last two LL spots - should probably abstract this...
         LL=NQ-NP+K

      WNF = 0.5*(WFALL(I,K)*DEN(I)*USOL(LL,I) + WFALL(I+1,K)*DEN(I+1)
     2  *USOL(LL,I+1))
      FLUXO(LL,I) = FLUXO(LL,I) - WNF
      endif


       enddo
      enddo
      endif
c add in molecular diffusion

c      print *, FLUX(
C
c         \phi = b*f*({1/over H_A}-{1/over H}) - b*df/dz
c     do k=1,NQ
c       do i=1,NZ1
c       flux(k,i) = flux(k,i)
c    5    - bXN2(i)*(usol(k,i+1) - usol(k,i))/dz
c    6    + bXN2(i)*usol(k,i)*(1./H_atm(i) - 1./scale_H(k,i))
c       enddo
c     enddo

c these should be flagged for IH2=-1
c         \phi = b*f*({1/over H_A}-{1/over H}) - b*df/dz
c the following gives fluxes from 1 to 80 km.

       do i=1,NZ-1
        do j=1,(NQ-NP)
         fluxo(j,i) = fluxo(j,i)
     5    - bX1X2(j,i)*(usol(j,i+1) - usol(j,i))/dz(i)
     6    + bX1X2(j,i)*0.5*(usol(j,i) + usol(j,i+1))
     7       *(1./H_atm(i) - 1./scale_H(j,i))
         enddo
       enddo

C      else
C       do i=1,NZ-1
C        fluxo(LH,i) = fluxo(LH,i)
C     5    - bHN2(i)*(usol(LH,i+1) - usol(LH,i))/dz(i)
C     6    + bHN2(i)*0.5*(usol(LH,i) + usol(LH,i+1))
C     7       *(1./H_atm(i) - 1./scale_H(LH,i))
C        fluxo(LH2,i) =fluxo(LH2,i)
C     5    - bH2N2(i)*(usol(LH2,i+1) - usol(LH2,i))/dz(i)
C     6    + bH2N2(i)*0.5*(usol(LH2,i) + usol(LH2,i+1))
C     7      *(1./H_atm(i) - 1./scale_H(LH2,i))
C       enddo
C      endif


C
c
c the following is self-consistent? molecular diffusion is in FLUX(K,1)
c  YP(1) and YL(1)*SL(1) are chemical production and loss in lowest kilometer (number densities)/sec
c  but its not entirely clear that the boundary condition knows this
      DO 15 K=1,NQ1
      FLOW(K) = FLUXO(K,1) - (YP(K,1) - YL(K,1)*SL(K,1))*DZ(1)
      FUP(K) = FLUXO(K,NZ1) + (YP(K,NZ) - YL(K,NZ)*SL(K,NZ))*DZ(NZ)
      CON(K) = TP(K) - TL(K) + FLOW(K) - FUP(K)
  15  CONTINUE
C

c water is not conserved.  should this be jtrop or jtrop + 1?
      FLOW(LH2O) = FLUXO(LH2O,jtrop)   ! jim had 11 hard-wired
      CON(LH2O) = TP(LH2O) - TL(LH2O) + FLOW(LH2O) - FUP(LH2O)
c-mab disabling FLUXO of water assignment to 0.0 for giants based on WT.
c-mab This is something I'm trying to move from a planet-based if loop.
c-mab Most solar system gas giants have a mwt of <3.0 and terrestrials > 20.0
c-mab Using wt of methane - 16.0 as a distinction here
c-mab (May change later after input soliciation)
      IF (WT.GT.16.0) THEN
      DO I=1,jtrop-1   ! jim had 10 hard wired... its not zero
      FLUXO(LH2O,I) = 0.
      END DO
      ENDIF
c there are issues here -


c  ok.  I've got the fluxes.

c  there is another i hope equal measure of H, H2 escape
c   experiment indicates that this is the same as what's in the model
      F_esc_H = usol(LH,NZ)*bHN2(NZ)*(1./H_atm(nz) - 1./scale_H(LH,nz))
      F_esc_H2 = usol(LH2,NZ)*bH2N2(NZ)
     $ *(1./H_atm(nz) - 1./scale_H(LH2,nz))
C
c  sum up the H fluxes at all heights except the ground
c  this works good for stratosphere but in the troposphere
c it fails -
c-mab abstracted by looping over species where atomsH ne 0
      flux_H = 0.0
      do i=1,NZ-1
        do j=1,NQ1   !note
          if(atomsH(j).gt.0.0) then
c-mab then only sum for species that contain H for a given layer
        	flux_H(i) = flux_H(i) + atomsH(j)*FLUXO(j,i)
          endif
        enddo
      enddo
        flux_H(nz) = FUP(LH) + 2.*FUP(LH2)

      write(24,974) time,
     $  usol(LO2,1), usol(LH2,1), usol(LCO,1),usol(LCH4,1),
     1  FLOW(LO2), FLOW(LCH4), FLOW(LH2), FLOW(LCO), distflux(LSO2),
     2  FLOW(LH2S), FUP(LH2), FUP(LH), TLOSS(LH2CO), TLOSS(LH2O2),
     3  TLOSS(LH2SO4)+TLOSS(LSO4aer), TLOSS(LSO2), TLOSS(LS8aer)
     4  ,O3COL,SL(NSP-1,1)/den(1), distflux(LH2), distflux(LCO),
     5  distflux(LO2), CO_column
  974 format(1x,1PE13.6,23E10.2)
c  I want the averages of these things




      IF(N.LT.NSTEPS) RETURN
C
C ***** write on LAST ITERATION ONLY *****



      write(14, 125)
 125  format(/1X,'NUMBER DENSITIES OF LONG-LIVED SPECIES'/)
      DO 9 L=1,LR
      K1 = 1 + (L-1)*IROW
      K2 = K1 + IROW - 1
      IF (L.EQ.LR) K2 = NQ
      write(14, 110) (ISPEC(K),K=K1,K2)
      DO I=1,NZ,ISKIP
      write(14, 120) Z(I),(SL(K,I),K=K1,K2)
      END DO
   9  CONTINUE
C
c-mc
c-mc want to write out number densities of SO2, O2,H20,CO2

         L1=LSO2
         L2=LS8
         if(L2.eq.0) L2=LS8AER



      DO I=1,NZ
      write(27, 224) SL(L1,I),SL(LO2,I),SL(LH2O,I),SL(LCO2,I),
     2 SL(L2,I)
      enddo

 224  format(6(1pe10.3,2X))


c-mc writing out final number densities to out.finalden

      write(42,112) (ISPEC(K),K=1,NQ)
       do I=1,nz
        write(42,114) (SL(K,I),K=1,NQ)
       enddo


 112  format(100(1X,A8,1X))   !ACK update if NQ ever goes > 100
 114  format(100(1pe10.3))    !ditto
c-mc

      write(14, 155)
 155  FORMAT(/1X,'FLUXES OF LONG-LIVED SPECIES'/)
      ZFL = 0.
      ZFT = ZF(NZ)
      DO 10 L=1,LR
      K1 = 1 + (L-1)*IROW
      K2 = K1 + IROW - 1
      IF (L.EQ.LR) K2 = NQ1
      write(14, 510) (ISPEC(K),K=K1,K2)
 510  format(/5X,'Z',8X,18(A8,2X))
 520  format(1X,1P1E9.2,1P18E10.2)
      write(14, 520) ZFL,(FLOW(K),K=K1,K2)
      DO I=1,NZ,ISKIP
      write(14, 520) ZF(I),(FLUXO(K,I),K=K1,K2)
      END DO
      write(14, 520) ZFT,(FUP(K),K=K1,K2)
  10  CONTINUE
C

      write(45,*) ZFL, FLOW(2)
       do I=1,NZ
        write(45, *) ZF(I),FLUXO(2,I)
       enddo
      write(45, *) ZFT,FUP(2)



c-mc write out total production and loss terms
c-mc this needs to be fixed for the generic case.  should just write out all sulfur gases, I suppose.
c-mc or nest some ifs somehow.  what a pain. for now, just hacking S8's out.  this will probably break the analysis scripts...

c$$$      write(27,533) (YP(LS8,I), YL(LS8,I),YP(LH2SO4,I),
c$$$     2    YP(LS,I),YP(LS2,I),YP(LSO4AER,I),YP(LS8AER,I),
c$$$     3    YP(LSO2,I),YL(LSO2,I), I=1,NZ)

c$$$ 533  format(9(1PE9.2,2X))
c above are origs.  below is hacked for my S8 gas removal

      write(27,533) (YP(LH2SO4,I),
     2    YP(LS,I),YP(LS2,I),YP(LSO4AER,I),YP(LS8AER,I),
     3    YP(LSO2,I),YL(LSO2,I), I=1,NZ)

 533  format(7(1PE9.2,2X))





c-copying tp/tpl loop here to help with print out. checking units..
c - find the sulfur production plots as a function of height...
      do j=1,NZ
       write(44,534) (YP(I,J), I=1,NQ1)
      enddo
      do j=1,NZ
       write(44,534) (YL(I,J)*SL(I,J), I=1,NQ1)  !'loss' units in .prod file in 1/cm^3/s to comapre with 'production'(yp)
      enddo
 534  format(1P100E10.3)    !ACK - update if NQ>100

c-when I use this format, IDL reads in <10^-72 as 0 which is fine by me)


c-mc


      write(14, 205)
 205  format(/1X,'AQUEOUS PHASE SPECIES'/)
      write(14, 206)
 206  format(5X,'Z',6X,'(SO2)g',3X,'(H2CO)g',2X,'(SO2)aq',1X,
     2  'CH2(OH)2',3X,'HCO3-',4X,'CO3=',5X,'HSO3-',4X,'SO3=',3X,
     3  'CH2OHSO3-',3X,'OH-',5X,'SO4=',6X,'PH')
      DO I=1,NH,ISKIP
      write(14, 120) Z(I),(XSAVE(K,I),K=1,10),SO4SAV(I),PH(I)
      END DO
C
      write(14, 210)
 210  format(/1X,'NORMAL HENRYS LAW COEFFICIENTS'/)
      DO 25 L=1,LR
      K1 = 1 + (L-1)*IROW
      K2 = K1 + IROW - 1
      IF (L.EQ.LR) K2 = NQ
      write(14, 110) (ISPEC(K),K=K1,K2)
      DO I=1,NH,ISKIP
      write(14, 120) Z(I),(H(K,I),K=K1,K2)
      END DO
  25  CONTINUE
C
      write(14, 215)
 215  format(/1X,'ENHANCEMENTS OF HENRYS LAW COEFFICIENTS'/)
      DO 27 L=1,LR
      K1 = 1 + (L-1)*IROW
      K2 = K1 + IROW - 1
      IF (L.EQ.LR) K2 = NQ
      write(14, 110) (ISPEC(K),K=K1,K2)
      DO I=1,NH,ISKIP
      write(14, 120) Z(I),(ENHAN(K,I),K=K1,K2)
      END DO
  27  CONTINUE
C
      write(14, 220)
 220  format(/1X,'GIORGI AND CHAMEIDES RAINOUT RATES'/)
      DO 29 L=1,LR
      K1 = 1 + (L-1)*IROW
      K2 = K1 + IROW - 1
      IF (L.EQ.LR) K2 = NQ
      write(14, 110) (ISPEC(K),K=K1,K2)
      DO I=1,NH,ISKIP
      write(14, 120) Z(I),(RAINGC(K,I),K=K1,K2)
      END DO
  29  CONTINUE
C
      write(14, 175)
  175 format(/1X,'RAINOUT RATE, PHIDEP, TLOSS, LOWER B.C. and vdep'/)
      write(14, 176)
 176  format(1X,'FOLLOWED BY TP, TL, FUP, FLOW, CON'/)
      DO 13 L=1,LR
      K1 = 1 + (L-1)*IROW
      K2 = K1 + IROW - 1
      IF (L.EQ.LR) K2 = NQ
      write(14, 110) (ISPEC(K),K=K1,K2)
      write(14, 145) (SR(K),K=K1,K2)
      write(14, 145) (PHIDEP(K),K=K1,K2)
      write(14, 145) (TLOSS(K),K=K1,K2)
      write(14, 146) (LBOUND(K),K=K1,K2)
      write(14, 145) (vdep(K),K=K1,K2)
 146  format(14X,18(I1,8X))
      write(14, 145)
      write(14, 145) (TP(K),K=K1,K2)
      write(14, 145) (TL(K),K=K1,K2)
      write(14, 145) (FUP(K),K=K1,K2)
      write(14, 145) (FLOW(K),K=K1,K2)
      write(14, 145) (CON(K),K=K1,K2)

c  terse output
      write(19, 110) (ISPEC(K),K=K1,K2)    ! terse
      write(19, 120) Z(1),(USOL(K,1),K=K1,K2)
      write(19, 145) (TLOSS(K),K=K1,K2)    ! terse
      write(19, 145) (FLOW(K),K=K1,K2)
      write(19, 145) (FUP(K),K=K1,K2)
  13  CONTINUE


      write(62,*) FLOW   !write out lower boundary fluxes (not mass weighted)



C   COMPUTE CONSERVATION OF SULFUR

       Sdep=0.0
       Srain=0.0
       Spro=0.0

       do i=1,nq1
          if (flow(i).le.0.0) then
c          print *, 'Sdep:   ', ISPEC(I),FLOW(I)*atomsS(i)
           Sdep=Sdep-flow(i)*atomsS(i)
          else
c          print *, 'Sprod:  ', ISPEC(I),FLOW(I)*atomsS(i)
           Spro=Spro+flow(i)*atomsS(i)  !this knows about the ground fluxes, but not distributed fluxes
          endif

          if(i.le.nq) then
          if(lbound(i).eq.3) then
            Spro = Spro+ distflux(i)*atomsS(i)  !add in any distributed fluxes
           endif
          endif

          Srain=Srain+SR(i)*atomsS(i)
       enddo

       Sloss=Sdep+Srain
       difference  = Spro - Sloss

      SO4LOS = TLOSS(LH2SO4) + TLOSS(LSO4AER)   ! these are redundant when taking accounts
      S8LOS = 8.*TLOSS(LS8AER)   !TLOSS is rainout+deposition

       ! OK SO4LOS/S8LOS are no longer useful I think. Check this at somepoint
       ! if the are, should but a test for LS8 here just in case anyone needs to bring back
       ! S8 in the gas phase (most of us just go S4+S4 -> S8AER without bothering with S8 gas
       ! but there is a heritage code for S8 gas phase condensation. This is only really needed for high temperatures
       ! anyway - if that is back and this is useful, might need to add 8.*TLOSS(LS8) to the above


      print 177,  Sloss,Spro,SO4LOS,S8LOS,Srain, difference
      write(14, 177) Sloss,Spro,SO4LOS,S8LOS,Srain, difference
      write(19, 177) Sloss,Spro,SO4LOS,S8LOS,Srain, difference
 177  format(/1X,'CONSERVATION OF SULFUR:',/5X,'S loss =',1PE10.3,
     2  2X,'S prod =',E13.6,2X,'SO4LOS =',E13.6,2X,'8 S8LOSS =',
     3  E13.6,2X, 'sulrain =',E10.3,3x,'Difference =',E10.3)

C compute net model redox

      H2ESC = 0.
      RedDist = 0.
      OxidDist = 0.
      do i=1,NQ1
      if (redoxstate(I) > 0.) then
        OxidIn=OxidIn + FLOW(I)*redoxstate(I)
        if(FLOW(I)<0) OxidOcean=OxidOcean+FLOW(I)*redoxstate(I)
        if(FLOW(I)>0) OxidAtm=OxidAtm+FLOW(I)*redoxstate(I)
        if(LBOUND(I)>=3) then   !This is the sgflux+distflux BC
          OxidIn = OxidIn + redoxstate(I)*distflux(I)
          OxidDist = OxidDist + redoxstate(I)*distflux(I)
        endif
        OxidOut=OxidOut + FUP(I)*redoxstate(I)
        if(FUP(I)<0) OxidDown=OxidDown+FUP(I)*redoxstate(I)
        if(FUP(I)>0) OxidUp=OxidUp+FUP(I)*redoxstate(I)
        OxidRain=OxidRain + SR(I)*redoxstate(I)
        else if (redoxstate(I) < 0) then
          RedIn=RedIn+ FLOW(I)*redoxstate(I)*(-1.0)
          if(FLOW(I)<0) RedOcean=RedOcean-FLOW(I)*redoxstate(I)
          if(FLOW(I)>0) RedAtm=RedAtm-FLOW(I)*redoxstate(I)
          if(LBOUND(I)>=3) then   !This is the sgflux+distflux BC
            RedIn = RedIn + redoxstate(I)*distflux(I)*(-1.0)
            RedDist=RedDist-redoxstate(I)*distflux(I)
          endif
          RedOut = RedOut + FUP(I) * redoxstate(I) * (-1.0)
          if(FUP(I)<0) RedDown=RedDown-FUP(I)*redoxstate(I)
          if(FUP(I)>0) RedUp=RedUp-FUP(I)*redoxstate(I)
          RedRain=RedRain + SR(I)*redoxstate(I)*(-1.0)
        endif
      enddo


      OxidFlow = OxidAtm + OxidOcean
      RedFlow = RedAtm + RedOcean

      redox = RedIn - RedOut - OxidIn + OxidOut - RedRain + OxidRain

      !CSH MeanFlux takes the log average of all the fluxes as a metric for how good redox balance is.
      MeanFlux = 10.**((L10A(OxidIn)+L10A(OxidOut)+L10A(RedIn)+
     &              L10A(RedOut)+L10A(RedRain)+L10A(OxidRain))/6.)

c      print 679, OxidIn,OxidOut,RedIn,RedOut,RedRain,OxidRain,redox
c      write(14, 679) OxidIn,OxidOut,RedIn,RedOut,RedRain,OxidRain,redox
c      write(19, 679) OxidIn,OxidOut,RedIn,RedOut,RedRain,OxidRain,redox
  679 format(/,'redox budget:',/5X,'OxidIn =',1PE13.6,
     &  2X,'OxidOut =',E13.6,2X,'RedIn =',E13.6,2X,'RedOut =',
     &  E13.6,2X,'RedRain =',E13.6,2X,'OxidRain =',E13.6,
     &  3X, 'redox =',E13.6)

c      print 667 ,redox, abs(redox/MeanFlux)
  667 format(/5X,'redox conservation = ',1PE10.3,
     &  ' mol/cm^2/s, a factor of ',1PE9.3,' vs. log-mean of fluxes' /)

      red1=RedIn-RedRain-RedDist-OxidIn+OxidRain-OxidDist !Includes wet and dry dep
      !Positive when more H2 going into atmosphere

      if(H2BAL>0) then
        H2ESC = RedOut - OxidOut !Need both to cancel influx of CO and O (dissociated CO2)
        print 936,H2ESC,RedDist
        mult = 3. !this multiplication factor greatly speeds up how fast global redox balance is achieved;
        ! it still works as 1, but takes about twice as long to reach a satisfactory value.
        red2 = -1.*redoxstate(LH2)*FLOW(LH2) - mult*red1 !Need to account for current flux of H2
        !Calculate the pertinent H2 boundary condition, H2BC
        if(red2>0) H2BC= red2 !return H2 to atmosphere, as more reducing power is going into the ocean
        if(red2<0) H2BC= red2/SL(LH2,1) !draw H2 down into the ocean, as more oxidizing power is going into the ocean
        if(H2BC>0.) then
          print 933,H2BC
          write(14,933) H2BC
          write(19,933) H2BC
          print 935, H2BC
        elseif(H2BC<0.) then
          print 934,abs(H2BC)
          write(14,934) abs(H2BC)
          write(19,934) abs(H2BC)
          print 935, H2BC
        endif
        write(*,935) red1
      endif
  933 format(1X,'Global redox adjustment sets H2 flux = ',1PE16.10)
  934 format(1X,'Global redox adjustment sets H2 vdep = ',1PE16.10)
  935 format(1X,'Global redox imbalance = ',1PE17.10,/)
  936 format(1X,'Hydrogen escape= ',1PE10.4,' Volcanism= ',1PE10.4,/)

!CSH glOcean represents the net reductant flux going into the atmosphere (+) or going into
!    the ocean (-); should be equivalent to the difference between volcanic outgassing
!    (distfluxes) and escape to space.
      glOcean = RedFlow - OxidFlow - RedRain + OxidRain
      glEscape = RedOut - OxidOut - RedDist + OxidDist

      write(14,937) OxidOut,RedOut,OxidUp,RedUp,OxidDown,
     &  RedDown,OxidDist,RedDist,OxidRain,OxidAtm,RedRain,
     &  RedAtm,OxidOcean,RedOcean,OxidFlow,RedFlow,glOcean,
     &  glEscape

  937 format(/20X,'OxidOut = ',1PE11.4,/20X,'RedOut  = ',1PE11.4,
     &  /23X,'/ \',/22X,'/\',/2X,'OxidUp = ',1PE10.4,' |',/2X,
     &  'RedUp  = ',1PE10.4,' | |',/,
     &  '-----------------------|-|----------------',/23X,
     &  '| | OxidDown = ',1PE11.4,/25X,'| RedDown  = ',1PE11.4,/25X,
     &  '\/',/,'/\',/,' | OxidDist = ',1PE10.4,/,' | RedDist  = ',
     &  1PE10.4,/,' |',29X,'/\',/,' | | OxidRain = ',1PE10.4,5X,
     &  '| OxidAtm = ',1PE10.4,/,' | | RedRain  = ',1PE10.4,5X,
     &  '| RedAtm  = ',1PE10.4,/,
     &  '---|--------------------------|-|---------',/2X,
     &  '\/  OxidOcean = ',1PE11.4,'| |',/6X,'RedOcean  = ',1PE11.4,
     &  '|',/29X,'\/',/30X,'\ /',/25X,'OxidFlow = ',1PE11.4,/25X,
     &  'RedFlow  = ',1PE11.4,//4X,'global imbalance (ocean)    = ',
     &  1PE11.4,/4X,'global imbalance (dist-fup) = ',1PE11.4)

c commenting out old code - will replace with Sonny's versions
      oxid_in=0.0
      oxid_out=0.0
      red_in=0.0
      red_out=0.0
      red_rain=0.0
      oxid_rain=0.0

      do i=1,NQ1
c         print *, ISPEC(I),redoxstate(I)
         if (redoxstate(I) .GT. 0.) then
            !print *, ISPEC(I),FLOW(I)
            oxid_in=oxid_in + FLOW(I)*redoxstate(I)
            oxid_out=oxid_out + FUP(I)*redoxstate(I)
            !print *, 'fup', fup(i)
            oxid_rain=oxid_rain + SR(I)*redoxstate(I)
           ! print *,'oxy i, ispec, sr',  i, ispec(i), sr(i)
            !print *, i , ispec(i), oxid_out_new
         else if (redoxstate(I) .LT. 0) then
            !print 888, ISPEC(I),redoxstate(I)
            red_in=red_in+ FLOW(I)*redoxstate(I)*(-1.0)
            red_out=red_out + FUP(I)*redoxstate(I)*(-1.0)
            red_rain=red_rain + SR(I)*redoxstate(I)*(-1.0)
           !  print *,'red i, ispec, sr',  i, ispec(i), sr(i)
         endif
      enddo
 888  format (A8,2X,F5.1)
c 889  format (A8,2X,1PE10.3)

!distributed fluxes
      do i=1,NQ1
        if (lbound(i).eq.3) then
          if (redoxstate(I).GT.0) then
            oxid_in=oxid_in + redoxstate(I)*distflux(I)
          else if (redoxstate(I) .LT. 0) then
            red_in=red_in + redoxstate(I)*distflux(I)*(-1.0) ! +1.5*distflux(LHCL)    !ACK
          end if
        end if
      end do
!ok this needs to be finished up.  I also need to make sure that the boundary conditions are actually working as intended.
!check in particular the distributed fluxes.  In general, this should be ready to go.
!I have updated the sulfur balance, so will probably be able to follow the same scheme.
      redox = 0.
      redox = red_in-red_out-oxid_in+oxid_out-red_rain+oxid_rain
      write(14,1679) oxid_in,oxid_out,red_in,red_out,red_RAIN,oxid_rain,
     2              redox
      write(19,1679) oxid_in,oxid_out,red_in,red_out,red_RAIN,oxid_rain,
     2              redox
      print 1667 ,redox, redox/oxid_in   !mc - for ease in debugging lightning print to screen

1667  format(/,'redox conservation = ',1PE10.3,
     2  ' mol/cm^2/s, a factor of ',1PE10.3, ' NEW METHOD' /)


 919  FORMAT(1X, F8.3, 5X, E10.5, 5X, F8.3, 6X, I4, 6X, E12.6,   !STB auxiliary coupling file for RH_surf calculation according to Ramirez et al. 2014
     &     6X, E12.6,6X, E12.6,6X, E12.6,6X, E12.6, 6X, E12.6)
 918  FORMAT(1X, 'WT_surf', 8X, 'Den0', 8X, 'T_surf',
     &    9X, 'NZ', 9X, 'FLOW(LCH4)', 9X, 'FLOW(LO2)',9x,'USOL(H2O)',
     &    9X, 'O3col_depth',9x 'CH4col_depth', 9X, 'O2col_depth')
      WRITE(118,918)
      WRITE(118,919) WTa(1), Den(1), T(1), NZ, FLOW(LCH4),
     $  FLOW(LO2), USOL(LH2O,1),O3COL ,CH4_column ,O2_column

      print 932,  FLOW(LCH4), FLOW(LO2), FLOW(LCH4)/FLOW(LO2)
 932  format('CH4 flux = ',1PE13.6,'  O2 flux = ',1PE13.6,
     2  ' a ratio of ', 1PE12.4,/)

      PRINT *, 'TAUO2 = ', TAUO2
      PRINT *, 'TAUCH4 = ', TAUCH4
      PRINT *, 'TAUSO2 = ', TAUSO2
      print *,''

      write(25, 680) usol(LO2,1),usol(LCH4,1),oxid_in,oxid_out,red_in,
     1  red_out, red_RAIN,oxy_rain, redox,SULPRO,SULDEP,SULRAIN,
     2  difference,SO4LOS,S8LOS

c 680  format(2(1PE10.2,1X) ,13(E10.3,1X))
 680  format(3(1PE12.6,1X),1e10.3,1X,4(1PE12.6,1X),1e10.3,1X,
     2  6(1PE12.6,1x))

 1679 format(/,'old redox budget numbers:',/5X,'OxidIn =',1PE13.6,
     &  2X,'OxidOut =',E13.6,2X,'RedIn =',E13.6,2X,'RedOut =',
     &  E13.6,2X,'RedRain =',E13.6,2X,'OxidRain =',E13.6,
     &  3X, 'redox =',E13.6)
c 679  format(/1X,'redox budget:',/5X,'oxid_in =',1PE13.6,
c     2  2X,'oxid_out =',E13.6,2X,'red_in =',E13.6,2X,'red_out =',
c     2  E13.6,2X,'red_RAIN =',E10.3,2X,'oxy_rain =',E10.3,
c     3  3X, 'redox =',E10.3)
c      red_out = 0.5*F_esc_H + F_esc_H2 + FUP(LCO)                         ! red leaving through upper boundary (I'm assuming KZ recomputed this to check F_esc_H versus FUP(LH2)?)
c      redox = red_in-red_out-oxid_in+oxid_out-red_Rain+oxid_rain
c      write(14,679)oxid_in,oxid_out,red_in,red_out, red_RAIN,oxid_rain,
c     2             redox

c last and most terse
      write(19,973) usol(LO2,1), usol(LH2,1), usol(LCO,1),usol(LCH4,1),
     1  FLOW(LO2), FLOW(LCH4), FLOW(LH2), FLOW(LCO), distflux(LSO2),
     2  FLOW(LH2S), FUP(LH2), FUP(LH), TLOSS(LH2CO), TLOSS(LH2O2),
     3  TLOSS(LH2SO4)+TLOSS(LSO4aer), TLOSS(LSO2), TLOSS(LS8aer)
     4  ,O3COL,SL(NSP-1,1)/den(1), distflux(LH2), distflux(LCO),
     5  distflux(LO2), CO_column

  973 format(1x,1P23E10.2)

! including  Tropopause O2 mixing ratio, min Tropospheric O2 mixing ratio
! also want to include height of maximum flux of S8, then J_SO2_MIF at
!that height

c-mc
      O2atTP= usol(LO2,JTROP)          !O2 mixing ratio at TP
      tempO2min=1.0
      do 273 I=1,JTROP-1
      tempo2min=min(usol(LO2,I),tempo2min)
 273  continue
      O2minT=tempo2min                 !minimum O2 mr in troposphere

!todo kill the S8 flux part and just use YP(LS8) as I've done in out.so2
!compute the J_SO2_MIF as I did in so2.pro
!add max(J_SO2) into outputfile along with a metric something like
!if S8depflux gt 1e4 then:  if J_SO2(z) > YP(LS8,Z) then MIF = 1 ELSE MIF = 0

!ok - flux killed and removed from so2.pro - still need to implement the rest...

c-mc - three lines below should be cut - perhaps usable for JSO2
      do j=1,NZ           !at each height
        do k=20,33        !ACK hardcoded WL grid - over each wl in range from 190.5-220nm
          PSO2MIF(j) = PSO2MIF(j)+PSO2MC(k,j)    !add up photolysis rate constant
        enddo
      enddo



!now compute photochemical production in the MIF producing regime (i.e. J_SO2_MIF)
      do j=1,NZ
        PROSO2MIF(j) = PSO2MIF(j) * SL(LSO2,j)
        JSO=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'SO      ')
        PROSO(j) = prates(JSO,j)*SL(LSO,j)   !assuming the existence of SO
      enddo

      write(27,522) (PROSO2MIF(I),PROSO(I),I=1,NZ)
 522  format(2(1PE9.2,2X))


!compute maximum production rate of MIF-inducing SO2 photolysis
         PROSO2MIFmax=0.0

      do i=1,NZ
         PROSO2MIFmax= max(PROSO2MIFmax, PROSO2MIF(i))
      enddo


         TESTMIFSO2=1
         TESTMIFSO=1
         TESTMIFS=1
         TESTMIFS2=1



      if (S8los.gt.1e4) then       !if there is a decent amount of S8 deposition

         do i=1,NZ

            If (YP(LS8AER,I).gt.1e-6) then !test only in regions where S8 is produced
             if (PROSO2MIF(I).lt.YP(LS8AER,I)) TESTMIFSO2=0  !SO2 photolysis
             if (TESTMIFSO2.eq.0) goto 998

            if (PROSO(I).lt.YP(LS8AER,I)) TESTMIFSO=0   !SO photolysis
            if (TESTMIFSO.eq.0) goto 998

            if (YP(LS,I).lt.YP(LS8AER,I)) TESTMIFS=0    !S production
            if (TESTMIFS.eq.0) goto 998

            if (YP(LS2,I).lt.YP(LS8AER,I)) TESTMIFS2=0  !S2 production
            if (TESTMIFS2.eq.0) goto 998
          endif
         enddo

      else

      !who cares if there is no S8dep - there are all set to 1

      endif

 998  continue

      write(21,975) usol(LO2,1), usol(LH2,1), usol(LCO,1),usol(LCH4,1),
     1  FLOW(LO2), FLOW(LCH4), FLOW(LH2), FLOW(LCO), distflux(LSO2),
     2  FLOW(LH2S), FUP(LH2), FUP(LH), TLOSS(LH2CO), TLOSS(LH2O2),
     3  TLOSS(LH2SO4)+TLOSS(LSO4aer), TLOSS(LSO2), TLOSS(LS8aer)
     4  ,O3COL,SL(NSP-1,1)/den(1), distflux(LH2), distflux(LCO),
     5  distflux(LO2), CO_column,
     5  oxid_in,oxid_out,red_in, red_out, red_RAIN,oxy_rain, redox,
     6  abs(redox)/oxid_in,SULPRO,SULDEP,SULRAIN, difference,
     7  SO4LOS,S8LOS,S8LOS/SULPRO,finaln,O2atTP,O2minT,
     8  PROSO2MIFmax,TESTMIFSO2,TESTMIFSO,TESTMIFS,TESTMIFS2

 975  format(46(1PE10.3,1X))



      write(14, 179)
 179  format(/1X,'INTEGRATED REACTION RATES'/)
      write(14, 191)
      IROW = 10
      LR = NR/IROW + 1
      RL = FLOAT(NR)/IROW + 1
      DIF = RL - LR
      IF (DIF.LT.0.001) LR = LR - 1

      DO 17 L=1,LR
      K1 = 1 + (L-1)*IROW
      K2 = K1 + IROW - 1
      IF (L.EQ.LR) THEN
        K2 = INT(NR)
        write(14, 186) K1,(RAT(K),K=K1,K2)
  186   FORMAT(I3,2X,1P10E10.3)     !this version is not hardcoded to NR
        GO TO 17
      ENDIF
      write(14, 193) K1,(RAT(K),K=K1,K2),K2
 193  FORMAT(I3,2X,1P10E10.3,2X,I3)
   17 CONTINUE

      write(14, 191)
  191 FORMAT(9X,"1",9X,"2",9X,"3",9X,"4",9X,"5",9X,"6",9X,"7",9X,
     2    "8",9X,"9",8X,"10")

C   - the following code is taken from the Kasting group's version to
C   - print out P&L tables with integrated rxn rates, "int.rates.out.dat"
      DO 702 I=1,NSP
         WRITE(15,703) ISPEC(I),TP(I)
 703     FORMAT(/A8,12X,'PRODUCTION RXS',14X,'INT RX RATE',4X,
     2      'TP = ',1PE9.2)
       DO 704 NN=1,NR
          IF(JCHEM(3,NN).EQ.I .OR. JCHEM(4,NN).EQ.I .OR.
     2       JCHEM(5,NN).EQ.I)THEN
        IF(RAT(NN).NE.0.) WRITE(15,705) NN,(CHEMJ(J,NN),J=1,5),RAT(NN)
 705       FORMAT(1X,I3,1H),1X,A7,3H + ,A7,3H = ,A7,3H + ,A6,2X,A4,
     2      1PE10.3)
          ENDIF
 704   CONTINUE
C
       WRITE(15,706) ISPEC(I),TL(I)
 706     FORMAT(/A8,15X,'LOSS RXS',16X,'INT RX RATE',4X,'TL = ',1PE9.2)
      DO 707 NN=1,NR
          IF(JCHEM(1,NN).EQ.I .OR. JCHEM(2,NN).EQ.I)THEN
        IF(RAT(NN).NE.0.) WRITE(15,705) NN,(CHEMJ(J,NN),J=1,5),RAT(NN)
          ENDIF
 707   CONTINUE
 702  CONTINUE


c-mc
c-mc write out loop of NR

c      string='(    (E10.3))'
c      write(string(2:6),'(I3)')NZ

!the above is a way to dynamically create a format string at runtime.replaces:  599  FORMAT(55(1PE9.2,2x))
c       write(28,string) (Z(I),(prates(J,I),J=1,KJ),I=1,NZ)

      fmtstr='(    (E10.3))'
      write(fmtstr(2:5),'(I4)')NZ


      do i=1,NR
         write(28,fmtstr) (REACRAT(i,j),j=1,NZ)
      enddo

      write(28,182) JCHEM
c 182  format(5(I2,1X)) !fails if NQ>99
 182  format(5(I3,1X))



      write(28, *) ISPEC

c maybe a write here of yp and yl to rates fle


C
      write(14, 160)
 160  format(/1X,'PHOTOCHEMICAL EQUILIBRIUM AND INERT SPECIES')
      NPE = NSP - NQ             !11
      LR = NPE/IROW + 1          !11/18 + 1    (1?)
      RL = FLOAT(NPE)/IROW + 1   !11./18 + 1   (~1.6)
      DIF = RL - LR              !0.6


      IF (DIF.LT.0.001) LR = LR - 1
C
      DO 12 L=1,LR
c      K1 = NQ1 + (L-1)*IROW  !orig
      K1 = NQ1 +1 + (L-1)*IROW   !gets the final long-lived species out of the PE loop
      K2 = K1 + IROW - 1
      IF (L.EQ.LR) K2 = NSP
      write(14, 110) (ISPEC(K),K=K1,K2)
      DO I=1,NZ,ISKIP
      write(14, 120) Z(I),(SL(K,I),K=K1,K2)
      END DO
  12  CONTINUE
C
      write(14, 190)
 190  FORMAT(/1X,'ATMOSPHERIC PARAMETERS')
      write(14, 195)
 195  FORMAT(/4X,'Z',9X,'T',9X,'EDD',7X,'DEN',8X,'P',8X,'H2OSAT',
     2  5X,'H2O',7X,'RELH',5X,'CONDEN',4x,'H FLUX')
      ISKIP = 2
      write(14, 200) (Z(I),T(I),EDD(I),DEN(I),P(I),H2OSAT(I),H2O(I),
     2  RELH(I),CONDEN(I), flux_H(i), I=1,NZ,ISKIP)
 200  FORMAT(1X,1P10E10.3)
C
      IF (NP.GT.O) THEN
c-mab: THE STUFF BELOW ARE RELEVANT TO AEROSOLS/PARTICLES ONLY (AGAIN)
c-mab: unhardcoded to now depend on NP, and print out HC type as well.
c-mab: the hardcoding in the species and expected order still remain...
       DO J=1,NP
       	IF (J.EQ.1) THEN
      write(14, 230)
 230  FORMAT(/1X,'SULFATE AEROSOL PARAMETERS')
      write(14, 235)
 235  FORMAT(/4X,'Z',8X,'AERSOL',5X,'RPAR',6X,'WFALL',5X,'FSULF',4X,
     2  'TAUSED',4X,'TAUEDD',4X,'TAUC',6X,'H2SO4S',4X,'H2SO4',5X,
     3  'CONSO4',4X,'CONVER')
      LP=LH2SO4
       write(14, 240) (Z(I),AERSOL(I,J),RPAR(I,J),WFALL(I,J),FSULF(I),     !ACK - harcoded particle numbers here
     2  TAUSED(I,J),TAUEDD(I),TAUC(I,J),H2SO4S(I),USOL(LP,I),
     3  CONSO4(I),CONVER(I,J),I=1,NZ,ISKIP)
 240   FORMAT(1X,1P12E10.3)
C
       	ELSE IF (J.EQ.2) THEN
      write(14, 250)
 250  FORMAT(/1X,'S8 AEROSOL PARAMETERS')
      write(14, 255)
 255  FORMAT(/4X,'Z',8X,'AERSOL',5X,'RPAR',6X,'WFALL',5X,'TAUSED',4X,
     2 'TAUEDD',4X,'TAUC',6X,'CONVER')

       write(14, 261) (Z(I),AERSOL(I,J),RPAR(I,J),WFALL(I,J),  !ACK - harcoded particle numbers here
     2  TAUSED(I,J),TAUEDD(I),TAUC(I,J),CONVER(I,J),I=1,NZ,ISKIP)

       	ELSE IF (J.EQ.3) THEN
      write(14,301)
 301  FORMAT(/1X,"HC TYPE 1 AEROSOL PARAMETERS")
      write(14,256)
      write(14,261) (Z(I),AERSOL(I,J),RPAR(I,J),WFALL(I,J),
     2  TAUSED(I,J),TAUEDD(I),TAUC(I,J),CONVER(I,J),I=1,NZ,ISKIP)
 261  FORMAT(1X,1P8E10.3)

        	ELSE IF (J.EQ.4) THEN
      write(14,302)
 302  FORMAT(/1X,"HC TYPE 2 AEROSOL PARAMETERS")
      write(14,256)
      write(14,262) (Z(I),AERSOL(I,J),RPAR(I,J),WFALL(I,J),
     2  TAUSED(I,J),TAUEDD(I),TAUC(I,J),CONVER(I,J),I=1,NZ,ISKIP)
 262  FORMAT(1X,1P8E10.3)

        ELSE
      write(14,*) "(Other particles (if any) after HC are not printed)"

       	ENDIF
       ENDDO
 256  FORMAT(/4X,'Z',8X,'AERSOL',5X,'RPAR',6X,'WFALL',5X,'TAUSED',4X,
     2  'TAUEDD',4X,'TAUC',6X,'CONVER')

      ENDIF
c GNA
c HCAER type 1
      if (NP .GT. 2) then
      J3 = 3
c photo-clima couple file
      write(90,8901)
 8901  FORMAT(1X,"HC AEROSOL PARAMETERS")
      write(90,8257)
 8257 FORMAT(/4X,'Z',8X,'AERSOL',5X,'RPARs',6X,
     2  'WFALL',5X,'TAUSED',4X,'TAUEDD',4X,'TAUC',6X,'CONVER')
      write(90,8262) (Z(I),AERSOL(I,J3),RPAR(I,J3),WFALL(I,J3),
     2  TAUSED(I,J3),TAUEDD(I),TAUC(I,J3),CONVER(I,J3),I=1,NZ,ISKIP)
 8262 FORMAT(1X,1P8E10.3)

      if (frak .eq. 1) then
      write(68,9901)
 9901  FORMAT(/1X,"HC AEROSOL PARAMETERS")
      write(68,9256)
      print *, 'fractals #1'
 9256  FORMAT(/4X,'Z',8X,'AERSOL',5X,'RPAR',6X,'RFRAC',6X,
     2    'WFALL',5X,'TAUSED',4X,'TAUEDD',4X,'TAUC',6X,'CONVER')
      write(68,9261) (Z(I),AERSOL(I,J3),RPAR(I,J3),RFRAC(I,J3),
     2    WFALL(I,J3), TAUSED(I,J3),TAUEDD(I),TAUC(I,J3),CONVER(I,J3),
     3    I=1,NZ,ISKIP)
 9261  FORMAT(1X,1P9E10.3)

      else
      write(68,9902)
 9902 FORMAT(/1X,"HC AEROSOL PARAMETERS")
      write(68,9257)
 9257 FORMAT(/4X,'Z',8X,'AERSOL',5X,'RPARs',6X,
     2  'WFALL',5X,'TAUSED',4X,'TAUEDD',4X,'TAUC',6X,'CONVER')
      write(68,9262) (Z(I),AERSOL(I,J3),RPAR(I,J3),WFALL(I,J3),
     2  TAUSED(I,J3),TAUEDD(I),TAUC(I,J3),CONVER(I,J3),I=1,NZ,ISKIP)
 9262 FORMAT(1X,1P8E10.3)
      endif
      endif

c HCAER type 2 (mostly spectrally irrelevant - Giada)
      if (NP .GT. 2) then
      J4 = 4
      if (frak .eq. 1) then
      write(69,7901)
 7901  FORMAT(/1X,"HC AEROSOL PARAMETERS")
      write(69,7256)
      print *, 'fractals #2 '
 7256  FORMAT(/4X,'Z',8X,'AERSOL',5X,'RPAR',6X,'RFRAC',6X,
     2  'WFALL',5X,'TAUSED',4X,'TAUEDD',4X,'TAUC',6X,'CONVER')
      write(69,7261) (Z(I),AERSOL(I,J4),RPAR(I,J4),RFRAC(I,J4),
     2  WFALL(I,J4),TAUSED(I,J4),TAUEDD(I),TAUC(I,J4),CONVER(I,J4),
     3  I=1,NZ,ISKIP)
 7261  FORMAT(1X,1P9E10.3)

      else
      write(69,7902)
 7902 FORMAT(/1X,"HC AEROSOL PARAMETERS")
      write(69,7257)
 7257 FORMAT(/4X,'Z',8X,'AERSOL',5X,'RPARs',6X,
     2  'WFALL',5X,'TAUSED',4X,'TAUEDD',4X,'TAUC',6X,'CONVER')
      write(69,7262) (Z(I),AERSOL(I,J4),RPAR(I,J4),WFALL(I,J4),
     2  TAUSED(I,J4),TAUEDD(I),TAUC(I,J4),CONVER(I,J4),I=1,NZ,ISKIP)
 7262 FORMAT(1X,1P8E10.3)
      endif
      endif




C compute numerical distance between input and output files
c Err1 is the NGE, which normalizes each species so includes changes in all species down to mr of 1e-30
c Err2 is the L2 norm, which is more sensitive to major species

      ERR1 = SUM(ABS((USOL-USOLORIG)/USOLORIG),mask=USOLORIG .gt.1e-30)/
     $        COUNT(mask=USOLORIG.gt.1e-30)


      ERR2 = SUM(SQRT(ABS(USOLORIG**2 - USOL**2)))
      print *, 'Comparison between input and output files'
      print *, 'Normalized Gross Error       L2 Norm'
      print *,ERR1, ERR2
      write (50,*) ERR1, ERR2



      if (SUM(atomsCL).gt.1) then
c- screen output for chlorine stuf
      print 810, ISPEC(LHCL),ISPEC(LCLO),ISPEC(LCLO3),ISPEC(LHCLO4),
     $ ISPEC(LHNO3),ISPEC(LHO2NO2),ISPEC(LSO4AER)
      print 845, TP(LHCL),TP(LCLO),TP(LCLO3),TP(LHCLO4),TP(LHNO3),
     $ TP(LHO2NO2),TP(LSO4AER)+TP(LH2SO4)
      print 846,  FLOW(LHCL),FLOW(LCLO),FLOW(LCLO3),FLOW(LHCLO4),
     $ FLOW(LHNO3),FLOW(LHO2NO2),SO4LOS,USOL(LHCL,1)
 810  format(/8X,7(A8,2X))
 845  format(1X,'TP',4X,7(1PE9.2,1X))
 846  format(1X,'FLOW',2X,8(1PE9.2,1X))
      write(51, 810) ISPEC(LHCL),ISPEC(LCLO),ISPEC(LCLO3),ISPEC(LHCLO4),
     $ ISPEC(LHNO3),ISPEC(LHO2NO2),ISPEC(LSO4AER)
      write(51, 845) TP(LHCL),TP(LCLO),TP(LCLO3),TP(LHCLO4),TP(LHNO3),
     $ TP(LHO2NO2),TP(LSO4AER)+TP(LH2SO4)
      write(51, 846) FLOW(LHCL),FLOW(LCLO),FLOW(LCLO3),FLOW(LHCLO4),
     $ FLOW(LHNO3),FLOW(LHO2NO2),SO4LOS,USOL(LHCL,1)
      endif

c-mab:"Human readable" input and output mixing ratio profiles + P,T,Z:
c-mab: Format below abstracted to work for all templates.
      fmtdatastr='(1X,1P000E10.3)'!3 0's reserved assuming NQ < 1000 always
      fmtheadstr='(4X,"PRESS",5X,"TEMP",6X,"ALT",7X,0  (A8,2X))'
       if(NQ.lt.97) then !97 + 3 = 100
        write(fmtdatastr(8:9),'(I2)')NQ+3 !3 for P, T & Z columns
        write(fmtheadstr(36:37),'(I2)')NQ
       else
        write(fmtdatastr(7:9),'(I3)')NQ+3 !allowing NQ with 3 digits
        write(fmtheadstr(35:37),'(I3)')NQ
       endif
        !print*,'fmtdatastr=',fmtdatastr
        !print*,'fmtheadstr =',fmtheadstr
       write(34,fmtheadstr) (ISPEC(k),K=1,NQ) !headers
       write(34,fmtdatastr)  !PTZ + input LL mixing ratios
     &          (P(i),T(i),Z(i),(USOLORIG(k,i),K=1,NQ),i=1,nz)
       write(35,fmtheadstr) (ISPEC(k),K=1,NQ) !headers
       write(35,fmtdatastr) !PTZ + final output LL mixing ratios
     &          (P(i),T(i),Z(i),(USOL(k,i),K=1,NQ),i=1,nz)
c-mab: Another file with all the final output species fluxes. (to be added)
       write(255,*) "Placeholder for master....will be added later."

      RETURN

      contains
      function L10A(val)
       implicit none
       real (kind=8), intent(IN) :: val
       real (kind=8) :: L10A
       L10A = log10(abs(val))
      END FUNCTION L10A

      END
