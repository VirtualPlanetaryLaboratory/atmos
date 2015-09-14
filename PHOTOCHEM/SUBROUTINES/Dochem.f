      SUBROUTINE DOCHEM(FVAL,N,JTROP,NSHORT,USETD)
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      real*8 mass
      CHARACTER*8 ISPEC,REACTYPE,PLANET,CHEMJ
      DIMENSION FVAL(NQ,NZ),XP(NZ),XL(NZ),D(NSP2,NZ)
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PHOTABLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/BBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/CBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/GBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/NBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/RBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/ZBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/LTBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/SATBLK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/SULBLK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/RRATS.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/ISOBLOK.inc'
      COMMON/LifeTime/TAUO2,TAUCH4,TAUSO2

C
C   THIS SUBROUTINE DOES THE CHEMISTRY BY CALLING CHEMPL.  PHOTO-
C   CHEMICAL EQUILIBRIUM SPECIES ARE DONE FIRST.  THESE MUST CON-
C   TAIN NO NONLINEARITIES (SUCH AS S4 REACTING WITH ITSELF TO FORM
C   S8) AND MUST BE DONE IN THE PROPER ORDER (I.E. IF SPECIES A
C   REACTS TO FORM B, THEN A MUST BE FOUND FIRST).  LONG-LIVED
C   SPECIES CAN BE DONE IN ANY ORDER.



!compute number densities for long-lived and particle species 
!these are re-computed below...
      do I=1,NQ1
       do J=1,NZ
        if (I.LE.NQ) then 
c         if(j.eq.1)print *, i, ISPEC(I),' long-lived densities'
          D(I,J) = USOL(I,J) * DEN(J)
c          print  *, usol(1,j)

        else 
c         if(j.eq.1)print *, i, ISPEC(I),' particle densities'
          D(I,J) = PARTICLES(J,I-NQ)* DEN(J)   !gets particles if using tri-diag
        endif   
       enddo
      enddo  

! compute densities for INERT species in ISOTOPE code
      if (ISOTOPE.EQ.1) then
         do I=NQ1+NSHORT+1, NSP-2   !loop over INERT species
c            print *, i, ISPEC(I),' inert densities'
          do J=1,NZ
             D(I,J) = UINERT(I-Loff,J)*DEN(J)
          enddo
         enddo   
      endif

!do the last 4 inert species that are the same in both codes
      do J=1,NZ

       if(ISPEC(NSP-1).eq.'CO2') D(LCO2,J) = FCO2 * DEN(J) !if CO2 is inert, place above N2 in list (hardcoded position for CO2/N2 (NSP-1/NSP)
       if(ISOTOPE.EQ.0) D(NSP,J) = (1. - USOL(LO2,J) - FCO2 - FAR 
     2  - FCO)* DEN(J)    ! N2 -defined as the "rest" of the density - EWS - now includes Argon and CO. 
       if(ISOTOPE.EQ.1) D(NSP,J)=(1.-UINERT(LO2-Loff,J)-FCO2 - FAR 
     2  - FCO) * DEN(J)  ! N2 -defined as the "rest" of the density - EWS - now includes argon, CO
       D(NSP1,J) = 1.                             ! HV has density of 1 for photorate calculations
       D(NSP2,J) = DEN(J)                         ! M - background density for three body reactions
      enddo

      
      if (N.GE.0) then  !normal operation mode (-1 just fills up D and SL for first timestep)
C
C ***** SOLVE FOR THE PHOTOCHEMICAL EQUILIBRIUM SPECIES *****
C

      SLS4=0  !checking if S4 is the the short-lived loop

      DO 3 I=NQ1+1,NQ1+NSHORT  !loop through the equilibrium species...
c         print *, I,ISPEC(I),' short-lived'
         if (ISPEC(I).EQ.'S4') ISS4SL=I

      CALL CHEMPL(D,XP,XL,I)
      DO 3 J=1,NZ
   3  D(I,J) = XP(J)/XL(J)


C   SOLVE QUADRATIC FOR S4, if S4 is in the SL lived loop
c equation is production=loss, which turns into a quadratic equation in S4 density
c production (C) is from S+S3 and S2+S2
c losses are S4+S4 ->S8AER (A) and S4 + Hv -> S2+S2 (B)
c so we have C=A*DenS4^2 + B*DENS4*DENHv, where DENHv=1 by definition.

      if(ISS4SL.gt.0) then

c add a loop over all reactions here, just dipping into to do loop if S4 is invloved
can use the chempl stuff, but need to figure out K, the species number of S4
         NPS4=NUMP(ISS4SL)
         NLS4=NUML(ISS4SL)
       DO J=1,NZ

        AQ = 2.*A(148,J)    !S4+S4 -> S8Aer loss term
        BQ = A(149,J)       ! S4+ Hv -> S2 + S2 loss term 
        CQ = A(146,J)*D(LS2,J)*D(LS2,J) + A(147,J)*D(LS,J)*D(LS3,J)   ! production terms

c        CQ=0.0

!ack need to abstract the above - below is a start, altough something is wrong.  check carefully against code in CHEMPL and shawn's code
!the problem is the N - it is overwriting one of the inputs... anyway.

!ACK - also down at the bottom of the file is an s8col printout that is suppressed.  Also a diagnostic printout in out.so2 in output

        !do the production terms generically
c        do I=1,nps4
c          L=IPROD(ISS4SL,I)   !reaction number
c          M = JCHEM(1,L)   !reactant 1 for reaction number J
c          N = JCHEM(2,L)   !reactant 2 for reaction number J
c          CQ=CQ + A(L,J)*D(M,J)*D(N,J) !rate*density1*density2 
c        enddo   

       DLS4 = (SQRT(BQ*BQ + 4.*AQ*CQ) - BQ)/(2.*AQ)
       D(LS4,J) = MAX(DLS4,1.e-99)
       enddo

      endif

C
C ***** LONG-LIVED SPECIES CHEMISTRY *****
      DO I=1,NQ
c         print *, I,ISPEC(I),' long-lived'
      CALL CHEMPL(D,XP,XL,I)

      if (ISOTOPE.EQ.0) then
       if (ISPEC(I).EQ.'O2')   TAUO2 = 1/XL(1)
       if (ISPEC(I).EQ.'CH4') TAUCH4 = 1/XL(1)
       if (ISPEC(I).EQ.'SO2') TAUSO2 = 1/XL(1)
      endif


      DO  J=1,NZ
      XLJ = XL(J) + RAINGC(I,J)
      FVAL(I,J) = XP(J)/DEN(J) - XLJ*USOL(I,J)

c      print *, usol(i,j)
      YP(I,J) = XP(J)
      YL(I,J) = XLJ
c      IF (ISPEC(I).EQ.'H2SO4') print *, J,XL(J),RAINGC(I,J),XP(j),XLJ
      ENDDO
      ENDDO


      if(USETD.EQ.1) then
C ***** TRIDIAGONAL SPECIES  *****
      do I = NQ+1, NQ1
c         print *, I,ISPEC(I),' tri-diag'
      CALL CHEMPL(D,XP,XL,I)
       do J = 1, NZ
c        YL(I,J) = XL(J) + RAINGC(LH2SO4,J)   !this seems wrong for S8 but is status quo
c        YL(I,J) = XL(J) + RAINGC(LSO4AER,J)   !this seems wrong for S8 but is status quo


c - these two below are the behavior that I am abstracting, but may be wrong
c         if (ISPEC(I).EQ.'SO4AER') YL(I,J) = XL(J) + RAINGC(LSO4AER,J)   
c         if (ISPEC(I).EQ.'S8AER') YL(I,J) = XL(J) + RAINGC(LS8AER,J)   
c - Jim/Kevin's code originally had the behavoir where H2SO4 was used for all particles.  Is this the source of the sulfur redox balance issues?

C-mc - OK, if I come back to this later I should remember two things here:
c- there is a mix of two different codes here. In one, I took RAINGC , HEFF, and some of the others to NQ+NP, and computed them directly.  If we are going to do the tridiagonal, they should be removed completely and I should go back to the original behavior, which is jsut uing H2SO4 everywhere.  At the same time, Kevin had S8 not being rained out at all because it isn't soluble. Either way, what we have in this section is both wrong and confusing...


c         YL(I,J) = XL(J) + RAINGC(I,J)
c         YL(I,J) = XL(J) + RAINGC(LH2SO4,J)   !back to original behavior (will need to tweak for isotope) !ISOHACK

         YL(I,J) = XL(J)+ RAINGC(LH2SO4,J)   !back to original behavior (will need to tweak for isotope) !ISOHACK
         YP(I,J) = XP(J)
c         if (ISPEC(I).EQ.'SO4AER') YL(I,J) = YL(I,J) + RAINGC(LH2SO4,J)   
       enddo   
      enddo   

      endif


      if (PLANET .EQ. 'EARTH') then 
       CONFAC = 1.6E-5          !condensation factor
      else if (PLANET .EQ. 'MARS') then
       CONFAC = 1.6E-5 * 10.    ! reduce supersaturation of stratosphere
      else if (PLANET .EQ. 'DRY') then !added by EWS 9/14/2015
       CONFAC = 1.0 ! atmosphere is dry everywhere
      endif   



C   ZERO OUT H2O TERMS IN THE TROPOSPHERE AND INCLUDE LIGHTNING
C   PRODUCTION OF NO AND O2

      if (ISOTOPE.EQ.0) then   !skip the next little bit for the ISOTOPE code as all this stuff is inert...
      changeL = 1        !mc temp var for testing lightning changes versus OLD JFK method
                         !changeL=1 uses new code, changeL=0 uses old code

      JT1 = JTROP + 1          ! same as NH1
      DO 5 J=1,JTROP
        FVAL(LH2O,J) = 0.
        YP(LH2O,J) = 0.
        YL(LH2O,J) = 0.
        SCALE = RAIN(J)/RAIN(1)

        ZAP = ZAPNO * SCALE
        FVAL(LNO,J) = FVAL(LNO,J) + ZAP/DEN(J)  
        YP(LNO,J) = YP(LNO,J) + ZAP


        if (changeL.eq.1) then
c-mc 4/28/06      making NO requires subtracting 1/2 O2   (1/2 N2 + 1/2 O2 <-> NO)
        FVAL(LO2,J) = FVAL(LO2,J) - 0.5*ZAP/DEN(J)
        YP(LO2,J) = YP(LO2,J) - 0.5*ZAP  
        else

c-mc   4/28/06 the following is no longer needed as the code computes how much of each reductant is produced
c-mc   i use these numbers to compute the amount of O2 produced. This is fine because CO2 and H2O reservoirs are effectivly infinte
c-mc  "un-commenting" these for test versus JFK's original code.  in the else statement
        ZAP = ZAPO2 * SCALE
        FVAL(LO2,J) = FVAL(LO2,J) + ZAP/DEN(J)  
        YP(LO2,J) = YP(LO2,J) + ZAP
        endif
        



        if (changeL.eq.1) then

c - mc adding lightning CO/H2 into chemistry for better redox conservation
c- kevin's addition of 3-20-06
        ZAP = ZAPH2 * SCALE
        FVAL(LH2,J) = FVAL(LH2,J) + ZAP/DEN(J)  
        YP(LH2,J) = YP(LH2,J) + ZAP

c-mc 4/28/06 adding H2 also requires adding 1/2 O2    (H2O <-> H2 + 1/2 O2)
        FVAL(LO2,J) = FVAL(LO2,J) + 0.5*ZAP/DEN(J)
        YP(LO2,J) = YP(LO2,J) + 0.5*ZAP
c-mc

        ZAP = ZAPCO * SCALE
        FVAL(LCO,J) = FVAL(LCO,J) + ZAP/DEN(J)  
        YP(LCO,J) = YP(LCO,J) + ZAP

c-mc 4/28/06  adding CO also requires adding 1/2 O2
        FVAL(LO2,J) = FVAL(LO2,J) + 0.5*ZAP/DEN(J)
        YP(LO2,J) = YP(LO2,J) + 0.5*ZAP
c-mc


c-mc now figure out how much O is produced: add to O, subtract 1/2 O2  (1/2 O2 <-> O)
         ZAP = ZAPO * SCALE
         FVAL(LO,J) = FVAL(LO,J) + ZAP/DEN(J)
         YP(LO,J) = YP(LO,J) + ZAP
         FVAL(LO2,J) = FVAL(LO2,J) - 0.5*ZAP/DEN(J)
         YP(LO2,J) = YP(LO2,J) - 0.5*ZAP
c-mc

         endif

c        stop
c-end 3-20-06 addition

   5  CONTINUE


! ACK - this may be part of the reason the time-dependent code is having trouble
C
C   H2O CONDENSATION IN THE STRATOSPHERE
C   (RHCOLD IS THE ASSUMED RELATIVE HUMIDITY AT THE COLD TRAP)
c dunno what to do here, I'll take it to be small
      if (PLANET .EQ. 'EARTH') then 
      RHCOLD = 0.1   ! Jim had 0.1 ; what needs to be here is something that will give the right stratospheric H2O 3ppm
      else if (PLANET .EQ. 'MARS') then
c      RHCOLD = 0.4  ! Jim had 0.1 ?  my standard is 0.4  <-Kevin words (mc - this seems wrong)
      RHCOLD = 0.17  ! from Kevin's Mars paper
c       RHCOLD = 0.10  ! Jim had 0.1 ?  my standard is 0.4  <-Kevin words (mc - this seems wrong)
      endif   
      DO 13 J=JT1,NZ
        H2OCRT = RHCOLD * H2OSAT(J)
        IF (USOL(LH2O,J) .LT. H2OCRT) GO TO 13
        CONDEN(J) = CONFAC * (USOL(LH2O,J) - H2OCRT)   !this is saved in SATBLK to be printed out in output file - wont be missed in ISOCODE
        FVAL(LH2O,J) = FVAL(LH2O,J) - CONDEN(J)
  13  CONTINUE

      endif  !end skip loop for isotope code

C
C   H2SO4 CONDENSATION

      if (ISOTOPE.EQ.0) then
       LL=LH2SO4
       LLA=LSO4AER
      else
       LL=LH2SXO4
       LLA=LSXO4AER
      endif


      DO 14 J=1,NZ
      CONSO4(J) = CONFAC * (USOL(LL,J) - H2SO4S(J))

      if (CONSO4(j).gt.0) then  !dont allow artifical evaporation
      FVAL(LL,J) = FVAL(LL,J) - CONSO4(J)
      if(USETD.EQ.0) FVAL(LLA,J) = FVAL(LLA,J) + CONSO4(J) !else handled in RHS of tri-diag in main code
      YL(LL,J) = YL(LL,J) + CONFAC
      YP(LL,J) = YP(LL,J) + CONFAC*H2SO4S(J)*DEN(J)  
      YP(LLA,J) = YP(LLA,J) + CONSO4(J)*DEN(J)
c      print *,j,CONSO4(J),CONFAC*H2SO4S(J)*DEN(J),CONSO4(J)*DEN(J)       
c      print *,j,USOL(LL,J),H2SO4S(J),USOL(LL,J) - H2SO4S(J)
      endif   
c      endif

  14  CONTINUE

c      stop

C    !what follows is not in Jim's code

C   S8 CONDENSATION (this is needed if we every want to deal with 'hot air' - s8 stays in the vapor phase       

      skipS8=1

      if (skipS8.eq.0) then 
  
      if(ISOTOPE.EQ.0) then
       LL=LS8
       LLA=LS8AER
      else
       LL=LSXS7
       LLA=LSXS7AER
      endif


      CONFC2 = 1.E-2   !whats this?  CONFAC = 1.6E-5   whatever that is

      DO  J=1,NZ  !(fixes an error in our earlier codes where S8 didn't condense in the troposphere)
c      DO 15 J=JT1,NZ  !(fixes an error in our earlier codes where S8 didn't condense in the troposphere) (temp return..)
c      EVAPS8(J) = 0.!what is this? it appears to do nothing,and is printed out as 0 in the output file...

      CONS8(J) = CONFC2 * (USOL(LL,J) - S8S(J))
!note that S8S is higher in the ISOTOPE CASE in order to make these budgets work out correctly

      IF (CONS8(J).gt.0.) then !dont allow artifical evaporation

      FVAL(LL,J) = FVAL(LL,J) - CONS8(J)
      if(USETD.EQ.0) FVAL(LLA,J) = FVAL(LLA,J) + CONS8(J)  !else do RHS in tri-diag
       YL(LL,J) = YL(LL,J) + CONFC2
       YP(LL,J) = YP(LL,J) + CONFC2*S8S(J)*DEN(J)
       YP(LLA,J) = YP(LLA,J) + CONS8(J)*DEN(J)

c       print *,j,CONS8(J),CONFC2*S8S(J)*DEN(J),CONS8(J)*DEN(J)       
c       print *, j, USOL(LL,J),S8S(J),(USOL(LL,J) - S8S(J))

      endif
      enddo

      endif  !end S8 skip loop

      endif  !end normal operation loop (i.e. if IDO = 0 or 1)


c      stop

C
c what is this?  its hardwired for case where O3 is trace
c     DO 7 J=1,NZ
c  7  O3(J) = D(LO3,J)/DEN(J)


C ***** SAVE THESE DENSITIES FOR PRINTOUT *****
c-mc and for allowing short-lived species to photolyze...?
c-mc modifying to contain all species.  Will this affect the isotope model?
! orig      DO 9 I=NQ1,NSP   
      DO 9 I=1,NSP   
      DO 9 J=1,NZ
   9  SL(I,J) = D(I,J)



      IF(N.LT.1) RETURN              !on final time step compute TP and TL

C
C ***** CALCULATE COLUMN-INTEGRATED PRODUCTION AND LOSS *****
      O3COL = 0.
      H2SCOL = 0.
      SO2COL = 0.
      S2COL = 0.
      S4COL = 0.
      S8COL = 0.
C
      DO 10 L=1,NR
       DO 420 J=1,NZ
 420    REACRAT(L,J)=0
  10  RAT(L) = 0.



      DO 11 K=1,NQ1    !num species
      TP(K) = 0.
      TL(K) = 0.
  11  CONTINUE
C

      DO J=1,NZ
      RELH(J) = D(LH2O,J)/DEN(J)/H2OSAT(J)  !this gets H2O mixing ratio in a more general way
      O3COL = O3COL + D(LO3,J)*DZ(J)

      if(ISOTOPE.EQ.0) then
       H2SCOL = H2SCOL + D(LH2S,J)*DZ(J)
       SO2COL = SO2COL + D(LSO2,J)*DZ(J)
       S2COL = S2COL + D(LS2,J)*DZ(J)
       S4COL = S4COL + D(LS4,J)*DZ(J)
c       S8COL = S8COL + D(LS8,J)*DZ(J)   !ACK need to deal if S8 in gas phase
      else
       H2SCOL = H2SCOL + D(LH2SX,J)*DZ(J)
       SO2COL = SO2COL + D(LSXO2,J)*DZ(J)
       S2COL = S2COL + D(LSXS,J)*DZ(J)
       S4COL = S4COL + D(LSXS3,J)*DZ(J)
c       S8COL = S8COL + D(LSXS7,J)*DZ(J)  !ACK need to deal if S8 in gas phase
      endif   


      enddo
C
      DO 12 L=1,NR
      M = JCHEM(1,L)         !identifies first reactant of equation L
      K = JCHEM(2,L)         !second reactant
      DO 12 J=1,NZ
         REACRAT(L,J) =  A(L,J)*D(M,J)*D(K,J)  !reaction rate*densities
c cm^3/mol/s * (mol/cm^3)^2 ->  mol/cm^3/s (i.e. rate units)

  12  RAT(L) = RAT(L) + REACRAT(L,J)*DZ(J)  
c      mol/cm^3/s * cm ->  mol/cm^2/s (i.e. RAT is in height integrated flux units)

      DO 8 I=1,NQ1
      XLG(I) = YL(I,1)
      DO 8 J=1,NZ
      TP(I) = TP(I) + YP(I,J)*DZ(J)
      TL(I) = TL(I) + YL(I,J)*D(I,J)*DZ(J)
   8  CONTINUE
C
      RETURN
      END
