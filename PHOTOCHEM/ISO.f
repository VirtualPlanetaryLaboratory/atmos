            PROGRAM ISO

c-mc first attempt at making an ISOTOPE code - note I am effectivly branching off TOTCtester at this point.
! a goal should be to re-implement this to TOTCs with ifs, but that is for down the road...
      INCLUDE 'INPUTFILES/parameters.inc'   !this file needs to have the same name because you can use IF for include statements.

      implicit real*8(A-H,O-Z)
      CHARACTER*8 ISPEC, CHEMJ,SPECIES,SPECTYPE,REACTYPE,PLANET
      CHARACTER*20 string,fmtstr,fmtstr2
      real*8 mass      
      DIMENSION FVAL(NQ,NZ),FV(NQ,NZ),DJAC(LDA,NEQ),RHS(NEQ),IPVT(NEQ)
      dimension SMFLUX(NQ),SGFLUX(NQ),VDEP(NQ),VEFF(NQ)
      dimension   DD(NQ1,NZ),DL(NQ1,NZ),DU(NQ1,NZ)
      dimension ADL(NQ,NZ), ADU(NQ,NZ), ADD(NQ,NZ)
      dimension USAVE(NQ,NZ),R(NZ),U(NQ),VDEP0(NQ),VEFF0(NQ)
      dimension USAVEOLD(NQ,NZ), USOLSAVE(NQ,NZ), USOLPREV(NQ,NZ)
      dimension USOLORIG(NQ,NZ)
!c-mc  USOLSAVE added for convergence testing 4/29/06, removed from code on 2/28/07
!c-mc  going to keep declaration here and reserve it in case I ever bring back the
!c-mc  the convergence-testing code, which is r37:36 in the /td branch
c-mc  USOLPREV(NQ,NZ) added for second order reverse Euler calculations
      DIMENSION DPU(NZ,NP),DPL(NZ,NP)
      DIMENSION TA(NZ),TB(NZ),TC(NZ),TY(NZ)

c      dimension atomsO(NSP2),atomsH(NSP2),atomsC(NSP2)   
c      dimension atomsN(NSP2),atomsCL(NSP2),atomsS(NSP2)

!temp
      dimension vdep1(nq),USOL2(NQ,NZ),testvec(NR),fixedmr(nq)
      dimension distheight(nq)

      CHARACTER*11 photolabel, AA
      CHARACTER*5 ISOPREFIX

      INCLUDE 'DATA/INCLUDE/ABLOK.inc'
      INCLUDE 'DATA/INCLUDE/BBLOK.inc'
      INCLUDE 'DATA/INCLUDE/CBLOK.inc'
      INCLUDE 'DATA/INCLUDE/DBLOK.inc'
      INCLUDE 'DATA/INCLUDE/EBLOK.inc'   !can go away when MSCAT does
      INCLUDE 'DATA/INCLUDE/FBLOK.inc'
      INCLUDE 'DATA/INCLUDE/GBLOK.inc'
      INCLUDE 'DATA/INCLUDE/JBLOK.inc'
      INCLUDE 'DATA/INCLUDE/NBLOK.inc'
      INCLUDE 'DATA/INCLUDE/RBLOK.inc'
      INCLUDE 'DATA/INCLUDE/SBLOK.inc'
      INCLUDE 'DATA/INCLUDE/ZBLOK.inc'
      INCLUDE 'DATA/INCLUDE/LTBLOK.inc'
      INCLUDE 'DATA/INCLUDE/AERBLK.inc'
      INCLUDE 'DATA/INCLUDE/SULBLK.inc'
      INCLUDE 'DATA/INCLUDE/PBLOK.inc'
      INCLUDE 'DATA/INCLUDE/ISOBLOK.inc'




C   CONSTANTS FOR 1ST EXPONENTIAL INTEGRAL  !can go away when MSCAT does
      DATA AI/.99999193, -.24991055, .05519968, -.00976004,
     2  .00107857, -.57721566/
      DATA BI/8.5733287401, 18.0590169730, 8.6347608925,
     2  .2677737343/
      DATA CI/9.5733223454, 25.6329561486, 21.0996530827,
     2  3.9584969228/
      DATA NUML,NUMP/NSP*0,NSP*0/


      ISOTOPE=1    !Isotope flag

      if (ISOTOPE.EQ.1) then
       print *, 'enter minor Sulfur Isotope (33 or 34)'
       read(*,'(I2)') ISOS
       ISOPREFIX='ISO  '
       write(ISOPREFIX(4:5),'(I2)')ISOS
      else 
       ISOPREFIX = ''   !of course needs testing to see if this works in main code...
      endif

c OPEN FILES  (ISOHACK - would just need some iff's in place for the opens...)

      open(2, file='DATA/aerosol.table',status='OLD') !,form='UNFORMATTED')  
      open(3, file='DATA/photo.dat',status='OLD')  
      open(4, file='INPUTFILES/ISOspecies.dat',status='OLD')    
      open(7, file='INPUTFILES/PLANET.dat',status='OLD')                !planet parameters (G, FSCALE, ALB,ZTROP,etc)
      open(231, file='INPUTFILES/input_photchem.dat',status='OLD')       !model parameters (AGL, IO2,INO, LGRID, etc)
      open(9, file='INPUTFILES/ISOreactions.rx',status='OLD')          ! reaction file

c      open(14, file=ISOPREFIX//'out.out',status='UNKNOWN')    ! output  !ISOHACK
      open(14, file='out'//ISOPREFIX//'.out',status='UNKNOWN')    ! output  !ISOHACK

      open(17, file='ISOin.dist',status='OLD')         ! formatted input (isotope species) !ISOHACK
      open(55, file='ISOinert.dist',status='OLD')         ! formatted input (all other "inert" species) !ISOHACK

c      open(18, file=ISOPREFIX//'out.dist',status='UNKNOWN')    ! formatted output !ISOHACK
      open(18, file='out'//ISOPREFIX//'.dist',status='UNKNOWN')    ! formatted output !ISOHACK


c      open(19, file='out.terse',status='UNKNOWN')    ! terse output
c      open(21, file='out.trs',status='UNKNOWN')    ! very terse output
c      open(23, file='out.time',status='UNKNOWN')    ! time output
c      open(24, file='out.tim',status='UNKNOWN')    ! time output
c      open(25, file='out.redox',status='UNKNOWN')    ! redox output - eventually combine into out.trs
c      open(26, file='out.converge',status='UNKNOWN')  ! temporary output file for looking at convergence
                                                     ! in the reverse Euler iteration
c      open(27, file=ISOPREFIX//'out.so2',status='UNKNOWN')      !printing out relevant SO2 photolysis pieces
                                               !between 190-220 as MIF signature
      open(27, file='out'//ISOPREFIX//'.so2',status='UNKNOWN')      !printing out relevant SO2 photolysis pieces


!revisit above - I probably want this (or something like it ) in isotope and fine grid codes...

c       open(28, file=ISOPREFIX//'out.rates',status='UNKNOWN') !reaction rates,reactions,species,densities,rate constants
       open(28, file='out'//ISOPREFIX//'.rates',status='UNKNOWN') !reaction rates,reactions,species,densities,rate constants



c       open(58, file='out.prates',status='OLD') !photolysis rates from main run !ISOHACK
c       open(59, file='out.raingc',status='OLD') !rainout rates from main run !ISOHACK

c      open(29, file='out.xsec',status='UNKNOWN')  !cross sections

c      open(30, file=ISOPREFIX//'out.gridz',status='UNKNOWN')  !height grid
      open(30, file='out'//ISOPREFIX//'.gridz',status='UNKNOWN')  !height grid


c      open(31, file='out.gridw',status='UNKNOWN')  !wavelength grid

c       open(41, file=ISOPREFIX//'out.rad',status='UNKNOWN')
       open(41, file='out'//ISOPREFIX//'.rad',status='UNKNOWN')



c      open(42, file='out.finalden',status='UNKNOWN') !this file should contain TATLTUAE
c      open(43, file='out.densities',status='UNKNOWN') !number densities at each timestep

c      open(44, file=ISOPREFIX//'out.prod',status='UNKNOWN') !total production and loss at steady state
      open(44, file='out'//ISOPREFIX//'.prod',status='UNKNOWN') !total production and loss at steady state



c      open(45, file='out.flux',status='UNKNOWN') !fluxes
c      open(48, file='out.tau',status='UNKNOWN') !tau=1 (at final step)


c       open(49, file=ISOPREFIX//'out.params',status='UNKNOWN') !some model params
       open(49, file='out'//ISOPREFIX//'.params',status='UNKNOWN') !some model params


c      open(50, file='out.error',status='UNKNOWN') !NGE and L2 norm between start and finish
c      open(51, file='out.cl',status='UNKNOWN') !TP/FLOW for chlorine species,nitrate, adn sulfate

       open(62, file='out'//ISOPREFIX//'.flow',status='UNKNOWN') !lower boundary fluxes (not including rainout)

       open(61, file='out'//ISOPREFIX//'.so2HR',status='UNKNOWN') !wavelength specific so2 photorates on HR grid

C - other model parameters read in from input_photochem.dat
      READ(231,555)
      if(IDEBUG.eq.1) print *, "input_photochem.dat data:"
      READ(231,*)AA,AGL
      IF(IDEBUG.eq.1) print *, "AGL =",AGL
      READ(231,*)AA,ISEASON
      IF(IDEBUG.eq.1) print *, "ISEASON =",ISEASON
      READ(231,*)AA,IZYO2
      IF(IDEBUG.eq.1) print *, "IZYO2 =",IZYO2
      READ(231,*)AA,LGRID
      IF(IDEBUG.eq.1) print *, "LGRID =",LGRID
      READ(231,*)AA,IO2
      IF(IDEBUG.eq.1) print *, "IO2 =",IO2
      READ(231,*)AA,INO
      IF(IDEBUG.eq.1) print *, "INO =",INO
      READ(231,*)AA,EPSJ
      IF(IDEBUG.eq.1) print *, "EPSJ =",EPSJ
      READ(231,*)AA,PRONO
      IF(IDEBUG.eq.1) print *, "PRONO =",PRONO
      READ(231,*)AA,frak
      IF(IDEBUG.eq.1) print *, "frak =",frak
 555  format(3/)
      close(231)

      


c-mc  reserving unit 33 for wavelength grids...


C - READ IN SPECIES NAMES, ATOMIC NUMBERS, AND BOUNDARY CONDITIONS

      iLL=0   !counter for long lived species
      iSL=0   !counter for short lived species
      iTD=0   ! counter for tridiagonal species
      iIN=0   !counter for inert species
      iSP=0   !counter for number of lines in species.dat file

      do while (I.LT.300)  !ACK this will crash if species.dat is longer than 300 lines.
         read (4,203, end=96) SPECIES,SPECTYPE
         if (scan(species,'*').LT.1) then    !ignore comments in species.dat file
            iSP=iSP+1
            ISPEC(iSP)=species
            print *, iSP, species
            call LNUM(ISPEC(isP),iSP)  !this loads the "Lnumbers" for ease of use later in the code

              backspace 4  !return to previous line in species.dat file
              read(4,207) LA,LB,LC,LD,LE,LF   !read in atmoic number data - NOTE: don't ever use LH,LN,LO,LS as placeholders as these mean something.

            if (SPECTYPE.EQ.'LL') then 
               iLL=iLL+1
               backspace 4  !return to previous line in species.dat file
               read(4,208) LBC, XX,YY,ZZ,XXX,LG,YYY,ZZZ  !read in boundary conditions
               LBOUND(iLL)=LBC
               VDEP0(iLL)=XX
               FIXEDMR(iLL)=YY
               if (LBOUND(iLL).eq.3) then
                distflux(iLL)=ZZ !distributed flux
               else 
                  SGFLUX(iLL)=ZZ  !lower boundary flux
               endif   
               distheight(iLL)=XXX
               MBOUND(iLL)=LG
               SMFLUX(iLL)=YYY
               VEFF0(iLL)=ZZZ
            endif   

            if (SPECTYPE.EQ.'IN') then
               iIN=iIN+1
c               backspace 4  !return to previous line in species.dat file   !ISOHACK
c               read(4,209) XX  !read in fixed mixing ratios                !ISOHACK
c               if (species.EQ.'CO2') FCO2=XX   !hardcoding woohoo!  !need to do N2 as well  !ISOHACK
            endif   

            if (SPECTYPE.EQ.'TD')iTD=iTD+1
            if (SPECTYPE.EQ.'SL') iSL=iSL+1
            if (SPECTYPE.EQ.'HV') iIN=iIN+1
            if (SPECTYPE.EQ.'M') iIN=iIN+1

               atomsO(iLL+iTD+iSL+iIN)=LA
               atomsH(iLL+iTD+iSL+iIN)=LB
               atomsC(iLL+iTD+iSL+iIN)=LC
               atomsS(iLL+iTD+iSL+iIN)=LD
               atomsN(iLL+iTD+iSL+iIN)=LE
               atomsCL(iLL+iTD+iSL+iIN)=LF

         endif
         I=I+1
      enddo   


 203  FORMAT(A8,3X,A2)  !for species name and type
 207  format(15X,6(I1,1X))      !for elemental counts
! 208  format(30X,I1,5X,4(E7.1,1X),I1,6X,2(E7.1,1X))  !for boundary conditions (original)
 208  format(30X,I1,5X,2(E7.1,1X),E9.3,1X,E7.1,1X,I1,6X,2(E7.1,1X))  !for boundary conditions
c 208  format(30X,I1,5X,2(E8.1),E9.3,1X,E7.1,1X,I1,6X,2(E7.1,1X))  !for boundary conditions
 209  format(30X,E7.1) !for INERT species boundary conditions

 96   CONTINUE

      mass=atomsO*16.+atomsH*1.+atomsC*12.+atomsS*32.+atomsN*14.
     $  + atomsCL*34. !molecular weights (NQ1)
!so far, the only use of mass in Difco,where it is summed over NQ, so would be OK to go higher


      redoxstate = atomsO*1.0 + atomsH*(-0.5) + atomsS*(-2.) + 
     $  atomsCL*(-1.0) + atomsC*(-2)  !N=0
!we are setting CLO as 'neutral"
!redoxstate goes from 1-NQ1 in Output

c      print *, redoxstate
c      print *, fixedmr
c      print *, sgflux
c      print *, distflux
c      print *, distheight
c      stop


      if (iTD.gt.0) then
         USETD=1
         print *, 'using tri-diagonal solver for particles'
         print *, 'redox and sulfur consevation diagnostics not '
         print *, 'guaranteed to work correctly'
      else
         USETD=0
      endif



C ***** READ THE CHEMISTRY DATA CARDS ***** 
corig      read (9,200) CHEMJ                 
      read (9,200) CHEMJ                 
corig 200  FORMAT(10X,A8,2X,A8,2X,A8,2X,A8,2X,A8)
 200  FORMAT(A8,2X,A8,2X,A8,2X,A8,2X,A8)
      write(14, 201) (J,(CHEMJ(M,J),M=1,5),J=1,NR)
 201  FORMAT(1X,I3,' ',5X,A8,' +  ',A8,'  =    ',A8,' +  ',A8,4X,A8)
      KJAC = LDA*NEQ
      write(14, 202) NQ,NZ,KJAC
 202  FORMAT(//1X,'NQ=',I2,5X,'NZ=',I3,5X,'KJAC=',I7)



c-mc this is bad code.  there is probably some way to do this with the read in above
c-mc or even if necessary, some way to redo the read without closing and opening the file again
c-mc whatever
      close(9)
      open(9, file='INPUTFILES/ISOreactions.rx',status='OLD')          ! chemical reaction file  !ISOHACK
      read(9,204) REACTYPE
 204  FORMAT(48X,A5)

      close(9)  !close this because Rates.f and Initphoto.f will re-open it later.



C  ****  JCHEM has species numbers; CHEMJ is corresponding characters
C ***** REPLACE HOLLERITH LABELS WITH SPECIES NUMBERS IN JCHEM *****
      DO 5 J=1,NR
      DO 5 M=1,5
      IF(CHEMJ(M,J).EQ.' ') GO TO 5
      DO 6 I=1,NSP2
      IF(CHEMJ(M,J).NE.ISPEC(I)) GO TO 6
      JCHEM(M,J) = I
      GO TO 5
   6  CONTINUE
      IERR = J
      print *, ISPEC
      print *, (CHEMJ(L,J),L=1,5)
      GOTO 25     ! quit; error in reactions
   5  CONTINUE
C

C ***** FILL UP CHEMICAL PRODUCTION AND LOSS MATRICES *****
      DO 7 M=1,2
      N = 3-M          !so N=2, then 1
      DO 7 J=1,NR
      I = JCHEM(M,J)   !so I = JCHEM(1,NR) then JCEHM(2,NR)
      IF(I.LT.1.OR.I.GT.NSP) GO TO 7    !skips 0 (i.e. nothing) and NSP1 (HV)
      NUML(I) = NUML(I) + 1             !counter of how many reactions species I is involved with
      IF(NUML(I).GT.NMAX) GOTO 20    ! quit; too many reactions  (seems unnecesary, but whatever)
      K = NUML(I)
      ILOSS(1,I,K) = J               !ILOSS(1,species in pos 1, species reac#) then ILOSS(1,spec in pos 2, reac#)= global reaction #
      ILOSS(2,I,K) = JCHEM(N,J)      !ILOSS(1,species in pos 1, species reac#) then ILOSS(1,spec in pos 2, reac#)= other species
   7  CONTINUE
C
      DO 8 M=3,5
      DO 8 J=1,NR
      I = JCHEM(M,J)
      IF(I.LT.1.OR.I.GT.NSP) GO TO 8
      NUMP(I) = NUMP(I) + 1
      IF(NUMP(I).GT.NMAX) GO TO 20
      K = NUMP(I)
      IPROD(I,K) = J
   8  CONTINUE

c-mc check mass balance of chemical reactions

      do i=1,nr

      rhOcount=0
      rhHcount=0
      rhCcount=0
      rhScount=0
      rhNcount=0
      rhCLcount=0

      numprod=3    !assume 3 products unless..
      if (JCHEM(5,i).EQ.0) numprod=2
      if (JCHEM(4,i).EQ.0) numprod=1

      do j=0,numprod-1   !this loop counts up the mass on the right hand side of the .rx
         rhOcount=rhOcount+atomsO(JCHEM(3+j,i))
         rhHcount=rhHcount+atomsH(JCHEM(3+j,i))
         rhCcount=rhCcount+atomsC(JCHEM(3+j,i))
         rhScount=rhScount+atomsS(JCHEM(3+j,i))
         rhNcount=rhNcount+atomsN(JCHEM(3+j,i))
         rhCLcount=rhCLcount+atomsCL(JCHEM(3+j,i))
      enddo  

      bad=0
      if (rhOcount.ne.atomsO(JCHEM(1,i))+atomsO(JCHEM(2,i))) bad=1
      if (rhHcount.ne.atomsH(JCHEM(1,i))+atomsH(JCHEM(2,i))) bad=1
      if (rhCcount.ne.atomsC(JCHEM(1,i))+atomsC(JCHEM(2,i))) bad=1
      if (rhScount.ne.atomsS(JCHEM(1,i))+atomsS(JCHEM(2,i))) bad=1
      if (rhNcount.ne.atomsN(JCHEM(1,i))+atomsN(JCHEM(2,i))) bad=1
      if (rhCLcount.ne.atomsCL(JCHEM(1,i))+atomsCL(JCHEM(2,i))) bad=1

         if (bad .eq. 1) then
         print *, 'bad mass balance in reaction',i
         print *, (CHEMJ(j,i),j=1,5)
         print *, numprod
         !the problem is either in the .rx file or the species.dat file
         print *, rhNcount,atomsN(JCHEM(1,i)),atomsN(JCHEM(2,i))
         stop
         endif
      enddo   !end mass balance check


C PROCESS the photolysis reactions
c
c this next little bit creates:
c photoreac(kj) - an array of species numbers for each photolysis reaction. 
c                 used in absorbers/columndepth computations
c photospec(ks) - the unique photolysis reaction numbers (i.e. unique elements of photoreac)
c                 used in Initphoto.f to fill up sq, the cross section vector 
c photonums(kj) - the reaction number of each photolysis reaction
c                 used in Photo.f to fill up the A vector of rates


       testvec=INDEX(REACTYPE,'PHOTO')       
       !a vector of length nr with 1's in the location where photolysis reactions are

        jcount=1
        jrcount=1
        juniq=1
       do i=1,nr
        if (testvec(i).eq.1.) then
           photoreac(jrcount)=JCHEM(1,i)   !capture the species number of each photo reaction
           photonums(jrcount)=i             !capture the reaction number of each photoreaction

c           print *, jrcount,i,JCHEM(1,i),(CHEMJ(m,i),m=1,5)

           jrcount=jrcount+1
           
           !capture the unique photolysis species in photospec
           if (juniq.eq.1) then
              photospec(juniq)=JCHEM(1,i)
              juniq=juniq+1
           else
              bad=0
              do m=1,ks
               if (JCHEM(1,i).eq.photospec(m)) bad=1
              enddo   
              if (bad.eq.0) then
                photospec(juniq)=JCHEM(1,i)
                juniq=juniq+1
              endif
           endif
c          print *, jrcount,juniq,photospec(juniq-1)
      endif
           jcount=jcount+1

       enddo

c       print *, photoreac
c       print *, ''
c       print *, INT(photospec)
c       print *, juniq
c       print *, photonums
c       stop


C - SOME CHECKS TO MAKE SURE THE INPUT FILES JIVE WITH PARAMATERS.INC

       if (juniq-1.ne.ks) then
          print *, 'discrepency between unique photolysis reactions/ks'
          print *, juniq-1, ks
          stop
       endif

       if (SUM(INDEX(REACTYPE,'PHOTO')) .NE. kj) then  
          print *,'discrepency between number of photo reactions and kj'
          print *, SUM(INDEX(REACTYPE,'PHOTO')), kj
       stop
       endif   

      IF (ILL.NE.NQ.OR.iLL+iTD.NE.NQ1.OR.iLL+iTD+iSL+iIN.NE.NSP2) then
      print *, 'discrepancy between INPUTFILES/species.dat and
     $    INPUTFILES/parameters.inc'
      PRINT *, ILL,NQ, iLL+iTD,NQ1,iLL+iTD+iSL+iIN,NSP2
      stop
      endif


c-mc the below could be made into a nice printout in the out.out file if someone felt like it.
c       print *, SUM(INDEX(REACTYPE,'PHOTO'))
c       print *, SUM(INDEX(REACTYPE,'2BODY'))
c       print *, SUM(INDEX(REACTYPE,'3BODY'))
c       print *, SUM(INDEX(REACTYPE,'WEIRD'))
c       print *, NR


C
C ***** READ THE INPUT DATAFILE *****

c read in formatted input data file
c ISOHACK - these are just the sulfur species(ISOin.dist from main code)

      IROW = 10  !num columns in .dist file
      LR = NQ/IROW + 1      !NQ=80  -> 9
      RL = FLOAT(NQ)/IROW + 1  !9
      DIF = RL - LR  !0
      IF (DIF.LT.0.001) LR = LR - 1  !so LR=8

      DO L=1,LR
       K1 = 1 + (L-1)*IROW
       K2 = K1 + IROW - 1
       IF (L.LT.LR) then
        read(17, 880) ((USOL(k,i),K=K1,K2),i=1,nz) 
       ELSE
          K2 = NQ
          if (K2-K1 .EQ. 9) then   !this only occurs if NQ is a multiple of 10
            read(17, 880) ((USOL(k,i),K=K1,K2),i=1,nz) 
          else  !if not, need to generate a dynamic format statement
           fmtstr='(  E17.8)'
           write(fmtstr(2:3),'(I1)')K2-K1+1   !OK for one character as this should always be <10
           read(17, fmtstr) ((USOL(k,i),K=K1,K2),i=1,nz)
          endif    
       ENDIF
      enddo 

        read(17,881) (T(i),EDD(i),DEN(i),O3(i), SL(NSP-1,i),i=1,nz) !this final one is CO2 number density

        fmtstr='(  E17.8)'
        write(fmtstr(2:3),'(I2)')NP*3

        do i=1,nz
         read (17,fmtstr) (AERSOL(i,j),j=1,NP),(WFALL(i,j),j=1,NP),
     $                  (RPAR(i,j),j=1,NP) 
        enddo  


      if(USETD.EQ.1) then !particles in tri-diag
        fmtstr='(  E17.8)'
        write(fmtstr(2:3),'(I2)')NP
       do i=1,nz
        read(17,fmtstr)  (PARTICLES(i,j),j=1,np) 
       enddo
      endif


!ISOHACK
!(this loop multiples isotopic S2 by 2, S4 by etc.)
!also changes lower boundary flux for any species with multiple S (CS2)

      do j=1,nq
       if (atomsS(j).gt.1) then
         SGFLUX(j)=SGFLUX(j)*atomsS(j)          
         print *, 'multiplying lower boundary flux/mixing ratios for ' 
     $             ,ISPEC(j),' by ',atomsS(j)
        do i=1,nz
         USOL(j,i)=USOL(j,i)*atomsS(j)
        enddo   
       endif
      enddo

      if (USETD.EQ.1) then  
       do i=1,nz
         PARTICLES(i,2)=PARTICLES(i,2)*8.  !ack - hardcoded particle number
         !particle 1 in isotope code is SO4AER (1 S), particle 2 is S8AER (8 S's)
       enddo   
      endif

!ISOHACKS


c      print *, USOL(LSXS,1),USOL(LSXS7,1)

c      stop

c read in formatted input data file
c ISOHACK - this should be all species from the main code

      iCOL=iIN-2 !number of inert minus hv and M

      IROW = 10  !num columns in .dist file
      LR = (iCOL)/IROW + 1      
      RL = FLOAT(iCOL)/IROW + 1  
      DIF = RL - LR  !0
      IF (DIF.LT.0.001) LR = LR - 1  

      DO L=1,LR
       K1 = 1 + (L-1)*IROW
       K2 = K1 + IROW - 1
       IF (L.LT.LR) then
        read(55, 880) ((UINERT(k,i),K=K1,K2),i=1,nz) 
       ELSE
          K2 = iCOL
          if (K2-K1 .EQ. 9) then   !this only occurs if NQ is a multiple of 10
            read(55, 880) ((UINERT(k,i),K=K1,K2),i=1,nz) 
          else  !if not, need to generate a dynamic format statement
           fmtstr='(  E17.8)'
           write(fmtstr(2:3),'(I1)')K2-K1+1   !OK for one character as this should always be <10
           read(55, fmtstr) ((UINERT(k,i),K=K1,K2),i=1,nz)
          endif    
       ENDIF
      enddo 

c      do k=1,iIN-4
c      print *, UINERT(k,1)
c      enddo
c      stop
      
      Loff=iLL+iSL+iTD  !(needed to match up positions in UINERT with USOL) - number of non-inert (i.e. isotopic) species...

      do K=1,NQ
       VDEP(K) = VDEP0(K)
       VEFF(K) = VEFF0(K)
        do  I=1,NZ
         USOL(K,I) = ABS(USOL(K,I))
         USOLORIG(K,I)=USOL(K,I)
        enddo
      enddo
c

c finish setting boundary conditions:
!sgflux, vdep, smflux, and distributed fluxes are already set
        do i=1,nq
         if (LBOUND(i).EQ.1) USOL(i,1)=fixedmr(i)    
c         print *, ISPEC(i),USOL(i,1)
        enddo
c        stop

C ***** READ THE PLANET PARAMETER DATAFILE *****
      READ(7,502) G,FSCALE,ALB,ZTROP,FAR,R0,P0,PLANET,TIMEGA
 502  FORMAT(F7.1/,F7.2/,F7.3/,E7.1/,F7.3/,E8.3/,F8.3/,A8/,F4.2)





C ***** SET MODEL PARAMETERS *****
C     ZY = SOLAR ZENITH ANGLE (IN DEGREES)
C     LTIMES = COUNTER FOR PHOTORATE SUBROUTINE
C     DT = INITIAL TIME STEP
C     TSTOP = TIME AT WHICH CALCULATION IS TO STOP
C     NSTEPS = NUMBER OF TIME STEPS TO RUN (IF TSTOP IS NOT REACHED)
C     FO2 = ground level O2 mixing ratio used in heritage calculations

      LTIMES = 0
      ZY = 50.
      FO2 = UINERT(LO2-Loff,1)  !ISOHACK
      FCO2= UINERT(LCO2-Loff,1) !ISOHACK



      CALL GRID

      JTROP=minloc(Z,1, Z .ge. ztrop)-1 !height index for the tropopause (-1 given the staggered grid)
                                        !the 1 in the second postion tells minloc to return a scalar
      poop3=1.0  !Kevin's pressure diddling factor - need to make sure this still works...

      CALL DENSTY(FO2,poop3)
      CALL RATES
      CALL DIFCO(FO2)  !computes diffusion coefficents (K*N) and binary diffusion coefficents for H and H2
      CALL SATRAT(JTROP,H2O)
      CALL DOCHEM(FVAL,-1,JTROP,iIN,iSL,USETD)      !IDO=-1, fill up SL for accurate calculation on first timestep !ISOHACK

      if (PLANET .EQ. 'EARTH') then 
       PRONO = PRONO/1.  ! current column integrated NO production rate on Earth
                        ! divide by 1000 turns off ltning  
c       PRONO = PRONO/1.e6 ! ATACAMA
      else if (PLANET .EQ. 'MARS') then
       PRONO = PRONO/1.E9 ! divide by 1e9 turns off lightning for dry mars
      endif   

c$$$      if (mbound(LH2) .gt. 0) then  !i.e if constant mr or constant flux UBC
c$$$        do i=1,nz
c$$$            bHN2(i) = 0.0   ! don't use molecular diffusion
c$$$            bH2N2(i) = 0.0
c$$$c           bXN2(i) = 0.0   ! don't use molecular diffusion
c$$$        enddo
c$$$      else !  use effusion velocity formulation of diffusion limited flux 
c$$$        Veff(LH) = 1.0*bhN2(nz)/DEN(NZ)
c$$$     $     *(1./Hscale(nz) - 1./scale_H(LH,nz))   ! diff lim flux
c$$$        Veff(LH2) = 1.0*bH2N2(nz)/DEN(NZ)
c$$$     $     *(1./Hscale(nz) - 1./scale_H(LH2,nz))
c$$$      endif
      
c$$$
c$$$      do J=1,JTROP
c$$$       USOL(LH2O,J) = H2O(J)   !sets H2O to relative humidity in troposphere
c$$$      enddo

      CALL LTNING(FO2)
      CALL AERTAB   !makes table of vapor pressures for H2O and H2SO4
      NZ1 = NZ - 1
      HA = HSCALE(NZ)
      NRAIN = 0   ! count of calls to rainout
c orig      DT = 1.E-15   !why the hell is Kevin starting at the Planck time?
      DT = 1.E-6
      DTINV = 1./DT
      TIME = 0.

      TSTOP = 1.E17    !as it was...
      NSTEPS = 10000    !as it was...


      TIME = 1.E16   !isohack
      DT = TIME
      DTINV = 1./DT
      TSTOP = 1.E27   !isohack

!     temporarily unisohacking
c      NSTEPS =1  !ISOHACK
c      DTINV= 1e-40   !ISOHACK



C ***** write OUT INITIAL DATA *****
      CALL OUTPUT(0,NSTEPS,0.D0,jtrop, vdep,USOLORIG,USETD)


C
C ***** STORE CONSTANT JACOBIAN COEFFICIENTS *****
c-mc      DZ2 = DZ*DZ
      do i=1,nq1
        DU(i,1) = DK(1)/DEN(1)/DZ(1)**2
        DL(i,NZ) = DK(NZ1)/DEN(NZ)/DZ(NZ)**2
        DD(i,1) = DU(i,1)
        DD(i,NZ) = DL(i,NZ)
         do J=2,NZ1
          DU(i,J) = DK(J)/DEN(J)/DZ(J)**2
          DL(i,J) = DK(J-1)/DEN(J)/DZ(J)**2
   	  DD(i,J) = DU(i,J) + DL(i,J)
         enddo
      enddo

      do i=1,nq   ! first order molecular diffusion terms 
        do j=1,nz
           ADU(i,j) = 0.0
           ADL(i,j) = 0.0
           ADD(i,j) = 0.0
        enddo
      enddo

!below is ISOHACK commenting
c$$$c  for H and H2
c$$$c  lower boundary condition
c$$$
c$$$      if (mbound(LH2).eq.0) then   ! diff limited flux implemented as effusion velocity  
c$$$        DU(LH,1) = DU(LH,1) + bHN2(1)/Den(1)/DZ(1)**2
c$$$        ADU(LH,1) = bHN2(1)/Den(1)/DZ(1)/2.*
c$$$     6      (1./scale_H(LH,1)-1./H_atm(1))
c$$$        DU(LH2,1) = DU(LH2,1) + bH2N2(1)/Den(1)/DZ(1)**2
c$$$        ADU(LH2,1) = bH2N2(1)/Den(1)/DZ(1)/2.*
c$$$     6      (1./scale_H(LH2,1)-1./H_atm(1))
c$$$c upper boundary condition
c$$$        DL(LH,NZ) = DL(LH,NZ) + bHN2(nz1)/Den(nz)/DZ(NZ)**2
c$$$        ADL(LH,NZ) = -bHN2(nz1)/Den(nz)/DZ(nz)/2.*
c$$$     6      (1./scale_H(LH,nz1)-1./H_atm(nz1))
c$$$        DL(LH2,NZ) = DL(LH2,NZ) + bH2N2(nz1)/Den(nz)/DZ(NZ)**2
c$$$        ADL(LH2,NZ) = -bH2N2(nz1)/Den(nz)/DZ(NZ)/2.*
c$$$     6      (1./scale_H(LH2,nz1)-1./H_atm(nz1))
c$$$c  unused...
c$$$        DD(LH,1) = DU(LH,1)
c$$$        ADD(LH,1) = -ADU(LH,1)
c$$$        DD(LH2,1) = DU(LH2,1)
c$$$        ADD(LH2,1) = -ADU(LH2,1)
c$$$				
c$$$c interior grid points   ?fixed 8-13-05
c$$$        do j=2,nz1
c$$$            DU(LH,j) = DU(LH,j) + bHN2(j)/Den(j)/DZ(j)**2
c$$$            ADU(LH,j) = bHN2(j)/Den(j)/DZ(j)/2.*
c$$$     6            (1./scale_H(LH,j)-1./H_atm(j))
c$$$            DL(LH,j) = DL(LH,j) + bHN2(j-1)/Den(j)/DZ(j)**2
c$$$            ADL(LH,j) = -bHN2(j-1)/Den(j)/DZ(j)/2.*
c$$$     6            (1./scale_H(LH,j-1)-1./H_atm(j-1))
c$$$            DU(LH2,j) = DU(LH2,j) + bH2N2(j)/Den(j)/DZ(j)**2
c$$$            ADU(LH2,j) = bH2N2(j)/Den(j)/DZ(j)/2.*
c$$$     6             (1./scale_H(LH2,j)-1./H_atm(j))
c$$$            DL(LH2,j) = DL(LH2,j) + bH2N2(j-1)/Den(j)/DZ(j)**2
c$$$            ADL(LH2,j) = -bH2N2(j-1)/Den(j)/DZ(j)/2.*
c$$$     6             (1./scale_H(LH2,j-1)-1./H_atm(j-1))
c$$$            DD(LH,j) = DU(LH,j) + DL(LH,j)
c$$$            ADD(LH,j) = -ADU(LH,j) - ADL(LH,j)
c$$$            DD(LH2,j) = DU(LH2,j) + DL(LH2,j)
c$$$            ADD(LH2,j) = -ADU(LH2,j) - ADL(LH2,j)
c$$$        enddo
c$$$      endif  !end molecular diffusion for H and H2 loop
c$$$
c$$$
c$$$
c$$$      do I=1,nz
c$$$c         write(10, 1210) Z(I),DU(1,i),DU(LH,i),DU(LH2,i),DL(1,i),
c$$$c    5     DL(LH,i),DL(LH2,i),DD(nq,i),DD(LH2,i)
c$$$c          write(10, 1210) Z(I),scale_H(1,i),scale_H(LH,i),
c$$$c     5     scale_H(LH2,i), Hscale(i)
c$$$      enddo
c$$$ 1210 FORMAT(1X,1P12E10.2)

c  end of not yet ready to test section
C
      KD = 2*NQ + 1
      KU = KD - NQ
      KL = KD + NQ
C
C   write OUT RESULTS EVERY NPR TIME STEPS
      NPR = 50
      PRN = NPR
      
C   DO PHOTORATES EVERY MP TIME STEPS
      NPHOT = 0
      MP = 1        !integer since starts with M
      PM = MP       !float since starts with P
      NN = 0
C
C ***** START THE TIME-STEPPING LOOP *****
      DO 1 N=1,NSTEPS     !calculation is stopped if it hasn't converged by NSTEPS...
      TIME = TIME + DT
      NN = NN + 1                  !counter for number of timesteps
      MS = (N-1)/MP                !both integers, so answer is modular (0 for N=1-3 , 1 for N=4-6, etc.)
      SM = (N-1)/PM                !SM=0,1/3,2/3,1,4/3,5/3,2......
      IF (NN.EQ.NSTEPS) SM = MS
      IF (SM-MS.GT.0.01) goto 18   ! skip PHOTO
      IF (N.GT.1 .AND. TIME.LT.1.E2) goto 18  !skip PHOTO    !do photo on first time step, not again until 100 seconds


! store mixing ratio of all species that take place in photolysis reactions
! these are the absorbers which block solar radiation (this exclude SL species - check this assumption at some poin)
!or does it? - 

!OK - stopping here for now.  check carefully below to see that photoreac, photospec
!are OK, and probably need to allow for absorbtion again...

c      print *, LSXO, LSO, Loff,LSO-Loff
c      print *, USOL(LSXO,1),UINERT(LSO-Loff,1)

c      print *, photoreac
c      print *, ISPEC(INT(photoreac))
c      print *, Loff, iLL
c      stop
                   isooffset=44    !ACK - hardcoding...  this should really be fixed at some point
                   !this is number of non-isotopic major species, and is not possible to calculate from what ISO.f currently knows
                   !it should be read from species.dat, or passed, or something
                   !I think NQ from parametersREGULAR-iLL is what is really wanted here (also minus number of non-iso particles)
                   !after the Sorgs main NQ was 42 and iLL=22
                   !after the HC's main NQ was 67 and iLL=22  (but what about HCAER? - do I need to subtract and extra 1 here?)
                   !after HC's without SORG: main NQ was 61 and ill=15 - so 46? NO - needs to be 44 - due to HCAERs?
                   !tried changing to 45 after adding CO2 to the main loop - doesn't work THIS REALLY NEEDS TO BE ABSTRACTED - ACK ACK ACK 

         lpolyscount=1
      do k=1,kj

       do i=1,nz

          if (photoreac(k).gt.nq) then   !this gets any SL/IN species that have photorates
!but, this isotope code is a bit different. All of the "major" species are in the inert loop.

               where= index(ISPEC(INT(photoreac(k))),'S')  !captures INERT species that don't contain S
               !INERT species containing S should not have PHOTO reactions in ISOreactions.rx, so this will nab all inert major species

               if (where.eq.0) then
                absorbers(k,i)=UINERT(INT(photoreac(k))-Loff,i)    !ISOHACK    
c                if(i.eq.1) print *, k,ISPEC(INT(photoreac(k))),'HELLO1'
!it's OK to use UINERT for the ISOTOPE code because the densitites don't change anyway.  This will get CO2 as well
!UINERT is being indexed by photoreac-Loff, which basically returns the order of ISPEC in the original code
               endif 

           else if (ISPEC(INT(photoreac(k))).eq.'S8      '    !be crafty about S8
     $          .or.ISPEC(INT(photoreac(k))).eq.'SXS7   '   ) then    !quasi-hardcoded S8 behavior


!              absorbers(k,i)=ABS(USOL(INT(photoreac(k)),i))     !this is 8 times as much as the original in the isotope code
!              if(i.eq.1)print *, absorbers(k,i)

              absorbers(k,i)=ABS(UINERT(INT(photoreac(k))+isooffset,i)) !these are the original unperterbed species                  
c              if(i.eq.1)print *, absorbers(k,i),INT(photoreac(k))
c              if(i.eq.1) print *, k,ISPEC(INT(photoreac(k))),'HELLO3'

              if (lpolyscount .eq. 2) absorbers(k,i)=0.0 !S8R doesn't really exist
              !S8L doesn't really either, but this is how we are rolling for now...
              !we are filling up absorbers for S8R and S8, but S8 doesn't have a cross section
              !so only S8R is actually absorbing photons here in the RT scheme.
              if (i.eq.nz) then
c            print *, k,photoreac(k),ISPEC(INT(photoreac(k))),lpolyscount
                 lpolyscount=lpolyscount+1
              endif   
          else     !fill in rest of long-livied species (which are only S bearing species)
c                absorbers(k,i)=ABS(USOL(INT(photoreac(k)),i))  !computed values with enhanced densities
!what I need is for absorbers S2, S3, S4, and S8 to use the UINERT values
c                if(i.eq.1) then
c                   print *, k,ISPEC(INT(photoreac(k))),'HELLO4'  
c                 endif    

               absorbers(k,i)=ABS(UINERT(INT(photoreac(k))+isooffset,i))  !using inert values for RT

          endif 
       enddo                                           
      enddo


c      do k=1,kj
c         print *, photoreac(k),ISPEC(INT(photoreac(k))),absorbers(k,1)
c      enddo   

c      stop

!fill up some heritage vectors rather than try to extract them from the code...
! I am assuming here that the code will always have these species in it so no if's are needed
c$$$
c$$$       JH2O=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'H2O    ')
c$$$       JO2_O1D=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'O2    ')
c$$$       JO3_O1D=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'O3    ')
c$$$!the above finds the first entry in photoreac for the given speices
c$$$
c$$$      do I=1,NZ
c$$$       H2O(I) = absorbers(JH2O,I)   
c$$$       O2(I) =  absorbers(JO2_O1D,I)
c$$$       O3(I) = absorbers(JO3_O1D,I)
c$$$       CO2(I) = FCO2
c$$$      enddo

      IDO = 0
      IF (NN.EQ.NSTEPS) IDO = 1
      CALL PHOTO(ZY,AGL,LTIMES,ISEASON,IZYO2,IO2,INO,IDO,timega,frak) 
      CALL RAINOUT(JTROP,NRAIN,USETD)
      CALL AERCON
      
c      print *, photoreac
c      print *, 'stopping after photo in Iso'
c      stop




C
C   TIME-DEPENDENT BOUNDARY CONDITIONS 
C     (NOTE THAT THE INCLUSION OF
C     A FALL VELOCITY FOR PARTICLES IS MATHEMATICALLY EQUIVALENT TO
C     INCREASING THE DEPOSITION VELOCITY AT THE SURFACE.  AT THE TOP
C     BOUNDARY, ZERO FLUX IS ACHIEVED BY SETTING THE EFFUSION VELO-
C     CITY EQUAL TO THE FALL VELOCITY)

      if(USETD.EQ.0) then
      VDEP(LSXO4AER) = VDEP0(LSXO4AER) + WFALL(1,1)   !ISOHACK
      VEFF(LSXO4AER) = VEFF0(LSXO4AER) + WFALL(NZ,1)  !ISOHACK

      VDEP(LSXS7AER) = VDEP0(LSXS7AER) + WFALL(1,2)   !ISOHACK
      VEFF(LSXS7AER) = VEFF0(LSXS7AER) + WFALL(NZ,2)  !ISOHACK
      endif

       JCO2=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'CO2    ')
       JCO2_O1D = JCO2+1
      VCO2 = (prates(JCO2,NZ) + prates(JCO2_O1D,NZ)) * HA
c      SMFLUX(LO) = - VCO2*CO2(NZ)*DEN(NZ) !ISOHACK
c      SMFLUX(LCO) = SMFLUX(LO)            !ISOHACK


      NMP = NSTEPS - MP
      if ((NN/NPR)*NPR.eq.NN .or. NN.eq.1) then

C  PHOTOLYSIS RATES FORMATTED I/O AND PRINTOUT

      IROW = 8  !num columns of photorates
      LR = KJ/IROW + 1      
      RL = FLOAT(KJ)/IROW + 1 
      DIF = RL - LR  
      IF (DIF.LT.0.001) LR = LR - 1 

      write(14,885)
      DO L=1,LR
       K1 = 1 + (L-1)*IROW
       K2 = K1 + IROW - 1
       IF (L.lt.LR) then
        write(14, 884) (photolabel(K),K=K1,K2) 
        write(14, 883) (Z(I),(prates(k,i),K=K1,K2),i=1,nz,3) 
       ELSE
          K2 = kj
           fmtstr="(/5X,'Z',6X,  A11)"
           write(fmtstr(13:14),'(I1)')K2-K1+1   !OK for one character as this should always be <10
           write(14, fmtstr) (photolabel(K),K=K1,K2) 

           fmtstr2='(  (1PE9.2,2X))'
           write(fmtstr2(2:3),'(I2)')K2-K1+2   !+2 here so it writes the Z a well   
           write(14, fmtstr2) (Z(I),(prates(k,i),K=K1,K2),i=1,nz,3) 
       ENDIF
      enddo 

 883  format(9(1PE9.2,2X))
 884  format(/5X,'Z',6X,8A11)
 885  FORMAT(//1X,'PHOTOLYSIS RATES')

      endif !end photorates printout

 18   continue   !start here if we are skipping PHOTO
C
      CALL SEDMNT(FSULF,USETD,frak)  !this subroutine updates RPAR and WFALL, which are skiped in the ISOTOPE code
! it also computes CONVER, which is needed


      if (USETD.EQ.0) then   !particles in main loop

      do J=1,NZ

!ISOHACK QUESTION - revisit this - perhaps this should not be touched...
c      AERSOL(J,1) = USOL(LSXO4AER,J)*DEN(J)/CONVER(J,1)  !ISOHACK
c      AERSOL(J,2) = USOL(LSXS7AER,J)*DEN(J)/CONVER(J,2)  !ISOHACK
c      print *, AERSOL(J,1),AERSOL(J,2),USOL(LSO4AER,J),CONVER(J,1),
c     &         USOL(LS8AER,J),CONVER(J,2)
      enddo


C
C   COMPUTE ADVECTION TERMS FOR PARTICLES  ! jim is using centered differences.
c   this makes sense for the inner points from the differential equation,
c   where he is taking a centered difference.
c   i am not so sure this applies at the boundaries
c   I need to find where he has done the boundaries
c   it may be that I should instead change flux estimates to use the "2" 
      DO 38 J=1,NP
      DPU(1,J) = WFALL(2,J)*DEN(2)/DEN(1)/(2.*DZ(1))
      DPL(NZ,J) = WFALL(NZ1,J)*DEN(NZ1)/DEN(NZ)/(2.*DZ(NZ))
      DO 38 I=2,NZ1
      DPU(I,J) = WFALL(I+1,J)*DEN(I+1)/DEN(I)/(2.*DZ(I))
  38  DPL(I,J) = WFALL(I-1,J)*DEN(I-1)/DEN(I)/(2.*DZ(I))
      
      else  !particles in tri-diag

      do J=1,NZ
       do k=1,NP
        AERSOL(J,K)=PARTICLES(J,K)*DEN(J)/CONVER(J,K)  !ISOHACK - new return after reviewing below
       enddo                                           !which has an if isotope.eq.)
      enddo


      endif   




C ***** SET UP THE JACOBIAN MATRIX AND RIGHT-HAND SIDE *****

c-mc gonna start here and try no to mess anything up - point is to feedback the solution on itself
c-mc first go is to try to just repeat a few times - if this works, then set up a check for convergence

c-mc OK seems to work.  4 iterations seems to provide a converged solution.  
c-mc test 1  - check the solution with 4 iterations against that of feeding back the code 4 times...
c-mc  this test looks good - see 'changetesting' sheet in ~/keff/redox.xls
c-mc next up is to save some intermediate output to see how much things are changing in subsquent iterations
c-mc  I should try to identify an epsilon that satisfies a similar global condition to that of just repeating
c-mc the loop four times...

c-MC - turning off iterated jacobians for chlorine testing
c-MC turning back on.../off


      DO 73 ITERATE=1,1

      DO 17 J=1,LDA
      DO 17 K=1,NEQ
  17  DJAC(J,K) = 0.
      DO 19 K=1,NEQ
  19  RHS(K) = 0.
C
C     (DJAC IS EQUAL TO (1/DT)*I - J, WHERE J IS THE JACOBIAN MATRIX)
c-mc  NOTE - DJAC is in band storage form.  see sgbfa header for details.
c-mACK expand me...

C - in ISOHACK mode DT-> INF so remove the 1/DT term and DJAC = -J


C
C   COMPUTE CHEMISTRY TERMS AT ALL GRID POINTS
      IDO = 0
      if ((NN/NPR)*NPR.eq.NN) IDO = 1  !computes TP, TL
      IF (NN.EQ.NSTEPS) IDO = 1

      CALL DOCHEM(FVAL,IDO,JTROP,iIN,iSL,USETD)      !IDO=1 happens on last step- computes total production and loss...

      DO 9 I=1,NQ
      DO 9 J=1,NZ
      K = I + (J-1)*NQ
      RHS(K) = FVAL(I,J)
      if (ITERATE.eq.1) USAVEOLD(I,J) = USOL(I,J)   
!mc testing 4/29/06  used to revert if timestep is too big.
!mc and also for USOLPREV once the timestep succeedes
   9  USAVE(I,J) = USOL(I,J)       !original code  - used as part of the reverse euler solver

C
      DO 3 I=1,NQ
      DO 11 J=1,NZ
c     R(J) = EPSJ * ABS(USOL(I,J))   !as it was
      R(J) = EPSJ * USOL(I,J)
  11  USOL(I,J) = USAVE(I,J) + R(J)
      CALL DOCHEM(FV,0,JTROP,iIN,iSL,USETD)

C
      DO 12 M=1,NQ
      MM = M - I + KD
      DO 12 J=1,NZ
      K = I + (J-1)*NQ
  12  DJAC(MM,K) = (FVAL(M,J) - FV(M,J))/R(J)      !-J since its orig - perturbed
C


      DO 10 J=1,NZ
  10  USOL(I,J) = USAVE(I,J)
   3  CONTINUE



C
C   COMPUTE TRANSPORT TERMS AT INTERIOR GRID POINTS
      DO 13 I = 1,NQ
      DO 14 J=2,NZ1
      K = I + (J-1)*NQ
      RHS(K) = RHS(K) - DD(i,J)*USOL(I,J) - ADD(i,j)*USOL(I,J)
     1  + DU(i,J)*USOL(I,J+1) + ADU(i,j)*USOL(i,j+1)
     2  + DL(i,J)*USOL(I,J-1) + ADL(i,J)*USOL(I,J-1)
      DJAC(KD,K) = DJAC(KD,K) + DTINV + DD(i,J) + ADD(i,j)
      DJAC(KU,K+NQ) = - DU(i,J) - ADU(i,j)
  14  DJAC(KL,K-NQ) = - DL(i,J) - ADL(i,j)
  13  CONTINUE
C

c-mc  ok, I need to verify that this is adding in -J in all DJAC calls (check signs above?)
c-mc ok - Jacobian for transport diagonals is: J~Chem-DD. We want -J and have already 
c-mc filled with -CHEM, so adding DD is appropriate. DTINV is the extra term in the main diagonal
c-mc J_upper=DU and J_lower=DL, so DJAC (which is -J)  uses -DU and -DL respectivly 


      if(USETD.EQ.0) then  !particles in main loop
!ack - these need to be abstracted...
C   ADD ADVECTION TERMS FOR PARTICLES
      L = 1
      I = LSXO4AER     !ISOHACK
C
      DO 24 J=2,NZ1
      K = I + (J-1)*NQ
      RHS(K) = RHS(K) + DPU(J,L)*USOL(I,J+1) - DPL(J,L)*USOL(I,J-1)
      DJAC(KU,K+NQ) = DJAC(KU,K+NQ) - DPU(J,L)
  24  DJAC(KL,K-NQ) = DJAC(KL,K-NQ) + DPL(J,L)

      L = 2
      I = LSXS7AER   !ISOHACK
      do J=2,NZ1
        K = I + (J-1)*NQ
        RHS(K) = RHS(K) + DPU(J,L)*USOL(I,J+1) - DPL(J,L)*USOL(I,J-1)
        DJAC(KU,K+NQ) = DJAC(KU,K+NQ) - DPU(J,L)
        DJAC(KL,K-NQ) = DJAC(KL,K-NQ) + DPL(J,L)
      enddo

      endif
C
C ***** LOWER BOUNDARY CONDITIONS *****
      DO 15 K=1,NQ
      U(K) = USOL(K,1)
      LB = LBOUND(K)   !ok as long as we don't model atmospheric Boron

      if (LB.eq.0 .OR. LB.eq.3) then
C       CONSTANT DEPOSITION VELOCITY/SPECIES WITH DISTRIBUTED FLUXES
        RHS(K) = RHS(K) + DU(k,1)*(USOL(K,2) - U(K)) 
     2   + ADU(k,1)*(USOL(K,2) + U(K)) - VDEP(K)*U(K)/DZ(1)
        DJAC(KD,K) = DJAC(KD,K) +DTINV +DU(k,1) -ADU(k,1) +VDEP(K)/DZ(1)
        DJAC(KU,K+NQ) = - DU(k,1) - ADU(k,1)  
c  is this right for particles??

      else if (LB .eq. 1) then
C       CONSTANT MIXING RATIO
        RHS(K) = 0.
        do M=1,NQ
          MM = KD + K - M
          DJAC(MM,M) = 0.
        enddo
        DJAC(KU,K+NQ) = 0.
        DJAC(KD,K) = DTINV + DU(k,1) - ADU(k,1)
      else  
C       CONSTANT UPWARD FLUX
        RHS(K) = RHS(K) + DU(k,1)*(USOL(K,2) - U(K)) 
     2   + ADU(k,1)*(USOL(K,2) + U(K)) + SGFLUX(K)/DEN(1)/DZ(1)
        DJAC(KD,K) = DJAC(KD,K) + DTINV + DU(k,1) - ADU(k,1)
        DJAC(KU,K+NQ) = - DU(k,1) - ADU(k,1)
      endif
 
 15   continue


C
C
C ***** UPPER BOUNDARY CONDITIONS *****
      DO 30 I=1,NQ
      U(I) = USOL(I,NZ)
      K = I + NZ1*NQ
      MB = MBOUND(I)
C
      if (MB.eq.0) then
C       CONSTANT EFFUSION VELOCITY
        RHS(K) = RHS(K) + DL(i,NZ)*(USOL(I,NZ1) - U(I))
     2  + ADL(i,NZ)*(USOL(I,NZ1) + U(I)) - VEFF(I)*U(I)/DZ(NZ)
       DJAC(KD,K) = DJAC(KD,K) +DTINV +DL(i,NZ) -ADL(i,NZ)
     2  + VEFF(I)/DZ(NZ)
        DJAC(KL,K-NQ) = - DL(i,NZ) -ADL(i,NZ)
      else if (MB.eq. 1) then
c       constant mixing ratio at the top. not debugged
        RHS(K) = 0.
        do M=1,NQ
          MM = KD + K - M
          DJAC(MM,M) = 0.
        enddo
        DJAC(KU,K+NQ) = 0.
        DJAC(KD,K) = DTINV + DL(i,NZ) -ADL(i,NZ)
      else
C   CONSTANT UPWARD FLUX
        RHS(K) = RHS(K) + DL(i,NZ)*(USOL(I,NZ1) - U(I)) 
     2   + ADL(i,NZ)*(USOL(I,NZ1) + U(I)) - SMFLUX(I)/DEN(NZ)/DZ(NZ)
        DJAC(KD,K) = DJAC(KD,K) + DTINV + DL(i,NZ) - ADL(i,NZ)
        DJAC(KL,K-NQ) = - DL(i,NZ) - ADL(i,NZ)
      endif
 30   continue
 
C   HOLD H2O AND S8 CONSTANT BELOW ZTROP
c   why am I doing this for S8??
c  turn it off for S8
      if(ISOTOPE.EQ.0) then    !ISOHACK - skip
      DO 34 I=1,1   ! Jim apparently was prepared to do this for many species
        L = LH2O
c        IF (I.EQ.2) L = LS8
        DO 33 J=1,JTROP
          K = L + (J-1)*NQ
          RHS(K) = 0.
          DO 32 M=1,NQ
          MM = M - L + KD
          DJAC(MM,K) = 0.
  32    continue
        DJAC(KD,K) = DTINV
        DJAC(KU,K+NQ) = 0.
        IF(J.EQ.1) GO TO 33
        DJAC(KL,K-NQ) = 0.
  33    continue
  34  CONTINUE
      endif   !ISOHACK - end skip loop

C distributed (volcanic) sources

      do i=1,nq
       if (LBOUND(i).eq.3) then
        disth=distheight(i)*1.e5  !convert to cm  
        jdisth=minloc(Z,1, Z .ge. disth)-1 !height index (-1 given the staggered grid)
                          !the 1 in the second postion tells minloc to return a scalar
        ZTOP=Z(jdisth)-Z(1)  
        ZTOP1=Z(jdisth)+0.5*DZ(jdistH)
c        print *, ISPEC(i),distH,jdistH,ZTOP,ZTOP1

        do j=2,jdisth !distribute from second level to distheight
          K = i + (J-1)*NQ
          rhs(k) = rhs(k) + 2.*distflux(i)*(ZTOP1-Z(j))/(Den(j)*ZTOP**2)
        enddo
       endif   
      enddo
      
C ***** FACTOR THE JACOBIAN AND SOLVE THE LINEAR SYSTEM *****

      CALL SGBFA(DJAC,LDA,NEQ,NQ,NQ,IPVT,INFO)
      IF(INFO.NE.0) then
c         print 103, N,INFO
         print *, N,INFO
         print *, 'ERROR in SGBFA'
         stop
      endif   
c 103  FORMAT(/1X,'N =',I3,5X,'INFO =',I3)


c-mc  we have set up RHS so that DJAC*X=RHS
c-mc DJAC is now upper triangular...

      CALL SGBSL(DJAC,LDA,NEQ,NQ,NQ,IPVT,RHS,0)

c-mc  after this point, RHS changes to solution vector for DJAC*X = RHS
c-mc  i.e. "RHS" = X = f_{n+1} - f_{n}, so that f_{n+1} = f{n} + RHS

      J15=minloc(Z,1, Z/1.e5 .ge. 15)-1 !height index (-1 given the staggered grid)
      J25=minloc(Z,1, Z/1.e5 .ge. 25)-1 !height index (-1 given the staggered grid)
      J70=minloc(Z,1, Z/1.e5 .ge. 70)-1 !height index (-1 given the staggered grid)
      J50=minloc(Z,1, Z/1.e5 .ge. 50)-1 !height index (-1 given the staggered grid)

C
C   COMPUTE NEW CONCENTRATIONS (IGNORE ERRORS IN SEVERAL SPECIES
C     THAT VIRTUALLY DISAPPEAR UP HIGH)
      EMAX = 0.
      DO 26 I=1,NQ
      DO 26 J=1,NZ
        K = I + (J-1)*NQ
c        IF (I.EQ.LH2S .AND. J.GT.J25) THEN  ! 30        
c          USOL(I,J) = USOL(I,J) + RHS(K)
c        ELSEIF (I.EQ.LS2 .AND. USOL(I,J).LT.1.E-20) THEN
c          USOL(I,J) = USOL(I,J) + RHS(K)
c        ELSEIF (I.EQ.LS4 .AND. USOL(I,J).LT.1.E-20) THEN
c          USOL(I,J) = USOL(I,J) + RHS(K)
c        ELSEIF (I.EQ.LS8 .AND. USOL(I,J).LT.1.E-20) THEN
c          USOL(I,J) = USOL(I,J) + RHS(K)
corig        ELSEIF(I.EQ.LSO4AER .AND. J.GT.J25) THEN  ! 50, this often causes problems causes the program to fail at 49.5 km!
c        IF (I.EQ.LSO4AER .AND. J.GT.J25) THEN  ! 50, this often causes problems causes the program to fail at 49.5 km!

!ACK - should return to this now that particles condensation is fixed up (ISOHACK AS well)
!should somehow enforce a check that this is the same as TOTCtester.f  (or just do the damn port...)

        IF (I.EQ.LSXO4AER.AND. J.GT.J25) THEN  ! a less drastic measure
          USOL(I,J) = USOL(I,J) + RHS(K)
        ELSEIF(I.EQ.LSXS7AER.AND. J.GT.J25) THEN  ! !!!
          USOL(I,J) = USOL(I,J) + RHS(K)

c  Jim set this at 50 km.  program fails at 49.5 km.
c  if I turn it off the program fails at 63.5 km immediately
c   so I may conclude that there is an issue with the UBC on SO4AER
c  anyway, I'll set it to say 35 km
c
c   more generally I want to ignore errors in anything with
c   mixing ratios less than say 1e-20
c       ELSEIF(I.EQ.LNO .AND. J.GT.70) THEN  ! 50, this causes the program to fail at 49.5 km!
c         USOL(I,J) = USOL(I,J) + RHS(K)
c       ELSEIF(I.EQ.LNO2 .AND. J.GT.70) THEN  ! 50, this causes the program to fail at 49.5 km!
c         USOL(I,J) = USOL(I,J) + RHS(K)
c       ELSEIF(I.EQ.LHNO .AND. J.GT.70) THEN  ! 50, this causes the program to fail at 49.5 km!
c         USOL(I,J) = USOL(I,J) + RHS(K)

        ELSEIF(I.EQ.LSXS3) THEN  !isohack
                     USOL(I,J) = USOL(I,J) + RHS(K)


        ELSEIF (USOL(I,J).LT. 1.E-20) THEN
c        IF (USOL(I,J).LT. 1.E-20) THEN
          USOL(I,J) = USOL(I,J) + RHS(K)
c  
        ELSE
          REL(I,J) = RHS(K)/USOL(I,J)
          EREL = ABS(REL(I,J))
          EMAX = max(EMAX,EREL)
          IF(EREL.LT.EMAX) THEN
            USOL(I,J) = USOL(I,J) + RHS(K)
          ELSE            !store info on species with largest error
            IS = I
            JS = J        !mc -this label is OK, because S will never have a photolysis reaction 
            UMAX = USOL(I,J)
            RMAX = RHS(K)
            USOL(I,J) = USOL(I,J) + RHS(K)
          ENDIF
        ENDIF
26    CONTINUE



C
C   RESET TROPOSPHERIC H2O TO ITS ORIGINAL VALUES, IN CASE IT CHANGED.
C     (IT SHOULDN'T HAVE, BUT IT MIGHT.)  ALSO, SULFATE AEROSOL HAS A
C     TENDENCY TO GO NEGATIVE NEAR THE UPPER BOUNDARY, SO MAKE SURE
C     IT DOESN'T STAY THAT WAY.)
c-mc I don't get this. If the code works right, it shouldn't change.
c-mc and conversly if the code doesn't work right, it should be fixed...
c-mc test this at some point down the road

      if(ISOTOPE.EQ.0) then
      DO 4 J=1,JTROP
        USOL(LH2O,J) = H2O(J)
c       USOL(LS8,J) = S8S(J)
   4  CONTINUE
      endif
c      do i=1,nq
c       do j=1,nz
c          print *,Z(j), USOL(LSO4AER,j),USOL(LS8AER,j)
c         USOL(i,j)=abs(USOL(i,j)) 
c       enddo
c      enddo 
      !temp

!ACK - needs to be something here to deal with names

!diddle particles
      smallest = 1.d-99
      lcountso4=0
      lcounts8=0

      DO J=1,NZ
      if(USOL(LSXO4AER,J).LT.0) lcountso4=lcountso4+1

      if(lcountso4.gt.0) then
       if (J.GT.1) then
       USOL(LSXO4AER,J)=USOL(LSXO4AER,J-1)*EXP(-WFALL(J,1)*DZ(J)/EDD(J))
       else
        USOL(LSXO4AER,J)=-1
       endif 
      endif 

      if(USOL(LSXS7AER,J).LT.0) lcounts8=lcounts8+1 
      if(lcounts8.gt.0) then
       if(J.GT.1) then
       USOL(LSXS7AER,J) = USOL(LSXS7AER,J-1) * 
     &   EXP(-WFALL(J,2)*DZ(J)/EDD(J))
       else
       USOL(LSXS7AER,J)=-1
       endif   
      endif

      USOL(LSXO4AER,J) = max (USOL(LSXO4AER,J),smallest)
      USOL(LSXS7AER,J) = max(USOL(LSXS7AER,J),smallest) 
      enddo  

      !test
      do i=1,nq
         do j=1,nz
c           USOL(i,j)=max(USOL(i,j),smallest) 
           USOL(i,j)=abs(USOL(i,j)) 
         enddo  
      enddo   


      if(USETD.EQ.1) then

*********TRIDIAG STARTS HERE  

C ***** SOLVE FOR S8 AND SO4 PARTICLES USING A TRIDIAGONAL INVERSION *****
c-mc I haven't done the work to port the hc aerosols to the tri-diag.  this would need to be done if desired.

      DO 58 L=1,NP  
      I = NQ + L   !tridiag particles must appear right after LL species in species.dat

      IF(I.EQ.LSO4AER) MZ = minloc(Z,1, Z/1.e5 .ge. 50)-1 
      IF(I.EQ.LSXO4AER) MZ = minloc(Z,1, Z/1.e5 .ge. 50)-1  !ISOHACK
      IF(I.EQ.LS8AER) MZ = minloc(Z,1, Z/1.e5 .ge. 40)-1 
      IF(I.EQ.LSXS7AER) MZ = minloc(Z,1, Z/1.e5 .ge. 40)-1  !ISOHACK
      !at some point check/abstract these 40/50km assumptions - this could be easily shunted to the species.dat file...
      !height indexes are -1 given the staggered grid

!the above is OK for the isotope loop because the main and isotopic species shouldn't be listed as TD in a given species.dat file
!I might need to rethink this if the densities really are enhanced (i.e. SXS7 = 8*S8) - although these all seem to be ratioed...

      MZ1 = MZ - 1
      MZP1 = MZ + 1

C   COMPUTE ADVECTION TERMS FOR PARTICLES
      DPU(1,L) = WFALL(2,L)*DEN(2)/DEN(1)/(2.*DZ(1))
      DPL(NZ,L) = WFALL(NZ1,L)*DEN(NZ1)/DEN(NZ)/(2.*DZ(NZ))
      DO 381 J=2,NZ1
      DPU(J,L) = WFALL(J+1,L)*DEN(J+1)/DEN(J)/(2.*DZ(J))
 381  DPL(J,L) = WFALL(J-1,L)*DEN(J-1)/DEN(J)/(2.*DZ(J))
    ! jim is using centered differences.
c   this makes sense for the inner points from the differential equation,
c   where he is taking a centered difference.
c   i am not so sure this applies at the boundaries
c   I need to find where he has done the boundaries
c   it may be that I should instead change flux estimates to use the "2" 

c      print *, DPU

C
C   TA = LOWER DIAGONAL, TB = DIAGONAL, TC = UPPER DIAGONAL, TY =
C   RIGHT-HAND SIDE
      DO 70 J=1,NZ
      TA(J) = 0.
      TB(J) = 0.
      TC(J) = 0.
  70  TY(J) = 0.
C
      DO  J=1,MZ
      TB(J) = YL(I,J)
      TY(J) = YP(I,J)/DEN(J)
c      print *,ISPEC(I),J,YL(I,J),YP(I,J)/DEN(J)
      enddo
C
      DO 45 J=2,MZ1
      TA(J) = - DL(I,J) + DPL(J,L)
      TB(J) = TB(J) + DD(I,J)
  45  TC(J) = - DU(I,J) - DPU(J,L)
C

! why are there no dl*PARTICLES() in here?  all the other dl's are multiplied by USOL.
!just on RHS...

      vturb=0.01
c      vturb=0.000005
c      vturb=0.0

C
C   BOUNDARY CONDITIONS
      TA(MZ) = - DL(I,MZ) + DPL(MZ,L)
      TB(MZ) = TB(MZ) + DL(I,MZ) + 0.5*WFALL(MZ,L)/DZ(MZ)
c      TB(1) = TB(1) + DU(I,1) + (0. - 0.5*WFALL(1,L))/DZ(1) !orig from Jim's code, as it was..
c      TB(1) = TB(1) + DU(I,1) + (vturb - 0.5*WFALL(1,L))/DZ(1) !orig from Jim's code, as it was..
      TB(1) = TB(1) + DU(I,1) - (vturb+0.5*WFALL(1,L))/DZ(1) !testing...
!0.01 is turbulent deposition velocity (see Pavlov 02 appendix)
      TC(1) = - DU(I,1) - DPU(1,L)
C
      NFLAG=0
      CALL SGTSL(MZ,TA,TB,TC,TY,NFLAG)
      IF (NFLAG.NE.0) PRINT 400, N,NFLAG,I
 400  FORMAT(//1X,"TRIDIAGONAL SOLVER FAILED AT N =",I3,2X,
     2  "NFLAG =",I2,2X,"SPECIES #",I2)
C

      do J=1,mz
       PARTICLES(J,L)=ABS(TY(J)) !is the abs really needed here?   
      enddo   

      smallest = 1d-99

C   FILL UP UPPER PORTION WITH APPROXIMATE ANALYTIC SOLUTION

      do J=MZP1,NZ
        PARTICLES(J,L)=PARTICLES(J-1,L) * EXP(-WFALL(J,L)*DZ(J)/EDD(J))
        PARTICLES(J,L) = MAX(PARTICLES(J,L),smallest)   !needed?
      enddo   



   58 CONTINUE
C
      endif  !end tri-diag loop
c      print *, 'stopping below tri-diag'
c      stop


 73   continue           !continue doing newton steps (4 seems to work best)
                         !someday I should see if this is justified 
                         !(r37:36 in /td branch has first attempt at convergence checking)



c      print *, EMAX
c      print *, USOL
c      print *, 'stopping at automatic time step control'



C   AUTOMATIC TIME STEP CONTROL
      DTSAVE = DT
c-mc      these are the ones that kevin originally used...
c      IF(EMAX.LT.0.2)  DT = 1.1*DTSAVE
c      IF(EMAX.LT.0.1)  DT = 1.2*DTSAVE
c      IF(EMAX.LT.0.04)  DT = 1.4*DTSAVE
c      IF(EMAX.LT.0.02)  DT = 2.0*DTSAVE
c      IF(EMAX.LT.0.01)  DT = 3.0*DTSAVE
c      IF(EMAX.LT.0.003) DT = 4.0*DTSAVE
c      IF(EMAX.LT.0.001) DT = 5.*DTSAVE

c-mc      these are even stricter...
      IF(EMAX.LT.0.15)  DT = 1.1*DTSAVE
      IF(EMAX.LT.0.07)  DT = 1.2*DTSAVE
      IF(EMAX.LT.0.01)  DT = 1.4*DTSAVE
      IF(EMAX.LT.0.008)  DT = 2.0*DTSAVE
      IF(EMAX.LT.0.004)  DT = 3.0*DTSAVE
      IF(EMAX.LT.0.001) DT = 4.0*DTSAVE
      IF(EMAX.LT.0.0005) DT = 5.*DTSAVE


c-mc      these are even stricter...
c      IF(EMAX.LT.0.015)  DT = 1.1*DTSAVE
c      IF(EMAX.LT.0.007)  DT = 1.2*DTSAVE
c      IF(EMAX.LT.0.001)  DT = 1.4*DTSAVE
c      IF(EMAX.LT.0.0008)  DT = 2.0*DTSAVE
c      IF(EMAX.LT.0.0004)  DT = 3.0*DTSAVE
c      IF(EMAX.LT.0.0001) DT = 4.0*DTSAVE
c      IF(EMAX.LT.0.00005) DT = 5.*DTSAVE


      DTINV = 1./DT
      ZMAX = Z(JS)
      if(ISOTOPE.EQ.0) then
      print 373, n, TIME, DT,EMAX,ISPEC(IS),ZMAX, USOL(is,js),
     $ USOL(LO2,1), USOL(LCH4,1), USOL(LH2,1), USOL(LCO,1)
      else
      print 373, n, TIME, DT,EMAX,ISPEC(IS),ZMAX, USOL(is,js),
     $ UINERT(LO2-Loff,1), UINERT(LCH4-Loff,1), UINERT(LH2-Loff,1), 
     $ UINERT(LCO-Loff,1)
      endif

 373  format (1x, I6, 1P3E12.3,2x,A8,1P1E12.3,4x,1P5E12.3)
C
      IF (SM-MS.GT.0.01) GOTO 317 ! skip oxidation state and sulfur budget

!following section run every three timesteps and printed to out.out
!oxidation state stuff commented out for now.

      write(14, 100) N,EMAX,ISPEC(IS),ZMAX,UMAX,RMAX,DT,TIME
 100  FORMAT(1X,'N =',I4,2X,'EMAX =',1PE9.2,' FOR ',A8,
     2  'AT Z =',E9.2,1X,'U =',E9.2,1X,'RHS =',E9.2,
     3  2X,'DT =',E9.2,2X,'TIME =',E9.2)
C
C   COMPUTE ATMOSPHERIC OXIDATION STATE
c   what follows needs work - 
      DO 42 I=1,NQ
      SR(I) = 0.
      DO 43 J=1,JTROP
  43  SR(I) = SR(I) + RAINGC(I,J)*USOL(I,J)*DEN(J)*DZ(J)
      PHIDEP(I) = VDEP(I)*USOL(I,1)*DEN(1)
  42  TLOSS(I) = SR(I) + PHIDEP(I)
c nb that vdep for particles is defined to include wfall when particles are in the main loop

      if(USETD.EQ.1) then
!this mimics the code Jim has, but is more general.
!I don't think I ever got this working 100% correctly to where I could balance the sulfur budget when using the tri-diag
       LP=1
      do i=NQ+1,NQ1 !need to fill up rainout and depostion vectors for the triadiagonal species to make the budgets work out. 
       SR(I)=0.
       if (ISPEC(I).EQ.'SXO4AER') then   !ISOHACK
       do j=1,JTROP
        SR(I) = SR(I)+ RAINGC(LH2SO4,J)*PARTICLES(J,LP)*DEN(J)*DZ(J)   !ACK - all particles raining out like H2SO4 !ISOHACK
       enddo
       endif

       PHIDEP(I)=(WFALL(1,LP)+vturb)* PARTICLES(1,LP)*DEN(1)  !ACK - hardcoded turbulent diffusion velocity !ISOHACK - will need to change in ISO
       TLOSS(I) = SR(I) + PHIDEP(I)
c       print *, ISPEC(I),SR(I),PHIDEP(I),SR(I)+PHIDEP(I)  !in general SR>>PHIDEP for particles,by about 100X
       LP=LP+1
      enddo   

c      stop
      endif  !end tri-diag budgeting loop

C
c the following is obsolete I hope
c      PHIESC = 2.5E13 * (0.5*USOL(LH,NZ) + USOL(LH2,NZ))
c    2  + USOL(LH2O,NZ) + 2.*USOL(LCH4,NZ) + 3.*USOL(LC2H6,NZ))
C


C - MC check the below with respect to the new boundary conditions.  Are some SGFLUXes that are not activly used in the code being counted here?

c - these could be done better in the manner of the redox computation
c - also how will these do if the L number doesn't exist, or if the sl/ll thing.
!tloss runs from 1-nq1 so cant have any SL's in here...

c      H2PROD = SGFLUX(LH2) + SGFLUX(LCO) + 4.*SGFLUX(LCH4) + 3.*
c     2  SGFLUX(LH2S) + 2.*TLOSS(LO2) + 0.5*TLOSS(LOH) + 1.5*TLOSS(LHO2)
c     3  + TLOSS(LO) + TLOSS(LNO) + TLOSS(LSO3)
c     3  + TLOSS(LH2O2) + 0.5*TLOSS(LHNO) + TLOSS(LH2SO4)
c     4  + TLOSS(LSO4AER) + 2.5*TLOSS(LHNO3) + 3.*TLOSS(LO3)



c      H2LOSS = PHIESC + 2.*TLOSS(LH2CO) + 3.*TLOSS(LH2S) 
c     2  + 2.5*TLOSS(LHS) + 1.5*TLOSS(LHSO) + 16.*TLOSS(LS8AER)
c     3  + 16.*TLOSS(LS8) + 8.*TLOSS(LS4) + 3.*TLOSS(LOCS)
c     4  + 6.*TLOSS(LS3) + 4.*TLOSS(LS2) + 1.5*TLOSS(LHCO)
c     5  + 2.*SGFLUX(LO2) + TLOSS(LSO)
c     6  + 7.*TLOSS(LC2H6) + 3.5*TLOSS(LCH3)
c  this looks ok provided that these terms are propoerly computed
C

c - note that these aren't really used for anything.  I am tempted to kill to allow for greater ease in switching SL and LL.  If I want to have a redox printout at every time step, I should find a way to make these generic like I did with redox in Output.


c-mc  the following is obsolete I hope
C   COMPUTE SULFUR BUDGET AND READJUST SO2 (H2S) OUTGASSING RATE IF SO
C   DESIRED (PROGRAM IS SET UP FOR PURE SO2 OUTGASSING)

c      SLOSS = TLOSS(LH2S) + TLOSS(LHS) + TLOSS(LS) + TLOSS(LSO) +
c     2  TLOSS(LSO2) + TLOSS(LH2SO4) + TLOSS(LHSO) + 2.*TLOSS(LS2) +
c     3  TLOSS(LSO4AER) + 4.*TLOSS(LS4) + 8.*(TLOSS(LS8) +
c     4  TLOSS(LS8AER)) + TLOSS(LOCS) + TLOSS(LSO3) + 3.*TLOSS(LS3)

c      SLOSSP = SLOSS - TLOSS(LSO2)
C
c      SFLUX = SGFLUX(LH2S) + SGFLUX(LSO2)    !check me at some point...

c-mc these sloss/slossP/sgflux printouts are meaningless now. 
c-mc should think about what (if anything) would be actually useful here to print out on every timestep


c      write(14,101)H2PROD,H2LOSS,SGFLUX(LSO2),SLOSS,SLOSSP,TLOSS(LS8AER)
c 101  FORMAT(10X,'H2PROD =',1PE10.3,2X,'H2LOSS =',E10.3,2X,'SO2FLX =',
c     2  E10.3,2X,'SLOSS =',E10.3,2X,'SLOSSP =',E10.3,2X,'S8LOSS =',
c     3  E10.3/)

 317  continue
C
C   RETRY TIME STEP IF EMAX EXCEEDS 25 PERCENT
      IF (EMAX.GT.0.25) THEN 
        DT = 0.7*DTSAVE
        TIME = TIME - DTSAVE
        do I=1,NQ
         do J=1,NZ
          USOL(I,J) = USAVEOLD(I,J)
         enddo
        enddo
      ELSE  !valid timestep, so update USOLPREV vector
       do i=1,nq
        do j=1,nz
         USOLPREV(I,J)=USAVEOLD(I,J)
        enddo
       enddo 
      ENDIF

c-should do some testing here to verify this is working as intended.
c-seems to be working
c$$$      IF (N .EQ. 50) THEN
c$$$         print *, 'old'
c$$$         print *, (USOLPREV(K,1),K=1,NQ)
c$$$         print *, ''
c$$$         print *, 'new'
c$$$         print *, (USOL(K,1),K=1,NQ)
c$$$         print *,''
c$$$         print *, 'diff'
c$$$         print *, (USOLPREV(K,1)-USOL(K,1),K=1,NQ)
c$$$         stop
c$$$      ENDIF



C
      NS = N/NPR    !NPR=PRN set above to 50 - i.e. write out every 50 steps
      SN = N/PRN
      IF(NN.EQ.NSTEPS) SN = NS
      IF (SN-NS .LT. 1.E-3) THEN 
        CALL OUTPUT(NN,NSTEPS,TIME,jtrop, vdep,USOLORIG,USETD)
      ENDIF 

c ISOHACK
c      write(23, 374) n, TIME, USOL(LO2,1), USOL(LH2,1), USOL(LCO,1),
c     $ USOL(LCH4,1)!, USOL(LS8aer,1)   !tri-diag

 374  format (1x, I6, 1P6E13.4)



c-mc writing out full number densities at each timestep
!this eventually should be an option as it is a large output file (1MB per 50 steps)

      
      if (tfac.eq.0.5) then 
         NSKIP=2
      else if (tfac.eq.0.25) then 
         NSKIP=4
      else if (tfac.eq.0.125) then
         NSKIP=8
      else if (tfac.eq.0.0625) then
         NSKIP=16
      else 
         NSKIP=1
      endif   
      SKIPN=NSKIP
      
      NS=N/NSKIP
      SN=N/SKIPN
      if (NN.EQ.NSTEPS) SN=NS
      if (SN-NS .lt. 1.e-3) then
c       write(43,115) N,TIME,DT,EMAX   !ISOHACK
        do I=1,nz
c not needed in whiff testing         write(43,114), (USOL(K,I)*DEN(I),K=1,NQ)   
c         write(43,114), (USOL(K,I)*DEN(I),K=1,NQ)   
        enddo
      endif

 114  format(100(1pe10.3))   !ACK hardcoded NQ - update if NQ>100
 115  format(I5, 3(1pe14.6)) 


      IF (INFO.NE.0) STOP
      IF (NN.EQ.NSTEPS) then
         finaln=NN
         GO TO 22
      endif

      IF (TIME.GT.TSTOP) then
        finaln=NN+1
        NN = NSTEPS - 1
      endif

   1  CONTINUE
C ***** END THE TIME-STEPPING LOOP *****
C
  22  CONTINUE   ! successful completion


C write out formatted out.dist file
c this new format works automatically even if NQ changes

      IROW = 10  !num columns in .dist file
      LR = NQ/IROW + 1      
      RL = FLOAT(NQ)/IROW + 1  
      DIF = RL - LR  
      IF (DIF.LT.0.001) LR = LR - 1 
C

      DO L=1,LR
       K1 = 1 + (L-1)*IROW
       K2 = K1 + IROW - 1
       IF (L.lt.LR) then
        write(18, 880) ((USOL(k,i),K=K1,K2),i=1,nz) 
       ELSE
          K2 = NQ
          if (K2-K1 .EQ. 9) then   !this only occurs if NQ is a multiple of 10
            write(18, 880) ((USOL(k,i),K=K1,K2),i=1,nz) 
          else  !if not, need to generate a dynamic format statement
           fmtstr='(  E17.8)'
           write(fmtstr(2:3),'(I1)')K2-K1+1   !OK for one character as this should always be <10
           write(18, fmtstr) ((USOL(k,i),K=K1,K2),i=1,nz)
          endif
       ENDIF
      enddo 

 880  format(10E17.8)
 881  format(5E17.8)


        write (18,881) (T(i),EDD(i),DEN(i),O3(i), SL(NSP-1,i),i=1,nz) !this final one is CO2 number density


        fmtstr='(  E17.8)'
        write(fmtstr(2:3),'(I2)')NP*3   
        do i=1,nz
         write(18,fmtstr) (AERSOL(i,j),j=1,NP),(WFALL(i,j),j=1,NP),
     $                  (RPAR(i,j),j=1,NP) 
        enddo  

        fmtstr='(  E17.8)'
        write(fmtstr(2:3),'(I2)')NP   
        
      if(USETD.EQ.1) then
      do i=1,nz
       write(18,fmtstr)  (PARTICLES(i,j),j=1,np)  !ordering is that in species.dat
      enddo
      endif

c new abstracted photorates printout

      IROW = 8  !num columns of photorates
      LR = KJ/IROW + 1      
      RL = FLOAT(KJ)/IROW + 1 
      DIF = RL - LR  
      IF (DIF.LT.0.001) LR = LR - 1 

      write(14,885)
      DO L=1,LR
       K1 = 1 + (L-1)*IROW
       K2 = K1 + IROW - 1
       IF (L.lt.LR) then
        write(14, 884) (photolabel(K),K=K1,K2) 
        write(14, 883) (Z(I),(prates(k,i),K=K1,K2),i=1,nz,3) 
       ELSE
          K2 = kj
           fmtstr="(/5X,'Z',6X,  A11)"
           write(fmtstr(13:14),'(I1)')K2-K1+1   
           write(14, fmtstr) (photolabel(K),K=K1,K2) 

           fmtstr2='(  (1PE9.2,2X))'
           write(fmtstr2(2:3),'(I2)')K2-K1+2   !+2 here so it writes the Z as well   
           write(14, fmtstr2) (Z(I),(prates(k,i),K=K1,K2),i=1,nz,3) 
       ENDIF
      enddo 

       fmtstr="(  A12)"
       write(fmtstr(2:3),'(I2)')kj    !format for photorate labels, with extra space !will crash if kj>99


c-mc  write out important parameters:
      write(49,*) NZ,NQ,NQ1,NSP2,NR,KJ,NP
      write(49,fmtstr) photolabel
      write(49, *) ISPEC

c-mc write out photorates at all heights in the .so2 file
      if(ISOTOPE.EQ.0) then
       JSO2=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'SO2    ')
       JO2=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'O2     ')+1
! note this is for the O2 + Hv -> O + O reaction,which is the second O2 reaction
! JH2O and JCO2 defined above...
      write(27,299) (Z(I),prates(JSO2,I),prates(JO2,I),prates(JH2O,I),
     $     prates(JCO2,I),I=1,NZ)
      else
       JSO2=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'SXO2    ')
      write(27,299) (Z(I),prates(JSO2,I),prates(JSO2,I),prates(JSO2,I),
     $     prates(JSO2,I),I=1,NZ)
!hack to include same number of lines in the ISOout.so2 - I will just have to remember not to do anything
!with the O2,H2O,and CO2 "photorates" when looking at this file
      endif
      


 299  FORMAT(5(1PE9.2,2x))
c-mc

c-mc write out all photorates at all heights in the .rates file
      string='(  (1PE9.2,2X))'
      write(string(2:3),'(I2)')KJ+1   !ack - hardcoded kj+1 here. will break if kj>99

!the above is a way to dynamically create a format string at runtime.replaces:  599  FORMAT(55(1PE9.2,2x)) 
       write(28,string) (Z(I),(prates(J,I),J=1,KJ),I=1,NZ)
c-mc

        STOP 


  20  CONTINUE   ! error in reactions
        PRINT 300,I
 300    FORMAT(//1X,'NMAX EXCEEDED FOR SPECIES ',I3)
        STOP 
  25  CONTINUE   ! error in reactions
        PRINT 301,IERR
 301    FORMAT(//1X,'ERROR IN REACTION ',I3)
        STOP
      END

