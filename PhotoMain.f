            PROGRAM PhotoMain
c
c- at some point go through and clean up all comments
c
c
c - this code contains the variable grid size changes used to compute
c - the suite of models for the whiff paper.
c - this code could/should be modified to use a variable grid size at some point
c - to make it faster all common blocks abstracted to DATA/INCLUDE


c - g77 will no longer work as there are some F90 built-ins being used
c see makefile for compilation syntax

cc

c-mc
c this code has iterated jacobian and batch multipliers
c-072406 - has 195K CO2 cross-section ! turned off for photo development

c-mc 080206 - iterated jacobion, new co2 x-section are turned off,TSTOP is 1e16
c-mc 112107 - iterated jacobians back in play for reruns of kevin's models for
c-            the whiff paper
c-mc 111308 - iterated jacbians off for chlorine testing...
c-mc 010609 - iterated jacobians back on during chlorine testing...

c-mc  THIS VERSION OF THE CODE WILL BE TIME-INDEPENDENT
c-    (i.e. P,T,zenith angle remain constant) (or maybe not...)

c-mc things to research at some point
c CO photolysis
c - whats up with HO2 photolysis?
c- where is HNO2?

c   the code reads in 4 files:
c     aerosol.table
c     photo.dat
c     sulfur.rx
c     in.dist

c  fixed typo in molecular transport terms 8-13-05

c  this version uses N, O3, and HNO3 as long-lived
c  it does not use CO2 as long-lived.  there problems be.

c  I need to do some CO expts for MC and DC.
c  but for giggles I am doing an experiment wiht big O2 and CH4 fluxes
c  to see what happens
c  this is hard to do if one uses O2 flux; it works better by setting
c  the O2 mixing ratio

c   some observations
c      1. problems with NOx at the upper boundary often means a problem with
c         H escape; NOx problems often resolve themselves with patience...
c      2. getting lost often means a problem with H escape
c      3. problems with S also means that there is a problem with redox balance.


c  update 6-23-05 inserts molecular diffusion of H
c  looks good

c   I added molecular diffusion for H2 and H.
c
c   the more general equation would include the molecular diffusive flux
c         \phi = b*f*({1/over H_A}-{1/over H}) - b*df/dz
c
c    [ w_1 - w_2 = {b_{12}\over n} \left( (1/n_1)(dn_1/dz) - (1/n_2)(dn_2/dz)
c                  - (1/f)(df/dz) \right)  ]
c  (I've not been careful. here f=n_1/n_2 <<1, so n=n_2.
c     I should look this up for the general case)
c
c     where H would be the scale height in vacuum of the species corresponding
c     to f and H_A is the scale height of the atmosphere.

c         \phi = b*f*({1/over H_A}-{1/over H}) - (b+KN)*df/dz
c          N df/dt = P - LNf - d\phi/dz - f*dN/dt
c  I have not done the work to allow for N to change with time - so let dN/dt=0
c
c    From this point hence one must be careful
c
C         df/dt  =  P/N - L*f + (1/N)*d/dz[(K*N+b)*(df/dz)]
c                   + (b/N)*d/dz({1/over H}-{1/over H_A})f
c
c    the finite difference form is
c
C     df(j)/dt  =  P(j)/N(j) + (DU(j) +ADU(j))*f(j+1)
c                 + (DL(j) +ADL(j))*f(j-1) - (DD(j) + ADD(j) + L(j)) *f(j)
c
c     where
c        DU(j) = [KN(j+) + b(j+)]/[N(j)*dz*dz]
c       ADU(j) = +b(j+)/[2N(j)*dz] *({1/over H} - {1/over H_A})
c        DL(j) = [KN(j-) + b(j-)]/[N(j)*dz*dz]
c       ADL(j) = -b(j-)/[2N(j)*dz] *({1/over H} - {1/over H_A})
c        DD(j) = DU(j) + DL(j)
c       ADD(j) = -ADU(j) - ADL(j)

c  at the lower boundary
c
C     df(0)/dt  =  P(0)/N(0) + DU(0)*(f(1)-f(0))  + ADU(0)*(f(1)+f(0))
c                 + FLUX(0)/N(0)/dz
c     where
c        DU(0) = [KN(0.5) + b(0.5)]/[N(0)*dz*dz]
c       ADU(0) = +b(0.5)/[2N(0)*dz] *({1/over H} - {1/over H_A})

c     in constant flux
c        djac(0) = chemistry(0) + 1/dt  + DU(0)  - ADU(0)
c        djac(1) = - DU(0)  - ADU(0)   (off diagonal)
c
c     in constant vdep
c        djac(0) = chemistry(0) + 1/dt  + DU(0)  - ADU(0) + vdep/dz
c        djac(1) = - DU(0)  - ADU(0)   (off diagonal)
c
c     in constant mixing ratio - i've left the chemistry out...
c        djac(0) = 1/dt  + DU(0)  - ADU(0)
c        djac(1) =  0.0   (off diagonal - is this correct??
c       - it may not be, but it should not matter given that it can't change f)
c
c
c     at the upper boundary:
C     df(nz)/dt  =  P(nz)/N(nz) + DL(nz)*(f(nz)-f(nz-1))
c                   + ADL(nz)*(f(nz)+f(nz-1)) - FLUX(nz)/N(nz)/dz
c     where
c        DL(nz) = [KN(nz-0.5) + b(nz-0.5)]/[N(nz)*dz*dz]
c       ADL(nz) = +b(nz-0.5)/[2N(nz)*dz] *({1/over H} - {1/over H_A})
c
c     in constant flux
c        djac(nz) = chemistry(nz) + 1/dt  + DL(nz)  - ADL(nz)
c        djac(nz-1) = - DL(nz)  - ADL(nz)   (off diagonal)
c
c     in constant veff
c        djac(nz) = chemistry(nz) + 1/dt  + DL(nz)  - ADL(nz) + veff/dz
c        djac(nz-1) = - DL(nz)  - ADL(nz)   (off diagonal)
c
c     in zero mixing ratio gradient   (df/dz=0)  (I've never used this...)
c       this turns out to be exactly the same as constant veff with
c        veff == +b(nz+0.5)/[ N(nz)*dz] *({1/over H} - {1/over H_A})
c
c

c  - what Jim did was to set the gradient of the H, H2 mixing ratios
c   to zero at the top (no flux upper boundary).
c   This does a better job with the chemistry.
c   i've left this intact for ch4co2.for
c
c
c- 4-25-05  jfk suggests distributing SO2 source over the troposphere
c           this way i can use a surface BC
c    This was implemented successfully
c
c  4-27-05  FIXED the accounting errors with fluxes.
c  the tridiagonal was not being included in the accounting
c  I wasn't sure how to do this save by putting S8aer into long lived species
c  This has now been done.  Fluxes are now under control.
c  redox balance and S conservation are both enforced. Fundamentally
c  nothing has really changed, but it has the advantage of being done right.

c  checked reactions 4/11/05
c  HSO is problematic - very few known reactions.
c  It does not react with O2

c  I have also likely over-estimated abstraction by S

c  CH3OH, C2H4, and C2H2 are missed in the chemical scheme
c   and would need to return if CH4 is made big

c to use this for Mars I need to make SO3 a long-lived species
c  while making S, S2, S4, and S8 short-lived if possible
c  program does not work well.
c  it gets lost on CO, which wanders all over the place...
c  jim suggests that I need to do H2 better
c   I think he is right.
c  I need to add OCS and treat CO2 as free
c  I have a lot of OCS reactions...
c  for giggles I added fictional reaction
c     SO3     CO   -> CO2     SO2             1.00E-19
c
c  to add OCS, I'd want these (taken from hncos.f)
c  I'd also want to find out if S + CO can go...
c     H       OCS     CO      HS         9.10E-12  0.00   1940.0   0.0      0.0
c     HS      CO      OCS     H          5.95E-14  1.12   8330.0   0.0      0.0
c     O       OCS     CO      SO         7.83E-11  0.00   2620.0   0.0      0.0
c     O       OCS     S       CO2        8.33E-11  0.0   5.53E+03  0.0      0.0
c     OCS     S       CO      S2         1.00E-10  0.0   4.56E+03  0.0      0.0
c     OCS     OH      CO2     HS         1.10E-13  0.0   1.20E+03  0.0      0.0
c     S       HCO     OCS     H          6.00E-11  0.0     0.0     0.0      0.0
c     OCS     HV      CO      S
c and...
c     S        CO     OCS                1.70E-33  0.0   1.51E+03  1.0      0.0
c  I'm guessing that S + CO goes at the same rate as O + CO

c  the current version has the fictional
c     SO3     CO      CO2     SO2        1.00E-16  0.0     0.0     0.0      0.0


C         THIS PROGRAM DIFFERS FROM PRIMS2 BY EXPLICITLY INCLUDING VAPOR
C     PHASE S8.  IT WAS USED TO GENERATE THE RESULTS FOR THE UV SCREEN
C     PAPER AT 2 BARS OF CO2.

c-mc the above needs to be exorcised...

C         THIS PROGRAM IS DESIGNED FOR EXTREMELY LOW (PRE-PHOTOSYNTHETIC
C     O2 LEVELS.  THE PHOTOLYSIS OF O2 IN THE SCHUMANN-RUNGE BANDS CAN
C     BE CALCULATED EITHER BY THE METHOD OF ALLEN AND FREDERICK OR BY
C     EXPONENTIAL SUMS.  THE LATTER METHOD SHOULD BE USED FOR LOW-O2
C     ATMOSPHERES (I.E. USE IO2 = 1).  LIKEWISE, NO PHOTOLYSIS SHOULD
C     BE CALCULATED USING THE CIESLIK AND NICOLET METHOD (INO = 1),
C     WHICH HAS BEEN MODIFIED TO AGREE WITH THE BAND INTENSITIES OF
C     FREDERICK AND HUDSON (1979).
C          THIS VERSION OF THE PROGRAM HAS A NEW OPTION FOR INCLUDING
C     TRANSPORT OF SPECIES THAT ONE DOES NOT WISH TO INCLUDE IN THE BIG
C     REVERSE EULER MATRIX.  THEY CAN BE SOLVED SEPARATELY USING A TRI-
C     DIAGONAL INVERSION METHOD.  S8 PARTICLES ARE TREATED THIS WAY IN
C     THIS PROGRAM BECAUSE THEY ARE VIRTUALLY NON-EXISTENT UP HIGH.
C     THE TOP OF THE TRIDIAGONAL MATRIX IS AT GRID POINT MZ, WHICH CAN
C     BE SET EQUAL TO OR LESS THAN NZ.
C
C       THIS PROGRAM IS A ONE-DIMENSIONAL MODEL OF THE PRIMORDIAL
C     ATMOSPHERE.  THE MIXING RATIOS OF THE LONG-LIVED SPECIES
C     ARE CALCULATED FROM THE EQUATION
C
C         DF/DT  =  (1/N)*D/DZ(KN*DF/DZ + WNF) + P/N - LF
C
C     WHERE
C     F = MIXING RATIO (USOL)
C     K = EDDY DIFFUSION COEFFICIENT (EDD)
C     N = TOTAL NUMBER DENSITY (DEN)
C     L = CHEMICAL LOSS FREQUENCY (XL)
C     P = CHEMICAL PRODUCTION RATE (XP)
C     W = FALL VELOCITY (WFALL) FOR PARTICLES, POSITIVE DOWNWARD (ALL
C         FLUXES, BY CONTRAST, ARE POSITIVE UPWARD)
C
C          TWO TYPES OF PARTICLES ARE INCLUDED: SULFATE (SO4AER) AND
C     ELEMENTAL SULFUR (S8AER).  SULFATE IS INCLUDED IN THE BIG MATRIX,
C     WHEREAS S8 IS TREATED SEPARATELY.  THE PARTICLES ARE ASSUMED TO
C     BE 0.1 UM IN RADIUS UP HIGH; THEY GET BIGGER AT LOW ALTITUDES BE-
C     CAUSE OF COAGULATION. THIS IS DONE CRUDELY, BASED ON A COMPARISON
C     OF THE RELATIVE LIFETIMES AGAINST DIFFUSION, SEDIMENTATION, AND
C     RAINOUT.
C
C          THE SYSTEM OF PDES IS SOLVED USING THE REVERSE EULER
C     METHOD.  LETTING THE TIME STEP GO TO INFINITY GIVES YOU NEWTONS
C     METHOD, E.G. IT REVERTS TO AN EFFICIENT STEADY-STATE SOLVER.
C
C          THE LIST OF SUBROUTINES IS AS FOLLOWS:
C     (1) GRID   -  SETS UP THE ALTITUDE GRID
C     (2) RATES  -  DEFINES CHEMICAL REACTION RATES AND RAINOUT RATE
C     (3)   RAINOUT - COMPUTES RAINOUT RATES USING THE METHOD OF GIORGI
C                     AND CHAMEIDES (1985)
C     (3.1) AQUEOUS - DOES RAINWATER CHEMISTRY (CALLED BY RAINOUT)

c-mc  Initphoto sets wl grid and cross-sections
c-mc Initmie read scattering params

c-mc should fill this out at some point

C     (3.2) PHOTO   - COMPUTES PHOTOLYSIS RATES (CALLS MSCAT)
C     (3.3) O3PHOT  - COMPUTES COEFFICIENTS USED TO FIND O(1D) QUANTUM
C                     YIELDS IN O3 PHOTOLYSIS
c-mc the above is no longer used...
C     (4)   DENSTY  - COMPUTES ATMOSPHERIC NUMBER DENSITIES FROM HYDRO-
C                     STATIC EQUILIBRIUM
C     (5)   DIFCO   - COMPUTES DK = K*N BETWEEN GRID POINTS; ALSO FINDS
C                     DIFFUSION LIFETIMES H*H/K
C     (5.1) SEDMNT  - CALCULATES FALL VELOCITIES AND ESTIMATES PARTICLE
C                     SIZES
C     (5.2) SATRAT  - COMPUTES SATURATION H2O MIXING RATIOS AT ALL HEIGHTS
C                     FINDS MANABE/WETHERALD RH DISTRIBUTION IN TROPOSPHERE
C     (5.3) AERTAB  - READS PAT HAMILL'S H2SO4 TABLE AND CALCULATES VAPOR
C                     PRESSURES OF H2O AND H2SO4 AT EACH HEIGHT AS A FUNC-
C                     TION OF SULFATE CONTENT OF THE PARTICLES
C     (5.4) AERCON  - FINDS WEIGHT PERCENT OF H2SO4 IN SULFATE PARTICLES
C                     ALONG WITH H2SO4 VAPOR PRESSURE, GIVEN T AND H2O AT
C                     EACH HEIGHT
C     (6) OUTPUT -  PRINTS OUT RESULTS
C     (7) DOCHEM - DOES CHEMISTRY FOR ALL SPECIES AT ALL GRID
C                  POINTS BY CALLING CHEMPL
C     (8) CHEMPL - COMPUTES CHEMICAL PRODUCTION AND LOSS RATES
C                  FOR ONE SPECIES AT ALL GRID POINTS
C     (9) LTNING -  COMPUTES LIGHTNING PRODUCTION RATES FOR O2 AND
C                   N2 BASED ON CHAMEIDES' RESULTS
C    (11) MSCAT  -  DOES RAYLEIGH SCATTERING USING YUK YUNG'S TECHNIQUE
C
C          OTHER DEFINED FUNCTIONS INCLUDE:
C     (1) TBDY   -  COMPUTES 3-BODY REACTION RATES
C     (2) E1     - EXPONENTIAL INTEGRAL OF ORDER ONE

c - some of Kevin's notes on possible reactions
c-keeping around for the time when more research is done on the chemical species
C    9.6e-12. HSO2 is not a species in this system - it reacts quickly with O2
C     ) HSO + NO2 -> NO + HSO2
C    this has measured but very slow reaction
c     )  H + SO2 + M -> HSO2 + M
C    3e-13 measured HO2 as product.  not a clear situation;
c         - my temptation is to assign OH, H rates measured at high T
c     ) HSO2 + O2 -> HO2 + SO2
c        RATE: 3e-12
c     ) HSO2 + H -> H2 + SO2
c        RATE: 8e-12
c     ) HSO2 + OH -> H2O + SO2
C        JPL has this!  but what do you do with it??
c     ) HS + NO + M -> HSNO + M
C        NIST 2e-12*exp(-1000/T)
c     ) HCO + HNO  -> H2CO + NO
C
C        NIST - v. fast: 2.5e-27*(298/T)^-3.8, 1.5e-10 limit
C     )  CH3 + OH + M ->  CH3OH + M
C     )  CH3 + OH  ->  H2O + CH2    ! NIST 1.2e-10 *exp(-1400/T)
C     )  CH3 + OH  ->  HCOH + H2    ! NIST 9e-11 *exp(-1500/T)
C     )  CH3 + OH  ->  H2COH + H    ! NIST 1.3e-11
c     )  S2  + H  + M ->  HS2 + M
c     )  S2  + OH  ->               ! SO +OH -> H + SO2 8.6e-11
c     )  S3  + OH  ->
c     )  S3  + O  ->  SO + S2
c     )  S3  + H  + M ->  HS3 + M
c     )  S4  + OH  ->
c     )  S4  + O  -> SO + S3
c     )  S4  + H  + M -> HS4 + M
c     )  HS2 + OH -> H2O + S2
c     )  HS2 + O  -> OH + S2
c     )  HS2 + H  -> H2 + S2
c     )  HS2 + S  -> H + S3
c  etc.
C
C        THIS PROGRAM DOES THE CHEMISTRY AUTOMATICALLY.
C     THE CHEMICAL
C     REACTIONS ARE ENTERED ON DATA CARDS IN FIVE 10-DIGIT COLUMNS
C     STARTING IN COLUMN 11, I.E.
C
C         REAC1     REAC2     PROD1     PROD2     PROD3
C
C     THE IMPORTANT PARAMETERS DESCRIBING THE CHEMISTRY ARE
C        NR   = NUMBER OF REACTIONS
C        NSP  = NUMBER OF CHEMICAL SPECIES
C        NSP1 = NSP + 1 (INCLUDES HV)
C        NSP2 = NSP + 2 (INCLUDES M)
C        NQ   = NUMBER OF SPECIES FOR WHICH A DIFFUSION EQUATION
C                 IS SOLVED AND WHICH ARE IN THE BIG, REVERSE EULER MATRIX
C        NQ1 = TOTAL NUMBER OF SPECIES FOR WHICH TRANSPORT IS INCLUDED
C              (THOSE NOT IN THE BIG MATRIX ARE SOLVED WITH A STEADY-STATE
C              TRIDIAGONAL INVERSION METHOD)
C        NMAX = MAXIMUM NUMBER OF REACTIONS IN WHICH AN INDIVIDUAL
C               SPECIES PARTICIPATES
C
C        PHOTOLYSIS REACTIONS ARE IDENTIFIED BY THE SYMBOL HV (NOT
C     COUNTED IN EVALUATING NSP).  THREE-BODY REACTIONS ARE WRITTEN
C     IN TWO-BODY FORM, SO THE DENSITY FACTOR MUST BE INCLUDED IN
C     THE RATE CONSTANT.
C        THE CHEMICAL REACTION SCHEME IS STORED IN THE FOLLOWING MATRICE
C
C     ISPEC(NSP2) = VECTOR CONTAINING THE character NAMES OF THE
C                  CHEMICAL SPECIES.  THE LAST TWO ENTRIES ARE HV AND M.
C     CHEMJ(5,NR) = MATRIX OF CHEMICAL REACTIONS (chars). THE FIRST TWO ARE
C                   REACTANTS, THE LAST THREE ARE PRODUCTS.
C     JCHEM(5,NR) = MATRIX OF CHEMICAL REACTIONS (indexed).  THE FIRST TWO ARE
C                   REACTANTS, THE LAST THREE ARE PRODUCTS.
C     ILOSS(2,NSP,NMAX) = MATRIX OF LOSS PROCESSES.  ILOSS(1,I,L)
C                         HOLDS REACTION NUMBER J, ILOSS(2,I,L) HOLDS
C                         REACTANT NUMBER.
C     IPROD(NSP,NMAX) = MATRIX OF PRODUCTION PROCESSES.  IPROD(I,L)
C                       HOLDS REACTION NUMBER J.
C     NUML(NSP) = NUMBER OF NON-ZERO ELEMENTS FOR EACH ROW OF ILOSS
C     NUMP(NSP) = NUMBER OF NON-ZERO ELEMENTS FOR EACH ROW OF IPROD
C
      INCLUDE 'PHOTOCHEM/INPUTFILES/parameters.inc'
      implicit real*8(A-H,O-Z)
      character*10 buffer
      character*8 pstar
      CHARACTER*8 ISPEC, CHEMJ,SPECIES,SPECTYPE,REACTYPE,PLANET
      CHARACTER*20 string,fmtstr,fmtstr2
      real*8 mass
      DIMENSION FVAL(NQ,NZ),FV(NQ,NZ),DJAC(LDA,NEQ),RHS(NEQ),IPVT(NEQ)
      dimension SMFLUX(NQ),SGFLUX(NQ),VDEP(NQ),VEFF(NQ)
      dimension   DD(NQ1,NZ),DL(NQ1,NZ),DU(NQ1,NZ)
      dimension ADL(NQ,NZ), ADU(NQ,NZ), ADD(NQ,NZ)
      dimension USAVE(NQ,NZ),R(NZ),U(NQ),VDEP0(NQ),VEFF0(NQ)
      dimension USAVEOLD(NQ,NZ), USOLPREV(NQ,NZ)
      dimension USOLORIG(NQ,NZ)
      dimension alt_new(NZ), T_new(NZ), water(NZ)
      dimension alt_dontuse(NZ), T_dontuse(NZ), water_fix(NZ)
Cc-mc USOLSAVE added for convergence testing 4/29/06,
Cc-mc removed from code on 2/28/07
Cc-mc going to keep declaration here and reserve it in case I ever bring back
Cc-mc the convergence-testing code, which is r37:36 in the /td branch
Cc-mc USOLPREV(NQ,NZ) added for second order reverse Euler calculations
      DIMENSION DPU(NZ,NP),DPL(NZ,NP)
      DIMENSION TA(NZ),TB(NZ),TC(NZ),TY(NZ)
      dimension PRES_bar(NZ)

      integer rhOcount,rhHcount,rhCcount,rhScount,rhNcount,rhCLcount
c      dimension atomsO(NSP2),atomsH(NSP2),atomsC(NSP2)
c      dimension atomsN(NSP2),atomsCL(NSP2),atomsS(NSP2)

!temp
      dimension testvec(NR),fixedmr(nq)
      dimension distheight(nq)

      CHARACTER*11 photolabel, AA
      CHARACTER*8  AX



c for sgbco testing (to get at Jacobian condition number)
c      real*8 RCOND
c      real*8 WORK(neq)


      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PHOTABLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/BBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/CBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc'
C     !can go away when MSCAT does WARNING DO WE NEED THIS COMMENT
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/EBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/FBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/GBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/JBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/NBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/RBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/SBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/ZBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/LTBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/AERBLK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/SULBLK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/SATBLK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/comPRESS1.inc'


C   CONSTANTS FOR 1ST EXPONENTIAL INTEGRAL  !can go away when MSCAT does
C                       -WARNING do we need can go away comment?
      DATA AI/.99999193, -.24991055, .05519968, -.00976004,
     2  .00107857, -.57721566/
      DATA BI/8.5733287401, 18.0590169730, 8.6347608925,
     2  .2677737343/
      DATA CI/9.5733223454, 25.6329561486, 21.0996530827,
     2  3.9584969228/
      DATA NUML,NUMP/NSP*0,NSP*0/

c OPEN FILES


C      ,form='UNFORMATTED')
      open(2, file='PHOTOCHEM/DATA/aerosol.table',status='OLD')
      open(3, file='PHOTOCHEM/DATA/photo.dat',status='OLD')
      open(4, file='PHOTOCHEM/INPUTFILES/species.dat',status='OLD')
      !planet parameters (G, FSCALE, ALB,ZTROP,etc)
      open(7, file='PHOTOCHEM/INPUTFILES/PLANET.dat',status='OLD')
C      Model Parameters (AGL, IO2,INO, LGRID, etc)
      open(231, file='PHOTOCHEM/INPUTFILES/input_photchem.dat',
     &          status='OLD')
C      REACTION FILE
      open(9, file='PHOTOCHEM/INPUTFILES/reactions.rx',status='OLD')
C      CODE OUTPUT
      open(14, file='PHOTOCHEM/OUTPUT/out.out',status='UNKNOWN')
C      FORMATTED INPUT
      open(17, file='PHOTOCHEM/in.dist',status='OLD')
C      FORMATTED OUTPUT
      open(18, file='PHOTOCHEM/OUTPUT/out.dist',status='UNKNOWN')
      open(34, file='PHOTOCHEM/PTZ_mixingratios_in.dist',
     &          status='UNKNOWN')
      open(35, file='PHOTOCHEM/OUTPUT/PTZ_mixingratios_out.dist',
     &          status='UNKNOWN')
C      TERSE OUTPUT
      open(19, file='PHOTOCHEM/OUTPUT/out.terse',status='UNKNOWN')
C      PARTIAL OUTPUT OF MIXING RATIOS
      open(67, file='PHOTOCHEM/OUTPUT/profile.pt',status='UNKNOWN')
C      HYDROCARBONS
      open(68, file='PHOTOCHEM/OUTPUT/hcaer.out',status='UNKNOWN')
C      OTHER HYDROCARBONS
      open(69, file='PHOTOCHEM/OUTPUT/hcaer2.out',status='UNKNOWN')
C      VERY TERSE OUTPUT
      open(21, file='PHOTOCHEM/OUTPUT/out.trs',status='UNKNOWN')
C      TIME OUTPUTS
      open(23, file='PHOTOCHEM/OUTPUT/out.time',status='UNKNOWN')
      open(24, file='PHOTOCHEM/OUTPUT/out.tim',status='UNKNOWN')
C-AVB,MAB  Final Output for P, T, Z and species fluxes (to be added)
      open(255, file='PHOTOCHEM/OUTPUT/PTZ_out.flux')
c-mc
C     redox output - eventually combine into out.trs
      open(25, file='PHOTOCHEM/OUTPUT/out.redox',status='UNKNOWN')
C     Temporary output file for looking at convergence
C         in the reverse Euler iteration
      open(26, file='PHOTOCHEM/OUTPUT/out.converge',status='UNKNOWN')
C     Printing out relevant SO2 photolysis pieces
C           between 190-220 as MIF signature
      open(27, file='PHOTOCHEM/OUTPUT/out.so2',status='UNKNOWN')
C     Reaction rates,reactions,species,densities,rate constants
      open(28, file='PHOTOCHEM/OUTPUT/out.rates',status='UNKNOWN')
C     FOLLOWING ARE CROSS SECTIONS, HEIGHT GRID, WAVELENGTH GRID:
      open(29, file='PHOTOCHEM/OUTPUT/out.xsec',status='UNKNOWN')
      open(30, file='PHOTOCHEM/OUTPUT/out.gridz',status='UNKNOWN')
      open(31, file='PHOTOCHEM/OUTPUT/out.gridw',status='UNKNOWN')

c-mc  unit 33 reserved for wavelength grids...

      open(41, file='PHOTOCHEM/OUTPUT/out.rad',status='UNKNOWN')
C     File below should contain TATLTUAE
      open(42, file='PHOTOCHEM/OUTPUT/out.finalden',status='UNKNOWN')
C     Number densities at each timestep
      open(43, file='PHOTOCHEM/OUTPUT/out.densities',status='UNKNOWN')
C     Total production and loss at steady state
      open(44, file='PHOTOCHEM/OUTPUT/out.prod',status='UNKNOWN')
C     Fluxes
      open(45, file='PHOTOCHEM/OUTPUT/out.flux',status='UNKNOWN')
C     tau=1 (at final step) in out.tau
      open(48, file='PHOTOCHEM/OUTPUT/out.tau',status='UNKNOWN')
C      Some model parameters
      open(49, file='PHOTOCHEM/OUTPUT/out.params',status='UNKNOWN')
C      NGE and L2 are normally between start and finish
      open(50, file='PHOTOCHEM/OUTPUT/out.error',status='UNKNOWN')
C      TP/FLOW for chlorine species,nitrate, adn sulfate
      open(51, file='PHOTOCHEM/OUTPUT/out.cl',status='UNKNOWN')

C     For testing O2 prates with various grid sizes
      open(58, file='PHOTOCHEM/OUTPUT/out.O2prates',status='UNKNOWN')
C        rainggc out - ISOHACK
      open(59, file='PHOTOCHEM/OUTPUT/out.raingc',status='UNKNOWN')

!c-mc 60 and 61 are opened below after LGRID is read in

C      lower boundary fluxes (not including rainout)
       open(62, file='PHOTOCHEM/OUTPUT/out.flow',status='UNKNOWN')
C      Haze optical depths
       open(63, file='PHOTOCHEM/OUTPUT/out.od',status='UNKNOWN')
C      Haze TOA and sur rpar
       open(73, file='PHOTOCHEM/OUTPUT/out.rp',status='UNKNOWN')


C - This file gives the input needed for SMART - radiative transfer code
c-mab: Commenting out smart inputs for now. (11/4/2016)
c       open(159, file='photochem_smart/photchem/photfile.pt',
c     &          status='UNKNOWN')
c       open(164, file='photochem_smart/atm/tag.atm',
c     &          status='UNKNOWN')

C - Seperating the out.dist file into seperate files

       open(70, file='PHOTOCHEM/OUTPUT/out.chem', status='UNKNOWN')
       open(66, file='PHOTOCHEM/OUTPUT/out.strctr', status='UNKNOWN')
       open(71, file='PHOTOCHEM/OUTPUT/out.aersol', status='UNKNOWN')
       open(72, file='PHOTOCHEM/OUTPUT/out.tridag', status='UNKNOWN')


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
      READ(231,*)AA,HCDENS
      IF(IDEBUG.eq.1) print *, "HCDENS =",HCDENS
      READ(231,*)AA, ICOUPLE
      IF(IDEBUG.eq.1) print *, "ICOUPLE =",ICOUPLE
      READ(231,*)AA, ihztype
      IF(IDEBUG.eq.1) print *, "IHZTYPE =",ihztype
      READ(231,*)AA, ZY
      IF(IDEBUG.eq.1) print *, "ZY =",ZY
 555  format(3/)
      close(231)

C     NO photolysis rates output
      if (LGRID.EQ.0) open(60, file='PHOTOCHEM/OUTPUT/out.NOprates',
     &                         status='UNKNOWN')
C    Wavelength specific SO2 photorates on HR grid
      if (LGRID.EQ.1) open(61, file='PHOTOCHEM/OUTPUT/out.so2HR',
     &                         status='UNKNOWN')



c The next four files are used when this model is coupled
C           with the climate model (ICOUPLE=1)
C   To be used as input for the climate model
C they are created and updated regardless of whether ICOUPLE=1 in input_photochem
       open(90, file='COUPLE/hcaer.photoout.out',status='UNKNOWN')
       open(84, file='COUPLE/fromPhoto2Clima.dat', status='UNKNOWN')
       open(116, file='COUPLE/fromClima2Photo.dat', status='UNKNOWN')
       open(117, file='COUPLE/mixing_ratios.dat', status='UNKNOWN')


C - READ IN SPECIES NAMES, ATOMIC NUMBERS, AND BOUNDARY CONDITIONS

C   counter for long lived species
      iLL=0
C   counter for short lived species
      iSL=0
C   counter for tridiagonal species
      iTD=0
C   counter for inert species
      iIN=0
C   counter for number of lines in species.dat file
      iSP=0
C   So the species.dat statements print only once
      iprint = 0

      do while (I.LT.300)
C     Note: Below will crash if species.dat is longer than 300 lines.
         read(4,*, end=96) SPECIES,SPECTYPE
         if (scan(species,'*').LT.1) then  ! else ignore comments in species.dat file (lines that start with *)
            iSP=iSP+1
            ISPEC(iSP)=species
         !   print *, iSP, species
C         This loads the "Lnumbers" for ease of use later in the code
            call LNUM(ISPEC(isP),iSP)
C             Return to previous line in species.dat file
              backspace 4
C  read in atmoic number data, NEVER use LC,LH,LN,LO,LS as placeholders
C  as they mean something else...
              read(4,*) AX,AX,LA,LB,LD,LE,LF,LM

            if (SPECTYPE.EQ.'LL') then
               iLL=iLL+1
C             Return to previous line in species.dat file
               backspace 4

C This section reads in the boundary conditions from species.dat.
C Note this now uses non-fixed formatting! :-D
                  read(4,*) AX,AX,LX,LX,LX,LX,LX,LX,LBC,XX,YY,ZZ,XXX,LG
     &   ,YYY,ZZZ

            !   print *, LBC
               LBOUND(iLL)=LBC
               VDEP0(iLL)=XX
               FIXEDMR(iLL)=YY
               if (LBOUND(iLL).eq.3) then
C             distributed flux
                distflux(iLL)=ZZ
               else
C             lower boundary flux
                  SGFLUX(iLL)=ZZ
               endif
               distheight(iLL)=XXX
               MBOUND(iLL)=LG
               SMFLUX(iLL)=YYY
               VEFF0(iLL)=ZZZ
C      CO2 only works as fixed mixing ratio. This could be handled better.
C-mab: What does above comment mean? It is NOT fixed/inert in my giant templates...
               if (species.EQ.'CO2') FCO2=YY

            endif

            if (SPECTYPE.EQ.'IN') then
               iIN=iIN+1
C              Returns to previous line in species.dat file
               backspace 4
C              Reads in fixed mixing ratios
               read(4,*) AX,AX,LX,LX,LX,LX,LX,LX,XX   !read in fixed mixing ratios
C            Hardcoding woohoo! need to do N2 as well WARNING
               if (species.EQ.'HE') FHE=XX
               if (species.EQ.'CO2') FCO2=XX
               if (species.EQ.'N2') FN2=XX
            endif

            if (SPECTYPE.EQ.'TD')iTD=iTD+1
            if (SPECTYPE.EQ.'SL') iSL=iSL+1
            if (SPECTYPE.EQ.'HV') iIN=iIN+1
            if (SPECTYPE.EQ.'M') iIN=iIN+1

               atomsO(iLL+iTD+iSL+iIN)=LA
               atomsH(iLL+iTD+iSL+iIN)=LB
               atomsC(iLL+iTD+iSL+iIN)=LD
               atomsS(iLL+iTD+iSL+iIN)=LE
               atomsN(iLL+iTD+iSL+iIN)=LF
               atomsCL(iLL+iTD+iSL+iIN)=LM

         endif
         I=I+1
      enddo

 96   CONTINUE

c      stop

C so far, the only use of mass in Difco,where it is summed over NQ,
C           -so it would be OK to go higher
C           Below are molecular weights (NQ1)
      mass=atomsO*16.+atomsH*1.+atomsC*12.+atomsS*32.+atomsN*14.
     $  + atomsCL*34.

C     !we are setting CLO as 'neutral"
C     !redoxstate goes from 1-NQ1 in Output
      redoxstate = atomsO*1.0 + atomsH*(-0.5) + atomsS*(-2.) +
     $  atomsCL*(-1.0) + atomsC*(-2)  !N=0

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



C Reading in the temperature and water profiles from the climate code

      IF (ICOUPLE.eq.1) then
        print *, 'alt_new, T_new, water'
        DO J=1, NZ

C Reading the altitude, temperature, and water profiles from the climate code
           READ(116,*) alt_new(J), T_new(J), water(J)
           print 351, alt_new(J), T_new(J), water(J)
        END DO
      close(116)
      endif
 351  FORMAT (1PE10.3, 1PE12.3, 1PE12.3)




C ***** READ THE CHEMISTRY DATA CARDS *****
corig      read (9,200) CHEMJ
 !     print *, 'this is nr', NR
      read (9,200) CHEMJ
corig 200  FORMAT(10X,A8,2X,A8,2X,A8,2X,A8,2X,A8)
 200  FORMAT(A8,2X,A8,2X,A8,2X,A8,2X,A8)
      write(14, 201) (J,(CHEMJ(M,J),M=1,5),J=1,NR)
 201  FORMAT(1X,I3,' ',5X,A8,' +  ',A8,'  =    ',A8,' +  ',A8,4X,A8)
      KJAC = LDA*NEQ
      write(14, 202) NQ,NZ,KJAC
 202  FORMAT(//1X,'NQ=',I2,5X,'NZ=',I3,5X,'KJAC=',I7)



c-mc this is bad code. there is probably some way to do this
c-mc    -with the read in above
c-mc or even if necessary, some way to redo the read without
c-mc closing and opening the file again     WARNING
      close(9)
      ! chemical reaction file
      open(9, file='PHOTOCHEM/INPUTFILES/reactions.rx',status='OLD')
      read(9,204) REACTYPE
 204  FORMAT(48X,A5)

C    close this because Rates.f and Initphoto.f will re-open it later.
      close(9)



C  ****  JCHEM has species numbers; CHEMJ is corresponding characters
C ***** REPLACE HOLLERITH LABELS WITH SPECIES NUMBERS IN JCHEM *****
      DO 5 J=1,NR
      DO 5 M=1,5
      IF(CHEMJ(M,J).EQ.' ') GO TO 5
      DO 6 I=1,NSP2
  !     print *, CHEMJ(M,J),ISPEC(I)
      IF(CHEMJ(M,J).NE.ISPEC(I)) GO TO 6
      JCHEM(M,J) = I
      GO TO 5
   6  CONTINUE
      IERR = J
      print *, ISPEC
      print *, 'ispec(i)', ISPEC(i)
      print *, (CHEMJ(L,J),L=1,5)
      ! quit; error in reactions
      GOTO 25
   5  CONTINUE
C

C ***** FILL UP CHEMICAL PRODUCTION AND LOSS MATRICES *****
      DO 7 M=1,2
C   !so N=2, then 1
      N = 3-M
      DO 7 J=1,NR
C   !so I = JCHEM(1,NR) then JCEHM(2,NR)
      I = JCHEM(M,J)
C   !skips 0 (i.e. nothing) and NSP1 (HV)
      IF(I.LT.1.OR.I.GT.NSP) GO TO 7
C   !counter for how many reactions species I is involved with
      NUML(I) = NUML(I) + 1
C   !quit; too many reactions  (seems unnecesary, but whatever)
      IF(NUML(I).GT.NMAX) GOTO 20
      K = NUML(I)
C   !ILOSS(1,species in pos 1, species reac#)
C     -then ILOSS(1,spec in pos 2, reac#)= global reaction #
      ILOSS(1,I,K) = J
C   !ILOSS(1,species in pos 1, species reac#)
C     -then ILOSS(1,spec in pos 2, reac#)= other species
      ILOSS(2,I,K) = JCHEM(N,J)
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

C     !assume 3 products unless..
      numprod=3
      if (JCHEM(5,i).EQ.0) numprod=2
      if (JCHEM(4,i).EQ.0) numprod=1

C   This loop counts up the mass on the right hand side of the .rx
      do j=0,numprod-1
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
c photospec(ks) - the unique photolysis reaction numbers
C                 (i.e. unique elements of photoreac)
c                 used in Initphoto.f to fill up sq, the cross section vector
c photonums(kj) - the reaction number of each photolysis reaction
c                 used in Photo.f to fill up the A vector of rates


C A vector of length nr with 1's in the location where photolysis reactions are
       testvec=INDEX(REACTYPE,'PHOTO')



        jcount=1
        jrcount=1
        juniq=1
       do i=1,nr
        if (testvec(i).eq.1.) then
C    Captures the species number of each photo reaction
           photoreac(jrcount)=JCHEM(1,i)
C    Captures the reaction number of each photoreaction
           photonums(jrcount)=i

c           print *, jrcount,i,JCHEM(1,i),(CHEMJ(m,i),m=1,5)

           jrcount=jrcount+1

C    Captures the unique photolysis species in photospec
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

c       print *, jnums
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


c-mc the below could be made into a nice printout in the out.out file
c     WARNING was this dealt with
c       print *, SUM(INDEX(REACTYPE,'PHOTO'))
c       print *, SUM(INDEX(REACTYPE,'2BODY'))
c       print *, SUM(INDEX(REACTYPE,'3BODY'))
c       print *, SUM(INDEX(REACTYPE,'2BACK'))
c       print *, SUM(INDEX(REACTYPE,'3BACK'))
c       print *, SUM(INDEX(REACTYPE,'WEIRD'))
c       print *, NR


C ***** READ THE PLANET PARAMETER DATAFILE *****
      READ(7,502) G,FSCALE,ALB,ZTROP,FAR,R0,P0,PLANET,TIMEGA,IRESET,
     &   pstar, ihzscale, uvscale
c-mab: Uncomment below for debugging with this part
C      	print*,'G,FSCALE,ALB,ZTROP,FAR,R0 = ',G,FSCALE,ALB,ZTROP,FAR,R0
C      	print*,'P0,PLANET,TIMEGA,IRESET = ',P0,PLANET,TIMEGA,IRESET
 502  FORMAT(F7.1/,F7.2/,F7.3/,E7.1/,F7.3/,E8.3/,F8.3/,A8/,F4.2/,I1/
     &  ,A8/,I1/,F8.3)
C     adding IRESET to create the atmospheric profile for modern earth
C     gna - added pstar keyword to change the star from planet.dat
C     print *, 'pstar is ', pstar
C     print *, 'ihzscale is ', ihzscale
C gna-scale stellar flux based on spectral type:
C uses Kopparapu et al 2012 scalings for earth-equivalent distance
      IF (ihzscale.eq.1) then
         IF (pstar.eq.'adleo') FSCALE = FSCALE * 0.870
         IF (pstar.eq.'t3200') FSCALE = FSCALE * 0.859
         IF (pstar.eq.'k2v') FSCALE = FSCALE * 0.950
         IF (pstar.eq.'f2v') FSCALE = FSCALE * 1.110
         IF (pstar.eq.'gj876') FSCALE = FSCALE * 0.866
         IF (pstar.eq.'proxima') FSCALE = FSCALE * 0.920
         IF (pstar.eq.'m8v') FSCALE = FSCALE * 0.911
         IF (pstar.eq.'wasp12') FSCALE = FSCALE * 1890.40
         print *, 'fscale is ', fscale
      ENDIF



C
      IF (IRESET.eq.1) CALL RESET
C Defines mixing ratio values for modern earth; don't need next section that
C       pulls ancient earth mixing ratios from in.dist

      IF (IRESET.eq.0) then
C ***** READ THE INPUT DATAFILE *****

c read in formatted input data file

C    num columns in .dist file
      IROW = 10
C     NQ=80  -> 9 for below
      LR = NQ/IROW + 1
C     RL is equal to 9
      RL = FLOAT(NQ)/IROW + 1
C    DIF = 0
      DIF = RL - LR
C    LR is equal to 8
      IF (DIF.LT.0.001) LR = LR - 1

      DO L=1,LR
       K1 = 1 + (L-1)*IROW
       K2 = K1 + IROW - 1
       IF (L.LT.LR) then
        read(17, 880) ((USOL(k,i),K=K1,K2),i=1,nz)
       ELSE
          K2 = NQ
C       this only occurs if NQ is a multiple of 10
C         if not, code needs to generate a dynamic format statement WARNING?
          if (K2-K1 .EQ. 9) then
            read(17, 880) ((USOL(k,i),K=K1,K2),i=1,nz)
          else
           fmtstr='(  E17.8)'
C           OK for one character as this should always be <10
           write(fmtstr(2:3),'(I1)')K2-K1+1
           read(17, fmtstr) ((USOL(k,i),K=K1,K2),i=1,nz)
          endif
       ENDIF
      enddo

C   This final one is CO2 number density
      read(17,881) (T(i),EDD(i),DEN(i),O3(i), SL(LCO2,i),i=1,nz)
C-SL(NSP-1) will be the density of the final short-lived species
C if CO2 is removed from the list. I dont THINK it will matter now
C that I call DOCHEM with the -1 before starting.... WARNING
! Use LCO2 instead of NSP-1 if CO2 specifically is desired...

C gna - we need to make it so that T = T_new
      IF(ICOUPLE.EQ.1) THEN
         DO I=1, NZ
               T(I) = T_new(I)
         END DO
         print *, 'T0 is'
         print *, T(1)
      ENDIF
      if (NP.gt.0) then
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
      else
        print*,'Warning: NP = 0, so no particles are being'
        print*,'assumed in this model. Proceed with caution...'
        print*,'(NP=0 to be used in giant planet templates only.)'
      endif


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
C sgflux, vdep, smflux, and distributed fluxes are already set
        do i=1,nq
         if (LBOUND(i).EQ.1) USOL(i,1)=fixedmr(i)
        enddo
      endif


C added by giada
      OPEN(unit=999, file='COUPLE/coupling_params.out')
 909  FORMAT(1X, F4.2, 5X, F8.3, 5X, F3.1, 5X, A8, 5X, I2,
     &     9X, I4, 6X, F4.2, 6X, F7.3)
 908  FORMAT(1X, 'timega', 6X, 'P0', 8X, 'frak', 3X, 'pstar', 8X,
     &   'ihztype', 6X, 'NZ', 6X, 'FSCALE', 6X, 'G')
      print *, 'FRAK = ', frak
      print *, 'P0 = ', P0
      print *, 'pstar = ', pstar
      if(NP.eq.4) print *, 'ihztype = ', ihztype
      print *, 'NZ =', NZ
      print *, 'FSCALE = ', FSCALE
      print *, 'G = ', G
      if(NP.eq.4) ihz = ihztype
      if(NP.lt.4) ihz = 99
      print *, 'ihz = ', ihz
      print *, 'NP =', NP
      WRITE(999,908)
      WRITE(999,909) timega, P0, frak, pstar, ihz, NZ, FSCALE,G

C
C
c-mc
c-mc  this block reads in a parameter file which, for the moment, is just
c-mc  a multulplier for the O2 or CH4 concentration
c-mc

       open(10, file='PHOTOCHEM/INPUTFILES/params.dat',status='OLD')
       read(10,303) o2mult, ch4mult
 303   format(F8.1,1x,F8.1)


c       do i=1,nq
c          print *, i, ISPEC(i)
c       enddo
c       stop

c       print *,USOL(LO2,1),2**ch4mult,2.0**ch4mult,
c     $       USOL(LO2,1)/(2.0**ch4mult)

c        USOL(LO2,1) = USOL(LO2,1)/(2.0**ch4mult)

c        print *, USOL(LO2,1)
c        stop


      do I=1,NZ
c-mc fix this stuff when all done (or maybe think of a better way to perturb)
C       WARNING ABOVE

!commenting these out for the HCL run

!        USOL(LCH4,I) = USOL(LCH4,I)*o2mult
c        USOL(LO2,I) = USOL(LO2,I)/2**ch4mult


c - fix this for new lbc thing...
c      HCLFLUX=1.e8*o2mult


        Den(i) = 1.0*Den(i)
      enddo

C ***** SET MODEL PARAMETERS *****
C     LTIMES = COUNTER FOR PHOTORATE SUBROUTINE
C     DT = INITIAL TIME STEP
C     TSTOP = TIME AT WHICH CALCULATION IS TO STOP
C     NSTEPS = NUMBER OF TIME STEPS TO RUN (IF TSTOP IS NOT REACHED)
C     FO2 = ground level O2 mixing ratio used in heritage calculations

      LTIMES = 0
C   fill up a heritage constant, eventually this should be purged. WARNING
c-mab: These values are needed for molecular weight computation in DIFCO, DENSTY, PHOTO.
      FO2 = USOL(LO2,1)
c-mab: Initializing these additional things below for use in giant planet template.
      FCO = USOL(LCO,1)
      FH2O = USOL(LH2O,1)
      FH2 = USOL(LH2,1)
      if(FH2.gt.0.50)FH2 = (1.0-FHE-FH2O-FCO-FCO2-FCH4) !do (1-everything)
c-mab: Note: this last adjustment only necessary if LBC = 1 for H2.
      FH = USOL(LH,1)
      FOH = USOL(LOH,1)
      FCH4 = USOL(LCH4,1)


C  CALL below sets up vertical grid
      CALL PHOTGRID

C   height index for the tropopause (-1 given the staggered grid)
C      the 1 in the second postion tells minloc to return a scalar
c-mab: JTROP is set to a constant value in planet.dat for the giant templates...
      IF (PLANET.NE.'WASP12B')JTROP=minloc(Z,1, Z .ge. ztrop)-1



c       print *, 'enter surface T'
c       read(*,'(F3.0)') TINC
c       print *, TINC

C    Swtich Around Tempterature
      do i=1,NZ
       if (i.le. JTROP) then
          zkm=Z(i)/1.e5
C        in Jim's original model BELOW WARNING
c         T(i) = 273. - 6.5*i - 0.25*i*i
C        what we are used for the sulfur paper
c         T(i) = 288. - 6.5*i - 0.25*i*i
C        what we are using for now - 286K at 0.5 km.
C             - i dont think this was the intended effect
c         T(i) = 293. - 6.5*i - 0.25*i*i

C This somewhat mimics what's above - where T was computed on a 1km grid that
C      didn't match up directly with the Z grid at 0.5,1.5, etc.
c         T(i) = 290. - 6.5*zkm - 0.25*zkm**2

C     testing a return to cold Archean
c          T(i) = Tinc - 6.5*zkm - 0.25*zkm**2
C     testing a return to cold Archean
c          print *, 289 - 6.5*zkm - 0.25*zkm**2

C   this somewhat mimics what's above - where T was computed on a 1km grid that
C          didn't match up directly with the Z grid at 0.5,1.5, etc.
c         T(i) = 293. - 6.5*zkm - 0.25*zkm**2
C   ok, this is close to the T profile that we were using in the NZ=80 case
c         print *,zkm, T(i)
C         T(i) = 290. - 6.5*i - 0.25*i*i  ! TD test - uncommenting for now

c        T(i) = 233. - 6.2*i - 0.25*i*i  ! cold earth ???
c        T(i) = 260. - 7.0*i !- 0.25*i*i  ! cold earth ???
c         Edd(i) = 1.E+5
c      else if (i.le. 3*jt) then
c        T(i) = 160. + 3.0*(i-jt)
       else
c        T(i) = 226. - 1.0*(i-3*jt)
c         T(i) = T(JTROP)
         Edd(i) = 1.0*Edd(i)
c        Edd(i) = 1.E+5
       endif
c       print *, Z(I),T(I)
      enddo
c      stop

ctemp
c      do i=jtrop,NZ
c      do i=136,NZ
c       edd(i)=1.3*edd(i)
c      enddo

      CALL DENSTY

c-mab: Note: Loop below is temporary....
      IF (PLANET.EQ.'WASP12B') THEN
	CALL RATESHJS
        PRINT*, "CALLING RATESHJS..."
      ELSE
        CALL RATES
        PRINT*, "CALLING RATES..."
      ENDIF

C    computes diffusion coefficents (K*N) and binary
C          diffusion coefficents for H and H2
      CALL DIFCO
      IF(PLANET.NE.'WASP12B')CALL PHOTSATRAT(JTROP,H2O)
      IF(PLANET.EQ.'WASP12B')THEN
      	DO I=1,NZ
      		P(I)=(1e-6)*PRESS(I)
      	ENDDO
      ENDIF
C    IDO=-1, fill up SL for accurate calculation on first timestep
      CALL DOCHEM(FVAL,-1,JTROP,iSL,USETD)

      if (PLANET .EQ. 'EARTH') then
C    current column integrated NO production rate on Earth
       PRONO = PRONO/1.


C                   ! divide by 1000 turns off ltning
c       PRONO = PRONO/1.e6 ! ATACAMA
      else if (PLANET .EQ. 'MARS') then
C      !divide by 1e9 turns off lightning for dry mars
       PRONO = PRONO/1.E9
      endif

C    !i.e if constant mr or constant flux UBC
      if (mbound(LH2) .gt. 0) then
        do i=1,nz
C        !don't use molecular diffusion
            bHN2(i) = 0.0
            bH2N2(i) = 0.0
            do j=1,NSP
              bX1X2(j,i) = 0.0 !Generalized form for giant atmospheres
            enddo
C           !don't use molecular diffusion
c           bXN2(i) = 0.0
        enddo
      else
C      !use effusion velocity formulation of diffusion limited flux
        Veff(LH) = 1.0*bhN2(nz)/DEN(NZ)
C      !diff lim flux
     $     *(1./Hscale(nz) - 1./scale_H(LH,nz))
        if(PLANET.EQ.'WASP12B')Veff(LH) = 0.0
        Veff(LH2) = 1.0*bH2N2(nz)/DEN(NZ)
     $     *(1./Hscale(nz) - 1./scale_H(LH2,nz))
        if(PLANET.EQ.'WASP12B')Veff(LH2) = 0.0
      endif

!gna - added coupling stuff for water here (just below tropopause)
c-mab: Executing this only for terrestrial planets based on FH2 prevalance.
      IF (FH2.LT.0.50) THEN
       do J=1,JTROP
        IF(ICOUPLE.eq.0) THEN
C sets H2O to relative humidity in troposphere
         USOL(LH2O,J) = H2O(J)
        ELSE
C     !set to h2o from clima if coupling on
         USOL(LH2O,J) = water(J)
        ENDIF
       enddo
      ENDIF


      IF (PLANET .eq. 'EARTH') CALL LTNING
c-mc: Lightning routine works for N2/CO2 atmospheres with variable O2
C     !makes table of vapor pressures for H2O and H2SO4
      CALL AERTAB
      NZ1 = NZ - 1
      HA = HSCALE(NZ)
C     The count of calls to rainout
      NRAIN = 0
      DT = 1.E-6
      DTINV = 1./DT
      TIME = 0.

c      NSTEPS = 681
C    for runs that are unstable...
c      TSTOP = 1.E14


      TSTOP = 1.E17
C     AVB changed to 1 for Earth+CL debugging
c      Commenting this out for now. We may use this later - but likely not in the main FORTRAN code.
c-mc - I'd be happy for this to be turned on when IDEBUG=1
c      PRINT*,' '
c      PRINT*,'Do you want to run for a single time step only (NSTEPS=1)'
c      PRINT*,'instead of till convergence (NSTEPS=10000)? [Note:'
c      CALL prompt('Default=n, i.e. use if debugging only.] (y/n?):')
c      READ(5,2)buffer
c2     FORMAT(a)
c      if(buffer(1:3).eq.'yes'.or.buffer(1:3).eq.'YES'.
c     &      or.buffer(1:1).eq.'y'.or.buffer(1:1).eq.'Y') then
c       NSTEPS = 1 !To run a single time step only without convergence
c      else
c       NSTEPS = 10000 !the default, to allow converging runs
c      endif
C      Default number of steps is 50,000. The code shouldn't take nearly this long to run except hot planets.
       NSTEPS = 10000

c-mab: nsteps = 1 recommended for initial model debugging
c      NSTEPS = 1
C     for standalone mode this should probably be the default
C      ICOUPLE = 1


C ***** write OUT INITIAL DATA *****
      CALL OUTPUT(0,NSTEPS,0.D0,jtrop, vdep,USOLORIG,USETD, frak)


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

C  First order molecular diffusion terms below
      do i=1,nq
        do j=1,nz
           ADU(i,j) = 0.0
           ADL(i,j) = 0.0
           ADD(i,j) = 0.0
        enddo
      enddo

c the following implements molecular diffusion - under construction

c  for everything other than H2, H
c     do i=1,nq
c       DU(i,1) = DU(i,1) + bXN2(1)/Den(1)/DZ/DZ
c       ADU(i,1) = bXN2(1)/Den(1)/DZ/2.*
c    6      (1./scale_H(i,1)-1./H_atm(1))
c       DL(i,NZ) = DL(i,NZ) + bXN2(nz1)/Den(nz)/DZ/DZ

c       ADL(i,NZ) = -bXN2(nz1)/Den(nz)/DZ/2.*
c    6      (1./scale_H(i,nz1)-1./H_atm(nz1))
c       DD(i,1) = DU(i,1)
c       ADD(i,1) = -ADU(i,1)
c why not DD, ADD for NZ?  NOT USED!
c     enddo

c     do i=1,nq
c interior grid points
c       do j=2,nz1
c           DU(i,j) = DU(i,j) + bXN2(j)/Den(j)/DZ/DZ
c           ADU(i,j) = bXN2(j)/Den(j)/DZ/2.*
c    6            (1./scale_H(i,j)-1./H_atm(j))
c           DL(i,j) = DL(i,j) + bXN2(j-1)/Den(j)/DZ/DZ
c           ADL(i,j) = -bXN2(j-1)/Den(j)/DZ/2.*
c    6            (1./scale_H(i,j-1)-1./H_atm(j-1))
c           DD(i,j) = DU(i,j) + DL(i,j)
c           ADD(i,j) = -ADU(i,j) - ADL(i,j)
c       enddo
c     enddo

c  for H and H2
c  lower boundary condition

      if(PLANET.EQ.'WASP12B'.or.PLANET.EQ.'EARTH') then
       do k=1,(NQ-NP)
       if(mbound(k).eq.0) then
        DU(k,1) = DU(k,1) + bX1X2(k,1)/Den(1)/DZ(1)**2
        ADU(k,1) = bX1X2(k,1)/Den(1)/DZ(1)/2.*
     6      (1./scale_H(k,1)-1./H_atm(1))
c upper boundary condition
        DL(k,NZ) = DL(k,NZ) + bX1X2(k,nz1)/Den(nz)/DZ(NZ)**2
        ADL(k,NZ) = -bX1X2(k,nz1)/Den(nz)/DZ(nz)/2.*
     6      (1./scale_H(k,nz1)-1./H_atm(nz1))
c  unused...
        DD(k,1) = DU(k,1)
        ADD(k,1) = -ADU(k,1)

c interior grid points   ?fixed 8-13-05
        do j=2,nz1
            DU(k,j) = DU(k,j) + bX1X2(k,j)/Den(j)/DZ(j)**2
            ADU(k,j) = bX1X2(k,j)/Den(j)/DZ(j)/2.*
     6            (1./scale_H(k,j)-1./H_atm(j))
            DL(k,j) = DL(k,j) + bX1X2(k,j-1)/Den(j)/DZ(j)**2
            ADL(k,j) = -bX1X2(k,j-1)/Den(j)/DZ(j)/2.*
     6            (1./scale_H(k,j-1)-1./H_atm(j-1))
            DD(k,j) = DU(k,j) + DL(k,j)
            ADD(k,j) = -ADU(k,j) - ADL(k,j)
        enddo
        endif
       enddo

      else
       if (mbound(LH2).eq.0) then
c diff limited flux implemented as effusion velocity
        DU(LH,1) = DU(LH,1) + bHN2(1)/Den(1)/DZ(1)**2
        ADU(LH,1) = bHN2(1)/Den(1)/DZ(1)/2.*
     6      (1./scale_H(LH,1)-1./H_atm(1))
        DU(LH2,1) = DU(LH2,1) + bH2N2(1)/Den(1)/DZ(1)**2
        ADU(LH2,1) = bH2N2(1)/Den(1)/DZ(1)/2.*
     6      (1./scale_H(LH2,1)-1./H_atm(1))
c upper boundary condition
        DL(LH,NZ) = DL(LH,NZ) + bHN2(nz1)/Den(nz)/DZ(NZ)**2
        ADL(LH,NZ) = -bHN2(nz1)/Den(nz)/DZ(nz)/2.*
     6      (1./scale_H(LH,nz1)-1./H_atm(nz1))
        DL(LH2,NZ) = DL(LH2,NZ) + bH2N2(nz1)/Den(nz)/DZ(NZ)**2
        ADL(LH2,NZ) = -bH2N2(nz1)/Den(nz)/DZ(NZ)/2.*
     6      (1./scale_H(LH2,nz1)-1./H_atm(nz1))
c  unused....
        DD(LH,1) = DU(LH,1)
        ADD(LH,1) = -ADU(LH,1)
        DD(LH2,1) = DU(LH2,1)
        ADD(LH2,1) = -ADU(LH2,1)

c interior grid points   fixed 8-13-05
        do j=2,nz1
            DU(LH,j) = DU(LH,j) + bHN2(j)/Den(j)/DZ(j)**2
            ADU(LH,j) = bHN2(j)/Den(j)/DZ(j)/2.*
     6            (1./scale_H(LH,j)-1./H_atm(j))
            DL(LH,j) = DL(LH,j) + bHN2(j-1)/Den(j)/DZ(j)**2
            ADL(LH,j) = -bHN2(j-1)/Den(j)/DZ(j)/2.*
     6            (1./scale_H(LH,j-1)-1./H_atm(j-1))
            DU(LH2,j) = DU(LH2,j) + bH2N2(j)/Den(j)/DZ(j)**2
            ADU(LH2,j) = bH2N2(j)/Den(j)/DZ(j)/2.*
     6             (1./scale_H(LH2,j)-1./H_atm(j))
            DL(LH2,j) = DL(LH2,j) + bH2N2(j-1)/Den(j)/DZ(j)**2
            ADL(LH2,j) = -bH2N2(j-1)/Den(j)/DZ(j)/2.*
     6             (1./scale_H(LH2,j-1)-1./H_atm(j-1))
            DD(LH,j) = DU(LH,j) + DL(LH,j)
            ADD(LH,j) = -ADU(LH,j) - ADL(LH,j)
            DD(LH2,j) = DU(LH2,j) + DL(LH2,j)
            ADD(LH2,j) = -ADU(LH2,j) - ADL(LH2,j)
        enddo
       endif  !end molecular diffusion for H and H2 loop
      endif !End loop with planet-based distinction



      do I=1,nz
c         write(10, 1210) Z(I),DU(1,i),DU(LH,i),DU(LH2,i),DL(1,i),
c    5     DL(LH,i),DL(LH2,i),DD(nq,i),DD(LH2,i)
c          write(10, 1210) Z(I),scale_H(1,i),scale_H(LH,i),
c     5     scale_H(LH2,i), Hscale(i)
      enddo
c 1210 FORMAT(1X,1P12E10.2)

c  end of not yet ready to test section
C
      KD = 2*NQ + 1
      KU = KD - NQ
      KL = KD + NQ
C
C   write OUT RESULTS EVERY NPR TIME STEPS
      NPR = 1000
      PRN = NPR

C   DO PHOTORATES EVERY MP TIME STEPS
      NPHOT = 0
C    integer since starts with M
      MP = 1
C    float since starts with P
      PM = MP
      NN = 0
C
C ***** START THE TIME-STEPPING LOOP *****
C   calculation is stopped if it hasn't converged by NSTEPS
      DO 1 N=1,NSTEPS
      TIME = TIME + DT
C   counter for number of timesteps
      NN = NN + 1
C   both are integers, so answer is modular(0 for N=1-3 , 1 for N=4-6, etc.)
      MS = (N-1)/MP
C   SM=0,1/3,2/3,1,4/3,5/3,2......
      SM = (N-1)/PM
      IF (NN.EQ.NSTEPS) SM = MS
C   To skip PHOTO
      IF (SM-MS.GT.0.01) goto 18
C   Doing photo on first time step, not again until 100 seconds
      IF (N.GT.1 .AND. TIME.LT.1.E2) goto 18



C store mixing ratio of all species that take place in photolysis reactions
C these are the absorbers which block solar radiation
C If S8 occurs in the gas phase, we use Andy Youngs method to calculate S8
C  photorates which requires care

         lpolyscount=1
      do k=1,kj

       do i=1,nz

C     this gets any SL/IN species that have photorates
          if (photoreac(k).gt.nq) then

             absorbers(k,i)=ABS(SL(INT(photoreac(k)),i))/DEN(I)

c           if (i.eq.1) print *, k,photoreac(k),ISPEC(INT(photoreac(k))),
c     $       absorbers(k,i)

C     quasi-hardcoded S8 behavior WARNING
           else if (ISPEC(INT(photoreac(k))).eq.'S8      ') then

              absorbers(k,i)=ABS(USOL(INT(photoreac(k)),i))
              if (lpolyscount .eq. 2) absorbers(k,i)=0.0
C           S8R doesn't really exist
C           S8L doesn't really either, but this is how we are rolling for now...
C           we are filling up absorbers for S8R and S8, but S8 doesn't
C           have a cross section
C           so only S8R is actually absorbing photons here in the RT scheme.
              if (i.eq.nz) then
c            print *, k,photoreac(k),ISPEC(INT(photoreac(k))),lpolyscount
                 lpolyscount=lpolyscount+1
              endif
          else
           absorbers(k,i)=ABS(USOL(INT(photoreac(k)),i))
          endif
       enddo
      enddo


c      print *, 'stopping in main'
c      stop

C Fill up some heritage vectors rather than try to extract them from the code
C I am assuming here that the code will always have these species
C  in it so no if's are needed

       JH2O=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'H2O    ')
       JO2_O1D=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'O2    ')
       JO3_O1D=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'O3    ')
C    Finds the first entry in photoreac for the given species
       JCO2=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'CO2    ')


      do I=1,NZ
       H2O(I) = absorbers(JH2O,I)
       O2(I) =  absorbers(JO2_O1D,I)
       O3(I) = absorbers(JO3_O1D,I)
C orig       CO2(I) = FCO2
       CO2(I) = absorbers(JCO2,I)
      enddo



      IDO = 0
      IF (NN.EQ.NSTEPS) IDO = 1
      CALL PHOTO(ZY,AGL,LTIMES,ISEASON,IZYO2,IO2,INO,IDO,timega,frak,
     &      pstar,ihztype,uvscale)

      IF (NAQ.GT.0) CALL RAINOUT(JTROP,NRAIN,USETD)  !ok




      CALL AERCON



c      print *, photoreac
c      stop

       JO1D=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'O1D    ')
       JO2=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'O2     ')

c      print *, 'stopping after photo'
c      stop

C
C   TIME-DEPENDENT BOUNDARY CONDITIONS
C     (NOTICE THAT THE INCLUSION OF
C     A FALL VELOCITY FOR PARTICLES IS MATHEMATICALLY EQUIVALENT TO
C     INCREASING THE DEPOSITION VELOCITY AT THE SURFACE.  AT THE TOP
C     BOUNDARY, ZERO FLUX IS ACHIEVED BY SETTING THE EFFUSION VELO-
C     CITY EQUAL TO THE FALL VELOCITY)

c Shawn's code has some changes to h escape terms here. Consider porting.
C -Mark C
c Using variables for the particle indeces to kill warnings. We should
c eventually migrate this to the read-in. But that will require more
c testing, or someone to figure out how to do that in Mark's input
c framework. -Shawn D-G       WARNING


      if(USETD.EQ.0.AND.NP.EQ.0) then
c-mab Time-dependent boundary conditions for particles set to 0
        NPSO4 = 0
        NPS8 = 0
        NPHC = 0
        NPHC2 = 0
      elseif(USETD.EQ.0.AND.NP.GT.0) then
        NPSO4 = LSO4AER - NQ + NP
        NPS8 = LS8AER - NQ + NP
        NPHC = LHCAER - NQ + NP
        NPHC2 = LHCAER2 - NQ + NP

        VDEP(LSO4AER) = VDEP0(LSO4AER) + WFALL(1,NPSO4)
        VEFF(LSO4AER) = VEFF0(LSO4AER) + WFALL(NZ,NPSO4)

        VDEP(LS8AER) = VDEP0(LS8AER) + WFALL(1,NPS8)
        VEFF(LS8AER) = VEFF0(LS8AER) + WFALL(NZ,NPS8)

        if(NP.ge.3) then
          VDEP(LHCAER) = VDEP0(LHCAER) + WFALL(1,NPHC)
          VEFF(LHCAER) = VEFF0(LHCAER) + WFALL(NZ,NPHC)

          if (NP.eq.4) then
            VDEP(LHCAER2) = VDEP0(LHCAER2) + WFALL(1,NPHC2)
            VEFF(LHCAER2) = VEFF0(LHCAER2) + WFALL(NZ,NPHC2)
          endif
        endif
      else
        NPSO4 = LSO4AER - NQ
        NPS8 = LS8AER - NQ
        NPHC = LHCAER - NQ
        NPHC2 = LHCAER2 - NQ
      endif



c estimate CO2 photolysis above the top of the grid
C and return CO + O to the upper grid point
c NOTE: this behavior is turned on and off by setting MBOUND=2 for CO2 in species.dat

c   JCO2 already defined above
      JCO2_O1D = JCO2+1
      VCO2 = (prates(JCO2,NZ) + prates(JCO2_O1D,NZ)) * HA
      SMFLUX(LO) = - VCO2*CO2(NZ)*DEN(NZ)
      SMFLUX(LCO) = SMFLUX(LO)


      NMP = NSTEPS - MP
      if ((NN/NPR)*NPR.eq.NN .or. NN.eq.1) then

C  PHOTOLYSIS RATES FORMATTED I/O AND PRINTOUT

C  Number of columns of photorates
      IROW = 8
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
C       OK for one character as this should always be <10
           write(fmtstr(13:14),'(I1)')K2-K1+1
           write(14, fmtstr) (photolabel(K),K=K1,K2)

           fmtstr2='(  (1PE9.2,2X))'
C       +2 is here so that it writes the Z a well
           write(fmtstr2(2:3),'(I2)')K2-K1+2
           write(14, fmtstr2) (Z(I),(prates(k,i),K=K1,K2),i=1,nz,3)
       ENDIF
      enddo



 883  format(9(1PE9.2,2X))
 884  format(/5X,'Z',6X,8A11)
 885  FORMAT(//1X,'PHOTOLYSIS RATES')

C   end photorates printout
      endif

C   start here if we are skipping PHOTO
 18   continue


C
      if (NP.GT.0) then   !particles in main loop
       CALL SEDMNT(FSULF,USETD,frak,HCDENS,ihztype)

C    particles in main loop
       if (USETD.EQ.0) then

      do J=1,NZ
      do JJ=1, NP

            if (JJ.eq.1)  nparti = LSO4AER
            if (JJ.eq.2)  nparti = LS8AER
            if (JJ.eq.3)  nparti = LHCAER
            if (JJ.eq.4)  nparti = LHCAER2


            AERSOL(J,JJ) = USOL(nparti,J)*DEN(J)/(CONVER(J,JJ))


c O2 CODE CHANGES
c gna: the way it's done below, you have to COMMENT OUT these lines to remove
c particles from the code.  Aaak!  The lines above are maybe not perfect
c but better than commenting out stuff when the input files change...

C-ACK Hardcoded particle numbers below
c      AERSOL(J,2) = USOL(LS8AER,J)*DEN(J)/CONVER(J,2)
c      AERSOL(J,3) = USOL(LHCAER,J)*DEN(J)/CONVER(J,3)
c      AERSOL(J,4) = USOL(LHCAER2,J)*DEN(J)/CONVER(J,4)
c      print *, AERSOL(J,1),AERSOL(J,2),USOL(LSO4AER,J),CONVER(J,1),
c     &         USOL(LS8AER,J),CONVER(J,2)

C   conver is the number of molecules/particle, so that main calcuations
C     are done in molecule space
C      molecules/cm3 *#particles/molecule - > AERSOL [# particles/cm3]
      enddo
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
       elseif (USETD.GT.0) then
C   particles in tri-diag
        do J=1,NZ
         do k=1,NP
           AERSOL(J,K)=PARTICLES(J,K)*DEN(J)/CONVER(J,K)
         enddo
        enddo
       endif

       else
       print*,'Note: Since NP = 0, did not call SDMNT...'

       endif





C ***** SET UP THE JACOBIAN MATRIX AND RIGHT-HAND SIDE *****

c-mc gonna start here and try no to mess anything up - point is to
c-mc feedback the solution on itself
c-mc first go is to try to just repeat a few times - if this works,
c-mc then set up a check for convergence

c-mc OK seems to work.  4 iterations seems to provide a converged solution.
c-mc test 1  - check the solution with 4 iterations against that of
c-mc feeding back the code 4 times...
c-mc  this test looks good - see 'changetesting' sheet in ~/keff/redox.xls
c-mc next up is to save some intermediate output to see how much things are
c-mc changing in subsquent iterations
c-mc  I should try to identify an epsilon that satisfies a similar global
c-mc condition to that of just repeating the loop four times...

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
c-mACK expand me... WARNING

C



C   COMPUTE CHEMISTRY TERMS AT ALL GRID POINTS
      IDO = 0
C   computes TP, TL
      if ((NN/NPR)*NPR.eq.NN) IDO = 1
      IF (NN.EQ.NSTEPS) IDO = 1

C   IDO=1 happens on last step- computes total production and loss.
      CALL DOCHEM(FVAL,IDO,JTROP,iSL,USETD)

      DO 9 I=1,NQ
      DO 9 J=1,NZ
      K = I + (J-1)*NQ
      RHS(K) = FVAL(I,J)
      if (ITERATE.eq.1) USAVEOLD(I,J) = USOL(I,J)
C-mc testing 4/29/06  used to revert if timestep is too big.
C-mc and also for USOLPREV once the timestep succeedes
C  original code  - used as part of the reverse euler solver
   9  USAVE(I,J) = USOL(I,J)

C


C     NEW CODE FROM EDDIE
C   Loop through all chemical species
      DO 3 I=1,NQ
C   Loop through all vertical atmospheric layers
      DO 11 J=1,NZ
C   as it was - USOL should be positive here anyway,
C       can probably remove this (Eddie)
c     R(J) = EPSJ * ABS(USOL(I,J))
C   R(J) is value to perturb USOL(I,J) by in Jacobian calculation,
C       EPSJ much less than 1
      R(J) = EPSJ * USOL(I,J)
C   This is my debug in other version of code - Eddie !!!  WARNING
C           ! PERTURB DEBUG !!!
C   Below ensures no USOL(I,J) falls below double precision limit !!!
      IF(R(J).LT.1.e-100) R(J) = 1.e-100
c   Add perturbing quantity to mixing ratio
  11  USOL(I,J) = USAVE(I,J) + R(J)
C   Call the photochemistry routine
C      FV has dimension (NQ,NZ) and holds gas densities
      CALL DOCHEM(FV,0,JTROP,iSL,USETD)
c
c end new code from eddie

c old code
c      DO 3 I=1,NQ
c      DO 11 J=1,NZ
c     R(J) = EPSJ * ABS(USOL(I,J))   !as it was
c      R(J) = EPSJ * USOL(I,J)
c  11  USOL(I,J) = USAVE(I,J) + R(J)
c      CALL DOCHEM(FV,0,JTROP,iSL,USETD)


C
      DO 12 M=1,NQ
      MM = M - I + KD
      DO 12 J=1,NZ
      K = I + (J-1)*NQ
C  -J since its orig - perturbed
  12  DJAC(MM,K) = (FVAL(M,J) - FV(M,J))/R(J)
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

c-mc ok, I need to verify that this is adding in -J in all DJAC calls
c-mc (check signs above?)     WARNING?
c-mc ok - Jacobian for transport diagonals is:
c-mc J~Chem-DD. We want -J and have already
c-mc filled with -CHEM, so adding DD is appropriate.
c-mc DTINV is the extra term in the main diagonal
c-mc J_upper=DU and J_lower=DL, so DJAC (which is -J)
c-mc uses -DU and -DL respectivly


      if(USETD.EQ.0.and.NP.gt.0) then  !particles in main loop
C ack - these need to be abstracted... WARNING
C   ADD ADVECTION TERMS FOR PARTICLES

      do L=1,NP
      if(L.EQ.1) I = LSO4AER
      if(L.EQ.2) I = LS8AER
      if(L.EQ.3) I = LHCAER
      if(L.EQ.4) I = LHCAER2
       do J=2,NZ1
        K = I + (J-1)*NQ
        RHS(K) = RHS(K) + DPU(J,L)*USOL(I,J+1) - DPL(J,L)*USOL(I,J-1)
        DJAC(KU,K+NQ) = DJAC(KU,K+NQ) - DPU(J,L)
        DJAC(KL,K-NQ) = DJAC(KL,K-NQ) + DPL(J,L)
       enddo
      enddo


      endif
C
C ***** LOWER BOUNDARY CONDITIONS *****
      DO 15 K=1,NQ
      U(K) = USOL(K,1)
C   OK as long as we don't model atmospheric Boron (WARNING?)
      LB = LBOUND(K)

      if (LB.eq.0 .OR. LB.eq.3) then
C       CONSTANT DEPOSITION VELOCITY/SPECIES WITH DISTRIBUTED FLUXES
        RHS(K) = RHS(K) + DU(k,1)*(USOL(K,2) - U(K))
     2   + ADU(k,1)*(USOL(K,2) + U(K)) - VDEP(K)*U(K)/DZ(1)
        DJAC(KD,K) = DJAC(KD,K) +DTINV +DU(k,1) -ADU(k,1) +VDEP(K)/DZ(1)
        DJAC(KU,K+NQ) = - DU(k,1) - ADU(k,1)
c  is this right for particles?? WARNING

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
c       Constant mixing ratio at the top. not debugged WARNING
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
c   why am I doing this for S8?? WARNING
c  turn it off for S8
C  Jim apparently was prepared to do this for many species
      DO 34 I=1,1
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



C distributed (volcanic) sources

      do i=1,nq
       if (LBOUND(i).eq.3) then
C     convert to cm
        disth=distheight(i)*1.e5
C     height index (-1 given the staggered grid)
C     the 1 in the second postion tells minloc to return a scalar
        jdisth=minloc(Z,1, Z .ge. disth)-1


        ZTOP=Z(jdisth)-Z(1)
        ZTOP1=Z(jdisth)+0.5*DZ(jdistH)
c        print *, ISPEC(i),distH,jdistH,ZTOP,ZTOP1

C     distribute from second level to distheight
        do j=2,jdisth
          K = i + (J-1)*NQ
          rhs(k) = rhs(k) + 2.*distflux(i)*(ZTOP1-Z(j))/(Den(j)*ZTOP**2)
        enddo
       endif
      enddo


C below is Kevin's test of a distributed sink.
C Keeping in case this is ever desired

C distribute H2SO4 sink uniformly between 0 and 20 km
c   I don't know how to keep track of the fluxes
c      jdistSO4 = 21
c      ZSO4 = Z(jdistSO4-1)
c        L = LH2SO4
c        do j=1,jdistso4
c          K = L + (J-1)*NQ
c          rhs(k) = rhs(k) - 0.01  ! per second
c        enddo



C testing device to toggle 1st and 2nd order time stepping
      revEu2=0

C Do 2nd order (first step always 1st order R.E.)
      if (N .gt. 1 .and. revEu2 .eq.1) then

c-mc first attempt at second order euler
c  [1/DT*I - a*J]*X = 1/DT*[b*f_{n} - c*f_{n-1} - ?fn] + e*d/dt{f_n}
c
c where a=2/3, b=4/3, c=1/3, d=?,e=2/3, X=f_{n+1} - f_{n}

c DJAC has alread been computed as 1/DT*I - J,
C  and d/dt(f_{n}) has been computed as RHS


      do j=1,lda
        do k=1,neq
          DJAC(j,k) = 2.0/3.0*DJAC(j,k)
        enddo
      enddo

c-mc could make this shorter by not multiplying the upper NQ rows
c-mc of djac (which are 0 by definition) WARNING

c renormalize DTINV
c the previous loop took the main diagonal
C (which is row KD in band form) to 2/3*1/DT - 2/3J
c so the following returns us to 1/DT - 2/3J

      do i=1,nq
       do j=1,nz
        K = I + (J-1)*NQ
        DJAC(KD,K) = DJAC(KD,K) + 1.0/3.0*DTINV
       enddo
      enddo

      do K=1,NEQ
         RHS(K) = 2.0/3.0*RHS(K)
      enddo

c RHS(K) where K=1,NEQ where NEQ=NQ*NZ, so need to make sure the right
C   usol/height combination is being added at each step:
      do I=1,NQ
       do J=1,NZ
        K = I + (J-1)*NQ
        RHS(K) = RHS(K) +  DTINV*1.0/3.0*(USOL(I,J)-USOLPREV(I,J))
       enddo
      enddo
C  End 2nd order reverse Euler correction loop
      endif



C ***** FACTOR THE JACOBIAN AND SOLVE THE LINEAR SYSTEM *****

c-mc  comment out original to get at condition number of the matrix
      CALL SGBFA(DJAC,LDA,NEQ,NQ,NQ,IPVT,INFO)
      IF(INFO.NE.0) then
         print *, N,INFO
         print *, 'ERROR in SGBFA'
         stop
      endif


c      CALL SGBCO(DJAC,LDA,NEQ,NQ,NQ,IPVT,RCOND,WORK)
c      IF(INFO.NE.0) then
c         print *, N,INFO
c         print *, 'ERROR in SGBFA'
c         stop
c      endif
c      print *, RCOND


c-mc  we have set up RHS so that DJAC*X=RHS
c-mc DJAC is now upper triangular...

      CALL SGBSL(DJAC,LDA,NEQ,NQ,NQ,IPVT,RHS,0)

c-mc  after this point, RHS changes to solution vector for DJAC*X = RHS
c-mc  i.e. "RHS" = X = f_{n+1} - f_{n}, so that f_{n+1} = f{n} + RHS

C   Height Index (-1 given the staggered grid) Below
      J15=minloc(Z,1, Z/1.e5 .ge. 15)-1
      J25=minloc(Z,1, Z/1.e5 .ge. 25)-1
      J70=minloc(Z,1, Z/1.e5 .ge. 70)-1
      J50=minloc(Z,1, Z/1.e5 .ge. 50)-1

C
C   COMPUTE NEW CONCENTRATIONS (IGNORE ERRORS IN SEVERAL SPECIES
C    THAT VIRTUALLY DISAPPEAR UP HIGH)
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
c    50, this often causes problems causes the program to fail at 49.5 km!
c-orig        ELSEIF(I.EQ.LSO4AER .AND. J.GT.J25) THEN
c        IF (I.EQ.LSO4AER .AND. J.GT.J25) THEN

!ACK - should return to this now that
C      particles condensation is fixed up 

c      A less drastic measure
        IF (I.EQ.LSO4AER.AND. J.GT.J25) THEN
          USOL(I,J) = USOL(I,J) + RHS(K)
        ELSEIF(I.EQ.LS8AER.AND. J.GT.J25) THEN
          USOL(I,J) = USOL(I,J) + RHS(K)
        ELSEIF(I.EQ.LHCAER.AND. J.GT.J25) THEN
c        ELSEIF(I.EQ.LHCAER) THEN  ! !!!
          USOL(I,J) = USOL(I,J) + RHS(K)
        ELSEIF(I.EQ.LHCAER2.AND. J.GT.J25) THEN
c        ELSEIF(I.EQ.LHCAER2) THEN  ! !!!
          USOL(I,J) = USOL(I,J) + RHS(K)



c  Jim set this at 50 km.  program fails at 49.5 km.
c  if I turn it off the program fails at 63.5 km immediately
c   so I may conclude that there is an issue with the UBC on SO4AER
c  anyway, I'll set it to say 35 km
c
c   more generally I want to ignore errors in anything with
c   mixing ratios less than say 1e-20
c    50, this causes the program to fail at 49.5 km!
c       ELSEIF(I.EQ.LNO .AND. J.GT.70) THEN
c         USOL(I,J) = USOL(I,J) + RHS(K)
c       ELSEIF(I.EQ.LNO2 .AND. J.GT.70) THEN
c         USOL(I,J) = USOL(I,J) + RHS(K)
c       ELSEIF(I.EQ.LHNO .AND. J.GT.70) THEN
c         USOL(I,J) = USOL(I,J) + RHS(K)

        ELSEIF(I.EQ.LS4) THEN
                     USOL(I,J) = USOL(I,J) + RHS(K)

        ELSEIF(PLANET.EQ.'WASP12B'.AND.I.EQ.LO2) THEN !!!!!! Ignoring O2 entirely right now for hot jupi to get it to converge quickly
!!!! btw this was kept disabled in Ravi's CURRENT version of main
                     USOL(I,J) = USOL(I,J) + RHS(K)
C        ELSEIF(PLANET.EQ.'WASP12B'.AND.I.EQ.LO2) THEN
C        	IF(J.GT.21.AND.J.LT.71) THEN !!! Using this to ignore only PARTS of O2, i.e. more carefully, than entiety of O2 (above this)
C                     USOL(I,J) = USOL(I,J) + RHS(K)
C            ENDIF
!!!!!!!!!!!! The conditions below are necessary to make WASP12B converge !!!!  
c-mab: These from trial/error + Kopparapu et al. 2012 code.
c-mab: For the solar system templates, all species have the same USOL = 1e-20 limit.
c-mab: Similar species-specific conditions may be necessary to get future templates to converge.
c-mab: This is temporary till we can find a "cleaner" solution.                   
        ELSEIF(PLANET.EQ.'WASP12B'.AND.I.EQ.LH2CO .AND. J.GT. 61) THEN
                     USOL(I,J) = USOL(I,J) + RHS(K)
        ELSEIF(PLANET.EQ.'WASP12B'.AND.I.EQ.LH2COH .AND. J.GT. 35) THEN
                     USOL(I,J) = USOL(I,J) + RHS(K)
        ELSEIF(PLANET.EQ.'WASP12B'.AND.I.EQ.LCH .AND. J.GT. 31) THEN
                     USOL(I,J) = USOL(I,J) + RHS(K)
        ELSEIF(PLANET.EQ.'WASP12B'.AND.I.EQ.LCH23 .AND. J.GT. 31) THEN
                     USOL(I,J) = USOL(I,J) + RHS(K)
        ELSEIF(PLANET.EQ.'WASP12B'.AND.I.EQ.LCH3 .AND. J.GT. 50) THEN
                     USOL(I,J) = USOL(I,J) + RHS(K)
        ELSEIF(PLANET.EQ.'WASP12B'.AND.I.EQ.LCH4 .AND. J.GT. 56) THEN
                     USOL(I,J) = USOL(I,J) + RHS(K)
C        ELSEIF(PLANET.EQ.'WASP12B'.AND.I.EQ.LH2O .AND. J.GT. 75) THEN  !!!! btw this was kept disabled in Ravi's CURRENT version of main
C                     USOL(I,J) = USOL(I,J) + RHS(K) !!! Incorporating this in this version actually slightly increased convergence time
C        ELSEIF(PLANET.EQ.'WASP12B'.AND.I.EQ.LCO2) THEN  !!!! btw this was kept disabled in Ravi's CURRENT version of main
C                     USOL(I,J) = USOL(I,J) + RHS(K)  !(PS - woah ignoring CO2 entirely? Really? - disabling this actually slightly increased convergence as well, more than H2O alone did)
        ELSEIF(PLANET.EQ.'WASP12B'.AND.I.EQ.LCH3OH .AND. J.GT. 10) THEN
                     USOL(I,J) = USOL(I,J) + RHS(K)
C        ELSEIF(PLANET.EQ.'WASP12B'.AND.I.EQ.LO .AND. J.GT. 60) THEN !!!! btw this was kept disabled in Ravi's CURRENT version of main
C                     USOL(I,J) = USOL(I,J) + RHS(K)
        ELSEIF(PLANET.EQ.'WASP12B'.AND.I.EQ.LC .AND. J.GT. 60) THEN
                     USOL(I,J) = USOL(I,J) + RHS(K)
        ELSEIF(PLANET.EQ.'WASP12B'.AND.I.EQ.LCH3O .AND. J.GT. 50) THEN
                     USOL(I,J) = USOL(I,J) + RHS(K)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !!!!ELSEIF (USOL(I,J).LT. 1.E-15) THEN !TEMPORARILY DISABLING THIS ALTOGETHER TO USE RAVI'S CONDITIONS DIRECTLY
        ELSEIF(PLANET.NE.'WASP12B'.AND.USOL(I,J).LT. 1.E-19) THEN !DEFAULT ATMOS CONDITION -- USING FOR NON HOT JUP PLANETS

c-orig        IF (USOL(I,J).LT. 1.E-20) THEN
c        IF (USOL(I,J).LT. 1.E-22) THEN
c        ELSEIF (USOL(I,J).LT. 1.E-22) THEN
          USOL(I,J) = USOL(I,J) + RHS(K)
c
        ELSE
          REL(I,J) = RHS(K)/USOL(I,J)
          EREL = ABS(REL(I,J))
          EMAX = max(EMAX,EREL)
          IF(EREL.LT.EMAX) THEN
            USOL(I,J) = USOL(I,J) + RHS(K)
C       store info on species with largest error
          ELSE
            IS = I
C       mc -this label is OK, because S will never have a photolysis reaction
            JS = J
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
c-mc test this at some point down the road WARNING

C-mab: Foregoing this for giant planets since we don't fudge the values above...
      IF(FH2.LT.0.50) THEN
      DO 4 J=1,JTROP
        IF(ICOUPLE.eq.0)THEN
        USOL(LH2O,J) = H2O(J)
c       USOL(LS8,J) = S8S(J)
        ELSE
         USOL(LH2O,J) = water(J)
        ENDIF
   4  CONTINUE
      ENDIF

c      do i=1,nq
c       do j=1,nz
c          print *,Z(j), USOL(LSO4AER,j),USOL(LS8AER,j)
c         USOL(i,j)=abs(USOL(i,j))
c       enddo
c      enddo
C      temp


      if (USETD.EQ.0.AND.NP.GT.0) then
C switch around main loop particles
c 1.e-38 is the smallest number for single precision. We should upgrade
c this to double precision at some point. -Shawn D-G  WARNING
      smallest = 1.e-38
      lcountso4=0
      lcounts8=0
      lcountHC=0
      lcountHC2=0


      DO J=1,NZ
      if(USOL(LSO4AER,J).LT.0) lcountso4=lcountso4+1

      if(lcountso4.gt.0) then
       if (J.GT.1) then
        USOL(LSO4AER,J)=
     &  USOL(LSO4AER,J-1)*EXP(-WFALL(J,NPSO4)*DZ(J)/EDD(J))
       else
        USOL(LSO4AER,J)=-1
       endif
      endif

      if(USOL(LS8AER,J).LT.0) lcounts8=lcounts8+1
      if(lcounts8.gt.0) then
       if(J.GT.1) then
        USOL(LS8AER,J) =
     &  USOL(LS8AER,J-1) * EXP(-WFALL(J,NPS8)*DZ(J)/EDD(J))
       else
        USOL(LS8AER,J)=-1
       endif
      endif

      if (NP.GE.3) then
      if(USOL(LHCAER,J).LT.0) lcountHC=lcountHC+1

      if(lcountHC.gt.0) then
       if (J.GT.1) then
        USOL(LHCAER,J) =
     &  USOL(LHCAER,J-1)*EXP(-WFALL(J,NPHC)*DZ(J)/EDD(J))
       else
         USOL(LHCAER,J)=-1
       endif
      endif

      if (NP.EQ.4) then
      if(USOL(LHCAER2,J).LT.0) lcountHC2=lcountHC2+1

      if(lcountHC2.gt.0) then
        if (J.GT.1) then
         USOL(LHCAER2,J)=
     &    USOL(LHCAER2,J-1)*EXP(-WFALL(J,NPHC2)*DZ(J)/EDD(J))
        else
         USOL(LHCAER2,J)=-1
        endif
      endif
      endif
      endif

      USOL(LSO4AER,J) = max (USOL(LSO4AER,J),smallest)
      USOL(LS8AER,J) = max(USOL(LS8AER,J),smallest)
      if (NP.GE.3) then
      USOL(LHCAER,J) = max (USOL(LHCAER,J),smallest)
      USOL(LHCAER2,J) = max (USOL(LHCAER2,J),smallest)
      endif
      enddo

C end particle switching
      endif


      !test
      do i=1,nq
         do j=1,nz
c           USOL(i,j)=max(USOL(i,j),smallest)
           USOL(i,j)=abs(USOL(i,j))
         enddo
      enddo


      if(USETD.EQ.1.AND.NP.GT.0) then



*********TRIDIAG STARTS HERE

C ***** SOLVE FOR S8 AND SO4 PARTICLES USING A TRIDIAGONAL INVERSION *****
c-mc I haven't done the work to port the hc aerosols to the tri-diag.
C-mc this would need to be done if desired.

      DO 58 L=1,NP
C tridiag particles must appear right after LL species in species.dat
      I = NQ + L
C     Below is height index (-1 given the staggered grid)
      IF(I.EQ.LSO4AER) MZ = minloc(Z,1, Z/1.e5 .ge. 50)-1
      IF(I.EQ.LS8AER) MZ = minloc(Z,1, Z/1.e5 .ge. 40)-1
      IF(I.EQ.LHCAER) MZ = minloc(Z,1, Z/1.e5 .ge. 70)-1
      IF(I.EQ.LHCAER2) MZ = minloc(Z,1, Z/1.e5 .ge. 70)-1
C    at some point check/abstract these 40/50km assumptions
C      -this could be easily shunted to the species.dat file. WARNING
C      height index (-1 given the staggered grid)
c      IF(I.EQ.LSO4AER) MZ = minloc(Z,1, Z/1.e5 .ge. 79.5)-1

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

C why are there no dl*PARTICLES() in here?  WARNING
C all the other dl's are multiplied by USOL.
C just on RHS...

      vturb=0.01
c      vturb=0.000005
c      vturb=0.0

C
C   BOUNDARY CONDITIONS
      TA(MZ) = - DL(I,MZ) + DPL(MZ,L)
      TB(MZ) = TB(MZ) + DL(I,MZ) + 0.5*WFALL(MZ,L)/DZ(MZ)
C    Orig from Jim's code, as it was.
c      TB(1) = TB(1) + DU(I,1) + (0. - 0.5*WFALL(1,L))/DZ(1)
C    Orig from Jim's code, as it was.
c      TB(1) = TB(1) + DU(I,1) + (vturb - 0.5*WFALL(1,L))/DZ(1)
C    Testing Below
      TB(1) = TB(1) + DU(I,1) - (vturb+0.5*WFALL(1,L))/DZ(1)
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
c 1.e-38 is the smallest number for single precision. We should upgrade
c this to double precision at some point. -Shawn D-G
      smallest = 1e-38

C   FILL UP UPPER PORTION WITH APPROXIMATE ANALYTIC SOLUTION

      do J=MZP1,NZ
        PARTICLES(J,L)=PARTICLES(J-1,L) * EXP(-WFALL(J,L)*DZ(J)/EDD(J))
C   needed?  WARNING
        PARTICLES(J,L) = MAX(PARTICLES(J,L),smallest)
      enddo



   58 CONTINUE
C

      endif
*****END TRIDIAG
C continue doing newton steps (4 seems to work best)
C someday I should see if this is justified
C (r37:36 in /td branch has first attempt at convergence checking) WARNING
 73   continue






C   AUTOMATIC TIME STEP CONTROL
      DTSAVE = DT
c-mc these are the ones that kevin originally used...
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

Ctesting changing to these
c-mc      these are even stricter...
c      IF(EMAX.LT.0.015)  DT = 1.1*DTSAVE
c      IF(EMAX.LT.0.007)  DT = 1.2*DTSAVE
c      IF(EMAX.LT.0.001)  DT = 1.4*DTSAVE
c      IF(EMAX.LT.0.0008)  DT = 2.0*DTSAVE
c      IF(EMAX.LT.0.0004)  DT = 3.0*DTSAVE
c      IF(EMAX.LT.0.0001) DT = 4.0*DTSAVE
c     IF(EMAX.LT.0.00005) DT = 5.*DTSAVE


      DTINV = 1./DT
      ZMAX = Z(JS)
C      print 373, n, TIME, DT,EMAX,ISPEC(IS),ZMAX, USOL(is,js),
C     $ USOL(LO2,1), USOL(LCH4,1), USOL(LH2,1), USOL(LCO,1)
C 373  format (1x, I6, 1P3E12.3,2x,A8,1P1E12.3,4x,1P5E12.3)
C
C   skip oxidation state and sulfur budget
      IF (SM-MS.GT.0.01) GOTO 317

C following section run every three timesteps and printed to out.out
C Oxidation state stuff commented out for now.

      write(14, 100) N,EMAX,ISPEC(IS),ZMAX,UMAX,RMAX,DT,TIME
 100  FORMAT(1X,'N =',I5,2X,'EMAX =',1PE9.2,' FOR ',A8,
     2  'AT Z =',E9.2,1X,'U =',E9.2,1X,'RHS =',E9.2,
     3  2X,'DT =',E9.2,2X,'TIME =',E9.2)
C     print this to terminal
       print 100, N, EMAX, ISPEC(IS), ZMAX, UMAX, RMAX, DT, TIME

C
C   COMPUTE ATMOSPHERIC OXIDATION STATE
c   what follows needs work - WARNING
      DO 42 I=1,NQ
      SR(I) = 0.
      DO 43 J=1,JTROP
  43  SR(I) = SR(I) + RAINGC(I,J)*USOL(I,J)*DEN(J)*DZ(J)
      PHIDEP(I) = VDEP(I)*USOL(I,1)*DEN(1)
  42  TLOSS(I) = SR(I) + PHIDEP(I)


c nb that vdep for particles is defined to include wfall
C    when particles are in the main loop

      if(USETD.EQ.1.AND.NP.GT.0) then
C This mimics the code Jim has, but is more general.
C I don't think I ever got this working 100% correctly to where
C  I could balance the sulfur budget when using the tri-diag WARNGING

CTHIS STILL NEEDS WORK
C need to fill up rainout and depostion vectors for
C the triadiagonal species to make the budgets work out.
        do i=NQ+1,NQ1
          SR(I)=0.
            do j=1,JTROP
C      ACK - all particles raining out like H2SO4 
              SR(I) = SR(I)+ RAINGC(LH2SO4,J)*PARTICLES(J,i-nq)
     $                *DEN(J)*DZ(J)
            enddo
C-ACK hardcoded turbulent diffusion velocity
            PHIDEP(I)=(WFALL(1,i-nq)+vturb)* PARTICLES(1,i-nq)*DEN(1)
            TLOSS(I) = SR(I) + PHIDEP(I)
C     in general SR>>PHIDEP for particles,by about 100X
c       print *, ISPEC(I),SR(I),PHIDEP(I),SR(I)+PHIDEP(I)
        enddo

c      stop
C  End of the tri-diag budgeting loop
      endif

C
c the following is obsolete I hope WARNING
c      PHIESC = 2.5E13 * (0.5*USOL(LH,NZ) + USOL(LH2,NZ))
c    2  + USOL(LH2O,NZ) + 2.*USOL(LCH4,NZ) + 3.*USOL(LC2H6,NZ))
C


C - MC check the below with respect to the new boundary conditions.
C  Are some SGFLUXes that are not activly used in the code being counted here?

c - these could be done better in the manner of the redox computation
c - also how will these do if the L number doesn't exist, or if the sl/ll thing.
Ctloss runs from 1-nq1 so cant have any SL's in here...

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

c - note that these aren't really used for anything.
C  I am tempted to kill to allow for greater ease in switching SL and LL.
C  If I want to have a redox printout at every time step,
C  I should find a way to make these generic like I did with redox in Output.


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
c-mc should think about what (if anything) would be actually useful
C-MC here to print out on every timestep WARNING


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
C   valid timestep, so update USOLPREV vector
      ELSE
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

C     NAN HAPPENING SOMEWHERE BELOW HERE...MUST BE IN OUTPUT.F


C
C
C NPR=PRN set above to 50 - i.e. write out every 50 steps
      NS = N/NPR
      SN = N/PRN
      IF(NN.EQ.NSTEPS) SN = NS
      IF (SN-NS .LT. 1.E-3) THEN
        CALL OUTPUT(NN,NSTEPS,TIME,jtrop, vdep,USOLORIG,USETD,frak)
      ENDIF
C    NAN HAPPENING SOMEWHERE ABOVE HERE...MUST BE IN OUTPUT.F WARNING


      write(23, 374) n, TIME, USOL(LO2,1), USOL(LH2,1), USOL(LCO,1),
C   tri-diag (could fix if I wanted...) WARNING
     $ USOL(LCH4,1)!, USOL(LS8aer,1)
 374  format (1x, I6, 1P6E13.4)



c-mc writing out full number densities at each timestep
C this eventually should be an option as it is
C a large output file (1MB per 50 steps)


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
       write(43,115) N,TIME,DT,EMAX
        do I=1,nz
c Not needed in whiff testing   WARNING
C         write(43,114), (USOL(K,I)*DEN(I),K=1,NQ)
c         write(43,114), (USOL(K,I)*DEN(I),K=1,NQ)
        enddo
      endif

c      print *, RPAR(1,3),RPAR(NZ,3)

C ACK hardcoded NQ - update if NQ>100
c 114  format(100(1pe10.3))
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

C Successful completion
  22  CONTINUE


C write out formatted CONTINUE
C-PK Write to file used for spectrum (VPL-SMART)
        BOLTZ = 1.38E-16 !cgs
       do J=1,NZ
        PRES_bar(J) = DEN(J)*BOLTZ*T(J)*1e-6
       enddo
      WRITE(159,938)
      WRITE(159,937) (Z(I),T(I),DEN(I),PRES_bar(I),USOL(LH2O,I),
     2 SL(LCH4,I)/DEN(I),SL(LC2H6,I)/DEN(I),SL(LCO2,I)/DEN(I),
     & SL(LO2,I)/DEN(I),O3(I),
     3 USOL(LCO,I),USOL(LH2CO,I),SL(LHNO3,I)/DEN(I),
     4 USOL(LNO2,I),USOL(LSO2,I),USOL(LN2O,I),I=1,NZ)
 938  FORMAT(1x,'    Alt      Temp       Den      Press      H2O ',
     2        '      CH4      C2H6       CO2      O2         O3  ',
     3        '      CO       H2CO       HNO3     NO2        SO2',
     4        '      N2O')

      WRITE(67,738)
      WRITE(67,941) (Z(I),T(I),DEN(I),PRES_bar(I),USOL(LH2O,I),
     2 SL(LCH4,I)/DEN(I),SL(LC2H6,I)/DEN(I),SL(LCO2,I)/DEN(I),
     & SL(LO2,I)/DEN(I),O3(I),
     3 USOL(LCO,I),USOL(LH2CO,I),SL(LHNO3,I)/DEN(I),
     4 USOL(LNO2,I),USOL(LSO2,I),USOL(LN2O,I),I=1,NZ)
 738  FORMAT(1x,'    Alt      Temp       Den      Press      H2O ',
     2        '      CH4      C2H6       CO2      O2         O3  ',
     3        '      CO       H2CO       HNO3     NO2        SO2',
     4        '      N2O')


 937  FORMAT(1X,1P16E10.3)
 941  FORMAT(1X,1P16E10.3)
      WRITE(164,939)
      WRITE(164,940) (Z(I),T(I),DEN(I),PRES_bar(I),
     & SL(LCO2,I)/DEN(I),SL(LCO,I)/DEN(I),SL(LC2H6,I)/DEN(I),
     & SL(LH2CO,I)/DEN(I),
     & USOL(LCH4,I),USOL(LHNO3,I),USOL(LNO2,I),USOL(LO2,I),O3(I),
     & USOL(LSO2,I),USOL(LH2O,I),USOL(LN2O,I),I=1,NZ)
 939  FORMAT(1X,'    Alt      Temp       Den      Press      CO2 ',
     2       '       CO       C2H6       H2CO     CH4        HNO3',
     3       '       NO2      O2         O3       SO2        H2O ',
     4        '      N2O')
 940  FORMAT(1X,1P16E10.3)




C write out formatted out.dist file
c this new format works automatically even if NQ changes

C   Number columns in .dist file
      IROW = 10
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
C    Print USOLS into another file .chem
        write(70, 880) ((USOL(k,i),K=K1,K2),i=1,nz)
       ELSE
          K2 = NQ
C    this only occurs if NQ is a multiple of 10
          if (K2-K1 .EQ. 9) then
            write(18, 880) ((USOL(k,i),K=K1,K2),i=1,nz)
            write(70, 880) ((USOL(k,i),K=K1,K2),i=1,nz)
C    if not, need to generate a dynamic format statement
          else
           fmtstr='(  E17.8)'
C    OK for one character as this should always be <10
           write(fmtstr(2:3),'(I1)')K2-K1+1
           write(18, fmtstr) ((USOL(k,i),K=K1,K2),i=1,nz)
           write(70, fmtstr) ((USOL(k,i),K=K1,K2),i=1,nz)
          endif
       ENDIF
      enddo

 880  format(10E17.8)
 881  format(5E17.8)
C     (NSP - 1), indicates (inert) CO2 density for terrestrial planets ONLY
        write (18,881) (T(i),EDD(i),DEN(i),O3(i), SL(LCO2,i),i=1,nz)
C     Print into another file .strctr
        write (66,881) (T(i),EDD(i),DEN(i),O3(i), SL(LCO2,i),i=1,nz)

      IF (NP.GT.0) THEN
        fmtstr='(  E17.8)'
        write(fmtstr(2:3),'(I2)')NP*3
        do i=1,nz
         write(18,fmtstr) (AERSOL(i,j),j=1,NP),(WFALL(i,j),j=1,NP),
     $                  (RPAR(i,j),j=1,NP)
         write(71,fmtstr) (AERSOL(i,j),j=1,NP),(WFALL(i,j),j=1,NP),
C      Print aerosols into another file .aersol
     $                  (RPAR(i,j),j=1,NP)
        enddo



        fmtstr='(  E17.8)'
        write(fmtstr(2:3),'(I2)')NP

       if(USETD.EQ.1) then
        do i=1,nz
C   Ordering is that in species.dat
         write(18,fmtstr)  (PARTICLES(i,j),j=1,np)
C   Print tridag into another file .tridag
         write(72,fmtstr)  (PARTICLES(i,j),j=1,np)
        enddo
       endif
      ENDIF

c New abstracted photorates printout

C Number columns of photorates
      IROW = 8
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
C       +2 here so it writes the Z as well
           write(fmtstr2(2:3),'(I2)')K2-K1+2
           write(14, fmtstr2) (Z(I),(prates(k,i),K=K1,K2),i=1,nz,3)
       ENDIF
      enddo

       fmtstr="(  A12)"
C  format for photorate labels, with extra space will crash if kj>99
       write(fmtstr(2:3),'(I2)')kj


c-mc  write out important parameters:
      write(49,*) NZ,NQ,NQ1,NSP2,NR,KJ,NP
      write(49,fmtstr) photolabel
      write(49, *) ISPEC

c-mc write out photorates at all heights in the .so2 file

       JSO2=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'SO2    ')
       JO2=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'O2     ')+1

C note this is for the O2 + Hv -> O + O reaction,which is the second O2 reaction
C JH2O and JCO2 defined above...
      write(27,299) (Z(I),prates(JSO2,I),prates(JO2,I),prates(JH2O,I),
     $     prates(JCO2,I),I=1,NZ)
 299  FORMAT(5(1PE9.2,2x))
c-mc

c-mc write out all photorates at all heights in the .rates file
      string='(  (1PE9.2,2X))'
C   ack - hardcoded kj+1 here. will break if kj>99
      write(string(2:3),'(I2)')KJ+1

C The above is a way to dynamically create a format string at runtime.
C   replaces:  599  FORMAT(55(1PE9.2,2x))
       write(28,string) (Z(I),(prates(J,I),J=1,KJ),I=1,NZ)
c-mc

C Write out O2 photolysis rates, for testing the S-R bands using various methods
       do i=1,nz
C ack - hardcoded O1D rate number (OK as long as O+O follows O+O1D, which is OK)
       write(58,*) Z(I),prates(JO2-1,I),prates(JO2,I)
       enddo

C print out rainout rates
       IF (NAQ.GT.0) THEN
       do i=1,nq
           write(59, *) (RAINGC(i,j),j=1,jtrop)
       enddo
       ENDIF

       JNO=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'NO    ')

C  If using old grid, write out NOprates so they can be used for new grid
       if(LGRID.EQ.0) write(60,*) (prates(JNO,I),I=1,NZ)

c      do i=1,nsp
c         print *, i, ISPEC(i)
c      enddo


c       if (ICOUPLE.eq.1) then
c GNA - I think it should be writing this file regardless
c if e.g. you run photo uncoupled to get initial conditions
c and then want to feed those into clima, you need to have
c this file written even if ICOUPLE = 0
C Transfer results to the climate model (COUPLING)
        DO 255 I=1,NZ
C Transfer O3 and H2O
         WRITE(84,254) Z(I),PRESS(I),USOL(LO3,I),USOL(LH2O,I),
     &                 max(USOL(LCH4,I),1.e-60),SL(LCO2,I)/DEN(I),
C  EWS debug to prevent floating point errors WARNING
     &                 max(USOL(LC2H6,I),1.e-60)
  254    FORMAT(1PE9.3,6(E10.2))
  255   CONTINUE
        close(84)
c       endif

C      Because argon is not normally calculated in photochem
         IF (IRESET.eq.0) FAR=1.000E-02
C-EWS  debug to prevent floating point errors
         FCH4=max(USOL(LCH4,1),1.e-60)
C-EWS  debug to prevent floating point errors
         FC2H6=max(USOL(LC2H6,1),1.e-60)
         FCO2=SL(LCO2,1)/DEN(1)
         FN2=SL(LN2,1)/DEN(1) + FCO2
C-EWS Note that CLIMA/mixing_ratios.dat treats the condensible (CO2) and
C     non-consibles mixing ratios differently
C     e.g., if the atmosphere is 99% N2 and 1% CO2, then the CO2 fraction is
C     0.01 and N2 should be set to 1, because it is 100% of noncondensibles.
C     In practice, N2 should be = (1 - [everything but CO2]).
         FO2=USOL(LO2,1)
C-mab: Using presence of HE in template to decide if H2 should be reset as below.
         IF(FHE.LE.0.0)FH2=USOL(LH2,1)
C-mab: Since H2 is major for giant planets, we do (1 - everything) for them.
C     But in other parts of photochem, N2 is of the total,
C     not excluding CO2, so we only change it here. 9/8/2015
C-gna clima can't currently cope with NO2 and having it is screwing it up
C WARNING
         FNO2=USOL(LNO2,1)*1.0e-60
         JCOLD=JTROP

         WRITE(117,102) FAR, FCH4, FC2H6, FCO2,
     &                  FN2, FO2, FH2, FNO2, JCOLD
  102    FORMAT(1PE10.3,10x,'!Argon'/,1PE10.3,10x,'!Methane'/,
     &          1PE10.3,10x,'!Ethane'/,1PE10.3,10x,'!Carbon Dioxide'/,
     &          1PE10.3,10x,'!Nitrogen'/,1PE10.3,10x,'!Oxygen'/,
     &          1PE10.3,10x,'!Hydrogen'/,1PE10.3,10x,
     &          '!Nitrogen Dioxide'/,I3,17x,'!Tropopause layer')

C MOVING HIGHER UP Reading in the temperature and water profiles
C from the climate code
C
C      IF (ICOUPLE.eq.1) then
C        print *, 'alt_new, T_new, water'
C        DO J=1, NZ
C
C Reading the altitude, temperature, and water profiles from the climate code
C           READ(116,*) alt_new(J), T_new(J), water(J)
C           print 351, alt_new(J), T_new(J), water(J)
C        END DO
C      close(116)
C      endif
C 351  FORMAT (1PE10.3, 1PE12.3, 1PE12.3)


      print*,"PhotoMain run completed...."
      print*,"See PHOTOCHEM/OUTPUT/ directory for the output files."

        STOP
C    error in reactions
  20  CONTINUE
        PRINT 300,I
 300    FORMAT(//1X,'NMAX EXCEEDED FOR SPECIES ',I3)
        STOP
C    error in reactions
  25  CONTINUE
        PRINT 301,IERR
 301    FORMAT(//1X,'ERROR IN REACTION ',I3)
        STOP
      END

      SUBROUTINE PROMPT(CHAR)
c-mab: PROMPT: sends a prompt input (CHAR) to terminal from this routine.
c-mab: Use this to be able to pass inputs mid-routine.

      INTEGER L
      CHARACTER*(*) CHAR
      DO 100 L=LEN(CHAR),1,-1
      IF(CHAR(L:L).NE.' ')THEN
        IF(CHAR(L:L).EQ.':')THEN
          WRITE(*,12) CHAR(1:L)
 12       format(' ',a,$)
         ELSE
          WRITE(*,13) CHAR(1:L)
 13       format(' ',a,' :',$)
          END IF
        RETURN
       END IF
 100  CONTINUE
      RETURN
      END
