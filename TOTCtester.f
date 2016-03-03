            PROGRAM TOTCtester
c
ccc test test
c I am attempting to abtract this so we can have an Earth/Mars switch
c I should also think about adding extra abstraction for the humidty stuff
c
c - this code is undergoing active development for chlorine speciesS
c - as such it is not really the best code to be used as a base for TOTC, but such it is.
c - it is starting from one of the more up-to-date branches of Mark's code, but doesn't have any time-dependent gear stuff in it

c- at some point go through and clean up all comments
c
c
c - this code contains the variable grid e changes used to compute the suite of models for the whiff paper.
c - this code could/should be modified to use a variable grid size at some point to make it faster
c - all common blocks abstracted to DATA/INCLUDE

c - re-introducing the tridiagonal...
c - and attempting to take it away...

c - g77 will no longer work as there are some F90 built-ins being used
c see makefile for compilation syntax

cc

c-mc
c this code has iterated jacobian and batch multipliers
c-072406 - has 195K CO2 cross-section ! turned off for photo development

c-mc 080206 - iterated jacobion, new co2 x-section are turned off, TSTOP is 1e16
c-mc 112107 - iterated jacobians back in play for reruns of kevin's models for the whiff paper
c-mc 111308 - iterated jacbians off for chlorine testing...
c-mc 010609 - iterated jacobians back on during chlorine testing...

c-mc  THIS VERSION OF THE CODE WILL BE TIME-INDEPENDENT (i.e. P,T,zenith angle remain constant) (or maybe not...)

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
c  but for giggles I am doing an experiment wiht big O2 and CH4 fluxes to see what happens
c  this is hard to do if one uses O2 flux; it works better by setting the O2 mixing ratio

c   some observations
c      1. problems with NOx at the upper boundary often means a problem with H escape
c          NOx problems often resolve themselves with patience...
c      2.  getting lost often means a problem with H escape
c      3.  problems with S also means that there is a problem with redox balance.


c  update 6-23-05 inserts molecular diffusion of H
c  looks good

c   I added molecular diffusion for H2 and H.
c
c   the more general equation would include the molecular diffusive flux
c         \phi = b*f*({1/over H_A}-{1/over H}) - b*df/dz
c
c     [ w_1 - w_2 = {b_{12}\over n} \left( (1/n_1)(dn_1/dz) - (1/n_2)(dn_2/dz) - (1/f)(df/dz) \right)  ]
c  (I've not been careful. here f=n_1/n_2 <<1, so n=n_2.  I should look this up for the
c    general case)
c
c    where H would be the scale height in vacuum of the species corresponding to f
c     and H_A is the scale height of the atmosphere.

c         \phi = b*f*({1/over H_A}-{1/over H}) - (b+KN)*df/dz
c          N df/dt = P - LNf - d\phi/dz - f*dN/dt
c    I have not done the work to allow for N to change with time  - so let dN/dt=0
c
c    From this point hence one must be careful
c
C         df/dt  =  P/N - L*f + (1/N)*d/dz[(K*N+b)*(df/dz)] + (b/N)*d/dz({1/over H}-{1/over H_A})f
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
C     df(0)/dt  =  P(0)/N(0) + DU(0)*(f(1)-f(0))  + ADU(0)*(f(1)+f(0)) + FLUX(0)/N(0)/dz
c
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
c        djac(1) =  0.0   (off diagonal - is this correct??  - it may not be, but it should not matter given that it can't change f)
c
c  at the upper boundary
c
C     df(nz)/dt  =  P(nz)/N(nz) + DL(nz)*(f(nz)-f(nz-1))  + ADL(nz)*(f(nz)+f(nz-1)) - FLUX(nz)/N(nz)/dz
c
c     where
c        DL(nz) = [KN(nz-0.5) + b(nz-0.5)]/[N(nz)*dz*dz]
c       ADL(nz) = +b(nz-0.5)/[2N(nz)*dz] *({1/over H} - {1/over H_A})

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
c   to zero at the top (no flux upper boundary).  This does a better job with the chemistry.
c   i've left this intact for ch4co2.for
c

c  4-25-05  jfk suggests distributing SO2 source over the troposphere
c           this way i can use a surface BC
c    This was implemented successfully

c  4-27-05  FIXED the accounting errors with fluxes.
c  the tridiagonal was not being included in the accounting
c  I wasn't sure how to do this save by putting S8aer into long lived species
c  This has now been done.  Fluxes are now under control.
c  redox balance and S conservation are both enforced.  Fundamentally nothing has really
c   changed, but it has the advantage of being done right.

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
c     H       OCS     CO      HS              9.10E-12  0.00   1940.0   0.0      0.0
c     HS      CO      OCS     H               5.95E-14  1.12   8330.0   0.0      0.0
c     O       OCS     CO      SO              7.83E-11  0.00   2620.0   0.0      0.0
c     O       OCS     S       CO2             8.33E-11  0.0   5.53E+03  0.0      0.0
c     OCS     S       CO      S2              1.00E-10  0.0   4.56E+03  0.0      0.0
c     OCS     OH      CO2     HS              1.10E-13  0.0   1.20E+03  0.0      0.0
c     S       HCO     OCS     H               6.00E-11  0.0     0.0     0.0      0.0
c     OCS     HV      CO      S
c and...
c     S        CO     OCS                     1.70E-33  0.0   1.51E+03  1.0      0.0
c  I'm guessing that S + CO goes at the same rate as O + CO

c  the current version has the fictional
c     SO3     CO      CO2     SO2             1.00E-16  0.0     0.0     0.0      0.0


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
c - keeping around for the time when more research is done on the chemical species
C
C       )  HSO + NO2 -> NO + HSO2   ! 9.6e-12.  HSO2 is not a species in this system - it reacts quickly with O2
c       )  H + SO2 + M -> HSO2 + M  !  this has measured but very slow reaction
c       )  HSO2 + O2 -> HO2 + SO2   ! 3e-13 measured HO2 as product.  not a clear situation;
c                                     my temptation is to assign OH, H rates measured at high T
c       )  HSO2 + H -> H2 + SO2     ! 3e-12
c       )  HSO2 + OH -> H2O + SO2   ! 8e-12
c       )  HS + NO + M -> HSNO + M   ! JPL has this!  but what do you do with it??
c       )  HCO + HNO  -> H2CO + NO   ! NIST 2e-12*exp(-1000/T)
C       )  CH3 + OH + M ->  CH3OH + M  ! NIST - v. fast: 2.5e-27*(298/T)^-3.8, 1.5e-10 limit
C       )  CH3 + OH  ->  H2O + CH2    ! NIST 1.2e-10 *exp(-1400/T)
C       )  CH3 + OH  ->  HCOH + H2    ! NIST 9e-11 *exp(-1500/T)
C       )  CH3 + OH  ->  H2COH + H    ! NIST 1.3e-11
c       )  S2  + H  + M ->  HS2 + M
c       )  S2  + OH  ->               ! SO +OH -> H + SO2 8.6e-11
c       )  S3  + OH  ->
c       )  S3  + O  ->  SO + S2
c       )  S3  + H  + M ->  HS3 + M
c       )  S4  + OH  ->
c       )  S4  + O  -> SO + S3
c       )  S4  + H  + M -> HS4 + M
c       )  HS2 + OH -> H2O + S2
c       )  HS2 + O  -> OH + S2
c       )  HS2 + H  -> H2 + S2
c       )  HS2 + S  -> H + S3
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
!c-mc  USOLSAVE added for convergence testing 4/29/06, removed from code on 2/28/07
!c-mc  going to keep declaration here and reserve it in case I ever bring back the
!c-mc  the convergence-testing code, which is r37:36 in the /td branch
c-mc  USOLPREV(NQ,NZ) added for second order reverse Euler calculations
      DIMENSION DPU(NZ,NP),DPL(NZ,NP)
      DIMENSION TA(NZ),TB(NZ),TC(NZ),TY(NZ)
      dimension PRES_bar(NZ)

c      dimension atomsO(NSP2),atomsH(NSP2),atomsC(NSP2)
c      dimension atomsN(NSP2),atomsCL(NSP2),atomsS(NSP2)

!temp
      dimension testvec(NR),fixedmr(nq)
      dimension distheight(nq)

      CHARACTER*11 photolabel, AA



C ISOHACK
c first one is for standard model
      parameter(NISO=15,NQI=NQ+NISO,NPISO=0,NISOSL=7)  !ISOHACK  - OK this is being hardcoded to also include the HCAER species with the particles
c this is for biogenic sulfur...
c      parameter(NISO=19,NQI=NQ+NISO,NPISO=0,NISOSL=7)  !ISOHACK  - OK this is being hardcoded to also include the HCAER species with the particles

      parameter(NPN=NISO+NPISO+NISOSL)           !number of sulfur species in the model  !ISOHACK
c      dimension USOLISO(NSP-NPN-2,NZ)  !OK - this was used in the initial testing of the ISOTOPE code (where S* = S)
      dimension USOLISO(NSP,NZ)

c for sgbco testing (to get at Jacobian condition number)
c      real*8 RCOND
c      real*8 WORK(neq)


      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PHOTABLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/BBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/CBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/DBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/EBLOK.inc'   !can go away when MSCAT does
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
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/PBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/ISOBLOK.inc'
      INCLUDE 'PHOTOCHEM/DATA/INCLUDE/comPRESS1.inc'


C   CONSTANTS FOR 1ST EXPONENTIAL INTEGRAL  !can go away when MSCAT does
      DATA AI/.99999193, -.24991055, .05519968, -.00976004,
     2  .00107857, -.57721566/
      DATA BI/8.5733287401, 18.0590169730, 8.6347608925,
     2  .2677737343/
      DATA CI/9.5733223454, 25.6329561486, 21.0996530827,
     2  3.9584969228/
      DATA NUML,NUMP/NSP*0,NSP*0/

      ISOTOPE=0   !Isotope flag
      ISOS=32     !use 32S (i.e. the dominant isotope)

c OPEN FILES



      open(2, file='PHOTOCHEM/DATA/aerosol.table',status='OLD') !,form='UNFORMATTED')
      open(3, file='PHOTOCHEM/DATA/photo.dat',status='OLD')
      open(4, file='PHOTOCHEM/INPUTFILES/species.dat',status='OLD')
      !planet parameters (G, FSCALE, ALB,ZTROP,etc)
      open(7, file='PHOTOCHEM/INPUTFILES/PLANET.dat',status='OLD')
      open(231, file='PHOTOCHEM/INPUTFILES/input_photchem.dat',
     &          status='OLD')       !model parameters (AGL, IO2,INO, LGRID, etc)
      ! reaction file
      open(9, file='PHOTOCHEM/INPUTFILES/reactions.rx',status='OLD')
      open(14, file='PHOTOCHEM/out.out',status='UNKNOWN')    ! output
      open(17, file='PHOTOCHEM/in.dist',status='OLD')         ! formatted input
      open(18, file='PHOTOCHEM/out.dist',status='UNKNOWN')    ! formatted output
      open(19, file='PHOTOCHEM/out.terse',status='UNKNOWN')    ! terse output
      open(67, file='PHOTOCHEM/profile.pt',status='UNKNOWN')   !mixing ratios
      open(68, file='PHOTOCHEM/hcaer.out',status='UNKNOWN')   !hydrocarbons
      open(69, file='PHOTOCHEM/hcaer2.out',status='UNKNOWN')   !other hydrocarbons
      open(21, file='PHOTOCHEM/out.trs',status='UNKNOWN')    ! very terse output
      open(23, file='PHOTOCHEM/out.time',status='UNKNOWN')    ! time output
      open(24, file='PHOTOCHEM/out.tim',status='UNKNOWN')    ! time output

c The next four files are used only when this model is coupled with the climate model (ICOUPLE=1)
      open(90, file='COUPLE/hcaer.photoout.out',status='UNKNOWN')   !hcaer for clima - GNA
      open(84, file='COUPLE/fromPhoto2Clima.dat',
     &         status='OLD')          ! To be used as input for the climate model (coupling)
      open(116, file='COUPLE/fromClima2Photo.dat')
      open(117, file='COUPLE/mixing_ratios.dat')
      open(118, file='COUPLE/fromClima2Photo_works.dat')

c-mc
      open(25, file='PHOTOCHEM/out.redox',status='UNKNOWN')    ! redox output - eventually combine into out.trs
      open(26, file='PHOTOCHEM/out.converge',status='UNKNOWN')  ! temporary output file for looking at convergence
                                                     ! in the reverse Euler iteration
      open(27, file='PHOTOCHEM/out.so2',status='UNKNOWN')      !printing out relevant SO2 photolysis pieces
                                               !between 190-220 as MIF signature
      open(28, file='PHOTOCHEM/out.rates',status='UNKNOWN')  !reaction rates,reactions,species,densities,rate constants


      open(29, file='PHOTOCHEM/out.xsec',status='UNKNOWN')  !cross sections
      open(30, file='PHOTOCHEM/out.gridz',status='UNKNOWN')  !height grid
      open(31, file='PHOTOCHEM/out.gridw',status='UNKNOWN')  !wavelength grid

c-mc  unit 33 reserved for wavelength grids...

      open(41, file='PHOTOCHEM/out.rad',status='UNKNOWN')
      open(42, file='PHOTOCHEM/out.finalden',status='UNKNOWN') !this file should contain TATLTUAE
      open(43, file='PHOTOCHEM/out.densities',status='UNKNOWN') !number densities at each timestep
c not creating this file for whiff testing, but needs to be in real td versions...
      open(44, file='PHOTOCHEM/out.prod',status='UNKNOWN') !total production and loss at steady state
      open(45, file='PHOTOCHEM/out.flux',status='UNKNOWN') !fluxes
      open(48, file='PHOTOCHEM/out.tau',status='UNKNOWN') !tau=1 (at final step)
      open(49, file='PHOTOCHEM/out.params',status='UNKNOWN') !some model params
      open(50, file='PHOTOCHEM/out.error',status='UNKNOWN') !NGE and L2 norm between start and finish
      open(51, file='PHOTOCHEM/out.cl',status='UNKNOWN') !TP/FLOW for chlorine species,nitrate, adn sulfate


      open(52, file='PHOTOCHEM/ISOin.dist',status='UNKNOWN')    ! formatted output - ISOHACK - for ISO model...
      open(53, file='PHOTOCHEM/ISOinert.dist',status='UNKNOWN')    ! formatted output - ISOHACK - for ISO model...

      open(58, file='PHOTOCHEM/out.O2prates',status='UNKNOWN')  !for testing O2 prates with various grid sizes
      open(59, file='PHOTOCHEM/out.raingc',status='UNKNOWN')  !rainggc out !ISOHACK

!c-mc 60 and 61 are opened below after LGRID is read in

       open(62, file='PHOTOCHEM/out.flow',status='UNKNOWN') !lower boundary fluxes (not including rainout)
       open(63, file='PHOTOCHEM/out.od',status='UNKNOWN')   !haze optical depths
       open(73, file='PHOTOCHEM/out.rp',status='UNKNOWN')   !haze toa and sur rpar
 

C - This file gives the input needed for SMART - radiative transfer code

       open(159, file='photochem_smart/photchem/photfile.pt',
     &          status='UNKNOWN')
       open(164, file='photochem_smart/atm/tag.atm',
     &          status='UNKNOWN')

C - Seperating the out.dist file into seperate files

       open(70, file='PHOTOCHEM/out.chem', status='UNKNOWN')
       open(66, file='PHOTOCHEM/out.strctr', status='UNKNOWN')
       open(71, file='PHOTOCHEM/out.aersol', status='UNKNOWN')
       open(72, file='PHOTOCHEM/out.tridag', status='UNKNOWN')


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
      READ(231,*)AA, NEWSPEC
      IF(IDEBUG.eq.1) print *, "NEWSPEC =",NEWSPEC
      READ(231,*)AA, ihztype
      IF(IDEBUG.eq.1) print *, "IHZTYPE =",ihztype
      READ(231,*)AA, ZY
      IF(IDEBUG.eq.1) print *, "ZY =",ZY
 555  format(3/)
      close(231)

      if (LGRID.EQ.0) open(60, file='PHOTOCHEM/out.NOprates',
     &                         status='UNKNOWN')  !NO photolysis rates output
      if (LGRID.EQ.1) open(61, file='PHOTOCHEM/out.so2HR',
     &                         status='UNKNOWN') !wavelength specific so2 photorates on HR grid

C - READ IN SPECIES NAMES, ATOMIC NUMBERS, AND BOUNDARY CONDITIONS

      iLL=0   !counter for long lived species
      iSL=0   !counter for short lived species
      iTD=0   ! counter for tridiagonal species
      iIN=0   !counter for inert species
      iSP=0   !counter for number of lines in species.dat file
      iprint = 0                !just so the species.dat statements just print once

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



               if (NEWSPEC.eq.1) then !new species.dat formatting
                  if (iprint.eq.0) then
                     print *, 'species.dat should have new formatting'
                     print *, "for VDEP and FIXEDMR (E8.2)"
                     iprint = 1
                  endif
                  read(4,208) LBC, XX,YY,ZZ,XXX,LG,YYY,ZZZ !read in boundary conditions
               endif

               if (NEWSPEC.eq.0) then !old species.dat formatting
                   if (iprint.eq.0) then
                     print *, 'species.dat should have old formatting'
                     print *, "for VDEP and FIXEDMR (E7.2)"
                     iprint = 1
                  endif
                  read(4,210) LBC, XX,YY,ZZ,XXX,LG,YYY,ZZZ !read in boundary conditions
               endif
               print *, LBC
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

               if (species.EQ.'CO2') FCO2=YY !hardcoding woohoo! CO2 only works as fixed mixing ratio. This could be handled better...

            endif

            if (SPECTYPE.EQ.'IN') then
               iIN=iIN+1
               backspace 4  !return to previous line in species.dat file
               read(4,209) XX  !read in fixed mixing ratios
               if (species.EQ.'CO2') FCO2=XX   !hardcoding woohoo!  !need to do N2 as well
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
 208  format(30X,I1,5X,2(E8.2,1X),E9.3,1X,E7.1,1X,I1,6X,2(E7.1,1X))  !for boundary conditions
 210  format(30X,I1,5X,2(E7.1,1X),E9.3,1X,E7.1,1X,I1,6X,2(E7.1,1X)) !for boundary conditions
c 208  format(30X,I1,5X,2(E8.1),E9.3,1X,E7.1,1X,I1,6X,2(E7.1,1X))  !for boundary conditions
 209  format(30X,E7.1) !for INERT species boundary conditions

 96   CONTINUE

c      stop
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
      ! chemical reaction file
      open(9, file='PHOTOCHEM/INPUTFILES/reactions.rx',status='OLD') 
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


c-mc the below could be made into a nice printout in the out.out file if someone felt like it.
c       print *, SUM(INDEX(REACTYPE,'PHOTO'))
c       print *, SUM(INDEX(REACTYPE,'2BODY'))
c       print *, SUM(INDEX(REACTYPE,'3BODY'))
c       print *, SUM(INDEX(REACTYPE,'WEIRD'))
c       print *, NR


C ***** READ THE PLANET PARAMETER DATAFILE *****
      READ(7,502) G,FSCALE,ALB,ZTROP,FAR,R0,P0,PLANET,TIMEGA,IRESET,
     &   msun, ihzscale
 502  FORMAT(F7.1/,F7.2/,F7.3/,E7.1/,F7.3/,E8.3/,F8.3/,A8/,F4.2/,I1/
     &  ,I2/,I1)
C     adding IRESET to create the atmospheric profile for modern earth
C     gna - added msun keyword to change the star from planet.dat
      print *, 'msun is ', msun
      print *, 'ihzscale is ', ihzscale
C gna-scale stellar flux based on spectral type:
C uses Kopparapu et al 2012 scalings for earth-equivalent distance
      IF (ihzscale.eq.1) then
         IF (msun.eq.16) FSCALE = FSCALE * 0.870
         IF (msun.eq.17) FSCALE = FSCALE * 0.859
         IF (msun.eq.18) FSCALE = FSCALE * 0.950
         IF (msun.eq.19) FSCALE = FSCALE * 1.110
         IF (msun.eq.76) FSCALE = FSCALE * 0.866
         print *, 'fscale is ', fscale
      ENDIF



C
      IF (IRESET.eq.1) CALL RESET
C Defines mixing ratio values for modern earth; don't need next section that pulls ancient earth mixing ratios from in.dist

      IF (IRESET.eq.0) then
C ***** READ THE INPUT DATAFILE *****

c read in formatted input data file


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

      !EWS - there was an error in this loop fixed 8/12/2015
      DO K=1,LR
         DO I=1,NZ
           IF(USOL(K,I).lt.1.e-30) USOL(K,I)=1.e-30  !EWS - fixed an error here where USOL=1.e-30 instead of USOL(K,I)=1.e-30
         END DO
      END DO

      read(17,881) (T(i),EDD(i),DEN(i),O3(i), SL(LCO2,i),i=1,nz) !this final one is CO2 number density
!SL(NSP-1) will be the density of the final short-lived species if CO2 is removed from the list.  I dont THINK it will matter now that I call DOCHEM with the -1 before starting....
! Changed NSP-1 to LCO2 to avoid hard coding

!gna - we need to make it so that T = T_new
      IF(ICOUPLE.EQ.1) THEN
         DO I=1, NZ
            IF(T_new(I).gt.100) THEN !gna - avoid some crazy unconverged clima solutions that make photo not converge -- better ideas to check for this?
               IF(T_new(I).lt.400) THEN
               T(I) = T_new(I)
               ENDIF
            ENDIF
         END DO
         print *, 'T0 is'
         print *, T(1)
      ENDIF


 
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
        enddo
      endif


C added by giada
      OPEN(unit=999, file='COUPLE/time_frak_photo.out')
 909  FORMAT(1X, F4.2, 7X, F8.3, 5X, F18.16, 5X, I2, 5X, I2, 
     &     9X, I4, 6X, F4.2)
 908  FORMAT(1X, 'timega', 5X, 'P0', 10X, 'frak', 18X, 'msun', 4X, 
     &   'ihztype', 6X, 'NZ', 6X, 'FSCALE')
      print *, frak
      print *, P0
      print *, msun
      print *, ihztype
      print *, NZ     
      print *, FSCALE
      WRITE(999,908)
      WRITE(999,909) timega, P0, frak, msun, ihztype, NZ, FSCALE

C
C
c-mc
c-mc  this block reads in a parameter file which, for the moment, is just
c-mc  a multulplier for the o2 or ch4 concentration
c-mc

       open(10, file='PHOTOCHEM/INPUTFILES/params.dat',status='OLD')
       read(10,303) o2mult, ch4mult
 303   format(F8.1,1x,F8.1)


c       do i=1,nq
c          print *, i, ISPEC(i)
c       enddo
c       stop

c make any changes
       poop3 = 1./1.!*10.**1.94 ! use this to diddle pressure  (from KZ mars code)

c       print *,USOL(LO2,1),2**ch4mult,2.0**ch4mult,
c     $       USOL(LO2,1)/(2.0**ch4mult)

c        USOL(LO2,1) = USOL(LO2,1)/(2.0**ch4mult)

c        print *, USOL(LO2,1)
c        stop


      do I=1,NZ
c-mc fix this stuff when all done (or maybe think of a better way to perturb)

!commenting these out for the HCL run

!        USOL(LCH4,I) = USOL(LCH4,I)*o2mult
c        USOL(LO2,I) = USOL(LO2,I)/2**ch4mult


c - fix this for new lbc thing...
c      HCLFLUX=1.e8*o2mult


        Den(i) = 1.0*Den(i)
      enddo

C ***** SET MODEL PARAMETERS *****
C     ZY = SOLAR ZENITH ANGLE (IN DEGREES)
C     LTIMES = COUNTER FOR PHOTORATE SUBROUTINE
C     DT = INITIAL TIME STEP
C     TSTOP = TIME AT WHICH CALCULATION IS TO STOP
C     NSTEPS = NUMBER OF TIME STEPS TO RUN (IF TSTOP IS NOT REACHED)
C     FO2 = ground level O2 mixing ratio used in heritage calculations

      LTIMES = 0
C      ZY = 50. why was this ever hardcoded? :(  In input_photochem.dat now
      FO2 = USOL(LO2,1)  ! fill up a heritage constant, eventually this should be purged.


      CALL PHOTGRID  !sets up vertical grid

      JTROP=minloc(Z,1, Z .ge. ztrop)-1 !height index for the tropopause (-1 given the staggered grid)
                                        !the 1 in the second postion tells minloc to return a scalar


c       print *, 'enter surface T'
c       read(*,'(F3.0)') TINC
c       print *, TINC

c diddle Tempterature
      do i=1,NZ
       if (i.le. JTROP) then
          zkm=Z(i)/1.e5
c         T(i) = 273. - 6.5*i - 0.25*i*i  ! in Jim's original model
c         T(i) = 288. - 6.5*i - 0.25*i*i  ! what we are used for the sulfur paper
c         T(i) = 293. - 6.5*i - 0.25*i*i  ! what we are using for now - 286K at 0.5 km. - i dont think this was the intended effect

c         T(i) = 290. - 6.5*zkm - 0.25*zkm**2  !this somewhat mimics what's above - where T was computed on a 1km grid that didn't match up directly with the Z grid at 0.5,1.5, etc.


c          T(i) = Tinc - 6.5*zkm - 0.25*zkm**2   !testing a return to cold Archean
c          print *, 289 - 6.5*zkm - 0.25*zkm**2   !testing a return to cold Archean

c         T(i) = 293. - 6.5*zkm - 0.25*zkm**2  !this somewhat mimics what's above - where T was computed on a 1km grid that didn't match up directly with the Z grid at 0.5,1.5, etc.
c         print *,zkm, T(i)  !ok, this is close to the T profile that we were using in the NZ=80 case
!         T(i) = 290. - 6.5*i - 0.25*i*i  ! TD test - uncommenting for now

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

      CALL DENSTY(FO2,poop3)
      CALL RATES
      CALL DIFCO(FO2)  !computes diffusion coefficents (K*N) and binary diffusion coefficents for H and H2
      CALL PHOTSATRAT(JTROP,H2O)
      CALL DOCHEM(FVAL,-1,JTROP,iIN,iSL,USETD)      !IDO=-1, fill up SL for accurate calculation on first timestep

      if (PLANET .EQ. 'EARTH') then
       PRONO = PRONO/1.  ! current column integrated NO production rate on Earth

 
                      ! divide by 1000 turns off ltning
c       PRONO = PRONO/1.e6 ! ATACAMA
      else if (PLANET .EQ. 'MARS') then
       PRONO = PRONO/1.E9 ! divide by 1e9 turns off lightning for dry mars
      endif

      if (mbound(LH2) .gt. 0) then  !i.e if constant mr or constant flux UBC
        do i=1,nz
            bHN2(i) = 0.0   ! don't use molecular diffusion
            bH2N2(i) = 0.0
c           bXN2(i) = 0.0   ! don't use molecular diffusion
        enddo
      else !  use effusion velocity formulation of diffusion limited flux
        Veff(LH) = 1.0*bhN2(nz)/DEN(NZ)
     $     *(1./Hscale(nz) - 1./scale_H(LH,nz))   ! diff lim flux
        Veff(LH2) = 1.0*bH2N2(nz)/DEN(NZ)
     $     *(1./Hscale(nz) - 1./scale_H(LH2,nz))
      endif

!gna - added coupling stuff for water here (just below tropopause)
      do J=1,JTROP
       IF(ICOUPLE.eq.0) THEN
       USOL(LH2O,J) = H2O(J)   !sets H2O to relative humidity in troposphere
       ELSE
       USOL(LH2O,J) = water(J)  !set to h2o from clima if coupling on
       ENDIF
      enddo

       !it's having mega problems for low h2o profiles for cold surface environments
       !here is a fix for now
       IF(T(1).lt.260) THEN
          print *, 'scaling water'

        DO J=1, NZ
        READ(118,*) alt_dontuse(J), T_dontuse(J), water_fix(J)
        END DO
        close(118)
        do J=1,JTROP
        USOL(LH2O,J) = water_fix(J)
        enddo
        endif




      IF (PLANET .eq. 'EARTH') CALL LTNING(FO2)
      CALL AERTAB   !makes table of vapor pressures for H2O and H2SO4
      NZ1 = NZ - 1
      HA = HSCALE(NZ)
      NRAIN = 0   ! count of calls to rainout
c orig      DT = 1.E-15   !why the hell is Kevin starting at the Planck time?
      DT = 1.E-6
      DTINV = 1./DT
      TIME = 0.

c      NSTEPS = 681     !...
c      TSTOP = 1.E14    !for runs that are unstable...


      TSTOP = 1.E17    !as it was...
      NSTEPS = 100000
C      ICOUPLE = 1      ! for standalone mode this should probably be the default


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

      do i=1,nq   ! first order molecular diffusion terms
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

      if (mbound(LH2).eq.0) then   ! diff limited flux implemented as effusion velocity
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
c  unused...
        DD(LH,1) = DU(LH,1)
        ADD(LH,1) = -ADU(LH,1)
        DD(LH2,1) = DU(LH2,1)
        ADD(LH2,1) = -ADU(LH2,1)

c interior grid points   ?fixed 8-13-05
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
! these are the absorbers which block solar radiation
! If S8 occurs in the gas phase, we use Andy Youngs method to calculate S8 photorates which requires care

         lpolyscount=1
      do k=1,kj

       do i=1,nz

          if (photoreac(k).gt.nq) then  !this gets any SL/IN species that have photorates

             absorbers(k,i)=ABS(SL(INT(photoreac(k)),i))/DEN(I)

c           if (i.eq.1) print *, k,photoreac(k),ISPEC(INT(photoreac(k))),
c     $       absorbers(k,i)

           else if (ISPEC(INT(photoreac(k))).eq.'S8      ') then    !quasi-hardcoded S8 behavior

              absorbers(k,i)=ABS(USOL(INT(photoreac(k)),i))
              if (lpolyscount .eq. 2) absorbers(k,i)=0.0 !S8R doesn't really exist
              !S8L doesn't really either, but this is how we are rolling for now...
              !we are filling up absorbers for S8R and S8, but S8 doesn't have a cross section
              !so only S8R is actually absorbing photons here in the RT scheme.
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

!fill up some heritage vectors rather than try to extract them from the code...
! I am assuming here that the code will always have these species in it so no if's are needed

       JH2O=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'H2O    ')
       JO2_O1D=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'O2    ')
       JO3_O1D=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'O3    ')
       JCO2=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'CO2    ')
!the above finds the first entry in photoreac for the given speices

      do I=1,NZ
       H2O(I) = absorbers(JH2O,I)
       O2(I) =  absorbers(JO2_O1D,I)
       O3(I) = absorbers(JO3_O1D,I)
corig       CO2(I) = FCO2
       CO2(I) = absorbers(JCO2,I)
      enddo



      IDO = 0
      IF (NN.EQ.NSTEPS) IDO = 1
      CALL PHOTO(ZY,AGL,LTIMES,ISEASON,IZYO2,IO2,INO,IDO,timega,frak,
     &      msun,ihztype)
  
      CALL RAINOUT(JTROP,NRAIN,USETD)  !ok



   
      CALL AERCON
 

      
c      print *, photoreac
c      stop

       JO1D=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'O1D    ')
       JO2=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'O2     ')

c      print *, 'stopping after photo'
c      stop

C
C   TIME-DEPENDENT BOUNDARY CONDITIONS
C     (NOTE THAT THE INCLUSION OF
C     A FALL VELOCITY FOR PARTICLES IS MATHEMATICALLY EQUIVALENT TO
C     INCREASING THE DEPOSITION VELOCITY AT THE SURFACE.  AT THE TOP
C     BOUNDARY, ZERO FLUX IS ACHIEVED BY SETTING THE EFFUSION VELO-
C     CITY EQUAL TO THE FALL VELOCITY)

c Shawn's code has some changes to h escape terms here. Consider porting.
C -Mark C
c Using variables for the particle indeces to kill warnings. We should
c eventually migrate this to the read-in. But that will require more
c testing, or someone to figure out how to do that in Mark's input
c framework. -Shawn D-G


      if(USETD.EQ.0) then
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

 
 
c estimate CO2 photolysis above the top of the grid and return CO + O to the upper grid point
c NOTE: this behavior is turned on and off by setting MBOUND=2 in species.dat

c (now above)     JCO2=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'CO2    ')
      JCO2_O1D = JCO2+1
      VCO2 = (prates(JCO2,NZ) + prates(JCO2_O1D,NZ)) * HA
      SMFLUX(LO) = - VCO2*CO2(NZ)*DEN(NZ)
      SMFLUX(LCO) = SMFLUX(LO)


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
      CALL SEDMNT(FSULF,USETD,frak,HCDENS,ihztype)


      if (USETD.EQ.0) then   !particles in main loop

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

c      AERSOL(J,2) = USOL(LS8AER,J)*DEN(J)/CONVER(J,2)     !ACK hardcoded particle numbers
c      AERSOL(J,3) = USOL(LHCAER,J)*DEN(J)/CONVER(J,3)     !ACK hardcoded particle numbers
c      AERSOL(J,4) = USOL(LHCAER2,J)*DEN(J)/CONVER(J,4)     !ACK hardcoded particle numbers
c      print *, AERSOL(J,1),AERSOL(J,2),USOL(LSO4AER,J),CONVER(J,1),
c     &         USOL(LS8AER,J),CONVER(J,2)

!conver is the number of molecules/particle, so that main calcuations are done in molecule space
      !molecules/cm3 *#particles/molecule - > AERSOL [# particles/cm3]
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

      else  !particles in tri-diag

      do J=1,NZ
       do k=1,NP
        AERSOL(J,K)=PARTICLES(J,K)*DEN(J)/CONVER(J,K)
       enddo
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


c new code from eddie
      DO 3 I=1,NQ                            ! Loop through all chemical species
      DO 11 J=1,NZ                           ! Loop through all vertical atmospheric layers
c     R(J) = EPSJ * ABS(USOL(I,J))           ! as it was - USOL should be positive here anyway, can probably remove this (Eddie)
      R(J) = EPSJ * USOL(I,J)                ! R(J) is value to perturb USOL(I,J) by in Jacobian calculation, EPSJ much less than 1
      !!! This is my debug in other version of code - Eddie !!!
      IF(R(J).LT.1.e-100) R(J) = 1.e-100 ! PERTURB DEBUG !!!
      !!! Above ensures no USOL(I,J) falls below double precision limit !!!
  11  USOL(I,J) = USAVE(I,J) + R(J)          ! Add perturbing quantity to mixing ratio
      CALL DOCHEM(FV,0,JTROP,iIN,iSL,USETD)  ! Call the photochemistry routine
c                                           ! FV has dimension (NQ,NZ) and holds gas densities
c end new code from eddie

c old code
c      DO 3 I=1,NQ
c      DO 11 J=1,NZ
c     R(J) = EPSJ * ABS(USOL(I,J))   !as it was
c      R(J) = EPSJ * USOL(I,J)
c  11  USOL(I,J) = USAVE(I,J) + R(J)
c      CALL DOCHEM(FV,0,JTROP,iIN,iSL,USETD)


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


C below is Kevin's test of a distributed sink.  Keeping in case this is ever desired

C distribute H2SO4 sink uniformly between 0 and 20 km
c   I don't know how to keep track of the fluxes
c      jdistSO4 = 21
c      ZSO4 = Z(jdistSO4-1)
c        L = LH2SO4
c        do j=1,jdistso4
c          K = L + (J-1)*NQ
c          rhs(k) = rhs(k) - 0.01  ! per second
c        enddo




      revEu2=0    !testing device to toggle 1st and 2nd order time stepping


      if (N .gt. 1 .and. revEu2 .eq.1) then   !do 2nd order (first step always 1st order R.E.)

c-mc first attempt at second order euler
c  [1/DT*I - a*J]*X = 1/DT*[b*f_{n} - c*f_{n-1} - ?fn] + e*d/dt{f_n}
c
c where a=2/3, b=4/3, c=1/3, d=?,e=2/3, X=f_{n+1} - f_{n}

c DJAC has alread been computed as 1/DT*I - J  , and d/dt(f_{n}) has been computed as RHS


      do j=1,lda
        do k=1,neq
          DJAC(j,k) = 2.0/3.0*DJAC(j,k)
        enddo
      enddo

c-mc could make this shorter by not multiply the upper NQ rows of djac (which are 0 by definition)

c renormalize DTINV
c the previous loop took the main diagonal(which is row KD in band form) to 2/3*1/DT - 2/3J
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

c RHS(K) where K=1,NEQ where NEQ=NQ*NZ, so need to make sure the right usol/height combination is being added at each step:
      do I=1,NQ
       do J=1,NZ
        K = I + (J-1)*NQ
        RHS(K) = RHS(K) +  DTINV*1.0/3.0*(USOL(I,J)-USOLPREV(I,J))
       enddo
      enddo

      endif    !end 2nd order reverse Euler correction loop



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

        IF (I.EQ.LSO4AER.AND. J.GT.J25) THEN  ! a less drastic measure
          USOL(I,J) = USOL(I,J) + RHS(K)
        ELSEIF(I.EQ.LS8AER.AND. J.GT.J25) THEN  ! !!!
          USOL(I,J) = USOL(I,J) + RHS(K)
        ELSEIF(I.EQ.LHCAER.AND. J.GT.J25) THEN  ! !!!
c        ELSEIF(I.EQ.LHCAER) THEN  ! !!!
          USOL(I,J) = USOL(I,J) + RHS(K)
        ELSEIF(I.EQ.LHCAER2.AND. J.GT.J25) THEN  ! !!!
c        ELSEIF(I.EQ.LHCAER2) THEN  ! !!!
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

        ELSEIF(I.EQ.LS4) THEN
                     USOL(I,J) = USOL(I,J) + RHS(K)


        ELSEIF (USOL(I,J).LT. 1.E-20) THEN
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

      DO 4 J=1,JTROP
        USOL(LH2O,J) = H2O(J)
c       USOL(LS8,J) = S8S(J)
   4  CONTINUE


c      do i=1,nq
c       do j=1,nz
c          print *,Z(j), USOL(LSO4AER,j),USOL(LS8AER,j)
c         USOL(i,j)=abs(USOL(i,j))
c       enddo
c      enddo
      !temp


      if (USETD.EQ.0) then
!diddle main loop particles
c 1.e-38 is the smallest number for single precision. We should upgrade
c this to double precision at some point. -Shawn D-G
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

      endif ! end particle diddlation


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
      IF(I.EQ.LSO4AER) MZ = minloc(Z,1, Z/1.e5 .ge. 50)-1 !height index (-1 given the staggered grid)
      IF(I.EQ.LS8AER) MZ = minloc(Z,1, Z/1.e5 .ge. 40)-1 !height index (-1 given the staggered grid)
      IF(I.EQ.LHCAER) MZ = minloc(Z,1, Z/1.e5 .ge. 70)-1 !height index (-1 given the staggered grid)
      IF(I.EQ.LHCAER2) MZ = minloc(Z,1, Z/1.e5 .ge. 70)-1 !height index (-1 given the staggered grid)
      !at some point check/abstract these 40/50km assumptions - this could be easily shunted to the species.dat file...

c      IF(I.EQ.LSO4AER) MZ = minloc(Z,1, Z/1.e5 .ge. 79.5)-1 !height index (-1 given the staggered grid)

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
c 1.e-38 is the smallest number for single precision. We should upgrade
c this to double precision at some point. -Shawn D-G
      smallest = 1e-38

C   FILL UP UPPER PORTION WITH APPROXIMATE ANALYTIC SOLUTION

      do J=MZP1,NZ
        PARTICLES(J,L)=PARTICLES(J-1,L) * EXP(-WFALL(J,L)*DZ(J)/EDD(J))
        PARTICLES(J,L) = MAX(PARTICLES(J,L),smallest)   !needed?
      enddo



   58 CONTINUE
C
      endif  !end tri-diag loop
*****END TRIDIAG

 73   continue           !continue doing newton steps (4 seems to work best)
                         !someday I should see if this is justified
                         !(r37:36 in /td branch has first attempt at convergence checking)



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

!was using these before
c-mc      these are even stricter...
       IF(EMAX.LT.0.15)  DT = 1.1*DTSAVE
       IF(EMAX.LT.0.07)  DT = 1.2*DTSAVE
       IF(EMAX.LT.0.01)  DT = 1.4*DTSAVE
       IF(EMAX.LT.0.008)  DT = 2.0*DTSAVE
       IF(EMAX.LT.0.004)  DT = 3.0*DTSAVE
       IF(EMAX.LT.0.001) DT = 4.0*DTSAVE
       IF(EMAX.LT.0.0005) DT = 5.*DTSAVE

!testing changing to these
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
      print 373, n, TIME, DT,EMAX,ISPEC(IS),ZMAX, USOL(is,js),
     $ USOL(LO2,1), USOL(LCH4,1), USOL(LH2,1), USOL(LCO,1)
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

!THIS STILL NEEDS WORK

        do i=NQ+1,NQ1 !need to fill up rainout and depostion vectors for the triadiagonal species to make the budgets work out.
          SR(I)=0.
            do j=1,JTROP
              !ACK - all particles raining out like H2SO4 !ISOHACK
              SR(I) = SR(I)+ RAINGC(LH2SO4,J)*PARTICLES(J,i-nq)
     $                *DEN(J)*DZ(J)
            enddo

            PHIDEP(I)=(WFALL(1,i-nq)+vturb)* PARTICLES(1,i-nq)*DEN(1)  !ACK - hardcoded turbulent diffusion velocity !ISOHACK - will need to change in ISO
            TLOSS(I) = SR(I) + PHIDEP(I)
c       print *, ISPEC(I),SR(I),PHIDEP(I),SR(I)+PHIDEP(I)  !in general SR>>PHIDEP for particles,by about 100X
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

      !NAN HAPPENING SOMEWHERE BELOW HERE...MUST BE IN OUTPUT.F
  

C
      NS = N/NPR    !NPR=PRN set above to 50 - i.e. write out every 50 steps
      SN = N/PRN
      IF(NN.EQ.NSTEPS) SN = NS
      IF (SN-NS .LT. 1.E-3) THEN
        CALL OUTPUT(NN,NSTEPS,TIME,jtrop, vdep,USOLORIG,USETD,frak)
      ENDIF
          !NAN HAPPENING SOMEWHERE ABOVE HERE...MUST BE IN OUTPUT.F

 
      write(23, 374) n, TIME, USOL(LO2,1), USOL(LH2,1), USOL(LCO,1),
     $ USOL(LCH4,1)!, USOL(LS8aer,1)   !tri-diag (could fix if I wanted...)
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
       write(43,115) N,TIME,DT,EMAX
        do I=1,nz
c not needed in whiff testing         write(43,114), (USOL(K,I)*DEN(I),K=1,NQ)
c         write(43,114), (USOL(K,I)*DEN(I),K=1,NQ)
        enddo
      endif

      print *, RPAR(1,3),RPAR(NZ,3)   

c 114  format(100(1pe10.3))   !ACK hardcoded NQ - update if NQ>100
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
   
  22  CONTINUE   ! successful completion


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
     4 USOL(LNO2,I),USOL(LSO2,I),USOL(LOCS,I),I=1,NZ)
 938  FORMAT(1x,'    Alt      Temp       Den      Press      H2O ',
     2        '      CH4      C2H6       CO2      O2         O3  ',
     3        '      CO       H2CO       HNO3     NO2        SO2',
     4        '      OCS ')

      WRITE(67,738)
      WRITE(67,941) (Z(I),T(I),DEN(I),PRES_bar(I),USOL(LH2O,I),
     2 SL(LCH4,I)/DEN(I),SL(LC2H6,I)/DEN(I),SL(LCO2,I)/DEN(I),
     & SL(LO2,I)/DEN(I),O3(I),
     3 USOL(LCO,I),USOL(LH2CO,I),SL(LHNO3,I)/DEN(I),
     4 USOL(LNO2,I),USOL(LSO2,I),USOL(LOCS,I),I=1,NZ)
 738  FORMAT(1x,'    Alt      Temp       Den      Press      H2O ',
     2        '      CH4      C2H6       CO2      O2         O3  ',
     3        '      CO       H2CO       HNO3     NO2        SO2',
     4        '      OCS ')


 937  FORMAT(1X,1P16E10.3)
 941  FORMAT(1X,1P16E10.3)
      WRITE(164,939)
      WRITE(164,940) (Z(I),T(I),DEN(I),PRES_bar(I),
     & SL(LCO2,I)/DEN(I),SL(LCO,I)/DEN(I),SL(LC2H6,I)/DEN(I),
     & SL(LH2CO,I)/DEN(I),
     & USOL(LCH4,I),USOL(LHNO3,I),USOL(LNO2,I),USOL(LO2,I),O3(I),
     & USOL(LSO2,I),USOL(LH2O,I),USOL(LOCS,I),I=1,NZ)
 939  FORMAT(1X,'    Alt      Temp       Den      Press      CO2 ',
     2       '       CO       C2H6       H2CO     CH4        HNO3',
     3       '       NO2      O2         O3       SO2        H2O ',
     4        '      OCS ')
 940  FORMAT(1X,1P16E10.3)




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
        write(70, 880) ((USOL(k,i),K=K1,K2),i=1,nz)    ! print USOLS into another file .chem
       ELSE
          K2 = NQ
          if (K2-K1 .EQ. 9) then   !this only occurs if NQ is a multiple of 10
            write(18, 880) ((USOL(k,i),K=K1,K2),i=1,nz)
            write(70, 880) ((USOL(k,i),K=K1,K2),i=1,nz)
          else  !if not, need to generate a dynamic format statement
           fmtstr='(  E17.8)'
           write(fmtstr(2:3),'(I1)')K2-K1+1   !OK for one character as this should always be <10
           write(18, fmtstr) ((USOL(k,i),K=K1,K2),i=1,nz)
           write(70, fmtstr) ((USOL(k,i),K=K1,K2),i=1,nz)
          endif
       ENDIF
      enddo

 880  format(10E17.8)
 881  format(5E17.8)

        write (18,881) (T(i),EDD(i),DEN(i),O3(i), SL(NSP-1,i),i=1,nz) !this final one is CO2 number density
        write (66,881) (T(i),EDD(i),DEN(i),O3(i), SL(NSP-1,i),i=1,nz) !print into another file .strctr

        fmtstr='(  E17.8)'
        write(fmtstr(2:3),'(I2)')NP*3
        do i=1,nz
         write(18,fmtstr) (AERSOL(i,j),j=1,NP),(WFALL(i,j),j=1,NP),
     $                  (RPAR(i,j),j=1,NP)
         write(71,fmtstr) (AERSOL(i,j),j=1,NP),(WFALL(i,j),j=1,NP),
     $                  (RPAR(i,j),j=1,NP)           ! print aerosols into another file .aersol
        enddo

        

        fmtstr='(  E17.8)'
        write(fmtstr(2:3),'(I2)')NP

      if(USETD.EQ.1) then
      do i=1,nz
       write(18,fmtstr)  (PARTICLES(i,j),j=1,np)  !ordering is that in species.dat
       write(72,fmtstr)  (PARTICLES(i,j),j=1,np)  !print tridag into another file .tridag
      enddo
      endif

C write out formatted ISOin.dist file
C - ISOHACK!
      fac=1.0

      IROW = 10  !num columns in .dist file
      LR = NISO/IROW + 1
      RL = FLOAT(NISO)/IROW + 1
      DIF = RL - LR
      IF (DIF.LT.0.001) LR = LR - 1
C

      DO L=1,LR
       K1 = 1 + (L-1)*IROW +NQ-NISO-2   !requires isotopic species to be at the end of species.dat list
       !this needs to be re-written to just grab 'S' species or something
       !right now it works if the all the S are at the end AND two hydrocarbon particles are in the last two spots
       K2 = K1 + IROW - 1
c       print *, k1,k2, ISPEC(k1),ISPEC(k2)
       IF (L.lt.LR) then
        write(52, 880) ((USOL(k,i),K=K1,K2),i=1,nz)
       ELSE
          K2 = NQ-2       !temp hack - this needs to be re-written somehow
c       print *, 'final k2',k2,ISPEC(k2)
          if (K2-K1 .EQ. 9) then   !this only occurs if NQ is a multiple of 10
            write(52, 880) ((USOL(k,i),K=K1,K2),i=1,nz)
          else  !if not, need to generate a dynamic format statement
           fmtstr='(  E17.8)'
           write(fmtstr(2:3),'(I1)')K2-K1+1   !OK for one character as this should always be <10
           write(52, fmtstr) ((USOL(k,i),K=K1,K2),i=1,nz)
          endif
       ENDIF
      enddo


        write (52,881) (T(i),EDD(i),DEN(i),O3(i), SL(NSP-1,i),i=1,nz) !this final one is CO2 number density
!the final one is not the CO2 numberdensity in the case of CO2 in the main loop - CHECK the read in here.

        fmtstr='(  E17.8)'
        write(fmtstr(2:3),'(I2)')NP*3
        do i=1,nz
         write(52,fmtstr) (AERSOL(i,j)*fac,j=1,NP),
     $    (WFALL(i,j)*fac,j=1,NP),(RPAR(i,j)*fac,j=1,NP)
        enddo

      if(USETD.EQ.1) then
        fmtstr='(  E17.8)'
        write(fmtstr(2:3),'(I2)')NP
       do i=1,nz
        write(52,fmtstr) (PARTICLES(i,j)*fac,j=1,np)
       enddo
      endif


C - write out file with all the inert species...
!the commented version was for the TESTING code where (S*=S)
c$$$      fac=1.0
c$$$      do i=1,nsp
c$$$       do j=1,nz
c$$$         if (i.le.NQ) then
c$$$           skip= index(ISPEC(i),'S') !used in testing ISOTOPE script (where S*=S)
c$$$           if (skip.ge.1) then
c$$$             if(j.eq.1)print *, 'tempL skipping ISO inert loop',ISPEC(i)
c$$$           else
c$$$              if(j.eq.1)print *, 'loading',i,ISPEC(I)
c$$$            USOLISO(i,j)=USOL(i,j) !load up major species
c$$$           endif
c$$$         endif
c$$$c         if (i.gt.NQ .and. i.le.NQ1) USOLISO(i,j)=PARTICLES(j,i-NQ) ! load up particles   (skipped in ISOTOPE testing mode)
c$$$         if (i.gt.NQ1 .and. i.le.NSP-2) then
c$$$            ind=i-NISO-NP                            !temp
c$$$            if(i.eq.NQ1+1)idex=ind
c$$$            if(j.eq.1)print *, i,ind, ISPEC(i)
c$$$           skip = index(ISPEC(i),'S')
c$$$           if (skip.ge.1) then
c$$$              if(j.eq.1) print *, 'temp skipping SL  ',ISPEC(i)
c$$$            else
c$$$            if(j.eq.1)print *, 'loading',idex,ISPEC(i)
c$$$c            USOLISO(i,j)= SL(i,j)/den(j) ! load up short-lived and inert
c$$$             USOLISO(idex,j)= SL(i,j)/den(j) ! temp
c$$$             if(j.eq.NZ)idex=idex+1
c$$$            endif
c$$$
c$$$
c$$$         endif
c$$$       enddo
c$$$      enddo


      fac=1.0
      do i=1,nsp
       do j=1,nz
         if (i.gt.NQ .and. i.le.NQ1) then  !if there are tri-diagonal species  (causes a gfortran warning about loop over nothing...)
          USOLISO(i,j)=PARTICLES(j,i-NQ) ! load up particles   (skipped in ISOTOPE testing mode)
         else
          USOLISO(i,j)= SL(i,j)/den(j)
         endif
c         if(j.eq.1) print *, ISPEC(i),USOLISO(i,j)
       enddo
      enddo
!the above should load mixing ratios of all species,particles, and short-lived species (in the order of species.dat) into USOLISO



!original versions
c      NXXX=NSP-NPN-2
c      print *, NXXX,NSP,NPN
!full up version
c      NXXX=NSP-2
      NXXX=NSP   !also send CO2 and N2


      IROW = 10  !num columns in .dist file
      LR = (NXXX)/IROW + 1
      RL = FLOAT(NXXX)/IROW + 1
      DIF = RL - LR
      IF (DIF.LT.0.001) LR = LR - 1
C

      DO L=1,LR
       K1 = 1 + (L-1)*IROW
       K2 = K1 + IROW - 1
       IF (L.lt.LR) then
        write(53, 880) ((USOLISO(k,i),K=K1,K2),i=1,nz)
       ELSE
          K2 = NXXX
          if (K2-K1 .EQ. 9) then   !this only occurs if NQ is a multiple of 10
            write(53, 880) ((USOLISO(k,i),K=K1,K2),i=1,nz)
          else  !if not, need to generate a dynamic format statement
           fmtstr='(  E17.8)'
           write(fmtstr(2:3),'(I1)')K2-K1+1   !OK for one character as this should always be <10
           write(53, fmtstr) ((USOLISO(k,i),K=K1,K2),i=1,nz)
          endif
       ENDIF
      enddo

!end ISO HACK

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

       JSO2=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'SO2    ')
       JO2=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'O2     ')+1

! note this is for the O2 + Hv -> O + O reaction,which is the second O2 reaction
! JH2O and JCO2 defined above...
      write(27,299) (Z(I),prates(JSO2,I),prates(JO2,I),prates(JH2O,I),
     $     prates(JCO2,I),I=1,NZ)
 299  FORMAT(5(1PE9.2,2x))
c-mc

c-mc write out all photorates at all heights in the .rates file
      string='(  (1PE9.2,2X))'
      write(string(2:3),'(I2)')KJ+1   !ack - hardcoded kj+1 here. will break if kj>99

!the above is a way to dynamically create a format string at runtime.replaces:  599  FORMAT(55(1PE9.2,2x))
       write(28,string) (Z(I),(prates(J,I),J=1,KJ),I=1,NZ)
c-mc

! write out O2 photolysis rates, for testing the S-R bands using various methods
       do i=1,nz
       write(58,*) Z(I),prates(JO2-1,I),prates(JO2,I) !ack - hardcoded O1D rate number (OK as long as O+O follows O+O1D, which is OK)
       enddo

! print out rainout rates
       do i=1,nq
           write(59, *) (RAINGC(i,j),j=1,jtrop)
       enddo


       JNO=minloc(photoreac,1,ISPEC(INT(photoreac)).eq.'NO    ')

       if(LGRID.EQ.0) write(60,*) (prates(JNO,I),I=1,NZ)  !if using old grid, write out NOprates so they can be used for new grid

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
     &                 max(USOL(LC2H6,I),1.e-60) !EWS debug to prevent floating point errors
  254    FORMAT(1PE9.3,6(E10.2))
  255   CONTINUE
        close(84)
c       endif

         IF (IRESET.eq.0) FAR=1.000E-02     ! Because argon is not normally calculated in photochem
         FCH4=max(USOL(LCH4,1),1.e-60) !EWS - debug to prevent floating point errors
         FC2H6=max(USOL(LC2H6,1),1.e-60) !EWS - debug to prevent floating point errors
         FCO2=SL(LCO2,1)/DEN(1)
         FN2=SL(LN2,1)/DEN(1) + FCO2 !EWS - note that CLIMA/mixing_ratios.dat treats the condensible (CO2) and non-consibles mixing ratios differently
         FO2=USOL(LO2,1)             !      e.g., if the atmosphere is 99% N2 and 1% CO2, then the CO2 fraction is 0.01 and N2 should be set to 1,
         FH2=USOL(LH2,1)             !      because it is 100% of noncondensibles. In practice, N2 should be = (1 - [everything but CO2]). 
c                                    !      But in other parts of photochem, N2 is of the total, not excluding CO2, so we only change it here. 9/8/2015
         FNO2=USOL(LNO2,1)/1.0E60 !gna - clima can't currently cope with NO2 and having it is screwing it up
         JCOLD=JTROP

         WRITE(117,102) FAR, FCH4, FC2H6, FCO2,
     &                  FN2, FO2, FH2, FNO2, JCOLD
  102    FORMAT(1PE10.3,10x,'!Argon'/,1PE10.3,10x,'!Methane'/,
     &          1PE10.3,10x,'!Ethane'/,1PE10.3,10x,'!Carbon Dioxide'/,
     &          1PE10.3,10x,'!Nitrogen'/,1PE10.3,10x,'!Oxygen'/,
     &          1PE10.3,10x,'!Hydrogen'/,1PE10.3,10x,
     &          '!Nitrogen Dioxide'/,I3,17x,'!Tropopause layer')

C MOVING HIGHER UP Reading in the temperature and water profiles from the climate code
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
