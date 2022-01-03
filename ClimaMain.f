cC 1    C     PROGRAM SURFT(INPUT,OUTPUT,TAPE1,TAPE2,TAPE3)
c
C  This program is a modified version of the climate model SURFTEM made
c  by James Kasting. The program has been modified by Michael Mischna (mm),
c  Alex Pavlov (AP), Kara Krelove (KK), Hilary Justh (HJ), Ravi Kopparapu (RK),
c  Ramses Ramirez(RR), and Antigona Segura (AS). Some changes are identified
c  with the initials of the author.

c  The code is mostly written in f77 but is compiled in f90 and it
c  contains some f90 features.

c  This code is a 1-D, cloud-free, radiative-convective climate model.R
c  The calculation of temperature profiles begins with an initial
c  temperature-pressure profile and a solar constant.

c  The net absorbed solar radiation is calculated by a delta two-stream
c  approximation (Toon, et al. JGR Vol. 94, 16287-16301, 1989). It uses
c  4-term correlated k coefficients to parameterize absorption by O3,
c  CO2, H2O, O2 and CH4 in 38 spectral intervals.

c  The IR is calculated by the RRTM routine developed by Mlawer et. al
c  (JGR, Vol.102 (D14), 16663-16682, 1997). It uses 16 term sums in
c  each of its spectral bands in which the k-coefficients are concentrated
c  in areas of most rapidly changing absorption. The version 3.0 of RRTM
c  was implemeted on August/2003 (www.rtweb.aer.com).

c  When the mixing ratio of CO2 is greater than CO2MAX, the maximum
c  level of CO2 that RRTM can manage, the former IR subroutine is used.
c  (Pavlov et al. J. Geophys. Res. 105: 11,981-11,990, 2000).

c  Units in cgs unless otherwise is stated.

c  Temperature in each layer is calculated from:
c              dT/dt = - (g/c_p) dF/dp
c  in this case the derivates are partial. T= temperature, t= time,
c  g= gravitational constant, F=Flux, c_p= Heat capacity, p=pressure.

c  Two types of reach convergence have been set up. One uses a non-strict
c  time stepping mode which is faster and better for high O2-low CO2 runs,
c  like present Earth. The other one is slower but needed on high CO2
c  atmospheres.

c  This model can work alone or coupled to a photochemical model.
c  Modifications for the coupled mode were made by Kara Krelove.

c Input data files required by the program are:
C     Unit   File
C      3     H2O_tables.pdatC      4     solar_data_38.pdat (Read by 2-stream code)
C      8     nearir_expsums.pdat
c      9     CO2_tables.pdat
c     20     ir_expsums.pdat
c        21         BIG_DATAFILE.dat

C
C   THE VERTICAL GRID IS STAGGERED BETWEEN TEMPERATURE AND FLUX
C   GRID POINTS.  THE FLUX GRID IS DEFINED FROM THE VERY TOP OF THE
C   ATMOSPHERE (J=1) TO THE GROUND (J=ND).  THE TEMPERATURE GRID POINTS
C   ARE HALFWAY BETWEEN THE FLUX POINTS, EXCEPT FOR T(ND) WHICH IS
C   LOCATED AT THE GROUND.
C
C   PARAMETERS:
C   ND = # OF ALTITUDE POINTS  (J)
C   NF = # OF FREQUENCIES  (N)
C   NA = # OF ANGLES  (M)
C   NS = # OF CHEMICAL SPECIES  (I)
C   NT = # OF TEMPERATURES IN THE STEAM TABLE
C   NSOL = # OF SOLAR FREQUENCIES
C
C   T = TEMPERATURE (K)
C   P = PRESSURE (bar)
C   Z = LOG PRESSURE + A CONSTANT (ZCON)
C   PF = PRESSURE AT FLUX GRID POINTS
C   ZF = LOG P AT FLUX POINTS
C   ALT = ALTITUDE (KM)
C   GAM = DTDZ
C   BVK = PLANCK FUNCTION
C   LAM = WAVELENGTHS (MICRONS)
C   AV = FREQUENCIES (1/S)
C   TAU = SLANT OPTICAL DEPTH TO OTHER PRESSURE LEVELS
C   F = INTEGRATED NET FLUX
C   FS = INTEGRATED SOLAR FLUX
C   FI = SPECIES MIXING RATIOS   1 = water, 2 = co2, 3 = ch4, 4 = o3, 5 = ethane
C   FH2O - H2O MIXING RATIO
C   T,TN - TEMPERATURES
C   FLAGCONVEC - Tags for the type of convection
c                1. = Water moist adiabat
c                2. = Water dry diabat
c                3. = CO2 adiabat
c                0. = Non convective layer

C-KK        NLAYERS is a translation parameter between this climate model
C-KK    and Mlawer's RRTM code.
C_KK    SurfTem indexes from 1 at the top to ND at the ground, while
C_KK    RRTM indexes from 0 at the ground to NLAYERS at the top.
C-KK        NZ is the number of layers being carried in atm_chem.
      INCLUDE 'CLIMA/INCLUDE/header.inc'
c     PARAMETER(ND=52)
      PARAMETER(NF=55, NA=1, NLAYERS=ND-1, NZ=200)
c     PARAMETER(NF=55, NA=1, NLAYERS=51, NZ=64)
      PARAMETER(NS=3, NS1=NS+2, NS4=NS+5) !gna: changed NS1 from NS+1 to NS+2 to add ethane
      PARAMETER(NT=76, MT=36)
      PARAMETER(NSOL=38, NGS=8, IK=8)  ! Added IK=8 parameter and NGS is 7 now, 3/26/2012
      !gna - changed ngs to 8 (ethane)
      parameter(nrow=11)
      character*8 pstar
C      CHARACTER*5 :: ICH4A   !Changed to make STARR hold up to 5 characters
C      CHARACTER*5 :: ICO2A   !Changed to make STARR hold up to 5 characters
      CHARACTER*5 :: STARR   !Changed to make STARR hold up to 5 characters
      CHARACTER*11 :: AA
      CHARACTER :: DIRINOUT*8,DIRDATA*10
      LOGICAL :: file_e

      DIMENSION TRAD(ND),DZ(ND),Z(ND),ZF(ND)
      DIMENSION temp_alt(NZ), water(NZ), O3(NZ), PRESS(NZ), !EWS - temp_t(NZ) removed because it wasn't used
     &  CH4(NZ), CO2(NZ), ethane(NZ)
      DIMENSION T(ND),TOLD(ND),FTOTAL(ND),FTIR(ND),
     &  FTSO(ND),PF1(ND),DELT(ND),DELTRAD(ND),TN(ND),
     &  DIVF(ND),TCOOL(ND),THEAT(ND),FLAGCONVEC(ND)
      DIMENSION FSATURATION(ND),FSATUR(ND),FSAVE(ND) !EWS dt(ND) removed, not used 8/18/2015
C      DIMENSION newalt(ND),HEATNET(ND),BETA(ND),FCO2V(ND),FH2O(ND)
      DIMENSION HEATNET(ND),BETA(ND),FCO2V(ND),FH2O(ND)
      DIMENSION AVOLD(NF) !EWS - ALAM not used
C-jdh DIMENSION LAM(NF),ALAM(NF),AVOLD(NF)
      DIMENSION PSATCO2(ND) !EWS - PML(ND) removed because it wasn't used
      DIMENSION FNC(ND)        ! Added FNC array c-rr 6/7/2012


c     vectors for gaussian zenith angles
      dimension fdnsoltot(nd), fupsoltot(nd)
      dimension xi(nrow,20), wi(nrow,20), ngauss(nrow)

      REAL*8 newalt ! removed extraneous kappa, kappa_ir, and FLAGCONVE(ND) declarations

      COMMON/DIR/DIRINOUT,DIRDATA
      COMMON/WAVE/AV(NF),LAM(NF),W(NF)
      COMMON/ABLOK/LTYPE(NF,3),XH2O(NF),YH2O(NF),XCO2(NF),YCO2(NF),
     2  AXH(NF),AYH(NF),BXH(NF),BYH(NF),AXC(NF),AYC(NF),BXC(NF),
     3  BYC(NF),PDOP(NF),CPR(NF),TPR(NF),PATH(NS4),PATHP(NS4),
     4  PATHT(NS4),P1,P2,T1,T2,TAU2,TAUP2,
     5  ALPHA(4),BETH2O(4,5,NF),
     6  BETCO2(4,5,NF),CA(19),CB(19),CC(19),CD(19),CE(19),CH(19),
     7  CK(19)
c     COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4
      COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4,FC2H6,FNO2, FI(NS1,ND),FH22
       ! Added FH2 to CBLOK 5/29/2012 c-rr
      COMMON/ALTBLOK/DALT(ND-1),RADIUS(ND-1),PARTICLES(ND),RAER(ND),
     2 ALT(ND)
      COMMON/EBLOK/PG,TG,PG0,IMW,RSURF,OMEGA,POCEAN,IMOIST,
     2  BETA1,BETA2,FVDRY,PDRY
      COMMON/FBLOK/TTAB(NT),PVAP(NT),DPVAP(NT),SVC(NT),DSV(NT),DSC(NT)
     2  ,RHOV(NT),DRHOV(NT),BETAM(70,75),TCP(75),PCP(70),DPDTL(70,75),
     3  DRDTL(70,75)
      COMMON/GBLOK/TCTAB(MT),PCVAP(MT),BETASC(MT),DPCVAP(MT),
     &  DRCVAP(MT),SVSC(MT),DSCC(MT),TKTAB(MT),TCC(25),PCC(36),
     &  BETAMC(25,36),CPC(25,36),DVDTC(25,36),DVDPC(25,36),
     &  DSVC(MT)
      COMMON/SBLOK/P0P,T0P,R,SUBL
      COMMON/PRESSURE/P(ND),PLOG(ND)
       COMMON/PRESS/BETIR1(4,5,NSOL),BETIR2(4,5,NSOL),
     &  kappa_solh2o(NSOL,8,8,IK), kappa_solco2(NSOL,8,8,IK) ! Added new kappa matricies for each of CO2 and H2O coefficients. 8/26/2012

      COMMON/AOZONE/BETAO3(nsol), BETAO2(2),
     &  WGHTO2(NSOL,2)
      COMMON/RSOL/ALPHAZ(4,2),BETAZ(4,2),NPROB(2),
     &  NG(2),SIGG(4,2,NSOL)

!       COMMON/SOLARDATA/weightco2_h2oSOL(16), weights(2,NSOL,IK)
!     & ,kmatrix_sol(NSOL,IK) ! Added SOLARDATA block here 3/26/2012

      COMMON/SOLARBLK/AMU0,SRFALB,OMG0A(NSOL,ND-1),
     &  ASYA(NSOL,ND-1),TAUAER(NSOL),SIGERT(NSOL),FMA(NSOL),PF(ND),
     &  ALAMBDA(NSOL),CGAS(ND,NGS),FUPSOL(ND),FDNSOL(ND),
     &  NGAS(2,NSOL),WGHT(4,2,NSOL),NPR(2,NSOL),SOLINT(NSOL),
     &  TAULAM(ND-1),ASY(ND-1),OMG0(ND-1),FMT(ND-1),QEXT(NSOL,ND-1)


C new common block, von Paris, 21/04/2006
      COMMON/IRDATA/WEIGHTCH4(6),xkappa(3,12,55,8),
     & CIA(7,NF), CPRW(ND,NF)
c-rr !3/23/11 put CIA matrix in IRDATA
      COMMON/VARIR/kappa_irh2o(NF,8,8,IK), kappa_irco2(NF,8,8,IK)! Added kappa matrix in IR for kpsectrum Co2 and H2O coefficients 8/26/2012
      COMMON/weightsIR/ weightco2_h2oIR(IK)



      COMMON/IRBLK/FUPIR(ND),FDNIR(ND),SRFALBIR,OMG0AIR(NF,ND-1),
     & ASYAIR(NF,ND-1),IO3,QEXTIR(NF,ND-1)
      COMMON/HYDROCARB/Qextirst(73,55),w0irst(73,55),
     &  girst(73,55),Qextsolst(73,38),w0solst(73,38),gsolst(73,38),
     &  radstand(73)
      COMMON/CH4BLOCK/ALPHACH4T188(4,17),BETACH4T188(4,17),
     & ALPHACH4T295(4,17),BETACH4T295(4,17),ALPHACH4Kark(4,21),
     & BETACH4Kark(4,21),GAMMAEXP188(17),GAMMAEXP295(17),
     & ALPHACH4NEW(6),BETACH4NEW(17,3,5,6),ALCH4(6,38)
      COMMON/CO2BLOK/betac1,betac2,PC0,TC0,VAPCL0,SUBCL0,DLVCDT
     & ,DLSCDT,CCL,CCS
      COMMON/NO2BLOK/SIGNO2(NSOL)
C
C-KK  Added 6/15/01 to integrate Mlawer RRTM.
      COMMON/ MLAWERI/  layers, numspec, newalt(ND), TempT(0:NLAYERS),
     &         Pres(0:NLAYERS), gasses(7, 0:NLAYERS), COLDEP(ND)
      COMMON/CONSS/C,BK,G,GNEW(ND),PI,SM,DM,DM2   ! Adding DM2 to common block entry 5/3/2011. DM and DM2 are AMN and AMN2 respectively in CONVEC
c-rr  3/29/11
      COMMON/CPHEAT/CPO2(ND),CPCO2(ND), CPN2(ND), CPH2O(ND),
     &  CPN(ND), CPNT(ND), CPH2(ND)  ! Added CPH2 5/31/2012 c-rr
C
      DATA BETA/ND*1./
      DATA BETH2O/1100*0./
      DATA BETCO2/1100*0./
      DATA SIGNO2/1.1E-20, 5.0E-20, 9.5E-20, 2.23E-19, 3.36E-19,
     2   5.1E-19, 5.36E-19, 2.58E-19, 1.07E-19, 8.0E-20, 4.75E-20,
     3   2.65E-20, 1.25E-20, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
     4   0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0./


C-KK        Change to allow differing T-profiles in differing atmospheres.
c-jdh commented out to allow ND!=52
c     DATA ALT /  69., 67.9, 66.8, 65.7, 64.5, 63.1, 61.7, 60.3,
c    2                58.4, 56.4, 53.9, 50.5, 47.2, 44.9, 42.9, 40.9,
c    3                38.9, 36.9, 34.9, 32.9, 30.8, 29.1, 27.6, 26.2,
c    4                24.9, 23.7, 22.5, 21.3, 20.2, 19.1, 18., 16.9,
c    5                15.8, 14.7, 13.6, 12.4, 11.7, 11., 10.3, 9.6,
c    6                8.9, 8.2, 7.5, 6.8, 6.1, 5.3, 4.5, 3.7, 2.9, 2.1,
c    7                1.1, 0.0/

C
C   FREQUENCIES AT ENDS OF SPECTRAL INTERVALS (1/CM)
      DATA AV/40., 100., 160., 220., 280., 330., 380., 440., 495.,
     2  545., 617., 667., 720., 800., 875., 940., 1000., 1065.,
     3  1108., 1200., 1275., 1350., 1450., 1550., 1650., 1750., 1850.,
     4  1950., 2050., 2200., 2397., 2494., 2796., 3087., 3425., 3760.,
     5  4030., 4540., 4950., 5370., 5925., 6390., 6990., 7650., 8315.,
     6  8850., 9350., 9650., 10400., 11220., 11870., 12790., 13300.,
     7  14470., 15000./

      DATA C,HP,BK,SIGMA,PI,SM/3.E10, 6.63E-27, 1.38E-16, 5.67E-5,
     2  3.14159, 1.67E-24/

c Names of the subdirectories for the data, inputs and outputs
      DIRINOUT = 'CLIMA/IO'
      DIRDATA =  'CLIMA/DATA'


c   =============    FILE SECTION ==================

c      print *, 'running'

C  INPUT FILES
      OPEN (unit=1,file= DIRINOUT//'/input_clima.dat')
      OPEN (unit=3,file= DIRDATA//'/H2O_tables.pdat',status='old')
      OPEN (unit=4,file= DIRDATA//'/solar_data_38.pdat',status='old')
      OPEN (unit=8,file= DIRDATA//'/nearIR_expsums.pdat',status='old')

!====================================================================
C New k-coefficients for H2O and CO2 were calculated by Eric Wolf
C using HELIOS-K (https://github.com/exoclime/HELIOS-K), an
C ultrafast GPU-driven correlated-k sorting program
C (Grimm et al. 2015, doi.org/10.1088/0004-637X/808/2/182).
C For H2O we use the HITRAN2016 line-list, assuming 25 cm-1
C line cut-offs using Lorentz profiles and with the plinth removed.
C For CO2 we also use the HITRAN2016 database, but we assume
C 500 cm-1 line cut-offs using the Perrin and Hartman
C (1989, doi.org/10.1016/0022-4073(89)90077-0) sub-Lorentzian
C line profiles.  These conventions represent the current standard
C practices for the treatment of H2O and CO2 lines within coarse
C spectral resolution climate model radiation schemes.
C It is assumed that the H2O self and foreign broadening components,
C and CO2-CO2 CIA, are included elsewhere in the code,
C both of which are independent of the line treatment.
C For further discussions contact eric.wolf@colorado.edu.

      OPEN (unit=15, file=DIRDATA//'/Wolf_HITRAN2016_solar_38_H2O.dat',
     &              status='old')  ! H2O solar coefficients
C
      OPEN (unit=16, file=DIRDATA//'/Wolf_HITRAN2016_solar_38_CO2.dat',
     &              status='old')  ! CO2 solar coefficients

C
      OPEN (unit=17, file=DIRDATA//'/Wolf_HITRAN2016_ir_55_H2O.dat',
     &              status='old')     ! H2O ir coefficients

      OPEN (unit=18, file=DIRDATA//'/Wolf_HITRAN2016_ir_55_CO2.dat',
     &              status='old')     ! CO2 ir coefficients



!====================================================================


      OPEN (unit=9,file= DIRDATA//'/CO2_tables.pdat',status='old')
      OPEN (unit=21,file= DIRDATA//'/BIG_DATAFILE.DAT',status='old')
      OPEN(unit=66, file = DIRINOUT//'/weight_factors.txt')
c      OPEN(unit=10, file=DIRDATA//'/GJ581_1AU.dat',status='old')
c      OPEN(unit=10,
c     .        file=DIRDATA//'/STELLAR_SPECTRA_new.pdat',status='old')
      OPEN(unit=10,
     .        file=DIRDATA//'/STELLAR_SPECTRA_update.pdat',status='old')
!EWS - new stellar spectra updated
      OPEN(unit=30, file=DIRDATA//'/FinalCIAcoeffs2.dat', status='old')
c      OPEN(unit=30, file=DIRDATA//'/FinalCIAcoeffs.dat', status='old')
c      OPEN(unit=89, file=DIRINOUT//'/TPRIND.dat')   !write tau sums experiment
       open(unit=90, file=DIRINOUT//'/FTIR.dat')
       open(unit=91, file=DIRINOUT//'/FTSO.dat')
c  Starting temperature profile
      OPEN (unit=11,file= DIRINOUT//'/TempIn.dat')
c  US standard atmosphere O3 profile used when the climate model is not
c  coupled to the photochemical model - no - it calls SUBROUTINE OZONE for this...
c      OPEN (unit=22,file= DIRINOUT//'/Ozone_standard.dat')

c-mc commenting these out as they are not currently being written to
c-mc If you turn them back on, please direct them to the IO directory and
c-mc add them to the .gitignore files.  We don't want output files clogging up
c-mc version contril
c      OPEN(UNIT=244,file='IHZ.dat')
c      OPEN(UNIT=24,file='waterloss_IHZ.dat')
c  Ozone and water profiles from the photochemical model
c    File formerly called Pass2SurfMP.dat
      OPEN (unit=113,file= 'COUPLE/fromPhoto2Clima.dat')
c  Surface mixing rations to set the chemical composition of the atmosphere.

c gna - eek!  was choosing which mixing_ratios.dat file to read BEFORE
c reading in whether ICOUPLE = 1 or 0!  Moving this block of code down below...

C These INPUT files are open along the program
c  Subroutine IR
c        UNIT           NAME
c         20        DIRDATA/ir_expsums.pdat

c  IMPORTANT Files these are read in the subroutine CHOOSE_STAR
c  IF the character variable STARR is different than "Sun"
c        80     DIRDATA/fluxesKGF_surf.pdat
c        81     DIRINOUT/M star flux (name it as you like)

c  Next files are used for the subroutine AERABSDATA
c         40         DIRDATA/irtotal.DAT
c         41        DIRDATA/soltotal.DAT

C   OUTPUT FILES
c-as Next file has the same structure as TempIn.dat, and should be copied
c-as to TempIn.dat in order to start from the last solution, if IUP=0
      OPEN (unit=12,file= DIRINOUT//'/TempOut.dat')
      OPEN(unit=116,file= 'COUPLE/fromClima2Photo.dat')

      OPEN(UNIT=98,FILE= DIRINOUT//'/clima_allout.tab')
      OPEN(UNIT=96,FILE= DIRINOUT//'/SolarHeating.tab')
      OPEN(UNIT=97,FILE= DIRINOUT//'/clima_last.tab')
      OPEN(UNIT=80,FILE= DIRINOUT//'/IR_wavelength_grid.tab')


c======================================================
c             VARIABLE INPUT PARAMETERS
c======================================================
C      NSTEPS - NUMBER OF ITERATIONS
C         IMW - 0 FOR SATURATED TROPOSPHERE, 1 FOR MANABE/WETHERALD
C               RELATIVE HUMIDITY, 2 FOR M/W WITH CONSTANT
C               STRATOSPHERIC H2O CONTENT
C        RSURF - SURFACE RELATIVE HUMIDITY
C           ZY - SOLAR ZENITH ANGLE (DEGREES)
C        DTAU0 - OPTICAL DEPTH STEP IN SUBLEVEL INTEGRATION
C         ZCON - ARBITRARY CONSTANT ADDED TO Z TO KEEP IT POSITIVE
C           P0 - PRESSURE AT TOP OF GRID
C          PG0 - DRY PRESSURE AT BOTTOM OF GRID (atm)
c            G - Gravity aceleration (cgs)
C          FAC - RATIO OF GRID SPACING AT TOP TO SPACING AT
C                BOTTOM
C          IO3 - 1 TO INCLUDE O3, 0 TO LEAVE IT OUT
C          IUP - SPECIFIES TYPE OF INITIALIZATION (0 IF YOU WISH TO
C                START FROM AN EXISTING SOLUTION, 1 IF YOU WISH TO
C                SPECIFY A NEW SURFACE TEMPERATURE)
C               IF OPTION 1 IS SELECTED YOU MUST MAKE SURE THAT
C               THE STARTING TEMPERATURES ABOVE GROUND LEVEL ARE LESS
C               THAN TG0, SINCE THE TROPOSPHERIC LAPSE RATE IS INTEGRA-
C               TED UPWARDS IN THIS CASE.
C         TG0 - INITIAL SURFACE TEMPERATURE (FOR IUP = 1 CASE)
C      TSTRAT - Stratospheric temperature for IUP=1
C       STARR - Character variable to choose a star, it can be:
c               Sun, F2V, K2V, G2V
c               NOTES: G2V is NOT the Sun.
c               Write it exactly as it is listed.
c               DO NOT FORGET quotation marks.
c    ICONSERV - O = Non strict time-stepping method (faster)
c               1 = Each time step conservs energy (better for high CO2)
c     ICOUPLE - 1 = Coupled to the photochemical model
c               0 = Not coupled
c      SRFALB - Planetary albedo (0.2 for Present Earth)
c      SOLCON - Solar constant (S/So)
c       dtmax - Maximum time step in seconds
c      CO2MAX - Maximum CO2 mixing ratio that RRTM can manage with accuracy,
c               for greater values of CO2 the former IR subroutine is used.
c                ***This version always uses the old IR
C JK   Idry   - If Idry = 0, use the moist adiabat. If Idry = 1, use a dry adiabat
      Idry = 0



      READ(1,51)
!      print *,'reading in'
      READ(1,*) AA,NSTEPS       !step number
      READ(1,*) AA,IMW
      READ(1,*) AA,RSURF
      READ(1,*) AA,zy
      READ(1,*) AA,DTAU0
      READ(1,*) AA,ZCON
      READ(1,*) AA,P0           !Pressure at the top
c      READ(1,*) AA,PG0          !Surface pressure (bar)
c*******Changed for now*********
      READ(1,*) AA,PG0
      READ(1,*) AA,G            !Gravity (Mars=373., Earth=980.)
      READ(1,*) AA,FAC
      READ(1,*) AA,IO3                !Ozone?
      READ(1,*) AA,IUP
      READ(1,*) AA,TG0                !Surface temperature for IUP=1
      READ(1,*) AA,TSTRAT       !Stratospheric temperature for  IUP=1
      READ(1,*) AA,STARR        !What star?
      READ(1,*) AA,ICONSERV     !Type of energy conservation
      READ(1,*) AA,ICOUPLE      !Coupled(1) or not(0)
      READ(1,*) AA,SRFALB       !fixed planetary albedo (0.2)
      READ(1,*) AA,SOLCON       !SOLCON=S/So
      READ(1,*) AA,dtmax        !maximum time step allowed (seconds)
      READ(1,*) AA,CO2MAX
      READ(1,*) AA, IMET        ! IMET (flag 0 or 1)
      READ(1,*) AA, IMETETH     ! IMETETH (flag 0 or 1)
      READ(1,*) AA, nga
      READ(1,*) AA, IHAZE       ! IHAZE (flag 0 or 1)
      READ(1,*) AA, ihztype
      READ(1,*) AA, icealbedo
      READ(1,*) AA, INVERSE
      READ(1,*) AA, FRAK        !can get a fractal haze without being coupled now
!      print *, 'inverse', inverse


!gna - moved this part here so now we know what ICOUPLE is supposed to be
!(before this piece of code was before ICOUPLE was read in)
      IF (ICOUPLE.eq.0) THEN
         OPEN (unit=114,file= DIRINOUT//'/mixing_ratios.dat')
      ELSE
         OPEN (unit=114,file= 'COUPLE/mixing_ratios.dat')
      END IF

      print *, 'icouple is', ICOUPLE

!gna - read more inputs from photo for coupling
      IF (ICOUPLE.eq.1) THEN
      OPEN(unit=999,FILE= 'COUPLE/coupling_params.out')
! 107  FORMAT(1X, F4.2, 5X, F8.3, 5X, F3.1, 5X, I2, 5X, I2,
!     &     9X, I4, 6X, F4.2, 6X, F8.2)

! 107  FORMAT(1X, F4.2, 5X, F8.3, 5X, F3.1, 5X, A8, 5X, I2,
!     &     9X, I4, 6X, F5.3, 6X, F7.3)
      READ(999,*)
      READ(999,*) timega, P0ground, frak, pstar, ihztype, nzp, fscale,
     & G
      print *, 'COUPLING PARAMETERS ARE:'
      print *, 'TIMEGA = ', timega
      print *, 'P0ground = ', P0ground
      print *, 'frak = ', frak
      print *, 'pstar = ', pstar
      print *, 'ihztype = ', ihztype
      print *, 'nzp = ', nzp
      print *, 'fscale = ', fscale
      print *, 'G = ', G
      print *, 'RSURF = ', RSURF

      !remove haze in input file if ihztype = 99 (this means no hcaer was run in PHOTO so nonsensical to include it)
      if(ihztype.eq.99) IHAZE = 0

      IF (pstar == '13') STARR = "Sun"
      IF (pstar == '14') STARR = "Sun"
      IF (pstar == '15') STARR = "ADLEO"
      !using the stellar parameterization implemented by Ramses for these next few
      !see pickstar.f for details
      IF (pstar == '16') STARR = "B5034" !adleo ("adleo" isn't working well) B5035
      IF (pstar == '17') STARR = "B5032" !T3200
      IF (pstar == '18') STARR = "B5050" !K2V
      IF (pstar == '19') STARR = "B4070" !F2V
      IF (pstar == '76') STARR = "B5034" !GJ876
      IF (pstar == '20') STARR = "B5026" !M8V
!      IF (pstar == 21) STARR = "B5030" !M5V
      IF (pstar == '21') STARR = "M5V" !EWS - testing this implementation
      IF (pstar == '24') STARR = "B5050" !K2.5V
      IF (pstar == '25') STARR = "B5042" !K6V
      IF (pstar == '26') STARR = "B5052" !K1V
      IF (pstar == '27') STARR = "B5040" !M0V
      ! Below are Teal's MUSCLES stars, see the end of the pickstar.f
      ! file
      IF (pstar == '80') STARR = "G80" ! GJ876
      IF (pstar == '81') STARR = "G81" ! GJ551 (Proxima)
      IF (pstar == '82') STARR = "G82" ! GJ581
      IF (pstar == '83') STARR = "G83" ! GJ667c
      IF (pstar == '84') STARR = "G84" ! GJ1214b
      IF (pstar == '85') STARR = "G85" ! GJ176
      IF (pstar == '86') STARR = "G86" ! GJ436
      IF (pstar == '87') STARR = "G87" ! GJ832
      IF (pstar == '88') STARR = "G88" ! HD40307
      IF (pstar == '89') STARR = "G89" ! HD40307
      IF (pstar == '90') STARR = "G90" ! HD40307
      IF (pstar == '91') STARR = "G91" ! HD40307

         age = 4.7
         time = age-timega
         SOLCON = (1+0.4*(1-time/4.7))**(-1)*FSCALE
         PG0 = P0ground
         print *, 'STAR = ', STARR

      !correction to SOLCON based on kopparapu HZ (earth distance ~> moist IHZ)
      !IF (msun.eq.16) SOLCON = SOLCON * 0.870
      !IF (msun.eq.17) SOLCON = SOLCON * 0.859
      !IF (msun.eq.18) SOLCON = SOLCON * 0.950
      !IF (msun.eq.19) SOLCON = SOLCON * 1.110
      !IF (msun.eq.76) SOLCON = SOLCON * 0.866

      ENDIF !icouple = 1


  51  FORMAT(4/)

      DO I = 1, ND
        GNEW(I) = 0.0D0 ! initialize GNEW
       ENDDO


c      print *, 'Hello1'
      WRITE(80,3000)
 3000 FORMAT(1X,'Wavelength intervals in 1/cm')
      write(80,3001)
 3001 FORMAT(/1x,'Int',2x,'Wavel',2x,'Waveu')
      wavel = 0.
      waveu = av(1)
      i = 1
      write(80,3002) i,wavel, waveu
 3002 format(1x,i2,2f7.0)
      DO I=2,NF
      write(80,3002) i,av(i-1),av(i)
      enddo
C
C **** Read the Gauss points and weights for the solar zenith angle int
      call data_grabber(xi,wi,ngauss)

c*********Calculate PGO*************
c sk      PG0 = .8 + PCO2
c      print 999, PG0
c999   FORMAT(1x,'PG0 =',1PE12.5)
c      print 999, PCO2
c===================================================================

c Reading the atmospheric composition from mixing_ratios.dat
         READ(114,*) FAR                  !Argon
         READ(114,*) FCH4                        !Methane
         READ(114,*) FC2H6                !Ethane
         READ(114,*) FCO2                        !Carbon dioxide
         READ(114,*) FN2                        !Nitrogen - added Nitrogen mixing ratio c-rr 6/5/2012
         READ(114,*) FO2                        !Oxygen
         READ(114,*) FH22                        ! c-rr 5/29/2012 added H2 mixing ratio
         READ(114,*) FNO2                        !Nitrogen dioxide
         READ(114,*) Jcold                !Tropopause layer read in, Photochem gives out the Tropopause layer on it's own grid, so that is not compatible

         INQUIRE(FILE="COUPLE/aux_c_couple.dat", EXIST=file_e)
        IF (IUP.EQ.0.AND.file_e) THEN
          OPEN (unit=333,file= 'COUPLE/aux_c_couple.dat')
          READ(333,*)
          READ(333,*) aaa, aaa, aaa, aaa, Jcold, aaa
          print *, 'Jcold (for coupling) is =', Jcold
          close(333)
        END IF


c***********Calculate new FCO2**************
c sk        FCO2 = PCO2/((.8/28.+PCO2/44.)*44.)
c        print 997, FCO2
c997   FORMAT(1x,'FCO2 =',1PE12.5)
c WJL- changed mixing ratios below to include NO2

c-rr Added methane flag. Ensures FCH4 is 0 when methane flag is turned off 5/2/11
        IF ((IMET.eq.0).and.(IMETETH.eq.0)) THEN
           FCH4=1.e-60
           FC2H6=1.e-60
        ENDIF

        IF ((IMET.eq.1).and.(IMETETH.eq.0)) THEN
           FC2H6=1.e-60
        ENDIF



        IO2 = 0
        IF (FO2.ge.1e-40)IO2 = 1


c Nitrogen mixing ratio
!      FN2 = 1. - FO2 - FAR - FCH4 - FNO2 - FH22 ! c-rr 5/29/2012 added H2 mixing ratio


c-rr Noncondensible molecular weight of the atmosphere when CO2 is condensing (for a colder planet)          5/3/2011
      DM2 = 28.*FN2 + 32.*FO2 + 40.*FAR + 16.*FCH4
     & + 46.*FNO2+ 2.*FH22 ! c-rr 5/29/2012 added H2 mixing ratio

c jfk DM is the noncondensible molecular weight when CO2 is not condensing
      DM = 44.*FCO2 + (1.-FCO2)*DM2

c        print*,FCH4, FCO2, FO2, JCOLD

c      IF(FCO2.gt.CO2MAX) print 550

      LAST = 0
      AMU0 = COS(ZY * PI/180.)

C   CONSTANT FACTORS (cgs)
      BCON = 2.*HP/C/C
      HK = HP/BK
      BKM = BK/(SM*G)
      ND1 = ND - 1
c      print *, 'Hello2'
      R = 1.9872
      P0P = 6.103E-3
      T0P = 273.15
      SUBL = 677.
c      print *, 'Hello10'
c  TRIPLE POINT PARAMETERS FOR CO2
      PC0 = 5.179
      TC0 = 216.56
      VAPCL0 = 83.2765
      SUBCL0 = 130.893
      DLVCDT = - 0.4817
      DLSCDT = - 0.1732
      CCL = 0.5
      CCS = 0.3



C Read Solar Data
      CALL READSOL

c Choosing a star
c-rr      CALL CHOOSE_STAR(STARR,SOLINT)  3/29/11
       CALL pickstar(STARR,SOLINT)
C try to accelerate ir.f, von Paris, 21/04/2006
c     CALL IREXPSUMS(WEIGHT,XKAPPA)
      CALL IREXPSUMS
c Reading an initial temperature and water profile
  998 FORMAT(3x,F16.12,7x,E22.15)
c-jdh format statement to read in TempIn_Standard.tab
c 998 FORMAT(3x,F7.3,7x,E9.3)
      IF(IUP.EQ.0) THEN
        DO J = 1,ND
          READ(11,998) T(J), FSAVE(J)
        fi(1,j)=FSAVE(J)
        END DO
        TG=T(ND)
      ENDIF
c      print *, 'Hello11'



c Reading the ozone and water from the photochemical model

c      IF(ICOUPLE.EQ.1) THEN
c        DO JREAD=1,NZ  !number of layers in photochem code
c         READ(13,*) temp_alt(JREAD),PRESS(JREAD),O3(JREAD),water(JREAD) !want to put in methane there
c         temp_alt(JREAD)=temp_alt(JREAD)/1.0e5
C         print *, 'temp_alt, press, o3, water'
C         print *, temp_alt(JREAD),PRESS(JREAD),O3(JREAD),water(JREAD)
c        END DO
c  352   FORMAT("Alt = ",1PE12.3," H20=",1PE12.3)

c  Interpolate the grid from the photochemical model to the grid of the
c  climate model
c        CALL INPUT_INTERP(temp_alt, water, O3, Jcold, T, FI)
c        DO J=1,ND
c         print *, alt(J), (FI(I,J),I=1,4)
c        ENDDO
c       ENDIF

C  Initialize pressure grid
      IF(IUP.EQ.1) TG = TG0
c      PRINT *, "Calling Grid()..."
      CALL GRID(P0,FAC,ZCON,Z,ZF,DZ)

c  Reading the US Standard Atmosphere ozone profile
        if(IO3.eq.1.and.ICOUPLE.eq.0) then
           CALL OZONE(FI,P)
c          do i=1,ND
c             read(22,*) x, FI(4,i)
c          enddo
        endif

C   CONVERT FREQUENCIES TO UNITS OF 1/SEC AND COMPUTE WEIGHTING FACTORS
      DO N=1,NF
        AV(N) = C*AV(N)
      END DO
      W(1) = AV(1)
      DO N=2,NF
        W(N) = AV(N) - AV(N-1)
      END DO
C
C   CENTER FREQUENCIES IN MIDDLE OF INTERVALS AND COMPUTE WAVELENGTHS
      SAV = 0.
      DO N=1,NF
        AVOLD(N) = AV(N)/C
        SAV2 = AV(N)
        AV(N) = 0.5*(AV(N) + SAV)
        LAM(N) = 3.E14/AV(N)
        SAV = SAV2
      END DO

c -rr        Tstratospheric iteration loop 4/22/2011
c        DO KK = 1,4

c Constructing temperature and water profiles in case they are not provided
       IF(IUP.EQ.1) THEN
          JCOLD = 1
          CALL PROFILE(TSTRAT,P,T,DZ,FSAVE,FCO2V,BETA,JCOLD,
     &    IDRY,FLAGCONVEC)
       ENDIF
c      print *,' JCOLD =',JCOLD
c Building the water profile
      if(ICOUPLE.eq.0)then
        DO J = 1,ND
         FI(1,J)=FSAVE(J)
        END DO

c jfk 6/25/08 Added four lines below
       IF(IMW.EQ.2) THEN
        DO J=1,JCOLD
        FI(1,J) = 4.E-6
c        print *,'j =',j,'  fi(1,j)=',fi(1,j)
        END DO
       END IF

      else !icouple.eq.1

       DO J=1,ND
         CALL SATRAT(T(J),PSAT)
         FSATURATION(J) = (PSAT/P(J))*RELHUM(P(J))

C-KK   The following line was modified to finish filling H2O grid.
         END DO

!       print *, 'JCOLD original is', JCOLD

       !first try to make JCOLD more sensible - giada
       !it's looking for where FSATURATION goes above 1
       !testing needed for all situations when it needs to change JCOLD?

       if (IUP.EQ.0.and.P(ND).GE.0.93.AND. .NOT.file_e) then !it's starting w/ fresh JCOLD otherwise ! EWS- and pressure is high enough
       IF (FSATURATION(1).GT.1) THEN !test if > 1 at start of grid
          JCOLD_NEW = -1
          DO j =1, ND
             IF ((JCOLD_NEW .EQ. -1).and.(FSATURATION(J).LT.1)) THEN
                JCOLD_NEW = J
             END IF
          END DO
          !Update JCOLD if needed
          IF (JCOLD_NEW.NE.-1) THEN
             JCOLD = JCOLD_NEW
          END IF
       END IF
            JCOLD = max(JCOLD,13) !EWS ensure JCOLD isn't too small.
       end if     !end cold trap search


C       print *, 'JCOLD updated is ', JCOLD
       DO J=1, ND
          IF (J .GE. JCOLD) FI(1,J)=FSATURATION(J)
       END DO
       FI(1,JCOLD)=(3.*FI(1,JCOLD+1)+FI(1,JCOLD)+3.*FI(1,JCOLD-1))/7.
      print*, 'it is in that loop that is icouple=0'

c
c jkf 6/26/08 Change H2O initialization in the stratosphere
       do j=1,jcold
       if (imw.eq.2) FI(1,J) = 4.e-6
       end do
      endif !end if icouple.eq.0



      DO 2 J=1,ND
      PF1(J) = PF(J)*1.E6      !PF1 in dyn/cm^2
      TOLD(J) = T(J)
      FI(2,J) = FCO2
      IF(IUP.EQ.1) FI(2,J)=FCO2V(J)
      FI(5,J) = max(FC2H6,1.e-60) !ethane !EWS - debug for low mixing ratios
   2  FI(3,J) = max(FCH4,1.e-60)  !methane !EWS - debug for low mixing ratios
c
c jfk 6/27/08
      do j=1,nd
      fsave(j) = fi(1,j)


      end do







C *** Initial time step
c-rr         3/29/11 Changed time step from 5.e3 to 2.5e3
c       dt0 = 5.e3
        dt0=5.e3
       IFLAGTIME = 0
       TIME = 0.

c  Altitude calculation
      CALL ALTITUDE(NST,T,FI,DZ)


c Reading the ozone and water from the photochemical model

      IF(ICOUPLE.EQ.1) THEN
c        print *, 'temp_alt, press, o3, water, ch4, co2'
        DO JREAD=1,NZ  !number of layers in photochem code
         READ(113,*) temp_alt(JREAD),PRESS(JREAD),O3(JREAD),
     &                 water(JREAD),CH4(JREAD), CO2(JREAD),
     &                 ethane(JREAD)
         temp_alt(JREAD)=temp_alt(JREAD)/1.0e5
c         print 353, temp_alt(JREAD),PRESS(JREAD),O3(JREAD),water(JREAD),
c     &         CH4(JREAD), CO2(JREAD)
        END DO
        FH2O=water(1)
c        FI(1,ND)=FH2O
        FCH4=CH4(1)
c        FI(3,ND)=FCH4
        FO3=O3(1)
c        FI(4,ND)=FO3
        FCO2=CO2(1)
c        FI(2,ND)=FCO2
        FC2H6 = ethane(1)
        IF(FC2H6.LT.1.e-60) FC2H6 = 1.e-60 !!! Debug to prevent memory underflow issues - Eddie (8/3/2015)
        IF(FCH4.LT.1.e-60) FCH4 = 1.e-60   !!! Note that these are read whether or not IMETH or IEMETH flags are set
        print *, 'FC2H6 is ', FC2H6
c  352   FORMAT("Alt = ",1PE12.3," H20=",1PE12.3)
c  353   FORMAT(6(1PE9.2,1x))
c  Interpolate the grid from the photochemical model to the grid of the
c  climate model

        CALL INPUT_INTERP(temp_alt, water, O3, CH4, CO2, ethane, Jcold,
     &   T, FI)
!        print *,'called input_interp'
c        print *, 'temp_alt,water,co2,ch4,o3(after input_interp)'
c        DO J=1,ND
c         print 353, alt(J), (FI(I,J),I=1,5)
c        ENDDO
       DO J=1, ND
          IF (J .GE. JCOLD) FI(1,J)=FSATURATION(J)
       END DO
       ENDIF

!     Initializing FNC c-rr 6/7/2012
      do J = 1,ND
      FNC(J) = 0.0
      enddo
      do J = 1, ND
      FNC(J) = 1. - FI(1,J) - FI(2,J)   ! Added initial FNC array c-rr 6/7/2012
c     print *, 'IN CLIMA_FI(1,J)=', FI(1,J), J
c     print *, 'IN CLIMA_FI(2,J)=', FI(2,J), J
c     print *, 'IN CLIMA_FNC=', FNC(J), J
      enddo

       close(113)

c Initial non-condensible mixing ratio at surface (used for write statement in output file) 6/7/2011
      FNCI = FNC(ND) ! c-rr 5/29/2012 added H2 mixing ratio

c Aerosol calculation (commented when not used)
      CALL AERABSDATA(FRAK, ihztype)
      CALL GRIDAER(ICOUPLE, IHAZE)
      CALL INTERPAR1(RAER)
C***********************************************************
C ****************** START ITERATIVE LOOP *******************
      DO 40 NST=1,NSTEPS
C************************************************************
      print *, 'TIME STEP = ', NST
      ITROP = 1
c      PRINT 161,NST
c 160  FORMAT(/1X,"---------------------------------------------",
c     2  //1X,"NST =",I6)
c 161  format(1x,"NST =", I6)
      TIME = Time + dt0
!      print *, 'found time'

C Set up gas concentrations for Solar code and former IR code



      CALL GASCON(T,PF,FO2,FH22,FI,FNC,CGAS,NST)  ! Added FH2 to GASCON input argument 5/30/2012
!      print *, 'called gascon'
C SWITCH FOR OLD IR CODE
c      if(FCO2.gt.CO2MAX) then
c-as Old subroutine to calculate IR flux
C PLANCK FUNCTION WAS CHANGED

c-rr gna  Created IRME.F (IR clone with methane and ethane loops turned on). When there is methane call IRM instead of IR. 5/2/2011

      IF (IMET.eq.0) THEN
!       print *, 'calling ir.f'
       CALL IR(T,PF,P,FNC,CGAS)  ! Passes FNC to IR c-rr 6/7/2012
      ENDIF

       IF ((IMET.eq.1).and.(IMETETH.eq.0)) THEN
!       print *, 'calling IRM'
       CALL IRM(T,PF,P,FNC,CGAS) ! Passes FNC to IRM c-rr 6/7/2012
!       print *, 'called IRM'
       ENDIF

      IF (IMETETH.eq.1) THEN
!       print *, 'calling IRME'
       CALL IRME(T,PF,P,FNC,CGAS) ! Passes FNC to IRM c-rr 6/7/2012
!       print *, 'called IRME'
      ENDIF

c      else
C    Code modified 6/15/01 to integrate Mlawer's RRTM
c      CALL TRANSLATEM(G,FI,T,PF,ND1,DM,BKM)
C IR subroutine v3.0 loaded August/2003 (www.rtweb.aer.com)
c      CALL RRTM
c      endif

      IF (NST .EQ. NSTEPS) LAST = 1
C  Solar code
c =================================================================
c  approximating the solar zenith angle with a gaussian summation
c      print *, 'Hello4'
      do j = 1, nd
      fdnsoltot(j) = 0.
      fupsoltot(j) = 0.
      enddo

C Find the right row in the matrix
        do i=1,11
        isave = i
        if (ngauss(i).eq.nga) exit
        enddo
C  isave holds the correct row number for the matrix


      do k=1,nga
        amu0 = xi(isave,k)
        zy = acos(amu0)*180./3.14159
c-rr    setting zenith angle to 60 degrees when nga = 1.
        if (nga.eq.1)then
        amu0 = 0.5
        zy=60.
        endif

        weightt = wi(isave,k)


C Heat capacity calculation
      DO J=1,ND-1
c-rr 3/30/11 The new CPCO2 and CPN2 curve fit equations
      CPCO2(J) = 5.89 + 6.06E-3*T(J) + 2.39E-5*T(J)*T(J)
     &           -3.44E-8*T(J)*T(J)*T(J)
c        if(j.eq.1)print *, 'CPCO2new=', CPCO2(J), T(J)

      CPN2(J) = 6.76 + 6.06E-4*T(J) + 1.3E-7*T(J)*T(J)
      CPO2(J) = 7.47 -4.84E-3*T(J) + 1.38E-5*T(J)*T(J)
     &          -8.73E-9*T(J)*T(J)*T(J) - 1.76E-9/T(J)/T(J)
      CPH2(J) = 7.17e-11*T(J)*T(J)*T(J)*T(J)
     & -1.0e-07*T(J)*T(J)*T(J) + 4.77E-05*T(J)*T(J)
     & -8.10E-03*T(J) + 7.17
      CPH2O(J) = 7.46 +4.52E-3*T(J)-1.38E-5*T(J)*T(J)
     &      + 1.74E-08*T(J)*T(J)*T(J)

c-rr   old curve fit equations
c      CPCO2(J) = 7.7 + 5.3E-3*T(J) - 8.3E-7*T(J)*T(J)
c        if(j.eq.1) print *, 'CPCO2old=', CPCO2(J), T(J)
c       CPN2(J) = 6.76 + 6.06E-4*T(J) + 1.3E-7*T(J)*T(J)


c      CPO2(J) = 8.27 + 2.58E-4*T(J) - 1.877E5/T(J)/T(J)


      CPO2(J) = AMAX1(CPO2(J),CPN2(J))


c-rr Recalculation of mixing rations for the noncondensibles
c FI(1,J)= water
c FI(2,J)= carbon dioxide
c FI(3,J) = methane

c The condensibles are water and carbon dioxide. Water convects for planets closer in and CO2 condenses for planets further out



       CpNC = FN2*CPN2(J) + FO2*CPO2(J) + FAR*4.97 +FCH4*8.3

c Total heat capacity
       CPN(J) = FI(1,J)*CPH2O(J)+FI(2,J)*CPCO2(J) + FNC(J)*CPNC

c-rr This is the total cp of all the gases (condensible + noncondensible)
c      CPN(J) = FI(2,J)*CPCO2(J) + FN2*CPN2(J) + FO2NC*CPO2(J) +   ! CPN Modified to reflect above mixing ratio changes 5/3/2011
c     &         FARNC*4.97 +FCH4NC*8.3

C since CPN is in calories/mol/K we should convert them to erg/g/K
      CPNT(J) = CPN(J)*4.18*1.E7/DM
!     CPNT(J) = CPN(J)*4.18*1.E7/DM2 To calculate total CP when CO2 is condensing right????  5/3/2011
      ENDDO

C   Surface heat capacity (assumes a 50 cm deep ocean mixed layer)
c   Units erg/K/cm^2
      CPNT(ND) = 50.* 4.18*1.E7

C
c        print 1300,nga,ng2,k,amu0,zy,weightt
c1300    format(1x,'nga=',i2,' ng2=',i2,' k=',i2,' amu0=',f8.5,
c     2  ' zy =',f5.2,' weightt=',f8.5)



        IF ((IMET.eq.1).and.(IO2.eq.1))THEN
        CALL SOLARMOX(T,LAST,FNC,NST)
        ELSEIF((IMET.eq.1).and.(IO2.eq.0))THEN
        CALL SOLARM(T,LAST,FNC,NST)
        ELSEIF((IMET.eq.0).and.(IO2.eq.1))THEN
        CALL SOLAROX(T,LAST,FNC,NST)
        ELSEIF((IMET.eq.0).and.(IO2.eq.0))THEN
        CALL SOLAR(T,LAST,FNC,NST)
        ENDIF




          do j=1,nd
c          if (j.eq.1)print *, 'fdnsoltot', fdnsoltot(j)

          fdnsoltot(j) = fdnsoltot(j) + fdnsol(j)*weightt

          fupsoltot(j) = fupsoltot(j) + fupsol(j)*weightt
          enddo
c          print *, 'FCH4=', FI(3,1)
c          print *, 'FC2H6=', FI(5,1)
c        print *, 'fdnsol',f
c      print 1301,fdnsol
c1301  format(1x,1p8e9.2)
c        print *
c      print *, 'fupsol'
c      print 1301,fupsol
      enddo
c
c      print *
c      print *,'Integrated fluxes'
c        print *, 'fdnsoltot'
c      print 1301,fdnsoltot
c        print *
c      print *, 'fupsoltot'
c      print 1301,fupsoltot


c =================================================================


c IR and SOLAR fluxes (erg/cm^2/s)
      DO 31 J=1,ND
c      if (j.eq.1) print *, 'FDNSOL=', FDNSOL(1)
      FDNSOL(J) = SOLCON*0.5 * FDNSOLTOT(J)
      FUPSOL(J) = SOLCON*0.5 * FUPSOLTOT(J)
      FTOTAL(J) = FDNSOL(J)-FUPSOL(J)+FDNIR(J)-FUPIR(J)
      FTIR(J) = FDNIR(J)-FUPIR(J)
      FTSO(J) = FDNSOL(J)-FUPSOL(J)

 !     print *, 'FDNSOL(J)', FDNSOL(J)
 !     print *, 'FUPSOL(J)', FUPSOL(J)
 !     print *, 'FTOTAL(J)', FTOTAL(J)
 !     print *,  'FTIR(J)',  FTIR(J)
 !     print *, 'FTSO(J)', FTSO(J)
 !     print *, 'FDNIR(J)', FDNIR(J)
 !     print *, 'FUPIR(J)', FUPIR(J)

 !     call sleep(1)

  31  CONTINUE
      ALBP = FUPSOL(1)/FDNSOL(1)
      SEFF = abs(FTIR(1)/FTSO(1))      !to print out Seff, c-rr 4/21/2011
      print *, 'FTIR= ', FTIR(1)
      print *, 'FTSO= ', FTSO(1)
      PRINT *, 'Seff=',SEFF
      PRINT *, 'JCOLD, alt(JCOLD), P(JCOLD)',JCOLD, ALT(JCOLD),P(JCOLD)
      PRINT 166,ALBP
 166  FORMAT(/1X,"PLANETARY ALBEDO:  ALBP = ",F6.4)

C




C BEGIN INVERSE SKIPS




      IF(INVERSE.EQ.0) THEN !only do if not wanting inverse calculations
C New temperature calculation for all layers from radiative equilibrum
      DO 41 J=1,ND-1
      TN(J)=T(J)-(FTOTAL(J+1)-FTOTAL(J))*dt0*GNEW(J)/CPNT(J)
     &        /(PF1(J+1)-PF1(J))
c        if(j.eq.1) print *, 'CPNT=', CPNT(J), 'FTOTAL(J+1)=',
c     &  FTOTAL(J+1),'FTOTAL(J)=',FTOTAL(J), 'dt0=', dt0,
c     &  'G=', G, 'PF1(J+1)', PF1(J+1), 'PF1(J)', PF1(J),
c     & 'TN(J)=', TN(J), 'T(J)=', T(J)
      TCOOL(J)=-(FTIR(J+1)-FTIR(J))*GNEW(J)/CPNT(J)
     &          /(PF1(J+1)-PF1(J))*86400.
      THEAT(J)=-(FTSO(J+1)-FTSO(J))*GNEW(J)/CPNT(J)
     &         /(PF1(J+1)-PF1(J))*86400.
  41  CONTINUE

c New surface temperature from radiative equilibrum
      select case (ICONSERV)
      case(1)
      TN(ND)=T(ND)+FTOTAL(ND)*dt0/CPNT(ND)
      TCOOL(ND)= FTIR(ND)*86400./CPNT(ND)
      THEAT(ND)= FTSO(ND)*86400./CPNT(ND)
      case(0)
C Lower atmospheric layer temperature calculated from the total flux at
C the TOP of the atmosphere
      TN(ND-1) = T(ND-1)+FTOTAL(1)/(PF1(ND)-PF1(ND-1))*
     &           GNEW(ND-1)/CPNT(ND-1)*dt0
      CALL SATRAT(TN(ND-1),PSAT)
      FI(1,ND-1) = RELHUM(P(ND-1))*PSAT/P(ND-1)
      end select
c      print *, TN(ND), FTOTAL(ND)

* Total heating rate
      do j=1,ND
       HEATNET(j)=THEAT(j)+TCOOL(j)
      enddo
c      print *, 'Hello5'
c-as TRAD is defined for printing and diagnostic purposes
      DO J=1,ND
        TRAD(J)=TN(J)
      ENDDO

cTEMPORARY DEBUGGING STATEMENT**************
C      print *,'Calling output early'
C     GOTO 571
c Calculating tropospheric temperatures
      select case(ICONSERV)
*** Non strict time-stepping model
      case(0)
c Calculation of the ground temperature
c jfk 7/14/08 Redo this! This is NOT how the ground temperature ought
c     to be calculated in this method! The ground temperature should be
c     adjusted in such a way as to balance the fluxes at the top of the
c     atmosphere, i.e., such that DIVF(1)=0. It has been recoded this
c     way below.
      DIVF(1) = FTOTAL(1)/FUPIR(1)
      TN(ND) = T(ND) * (1. + 0.1*DIVF(1))
c      print *,'ftotal(1)=',ftotal(1),' fupir(1)=',fupir(1),
c     2  'x(1)=',divf(1)
c      print *,'T(ND)=',t(nd),' TN(ND)=',tn(nd)
C
c jfk 7/16/08 One needs different logic, depending on whether the
c     surface temperature is increasing or decreasing.
      IF (TN(ND) .LT. T(ND)) GO TO 1400
c
c   Surface temperature is increasing, so do a normal convection
c   calculation, adjusting those layers that are unstable.
      JCONV=ND
      ITROP=1
      DO J1=ND, 2, -1
        FLAGCONVEC(J1) = 0.
        T1=TN(J1)
        DZP = DZ(J1)
        P1 = P(J1)
        P2 = P(J1-1)
        FC1 =FI(2,J1)
        FH1 =FI(1,J1)


        CALL CONVEC(T1,T2,P1,P2,FH1,FH2,FC1,FC2,DZP,ITROP,cflag,
     & Idry, imco2)

c         IF(IO3.EQ.1 .AND. ALT(J1) .GT. 40.) GOTO 1401   ! Skip convection if ozone beyond 40km
c
c  jkf 7/15/08 I am replacing the following logic with simpler logic.
           IF (TN(J1-1) .LE. T2) THEN
                   TN(J1-1) = T2
c                   FI(2,J1-1) = FC2
c jfk 4/15/11 The program is not doing CO2 condensation well. It is also jiggered in
c             CONVEC to stay on the saturation vapor pressure curve for CO2 mixing
c             ratios >0.9. Right now, the code only work for nearly pure CO2 atmospheres.
c             So, let's keep the CO2 vertical profile fixed at the surface CO2
c             mixing ratio.
                   FI(2,J1-1) = FCO2
                     FLAGCONVEC(J1) = cflag
                   JCONV = J1
           END IF
      END DO
      FLAGCONVEC(ND) = 1.
      ND1 = ND-1
      DO J=ND1,1,-1
      FI(2,J) = AMIN1(FI(2,J),FI(2,J+1))
      END DO
      GO TO 1401
c
c   If surface temperature is decreasing, then adjust all temperatures
c   below the cold trap downward by the same amount. This ensures that
c   the upward IR flux will decrease as surface temperature decreases.
 1400 CONTINUE
      DTSURF = T(ND) - TN(ND)
      DO J=JCOLD, ND-1
      TN(J) = T(J) - DTSURF
      ENDDO
      JCONV = JCOLD    ! 5/23/2011 So that it knows what JCONV is when it is not convecting
 1401 CONTINUE
** End of the non strict time-stepping model

*** Here the temperatures are calculated conserving energy on each
*** layer
      case(1)

      DO ITER=1,20        !starting convective adjustment
       ITROP = 1
       imco2=0
       JCONV=ND
c-as   Adjusting the temperature on the surface and the layer above
c-as   it (ND and ND-1) (sept-2004)
       HC1=CPNT(ND-1)*(PF1(ND)-PF1(ND-1))/g
       DZP = DZ(ND)
       T1 = TN(ND)
       P1 = P(ND)
       P2 = P(ND-1)
       FH1 = FI(1,ND)
       FC1 = FI(2,ND)

       CALL CONVEC(T1,TadND1,P1,P2,FH1,FH2,FC1,FC2,DZP,1,cflag,
     & Idry, imco2)


       FI(1,ND-1)=FH2*RELHUM(P(ND-1))
       TnewND=(CPNT(ND)*TN(ND)-HC1*(TadND1-TN(ND)-TN(ND-1)))/
     & (HC1+CPNT(ND))
       TnewND1=TadND1-TN(ND)+TnewND
       TN(ND-1)=TnewND1
       TN(ND)=TnewND
       FLAGCONVEC(ND)= 1.
       FLAGCONVEC(1) = 0.
       imco2=0
c       print *, 'Hello6'
c-as This part has been modified to consider energy balance in each
c-as  layer, as Hilary Justh did it (oct-2003)
***** CONVECTIVE ADJUSTMENT (considering energy balance and convection)

       DO J1=ND-1,2,-1
        T1 = TN(J1)
        DZP = DZ(J1)
        P1 = P(J1)
        P2 = P(J1-1)
        FH = FI(1,J1)
        FC1 = FI(2,J1)
       CALL CONVEC(T1,T2,P1,P2,FH,FH2,FC1,FC2,DZP,ITROP,cflag,
     & Idry, imco2)

c        IF(IO3 .EQ. 1 .AND. ALT(J1) .GT. 40.) GOTO 1403

        IF (TN(J1-1).LE.T2) THEN
          FLAGCONVEC(J1) =  cflag
          IF(cflag.eq.1.or.cflag.eq.3) JCONV=J1
          DELPCP1=(PF1(J1+1)-PF1(J1))*CPNT(J1)
          DELPCP2=(PF1(J1)-PF1(J1-1))*CPNT(J1-1)
          T2P=TN(J1-1)*(DELPCP2/(DELPCP1+DELPCP2))+(T2)*
     &     (DELPCP1/(DELPCP1+DELPCP2))
          T1P=T1-T2+T2P
          TN(J1-1)=T2P
          TN(J1)=T1P
c          FI(2,J1-1)=FC2
          FI(2,J1-1) = FCO2  ! jiggered again
        ELSE
            ITROP=0
            FLAGCONVEC(J1)=0.
        ENDIF
       ENDDO
1403   CONTINUE
      ENDDO          ! End of convection adjustment loop
      ND1 = ND-1
      DO J=ND1,1,-1
      FI(2,J) = AMIN1(FI(2,J),FI(2,J+1))
      END DO
      end select

c      print *, 'Hello7'
c Water recalculation
      DO J=1,ND
       CALL SATRAT(TN(J),PSAT)
       FSATUR(J) = (PSAT/P(J))
      ENDDO
      FCT=FSATUR(ND)
      DO J=ND-1,1, -1
       FCT=AMIN1(FCT,FSATUR(J+1))
      ENDDO
c      print *
c      print *,'FSATUR'
c      print 4309,fsatur
c Finding the cold trap (JCOLD)

        JCOLD = 1
c
c jfk 7/15/08 Simplify the logic for finding the cold trap
c        DO J = ND-1, 2, -1
c           IF (JCOLD .EQ. 1) THEN
c                IF (FSATUR(J) .LT. FSATUR(J-1)) THEN
c                        JCOLD = J
c                  END IF
c           END IF
c        END DO
c
        DO J = ND-1,2,-1
          IF (P(J) .LT. 1.E-01) GO TO 3100
          JCOLD = J
          IF (FSATUR(J-1) .GT. FSATUR(J)) GO TO 3100
        END DO
 3100   CONTINUE
C
c        print*,"JCOLD = ", JCOLD
c        print*,"CONVEC = ",FLAGCONVEC(ND)

c Water from the cold trap to the ground
      DO J = JCOLD, ND
       FI(1,J) = FSATUR(J)*RELHUM(P(J))
       if(imw.eq.2) FI(1,J) = amax1(FI(1,J),4.e-6)
      END DO

c Water from the cold trap to the top (if it is used in the coupled
c mode these values are given by the photochemical code)
       if(ICOUPLE.eq.0)then
         DO J = JCOLD-1, 1, -1
           FI(1,J)= FI(1,JCOLD)
c          if(imw.eq.2) FI(1,J) = 4.e-6
          END DO
        endif

C-KK To smooth over the profile around JCOLD.
      sum = 2*FI(1,(JCOLD-1)) + FI(1,(JCOLD+1)) + 2*FI(1,JCOLD)
      FI(1,JCOLD) = sum/5.

c Saving the former temperature profile
      DO J=1,ND
       TOLD(J) = T(J)
      ENDDO
c      print *, 'Hello8'
C Smoothing of temperature profile conserving energy
       if(ICONSERV.eq.1) then
c jfk Replace the DO logic below to make sure that the smoothing does
c  not occur when CO2 is condensing (FLAGCONVEC=3)
          DO J=2,JCOLD
c         DO J=2,ND-2
c         IF (FLAGCONVEC(J).LT.1.E-3) THEN
          Tj1 = 0.5*TN(J) + 0.25*(TN(J-1) + TN(J+1))
          CPP0=(PF1(J)-PF1(J-1))*CPNT(J-1)
          CPP1=(PF1(J+1)-PF1(J))*CPNT(J)
          CPP2=(PF1(J+2)-PF1(J+1))*CPNT(J+1)
          En1=CPP0*TN(J-1)+CPP1*TN(J)+CPP2*TN(J+1)
          DELT1=Tj1-TN(J)
          DELT2=-(CPP1/(CPP0+CPP2))*DELT1
          TN(J+1)=TN(J+1)+DELT2
          TN(J-1)=TN(J-1)+DELT2
          TN(J)=Tj1
c         ENDIF
         END DO
       endif

C  Diagnostics parameters
      DO J=1,ND
      DELT(J) = (TN(J)-TOLD(J))
      DELTRAD(J) =TRAD(J)-TOLD(J)
      T(J) = TN(J)
      DIVF(J) = FTOTAL(J)/FUPIR(J)
      ENDDO
c      print *,'T after DIVF calculation'
c      print 4309,T
c      print *
c      print *,'JCOLD =',JCOLD
c 4309 format(1p8e9.2)

c Smoothing the temperature profile in the non-strict time step case
c-jdh **check on this when comparing "conserving" vs "non-conserving"**
      if(ICONSERV.eq.0) then
c  Replace the DO logic below to make sure that the smoothing does
c  not occur when CO2 is condensing (FLAGCONVEC=3)
          DO J=2,JCOLD
C         DO J=2,ND-2
C         IF (FLAGCONVEC(J).LT.1.E-3) THEN
          T(J) = 0.5*TN(J) + 0.25*(TN(J-1) + TN(J+1))
C         ENDIF
        END DO
      endif

      print*,'Surface temperature=',T(ND)

c adjust albedo based on ice-albedo feedback
c parameterization added by Giada based on Charnay et al 2014
      if (icealbedo.eq.1) then
          IF (T(ND).LT.240.) SRFALB = 0.65

          IF (T(ND).GT.290.) SRFALB = 0.30

          IF (T(ND).GE.240. .AND. T(ND).LE.290.) then

              SRFALB=0.65+(0.3-0.65)*( (T(ND)-240)/(290-240) )**0.37

         end if
       print *, 'Surface albedo=', SRFALB

       end if


c Adjusting the time stepper
       DTS = dt0
       CHG = 0.
       DO J=2,ND-1
         REL = ABS(DELT(J)/TOLD(J))
         CHG = AMAX1(CHG,REL)
       END DO
       IF (CHG.LT.0.01) dt0 = DTS*1.5
       IF (CHG.LT.0.001) dt0 = DTS*5.
       IF (CHG.GT.0.02) dt0 = DTS/2.
       IF (dt0.GE.dtmax) dt0 = dtmax
c       print *, 'Hello9'
       CALL ALTITUDE(NST,T,FI,DZ)
       END IF !end skipping for inverse model



C /INVERSE


C***********************************************************
c***  WRITING OUTPUT FILES
************************************************************
C 571   CONTINUE
      IF(NST.EQ.1) THEN
      WRITE(98,*)
      WRITE(98,*) "   OUTPUT FILES FOR THE ",STARR
      WRITE(98,*)
       WRITE(98,555) SOLCON,FCH4*FNCI,FCO2,FO2*FNCI,FN2*FNCI,
     &   FH22*FNCI,FAR*FNCI,IO3,IUP

 555  format(1x,"Solar Constant= ",F5.3,3x,"F_CH4= ",1pe10.4,2x,
     & "F_CO2= ",1pe10.4,2x,"F_O2= ",1pe10.4,2x,"F_N2= ",1pe10.4,2x,
     & "F_H2= ", 1pe10.4,2x,"F_AR= ", 1pe10.4,2x,
     & "IO3 = ",I2,3x,"IUP= ",I2)

      WRITE(98,556) ICONSERV,FAC,ND,SRFALB,G,IMW, INVERSE
 556  format(1x,'ICONSERV=',I2,2X,'FAC=',F4.1,2X,'ND=',I3,2X,'SRFALB='
     & ,F6.3,2X,'G=',F6.1,2X,'IMW=',I2, 2X, 'INVERSE=', I2)
      WRITE(98,557) FNO2
 557  FORMAT(/1x,'FNO2 =',1pe10.3)
      WRITE(98,*)
c      if(FCO2.gt.CO2MAX)write(98,550)
      ENDIF
c 550  format(10x,'*** fCO2 > CO2MAX The old IR subroutine is used')
      nsteps2 = nsteps-2
      nsteps3 = nsteps-3
C
C JK Calculate root-mean-square flux divergence above the tropopause
       DIVFrms = 0.
       JC1 = JCONV-1
       DO J=1,JC1
       DIVFrms = DIVFrms + DIVF(J)*DIVF(J)
       END DO
       DIVFrms = SQRT(DIVFrms/JC1)
C
c      if(nst.gt.2 .and. nst.lt.nsteps2) then
       WRITE(98,966) NST,JCONV,CHG,dt0,DIVF(1),DIVFrms,DELT(ND),T(ND)
 966   FORMAT(1x,"NST=",I6,1X,'JCONV=',I3,1x,'CHG=',1pe8.2,1x,"dt0=",
     & 1pe8.2,1X,"DIVF(1)=",1PE9.2,1X,"DIVFrms=",1PE9.2,1x,
     & "DT(ND)=",1PE9.2,1x,"T(ND)=",1PE10.4)
c       ENDIF
c      if(nst.eq.nsteps3) write(98,*)
      if(nst.eq.1 .or. nst.eq.nsteps) then
c       WRITE(98,965) NST,dt0,DIVF(1),FTOTAL(ND-1),FTIR(ND-1),
c    & FTSO(ND-1),DELT(ND),T(ND)
c 965   FORMAT(1x,"NST=",I3,2X,"dt0=",1PE9.3,
c     & 2X,"DIVF(1)=",1PE12.5,2X,"Ftot(ND-1)=",1pe11.4,2x,"FtIR(ND-1)="
c     & ,1pe11.4,2x,"FtSol(ND-1)=",1pe11.4,/,1x,"DT(ND)=",1PE10.3,2x,
c     & "T(ND)=",1PE10.4)
c
      TIMEDAYS = TIME/24./3600.
      WRITE(98,567) TIME,TIMEDAYS
 567  FORMAT(/1X,'TIME=',1PE10.3,2X,'TIME IN DAYS =',E10.3)
      WRITE (98,166) ALBP
      WRITE(98,683)
        DO J=1,ND
c      WRITE(98,680) J,P(J),ALT(J),T(J),FLAGCONVEC(J),
c     & DELT(J),TOLD(J),FI(1,J),HEATNET(J),TCOOL(J),THEAT(J)
      WRITE(98,680) J,P(J),ALT(J),T(J),FLAGCONVEC(J),
     & DELT(J),TOLD(J),FI(1,J),FSAVE(J),FI(4,J),TCOOL(J),THEAT(J)
       ENDDO
      WRITE(98,1683)
      Write (98,1112)
 1112 FORMAT(/1x,'FCO2')
      WRITE(98,1111) (FI(2,J),J=1,ND)
 1111 FORMAT(1x,1p10e9.2)
      WRITE(98,685)
         DO J=1,ND
      WRITE(98,684) J,PF(J),ALT(J),FTOTAL(J),FTIR(J),FDNIR(J),
     & FUPIR(J),FTSO(J),FDNSOL(J),FUPSOL(J),DIVF(J)
         ENDDO
      WRITE(98,*)
      END IF
  683  FORMAT(/2x,"J",5X,"P",9X,"ALT",9X,"T",8X,"CONVEC",
     & 7X,"DT",10X,"TOLD",8x,"FH20",
     &  7x,'FSAVE',8x,'FO3',8x,'TCOOL',7x,'THEAT') ! top of file
 1683  FORMAT(2x,"J",5X,"P",9X,"ALT",9X,"T",8X,"CONVEC",
     & 7X,"DT",10X,"TOLD",8x,"FH20",
     &  7x,'FSAVE',8x,'FO3',8x,'TCOOL',7x,'THEAT')  !bottom of file
 680  FORMAT(I3,3(1x,1PE10.4),1X,1PE9.2,2X,1PE11.4,7(1X,1PE11.4))
 685  FORMAT(/2x,"J",4X,"PF",9X,"ALT",7X,"FTOTAL",7X,"FTIR",7X,"FDNIR",
     & 7X,"FUPIR",7X,"FTSOL",7X,"FDNSOL",7X,"FUPSOL",7X,"DIVF")
 684  FORMAT(I3,2(1x,1PE10.4),8(1X,1PE11.4))
C
      DO J=1,ND
      TSAT = T(J)
      CALL SATCO2(TSAT,PSAT)
      PSATCO2(J) = PSAT
      END DO

      do j=1,nd   ! FH2O becomes old FSAVE for next time step 4/23/2012
      fsave(j) = fi(1,j)
      end do

      do j = 1,nd  ! redefines FNC for next time step c-rr 6/7/2012
      FNC(J) = 1. - FI(1,J) - FI(2,J)   ! Added initial FNC array c-rr 6/7/2012
      enddo


***************************************************************
C   End of iterative loop
  40  CONTINUE
***************************************************************

      if(ICOUPLE.eq.1) then
       print *, 'output photo'
       CALL OUTPUT_PHOTO(T, FI, water, ALT, nzp)
      endif

        WRITE(97,466)
 466  FORMAT(5X,'ALT',10X,'P',10X,'T',10X,'FH2O',11X,'O3',11X,
     2  'THEAT',8X,'TCOOL',8X,'PSATCO2',8x,'FCO2')
      DO J=1,ND
c jkl 6/27/08 Print out H2O from the initial profile
        WRITE(12,998) T(J),FI(1,J)
c        WRITE(97,467) ALT(J),P(J),T(J),FI(1,J),FI(4,J),THEAT(J),
c     &  TCOOL(J)
        WRITE(97,467) ALT(J),P(J),T(J),FSAVE(J),FI(4,J),THEAT(J),
     &  TCOOL(J),PSATCO2(J),FI(2,J)
 467  FORMAT(1PE10.4,2X,1PE10.4,1X,1PE10.4,3X,1PE11.4,3X,1PE11.4,
     & 3X,1PE11.4,3X,1PE11.4,3x,1pe11.4,3x,1pe11.4)
      END DO
c       close(89)
       OPEN (unit=333,file= 'COUPLE/aux_c_couple.dat') !stb

 919  FORMAT(1X,  E12.5, 9x, E10.5,9x, E10.5,9x, E10.5,9x,I3,9x, E10.5)
 918  FORMAT(1X, 'FTOTAL(101)',9x,'T(101)',9x, 'P(JCOLD)',9x, 'P_surf',
     2 9x, 'jcold',9x, 'DIVFrms')
        WRITE(333,918)
        WRITE(333,919) FTOTAL(101),T(101),P(JCOLD), P(101),JCOLD,
     & DIVFrms
        close(333)


c        Iterative Tstrat procedure   c-rr 4/22/2011
c        174K is the stratospheric temperature computed for present Mars given from the
c        point of first CO2 condensation (CONVEC = 3 or when PSAT/P=1).
c       TCONST is the Seff(1-ALBP)*SOLCON (SOLCON is .43 for Mars)
c        TCONST= .7206*.43
c        TSTRAT = 174*((SEFF*(1.0-ALBP))/TCONST)**0.25
c        print *, 'TSTRAT=', TSTRAT
c        ENDDO

       STOP


      END                 !end of the main program

*********************************************************************
      SUBROUTINE ALTITUDE(NST,T,FI,DZ)
      INCLUDE 'CLIMA/INCLUDE/header.inc'
      PARAMETER(NS1=5)       !gna: changed ns1 from 4 to 5
      COMMON/CONSS/C,BK,G,GNEW(ND),PI,SM,DM,DM2
      COMMON/ALTBLOK/DALT(ND-1),RADIUS(ND-1),PARTICLES(ND),RAER(ND),
     & ALT(ND)
      DIMENSION T(ND),FI(NS1,ND),DZ(ND)

c-as  This subroutine calculates the altitude.
c-as  The water vapor was eliminated the first time this subroutine
c-as  is called (before the NSTEPS DO loop) in order to make easier
c-as  the parameter translation to the photochemical model

      ALT(ND) = 0.

      DO J=ND-1,1,-1
       GNEW(J) = G*(6378.**2)/(6378. + ALT(J+1))**2  ! Turned Gravity into arrays 10/12/2012
       BKM = BK/(SM*GNEW(J))
       TA = 0.5*(T(J) + T(J+1))
       FH2O = 0.5 * (FI(1,J) + FI(1,J+1))
       FCO2J = 0.5* (FI(2,J) + FI(2,J+1))

       FNCA = 1. - FH2O - FCO2J
       AM = 18.*FH2O + 44.*FCO2J + DM2*FNCA  ! AM is the weight of entire parcel (noncondensible + condensible)


c       IF(NST.lt.1) AM = DM2
       BMG = BKM/AM
       ALT(J) = ALT(J+1) + BMG*TA*DZ(J+1)*1.E-5

       DALT(J) = ALT(J) - ALT(J+1)
      ENDDO

      RETURN
      END

C=============================================================================
        subroutine data_grabber(xi,wi,n)
        parameter(nrow=11)
        dimension xi(nrow,20), wi(nrow,20), n(nrow)
100     format(2x, I2)
200     format(F7.5,1x,F7.5)
c300     format (20(f7.5,1x))
c400     format (/)
c500     format (11i3)
        do i = 1,nrow
        read (66,100) n(i)
          do j=1,n(i)
          read(66,200) xi(i,j), wi(i,j)
          enddo
        enddo
C        print*, 'n='
C        print 500,n
        !print 400
        !print*,'the Abscissas(xi)'
        !print 300,((xi(i,j),j=1,n2), i=1,nrow)
        !print 400
        !print*, 'the weight factors(wi)'
        !print 300,((wi(i,j),j=1,n2), i=1,nrow)


        do i=1,nrow
        sum = 0.
           do j=1,n(i)
        !print 20,j,xi(i,j),wi(i,j)
c20      format(5x,'j=',i2,1x,2f8.5)
          sum = sum + xi(i,j)*wi(i,j)
          enddo
C        print*,'n=',n(i),' sum =', sum
        enddo

        end
