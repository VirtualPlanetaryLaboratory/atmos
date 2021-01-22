Converged templates can be found in the TEMPLATES folder.
Some informmation on each template can be found in the README
in the TEMPLATES folder. 

Variable definitions of the input files are as follows: 
----------------------------------------------
----------------PLANET.dat--------------------
This has basic info on the planet-star configuration you're dealing with.

* G: gravitational accelleration on your planet [cm/s]
* FSCALE: solar constant scaling. 1 = modern Earth
    Note that this parameter should be used to scale
    your solar flux based on the planet's DISTANCE
    (IMPORTANT CAVEAT: if you wish to adjust the solar constant
    based on the sun's age -- and this applies ONLY to the sun, not other stars --
    it is best to use the TIMEGA variable described below)   
* ALB: Planet's surface albedo.
* ZTROP: Height of tropopause. [cm]
* FAR: Argon partial pressure 
* R0: planet's radius in cm
* P0: surface pressure in atm
* PLANET: string naming the planet. "EARTH" or "MARS" are possible
* TIMEGA: Time in billions of years ago for scaling the solar constant
     based on the time you're simulating based on Claire et al. 2012
     IMPORTANT CAVEAT: This applies ONLY to the sun. 
* IRESET: Depreciated? 
* msun: Choice of star you want to run (apologies for the convoluted legacy numbering system)
choices for msun (these numbers are from photgrid.f)
      pstar = 12    !Gj 581 from Lucianne Walkowicz  
      pstar = 13    !high resolution solar data from ATLAS1/3 (Thullier et al 2004) - best choice for the sun 
      pstar = 14    !Kevin Zahnle's data from photo.dat
      pstar = 15    !AD Leo from VPL climate DB
      pstar = 16    !AD LEO from VPL website
      pstar = 17    !T3200.dat
      pstar = 18    !K2V
      pstar = 19    !F2V
      pstar = 76    !GJ876
      pstar = 20    !M8V
      pstar = 21    !M5V/Proxima Centauri
      pstar = 22    !G0V star for wasp12b (first Hot Jupiter model)
      pstar = 23    !Phoenix 3300K star (for GCM Earth Project)
      pstar = 24    !HD40307/K2.5V
      pstar = 25    !HD85512/K6V
      pstar = 26    !HD97658/K1V
      pstar = 80    !GJ876 MUSCLES data adaptive constant res panchromatic
      pstar = 81    !GJ551 """""" (Proxima Centauri)
      pstar = 82    !GJ581 """"""
      pstar = 83    !GJ667c """""" (note: GJ667 is a triple star system)
      pstar = 84    !GJ1214 """"""
* IHZSCALE: depreciated
* uvscale: apply constant UV scaling factor

---------------------------------------------------
---------------input_photochem.dat-----------------    
This file has the basic info on certain photochem model parameters 
you might want to set and also can turn various model functionalities on or off.

C     AGL = DIURNAL AVERAGING FACTOR FOR PHOTORATES
C     ISEASON = TELLS WHETHER P AND T VARY WITH TIME (THEY DON'T FOR
C               ISEASON < 3)
C     IZYO2 = TELLS WHETHER SOLAR ZENITH ANGLE VARIES WITH TIME (0 SAYS
C             IT DOESN'T; 1 SAYS IT DOES)
C     LGRID = 0  USE OLD JPL WAVELENGTH GRID
C     LGRID = 1  USE MARK'S NEW HIGH RESOLUTION GRID (SEE ALSO GRID.F)
C                IF LGRID=1, IO2 and INO should be set to 2
C     IO2 = 0 FOR ALLEN AND FREDERICK O2 SCHUMANN-RUNGE COEFFICIENTS
C         = 1 FOR EXPONENTIAL SUM FITS (FOR LOW-O2 ATMOSPHERES)
C         = 2 FOR USING HIGH RESOLUTION CROSS SECTION
C     INO = 0 FOR ALLEN AND FREDERICK NO PREDISSOCIATION COEFFICIENTS
C         = 1 FOR MODIFIED CIESLIK AND NICOLET FORMULATION
C         = 2 USE OUT.NOPRATES FILE
C     EPSJ = AMOUNT BY WHICH TO PERTURB THINGS FOR JACOBIAN CALCULATION
C     PRONO  = COLUMN-INTEGRATED NO PRODUCTION RATE FROM LIGHTNING IN
C              EARTH'S PRESENT ATMOSPHERE
C     FRAK = 0  use Mie scattering for hydrocarbon particles
C          = 1  use Fractal scattering for hydrocarbon particles
C     HCDENS = 0.63 g/cm^3 (Archean hydrocarbon density - Trainer et al 2006)
C            = 0.8 g/cm^3 (titan tholins - Trainer et al 2006)
C 	     = 1.0 g/cm^3 (old suspicious default value)
C     ICOUPLE = 0 to turn off coupling / 1 to turn on coupling
C             = 0 to use old formatting (E7.2 for VDEP and FIXEDMR)
C     MONSIZE = 0 for 0.05 um monomers for fractal particles
C      	      = 1 for 0.01 um monomers for fractal particles
C	      = 2 for 0.02 um monomers for fractal particles
C	      = 3 for 0.07 um monomers for fractal particles
C	      = 4 for 0.10 um monomers for fractal particles
C     ZY      = solar zenith angle
C     USOLMIN = Value under which the time stepper no longer considers
C               a species when calculating "max relative change."
C               This should be set to 1.e-20 or lower, unless you are
C               moving the model to a new regime or creating a new
C               template and require looser convergence criteria.
C               Even in that case, you should lower/tighten this again
C               before production runs.
C     KIDA = 1 if using "kida-style" reactions.rx format
C          = 0 if using the code's "old-style" reactions.rx format


-------------------------------------------------
-------------------species.dat-------------------
This is the boudnary conditions file tht tells the model what's
in the atmosphere and what fluxes are entering/leaving the atmosphere
at the top and bottom.


LBOUND: type of lower boundary condition
* 0 = constant deposition velocity (VDEP) [cm/s]
* 1 = constant mixing ratio (FIXEDMR)
* 2 = constant upward flux (SGFLUX) [molecules/cm**2/s]
* 3 = constant vdep + vertically distributed upward flux  (uses both SGFLUX and DISTH)
-- If you set 1 (constant mixing ratio) the model will calculate the flux necessary
   at the surface to sustain that mixing ratio.  That flux will be printed in out.out and
   can be found for searching for the string "FLUXES OF LONG-LIVED SPECIES"
-- If you set the surface LBOUND to a type of flux, the model will, conversely,
   calculate the mixing ratio that results from that flux. To know what this is,
   search in out.out for the last time " MIXING RATIOS OF LONG-LIVED SPECIES" appears.

*MBOUND - Upper boundary conditions
* 0 = CONSTANT EFFUSION VELOCITY (VEFF)  - (H and H2 set in code for molecular diffusion/diffusion limited flux)
* 1 = constant mixing ratio - never been used so would need testing
* 2 = CONSTANT FLUX (SMFLUX) (option for CO2/CO/0 in code)

Notes:
The model treats long lived, short lived, etc species differently.
You cannot move species between long lived and short lived, etc, without
generating a new in.dist file (inidial conditions) 

-------------------------------------------------
-----------------parameters.inc------------------
This file is actually part of the code and so, if you edit it, you will need to
recompile the code. 


NZ: number of layers in the model atmosphere
NQ: number of species you need to solve a diffusion equation for. This
    is the number of long lived (LL) species in species.dat
NP: number of types of particles. (=1 for SO4 aerosols, =2 for SO4 and S8
    = 3 for sulfur aerosols and C4H2 (HCAER) hydrocarbons, 
    = 4 for all of the above and C6H4 (HCAER2) hydrocarbons
    In general, set NP to 4 if you're doing organic hazes.
NR: number of reactions (see reactions.rx)
NSP: total number of chemical species. NSP1 = NSP + 1 (includes HV)
     and NSP2 = NSP +2 (includes M and HV)
NMAX: Max # of reactions that a single species participates in. This value
      is used in the code to set the size of the matrix of loss processes
      and the matrix of production processes
KW: Max number of wl intervals. Used to set the size of lots of arrays
    in the model.
KJ: Number of photolysis reactions in reactions.rx
KS: Number of species in photolysis reactions
NAQ: Number of aqueous phase species (see aqueous.f)
     I do not think you can edit this parameter without editing aqueous.f
NT: number of temperatures in sulfate/H2O vapor pressure file
    Refers to DATA/aerosol.table, which is a table of sulfuric acid
    and water vapor pressures as a function of temperature and H2SO4
    concentration in the particles.
NF: Number of pressures in DATA/aerosol.table
ML: Honestly not sure...but it's used in Mscat.f. My guess is it's the 
    number of moments computed for the scattering phase function
KIN: File unit for cross sections.  Probably no need to change this?

-------------------------------------------------
-----------------reactions.rx--------------------
All of the reactions in the model. "WEIRD" reactions are hardcoded
in Rates.f
