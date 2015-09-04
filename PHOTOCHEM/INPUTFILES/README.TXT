Hi!  You are probably wondering what files to start with and what files are 
necessary as inputs to the model.  Read on...


~~~~~~~~~~~~TEMPLATES FOLDER~~~~~~~~~~~~~~~~~~
I suggest starting with the files in:
 /templates/Archean+haze 
These were designed to work with no edits necessary with this version of 
the model. Copy everything except in.dist to:
 /PHOTOCHEM/INPUTFILES   
Copy in.dist to:
 /PHOTOCHEM   
If the model does not converge in less than ~ 100 steps, check to make sure
you correctly copied in.dist to the PHOTOCHEM directory. If you still have 
problems, talk to Giada. The other templates SHOULD work, too, but 
your mileage may vary....

----------------------------------------------
----------------planet.dat--------------------
This has basic info on the planet-star configuration you're dealing with.

* G: gravitational accelleration on your planet [cm/s]
* FSCALE: solar constant scaling. 1 = modern Earth
    Note that this parameter should be used to scale
    your solar flux based on the planet's DISTANCE
    but typically NOT for the sun's age. That is done
    with the TIMEGA parameter described below.
* ALB: Planet's surface albedo.
* ZTROP: Height of tropopause.
* FAR: Argon partial pressure 
* R0: planet's radius in cm
* P0: surface pressure in atm
* PLANET: string naming the planet. "EARTH" or "MARS" are possible
* TIMEGA: Time in billions of years ago for scaling the solar constant
     based on the time you're simulating. This includes Mark C.'s 
     wavelength-dependent solar activity scaling
* IRESET: Should reset all parameters to defaults regardless
     of how they are set in the inputfiles (I think...)
* msun: Choice of star you want to run (sorry for the terrible
     numbering system below...partly due to legacy code)
      msun = 13    !high resolution solar data from ATLAS1/3 
      msun = 14    !kevin's data from photo.dat
      msun = 15    !AD Leo from VPL climate DB (broken...)
      msun = 16    !AD LEO from VPL website
      msun = 17    !T3200.dat (M dwarf with no chromosphere)
      msun = 18    !K2V
      msun = 19    !F2V
      msun = 76    !GJ876

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
C                NOTE: the model will output coupling files for clima regardless
C		 of whether this is on or off!  This was done intentionally because
C		 you need to use photochem to generate an initial atmospheric state
C		 uncoupled from clima.
C     NEWSPEC = 1 to use species.dat that has E8.2 format for VDEP and FIXEDMR
C             = 0 to use old formatting (E7.2 for VDEP and FIXEDMR)
C     MONSIZE = 0 for 0.05 um monomers for fractal particles
C      	      = 1 for 0.01 um monomers for fractal particles
C	      = 2 for 0.02 um monomers for fractal particles
C	      = 3 for 0.07 um monomers for fractal particles
C	      = 4 for 0.10 um monomers for fractal particles

-------------------------------------------------
-------------------species.dat-------------------
This file tells the model what's in the atmosphere and what fluxes
are entering/leaving the atmosphere at the top and bottom.

There are two possible ways to format this file depending on what value you set
NEWSPEC to in input_photochem.dat. 

NEWSPEC = 1 means you have an extra sig fig for some of the columns:
format(30X,I1,5X,2(E8.2,1X),E9.3,1X,E7.1,1X,I1,6X,2(E7.1,1X))

NEWSPEC = 0 means you don't have that extra sig fig:
format(30X,I1,5X,2(E7.1,1X),E9.3,1X,E7.1,1X,I1,6X,2(E7.1,1X))

LBOUND: type of lower boundary condition
* 0 = constant deposition velocity (VDEP)
* 1 = constant mixing ratio (FIXEDMR)
* 2 = constant upward flux (SGFLUX)
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
generating a new in.dist file (ugh...). (in.dist has model initial conditions).
You also cannot just add species to this file without also editing 
in.dist, reactions.rx, and parameters.inc.

!!!NOTE!!!: There's an annoying bit of hardcoding in Sedmnt.f that assumes
particles are the last entries in the long lived (LL) species.

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
All the reactions in your model!  Adding to this or removing from it 
means you need to also edit parameters.inc 
Format of reactions is:
first column + second column ---> third column + fourth column + fifth column

NOTES ON THE OTHER COLUMNS: 
2BODY: TWO BODY REACTIONS: Arhenius and Temperature factor
 A(J,I) = ARH * EXP (TFAC/T(I))   !two body reaction rates

3BODY: THREE BODY REACTIONS: K0, KINF, T0exp, Tinfexp

WEIRD: reactions are hardcoded into rates.f

PHOTO: photolysis reactions

-------------------------------------------------
-----------------params.dat----------------------
multipliers for o2 and ch4 concentrations.
They are read into the model but not actually used anywhere in the code...
Mark was doing some kind of test, I guess.



~~~~~~~~~~~~~~Note on the DEFAULT* and INITIAL* files~~~~~~~~~~~~
These were generated to work with my coupling wrapper script are not designed
to work "out of the box" because they're meant to be externally tweaked by 
the  wrapper script based on user inputs. Use them if you want, but know that
they may require (minor) editing before they work.  In particular, if you use
these, pay attention to the parameter "NEWSPEC" in input_photochem.  This is
set to 0 by default (old formatting of species.dat), but we typically now
use the new formatting (NEWSPEC = 1).  This is the biggest trip up people
have had using these DEFAULT* and INITIAL* files.  Also, just keep in mind 
these files were generated for hazy Archean atmospheres except the ones with
the OXY_* prefix, so if you aren't interested in that, they will require
editing.


Have fun!
- Giada