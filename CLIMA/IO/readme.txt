This subdirectory contains the input and output files for the 
climate model.

The initial files are the input files to run Earth's atmposphere 
around the Sun.	

** Input files

- mixing_ratios.dat
It has the composition of the planetary atmosphere 

- input_clima.dat
Variable input parameters needed to set initialize the program

- TempIn.dat
Initial temperature and water profiles. Before you run the program for 
the first time this file is identical to TempOut.dat

- fromPhoto2Clima.dat
Water and ozone profiles from photochemical model. 
It is not necesary if you are running the climate code alone.
The file that is in this subdirectory, before running the model for the 
first time in coupled mode corresponds to a coupled converged run for 
present Earth around the Sun.


** Output files 
NOTE:Before running the program for the first time only TempOut.dat
will be present.

- OUTPUT_RRTM
Output file from the RRTM code contains the IR fluxes and heating rate 
calculated by this program 

- SolarHeating.tab
Heating rate in each wavelength interval from the solar subroutine
(visible/near IR  wavelengths).

- TempOut.dat 
Temperature and water profiles from the last solution obtained by the code. 
It has the same structure as TempIn.dat To restart from the last solution 
copy this file to TempIn.dat
Before you run the program for the first time this file is identical to
TempIn.dat

- clima_last.tab 
Profiles from the last solution (from right to left):
Altitude (km)
Pressure (atm)
Temperature (K)
Water mixing ratio
Ozone mixing ratio
heating rate (ergs/day)
cooling rate (ergs/day)

- clima_allout.tab
This file is very important to check the convergence of the program.
It contains the fluxes, temperature, heating rates and water profiles 
calculated by the program form the first to the last step. 
The first 3 and the last 3 runs are written in detail. 
J = Number of the layer (1 is the top of the atmosphere)
P = Pressure (atm)
ALT = Altitude (km)
T = Temperature calculated at the end of this step (K)
CONV = Tags for convective layers (non cero are convective layers)
DT = T-TOLD
TOLD = Temperature from the previous iteration
FH2O = Water mixing ratio
TCOOL = Cooling rate (ergs/day)
THEAT = Heating rate (ergs/day)
PF = Pressure at the 
FTOTAL = Total flux (ergs/cm^2/s)
FTIR = Total IR flux
FDNIR = Downward IR flux
FUPIR = Upward IR flux 
FTSOL = Total solar flux (visible/nearIR)
FDNSOL = Downward solar flux
FUPSOL = Upward solar flux 
DIVF = total flux / IR flux at the top of the atmosphere
 
For the rest of the runs the next quantities are listed:
NST = iteration number
JCONV = Number of the last convective layer
CHG = Maximum relative change between the temperature caltulated at the step
      NST and the step NST-1. It sets the time step.
dt0 = time step in seconds
DIVF(1) = total flux / IR flux at the top of the atmosphere. Smaller the number
         better the convergence. DIVF(1) < 1e-2 is acceptable, <1e-3 is 
         very good.
DT(1) = Difference between the temperatures at the rop of the atmosphere
        calculated at the step NST-1 and NST. Smaller the number, better the 
        solution.
TN(ND)= Surface temperature calculated at the end of the step NST


- fromClima2Photo.dat 
Temperature and water profiles for the photochemical model