Running the photochemical model for a high-O2 (95% O2) Earth with the Sun remaining as the host star.

Assumption is a post-runaway abiotic O2 atmosphere (e.g., Luger & Barnes 2015), but recognizing this can
only happen inside the inner edge of the habitable zone for a sun-type star. Mantle/surface oxidized. 
No reduced gases should be present. 

Changes to the species.dat file from the "ModernEarthThatWorks" template:

H2 --> LBound 1--> 0 (Vdep remains 2.4E-04)
CO --> Lbound 1--> 0 (Vdep remaints at 2.4E-04)
CH4 --> Lbound still 1; mixing ratio fixed at 0.0 (Vdep changed to 1.e00)
N2O --> Mixing ratio goes from 3.7e-07 --> 0.0
N2 --> Rmix 0.78 --> 0.04

in.dist close to convergence. Note H2O still present. 

-EWS 07/29/2015 (eschwiet@uw.edu)
