Running the photochemical model for a high-O2 (95% O2) Earth with GJ 876 as the host star.

Assumption is a post-runaway abiotic O2 atmosphere (e.g., Luger & Barnes 2015), with the planet residing in the CHZ, but
experiencing a runaway during the pre-main sequence super-luminous phase of the star. Mantle/surface oxidized. 
No reduced gases should be present. 

Changes to the species.dat file from the "ModernEarthThatWorks" template:

H2 --> LBound 1--> 0 (Vdep remains 2.4E-04)
CO --> Lbound 1--> 0 (Vdep remaints at 2.4E-04)
CH4 --> Lbound still 1; mixing ratio fixed at 0.0 (Vdep changed to 1.e00)
N2O --> Mixing ratio goes from 3.7e-07 --> 0.0
N2 --> Rmix 0.78 --> 0.02

in.dist close to convergence. Note H2O still present. 

-EWS 09/08/2015 (eschwiet@uw.edu)
