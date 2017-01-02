Running the photochemical model for a high-O2 (95% O2) Earth with the Sun remaining as the host star.

Assumption is a post-runaway abiotic O2 atmosphere (e.g., Luger & Barnes 2015), but recognizing this can
only happen inside the inner edge of the habitable zone for a sun-type star. Mantle/surface oxidized. 
No reduced gases should be present. 
In this "dry" case, almost all the hydrogen has escaped, so the total water content is extremely low. 

Updated "species.dat" for oxidized mantle case:

H2 --> LBound 1--> 0 (Vdep = 2.4E-02)
CO --> Lbound 1--> 0 (Vdep = 1.2E-02)
CH4 --> Lbound still 1; mixing ratio fixed at 0.0 (Vdep changed to 1.e00)
N2O --> Mixing ratio goes from 3.7e-07 --> 0.0
N2 --> Rmix 0.78 --> 0.02

Other reduced species have a large deposition velocity because they would be consumed by an oxidized surface. 
Oxidized species have a 0 deposition velocity. 

in.dist close to convergence.

-EWS - updated 09/21/2015 (eschwiet@uw.edu)

03/04/2016: Updated in.dist to correct for earlier code error - EWS (eschwiet@uw.edu)
03/14/2016: Updated species.dat, parameters.inc, and in.dist to correct for NR (no significant effect),
            and redox calc error from Tri-diag species. 
