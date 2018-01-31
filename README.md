# README #

Atmos is a package containing two atmospheric models, along with scripts to couple them together. One of the atmospheric models calculates the profiles of chemical species, including both gaseous and aerosol phases. The second model calculated the temperature profile. Because these profiles depend on each other - kinetic reaction rates are temperature-dependent and radiative transfer is subject to radiatively active gases - we have set up atmos to run alternate the running of these two models until both models have solutions consistent with the other one. While either of these models can be run with time-dependence, most applications of these models are to find steady-state solutions for the atmosphere that would be stable over long (geological/astronomical) time periods, given constant inputs to the atmosphere.

#### What is this repository for? ####

This is a coupled photochemistry-climate code. This model leverages
the work by the Kasting, Zahnle, and Catling groups and represents an
an effort to merge the various versions and features of the codes
that have been developed over the years by the groups using the code. Ashley Horan and Shawn Domagal-Goldman, using the past work of Antigona Segura as a guide, coupled the modern versions of photochem and clima together.

The model development was supported by NASA Astrobiology Institute's  Virtual Planetary Laboratory lead team, supported by NASA under cooperative agreement NNH05ZDA001C.

### [Contributors and Development History](https://bitbucket.org/ravikopparapu/contributors/wiki/Contributors) ###

### [How to run Atmos ?](https://bitbucket.org/ravikopparapu/how-to-run-atmos/wiki/How%20to%20Run%20Atmos) ###

### [What does the converged output look like?](https://bitbucket.org/ravikopparapu/output-examples/wiki/Converged%20Output%20Examples) ###