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

~~~~~~~~~~Atmos Basics/How-tos/Frequent Problems~~~~~~~~~~~

This file represents a high-level "getting started guide" for Atmos
+ FAQs for frequently asked problems/issues.

"Atmos" refers to the coupled climate-photochemical code.
"CLIMA" is the climate code.
"PHOTO" is the photochemical code. (sometimes called "PHOTOCHEM" or "PHOTCHEM")

In general, PHOTO is more developed than CLIMA.
PHOTO works for a wider range of planet types than CLIMA, which is
only stable for "habitable" temperature planets. So, you can use
CLIMA for "Earth-like" planets but exercise caution for anything else.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
----------------HOW TO RUN THE MODEL-----------------------

***
QUESTION: What's the quickest way to get started with the code?
ANSWER: The quickest and easiest way to initialize Atmos is by typing:
> ./RunModels.sh
This will allow you to select a template, compile PHOTO, run photo,
compile CLIMA, and run CLIMA.
When initializing a new template you should recompile PHOTO.
Only ModernEarth, Archean+haze, and ArcheanSORG+haze can properly
couple to CLIMA.


QUESTION: How do I compile these codes?
ANSWER: There are two makefilesin the atmos top level directory:
makefileClima and makefilePhotochem.
To make these two codes in the top level directory you do:
> make -f ClimaMake
> make -f PhotoMake

***

QUESTION: Ok, now how do I run these codes?
ANSWER: To run the climate model type:
> ./Clima.run
To run the photochemical model type:
> ./Photo.run

***

QUESTION: Where are the inputfiles for these codes?
ANSWER: The PHOTO inputfiles live in PHOTOCHEM/INPUTFILES
The CLIMA input files live in CLIMA/IO.

***

QUESTION: I want to edit the code.  There are too many files here.  What
ones do I care about?
ANSWER: If you want to edit the photochemical code, you care about
the file PhotoMain.f in the top-level directory.  Likewise, for the climate
model, it's ClimaMain.f.
Subroutines for Clima are in various CLIMA/ subfolders.  Subroutines for
PHOTO are easier to find because they are all in one place in
PHOTOCHEM/SUBROUTINES.

***


QUESTION: How do I know if the models have converged?
ANSWER: The photochemical model should stop in < 10000 steps if it converges
(ideally much fewer than that...the number of steps is set by the variable
"NSTEPS" in PhotoMain.f).
For the climate model, go to CLIMA/IO/clima_allout.tab
DT, DIVF, DIVFrms should be small-valued.  How small?  Values on the order of
1E-5ish and smaller and you should be golden.  Values 1E-3...probably ok.
Any bigger?  Use at your own risk. Or don't.  Because you're probably not converged.

***

QUESTION: Ok, the model isn't converging. Help?
ANSWER: See below under the Problems & Errors section...

***

QUESTION: How does this whole coupling thing work?
ANSWER:  The power of the coupled model is interactive chemistry:
the climate responds to the photochemistry, and the photochemistry, in turn,
responds to the atmosphere's temperature. Awesome, right?!

You will initially need to run ONE of the models uncoupled (ICOUPLE=0) to generate an initial
atmospheric state. In general, we run the photochemical model first, and this seems to work well.

Here are your steps for "manual" coupling:
1. Run PHOTO uncoupled (set ICOUPLE = 0 in input_photochem). PHOTO will automatically generate coupling files
   in the /COUPLE folder for clima to use (or not). These files are:
   * fromPhoto2Clima.dat: contains altitude, pressure, O3, water, Ch4, Co2, c2h6
   * hcaer.photoout.out: hydrocarbons in the atmosphere
   * mixing_ratios.dat: mixing ratio file for clima
   * output_couple.dat: needed coupling parameters passed from PHOTO to CLIMA

2. Turn on ICOUPLE (=1) in CLIMA.  Turning on coupling means CLIMA will ignore certain
   inputs in input_clima.dat because it now gets these from the photochemical model outputs,
   and it will ignore the mixing_ratios.dat in CLIMA/IO, instead
   reading this in from /COUPLE/mixing_ratios.dat. You do not need to change PG0 (surface pressure)
   in input_clima when coupling is on...it gets this from time_frak_photo.out.  You also do not
   need to change SOLCON in input_clima...it will scale get this based on TIMEGA in time_frak_photo.

3. Run CLIMA. When it finishes, it will generate certain files for PHOTO:
   * fromClima2Photo.dat: file with altitude, temperature, water mixing ratio
   * hcaer.climaout.out: sanity check that clima is dealing with hydrocarbons properly.
                         you can ignore this.

4. Set ICOUPLE = 1 in input_photochem and run PHOTO...

5. Run CLIMA...

6. Run PHOTO...

   When do you stop? Up to whatever convergence criteria you are happy with.

Think manual coupling sounds like a bore?  We agree. There are numerous coupling scripts
that can be provided on request.


***


----------------PROBLEMS & ERRORS :(-----------------------
Please help update this section as you find and fix problems...

***

PROBLEM: This happens at the end of a CLIMA run:
"IEEE_INVALID_FLAG IEEE_DIVIDE_BY_ZERO IEEE_UNDERFLOW_FLAG IEEE_DENORMAL"

SOLUTION: This happens when you compile the model using a newer version of the
gfortran compiler.  It means, to quote Shawn: "somewhere in a ginormous matrix
there are some tiny numbers."

It's probably not something you need to worry about...

***

PROBLEM: You ran RunModels.sh to replace your input files with a different template, but
cannot run without re-compiling code first... You get an error regarding reactions.rx:

Would you like to compile PHOTOCHEM (y/n)?:n

Would you like to run PHOTOCHEM (y/n)?:y
At line 783 of file PhotoMain.f (unit = 9, file = 'PHOTOCHEM/INPUTFILES/reactions.rx')
Fortran runtime error: End of file

Error termination. Backtrace:
#0  0x12081e729
#1  0x12081f3f5
#2  0x12081fb59
#3  0x1208e7f8b
#4  0x1208e8527
#5  0x1208e55c3
#6  0x1208ea5b4
#7  0x1208e7b19
#8  0x109fb5482
#9  0x10a011ade

SOLUTION: Since some of the output files we're writing have dynamic format
statements that depend on the particular template's parameters.inc values,
always re-compile (select y) when using a different template! (12/16 MAB)

***

PROBLEM: The photochemical and/or climate model is taking FOREVER to run and/or
not converging EVER!

SOLUTION: Well...there might not be a solution.  It might really need to take that long.
However, here are some reasons why it might be taking forever:

Photochemical Model:

- Are you using an appropriate in.dist file?  in.dist lives in the PHOTOCHEM folder and
  contains the initial conditions for your photochemical run. If you perturb your input
  files too much compared to in.dist, convergence may never happen.  You may need to
  slowly step from one part of parameter space to another...(which is annoying, but sometimes
  necessary).  To do this, get an in.dist that works for a starting part of parameter space (you can do
  this by starting from one of the templates).
  How do you know you're converging? If the photochemical model converges in less than 10000 steps
  (hopefully a lot less!)
  Perturb a quantity in the input files a little bit -> run the model to convergence -> copy
  the out.dist into in.dist -> Perturb the quantity some more -> run the model -> out.dist
  to in.dist -> again and again until you reach the target region of parameter space.

- You can now plot from PHOTOCHEM/PTZ_mixingratios_in.dat & PHOTOCHEM/OUTPUT/PTZ_mixingratios_out.dat,
  where the pressure, temperature, altitude and mixing ratios (see headers) of all long-lived species,
  in your most recently ran template, are tabulated in a human-readable style. (Added by MAB, 12/16)

Climate Model:
- What's NSTEPS in CLIMA/IO/input_clima.dat set to?  Is it unreasonably high?
  (~hundreds of steps is probably reasonable)

- In input_clima.dat, what is IMET set to?  If it's 1, methane is turned on.  If you
  don't have a lot of methane in your atmosphere, turn this off and things will be a lot
  faster. Clima calculates gas absorption in a series of nested for loops so the run time
  does not scale linearly with the number of gases in the atmosphere.

- In input_clima.dat, what is IMETETH set to?  If it's 1, methane and ethane are turned on.  If you
  don't have a lot of methane and ethane in your atmosphere, turn this off and things will be a LOOOOT
  faster. Clima calculates gas absorption in a series of nested for loops so the run time
  does not scale linearly with the number of gases in the atmosphere.

- Did you run Clima for an obscene number of steps and it's still not converged? Try
  copying CLIMA/IO/TempOut.dat to TempIn.dat.  Then edit input_clima.dat and change
  IUP to 0. Re-run clima, and it will use the previous run's ending state as its initial
  conditions. This will hopefully help it achieve convergence. You can also try turning off
  ICONSERV in input_clima.dat. That sometimes helps.

- Are you running CLIMA for a planet that is in a "non-habitable" temperature regime?
  The model typically cannot converge for very hot or very cold planets.
