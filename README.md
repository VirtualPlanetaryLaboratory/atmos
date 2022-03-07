# README #

Atmos is a package containing two atmospheric models, along with scripts to couple them together. One of the atmospheric models calculates the profiles of chemical species, including both gaseous and aerosol phases. The second model calculated the temperature profile. Because these profiles depend on each other - kinetic reaction rates are temperature-dependent and radiative transfer is subject to radiatively active gases - we have set up atmos to run alternate the running of these two models until both models have solutions consistent with the other one. While either of these models can be run with time-dependence, most applications of these models are to find steady-state solutions for the atmosphere that would be stable over long (geological/astronomical) time periods, given constant inputs to the atmosphere.

### [TL;DR] ###
If this all is too much or if you're a visual learner, you can go to video explanations of how to run the photochemical model, how to set up photochemical model inputs, and how to read photochemical model outputs. We hope to also provide similar video tutorials on the climate model in the future.

For a video explanation on how to run the photochemical model, go here:
https://youtu.be/-9rFn_sCnmo

And for a second video, on photochemical inputs, go here:
https://youtu.be/pcpAddcsmds

Finally, for understanding photochemical model outputs, go here:
https://youtu.be/s4JKzr-6ugk




### [Atmos Jargon/Terminology] ###

This file represents a high-level "getting started guide" for Atmos
+ FAQs for frequently asked problems/issues.

"Atmos" refers to the coupled climate-photochemical code.
"CLIMA" is the climate code.
"PHOTO" is the photochemical code. (sometimes called "PHOTOCHEM" or "PHOTCHEM")

In general, PHOTO is more developed than CLIMA. PHOTO works for a wider range of planet types than CLIMA, which is
only stable for "habitable" temperature planets. So, you can use CLIMA for "Earth-like" planets but exercise extreme caution for anything else.



### [How to run the model] ###

For both models, you need to do this in three steps. First, put the input files into place. Second, compile the model. Third, run the model executables. There is a script that automates this process (RunModels.sh), which we recommend for people that are new to the model. Advanced users will also benefit greatly from a script, but might choose to write their own.

The input files are contained in template folders. For the photochemical model, these are located in: PHOTOCHEM/INPUTFILES/TEMPLATES/. Each sub-folder there contains the inputfiles you need to run the code for a particular set of default conditions. A few initial conditions input files for the climate model are in CLIMA/IO/TEMPLATES. **IMPORTANT NOTE**: Clima has fewer templates set up than Photo because it performs reliably across a narrower range of planetary conditions, and because Clima's convergence performance is less sensitive to initial conditions compared to Photo.

We strongly recommend using our packaged script to make sure all the input files are in the proper place before running. You can execute this with:
./RunModels.sh

This will present you with a choice of templates. Once you choose a template, the script will copy all the input files into the right place. Then is will ask if you want to compile - say yes if this is your first time running the model. Then it will ask if you want to run the model. After that, it will ask if you want to compile and/or run the climate model.  **IMPORTANT NOTE**: Running the Climate model from this script will not work if Clima does not contain a template for the type of atmosphere you are running. Check CLIMA/IO/TEMPLATES to be sure.

If you have the correct template for Clima, it will then bounce back and forth - asking if you want to run photochem and clima in succession until you quit the executable or tell the script "no." This is also the easiest way to couple the models together, or to run the climate model using inputs from the photochemical model.



### [Compiling and Running Photo Manually] ###

If you want to compile the photochemical model manually - for example after making changes to the fortran code or in your own script - you can do so with:
make -f Photomake

You can then manually run the model from the command line with:
./Photo.run

Once you run the model one time with RunModels.sh, you do not need to run it again. (But you can use it if you'd like!) After that first time, the input files for your chosen template will be in the correct place. So you can change the input files as needed and then re-run the model directly:
./Photo.run

If you change the fortran files, you'll need to recompile the code before re-running it:
make -f Photomake
./Photo.run

One particular thing to pay attention to is if you change any "include files" that have a *.inc suffix. The most common reason for this is that you have changed the size of one or more vector size declarations contained in PHOTOCHEM/INPUTFILES/parameters.inc. If you do change any of these *.inc files it is safest to recompile the entire photochem model. you can do this with:
make -f PhotoMake clean; make -f PhotoMake

This will erase all the model executables, and recompile everything.



### [Compiling and Running Clima Manually] ###

If you want to run the 'Climate' portion of 'Atmos' in standalone mode (without using inputs from 'Photochem') do this:
Before you start compiling and running the code, make sure your input files make sense to the problem you are addressing. The input files that are needed for the climate code part of 'Atmos' are: (1) input_clima.dat, which contains several parameters like the number of steps to run the model, pressure at the surface, pressure at the top of the atmosphere, surface temperature, surface albedo, solar constant and model switches. (2) mixing_ratios.dat - which is only needed when running Clima in standalone mode (no inputs from photochem) - contains the species mixing ratios are adjusted. Note that you only need to change these files once before the start of the run. There is a standard set of input files already included with the release. They are located in atmos/CLIMA/IO/TEMPLATES/. The three templates that are available currently are (1) Archean Earth + Haze (2) Archean Earth + Haze + Sulfur species (3) Modern Earth.

To put these models into the proper location, copy one of these template files (i.e, both input_clima.dat and mixing_ratios.dat) to the input directory 'atmos/CLIMA/IO'. Then you can change the 'input_clima.dat' and 'mixing_ratios.dat' according to your needs. We recommend testing the model first, before changing these default values.

Once the input file parameters are properly chosen, it is time to compile the code. There is a makefile in the top level directory (atmos) called 'ClimaMake'. This make fill will compile each individual subroutine. (These subroutines are located in CLIMA/CONVEC, CLIMA/COUPLE, CLIMA/PLUME, CLIMA/PRTCL, CLIMA/RADTRANS, CLIMA/RRTM, and CLIMA/SETUP). The makefile then  stitches the individual runfiles (located in CLIMA/OBJECT_CLIMA/) into a single executable file, Clima.run. That file is what you use to run the climate model.

When compiling the code for the first time, it is best to start with 'cleaning' the code using:
make -f ClimaMake clean

This command will remove the executable files, and force the model to recompile everything.

To compile the Climate code portion, simply use the command:
make –f ClimaMake

After compilation, time to run the code:
./Clima.run

Once you have compiled and run the code, you can run it again without re-compiling so long as you do not change the fortran code itself. In other words, you can change the input files, and run the model again without recompiling the code. If you do change the fortran files, then recompile before running again:
make –f ClimaMake

If for any reason the code is not compiling or running properly, you can check the main Climate code to address the issue. The main Climate code is called ClimaMain.f, try saving your input files with other names and running for default versions of the input files. If that gets rid of the error, then you had a problem in your input files. If it doesn't, the model itself is having trouble and you should get in touch with us for technical support. See the bottom of this read.me for contact information.



### [Coupling Photo and Clima] ###

If you want to manually run the coupled Clima-Photochem model, do this:
Put the input files in place - ideally with the ./RunModels.sh script. Then run Photo uncoupled (set ICOUPLE = 0 in input_photochem), and then compile and execute Photochem. Photochem will automatically generate coupling files in the /COUPLE folder for Clima to use. These files are:
fromPhoto2Clima.dat: contains altitude, pressure, O3, H2O, CH4, CO2, C2H6
hcaer.photoout.out: hydrocarbons in the atmosphere
mixing_ratios.dat: mixing ratio file for clima
coupling_params.out: time ago (billions of years), pressure at bottom of atmosphere, fractal particles (=1 if yes), type of star, haze type, number of atmospheric layers, solar constant scaling, surface gravity

Turn on ICOUPLE (=1) in CLIMA/IO/input_clima.dat. Turning on coupling means Clima will ignore certain variables in CLIMA/IO/input_clima.dat, because it now gets these from the photochemical model outputs. Clima will also ignore the mixing_ratios.dat file in CLIMA/IO, instead reading this information from /COUPLE/mixing_ratios.dat. You do not need to change PG0 (surface pressure) in input_clima when coupling is on. Clima gets this from coupling_params.out. You also do not need to change SOLCON in input_clima. It will scaled this based on TIMEGA in coupling_params. In general, any parameter it reads in from a coupling file will OVERWRITE any equivalent parameter set in one of the Clima input files.

Run Clima (./Clima.run). When it finishes, it will generate certain files for Photo:
fromClima2Photo.dat: file with altitude, temperature, water mixing ratio
hcaer.climaout.out: sanity check that clima is dealing with hydrocarbons properly. You can ignore this unless you're debugging hazes.

Next, Set ICOUPLE = 1 in input_photochem and run photo...

Run Clima.
Run Photo.
Run Clima.
Run Photo... And so on.
When do you stop? Up to whatever convergence criteria you are happy with. There are coupling wrapper scripts that may be provided upon request.



### [Other Help Sources] ###
There are other sources of help and information contained in this repository. First, PHOTOCHEM/INPUTFILES/README.TXT is a file containing tips, tricks and common errors, created originally by Giada Arney and added to by a few other folks over the years. Additionally there are some informational files in some of the data directories. For example, see PHOTOCHEM/INPUTFILES/TEMPLATES/README.txt for information on the templates. There are notes and references on most absorption cross sections inside the folder for each individual molecule inside PHOTOCHEM/DATA/XSECTIONS. For further help, save everything, go back to default input files, and the original code base. If none of that works, or if you still need help, contact us via the information below.



### [Contributors and Development History] ###
This version of the code was started by Mark Claire in an effort to merge the best features from photochemical models that had been previously developed by Kevin Zahnle, James Kasting, and David Catling. In the process, Mark generalized the model to be less Earth-specific, and "modernized" the FORTRAN coding by removing most of the vestiges of F66 code. This made it into a F77/F90 hybrid - ugly but MUCH more could be done. Much of the effort was in removing hardcoding, and in allowing flexibility in the choice of wavelength grid, stellar flux, atmospheric species, and reactions used. The generalization of the code allowed us to develop “templates” for different types of planets. These templates contain different input files to get the code to run for different planet types - without having to change the Fortran code itself. This prevents the need to maintain different “versions” of the code for these different cases. Instead, all one has to do is move these input files to the correct directory (the RunModels.sh script automates this for you), (re)compile the code, and then run it. Will Sluder continued this effort, working with Eric Hebrard and Ryan Felton to add more flexible input file formats, streamline the cross sections subroutine, and allow for the use of reaction rate constants from the KInetic Database for Astrophysics (KIDA). And Eddie Schwieterman upgraded the numerical solver to work under highly oxidized conditions. D.J. Teal has led further modernization of the code, removing redundancies in many places, cleaning up the comments in others, and addressing many (hundreds?) of compilations errors when compiling with gfortran. Andrew Lincowski increased the resolution of the photolysis grid and updated many cross sections and reaction rates. He further developed the ModernEarthComplex template. Sandra Bastelberger has led major bug fixes and model upgrades throughout most of the subroutines, with critical input and assistance from Jaime Crouse. Multiple people have also contributed to specific "templates" over time - for example Amber Young for ModernEarth and Mars, Eddie Schwieterman for high-O2 cases, Giada Arney for Archean Earth. Sandra Bastelberger has further improved and standardized many of these templates, as part of her efforts to improve the overall quality and consistency of the model.

The climate code was most recently advanced by Sandra Bastelberger and Ravi Kopparapu, who added k-coefficients derived by Eric Wolf, and implemented some suggested improvements from Ben Hayworth. They improved the model's k-coefficients, and upgraded many of the numerical approaches. Ashley Horan and Shawn Domagal-Goldman. The coupling of the two models used logic originally developed by Antigona Segura. She also introduced much of the multi-star functionality used here. Giada Arney added to the model's ability to self-consistently treat aerosol particles. And Giada Arney upgraded the treatment of aerosols in both models, as she also ensured they were treated more self-consistently between the two codes. D.J. Teal, Sandra Bastelberger, and Shawn Domagal-Goldman developed run scripts to automate input file generation, code compiling, and code convergence criteria.

The model development was supported by NASA Astrobiology Institute's Virtual Planetary Laboratory lead team, supported by NASA under cooperative agreement NNH05ZDA001C. The model development is also supported under the Sellers Exoplanets Environments Collaboration, which is part of the NASA Internal Scientist Funding Model.

Contact email - we suggest emailing all of us as our work ebbs and flows and at any given time its almost certain that at least one of us is over-committed. If you email us all, you'll have a better chance that at least one of us is not.
Giada Arney (giada.n.arney@nasa.gov)
Mark Claire (mc229@st-andrews.ac.uk)
Shawn Domagal-Goldman (shawn.goldman@nasa.gov)
Ravi Kopparapu (ravikumar.kopparapu@nasa.gov)




### [FAQs] ###

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
"NSTEPS" in PhotoMain.f). The convergence criterion for PHOTO is when the timestep
length is > the age of the universe in seconds.
For the climate model, go to CLIMA/IO/clima_allout.tab
DT, DIVF, DIVFrms should be small-valued.  How small?  Values on the order of
1E-5ish and smaller and you should be golden.  Values 1E-3...probably ok.
Any bigger?  Use at your own risk. Or don't.  Because you're probably not converged.

***

QUESTION: Ok, the model isn't converging. Help?
ANSWER: See below under the Problems, Bugs, and Debugs section...

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
   * coupling_parsms.out: needed coupling parameters passed from PHOTO to CLIMA

2. Turn on ICOUPLE (=1) in CLIMA.  Turning on coupling means CLIMA will ignore certain
   inputs in input_clima.dat because it now gets these from the photochemical model outputs,
   and it will ignore the mixing_ratios.dat in CLIMA/IO, instead
   reading this in from /COUPLE/mixing_ratios.dat. You do not need to change PG0 (surface pressure)
   in input_clima when coupling is on...it gets this from coupling_parms.out.  You also do not
   need to change SOLCON in input_clima...it will scale get this based on TIMEGA in coupling_parms.

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

### [Problems, Bugs and De-bugs] ###

Please help update this section as you find and fix problems...

***

PROBLEM: This happens at the end of a run:
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
