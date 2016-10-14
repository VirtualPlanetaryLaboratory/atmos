# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

* Quick summary
This is a coupled photochemistry-climate code. This model leverages
the work by the Kasting, Zahnle, and Catling groups and represents an
an effort to merge the various versions and features of the codes
that have been developed over the years by the groups using the code.

### Contributors ###

The coupling routines used here were developed by Antigona Segura.
She also introduced much of the multi-star functionality used here.

The photochemical code was significantly advanced by Mark Claire,
who combined features, generalized the model to be less Earth-
specific, and modernized the FORTRAN coding.

The climate code was updated by Ravi Kopparapu and
Ramses Ramirez, working for James Kasting. They improved the model's water-vapor and CO2 absorption coefficients, and upgraded many of the numerical approaches.

Eddie Schwieterman improved the numerical solver so the code would run for a wider range (more oxidizing) redox conditions. He also added templates for very high O2 conditions that could be caused by “Luger-Barnes” atmospheres where H is stripped away by high-energy radiation, leaving behind O. Along with Giada Arney, he developed tools to couple these codes to SMART, a line-by-line radiative transfer tool that we can use to predict the transit and reflected light spectra of the atmospheres we simulate.

Giada Arney added fractal hazes to the climate code and ensured coupling between photochem and clima, including self-cosnistent treatment of hazes. She also added (or repaired) all of the Archean Earth templates.

Eric Hébrard, Will Sluder, and Ryan Felton generalized the subroutine that calculated absorption crosss-sections (Xsections.f) and modernized many of the reaction rate constants. Ryan Felton also added a Titan template to the model (not public yet), under the guidance of Eric Hébrard.

Ashley Horan and Shawn Domagal-Goldman, using the past work of Antigona Segura as a guide, coupled the modern versions of photochem and clima together.

Dillon Teal developed run scripts to automate input file generation, code compiling, and code convergence criteria.

Mahmuda Afrin Badhan added a hot Jupiters template, under the guidance of and based on the prior work of Ravi Kopparapu.

Amber Britt fixed the modernEarth+Chlorine template, and incorporated the Mars template that was originally developed by Mark Claire and Meg Smith.

Dillon Teal, Will Sluder, Amber Britt, and Shawn Domagal-Goldman did extensive edits to the code formatting to standardize things and make everything more readable.


* Version

* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Summary of set up

* Configuration

* Dependencies

* Database configuration

* How to run tests

* Deployment instructions

### Contribution guidelines ###

* Writing tests

* Code review

* Other guidelines

### Who do I talk to? ###

* Repo owner or admin

* Other community or team contact