PRO run_coupled_model, spherical=spherical, dirpath=dirpath, numclima=numclima, doskips=doskips,  retainindist = retainindist,  climatwice=climatwice 

;-------------------------------------
;-------WHAT IT DOES:-----------------
;This IDL code runs the coupled climate-photochemical model for a variety
;of conditions set under the "EDIT THIS STUFF" label found below.
;It considers the model converged when change in temperature between runs
;is less than 1 K.  If it takes > 10 iterations to get to that,
;it will give up and move to the next run. It also checks for whether
;the model is producing NaNs or oscillating without converging. 
;
;Note that you need to define a file in:
;PHOTOCHEM/default_indist/DEFAULT_in.dist
;This file should contain a reasonable starting point for your simulations.
;--------------------------------------
;--------HOW TO RUN IT:----------------
;0. If you're a non-IDL person, here's a bit of help.
;   Put this program in a folder where IDL can see
;   it (typically something like ~/idl)
;   Then open idl on the command line by typing "idl" (without the quotes)
;   You compile a program in idl by typing:
; IDL> .r run_coupled_model
; Each time you edit this program, you will need to recompile. 
; Then, to run it:
; IDL> run_coupled_model 
; Although, see the note below about turning on optional keyword flags
; 
;
;1. Edit the items under the "Edit this stuff" label below to 
;   choose what kind of photochemical-climate runs you want to do.
;   You can change these parameters there:
;  pressures (atmos pressure in bars)
;  oxygen    (oxygen mixing ratio)
;  timeago   (time ago to run model in billions of years [ga])
;  dens      (density of haze particle material in g/cm^3)
;  methane   (methane mixing ratio)
;  co2       (CO2 mixing ratio)
;  co        (CO flux)
;  startype  (index of star you're using)
;  monsize   (size of haze monomers) 
;  iconserv  (iconserv = 1 or 2)
;  imw       (type of humidity parameterization) 
;  fscale    (solar flux scaling)
;  Don't forget to recompile after editing the code!
;
;2. Edit dirpathdefault as needed below (can also change working
;directory from program call or it'll prompt you)
;CHANGE THIS TO YOUR DIRECTORY!
  dirpathdefault = '/astro/users/giada/atmos/' ;top level directory of the coupled model
;
;
;3. Choose optional keyword flags.  They are:
;   spherical - True/False. set flag to make hydrocarbon aerosols spheres (fractals set by
;               default)
;   dirpath - String. set as string of directory to work in if you
;             want a different path from the default (should be TOP LEVEL directory
;             of atmos) e.g. '/astro/users/giada/atmos/'
;   numclima - String (yes, it's a string!). Set number of steps clima goes through
;              It will do 400 by default, which takes forever, but
;              a huge number may be necessary if you find the model
;              is having a hard time converging.
;   doskips - True/False. set flag to skip certain combinations in the "edit this stuff" 
;             arrays below.  Set which combinations to skip under the section
;             labeled "!!!!!!!!skips!!!!!!!!" down below
;   retainindist - True/False. by default, it restarts with a default in.dist file
;                  it copies over at the start of each iteration
;                  through the nested for loops.  If you want to use
;                  the previous loop iteration's in.dist, turn
;                  this on. I recommend turning this on! It helps
;                  convergence! 
;   climatwice - True/False. run clima twice each time. This can help convergence
;                sometimes... though you could also just have clima
;                run double the number of steps.
;   
; You can turn on keyword flags that are true/false by typing
; "/keywordflag". You can set values for keyword flags that need them
; by doing "keywordflag = value" 

;EXAMPLES:

; e.g. run it with no optional keyword flags:
; IDL> run_coupled

; e.g. run it with spherical particles, 1000 clima steps, retain
; in.dist
; IDL> run_coupled, /spherical, numclima = '1000', /retainindist 
;---------------------------------------------------------------

;------prompt user for input-------------

  IF not keyword_set(dirpath) THEN BEGIN
     print, ''
     print, '-----'
     print, 'ATMOS directory path not set.'
     print, "I'm going to use this default one unless I'm told otherwise: "+dirpathdefault
     print, 'do you want to type in a new directory path?'
     print, ''
     inputyn = ''
     print, 'Type "y" to provide a new directory manually, "n" to keep it as the default, "here" to use the directory IDL is currently open in: '
     READ, inputyn, prompt = 'Type "y/n/here" (all lowercase): '
     IF inputyn EQ 'y' THEN BEGIN
        print, 'what directory do you want to use instead?'
        print, 'example: /astro/users/giada/FOLDER/'
        inputdir = ''
        READ, inputdir, prompt = 'Type new directory path here: '
        print, 'you typed: ' + inputdir
        dirpath = inputdir
     ENDIF
     IF inputyn EQ 'n' THEN dirpath = dirpathdefault
     IF inputyn EQ 'here' THEN BEGIN
        print, 'ok, I will use the current working directory that IDL is open in'
        CD, current=TheDirectory
        print, 'the current directory is : ' + thedirectory
        dirpath = thedirectory+"/"
     ENDIF
     
     IF inputyn NE 'y' and inputyn NE 'n' and inputyn NE 'here' THEN BEGIN
        print, 'Invalid input.  You need to type either "y" or "n" or "here" (lowercase, no quotes when you enter it)'
        stop
     ENDIF

                                ;check to make sure the directory is formatted correctly:
     lastchar = strmid(dirpath, 0,1,/reverse_offset)
     firstpart = strmid(dirpath, 0,1)
     IF lastchar NE '/' THEN dirpath = dirpath+'/'
     IF firstpart EQ '/' THEN dirpath = '/'+dirpath

  ENDIF
;----------!!!!!!EDIT THIS STUFF!!!!!!!!!-----------------
;Because of the way this code edits the input files, these all need to
;be STRINGS!  This is because the photochemical and climate input
;files are finicky about formatting and it's easiest to capture this
;with string inputs. Please format all of these values with the same
;style and number of digits as the examples. 

  IMW = '2'
; CLIMA -> Integer (disguised as a string). type of humidity parameterization:
;      0 FOR SATURATED TROPOSPHERE, 1 FOR MANABE/WETHERALD
;      FIXED RELATIVE HUMIDITY, 2 FOR M/W WITH CONSTANT
;      STRATOSPHERIC H2O CONTENT AND EMPIRICAL TROPOSPHERIC
;      H2O, 3 for DRY ADIABAT
;      e.g. IMW = '2'

  ICONSERV = '0'
; CLIMA -> Integer (disguised as a string). type of energy conservation:
;     O = Non strict time-stepping method (faster)
;     1 = Each time step conservs energy (better for high CO2)  
;     e.g. ICONSERV = '0'

  pressures =  ['1.000']
; CLIMA & PHOTO -> pressure at the bottom of the atmosphere (bars)
;     This is an array of strings. e.g. pressures=['1.000', '0.500']


  oxygen = ['1.00E-08']
; CLIMA & PHOTO -> mixing ratio of oxygen
;     Array of strings. e.g. oxygen = ['1.00E-08']

  times = ['2.7']
; CLIMA & PHOTO -> time ago for scaling solar flux (billions of years)
;     Array of strings. e.g. times = ['3.8', '2.7', '0.0']


  dens = ['0.63']
; PHOTO -> density of hydrocarbon haze particle material 
;      Array of strings
;      = 0.63 g/cm^3 (Archean hydrocarbon density - Trainer et al 2006)
;      = 0.8 g/cm^3 (titan tholins - Trainer et al 2006)
;      = 1.0 g/cm^3 (old suspicious default value)

  monsize = ['0'] 
; PHOTO & CLIMA -> which sized monomers to use for fractal haze
; particles
;      Array of strings
;      0 = 0.05um; 1 = 0.01um; 2 = 0.02um; 3 = 0.07um; 4 = 0.10um
; **NOTE** Only option '0' (0.05 um monomers) is guaranteed to work.
;  The other options are still not 100% implemented yet. 

  methane = ['1.70e-03']
; PHOTO & CLIMA -> mixing ratio of methane
;     Array of strings
;     e.g. methane = ['1.00e-03', '2.00e-03', '3.00e-03']


  co2 = ['2.0E-2']
; PHOTO & CLIMA -> mixing ratio of CO2 
;     Array of strings
;     e.g. CO2 = ['2.0E-2']

  fscale = ['1.00']      
; PHOTO & CLIMA -> Solar constant scaling 
;     Array of strings
;     e.g. fscale = ['1.00']
; *NOTE*: If you set the 'times' array to be some time in the past
; and are using the SUN, you don't need to mess with this to scale
; the solar constant since Atmos knows how to scale the solar constant
; based on the time ago. Leave this at 1 if you are doing the sun and
; have set the times array. If you run Atmos with the sun and set both
; this parameter AND the times parameter, it will scale the solar
; constant twice based on this and based on the time ago. 
; If you are doing a DIFFERENT STAR THAN THE
; SUN, then set the times array to 0 ga and use fscale to alter the
; stellar insolation. 

  startype = ['13']
; PHOTO & CLIMA -> type of star
;     Array of strings
;     e.g. startype = ['13']
;key for stars (13 is default sun):
;msun = 12    !Gj 581 from Lucianne
;msun = 13    !high resolution solar data from ATLAS1/3 (Thullier et al 2004)
;msun = 14    !kevin's data from photo.dat
;msun = 15    !AD Leo from VPL climate DB
;msun = 16    !AD LEO from VPL website
;msun = 17    !T3200.dat
;msun = 18    !K2V
;msun = 19    !F2V
;msun = 76    !GJ876

  steplimit = 10                ; number of coupling iterations you will allow in each step in for loop before forcing it to move on

;-----------------STOP EDITING THIS STUFF------------------------------------------

  print, '-------------------------------------'
  print, '-------------------------------------'
  print, 'DOUBLE CHECK!'
  print, "I'm using this as the directory path: " +dirpath
  IF keyword_set(spherical) then print, 'The haze particles are spheres.'
  IF not keyword_set(spherical) then print, 'The haze particles are fractals.'
  IF keyword_set(numclima) then print, 'The number of clima iterations is: ' + numclima
  If not keyword_set(numclima) then print, 'The number of clima iterations is: 400 (default)'
  if keyword_set(doskips) then print, 'I will skip certain model combinations defined in the code by YOU. (you need to edit this code accordingly)'
  print, 'the stellar type array is :'
  print, startype
  print, 'the pressure array is: '
  print, pressures
  print, 'the oxygen array is: ' 
  print, oxygen
  print, 'the times array is: '
  print, times
  print, 'the particle density array is: '
  print, dens
  print, 'the methane mixing ratios and are: ' 
  print, methane
  print, 'the co2 mixing ratios array are: ' 
  print, co2
  print, 'the monomer size array is: ' 
  print, monsize
  print, 'fscale is:'
  print, fscale
  print, 'IMW is: ', imw
  print, 'ICONSERV is: ', iconserv  
  print, 'if not correct, type CTRL + C to cancel '

  wait, 7

;------****!!!!do not edit this section UNLESS you also edit the default input files!!!!*****----
; DEFINE DEFAULT (INITIAL) VALUES FROM DEFAULT (INITIAL) INPUT FILES
  defaultppmv = '3.00E-03'      ;default amount of methane
  defaultpress = '1.013'
  defaultga = '0.0'
  defaultflux = '3.410'
  defaultdens = '0.63'
  defaultoxy = '1.00E-08'
  defaultco2 = '1.0E-2'
  defaultstar = '13'
  defaultmonsize = '0'
  defaultfscale = '1.00'
  defaultimw = '2'              ;only needs to change at beginning
  defaulticonserv = '0'         ;only needs to change at beginning
;---------------------------------------------------------------------------------------------------
;initialize some things 
  imacounter = 0                ;I'm a counter!
  madefile = 0 

  lastmethane = defaultppmv
  lastpress = defaultpress
  lasttime = defaultga
  lastdens = defaultdens
  lastoxy  = defaultoxy
  lastco2 = defaultco2
  laststar = defaultstar
  lastmonsize = defaultmonsize
  lastfscale = defaultfscale
  
  file = dirpath+'edit_inputs.csh'
  fileuncouple = dirpath+'make_uncoupled.csh'
  filecouple = dirpath+'make_coupled.csh'
  fileiup0 = dirpath+'make_iup0.csh'
  fileiup1 = dirpath+'make_iup1.csh'
  filedefaults = dirpath+'make_defaults.csh'

  firsttime = 1
  skipthis = 0
  thiskind = 0
  NaN = 0

;the big loop
  for ss=0, n_elements(startype)-1 do begin
     for gg =0, n_elements(times) -1 do begin
        for ff = 0, n_elements(fscale) -1 do begin
           for mm=0, n_elements(methane) -1 do begin
              for pp=0, n_elements(pressures)-1 do begin
                 for dd=0, n_elements(dens)-1 do begin
                    for ox=0, n_elements(oxygen)-1 do begin
                       for cc=0, n_elements(co2)-1 do begin
                          for mon=0, n_elements(monsize)-1 do begin
                             
                             arraytemps = fltarr(100)
                             converged = 0  
                             skipthis = 0
                             prefix=''
                             prefix = prefix+'msun'+startype[ss]
                             IF keyword_set(spherical) then prefix = prefix+'sphere'

                                ;for every one of these we need to run the coupled model
                                ;------!!!!!!!!skips!!!!!!!!------------
                                ;you might not want to do every
                                ;combination...define things to skip here
                             if keyword_set(doskips) then begin
                                ;  IF startype[ss] eq '16' and float(methane[mm]) lt 4.00e-3 then skipthis=1
                                ;  IF startype[ss] eq '76' and float(methane[mm]) gt 6.00e-3 then skipthis=1
                                ;  IF startype[ss] eq '17' and float(methane[mm]) lt 6.00e-3 then skipthis=1
                                ;  IF startype[ss] eq '18' and float(methane[mm]) gt 3.0E-3 then skipthis=1
                                ;  IF startype[ss] eq '19' and float(methane[mm]) lt 3.0E-3 then skipthis=1
                             endif

                             if skipthis eq 0 then begin
                                outcounter = 0 ;count number of interations for each step in for loop
                                ;------------------------------------------------------------------------  
                                ;here we define a lot of strings that
                                ;will be used to edit input files : 
                                
                                ;IMW change - only need to do on initial step
                                imwchange = "perl -pi -e 's/IMW=       2/IMW=       "+imw+"/g' "+dirpath+"CLIMA/IO/input_clima.dat"

                                ;Ptop - top of clima pressure grid - only need to do on initial step (note this sometimes has
                                ;       troubles so you might need to check and test things) 
                                if float(pressures[pp]) le 0.500 then ptopchange = "perl -pi -e 's/P0=        7.e-6/P0=        7.e-7/g' "+dirpath+"CLIMA/IO/input_clima.dat"
                                if float(pressures[pp]) gt 0.500 then ptopchange = ''

                                ;ICONSERV change - only need to do on initial step
                                iconservchange = "perl -pi -e 's/ICONSERV=  0/ICONSERV=  "+iconserv+"/g' "+dirpath+"CLIMA/IO/input_clima.dat"

                                ;billion years ago
                                timechange = "perl -pi -e 's/"+lasttime+"      = TIMEGA - time in Ga, for modifying the solar flux/"+times[gg]+"      = TIMEGA - time in Ga, for modifying the solar flux/g' "+dirpath+"PHOTOCHEM/INPUTFILES/PLANET.dat"

                                ;hcaer density
                                denschange =  "perl -pi -e 's/HCDENS=    "+lastdens+"/HCDENS=    "+dens[dd]+"/g' "+dirpath+"PHOTOCHEM/INPUTFILES/input_photchem.dat"

                                ;monomer size
                                monchange = "perl -pi -e 's/monsize=   "+lastmonsize+"/monsize=   "+monsize[mon]+"/g' "+dirpath+"PHOTOCHEM/INPUTFILES/input_photchem.dat"

                                ;oxygen amount - mixing ratio                                
                                oxychange = "perl -pi -e 's/O2         LL  2 0 0 0 0 0    1     0.       "+lastoxy+"/O2         LL  2 0 0 0 0 0    1     0.       "+oxygen[ox]+"/g' "+dirpath+"PHOTOCHEM/INPUTFILES/species.dat"

                                ;atmospheric pressure
                                pchange = "perl -pi -e 's/"+lastpress+"    = P0  = surface pressure of modern Earth in atm/"+pressures[pp]+"    = P0  = surface pressure of modern Earth in atm/g' "+dirpath+"PHOTOCHEM/INPUTFILES/PLANET.dat"

                                ;coupling changes
                                makecoupled = "perl -pi -e 's/ICOUPLE=   0/ICOUPLE=   1/g' "+dirpath+"PHOTOCHEM/INPUTFILES/input_photchem.dat"
                                makeuncoupled = "perl -pi -e 's/ICOUPLE=   1/ICOUPLE=   0/g' "+dirpath+"PHOTOCHEM/INPUTFILES/input_photchem.dat"

                                ;co2 change
                                co2change = "perl -pi -e 's/CO2        IN  2 0 1 0 0 0    "+lastco2+"/CO2        IN  2 0 1 0 0 0    "+co2[cc]+"/g' "+dirpath+"PHOTOCHEM/INPUTFILES/species.dat"

                                ;iup (0 starts from existing solution; 1 for new solution)
                                iupmake0 = "perl -pi -e 's/IUP=       1/IUP=       0/g' "+dirpath+"CLIMA/IO/input_clima.dat"
                                iupmake1 = "perl -pi -e 's/IUP=       0/IUP=       1/g' "+dirpath+"CLIMA/IO/input_clima.dat"
                                

                                ;amount of methane ppm
                                ch4change = "perl -pi -e 's/CH4        LL  0 4 1 0 0 0    1     0.       "+lastmethane+"/CH4        LL  0 4 1 0 0 0    1     0.       "+methane[mm]+"/g' "+dirpath+"PHOTOCHEM/INPUTFILES/species.dat"


                                ;stellar type
                                starchange = "perl -pi -e 's/"+laststar+"       = msun/"+startype[ss]+"       = msun/g' "+dirpath+"PHOTOCHEM/INPUTFILES/PLANET.dat" 


                                ;fscale change
                                fscalechange = "perl -pi -e 's/"+lastfscale+"     = FSCALE/"+fscale[ff]+"     = FSCALE/g' "+dirpath+"PHOTOCHEM/INPUTFILES/PLANET.dat" 
                                ;------------------------------------------------------------------------  


                                copy = 'cp '+dirpath+'PHOTOCHEM/INPUTFILES/defaults/DEFAULT_planet.dat '+dirpath+'PHOTOCHEM/INPUTFILES/PLANET.dat'
                                copy2 = 'cp '+dirpath+'PHOTOCHEM/INPUTFILES/defaults/DEFAULT_species_newspec.dat '+dirpath+'PHOTOCHEM/INPUTFILES/species.dat'

                                IF not keyword_set(spherical) THEN copy3 = 'cp '+dirpath+'PHOTOCHEM/INPUTFILES/defaults/INITIAL_input_photchem.dat '+dirpath+'PHOTOCHEM/INPUTFILES/input_photchem.dat'
                                IF keyword_set(spherical) THEN copy3 = 'cp '+dirpath+'PHOTOCHEM/INPUTFILES/defaults/INITIAL_SPHERES_input_photchem.dat '+dirpath+'PHOTOCHEM/INPUTFILES/input_photchem.dat'
                                if float(methane[mm]) ge 1.0e-4 then copy4 = 'cp '+dirpath+'CLIMA/IO/DEFAULT_input_clima.dat '+dirpath+'CLIMA/IO/input_clima.dat'
                                if float(methane[mm]) lt 1.0e-4 then copy4 = 'cp '+dirpath+'CLIMA/IO/DEFAULT_input_clima_NO_METHANE.dat '+dirpath+'CLIMA/IO/input_clima.dat'
                                ;in.dist
                                if firsttime eq 1 or NaN eq 1 then begin ;copy new in.dist regardless if NaNs appear since old one will be contaminated 
                                   copyindist = 'cp '+dirpath+'PHOTOCHEM/default_indist/DEFAULT_in.dist '+dirpath+'PHOTOCHEM/in.dist'
                                endif
                                if not keyword_set(retainindist) then copyindist = 'cp '+dirpath+'PHOTOCHEM/default_indist/DEFAULT_in.dist '+dirpath+'PHOTOCHEM/in.dist'
                                if keyword_set(retainindist) and firsttime eq 0 and NaN eq 0 then copyindist = ''

                                ;need to change the number of clima iterations?
                                IF keyword_set(numclima) THEN $
                                   numchange = "perl -pi -e 's/NSTEPS=    400/NSTEPS=    "+numclima+"/g' "+dirpath+"CLIMA/IO/input_clima.dat"
                                IF not keyword_set(numclima) THEN $
                                   numchange = ""
                                
                                
                                IF firsttime eq 1 then perl_lines = [copyindist, copy, copy2, copy3, copy4, timechange, denschange, oxychange, $ 
                                                                     co2change, pchange, ch4change, starchange, fscalechange, monchange, $ 
                                                                     numchange, imwchange, iconservchange, ptopchange] ;by default starts off uncoupled (copy3)
                                IF firsttime ne 1 then perl_lines = [copyindist, makeuncoupled, iupmake1, timechange, denschange, $ 
                                                                     oxychange, co2change, pchange, ch4change, starchange, fscalechange, monchange]
                                
                                run = strarr(29)
                                firsttime = 0
                                ;write csh script to edit input files
                                openw, lun, file, /get_lun
                                for j =0, n_elements(perl_lines) -1 do begin 
                                   size = strlen(perl_lines[j])
                                   printf, lun, perl_lines[j], format='(a'+strtrim(size,1)+')'
                                endfor
                                close,lun

                                openw, lun, fileuncouple, /get_lun
                                size = strlen(makeuncoupled)
                                printf, lun, makeuncoupled, format='(a'+strtrim(size,1)+')'
                                close,lun

                                openw, lun, filecouple, /get_lun
                                size = strlen(makecoupled)
                                printf, lun, makecoupled, format='(a'+strtrim(size,1)+')'
                                close,lun

                                openw, lun, fileiup0, /get_lun
                                size = strlen(iupmake0)
                                printf, lun, iupmake0, format='(a'+strtrim(size,1)+')'                    
                                close,lun

                                openw, lun, fileiup1, /get_lun
                                size = strlen(iupmake1)
                                printf, lun, iupmake1, format='(a'+strtrim(size,1)+')'
                                close,lun


                                close, /all

                                ;long yet descriptive folder name! 
                                runfolder = prefix+'_'+times[gg]+"Ga_"+methane[mm]+"CH4_"+co2[cc]+"CO2_"+monsize[mon]+'rmon_'+pressures[pp]+"bar_"+oxygen[ox]+"o2_"+dens[dd]+"gcm3"+fscale[ff]+"fscale"

                                ;initial run
                                CD, dirpath
                                spawn, "chmod +x edit_inputs.csh"
                                spawn, "./edit_inputs.csh"


                                spawn, "mkdir "+dirpath+"IO/"+runfolder+"/"
                                spawn, "cp "+dirpath+"PHOTOCHEM/INPUTFILES/species.dat "+dirpath+"IO/"+runfolder+"/"

                                spawn, "./TOTCdev" ; runs photochemical model

                                ; these are the photochemical
                                ; outputs I'm choosing to save. 
                                spawn, "cp PHOTOCHEM/OUTPUT/out.out "+dirpath+"IO/"+runfolder+"/out"+strtrim(outcounter,1)+".out"
                                spawn, "cp PHOTOCHEM/OUTPUT/out.od "+dirpath+"IO/"+runfolder+"/out"+strtrim(outcounter,1)+".od"
                                spawn, "cp PHOTOCHEM/OUTPUT/out.dist PHOTOCHEM/in.dist"
                                spawn, "cp PHOTOCHEM/profile.pt "+dirpath+"IO/"+runfolder+"/profile"+strtrim(outcounter,1)+".pt"
                                spawn, "cp PHOTOCHEM/OUTPUT/out.rates "+dirpath+"IO/"+runfolder+"/out.rates"+strtrim(outcounter,1)+".pt"
                                spawn, "cp PHOTOCHEM/OUTPUT/out.params "+dirpath+"IO/"+runfolder+"/out.params"+strtrim(outcounter,1)+".pt"
                                spawn, "cp PHOTOCHEM/OUTPUT/out.prod "+dirpath+"IO/"+runfolder+"/out.prod"+strtrim(outcounter,1)+".pt"
                                spawn, "cp PHOTOCHEM/OUTPUT/out.gridz "+dirpath+"IO/"+runfolder+"/out.gridz"+strtrim(outcounter,1)+".pt"
                                spawn, "cp PHOTOCHEM/hcaer.out "+dirpath+"IO/"+runfolder+"/hcaer"+strtrim(outcounter,1)+".pt"
                                spawn, "cp PHOTOCHEM/hcaer.out "+dirpath+"IO/"+runfolder+"/hcaer_2_"+strtrim(outcounter,1)+".pt"
                                spawn, "cp PHOTOCHEM/OUTPUT/out.dist "+dirpath+"IO/"+runfolder+"/out"+strtrim(outcounter,1)+".dist"
                                spawn, "cp PHOTOCHEM/OUTPUT/out.rates "+dirpath+"IO/"+runfolder+"/out"+strtrim(outcounter,1)+".rates"
                                spawn, "cp COUPLE/fromPhoto2Clima.dat "+dirpath+"IO/"+runfolder+"/fromPhoto2Clima.dat"+strtrim(outcounter,1)+".dat"
                                spawn, "cp COUPLE/time_frak_photo.out "+dirpath+"IO/"+runfolder+"/time_frak_photo"+strtrim(outcounter,1)+".out"
                                spawn, "cp COUPLE/mixing_ratios.dat "+dirpath+"IO/"+runfolder+"/mixing_ratios"+strtrim(outcounter,1)+".dat"
                                spawn, "cp COUPLE/hcaer.photoout.out "+dirpath+"IO/"+runfolder+"/hcaer.photoout"+strtrim(outcounter,1)+".out"

                                spawn, "chmod +x make_coupled.csh" ;make it coupled now
                                spawn, "./make_coupled.csh"

                                spawn, "./runclima" ; runs climate model 

                                spawn, "cp CLIMA/IO/clima_allout.tab "+dirpath+"IO/"+runfolder+"/clima_allout"+strtrim(outcounter,1)+".tab"
                                spawn, "cp COUPLE/fromClima2Photo.dat "+dirpath+"IO/"+runfolder+"/fromClima2Photo"+strtrim(outcounter,1)+".dat"
                                spawn, "cp CLIMA/IO/TempOut.dat "+dirpath+"IO/"+runfolder+"/TempOut"+strtrim(outcounter,1)+".dat"
                                spawn, "cp CLIMA/IO/TempOut.dat CLIMA/IO/TempIn.dat"
                                spawn, "chmod +x make_iup0.csh"
                                spawn, "./make_iup0.csh" ;for this iteration of for loop -> perl lines should undo it for next step through for loop

                                ;read in initial temperature for
                                ;storage in this program
                                readcol, 'CLIMA/IO/TempOut.dat', temperature, wateramt
                                arraytemps[outcounter] = temperature[n_elements(temperature)-1]

                                if keyword_set(onestep) then converged = 1
                                
;now checks whether clima converged and if not then repeat
                                WHILE converged eq 0 DO BEGIN
                                   outcounter = outcounter + 1    ; count number of interations for each step in for loop
                                   imacounter = imacounter + 1    ; count total number of interations we do 

                                   spawn, "./TOTCdev" ; run photo

                                   spawn, "cp PHOTOCHEM/OUTPUT/out.out "+dirpath+"IO/"+runfolder+"/out"+strtrim(outcounter,1)+".out"
                                   spawn, "cp PHOTOCHEM/OUTPUT/out.od "+dirpath+"IO/"+runfolder+"/out"+strtrim(outcounter,1)+".od"
                                   spawn, "cp PHOTOCHEM/OUTPUT/out.dist PHOTOCHEM/in.dist"
                                   spawn, "cp PHOTOCHEM/profile.pt "+dirpath+"IO/"+runfolder+"/profile"+strtrim(outcounter,1)+".pt"
                                   spawn, "cp PHOTOCHEM/OUTPUT/out.dist "+dirpath+"IO/"+runfolder+"/out"+strtrim(outcounter,1)+".dist"
                                   spawn, "cp PHOTOCHEM/OUTPUT/out.rates "+dirpath+"IO/"+runfolder+"/out"+strtrim(outcounter,1)+".rates"                           
                                   spawn, "cp PHOTOCHEM/hcaer.out "+dirpath+"IO/"+runfolder+"/hcaer"+strtrim(outcounter,1)+".out"
                                   spawn, "cp PHOTOCHEM/hcaer2.out "+dirpath+"IO/"+runfolder+"/hcaer_2_"+strtrim(outcounter,1)+".out"
                                   spawn, "cp PHOTOCHEM/OUTPUT/out.params "+dirpath+"IO/"+runfolder+"/out.params"+strtrim(outcounter,1)+".pt"
                                   spawn, "cp PHOTOCHEM/OUTPUT/out.prod "+dirpath+"IO/"+runfolder+"/out.prod"+strtrim(outcounter,1)+".pt"
                                   spawn, "cp PHOTOCHEM/OUTPUT/out.gridz "+dirpath+"IO/"+runfolder+"/out.gridz"+strtrim(outcounter,1)+".pt"
                                   spawn, "cp COUPLE/fromPhoto2Clima.dat "+dirpath+"IO/"+runfolder+"/fromPhoto2Clima.dat"+strtrim(outcounter,1)+".dat"
                                   spawn, "cp COUPLE/time_frak_photo.out "+dirpath+"IO/"+runfolder+"/time_frak_photo"+strtrim(outcounter,1)+".out"
                                   spawn, "cp COUPLE/mixing_ratios.dat "+dirpath+"IO/"+runfolder+"/mixing_ratios"+strtrim(outcounter,1)+".dat"
                                   spawn, "cp COUPLE/hcaer.photoout.out "+dirpath+"IO/"+runfolder+"/hcaer.photoout"+strtrim(outcounter,1)+".out"

                                   spawn, "./runclima" ;run clima

                                   spawn, "cp CLIMA/IO/TempOut.dat CLIMA/IO/TempIn.dat"
                                   spawn, "cp CLIMA/IO/clima_allout.tab "+dirpath+"IO/"+runfolder+"/clima_allout"+strtrim(outcounter,1)+"_first.tab"
                                   spawn, "cp CLIMA/IO/TempOut.dat "+dirpath+"IO/"+runfolder+"/TempOut"+strtrim(outcounter,1)+"_first.dat"
                                   if keyword_set(climatwice) then begin
                                      spawn, "./runclima" ;run again to make sure clima converges 
                                      spawn, "cp CLIMA/IO/TempOut.dat CLIMA/IO/TempIn.dat"
                                   endif

                                   spawn, "cp CLIMA/IO/clima_allout.tab "+dirpath+"IO/"+runfolder+"/clima_allout"+strtrim(outcounter,1)+".tab"
                                   spawn, "cp COUPLE/fromClima2Photo.dat "+dirpath+"IO/"+runfolder+"/fromClima2Photo"+strtrim(outcounter,1)+".dat"
                                   spawn, "cp CLIMA/IO/TempOut.dat "+dirpath+"IO/"+runfolder+"/TempOut"+strtrim(outcounter,1)+".dat"
                                   readcol, 'CLIMA/IO/TempOut.dat', temperature, wateramt
                                   arraytemps[outcounter] = temperature[n_elements(temperature)-1]
                                   deltaT = abs(arraytemps[outcounter] - arraytemps[outcounter-1])
                                   print, 'delta T is ', strtrim(deltaT,1)
                                   
                                ;sometimes it hits some kinda funky
                                ;oscillating equlibrium
                                ;point...won't go on unless we
                                ;force it:
                                   IF outcounter gt 4 then begin ;check for oscillation 
                                      if abs(arraytemps[outcounter] - arraytemps[outcounter-2]) lt 1.5 then begin
                                         if abs(arraytemps[outcounter-1] - arraytemps[outcounter-3]) lt 1.5 then begin 
                                            converged =1
                                            spawn, "cp CLIMA/IO/TempOut.dat "+dirpath+"IO/"+runfolder+"/OSCILLATION_DETECTED.TXT"
                                            NaN = 0
                                         endif
                                      endif
                                   endif

                                   
                                   
                                   IF finite(temperature[outcounter]) eq 0 THEN BEGIN ;check for NaNs 
                                      converged =1
                                      NaN = 1
                                      print, 'model is producing NaNs.  moving on to next step in for loop.'
                                      spawn, "cp CLIMA/IO/TempOut.dat "+dirpath+"IO/"+runfolder+"/KILLED_FROM_NANS.TXT"
                                      
                                   ENDIF

                                   IF deltaT LT 1.0 THEN BEGIN ; check for convergence (note this doesn't mean you shouldn't check the files manually!)
                                      converged = 1
                                      print, 'model is converged (delta T < 1) for this step in for loop'
                                      NaN = 0
                                   ENDIF ELSE BEGIN
                                      print, 'model is not converged yet (delta T > 1).  Trying again.'
                                      NaN = 0
                                   ENDELSE

                                   IF outcounter GT steplimit THEN BEGIN ;check for exceeding step limit
                                      converged = 1
                                      print, 'model did not converge after '+strtrim(steplimit,1)+' iterations.  moving on to next step in for loop.'
                                      spawn, "cp CLIMA/IO/TempOut.dat "+dirpath+"IO/"+runfolder+"/DID_NOT_CONVERGE.TXT"
                                      NaN = 0
                                   ENDIF
                                   
                                ENDWHILE ;while not converged 


                                ;to finish up copy outputs to photochem_smart:
                                spawn, "cp PHOTOCHEM/profile.pt "+dirpath+"photochem_smart/photchem/"+prefix+'_'+times[gg]+"Ga_"+methane[mm]+"ch4_"+co2[cc]+"_"+monsize[mon]+'rmon_'+pressures[pp]+"bar_"+oxygen[ox]+"o2_"+dens[dd]+"gcm3"+fscale[ff]+"fscale"
                                spawn, "cp PHOTOCHEM/hcaer.out "+dirpath+"photochem_smart/hcaer/hcaer_"+prefix+'_'+times[gg]+"Ga_"+methane[mm]+"ch4_"+co2[cc]+"_"+monsize[mon]+'rmon_'+pressures[pp]+"bar_"+oxygen[ox]+"o2_"+dens[dd]+"gcm3"+fscale[ff]+"fscale"
                                spawn, "cp PHOTOCHEM/hcaer2.out "+dirpath+"photochem_smart/hcaer/hcaer2_"+prefix+'_'+times[gg]+"Ga_"+methane[mm]+"ch4_"+co2[cc]+"_"+monsize[mon]+'rmon_'+pressures[pp]+"bar_"+oxygen[ox]+"o2_"+dens[dd]+"gcm3"+fscale[ff]+"fscale"

                                
                                imacounter = imacounter + 1

                                lastmethane = methane[mm]
                                lastpress = pressures[pp]
                                lasttime = times[gg]
                                lastdens = dens[dd]
                                lastoxy = oxygen[ox]
                                laststar = startype[ss]
                                lastmonsize = monsize[mon]
                                lastfscale = fscale[ff]
                             endif

                          endfor

                       endfor
                    endfor
                 endfor
              endfor
           endfor
        endfor
     endfor


  endfor

  
  print, 'done!!!  I did a total of ' +strtrim(imacounter,1)+' iterations!  wow!'
;done!

END

