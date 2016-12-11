PRO run_coupled, rmix=rmix, flux=flux, ADLEO=ADLEO, spherical=spherical, dirpath=dirpath, timeago=timeago, numclima=numclima
;-------------------------------------
;-------WHAT IT DOES:-----------------
;This IDL code runs the coupled climate-photochemical model for a variety
;of conditions set under the "EDIT THIS STUFF" label found below.
;It considers the model converged when change in temperature between runs
;is less than 1 K.  If it takes > 10 iterations to get to that,
;it will give up and move to the next run.
;--------------------------------------
;--------HOW TO RUN IT:----------------
;1. Edit the items under the "Edit this stuff" label below to 
;   choose what kind of photochemical-climate runs you want to do.
;   You can change these parameters there:
;  pressures (atmos pressure in bars)
;  oxygen    (oxygen mixing ratio)
;  timeago   (time ago to run model in billions of years [ga] - note you can also set from program call)
;  dens      (density of haze particle material in g/cm^3)
;  methane   (methane mixing ratio OR flux -- choose whether mixing ratio (rmix) or flux)
;  co2       (CO2 mixing ratio)
;
;2. Edit dirpathdefault as needed below (can also change working
;directory from program call or it'll prompt you)
;
;3. Set optional keyword flags.  They are:
;   rmix - set flag to make methane amount mixing ratio
;   flux - set flag to make methane amount a flux 
;   ADLEO - set flag to make star planet is orbiting as AD Leo
;   spherical - set flag to make hydrocarbon aerosols spheres (fractals set by
;               default)
;   dirpath - string. set as string of directory to work in (should be TOP LEVEL directory
;             of atmos) e.g. '/astro/users/giada/atmos_testing/'
;   timeago - string (really). you can either set the time ago the model runs here or
;             in the "edit this stuff" arrays below.  It was
;             convenient for me to allow both in the way I run this
;             program.  e.g. '0.0', '2.7' for 0 or 2.7 Ga (default is
;             what's in the times array below)
;   numclima - string integer. set number of iterations clima goes through
;              It will do 400 by default, which takes forever, but
;              that huge number may be necessary if you find the model
;              is having a hard time converging.
;EXAMPLES:
; e.g. run it with methane mixing ratios, around ad leo, and 200 clima
; iterations:
; IDL> run_coupled, /rmix, /adleo, numclima = '200'
;
; e.g. run it with methane fluxes, spherical particles, 2.5 billion
; years ago:
; IDL> run_coupled, /flux, /spherical, timeago='2.5' 
;---------------------------------------------------------------

;------prompt user for input--------------

  IF keyword_set(rmix) THEN type = 'rmix'
  IF keyword_set(flux) THEN type = 'flux'

  IF not keyword_set(rmix) and not keyword_set(flux) THEN BEGIN
     print, 'do you want methane mixing ratio or flux?'
     input = ''
     READ, input, prompt = 'Type rmix/flux: '
     IF input EQ 'rmix' THEN type = 'rmix'
     IF input EQ 'flux' THEN type = 'flux'
     IF input NE 'rmix' and input NE 'flux' THEN BEGIN
        print, 'invalid input.  you need to type "rmix" or "flux" (lowercase, no quotes when you type it)'
        stop
     ENDIF
  ENDIF

  dirpathdefault = '/astro/users/giada/atmos_testing/' ;path of the coupled model

  IF not keyword_set(dirpath) THEN BEGIN
     print, ''
     print, '-----'
     print, 'optional directory path not set.'
     print, "going to use this default one unless I'm told otherwise: "+dirpathdefault
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
        print, 'invalid input.  you need to type either "y" or "n" or "here" (lowercase, no quotes when you enter it)'
        stop
     ENDIF

     ;check to make sure the directory is formatted correctly:
     lastchar = strmid(dirpath, 0,1,/reverse_offset)
     firstpart = strmid(dirpath, 0,11)
     IF lastchar NE '/' THEN dirpath = dirpath+'/'
     IF firstpart EQ 'astro/users' THEN dirpath = '/'+dirpath

  ENDIF
;----------!!!!!!EDIT THIS STUFF!!!!!!!!!-----------------
  
  pressures = ['1.013']
  oxygen = ['1.0E-08']
  IF not keyword_set(timeago) THEN times = ['2.7']  ;in Ga
  IF keyword_set(timeago) THEN times = [timeago]
  dens = ['0.63']
                                ;IF type eq 'rmix' THEN methane = ['1.0E-02', '2.0E-03', '2.5E-03', '3.0E-03', '3.5E-03'] ;rmix
                                ;IF type eq 'rmix' THEN methane = ['2.0E-03', '2.1E-03', '2.2E-03', '2.3E-03','2.4E-03', '2.5E-03', '2.6E-03', '2.7E-03', '2.8E-03', '2.9E-03', '3.0E-03'] ;rmix
  IF type eq 'rmix' THEN methane = ['2.0E-03', '2.1E-03', '2.2E-03', '2.3E-03', '2.4E-03', '2.5E-03', '2.6.0E-03', '2.7E-03', '2.8E-03', '2.9E-03'] ;rmix
;  IF type eq 'rmix' THEN methane = ['2.4E-03', '2.5E-03']

;for fluxes they are x10^11, but don't write that part
  IF type eq 'flux' THEN methane = ['3.410', '1.000', '0.300', '5.000'] ;fluxes times E+11

  co2 = ['1.0E-2']
  
;----------------------------------------------

  print, '-------------------------------------'
  print, '-------------------------------------'
  print, 'DOUBLE CHECK'
  IF keyword_set(ADLEO) THEN print, "I'm using star AD Leo"
  IF not keyword_set(ADLEO) then print, "I'm using star sun"
  print, "I'm using this as the directory path: " +dirpath
  IF keyword_set(spherical) then print, 'The particles are spheres.'
  IF not keyword_set(spherical) then print, 'The particles are fractals.'
  IF keyword_set(numclima) then print, 'The number of clima iterations is: ' + numclima
  If not keyword_set(numclima) then print, 'The number of clima iterations is: 400 (default)'
  print, 'the pressure array is: '
  print, pressures
  print, 'the oxygen array is: ' 
  print, oxygen
  print, 'the times array is: '
  print, times
  print, 'the particle density array is: '
  print, dens



  If type eq 'rmix' then begin
     print, 'the methane is mixing ratios and the array is: ' 
     print, methane
  ENDIF

  If type eq' flux' then begin
     print, 'the methane is fluxes and the array is: ' 
     print, methane
  ENDIF

  print, 'the co2 array is: ' 
  print, co2
  
  print, 'if not correct, type CTRL + C'

  wait, 5

;------do not edit UNLESS you also edit the default input files----
  defaultppmv = '3.0E-03'
  defaultpress = '1.013'
  defaultga = '0.0'
  defaultflux = '3.410'
  defaultdens = '0.63'
  defaultoxy = '1.0E-08'
  defaultco2 = '1.0E-2'
;--------------------------------------------------------------------
  prefix=''
  IF keyword_set(ADLEO) then prefix = prefix+'adleo'
  IF keyword_set(spherical) then prefix = prefix+'sphere'


  imacounter= 0
  madefile=0

  IF type eq 'rmix'  then lastmethane = defaultppmv
  IF type eq 'flux'  then lastmethane = defaultflux
  lastpress = defaultpress
  lasttime = defaultga
  lastdens = defaultdens
  lastoxy  = defaultoxy
  lastco2 = defaultco2
  
  file = dirpath+'edit_inputs.csh'
  fileuncouple = dirpath+'make_uncoupled.csh'
  filecouple = dirpath+'make_coupled.csh'
  fileiup0 = dirpath+'make_iup0.csh'
  fileiup1 = dirpath+'make_iup1.csh'
  filedefaults = dirpath+'make_defaults.csh'



  firsttime = 1

  for gg =0, n_elements(times) -1 do begin
     for mm=0, n_elements(methane) -1 do begin
        for pp=0, n_elements(pressures)-1 do begin
           for dd=0, n_elements(dens)-1 do begin
              for ox=0, n_elements(oxygen)-1 do begin
                 for cc=0, n_elements(co2)-1 do begin
                    arraytemps = fltarr(100)
                    converged = 0  
;for every one of these we need to run the coupled model
                    
                    outcounter = 0

;billion years ago
                    timechange = "perl -pi -e 's/"+lasttime+"      = TIMEGA - time in Ga, for modifying the solar flux/"+times[gg]+"      = TIMEGA - time in Ga, for modifying the solar flux/g' "+dirpath+"PHOTOCHEM/INPUTFILES/PLANET.dat"

;hcaer density
                    denschange =  "perl -pi -e 's/HCDENS=    "+lastdens+"/HCDENS=    "+dens[dd]+"/g' "+dirpath+"PHOTOCHEM/INPUTFILES/input_photchem.dat"

;oxygen amount
                    oxychange = "perl -pi -e 's/O2         LL  2 0 0 0 0 0    1     0.      "+lastoxy+"/O2         LL  2 0 0 0 0 0    1     0.      "+oxygen[ox]+"/g' "+dirpath+"PHOTOCHEM/INPUTFILES/species.dat"

;atmos pressure
                    pchange = "perl -pi -e 's/"+lastpress+"    = P0  = surface pressure of modern Earth in atm/"+pressures[pp]+"    = P0  = surface pressure of modern Earth in atm/g' "+dirpath+"PHOTOCHEM/INPUTFILES/PLANET.dat"

;couple change 
                    makecoupled = "perl -pi -e 's/ICOUPLE=   0/ICOUPLE=   1/g' "+dirpath+"PHOTOCHEM/INPUTFILES/input_photchem.dat"
                    makeuncoupled = "perl -pi -e 's/ICOUPLE=   1/ICOUPLE=   0/g' "+dirpath+"PHOTOCHEM/INPUTFILES/input_photchem.dat"

;co2 change
                    co2change = "perl -pi -e 's/CO2        IN  2 0 1 0 0 0    "+lastco2+"/CO2        IN  2 0 1 0 0 0    "+co2[cc]+"/g' "+dirpath+"PHOTOCHEM/INPUTFILES/species.dat"

;iup (0 starts from existing solution; 1 for new solution)
                    iupmake0 = "perl -pi -e 's/IUP=       1/IUP=       0/g' "+dirpath+"CLIMA/IO/input_clima.dat"
                    iupmake1 = "perl -pi -e 's/IUP=       0/IUP=       1/g' "+dirpath+"CLIMA/IO/input_clima.dat"


;amount of methane ppm
                    IF type eq 'rmix' THEN BEGIN
                       ch4change = "perl -pi -e 's/CH4        LL  0 4 1 0 0 0    1     0.      "+lastmethane+"/CH4        LL  0 4 1 0 0 0    1     0.      "+methane[mm]+"/g' "+dirpath+"PHOTOCHEM/INPUTFILES/species.dat"
                    ENDIF

;amount of methane flux
                    IF type eq 'flux' THEN BEGIN
                       IF gg eq 0 and mm eq 0 and pp eq 0 then begin                
                          ch4change = "perl -pi -e 's/CH4        LL  0 4 1 0 0 0    1     0.      1.0E-02 "+lastmethane+"/CH4        LL  0 4 1 0 0 0    2     0.      1.0E-02 "+methane[mm]+"/g' "+dirpath+"PHOTOCHEM/INPUTFILES/species.dat"            
                       ENDIF ELSE BEGIN             
                          ch4change = "perl -pi -e 's/CH4        LL  0 4 1 0 0 0    2     0.      1.0E-02 "+lastmethane+"/CH4        LL  0 4 1 0 0 0    2     0.      1.0E-02 "+methane[mm]+"/g' "+dirpath+"PHOTOCHEM/INPUTFILES/species.dat"
                       ENDELSE

                    ENDIF

                    IF not keyword_set(ADLEO) THEN copy = 'cp '+dirpath+'PHOTOCHEM/INPUTFILES/DEFAULT_planet.dat '+dirpath+'PHOTOCHEM/INPUTFILES/PLANET.dat'
                    IF keyword_set(ADLEO) THEN copy = 'cp '+dirpath+'PHOTOCHEM/INPUTFILES/DEFAULT_planet_adleo.dat '+dirpath+'PHOTOCHEM/INPUTFILES/PLANET.dat'
                    copy2 = 'cp '+dirpath+'PHOTOCHEM/INPUTFILES/DEFAULT_species.dat '+dirpath+'PHOTOCHEM/INPUTFILES/species.dat'
                    IF not keyword_set(spherical) THEN copy3 = 'cp '+dirpath+'PHOTOCHEM/INPUTFILES/INITIAL_input_photchem.dat '+dirpath+'PHOTOCHEM/INPUTFILES/input_photchem.dat'
                    IF keyword_set(spherical) THEN copy3 = 'cp '+dirpath+'PHOTOCHEM/INPUTFILES/INITIAL_SPHERES_input_photchem.dat '+dirpath+'PHOTOCHEM/INPUTFILES/input_photchem.dat'
                    copy4 = 'cp '+dirpath+'CLIMA/IO/DEFAULT_input_clima.dat '+dirpath+'CLIMA/IO/input_clima.dat'
                    copyindist = 'cp '+dirpath+'PHOTOCHEM/DEFAULT_in.dist '+dirpath+'PHOTOCHEM/in.dist'

                    ;need to change the number of clima iterations?
                    IF keyword_set(numclima) THEN $
                       numchange = "perl -pi -e 's/NSTEPS=    400/NSTEPS=    "+numclima+"/g' "+dirpath+"CLIMA/IO/input_clima.dat"
                    IF not keyword_set(numclima) THEN $
                       numchange = ""
                    
                    perl_lines = [copyindist, makeuncoupled, iupmake1, timechange, denschange, oxychange, co2change, pchange, ch4change]
                    IF firsttime eq 1 then perl_lines = [copyindist, copy, copy2, copy3, copy4, timechange, denschange, oxychange, $ 
                                                         co2change, pchange, ch4change,numchange] ;by default starts off uncoupled (copy3)
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

                    

                    runfolder = prefix+'_'+times[gg]+"Ga_"+methane[mm]+"ch4_"+type+"_"+co2[cc]+"_"+pressures[pp]+"bar_"+oxygen[ox]+"o2_"+dens[dd]+"gcm3"


                                ;initial run
                    CD, dirpath
                    spawn, "chmod +x edit_inputs.csh"
                    spawn, "./edit_inputs.csh"

                    spawn, "mkdir "+dirpath+"IO/"+runfolder+"/"                    
                    spawn, "./TOTCdev"                   
                    spawn, "cp PHOTOCHEM/OUTPUT/out.out "+dirpath+"IO/"+runfolder+"/out"+strtrim(outcounter,1)+".out"
                    spawn, "cp PHOTOCHEM/OUTPUT/out.dist PHOTOCHEM/in.dist"
                    spawn, "cp PHOTOCHEM/profile.pt "+dirpath+"IO/"+runfolder+"/profile"+strtrim(outcounter,1)+".pt"
                    spawn, "cp PHOTOCHEM/hcaer.out "+dirpath+"IO/"+runfolder+"/hcaer"+strtrim(outcounter,1)+".pt"
                    spawn, "cp COUPLE/fromPhoto2Clima.dat "+dirpath+"IO/"+runfolder+"/fromPhoto2Clima.dat"+strtrim(outcounter,1)+".dat"
                    spawn, "cp COUPLE/time_frak_photo.out "+dirpath+"IO/"+runfolder+"/time_frak_photo"+strtrim(outcounter,1)+".out"
                    spawn, "cp COUPLE/mixing_ratios.dat "+dirpath+"IO/"+runfolder+"/mixing_ratios"+strtrim(outcounter,1)+".dat"
                    spawn, "cp COUPLE/hcaer.photoout.out "+dirpath+"IO/"+runfolder+"/hcaer.photoout"+strtrim(outcounter,1)+".out"

                    spawn, "chmod +x make_coupled.csh"
                    spawn, "./make_coupled.csh"

                    spawn, "./runclima"
                    spawn, "cp CLIMA/IO/clima_allout.tab "+dirpath+"IO/"+runfolder+"/clima_allout"+strtrim(outcounter,1)+".tab"
                    spawn, "cp COUPLE/fromClima2Photo.dat "+dirpath+"IO/"+runfolder+"/fromClima2Photo"+strtrim(outcounter,1)+".dat"
                    spawn, "cp CLIMA/IO/TempOut.dat "+dirpath+"IO/"+runfolder+"/TempOut"+strtrim(outcounter,1)+".dat"
                    spawn, "cp CLIMA/IO/TempOut.dat CLIMA/IO/TempIn.dat"
                    spawn, "chmod +x make_iup0.csh"
                    spawn, "./make_iup0.csh" ;for this iteration of for loop -> perl lines should undo it for next step through for loop

                                ;read in initial temperature
                    readcol, 'CLIMA/IO/TempOut.dat', temperature, wateramt
                    arraytemps[outcounter] = temperature[n_elements(temperature)-1]

                    
;now checks whether clima converged and if not then repeat
                    WHILE converged eq 0 DO BEGIN
                       outcounter = outcounter + 1
                       imacounter = imacounter + 1
                       spawn, "./TOTCdev"
                       spawn, "cp PHOTOCHEM/OUTPUT/out.out "+dirpath+"IO/"+runfolder+"/out"+strtrim(outcounter,1)+".out"
                       spawn, "cp PHOTOCHEM/OUTPUT/out.dist PHOTOCHEM/in.dist"
                       spawn, "cp PHOTOCHEM/profile.pt "+dirpath+"IO/"+runfolder+"/profile"+strtrim(outcounter,1)+".pt"
                       spawn, "cp PHOTOCHEM/hcaer.out "+dirpath+"IO/"+runfolder+"/hcaer"+strtrim(outcounter,1)+".out"
                       spawn, "cp COUPLE/fromPhoto2Clima.dat "+dirpath+"IO/"+runfolder+"/fromPhoto2Clima.dat"+strtrim(outcounter,1)+".dat"
                       spawn, "cp COUPLE/time_frak_photo.out "+dirpath+"IO/"+runfolder+"/time_frak_photo"+strtrim(outcounter,1)+".out"
                       spawn, "cp COUPLE/mixing_ratios.dat "+dirpath+"IO/"+runfolder+"/mixing_ratios"+strtrim(outcounter,1)+".dat"
                       spawn, "cp COUPLE/hcaer.photoout.out "+dirpath+"IO/"+runfolder+"/hcaer.photoout"+strtrim(outcounter,1)+".out"
                       spawn, "./runclima"
                       spawn, "cp CLIMA/IO/TempOut.dat CLIMA/IO/TempIn.dat"
                       spawn, "./runclima" ;run again to make sure clima converges because ARGH CLIMA
                       spawn, "cp CLIMA/IO/TempOut.dat CLIMA/IO/TempIn.dat"
                       spawn, "cp CLIMA/IO/clima_allout.tab "+dirpath+"IO/"+runfolder+"/clima_allout"+strtrim(outcounter,1)+".tab"
                       spawn, "cp COUPLE/fromClima2Photo.dat "+dirpath+"IO/"+runfolder+"/fromClima2Photo"+strtrim(outcounter,1)+".dat"
                       spawn, "cp CLIMA/IO/TempOut.dat "+dirpath+"IO/"+runfolder+"/TempOut"+strtrim(outcounter,1)+".dat"
                       readcol, 'CLIMA/IO/TempOut.dat', temperature, wateramt
                       arraytemps[outcounter] = temperature[n_elements(temperature)-1]
                       deltaT = abs(arraytemps[outcounter] - arraytemps[outcounter-1])
                       print, 'delta T is ', strtrim(deltaT,1)

                       IF finite(temperature[outcounter]) eq 0 THEN BEGIN
                          converged =1
                          print, 'model is producing NaNs.  moving on to next step in for loop.'
                          spawn, "cp CLIMA/IO/TempOut.dat "+dirpath+"IO/"+runfolder+"/KILLED_FROM_NANS.TXT"
                       ENDIF

                       IF deltaT LT 1.0 THEN BEGIN
                          converged = 1
                          print, 'model is converged (delta T < 1) for this step in for loop'
                       ENDIF ELSE BEGIN
                          print, 'model is not converged yet (delta T > 1).  Trying again.'
                       ENDELSE

                       IF outcounter GT 10 THEN BEGIN
                          converged = 1
                          print, 'model did not converge after 10 iterations.  moving on to next step in for loop.'
                          spawn, "cp CLIMA/IO/TempOut.dat "+dirpath+"IO/"+runfolder+"/DID_NOT_CONVERGE.TXT"
                          
                       ENDIF
                       
                    ENDWHILE



;to finish up:
                    spawn, "cp PHOTOCHEM/profile.pt "+dirpath+"photochem_smart/photchem/"+prefix+times[gg]+"Ga_"+methane[mm]+"ch4_"+type+"_"+co2[cc]+"_"+pressures[pp]+"bar_"+oxygen[ox]+"o2_"+dens[dd]+"gcm3.pt"
                    spawn, "cp PHOTOCHEM/hcaer.out "+dirpath+"photochem_smart/hcaer/hcaer_"+prefix+times[gg]+"Ga_"+methane[mm]+"ch4_"+type+"_"+co2[cc]+"_"+pressures[pp]+"bar_"+oxygen[ox]+"o2_"+dens[dd]+"gcm3_.out"
                    spawn, "cp PHOTOCHEM/hcaer2.out "+dirpath+"photochem_smart/hcaer/hcaer2_"+prefix+times[gg]+"Ga_"+methane[mm]+"ch4_"+type+"_"+co2[cc]+"_"+pressures[pp]+"bar_"+oxygen[ox]+"o2_"+dens[dd]+"gcm3_.out"
                                ; spawn, "./make_uncoupled.txt" ;not
                                ; necessary b/c this is in perl_lines

                    imacounter = imacounter + 1

                    lastmethane = methane[mm]
                    lastpress = pressures[pp]
                    lasttime = times[gg]
                    lastdens = dens[dd]
                    lastoxy = oxygen[ox]
                    lastco2 = co2[cc]

                 endfor
              endfor
           endfor
        endfor
     endfor
  endfor


;restore defaults:
;billion years ago
  defaulttime = "perl -pi -e 's/"+lasttime+"      = TIMEGA - time in Ga, for modifying the solar flux/"+defaultga+"      = TIMEGA - time in Ga, for modifying the solar flux/g' "+dirpath+"PHOTOCHEM/INPUTFILES/PLANET.dat"

;atmos pressure
  defaultp = "perl -pi -e 's/"+lastpress+"    = P0  = surface pressure of modern Earth in atm/"+defaultpress+"    = P0  = surface pressure of modern Earth in atm/g' "+dirpath+"PHOTOCHEM/INPUTFILES/PLANET.dat"

;amount of methane ppm
  IF type eq 'rmix' THEN BEGIN
     defaultch4 = "perl -pi -e 's/CH4        LL  0 4 1 0 0 0    1     0.      "+lastmethane+"/CH4        LL  0 4 1 0 0 0    1     0.      "+defaultppmv+"/g' "+dirpath+"PHOTOCHEM/INPUTFILES/species.dat"
  ENDIF

;amount of methane flux
  IF type eq 'flux' THEN BEGIN
     defaultch4 = "perl -pi -e 's/CH4        LL  0 4 1 0 0 0    2     0.      1.0E-02 "+lastmethane+"/CH4        LL  0 4 1 0 0 0    1     0.      1.0E-02 "+defaultflux+"/g' "+dirpath+"PHOTOCHEM/INPUTFILES/species.dat"
  ENDIF

;co2
  co2change = "perl -pi -e 's/CO2        IN  2 0 1 0 0 0    "+lastco2+"/CO2        IN  2 0 1 0 0 0    "+defaultco2+"+/g' "+dirpath+"PHOTOCHEM/INPUTFILES/species.dat"

;hcaer density
  defaultdens =  "perl -pi -e 's/HCDENS=    "+lastdens+"/HCDENS=    "+defaultdens+"/g' "+dirpath+"PHOTOCHEM/INPUTFILES/input_photchem.dat"

;oxygen amount
  defaultoxy = "perl -pi -e 's/O2         LL  2 0 0 0 0 0    1     0.      "+lastoxy+"/O2         LL  2 0 0 0 0 0    1     0.      "+defaultoxy+"/g' "+dirpath+"PHOTOCHEM/INPUTFILES/species.dat"


  
  defaults = [defaultp, defaultch4, defaulttime, defaultdens, defaultoxy, defaultco2]

  openw, lun, filedefaults, /get_lun
  for j =0, n_elements(defaults) -1 do begin 
     size = strlen(defaults[j])
     printf, lun, defaults[j], format='(a'+strtrim(size,1)+')'
  endfor
  close,lun
  close, /all

  spawn, "chmod +x make_defaults.csh"
  spawn, "./make_defaults.csh"

  
  print, 'done!!!  I did a total of ' +strtrim(imacounter,1)+' iterations!  wow!'
;done!

END

