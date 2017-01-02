#!/bin/bash
#############################################################################
# TestModels.scr - Model Testing Script                                     #
# This script will test all of the templates for PHOTOCHEM and CLIMA. This  #
# happens in three stages. First, it copies the files in place, then it     #
# comptiles the model with the files in place. Finally, it runs the code.   #
# This process is repeated for each template. There are traps set up at     #
# each stage - for when files are not present, the code will not compile,   #
# or the code will not run. The script will stop if any of these conditions #
# are true. The script has a trap but will NOT stop if the model runs       #
# "successfully" but does not converge.                                     #
#                                                                           #
# We plan to add flags for when the model converges but to a different      #
# position than it was in before.                                           #
#############################################################################
# Set paths and constants
temp_path='PHOTOCHEM/INPUTFILES/TEMPLATES/'
temp_path2='CLIMA/IO/TEMPLATES/'
temp_path3='CLIMA/IO/'
RED='\033[0;31m' #Red
NC='\033[0m' # No Color
GREEN='\033[0;32m' #Green
YELLOW='\033[0;33m'
true_path=$(pwd)

# *** Initial messages to user ***
clear
echo "***************** THIS IS A TEST. THIS IS ONLY A TEST ****************"
echo "This script will automatically run your code for each of the *working*"
echo "templates. It will run the code once for each of these templates, and"
echo "check the results from prior runs of each template. It will flag cases"
echo "where the code doesn't compile, crashes, or doesn't converge."
echo "**********************************************************************"
echo " "
echo "You will need to be in your top folder for atmos, not a subfolder."
echo "Your current folder is $(pwd)."
echo ""
echo "Your PHOTOCHEM has these templates:"
ls $temp_path
echo ""
# *** Initiate loop for testing ***
# Find all the templates and copy/compile/run/test each one
cd $temp_path
# This command finds and loops on all the top-level template folders
for folder in `find . -type d -maxdepth 1 -mindepth 1 | sed 's|./||'`; do
  cd $true_path
  cd $temp_path
  # copy files
  cd $folder
  echo "Copying files from ${temp_path}${folder}."
  # *** Copy files and exit script for missing files ***
  if cp 'in.dist' '../../..'; then
    if cp 'input_photchem.dat' '../..'; then
      if cp 'reactions.rx' '../..'; then
        if cp 'parameters.inc' '../..'; then
          if cp 'species.dat' '../..'; then
            if cp 'PLANET.dat' '../..'; then
              printf "${GREEN}${folder}Finished copying photochem templates over for ${folder} run!${NC}"
              echo ""
            else
              printf "${RED}${temp_path}${folder}/PLANET.dat appears to be missing! Quitting.${NC}"
              echo ""
              exit 1
            fi
          else
            printf "${RED}${temp_path}${folder}/species.dat appears to be missing! Quitting.${NC}"
            echo ""
            exit 1
          fi
        else
          printf "${RED}${temp_path}${folder}/parameters.inc appears to be missing! Quitting.${NC}"
          echo ""
          exit 1
        fi
      else
        printf "${RED}${temp_path}${folder}/reactions.rx appears to be missing! Quitting.${NC}"
        echo ""
        exit 1
      fi
    else
      printf "${RED}${temp_path}${folder}/input_photochem.dat appears to be missing! Quitting.${NC}"
      echo ""
      exit 1
    fi
  else
    printf "${RED}${temp_path}${folder}/in.dist appears to be missing! Quitting.${NC}"
    echo ""
    exit 1
  fi
  cd $true_path

  # *** PHOTOCHEM tests ***
  # Compile photochemistry code
  echo "Compiling ${folder}"
  make '-f' 'PhotoMake' 'clean' >> ${folder}.log 2>&1
  if make '-f' 'PhotoMake' >> ${folder}.log 2>&1; then
    printf "${GREEN}${folder} compiled!${NC}"
    echo ""
  else
    # exit if compliation fails
    printf "${RED}${folder} did not compile! Something is wrong. Quitting.${NC}"
    echo "See ${folder}.log for compliation output and error messages."
    echo ""
    exit 1
  fi

  # run photochemistry code
  echo ${folder} > RunModel.txt
  echo "y" >> RunModel.txt
  echo "y" >> RunModel.txt
  echo "n" >> RunModel.txt
  echo "n" >> RunModel.txt
  echo "Running ${folder}....... Please wait. This might take a while."
  echo "If you need to peek take a look at ${folder}.out.out."
  if ./'Photo.run' >& ${folder}.out.out; then
    # Tests for non-convergence
    echo "Here is information from your last time step:"
    echo $(grep "N =" ${folder}.out.out | tail -n 1)
    NSTEPS=$(grep "NSTEPS =" PhotoMain.f | grep -vi c | sed 's/NSTEPS = //' | sed 's/ //g')
    if grep "N =" ${folder}.out.out | grep ${NSTEPS} >/dev/null 2>&1; then
      # In this case, the code took too long to run.
      # We flag this case, but do not stop the script.
      printf "${RED}${folder} did not converge after ${NSTEPS} steps!${NC}"
      echo ""
      echo "See ${folder}.out.out for model output to diagnose this."
      echo "Please DO NOT submit your code changes until resolving this issue."
    elif grep "N =" ${folder}.out.out | grep 500 >/dev/null 2>&1; then
      printf "${YELLOW}${folder} run complete in > 500 steps... This means your model did not converge very quicly.${NC}"
      echo ""
      echo "See ${folder}.out.out for model output to diagnose this."
      echo "Please DO NOT submit your code changes until resolving this issue."
    else
      printf "${GREEN}${folder} run complete in < 500 steps! Output stored in ${folder}.out.out${NC}"
      echo ""
    fi
  else
    # In this case, the model crashed before completion.
    printf "${RED}${folder} run did NOT complete! Something is worng. Quitting.${NC}"
    echo "See ${folder}.out.out for model output and error messages."
    echo "Please DO NOT submit your code changes until resolving this issue."
    echo ""
    exit 1
  fi
  echo ""

  # *** Test photochem outputs ***
  # FOR LATER - CAN PROBABLY USE NSTEPS and L2 ERROR FOR THIS.
  # GIADA PROBABLY KNOWS SOME TESTS WE CAN USE FROM HER COUPLING EXPERTISE.

#  # *** CLIMA tests ***

#  # copy files
#  cd $true_path
#  if cd $folder_path2; then
#    if cp 'input_clima.dat' '../..'; then
#      cd ${true_path}
#  fi
#  cd $true_path
#  if [ ! -d "$folder_path2" ]; then
#    cd $temp_path3
#    echo "!!WARNING!!: no clima template exists for $folder. Copying GENERIC version."
#    echo "This will *probably* cause problems for you. Check the file manually to see if you want the options set in it."
#    cp 'DEFAULT_input_clima.dat' 'input_clima.dat' && echo "Copied GENERIC input_clima.dat to  $(pwd ../..)"
#  fi
#  echo 'Finished copying clima templates over'
#  cd $true_path
#  pwd

#  # compile climate code
#  make '-f' 'ClimaMake' 'clean'

#  # run climate code
#  ./'Clima.run'

done

printf "${GREEN}Everything compiled and ran!${NC}"
echo ""
printf "Look above for ${RED}red text${NC} and ${YELLOW}yellow text${NC}, to make sure all models converged."
echo ""
echo "This script does not test for accuracy."
echo "But... at least your templates compiled and ran! (And converged?)"
echo ""
