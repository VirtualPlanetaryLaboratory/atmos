#!/bin/bash

temp_path='PHOTOCHEM/INPUTFILES/TEMPLATES/'

clear
echo "This script imports a PHOTOCHEM template from a template file."
echo "You will need to be in your top folder for atmos, not a subfolder."
echo "Your current folder is $(pwd)."
true_path=$(pwd)
echo -n "Enter the folder title (NOT THE PATH) and press [ENTER]: "
read folder

folder_path=$temp_path$folder

cd $folder_path
cp 'in.dist' '../../..' && echo "Copied in.dist to $(pwd ../../..)"
cp 'input_photchem.dat' '../..' && echo "Copied input_photchem.dat to $(pwd ../..)"
cp 'ISOreactions.rx' '../..' && echo "Copied ISOreactions.rx to $(pwd ../..)"
cp 'reactions.rx' '../..' && echo "Copied reactions.rx to $(pwd ../..)"
cp 'parameters.inc' '../..' && echo "Copied parameters.inc to $(pwd ../..)"
cp 'species.dat' '../..' && echo "Copied species.dat to $(pwd ../..)"
cp 'PLANET.dat' '../..' && echo "Copied PLANET.dat to $(pwd ../..)"
cp 'ISOparameters.inc' '../..' && echo "Copied ISOparameters.inc to $(pwd ../..)"
cp 'parametersREGULAR.inc' '../..' && echo "Copied parametersREGULAR.in to $(pwd ../..)"
echo 'Finished copying templates over'
cd $true_path
pwd

echo -n 'Would you like to run PHOTOCHEM (y/n)? '
read run_photo
if [ "$run_photo" == "y" -o "$run_photo" == 'Y' ]
   then
       echo -n 'Would you like to recompile the model? (y/n?)'
       read recompile
       if [ "$recompile" == "y" -o "$recompile" == 'Y' ]
       then
	   make '-f' 'makefilePhotochem' 'clean'
	   make '-f' 'makefilePhotochem'
       fi
  ./'TOTCdev'
fi
