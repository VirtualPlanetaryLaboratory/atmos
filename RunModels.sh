#!/bin/bash

temp_path='PHOTOCHEM/INPUTFILES/TEMPLATES/'
temp_path2='CLIMA/IO/TEMPLATES/'
temp_path3='CLIMA/IO/'

clear
echo "This script imports PHOTOCHEM & CLIMA templates from a template file."
echo "You will need to be in your top folder for atmos, not a subfolder."
echo "Your current folder is $(pwd)."
echo "Your PHOTOCHEM/INPUTFILES/TEMPLATE folder has these templates:"
ls $temp_path --hide=*.*
true_path=$(pwd)
echo -n "Enter the folder title (NOT THE PATH) and press [ENTER]: "
read folder

folder_path=$temp_path$folder
folder_path2=$temp_path2$folder

photo_ins='PHOTOCHEM/INPUTFILES/'
clima_ins='CLIMA/IO'

#cd $folder_path

# If the / is at the end, get rid of it
if [ "${folder_path:$length:-1}" == '/' ]
    then
    folder_path="${folder_path:-1}"
fi


# Copies over the newer dist file if it's present, 
# otherwise copies over the normal in.dist
cp "$folder_path/in.dist" "$photo_ins/../in.dist" && echo "Copied in.dist from $folder_path"

cp "$folder_path/input_photchem.dat" $photo_ins && echo "Copied input_photchem.dat from $folder_path"
cp "$folder_path/reactions.rx" $photo_ins && echo "Copied reactions.rx from $folder_path"
cp "$folder_path/parameters.inc" $photo_ins && echo "Copied parameters.inc from $folder_path"
cp "$folder_path/species.dat" $photo_ins && echo "Copied species.dat from $folder_path"
cp "$folder_path/PLANET.dat" $photo_ins && echo "Copied PLANET.dat from $folder_path"
echo 'Finished copying photochem templates over'
cd $true_path
pwd

echo -e '\n'
echo -n 'Would you like to compile PHOTOCHEM (y/n)?:'
 read recompile
  if [ "$recompile" == "y" -o "$recompile" == 'Y' ]
   then
       make '-f' 'PhotoMake' 'clean'
       make '-f' 'PhotoMake'
  fi
echo -e '\n'
echo -n 'Would you like to run PHOTOCHEM (y/n)?:'
  read run_photo
  if [ "$run_photo" == "y" -o "$run_photo" == 'Y' ]
   then             
       ./'Photo.run'
  fi
echo -e '\n'
echo -n 'Would you like to recompile CLIMA (y/n)?:'
read recompile
 if [ "$recompile" == "y" -o "$recompile" == 'Y' ]
 then
     make '-f' 'ClimaMake' 'clean'
     make '-f' 'ClimaMake'
 fi 
echo -e '\n'
echo -n 'Would you like to run CLIMA (y/n)?:'
read run_clima
if [ "$run_clima" == "y" -o "$run_clima" == 'Y' ]
then
    cd $true_path
    echo "folder path is $folder_path2"
	cp "$folder_path2/input_clima.dat" $clima_ins && echo "Copied input_clima.dat from $folder_path2"
	if [ "$run_photo" == "n" -o "$run_photo" == "N" ]
	then
         echo "NOTE: You did not run photo but you want to run Clima. Copying stored Clima coupling files."
         cp "$folder_path2/coupling_params.out" "COUPLE/" && echo "Copied coupling_params.out from $folder_path2"
         cp "$folder_path2/fromPhoto2Clima.dat" "COUPLE/" && echo "Copied fromPhoto2Clima.dat from $folder_path2"
         cp "$folder_path2/hcaer.photoout.out" "COUPLE/" && echo "Copied hcaer.photoout.out from $folder_path2"
         cp "$folder_path2/mixing_ratios.dat" "COUPLE/" && echo "Copied mixing_ratios.dat from $folder_path2"
         cp "$folder_path2/time_frak_photo.out" "COUPLE/" && echo "Copied time_frak_photo.out from $folder_path2"
    fi
 
     
        cd $true_path
    if [ ! -d "$folder_path2" ]; then
	echo "folder path is now $(pwd)"
	cd $temp_path3
	echo "!!WARNING!!: no clima template exists for $folder. Copying GENERIC version."
	echo "This will *probably* cause problems for you. Check the file manually to see if you want the options set in it."
	cp 'DEFAULT_input_clima.dat' 'input_clima.dat' && echo "Copied GENERIC input_clima.dat to  $(pwd ../..)"
    fi
    #cp 'mixing_ratios.dat' '../../../../COUPLE/' && echo "Copied mixing_ratios.dat $(pwd ../../../../COUPLE/)"
    #cp 'mixing_ratios.dat' '../..' && echo "Copied mixing_ratios.dat $(pwd ../..)"
    #cp 'time_frak_photo.out' '../../../../COUPLE/' && echo "Copied time_frak_photo.out to $(pwd ../../../../COUPLE/)"
    #cp 'hcaer.photoout.out' '../../../../COUPLE/' && echo "Copied hcaer.photoout.out to $(pwd ../../../../COUPLE/)"
    echo 'Finished copying clima templates over'
    cd $true_path
    pwd

  ./'Clima.run'
fi
