#!/bin/bash

# call literlly like ./particle_format.sh 100 200

#----------------------------------------------------------------------
# Extract particle data and format for use with the climate model
#----------------------------------------------------------------------

if [ $# -ne 2 ]; then
	echo 'Usage:', $0, '<fCH4 ppmv> <fCO2 ppmv>'
	exit 127
fi

wd='/astro/users/giada/archean_earth/cloud_files/hcaer/'
#filename=$wd'/'$1'ppmCH4_'$2'ppmCO2.hcaer.dat'
filename=$wd'hcaer_new_fractal_'$1'CH4.out'
#if [ -f $filename ]; then
#  echo -e $filename'\n\n'
#else
#	echo -e "Errah: $filename does not exist"#
#	exit 127
#fi

indices=( 4 8 12 16 20 24 28 32 36 40 44 48 52 56 60 64 68 72)
npart='      DATA NPART/'
rpart='      DATA RPART/'

for entry in ${indices[@]}
do
	numcmd=$entry'p'
	nval=`sed -n $numcmd $filename | awk '{print $2}'`
	rval=`sed -n $numcmd $filename | awk '{print $3}'`
	if [ $entry == 20 ] || [ $entry == 40 ] || [ $entry == 60 ]; then
		npart=$npart'\n     &  '$nval', '
		rpart=$rpart'\n     &  '$rval', '
	elif [ $entry == 72 ]; then
		npart=$npart$nval'/'
		rpart=$rpart$rval'/'
	else
		npart=$npart$nval', '
		rpart=$rpart$rval', '
	fi
done

echo -e "C*** FCH4=$1ppm, FCO2=$2ppm ***************************"
echo -e "      print *,'using particle CH4="$1"ppm, CO2="$2"ppm'"
echo -e "$npart"
echo -e ''
echo -e "$rpart"
echo 'C****************************************************************'
