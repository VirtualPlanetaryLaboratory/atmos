# Unix commands to save the output of a flare simulation 

set SDIR = SAVEOUT/GJ581d-2.45bar-CO290-SolCon1.3
mkdir $SDIR
#cp inputGl581-1AUeq.dat  $SDIR

cp IO/*.tab $SDIR
cp IO/*.dat $SDIR

cd CLIMA
cp IO/TempIn.dat ../$SDIR

cd ../PHOTCHEM
cp IO/*.dat ../$SDIR
cp photo_CO2.f ../$SDIR

cd ..

