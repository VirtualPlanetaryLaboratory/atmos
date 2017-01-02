# Unix commands to save the output of a flare simulation 

set SDIR = SAVEOUT/GJ581d-1bar-CO2only-SolCon1.0
mkdir $SDIR
#cp inputGl581-1AUeq.dat  $SDIR

cp IO/clima_allout.tab $SDIR
cp IO/input_clima.dat $SDIR
cp IO/mixing_ratios.dat $SDIR
cp IO/PLANET.dat $SDIR

cd CLIMA
cp IO/TempIn.dat ../$SDIR
cd ..

