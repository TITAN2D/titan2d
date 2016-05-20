#!/bin/bash

wd=`pwd`
titan=$wd/../../bin/titan

#source enviroment variables if needed
if [ -a "$wd/../../bin/titanvar" ]; then
    echo "loading titan2d enviroment variables"
    source "$wd/../../bin/titanvar"
fi

for d in Coulomb  Coulomb_MatMap
do 
echo $d
cd inclined/$d
echo "In directory: " `pwd`
echo "Running:" $titan -nt 4 input.py \>\& out_reg
$titan -nt 4 input.py >& out_reg
cd $wd
done

for d in Coulomb Pouliquen_Forterre Voellmy_Salm Coulomb_with_stopping Coulomb_all_param_explicit Coulomb_FluxSrc_DischPlane Coulomb_GeoTIFF TwoPhases_Pitman_Le
do 
echo $d
cd colimafinemini/$d
echo "In directory: " `pwd`
echo "Running:" $titan -nt 4 input.py \>\& out_reg
$titan -nt 4 input.py >& out_reg
cd $wd
done