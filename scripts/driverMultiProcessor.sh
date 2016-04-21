#!/bin/bash

# Args:
# $1 - Run Directory
# $2 - TOOLDIR

cp $1/Stdout.txt $1/vizout/Stdout.txt

echo "Running h5dump"
python $2/bin/h5dump.py $1/vizout
echo "Parsing XML"
python $2/bin/parseXML.py $1/vizout
echo "Combining Files"
python $2/bin/combineFiles.py $1/vizout
echo "Running boundariesMultiProcessor"
python $2/bin/boundariesMultiProcessor.py $1/vizout
echo "Running lava"
python $2/bin/lava.py $1/vizout
echo "Running lava_points"
python $2/bin/lava_points.py $1/vizout

echo "Running gnuplot.py"
if [ $# == 2 ] ; then
    python $2/bin/gnuplot.py $1/vizout
else
    python $2/bin/gnuplot.py $1/vizout $3
fi
