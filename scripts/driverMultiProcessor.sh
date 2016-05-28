#!/bin/bash

# Args:
# $1 - Run Directory
# $2 - TOOLDIR (not used)

fullpath="$(readlink -f $0)"
TITAN2D_HOME="$(echo "$fullpath" | sed "s?/bin/driverMultiProcessor.sh??")"
TITAN2D_PYUTIL_DIR=$TITAN2D_HOME/lib/python2.7/site-packages

cp $1/Stdout.txt $1/vizout/Stdout.txt

echo "Running h5dump"
python $TITAN2D_PYUTIL_DIR/h5dump.py $1/vizout
echo "Parsing XML"
python $TITAN2D_PYUTIL_DIR/parseXML.py $1/vizout
echo "Combining Files"
python $TITAN2D_PYUTIL_DIR/combineFiles.py $1/vizout
echo "Running boundariesMultiProcessor"
python $TITAN2D_PYUTIL_DIR/boundariesMultiProcessor.py $1/vizout
echo "Running lava"
python $TITAN2D_PYUTIL_DIR/lava.py $1/vizout
echo "Running lava_points"
python $TITAN2D_PYUTIL_DIR/lava_points.py $1/vizout

echo "Running gnuplot.py"
if [ $# == 2 ] ; then
    python $TITAN2D_PYUTIL_DIR/gnuplot.py $1/vizout
else
    python $TITAN2D_PYUTIL_DIR/gnuplot.py $1/vizout $3
fi
