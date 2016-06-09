#!/bin/bash

# Args:
# $1 - Run Directory
# $2 - TOOLDIR (not used)

fullpath="$(readlink -f $0)"
TITAN2D_HOME="$(echo "$fullpath" | sed "s?/bin/driverMultiProcessor.sh??")"
TITAN2D_PYUTIL_DIR=$TITAN2D_HOME/lib/python2.7/site-packages
TITAN2D_PYTHON=$TITAN2D_HOME/bin/titan

RUN_DIR=$1/vizout

cp $1/Stdout.txt $1/vizout/Stdout.txt

echo "Running h5dump"
$TITAN2D_PYTHON $TITAN2D_PYUTIL_DIR/h5dump.py $RUN_DIR
echo "Parsing XML"
$TITAN2D_PYTHON $TITAN2D_PYUTIL_DIR/parseXML.py $RUN_DIR
echo "Combining Files"
$TITAN2D_PYTHON $TITAN2D_PYUTIL_DIR/combineFiles.py $RUN_DIR
echo "Running boundariesMultiProcessor"
$TITAN2D_PYTHON $TITAN2D_PYUTIL_DIR/boundariesMultiProcessor.py $RUN_DIR
echo "Running lava"
$TITAN2D_PYTHON $TITAN2D_PYUTIL_DIR/lava.py $RUN_DIR
echo "Running lava_points"
$TITAN2D_PYTHON $TITAN2D_PYUTIL_DIR/lava_points.py $RUN_DIR

echo "Running gnuplot.py"
if [ $# == 2 ] ; then
    $TITAN2D_PYTHON $TITAN2D_PYUTIL_DIR/gnuplot.py $RUN_DIR
else
    $TITAN2D_PYTHON $TITAN2D_PYUTIL_DIR/gnuplot.py $RUN_DIR $3
fi
