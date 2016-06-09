#!/bin/bash

# Args:
# $1 - Run Directory
# $2 - TOOLDIR (not used)

fullpath="$(readlink -f $0)"
TITAN2D_HOME="$(echo "$fullpath" | sed "s?/bin/driverKML.sh??")"
TITAN2D_PYUTIL_DIR=$TITAN2D_HOME/lib/python2.7/site-packages
TITAN2D_PYTHON=$TITAN2D_HOME/bin/titan


cp $1/zone.txt $1/vizout/zone.txt
cp $1/height_scale_for_KML.data $1/vizout/height_scale_for_KML.data
 
echo "Running h5dump"
$TITAN2D_PYTHON $TITAN2D_PYUTIL_DIR/h5dump.py $1/vizout
echo "Parsing XML"
$TITAN2D_PYTHON $TITAN2D_PYUTIL_DIR/parseXML.py $1/vizout
echo "Combining Files"
$TITAN2D_PYTHON $TITAN2D_PYUTIL_DIR/combineFiles.py $1/vizout
# On VHub, FreeSerif.ttf font file is in the 
# /usr/share/fonts/truetype/freefont directory.
# Cannot assume this is true for all installations of the Titan2D GUI.
# Look for FreeSerif.ttf font file in $TOOLDIR/bin directory.
echo "Creating animated KML"
$TITAN2D_PYTHON $TITAN2D_PYUTIL_DIR/createVolcanoKML_POLYGONS_ANIMATION.py $1/vizout $TITAN2D_PYUTIL_DIR
echo "Creating flat KML"
$TITAN2D_PYTHON $TITAN2D_PYUTIL_DIR/createVolcanoKML_POLYGONS_LASTFRAME.py $1/vizout $TITAN2D_PYUTIL_DIR
