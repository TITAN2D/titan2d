#!/bin/bash

# Args:
# $1 - Run Directory
# $2 - TOOLDIR

cp $1/zone.txt $1/vizout/zone.txt
cp $1/height_scale_for_KML.data $1/vizout/height_scale_for_KML.data
 
echo "Running h5dump"
python $2/bin/h5dump.py $1/vizout
echo "Parsing XML"
python $2/bin/parseXML.py $1/vizout
echo "Combining Files"
python $2/bin/combineFiles.py $1/vizout
# On VHub, FreeSerif.ttf font file is in the 
# /usr/share/fonts/truetype/freefont directory.
# Cannot assume this is true for all installations of the Titan2D GUI.
# Look for FreeSerif.ttf font file in $TOOLDIR/bin directory.
echo "Creating animated KML"
python $2/bin/createVolcanoKML_POLYGONS_ANIMATION.py $1/vizout $2/bin
echo "Creating flat KML"
python $2/bin/createVolcanoKML_POLYGONS_LASTFRAME.py $1/vizout $2/bin
