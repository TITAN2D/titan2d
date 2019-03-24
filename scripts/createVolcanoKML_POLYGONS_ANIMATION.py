
# File:		createVolcanoKML.py
# Author: 	Kyle Marcus
# Date:   	Jan 27 2012

# Takes the output of Titan and creates the kml file that will plot the 
# output points on Google Earth.

# This file will create the kml file that will show the animation of the flow
# of a volcano from the output of the Titan simulation. It uses the time-step
# animation in GOOGLE EARTH to display the animation (NOTE: this will not work
# in GOOGLE MAPS).  The file output will be very large and to lower the file 
# size, you can compress it to a kmz file. To do this, change the output of the
# file created by this program to doc.kml and then zip this file up and change 
# the name to foo.kmz. The you can open this in GOOGLE EARTH and it will have a
# lower file size.

# HOW TO CREATE KMZ FILE:
# $mv merapi_polygons_animation.kml doc.kml
# $zip merapi.kmz doc.kml
# $google-earth merapi.kmz &

# HOW TO RUN:
# $python h5dump.py /Users/kyle/CCR/Spring2012/google_earth/lasttimestep/
# $python parseXML.py /Users/kyle/CCR/Spring2012/google_earth/lasttimestep/
# $python combineFiles.py /Users/kyle/CCR/Spring2012/google_earth/lasttimestep/
# $python createVolcanoKML_POLYGONS_ANIMATION.py /Users/kyle/CCR/Spring2012/google_earth/lasttimestep/ $TOOLDIR/bin

# TODO:
# - (DONE) Automatic way to get Points and PILE_HEIGHT from h5 file
# - (DONE) Take away points that PILE_HEIGHT !> 0
# - (DONE) Add time stepped animation

from utmToLatLng import utmToLatLng
from createLabel import makeLabel
import string
import datetime
import os
import sys
import shutil
import zipfile

zone = 13
northernHemisphere = True

header = '''<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
	
	<Document>
	
		<name>Pile Height</name>
		<open>0</open>
		
		<Style id="volcanoPolyStyle0">
			<PolyStyle>
				<color>a0E7551C</color>
				<outline>0</outline>
			</PolyStyle>
		</Style>
		
		<Style id="volcanoPolyStyle1">
			<PolyStyle>
				<color>a0B67950</color>
				<outline>0</outline>
			</PolyStyle>
		</Style>
		
		<Style id="volcanoPolyStyle2">
			<PolyStyle>
				<color>a0869D85</color>
				<outline>0</outline>
			</PolyStyle>
		</Style>
		
		<Style id="volcanoPolyStyle3">
			<PolyStyle>
				<color>a055C0B9</color>
				<outline>0</outline>
			</PolyStyle>
		</Style>
		
		<Style id="volcanoPolyStyle4">
			<PolyStyle>
				<color>a02CDEE4</color>
				<outline>0</outline>
			</PolyStyle>
		</Style>
		
		<Style id="volcanoPolyStyle5">
			<PolyStyle>
				<color>a014D5FF</color>
				<outline>0</outline>
			</PolyStyle>
		</Style>
		
		<Style id="volcanoPolyStyle6">
			<PolyStyle>
				<color>a014A8FF</color>
				<outline>0</outline>
			</PolyStyle>
		</Style>
		
		<Style id="volcanoPolyStyle7">
			<PolyStyle>
				<color>a01473FF</color>
				<outline>0</outline>
			</PolyStyle>
		</Style>
		
		<Style id="volcanoPolyStyle8">
			<PolyStyle>
				<color>a0143EFF</color>
				<outline>0</outline>
			</PolyStyle>
		</Style>
		
		<Style id="volcanoPolyStyle9">
			<PolyStyle>
				<color>a01408FF</color>
				<outline>0</outline>
			</PolyStyle>
		</Style>

		<Style id="volcanoPolyStyle10">
			<PolyStyle>
				<color>a01400A0</color>
				<outline>0</outline>
			</PolyStyle>
		</Style>

		<ScreenOverlay>
			<name>Pile Height Legend</name>
			<Icon> 
				<href>pileheight.png</href>
			</Icon>
			<overlayXY x="0" y="1" xunits="fraction" yunits="fraction" />
			<screenXY x="0" y="1" xunits="fraction" yunits="fraction" />
			<rotationXY x="0" y="0" xunits="fraction" yunits="fraction" />
			<size x="0" y="0" xunits="fraction" yunits="fraction" />
		</ScreenOverlay>
		
'''

footer = '''
	
	</Document>
	
</kml>
'''

placemarkHeader = '''
			<Polygon>
				<extrude>1</extrude>
				<altitudeMode>clampToGround</altitudeMode>
				<outerBoundaryIs>
					<LinearRing>
						<coordinates>
'''

placemarkFooter = '''
						</coordinates>
					</LinearRing>
				</outerBoundaryIs>
			</Polygon>
		</Placemark>
'''

# custom list dir that sorts directory correctly
def custom_listdir(path):
    """
    Returns the content of a directory by showing directories first
    and then files by ordering the names alphabetically
    """
    dirs = sorted([d for d in os.listdir(path) if os.path.isdir(path + os.path.sep + d)])
    dirs.extend(sorted([f for f in os.listdir(path) if os.path.isfile(path + os.path.sep + f)]))
    return dirs

# base directory
outDir = sys.argv[1]
fontDir = sys.argv[2]+"/"
baseDir = outDir+"/script_output"

print ("fontDir: %s" % fontDir)
print ("baseDir: %s" % baseDir)

# Output kml (Google Earth) file
merapi = open(outDir+'/pileheight_animation.kml', 'w')

#get Zone
zonefile = open(outDir+"/zone.txt","r")
zone = zonefile.readline()
zone = zone.strip()
zone = int(zone)
hemi = zonefile.readline()
zonefile.close()

hemi = hemi.strip()
if(hemi == "North"):
   northernHemisphere = True
if(hemi == "South"):
   northernHemisphere = False


#get the maximum height
maxHeightFile = open(outDir+"/maxHeight_for_KML.data","r")
maxHeight = maxHeightFile.readline()
maxHeight = maxHeight.strip()
#maxHeight = int(float(maxHeight))
maxHeight = float(maxHeight)

#get Image location
#imagefile = os.path.realpath(__file__)
#imagefile = string.replace(imagefile, "createVolcanoKML_POLYGONS_ANIMATION.py", "pileheight.png")
#print imagefile, " to ", outDir
#shutil.copy(imagefile, outDir)

#Make label
blankimage = os.path.realpath(__file__)
blankimage = string.replace(blankimage, "createVolcanoKML_POLYGONS_ANIMATION.py", "blankheight.png")
makeLabel(blankimage, fontDir, outDir + "/pileheightlabel.png", maxHeight)

# ----- MAIN ------

date = datetime.date(2010,1,1)

merapi.write(header)

dirList = custom_listdir(baseDir)

if (len(dirList) > 0):
    print ("file(s) found:")

for fname in dirList:
    fileSplit = fname.split('.')
    # check to see if its a file
    if len(fileSplit) > 1:
		# check to see if its a 'Points' file
		if fileSplit[1] == "Points":
			print fname
			points = open(baseDir+"/"+fname, 'r')
			ph = fileSplit[0]+".PILE_HEIGHT"
			heights = open(baseDir+"/"+ph, 'r')

			# every 1 line in PILE_HEIGHT corresponds to 4 lines in Points.txt
			height = heights.readline()
			begin = date
			end = date + datetime.timedelta(days=1)
			firstCordX = 0
			firstCordY = 0
			firstCordZ = 0
			while 1:
				# Convert UTM to lat/lon and then print <Placemark> tag
				if ( float(height) > 0 ):
					merapi.write('\t\t<Placemark>\n')
					merapi.write('\t\t\t<name>Height '+str(float(height)) +' meters</name>\n')
					style = int(float(height)*10/maxHeight)
					if(style > 10):
						style = 10
					merapi.write('\t\t\t<styleUrl>#volcanoPolyStyle'+str(style)+'</styleUrl>\n')
					merapi.write('\t\t\t<TimeSpan>\n')
					merapi.write('\t\t\t\t<begin>'+str(begin)+'</begin>\n')
					merapi.write('\t\t\t\t<end>'+str(end)+'</end>\n')
					merapi.write('\t\t\t</TimeSpan>\n')
					merapi.write(placemarkHeader)
                                        tryCount=0
                                        exceptCount=0
					for x in range(5):
						if (x != 4):
							point = points.readline()
							xyz = point.split('\t')
                                                        try:
                                                           tryCount+=1
							   (lat,lng) = utmToLatLng(zone, float(xyz[0]), float(xyz[1]), northernHemisphere)
							   merapi.write('\t\t\t\t\t\t\t'+str(lng)+','+str(lat)+','+str(float(height))+'\n')
                                                        except:
                                                           exceptCount+=1
                                                           break
							if (x == 0):
								firstCordX = str(lng)
								firstCordY = str(lat)
								firstCordZ = str(float(height))
						if (x == 4):
							merapi.write('\t\t\t\t\t\t\t'+firstCordX+','+firstCordY+','+firstCordZ)
					merapi.write(placemarkFooter)
				else:
					for x in range(4):
						point = points.readline()
				height = heights.readline()
				if not height:
					break
                        # end while
			date += datetime.timedelta(days=1)
			points.close()
			heights.close()

                        if exceptCount != 0:
                            print("try count: %d\n", tryCount)
                            print("except count: %d\n", exceptCount)
                            print("percent: %d\n", exceptCount/tryCount)

merapi.write(footer)
merapi.close()

zf = zipfile.ZipFile(outDir+"/pileheightAnimation.kmz", mode="w")
zf.write(outDir+"/pileheightlabel.png", "pileheight.png")
zf.write(outDir+"/pileheight_animation.kml", "pileheight_animation.kml")
zf.close()
os.remove(outDir+"/pileheight_animation.kml")
os.remove(outDir+"/pileheightlabel.png")
#os.remove(outDir+"/pileheight.png")
