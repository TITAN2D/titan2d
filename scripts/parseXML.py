import xml.dom.minidom
from xml.dom.minidom import parse
import os
import sys

baseDir = sys.argv[1]

#output directroy that everything will go in to make the animation and images
baseDir = baseDir+"/script_output"

dirList = os.listdir(baseDir)
for fname in dirList:
    basename, extension = fname.split('.')
    if extension == "xml":

        dom = parse(baseDir+"/"+fname)

        #---------------- POINTS --------------------

        pointsfn = basename + ".Points"
        pointsf = open(baseDir+"/"+pointsfn, "w")

        #get all of the connections in the xml file
        connectionsList = dom.getElementsByTagName('hdf5:DataFromFile')[0]
        connections = connectionsList.childNodes[0]
        arrayConnections = connections.nodeValue.split()
        # The connections are indices into the pointsList
        # print('len(arrayConnections): ', len(arrayConnections))
        count = len(arrayConnections)  #the number of connections is 4.
        # print('count: ', count)


        #get all of the points in the xml file
        pointsList = dom.getElementsByTagName('hdf5:DataFromFile')[1]
        points = pointsList.childNodes[0]
        arrayPoints = points.nodeValue.split()
        # print('len(arrayPoints): ', len(arrayPoints))

        #debug = 0

        #for connections in connectionsList.childNodes:

        c = 0  #keep count of how many ponits have been written for the current connection
        j = 0
        while c < count:
            i = int(arrayConnections[c])*3
            #if (debug < 5):
                #print('i: ', i)

            pointsf.write(arrayPoints[i]+"\t"+arrayPoints[i+1]+"\t"+arrayPoints[i+2]+"\n")
            j+=1
            c+=1  #increase connection count
            #debug = debug+1

        pointsf.close()
        #print('j: ', j)

        #for points in pointsList.childNodes:
        #    arrayPoints = points.nodeValue.split()
        #    count = len(arrayPoints)/3  #the number of points, 1 points has x, y, z
        #    i = 0  #iterate through x, y, z for 1 point
        #    c = 0  #keep count of how many points have been written
        #    while c < count:
        #        pointsf.write(arrayPoints[i]+"\t"+arrayPoints[i+1]+"\t"+arrayPoints[i+2]+"\n")
        #        i+=3
        #        c+=1  #increase point count

        pointsf.close()



        #---------------- PILE_HEIGHT --------------------

        #get all of the pile heights in the xml file
        pileHeightList = dom.getElementsByTagName('hdf5:DataFromFile')[2]

        pileHeightfn = basename + ".PILE_HEIGHT"
        pileHeightf = open(baseDir+"/"+pileHeightfn, "w")

        for pileHeight in pileHeightList.childNodes:
            arrayPileHeight = pileHeight.nodeValue.split()
            count = len(arrayPileHeight)  #the number of pile heights
            c = 0  #keep count of how pile heights have been written
            while c < count:
                pileHeightf.write(arrayPileHeight[c]+"\n")
                c+=1  #increase pile height count

        pileHeightf.close()
