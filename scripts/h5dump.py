#This will dump Points and PILE_HEIGHT from h5 file to text file

import os
import sys

baseDir = sys.argv[1]
directory = baseDir + "/script_output"

#output directroy that everything will go in to make the animation and images
if not os.path.exists(directory):
   os.makedirs(directory)
else:
   os.system("rm "+directory+"/*.*")

dirList = os.listdir(baseDir)
for fname in dirList:
    fileSplit = fname.split('.')
    #print ("fileSplit: ", fileSplit)
    #check to see if its a file
    if len(fileSplit) > 1:
        #check to see if its a h5 file
        if fileSplit[1] == "h5":
            #create xml file in new output directory script_ouput/
            os.system("h5dump -x "+baseDir+"/"+fname+" > "+baseDir+"/script_output/"+fileSplit[0]+".xml")
