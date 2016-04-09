#combine files for the multiprocessor run

import os
import sys
import shutil

#base directory
baseDir = sys.argv[1]

#what files to combine, Points or PILE_HEIGHT
combineFileTypes = ["Points", "PILE_HEIGHT"]

#output directory
baseDir = baseDir+"/script_output"

#keep list of all the time stamps
timeStampArray = list()

#custom list dir that sorts directory correctly
def custom_listdir(path):
    """
    Returns the content of a directory by showing directories first
    and then files by ordering the names alphabetically
    """
    dirs = sorted([d for d in os.listdir(path) if os.path.isdir(path + os.path.sep + d)])
    dirs.extend(sorted([f for f in os.listdir(path) if os.path.isfile(path + os.path.sep + f)]))
    return dirs


for combine in combineFileTypes:
    del timeStampArray[:]

    #dirList = os.listdir(baseDir)
    dirList = custom_listdir(baseDir)
    for fname in dirList:
        fileSplit = fname.split('.')
        #check to see if its a file
        if len(fileSplit) > 1:
            #check to see if its a 'Points' file
            if fileSplit[1] == combine:
                #get timestamp of current file
                fileTimeStamp = fileSplit[0][10:]
                #check if there is a new time stamp
                if fileTimeStamp not in timeStampArray:
                    #print "adding: " + fileTimeStamp + " " + combine
                    #update time stamp
                    timeStampArray.append(fileTimeStamp)
                    #file to write to
                    combinef = open(baseDir+"/xdmf000000"+fileTimeStamp+"."+combine+"C", "w")
                    for fname2 in dirList:
                        fileSplit2 = fname2.split('.')
                        #check to see if its a file
                        if len(fileSplit2) > 1:
                            #check to see if its a 'Points' file
                            if fileSplit2[1] == combine:
                                fileTimeStamp2 = fileSplit2[0][10:]
                                if fileTimeStamp2 == fileTimeStamp:
                                    #print fileSplit2[0] + combine
                                    tempf = open(baseDir+"/"+fileSplit2[0]+"."+combine, "r")
                                    for line in tempf:
                                        combinef.write(line)

os.system("rm "+baseDir+"/*.Points "+baseDir+"/*.PILE_HEIGHT")

# This form of rename not working consistently across Linux distributions
#os.system("rename -v 's/\.PointsC$/\.Points/' "+baseDir+"/*.PointsC")
#os.system("rename -v 's/\.PILE_HEIGHTC$/\.PILE_HEIGHT/' "+baseDir+"/*.PILE_HEIGHTC")

for fname in os.listdir(baseDir):
    fileSplit = fname.split('.')
    #print fileSplit
    if fileSplit[1] == combineFileTypes[0]+"C":
       shutil.move(baseDir+"/"+fname,baseDir+"/"+fileSplit[0]+"."+combineFileTypes[0])
    if fileSplit[1] == combineFileTypes[1]+"C":
       shutil.move(baseDir+"/"+fname,baseDir+"/"+fileSplit[0]+"."+combineFileTypes[1])
      
