# this program will extract the lava flow info from Points.txt
# this will show just the lava flow

# this will need 3 points added to fit with base map


import os
import sys

baseDir = sys.argv[1]

#output directroy
baseDir = baseDir+"/script_output"



dirList = os.listdir(baseDir)
for fname in dirList:
    fileSplit = fname.split('.')
    #check to see if its a file
    if len(fileSplit) > 1:
        #check to see if its a PILE_HEIGHT
        if fileSplit[1] == "PILE_HEIGHT":
        
        
            #open points
            pointsf = open(baseDir+"/"+fileSplit[0]+'.Points', 'r')
            #open pile_height
            pilef = open(baseDir+"/"+fileSplit[0]+'.PILE_HEIGHT', 'r')
            #open lava file to write
            lavaf = open(baseDir+"/"+fileSplit[0]+'.lava', 'w')

            #open boundaries file
            boundaries = open(baseDir+"/boundaries.txt", 'r')
        
            #read all line boundaries into lava file
            lavaf.write(boundaries.read())
        
            boundaries.close()
            
            pilel = pilef.readline()

            while pilel:
                if (pilel == '0\n'):
                    pointsf.readline()
                    pointsf.readline()
                    pointsf.readline()
                    pointsf.readline()
                else:
                    line1 = pointsf.readline().split()
                    x1 = float(line1[0])
                    y1 = float(line1[1])
                    z1 = float(line1[2])
                    line2 = pointsf.readline().split()
                    x2 = float(line2[0])
                    y2 = float(line2[1])
                    z2 = float(line2[2])
                    line3 = pointsf.readline().split()
                    x3 = float(line3[0])
                    y3 = float(line3[1])
                    z3 = float(line3[2])
                    line4 = pointsf.readline().split()
                    x4 = float(line4[0])
                    y4 = float(line4[1])
                    z4 = float(line4[2])
                    #get center points
                    center_x = ((x2-x1)/2)+x1
                    center_y = ((y3-y2)/2)+y2
                    #get average height
                    center_z = ((z1+z2+z3+z4)/4)+float(pilel)
                    #print
                    lavaf.write(str(center_x)+"\t"+str(center_y)+"\t"+str(center_z)+"\n\n")
            
                pilel = pilef.readline()

            pointsf.close()
            pilef.close()
            lavaf.close()