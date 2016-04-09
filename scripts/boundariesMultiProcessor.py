# this will write the boundaries of map

import os
import sys

baseDir = sys.argv[1]

#output directroy
baseDir = baseDir+"/script_output"

#open first points file
pointsf = open(baseDir+"/"+"xdmf000000_i00000000.Points", 'r')

save_x = ""
save_y = ""
save_z = ""

greatest_x = 0
greatest_y = 0
greatest_z = 0

least_x = 0
least_y = 0

for line in pointsf:
    arr = line.split()
    if (least_x == 0):
        least_x = arr[0]
        least_y = arr[1]
    #find least x
    if arr[0] < least_x:
        least_x = arr[0]
    #find greatest x
    if arr[0] > greatest_x:
        greatest_x = arr[0]
    #find y in scientific notation
    split = arr[1].split("e+")
    numberY = float(split[0]) * (10**int(split[1]))
    #print split[0] + " * " + split[1] + " = " + str(numberY)
    #find least y
    if numberY < least_y:
        least_y = numberY
    #find greatest y
    if numberY > greatest_y:
        greatest_y = numberY
    #find greatest z
    if float(arr[2]) > float(greatest_z):
        greatest_z = arr[2]

pointsf.close()

boundaries = open(baseDir+"/boundaries.txt", 'w')

least_y = "%e" % least_y
greatest_y = "%e" % greatest_y

boundaries.write(str(least_x)+"\t"+str(least_y)+"\t"+"0"+"\n\n")
boundaries.write(str(greatest_x)+"\t"+str(least_y)+"\t"+"0"+"\n\n")
boundaries.write(str(greatest_x)+"\t"+str(greatest_y)+"\t"+"0"+"\n\n")
boundaries.write(str(least_x)+"\t"+str(greatest_y)+"\t"+"0"+"\n\n")

boundaries.close()

zrange = open(baseDir+"/zrange.txt", 'w')

#round up to nearest 1000
greatest_z = int( float(greatest_z) / 1000 )
greatest_z += 1
greatest_z *= 1000

zrange.write(str(greatest_z))

zrange.close()
