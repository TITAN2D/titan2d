# this will write the boundaries of map

import os
import sys

baseDir = sys.argv[1]

#output directroy
baseDir = baseDir+"/script_output"

#open first points file
pointsf = open(baseDir+"/"+"xdmf000000_i00000000.Points", 'r')

save_x = ""
save_first = ""
save_last = ""

greatest_x = 0
greatest_z = 0

for line in pointsf:
    #save first line
    if save_first == "":
        save_first = line
    #save last line
    save_last = line
    #find greatest x
    arr = line.split()
    if arr[0] > greatest_x:
        greatest_x = arr[0]
        save_x = line
    if float(arr[2]) > float(greatest_z):
        greatest_z = arr[2]

pointsf.close()

boundaries = open(baseDir+"/boundaries.txt", 'w')

boundaries.write(save_first+"\n\n"+save_last+"\n\n"+save_x+"\n\n")

boundaries.close()

zrange = open(baseDir+"/zrange.txt", 'w')

#round up to nearest 1000
greatest_z = int( float(greatest_z) / 1000 )
greatest_z += 1
greatest_z *= 1000

zrange.write(str(greatest_z))

zrange.close()
