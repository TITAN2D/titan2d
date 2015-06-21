"""
Compares two hdf5 files generated for Paraview
"""

import sys
import h5py
import numpy
import math

if len(sys.argv)!=3:
    print "Usage:"
    print "\tpython compare_hdf5.py file1 file2"

filename1=sys.argv[1]
filename2=sys.argv[2]

print "Comparing "+filename1+" and "+filename2

f1 = h5py.File(filename1, "r")
f2 = h5py.File(filename2, "r")

#compare Mash
print "Comparing Mesh/Points"
ds1=numpy.array(f1['Mesh']['Points'])
ds2=numpy.array(f2['Mesh']['Points'])
if ds1.shape!=ds2.shape:
    print "\tDatasets have different size!"
else:
    ds=ds1-ds2
    rmsd=math.sqrt(ds.mean())
    print "\tRMSD=%.3e"%(rmsd,)

#compare Mash
print "Comparing Mesh/Connections"
ds1=numpy.array(f1['Mesh']['Connections'])
ds2=numpy.array(f2['Mesh']['Connections'])
if ds1.shape!=ds2.shape:
    print "\tDatasets have different size!"
else:
    ds=ds1-ds2
    sum=ds.sum()
    if sum>0:
        print "\tConnections are different!"
    else:
        print "\tConnections are identical"

#compare Properties
for prop in ('PILE_HEIGHT', 'XMOMENTUM', 'YMOMENTUM'):
    print "Comparing Properties/"+prop
    
    ds1=numpy.array(f1['Properties'][prop])
    ds2=numpy.array(f2['Properties'][prop])
    if ds1.shape!=ds2.shape:
        print "\tDatasets have different size!"
        continue
    ds=(ds1-ds2)*(ds1-ds2)
    rmsd=math.sqrt(ds.mean())
    print "\tRMSD=%.3e"%(rmsd,)
    


f2.close()
f1.close()