#this will create the gnuplot scripts and then run them

import os
import sys

if len(sys.argv) == 2:
    baseDir = sys.argv[1]
else:
    baseDir = sys.argv[1]
    localDir = sys.argv[2]

#output directroy that everything will go in to make the animation and images
baseDirOld = baseDir
baseDir = baseDir+"/script_output"



#this will read Stdout.txt and put the timestamps in a dictionary, Key:timestamp | Value:time
def timestamp():
    stdoutf = open(baseDirOld+"/Stdout.txt")
    line = stdoutf.readline()
    while line != "":
        if len(line) > 25:
            if line[0:23] == "At the end of time step":
                arr = line[24:].split(" ")
                global timestamps
                simTime = arr[4].split(":")
                #FORMAT the time hh:mm:ss.xx
                if len(simTime[0]) == 1:         #check to make sure there are 2 numbers in the hours place
                    simTime[0] = "0"+simTime[0]
                seconds = simTime[2]             #check seconds
                seconds = seconds.split(".")
                if len(seconds) == 2:  #if there is a period
                    if len(seconds[0]) != 2:
                        seconds[0] = "0"+seconds[0]
                    tmp = seconds[1]
                    seconds[1] = tmp[:2]
                    if len(seconds[1]) != 2:
                        seconds[1] = seconds[1]+"0"
                    simTime[2] = seconds[0]+"."+seconds[1]
                else:                  #if there is no period
                    if len(seconds[0]) == 1:
                        seconds[0] = "0"+seconds[0]
                    simTime[2] = seconds[0]+".00"

                timestamps[arr[0]] = simTime[0]+":"+simTime[1]+":"+simTime[2]
        line = stdoutf.readline()

timestamps = {} #dictionary that holds the time stamps
timestamp() #populate the dictionary




dirList = os.listdir(baseDir)
for fname in dirList:
    basename, extension = fname.split('.')
    if extension == "lava":
        
        # contour line images
        
        gnuplotf = open(baseDir+"/"+basename+".gnuplot", 'w')
        
        zrange = open(baseDir+"/zrange.txt", 'r')
        z_range = zrange.readline()
        zrange.close()

        zMap = "set cbrange [0:"+z_range+"]\n"
        z3D  = "set zrange [0:"+z_range+"]\n"
        
        gnuplotf.write("set term gif\n")
        gnuplotf.write("set output '"+baseDir+"/"+basename+".gif'\n")
        gnuplotf.write("set dgrid 30,30,4\n")
        gnuplotf.write("set view map\n")
        gnuplotf.write("set palette gray\n")
        gnuplotf.write("set multiplot\n")
        title = "set title '"+basename
        if str(int(basename[12:])) in timestamps:
            title = title + " ~ " + timestamps[str(int(basename[12:]))] + "'\n"
        else:
            title = title + "'\n"
        gnuplotf.write(title)
        #gnuplotf.write("set title '"+basename+"'\n")
        #gnuplotf.write("set cbrange [0:4000]\n")
        gnuplotf.write(zMap)
        gnuplotf.write("splot '"+baseDir+"/"+"xdmf000000_i00000000.Points' w pm3d notitle\n")
        gnuplotf.write("unset dgrid\n")
        gnuplotf.write("splot '"+baseDir+"/"+basename+".lava' w p notitle\n")
        gnuplotf.write("set dgrid 30,30,4\n")
        gnuplotf.write("unset surface\n")
        gnuplotf.write("set contour\n")
        gnuplotf.write("set cntrparam levels 20\n")
        gnuplotf.write("set key\n")
        gnuplotf.write("splot '"+baseDir+"/"+basename+".lava_points' w pm3d notitle\n")
        
        gnuplotf.close()
        
        # 3D images
        
        gnuplotf = open(baseDir+"/"+basename+".gnuplot3D", 'w')
        
        gnuplotf.write("set term gif\n")
        gnuplotf.write("set output '"+baseDir+"/"+basename+".gif3D'\n")
        gnuplotf.write("set dgrid 30,30,4\n")
        gnuplotf.write("set style data pm3d\n")
        gnuplotf.write("unset key\n")
        gnuplotf.write("set multiplot\n")
        title = "set title '"+basename
        if str(int(basename[12:])) in timestamps:
            title = title + " ~ " + timestamps[str(int(basename[12:]))] + "'\n"
        else:
            title = title + "'\n"
        gnuplotf.write(title)
        #gnuplotf.write("set title '"+basename+"'\n")
        #gnuplotf.write("set zrange [0:4000]\n")
        gnuplotf.write(z3D)
        gnuplotf.write("set palette gray\n")
        gnuplotf.write("splot '"+baseDir+"/"+"xdmf000000_i00000000.Points'\n")
        gnuplotf.write("set palette color\n")
        gnuplotf.write("unset dgrid\n")
        gnuplotf.write("splot '"+baseDir+"/"+basename+".lava' w p\n")
        
        gnuplotf.close()
        
        os.system("gnuplot "+baseDir+"/"+basename+".gnuplot")
        os.system("gnuplot "+baseDir+"/"+basename+".gnuplot3D")

os.system("convert -delay 100 "+baseDir+"/"+"*.gif -loop 0 "+baseDirOld+"/"+"animation.gif")
os.system("convert -delay 100 "+baseDir+"/"+"*.gif3D -loop 0 "+baseDirOld+"/"+"animation3D.gif")

if len(sys.argv) == 3:
    os.system("cp "+baseDirOld+"/"+"animation.gif "+localDir)
    os.system("cp "+baseDirOld+"/"+"animation3D.gif "+localDir)

#os.system("rm -r "+baseDir)
