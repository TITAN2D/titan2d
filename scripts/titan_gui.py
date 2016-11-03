#!/usr/bin/env python
#*******************************************************************
#* Copyright (C) 2003 University at Buffalo
#*
#* This software can be redistributed free of charge.  See COPYING
#* file in the top distribution directory for more details.
#*
#* This software is distributed in the hope that it will be useful,
#* but WITHOUT ANY WARRANTY; without even the implied warranty of
#* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#*
#* Author: 
#* Description: this is an executable python script
#*
#*******************************************************************
#* $Id: titan_gui.py 211 2009-06-16 20:02:10Z dkumar $ 
#*

import sys,os,math,string,re,socket

import Tkinter
from Tkinter import *
import SimpleDialog,tkMessageBox

class FinalInstructions:

    def __init__(self,master,instruction_string):

        Label(master, text=instruction_string).grid(row=0,column=0,sticky=W,columnspan=3)

        Button(master, text="Quit", command=master.quit).grid(row=3,column=1)

class MixtureProps:
     def __init__ (self):
         self.solid_rho = 2800.
         self.fluid_rho = 1200.
         self.viscosity = 0.1
         self.diameter = 0.005

class QuestionTemplate6:
    def __init__ (self, master, mixprops):
        Label(master, text="Solid-fluid mixture properties.").grid(row=0, column=0, sticky=W, columnspan=2)
        Label(master, text="Solid Density (Kg/m^3)").grid(row=1, column=0, sticky=W)
        Label(master, text="Fluid Density (Kg/m^3)").grid(row=2, column=0, sticky=W)
        Label(master, text="Fluid Viscosity (N-s/m^2)").grid(row=3, column=0, sticky=W)
        Label(master, text="Mean grain diameter (m)").grid(row=4, column=0, sticky=W)

        self.solid_rho = Entry(master)
        self.fluid_rho = Entry(master)
        self.viscosity = Entry(master)
        self.diameter  = Entry(master)

        self.solid_rho.grid(row=1, column=1)
        self.fluid_rho.grid(row=2, column=1)
        self.viscosity.grid(row=3, column=1)
        self.diameter.grid(row=4, column=1)

        Button(master, text="Done", command=self.done).grid(row=5,column=0)
        Button(master, text="Quit", command=master.quit).grid(row=5,column=1)

    def done(self):
        if ( self.solid_rho.get() != '' ):
            mixprops.solid_rho = float(self.solid_rho.get())

        if ( self.fluid_rho.get() != '' ):
            mixprops.fluid_rho = float(self.fluid_rho.get())

        if ( self.viscosity.get() != '' ):
            mixprops.viscosity = float(self.viscosity.get())

        if ( self.diameter.get() != '' ):
            mixprops.diameter = float(self.diameter.get())


class QuestionTemplate5:
    def __init__(self,master,src_number,filename,directory,heightscale):
        Label(master, text="Information for Flux Source Number "+str(src_number+1)).grid(row=0,column=0,sticky=W,columnspan=2)
        Label(master, text="Extrusion flux rate [m/s]:").grid(row=1,column=0,sticky=W)
        Label(master, text="Active Time [s], start, end:").grid(row=2,column=0,sticky=W)
        Label(master, text="Center of the source, xc, yc (UTM E, UTM N):").grid(row=3,column=0,sticky=W)
        Label(master, text="Major and Minor Extent, majorR, minorR (m, m):").grid(row=4,column=0,sticky=W)        
        Label(master, text="Orientation (angle [degrees] from X axis to major axis):").grid(row=5,column=0,sticky=W)
        Label(master, text="Initial speed [m/s]:").grid(row=6,column=0,sticky=W)
        Label(master, text="Initial direction ([degrees] from X axis):").grid(row=7,column=0,sticky=W)

        self.influx = Entry(master)
        self.start_time = Entry(master)
        self.end_time   = Entry(master)
        self.xcenter = Entry(master)
	self.ycenter = Entry(master)
        self.majradius = Entry(master)
        self.minradius = Entry(master)
        self.orientation = Entry(master)
        self.Vmagnitude = Entry(master)
        self.Vdirection = Entry(master)

        self.influx.grid(row=1,column=1)
        self.start_time.grid(row=2,column=1)
        self.end_time.grid(row=2,column=2)
        self.xcenter.grid(row=3,column=1)
	self.ycenter.grid(row=3,column=2)
        self.majradius.grid(row=4,column=1)
        self.minradius.grid(row=4,column=2)
        self.orientation.grid(row=5,column=1)
        self.Vmagnitude.grid(row=6,column=1)
        self.Vdirection.grid(row=7,column=1)

        Button(master, text="Done", command=self.done).grid(row=8,column=0)
        Button(master, text="Quit", command=master.quit).grid(row=8,column=1)

        #passed in variables
        self.master = master
        self.filename = filename
        self.heightscale = heightscale
        self.write_flag = 0 
        self.directory = directory
        self.input_flag = 0 
        self.src_number = src_number

    def done(self):
	# check values
	if self.influx.get() == '':
	    influx = 0.0
	else:
	    influx = float(self.influx.get())
	    if influx <= float(0.):
		influx = 0.

        if self.start_time.get() == '':
            start_time = 0.0
        else:
            start_time = float(self.start_time.get())
            if start_time < 0.0:
                start_time = 0

        if self.end_time.get() == '':
            end_time = 0.0
        else:
            end_time = float(self.end_time.get())
            if end_time < 0:
                end_time = 0

	if self.xcenter.get() == '':
	    xcenter = 1.0
	else:
	    xcenter = float(self.xcenter.get())

	if self.ycenter.get() == '':
	    ycenter = 1.0
	else:
	    ycenter = float(self.ycenter.get())

	if self.majradius.get() == '':
	    majradius = 1.0
	else:
	    majradius = float(self.majradius.get())
	    if majradius <= 0:
		majradius = 1

	if self.minradius.get() == '':
	    minradius = 1.0
	else:
	    minradius = float(self.minradius.get())
	    if minradius <= 0:
		minradius = 1

        if self.orientation.get() == '':
            orientation = 0.0
        else:
            orientation = float(self.orientation.get())

        if self.Vmagnitude.get() == '':
            Vmagnitude = 0.0
            Vdirection = 0.0
        else:
            Vmagnitude=float(self.Vmagnitude.get())
            if self.Vdirection.get() == '':
                Vdirection = 0.0
            else:
                Vdirection=float(self.Vdirection.get())
                            
        #print 'Vmagnitude= ' + str(Vmagnitude) + '\nVdirection= ' + str(Vdirection) + '\n'

        if self.write_flag == 0:
            self.write_flag = 1
            file = open(self.filename, "a+", 0)
            file.write( str(influx) + '\n' + str(start_time) + '\n' + str(end_time) + '\n' + str(xcenter) + '\n' + str(ycenter) + '\n' + str(majradius) + '\n' +str(minradius) + '\n' + str(orientation) + '\n' + str(Vmagnitude) + '\n' + str(Vdirection) + '\n')
            file.close

            #approx: h=influx*t-0.5*a*t^2
            #if no s => t1=N*(2*h/g)^0.5  N is a empirical constant,
            #for cylindrical piles of aspect ratio (height/radius) of approx 1
            #2<=N<=3 (closer to 2) but there are 3 reasons we should increase N
            #(1) cylindrical pile does not collapse the whole way, shorter
            #distance means decreased acceleration means increased time, N
            #(2) paraboloid piles are closer to conical than cylinder so it
            #should collapse even less, so increase N
            #(3) "influx" is a constant source "velocity" not an initial
            #velocity which should increase h in "approx: h=..." equation, so
            #as a fudge factor increase N some more
            #calibrated on a single starting condition at tungaruhau says
            #N=3.21   N=X
            #anyway a=2*h/t1^2 = g/N^2
            #approx: v=influx-a*t2 at hmax v=0 => t2=influx/a = N^2*influx/g
            #t3=min(t2,end_time-start_time)
            #plug int first equation
            #approx hmax=influx*t3-0.5*a*t3^2
            #if t3==t2=> hmax= N^2/2*s^2/g
            #DEM: tungfla2
            #influx 12 m/s (vel 50 m/s at +35 degrees cc from +x direction
            #starts in filled crater which means this velocity points up hill
            #so pile is mostly stationary while flux source is active, 50 m/s
            #is just short of what is needed to top the crater)
            #end_time-start_time=20 gives actual hmax=75.6 m
            #g=9.8 m/s^2, N=3.21, t3=t2=12.62<20 s => computed hmax=75.7 m
            X = 3.21
            g = 9.8
            a = g/X/X
            t3 = X*X*influx/g
            if t3 > (end_time-start_time):
                t3=(end_time-start_time)
            effectheight=influx*t3 - 0.5*a*t3*t3
            if effectheight>self.heightscale:
                self.heightscale = effectheight

            self.input_flag = 1

            #write to the python_input.data file
            file = open(self.directory+'python_input.data', "a+", 0)
            python_data = """ Mean Flux  (kg/(m^2-s)): """+ str(influx) +"""
Active duration of source (start time, end time) (s) " """ + str(start_time) + """ """+ str(end_time) +"""
Center of the Source, xc, yc (UTM E, UTM N): """ + str(xcenter) + """ """+str(ycenter)+"""
Major and Minor Extent, majorR, minorR (m, m): """ + str(majradius)+ """ """ +str(minradius)+"""
Angle from X axis to major axis (degrees): """ +str(orientation)+"""
Initial speed [m/s]: """ + str(Vmagnitude) + """
Initial direction ([degrees] from X axis): """ + str(Vdirection)+"""
"""
            file.write(python_data)
            file.close

class QuestionTemplate4:

    def __init__(self,master,discharge_indice,discharge_planes):

        Label(master, text="Enter Discharge Plane "+str(discharge_indice+1)+" of "+str(discharge_planes)+": ").grid(row=0,column=0,sticky=W,columnspan=2)
        Label(master, text="Point A (UTM E, UTM N): ").grid(row=1,column=0,sticky=W)
        Label(master, text="Point B (UTM E, UTM N): ").grid(row=2,column=0,sticky=W)

        self.x_a = Entry(master)
        self.y_a = Entry(master)
        self.x_b = Entry(master)
        self.y_b = Entry(master)

        self.x_a.grid(row=1,column=1)
        self.y_a.grid(row=1,column=2)
        self.x_b.grid(row=2,column=1)
        self.y_b.grid(row=2,column=2)

        Button(master, text="Done", command=self.done).grid(row=4,column=0)
        Button(master, text="Quit", command=master.quit).grid(row=4,column=1)

    def done(self):
	# check values and store them
        if self.x_a.get() != '':
            self.xa = float(self.x_a.get())
        else:
            self.xa = ''
        if self.x_b.get() != '':
            self.xb = float(self.x_b.get())
        else:
            self.xb = ''
        if self.y_a.get() != '':
            self.ya = float(self.y_a.get())
        else:
            self.ya = ''
        if self.y_b.get() != '':
            self.yb = float(self.y_b.get())
        else:
            self.yb = ''


class QuestionTemplate3:

    def __init__(self,master,nummat,matindice,matname,previntfrict,prevbedfrict):

        self.previntfrict=previntfrict
        self.prevbedfrict=prevbedfrict

        Label(master, text="Material "+str(matindice+1)+" of "+str(nummat)+": "+matname).grid(row=0,column=0,sticky=W,columnspan=2)
        Label(master, text="Internal Friction Angle (deg):").grid(row=1,column=0,sticky=W)
        Label(master, text="Bed Friction Angle (deg):").grid(row=2,column=0,sticky=W)

        self.intfriction = Entry(master)
        self.bedfriction = Entry(master)

        self.intfriction.grid(row=1,column=1)
        self.bedfriction.grid(row=2,column=1)

        Button(master, text="Done", command=self.done).grid(row=4,column=0)
        Button(master, text="Quit", command=master.quit).grid(row=4,column=1)

    def done(self):
	# check values and store them
	if self.intfriction.get() == '':
            self.intfrict = self.previntfrict
	else:
	    self.intfrict = float(self.intfriction.get())
        if self.intfrict <= 0.0 :
            self.intfrict = 30.0

        if self.bedfriction.get() == '':
            self.bedfrict = self.prevbedfrict
	else:
	    self.bedfrict = float(self.bedfriction.get())
        if self.bedfrict <= 0.0 :
            self.bedfrict = 15.0


class QuestionTemplate2:

    def __init__(self,master,pile_number,filename,directory,max_height,topomap):
        
        Label(master, text="Information for Pile Number "+str(pile_number+1)).grid(row=0,column=0,sticky=W,columnspan=2)
        Label(master, text="Thickness of Initial Volume, h(x,y):").grid(row=1,column=0,sticky=W)
        Label(master, text="P*(1-((x-xc)/xr)^2 - ((y-yc)/yr)^2)").grid(row=1,column=1,sticky=W)
        Label(master, text="Maximum Initial Thickness, P (m):").grid(row=2,column=0,sticky=W)
        Label(master, text="Initial solid-volume fraction,(0:1.):").grid(row=3,column=0,sticky=W)
        Label(master, text="Center of Initial Volume, xc, yc (UTM E, UTM N):").grid(row=4,column=0,sticky=W)
        Label(master, text="Major and Minor Extent, majorR, minorR (m, m):").grid(row=5,column=0,sticky=W)        
        Label(master, text="Orientation (angle [degrees] from X axis to major axis):").grid(row=6,column=0,sticky=W)
        Label(master, text="Initial speed [m/s]:").grid(row=7,column=0,sticky=W)
        Label(master, text="Initial direction ([degrees] from X axis):").grid(row=8,column=0,sticky=W)

        self.pileheight = Entry(master)
        self.volfract   = Entry(master)
        self.xpilecenter = Entry(master)
	self.ypilecenter = Entry(master)
        self.majradius = Entry(master)
        self.minradius = Entry(master)
        self.orientation = Entry(master)
        self.Vmagnitude = Entry(master)
        self.Vdirection = Entry(master)

        self.pileheight.grid(row=2,column=1)
        self.volfract.grid(row=3,column=1)
        self.xpilecenter.grid(row=4,column=1)
	self.ypilecenter.grid(row=4,column=2)
        self.majradius.grid(row=5,column=1)
        self.minradius.grid(row=5,column=2)
        self.orientation.grid(row=6,column=1)
        self.Vmagnitude.grid(row=7,column=1)
        self.Vdirection.grid(row=8,column=1)

        Button(master, text="Done", command=self.done).grid(row=9,column=0)
        Button(master, text="Quit", command=master.quit).grid(row=9,column=1)
        Button(master, text="Calculate\n Volume", command=self.showVolume).grid(row=9,column=2)

        # create Map button if started from within grass
        mapbutt = Button(master, text="Map", command=self.showMap )
        mapbutt.grid(row=1,column=2)
        self.newmon = IntVar()
        newmonbutt = Checkbutton(master, text="New Monitor", command=self.manual_toggle)
        newmonbutt.grid(row=2, column=2)
        newmonbutt.toggle()
        self.newmon.set( 1 ) # workaround, bug in TkInter?

        if not os.environ.has_key('GISBASE'):
            mapbutt.config(state=DISABLED)
            newmonbutt.config(state=DISABLED)
        
        #passed in variables
        self.master = master
        self.top = master.winfo_toplevel()
        self.filename = filename
        self.max_height = max_height
        self.write_flag = 0 # used to make sure that the information is written out only once per call
        self.directory = directory
        self.input_flag = 0 #used to make sure that the values are entered
        self.topomap = topomap
        self.pile_number = pile_number

        # initialize entries from environment variables
        if pile_number == 0: # but only for the 1st pile
            if os.environ.has_key('GMFG_PILEH'):
                self.pileheight.insert(INSERT, os.environ['GMFG_PILEH'])
            if os.environ.has_key('GMFG_PILEX'):
                self.xpilecenter.insert(INSERT, os.environ['GMFG_PILEX'])
            if os.environ.has_key('GMFG_PILEY'):
                self.ypilecenter.insert(INSERT, os.environ['GMFG_PILEY'])
            if os.environ.has_key('GMFG_PILEMAJE'):
                self.majradius.insert(INSERT, os.environ['GMFG_PILEMAJE'])
            if os.environ.has_key('GMFG_PILEMINE'):
                self.minradius.insert(INSERT, os.environ['GMFG_PILEMINE'])
            if os.environ.has_key('GMFG_PILEOR'):
                self.orientation.inser(INSERT, os.environ['GMFG_PILEOR'])

    def manual_toggle(self):
        # this function is a workaround a strange problem in TkInter
        # the same piece of code works fine in the main window
        if self.newmon.get() == 1:
            self.newmon.set( 0 )
        else :
            self.newmon.set( 1 )
        #print "variable is ", self.newmon.get()
                  
    def showMap(self):
        # calls the _pilehelper with the map name and output file name
        # the takes the output file and parses it into coordinates
        if self.newmon.get() == 1:
            helper = os.path.dirname(os.path.abspath(sys.argv[0])) + '/titan_pilehelper'
            os.system('sh ' + helper + ' ' + self.topomap + ' ' + self.directory + '/pile.coord')
        else:
            os.system('d.where -1 > ' + self.directory + '/pile.coord')
        m = re.compile('^\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)').match( open(self.directory + '/pile.coord').readline() )
        if m:
           self.xpilecenter.delete(0, END) 
           self.xpilecenter.insert(END, m.group(1)) 
           self.ypilecenter.delete(0, END) 
           self.ypilecenter.insert(END, m.group(2)) 
        else:
            tkMessageBox.showerror(
                "Pile location",
                "Unable to parse d.where output,\ne-mail sorokine@buffalo.edu"
                )

    def showVolume(self):
        if self.pileheight.get() == '':
            self.ph = 0.
        else:
            self.ph=float(self.pileheight.get())

        if self.majradius.get() == '':
            self.majr = 0.
        else:
            self.majr=float(self.majradius.get())
            
        if self.minradius.get() == '':
            self.minr = 0
        else:
            self.minr=float(self.minradius.get())
            
        self.volume=math.pi*self.ph*self.majr*self.minr/2.
        tkMessageBox.showinfo(self,message=self.volume)

    def done(self):
	# check values
	if self.pileheight.get() == '':
	    pileheight = 0.0
	else:
	    pileheight = float(self.pileheight.get())
	    if pileheight <= float(0.):
		pileheight = 0.

        if self.volfract.get() =='':
            volfract = 1.0
        else:
            volfract = float(self.volfract.get())
            if volfract < 0:
                volfract = 0.
            elif volfract > 1.:
                volfract = 1. 

	if self.xpilecenter.get() == '':
	    xpilecenter = 1.0
	else:
	    xpilecenter = float(self.xpilecenter.get())

	if self.ypilecenter.get() == '':
	    ypilecenter = 1.0
	else:
	    ypilecenter = float(self.ypilecenter.get())

	if self.majradius.get() == '':
	    majradius = 1.0
	else:
	    majradius = float(self.majradius.get())
	    if majradius <= 0:
		majradius = 1

	if self.minradius.get() == '':
	    minradius = 1.0
	else:
	    minradius = float(self.minradius.get())
	    if minradius <= 0:
		minradius = 1

        if self.orientation.get() == '':
            orientation = 0.0
        else:
            orientation = float(self.orientation.get())

        if self.Vmagnitude.get() == '':
            Vmagnitude = 0.0
            Vdirection = 0.0
        else:
            Vmagnitude=float(self.Vmagnitude.get())
            if self.Vdirection.get() == '':
                Vdirection = 0.0
            else:
                Vdirection=float(self.Vdirection.get())
                            
        #print 'Vmagnitude= ' + str(Vmagnitude) + '\nVdirection= ' + str(Vdirection) + '\n'

        if self.write_flag == 0:
            self.write_flag = 1
            file = open(self.filename, "a+", 0)
            file.write( str(pileheight) + '\n' +  str(volfract) + '\n' + str(xpilecenter) + '\n' + str(ypilecenter) + '\n' + str(majradius) + '\n' +str(minradius) + '\n' + str(orientation) + '\n' + str(Vmagnitude) + '\n' + str(Vdirection) + '\n')
            file.close
            if pileheight > self.max_height:
                self.max_height = pileheight

            self.input_flag = 1
            #write to the python_input.data file
            file = open(self.directory+'python_input.data', "a+", 0)
            python_data = """Maximum Initial Thickness, P (m): """ + str(pileheight) + """Initial solid-volume fraction: """ + str(volfract) + """ Center of Initial Volume, xc, yc (UTM E, UTM N): """ + str(xpilecenter) + """ """+str(ypilecenter)+"""
Major and Minor Extent, majorR, minorR (m, m): """ + str(majradius)+ """ """ +str(minradius)+"""
Angle from X axis to major axis (degrees): """ +str(orientation)+"""
Initial speed [m/s]: """ + str(Vmagnitude) + """
Initial direction ([degrees] from X axis): """ + str(Vdirection)+"""
"""
            file.write(python_data)
            file.close

        # now feed the data back to the environment
        if self.pile_number == 1: # but only for the 1st pile
            os.environ['GMFG_PILEH']  = self.pileheight.get()
            os.environ['GMFG_VOLFRAC']  = self.volfract.get()
            os.environ['GMFG_PILEX']  = self.xpilecenter.get()
            os.environ['GMFG_PILEY']  = self.ypilecenter.get()
            os.environ['GMFG_PILEMAJE'] = self.majradius.get()
            os.environ['GMFG_PILEMINE'] = self.minradius.get()
            os.environ['GMFG_PILEOR'] = self.orientation.get()
        
class QuestionTemplate:

    def __init__(self,master):
        self.master=master

        Label(master, text="GIS Information Main Directory:").grid(row=0,column=0,sticky=W)
        Label(master, text="GIS Sub-Directory:").grid(row=1,column=0,sticky=W)
        Label(master, text="GIS Map Set:").grid(row=2,column=0,sticky=W)
        Label(master, text="GIS Map:").grid(row=3,column=0,sticky=W)
        Label(master, text="GIS Vector:").grid(row=4,column=0,sticky=W)
        Label(master, text="Use GIS Material Map?").grid(row=5,column=0,sticky=W)
        Label(master, text="Simulation Directory Location:").grid(row=6,column=0,sticky=W)
        Label(master, text="Number of Processors:").grid(row=7,column=0,sticky=W)
        Label(master, text="Number of Computational Cells Across\n Smallest Pile/Flux-Source Diameter:").grid(row=8,column=0,sticky=W)
        Label(master, text="Number of Piles:").grid(row=9,column=0,sticky=W)
        Label(master, text="Number of Flux Sources:").grid(row=10,column=0,sticky=W)
        Label(master, text="Number of Discharge Planes: ").grid(row=11,column=0,sticky=W)
        Label(master, text="Scale Simulation? ").grid(row=12,column=0,sticky=W)
        Label(master, text="If Scaled, Length Scale [m]:").grid(row=13,column=0,sticky=W)
        Label(master, text="Maximum Number of Time Steps:").grid(row=14,column=0,sticky=W)
        Label(master, text="Maximum Time [sec]:").grid(row=15,column=0,sticky=W)
        Label(master, text="Time [sec] between Results Output:").grid(row=16,column=0,sticky=W)
        Label(master, text="Time [sec] between Saves: ").grid(row=17,column=0,sticky=W)
        Label(master, text="Adapt the Grid?").grid(row=18,column=0,sticky=W)	
        Label(master, text="Visualization Output:").grid(row=19,column=0,sticky=W)
        Label(master, text="First/Second Order Method:").grid(row=20,column=0,sticky=W)
        Label(master, text="Minimum x and y location (UTM E, UTM N):").grid(row=21,column=0,sticky=W)
        Label(master, text="Maximum x and y location (UTM E, UTM N):").grid(row=22,column=0,sticky=W)
        Label(master, text="Height used to define flow outline (>0) [m]: ").grid(row=23,column=0,sticky=W)
        Label(master, text="Test if flow reaches height [m] ...:").grid(row=24,column=0,sticky=W)
        Label(master, text="... at test point (x and y location):").grid(row=25,column=0,sticky=W)

        Label(master,  text="Email Address:").grid(row=26,column=0,sticky=W)

        self.topomain = Entry(master)
        self.toposub = Entry(master)
        self.topomapset = Entry(master)
        self.topomap = Entry(master)
        self.vector =  Entry(master)
        self.matmap = IntVar()
        self.matmap_fn = Checkbutton(master,variable=self.matmap,text="Yes")
        self.directory = Entry(master)
        self.numprocs = Entry(master)
        self.numcellsacrosspile = Entry(master)
        self.numpiles = Entry(master)
        self.numsrcs = Entry(master)
        self.numdischarge = Entry(master)
        self.scale = IntVar()
        self.scale_fn = Checkbutton(master,variable=self.scale,text="Yes")
        self.lengthscale = Entry(master)
        self.steps = Entry(master)
        self.maxtime = Entry(master)
        self.timeoutput = Entry(master)
        self.timesave = Entry(master)
        self.adapt = IntVar()
        self.adapt_fn = Checkbutton(master,variable=self.adapt,text="Yes")
        self.vizoutput = Menubutton(master, text="Choose Formats",relief=RAISED)
        self.order = IntVar()
        self.order_fn = Checkbutton(master,variable=self.order,text="Second")
        self.min_location_x = Entry(master)
        self.min_location_y = Entry(master)
        self.max_location_x = Entry(master)
        self.max_location_y = Entry(master)
        self.edge_height = Entry(master)
        self.test_height = Entry(master)
        self.test_location_x = Entry(master)
        self.test_location_y = Entry(master)
        self.emailaddress = Entry(master) 
	
        self.topomain.grid(row=0,column=1)
        self.toposub.grid(row=1,column=1)
        self.topomapset.grid(row=2,column=1)
        self.topomap.grid(row=3,column=1)
        self.vector.grid(row=4,column=1)
        self.matmap_fn.grid(row=5,column=1)
        self.directory.grid(row=6,column=1)
        self.numprocs.grid(row=7,column=1)
        self.numcellsacrosspile.grid(row=8,column=1)
        self.numpiles.grid(row=9,column=1)
        self.numsrcs.grid(row=10,column=1)
        self.numdischarge.grid(row=11,column=1)
        self.scale_fn.grid(row=12,column=1)
        self.lengthscale.grid(row=13,column=1)
        self.steps.grid(row=14,column=1)
        self.maxtime.grid(row=15,column=1)
        self.timeoutput.grid(row=16,column=1)
        self.timesave.grid(row=17,column=1)
        self.adapt_fn.grid(row=18,column=1)
        self.vizoutput.grid(row=19,column=1)
        self.order_fn.grid(row=20,column=1)
        self.min_location_x.grid(row=21,column=1)
        self.min_location_y.grid(row=21,column=2)
        self.max_location_x.grid(row=22,column=1)
        self.max_location_y.grid(row=22,column=2)
        self.edge_height.grid(row=23,column=1)
        self.test_height.grid(row=24,column=1)
        self.test_location_x.grid(row=25,column=1)
        self.test_location_y.grid(row=25,column=2)
        #self.vizoutput.grid()
        self.emailaddress.grid(row=26,column=1)

        Button(master, text="Run", command=self.run).grid(row=27,column=0)
        Button(master, text="Quit", command=master.quit).grid(row=27,column=1)
        Button(master, bitmap="question", command=self.getHelp).grid(row=27,column=2)

        self.vizoutput.menu = Menu(self.vizoutput, tearoff=0)
        self.vizoutput["menu"] = self.vizoutput.menu
        self.tecplotVar = IntVar()
        self.mshplotVar = IntVar()
        self.xdmfVar = IntVar()
        self.padyVar = IntVar()
        self.grasssitesVar = IntVar()
        self.quickviewVar = IntVar()
        
        self.vizoutput.menu.add_checkbutton(label="tecplotxxxx.plt",variable=self.tecplotVar)
        self.vizoutput.menu.add_checkbutton(label="mshplotxxxx.plt",variable=self.mshplotVar)
        self.vizoutput.menu.add_checkbutton(label="XDMF/Paraview",variable=self.xdmfVar)
        self.vizoutput.menu.add_checkbutton(label="grass_sites",variable=self.grasssitesVar)
        

        # initialize entries from environment variables
        if os.environ.has_key('GISDBASE'):
            self.topomain.insert(INSERT, os.environ['GISDBASE'])
        if os.environ.has_key('LOCATION_NAME'):
            self.toposub.insert(INSERT, os.environ['LOCATION_NAME'])
        if os.environ.has_key('MAPSET'):
            self.topomapset.insert(INSERT, os.environ['MAPSET'])
        if os.environ.has_key('GMFG_MAP'):
            self.topomap.insert(INSERT, os.environ['GMFG_MAP'])
        if os.environ.has_key('GMFG_MAT_MAP'):
            if os.environ['GMFG_MAT_MAP'] == "1":
                self.matmap_fn.toggle()
        if os.environ.has_key('GMFG_DIR'):
            self.directory.insert(INSERT, os.environ['GMFG_DIR'])
        if os.environ.has_key('GMFG_MP'):
            self.numprocs.insert(INSERT, os.environ['GMFG_MP'])
        if os.environ.has_key('GMFG_MESH'):
            self.numcellsacrosspile.insert(INSERT, os.environ['GMFG_MESH'])
#        if os.environ.has_key('GMFG_IANG'):
#            self.intfrict.insert(INSERT, os.environ['GMFG_IANG'])
#        if os.environ.has_key('GMFG_BANG'):
#            self.bedfrict.insert(INSERT, os.environ['GMFG_BANG'])
        if os.environ.has_key('GMFG_PILES'):
            self.numpiles.insert(INSERT, os.environ['GMFG_PILES'])
        if os.environ.has_key('GMFG_SCALE'):
            if os.environ['GMFG_SCALE'] == "1":
                self.scale_fn.toggle()
        if os.environ.has_key('GMFG_LENGTH'):
            self.lengthscale.insert(INSERT, os.environ['GMFG_LENGTH'])
        if os.environ.has_key('GMFG_MAXTS'):
            self.steps.insert(INSERT, os.environ['GMFG_MAXTS'])
        if os.environ.has_key('GMFG_MAXTIME'):
            self.maxtime.insert(INSERT, os.environ['GMFG_MAXTIME'])
        if os.environ.has_key('GMFG_OUTTIME'):
            self.timeoutput.insert(INSERT, os.environ['GMFG_OUTTIME'])
        if os.environ.has_key('GMFG_ADAPT'):
            if os.environ['GMFG_ADAPT'] == "1":
                self.adapt_fn.toggle()
        if os.environ.has_key('GMFG_OUTFMT'):
            fmt = os.environ['GMFG_OUTFMT']
            if re.compile("tecplot").search(fmt):
                self.tecplotVar.set(1)
            if re.compile("mshplot").search(fmt):
                self.mshplotVar.set(1)
            if re.compile("gmf").search(fmt):
                self.padyVar.set(1)
            if re.compile("XDMF/Paraview").search(fmt):
                self.xdmfVar.set(1)
            if re.compile("grasssites").search(fmt):
                self.grasssitesVar.set(1)
            if re.compile("WebViz").search(fmt):
                self.quickviewVar.set(1)
        #if os.environ.has_key(''):
        #    self.vizoutput.insert(INSERT, os.environ['GMFG_'])
        if os.environ.has_key('GMFG_EMAIL'):
            self.emailaddress.insert(INSERT, os.environ['GMFG_EMAIL'])


        # create widget that are activated  if started from within grass
        regionbutt = Button(master, text="Region", command=self.selRegion)
        regionbutt.grid(row=0,column=2)
        self.newmon = IntVar()
        newmonbutt = Checkbutton(master, text="New Monitor", variable = self.newmon)
        newmonbutt.grid(row=1, column=2)
        newmonbutt.toggle()

        mapsetmenubutt = Menubutton(master, text="Mapsets",relief=RAISED)
        mapsetmenubutt.grid(row=2, column=2)
        mapmenubutt = Menubutton(master, text="  Maps  ",relief=RAISED)
        mapmenubutt.grid(row=3, column=2)

        testpointbutt = Button(master, text="test point", command=self.selTestPoint)
        testpointbutt.grid(row=24, column=2)
        
        if os.environ.has_key('GISBASE'):

            # disable these fields because they cannot be changed
            self.topomain.config(state=DISABLED, relief='flat')
            self.topomain.xview(END)
            self.toposub.config(state=DISABLED, relief='flat')
            
            # setup menu of mapsets
            self.mapsetname = StringVar()
            self.mapsetmenu = Menu(mapsetmenubutt,
                                   tearoff=0,
                                   postcommand=self.updateMapsetMenu
                                   )
            self.topomapset.bind("<Button-3>", self.mapsetPopup)
            mapsetmenubutt["menu"] = self.mapsetmenu

            # create a list of maps
            self.mapname = StringVar()
            self.mapmenu = Menu(mapmenubutt,
                                tearoff=0,
                                postcommand=self.updateMapMenu
                                )
            self.topomap.bind("<Button-3>", self.mapPopup)    
            mapmenubutt["menu"] = self.mapmenu

        else:
            regionbutt.config(state=DISABLED)
            mapsetmenubutt.config(state=DISABLED)
            mapmenubutt.config(state=DISABLED)
            newmonbutt.config(state=DISABLED)
            testpointbutt.config(state=DISABLED)
            
    # END OF __init__      
 
    def updateMapsetMenu(self):
        self.mapsetmenu.delete(0,END)
        for mapsetname in os.popen("g.mapsets -l | tr ' ' '\012' | sed -e '/^$/d' ").readlines():
            self.mapsetmenu.add_radiobutton(
                label=string.rstrip(mapsetname),
                command=self.mapsetSelected,
                variable=self.mapsetname,
                indicatoron=0
                )
        self.mapsetmenu.add_radiobutton(label='-<CANCEL>-',indicatoron=0)

    def mapsetSelected(self):
        self.topomapset.delete(0, END)
        self.topomapset.insert(0, self.mapsetname.get())

    def mapsetPopup(self, event):
        self.mapsetmenu.post(event.x_root, event.y_root)


    def updateMapMenu(self):
        self.mapmenu.delete(0,END)
        mapsetname = ''
        if self.topomapset.get() != '':
            mapsetname = 'mapset='+self.topomapset.get()
        for mapname in os.popen("g.list rast " + mapsetname + " | sed -e '/-/d; /:/d; /^\s*$/d; s/\t/ /g; s/ \+/ /g' | tr ' ' '\012'").readlines():
            self.mapmenu.add_radiobutton(
                label=string.rstrip(mapname),
                command=self.mapSelected,
                variable=self.mapname,
                indicatoron=0
                )
        self.mapmenu.add_radiobutton(label='-<CANCEL>-',indicatoron=0)

    def mapSelected(self):
        self.topomap.delete(0, END)
        self.topomap.insert(0, self.mapname.get())

    def mapPopup(self, event):
        self.mapmenu.post(event.x_root, event.y_root)

    def selRegion(self):
        if self.topomap.get() != '':
            helper = os.path.dirname(os.path.abspath(sys.argv[0])) + '/titan_regionhelper'
            os.system('sh ' + helper + ' ' +
                      self.topomap.get() + ' ' +
                      str(self.newmon.get()) + ' ' +
                      '/tmp/region.coord' )
            
            # parse coordinate file
            m = re.compile('^(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)').match( open('/tmp/region.coord').readline() )
            if m:
                self.min_location_x.delete(0, END)
                self.min_location_y.delete(0, END)
                self.max_location_x.delete(0, END)
                self.max_location_y.delete(0, END)

                self.min_location_x.insert(END, m.group(3)) 
                self.min_location_y.insert(END, m.group(2)) 
                self.max_location_x.insert(END, m.group(4)) 
                self.max_location_y.insert(END, m.group(1)) 
            else:
                tkMessageBox.showerror(
                    "Region Selection",
                    "Unable to parse region,\ne-mail sorokine@buffalo.edu"
                    )
        else:
            tkMessageBox.showwarning(
                "GRASS Region",
                "Specify GIS map name first!"
                )
            
    def selTestPoint(self):
        # just copied Alex's pile centroid point selector and changed
        # where the values were stored
        # calls the _pilehelper with the map name and output file name
        # the takes the output file and parses it into coordinates
        if self.newmon.get() == 1:
            helper = os.path.dirname(os.path.abspath(sys.argv[0])) + '/titan_pilehelper'
            os.system('sh ' + helper + ' ' + self.topomap.get() + ' ' + '/tmp/pile.coord')
            #os.system('sh ' + helper + ' ' + self.topomap.get() + ' ' + self.directory + '/pile.coord')
        else:
            os.system('d.where -1 > ' + '/tmp/pile.coord')
        m = re.compile('^\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)').match( open('/tmp/pile.coord').readline() )
        if m:
           self.test_location_x.delete(0, END) 
           self.test_location_x.insert(END, m.group(1)) 
           self.test_location_y.delete(0, END) 
           self.test_location_y.insert(END, m.group(2)) 
        else:
            tkMessageBox.showerror(
                "Test point location",
                "Unable to parse d.where output,\ne-mail sorokine@buffalo.edu"
                )

    def getHelp(self):
        a=open('../README').read()
        tk = Tkinter.Tk()
        frame = Tkinter.Frame(tk, relief=FLAT, borderwidth=0)
        frame.pack(fill=BOTH,expand=1)
        label = Tkinter.Label(frame, text="Instructions")
        label.pack(fill=X, expand=1)
        text=Tkinter.Text(frame,relief=RIDGE,width=75)
        text.insert(END,a)
        text.pack(side=LEFT)
        button=Tkinter.Button(frame,text="OK",command=tk.destroy)
        button.pack(side=RIGHT,anchor=S)
        scroll = Tkinter.Scrollbar(frame,command=text.yview,width=20)
        scroll.pack(side=RIGHT,ipady=80,anchor=W)
        text["yscrollcommand"] = scroll.set
        tk.mainloop()


    def run(self):
	os.system('echo "Please email abani@eng.buffalo.edu with any problems"')
        #get system information so it is known which system the script is running on
        machine = socket.gethostbyaddr(socket.gethostname())
        print 'Trying to run a job on ' + machine[0]
	# if there is no topo file, quit
	if self.topomain.get() == '' or self.toposub.get() == '' or self.topomapset.get() == '' or self.topomap.get() == '':
            print 'Missing GIS information.  No job will be run.'
            return
	#elif os.access(self.topofile.get(), os.F_OK) == 0:
	#    print 'cannot find '+self.topofile.get()
	#    return
		
	str_numprocs = self.numprocs.get()
	if str_numprocs == '':
	    numprocs = 1
	else:
	    numprocs = int(str_numprocs)
	    
	if numprocs <= 0:
	    print 'numprocs must be greater than 0, it is ' + str(numprocs)
	    numprocs = 1

        if numprocs != 1 and numprocs != 2 and numprocs != 4 and numprocs != 8 and numprocs != 16 and numprocs != 32 and numprocs != 64 and numprocs != 128 and numprocs != 256 and numprocs != 512:
            print 'wrong amount of processors!'
            return
            
        if numprocs > 8 and machine[0] == 'crosby.ccr.buffalo.edu':
            print 'you need to have access to the mp queues for that many processors on crosby, email ccr-help@ccr.buffalo.edu for help'
            return

        if numprocs > 32 and machine[0] == 'stills.ccr.buffalo.edu':
            print 'must use 32 or less processors on stills'
            return
        
	#create proper directory
	if os.access(self.directory.get(), os.F_OK) == 1:
	    print self.directory.get() + ' already exists, will not overwrite it'
	    return

	dir1 = self.directory.get()
        if dir1 == '':
            print 'Must specify a simulation directory.'
            return
        
	dirlen = len(dir1)
	if dir1[dirlen-1] != '/':
	    directory = dir1 + '/'
	else:
	    directory = dir1 
            
        instruction_string= "Completed setup ot Titan input data.\ncd to simulation directory, \"" + directory + "\"\nand run Titan (on linux PC's type \"./titan\")"

	os.system('mkdir ' + directory)

        #get (then write) list of material names and properties

        #if you don't enter a value for the current material it
        #defaults to the value for the previous material, so must
        #provide some previous value for the first material
        #since there is value checking in QuestionTemplate3, any
        #values can be used here
        previntfrict = 0.
        prevbedfrict = 0.

        fout=open("frict.data","w",0)
        counter=0

        #if they didn't want to use a GIS material map, get the material
        #properties once up front
        if self.matmap.get() == 0:
            nummat=1
            fout.write(str(nummat)+'\n')
            matname = 'all materials'
            root3=Tk()
            app=QuestionTemplate3(root3,nummat,counter,matname,previntfrict,prevbedfrict)
            #root3.pack()
            root3.mainloop()
            root3.withdraw()
            fout.write(matname+'\n')
            fout.write(str(app.intfrict)+' '+str(app.bedfrict)+'\n')
            del app

        else:  #if they did want to use a GIS material map...
            
            #run C++ program to generate a list of material names
            os.system('./titan_materialnames '+self.topomain.get() +' '+self.toposub.get()+' '+self.topomapset.get()+' '+self.topomap.get());
            
            fin=open("materialnames.dat", "r", 0);

            #read in the number of materials (the integer on the first 
            #line) from the file "materialnames.dat"
            str_nummat=''
            tempchar=fin.read(1);
            while tempchar != '\n':
                str_nummat = str_nummat + tempchar
                tempchar=fin.read(1);
            nummat=int(str_nummat)
            #print 'nummat=' + str(nummat)
            del str_nummat

            fout.write(str(nummat)+'\n')
            
            while counter < nummat:
                matname=''
                tempchar=fin.read(1);
                while tempchar != '\n':
                    matname = matname + tempchar
                    tempchar=fin.read(1);
                    
                #print 'matname= [' + matname + ']'
                root3=Tk()
                app=QuestionTemplate3(root3,nummat,counter,matname,previntfrict,prevbedfrict)
                root3.mainloop()
                root3.withdraw()
                fout.write(matname+'\n')
                fout.write(str(app.intfrict)+' '+str(app.bedfrict)+'\n')
                previntfrict=app.intfrict
                prevbedfrict=app.bedfrict
                del app
                counter = counter + 1


            fin.close
            del tempchar
            del counter
            #remove this file so it can't accidentally be reused next time
            os.system('rm -f materialnames.dat');


        del matname
        fout.close

       
        #check values
        srctype = 0
        if self.numpiles.get() == '':
            numpiles = 0
	    print 'Number of piles not set.  Setting it to ' + str(numpiles)
        else:
            numpiles = int(self.numpiles.get())
            if numpiles < 0:
                print 'Number of piles cannot be a negative number.  Setting it to 0'
                numpiles = 0

        if self.numsrcs.get() == '':
            numsrcs = 0
	    print 'Number of Flux not set.  Setting it to ' + str(numsrcs)
        else:
            numsrcs = int(self.numsrcs.get())
            if numsrcs < 0:
                print 'Number of Flux Sources cannot be a negative number.  Setting it to 0'
                numsrcs = 0

        if numpiles > 0:
            srctype = srctype +1

        if numsrcs  > 0:
            srctype = srctype +2

        coord_flag = 0
        min_location_x = 0
        max_location_x = 0
        min_location_y = 0
        max_location_y = 0
        if self.min_location_x.get() != '':
            coord_flag = 1
            min_location_x = float(self.min_location_x.get())

        if self.min_location_y.get() != '':
            coord_flag = 1+coord_flag
            min_location_y = float(self.min_location_y.get())

        if self.max_location_x.get() != '':
            coord_flag = 1+coord_flag
            max_location_x = float(self.max_location_x.get())

        if self.max_location_y.get() != '':
            coord_flag = 1+coord_flag
            max_location_y = float(self.max_location_y.get())

        if min_location_x > max_location_x:
            temp = max_location_x
            max_location_x = min_location_x
            min_location_x = temp

        if min_location_y > max_location_y:
            temp = max_location_y
            max_location_y = min_location_y
            min_location_y = temp

        if coord_flag != 0 and coord_flag != 4:
            print 'Must either fill in none or all minimum and maximum coordinate values'
            return

	numcellsacrosspile = 20
	if self.numcellsacrosspile.get() != '':
	    numcellsacrosspile = int(self.numcellsacrosspile.get())
	    if numcellsacrosspile <= 0:
		numcellsacrosspile = 20

	if self.steps.get() == '':
	    steps = 100
	else:
	    steps = self.steps.get()
	    if steps < 1:
		steps = 100

	if self.timeoutput.get() == '':
	    timeoutput = 10
	else:
	    timeoutput = self.timeoutput.get()
	    if timeoutput < 0:
		timeoutput = 10

	if self.maxtime.get() == '':
	    maxtime = 1.5
	else:
	    maxtime = float(self.maxtime.get())
	    if maxtime <= 0:
		maxtime = 1.5

        if self.timesave.get() == '':
            timesave = -1
        else:
            timesave = float(self.timesave.get())
            if timesave <= 0:
                timesave = -1
                if timesave > maxtime:
                    timesave = -1

        if self.edge_height.get() == '':
            edge_height = -1
        else:
            edge_height = float(self.edge_height.get())
            if edge_height <= 0:
                print 'you entered an edge height <= 0!\nusing default value hardcoded in titan instead\n'
                edge_height = -1
                            
        if self.test_height.get() == '':
            test_height = -1
        else:
            test_height =float(self.test_height.get())
            if test_height <= 0:
                 print 'you entered an point test height <= 0!\nusing default value hardcoded in titan instead\n'
                 test_height = -1

        if self.test_location_x.get() == '' or self.test_location_y.get() == '':
            test_height = -2
            test_location_x = 'none'
            test_location_y = 'none'
            
        else:
            test_location_x = str(float(self.test_location_x.get()))
            test_location_y = str(float(self.test_location_y.get()))

            
        # make a copy of the input file for reference
        if self.matmap.get() == 1:
            str_matmap = 'yes'
        else:
            str_matmap = 'no'

        if self.scale.get() == 1:
            str_scale = 'yes'
        else:
            str_scale = 'no'

        if self.adapt.get() == 1:
            str_adapt = 'yes'
        else:
            str_adapt = 'no'

        run_file = directory+'python_input.data'
        run_info = """NOTE:  This file outputs the data used in the runs, NOT the data input into the Python script
GIS Information Main Directory:  """ + self.topomain.get() + """
GIS Sub-Directory:  """ + self.toposub.get() + """
GIS Map Set:  """ + self.topomapset.get() + """
GIS Map:  """ + self.topomap.get() + """
GIS Vector:  """ + self.vector.get() + """
Use GIS Material Map?  """ + str_matmap + """
(there are """ + str(nummat) + """ different materials)
Simulation Directory Location:  """ + self.directory.get() + """
Number of Processors:  """ + str(numprocs) + """
Number of Computational Cells Across Smallest Pile Diameter:  """ + str(numcellsacrosspile) + """
Number of Piles:  """ + str(numpiles) + """
Number of Flux Sources: """ +str(numsrcs) + """
Scale Simulation? """ + str_scale + """
If Scaled, Length Scale [m]:  """ + self.lengthscale.get() + """
Maximum Number of Time Steps:  """ + str(steps) + """
Maximum Time [sec]:  """ + str(maxtime) +"""
Time [sec] between Results Output:  """ + str(timeoutput) + """
Time [sec] between Saves: """ + str(timesave) + """
Adapt the Grid?  """ + str_adapt + """
First/Second Order Method:  """ + str(self.order.get()+1) + """
Minimum x and y location (UTM E, UTM N):  """ + self.min_location_x.get() + """  """ +self.min_location_y.get() + """
Maximum x and y location (UTM E, UTM N):  """ + self.max_location_x.get() + """  """ +self.max_location_y.get() + """
Flow edge defined to have pile height (>0) [m]: """ + str(edge_height) + """
Test if flow reaches height [m] ...: """ + str(test_height) + """
...at test point (x and y location (UTM E, UTM N)):  """ + test_location_x + """ """ + test_location_y + """
"""
	f=open(run_file, "w", 0)
        f.write(run_info)
	f.close

        counter = 0
        max_height = 0.0
        #pile geometry stuff, friction coefficients, etc.
        # get all of the pile information
        f_p = directory+'simulation.data'
        f_p2=open(f_p, "w", 0)
        f_p2.write(str(srctype) + '\n')
        if int(numpiles) > 0:
            f_p2.write(str(numpiles) + '\n')
        if int(numsrcs) > 0:
            f_p2.write(str(numsrcs) + '\n')
        f_p2.close
        if numpiles > 0:
            while counter < numpiles:
                root2=Tk()
                app=QuestionTemplate2(root2,counter,f_p,directory,max_height,self.topomap.get())
                root2.mainloop()
                #print 'counter is ' + str(counter) + ' and numpiles is ' + str(numpiles)
                counter = counter + app.input_flag
                max_height = app.max_height
                root2.withdraw()
                del app

        counter = 0
        heightscale = max_height
        if numsrcs > 0:
            while counter < numsrcs:
                root5=Tk()
                app=QuestionTemplate5(root5,counter,f_p,directory, heightscale)
                root5.mainloop()
                counter = counter + app.input_flag
                heightscale = app.heightscale;
                root5.withdraw()
                del app
            
	output2 = str(numcellsacrosspile) + '\n' + str(steps) + '\n' + str(maxtime) + '\n' + str(timeoutput) + '\n' + str(timesave) + '\n' +  str(self.adapt.get())
        f_p2=open(f_p, "a+", 0)        
	f_p2.write(output2)

        #scaling stuff
        scale = self.scale.get()
        lengthscale = 1.
        if scale == 1:
	    if self.lengthscale.get() == '':
		lengthscale = 1.
	    else:
		lengthscale = float(self.lengthscale.get())
		if lengthscale <= 0.:
		    lengthscale = 1.
            if max_height <= 0.:
                max_height = 1.
            if heightscale <= 0.:
                heightscale = 1.
                
	    output1 = str(lengthscale) + "\n" + str(heightscale) + "\n9.80" 
	    
	else:
	    output1 = "1 \n1 \n1"
        # feed the parameters back to environment
        os.environ['GMFG_MAP'] = self.topomap.get()
        os.environ['GMFG_PROP_MAP'] = str(self.matmap.get())
        os.environ['GMFG_DIR'] = self.directory.get()
        os.environ['GMFG_MP'] = self.numprocs.get()
        os.environ['GMFG_MESH'] = self.numcellsacrosspile.get()
        os.environ['GMFG_PILES'] = self.numpiles.get()
        os.environ['GMFG_LENGTH'] = self.lengthscale.get()
        os.environ['GMFG_MAXTS'] = self.steps.get()
        os.environ['GMFG_MAXTIME'] = self.maxtime.get()
        os.environ['GMFG_OUTTIME'] = self.timeoutput.get()
        os.environ['GMFG_EMAIL'] = self.emailaddress.get()
        os.environ['GMFG_SCALE'] = str(self.scale.get())
        os.environ['GMFG_ADAPT'] = str(self.adapt.get())
        fmt = []
        if self.tecplotVar.get() == 1:
            fmt.append('tecplot')
        if self.mshplotVar.get() == 1:
            fmt.append('meshlot')
        if self.padyVar.get() == 1:
            fmt.append("gmf")
        if self.xdmfVar.get() == 1:
            fmt.append("XDMF/Paraview")
        if self.grasssitesVar.get() == 1:
            fmt.append("grasssites")
        if self.quickviewVar.get() == 1:
           fmt.append("quickview")
        os.environ['GMFG_OUTFMT'] = string.join(fmt, ',')

	f2 = directory+ 'scale.data'
	f=open(f2, "w", 0)
	f.write(output1)
	f.close

        #put in the stuff for viz the idea is to use prime numbers and the remainder function to determine which formats to output in
        viz_num = 0
        #print 'tecplot var is ' + str(self.tecplotVar.get())
        if self.tecplotVar.get() == 1:
            viz_num = viz_num +1 #first bit flag

        #print 'meshplot var is ' + str(self.mshplotVar.get())
        if self.mshplotVar.get() == 1:
            viz_num = viz_num +2 #second bit flag

        #print 'gmfg viz ' + str(self.padyVar.get())
        if self.padyVar.get() == 1:
            viz_num = viz_num +4 #third bit flag

        #print 'xdmf var is ' + str(self.xdmfVar.get())
        if self.xdmfVar.get() == 1:
            viz_num = viz_num +8 #fourth bit flag

        #print 'grasssites var is ' + str(self.grasssitesVar.get())
        if self.grasssitesVar.get()== 1:
            viz_num = viz_num + 16 #fifth bit flag

        #print 'quickview var is ' + str(self.quickviewVar.get())
        if self.quickviewVar.get() == 1:
            viz_num = viz_num +32 #sixth bit flag
            os.system('mkdir webviz')
            os.system('mkdir webviz/topo')
            os.system('mkdir webviz/data')
            os.system('mkdir webviz/textures')
            os.system('mv webviz ' + directory)

        f_p2.write('\n' + str(viz_num) + '\n' + str(self.order.get()+1))
        #GIS stuff
        f_p2.write('\n' + self.topomain.get() + '\n' + self.toposub.get() + '\n' + self.topomapset.get() + '\n' + self.topomap.get() +'\n' + str(self.matmap.get()))
        f_p2.write('\n' + str(edge_height) + '\n' + str(test_height) + '\n' + test_location_x + ' ' + test_location_y)
	f_p2.close

        #----------------------------------------
        #-----Number of Discharge Planes---------
        #----------------------------------------
        #check values
        if self.numdischarge.get() == '':
            numdischarge = 0
            print 'Number of Discharge Planes not set.  Setting it to ' + str(numdischarge)
        else:
            numdischarge = int(self.numdischarge.get())
            if numdischarge < 0:
                print 'Number of Discharge Planes cannot be a negative number.  Setting it to 0'
                numdischarge = 0

        #loop to allow user input of coordinates for variable # of discharge planes
        counter=0
        fout = open(directory+'simulation.data',"a+",0)
        fout.write('\n'+str(numdischarge)+'\n')
        while counter < numdischarge:
            root4=Tk()
            app=QuestionTemplate4(root4,counter,numdischarge)
            root4.mainloop()
            root4.withdraw()
            if app.xa != '' and app.xb != '' and app.ya != '' and app.yb != '':
                fout.write(str(app.xa)+' '+str(app.ya)+' '+str(app.xb)+' '+str(app.yb)+'\n')
                del app
                counter = counter + 1
            else:
                print 'Missing data.  Must re-enter.'
        fout.close
        del counter
        #----------------------------------------------
        #----------------------------------------------


        print 'max height is ' + str(max_height)
        print 'heightscale is ' + str(heightscale)
	#run the makefile (clean out all c and c++ object files first)
        if os.access('titan',os.X_OK) == 0:
            if machine[0] == 'crosby.ccr.buffalo.edu' or machine[0] == 'nash.ccr.buffalo.edu' or machine[0] == 'joplin.ccr.buffalo.edu' or  machine[0] == 'popo.eng.buffalo.edu' or machine[0] == 'colima.eng.buffalo.edu' or machine[0] == 'rainier.eng.buffalo.edu' or machine[0] == 'elchichon.eng.buffalo.edu' or machine[0] == 'casita.eng.buffalo.edu' or machine[0] == 'sthelens.eng.buffalo.edu':
                print 'Need to install titan with ub-compile.sh script'
            else:
                print 'Need to install titan'



	# run preproc.x to create the fem grid, if it is not already there
	#if os.access('PRE/preproc.x',os.X_OK)==0:
        #    os.system('cd PRE;gmake')
        # locate titan_preprocess: 
        # first check the local-diretory,
        if ( os.path.isfile('titan_preprocess') ):
            preproc = './titan_preprocess'
        else:
            preproc = 'titan_preprocess'
	    
        if coord_flag == 0:
            os.system(preproc+' '+str(numprocs)+' '+str(numcellsacrosspile)+' '+self.topomain.get() +' '+self.toposub.get()+' '+self.topomapset.get()+' '+self.topomap.get())
        else:
            print 'window is  '+str(min_location_x)+' '+str(min_location_y)+' '+str(max_location_x)+' '+str(max_location_y)
#            print '|'+str(self.topomain.get())+'|'+str(self.toposub.get())+'|'+str(self.topomapset.get())+'|'+str(self.topomap.get())+'|'
            os.system(preproc+' '+str(numprocs)+' '+str(numcellsacrosspile)+' '+self.topomain.get() +' '+self.toposub.get()+' '+self.topomapset.get()+' '+self.topomap.get()+' '+str(min_location_x)+' '+str(min_location_y)+' '+str(max_location_x)+' '+str(max_location_y))

        if self.vector.get() != "":
            os.system('./VecDataPreproc ' + self.topomain.get() +' '+self.toposub.get()+' '+self.topomapset.get()+' '+self.topomap.get()+' '+self.vector.get())
            os.system('mv VectorDataOutput.data ' + directory)
            
	#if os.access('PRE/preprocess',os.X_OK) != 1:
	#    os.system('cd PRE;gmake -f Makefile_C')

	#os.system('./titan_preprocess '+str(numprocs))
	#os.system('cp -f ' + self.topofile.get() +' '+ directory +'/topo.data')
	os.system('mv funky*inp '+directory)
	os.system('cp titan ' + directory)
        if os.access('Viewer',os.X_OK)==1:
            os.system('cp Viewer ' + directory)
        os.system('mv frict.data '+directory)
        os.system('cp ../src/stochastic/*.pl ../src/stochastic/pbslhs  lhsbed lhsvol  lhstitanstats  ../src/stochastic/lhs.readme '+directory)

        #	os.system('mpirun -np ' + str(numprocs) + directory+'titan')
############################################################
############################################################
###  machine specific stuff here...
############################################################
############################################################
        if machine[0] == 'crosby.ccr.buffalo.edu' or machine[0] == 'nash.ccr.buffalo.edu' or machine[0] == 'joplin.ccr.buffalo.edu':
            iam = self.emailaddress.get()
            print 'sending email to : '+iam
        elif machine[0] == 'stills.ccr.buffalo.edu':
            print 'no email notification on stills'
        elif machine[0] == 'popo.eng.buffalo.edu' or machine[0] == 'colima.eng.buffalo.edu' or machine[0] == 'rainier.eng.buffalo.edu' or machine[0] == 'elchichon.eng.buffalo.edu' or machine[0] == 'casita.eng.buffalo.edu' or machine[0] == 'sthelens.eng.buffalo.edu':
            print 'no email notification on Linux machines'
            
        if machine[0] == 'crosby.ccr.buffalo.edu':            
            batchscript = """#!/bin/csh -f 
#PBS -m e
#PBS -l ncpus="""+ str(numprocs)+ """
#PBS -l walltime=50:00:00
#PBS -l mem=128mb
#PBS -M """+iam +"""
#PBS -o volcano.out.crosby
#PBS -j oe
#PBS -N volc."""+str(numprocs)+ """.crosby
#PBS
limit coredumpsize 0
#
# stage input files and executable to scratch - $PBS_O_WORKDIR is
#  the directory from which the job was submitted.
#
cp $PBS_O_WORKDIR/*.inp $PBSTMPDIR
cp $PBS_O_WORKDIR/titan $PBSTMPDIR
cp $PBS_O_WORKDIR/*.data $PBSTMPDIR
# goto scratch dir and run job
cd $PBSTMPDIR
date
mpirun -cpr -np """+str(numprocs)+"""  ./titan
date
#
# stage output files back
## cp *.out $PBS_O_WORKDIR/
gzip *out *gz *h5
# remove scratch directory
#rm -r $PBSTMPDIR
"""
            f6 = directory+'pbs_script'
            f=open(f6, "w", 0)
            f.write(batchscript)
            f.close
            #os.system('cd '+directory+';qsub pbs_script')


#################################
# nash.ccr.buffalo.edu
#################################
        if machine[0] == 'nash.ccr.buffalo.edu':
            if numprocs ==1:
                b2 = "#PBS -l nodes=1:ppn=1"
            else:
                b2 = "#PBS -l nodes="+str(numprocs/2)+":ppn=2"
                
            batchscript = """#!/bin/csh
#PBS -m e
""" + b2 + """
#PBS -l walltime=30:00:00
#PBS -M """+iam +"""
#PBS -o volcano.out.nash
#PBS -j oe
#PBS -N volcano."""+str(numprocs)+ """.nash
#PBS
limit coredumpsize 0
#
cd $PBS_O_WORKDIR
set NP = `cat $PBS_NODEFILE | wc -l`
set NN = `cat $PBS_NODEFILE | wc -l`
echo "NN = "$NN
set plist = `cat $PBS_NODEFILE | sed "s/\.ccr\.buffalo\.edu//"`
set uplist = `cat $PBS_NODEFILE | uniq | sed "s/\.ccr\.buffalo\.edu//"`
echo "unique nodes list = "$uplist
echo "full process list = "$plist

@ NN = -1
foreach p ($plist)
 @ NN ++
  set NUM = `echo $NN | awk '{printf "%.4d", $1}'`
  rsh $p "cp $PBS_O_WORKDIR/funky$NUM.inp $PBSTMPDIR/."
end

foreach p ($uplist)
  rsh $p "cp $PBS_O_WORKDIR/titan $PBSTMPDIR/"
  rsh $p "cp $PBS_O_WORKDIR/*.data $PBSTMPDIR/"
  echo "Contents of "$PBSTMPDIR" on node "$p
  rsh $p "ls -l $PBSTMPDIR/"
end

source /util/tag-gm-gnu.csh
cd $PBSTMPDIR
date
/util/mpich-gm/gnu/ch_gm/bin/mpirun.ch_gm -v -machinefile $PBS_NODEFILE  -np """+str(numprocs)+"""  ./titan
date
#
foreach p ($uplist)
  rsh $p "cd $PBSTMPDIR; gzip *out *plt *h5"
  rsh $p "cp $PBSTMPDIR/*gz $PBS_O_WORKDIR/"
end
"""
            f6 = directory+'pbs_script'
            f=open(f6, "w", 0)
            f.write(batchscript)
            f.close
            #os.system('cd '+directory+';qsub pbs_script')
#################################
# joplin.ccr.buffalo.edu -- added by abani 02/27
################################
        elif machine[0] == 'joplin.ccr.buffalo.edu':
            if numprocs ==1:
                b2 = "#PBS -l nodes=1:ppn=1"
            else:
                b2 = "#PBS -l nodes="+str(numprocs/2)+":ppn=2"
                
            batchscript = """#!/bin/csh -f 
""" + b2 + """
#PBS -l walltime=00:59:00
#PBS -M """ + iam + """
#PBS -m e
#PBS -j oe
#PBS -o volcano."""+str(numprocs)+""".joplin
#
# Set executable name & mpiexec path
#  EXE = Executable name
#  INITDIR = dir holding input & executable
#  SAVEDIR = dir to place output files
#
source /util/tag-gm-gnu.csh
set EXE = titan
set INITDIR = $PBS_O_WORKDIR
set SAVEDIR = $PBS_O_WORKDIR
#
# Get list of processors from PBS

set plist = `cat $PBS_NODEFILE | sed "s/\.ccr\.buffalo\.edu//"`
set uplist = `cat $PBS_NODEFILE | uniq | sed "s/\.ccr\.buffalo\.edu//"`
echo "unique nodes list = "$uplist
echo "full process list = "$plist
#
# Stage out indexed input files (and common ones too) to each
#  processor
# (also build an mpiexec config file with the ranks to correspond
#  to the staged inputs, or so we hope - I still don't see that
#  as guaranteed!)
#
cd $PBS_O_WORKDIR
set conf = $PBSTMPDIR/conf.$$
if (-e $conf) rm -f $conf
@ NN = -1
foreach p ($plist)
 @ NN ++
  set NUM = `echo $NN | awk '{printf "%.4d", $1}'`
  rsh $p "cp $INITDIR/funky$NUM.inp $PBSTMPDIR"
  echo "$p : $EXE" >> $conf
end
#
#  files that only need to appear once on a node, no matter
#  how many procs
#
foreach p ($uplist)
  rsh $p "cp $INITDIR/$EXE $PBSTMPDIR/"
  rsh $p "cp $INITDIR/*.data $PBSTMPDIR/"
  echo "Contents of "$PBSTMPDIR" on node "$p
  rsh $p "ls -l $PBSTMPDIR/"
end
echo "Contents of mpiexec config file:"
cat $conf
#
# Run code from $PBSTMPDIR
#
cd $PBSTMPDIR
/usr/bin/time mpiexec -config $conf
#
# Stage back output files
#
foreach p ($uplist)
  rsh $p "gzip $PBSTMPDIR/*out; gzip $PBSTMPDIR/*plt; gzip $PBSTMPDIR/*h5; cp $PBSTMPDIR/*gz $SAVEDIR/"
end
rm -f $conf
"""
            f6 = directory+'pbs_script'
            f=open(f6, "w", 0)
            f.write(batchscript)
            f.close
            #os.system('cd '+directory+';qsub pbs_script')

#################################
# stills.ccr.buffalo.edu - not working as of 5/7/03 but left in for future use/installation of proper libraries
#################################
        elif machine[0] == 'stills.ccr.buffalo.edu':
            #figure out the distribution of processors
            if numprocs <=4:
                nodes = 1
                procs_per_node = numprocs
                requirements = '(Memory>256)'
            else:
                nodes = numprocs/4
                procs_per_node = 4
                requirements = '(Pool==2)'

            batchscript = """###########################################################################
## A sample Loadl file. 
## Only run jobs in $LOADL_SCRATCH or /gpfs.
###########################################################################
## classes are:  Short, Medium, Long, V.Long
#!/bin/csh -xf 
# @ class =  V.Long
# @ environment = COPY_ALL
# @ error       = volcano"""+str(numprocs)+""".err..$(jobid)
# @ output      = volcano"""+str(numprocs)+""".out.$(jobid)
# @ network.MPI = css0,not_shared,US
# @ job_type = parallel
# @ requirements = """+requirements+"""
# @ tasks_per_node="""+str(procs_per_node)+"""
# @ node ="""+str(nodes)+"""
# @ notification = never
# @ queue
#
##########################################################################
##
##  Change to suit your needs

set EXECDIR=`pwd`
set EXECNAME=titan
set WORKDIR=`pwd`

if (-e $EXECDIR/EXECNAME) then
   echo ' Executible not found'
   echo ' Directory checked was $EXECDIR'
   echo ' EXE filename checked was $EXECNAME'
   exit1
else
   echo ' Found the EXE file $EXECDIR/$EXECNAME'
endif
set SCRATCH = $LOADL_SCRATCHDIR

############################################################################
##
## No need to change this. 
## Automatically sets up a proper scratch directory on all nodes. 
##

set SCRATCH = $LOADL_SCRATCHDIR
cd $SCRATCH

set LL_hosts = ($LOADL_PROCESSOR_LIST)
echo 'Number of processes is' $#LL_hosts
echo 'The scheduled processors are on nodes ' $LL_hosts
set noclobber


#############################################################################
## Change to suit
## Add to the copy list as necessary - never remove $EXECNAME copy
##      
mcp $EXECDIR/$EXECNAME $SCRATCH/.
mcp $EXECDIR/scale.data $SCRATCH/.
mcp $EXECDIR/topo.data $SCRATCH/.
mcp $EXECDIR/simulation.data $SCRATCH/.
"""
            f6 = directory+'loadleveler_script'
            f=open(f6, "w", 0)
            f.write(batchscript)
            for i in range(numprocs):
                if i < 10:
                    f.write("mcp $EXECDIR/funky0"+str(i)+".inp $SCRATCH/.\n")
                else:
                    f.write("mcp $EXECDIR/funky"+str(i)+".inp $SCRATCH/.\n")

            f.write("""############################################################################
## Change to suit
## Invoke POE
time poe << EOF
    $SCRATCH/$EXECNAME 
    quit
EOF

############################################################################
## Change to suit
## Copy output files back to $WORKDIR
## need to do some tricky stuff here because the output files are
## in scratch space on all of the nodes and only 1 node runs this script

foreach node ($LL_hosts)
	rcp "$node":"$SCRATCH/viz_output*.out" $WORKDIR/.
end
gzip viz*out

# ----------------------------------------------------------------
# Last second information before giving up the nodes/disks
#
date
""")
            f.close
            os.system('cd '+directory+';llsubmit loadleveler_script')
            
        elif machine[0] == 'popo.eng.buffalo.edu' or machine[0] == 'colima.eng.buffalo.edu' or machine[0] == 'elchichon.eng.buffalo.edu' or machine[0] == 'sthelens.eng.buffalo.edu':
            machine2 = os.uname()
            f8 = directory+'machines'
            f=open(f8, "w",0)
            f.write(machine2[1])
            f.close
            i = 0
            while i == 0:
                i =  os.access(f_p, os.F_OK)
                                        
            os.system('cd '+directory+';/usr/local/mpich-1.2.4/ch_p4/bin/mpirun -machinefile machines -np '+str(numprocs)+' titan')

        elif machine[0] == 'rainier.eng.buffalo.edu':
            root3=Tk()
            app=FinalInstructions(root3,instruction_string)
            root3.mainloop()
            root3.withdraw()
            exit
#            self.master.quit
#            master.quit

            machine2 = os.uname()
            f8 = directory+'machines'
            f=open(f8, "w",0)
            f.write(machine2[1])
            f.close
            i = 0
            while i == 0:
                i =  os.access(f_p, os.F_OK)

#            os.system('cd '+directory+';/usr/local/mpich-1.2.6/bin/mpirun -machinefile machines -np '+str(numprocs)+' titan')
                
        elif machine[0] == 'casita.eng.buffalo.edu':

            machine2 = os.uname()
            f8 = directory+'machines'
            f=open(f8, "w",0)
            f.write(machine2[1])
            f.close
            i = 0
            while i == 0:
                i =  os.access(f_p, os.F_OK)
                                        
            os.system('cd '+directory+';/usr/local/mpich-1.2.5/bin/mpirun -machinefile machines -np '+str(numprocs)+' titan')

        else:
            print 'TITAN2D is ready to run in  ' +directory
            print """The files that the executable 'titan' needs to run are funky*.inp, scale.data, and simulation.data"""
            print 'The output files will be:'
            print 'output_summary.readme,'
            if self.tecplotVar.get() == 1:
                print 'tecplot*plt'

            if self.mshplotVar.get() == 1:
                print 'mshplot*plt'

            if self.padyVar.get() == 1:
                print 'tri_output*.out, viz_filenames.out, and viz_output*out'

            if self.xdmfVar.get() == 1:
                print 'xdmf*.xdf'
                
            if self.grasssitesVar.get() == 1:
                print 'grass_sites*.*'

            if self.quickviewVar.get() == 1:
                print 'l*r*.asc'
                

root=Tk()
app=QuestionTemplate(root)
root.mainloop()
            




