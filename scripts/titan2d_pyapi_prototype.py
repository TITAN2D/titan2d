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
#* $Id: titan_gui.py 235 2012-03-29 07:04:35Z dkumar $ 
#*

import sys,os,math,string,re,socket

TITAN2d_HOME='/home/mikola/titan_wsp/titan2d_bld/iccopt'


class TitanFluxSource:
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

class TitanDischargePlane:

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


class TitanPile:

    def __init__(self,pileheight=0.0,pilecenter=None,radii=None,
                 orientation=0.0,
                 Vmagnitude=0.0,
                 Vdirection=0.0):
        #def __init__(self,master,pile_number,filename,directory,max_height,topomap):
        
        
        #Information for Pile Number 
        #Thickness of Initial Volume, h(x,y)
        #P*(1-((x-xc)/xr)^2 - ((y-yc)/yr)^2)
        #Maximum Initial Thickness, P (m)
        self.pileheight = float(pileheight)
        #Center of Initial Volume, xc, yc (UTM E, UTM N)
        if pilecenter!=None:
            self.xpilecenter = float(pilecenter[0])
            self.ypilecenter = float(pilecenter[1])
        else:
            self.xpilecenter = 1.0
            self.ypilecenter = 1.0
        #Major and Minor Extent, majorR, minorR (m, m)
        if radii!=None:
            self.majradius = float(radii[0])
            self.minradius = float(radii[1])
        else:
            self.majradius = 1.0
            self.minradius = 1.0
        #Orientation (angle [degrees] from X axis to major axis)
        self.orientation = float(orientation)
        #Initial speed [m/s]
        self.Vmagnitude = float(Vmagnitude)
        #Initial direction ([degrees] from X axis)
        self.Vdirection = float(Vdirection)
        
        # create Map button if started from within grass
        #mapbutt = Button(master, text="Map", command=self.showMap )
        
        #passed in variables
        if 0:
            self.master = master
            self.top = master.winfo_toplevel()
            self.filename = filename
            self.max_height = max_height
            self.write_flag = 0 # used to make sure that the information is written out only once per call
            self.directory = directory
            self.input_flag = 0 #used to make sure that the values are entered
            self.topomap = topomap
            self.pile_number = pile_number
        
        self.validateValues()
        
    def validateValues(self):
        if self.pileheight < 0.0:
            raise ValueError('TitanPile::pileheight should be non negative')
    def showVolume(self): 
        self.volume=math.pi*self.pileheight*self.majradius*self.minradius/2.0
        print "Pile volume:",self.volume

    def done(self,filename):
        pileheight = self.pileheight
        xpilecenter = self.xpilecenter
        ypilecenter = self.ypilecenter
        
        majradius = self.majradius
        minradius = self.minradius
        
        orientation = self.orientation
        Vmagnitude=self.Vmagnitude
        Vdirection=self.Vdirection
        
        fout = open(filename, "a+", 0)
        fout.write( str(pileheight) + '\n' + str(xpilecenter) + '\n' + str(ypilecenter) + '\n' + str(majradius) + '\n' +str(minradius) + '\n' + str(orientation) + '\n' + str(Vmagnitude) + '\n' + str(Vdirection) + '\n')
        fout.close
        
class TitanSimulation:
    possible_vizoutputs={
        'tecplotxxxx.tec':1, # first bit flag
        'mshplotxxxx.tec':2, # second bit flag
        'XDMF/Paraview':4, # third bit flag
        'grass_sites':8, # fourth bit flag
        'tecplot':1, # first bit flag
        'meshlot':2, # second bit flag
        'grasssites':8 # fourth bit flag
    }
    possible_orders=('PlaceHolder','First','Second')
    def __init__(self,
                 
                 numcellsacrosspile,
                 
                 numsrcs,
                 numdischarge,
                 steps,
                 maxtime,
                 timeoutput,
                 timesave,
                 numprocs=1,
                 lengthscale=None,
                 adapt=False,
                 vizoutput="tecplotxxxx.tec",
                 order='First',
                 min_location_x=None,
                 min_location_y=None,
                 max_location_x=None,
                 max_location_y=None,
                 edge_height=None,
                 test_height=None,
                 test_location_x=None,
                 test_location_y=None
                 ):
        
        
        
        
        #init values
        self.gisformat = None
        self.topomain = None
        self.toposub = None
        self.topomapset = None
        self.topomap = None
        self.vector=None
        self.matmap=None
        
        #Number of Processors
        self.numprocs = numprocs
        if numprocs <= 0:
            raise ValueError('numprocs must be greater than 0, it is ' + str(numprocs))
        if numprocs not in (1,2,4,8,12,128,256,512):
            raise ValueError('wrong amount of processors!')
        
        #Number of Computational Cells Across Smallest Pile/Flux-Source Diameter
        self.numcellsacrosspile = numcellsacrosspile
        #Number of Piles
        self.numpiles = None #int(numpiles)
        #Number of Flux Sources
        self.numsrcs = int(numsrcs)
        #Number of Discharge Planes
        self.numdischarge = int(numdischarge)
        
        #Scale Simulation?
        if lengthscale!=None:
            self.scale = True
            #If Scaled, Length Scale [m]
            self.lengthscale = float(lengthscale)
            if self.lengthscale<=0.0:
                raise ValueError("TitanSimulation::lengthscale should be positive")
        else:
            self.scale = False
            self.lengthscale = 1.0
        
        #Maximum Number of Time Steps
        self.steps = steps
        #Maximum Time [sec]
        self.maxtime = maxtime
        #Time [sec] between Results Output
        self.timeoutput = timeoutput
        #Time [sec] between Saves
        self.timesave = timesave
        #Adapt the Grid?
        self.adapt = adapt
        #Visualization Output
        if vizoutput in TitanSimulation.possible_vizoutputs:
            self.vizoutput = TitanSimulation.possible_vizoutputs[vizoutput]
        else:
            raise ValueError("Unknown vizoutput "+str(vizoutput)+". Possible formats: "+str(possible_vizoutputs.keys()))
        
        #First/Second Order Method
        if order in TitanSimulation.possible_orders:
            self.order = TitanSimulation.possible_orders.index(order)
        else:
            raise ValueError("Unknown order "+str(order)+". Possible formats: "+str(possible_orders[1:]))
        
        
        #Minimum x and y location (UTM E, UTM N)
        self.min_location_x = min_location_x
        self.min_location_y = min_location_y
        #Maximum x and y location (UTM E, UTM N)
        self.max_location_x = max_location_x
        self.max_location_y = max_location_y
        #Height used to define flow outline (>0) [m]
        self.edge_height = edge_height
        #Test if flow reaches height [m] ...
        self.test_height = test_height
        #... at test point (x and y location)
        self.test_location_x = test_location_x
        self.test_location_y = test_location_y
        
        #other inits
        self.piles=[]
        

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

    def setTopo(self,gisformat='GIS_GRASS',
                 topomain=None,
                 toposub=None,
                 topomapset=None,
                 topomap=None,
                 vector=None):
        possible_gisformats=('PlaceHolder','GIS_GRASS', 'GDAL')
        if gisformat in possible_gisformats:
            self.gisformat = possible_gisformats.index(gisformat)
        else:
            raise ValueError("Unknown gisformat "+str(gisformat)+". Possible formats: "+str(gisformats[1:]))
        #GIS Information Main Directory
        self.topomain = topomain
        #GIS Sub-Directory
        self.toposub = toposub
        #GIS Map Set
        self.topomapset = topomapset
        #GIS Map
        self.topomap = topomap
        #GIS Vector
        self.vector = vector
        
        #here should be validator
        # if there is no topo file, quit
        if self.gisformat == 1:
            errmsg='Missing GIS information.  No job will be run.'
            if self.topomain == '' or self.toposub == '' or self.topomapset == '' or self.topomap == '':
                raise ValueError(errmsg)
            if self.topomain == None or self.toposub == None or self.topomapset == None or self.topomap == None:
                raise ValueError(errmsg)
            if (not isinstance(self.topomain, basestring)) or (not isinstance(self.toposub, basestring)) or \
                    (not isinstance(self.topomapset, basestring)) or (not isinstance(self.topomap, basestring)):
                raise ValueError(errmsg)
            
            p=""
            for dirname in (self.topomain,self.toposub,self.topomapset):
                p=os.path.join(p,dirname)
                print p
                if not os.path.isdir(p):
                    raise ValueError(errmsg+". "+p+" does not exist!")
            #self.topomap?
        elif  self.gisformat == 2:
            if (self.topomap == '') or (self.topomap == None) or (not isinstance(self.topomap, basestring)) or \
                    (not os.path.isdir(self.topomap)):
                raise ValueError(errmsg)
            
    def setMatMap(self,
            useGIS_MatMap=False,
            matMap=None):
        #Use GIS Material Map?
        self.matmap = useGIS_MatMap
        if matMap==None:
            matMap=[{'intfrict':30.0,'bedfrict':15.0}]
        
        #get (then write) list of material names and properties

        #if you don't enter a value for the current material it
        #defaults to the value for the previous material, so must
        #provide some previous value for the first material
        #since there is value checking in QuestionTemplate3, any
        #values can be used here
        previntfrict = 0.
        prevbedfrict = 0.
        

        fout=open("frict.data","w",0)

        #if they didn't want to use a GIS material map, get the material
        #properties once up front
        if self.matmap == False:
            nummat=1
            fout.write(str(nummat)+'\n')
            matname = 'all materials'
            
            # check values and store them
            intfrict=float(matMap[0]['intfrict'])
            bedfrict=float(matMap[0]['bedfrict'])
            
            if intfrict <= 0.0 :
                raise ValueError('intfrict can not be negative')
            if bedfrict <= 0.0 :
                raise ValueError('bedfrict can not be negative')
            
            fout.write(matname+'\n')
            fout.write(str(intfrict)+' '+str(bedfrict)+'\n')
            
        else:  #if they did want to use a GIS material map...
            raise Exception("Not implemented as there were no suitable example")
        fout.close
    def setPiles(self,piles):
        #Information for Pile Number
        #Thickness of Initial Volume, h(x,y)
        #P*(1-((x-xc)/xr)^2 - ((y-yc)/yr)^2)
        #Maximum Initial Thickness, P (m)
        self.pileheight = None
        #Center of Initial Volume, xc, yc (UTM E, UTM N)
        self.xpilecenter = None
        self.ypilecenter = None
        #Major and Minor Extent, majorR, minorR (m, m) 
        self.majradius = None
        self.minradius = None
        #Orientation (angle [degrees] from X axis to major axis)
        self.orientation = None
        #Initial speed [m/s]
        self.Vmagnitude = None
        #Initial direction ([degrees] from X axis)
        self.Vdirection = None
    
    def addPile(self,**kwargs):
        pile=TitanPile(**kwargs)
        if pile!=None:
            self.piles.append(pile)
        
    
    def run(self):
        #get system information so it is known which system the script is running on
        machine = socket.gethostbyaddr(socket.gethostname())
        print 'Trying to run a job on ' + machine[0]
        
        #check values
        srctype = 0
        
        self.numpiles=len(self.piles)
        numpiles = self.numpiles
        if numpiles < 0:
            raise ValueError('Number of piles cannot be a negative number')
        
        numsrcs = self.numsrcs
        if numsrcs < 0:
            raise ValueError('Number of Flux Sources cannot be a negative number')
        

        if numpiles > 0:
            srctype = srctype +1

        if numsrcs  > 0:
            srctype = srctype +2
            
        coord_flag = 0
        min_location_x = 0
        max_location_x = 0
        min_location_y = 0
        max_location_y = 0
        if self.min_location_x != None:
            coord_flag = 1
            min_location_x = float(self.min_location_x)

        if self.min_location_y != None:
            coord_flag = 1+coord_flag
            min_location_y = float(self.min_location_y)

        if self.max_location_x != None:
            coord_flag = 1+coord_flag
            max_location_x = float(self.max_location_x)

        if self.max_location_y != None:
            coord_flag = 1+coord_flag
            max_location_y = float(self.max_location_y)

        if min_location_x > max_location_x:
            temp = max_location_x
            max_location_x = min_location_x
            min_location_x = temp

        if min_location_y > max_location_y:
            temp = max_location_y
            max_location_y = min_location_y
            min_location_y = temp

        if coord_flag != 0 and coord_flag != 4:
            raise ValueError('Must either fill in none or all minimum and maximum coordinate values')
        
        numcellsacrosspile = 20
        if self.numcellsacrosspile != None:
            numcellsacrosspile = int(self.numcellsacrosspile)
            if numcellsacrosspile <= 0:
                numcellsacrosspile = 20

        if self.steps == None:
            steps = 100
        else:
            steps = self.steps
            if steps < 1:
                steps = 100

        if self.timeoutput == None:
            timeoutput = 10
        else:
            timeoutput = self.timeoutput
            if timeoutput < 0:
                timeoutput = 10

        if self.maxtime == None:
            maxtime = 1.5
        else:
            maxtime = float(self.maxtime)
            if maxtime <= 0:
                maxtime = 1.5

        if self.timesave == None:
            timesave = -1
        else:
            timesave = float(self.timesave)
            if timesave <= 0:
                timesave = -1
                if timesave > maxtime:
                    timesave = -1

        if self.edge_height == None:
            edge_height = -1
        else:
            edge_height = float(self.edge_height)
            if edge_height <= 0:
                print 'you entered an edge height <= 0!\nusing default value hardcoded in titan instead\n'
                edge_height = -1
                            
        if self.test_height == None:
            test_height = -1
        else:
            test_height =float(self.test_height)
            if test_height <= 0:
                 print 'you entered an point test height <= 0!\nusing default value hardcoded in titan instead\n'
                 test_height = -1

        if self.test_location_x == None or self.test_location_y == None:
            test_height = -2
            test_location_x = 'none'
            test_location_y = 'none'
            
        else:
            test_location_x = str(float(self.test_location_x))
            test_location_y = str(float(self.test_location_y))

            
        
        
        #pile geometry stuff, friction coefficients, etc.
        # get all of the pile information
        max_height = 0.0
        f_p = 'simulation.data'
        f_p2=open(f_p, "w", 0)
        f_p2.write(str(srctype) + '\n')
        if int(numpiles) > 0:
            f_p2.write(str(numpiles) + '\n')
        if int(numsrcs) > 0:
            f_p2.write(str(numsrcs) + '\n')
        f_p2.close
        
        max_height=0.0
        for iPile in range(len(self.piles)):
            self.piles[iPile].showVolume()
            self.piles[iPile].done(f_p)
            if self.piles[iPile].pileheight > max_height:
                max_height = self.piles[iPile].pileheight
        

        counter = 0
        heightscale = max_height
        if numsrcs > 0:
            raise NotImplementedError("numsrcs > 0 not implemented yet!")
            #while counter < numsrcs:
                #app=TitanFluxSource/QuestionTemplate5(root5,counter,f_p,directory, heightscale)
                #heightscale = app.heightscale;
            
        output2 = str(numcellsacrosspile) + '\n' +\
            str(steps) + '\n' + str(maxtime) + '\n' +\
            str(timeoutput) + '\n' + str(timesave) + '\n' +\
            str(int(self.adapt))
        f_p2=open(f_p, "a+", 0)        
        f_p2.write(output2)
        
        #scaling stuff
        scale = self.scale
        lengthscale = 1.
        if scale:
            lengthscale = self.lengthscale
            if max_height <= 0.:max_height = 1.
            if heightscale <= 0.:heightscale = 1.
                
            output1 = str(lengthscale) + "\n" + str(heightscale) + "\n9.80" 
            
        else:
            output1 = "1 \n1 \n1"
        
        f2 = 'scale.data'
        f=open(f2, "w", 0)
        f.write(output1)
        f.close

        #put in the stuff for viz the idea is to use prime numbers and the remainder function to determine which formats to output in
        viz_num = self.vizoutput

        f_p2.write('\n' + str(viz_num) + '\n' + str(self.order))
        #GIS stuff
        f_p2.write('\n' + str(self.gisformat))
        if self.gisformat == 1:
            f_p2.write('\n' + self.topomain + '\n' + self.toposub + '\n' + self.topomapset + '\n' + self.topomap +'\n' + str(int(self.matmap)))
        else:
            f_p2.write('\n' + self.topomap)

        f_p2.write('\n' + str(edge_height) + '\n' + str(test_height) + '\n' + test_location_x + ' ' + test_location_y)
        f_p2.close

        #----------------------------------------
        #-----Number of Discharge Planes---------
        #----------------------------------------
        #check values
        if self.numdischarge == 0:
            numdischarge = 0
            print 'Number of Discharge Planes not set.  Setting it to ' + str(numdischarge)
        else:
            numdischarge = self.numdischarge

        #loop to allow user input of coordinates for variable # of discharge planes
        fout = open('simulation.data',"a+",0)
        fout.write('\n'+str(numdischarge)+'\n')
        
        if numdischarge>0:
            raise NotImplementedError("numsrcs > 0 not implemented yet!")
            #app=TitanDischargePlane/QuestionTemplate4(root4,counter,numdischarge)
#             if app.xa != '' and app.xb != '' and app.ya != '' and app.yb != '':
#                 fout.write(str(app.xa)+' '+str(app.ya)+' '+str(app.xb)+' '+str(app.yb)+'\n')
#                 del app
#                 counter = counter + 1
#             else:
#                 print 'Missing data.  Must re-enter.'
        fout.close
        #----------------------------------------------
        #----------------------------------------------


        print 'max height is ' + str(max_height)
        print 'heightscale is ' + str(heightscale)
        
        
        # run preproc.x to create the fem grid, if it is not already there
        #if os.access('PRE/preproc.x',os.X_OK)==0:
        #    os.system('cd PRE;gmake')
        # locate titan_preprocess: 
        # first check the local-diretory,
        
        numprocs=self.numprocs
        
        if ( os.path.isfile('titan_preprocess') ):
            preproc = './titan_preprocess'
        else:
            preproc = '/home/mikola/titan_wsp/titan2d_bld/iccopt/bin/titan_preprocess'
            
        if coord_flag == 0:
            if self.gisformat == 1:
                command=preproc+' '+str(numprocs)+' '+str(numcellsacrosspile)+' '+str(self.gisformat)+' '+self.topomain +' '+self.toposub+' '+self.topomapset+' '+self.topomap
                print "executing:",command
                os.system(command)
            else:
                command=preproc+' '+str(numprocs)+' '+str(numcellsacrosspile)+' '+str(self.gisformat)+' '+self.topomap
                print "executing:",command
                os.system(command)
        else:
            print 'window is  '+str(min_location_x)+' '+str(min_location_y)+' '+str(max_location_x)+' '+str(max_location_y)
            if self.gisformat == 1:
                command=preproc+' '+str(numprocs)+' '+str(numcellsacrosspile)+' '+str(self.gisformat)+' '+self.topomain +' '+self.toposub+' '+self.topomapset+' '+self.topomap+' '+str(min_location_x)+' '+str(min_location_y)+' '+str(max_location_x)+' '+str(max_location_y)
                print "executing:",command
                os.system(command)
            else:
                command=preproc+' '+str(numprocs)+' '+str(numcellsacrosspile)+' '+str(self.gisformat)+' '+self.topomap+' '+str(min_location_x)+' '+str(min_location_y)+' '+str(max_location_x)+' '+str(max_location_y)
                print "executing:",command
                os.system(command)

        if self.vector != None:
            raise NotImplementedError("GIS Vector is Not Implemented yet in py api!")
            #os.system('./VecDataPreproc ' + self.topomain +' '+self.toposub+' '+self.topomapset+' '+self.topomap+' '+self.vector)
            #os.system('mv VectorDataOutput.data ' + directory)
        
    
if __name__ == "__main__":
    if len(sys.argv)!=2:
        print "usage titan2d.py <simulation.py>"
    execfile(sys.argv[1])

