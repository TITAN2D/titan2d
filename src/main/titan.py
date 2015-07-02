
#*******************************************************************
#* Copyright (C) 2003-2015 University at Buffalo
#*
#* This software can be redistributed free of charge.  See COPYING
#* file in the top distribution directory for more details.
#*
#* This software is distributed in the hope that it will be useful,
#* but WITHOUT ANY WARRANTY; without even the implied warranty of
#* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#*
#* Author: NA Simakov 
#* Description: titan python api
#*
#******************************************************************* 
#*

import sys,os,math,string,re,socket

from cxxtitan import *
#cxxTitanSimulation,TitanPreproc,cxxTitanPile,cxxTitanFluxSource,MaterialMap,vectordatpreproc


class TitanFluxSource(cxxTitanFluxSource):
    def __init__(self,influx,start_time,end_time,center=None,radii=None,
                 orientation=0.0,
                 Vmagnitude=0.0,
                 Vdirection=0.0):
        super(TitanFluxSource, self).__init__()
        #Extrusion flux rate [m/s]
        self.influx = float(influx)
        
        #Active Time [s], start, end
        self.start_time = float(start_time)
        self.end_time   = float(end_time)
        #Center of the source, xc, yc (UTM E, UTM N):
        if center!=None:
            self.xcenter = float(center[0])
            self.ycenter = float(center[1])
        else:
            self.xcenter = 1.0
            self.ycenter = 1.0
        #Major and Minor Extent, majorR, minorR (m, m)
        if radii!=None:
            self.majradius = float(radii[0])
            self.minradius = float(radii[1])
        else:
            self.majradius = 1.0
            self.minradius = 1.0
        #Orientation (angle [degrees] from X axis to major axis):
        self.orientation = float(orientation)
        #Initial speed [m/s]:
        self.Vmagnitude = float(Vmagnitude)
        #Initial direction ([degrees] from X axis):
        self.Vdirection = float(Vdirection)
        
        self.validateValues()
    def validateValues(self):
        if self.influx < 0.0:
            raise ValueError('TitanFluxSource::influx should be non negative')
        if self.start_time < 0.0:
            raise ValueError('TitanFluxSource::start_time should be non negative')
        if self.end_time < 0.0:
            raise ValueError('TitanFluxSource::start_time should be non negative')
    def done(self,filename):
        file = open(filename, "a+", 0)
        file.write( str(self.influx) + '\n' + \
                    str(self.start_time) + '\n' + str(self.end_time) + '\n' + \
                    str(self.xcenter) + '\n' + str(self.ycenter) + '\n' + \
                    str(self.majradius) + '\n' +str(self.minradius) + '\n' + \
                    str(self.orientation) + '\n' + \
                    str(self.Vmagnitude) + '\n' + \
                    str(self.Vdirection) + '\n')
        file.close
        
       

        self.input_flag = 1
        #write to the python_input.data file
        file = open('python_input.data', "a+", 0)
        python_data = """ Mean Flux  (kg/(m^2-s)): """+ str(self.influx) +"""
Active duration of source (start time, end time) (s) " """ + str(self.start_time) + """ """+ str(self.end_time) +"""
Center of the Source, xc, yc (UTM E, UTM N): """ + str(self.xcenter) + """ """+str(self.ycenter)+"""
Major and Minor Extent, majorR, minorR (m, m): """ + str(self.majradius)+ """ """ +str(self.minradius)+"""
Angle from X axis to major axis (degrees): """ +str(self.orientation)+"""
Initial speed [m/s]: """ + str(self.Vmagnitude) + """
Initial direction ([degrees] from X axis): """ + str(self.Vdirection)+"""
"""
        file.write(python_data)
        file.close

class TitanDischargePlane(cxxTitanDischargePlane):

    def __init__(self,x_a,y_a,x_b,y_b):
        super(TitanDischargePlane, self).__init__(float(x_a), float(y_a), float(x_b), float(y_b))
        #Enter Discharge Plane
        #Point A (UTM E, UTM N):
        #self.x_a = float(x_a)
        #self.y_a = float(y_a)
        #Point B (UTM E, UTM N):
        #self.x_b = float(x_b)
        #self.y_b = float(y_b)


class TitanSimulation(cxxTitanSimulation):
    possible_vizoutputs={
        'tecplotxxxx.tec':1, # first bit flag
        'mshplotxxxx.tec':2, # second bit flag
        'XDMF/Paraview':4, # third bit flag
        'grass_sites':8, # fourth bit flag
        'tecplot':1, # first bit flag
        'meshlot':2, # second bit flag
        'grasssites':8 # fourth bit flag
    }
    possible_gis_formats={
        'GIS_GRASS':cxxTitanSimulation.GIS_GRASS,
        'GDAL':cxxTitanSimulation.GDAL
    }
    possible_orders={'First':1,'Second':2}
    def __init__():
        pass
class TitanSinglePhase(cxxTitanSinglePhase):
    def __init__(self,
                 maxiter=100,
                 maxtime=1.5,
                 timeoutput=10.0,
                 timesave=None,
                 number_of_cells_across_axis=20,
                 length_scale=1.0,
                 gravity_scale=9.8,
                 height_scale=None,
                 adapt=False,
                 vizoutput="tecplotxxxx.tec",
                 order='First',
                 edge_height=None,
                 test_height=None,
                 test_location=None
                 ):
        super(TitanSinglePhase, self).__init__()
        #init values
        
        #Number of Processors
        numprocs=self.numprocs
        if numprocs <= 0:
            raise ValueError('numprocs must be greater than 0, it is ' + str(numprocs))
        if numprocs not in (1,2,4,8,12,128,256,512):
            raise ValueError('wrong amount of processors!')
        
        #Number of Computational Cells Across Smallest Pile/Flux-Source Diameter
        self.number_of_cells_across_axis = int(number_of_cells_across_axis)
        if self.number_of_cells_across_axis<=0:
            raise ValueError("TitanSimulation::number_of_cells_across_axis should be positive")
        
        #Length Scale [m]
        self.length_scale = float(length_scale)
        if self.length_scale<=0.0:
            raise ValueError("TitanSimulation::length_scale should be positive")
        #gravity scaling factor [m/s^2]
        self.gravity_scale = float(gravity_scale)
        if self.gravity_scale<=0.0:
            raise ValueError("TitanSimulation::gravity_scale should be positive")
        #height scaling factor
        if height_scale==None:
            self.height_scale=0.0
        else:
            self.height_scale = float(height_scale)
            if self.height_scale<=0.0:
                raise ValueError("TitanSimulation::height_scale should be positive")
            
        
        #Maximum Number of Time Steps
        self.maxiter = int(maxiter)
        if self.maxiter<=0:
            raise ValueError("TitanSimulation::maxiter should be positive")
        #Maximum Time [sec]
        self.maxtime = float(maxtime)
        if self.maxtime<=0.0:
            raise ValueError("TitanSimulation::maxtime should be positive")
        #Time [sec] between Results Output
        self.timeoutput = float(timeoutput)
        if self.timeoutput<=0.0:
            raise ValueError("TitanSimulation::timeoutput should be positive")
        #Time [sec] between Saves
        if timesave == None:
            self.timesave = -1.0
        else:
            self.timesave = float(timesave)
            if self.timesave<=0.0:
                raise ValueError("TitanSimulation::timesave should be positive or None")
            
        #Adapt the Grid?
        if adapt:
            self.adapt = 1
        else:
            self.adapt = 0
        #Visualization Output
        if vizoutput in TitanSimulation.possible_vizoutputs:
            self.vizoutput = TitanSimulation.possible_vizoutputs[vizoutput]
        else:
            raise ValueError("Unknown vizoutput "+str(vizoutput)+". Possible formats: "+str(possible_vizoutputs.keys()))
        
        #First/Second Order Method
        if order in TitanSimulation.possible_orders:
            self.order = TitanSimulation.possible_orders[order]
        else:
            raise ValueError("Unknown order "+str(order)+". Possible formats: "+str(possible_orders.keys()))
        
        #Test if flow reaches height [m] ...
        if edge_height == None:
            self.edge_height = -1.0
        else:
            self.edge_height = float(edge_height)
            if self.edge_height <= 0:
                raise ValueError('TitanSimulation::edge_height should be positive or None\n')
        
        #Height used to define flow outline (>0) [m]
        
        if test_height == None:
            self.test_height = -2.0
        else:
            self.test_height =float(test_height)
            if test_height <= 0:
                raise ValueError('TitanSimulation::test_height should be positive or None\n') 
        
        #... at test point (x and y location)
        if test_location==None:
            if test_height!=None:
                raise ValueError('TitanSimulation::test_location should be set if test_height>0\n') 
            self.test_location_x = 0.0
            self.test_location_y = 0.0
        else:
            if not isinstance(test_location, (list, tuple)):
                raise ValueError("Unknown format for test_location ("+str(test_location)+"). Should be [float, float]")
            if len(test_location)!=2:
                raise ValueError("Unknown format for test_location ("+str(test_location)+"). Should be [float, float]")
            self.test_location_x = float(test_location[0])
            self.test_location_y = float(test_location[1])
        
        #other inits
        self.flux_sourcesHelper=[]
        self.discharge_planesHelper=[]

    def setTopo(self,gis_format='GIS_GRASS',
                 topomain=None,
                 toposub=None,
                 topomapset=None,
                 topomap=None,
                 topovector=None,
                 min_location_x=None,
                 min_location_y=None,
                 max_location_x=None,
                 max_location_y=None):
        
        if gis_format in TitanSimulation.possible_gis_formats:
            self.gis_format = TitanSimulation.possible_gis_formats[gis_format]
            print "gis_format",gis_format,self.gis_format,TitanSimulation.possible_gis_formats[gis_format]
        else:
            raise ValueError("Unknown gis format "+str(gis_format)+". Possible formats: "+str(gis_formats[1:]))
        #GIS Information Main Directory
        self.topomain = topomain if topomain!=None else ''
        #GIS Sub-Directory
        self.toposub = toposub if toposub!=None else ''
        #GIS Map Set
        self.topomapset = topomapset if topomapset!=None else ''
        #GIS Map
        self.topomap = topomap if topomap!=None else ''
        #GIS Vector
        self.topovector = topovector if topovector!=None else ''
        
        if min_location_x!=None and min_location_y!=None and \
            max_location_x!=None and max_location_y!=None:
            self.region_limits_set=True
            #Minimum x and y location (UTM E, UTM N)
            self.min_location_x = min_location_x
            self.min_location_y = min_location_y
            #Maximum x and y location (UTM E, UTM N)
            self.max_location_x = max_location_x
            self.max_location_y = max_location_y
        else:
            self.region_limits_set=False
            #Minimum x and y location (UTM E, UTM N)
            self.min_location_x = 0.0
            self.min_location_y = 0.0
            #Maximum x and y location (UTM E, UTM N)
            self.max_location_x = 0.0
            self.max_location_y = 0.0
        
        
        
        #here should be validator
        # if there is no topo file, quit
        if self.gis_format == 1:
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
        elif  self.gis_format == 2:
            if (self.topomap == '') or (self.topomap == None) or (not isinstance(self.topomap, basestring)) or \
                    (not os.path.isdir(self.topomap)):
                raise ValueError(errmsg)
            
    def setMatMap(self,
            use_gis_matmap=False,
            matmap=None):
        #Use GIS Material Map?        
        self.use_gis_matmap = use_gis_matmap
        if matmap==None:
            matmap=[{'intfrict':30.0,'bedfrict':15.0}]
        
        #get (then write) list of material names and properties

        #if you don't enter a value for the current material it
        #defaults to the value for the previous material, so must
        #provide some previous value for the first material
        #since there is value checking in QuestionTemplate3, any
        #values can be used here
        previntfrict = 0.
        prevbedfrict = 0.
        
        if self.use_gis_matmap == False:
            self.material_map.name.push_back("all materials")
            self.material_map.intfrict.push_back(matmap[0]['intfrict'])
            self.material_map.bedfrict.push_back(matmap[0]['bedfrict'])
            #self.material_map.print0()
        else:  #if they did want to use a GIS material map...
            raise Exception("GIS material map Not implemented yet")
        #m.print0()

        fout=open("frict.data","w",0)

        #if they didn't want to use a GIS material map, get the material
        #properties once up front
        if self.use_gis_matmap == False:
            nummat=1
            fout.write(str(nummat)+'\n')
            matname = 'all materials'
            
            # check values and store them
            intfrict=float(matmap[0]['intfrict'])
            bedfrict=float(matmap[0]['bedfrict'])
            
            if intfrict <= 0.0 :
                raise ValueError('intfrict can not be negative')
            if bedfrict <= 0.0 :
                raise ValueError('bedfrict can not be negative')
            
            fout.write(matname+'\n')
            fout.write(str(intfrict)+' '+str(bedfrict)+'\n')
            
        else:  #if they did want to use a GIS material map...
            raise Exception("Not implemented as there were no suitable example")
        fout.close
    def validatePile(self, **kwargs):
        out={}
        out['height']=float(kwargs['height'])
        if out['height'] < 0.0:
            raise ValueError('TitanPile::height should be non negative')
        
        if kwargs['center']!=None:
            out['xcenter'] = float(kwargs['center'][0])
            out['ycenter'] = float(kwargs['center'][1])
        else:
            out['xcenter'] = 1.0
            out['ycenter'] = 1.0
            
        if kwargs['radii']!=None:
            out['majradius'] = float(kwargs['radii'][0])
            out['minradius'] = float(kwargs['radii'][1])
        else:
            out['majradius'] = 1.0
            out['minradius'] = 1.0
        out['orientation'] = float(kwargs['orientation'])
        out['Vmagnitude'] = float(kwargs['Vmagnitude'])
        out['Vdirection'] = float(kwargs['Vdirection'])
        return out
        
    def addPile(self,**kwargs):
        """
        Information for Pile Number 
        Thickness of Initial Volume, h(x,y)
        P*(1-((x-xc)/xr)^2 - ((y-yc)/yr)^2)
        
        height=float - Maximum Initial Thickness, P (m)
        
        center=[float,float] - Center of Initial Volume, xc, yc (UTM E, UTM N)
        
        radii=[float,float] - Major and Minor Extent, majorR, minorR (m, m)
        
        orientation=float - Orientation (angle [degrees] from X axis to major axis)
        
        Vmagnitude=float - Initial speed [m/s]
        
        Vdirection = float - Initial direction ([degrees] from X axis)
        """
        
        pile=self.validatePile(**kwargs)
        if pile!=None:
            self.pileprops.addPile(pile['height'], pile['xcenter'], pile['ycenter'], pile['majradius'], 
                                   pile['minradius'], pile['orientation'], pile['Vmagnitude'], pile['Vdirection'])
            
        
    def addFluxSource(self,**kwargs):
        fluxSource=TitanFluxSource(**kwargs)
        if fluxSource!=None:
            self.flux_sources.push_back(fluxSource)
            self.flux_sourcesHelper.append(fluxSource)
    
    def addDischargePlane(self,*args):
        dischargePlane=TitanDischargePlane(*args)
        if dischargePlane!=None:
            self.discharge_planes.push_back(dischargePlane)
            self.discharge_planesHelper.append(dischargePlane)
        
    def preproc(self):
        if self.myid==0:
            preproc=TitanPreproc(self)
            preproc.validate();
            preproc.run();
        
    def run(self):
        max_height=0.0
        for iPile in range(self.pileprops.numpiles):
            if self.pileprops.pileheight[iPile] > max_height:
                max_height = self.pileprops.pileheight[iPile]
        heightscale = max_height
        for i in range(len(self.flux_sourcesHelper)):
            self.flux_sourcesHelper[i].done(f_p)
            effective_height=self.flux_sources[i].get_effective_height()
            if effective_height > max_height:
                max_height = effective_height
        #scaling stuff
        if max_height <= 0.:max_height = 1.
        if heightscale <= 0.:heightscale = 1.
        
        if self.myid==0:
            print 'max height is ' + str(max_height)
            print 'heightscale is ' + str(heightscale)
        
        if self.myid==0:
            # run preproc.x to create the fem grid, if it is not already there
            #if os.access('PRE/preproc.x',os.X_OK)==0:
            #    os.system('cd PRE;gmake')
            # locate titan_preprocess: 
            # first check the local-diretory,
            
            print "\npreproc..."
            self.preproc()
    
            if self.topovector != "":
                print "\nvectordatpreproc..."
                vectordatpreproc(self.topomain, self.toposub, self.topomapset, self.topomap, self.topovector)
            print
        
        super(TitanSinglePhase, self).run()

class TitanTwoPhases(TitanSinglePhase):
    def __init__(self):
        super(TitanSinglePhase, self).__init__()
        self.pileprops2=PilePropsTwoPhases()
    def validatePile(self, **kwargs):
        out=super(TitanTwoPhases, self).validatePile(**kwargs)
        out['vol_fract'] = float(kwargs['vol_fract'])
        return out
        
    def addPile(self,**kwargs):
        """
        Information for Pile Number 
        Thickness of Initial Volume, h(x,y)
        P*(1-((x-xc)/xr)^2 - ((y-yc)/yr)^2)
        
        height=float - Maximum Initial Thickness, P (m)
        
        center=[float,float] - Center of Initial Volume, xc, yc (UTM E, UTM N)
        
        radii=[float,float] - Major and Minor Extent, majorR, minorR (m, m)
        
        orientation=float - Orientation (angle [degrees] from X axis to major axis)
        
        Vmagnitude=float - Initial speed [m/s]
        
        Vdirection = float - Initial direction ([degrees] from X axis)
        
        vol_fract = float - Initial solid-volume fraction,(0:1.)
        """
        
        pile=self.validatePile(**kwargs)
        if pile!=None:
            self.pileprops2.addPile(pile['height'], pile['xcenter'], pile['ycenter'], pile['majradius'], 
                                   pile['minradius'], pile['orientation'], pile['Vmagnitude'], pile['Vdirection'],pile['vol_fract'])
    
