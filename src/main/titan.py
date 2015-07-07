
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

class TitanSimulation(object):
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
    def __init__(self):
        self.sim=None
    def input_summary(self):
        if self.sim!=None:
            self.sim.input_summary()
        
class TitanSinglePhase(TitanSimulation):
    def __init__(self,
                 maxiter=100,
                 maxtime=1.5,
                 timeoutput=10.0,
                 timesave=None,
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
        self.sim=cxxTitanSinglePhase()
        #init values
        
        #Number of Processors
        numprocs=self.sim.numprocs
        if numprocs <= 0:
            raise ValueError('numprocs must be greater than 0, it is ' + str(numprocs))
        if numprocs not in (1,2,4,8,12,128,256,512):
            raise ValueError('wrong amount of processors!')
        
        
        
        #Length Scale [m]
        self.sim.length_scale = float(length_scale)
        if self.sim.length_scale<=0.0:
            raise ValueError("TitanSimulation::length_scale should be positive")
        #gravity scaling factor [m/s^2]
        self.sim.gravity_scale = float(gravity_scale)
        if self.sim.gravity_scale<=0.0:
            raise ValueError("TitanSimulation::gravity_scale should be positive")
        #height scaling factor
        if height_scale==None:
            self.sim.height_scale=0.0
        else:
            self.sim.height_scale = float(height_scale)
            if self.sim.height_scale<=0.0:
                raise ValueError("TitanSimulation::height_scale should be positive")
            
        
        #Maximum Number of Time Steps
        self.sim.maxiter = int(maxiter)
        if self.sim.maxiter<=0:
            raise ValueError("TitanSimulation::maxiter should be positive")
        #Maximum Time [sec]
        self.sim.maxtime = float(maxtime)
        if self.sim.maxtime<=0.0:
            raise ValueError("TitanSimulation::maxtime should be positive")
        #Time [sec] between Results Output
        self.sim.timeoutput = float(timeoutput)
        if self.sim.timeoutput<=0.0:
            raise ValueError("TitanSimulation::timeoutput should be positive")
        #Time [sec] between Saves
        if timesave == None:
            self.sim.timesave = -1.0
        else:
            self.sim.timesave = float(timesave)
            if self.sim.timesave<=0.0:
                raise ValueError("TitanSimulation::timesave should be positive or None")
            
        #Adapt the Grid?
        if adapt:
            self.sim.adapt = 1
        else:
            self.sim.adapt = 0
        #Visualization Output
        if vizoutput in TitanSimulation.possible_vizoutputs:
            self.sim.vizoutput = TitanSimulation.possible_vizoutputs[vizoutput]
        else:
            raise ValueError("Unknown vizoutput "+str(vizoutput)+". Possible formats: "+str(possible_vizoutputs.keys()))
        
        #First/Second Order Method
        if order in TitanSimulation.possible_orders:
            self.sim.order = TitanSimulation.possible_orders[order]
        else:
            raise ValueError("Unknown order "+str(order)+". Possible formats: "+str(possible_orders.keys()))
        
        #Test if flow reaches height [m] ...
        if edge_height == None:
            self.sim.edge_height = -1.0
        else:
            self.sim.edge_height = float(edge_height)
            if self.sim.edge_height <= 0:
                raise ValueError('TitanSimulation::edge_height should be positive or None\n')
        
        #Height used to define flow outline (>0) [m]
        
        if test_height == None:
            self.sim.test_height = -2.0
        else:
            self.sim.test_height =float(test_height)
            if test_height <= 0:
                raise ValueError('TitanSimulation::test_height should be positive or None\n') 
        
        #... at test point (x and y location)
        if test_location==None:
            if test_height!=None:
                raise ValueError('TitanSimulation::test_location should be set if test_height>0\n') 
            self.sim.test_location_x = 0.0
            self.sim.test_location_y = 0.0
        else:
            if not isinstance(test_location, (list, tuple)):
                raise ValueError("Unknown format for test_location ("+str(test_location)+"). Should be [float, float]")
            if len(test_location)!=2:
                raise ValueError("Unknown format for test_location ("+str(test_location)+"). Should be [float, float]")
            self.sim.test_location_x = float(test_location[0])
            self.sim.test_location_y = float(test_location[1])
        
        #other inits
        self.sim.flux_sourcesHelper=[]
        self.sim.discharge_planesHelper=[]

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
            self.sim.gis_format = TitanSimulation.possible_gis_formats[gis_format]
        else:
            raise ValueError("Unknown gis format "+str(gis_format)+". Possible formats: "+str(gis_formats[1:]))
        #GIS Information Main Directory
        self.sim.topomain = topomain if topomain!=None else ''
        #GIS Sub-Directory
        self.sim.toposub = toposub if toposub!=None else ''
        #GIS Map Set
        self.sim.topomapset = topomapset if topomapset!=None else ''
        #GIS Map
        self.sim.topomap = topomap if topomap!=None else ''
        #GIS Vector
        self.sim.topovector = topovector if topovector!=None else ''
        
        if min_location_x!=None and min_location_y!=None and \
            max_location_x!=None and max_location_y!=None:
            self.sim.region_limits_set=True
            #Minimum x and y location (UTM E, UTM N)
            self.sim.min_location_x = min_location_x
            self.sim.min_location_y = min_location_y
            #Maximum x and y location (UTM E, UTM N)
            self.sim.max_location_x = max_location_x
            self.sim.max_location_y = max_location_y
        else:
            self.sim.region_limits_set=False
            #Minimum x and y location (UTM E, UTM N)
            self.sim.min_location_x = 0.0
            self.sim.min_location_y = 0.0
            #Maximum x and y location (UTM E, UTM N)
            self.sim.max_location_x = 0.0
            self.sim.max_location_y = 0.0
        
        
        
        #here should be validator
        # if there is no topo file, quit
        if self.sim.gis_format == 1:
            errmsg='Missing GIS information.  No job will be run.'
            if self.sim.topomain == '' or self.sim.toposub == '' or self.sim.topomapset == '' or self.sim.topomap == '':
                raise ValueError(errmsg)
            if self.sim.topomain == None or self.sim.toposub == None or self.sim.topomapset == None or self.sim.topomap == None:
                raise ValueError(errmsg)
            if (not isinstance(self.sim.topomain, basestring)) or (not isinstance(self.sim.toposub, basestring)) or \
                    (not isinstance(self.sim.topomapset, basestring)) or (not isinstance(self.sim.topomap, basestring)):
                raise ValueError(errmsg)
            
            p=""
            for dirname in (self.sim.topomain,self.sim.toposub,self.sim.topomapset):
                p=os.path.join(p,dirname)
                print p
                if not os.path.isdir(p):
                    raise ValueError(errmsg+". "+p+" does not exist!")
            #self.sim.topomap?
        elif  self.sim.gis_format == 2:
            if (self.sim.topomap == '') or (self.sim.topomap == None) or (not isinstance(self.sim.topomap, basestring)) or \
                    (not os.path.isdir(self.sim.topomap)):
                raise ValueError(errmsg)
            
    def setMatMap(self,
            use_gis_matmap=False,
            number_of_cells_across_axis=20,
            intfrict=30.0,
            bedfrict=15.0,
            mat_names=None):
        #Use GIS Material Map?        
        self.sim.use_gis_matmap = use_gis_matmap
        #Number of Computational Cells Across Smallest Pile/Flux-Source Diameter
        self.sim.matprops.number_of_cells_across_axis = int(number_of_cells_across_axis)
        if self.sim.matprops.number_of_cells_across_axis<=0:
            raise ValueError("TitanSimulation::number_of_cells_across_axis should be positive")
        
        
        if not isinstance(bedfrict, (list, tuple)):
            bedfrict=[bedfrict]
        
        self.sim.matprops.intfrict=float(intfrict)
        if self.sim.use_gis_matmap == False:
            self.sim.matprops.material_count=1
            self.sim.matprops.matnames.push_back("all materials")
            self.sim.matprops.bedfrict.push_back(float(bedfrict[0]))
            #self.sim.material_map.print0()
        else:  #if they did want to use a GIS material map...
            self.sim.matprops.material_count=len(mat_names)
            if len(bedfrict)!=len(mat_names):
                raise Exception("number of mat_names does not match number of bedfrict")
            for i in range(len(bedfrict)):
                self.sim.matprops.matnames.push_back(mat_names[i])
                self.sim.matprops.bedfrict.push_back(float(bedfrict[i]))
            raise Exception("GIS material map Not implemented yet")
        
    def validatePile(self, **kwargs):
        out={}
        out['height']=float(kwargs['height'])
        if out['height'] < 0.0:
            raise ValueError('TitanPile::height should be non negative')
        
        if 'center' in kwargs and kwargs['center']!=None:
            out['xcenter'] = float(kwargs['center'][0])
            out['ycenter'] = float(kwargs['center'][1])
        else:
            out['xcenter'] = 1.0
            out['ycenter'] = 1.0
            
        if 'radii' in kwargs and kwargs['radii']!=None:
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
            self.sim.pileprops.addPile(pile['height'], pile['xcenter'], pile['ycenter'], pile['majradius'], 
                                   pile['minradius'], pile['orientation'], pile['Vmagnitude'], pile['Vdirection'])
            
    
    def validateFluxSource(self,influx,start_time,end_time,center,radii,
                 orientation,
                 Vmagnitude,
                 Vdirection):
        out={}
        #Extrusion flux rate [m/s]
        out['influx'] = float(influx)
        
        #Active Time [s], start, end
        out['start_time'] = float(start_time)
        out['end_time']   = float(end_time)
        #Center of the source, xc, yc (UTM E, UTM N):
        if center!=None:
            out['xcenter'] = float(center[0])
            out['ycenter'] = float(center[1])
        else:
            out['xcenter'] = 1.0
            out['ycenter'] = 1.0
        #Major and Minor Extent, majorR, minorR (m, m)
        if radii!=None:
            out['majradius'] = float(radii[0])
            out['minradius'] = float(radii[1])
        else:
            out['majradius'] = 1.0
            out['minradius'] = 1.0
        #Orientation (angle [degrees] from X axis to major axis):
        out['orientation'] = float(orientation)
        #Initial speed [m/s]:
        out['Vmagnitude'] = float(Vmagnitude)
        #Initial direction ([degrees] from X axis):
        out['Vdirection'] = float(Vdirection)
        
        if out['influx'] < 0.0:
            raise ValueError('TitanFluxSource::influx should be non negative')
        if out['start_time'] < 0.0:
            raise ValueError('TitanFluxSource::start_time should be non negative')
        if out['end_time'] < 0.0:
            raise ValueError('TitanFluxSource::start_time should be non negative')    
        return out
    def addFluxSource(self,influx,start_time,end_time,center=None,radii=None,
                 orientation=0.0,
                 Vmagnitude=0.0,
                 Vdirection=0.0):
        out=self.validateFluxSource(influx,start_time,end_time,center,radii,
                 orientation, Vmagnitude, Vdirection)
        if out!=None:
            self.sim.fluxprops.addFluxSource(out['influx'],out['start_time'],out['end_time'], out['xcenter'],out['ycenter'],
                                             out['majradius'],out['minradius'],out['orientation'],out['Vmagnitude'],out['Vdirection'])
    
    def addDischargePlane(self,x_a,y_a,x_b,y_b):
        self.sim.discharge_planes.addDischargePlane(float(x_a), float(y_a), float(x_b), float(y_b))
        
    def preproc(self):
        if self.sim.myid==0:
            preproc=TitanPreproc(self.sim)
            preproc.validate();
            preproc.run();
        
    def run(self):
        max_height=0.0
        for iPile in range(self.sim.pileprops.numpiles):
            if self.sim.pileprops.pileheight[iPile] > max_height:
                max_height = self.sim.pileprops.pileheight[iPile]
        heightscale = max_height
        for i in range(len(self.sim.flux_sourcesHelper)):
            self.sim.flux_sourcesHelper[i].done(f_p)
            effective_height=self.sim.flux_sources[i].get_effective_height()
            if effective_height > max_height:
                max_height = effective_height
        #scaling stuff
        if max_height <= 0.:max_height = 1.
        if heightscale <= 0.:heightscale = 1.
        
        if self.sim.myid==0:
            print 'max height is ' + str(max_height)
            print 'heightscale is ' + str(heightscale)
        
        if self.sim.myid==0:
            # run preproc.x to create the fem grid, if it is not already there
            #if os.access('PRE/preproc.x',os.X_OK)==0:
            #    os.system('cd PRE;gmake')
            # locate titan_preprocess: 
            # first check the local-diretory,
            
            print "\npreproc..."
            self.preproc()
    
            if self.sim.topovector != "":
                print "\nvectordatpreproc..."
                vectordatpreproc(self.sim.topomain, self.sim.toposub, self.sim.topomapset, self.sim.topomap, self.sim.topovector)
            print
        
        self.sim.run()

class TitanTwoPhases(TitanSinglePhase):
    def __init__(self, **kwargs):
        super(TitanTwoPhases, self).__init__(**kwargs)
        self.sim=cxxTitanTwoPhases()
        
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
            self.sim.pileprops.addPile(pile['height'], pile['xcenter'], pile['ycenter'], pile['majradius'], 
                                   pile['minradius'], pile['orientation'], pile['Vmagnitude'], pile['Vdirection'],pile['vol_fract'])
    
