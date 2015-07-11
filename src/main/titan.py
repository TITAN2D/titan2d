
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
        'GIS_GRASS':MapNames.GIS_GRASS,
        'GDAL':MapNames.GDAL
    }
    possible_orders={'First':1,'Second':2}
    
    possible_pile_types={'PARABALOID':PileProps.PARABALOID,
                         'CYLINDER':PileProps.CYLINDER,
                         'PLANE':PileProps.PLANE,
                         'CASITA':PileProps.CASITA,
                         'POPO':PileProps.POPO,
                         'ID1':PileProps.ID1,
                         'ID2':PileProps.ID2
                         }
    
    def __init__(self):
        self.sim=None
    def input_summary(self):
        if self.sim!=None:
            self.sim.input_summary()
        
class TitanSinglePhase(TitanSimulation):
    def __init__(self,
                 max_iter=100,
                 max_time=1.5,
                 time_output=10.0,
                 time_save=None,
                 length_scale=1.0,
                 gravity_scale=9.8,
                 height_scale=None,
                 adapt=False,
                 vizoutput="tecplotxxxx.tec",
                 order='First',
                 edge_height=None,
                 test_height=None,
                 test_location=None,
                 sim_class=None
                 ):
        super(TitanSinglePhase, self).__init__()
        if sim_class==None:
            self.sim=cxxTitanSinglePhase()
        else:
            self.sim=sim_class
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
        max_iter = int(max_iter)
        if max_iter<=0:
            raise ValueError("TitanSimulation::max_iter should be positive")
        #Maximum Time [sec]
        max_time = float(max_time)
        if max_time<=0.0:
            raise ValueError("TitanSimulation::max_time should be positive")
        #Time [sec] between Results Output
        time_output = float(time_output)
        if time_output<=0.0:
            raise ValueError("TitanSimulation::time_output should be positive")
        #Time [sec] between Saves
        if time_save == None:
            time_save = -1.0
        else:
            time_save = float(time_save)
            if time_save<=0.0:
                raise ValueError("TitanSimulation::time_save should be positive or None")
        self.sim.get_timeprops().set_time(max_iter,max_time,time_output,time_save)
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
            edge_height = -1.0
        else:
            edge_height = float(edge_height)
            if edge_height <= 0:
                raise ValueError('TitanSimulation::edge_height should be positive or None\n')
        
        #Height used to define flow outline (>0) [m]
        
        if test_height == None:
            test_height = -2.0
            test_location_x = 0.0
            test_location_y = 0.0
        else:
            test_height =float(test_height)
            if test_height <= 0:
                raise ValueError('TitanSimulation::test_height should be positive or None\n')
            #... at test point (x and y location)
            if test_location==None:
                raise ValueError('TitanSimulation::test_location should be set if test_height>0\n') 
            else:
                if not isinstance(test_location, (list, tuple)):
                    raise ValueError("Unknown format for test_location ("+str(test_location)+"). Should be [float, float]")
                if len(test_location)!=2:
                    raise ValueError("Unknown format for test_location ("+str(test_location)+"). Should be [float, float]")
                test_location_x = float(test_location[0])
                test_location_y = float(test_location[1])
        
        self.sim.get_statprops().set(edge_height, test_height, test_location_x, test_location_y);
        #other inits

    def setGIS(self,gis_format='GIS_GRASS',
                 gis_main=None,
                 gis_sub=None,
                 gis_mapset=None,
                 gis_map=None,
                 gis_vector=None,
                 min_location_x=None,
                 min_location_y=None,
                 max_location_x=None,
                 max_location_y=None):
        
        if gis_format in TitanSimulation.possible_gis_formats:
            gis_format = TitanSimulation.possible_gis_formats[gis_format]
        else:
            raise ValueError("Unknown gis gis_format "+str(gis_format)+". Possible formats: "+str(gis_formats[1:]))
        #GIS Information Main Directory
        gis_main = gis_main if gis_main!=None else ''
        #GIS Sub-Directory
        gis_sub = gis_sub if gis_sub!=None else ''
        #GIS Map Set
        set = gis_mapset if gis_mapset!=None else ''
        #GIS Map
        gis_map = gis_map if gis_map!=None else ''
        #GIS Vector
        gis_vector = gis_vector if gis_vector!=None else ''
        
        self.sim.get_mapnames().set(gis_format,gis_main, gis_sub, gis_mapset, gis_map, gis_vector, 0)
        
        
        if min_location_x!=None and min_location_y!=None and \
            max_location_x!=None and max_location_y!=None:
            
            #Minimum x and y location (UTM E, UTM N)
            #Maximum x and y location (UTM E, UTM N)
            self.sim.get_mapnames().set_region_limits(float(min_location_x),float(max_location_x),float(min_location_y),float(max_location_y))
            
        
        
        
        #here should be validator
        # if there is no topo file, quit
        if gis_format == 1:
            errmsg='Missing GIS information.  No job will be run.'
            if gis_main == '' or gis_sub == '' or gis_mapset == '' or gis_map == '':
                raise ValueError(errmsg)
            if gis_main == None or gis_sub == None or gis_mapset == None or gis_map == None:
                raise ValueError(errmsg)
            if (not isinstance(gis_main, basestring)) or (not isinstance(gis_sub, basestring)) or \
                    (not isinstance(gis_mapset, basestring)) or (not isinstance(gis_map, basestring)):
                raise ValueError(errmsg)
            
            p=""
            for dirname in (gis_main,gis_sub,gis_mapset):
                p=os.path.join(p,dirname)
                if not os.path.isdir(p):
                    raise ValueError(errmsg+". "+p+" does not exist!")
            #gis_map?
        elif  self.sim.gis_format == 2:
            if (gis_map == '') or (gis_map == None) or (not isinstance(gis_map, basestring)) or \
                    (not os.path.isdir(gis_map)):
                raise ValueError(errmsg)
            
    def setMatMap(self,
            use_gis_matmap=False,
            number_of_cells_across_axis=20,
            int_frict=30.0,
            bed_frict=15.0,
            mat_names=None):
        #Use GIS Material Map?        
        self.sim.use_gis_matmap = use_gis_matmap
        #Number of Computational Cells Across Smallest Pile/Flux-Source Diameter
        self.sim.get_matprops().number_of_cells_across_axis = int(number_of_cells_across_axis)
        if self.sim.get_matprops().number_of_cells_across_axis<=0:
            raise ValueError("TitanSimulation::number_of_cells_across_axis should be positive")
        
        
        if not isinstance(bed_frict, (list, tuple)):
            bed_frict=[bed_frict]
        
        self.sim.get_matprops().intfrict=float(int_frict)
        if self.sim.use_gis_matmap == False:
            self.sim.get_matprops().material_count=1
            self.sim.get_matprops().matnames.push_back("all materials")
            self.sim.get_matprops().bedfrict.push_back(float(bed_frict[0]))
        else:  #if they did want to use a GIS material map...
            self.sim.get_matprops().material_count=len(mat_names)
            if len(bed_frict)!=len(mat_names):
                raise Exception("number of mat_names does not match number of bed_frict")
            for i in range(len(bed_frict)):
                self.sim.get_matprops().matnames.push_back(mat_names[i])
                self.sim.get_matprops().bedfrict.push_back(float(bed_frict[i]))
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
        
        if 'pile_type' not in kwargs or kwargs['pile_type']==None:
            kwargs['pile_type']='CYLINDER'
        
        if kwargs['pile_type'] in TitanSimulation.possible_pile_types:
            out['pile_type'] = TitanSimulation.possible_pile_types[kwargs['pile_type']]
        else:
            raise ValueError("Unknown pile_type "+str(pile_type)+". Possible formats: "+str(possible_pile_types.keys()))
            
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
            self.sim.pileprops_single_phase.addPile(pile['height'], pile['xcenter'], pile['ycenter'], pile['majradius'], 
                                   pile['minradius'], pile['orientation'], pile['Vmagnitude'], pile['Vdirection'],pile['pile_type'])
            
    
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
        for iPile in range(self.sim.get_pileprops().numpiles):
            if self.sim.get_pileprops().pileheight[iPile] > max_height:
                max_height = self.sim.get_pileprops().pileheight[iPile]
        heightscale = max_height
        for i in range(self.sim.fluxprops.no_of_sources):
            effective_height=self.sim.flux_sources[i].get_effective_height()
            if effective_height > max_height:
                max_height = effective_height
        #scaling stuff
        if max_height <= 0.:max_height = 1.
        if heightscale <= 0.:heightscale = 1.
        
        if self.sim.myid==0:
            print 'max height is ' + str(max_height)
            #print 'heightscale based on max height is ' + str(heightscale)
        
        if self.sim.myid==0:
            # run preproc.x to create the fem grid, if it is not already there
            #if os.access('PRE/preproc.x',os.X_OK)==0:
            #    os.system('cd PRE;gmake')
            # locate titan_preprocess: 
            # first check the local-diretory,
            
            print "\npreproc..."
            self.preproc()
    
            if self.sim.get_mapnames().gis_vector != "":
                print "\nvectordatpreproc..."
                vectordatpreproc(self.sim.get_mapnames().gis_main, self.sim.get_mapnames().gis_sub, self.sim.get_mapnames().gis_mapset, self.sim.get_mapnames().gis_map, self.sim.get_mapnames().gis_vector)
            print
        #self.sim.input_summary()
        self.sim.run()

class TitanTwoPhases(TitanSinglePhase):
    def __init__(self, **kwargs):
        kwargs['sim_class']=cxxTitanTwoPhases()
        super(TitanTwoPhases, self).__init__(**kwargs)
        
    def validatePile(self, **kwargs):
        if 'pile_type' not in kwargs or kwargs['pile_type']==None:
            kwargs['pile_type']='PARABALOID'
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
            self.sim.pileprops_two_phases.addPile(pile['height'], pile['xcenter'], pile['ycenter'], pile['majradius'], 
                                   pile['minradius'], pile['orientation'], pile['Vmagnitude'], pile['Vdirection'],pile['pile_type'],pile['vol_fract'])
    
