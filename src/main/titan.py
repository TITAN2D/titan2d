
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
import copy
import inspect

from cxxtitan import *


class TitanSimulationBase(object):
    """ Base class for TitanSimulation
    it defines class members and static variables  
    """
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
    
    possible_pile_types={
         'PARABALOID':PileProps.PARABALOID,
         'CYLINDER':PileProps.CYLINDER,
         'PLANE':PileProps.PLANE,
         'CASITA':PileProps.CASITA,
         'POPO':PileProps.POPO,
         'ID1':PileProps.ID1,
         'ID2':PileProps.ID2
    }
    
    possible_internal_mat_models={
        'Coulomb':{
            'allParameters':('order','int_frict'),
            'defaultParameters':{'order':'First','int_frict':37.0},
            'elementType':ElementType_SinglePhase,
            'integrators':[
                {
                    'conditions' :[lambda tsim,model_parameters: model_parameters['order']==1 or model_parameters['order']==2],
                    'constructor':Integrator_SinglePhase_Coulomb
                }
            ]
        },
        'Vollmey':{
            'allParameters':('order','mu','xi','int_frict'), 
            'defaultParameters':{
                'order':'First',
                'mu' : 0.5,
                'xi' : 120.0,
                'int_frict':37.0
            },
            'elementType':ElementType_SinglePhase,
            'integrators':[{
                    'conditions' :[lambda tsim,model_parameters: model_parameters['order']==1],
                    'constructor':Integrator_SinglePhase_Vollmey_FirstOrder
            }]
        },
        'Pouliquen':{
            'allParameters':('order','phi1','phi2','partdiam','I_O'), 
            'defaultParameters':{
                'order':'First',
                'phi1':24.0,
                'phi2':30.0,
                'partdiam':1.0E-4,
                'I_O':0.3
            },
            'elementType':ElementType_SinglePhase,
            'integrators':[{
                    'conditions' :[lambda tsim,model_parameters: model_parameters['order']==1],
                    'constructor':Integrator_SinglePhase_Pouliquen_FirstOrder
            }]
        },
        'Maeno':{
            'allParameters':('order','phis', 'phi2','partdiam', 'I_not'),
            'defaultParameters':{
                'order':'First',
                'phis':24.0,
                'phi2':30.0,
                'partdiam':1.0E-4,
                'I_not':0.3
            },
            'elementType':ElementType_SinglePhase,
            'integrators':[{
                    'conditions' :[lambda tsim,model_parameters: model_parameters['order']==1],
                    'constructor':Integrator_SinglePhase_Maeno_FirstOrder
            }],
        },
        'TwoPhases_Coulomb':{
            'allParameters':('order','int_frict',), 
            'defaultParameters':{'order':'First','int_frict':37.0},
            'elementType':ElementType_TwoPhases,
            'integrators':[{
                    'conditions' :[lambda tsim,model_parameters: model_parameters['order']==1],
                    'constructor':Integrator_TwoPhases_Coulomb
            }]
        }
    }
    
    
    
    def __init__(self):
        #initiate all class members
        
        #referece to cxxTitanSimulation, the class which handles calculations on c++ side 
        self.sim=None
        
        #flag that integrator and internal material model is initialized        
        self.integrator_initialized=False
        self.gis_initialized=False
        
        self.pileprops=None
        
    def input_summary(self):
        if self.sim!=None:
            self.sim.input_summary()
        
class TitanSimulation(TitanSimulationBase):
    def __init__(self,
                 max_iter=100,
                 max_time=1.5,
                 time_output=10.0,
                 time_save=None,
                 length_scale=1.0,
                 gravity_scale=9.8,
                 height_scale=None,
                 adapt=False,
                 short_speed=False,
                 vizoutput="tecplotxxxx.tec",
                 edge_height=None,
                 test_height=None,
                 test_location=None,
                 sim_class=None
                 ):
        super(TitanSimulation, self).__init__()
        if sim_class==None:
            self.sim=cxxTitanSimulation()
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
        self.sim.scale_.length = float(length_scale)
        if self.sim.scale_.length<=0.0:
            raise ValueError("TitanSimulation::length_scale should be positive")
        #gravity scaling factor [m/s^2]
        self.sim.scale_.gravity = float(gravity_scale)
        if self.sim.scale_.gravity<=0.0:
            raise ValueError("TitanSimulation::gravity_scale should be positive")
        #height scaling factor
        if height_scale==None:
            self.sim.scale_.height=0.0
        else:
            self.sim.scale_.height = float(height_scale)
            self.sim.scale_.auto_calc_height_scale=False
            if self.sim.scale_.height<=0.0:
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
            
        #short_speed
        self.sim.set_short_speed(short_speed)
        #Visualization Output
        if vizoutput in TitanSimulation.possible_vizoutputs:
            self.sim.vizoutput = TitanSimulation.possible_vizoutputs[vizoutput]
        else:
            raise ValueError("Unknown vizoutput "+str(vizoutput)+". Possible formats: "+str(possible_vizoutputs.keys()))
        
        
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
    def setIntMatModel(self,**karg):
        """
        Setup internal material model and respective integrator
        Coulomb is default model
        """
        #
        model=karg.pop('model','Coulomb');
        if model not in TitanSimulationBase.possible_internal_mat_models.keys():
            raise ValueError("Unknown internal material model "+str(vizoutput)+". Possible models: "+str(TitanSimulationBase.possible_internal_mat_models.keys()))
        integrator=karg.pop('integrator',None);
        
        #initiate model parameters with default values
        model_parameters=copy.deepcopy(TitanSimulationBase.possible_internal_mat_models[model]['defaultParameters'])
        #get the list of all parameters for that model
        model_all_parameters=TitanSimulationBase.possible_internal_mat_models[model]['allParameters']
        
        #set user specified parameters        
        while len(karg):
            (key,value)=karg.popitem()
            if key not in model_all_parameters:
                raise ValueError("Unknown parameter '"+str(key)+"' for "+str(model)+" model. Possible parameters: "+str(model_all_parameters) )
            model_parameters[key]=value
        
        #check that all parameters present
        for key in model_all_parameters:
            if key not in model_parameters:
                raise ValueError("Parameter '"+str(key)+"' is not set for "+str(model)+" model!")
        
        #First/Second Order Method
        if model_parameters['order'] in TitanSimulation.possible_orders:
            model_parameters['order'] = TitanSimulation.possible_orders[model_parameters['order']]
        else:
            raise ValueError("Unknown order "+str(model_parameters['order'])+". Possible formats: "+str(possible_orders.keys()))
        
        #find integrator        
        if integrator==None:
            integrators=TitanSimulationBase.possible_internal_mat_models[model]['integrators']
            msg=""
            
            for desc in integrators:
                does_feat=True
                msg+="Testing: "+desc['constructor'].__name__+"\n"
                for cond in desc['conditions']:
                    cond_result=cond(self,model_parameters)
                    msg+="  "+inspect.getsource(cond).strip()+" ===>"+str(cond_result)+"\n"
                    does_feat=does_feat and cond_result
                if does_feat:
                    integrator=desc['constructor']
                    break
            if integrator==None:
                raise ValueError("Can not find suitable integrator, here is the hint\n"+msg+"\ncheck the manual.")
        
        #set element type
        elementType=TitanSimulationBase.possible_internal_mat_models[model]['elementType']
        self.sim.set_element_type(elementType)
        #construct integrator
        integrator_obj=integrator(self.sim)
        #set parameters
        for k,v in model_parameters.iteritems():
            setattr(integrator_obj,k,v)
        
        
        self.sim.set_integrator(integrator_obj)
        self.integrator_initialized=True
        #uncouple the c++ from python proxy 
        integrator_obj.thisown=0
        
    def setMatMap(self,
            use_gis_matmap=False,
            number_of_cells_across_axis=20,
            bed_frict=15.0,
            mat_names=None):
        if not self.integrator_initialized:
            raise ValueError("Internal material model should be set before setMatMap!(see setIntMatModel")
        
        #Use GIS Material Map?        
        self.sim.use_gis_matmap = use_gis_matmap
        
        #here we need to set mat prop
        if self.sim.get_element_type()==ElementType_SinglePhase:
            matprops=MatProps(self.sim.scale_)
        elif self.sim.get_element_type()==ElementType_TwoPhases:
            matprops=MatPropsTwoPhases(self.sim.scale_)
        
        #Number of Computational Cells Across Smallest Pile/Flux-Source Diameter
        matprops.number_of_cells_across_axis = int(number_of_cells_across_axis)
        if matprops.number_of_cells_across_axis<=0:
            raise ValueError("TitanSimulation::number_of_cells_across_axis should be positive")
        
        
        if not isinstance(bed_frict, (list, tuple)):
            bed_frict=[bed_frict]
        
        if self.sim.use_gis_matmap == False:
            matprops.material_count=1
            matprops.matnames.push_back("all materials")
            matprops.bedfrict.push_back(float(bed_frict[0]))
        else:  #if they did want to use a GIS material map...
            matprops.material_count=len(mat_names)
            if len(bed_frict)!=len(mat_names):
                raise Exception("number of mat_names does not match number of bed_frict")
            for i in range(len(bed_frict)):
                matprops.matnames.push_back(mat_names[i])
                matprops.bedfrict.push_back(float(bed_frict[i]))
            raise Exception("GIS material map Not implemented yet")
        
        self.sim.set_matprops(matprops)
        #uncouple the c++ from python proxy 
        matprops.thisown=0 
        

        
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
        
        pile_type = string - ['PARABALOID', 'CYLINDER']
        
        vol_fract = float - Initial solid-volume fraction,(0:1.) [for two phases]
        """
        if not self.integrator_initialized:
            raise ValueError("Internal material model should be set before adding piles!(see setIntMatModel")
        
        
        #couple helping functions
        def validateSinglePhasePile(**kwargs):
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
            
            if kwargs['pile_type'] in TitanSimulationBase.possible_pile_types:
                out['pile_type'] = TitanSimulationBase.possible_pile_types[kwargs['pile_type']]
            else:
                raise ValueError("Unknown pile_type "+str(pile_type)+". Possible formats: "+str(possible_pile_types.keys()))
                
            return out
        
        def validateTwoPhasesPile(**kwargs):
            out=validateSinglePhasePile(**kwargs)
            #if 'pile_type' not in kwargs or kwargs['pile_type']==None:
            #    kwargs['pile_type']='PARABALOID'
            out['vol_fract'] = float(kwargs['vol_fract'])
            return out
        
        if self.sim.get_element_type()==ElementType_SinglePhase:
            if self.pileprops==None:
                self.pileprops=PileProps()
                self.pileprops.thisown=0
                self.sim.set_pileprops(self.pileprops)
            if not isinstance(self.pileprops,PileProps):
                raise ValueError("Can not mix element type for piles!")
            if isinstance(self.pileprops,PilePropsTwoPhases):
                raise ValueError("Can not mix element type for piles!")
            
            pile=validateSinglePhasePile(**kwargs)
            if pile!=None:
                self.pileprops.addPile(pile['height'], pile['xcenter'], pile['ycenter'], pile['majradius'], 
                                       pile['minradius'], pile['orientation'], pile['Vmagnitude'], pile['Vdirection'],pile['pile_type'])
        elif self.sim.get_element_type()==ElementType_TwoPhases:
            if self.pileprops==None:
                self.pileprops=PilePropsTwoPhases()
                self.pileprops.thisown=0
                self.sim.set_pileprops(self.pileprops)
            if not isinstance(self.pileprops,PileProps):
                raise ValueError("Can not mix element type for piles!")
            if not isinstance(self.pileprops,PilePropsTwoPhases):
                raise ValueError("Can not mix element type for piles!")
            
            pile=validateTwoPhasesPile(**kwargs)
            if pile!=None:
                self.pileprops.addPile(pile['height'], pile['xcenter'], pile['ycenter'], pile['majradius'], 
                                       pile['minradius'], pile['orientation'], pile['Vmagnitude'], pile['Vdirection'],pile['pile_type'],pile['vol_fract'])
        else:
            raise ValueError("Unknown element type")
    

            
    
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
            effective_height=self.sim.fluxprops.get_effective_height(i)
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
