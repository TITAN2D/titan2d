
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
import shutil
import copy
import inspect
import glob

from cxxtitan import *


def VarTypeFloat(sectionName,varName,value):
    try:
        float(value)
    except:
        raise ValueError(sectionName+": incorrect data type for argument "+varName+"="+str(value)+", should be float!")
    return float(value)
def VarTypeInt(sectionName,varName,value):
    try:
        int(value)
    except:
        raise ValueError(sectionName+": incorrect data type for argument "+varName+"="+str(value)+", should be int!")
    return int(value)
def VarTypeString(sectionName,varName,value):
    if not isinstance(value, basestring):
        raise ValueError(sectionName+": incorrect data type for argument "+varName+"="+str(value)+", should be string!")
    return value
def VarTypeStringConvNone(sectionName,varName,value):
    if value==None:
        return ''
    return VarTypeString(sectionName,varName,value)
def VarTypeTupleSimple(sectionName,varName,value):
    if not isinstance(value, (list, tuple)):
        raise ValueError(sectionName+": incorrect data type for argument "+varName+"="+str(value)+", should be tuple!")
    if isinstance(value, list):
        return tuple(value)
    return value

class VarType:
    def __init__(self,dataconv,conditions=[],CanBeNone=False):
        self.dataconv=dataconv
        self.CanBeNone=CanBeNone
        self.conditions=conditions
    def chk(self,sectionName,varName,value):
        if self.CanBeNone:
            if value is None:
                return None
        value_out=None
        try:
            value_out=self.dataconv(value)
        except:
            raise ValueError(sectionName+": incorrect data type for argument "+varName+"="+str(value)+\
                ", should be "+str(self.dataconv)+(" or None!" if self.CanBeNone else "!"))
        for c in self.conditions:
            if not c['f'](value_out):
                raise ValueError(sectionName+": incorrect value of argument "+varName+"="+str(value)+"! "+c['msg'])
        return value_out

class VarTypeTuple:
    def __init__(self,dataconv,N=None,CanBeNone=False):
        self.N=N
        self.dataconv=dataconv
        self.CanBeNone=CanBeNone
    def chk(self,sectionName,varName,value):
        if self.CanBeNone:
            if value is None:
                return None
        if not isinstance(value, (tuple,list)):
            raise ValueError(sectionName+": incorrect data type for argument "+varName+"="+str(value)+", should be tuple of "+str(self.dataconv)+"!")
        if self.N!=None:
            if len(value)!=self.N:
                raise ValueError(sectionName+": incorrect tuple/list length for argument "+varName+"="+str(value)+" ! "+\
                    "Should be tuple of "+str(self.dataconv)+" of length "+str(self.N)+"!")
        value_out=[]
        for v in value:
            try:
                value_out.append(self.dataconv(v))
            except:
                raise ValueError(sectionName+": incorrect data type for argument "+varName+"="+str(value)+\
                    ", should be tuple of "+str(self.dataconv)+" of length "+str(self.N)+"!")
        return tuple(value_out)

class VarTypeDictConvert:
    def __init__(self,dictConv):
        self.dictConv=dictConv
    def chk(self,sectionName,varName,value):
        if value not in self.dictConv.keys():
            raise ValueError(sectionName+": argument "+varName+" (="+str(value)+") has incorrect value!"+\
                    "Possible values: "+",".join(self.dictConv.keys()))
        return self.dictConv[value]
    
class TiArgCheckerAndSetter(object):
    def __init__(self,sectionName="",levelZeroParameters={},defaultParameters={},switchArguments={},optionalParameters=[]):
        self.sectionName=sectionName
        self.levelZeroParameters=levelZeroParameters
        self.defaultParameters=defaultParameters
        self.switchArguments=switchArguments
        self.optionalParameters=optionalParameters
        
        
    def __process(self,arg_in,arg_out,all_param,check_extra_args=False, validate=False):
        sw_args=self.switchArguments
        
        all_param.update(self.levelZeroParameters)
        
        defaultParameters=self.defaultParameters.keys()
        for arg in defaultParameters:
            if arg not in sw_args.keys():
                arg_out[arg]=self.defaultParameters[arg]
        
        #process level 0
        #process values with default values
        for arg in defaultParameters:
            if arg in arg_in and arg not in sw_args.keys():
                arg_out[arg]=arg_in.pop(arg)
                
        #process values without default values
        for arg in self.levelZeroParameters.keys():
            if arg in arg_in and arg not in defaultParameters:
                arg_out[arg]=arg_in.pop(arg)
                
        #process level 1, induced by switches
        for sw_arg in sw_args.keys():
            if sw_arg in arg_in:
                sw_arg_val=arg_in.pop(sw_arg)
                arg_out[sw_arg]=sw_arg_val
            else:
                if sw_arg in defaultParameters:
                    sw_arg_val=self.defaultParameters[sw_arg]
                    arg_out[sw_arg]=sw_arg_val
                else:
                    raise ValueError(self.sectionName+": argument "+sw_arg+" is not specified!")
            if sw_arg_val not in sw_args[sw_arg]['switches'].keys():
                raise ValueError(self.sectionName+": argument "+sw_arg+" (="+str(sw_arg_val)+") has incorrect value!"+\
                    "Possible values: "+",".join(sw_args[sw_arg]['switches'].keys()))
            all_param[sw_arg]={'validator':sw_args[sw_arg]['validator'],'desc':sw_args[sw_arg]['desc']}
            all_param.update(sw_args[sw_arg]['switches'][sw_arg_val].levelZeroParameters)
            
            chkr=sw_args[sw_arg]['switches'][sw_arg_val]
            if chkr.sectionName=="":
                chkr.sectionName=self.sectionName+':'+sw_arg+'='+str(sw_arg_val)
            chkr.__process(arg_in,arg_out,all_param)
            
        if check_extra_args:
            if len(arg_in)==1:
                raise ValueError(self.sectionName+": unknown argument: "+arg_in.keys()[0])
            if len(arg_in)>1:
                raise ValueError(self.sectionName+": unknown arguments: "+", ".join(arg_in.keys()))
        
        if validate:
            errcount=0
            errmsg=''
            for param in all_param:
                if param in self.optionalParameters and param not in arg_out:
                    continue
                if param not in arg_out:
                    errcount+=1
                    errmsg+=self.sectionName+": argument "+param+" is not set!\n"
                else:
                    arg_out[param]=all_param[param]['validator'](self.sectionName,param,arg_out[param])
                    
            if errcount>0:
                raise ValueError(errmsg)
            
        
        
    def process(self,kwarg):
        arg_in=copy.deepcopy(kwarg)
        arg_out={}
        all_param={}
        
        
        self.__process(arg_in,arg_out,all_param, True, True)
        
        return arg_out
    def addSwitchArguments(self,name,ArgChecker,):
        pass
    
class VarTypeMaskOr:
    def __init__(self,dictConv):
        self.dictConv=dictConv
    def chk(self,sectionName,varName,value):
        
        if not isinstance(value, (tuple,list)):
            value_in=(value,)
        else:
            value_in=value
        value_out=0
        for v in value_in:
            if v not in self.dictConv.keys():
                raise ValueError(sectionName+": argument "+varName+" (="+str(value)+") has incorrect value!"+\
                    "Possible values: "+",".join(self.dictConv.keys())+" or tuple/list of these values!")
            value_out=value_out|self.dictConv[v]
        return value_out
# TitanSimulationBase - process the user input, validate it and prepare
# to be dijested by cxxTitanSimulationBase
#


# In TitanSimulationBase validation is happining in two stages
# first each setter validates its own input
# second __validate method performe cross setter validation
#
# variable convention
# ui_<section> user input for <section>
# default_<section> default values for <section> if <section> section uses 
# abitrary keyword arguments
# possible_<variable/sectoon> passible values
#
class TitanSimulationBase(object):
    """ Base class for TitanSimulation
    it defines setters, validate user input
    """
    possible_vizoutputs={
        'tecplotxxxx.tec':1, # first bit flag
        'mshplotxxxx.tec':2, # second bit flag
        'XDMF/Paraview':4, # third bit flag
        'grass_sites':8, # fourth bit flag
        'tecplot':1, # first bit flag
        'meshlot':2, # second bit flag
        'xdmf':4, # third bit flag
        'grasssites':8 # fourth bit flag
    }
    possible_gis_formats={
        'GIS_GRASS':MapNames.GIS_GRASS,
        'GDAL':MapNames.GDAL
    }
    possible_orders={'First':1,'Second':2}
    
    possible_pile_types={
         'Paraboloid':PileProps.PARABALOID,
         'Cylinder':PileProps.CYLINDER,
         #'PLANE':PileProps.PLANE,
         #'CASITA':PileProps.CASITA,
         #'POPO':PileProps.POPO,
         #'ID1':PileProps.ID1,
         #'ID2':PileProps.ID2
    }
    #defaultParameters can be removed
    possible_internal_mat_models={
        'Coulomb':{
            'allParameters':('order','int_frict'),
            'defaultParameters':{'order':'First','int_frict':37.0},
            'elementType':ElementType_SinglePhase,
            'integrators':[
                {
                    'conditions' :[lambda tsim,modprop,numprop: numprop['order']==1 or numprop['order']==2],
                    'constructor':Integrator_SinglePhase_Coulomb
                }
            ]
        },
        'Voellmy':{
            'allParameters':('order','mu','xi','int_frict'), 
            'defaultParameters':{
                'order':'First',
                'mu' : 0.5,
                'xi' : 120.0,
                'int_frict':37.0
            },
            'elementType':ElementType_SinglePhase,
            'integrators':[{
                    'conditions' :[lambda tsim,modprop,numprop: numprop['order']==1],
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
                'I_O':0.3,
                'int_frict':37.0
            },
            'elementType':ElementType_SinglePhase,
            'integrators':[{
                    'conditions' :[lambda tsim,modprop,numprop: numprop['order']==1],
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
                'I_not':0.3,
                'int_frict':37.0
            },
            'elementType':ElementType_SinglePhase,
            'integrators':[{
                    'conditions' :[lambda tsim,modprop,numprop: numprop['order']==1],
                    'constructor':Integrator_SinglePhase_Maeno_FirstOrder
            }],
        },
        'TwoPhases_Coulomb':{
            'allParameters':('order','int_frict',), 
            'defaultParameters':{'order':'First','int_frict':37.0},
            'elementType':ElementType_TwoPhases,
            'integrators':[{
                    'conditions' :[lambda tsim,modprop,numprop: numprop['order']==1],
                    'constructor':Integrator_TwoPhases_Coulomb
            }]
        }
    }
    default_GIS={
                 'gis_format':possible_gis_formats['GIS_GRASS'],
                 'gis_main':None,
                 'gis_sub':None,
                 'gis_mapset':None,
                 'gis_map':None,
                 'gis_vector':None,
                 'region_limits':None
                 }
    def __init__(self,overwrite_output=False):
        #initiate all class members
        self.overwrite_output=VarType(bool).chk("TitanSimulation", "overwrite_output", overwrite_output)
        #GIS
        self.ui_GIS=None
        self.chk_GIS=TiArgCheckerAndSetter(
            sectionName="setGIS",
            levelZeroParameters={
                'region_limits':{'desc':'',
                    'validator':VarTypeTuple(float,N=4,CanBeNone=True).chk
                }
            },
            defaultParameters={'region_limits':None},
            switchArguments={'gis_format':
                {
                    'desc':'',
                    'validator':VarTypeDictConvert(TitanSimulationBase.possible_gis_formats).chk,
                    'switches':{
                        'GIS_GRASS':TiArgCheckerAndSetter(
                            levelZeroParameters={
                                'gis_main':{'validator':VarTypeString,'desc':''},
                                'gis_sub':{'validator':VarTypeString,'desc':''},
                                'gis_mapset':{'validator':VarTypeString,'desc':''},
                                'gis_map':{'validator':VarTypeString,'desc':''},
                                'gis_vector':{'validator':VarTypeStringConvNone,'desc':''}
                            },
                            defaultParameters={'gis_vector':None}
                        ),
                        'GDAL':TiArgCheckerAndSetter(
                            levelZeroParameters={
                                'gis_map':{'validator':VarTypeString,'desc':''}
                            },
                            defaultParameters={'gis_vector':None}
                        )
                    }
                }
            }
        )
        #scale
        self.ui_Scale=None
        self.chk_Scale=TiArgCheckerAndSetter(
            sectionName="setScale",
            levelZeroParameters={
                'length_scale':{'desc':'',
                    'validator':VarType(float,conditions=[{'f':lambda v: v > 0.0,'msg':'should be positive!'}]).chk
                },
                'gravity_scale':{'desc':'',
                    'validator':VarType(float,conditions=[{'f':lambda v: v > 0.0,'msg':'should be positive!'}]).chk
                },
                'height_scale':{'desc':'',
                    'validator':VarType(float,conditions=[{'f':lambda v: v > 0.0,'msg':'should be positive or None!'}],CanBeNone=True).chk
                },
            },
            defaultParameters={'height_scale':None}
        )
        #NumProp
        self.ui_NumProp=None
        self.chk_NumProp=TiArgCheckerAndSetter(
            sectionName="setNumProp",
            levelZeroParameters={
                'AMR':{'desc':'',
                    'validator':VarType(bool).chk
                },
                'number_of_cells_across_axis':{'desc':'',
                    'validator':VarType(int,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk
                },
                'order':{'desc':'',
                    'validator':VarTypeDictConvert(TitanSimulationBase.possible_orders).chk
                },
                'geoflow_tiny':{'desc':'',
                    'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk
                },
                'short_speed':{'desc':'',
                    'validator':VarType(bool).chk
                },    
            },
            defaultParameters={'short_speed':False,'geoflow_tiny':0.0001}
        )
        #setMatModel
        self.ui_MatModel=None
        self.chk_MatModel=TiArgCheckerAndSetter(
            sectionName="setMatModel",
            levelZeroParameters={                                      
            },
            defaultParameters={'use_gis_matmap':False},
            switchArguments={
                'use_gis_matmap':
                {
                    'desc':'',
                    'validator':VarType(bool).chk,
                    'switches':{
                        False:TiArgCheckerAndSetter(
                            levelZeroParameters={
                                'bed_frict':{'desc':'',
                                    'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk
                                }
                            }
                        ),
                        True:TiArgCheckerAndSetter(
                            levelZeroParameters={
                                'bed_frict':{'desc':'should be tuple of tuple with name and bed_frict',
                                    #@todo better validator for bed_frict in case of use_gis_matmap=True
                                    'validator':VarTypeTupleSimple
                                }
                            }
                        ),
                    }
                },
                'model':
                {
                    'desc':'',
                    'validator':VarTypeString,
                    'switches':{
                        'Coulomb':TiArgCheckerAndSetter(
                            levelZeroParameters={
                                'int_frict':{'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk,'desc':''},
                            }
                        ),
                        'Voellmy':TiArgCheckerAndSetter(
                            levelZeroParameters={
                                'mu':{'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk,'desc':''},
                                'xi':{'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk,'desc':''},
                                'int_frict':{'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk,'desc':''},
                            }
                        ),
                        'Pouliquen':TiArgCheckerAndSetter(
                            levelZeroParameters={
                                'phi1':{'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk,'desc':''},
                                'phi2':{'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk,'desc':''},
                                'partdiam':{'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk,'desc':''},
                                'I_O':{'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk,'desc':''},
                                'int_frict':{'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk,'desc':''}
                            }
                        ),
                        'Maeno':TiArgCheckerAndSetter(
                            levelZeroParameters={
                                'phis':{'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk,'desc':''},
                                'phi2':{'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk,'desc':''},
                                'partdiam':{'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk,'desc':''},
                                'I_not':{'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk,'desc':''},
                                'int_frict':{'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk,'desc':''}
                            }
                        ),
                        'TwoPhases_Coulomb':TiArgCheckerAndSetter(
                            levelZeroParameters={
                                'int_frict':{'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk,'desc':''},
                            }
                        ),
                    }
                }
            }
        )
        #setTimeProps
        self.ui_TimeProps=None
        self.chk_TimeProps=TiArgCheckerAndSetter(
            sectionName="setTimeProps",
            levelZeroParameters={
                'max_iter':{'desc':'',
                    'validator':VarType(int,conditions=[{'f':lambda v: v > 0,'msg':'should be positive or None!'}],CanBeNone=True).chk
                },
                'max_time':{'desc':'',
                    'validator':VarType(float,conditions=[{'f':lambda v: v > 0.0,'msg':'should be positive or None!'}],CanBeNone=True).chk
                }
            },  
            defaultParameters={'max_iter':None,'max_time':None},  
        )
        self.setTimeProps()
        #setRestartOutput
        self.ui_RestartOutput=None
        self.chk_RestartOutput=TiArgCheckerAndSetter(
            sectionName="setRestartOutput",
            levelZeroParameters={
                'diter':{'desc':'',
                    'validator':VarType(int,conditions=[{'f':lambda v: v > 0,'msg':'should be positive or None!'}],CanBeNone=True).chk
                },
                'dtime':{'desc':'',
                    'validator':VarType(float,conditions=[{'f':lambda v: v > 0.0,'msg':'should be positive or None!'}],CanBeNone=True).chk
                },
                'keep_all':{'desc':'',
                    'validator':VarType(bool).chk
                },
                'keep_redundant_data':{'desc':'',
                    'validator':VarType(bool).chk
                },
                'output_prefix':{'desc':'',
                    'validator':VarTypeString
                },             
            },
            defaultParameters={'dtime':None, 'diter':1000, 'keep_all':False, 'keep_redundant_data':False,'output_prefix':'restart'}
        )
        self.setRestartOutput()
        #setTimeSeriesOutput
        self.ui_TimeSeriesOutput=None
        self.chk_TimeSeriesOutput=TiArgCheckerAndSetter(
            sectionName="setTimeSeriesOutput",
            levelZeroParameters={
                'vizoutput':{'desc':'',
                    'validator':VarTypeMaskOr(TitanSimulationBase.possible_vizoutputs).chk
                },
                'diter':{'desc':'',
                    'validator':VarType(int,conditions=[{'f':lambda v: v > 0,'msg':'should be positive or None!'}],CanBeNone=True).chk
                },
                'dtime':{'desc':'',
                    'validator':VarType(float,conditions=[{'f':lambda v: v > 0.0,'msg':'should be positive or None!'}],CanBeNone=True).chk
                },
                'output_prefix':{'desc':'',
                    'validator':VarTypeString
                },             
            },
            defaultParameters={'dtime':None, 'diter':1000, 'keep_all':False, 'keep_redundant_data':False,'output_prefix':'vizout'}
        )
        #setStatProps
        self.ui_StatProps=None
        self.chk_StatProps=TiArgCheckerAndSetter(
            sectionName="setStatProps",
            levelZeroParameters={
                'enabled':{'desc':'',
                    'validator':VarType(bool).chk
                },
                'edge_height':{'desc':'',
                    'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}],CanBeNone=True).chk
                },
                'test_height':{'desc':'',
                    'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}],CanBeNone=True).chk
                },
                'test_location':{'desc':'',
                    'validator':VarTypeTuple(float,N=2,CanBeNone=True).chk
                },
            },
            defaultParameters={'enabled':True,'edge_height':None, 'test_height':None, 'test_location':None}
        )
        self.setStatProps(True)
        #setOutlineProps
        self.ui_OutlineProps=None
        self.chk_OutlineProps=TiArgCheckerAndSetter(
            sectionName="setOutlineProps",
            levelZeroParameters={
                'enabled':{'desc':'',
                    'validator':VarType(bool).chk
                },
                'max_linear_size':{'desc':'',
                    'validator':VarType(int,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk
                },
                'use_DEM_resolution':{'desc':'',
                    'validator':VarType(bool).chk
                },
            },
            defaultParameters={'enabled':True,'max_linear_size':1024, 'use_DEM_resolution':False}
        )
        self.setOutlineProps(True)
        #addPile
        self.ui_Pile=[]
        self.chk_Pile=TiArgCheckerAndSetter(
            sectionName="addPile",
            levelZeroParameters={
                'height':{'desc':'',
                    'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk
                },                
                'center':{'desc':'',
                    'validator':VarTypeTuple(float,N=2).chk
                },
                'radii':{'desc':'',
                    'validator':VarTypeTuple(float,N=2).chk
                },
                'orientation':{'desc':'',
                    'validator':VarType(float).chk
                }, 
                'Vmagnitude':{'desc':'',
                    'validator':VarType(float).chk
                }, 
                'Vdirection':{'desc':'',
                    'validator':VarType(float).chk
                },
                'pile_type':{'desc':'',
                    'validator':VarTypeDictConvert(TitanSimulationBase.possible_pile_types).chk
                },               
                'vol_fract':{'desc':'',
                    'validator':VarType(float,conditions=[
                        {'f':lambda v: v >= 0.0,'msg':'should be >=0.0 and <=1.0!'},
                        {'f':lambda v: v <=1.0,'msg':'should be >=0.0 and <=1.0!'}]).chk
                },           
                
            },
            defaultParameters={'orientation':0.0,'Vmagnitude':0.0,'Vdirection':0.0,'pile_type':'Cylinder',},
            optionalParameters=['vol_fract']
        )
        #addFluxSource
        self.ui_FluxSource=[]
        self.chk_FluxSource=TiArgCheckerAndSetter(
            sectionName="addFluxSource",
            levelZeroParameters={
                'influx':{'desc':'',
                    'validator':VarType(float,conditions=[{'f':lambda v: v > 0.0,'msg':'should be positive!'}]).chk
                },  
                'start_time':{'desc':'',
                    'validator':VarType(float,conditions=[{'f':lambda v: v >= 0.0,'msg':'should be >=0.0!'}]).chk
                },  
                'end_time':{'desc':'',
                    'validator':VarType(float,conditions=[{'f':lambda v: v >= 0.0,'msg':'should be >=0.0!'}]).chk
                },                            
                'center':{'desc':'',
                    'validator':VarTypeTuple(float,N=2).chk
                },
                'radii':{'desc':'',
                    'validator':VarTypeTuple(float,N=2).chk
                },
                'orientation':{'desc':'',
                    'validator':VarType(float).chk
                }, 
                'Vmagnitude':{'desc':'',
                    'validator':VarType(float).chk
                }, 
                'Vdirection':{'desc':'',
                    'validator':VarType(float).chk
                },     
                
            },
            defaultParameters={'orientation':0.0,'Vmagnitude':0.0,'Vdirection':0.0,}
        )
        #addDischargePlane
        self.ui_DischargePlane=[]
        self.chk_DischargePlane=TiArgCheckerAndSetter(
            sectionName="addDischargePlane",
            levelZeroParameters={
                'x_a':{'desc':'',
                    'validator':VarType(float).chk
                }, 
                'y_a':{'desc':'',
                    'validator':VarType(float).chk
                },
                'x_b':{'desc':'',
                    'validator':VarType(float).chk
                }, 
                'y_b':{'desc':'',
                    'validator':VarType(float).chk
                },
            }
        )
        #loadRestart
        self.ui_loadRestart=None
        self.chk_loadRestart=TiArgCheckerAndSetter(
            sectionName="loadRestart",
            levelZeroParameters={
                'filename':{'desc':'',
                    'validator':VarTypeString
                },
            }
        )
        
        
    def setGIS(self,**kwarg):
        """Set where and how to get GIS
        gis_format - format of GIS
            possible values:'GIS_GRASS' or 'GDAL'
        Parameters for GIS_GRASS format:
        gis_main=string,
        gis_sub=string,
        gis_mapset=string,
        gis_map=string,
        Parameters for GDAL format:
        gis_map=string - location of map on filesystem
        
        gis_vector=None or string,
        General Parameters:
        region_limits - set region for simulation 
            possible values:
                None - use whole map <default value>
                list/tuple if four floats: (min_x,min_y,max_x,max_y) - set simulation region to specified values
            default:None
        """
        #process user input
        ui_GIS =self.chk_GIS.process(kwarg)
        
        #validate map names
        # if there is no topo file, quit
        if ui_GIS['gis_format'] == MapNames.GIS_GRASS:
            errmsg='Missing GIS information.  No job will be run.'
            if ui_GIS['gis_main'] == '' or ui_GIS['gis_sub'] == '' or ui_GIS['gis_mapset'] == '' or ui_GIS['gis_map'] == '':
                raise ValueError(errmsg)
            
            p=""
            for dirname in (ui_GIS['gis_main'],ui_GIS['gis_sub'],ui_GIS['gis_mapset']):
                p=os.path.join(p,dirname)
                if not os.path.isdir(p):
                    raise ValueError(errmsg+". "+p+" does not exist!")
            #gis_map?
        elif ui_GIS['gis_format'] == MapNames.GDAL:
            if ui_GIS['gis_map'] == '':
                raise ValueError(errmsg)
        #@todo validate gis_vector
        
        #validate region_limits
        
        self.ui_GIS=ui_GIS
    
    def setScale(self,**kwarg):
        self.ui_Scale=self.chk_Scale.process(kwarg)
    def setNumProp(self,**kwarg):
        self.ui_NumProp=self.chk_NumProp.process(kwarg)
    def setMatModel(self,**kwarg):
        self.ui_MatModel=self.chk_MatModel.process(kwarg)
    def setTimeProps(self,**kwarg):
        self.ui_TimeProps=self.chk_TimeProps.process(kwarg)
        if self.ui_TimeProps['max_iter']==None:
            self.ui_TimeProps['max_iter']=-1
        if self.ui_TimeProps['max_time']==None:
            self.ui_TimeProps['max_time']=-1.0
    def setRestartOutput(self,**kwarg):
        self.ui_RestartOutput=self.chk_RestartOutput.process(kwarg)
        if self.ui_RestartOutput['dtime']==None:
            self.ui_RestartOutput['dtime']=-1.0
        if self.ui_RestartOutput['diter']==None:
            self.ui_RestartOutput['diter']=-1
    def setTimeSeriesOutput(self,**kwarg):
        self.ui_TimeSeriesOutput=self.chk_TimeSeriesOutput.process(kwarg)
        if self.ui_TimeSeriesOutput['dtime']==None:
            self.ui_TimeSeriesOutput['dtime']=-1.0
        if self.ui_TimeSeriesOutput['diter']==None:
            self.ui_TimeSeriesOutput['diter']=-1
    def setStatProps(self,enabled,**kwarg):
        ui=self.chk_StatProps.process(kwarg)
        
        #Test if flow reaches height [m] ...
        if ui['edge_height'] == None:
            ui['edge_height'] = -1.0
        else:
            ui['edge_height'] = float(ui['edge_height'])
            if ui['edge_height'] <= 0:
                raise ValueError('TitanSimulation::edge_height should be positive or None\n')
        
        #Height used to define flow outline (>0) [m]
        if ui['test_height'] == None:
            ui['test_height'] = -2.0
            if ui['test_location'] != None:
                raise ValueError('TitanSimulation::test_height set to None, test_location also should be None!')
            ui['test_location'] =(0.0, 0.0)
        else:
            if ui['test_height'] <= 0:
                raise ValueError('TitanSimulation::test_height should be positive or None\n')
            #... at test point (x and y location)
            if ui['test_location']==None:
                raise ValueError('TitanSimulation::test_location should be set if test_height>0\n') 
        self.ui_StatProps=ui
    def setOutlineProps(self,enabled,**kwarg):
        self.ui_OutlineProps=self.chk_OutlineProps.process(kwarg)
    def addPile(self,**kwarg):
        self.ui_Pile.append(self.chk_Pile.process(kwarg))
    def addFluxSource(self,**kwarg):
        ui=self.chk_FluxSource.process(kwarg)
        if ui['start_time'] >= ui['end_time']:
            raise ValueError('TitanSimulation:addFluxSource:start_time should be less than end_time !')
        
        self.ui_FluxSource.append(ui)
        
    def addDischargePlane(self,x_a,y_a,x_b,y_b):
        args={'x_a':x_a,'y_a':y_a,'x_b':x_b,'y_b':y_b}
        self.ui_DischargePlane.append(self.chk_DischargePlane.process(args))
        
    def loadRestart(self,**kwarg):
        self.ui_loadRestart=self.chk_loadRestart.process(kwarg)
        
    def _validate(self):
        #######################################################################
        #check cross section
        if self.ui_MatModel['model']=='TwoPhases_Coulomb':
            for pile in self.ui_Pile:
                if 'vol_fract' not in pile:
                    raise ValueError('TitanSimulation:addFluxSource: TwoPhases_Coulomb is set addPile should set vol_fract!')
            if len(self.ui_FluxSource)>0:
                raise NotImplementedError("TitanSimulation:addFluxSource: FluxSources are not imlemented for TwoPhases_Coulomb model!")
            if len(self.ui_DischargePlane)>0:
                raise NotImplementedError("TitanSimulation:addDischargePlane: DischargePlanes are not imlemented for TwoPhases_Coulomb model!")
        else:
            for pile in self.ui_Pile:
                if 'vol_fract' in pile:
                    raise ValueError('TitanSimulation:addFluxSource: Single phase model is set, addPile should not have vol_fract argument!')
        
        #######################################################################
        #check satisfaction of integrator
        integrator=None
        model=self.ui_MatModel['model']
        integrators=TitanSimulationBase.possible_internal_mat_models[model]['integrators']
        msg=""
        
        for desc in integrators:
            does_feat=True
            msg+="Testing: "+desc['constructor'].__name__+"\n"
            for cond in desc['conditions']:
                cond_result=cond(self,self.ui_MatModel,self.ui_NumProp)
                msg+="  "+inspect.getsource(cond).strip()+" ===>"+str(cond_result)+"\n"
                does_feat=does_feat and cond_result
            if does_feat:
                integrator=desc['constructor']
                break
        if integrator==None:
            raise ValueError("Can not find suitable integrator, here is the hint\n"+msg+"\ncheck the manual.")
        self.integratorConstructor=integrator
        
        
        
class TitanSimulation(TitanSimulationBase):
    def __init__(self,overwrite_output=False):
        super(TitanSimulation, self).__init__(overwrite_output=overwrite_output)
        #initiate all class members
        self.sim=None
        
    def _setCxxTitanSimulation(self):
        """initiate and set all parameters of cxxTitanSimulation"""
        if 0:
            print self.ui_GIS
            print self.ui_Scale
            print self.ui_NumProp
            print self.ui_MatModel
            print self.ui_TimeProps
            print 'ui_RestartOutput',self.ui_RestartOutput
            print 'ui_TimeSeriesOutput',self.ui_TimeSeriesOutput
            print self.ui_StatProps
            print self.ui_OutlineProps
            print self.ui_Pile
            print self.ui_FluxSource
            print self.ui_DischargePlane
            print self.ui_loadRestart
        #convinience references
        ui_GIS=self.ui_GIS
        ui_Scale=self.ui_Scale
        ui_NumProp=self.ui_NumProp
        ui_MatModel=self.ui_MatModel
        ui_TimeProps=self.ui_TimeProps
        ui_RestartOutput=self.ui_RestartOutput
        ui_TimeSeriesOutput=self.ui_TimeSeriesOutput
        ui_StatProps=self.ui_StatProps
        ui_OutlineProps=self.ui_OutlineProps
        ui_Pile=self.ui_Pile
        ui_FluxSource=self.ui_FluxSource
        ui_DischargePlane=self.ui_DischargePlane
        ui_loadRestart=self.ui_loadRestart
        
        model=self.ui_MatModel['model']
        
        self.sim=cxxTitanSimulation()
        self.sim.overwrite_output=self.overwrite_output
        #######################################################################
        #check the presence of output files and delete them if nessesary also create directory for some outputs
        if self.sim.myid==0:
            def check_and_remove_filedir(filename):
                if os.path.exists(filename):
                    if self.overwrite_output:
                        if os.path.isfile(filename):
                            os.remove(filename)
                        if os.path.isdir(filename):
                            shutil.rmtree(filename)
                    else:
                        raise IOError("Output file or directory exists ("+filename+"). Remove it manually or set overwrite_output to True or change output prefix.")
            def check_and_remove_filedir_by_wildcard(filename):
                files=glob.glob(filename)
                if not self.overwrite_output and len(files)>0:
                    raise IOError("Output files or directories exists ("+",".join(files)+"). Remove it manually or set overwrite_output to True or change output prefix.")
                for f in files:
                    check_and_remove_filedir(f)
            #restarts     
            output_prefix=self.ui_RestartOutput['output_prefix']
            check_and_remove_filedir(output_prefix)
            check_and_remove_filedir_by_wildcard("%s_Quad[49]_p[0-9][0-9][0-9][0-9].xmf"%(output_prefix,));
            
            os.mkdir(output_prefix)
            #visoutput
            output_prefix=self.ui_TimeSeriesOutput['output_prefix']
            check_and_remove_filedir(output_prefix)
            check_and_remove_filedir_by_wildcard("%s_xdmf_p[0-9][0-9][0-9][0-9].xmf"%(output_prefix,));
            
            os.mkdir(output_prefix)
        
        #######################################################################
        #Set GIS
        self.sim.get_mapnames().set(
            ui_GIS['gis_format'],ui_GIS['gis_main'], ui_GIS['gis_sub'], 
            ui_GIS['gis_mapset'],ui_GIS['gis_map'],  ui_GIS['gis_vector'], 0
        )
        
        if ui_GIS['region_limits']!=None:
            #Minimum x and y location (UTM E, UTM N)
            #Maximum x and y location (UTM E, UTM N)
            self.sim.get_mapnames().set_region_limits(ui_GIS['region_limits'][0],ui_GIS['region_limits'][2],ui_GIS['region_limits'][1],ui_GIS['region_limits'][3])
        
        #######################################################################
        #Set Scale
        #Length Scale [m]
        self.sim.scale_.length = ui_Scale['length_scale']
        #gravity scaling factor [m/s^2]
        self.sim.scale_.gravity = ui_Scale['gravity_scale']
        #height scaling factor
        if ui_Scale['height_scale']==None:
            self.sim.scale_.height=0.0
            self.sim.scale_.auto_calc_height_scale=True
        else:
            self.sim.scale_.height = ui_Scale['height_scale']
            self.sim.scale_.auto_calc_height_scale=False
        
        #######################################################################
        # NumProp and MatModel
        #set element type
        elementType=TitanSimulationBase.possible_internal_mat_models[model]['elementType']
        self.sim.set_element_type(elementType)
        
        self.sim.adapt=int(ui_NumProp['AMR'])
        self.sim.set_short_speed(ui_NumProp['short_speed'])
        #geoflow_tiny
        
        
        ##################
        # build integrator
        
        # get parameters
        ui_NumModmProp=copy.deepcopy(ui_NumProp)
        ui_NumModmProp.update(ui_MatModel)
        
        ui_Integrator={}
        integrator_all_parameters=TitanSimulationBase.possible_internal_mat_models[model]['allParameters']
        for param in integrator_all_parameters:
            ui_Integrator[param]=ui_NumModmProp[param]
     
        #construct integrator
        integrator_obj=self.integratorConstructor(self.sim)
        #uncouple the c++ from python proxy 
        integrator_obj.thisown=0
        #set parameters
        for k,v in ui_Integrator.iteritems():
            setattr(integrator_obj,k,v)
        
        self.sim.set_integrator(integrator_obj)
        self.integrator_initialized=True
        
        ##################
        #
        
        ##################
        # MatMap
        #Use GIS Material Map?        
        self.sim.use_gis_matmap = ui_MatModel['use_gis_matmap']
        
        #here we need to set mat prop
        if self.sim.get_element_type()==ElementType_SinglePhase:
            matprops=MatProps(self.sim.scale_)
        elif self.sim.get_element_type()==ElementType_TwoPhases:
            matprops=MatPropsTwoPhases(self.sim.scale_)
            
        matprops.number_of_cells_across_axis = int(ui_NumProp['number_of_cells_across_axis'])
        
        if self.sim.use_gis_matmap == False:
            matprops.material_count=1
            matprops.matnames.push_back("all materials")
            matprops.bedfrict.push_back(float(ui_MatModel['bed_frict']))
        else:  #if they did want to use a GIS material map...
            #matprops.material_count=len(mat_names)
            #if len(bed_frict)!=len(mat_names):
            #    raise Exception("number of mat_names does not match number of bed_frict")
            #for i in range(len(bed_frict)):
            #    matprops.matnames.push_back(mat_names[i])
            #    matprops.bedfrict.push_back(float(bed_frict[i]))
            raise NotImplementedError("GIS material map Not implemented yet")
        
        self.sim.set_matprops(matprops)
        #uncouple the c++ from python proxy 
        matprops.thisown=0 
        
        #######################################################################
        # TimeProps
        self.sim.get_timeprops().setTime(ui_TimeProps['max_iter'],ui_TimeProps['max_time'])
        
        #######################################################################
        # RestartOutput
        self.sim.get_timeprops().setRestartOutput(ui_RestartOutput['diter'],ui_RestartOutput['dtime'])
        self.sim.restart_prefix=ui_RestartOutput['output_prefix']
        self.sim.restart_keep_all=ui_RestartOutput['keep_all']
        self.sim.restart_keep_redundant_data=ui_RestartOutput['keep_redundant_data']
        
        #######################################################################
        # TimeSeriesOutput
        self.sim.get_timeprops().setTimeSeriesOutput(ui_TimeSeriesOutput['diter'],ui_TimeSeriesOutput['dtime'])
        self.sim.vizoutput = ui_TimeSeriesOutput['vizoutput']
        self.sim.vizoutput_prefix=ui_TimeSeriesOutput['output_prefix'];
        
        #######################################################################
        # ui_StatProps
        # @todo implement enabled ui_StatProps
        self.sim.get_statprops().set(
            ui_StatProps['edge_height'], ui_StatProps['test_height'],
            ui_StatProps['test_location'][0], ui_StatProps['test_location'][1])
        
        #######################################################################
        # ui_OutlineProps
        # @todo implement ui_OutlineProps
        #c++ is not implemented yet
        #######################################################################
        # ui_Pile
        self.pileprops=None
        for pile in ui_Pile:
            if self.sim.get_element_type()==ElementType_SinglePhase:
                if self.pileprops==None:
                    self.pileprops=PileProps()
                    self.pileprops.thisown=0
                    self.sim.set_pileprops(self.pileprops)
                if not isinstance(self.pileprops,PileProps):
                    raise ValueError("Can not mix element type for piles!")
                if isinstance(self.pileprops,PilePropsTwoPhases):
                    raise ValueError("Can not mix element type for piles!")
                
                if pile!=None:
                    self.pileprops.addPile(pile['height'], pile['center'][0], pile['center'][1], pile['radii'][0], 
                                           pile['radii'][1], pile['orientation'], pile['Vmagnitude'], pile['Vdirection'],pile['pile_type'])
            elif self.sim.get_element_type()==ElementType_TwoPhases:
                if self.pileprops==None:
                    self.pileprops=PilePropsTwoPhases()
                    self.pileprops.thisown=0
                    self.sim.set_pileprops(self.pileprops)
                if not isinstance(self.pileprops,PileProps):
                    raise ValueError("Can not mix element type for piles!")
                if not isinstance(self.pileprops,PilePropsTwoPhases):
                    raise ValueError("Can not mix element type for piles!")
                
                if pile!=None:
                    self.pileprops.addPile(pile['height'], pile['center'][0], pile['center'][1], pile['radii'][0], 
                        pile['radii'][1], pile['orientation'], pile['Vmagnitude'], pile['Vdirection'],pile['pile_type'],pile['vol_fract'])
            else:
                raise ValueError("Unknown element type")
        #######################################################################
        # ui_FluxSource
        for fluxSource in ui_FluxSource:
            self.sim.fluxprops.addFluxSource(fluxSource['influx'],fluxSource['start_time'],fluxSource['end_time'], fluxSource['center'][0], fluxSource['center'][1], fluxSource['radii'][0], 
                fluxSource['radii'][1],fluxSource['orientation'],fluxSource['Vmagnitude'],fluxSource['Vdirection'])
        #######################################################################
        # ui_DischargePlane
        for disPlane in ui_DischargePlane:
            self.sim.discharge_planes.addDischargePlane(disPlane['x_a'],disPlane['y_a'],disPlane['x_b'],disPlane['y_b'])
        #######################################################################
        # ui_loadRestart
        if ui_loadRestart !=None:
            raise NotImplementedError("loadRestart not implemented yet!")
    def _preproc(self):
        if self.sim.myid==0:
            preproc=TitanPreproc(self.sim)
            preproc.validate();
            preproc.run();
        
    def run(self):
        """Perform simulation"""
        
        self._validate()
        #by this time all values should be sanitized
        self._setCxxTitanSimulation()
        
        if self.sim.myid==0:
            # run preproc.x to create the fem grid, if it is not already there
            #if os.access('PRE/preproc.x',os.X_OK)==0:
            #    os.system('cd PRE;gmake')
            # locate titan_preprocess: 
            # first check the local-diretory,
            
            print "\npreproc..."
            self._preproc()
    
            if self.sim.get_mapnames().gis_vector != "":
                print "\nvectordatpreproc..."
                vectordatpreproc(self.sim.get_mapnames().gis_main, self.sim.get_mapnames().gis_sub, self.sim.get_mapnames().gis_mapset, self.sim.get_mapnames().gis_map, self.sim.get_mapnames().gis_vector)
            print
        #self.sim.input_summary()
        self.sim.run()
        

class TitanSimulation1(TitanSimulationBase):  
    def __init__(self):
        super(TitanSimulation, self).__init__()
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
 
class TitanSimulation2(TitanSimulationBase):
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
        gis_mapset = gis_mapset if gis_mapset!=None else ''
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
