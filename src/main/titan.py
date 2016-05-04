
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

def VarTypeCoulombMatMapBedFrict(sectionName,varName,value):
    print sectionName,varName,value
    value_out={}
    if not isinstance(value, (dict)):
        raise ValueError(sectionName+": incorrect data type for argument "+varName+"="+str(value)+", should be dict!")
    for k,v in value.iteritems():
        value_out[k]=VarTypeFloat(sectionName, "values of dict "+varName,v)
    return value_out

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
    def __init__(self,dictConv,NoneToStringNone=False):
        self.dictConv=dictConv
        self.NoneToStringNone=NoneToStringNone
    def chk(self,sectionName,varName,value):
        value2=value
        if self.NoneToStringNone and value==None:
            value2="None"
            
        if value2 not in self.dictConv.keys():
            raise ValueError(sectionName+": argument "+varName+" (="+str(value2)+") has incorrect value!"+\
                    "Possible values: "+",".join(self.dictConv.keys()))
        return self.dictConv[value2]
    
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
    def __init__(self,dictConv,CanBeNone=False):
        self.dictConv=dictConv
        self.CanBeNone=CanBeNone
    def chk(self,sectionName,varName,value):
        if self.CanBeNone and value==None:
            return 0
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
        'meshplot':2, # second bit flag
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
    possible_Interface_Capturing_Type={
        'Heuristic':Interface_Capturing_Type_Heuristic,
        'LevelSet':Interface_Capturing_Type_LevelSet,
        'PhaseField':Interface_Capturing_Type_PhaseField
    }
    #defaultParameters can be removed
    possible_internal_mat_models={
        'Coulomb':{
            'allParameters':('order','int_frict','stopping_criteria'),
            'defaultParameters':{'order':'First','int_frict':37.0,'stopping_criteria':None},
            'elementType':ElementType_SinglePhase,
            'interface_capturing_type':Interface_Capturing_Type_Heuristic,
            'integrators':[
                {
                    'conditions' :[lambda tsim,modprop,numprop: numprop['order']==1 or numprop['order']==2],
                    'constructor':Integrator_SinglePhase_Heuristic_Coulomb
                }
            ]
        },
        'Coulomb':{
            'allParameters':('order','int_frict','stopping_criteria'),
            'defaultParameters':{'order':'First','int_frict':37.0,'stopping_criteria':None},
            'elementType':ElementType_SinglePhase,
            'interface_capturing_type':Interface_Capturing_Type_LevelSet,
            'integrators':[
                {
                    'conditions' :[lambda tsim,modprop,numprop: numprop['order']==1],
                    'constructor':Integrator_SinglePhase_LevelSet_Coulomb
                }
            ]
        },
        'Voellmy-Salm':{
            'allParameters':('order','mu','xi'), 
            'defaultParameters':{
                'order':'First',
                'mu' : 0.5,
                'xi' : 120.0,
            },
            'elementType':ElementType_SinglePhase,
            'interface_capturing_type':Interface_Capturing_Type_Heuristic,
            'integrators':[{
                    'conditions' :[lambda tsim,modprop,numprop: numprop['order']==1],
                    'constructor':Integrator_SinglePhase_Heuristic_Voellmy_Salm
            }]
        },
        'Voellmy-Salm':{
            'allParameters':('order','mu','xi'), 
            'defaultParameters':{
                'order':'First',
                'mu' : 0.5,
                'xi' : 120.0,
            },
            'elementType':ElementType_SinglePhase,
            'interface_capturing_type':Interface_Capturing_Type_LevelSet,
            'integrators':[{
                    'conditions' :[lambda tsim,modprop,numprop: numprop['order']==1],
                    'constructor':Integrator_SinglePhase_LevelSet_Voellmy_Salm
            }]
        },
        'Pouliquen-Forterre':{
            'allParameters':('order','phi1','phi2','phi3','Beta','L_material'), 
            'defaultParameters':{
                'order':'First',
                'phi1':32.9,
                'phi2':42.0,
                'phi3':33.9,
                'Beta':0.65,
                'L_material':1.0E-3
            },
            'elementType':ElementType_SinglePhase,
            'interface_capturing_type':Interface_Capturing_Type_Heuristic,
            'integrators':[{
                    'conditions' :[lambda tsim,modprop,numprop: numprop['order']==1],
                    'constructor':Integrator_SinglePhase_Heuristic_Pouliquen_Forterre
            }]
        },
        'Pouliquen-Forterre':{
            'allParameters':('order','phi1','phi2','phi3','Beta','L_material'), 
            'defaultParameters':{
                'order':'First',
                'phi1':32.9,
                'phi2':42.0,
                'phi3':33.9,
                'Beta':0.65,
                'L_material':1.0E-3
            },
            'elementType':ElementType_SinglePhase,
            'interface_capturing_type':Interface_Capturing_Type_LevelSet,
            'integrators':[{
                    'conditions' :[lambda tsim,modprop,numprop: numprop['order']==1],
                    'constructor':Integrator_SinglePhase_LevelSet_Pouliquen_Forterre
            }]
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
    possible_stopping_criteria={
        'None':0,
        'DragBased':1
     }
    possible_outline_init_size={
        'AMR':OutLine.AMR,
        'DEM':OutLine.DEM,                
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
                'interface_capturing_type':{'desc':'',
                    'validator':VarTypeDictConvert(TitanSimulationBase.possible_Interface_Capturing_Type).chk
                },
            },
            defaultParameters={'short_speed':False,'geoflow_tiny':0.0001,'interface_capturing_type':'Heuristic'}
        )
        #setMatModel
        self.ui_MatModel=None
        self.chk_MatModel=TiArgCheckerAndSetter(
            sectionName="setMatModel",
            levelZeroParameters={                                      
            },
            switchArguments={
                'model':
                {
                    'desc':'',
                    'validator':VarTypeString,
                    'switches':{
                        'Coulomb':TiArgCheckerAndSetter(
                            levelZeroParameters={
                                'int_frict':{'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk,'desc':''},
                                'stopping_criteria':{'validator':VarTypeDictConvert(TitanSimulationBase.possible_stopping_criteria,NoneToStringNone=True).chk,'desc':''},
                            },
                            defaultParameters={'use_gis_matmap':False,'stopping_criteria':None},
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
                                                    'validator':VarTypeCoulombMatMapBedFrict
                                                }
                                            }
                                        ),
                                    }
                                }
                            }
                        ),
                        'Voellmy-Salm':TiArgCheckerAndSetter(
                            levelZeroParameters={
                                'mu':{'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk,'desc':''},
                                'xi':{'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk,'desc':''},
                            }
                        ),
                        'Pouliquen-Forterre':TiArgCheckerAndSetter(
                            levelZeroParameters={
                                'phi1':{'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk,'desc':''},
                                'phi2':{'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk,'desc':''},
                                'phi3':{'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk,'desc':''},
                                'Beta':{'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk,'desc':''},
                                'L_material':{'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk,'desc':''}
                            }
                        ),
                        'TwoPhases_Coulomb':TiArgCheckerAndSetter(
                            levelZeroParameters={
                                'int_frict':{'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk,'desc':''},
                                'bed_frict':{'desc':'',
                                                    'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}]).chk
                                                }
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
                'enabled':{'desc':'',
                    'validator':VarType(bool).chk
                },
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
            defaultParameters={'enabled':True,'dtime':None, 'diter':1000, 'keep_all':False, 'keep_redundant_data':False,'output_prefix':'restart'}
        )
        self.setRestartOutput()
        #setTimeSeriesOutput
        self.ui_TimeSeriesOutput=None
        self.chk_TimeSeriesOutput=TiArgCheckerAndSetter(
            sectionName="setTimeSeriesOutput",
            levelZeroParameters={
                'vizoutput':{'desc':'',
                    'validator':VarTypeMaskOr(TitanSimulationBase.possible_vizoutputs,CanBeNone=True).chk
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
            defaultParameters={'dtime':None, 'diter':1000, 'output_prefix':'vizout'}
        )
        self.setTimeSeriesOutput(vizoutput=None)
        #setStatProps
        self.ui_StatProps=None
        self.chk_StatProps=TiArgCheckerAndSetter(
            sectionName="setStatProps",
            levelZeroParameters={
                'edge_height':{'desc':'',
                    'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}],CanBeNone=True).chk
                },
                'test_height':{'desc':'',
                    'validator':VarType(float,conditions=[{'f':lambda v: v > 0,'msg':'should be positive!'}],CanBeNone=True).chk
                },
                'test_location':{'desc':'',
                    'validator':VarTypeTuple(float,N=2,CanBeNone=True).chk
                },
                'runid':{'desc':'',
                    'validator':VarType(int).chk
                },
                'output_prefix':{'desc':'',
                    'validator':VarTypeString
                },      
            },
            defaultParameters={'edge_height':None, 'test_height':None, 'test_location':None, 'runid':-1, 'output_prefix':''}
        )
        self.setStatProps()
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
                'init_size':{'desc':'',
                    'validator':VarTypeDictConvert(TitanSimulationBase.possible_outline_init_size).chk
                },
                'output_prefix':{'desc':'',
                    'validator':VarTypeString
                },
            },
            defaultParameters={'enabled':True,'max_linear_size':1024, 'init_size':'AMR', 'output_prefix':''}
        )
        self.setOutlineProps(enabled=True)
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
    def setStatProps(self,**kwarg):
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
    def setOutlineProps(self,**kwarg):
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
            if ui_RestartOutput['enabled']:
                output_prefix=self.ui_RestartOutput['output_prefix']
                check_and_remove_filedir(output_prefix)
                check_and_remove_filedir_by_wildcard("%s_Quad[49]_p[0-9][0-9][0-9][0-9].xmf"%(output_prefix,));
            
                os.mkdir(output_prefix)
            #visoutput
            if self.ui_TimeSeriesOutput['vizoutput']!=0:
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
        self.sim.set_geoflow_tiny(ui_NumProp['geoflow_tiny'])
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
        
        self.sim.set_interface_capturing_type(ui_NumModmProp['interface_capturing_type'])
        
        ##################
        #
        
        ##################
        # MatMap
        #Use GIS Material Map?        
        self.sim.use_gis_matmap = ui_MatModel.get('use_gis_matmap',False)
        
        #here we need to set mat prop
        if self.sim.get_element_type()==ElementType_SinglePhase:
            matprops=MatProps(self.sim.scale_)
        elif self.sim.get_element_type()==ElementType_TwoPhases:
            matprops=MatPropsTwoPhases(self.sim.scale_)
            
        matprops.number_of_cells_across_axis = int(ui_NumProp['number_of_cells_across_axis'])
        
        if self.sim.use_gis_matmap == False:
            if ui_MatModel['model']=='Coulomb' or ui_MatModel['model']=='TwoPhases_Coulomb':
                matprops.material_count=1
                matprops.matnames.push_back("all materials")
                matprops.bedfrict.push_back(float(ui_MatModel['bed_frict']))
        else:
            print ui_MatModel['bed_frict']
            matprops.material_count=len(ui_MatModel['bed_frict'])
            for k,v in ui_MatModel['bed_frict'].iteritems():
                print k,v
                matprops.matnames.push_back(k)
                matprops.bedfrict.push_back(float(v))
            #if they did want to use a GIS material map...
            #matprops.material_count=len(mat_names)
            #if len(bed_frict)!=len(mat_names):
            #    raise Exception("number of mat_names does not match number of bed_frict")
            #for i in range(len(bed_frict)):
            #    matprops.matnames.push_back(mat_names[i])
            #    matprops.bedfrict.push_back(float(bed_frict[i]))
            #raise NotImplementedError("GIS material map Not implemented yet")
        
        self.sim.set_matprops(matprops)
        #uncouple the c++ from python proxy 
        matprops.thisown=0 
        
        #######################################################################
        # TimeProps
        self.sim.get_timeprops().setTime(ui_TimeProps['max_iter'],ui_TimeProps['max_time'])
        
        #######################################################################
        # RestartOutput
        self.sim.get_timeprops().setRestartOutput(ui_RestartOutput['diter'],ui_RestartOutput['dtime'])
        self.sim.restart_enabled=ui_RestartOutput['enabled']
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
        self.sim.get_statprops().runid=ui_StatProps['runid']
        self.sim.get_statprops().output_prefix=ui_StatProps['output_prefix']
        
        #######################################################################
        # ui_OutlineProps
        self.sim.get_outline().enabled=ui_OutlineProps['enabled']
        self.sim.get_outline().max_linear_size=ui_OutlineProps['max_linear_size']
        self.sim.get_outline().init_size=ui_OutlineProps['init_size']
        self.sim.get_outline().output_prefix=ui_OutlineProps['output_prefix']
        
        #######################################################################
        # ui_Pile
        self.pileprops=None
        if self.sim.get_element_type()==ElementType_SinglePhase:
            if self.pileprops==None:
                if self.sim.get_interface_capturing_type()==Interface_Capturing_Type_Heuristic:
                    self.pileprops=PileProps()
                    self.pileprops.thisown=0
                    self.sim.set_pileprops(self.pileprops)
                elif Interface_Capturing_Type_LevelSet:
                    self.pileprops=PilePropsLevelSet()
                    self.pileprops.thisown=0
                    self.sim.set_pileprops(self.pileprops)
                else:
                    raise Exception("ERROR: this interface_capturing_type is not implemented yet!")
        elif self.sim.get_element_type()==ElementType_TwoPhases:
            if self.pileprops==None:
                self.pileprops=PilePropsTwoPhases()
                self.pileprops.thisown=0
                self.sim.set_pileprops(self.pileprops)
        else:
            raise ValueError("Unknown element type")
        
        #if len(ui_Pile)==0 :
        #    raise NotImplementedError("Simulations without piles are not implemented yet, If you need it contact the developers")
        
        for pile in ui_Pile:
            if self.sim.get_element_type()==ElementType_SinglePhase:
                if not isinstance(self.pileprops,PileProps):
                    raise ValueError("Can not mix element type for piles!")
                if isinstance(self.pileprops,PilePropsTwoPhases):
                    raise ValueError("Can not mix element type for piles!")
                
                if pile!=None:
                    self.pileprops.addPile(pile['height'], pile['center'][0], pile['center'][1], pile['radii'][0], 
                                           pile['radii'][1], pile['orientation'], pile['Vmagnitude'], pile['Vdirection'],pile['pile_type'])
            elif self.sim.get_element_type()==ElementType_TwoPhases:
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
        self.sim.run(False)
