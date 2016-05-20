#
# About this simulation:
#
# Similar to Coulomb example but most parameters were set explicitly
#


sim=TitanSimulation(
    overwrite_output=True
)

sim.setGIS(
    gis_format='GIS_GRASS', 
    gis_main='../../dem',
    gis_sub='colimafinemini',
    gis_mapset='PERMANENT',
    gis_map='colima',
    gis_vector=None, #default=None
    region_limits=None, #default=None
)
    

sim.setScale(
    length_scale=4000.0,
    gravity_scale=9.8,
    height_scale=None,  #default=None
)

sim.setNumProp( #setNumericalMethods
    AMR=True,
    number_of_cells_across_axis=16, #number_of_elements_across_axis? should go either to setAdapt or piles
    order='First',
    geoflow_tiny=0.0001,
    short_speed=False,
)
sim.setMatModel(
    model='Coulomb',

    int_frict=37.0,# other model parameter
    
    use_gis_matmap=False,
    bed_frict=27.0
)
sim.addPile(
    height=30.0,
    center=[644956.0, 2157970.0],
    radii=[55.0, 55.0],
    orientation=0.0,
    Vmagnitude=0.0,
    Vdirection=0.0,
    pile_type='Cylinder', #Paraboloid or Cylinder
)
if 0: #set to 1 to try with flux source
    sim.addFluxSource(
        influx=10.0,
        start_time=0.0,
        end_time=10000.0,
        center=[644956.0, 2157970.0],
        radii=[55.0, 55.0],
        orientation=0.0,
        Vmagnitude=0.0,
        Vdirection=0.0
    )
    
if 0: #set to 1 to try with dischargePlane
    sim.addDischargePlane(637380.0, 2145800.0, 664380.0, 2169800.0)

sim.setTimeProps(
    max_iter=1000, #one of them can be none, so the default value for both is none, if both none run till killed 
    max_time=None
)
#dump restart file every dtime simulated seconds and/or diter interations
sim.setRestartOutput(
    enabled=True,
    dtime=1.0,#float or None, default=None
    diter=None,#int or None, default=None
    keep_all=True, #True - keep all restart files or False - only last default=False?
    keep_redundant_data=True,
    output_prefix='restart', #directory name where restart h5 and individual xdmf are stored, prefix for time series xdmf, default='restart'
)
#write output for visualization every dtime simulated seconds and/or diter interations
sim.setTimeSeriesOutput(
    vizoutput=('xdmf','meshplot'), #possible formats tecplot, meshlot, xdmf, grasssites
    dtime=1.0,#float or None
    diter=None,#int or None
    output_prefix='vizout' #? directory name where frames are recorded if applicable and/or prefix for single time series output, default='vizout?'
)
#if not set default settings are used

sim.setStatProps(
    runid=-1,
    edge_height=None,
    test_height=None,
    test_location=None,
    output_prefix=''
)

#if not set default settings are used
sim.setOutlineProps(
    output_prefix='',
    enabled=True,
    max_linear_size=1024,
    init_size='DEM' #default AMR
)


#start simulation
sim.run()
