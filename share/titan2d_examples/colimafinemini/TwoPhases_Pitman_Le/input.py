#
# About this simulation:
#
# NOTE what this simulation uses a different DEM
# This simulation require bigger region and to keep 
# example sizes small, the fine map was coarsened
#
# it also require much more iterations and execution time
# because pile project itself further

sim=TitanSimulation(overwrite_output=True)

sim.setGIS(
    gis_format='GIS_GRASS', 
    gis_main='../../dem',
    gis_sub='colimacoarse',
    gis_mapset='PERMANENT',
    gis_map='colima'
)

sim.setScale(
    length_scale=4000.0,
    gravity_scale=9.8,
)

sim.setNumProp(
    AMR=True,
    number_of_cells_across_axis=16,
    order='First',
)
sim.setMatModel(
    model='TwoPhases-Pitman-Le',
    int_frict=37.0,
    bed_frict=27.0
)
sim.addPile(
    pile_type='Paraboloid',
    height=30.0,
    center=[644956.0, 2157970.0],
    radii=[55.0, 55.0],
    orientation=0.0,
    Vmagnitude=0.0,
    Vdirection=0.0,
    vol_fract=0.7,
)

sim.setTimeProps(
    max_iter=18000,
)

sim.setTimeSeriesOutput(
    vizoutput='xdmf',
    dtime=1.0,
)

#if not set default settings are used
sim.setOutlineProps(
    output_prefix='',
    enabled=True,
    max_linear_size=1024,
    init_size='AMR'
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

#start simulation
sim.run()
