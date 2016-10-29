sim=TitanSimulation(overwrite_output=True)

sim.setGIS(
    gis_format='GIS_GRASS', 
    gis_main='../../dem',
    gis_sub='inclined',
    gis_mapset='PERMANENT',
    gis_map='inclined',
)

sim.setScale(
    length_scale=1.0,
    gravity_scale=9.8,
)

sim.setNumProp(
    AMR=True,
    number_of_cells_across_axis=12,
    order='First',
)
sim.setMatModel(
    model='Coulomb',
    int_frict=37.0,
    
    use_gis_matmap=True,
    bed_frict={
        "Mat1":27.0,
        "Mat2":35.0,
        "Mat3":20.0,
    }
)
sim.addPile(
    pile_type='Cylinder',
    height=0.01,
    center=[0.2, 0.2],
    radii=[.02, 0.02],
    orientation=0.0,
    Vmagnitude=0.0,
    Vdirection=0.0,
)

sim.setTimeProps(
    max_iter=5000,
)

sim.setTimeSeriesOutput(
    vizoutput=('xdmf'),
    dtime=1.0,
)

sim.setRestartOutput(
    enabled=True,
    diter=200,
    keep_all=True,
    keep_redundant_data=True,
    output_prefix='restart',
)

#start simulation
sim.run()
