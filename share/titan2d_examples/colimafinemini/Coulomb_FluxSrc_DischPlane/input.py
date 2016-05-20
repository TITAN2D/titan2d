sim=TitanSimulation(overwrite_output=True)

sim.setGIS(
    gis_format='GIS_GRASS', 
    gis_main='../../dem',
    gis_sub='colimafinemini',
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
    model='Coulomb',
    int_frict=37.0,
    bed_frict=27.0
)
sim.addPile(
    pile_type='Cylinder',
    height=30.0,
    center=[644956.0, 2157970.0],
    radii=[55.0, 55.0],
    orientation=0.0,
    Vmagnitude=0.0,
    Vdirection=0.0,
)
sim.addFluxSource(
    influx=10.0,
    start_time=5.0,
    end_time=20.0,
    center=[644956.0, 2157970.0],
    radii=[55.0, 55.0],
    orientation=0.0,
    Vmagnitude=0.0,
    Vdirection=0.0
)
sim.addDischargePlane(637380.0, 2145800.0, 664380.0, 2169800.0)

sim.setTimeProps(
    max_iter=5000,
)

sim.setTimeSeriesOutput(
    vizoutput=('xdmf','meshplot'),
    dtime=1.0,
)

#start simulation
sim.run()
