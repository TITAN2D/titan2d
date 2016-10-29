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
    model='Pouliquen-Forterre',
    phi1=32.9,
    phi2=42.0,
    phi3=33.9,
    Beta=0.65,
    L_material=1.0E-3,
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

sim.setTimeProps(
    max_iter=5000,
)

sim.setTimeSeriesOutput(
    vizoutput=('xdmf','meshplot'),
    dtime=1.0,
)

#start simulation
sim.run()
