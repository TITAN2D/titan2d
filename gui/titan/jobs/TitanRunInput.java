package titan.jobs;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import java.util.Vector;

import titan.graphics.TextInput;
import titan.gui.TitanConstants;

// For Run Type Normal:

public class TitanRunInput {

    public static final String SIMULATION_DATA = new String("simulation.py");
    public static final String tab = new String("   ");

    private BufferedWriter writer;
    private String baseDir;

    // Send pile volume_fraction field for the TwoPhases-Pitman-Le model only
    private boolean TWOPHASES_PITMAN_LE_MODEL = false;

    public TitanRunInput(String dir) {

        baseDir = new String(dir);
    }

    public void writeSimulationData(TitanSimulationData data, TitanScaleData scaleData) {
        String fn = new String(baseDir + File.separator + SIMULATION_DATA);

        try {
            writer = new BufferedWriter(new FileWriter(fn));

            // Check flag

            String overwriteOutput = "";
            if (data.overwriteOutput == true) overwriteOutput = new String("True");
            else overwriteOutput = new String("False");

            writer.write("sim=TitanSimulation(overwrite_output=" + overwriteOutput + ")\n\n");

            writeSetGIS(writer, data);
            writeSetScale(writer, data, scaleData);
            writeSetNumProp(writer, data);
            // Send pile vol_fract field for the TwoPhases-Pitman-Le model only,
            // therefore call SetMatModel before writeAddPile
            writeSetMatModel(writer, data);
            writeSetTimeProps(writer, data);
            writeSetRestartOutput(writer, data);
            writeSetTimeSeriesOutput(writer, data);
            writeSetStatProps(writer, data);
            writeSetOutlineProps(writer, data);
            writeAddPile(writer, data);
            writeAddFluxSource(writer, data);
            writeAddDischargePlane(writer, data);

            writer.write("sim=sim.run()\n\n");

            writer.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    // sim.setGIS

    public void writeSetGIS(BufferedWriter writer, TitanSimulationData data) {

        try {

            writer.write("sim.setGIS(\n");

            // GIS format
            writer.write(tab + "gis_format='" + data.format + "',\n");

            if (data.format.compareTo(TitanConstants.GisFormats[TitanConstants.GIS_FORMAT_GIS_GRASS]) == 0) {

                // GIS_GRASS

                // GIS Main directory
                writer.write(tab + "gis_main='" + data.mainDirectory + "',\n");

                // GIS subdirectory
                writer.write(tab + "gis_sub='" + data.subDirectory + "',\n");

                // GIS Mapset
                writer.write(tab + "gis_mapset='" + data.mapset + "',\n");

                // GIS Map
                writer.write(tab + "gis_map='" + data.map + "',\n");

                // GIS Vector
                if (data.vector.compareTo("") == 0) {
                    writer.write(tab + "gis_vector=" + "None" + ",\n");
                } else {
                    writer.write(tab + "gis_vector='" + data.vector + "',\n");
                }
            } else {

                // GDAL

                // GIS Map
                // Need full path name for GDAL
                writer.write(tab + "gis_map='" + data.map + "',\n");
            }

            // Min and Max x and y locations

            if ((data.min_x_loc == TitanSimulationData.BLANK_FIELD) || (data.min_y_loc == TitanSimulationData.BLANK_FIELD)
                    || (data.max_x_loc == TitanSimulationData.BLANK_FIELD) || (data.max_y_loc == TitanSimulationData.BLANK_FIELD)) {
                writer.write(tab + "region_limits=" + "None" + "\n");
            } else {
                writer.write(tab + "region_limits=[" + data.min_x_loc + ", " + data.min_y_loc +
                        ", " + data.max_x_loc + ", " + data.max_y_loc + "]\n");
            }
            writer.write(")\n\n");
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    /**
     * @param data
     */
    // sim.setScale
    public void writeSetScale(BufferedWriter writer,
                              TitanSimulationData data,
                              TitanScaleData scaleData) {

        if (scaleData.scaleSim == true) {

            try {
                writer.write("sim.setScale(\n");
                writer.write(tab + "length_scale=" + scaleData.lengthScale + ",\n");
                writer.write(tab + "gravity_scale=" + scaleData.gravityScale + ",\n");

                if (scaleData.heightScale == TitanSimulationData.BLANK_FIELD) {
                    writer.write(tab + "height_scale=" + "None" + "\n");
                } else {
                    writer.write(tab + "height_scale=" + scaleData.heightScale + "\n");
                }

                writer.write(")\n\n");

            } catch (IOException ex) {
                ex.printStackTrace();
            }
        } else {
            try {
                writer.write("sim.setScale(\n");
                writer.write(tab + "length_scale=1.0,\n");
                writer.write(tab + "gravity_scale=1.0,\n");
                writer.write(tab + "height_scale=1.0\n");
                writer.write(")\n\n");
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
    }

    // sim.setNumProp

    public void writeSetNumProp(BufferedWriter writer, TitanSimulationData data) {

        try {

            String AMR = "";
            String methodOrder = "";
            String interfaceCapturingType = "";

            // Check flags

            if (data.AMR == true) AMR = new String("True");
            else AMR = new String("False");

            if (data.firstOrder) methodOrder = new String("First");
            else if (data.secondOrder) methodOrder = new String("Second");
            else {
                System.out.println("Method Order error.  Using first order");
                methodOrder = new String("First");
            }

            if (data.heuristicInterfaceCapturingType) interfaceCapturingType = new String("Heuristic");
            else if (data.levelSetInterfaceCapturingType) interfaceCapturingType = new String("LevelSet");
            else {
                System.out.println("Interface Capturing Type error.  Using Heuristic");
                interfaceCapturingType = new String("Heuristic");
            }

            writer.write("sim.setNumProp(\n");

            writer.write(tab + "AMR=" + AMR + ",\n");

            writer.write(tab + "number_of_cells_across_axis=" + data.cellsAcross + ",\n");

            writer.write(tab + "order='" + methodOrder + "',\n");

            writer.write(tab + "interface_capturing_type='" + interfaceCapturingType + "'\n");

            writer.write(")\n\n");
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    // sim.setMatModel

    public void writeSetMatModel(BufferedWriter writer, TitanSimulationData data) {

        try {

            boolean COULOMB_MODEL = false;
            int physicsModelIndex = 0;

            // Check flag

            String useMaterialMap = "";

            if (data.useMaterialMap) useMaterialMap = new String("True");
            else useMaterialMap = new String("False");

            if (data.physicsModel.compareTo(TitanConstants.PhysicsModels[TitanConstants.PHYSICS_MODEL_COULOMB]) == 0) {
                COULOMB_MODEL = true;
            }
            // Send pile vol_fract field for the TwoPhases-Pitman-Le physics model only
            if (data.physicsModel.compareTo(TitanConstants.PhysicsModels[TitanConstants.PHYSICS_MODEL_TWOPHASES_PITMAN_LE]) == 0) {
                TWOPHASES_PITMAN_LE_MODEL = true;
            }

            for (int i = 0; i < TitanConstants.PhysicsModels.length; i++) {
                if (data.physicsModel.compareTo(TitanConstants.PhysicsModels[i]) == 0) {
                    physicsModelIndex = i;
                    break;
                }
            }

            String[] matParmNames = TitanConstants.matParmNames[physicsModelIndex];
            String matParmName;

            writer.write("sim.setMatModel(\n");

             /* Material model */
            writer.write(tab + "model='" + data.physicsModel + "',\n");

            // Use a material map?
            // titan errors if use material map and not Coulomb
            if (COULOMB_MODEL == true) {
                writer.write(tab + "use_gis_matmap=" + useMaterialMap + ",\n");
            }

            // Send stopping_criteria for the Coulomb physics model only
            if (COULOMB_MODEL == true) {

                if (data.stoppingCriteria.compareTo("None") == 0) {
                    writer.write(tab + "stopping_criteria=None,\n");
                } else {
                    writer.write(tab + "stopping_criteria='" + data.stoppingCriteria + "',\n");
                }
            }

            // Material model parameters
            if ((COULOMB_MODEL == true) && (data.useMaterialMap == true)) {

                // int_frict
                writer.write(tab + matParmNames[TitanConstants.INT_FRICT_INDEX] + "=" +
                        data.matParmsVals[TitanConstants.INT_FRICT_INDEX] + ",\n");

                // bed_frict
                writer.write(tab + matParmNames[TitanConstants.BED_FRICT_INDEX] + "={\n");

                // The GIS cats directory _Mat file material names may not have the white spaces removed.
                // The binary file material names have the white spaces removed.
                // (See titan2d src/gisapi/Gispai.C and src/gisapi/GisCats.C).
                // Remove white spaces so that the titan2d main/titan_simulation.C check will pass
                String materialName;

                for (int i = 1; i < data.numPhysicsModelParameters; i++) {
                    materialName = TitanConstants.MATERIAL_MAP_NAMES[i - 1].replaceAll("\\s", "");
                    if (i < data.numPhysicsModelParameters - 1)
                        writer.write(tab + tab + "'" + materialName +
                                "':" + data.matParmsVals[i] + ",\n");
                    else writer.write(tab + tab + "'" + materialName +
                            "':" + data.matParmsVals[i] + "\n");
                }
                writer.write(tab + "}\n");
            } else {

                for (int i = 0; i < matParmNames.length; i++) {
                    matParmName = matParmNames[i];
                    if (i < matParmNames.length - 1)
                        writer.write(tab + matParmName + "=" + data.matParmsVals[i] + ",\n");
                    else writer.write(tab + matParmName + "=" + data.matParmsVals[i] + "\n");
                }
            }

            writer.write(")\n\n");
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    // sim.setTimeProps

    public void writeSetTimeProps(BufferedWriter writer, TitanSimulationData data) {

        try {

            writer.write("sim.setTimeProps(\n");

            // Maximum Number of Steps
            if (data.maxNumberTimeSteps == TitanSimulationData.BLANK_FIELD) {
                writer.write(tab + "max_iter=" + "None" + ",\n");
            } else {
                writer.write(tab + "max_iter=" + data.maxNumberTimeSteps + ",\n");
            }

            // Maximum Time
            if (data.maxTime == TitanSimulationData.BLANK_FIELD) {
                writer.write(tab + "max_time=" + "None" + "\n");
            } else {
                writer.write(tab + "max_time=" + data.maxTime + "\n");
            }

            writer.write(")\n\n");
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    // sim.setRestartOutput

    public void writeSetRestartOutput(BufferedWriter writer, TitanSimulationData data) {

        if (data.restartOutputEnabled == true) {

            try {

                // Check flags

                String keepAll = "";
                String keepRedundantData = "";

                if (data.keepAll == true) keepAll = new String("True");
                else keepAll = new String("False");
                if (data.keepRedundantData == true) keepRedundantData = new String("True");
                else keepRedundantData = new String("False");

                writer.write("sim.setRestartOutput(\n");

                // Time between result output
                if (data.resultOutputDelta1 == TitanSimulationData.BLANK_FIELD) {
                    writer.write(tab + "dtime=" + "None" + ",\n");
                } else {
                    writer.write(tab + "dtime=" + data.resultOutputDelta1 + ",\n");
                }

                // Time between saves
                if (data.saveDelta1 == TitanSimulationData.BLANK_FIELD) {
                    writer.write(tab + "diter=" + "None" + ",\n");
                } else {
                    writer.write(tab + "diter=" + data.saveDelta1 + ",\n");
                }

                writer.write(tab + "keep_all=" + keepAll + ",\n");
                writer.write(tab + "keep_redundant_data=" + keepRedundantData + ",\n");

                writer.write(tab + "output_prefix='" + data.outputPrefix1 + "'\n");

                writer.write(")\n\n");
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
    }

    // sim.setTimeSeriesOutput

    public void writeSetTimeSeriesOutput(BufferedWriter writer, TitanSimulationData data) {

        // 0 - 4 vizoutput
        String[] vizTypes = new String[4];
        int numVizTypes = 0;

        if (data.tecplot) {
            vizTypes[numVizTypes] = "tecplot";
            numVizTypes = numVizTypes + 1;
        }

        if (data.meshplot) {
            // meshplot?
            vizTypes[numVizTypes] = "meshplot";
            numVizTypes = numVizTypes + 1;
        }

        if (data.xdmf) {
            vizTypes[numVizTypes] = "xdmf";
            numVizTypes = numVizTypes + 1;
        }

        if (data.grass) {
            vizTypes[numVizTypes] = "grasssites";
            numVizTypes = numVizTypes + 1;
        }

        if (data.webviz) {
            // Previously handled by vhub version of titan2d only.
            // Ability to select webviz disabled.
            // Stub for adding handling in this version of titan
        }

        if (data.webviz) {
            // Previously handled by vhub version of titan2d only.
            // Ability to select gmfg disabled.
            // Stub for adding handling in this version of titan
        }

        if (numVizTypes > 0) {

            try {

                writer.write("sim.setTimeSeriesOutput(\n");

                writer.write(tab + "vizoutput=(");

                for (int i = 0; i < numVizTypes - 1; i++) {
                    writer.write("'" + vizTypes[i] + "',");
                }
                writer.write("'" + vizTypes[numVizTypes - 1] + "'),\n");

                // Time between result output
                if (data.resultOutputDelta2 == TitanSimulationData.BLANK_FIELD) {
                    writer.write(tab + "dtime=" + "None" + ",\n");
                } else {
                    writer.write(tab + "dtime=" + data.resultOutputDelta2 + ",\n");
                }

                // Time between saves
                if (data.saveDelta2 == TitanSimulationData.BLANK_FIELD) {
                    writer.write(tab + "diter=" + "None" + ",\n");
                } else {
                    writer.write(tab + "diter=" + data.saveDelta2 + ",\n");
                }

                writer.write(tab + "output_prefix='" + data.outputPrefix2 + "'\n");

                writer.write(")\n\n");
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
    }

    // sim.setStatProps

    public void writeSetStatProps(BufferedWriter writer, TitanSimulationData data) {

        try {

            writer.write("sim.setStatProps(\n");

            // Run ID
            if (data.runID == TitanSimulationData.BLANK_FIELD)
                writer.write(tab + "runid=-1,\n");
            else writer.write(tab + "runid=" + data.runID + ",\n");

            // Height used to define flow outline (-1 if blank)
            if (data.flowOutlineHeight == TitanSimulationData.BLANK_FIELD)
                writer.write(tab + "edge_height=None,\n");
            else writer.write(tab + "edge_height=" + data.flowOutlineHeight + ",\n");

            // Test if flow reaches height ... (-2 if blank)
            if (data.flowHeight == TitanSimulationData.BLANK_FIELD)
                writer.write(tab + "test_height=None,\n");
            else writer.write(tab + "test_height=" + data.flowHeight + ",\n");

            String testX;
            String testY;

            if ((data.testX == TitanSimulationData.BLANK_FIELD) ||
                    (data.testY == TitanSimulationData.BLANK_FIELD)) {
                writer.write(tab + "test_location=None,\n");
            } else {
                testX = new String(Float.toString(data.testX));
                testY = new String(Float.toString(data.testY));
                writer.write(tab + "test_location=(" + testX + ", " + testY + "),\n");
            }

            writer.write(tab + "output_prefix='" + data.statPropsPrefix + "'\n");

            writer.write(")\n\n");
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    // sim.setOutlineProps

    public void writeSetOutlineProps(BufferedWriter writer, TitanSimulationData data) {

        try {

            // Check flag

            String outlinePropsEnabled = "";

            if (data.outlinePropsEnabled == true) outlinePropsEnabled = new String("True");
            else outlinePropsEnabled = new String("False");

            writer.write("sim.setOutlineProps(\n");

            writer.write(tab + "enabled=" + outlinePropsEnabled + ",\n");

            // If blank, titan will use the default value for this
            if (data.maxLinearSize != TitanSimulationData.BLANK_FIELD) {
                writer.write(tab + "max_linear_size=" + data.maxLinearSize + ",\n");
            }

            writer.write(tab + "init_size='" + data.initSize + "',\n");

            writer.write(tab + "output_prefix='" + data.outlinePropsPrefix + "'\n");

            writer.write(")\n\n");
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    // sim.addPile

    public void writeAddPile(BufferedWriter writer, TitanSimulationData data) {

        try {

			/*
             * Pile information
			 */
            for (int pile = 0; pile < data.pile.length; pile++) {

                writer.write("sim.addPile(\n");

                if (data.pile[pile].type == TitanConstants.PILE_TYPE_PARABOLOID_SELECT) {
                    writer.write(tab + "pile_type='" + TitanConstants.PILE_TYPE_PARABOLOID + "',\n");
                } else {
                    writer.write(tab + "pile_type='" + TitanConstants.PILE_TYPE_CYLINDER + "',\n");
                }
                writer.write(tab + "height=" + data.pile[pile].maxThickness + ",\n");
                writer.write(tab + "center=[" + data.pile[pile].centerVolumeX + ", " + data.pile[pile].centerVolumeY + "],\n");
                writer.write(tab + "radii=[" + data.pile[pile].majorExtent + ", " + data.pile[pile].minorExtent + "],\n");
                writer.write(tab + "orientation=" + data.pile[pile].orientation + ",\n");
                writer.write(tab + "Vmagnitude=" + data.pile[pile].speed + ",\n");

                // Send pile vol_fract field for the TwoPhases-Pitman-Le model only
                if (TWOPHASES_PITMAN_LE_MODEL == true) {
                    writer.write(tab + "Vdirection=" + data.pile[pile].direction + ",\n");
                    writer.write(tab + "vol_fract" + "=" + data.pile[pile].volumeFraction + "\n");
                } else {
                    writer.write(tab + "Vdirection=" + data.pile[pile].direction + "\n");
                }

                writer.write(")\n\n");
            }
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    // sim.addFlux

    public void writeAddFluxSource(BufferedWriter writer, TitanSimulationData data) {

        try {

			/*
             * Flux Source specific information (grouped by source) :
             */
            for (int flux = 0; flux < data.fluxSource.length; flux++) {

                writer.write("sim.addFluxSource(\n");

                writer.write(tab + "influx=" + data.fluxSource[flux].extrusionFluxRate + ",\n");
                writer.write(tab + "start_time=" + data.fluxSource[flux].activeTimeStart + ",\n");
                writer.write(tab + "end_time=" + data.fluxSource[flux].activeTimeEnd + ",\n");
                writer.write(tab + "center=[" + data.fluxSource[flux].centerX + ", " + data.fluxSource[flux].centerY + "],\n");
                writer.write(tab + "radii=[" + data.fluxSource[flux].majorExtent + ", " + data.fluxSource[flux].minorExtent + "],\n");
                writer.write(tab + "orientation=" + data.fluxSource[flux].orientationAngle + ",\n");
                writer.write(tab + "Vmagnitude=" + data.fluxSource[flux].initSpeed + ",\n");
                writer.write(tab + "Vdirection=" + data.fluxSource[flux].initDirection + "\n");

                writer.write(")\n\n");
            }
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    // sim.addDischargePlane

    public void writeAddDischargePlane(BufferedWriter writer, TitanSimulationData data) {

        try {

            /*
             * Discharge Plane Information :
             */

            for (int plane = 0; plane < data.plane.length; plane++) {

                writer.write("sim.addDischargePlane(");

                writer.write(data.plane[plane].a_e + ", ");
                writer.write(data.plane[plane].a_n + ", ");
                writer.write(data.plane[plane].b_e + ", ");
                writer.write(data.plane[plane].b_n + "");

                writer.write(")\n\n");
            }
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    public class TitanSimulationData {

        // BLANK_FIELD if field is left blank
        public final static int BLANK_FIELD = -32767;

        public TitanSimulationPile[] pile = new TitanSimulationPile[0];
        public TitanSimulationFluxSource[] fluxSource = new TitanSimulationFluxSource[0];
        public TitanSimulationPlane[] plane = new TitanSimulationPlane[0];

        // GIS Tab
        public String format;
        public String mainDirectory;
        public String subDirectory;
        public String mapset;
        public String map;
        public String vector;
        // This GIS tab Zone Override and Hemisphere parameters are used when creating KML files
        // when the KML button is clicked on the JobDetailsDialog tab
        public int zoneOverride = 0;
        public String hemisphere;
        public float min_x_loc = BLANK_FIELD;
        public float min_y_loc = BLANK_FIELD;
        public float max_x_loc = BLANK_FIELD;
        public float max_y_loc = BLANK_FIELD;

        // Run Parameters Tab
        public int cellsAcross = 0;
        public int maxNumberTimeSteps = BLANK_FIELD;
        public float maxTime = BLANK_FIELD;
        public boolean AMR = false;
        // Method order
        public boolean firstOrder = false;
        public boolean secondOrder = false;
        // Interface capturing type
        public boolean heuristicInterfaceCapturingType = false;
        public boolean levelSetInterfaceCapturingType = false;

        // Stat Props
        public int runID = BLANK_FIELD;
        public float flowOutlineHeight = BLANK_FIELD;
        public float flowHeight = BLANK_FIELD;
        public float testX = BLANK_FIELD;
        public float testY = BLANK_FIELD;
        public String statPropsPrefix;

        // Outline Props
        public boolean outlinePropsEnabled = false;
        public float maxLinearSize = BLANK_FIELD;
        public String initSize;
        public String outlinePropsPrefix;

        // Job Options Tab

        public boolean overwriteOutput = false;

        public boolean restartOutputEnabled = false;
        public float resultOutputDelta1 = BLANK_FIELD;
        public int saveDelta1 = BLANK_FIELD;
        public boolean keepAll = false;
        public boolean keepRedundantData = false;
        public String outputPrefix1;

        // Output types
        public boolean tecplot = false;
        public boolean meshplot = false;
        public boolean xdmf = false;
        public boolean grass = false;
        public boolean webviz = false;
        public boolean gmfg = false;

        public float resultOutputDelta2 = BLANK_FIELD;
        public int saveDelta2 = BLANK_FIELD;
        public String outputPrefix2;

        // Material Model and Map tab

        public String physicsModel;

        //public String matParmsVals[] = new String[TitanConstants.MAX_PHYSICS_MODEL_PARAMETERS];

        public float matParmsVals[] = new float[TitanConstants.MAX_PHYSICS_MODEL_PARAMETERS];

        ///for (int i=0; i<TitanConstants.MAX_PHYSICS_MODEL_PARAMETERS; i++){
        //matParmsVals[i] = BLANK_FIELD;
        //}

        public boolean useMaterialMap;

        public int numPhysicsModelParameters;

        public String stoppingCriteria;
    }

    public class TitanScaleData {

        // BLANK_FIELD if field is left blank
        public final static int BLANK_FIELD = -32767;

        public boolean scaleSim = false;
        public float lengthScale = BLANK_FIELD;
        public float gravityScale = BLANK_FIELD;
        public float heightScale = BLANK_FIELD;
    }

    public class TitanSimulationPile {
        public float type = 0.0f;
        public float maxThickness = 0.0f;
        public float centerVolumeX = 0.0f;
        public float centerVolumeY = 0.0f;
        public float volumeFraction = 0.0f;
        public float majorExtent = 0.0f;
        public float minorExtent = 0.0f;
        public float orientation = 0.0f;
        public float speed = 0.0f;
        public float direction = 0.0f;

        public String toString() {
            return new String(type + "\n" +
                    maxThickness + "\n" +
                    centerVolumeX + "\n" +
                    centerVolumeY + "\n" +
                    volumeFraction + "\n" +
                    majorExtent + "\n" +
                    minorExtent + "\n" +
                    orientation + "\n" +
                    speed + "\n" +
                    direction);
        }
    }

    public class TitanSimulationFluxSource {
        public float centerX = 0.0f;
        public float centerY = 0.0f;
        public float majorExtent = 0.0f;
        public float minorExtent = 0.0f;
        public float orientationAngle = 0.0f;
        public float initSpeed = 0.0f;
        public float initDirection = 0.0f;
        public float extrusionFluxRate = 0.0f;
        public float activeTimeStart = 0.0f;
        public float activeTimeEnd = 0.0f;

        public String toString() {

            return new String(extrusionFluxRate + "\n" +
                    activeTimeStart + "\n" +
                    activeTimeEnd + "\n" +
                    centerX + "\n" +
                    centerY + "\n" +
                    majorExtent + "\n" +
                    minorExtent + "\n" +
                    orientationAngle + "\n" +
                    initSpeed + "\n" +
                    initDirection);

        }
    }

    public class TitanSimulationPlane {
        public float a_e = 0.0f;
        public float a_n = 0.0f;
        public float b_e = 0.0f;
        public float b_n = 0.0f;

        public String toString() {
            return new String(a_e + " " + a_n + " " + b_e + " " + b_n);
        }
    }
}
