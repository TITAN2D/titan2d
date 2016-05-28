package titan.gui;

public abstract class TitanConstants {

	// Lists
	public static final String LIST_SEPARATOR = "/=";

	// Parameter file extensions
	public static final String PARM_FILE_EXT = "ascprm";
	public static final String FLUX_FILE_EXT = "ascflux";
	public static final String PILE_FILE_EXT = "ascpile";
	public static final String PLANE_FILE_EXT = "ascplane";

	public static final String FALSE = "False";
	public static final String TRUE = "True";
	public static final String[] TrueFalse = { TRUE, FALSE };

	// Parameter Names from ascprm file

	// Load/Save Tab  Parameters
	// Parameter names for loading/saving runs
	public static final String LOAD_DIRECTORY = "LOAD_DIRECTORY";
	public static final String SAVE_DIRECTORY = "SAVE_DIRECTORY";
	public static final String RUN_NAME = "RUN_NAME";

	public static final String SIM_DIR = "SIM_DIR";

	// GIS Tab Parameters
	public static final String GIS_FORMAT = "GIS_FORMAT";
	public static final String[] GisFormats = { "GIS_GRASS", "GDAL" };
	public static final int GIS_FORMAT_GIS_GRASS = 0;
	public static final int GIS_FORMAT_GDAL = 1;

	public static final String GIS_INFO_DIRECTORY = "GIS_INFO_DIRECTORY";
	public static final String GIS_SUBDIR = "GIS_SUBDIR";
	public static final String GIS_MAPSET = "GIS_MAPSET";
	public static final String GIS_MAP = "GIS_MAP";
	public static final String GIS_VECTOR = "GIS_VECTOR";
	public static final String MIN_X_LOC = "MIN_X_LOC";
	public static final String MIN_Y_LOC = "MIN_Y_LOC";
	public static final String MAX_X_LOC = "MAX_X_LOC";
	public static final String MAX_Y_LOC = "MAX_Y_LOC";

	// Parameters used when creating KML files
	public static final String ZONEOVERRIDE = "ZONEOVERRIDE";
	public static final String HEMISPHERE = "HEMISPHERE";
	public static final String[] HemisphereHalves = { "North", "South" };
	public static final int HEMISPHERE_NORTH = 0;
	public static final int HEMISPHERE_SOUTH = 1;

	public static final String EMAIL_ADDRESS = "EMAIL_ADDRESS";

	// Run Parameters Tab Parameters
	public static final String RESTART_ENABLED = "RESTART_ENABLED";
	public static final String RESTART_FILE = "RESTART_FILE";
	public static final String RESTART_LOAD_DIRECTORY = "RESTART_LOAD_DIRECTORY";
	public static final String RESTART_SAVE_DIRECTORY = "RESTART_SAVE_DIRECTORY";
	public static final String RESTART_RUN_NAME = "RESTART_RUN_NAME";
	public static final String NUM_CELLS_ACROSS = "NUM_CELLS_ACROSS";
	public static final String MAX_NUM_TIME_STEPS = "MAX_NUM_TIME_STEPS";
	public static final String MAX_TIME = "MAX_TIME";
	public static final String AMR = "AMR";
	public static final String ORDER_METHOD = "ORDER_METHOD";
	public static final String SCALE_SIM = "SCALE_SIM";
	public static final String LENGTH_SCALE = "LENGTH_SCALE";
	public static final String GRAVITY_SCALE = "GRAVITY_SCALE";
	public static final String HEIGHT_SCALE = "HEIGHT_SCALE";
	public static final String RUN_ID = "RUN_ID";
	public static final String FLOW_OUTLINE_HGT = "FLOW_OUTLINE_HGT";
	public static final String TEST_FLOW_HEIGHT_MIN = "TEST_FLOW_HEIGHT_MIN";
	public static final String TEST_FLOW_X_LOC = "TEST_FLOW_X_LOC";
	public static final String TEST_FLOW_Y_LOC = "TEST_FLOW_Y_LOC";
	public static final String STAT_PROPS_PREFIX = "STAT_PROPS_PREFIX";
	public static final String OUTLINE_PROPS_ENABLED = "OUTLINE_PROPS_ENABLED";
	public static final String MAX_LINEAR_SIZE = "MAX_LINEAR_SIZE";
	public static final String INIT_SIZE = "INIT_SIZE";
	public static final String OUTLINE_PROPS_PREFIX = "OUTLINE_PROPS_PREFIX";

	public static final String[] OrderMethod = {"First", "Second"};
	public static final int ORDER_METHOD_FIRST = 0;
	public static final int ORDER_METHOD_SECOND = 1;

	public static final String[] InitSizes = {"AMR", "DEM"};
	public static final int INIT_SIZE_AMR = 0;
	public static final int INIT_SIZE_DEM = 1;

	// Output Options Tab Parameters
	public static final String OVERWRITE_OUTPUT = "OVERWRITE_OUTPUT";
	public static final String RESTART_OUTPUT_ENABLED = "RESTART_OUTPUT_ENABLED";
	public static final String RESULT_OUTPUT_TIME_DELTA1 = "RESULT_OUTPUT_TIME_DELTA1";
	public static final String SAVE_TIME_DELTA1 = "SAVE_TIME_DELTA1";
	public static final String KEEP_ALL = "KEEP_ALL";
	public static final String KEEP_REDUNDANT_DATA = "KEEP_REDUNDANT_DATA";
	public static final String OUTPUT_PREFIX1 = "OUTPUT_PREFIX1";
	public static final String VIS_OUTPUT = "VIS_OUTPUT";
	public static final String RESULT_OUTPUT_TIME_DELTA2 = "RESULT_OUTPUT_TIME_DELTA2";
	public static final String SAVE_TIME_DELTA2 = "SAVE_TIME_DELTA2";
	public static final String OUTPUT_PREFIX2 = "OUTPUT_PREFIX2";

	public static final String[] VizTypes = {"tecplotxxx.tec", "mshplotxxx.tec",
			"XDMF/Paraview", "grass_sites",
			"Web Viz", "GMFG Viz"};
	public static final int VIS_OUTPUT_TECPLOT = 0;
	public static final int VIS_OUTPUT_MSHPLOT = 1;
	public static final int VIS_OUTPUT_XDMF = 2;
	public static final int VIS_OUTPUT_GRASS = 3;
	public static final int VIS_OUTPUT_WEBVIZ = 4;
	public static final int VIS_OUTPUT_GMFGVIZ = 5;

	// Material Model and Map Tab Parameters
	public static final String PHYSICS_MODEL = "PHYSICS_MODEL";
	public static final String NUM_PHYSICS_MODEL_PARAMETERS = "NUM_PHYSICS_MODEL_PARAMETERS";
	public static final String USE_GIS_MAT_MAP = "USE_GIS_MAT_MAP";
	public static final String STOPPING_CRITERIA = "STOPPING_CRITERIA";

	public static final String[] PhysicsModels = {"Coulomb", "TwoPhases-Pitman-Le", "Voellmy-Salm", "Pouliquen-Forterre"};
	public static final int PHYSICS_MODEL_COULOMB = 0;
	public static final int PHYSICS_MODEL_TWOPHASES_PITMAN_LE = 1;
	public static final int PHYSICS_MODEL_VOELLMY_SALM = 2;
	public static final int PHYSICS_MODEL_POULIQUEN_FORTERRE = 3;

	// Allows up to 6 material map materials.
	public static final int MAX_PHYSICS_MODEL_PARAMETERS = 7;
	public static final int MAX_MAT_MAP_MATERIALS = MAX_PHYSICS_MODEL_PARAMETERS - 1;

	public static final String MATERIAL_MODEL_TEXT_INPUT_0 = "MATERIAL_MODEL_TEXT_INPUT_0";
	public static final String MATERIAL_MODEL_TEXT_INPUT_1 = "MATERIAL_MODEL_TEXT_INPUT_1";
	public static final String MATERIAL_MODEL_TEXT_INPUT_2 = "MATERIAL_MODEL_TEXT_INPUT_2";
	public static final String MATERIAL_MODEL_TEXT_INPUT_3 = "MATERIAL_MODEL_TEXT_INPUT_3";
	public static final String MATERIAL_MODEL_TEXT_INPUT_4 = "MATERIAL_MODEL_TEXT_INPUT_4";
	public static final String MATERIAL_MODEL_TEXT_INPUT_5 = "MATERIAL_MODEL_TEXT_INPUT_5";
	public static final String MATERIAL_MODEL_TEXT_INPUT_6 = "MATERIAL_MODEL_TEXT_INPUT_6";

	public static final int INT_FRICT_INDEX = 0;
	public static final int BED_FRICT_INDEX = 1;

	public static final String[] CoulombParmNames =
			{"int_frict", "bed_frict"};

	public static final String[] TwoPhasesPitmanLeParmNames =
			{"int_frict", "bed_frict"};

	public static final String[] VoellmySalmParmNames =
			{"mu", "xi"};

	public static final String[] PouliquenForterreParmNames =
			{"phi1", "phi2", "phi3", "Beta", "L_material"};

	public static final String[] CoulombParmUnits =
			{"[deg]", "[deg]"};
	public static final String[] TwoPhasesPitmanLeParmUnits =
			{"[deg]", "[deg]"};
	public static final String[] VoellmySalmParmUnits =
			{"[deg]", "[deg]"};
	public static final String[] PouliquenForterreParmUnits =
			{"[deg]", "[deg]", "[deg]", "", ""};

	public static final String[][] matParmNames =
			{CoulombParmNames, TwoPhasesPitmanLeParmNames, VoellmySalmParmNames, PouliquenForterreParmNames};

	public static final String[][] matParmUnits =
			{CoulombParmUnits, TwoPhasesPitmanLeParmUnits, VoellmySalmParmUnits, PouliquenForterreParmUnits};

	// Note: these are not constants;
	// These get read in from the material map file.
	// Allows the TabMaterialModelMap class to set and access and the TitanRunInput class to access
	public static String[] MATERIAL_MAP_NAMES = new String[MAX_MAT_MAP_MATERIALS];

	public static final String[] StoppingCriterias = {"None", "DragBased"};
	public static final int STOPPING_CRITERIA_NONE = 0;
	public static final int STOPPING_CRITERIA_DRAGBASED = 1;

	// Flux Source Tab Parameters
	public static final String SRC_EXTRUSION_FLUX_RATE = "SRC_EXTRUSION_FLUX_RATE";
	public static final String SRC_ACTIVE_TIME_START = "SRC_ACTIVE_TIME_START";
	public static final String SRC_ACTIVE_TIME_END = "SRC_ACTIVE_TIME_END";
	public static final String SRC_CENTER_XC = "SRC_CENTER_XC";
	public static final String SRC_CENTER_YC = "SRC_CENTER_YC";
	public static final String SRC_MAJOR_EXTENT = "SRC_MAJOR_EXTENT";
	public static final String SRC_MINOR_EXTENT = "SRC_MINOR_EXTENT";
	public static final String SRC_ORIENTATION_ANGLE = "SRC_ORIENTATION_ANGLE";
	public static final String SRC_INIT_SPEED = "SRC_INIT_SPEED";
	public static final String SRC_INIT_DIRECTION = "SRC_INIT_DIRECTION";
	
	// Piles Tab Parameters
	public static final String PILE_TYPE = "PILE_TYPE";
	public static final String PILE_MAX_INIT_THICKNESS = "PILE_MAX_INIT_THICKNESS";
	public static final String PILE_CENTER_INIT_VOLUME_XC = "PILE_CENTER_INIT_VOLUME_XC";
	public static final String PILE_CENTER_INIT_VOLUME_YC = "PILE_CENTER_INIT_VOLUME_YC";
	public static final String PILE_VOLUME_FRACTION = "VOLUME_FRACTION";
	public static final String PILE_MAJOR_EXTENT = "PILE_MAJOR_EXTENT";
	public static final String PILE_MINOR_EXTENT = "PILE_MINOR_EXTENT";
	public static final String PILE_ORIENTATION_ANGLE = "PILE_ORIENTATION_ANGLE";
	public static final String PILE_INIT_SPEED = "PILE_INIT_SPEED";
	public static final String PILE_INIT_DIRECTION = "PILE_INIT_DIRECTION";

	public static final String PILE_TYPE_PARABOLOID = "Paraboloid";
	public static final String PILE_TYPE_CYLINDER = "Cylinder";

	public static final int PILE_TYPE_PARABOLOID_SELECT = 1;
	public static final int PILE_TYPE_CYLINDER_SELECT = 2;

	// Discharge Planes Tab Parameters
	public static final String PLANE_A_E = "PLANE_A_E";
	public static final String PLANE_A_N = "PLANE_A_N";
	public static final String PLANE_B_E = "PLANE_B_E";
	public static final String PLANE_B_N = "PLANE_B_N";

	public static final String[] parameterKeyList = {
			SIM_DIR, ZONEOVERRIDE, HEMISPHERE, EMAIL_ADDRESS,
			GIS_FORMAT,GIS_INFO_DIRECTORY, GIS_SUBDIR, GIS_MAPSET, GIS_MAP, GIS_VECTOR,
			MIN_X_LOC, MIN_Y_LOC, MAX_X_LOC, MAX_Y_LOC,
	        RESTART_ENABLED, RESTART_FILE, MAX_NUM_TIME_STEPS, MAX_TIME, NUM_CELLS_ACROSS, AMR, ORDER_METHOD,
			SCALE_SIM, LENGTH_SCALE, GRAVITY_SCALE, HEIGHT_SCALE,
			RUN_ID, FLOW_OUTLINE_HGT, TEST_FLOW_HEIGHT_MIN, TEST_FLOW_X_LOC, TEST_FLOW_Y_LOC, STAT_PROPS_PREFIX,
			OUTLINE_PROPS_ENABLED, MAX_LINEAR_SIZE, INIT_SIZE, OUTLINE_PROPS_PREFIX,
			OVERWRITE_OUTPUT,
			RESTART_OUTPUT_ENABLED, RESULT_OUTPUT_TIME_DELTA1, SAVE_TIME_DELTA1, KEEP_ALL,  KEEP_REDUNDANT_DATA, OUTPUT_PREFIX1,
			VIS_OUTPUT, RESULT_OUTPUT_TIME_DELTA2, SAVE_TIME_DELTA2, OUTPUT_PREFIX2,
			MATERIAL_MODEL_TEXT_INPUT_0, MATERIAL_MODEL_TEXT_INPUT_1, MATERIAL_MODEL_TEXT_INPUT_2,
			MATERIAL_MODEL_TEXT_INPUT_3, MATERIAL_MODEL_TEXT_INPUT_4, MATERIAL_MODEL_TEXT_INPUT_5, MATERIAL_MODEL_TEXT_INPUT_6,
			PHYSICS_MODEL, USE_GIS_MAT_MAP, STOPPING_CRITERIA, NUM_PHYSICS_MODEL_PARAMETERS};

	public static final String[] parameterValueDefaults = {
			"","",HemisphereHalves[0],"",
		    GisFormats[GIS_FORMAT_GIS_GRASS], "", "", "", "", "",
			"", "", "", "",
		    FALSE, "", "", "", "20", TRUE, OrderMethod[ORDER_METHOD_FIRST],
			FALSE, "", "9.8", "",
			"", "", "", "", "","",
			TRUE, "1024", InitSizes[INIT_SIZE_AMR], "",
			TRUE,
			FALSE, "", "",  FALSE, FALSE, "restart",
			"", "", "", "vizout",
			"37.0", "27.0", "", "", "", "", "",
			PhysicsModels[PHYSICS_MODEL_COULOMB], FALSE, StoppingCriterias[STOPPING_CRITERIA_NONE], ""+CoulombParmNames.length};
}

