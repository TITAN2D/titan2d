package titan.jobs;

import titan.gui.TitanConstants;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public abstract class RunMethod {
	public static final String TITAN_SUCCESSFUL = new String(".titan_successful");
	public static final String TITAN_FAILED = new String(".titan_failed");
	public static final String STDOUT = new String("Stdout.txt");
	public static final String STDERR = new String("Stderr.txt");
	public static final String GIS_DATA_TAR = new String("gis_data.tar");
	public static String execName;
	protected String basedir;
	protected String jobName;
	protected BufferedWriter writer;
	protected String[] gisParms;
	protected boolean binMode;
	protected String binDir;
	protected String restartEnabled;
	protected String restartFile;
	protected int restartMaxNumberTimeSteps;
	protected float restartMaxTime;
	
	public RunMethod(String srestartEnabled, String dir, String date, String sexecName) {
		restartEnabled = srestartEnabled;
		basedir = new String(dir);
		jobName = new String("titan_" + date);
		execName = sexecName;
	}

	public boolean allowMultiProcessors() { return false; }

	public void initialize(GisParms gParms,
			               boolean useExternal,
						   String binaryDir,
                           String preCommands,
                           JobSubmissionParameters jobParms,
						   String srestartFile,
						   int irestartMaxNumberTimeSteps,
						   float frestartMaxTime) {
		
		int numProcessors =
			Integer.parseInt(jobParms.getJobParameter(JobSubmissionParameters.NODES)) *
			Integer.parseInt(jobParms.getJobParameter(JobSubmissionParameters.CPUS));

        gisParms = new String[12];
        gisParms[0] = new String(Integer.toString(numProcessors));
        gisParms[1] = new String(gParms.NUM_CELLS_ACROSS);
        gisParms[2] = new String(gParms.GIS_INFO_DIRECTORY);
        gisParms[3] = new String(gParms.GIS_SUBDIR);
        gisParms[4] = new String(gParms.GIS_MAPSET);
        gisParms[5] = new String(gParms.GIS_MAP);
        gisParms[6] = new String(gParms.MIN_X_LOC);
        gisParms[7] = new String(gParms.MIN_Y_LOC);
        gisParms[8] = new String(gParms.MAX_X_LOC);
        gisParms[9] = new String(gParms.MAX_Y_LOC);
		gisParms[10] = new String(gParms.ZONEOVERRIDE);
        gisParms[11] = new String(gParms.HEMISPHERE);

        binMode = useExternal;
        binDir = binaryDir;
		restartFile = srestartFile;
		restartMaxNumberTimeSteps = irestartMaxNumberTimeSteps;
		restartMaxTime = frestartMaxTime;

	}
	
	public JobDetails runJob() { return null; }

	/**
	 * 
	 * @param gisData Directory where GIS data resides
	 */
	protected void zipGISData(String gisData) {
		try {

			String cdDir = gisData.substring(0, gisData.lastIndexOf("/"));
			String tarDir = gisData.substring(gisData.lastIndexOf("/")+1);
			Process tarJob = Runtime.getRuntime().exec("tar cvf " +
		    		basedir + File.separator + GIS_DATA_TAR + " -C" + cdDir + " " + tarDir);
			tarJob.waitFor();
			Process zipJob = Runtime.getRuntime().exec("gzip " + basedir +
					File.separator + GIS_DATA_TAR);
			zipJob.waitFor();
		}
		catch (IOException ex) { ex.printStackTrace(); }
		catch (InterruptedException ex) { ex.printStackTrace(); }
	}
	
	public class GisParms {
		public String NUM_CELLS_ACROSS = null;
		public String GIS_INFO_DIRECTORY = null;
		public String GIS_SUBDIR = null;
		public String GIS_MAPSET = null;
		public String GIS_MAP = null;
		public String MIN_X_LOC = null;
		public String MIN_Y_LOC = null;
		public String MAX_X_LOC = null;
		public String MAX_Y_LOC = null;
		public String GIS_VECTOR = null;
		public String ZONEOVERRIDE = null;
		public String HEMISPHERE = null;
	}

}
