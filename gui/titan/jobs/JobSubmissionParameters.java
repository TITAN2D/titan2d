package titan.jobs;

import titan.io.NameValuePairGroup;
import titan.io.NameValuePairSimple;


public class JobSubmissionParameters {

	public static final String QUEUE = "queue";
	public static final String MAX_RUN_TIME = "max_run_time";
	public static final String NODES = "nodes";
	public static final String CPUS = "cpus";
	public static final String NODES_SUBMIT = "nodes_submit";
	public static final String CPUS_SUBMIT = "cpus_submit";
	public static final String MAX_RUN_TIME_SUBMIT = "max_run_time_submit";
	public static final String REQUIREMENTS = "reqts";
	public static final String RUN_MODE = "run_mode";
	public static final String RUN_ARGS = "run_args";
	//use this for a flag telling if the output of the job should go into irods
	public static final String OUTPUT_IRODS = "output_irods";
	public static final String DATE = "date";
	public static final String SAVE_DIR = "save_dir";
	
	protected NameValuePairGroup jobParms;
	
	public JobSubmissionParameters() {
		jobParms = new NameValuePairGroup();
	}
	
	public void addJobParameter(String key, String value) {
		NameValuePairSimple nvps = new NameValuePairSimple(key, value);
		jobParms.addNameValuePair(nvps);
	}
	
	public String getJobParameter(String key) {
		return jobParms.getValue(key);
	}
	
}
