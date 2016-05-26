package titan.jobs;

import java.io.File;

public abstract class JobStatus {

	public static final String STATUS_UNKNOWN    = "UNKNOWN";
	public static final String STATUS_RUNNING    = "Running";
	public static final String STATUS_EXITING    = "Exiting";
	public static final String STATUS_HELD       = "Held";
	public static final String STATUS_QUEUED     = "Queued";
	public static final String STATUS_MOVING     = "Being Moved";
	public static final String STATUS_WAITING    = "Waiting";
	public static final String STATUS_SUSPENDED  = "Suspended";
	public static final String STATUS_COMPLETED  = "Completed";
	public static final String STATUS_UNEXPANDED = "Unexpanded";
	public static final String STATUS_REMOVED    = "Removed";
	//public static final String STATUS_PREPROCESS_FAIL = "Preprocess Failed";
	//public static final String STATUS_VECDATA_FAIL    = "Vecdata Failed";
	public static final String STATUS_TITAN_FAIL      = "Titan Failed";
	public static final String STATUS_SUCCESS         = "Successful";
	
	/**
	 * 
	 * @param outputDir The output directory for the run on the original host
	 * @return Returns true if the output directory can be found and the
	 *         file .complete exists in that directory
	 */
	protected static String jobStatusFromScript(String outputDir) {
		String status = STATUS_UNKNOWN;
    	if(new File(outputDir).exists()) {
    		//if(new File(outputDir + File.separator + RunMethod.PREPROCESS_FAILED).exists()) {
    		//	status = STATUS_PREPROCESS_FAIL;
    		//}
    		//else if(new File(outputDir + File.separator + RunMethod.VECDATA_FAILED).exists()) {
    		//	status = STATUS_VECDATA_FAIL;
    		//}
    		if(new File(outputDir + File.separator + RunMethod.TITAN_FAILED).exists()) {
    			status = STATUS_TITAN_FAIL;
    		}
    		else if(new File(outputDir + File.separator + RunMethod.TITAN_SUCCESSFUL).exists()) {
    			status = STATUS_SUCCESS;
    		}
    		// Job is probably running, but not really sure
    		else {
    			status = STATUS_UNKNOWN;
    		}
    	}
    	return status;
	}
	
	public static String getJobStatus(String jobID, String outputDir) { return STATUS_UNKNOWN; }
	
}
