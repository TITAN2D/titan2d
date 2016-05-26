package titan.jobs.local;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

import titan.jobs.JobStatus;

public class JobStatusLocal extends JobStatus {

	public static String getJobStatus(String jobID, String outputDir) {
		String status = STATUS_UNKNOWN;
		try {

		    Process psJob = Runtime.getRuntime().exec("ps -u");
		    
		    BufferedReader rs = new BufferedReader(new InputStreamReader(psJob.getInputStream()));

		    String line;
		    boolean jobFound = false;
		    while((line = rs.readLine()) != null) {
	    	    line = line.substring(line.indexOf(" ") + 1);
	    	    while(line.startsWith(" ")) line = line.substring(1);
	    	    line = line.substring(0, line.indexOf(" "));
	    	    if(line.compareTo(jobID) == 0) {
	    	        status = STATUS_RUNNING;
	    	        jobFound = true;
	    	    }
		    }
		    // If the job is not running (on this host) see if the
		    // job is complete
		    if(!jobFound) {
		    	status = jobStatusFromScript(outputDir);
		    }

		}
		catch (IOException ex) { ex.printStackTrace(); }
		return status;
	}
	

}
