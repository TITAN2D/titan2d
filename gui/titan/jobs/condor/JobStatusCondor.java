package titan.jobs.condor;

import titan.jobs.JobStatus;
import titan.jobs.condor.CondorJobStatusProcessor.CondorQData;

public class JobStatusCondor extends JobStatus {
	
	public static String getJobStatus(String jobID, String outputDir) {
		String status = STATUS_UNKNOWN;
		
		// Check condor_q for job
		CondorQProcessor condorQ = new CondorQProcessor();
		CondorQData jobData = condorQ.getJob(jobID);
		boolean jobFound = false;
		if(jobData != null) {
			if(jobData.jobID.startsWith(jobID)) {
				status = convertStatus(jobData.status);
				jobFound = true;
			}
		}

		// Check condor_history for job
		if(!jobFound) {
			CondorHistoryProcessor condorHistory = new CondorHistoryProcessor();
			CondorQData jobHistory = condorHistory.getJob(jobID);
			if(jobHistory != null) {
				if(jobHistory.jobID.startsWith(jobID)) {
					status = convertStatus(jobHistory.status);
					jobFound = true;
				}
			}
		}
		
	    // If the job is not running (on this host) see if the
	    // job is complete
	    if(!jobFound) {
	    	status = jobStatusFromScript(outputDir);
	    }
		return status;
	}
	
	private static String convertStatus(char inStatus) {
		String status = STATUS_UNKNOWN;
		if(inStatus == 'U')
			status = STATUS_UNEXPANDED;
		else if(inStatus == 'H')
			status = STATUS_HELD;
		else if(inStatus == 'R')
			status = STATUS_RUNNING;
		else if(inStatus == 'I')
			status = STATUS_QUEUED;
		else if(inStatus == 'C')
			status = STATUS_COMPLETED;
		else if(inStatus == 'X')
			status = STATUS_REMOVED;
		
		return status;
	}
}
