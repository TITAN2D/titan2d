package titan.jobs.pbs;

import titan.jobs.JobStatus;
import titan.jobs.pbs.QstatProcessor.QstatData;

public class JobStatusPBS extends JobStatus {
	
	public static String getJobStatus(String jobID, String outputDir) {
		String status = STATUS_UNKNOWN;
		QstatProcessor qstat = new QstatProcessor();
		QstatData jobData = qstat.getJob(jobID);
		boolean jobFound = false;
		if(jobData != null) {
			if(jobID.startsWith(jobData.jobID)) {
				if(jobData.status == 'E')
					status = STATUS_EXITING;
				else if(jobData.status == 'H')
					status = STATUS_HELD;
				else if(jobData.status == 'Q')
					status = STATUS_QUEUED;
				else if(jobData.status == 'R')
					status = STATUS_RUNNING;
				else if(jobData.status == 'T')
					status = STATUS_MOVING;
				else if(jobData.status == 'W')
					status = STATUS_WAITING;
				else if(jobData.status == 'S')
					status = STATUS_SUSPENDED;
				jobFound = true;
			}
		}
	    // If the job is not running (on this host) see if the
	    // job is complete
	    if(!jobFound) {
	    	status = jobStatusFromScript(outputDir);
	    }
		return status;
	}
	

	
}
