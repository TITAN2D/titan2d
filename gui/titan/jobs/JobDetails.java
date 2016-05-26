package titan.jobs;

import java.util.Calendar;
import java.lang.Process;

public class JobDetails {

	public String jobName;
	public String jobID;
	public String errorMessage;
	public String outputDirectory;
	public Calendar submissionDate = Calendar.getInstance();
	public Process process;
	
	public void setJobName(String jn) {
		jobName = new String(jn);
	}
	
	public String getJobName() {
		return jobName;
	}
	
	public void setJobID(String jid) {
		jobID = new String(jid);
	}
	
	public String getJobID() {
		return jobID;
	}
	
	public void setProcess(Process pro) {
	    process = pro;
	}
	
	public Process getProcess() {
	    return process;
	}
	
	public void setErrorMessage(String err) {
		errorMessage = new String(err);
	}
	
	public String getErrorMessage() {
		return errorMessage;
	}
	
	public void setSubmissionDate(Calendar c) {
		submissionDate = c;
	}
	
	public Calendar getSubmissionDate() {
		return submissionDate;
	}

	public void setOutputDirectory(String dir) {
		outputDirectory = new String(dir);
	}
	
	public String getOutputDirectory() {
		return outputDirectory;
	}
}
