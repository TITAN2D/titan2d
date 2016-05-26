package titan.jobs;

import titan.jobs.condor.JobStatusCondor;
import titan.jobs.condor.RunMethodCondor;
import titan.jobs.local.JobStatusLocal;
import titan.jobs.local.RunMethodLocal;
import titan.jobs.pbs.JobStatusPBS;
import titan.jobs.pbs.RunMethodPBS;
import titan.jobs.submit.JobStatusSubmit;
import titan.jobs.submit.RunMethodSubmit;
import titan.jobs.RunMethod;

public class JobSubmissionContainer {

	// Reference http://stackoverflow.com/questions/25287981/mpiexec-vs-mpirun:
	// mpiexec is defined in the MPI standard (well, the recent versions at least)
	// and I refer you to those (your favorite search engine will find them for you) for details.
	// mpirun is a command implemented by many MPI implementations.
	// It has never, however, been standardized and there have always been, often subtle, differences between implementations.
	// For details see the documentation of the implementation(s) of your choice.
	// And yes, they are both used to launch MPI programs, these days mpiexec is generally preferable because it is standardised.
	public RunMethod runMethod = null;

	// These would go away when reflection compete
	public static final String RUN_STYLE_LOCAL = "Local";
	public static final String RUN_STYLE_PBS = "PBS";
	public static final String RUN_STYLE_CONDOR = "Condor";
	public static final String RUN_STYLE_SUBMIT = "Hub-Submit";
	private int numMethods = 4;

	//public static final String MPI_STYLE_MPIEXEC = "mpiexec";
	//public static final String MPI_STYLE_MPIRUN  = "mpirun";
	//public static final String MPI_STYLE_MPINONE = "Do Not Use MPI";
	public static final String RUN_MODE_OPENMP = "openmp";
	public static final String RUN_MODE_HYBRID = "hybrid";
	private int numRunModes = 2;
	
	public static final int UNKNOWN_INDEX = -1;
	private String methods[];
	private String runModes[];

	// Write data from the HUB submit output strean
	public static final String HUB_SUBMIT_INFO = new String("HUB_submit_info.txt");

	// Really want to use reflection, but this will do for now
	public JobSubmissionContainer() {
        methods = new String[numMethods];
        methods[0] = new String(RUN_STYLE_LOCAL);
        methods[1] = new String(RUN_STYLE_PBS);
        methods[2] = new String(RUN_STYLE_CONDOR);
        methods[3] = new String(RUN_STYLE_SUBMIT);
        
        runModes = new String[numRunModes];
        runModes[0] = new String(RUN_MODE_OPENMP);
        runModes[1] = new String(RUN_MODE_HYBRID);
 	}
	
	public String[] getMethodNames() {
		return methods;
	}

	public String[] getRunModes() {
		return runModes;
	}
	
	// This method needs to be re-written for reflection
	public void createRunMethod(String m, String srestart, String parmDir, String dateString) {
		if(m.compareTo(RUN_STYLE_LOCAL) == 0) runMethod = new RunMethodLocal(srestart, parmDir, dateString);
		else if(m.compareTo(RUN_STYLE_PBS) == 0) runMethod = new RunMethodPBS(srestart, parmDir, dateString);
		else if(m.compareTo(RUN_STYLE_CONDOR) == 0) runMethod = new RunMethodCondor(srestart, parmDir, dateString);
		else if(m.compareTo(RUN_STYLE_SUBMIT) == 0) runMethod = new RunMethodSubmit(srestart, parmDir, dateString);
	}

	// This method may need to be re-written for reflection
	public RunMethod getRunMethod() {
		return runMethod;
	}

	// This method needs to be re-written for reflection
	public String getJobStatus(String runStyle, String jobID, String outDir) {
		String status = null;
		if(runStyle.compareTo(RUN_STYLE_LOCAL) == 0)
			status = JobStatusLocal.getJobStatus(jobID, outDir);
		else if(runStyle.compareTo(RUN_STYLE_PBS) == 0)
			status = JobStatusPBS.getJobStatus(jobID, outDir);
		else if(runStyle.compareTo(RUN_STYLE_CONDOR) == 0)
			status = JobStatusCondor.getJobStatus(jobID, outDir);
		else if(runStyle.compareTo(RUN_STYLE_SUBMIT) == 0)
			status = JobStatusSubmit.getJobStatus(jobID, outDir);
		else
		    status = "Error: run style unknown";
		return status;
	}
	
	// This method needs to be re-written for reflection
	public boolean deleteJob(String runStyle, String jobID, String outDir) {
		String status = getJobStatus(runStyle, jobID, outDir);
		if(runStyle.compareTo(RUN_STYLE_LOCAL) == 0)
			return RunMethodLocal.removeJob(status, jobID);
		else if(runStyle.compareTo(RUN_STYLE_PBS) == 0)
			return RunMethodPBS.removeJob(status, jobID);
		else if(runStyle.compareTo(RUN_STYLE_CONDOR) == 0)
			return RunMethodCondor.removeJob(status, jobID);
		else if(runStyle.compareTo(RUN_STYLE_SUBMIT) == 0)
			return RunMethodSubmit.removeJob(status, jobID);
		else
			return false;
	}
	
}
