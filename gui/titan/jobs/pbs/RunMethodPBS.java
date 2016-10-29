package titan.jobs.pbs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

import titan.jobs.JobDetails;
import titan.jobs.JobStatus;
import titan.jobs.JobSubmissionParameters;
import titan.jobs.RunMethod;

public class RunMethodPBS extends RunMethod {
	
	public RunMethodPBS(String srunType, String dir, String date) {
		super(srunType, dir, date, "runtitan.sh");
	}

	public boolean allowMultiProcessors() { return true; }

	public void initialize(GisParms gParms,
						   boolean useExternal, String binaryDir,
						   String preCommands,
						   JobSubmissionParameters jobParms,
						   String srestartFile, int irestartMaxNumberTimeSteps, float frestartMaxTime) {
		super.initialize(gParms, useExternal, binaryDir,
				preCommands, jobParms, srestartFile, irestartMaxNumberTimeSteps, frestartMaxTime);

		//createScript(preprocParms, vecParms, jobParms, preCommands);
	}

	protected void createScript(String[] preprocParms, String[] vecParms,
			                    JobSubmissionParameters jobParms,
			                    String pre) {

		try {
			writer = new BufferedWriter(
					new FileWriter(basedir + File.separator + execName));

			writer.write("#!/bin/bash\n");
			// Name the job
			writer.write("#PBS -N " + jobName + "\n");
			// wall time
			if(jobParms.getJobParameter(JobSubmissionParameters.MAX_RUN_TIME).compareTo("") != 0)
				writer.write("#PBS -l walltime=" + jobParms.getJobParameter(JobSubmissionParameters.MAX_RUN_TIME) + "\n");
			// Number of processors
			writer.write("#PBS -l nodes=" + jobParms.getJobParameter(JobSubmissionParameters.NODES) +
					     ":ppn=" + jobParms.getJobParameter(JobSubmissionParameters.CPUS) + "\n");
			// file for standard output
			writer.write("#PBS -o " + basedir + File.separator + STDOUT + "\n");
			// file for standard error
			writer.write("#PBS -e " + basedir + File.separator + STDERR + "\n");
			// Queue name
			if(jobParms.getJobParameter(JobSubmissionParameters.QUEUE).compareTo("") != 0)
				writer.write("#PBS -q " + jobParms.getJobParameter(JobSubmissionParameters.QUEUE) + "\n");
			// export the users environment for the job
			writer.write("#PBS -V\n");
			writer.write("\n");
			
//			super.createScript(preprocParms, vecParms, jobParms, pre);
			
			// Change directory into directory to run in
			writer.write("cd " + basedir + "\n");
			writer.write("\n");
			
			writer.write("# User specified commands\n");
			writer.write(pre);
			writer.write("\n");

			/*
			// titan_preprocess
			writer.write("\n# Run the titan preprocessor\n");
			
			String preprocess = null;
			// Full path to exectuable
			if(binMode == true)
			    preprocess = new String(binDir + File.separator + "titan_preprocess ");
			// Depend on path env variable
			else
				preprocess = new String("titan_preprocess ");
			
			for(int i = 0; i < preprocParms.length; i++) {
				preprocess = preprocess.concat(" " + preprocParms[i]);
			}
			writer.write(preprocess + "\n");
			writer.write("\n");
			writer.write("if [ \"$?\" -eq 0 ]\n");
			writer.write("then\n");
			writer.write("  touch " + PREPROCESS_SUCCESSFUL + "\n");
			writer.write("else\n");
			writer.write("  touch " + PREPROCESS_FAILED + "\n");
			writer.write("fi\n");
			writer.write("\n");
			
			// VecDataPreprocessor
			if(vecParms.length > 0) {
				writer.write("# Run Vec Data Preprocessor\n");
				
				String vecdata = null;
				// Full path to exectuable
				if(binMode == true)
				    vecdata = new String(binDir + File.separator + "VecDataPreproc ");
				// Depend on path env variable
				else
					vecdata = new String("VecDataPreproc ");
				
				for(int i = 0; i < vecParms.length; i++) {
					vecdata = vecdata.concat(" " + vecParms[i]);
				}
				writer.write(vecdata + "\n");
				writer.write("\n");
				writer.write("if [ \"$?\" -eq 0 ]\n");
				writer.write("then\n");
				writer.write("  touch " + VECDATA_SUCCESSFUL + "\n");
				writer.write("else\n");
				writer.write("  touch " + VECDATA_FAILED + "\n");
				writer.write("fi\n");
				writer.write("\n");
			}
			*/

			// titan
			writer.write("# Run titan\n");
			String mpiCommand = new String("");
			//if(jobParms.getJobParameter(JobSubmissionParameters.MPI_STYLE) != null) {
				//mpiCommand = mpiCommand.concat(jobParms.getJobParameter(JobSubmissionParameters.MPI_STYLE) + " ");
			    //mpiCommand = mpiCommand.concat(" " + jobParms.getJobParameter(JobSubmissionParameters.MPI_ARGS) + " ");
			//}

			// Full path to exectuable
			if(binMode == true)
			    writer.write(mpiCommand + " " + binDir + File.separator + "titan");
			// Depend on path env variable
			else
				writer.write(mpiCommand + " titan");
			writer.write("\n");
			writer.write("if [ \"$?\" -eq 0 ]\n");
			writer.write("then\n");
			writer.write("  touch " + TITAN_SUCCESSFUL + "\n");
			writer.write("else\n");
			writer.write("  touch " + TITAN_FAILED + "\n");
			writer.write("fi\n");
			writer.write("\n");
			
			writer.close();

		}
		catch (IOException ex) { ex.printStackTrace(); }
	}
	
	/**
	 * Submit a job to the PBS job manager.
	 * @return Returns information regarding the job submitted to the
	 *         PBS queue.
	 */
	public JobDetails runJob() {

		JobDetails jobDetails = new JobDetails();
		jobDetails.setJobName(jobName);
		jobDetails.setOutputDirectory(basedir);
		
		try {
		    Process titanJob = Runtime.getRuntime().exec("qsub " + basedir + File.separator + execName);
		    BufferedReader err = new BufferedReader(new InputStreamReader(titanJob.getErrorStream()));
		    String line;
		    if((line = err.readLine()) == null) {
			    BufferedReader out = new BufferedReader(new InputStreamReader(titanJob.getInputStream()));
			    line = out.readLine();
		    	jobDetails.setJobID(line);
		    }
		    else {
		    	jobDetails.setErrorMessage(line);
		    	jobDetails.setJobID("JOB NOT SUBMITTED");
		    }
		}
		catch (IOException ex) { ex.printStackTrace(); }

		return jobDetails;
	}
	
	/**
	 * 
	 * @param jobStatus
	 * @param jobID
	 * @return Returns true if the job was removed successfully and false
	 *         otherwise.
	 */
	public static boolean removeJob(String jobStatus, String jobID) {

		boolean status = true;
		try {
			if((jobStatus.compareTo(JobStatus.STATUS_RUNNING) == 0)    ||
				(jobStatus.compareTo(JobStatus.STATUS_HELD) == 0)      ||
				(jobStatus.compareTo(JobStatus.STATUS_QUEUED) == 0)    ||
				(jobStatus.compareTo(JobStatus.STATUS_MOVING) == 0)    ||
				(jobStatus.compareTo(JobStatus.STATUS_WAITING) == 0)   ||
				(jobStatus.compareTo(JobStatus.STATUS_SUSPENDED) == 0) ||
				(jobStatus.compareTo(JobStatus.STATUS_UNEXPANDED) == 0)) {
		        Runtime.getRuntime().exec("qdel " + jobID);
			}
			else
				status = false;

		}
		catch (IOException ex) { ex.printStackTrace(); }
		return status;
	}
}
