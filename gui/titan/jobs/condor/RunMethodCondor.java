package titan.jobs.condor;

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

public class RunMethodCondor extends RunMethod {

	public static final String CONDOR_SCRIPT = "runtitan.condor";
	
	public RunMethodCondor(String srunType, String dir, String date) {
		super(srunType, dir, date, "runtitan.sh");
	}

	public boolean allowMultiProcessors() { return false; }

	public void initialize(GisParms gParms,
						   boolean useExternal, String binaryDir,
						   String preCommands,
						   JobSubmissionParameters jobParms,
						   String srestartFile, int irestartMaxNumberTimeSteps, float frestartMaxTime) {
		super.initialize(gParms, useExternal, binaryDir,
				preCommands, jobParms, srestartFile, irestartMaxNumberTimeSteps, frestartMaxTime);

		this.zipGISData(gisParms[2]);
		gisParms[2] = gisParms[2].substring(gisParms[2].lastIndexOf("/")+1);
		//createScript(preprocParms, jobParms, preCommands);
		createCondorScript(jobParms.getJobParameter(JobSubmissionParameters.REQUIREMENTS));
	}
	
	/**
	 * Submit a job to the Condor job manager.
	 * @return Returns information regarding the job submitted to the
	 *         Condor queue.
	 */
	public JobDetails runJob() {
		JobDetails jobDetails = new JobDetails();
		jobDetails.setJobName(jobName);
		jobDetails.setOutputDirectory(basedir);

		try {
		    Process titanJob = Runtime.getRuntime().exec("condor_submit " + basedir + File.separator + CONDOR_SCRIPT);
		    BufferedReader err = new BufferedReader(new InputStreamReader(titanJob.getErrorStream()));
		    String line;
		    if((line = err.readLine()) == null) {
			    BufferedReader out = new BufferedReader(new InputStreamReader(titanJob.getInputStream()));
			    while((line = out.readLine()) != null) {
			    	if(line.contains("submitted to cluster")) {
			    		line = line.substring(line.indexOf("submitted to cluster ") + "submitted to cluster ".length());
				    	jobDetails.setJobID(line);
			    	}
			    }
		    }
		    else {
		    	String errorMsg = new String();
		    	while((line = err.readLine()) != null) {
		    		errorMsg = errorMsg.concat(line);
		    	}
		    	jobDetails.setErrorMessage(errorMsg);
		    	jobDetails.setJobID("JOB NOT SUBMITTED");
		    }
		}
		catch (IOException ex) { ex.printStackTrace(); }
		
		return jobDetails;
	}
	
	/**
	 * Create the runtitan.sh bash shell script used to launch
	 * titan.
	 */
	protected void createScript(String[] preprocParms, String[] vecParms,
            JobSubmissionParameters jobParms, String pre) {

		try {
		    writer = new BufferedWriter(
		    new FileWriter(basedir + File.separator + execName));
		
			try {
			
				writer.write("#!/bin/bash\n");
				writer.newLine();
				
				writer.write("# User specified commands\n");
				writer.write(pre);
				writer.newLine();
				
				// tar and zip GIS data directory structure
				writer.write("gunzip gis_data.tar.gz\n");
				writer.write("tar xvf gis_data.tar\n");
				
				// titan_preprocess
				/*
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
//				writer.write(preprocess + "1> " + STDOUT + " 2> " + STDERR + "\n");
				writer.write(preprocess + "\n");
				writer.newLine();
				writer.write("if [ \"$?\" -eq 0 ]\n");
				writer.write("then\n");
				writer.write("  touch " + PREPROCESS_SUCCESSFUL + "\n");
				writer.write("else\n");
				writer.write("  touch " + PREPROCESS_FAILED + "\n");
				writer.write("fi\n");
				writer.newLine();
				
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
//				    writer.write(vecdata + "1>> " + STDOUT + " 2>>" + STDERR + "\n");
				    writer.write(vecdata + "\n");
					writer.newLine();
					writer.write("if [ \"$?\" -eq 0 ]\n");
					writer.write("then\n");
					writer.write("  touch " + VECDATA_SUCCESSFUL + "\n");
					writer.write("else\n");
					writer.write("  touch " + VECDATA_FAILED + "\n");
					writer.write("fi\n");
					writer.newLine();
				}
				*/

				// titan
				writer.write("# Run titan\n");
				// Full path to exectuable
				if(binMode == true)
				    writer.write(binDir + File.separator + "titan ");
				// Depend on path env variable
				else
				    writer.write("titan ");
				
//				writer.write("titan 1>> " + STDOUT + " 2>>" + STDERR + "\n");
				writer.newLine();
				writer.write("if [ \"$?\" -eq 0 ]\n");
				writer.write("then\n");
				writer.write("  touch " + TITAN_SUCCESSFUL + "\n");
				writer.write("else\n");
				writer.write("  touch " + TITAN_FAILED + "\n");
				writer.write("fi\n");
				writer.newLine();
			
				writer.write("# Cleanup\n");
				writer.write("rm " + this.GIS_DATA_TAR);
				
			}
			catch (IOException ex) { ex.printStackTrace(); }
		
		writer.close();
		
		}
		catch (IOException ex) { ex.printStackTrace(); }
	}
	
	/**
	 * Create the condor description file used to launch a condor run.
	 */
	protected void createCondorScript(String requirements) {
		try {
			writer = new BufferedWriter(
					new FileWriter(basedir + File.separator + CONDOR_SCRIPT));

			try {
				writer.write("Universe = vanilla\n");
				writer.write("Executable = " + basedir + File.separator + execName + "\n");
				writer.newLine();

				if(requirements.compareTo("") != 0) {
				    writer.write("Requirements = " + requirements + "\n");
				    writer.newLine();
				}
				
				writer.write("input = /dev/null\n");
				writer.write("output = " + STDOUT + "\n");
				writer.write("error = " + STDERR + "\n");
				writer.newLine();

				writer.write("initialdir = " + basedir + "\n");
				writer.newLine();
				
				writer.write("transfer_input_files = ./simulation.data, ./scale.data, " +
						"./frict.data, " + "./" + GIS_DATA_TAR + ".gz\n");
				writer.newLine();
				
				writer.write("getenv = true\n");
				writer.newLine();
				
				writer.write("should_transfer_files = YES\n");
				writer.write("when_to_transfer_output = ON_EXIT\n");
				writer.newLine();
				
				writer.write("notification = Never\n");
				writer.write("Queue\n");

				
			}
			catch (IOException ex) { ex.printStackTrace(); }
			
			writer.close();

		}
		catch (IOException ex) { ex.printStackTrace(); }
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
		        Runtime.getRuntime().exec("condor_rm " + jobID);
			}
			else
				status = false;

		}
		catch (IOException ex) { ex.printStackTrace(); }
		return status;
	}
	
}
