package titan.jobs.local;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

import titan.gui.TabRunParameters;

import javax.swing.JOptionPane;

import titan.gui.TitanConstants;
import titan.jobs.*;

public class RunMethodLocal extends RunMethod {

    public RunMethodLocal(String srestart, String dir, String date) {
        super(srestart, dir, date, "runtitan.sh");
    }

    public boolean allowMultiProcessors() {
        return true;
    }

    public void initialize(GisParms gParms,
                           boolean useExternal, String binaryDir,
                           String preCommands,
                           JobSubmissionParameters jobParms,
                           String srestartFile, int irestartMaxNumberTimeSteps, float frestartMaxTime) {
        super.initialize(gParms, useExternal, binaryDir,
                preCommands, jobParms, srestartFile, irestartMaxNumberTimeSteps, frestartMaxTime);

        createScript(gisParms, jobParms, preCommands);
    }

    protected void createScript(String[] preprocParms,
                                JobSubmissionParameters jobParms,
                                String pre) {

        try {
            writer = new BufferedWriter(
                    new FileWriter(basedir + File.separator + execName));

            try {

                writer.write("#!/bin/bash\n");
                writer.write("\n");

                // Change directory into directory to run in
                writer.write("cd " + basedir + "\n");
                writer.write("\n");

                writer.write("# User specified commands\n");
                writer.write(pre);
                writer.write("\n");

                // titan
                // this is where we need to run with the different run modes.
                // openMP mode
                // source <path_to_titan>/bin/titanvars.sh
                //<path_to_titan>/bin/titan -nt <number of threads> <input python script>

                // In order to run Titan in hybrid MPI/OpenMP mode:
                // source <path_to_titan>/bin/titanvars.sh
                // mpirun -np <number of processes> <path_to_titan>/bin/titan -nt <number of threads> <input python script>

                Integer numNodes = Integer.parseInt(jobParms.getJobParameter(JobSubmissionParameters.NODES));
                Integer numCPUs = Integer.parseInt(jobParms.getJobParameter(JobSubmissionParameters.CPUS));
                Integer numTotalCPUs = numNodes * numCPUs;

                writer.write("# Run titan\n");

                String mpiArg = new String("");
                String openmpArg = new String("");
                String otherArgs = new String("");

                // For Hybrid Mode
                if (jobParms.getJobParameter(JobSubmissionParameters.RUN_MODE).compareTo(JobSubmissionContainer.RUN_MODE_HYBRID) == 0) {
                    mpiArg = "mpirun -np " + numTotalCPUs.toString() + " ";
                }
                writer.write(mpiArg);

                if (binMode == true)
                    writer.write(binDir + File.separator + "titan ");
                    // Depend on path env variable
                else
                    writer.write("titan ");

                // For both Hybrid and Openmp Modes
                if ((jobParms.getJobParameter(JobSubmissionParameters.RUN_MODE).compareTo(JobSubmissionContainer.RUN_MODE_HYBRID) == 0)
                        || (jobParms.getJobParameter(JobSubmissionParameters.RUN_MODE).compareTo(JobSubmissionContainer.RUN_MODE_OPENMP) == 0)) {
                    openmpArg = "-nt " + numCPUs.toString() + " ";
                }
                writer.write(openmpArg);

                if (restartEnabled.compareTo(TitanConstants.TRUE) == 0) {

                    // Restart is enabled

                    otherArgs = "-restart ";
                    if (restartMaxNumberTimeSteps != 0) {
                        otherArgs = otherArgs + "-add-iter " + restartMaxNumberTimeSteps + " ";
                    }
                    if (restartMaxTime != 0.0f) {
                        otherArgs = otherArgs + "-add-time " + restartMaxTime + " ";
                    }
                    otherArgs = otherArgs + restartFile + " ";

                } else {

                    // Restart is not enabled
                    otherArgs = TitanRunInput.SIMULATION_DATA.concat(" " + jobParms.getJobParameter(JobSubmissionParameters.RUN_ARGS) + " ");

                }
                writer.write(otherArgs);

                // $? returns the results of the last executed command
                writer.write("1>> " + STDOUT + " 2>> " + STDERR + "\n");
                writer.write("\n");
                writer.write("if [ \"$?\" -eq 0 ]\n");
                writer.write("then\n");
                writer.write("  touch " + TITAN_SUCCESSFUL + "\n");
                writer.write("else\n");
                writer.write("  touch " + TITAN_FAILED + "\n");
                writer.write("fi\n");
                writer.write("\n");

            } catch (IOException ex) {
                ex.printStackTrace();
            }

            writer.close();

        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    /**
     * Submit a job to the PBS job manager.
     *
     * @return Returns information regarding the job submitted to the
     * PBS queue.
     */
    public JobDetails runJob() {

        JobDetails jobDetails = new JobDetails();
        jobDetails.setJobName(jobName);
        jobDetails.setOutputDirectory(basedir);

        try {
            File execFile = new File(basedir + File.separator + execName);
            execFile.setExecutable(true, false);
            Process jobProcess = Runtime.getRuntime().exec(basedir + File.separator + execName);
            jobDetails.setProcess(jobProcess);
            Process psJob = Runtime.getRuntime().exec("ps -aux");

            BufferedReader rs = new BufferedReader(new InputStreamReader(psJob.getInputStream()));

            String line;
            boolean jobFound = false;
            while ((line = rs.readLine()) != null) {
                if (line.contains(basedir)) {
                    line = line.substring(line.indexOf(" ") + 1);
                    while (line.startsWith(" ")) line = line.substring(1);
                    line = line.substring(0, line.indexOf(" "));
                    jobDetails.setJobID(line);
                    jobFound = true;
                }
            }

            if (jobFound == false) {
                jobDetails.setErrorMessage("Error submitting job");
                jobDetails.setJobID("JOB NOT SUBMITTED");
            }

        } catch (IOException ex) {
            ex.printStackTrace();
        }
        return jobDetails;
    }

    /**
     * @param jobStatus
     * @param jobID
     * @return Returns true if the job was removed successfully and false
     * otherwise.
     */
    public static boolean removeJob(String jobStatus, String jobID) {
        boolean status = false;
        try {
            if ((jobStatus.compareTo(JobStatus.STATUS_RUNNING) == 0) ||
                    (jobStatus.compareTo(JobStatus.STATUS_HELD) == 0) ||
                    (jobStatus.compareTo(JobStatus.STATUS_QUEUED) == 0) ||
                    (jobStatus.compareTo(JobStatus.STATUS_MOVING) == 0) ||
                    (jobStatus.compareTo(JobStatus.STATUS_WAITING) == 0) ||
                    (jobStatus.compareTo(JobStatus.STATUS_SUSPENDED) == 0) ||
                    (jobStatus.compareTo(JobStatus.STATUS_UNEXPANDED) == 0)) {
                Runtime.getRuntime().exec("kill -9 " + jobID);
            } else
                status = false;

        } catch (IOException ex) {
            ex.printStackTrace();
        }
        return status;
    }

}
