package titan.jobs.submit;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

import javax.swing.JOptionPane;

import titan.gui.TitanConstants;
import titan.jobs.*;

public class RunMethodSubmit extends RunMethod {

    public String db_format;
    public String db_path;
    public String db_dir;
    public String h5path = "";
    public String h5name = "";

    public int use_irods = 0;  // 0 - dont use irods for import of DEM
    // 1 - use irods for import of DEM
    public int output_to_irods = 0;  // 0 - dont output files to irods
    // 1 - use irods for output
    public int user_commands = 0; // 0 - no user commands are included
    // 1 - user commands are included


    public RunMethodSubmit(String srestart, String dir, String date) {
        super(srestart, dir, date, "submittitan.sh");
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

        // See RunMethod.java for indices
        int preprocParms_GIS_FORMAT_index = 2;
        int preprocParms_GIS_INFO_DIRECTORY_index = 3;
        int preprocParms_GIS_MAP_index = 6;

	    /*
        // Test if it should grab DEM from ~/irods directory
        // Note: Per legacy code, use_irods set to 0.
        db_path = preprocParms[preprocParms_GIS_INFO_DIRECTORY_index];
        if ( db_path.indexOf("/irods/") != -1 ) {
            //use irods
            use_irods = 1;
        }
        else {
            //do not use irods
            use_irods = 0;
        }
        */
        use_irods = 0; //This will just load from the mount before sending over to the hub.

        //script to run submit

        int index;

        // Variables needed to zip the GIS directory
        db_format = preprocParms[preprocParms_GIS_FORMAT_index];

        if (db_format.compareTo(TitanConstants.GisFormats[TitanConstants.GIS_FORMAT_GIS_GRASS]) == 0) {

            // GIS_GRASS
            // For HUB-Submits, the information main directory is zipped into a tar.gz file

            index = preprocParms[preprocParms_GIS_INFO_DIRECTORY_index].lastIndexOf(File.separator);
            db_path = preprocParms[preprocParms_GIS_INFO_DIRECTORY_index].substring(0, index);
            db_dir = preprocParms[preprocParms_GIS_INFO_DIRECTORY_index].substring(index + 1);

        } else {

            // GDAL
            // For GDAL, GIS Map is a full path name to the mapset file.
            // For HUB-Submits, the mapset file is zipped into a tar.gz file.
            // In this case, db_dir is actually a db_file

            index = preprocParms[preprocParms_GIS_MAP_index].lastIndexOf(File.separator);
            db_path = preprocParms[preprocParms_GIS_MAP_index].substring(0, index);
            db_dir = preprocParms[preprocParms_GIS_MAP_index].substring(index + 1);
        }

        //System.out.println("db_path: " + db_path);
        //System.out.println("db_dir: " + db_dir);

        if (restartEnabled.compareTo(TitanConstants.TRUE) == 0) {

            // Restart is enabled

            // Variables needed to zip the restart .h5 file
            index = restartFile.lastIndexOf(File.separator);
            h5path = restartFile.substring(0, index);
            h5name = restartFile.substring(index + 1);
            System.out.println(": " + h5path);
            System.out.println(": " + h5name);
        }

        //testing
        if (jobParms.getJobParameter(JobSubmissionParameters.OUTPUT_IRODS).equalsIgnoreCase("true"))
            output_to_irods = 1;
        else
            output_to_irods = 0;

        try {
            writer = new BufferedWriter(new FileWriter(basedir + File.separator + "argfile"));
            try {

                /*
                In order to run Titan in openMP mode (preferred at this time):

                source <path_to_titan>/bin/titanvars.sh
                <path_to_titan>/bin/titan -nt <number of threads> <input python script>

                In order to run Titan in hybrid MPI/OpenMP mode:

                source <path_to_titan>/bin/titanvars.sh
                // Note: Need to determine what user input is required for allocated the processes and threads
                mpirun -np <number of processes> <path_to_titan>/bin/titan -nt <number of threads> <input python script>"
                */

                Integer numNodes = Integer.parseInt(jobParms.getJobParameter(JobSubmissionParameters.NODES_SUBMIT));
                Integer numCPUs = Integer.parseInt(jobParms.getJobParameter(JobSubmissionParameters.CPUS_SUBMIT));
                Integer numTotalCPUs = numNodes * numCPUs;

                // For initial testing of titan
                //String runCommand = new String("titan -nt " + numCPUs.toString() + " simulation.py");
                //writer.write(runCommand);

                String mpiArg = new String("");
                String openmpArg = new String("");
                String otherArgs = new String("");

                // For Hybrid Mode
                if (jobParms.getJobParameter(JobSubmissionParameters.RUN_MODE).compareTo(JobSubmissionContainer.RUN_MODE_HYBRID) == 0) {
                    mpiArg = "-np " + numTotalCPUs.toString();
                }
                writer.write(mpiArg);
                writer.write("\n");

                // For both Hybrid and Openmp Modes
                if ((jobParms.getJobParameter(JobSubmissionParameters.RUN_MODE).compareTo(JobSubmissionContainer.RUN_MODE_HYBRID) == 0)
                        || (jobParms.getJobParameter(JobSubmissionParameters.RUN_MODE).compareTo(JobSubmissionContainer.RUN_MODE_OPENMP) == 0)) {
                    openmpArg = "-nt " + numCPUs.toString();
                }
                writer.write(openmpArg);
                writer.write("\n");

                if (restartEnabled.compareTo(TitanConstants.TRUE) == 0) {

                    // Restart is enabled

                    otherArgs = "-restart ";
                    if (restartMaxNumberTimeSteps != 0) {
                        otherArgs = otherArgs + "-add-iter " + restartMaxNumberTimeSteps + " ";
                    }
                    if (restartMaxTime != 0.0f) {
                        otherArgs = otherArgs + "-add-time " + restartMaxTime + " ";
                    }

                    // Need just the restart .h5 filename here
                    otherArgs = otherArgs + h5name + " ";

                } else {

                    // Restart is not enabled

                    otherArgs = TitanRunInput.SIMULATION_DATA.concat(" " + jobParms.getJobParameter(JobSubmissionParameters.RUN_ARGS) + " ");
                }
                writer.write(otherArgs);

            } catch (IOException ex) {
                ex.printStackTrace();
            }
            writer.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }

        try {
            BufferedWriter writer2 = new BufferedWriter(new FileWriter(basedir + File.separator + execName));
            try {
                writer2.write("#!/bin/bash\n");
                writer2.write("\n");

                //change to correct directory
                writer2.write("cd " + basedir);
                writer2.write("\n\n");

                //tar GIS database for normal run or restart file for restart run directory
                if (restartEnabled.compareTo(TitanConstants.TRUE) == 0) {

                    // Restart is enabled

                    if (use_irods == 0) {

                        // tar GIS directory and restart file
                        writer2.write("tar -zcf db.tar.gz -C " + db_path + " " + db_dir + " -C " + h5path + " " + h5name);
                    }
                } else {

                    // Restart is not enabled

                    if (use_irods == 0) {

                        // tar GIS directory
                        writer2.write("tar -zcf db.tar.gz -C " + db_path + " " + db_dir);
                    }
                }
                writer2.write("\n\n");

                //tar everything in directory
                writer2.write("tar -zcf submit.tar.gz --exclude runtitan.sh *");
                writer2.write("\n\n");

                //submit
                if ((jobParms.getJobParameter(JobSubmissionParameters.NODES_SUBMIT).compareTo("1") == 0) &&
                        (jobParms.getJobParameter(JobSubmissionParameters.CPUS_SUBMIT).compareTo("1") == 0)) {
                    writer2.write("submit -w ");
                    writer2.write(jobParms.getJobParameter(JobSubmissionParameters.MAX_RUN_TIME_SUBMIT));
                    if (output_to_irods == 1) {
                        writer2.write(" -q ");
                    }
                    writer2.write(" -v u2-grid submitTitanRun ./argfile ");
                    writer2.write(Integer.toString(output_to_irods));
                    writer2.write(" -i submit.tar.gz");
                    if (use_irods == 1) {
                        writer2.write(" -i preprocess_irods.sh");
                    }
                    if (output_to_irods == 1) {
                        writer2.write(" -i postprocess_irods.sh");
                    }
                } else {
                    Integer numNodes = Integer.parseInt(jobParms.getJobParameter(JobSubmissionParameters.NODES_SUBMIT));
                    Integer numCPUs = Integer.parseInt(jobParms.getJobParameter(JobSubmissionParameters.CPUS_SUBMIT));
                    Integer numTotalCPUs = numNodes * numCPUs;
                    writer2.write("submit -n ");
                    //writer2.write(jobParms.getJobParameter(JobSubmissionParameters.NODES_SUBMIT));
                    writer2.write(numTotalCPUs.toString());
                    writer2.write(" -N ");
                    writer2.write(jobParms.getJobParameter(JobSubmissionParameters.CPUS_SUBMIT));
                    writer2.write(" -w ");
                    writer2.write(jobParms.getJobParameter(JobSubmissionParameters.MAX_RUN_TIME_SUBMIT));
                    if (output_to_irods == 1) {
                        writer2.write(" -q ");
                    }
                    writer2.write(" -v u2-grid submitTitanRun ./argfile 0 -i submit.tar.gz");
                    if (use_irods == 1) {
                        writer2.write(" -i preprocess_irods.sh");
                    }
                    if (output_to_irods == 1) {
                        writer2.write(" -i postprocess_irods.sh");
                    }
                }
                writer2.write("\n\n");

                //remove tarballs when done - keep for initial testing
                if (use_irods == 0) {
                    writer2.write("rm db.tar.gz submit.tar.gz");
                } else {
                    writer2.write("rm submit.tar.gz");
                }

            } catch (IOException ex) {
                ex.printStackTrace();
            }
            writer2.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }

        //preprocess script if using irods.
        // Note: Per legacy code, use_irods set to 0.
        // If enabled, need to modify following for GDAL
        if (use_irods == 1) {
            try {
                BufferedWriter writer3 = new BufferedWriter(new FileWriter(basedir + File.separator + "preprocess_irods.sh"));
                try {
                    writer3.write("#!/bin/bash\n");
                    writer3.write("\n");
                    //int index = preprocParms[i].lastIndexOf(File.separator);
                    //String db_path = preprocParms[i].substring(0,index);
                    db_path = preprocParms[preprocParms_GIS_INFO_DIRECTORY_index];
                    writer3.write("cp -r ~" + db_path.substring(db_path.indexOf("/irods")) + " .\n");
                } catch (IOException ex) {
                    ex.printStackTrace();
                }
                writer3.close();
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
        // end of creating preprocess script

        //postprocess script if outputing to irods
        if (output_to_irods == 1) {
            try {
                BufferedWriter writer4 = new BufferedWriter(new FileWriter(basedir + File.separator + "postprocess_irods.sh"));
                try {
                    //System.out.println( "DATE: " + jobParms.getJobParameter(JobSubmissionParameters.DATE) );
                    //System.out.println( "SAVE NAME: " + jobParms.getJobParameter(JobSubmissionParameters.SAVE_DIR) );
                    //System.out.println( "USERNAME: " + java.lang.System.getProperty("user.name") );
                    String output_dir;
                    String date = jobParms.getJobParameter(JobSubmissionParameters.DATE);
                    String save_dir = jobParms.getJobParameter(JobSubmissionParameters.SAVE_DIR);
                    String username = java.lang.System.getProperty("user.name");
                    writer4.write("#!/bin/bash\n");

                    writer4.write("pwd\n\n");
                    writer4.write("ls -lart *\n\n");

                    //specific output directory
                    if (save_dir.isEmpty()) {
                        //no save dir specified
                        //writer4.write("imkdir -pV /vhub/home/public/TitanVhubResults/" + date + "\n\n");
                        output_dir = "/vhub/home/public/TitanVhubResults/" + date;
                    } else {
                        //writer4.write("imkdir -pV /vhub/home/public/TitanVhubResults/" + save_dir + File.separator + date + "\n\n");
                        output_dir = "/vhub/home/public/TitanVhubResults/" + save_dir + File.separator + date;
                    }

                    writer4.write("imkdir -pV " + output_dir + "\n\n");

                    writer4.write("#remove DEM\n");
                    //writer4.write("rm -r grass5\n\n");
                    //writer4.write("rm -r " + db_path.substring(db_path.lastIndexOf(File.separator) + 1) + "\n\n");
                    writer4.write("rm -r " + db_dir + "\n\n");

                    writer4.write("#put all output files/directories into irods\n");
                    writer4.write("iput -rfV * " + output_dir + "\n\n");

                    //!!!!!!! NEED TO CHANGE THIS SO THAT THE FOLLOWING FILES CAN STILL BE SEEN
                    // Stdout.txt, Stderr.txt, perform*, *.stdout and *.stderr
                    writer4.write("#remove the local files/directories\n");
                    //writer4.write("rm *\n");
                    // grep -v: invert the sense of matching to select non-matching lines;
                    // the '^' (begin) and '$' (end) are anchor characters
                    writer4.write("rm -r restart\n");
                    writer4.write("rm -r vizout\n");
                    writer4.write("rm -rf $(ls * | grep -v '^Stdout.txt$' | grep -v '^Stderr.txt$' | grep -v '^perform' | grep -v '^ERROR$' | grep -v '.stdout$' | grep -v '.stderr$' )\n");
                } catch (IOException ex) {
                    ex.printStackTrace();
                }
                writer4.close();

            } catch (IOException ex) {
                ex.printStackTrace();
            }

        }
        // end of creating postprocess script

        //User specified commands script
        if (pre != null && !pre.isEmpty()) {
            try {
                BufferedWriter writer5 = new BufferedWriter(new FileWriter(basedir + File.separator + "userCommands.sh"));
                try {
                    writer5.write("#!/bin/bash\n");
                    writer5.write("\n");
                    writer5.write(pre);
                } catch (IOException ex) {
                    ex.printStackTrace();
                }
                writer5.close();
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
        // end of user specified commands
    }

    /**
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
            if (use_irods == 1) {
                File execFile3 = new File(basedir + File.separator + "preprocess_irods.sh");
                execFile3.setExecutable(true, false);
            }
            if (output_to_irods == 1) {
                File execFile4 = new File(basedir + File.separator + "postprocess_irods.sh");
                execFile4.setExecutable(true, false);
            }
            //change to current directory
            //Runtime.getRuntime().exec("cd " + basedir);
            //run submit
            Process jobProcess = Runtime.getRuntime().exec(basedir + File.separator + execName);
            jobDetails.setProcess(jobProcess);
            //Runtime.getRuntime().exec(basedir + File.separator + EXEC_NAME);
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
