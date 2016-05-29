package titan.gui;


import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.FileNotFoundException;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.FileInputStream;
import java.io.InputStream;

import java.net.InetAddress;
import java.net.UnknownHostException;
import java.lang.Runtime;
import java.lang.Exception;
import java.lang.Process;

import javax.swing.*;

import titan.graphics.LabelDisplay;
import titan.io.GenericReader;
import titan.io.TitanDBAccess.TitanJobMonitor;
import titan.jobs.JobStatus;
import titan.jobs.JobSubmissionContainer;
import titan.jobs.TitanRunInput;
import titan.jobs.condor.RunMethodCondor;

public class JobDetailsDialog extends JDialog {

    public static JLabel _label1; //used for the imagel label ~ contour image
    public static JLabel _label2; //used for the imagel label ~ 3D image
    public JButton kmlButton; //The button for making kml

    private String outputFolder;

    public JobDetailsDialog(Component owner, TitanJobMonitor data, JobSubmissionContainer sc) {

        getContentPane().setLayout(new BoxLayout(getContentPane(), BoxLayout.Y_AXIS));
        setTitle("Job Details");
        setResizable(true);
        this.setLocation(100, 100);
        setModal(false);
        this.setDefaultCloseOperation(JDialog.HIDE_ON_CLOSE);

        kmlButton = new JButton("Make KML");
        kmlButton.setEnabled(true);

        JPanel detailsPanel = new JPanel();
        detailsPanel.setLayout(new BoxLayout(detailsPanel, BoxLayout.Y_AXIS));

        LabelDisplay jobName = new LabelDisplay("Job Name : ", data.jobName);
        LabelDisplay jobID = new LabelDisplay("Job ID : ", data.jobID);
        LabelDisplay submitMethod = new LabelDisplay("Job Submit Method : ", data.submitMethod);
        LabelDisplay submitHost = new LabelDisplay("Host Job Submitted On : ", data.submitHost);
        LabelDisplay dateSubmitted = new LabelDisplay("Date/Time Job Submitted : ", data.dateSubmitted);
        LabelDisplay status;

        data.status = sc.getJobStatus(data.submitMethod, data.jobID, data.outputDirectory);

        status = new LabelDisplay("JobStatus : ", data.status);

        LabelDisplay outputDirectory = new LabelDisplay("Output Directory : ", data.outputDirectory);

        JTabbedPane detailsPane = new JTabbedPane();

        // display of output file

        JScrollPane stdoutScroll = null;
        JTextArea stdout = createBrowseText();

        try {

            StringBuilder stringBuilder = new StringBuilder();
            BufferedReader reader = new BufferedReader(new FileReader(data.outputDirectory + File.separator + sc.runMethod.STDOUT));
            String line;

            while ((line = reader.readLine()) != null) {
                stringBuilder.append(line);
                stringBuilder.append("\n");
            }

            reader.close();
            stringBuilder.deleteCharAt(stringBuilder.length() - 1);
            stdout.setText(stringBuilder.toString());

        } catch (FileNotFoundException e) {
            stdout.setText("THE FILE APPPEARS TO BE MISSING.");
        } catch (Exception e) {
            e.printStackTrace();
        }

        //stdout.append(addTextFileContent(
        //		data.outputDirectory + File.separator + RunMethod.STDOUT));
        stdout.setCaretPosition(0);
        stdoutScroll = new JScrollPane(stdout);
        stdoutScroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);

        // display of error  file
        JScrollPane stderrScroll = null;
        JTextArea stderr = createBrowseText();
        stderr.append(addTextFileContent(
                data.outputDirectory + File.separator + sc.runMethod.STDERR));
        stderr.setCaretPosition(0);
        stderrScroll = new JScrollPane(stderr);
        stderrScroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);

        // display of simulation.data  file
        JScrollPane simulationScroll = null;
        JTextArea simulation = createBrowseText();
        simulation.append(addTextFileContent(
                data.outputDirectory + File.separator + TitanRunInput.SIMULATION_DATA));
        simulation.setCaretPosition(0);
        simulationScroll = new JScrollPane(simulation);
        simulationScroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);

        // display of HUB submit info  file
        JScrollPane hubSubmitInfoScroll = null;
        JTextArea hubSubmitInfo = createBrowseText();
        hubSubmitInfo.append(addTextFileContent(
                data.outputDirectory + File.separator + JobSubmissionContainer.HUB_SUBMIT_INFO));
        hubSubmitInfo.setCaretPosition(0);
        hubSubmitInfoScroll = new JScrollPane(hubSubmitInfo);
        hubSubmitInfoScroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);

        // display of executable shell script file
        JScrollPane runtitanScroll = null;
        JTextArea runtitan = createBrowseText();
        runtitan.append(addTextFileContent(
                data.outputDirectory + File.separator + "runtitan.sh"));
        runtitan.setCaretPosition(0);
        runtitanScroll = new JScrollPane(runtitan);
        runtitanScroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);

        // display of executable shell script file
        JScrollPane submittitanScroll = null;
        JTextArea submittitan = createBrowseText();
        submittitan.append(addTextFileContent(
                data.outputDirectory + File.separator + "submittitan.sh"));
        submittitan.setCaretPosition(0);
        submittitanScroll = new JScrollPane(submittitan);
        submittitanScroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);

        // display perform files
        JScrollPane performScroll = null;
        JTextArea perform = createBrowseText();
        //retrieve only the files that start with string 'perform'
        File dir = new File(data.outputDirectory);
        FilenameFilter filter = new FilenameFilter() {
            public boolean accept(File dir, String name) {
                return name.startsWith("perform");
            }
        };
        String[] children = dir.list(filter);
        for (String fileName : children) {
            perform.append("===> " + fileName + " <===\n");
            perform.append(addTextFileContent(
                    data.outputDirectory + File.separator + fileName) + "\n\n");
        }
        perform.setCaretPosition(0);
        performScroll = new JScrollPane(perform);
        performScroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);

        //Display images
        String filePath = dir.getAbsolutePath();
        outputFolder = filePath;
        File animationContour = new File(filePath + "/vizout/animation.gif");
        File animation3D = new File(filePath + "/vizout/animation3D.gif");
        JScrollPane imageCountourScroll = null;
        JScrollPane image3DScroll = null;
        String outputdir_1 = "";  //used for iRODS
        String outputdir_2 = "";  //used for iRODS
        String pathfile1 = "";
        String pathfile2 = "";
        //Check if images have already been created
        if (animationContour.exists() && animation3D.exists()) {
            // Display the images
            ImageIcon i1 = new ImageIcon(filePath + "/vizout/animation.gif");
            ImageIcon i2 = new ImageIcon(filePath + "/vizout/animation3D.gif");
            imageCountourScroll = new JScrollPane(new JLabel(i1));
            imageCountourScroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
            image3DScroll = new JScrollPane(new JLabel(i2));
            image3DScroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
        } else {
            //call driver.sh to create images
            //get the directory from the environment variable

            File h5file = new File(filePath + "/vizout/xdmf_p0000_i00000000.h5");
            File tempDir = new File(filePath + "/vizout/script_output");

            if (h5file.exists() && (data.status == "Successful")) {

                _label1 = new JLabel("Creating images, please wait ...");
                _label2 = new JLabel("Creating images, please wait ...");

                imageCountourScroll = new JScrollPane(_label1);
                imageCountourScroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
                image3DScroll = new JScrollPane(_label2);
                image3DScroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);

                if (!(tempDir.exists())) {  // if thread to create images has not already been created, create it
                    Thread thread = this.new ImageOutput(filePath);
                    thread.start();
                }

            } else {

                //if the files are not in the output dir, then check if they are in iRODS
                outputdir_1 = filePath.substring(filePath.lastIndexOf(File.separator) + 1);  //get final output dir
                String tmp_dir = filePath.substring(0, filePath.lastIndexOf(File.separator));
                outputdir_2 = tmp_dir.substring(tmp_dir.lastIndexOf(File.separator) + 1);  //get possible run name

                //System.out.println("outputdir_1: " + outputdir_1);
                //System.out.println("tmp_dir: " + tmp_dir);
                //System.out.println("outputdir_2: " + outputdir_2);

                pathfile1 = java.lang.System.getProperty("user.home") + "/irods/" +
                        java.lang.System.getProperty("user.name") + "/TitanVhubResults/" +
                        outputdir_1;

                pathfile2 = java.lang.System.getProperty("user.home") + "/irods/" +
                        java.lang.System.getProperty("user.name") + "/TitanVhubResults/" +
                        outputdir_2 + File.separator + outputdir_1;

                //System.out.println("path1: " + pathfile1);
                //System.out.println("path2: " + pathfile2);

                File h5file_irods1 = new File(pathfile1 + "/vizout/xdmf_p0000_i00000000.h5");
                File h5file_irods2 = new File(pathfile2 + "/vizout/xdmf_p0000_i00000000.h5");

                if (h5file_irods1.exists() && (data.status == "Successful")) {

                    //file exists in iRODS
                    //System.out.println("files in iRODS:");
                    //System.out.println( pathfile1 + "/vizout/xdmf_p0000_i00000000.h5" );

                    _label1 = new JLabel("Creating images, please wait ...");
                    _label2 = new JLabel("Creating images, please wait ...");

                    imageCountourScroll = new JScrollPane(_label1);
                    imageCountourScroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
                    image3DScroll = new JScrollPane(_label2);
                    image3DScroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);

                    File tempDir1 = new File(pathfile1 + "/vizout/script_output");
                    if (!(tempDir1.exists())) {  // if thread to create images has not already been created, create it
                        Thread thread = this.new ImageOutput(filePath, pathfile1);
                        thread.start();
                    }

                } else if (h5file_irods2.exists() && (data.status == "Successful")) {

                    //file exists in iRODS
                    //System.out.println("files in iRODS:");
                    //System.out.println( pathfile2 + "/vizout/xdmf_p0000_i00000000.h5" );

                    _label1 = new JLabel("Creating images, please wait ...");
                    _label2 = new JLabel("Creating images, please wait ...");

                    imageCountourScroll = new JScrollPane(_label1);
                    imageCountourScroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
                    image3DScroll = new JScrollPane(_label2);
                    image3DScroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);

                    File tempDir2 = new File(pathfile2 + "/vizout/script_output");
                    if (!(tempDir2.exists())) {  // if thread to create images has not already been created, create it
                        Thread thread = this.new ImageOutput(filePath, pathfile2);
                        thread.start();
                    }
                }
            }
        }

        // display of condor description file if necessary
        //  hardcoded string : hate it
        JScrollPane condorScroll = null;
        JTextArea condor;
        if (data.submitMethod.compareTo("Condor") == 0) {
            condorScroll = null;
            condor = createBrowseText();
            condor.append(addTextFileContent(
                    data.outputDirectory + File.separator + RunMethodCondor.CONDOR_SCRIPT));
            condor.setCaretPosition(0);
            condorScroll = new JScrollPane(condor);
            condorScroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
        }

        // display additional error details if necessary
        JScrollPane errorScroll = null;
        JTextArea error;
        if (data.status.compareTo(JobStatus.STATUS_UNKNOWN) == 0) {
            error = createBrowseText();
            error.append(addJobStatusDetails(data));
            error.setCaretPosition(0);
            errorScroll = new JScrollPane(error);
            errorScroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
        }

        detailsPanel.add(jobName.getPanel());
        detailsPanel.add(jobID.getPanel());
        detailsPanel.add(submitMethod.getPanel());
        detailsPanel.add(submitHost.getPanel());
        detailsPanel.add(dateSubmitted.getPanel());
        detailsPanel.add(status.getPanel());
        detailsPanel.add(outputDirectory.getPanel());

        detailsPane.addTab(sc.runMethod.STDOUT, stdoutScroll);
        detailsPane.addTab(sc.runMethod.STDERR, stderrScroll);
        detailsPane.addTab(TitanRunInput.SIMULATION_DATA, simulationScroll);
        detailsPane.addTab("runtitan.sh", runtitanScroll);
        detailsPane.addTab("submittitan.sh", submittitanScroll);
        detailsPane.addTab("HUB-Submit Info", hubSubmitInfoScroll);
        // Performance parameters are now written to Stdout.txt
        //detailsPane.addTab("perform", performScroll);
        //only display images if h5 files are present

        if (animationContour.exists() && animation3D.exists()) {
            detailsPane.addTab("Contour GIF", imageCountourScroll);
            detailsPane.addTab("3D GIF", image3DScroll);
        } else {

            File h5file = new File(filePath + "/vizout/xdmf_p0000_i00000000.h5");
            if (h5file.exists() && (data.status == "Successful")) {
                detailsPane.addTab("Contour GIF", imageCountourScroll);
                detailsPane.addTab("3D GIF", image3DScroll);
            }

            //iRODS
            if (!pathfile1.isEmpty()) {
                //System.out.println("outputdir is not empty: "+pathfile1);
                File h5file_irods1 = new File(pathfile1 + "/vizout/xdmf_p0000_i00000000.h5");
                File h5file_irods2 = new File(pathfile2 + "/vizout/xdmf_p0000_i00000000.h5");
                if (h5file_irods1.exists() && (data.status == "Successful")) {
                    //System.out.println("creating tabs ...");
                    detailsPane.addTab("Contour GIF", imageCountourScroll);
                    detailsPane.addTab("3D GIF", image3DScroll);
                } else if (h5file_irods2.exists() && (data.status == "Successful")) {
                    //System.out.println("creating tabs ...");
                    detailsPane.addTab("Contour GIF", imageCountourScroll);
                    detailsPane.addTab("3D GIF", image3DScroll);
                }

            } else {
                //System.out.println("outputdir is empty: "+outputdir_1);
            }
        }

        if (data.submitMethod.compareTo("Condor") == 0) {
            detailsPane.addTab(RunMethodCondor.CONDOR_SCRIPT, condorScroll);
        }
        if (data.status.compareTo(JobStatus.STATUS_UNKNOWN) == 0) {
            detailsPane.addTab("ERROR", errorScroll);
        }

        detailsPanel.add(detailsPane);

        getContentPane().add(detailsPanel);

        JButton closeButton = new JButton("Close");
        closeButton.addActionListener(new CloseButtonListener());


        kmlButton.addActionListener(new MakeKmlListener());

        JPanel buttonPanel = new JPanel();
        buttonPanel.add(closeButton);
        buttonPanel.add(kmlButton);
        getContentPane().add(buttonPanel);

        pack();

        this.addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent e) {
                JobDetailsDialog.this.dispose();
            }
        });
    }

    private JTextArea createBrowseText() {

        JTextArea details = new JTextArea();
        details.setLineWrap(true);
        details.setEditable(false);
        details.setRows(20);
        return details;
    }

    private String addJobStatusDetails(TitanJobMonitor data) {
        String summary = new String("Additional Status Details :\n");
        String localHost;
        try {
            localHost = InetAddress.getLocalHost().getHostName();
        } catch (UnknownHostException ex) {
            localHost = "";
        }

        // YES, I hate the hardcoded string here
        if (data.submitMethod.compareTo("Local") == 0) {
            summary = summary.concat("  The job is not active on the local machine.\n");
            if (data.submitHost.compareTo(localHost) != 0) {
                summary = summary.concat("    The job was not submitted from this host and");
                summary = summary.concat("the job list for the machine " + data.submitHost);
                summary = summary.concat(" cannot be queried.  The job may be running on ");
                summary = summary.concat(data.submitHost + ".\n");
            }
	} else if (data.submitMethod.compareTo("PBS") == 0) {
            summary = summary.concat("  The PBS job manager does not appear to hold the job.\n");
            if (data.submitHost.compareTo(localHost) != 0) {
                summary = summary.concat("    The job was not submitted from this host and");
                summary = summary.concat("the PBS job manager on the local machine may not");
                summary = summary.concat("be the same as the PBS job manager on ");
                summary = summary.concat(data.submitHost + ".\n");
            }
        } else if (data.submitMethod.compareTo("Condor") == 0) {
            summary = summary.concat("  The Condor job manager does not appear to hold the job.\n");
            if (data.submitHost.compareTo(localHost) != 0) {
                summary = summary.concat("    The job was not submitted from this host and");
                summary = summary.concat("the Condor job manager on the local machine may not");
                summary = summary.concat("be the same as the Condor job manager on ");
                summary = summary.concat(data.submitHost + ".\n");
            }
       } else if (data.submitMethod.compareTo("Hub-Submit") == 0) {
            summary = summary.concat("  The Hub-Submit job manager does not appear to hold the job.\n");
            summary = summary.concat("  The Hub-Submit Run Style option is valid on VHub only.\n"); 
            summary = summary.concat("  Please verify the setting for the E_VHUB environment variable in <path to root titan2d>/bin/titan_gui.sh.\n"); 
            if (data.submitHost.compareTo(localHost) != 0) {
                summary = summary.concat("    The job was not submitted from this host and");
                summary = summary.concat("the job list for the machine " + data.submitHost);
                summary = summary.concat(" cannot be queried.  The job may be running on ");
                summary = summary.concat(data.submitHost + ".\n");
            }
        }
        File outputDir = new File(data.outputDirectory);
        if (!outputDir.exists()) {
            summary = summary.concat("  The output directory for the job does not exist.\n");
            if (data.submitHost.compareTo(localHost) != 0) {
                summary = summary.concat("    The job was not submitted from this host and");
                summary = summary.concat("the output directory may not be visible ");
                summary = summary.concat("from the local machine.\n");
            } else {
                summary = summary.concat("    It appears the directory has been removed.\n");
            }
        }
        return summary;
    }

    private String addTextFileContent(String textFile) {
        String content = new String("");
        File inFile = new File(textFile);

        try {
            GenericReader reader = new GenericReader(inFile);
            content = content.concat(reader.read());
        } catch (FileNotFoundException e) {
            content = new String("THE FILE APPPEARS TO BE MISSING.");
        } catch (Exception e) {
            e.printStackTrace();
        }
        return content;
    }

    private class CloseButtonListener implements ActionListener {
        public void actionPerformed(ActionEvent ev) {
            dispose();
        }
    }

    private class MakeKmlListener implements ActionListener {
        public void actionPerformed(ActionEvent ev) {
            String filePath = outputFolder;
            System.out.println("filePath: " + filePath);
            Thread thread = new KmlOutput(filePath);
            thread.start();
        }
    }

    //used to create images
    private class ImageOutput extends Thread {

        public String filePath;
        public String filePath_irods = "";

        public ImageOutput(String f) {
            filePath = f;
        }

        public ImageOutput(String f, String i) {
            filePath = f;
            filePath_irods = i;
        }

        public void run() {

            String tooldir = System.getenv("TITAN2D_HOME");
            kmlButton.setEnabled(false);

            try {
                Process p;
                if (filePath_irods.isEmpty()) {
                    p = Runtime.getRuntime().exec(tooldir + "/bin/driverMultiProcessor.sh " + filePath + " " + tooldir);
                } else {
                    p = Runtime.getRuntime().exec(tooldir + "/bin/driverMultiProcessor.sh " + filePath_irods + " " + tooldir + " " + filePath);
                }
                p.waitFor();
                //Debug code (comment out later)
                BufferedReader input = new BufferedReader(new InputStreamReader(p.getInputStream()));
                String line;
                while ((line = input.readLine()) != null) {
                    System.err.println(line);
                }
                input.close();
                BufferedReader error = new BufferedReader(new InputStreamReader(p.getErrorStream()));
                while ((line = error.readLine()) != null) {
                    System.err.println(line);
                }
                error.close();
            } catch (Exception e) {
                System.out.println(e.toString() + "\n");
                e.printStackTrace();
            }

            //display images
            ImageIcon i1 = new ImageIcon(filePath + "/vizout/animation.gif");
            ImageIcon i2 = new ImageIcon(filePath + "/vizout/animation3D.gif");

            //update images
            _label1.setText(null);
            _label2.setText(null);

            _label1.removeAll();
            _label2.removeAll();

            _label1.setIcon(i1);
            _label2.setIcon(i2);

            _label1.setHorizontalAlignment(SwingConstants.CENTER);
            _label2.setHorizontalAlignment(SwingConstants.CENTER);

            _label1.updateUI();
            _label2.updateUI();

            kmlButton.setEnabled(true);
       }
    }

    //used to create kml files
    private class KmlOutput extends Thread {

        public String filePath;
        public String filePath_irods = "";

        public KmlOutput(String f) {
            filePath = f;
        }

        public KmlOutput(String f, String i) {
            filePath = f;
            filePath_irods = i;
        }

        public void run() {
            kmlButton.setEnabled(false);
            String tooldir = System.getenv("TITAN2D_HOME");
            System.err.println(tooldir);

            try {
                Process p;
                if (filePath_irods.isEmpty()) {
                    System.err.println(tooldir + "/bin/driverKML.sh " + filePath + " " + tooldir);
                    p = Runtime.getRuntime().exec(tooldir + "/bin/driverKML.sh " + filePath + " " + tooldir);
                } else {
                    System.err.println(tooldir + "/bin/driverKML.sh " + filePath_irods + " " + tooldir);
                    p = Runtime.getRuntime().exec(tooldir + "/bin/driverKML.sh " + filePath_irods + " " + tooldir);
                }
                //p.waitFor();
                BufferedReader input = new BufferedReader(new InputStreamReader(p.getInputStream()));
                String line;
                while ((line = input.readLine()) != null) {
                    System.err.println(line);
                }
                input.close();
                BufferedReader error = new BufferedReader(new InputStreamReader(p.getErrorStream()));
                while ((line = error.readLine()) != null) {
                    System.err.println(line);
                }
                error.close();
            } catch (Exception e) {
                System.out.println(e.toString() + "\n");
                e.printStackTrace();
            }

            //Post-run update
            kmlButton.setEnabled(true);
        }
    }
}
