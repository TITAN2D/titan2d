package titan.gui;

import java.io.*;

import java.lang.Thread;

import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.net.InetAddress;
import java.net.UnknownHostException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.sql.SQLException;
import java.util.Calendar;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextArea;
import javax.swing.border.TitledBorder;
import javax.swing.JLabel;

import titan.graphics.RadioButtonGroup;
import titan.graphics.TextArea;
import titan.graphics.TextInput;
import titan.gui.Titan.FetchData;
import titan.gui.Titan.SaveTitanRun;
import titan.gui.Titan.LoadTitanRun;
import titan.gui.TabMaterialModelMap;
import titan.io.NameValuePairGroup;
import titan.io.TitanDBAccess;
import titan.io.TitanDBAccess.TitanJobMonitor;
import titan.jobs.JobDetails;
import titan.jobs.JobSubmissionContainer;
import titan.jobs.JobSubmissionParameters;
import titan.jobs.RunMethod;
import titan.jobs.TitanRunInput;
import titan.jobs.TitanRunInput.TitanSimulationData;
import titan.jobs.TitanRunInput.TitanScaleData;
import titan.jobs.TitanRunInput.TitanSimulationFluxSource;
import titan.jobs.TitanRunInput.TitanSimulationPile;
import titan.jobs.TitanRunInput.TitanSimulationPlane;
import titan.options.OptionsManager;

public class TabJobSubmission extends JPanel {


    private JButton runButton;
    private FetchData dataFetch;
    private TitanRunInput simWriter;
    private JTextArea submitSummary;
    private JTabbedPane tp;
    private TextInput nodesLocal;         // Local
    private TextInput cpusLocal;          // Local
    private TextInput maxRuntimePBS;      // PBS
    private TextInput queueNamePBS;       // PBS
    private TextInput nodesPBS;           // PBS
    private TextInput cpusPBS;            // PBS
    private TextInput requirementsCondor; // Condor
    private TextInput maxRuntimeSubmit;   // Hub-Submit
    private TextInput nodesSubmit;        // Hub-Submit
    private TextInput cpusSubmit;         // Hub-Submit
    private JLabel suggestedCPU;          // Hub-Submit
    private JLabel suggestedCPU2;         // Hub-Submit
    private JLabel suggestedCPU3;         // Hub-Submit
    private TextArea additional;
    private TitanDBAccess dba;
    private RadioButtonGroup runStyle;
    private RadioButtonGroup runMode;
    private TextInput runArgs;
    private OptionsManager optionsManager;
    private SaveTitanRun save;
    private LoadTitanRun load;
    private TabMaterialModelMap matModelMap;
    private JobSubmissionContainer submitContainer;

    // Note: on vhub only 1 node and 1 cpu are allowed per job.
    // E_VHUB environment variable read in Titan.java
    public boolean eVHUB;

    public TabJobSubmission(FetchData create, TitanDBAccess db,
                            OptionsManager om, Titan.SaveTitanRun sv,
                            TabMaterialModelMap mmm,
                            JobSubmissionContainer sc, boolean eVHUBFlag) {

        dba = db;
        optionsManager = om;
        save = sv;
        matModelMap = mmm;
        submitContainer = sc;

        eVHUB = eVHUBFlag;

        setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));

        String methods[] = submitContainer.getMethodNames();

        runStyle = new RadioButtonGroup("Run Style", methods, RadioButtonGroup.SINGLE_SELECTION);
        runStyle.setEditable(submitContainer.RUN_STYLE_PBS, false);
        runStyle.setEditable(submitContainer.RUN_STYLE_CONDOR, false);
        if (eVHUB == false) {
            runStyle.setEditable(submitContainer.RUN_STYLE_SUBMIT, false);
        }


        for (int i = 0; i < methods.length; i++) {

            final String method = methods[i];

            runStyle.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    initializeButtons(method);
                }
            }, method);
        }

        add(runStyle.getPanel());

        // Local job submission options
        nodesLocal = new TextInput("Nodes");
        nodesLocal.setValue("1");
        cpusLocal = new TextInput("CPUs per node");
        cpusLocal.setValue("1");

        JPanel localPanel = new JPanel();
        localPanel.setLayout(new BoxLayout(localPanel, BoxLayout.Y_AXIS));
        localPanel.add(nodesLocal.getPanel());
        localPanel.add(cpusLocal.getPanel());

        // PBS job submission options (PBS job submissions currently disabled)
        maxRuntimePBS = new TextInput("Maximum Run Time");
        queueNamePBS = new TextInput("Queue");
        nodesPBS = new TextInput("Nodes");
        nodesPBS.setValue("1");
        cpusPBS = new TextInput("CPUs per node");
        cpusPBS.setValue("1");

        JPanel pbsPanel = new JPanel();
        pbsPanel.setLayout(new BoxLayout(pbsPanel, BoxLayout.Y_AXIS));
        pbsPanel.add(maxRuntimePBS.getPanel());
        pbsPanel.add(queueNamePBS.getPanel());
        pbsPanel.add(nodesPBS.getPanel());
        pbsPanel.add(cpusPBS.getPanel());

        // Condor job submission options (Condor job submissions currently disabled)
        requirementsCondor = new TextInput("Requirements");
        requirementsCondor.setValue("");
        JPanel condorPanel = new JPanel();
        condorPanel.setLayout(new BoxLayout(condorPanel, BoxLayout.Y_AXIS));
        condorPanel.add(requirementsCondor.getPanel());
        requirementsCondor.setEditable(false);

        // Hub-submit job submission options
        maxRuntimeSubmit = new TextInput("Max Run Time ( 72hrs Max )");
        maxRuntimeSubmit.setToolTipText("hh:mm:ss");
        maxRuntimeSubmit.setValue("01:00:00");
        nodesSubmit = new TextInput("Nodes");
        nodesSubmit.setValue("1");
        cpusSubmit = new TextInput("CPUs per node");
        cpusSubmit.setValue("1");

        suggestedCPU = new JLabel("NOTE: 1 or 2 CPUs are suggested for small DEM's (32MB/process), although there is a max of 8 CPUs allowed per job.", javax.swing.SwingConstants.LEFT);
        suggestedCPU2 = new JLabel("Also note that, larger processor counts can experience longer queue waiting times on the cluster.", javax.swing.SwingConstants.LEFT);
        suggestedCPU3 = new JLabel("Currently there is only 1 node allowed per job, however this restraint will be changed in the future.", javax.swing.SwingConstants.LEFT);
        //suggestedCPU.setVerticalAlignment(javax.swing.SwingConstants.CENTER);
        suggestedCPU.setBorder(javax.swing.BorderFactory.createLineBorder(java.awt.Color.black));
        suggestedCPU2.setBorder(javax.swing.BorderFactory.createLineBorder(java.awt.Color.black));
        suggestedCPU3.setBorder(javax.swing.BorderFactory.createLineBorder(java.awt.Color.black));

        JPanel submitJobPanel = new JPanel();
        submitJobPanel.setLayout(new BoxLayout(submitJobPanel, BoxLayout.Y_AXIS));
        submitJobPanel.add(maxRuntimeSubmit.getPanel());
        submitJobPanel.add(nodesSubmit.getPanel());
        submitJobPanel.add(cpusSubmit.getPanel());
        submitJobPanel.add(suggestedCPU);
        submitJobPanel.add(suggestedCPU2);
        submitJobPanel.add(suggestedCPU3);

        JPanel runStylePanel = new JPanel();
        runStylePanel.setLayout(new BoxLayout(runStylePanel, BoxLayout.Y_AXIS));
        runStylePanel.setBorder(new TitledBorder("Run Style Specific Options"));
        runStylePanel.setMaximumSize(new Dimension(2000, 200)); // (width, height) 2000 -> a large number

        tp = new JTabbedPane();
        // index 0
        tp.add("Local Options", localPanel);
        // index 1
        tp.add("PBS Options", pbsPanel);
        // index 2
        tp.add("Condor Options", condorPanel);
        // index 3
        tp.add("Submit Options", submitJobPanel);

        runStylePanel.add(tp);

        add(runStylePanel);

        additional = new TextArea("Pre-execution Commands");
        add(additional.getPanel());

        // Run mode panel
        JPanel runModePanel = new JPanel();
        runModePanel.setLayout(new BoxLayout(runModePanel, BoxLayout.Y_AXIS));
        runModePanel.setBorder(new TitledBorder("Run Mode"));

        // Run mode
        String runModes[] = submitContainer.getRunModes();

        runMode = new RadioButtonGroup("Run Command", runModes, RadioButtonGroup.SINGLE_SELECTION);
        runMode.setEditable(submitContainer.RUN_MODE_HYBRID, false);
        runMode.setValue(runModes[0]);

        runArgs = new TextInput("Run Submission Arguments");
        runMode.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ev) {
                runArgs.setEditable(false);
            }
        }, runModes[0]);
        runMode.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ev) {
                runArgs.setEditable(false);
            }
        }, runModes[1]);

        // rlj 1/12/2016
        //
        // Only allow openmp mode for now
        //runModePanel.add(runMode.getPanel());
        //runModePanel.add(runArgs.getPanel());
        //add(runModePanel);

        runButton = new JButton("Run Job");
        runButton.addActionListener(TabJobSubmission.this.new SaveAction());
        add(runButton);

        JPanel submitPanel = new JPanel();
        submitPanel.setMaximumSize(new Dimension(2000, 185));
        submitPanel.setLayout(new BoxLayout(submitPanel, BoxLayout.Y_AXIS));
        submitPanel.setBorder(new TitledBorder("Job Execution Summary"));

        submitSummary = new JTextArea();
        submitSummary.setLineWrap(true);
        submitSummary.setEditable(false);

        JScrollPane textPain = new JScrollPane(submitSummary);
        textPain.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
        submitPanel.add(textPain);

        add(submitPanel);

        dataFetch = create;

        //initialize to local
        runStyle.setValue(methods[0]);
        initializeButtons(methods[0]);
    }

    private void initializeButtons(String selected) {

        nodesLocal.setEditable(false);         // Local
        cpusLocal.setEditable(false);          // Local
        maxRuntimePBS.setEditable(false);      // PBS
        queueNamePBS.setEditable(false);       // PBS
        nodesPBS.setEditable(false);           // PBS
        cpusPBS.setEditable(false);            // PBS
        requirementsCondor.setEditable(false); // Condor
        maxRuntimeSubmit.setEditable(false);   // Hub-Submit
        nodesSubmit.setEditable(false);        // Hub-Submit
        cpusSubmit.setEditable(false);         // Hub-Submit

        // For first update of the tool only enabling the openmp run mode
        if (selected.compareTo(JobSubmissionContainer.RUN_STYLE_LOCAL) == 0) {
            // Note: on vhub only 1 node and 1 cpu are allowed per job.
            // Otherwise, let the user pick the number of cpus
            if (eVHUB == false) {
                cpusLocal.setEditable(true);
            }
            //runMode.setEditable(true);
            String runmodes[] = submitContainer.getRunModes();
            runMode.setValue(runmodes[0]);
            runMode.setEditable(false);
            runArgs.setEditable(false);
            tp.setSelectedIndex(0);
        } else if (selected.compareTo(JobSubmissionContainer.RUN_STYLE_PBS) == 0) {
            maxRuntimePBS.setEditable(true);
            queueNamePBS.setEditable(true);
            nodesPBS.setEditable(true);
            cpusPBS.setEditable(true);
            runMode.setEditable(true);
            if (runMode.getValue().compareTo(submitContainer.getRunModes()[0]) == 0 ||
                    runMode.getValue().compareTo(submitContainer.getRunModes()[1]) == 0)
                runArgs.setEditable(true);
            else
                runArgs.setEditable(false);
            tp.setSelectedIndex(1);
        } else if (selected.compareTo(JobSubmissionContainer.RUN_STYLE_CONDOR) == 0) {
            requirementsCondor.setEditable(true);
            runMode.setEditable(false);
            runArgs.setEditable(false);
            tp.setSelectedIndex(2);
        } else { // Submit

            if (eVHUB == true) {
                maxRuntimeSubmit.setEditable(true);
                cpusSubmit.setEditable(true);
                String runmodes[] = submitContainer.getRunModes();
                runMode.setValue(runmodes[0]);
                runMode.setEditable(false);
                runArgs.setEditable(false);
                tp.setSelectedIndex(3);
            }
        }
    }

    private class SaveAction implements ActionListener {
        public void actionPerformed(ActionEvent e) {

            // Send pile volume_fraction field for the TwoPhases-Pitman-Le model only,
            boolean TWOPHASES_PITMAN_LE_MODEL = false;

            Calendar currentDate = Calendar.getInstance();

            String dateString = new String(currentDate.get(Calendar.YEAR) + "_");
            if ((currentDate.get(Calendar.MONTH) + 1) > 9)
                dateString = dateString.concat((currentDate.get(Calendar.MONTH) + 1) + "_");
            else
                dateString = dateString.concat("0" + (currentDate.get(Calendar.MONTH) + 1) + "_");
            if (currentDate.get(Calendar.DAY_OF_MONTH) > 9)
                dateString = dateString.concat(currentDate.get(Calendar.DAY_OF_MONTH) + "_");
            else
                dateString = dateString.concat("0" + currentDate.get(Calendar.DAY_OF_MONTH) + "_");
            if (currentDate.get(Calendar.HOUR_OF_DAY) > 9)
                dateString = dateString.concat(currentDate.get(Calendar.HOUR_OF_DAY) + "_");
            else
                dateString = dateString.concat("0" + currentDate.get(Calendar.HOUR_OF_DAY) + "_");
            if (currentDate.get(Calendar.MINUTE) > 9)
                dateString = dateString.concat(currentDate.get(Calendar.MINUTE) + "_");
            else
                dateString = dateString.concat("0" + currentDate.get(Calendar.MINUTE) + "_");
            if (currentDate.get(Calendar.SECOND) > 9)
                dateString = dateString.concat(currentDate.get(Calendar.SECOND) + ".");
            else
                dateString = dateString.concat("0" + currentDate.get(Calendar.SECOND) + ".");
            if (currentDate.get(Calendar.MILLISECOND) > 99)
                dateString = dateString.concat(Integer.toString(currentDate.get(Calendar.MILLISECOND)));
            else if (currentDate.get(Calendar.MILLISECOND) > 9)
                dateString = dateString.concat("0" + currentDate.get(Calendar.MILLISECOND));
            else
                dateString = dateString.concat("00" + currentDate.get(Calendar.MILLISECOND));

            // Read data from the tabs
            NameValuePairGroup mainData;
            NameValuePairGroup[] pileData;
            NameValuePairGroup[] fluxData;
            NameValuePairGroup[] planeData;

            try {
                mainData = dataFetch.getData();
            } catch (NullPointerException ex) {
                JOptionPane.showMessageDialog(TabJobSubmission.this,
                        "Problem reading data from the main tabs.\n" +
                                "Unable to submit the job.",
                        "Job Submission Error",
                        JOptionPane.ERROR_MESSAGE);
                return;
            }

            try {
                pileData = dataFetch.getPileData();
            } catch (NullPointerException ex) {
                JOptionPane.showMessageDialog(TabJobSubmission.this,
                        "Problem reading data from the Piles tab information stored in the tool's database.\n" +
                                "Please verify no blank fields on the Piles tab;\n" +
                                "when modifying a field, click enter <Enter> or another field.\n\n" +
                                "Note: this will occur when loading an Input Directory from a previous version of the Titan2D tool;\n" +
                                "to fix, set the Pile(s) Type(s) field(s) on the Piles Tab.  Please see Help for more information.\n\n" +
                                "Note: this will occur after changing the Material Model of a loaded test to TwoPhases-Pitman-Le;\n" +
                                "to fix, set the Volume Fraction field(s) on the Piles tab.\n\n" +
                                "Unable to submit the job.",
                        "Job Submission Error",
                        JOptionPane.ERROR_MESSAGE);
                return;
            }
            try {
                fluxData = dataFetch.getFluxData();
            } catch (NullPointerException ex) {
                JOptionPane.showMessageDialog(TabJobSubmission.this,
                        "Problem reading data from the Flux Sources tab information stored in the tool's database.\n" +
                                "Please verify no blank fields on the Flux Sources tab;\n" +
                                "when modifying a field, click enter <Enter> or another field.\n" +
                                "Unable to submit the job.",
                        "Job Submission Error",
                        JOptionPane.ERROR_MESSAGE);
                return;
            }
            try {
                planeData = dataFetch.getPlaneData();
            } catch (NullPointerException ex) {
                JOptionPane.showMessageDialog(TabJobSubmission.this,
                        "Problem reading data from the Discharge Planes tab information stored in the tool's database.\n" +
                                "Please verify no blank fields on the Discharge Planes tab;\n" +
                                "when modifying a field, click enter <Enter> or another field.\n" +
                                "Unable to submit the job.",
                        "Job Submission Error",
                        JOptionPane.ERROR_MESSAGE);
                return;
            }

            // Restart enabled?

            String restartEnabled = mainData.getValue(TitanConstants.RESTART_ENABLED);

            // Run Type Restart parameters
            String restartFile = mainData.getValue(TitanConstants.RESTART_FILE);
            String restartRunDir = "";
            int restartMaxNumberTimeSteps = 0;
            float restartMaxTime = 0.0f;

            String parmDir = "";

            if (restartEnabled.compareTo(TitanConstants.TRUE) == 0) {

                // Restart is enabled

                if (restartFile.compareTo("") == 0) {
                    JOptionPane.showMessageDialog(TabJobSubmission.this,
                            "The Restart File field on the Run Parameters tab is null.\n" +
                                    "Unable to submit the job.",
                            "Job Submission Error",
                            JOptionPane.ERROR_MESSAGE);
                    return;
                }

                String restartDir = "";
                String restartInputDir = "";
                String restartSaveDir = "";
                String restartRunName = "";

                System.out.println("Restart Directory: " + restartFile);
                restartDir = restartFile.substring(0, restartFile.lastIndexOf(File.separator));
                System.out.println("restartDir: " + restartDir);
                restartRunDir = restartDir.substring(0, restartDir.lastIndexOf(File.separator));
                System.out.println("restartRunDir: " + restartRunDir);
                restartInputDir = restartRunDir.substring(0, restartRunDir.lastIndexOf(File.separator));
                System.out.println("restartInputDir: " + restartInputDir);
                restartRunName = restartInputDir.substring(restartInputDir.lastIndexOf(File.separator) + 1);
                System.out.println("restartRunName: " + restartRunName);
                restartSaveDir = restartInputDir.substring(0, restartInputDir.lastIndexOf(File.separator));
                System.out.println("restartSaveDir: " + restartSaveDir);

                if (restartSaveDir.compareTo(mainData.getValue(TitanConstants.SAVE_DIRECTORY)) == 0 &&
                        restartRunName.compareTo(mainData.getValue(TitanConstants.RUN_NAME)) == 0) {

                    // Restart file is in currently loaded input directory.
                    // Allow the restart job to run, do not save

                    parmDir = new String(mainData.getValue(TitanConstants.SAVE_DIRECTORY) + File.separator +
                            mainData.getValue(TitanConstants.RUN_NAME) + File.separator + dateString);
                    File file = new File(parmDir);
                    file.mkdir();

                } else {

                    // No or different input directory loaded
                    int selectionOption = JOptionPane.showConfirmDialog(TabJobSubmission.this,
                            "The Base Save Directory and/or Run Name field on the Load/Save tab are null or different from \n" +
                                    "the Base Save Directory and/or Run Name derived from the Restart File field on the Run Parameters tab.\n\n" +
                                    "Select the OK Option to save results to the derived Base Save Directory: " + restartSaveDir + " and Run Name: " + restartRunName + ".\n\n" +
                                    "Note: The OK option is only valid for Restart Files created locally and for jobs submitted locally.\n\n" +
                                    "Recommended practice: Select the CANCEL option and load the Input Directory: " + restartInputDir + " on the Load/Save tab.\n" +
                                    "After loading: " + restartInputDir + ", verify the Restart settings on the Run Parameters tab.\n",
                            "Job Submission Option",
                            JOptionPane.OK_CANCEL_OPTION);

                    if (selectionOption == JOptionPane.OK_OPTION) {

                        // Allow the restart job to run, do not save
                        parmDir = new String(restartSaveDir + File.separator + restartRunName + File.separator + dateString);
                        File file = new File(parmDir);
                        file.mkdir();

                    } else {
                        return;
                    }
                }
            } else {

                // Restart is not enabled

                if (mainData.getValue(TitanConstants.SAVE_DIRECTORY).compareTo("") == 0 ||
                        mainData.getValue(TitanConstants.RUN_NAME).compareTo("") == 0) {
                    JOptionPane.showMessageDialog(TabJobSubmission.this,
                            "The Base Save Directory and/or Run Name fields on the Load/Save tab are null.\n" +
                                    "Unable to submit the job.",
                            "Job Submission Error",
                            JOptionPane.ERROR_MESSAGE);
                    return;
                }

                try {
                    if (save.Save(mainData.getValue(TitanConstants.SAVE_DIRECTORY),
                            mainData.getValue(TitanConstants.RUN_NAME)) == false) {
                        JOptionPane.showMessageDialog(TabJobSubmission.this,
                                "Titan.Save error encountered.\n" +
                                        "Unable to submit the job.",
                                "Job Submission Error",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                } catch (NullPointerException ex) {
                    JOptionPane.showMessageDialog(TabJobSubmission.this,
                            "The Base Save Directory and/or Run Name field on the Load/Save tab are invalid.\n" +
                                    "Unable to submit the job.",
                            "Job Submission Error",
                            JOptionPane.ERROR_MESSAGE);
                    return;
                }

                parmDir = new String(mainData.getValue(TitanConstants.SAVE_DIRECTORY) + File.separator +
                        mainData.getValue(TitanConstants.RUN_NAME) + File.separator + dateString);
                File file = new File(parmDir);
                file.mkdir();
            }

            // parmDir should not be null
            if (parmDir.compareTo("") == 0) {
                JOptionPane.showMessageDialog(TabJobSubmission.this,
                        "Job Submission Check parmDir is null",
                        "Job Submission Error",
                        JOptionPane.ERROR_MESSAGE);
                return;
            }

            // Make directories for web viz output
            //if(mainData.getValue(TitanConstants.VIS_OUTPUT).contains(TitanConstants.VizTypes[TitanConstants.VIS_OUTPUT_WEBVIZ])) {
            //	String webvizDir = new String(parmDir + File.separator + "webviz");
            //    File webvizFile = new File(webvizDir);
            //    webvizFile.mkdir();
            //    webvizFile = new File(webvizDir + File.separator + "topo");
            //    webvizFile.mkdir();
            //   webvizFile = new File(webvizDir + File.separator + "data");
            //    webvizFile.mkdir();
            //   webvizFile = new File(webvizDir + File.separator + "textures");
            //    webvizFile.mkdir();
            //}

            if (restartEnabled.compareTo(TitanConstants.TRUE) == 0) {

                // Restart is enabled

                if (mainData.getValue(TitanConstants.MAX_NUM_TIME_STEPS).compareTo("") != 0) {
                    try {
                        restartMaxNumberTimeSteps = Integer.parseInt(mainData.getValue(TitanConstants.MAX_NUM_TIME_STEPS));
                    } catch (NumberFormatException ex) {
                        JOptionPane.showMessageDialog(TabJobSubmission.this,
                                "The Max Number of Time Steps field on the Run Parameters tab is invalid.\n" +
                                        "Unable to submit the restart job.",
                                "Job Submission Error",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                }

                if (mainData.getValue(TitanConstants.MAX_TIME).compareTo("") != 0) {
                    try {
                        restartMaxTime = Float.parseFloat(mainData.getValue(TitanConstants.MAX_TIME));
                    } catch (NumberFormatException ex) {
                        JOptionPane.showMessageDialog(TabJobSubmission.this,
                                "The Max Time field on the Run Parameters tab is invalid.\n" +
                                        "Unable to submit the restart job.",
                                "Job Submission Error",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                }

                // Copy files required for creating KML files

                File fn1;
                File fn2;
                fn1 = new File(restartRunDir + File.separator + "zone.txt");
                fn2 = new File(parmDir + File.separator + "zone.txt");
                try {
                    Files.copy(fn1.toPath(), fn2.toPath(), StandardCopyOption.REPLACE_EXISTING);
                } catch (IOException ex) {
                    int selectionOption = JOptionPane.showConfirmDialog(TabJobSubmission.this,
                            "Unable to copy the zone.txt file, required for creating KML files, from\n" +
                                    restartRunDir + " to\n" +
                                    parmDir + "\n\n" +
                                    "Select the OK Option to continue job submission.",
                            "Job Submission Option",
                            JOptionPane.OK_CANCEL_OPTION);

                    if (selectionOption != JOptionPane.OK_OPTION) {

                        return;

                    }
                }

                fn1 = new File(restartRunDir + File.separator + "height_scale_for_KML.data");
                fn2 = new File(parmDir + File.separator + "height_scale_for_KML.data");
                try {
                    Files.copy(fn1.toPath(), fn2.toPath(), StandardCopyOption.REPLACE_EXISTING);
                } catch (IOException ex) {
                    int selectionOption = JOptionPane.showConfirmDialog(TabJobSubmission.this,
                            "Unable to copy the height_scale_for_KML.data file, required for creating KML files, from\n" +
                                    restartRunDir + " to\n" +
                                    parmDir + "\n\n" +
                                    "Select the OK Option to continue job submission.",
                            "Job Submission Option",
                            JOptionPane.OK_CANCEL_OPTION);

                    if (selectionOption != JOptionPane.OK_OPTION) {

                        return;

                    }
                }
            } else {

                // Restart is not enabled

                // Create the titan input file simulation.py

                simWriter = new TitanRunInput(parmDir);

                TitanSimulationData data = simWriter.new TitanSimulationData();

                // Get and verify input

                // GIS tab

                data.format = mainData.getValue(TitanConstants.GIS_FORMAT);

                if (data.format.compareTo(TitanConstants.GisFormats[TitanConstants.GIS_FORMAT_GIS_GRASS]) == 0) {

                    // GIS_GRASS

                    //For HUB-Submits, the information main directory is zipped into a tar.gz file
                    if (runStyle.getValue().compareTo(JobSubmissionContainer.RUN_STYLE_SUBMIT) == 0) {
                        String fullPathTemp = mainData.getValue(TitanConstants.GIS_INFO_DIRECTORY);
                        int index = fullPathTemp.lastIndexOf(File.separator);
                        data.mainDirectory = fullPathTemp.substring(index + 1);
                    } else {
                        data.mainDirectory = mainData.getValue(TitanConstants.GIS_INFO_DIRECTORY);
                    }

                    data.subDirectory = mainData.getValue(TitanConstants.GIS_SUBDIR);
                    data.mapset = mainData.getValue(TitanConstants.GIS_MAPSET);
                    data.map = mainData.getValue(TitanConstants.GIS_MAP);
                    data.vector = mainData.getValue(TitanConstants.GIS_VECTOR);

                    if (data.mainDirectory.compareTo("") == 0 ||
                            data.subDirectory.compareTo("") == 0 ||
                            data.mapset.compareTo("") == 0 ||
                            data.map.compareTo("") == 0) {
                        JOptionPane.showMessageDialog(TabJobSubmission.this,
                                "The Information Main Directory, Sub-Directory, Map Set and/or Map field(s) on the GIS tab are blank.\n" +
                                        "Unable to submit the job.",
                                "Job Submission Error",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                }else {

                    // GDAL

                    //For HUB-Submits, the mapset directory is zipped into a tar.gz file
                    if (runStyle.getValue().compareTo(JobSubmissionContainer.RUN_STYLE_SUBMIT) == 0) {

                        String fullPathTemp;
                        String map;
                        String cellhddir;
                        String cellhd;
                        String mapsetdir;
                        String mapset;

                        fullPathTemp = mainData.getValue(TitanConstants.GIS_MAP);
                        map = fullPathTemp.substring(fullPathTemp.lastIndexOf(File.separator) + 1);
                        cellhddir = fullPathTemp.substring(0, fullPathTemp.lastIndexOf(File.separator));
                        cellhd = cellhddir.substring(cellhddir.lastIndexOf(File.separator) + 1);
                        mapsetdir = cellhddir.substring(0, cellhddir.lastIndexOf(File.separator));
                        mapset = mapsetdir.substring(mapsetdir.lastIndexOf(File.separator) + 1);
                        data.map = mapset + File.separator + cellhd + File.separator + map;
                    } else {
                        data.map = mainData.getValue(TitanConstants.GIS_MAP);
                    }

                    if (data.map.compareTo("") == 0) {
                       JOptionPane.showMessageDialog(TabJobSubmission.this,
                                "The Map field on the GIS tab is blank.\n" +
                                        "Unable to submit the job.",
                                "Job Submission Error",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                }

                if (mainData.getValue(TitanConstants.ZONEOVERRIDE).compareTo("") != 0) {
                    try {
                        data.zoneOverride = Integer.parseInt(mainData.getValue(TitanConstants.ZONEOVERRIDE));
                    } catch (NumberFormatException ex) {
                        JOptionPane.showMessageDialog(TabJobSubmission.this,
                                "The Zone Override field on the Gis tab is invalid.\n" +
                                        "Unable to submit the job.",
                                "Job Submission Error",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                }

                data.hemisphere = mainData.getValue(TitanConstants.HEMISPHERE);

                if (mainData.getValue(TitanConstants.MIN_X_LOC).compareTo("") != 0) {
                    try {
                        data.min_x_loc = Float.parseFloat(mainData.getValue(TitanConstants.MIN_X_LOC));
                    } catch (NumberFormatException ex) {
                        JOptionPane.showMessageDialog(TabJobSubmission.this,
                                "The Minimum X Location field on the Gis tab is invalid.\n" +
                                        "Unable to submit the job.",
                                "Job Submission Error",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                }

                if (mainData.getValue(TitanConstants.MIN_Y_LOC).compareTo("") != 0) {
                    try {
                        data.min_y_loc = Float.parseFloat(mainData.getValue(TitanConstants.MIN_Y_LOC));
                    } catch (NumberFormatException ex) {
                        JOptionPane.showMessageDialog(TabJobSubmission.this,
                                "The Minimum Y Location field on the Gis tab is invalid.\n" +
                                        "Unable to submit the job.",
                                "Job Submission Error",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                }

                if (mainData.getValue(TitanConstants.MAX_X_LOC).compareTo("") != 0) {
                    try {
                        data.max_x_loc = Float.parseFloat(mainData.getValue(TitanConstants.MAX_X_LOC));
                    } catch (NumberFormatException ex) {
                        JOptionPane.showMessageDialog(TabJobSubmission.this,
                                "The Maximum X Location field on the Gis tab is invalid.\n" +
                                        "Unable to submit the job.",
                                "Job Submission Error",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                }

                if (mainData.getValue(TitanConstants.MAX_Y_LOC).compareTo("") != 0) {
                    try {
                        data.max_y_loc = Float.parseFloat(mainData.getValue(TitanConstants.MAX_Y_LOC));
                    } catch (NumberFormatException ex) {
                        JOptionPane.showMessageDialog(TabJobSubmission.this,
                                "The Maximum Y Location field on the Gis tab is invalid.\n" +
                                        "Unable to submit the job.",
                                "Job Submission Error",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                }

			/*if(data.mainDirectory.startsWith(File.separator) == false) {
                JOptionPane.showMessageDialog(TabJobSubmission.this,
						"The GIS Main Directory on the GIS tab must be an\n" +
                        "absolute path.  Unable to submit the job.",
                        "Job Submission Error",
                        JOptionPane.ERROR_MESSAGE);
				return;
			}*/

                // Run Parameters Tab

                try {
                    data.cellsAcross = Integer.parseInt(mainData.getValue(TitanConstants.NUM_CELLS_ACROSS));
                } catch (NumberFormatException ex) {
                    JOptionPane.showMessageDialog(TabJobSubmission.this,
                            "The Cells Across field on the General tab is invalid.\n" +
                                    "Unable to submit the job.",
                            "Job Submission Error",
                            JOptionPane.ERROR_MESSAGE);
                    return;
                }

                if ((mainData.getValue(TitanConstants.MAX_NUM_TIME_STEPS).compareTo("") == 0) &&
                        (mainData.getValue(TitanConstants.MAX_TIME).compareTo("") == 0)) {
                    JOptionPane.showMessageDialog(TabJobSubmission.this,
                            "Both the Max Number of Time Steps and Max Time fields on the Run Parameters tab are null.\n" +
                                    "Unable to submit the job.",
                            "Job Submission Error",
                            JOptionPane.ERROR_MESSAGE);
                    return;
                }

                if (mainData.getValue(TitanConstants.MAX_NUM_TIME_STEPS).compareTo("") != 0) {
                    try {
                        data.maxNumberTimeSteps = Integer.parseInt(mainData.getValue(TitanConstants.MAX_NUM_TIME_STEPS));
                    } catch (NumberFormatException ex) {
                        JOptionPane.showMessageDialog(TabJobSubmission.this,
                                "The Max Number of Time Steps field on the Run Parameters tab is invalid.\n" +
                                        "Unable to submit the job.",
                                "Job Submission Error",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                }

                if (mainData.getValue(TitanConstants.MAX_TIME).compareTo("") != 0) {
                    try {
                        data.maxTime = Float.parseFloat(mainData.getValue(TitanConstants.MAX_TIME));
                    } catch (NumberFormatException ex) {
                        JOptionPane.showMessageDialog(TabJobSubmission.this,
                                "The Max Time field on the Run Parameters tab is invalid.\n" +
                                        "Unable to submit the job.",
                                "Job Submission Error",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                }

                data.AMR = Boolean.valueOf(mainData.getValue(TitanConstants.AMR)).booleanValue();

                // Method order
                String methodOrder = mainData.getValue(TitanConstants.ORDER_METHOD);
                if (methodOrder.compareTo(TitanConstants.OrderMethod[TitanConstants.ORDER_METHOD_FIRST]) == 0)
                    data.firstOrder = true;
                else if (methodOrder.compareTo(TitanConstants.OrderMethod[TitanConstants.ORDER_METHOD_SECOND]) == 0)
                    data.secondOrder = true;

                TitanScaleData scaleData = simWriter.new TitanScaleData();

                scaleData.scaleSim = Boolean.valueOf(mainData.getValue(TitanConstants.SCALE_SIM)).booleanValue();
                if (scaleData.scaleSim) {
                    try {
                        scaleData.lengthScale = Float.parseFloat(mainData.getValue(TitanConstants.LENGTH_SCALE));

                        if (scaleData.lengthScale <= 0.0) {
                            JOptionPane.showMessageDialog(TabJobSubmission.this,
                                    "The Length Scale field on the Run Parameters tab is less than or equal to 0.0.\n" +
                                            "Unable to submit the job.",
                                    "Job Submission Error",
                                    JOptionPane.ERROR_MESSAGE);
                            return;
                        }
                    } catch (NumberFormatException ex) {
                        JOptionPane.showMessageDialog(TabJobSubmission.this,
                                "The Length Scale field on the Run Parameters tab is invalid.\n" +
                                        "Unable to submit the job.",
                                "Job Submission Error",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }

                    try {
                        scaleData.gravityScale = Float.parseFloat(mainData.getValue(TitanConstants.GRAVITY_SCALE));

                        if (scaleData.gravityScale == 0.0) {
                            JOptionPane.showMessageDialog(TabJobSubmission.this,
                                    "The Gravity Scale field on the Run Parameters tab is equal to 0.0.\n" +
                                            "Unable to submit the job.",
                                    "Job Submission Error",
                                    JOptionPane.ERROR_MESSAGE);
                            return;
                        }
                    } catch (NumberFormatException ex) {
                        JOptionPane.showMessageDialog(TabJobSubmission.this,
                                "The Gravity Scale field on the Run Parameters tab is invalid.\n" +
                                        "Unable to submit the job.",
                                "Job Submission Error",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }

                    if (mainData.getValue(TitanConstants.HEIGHT_SCALE).compareTo("") != 0) {

                        try {
                            scaleData.heightScale = Float.parseFloat(mainData.getValue(TitanConstants.HEIGHT_SCALE));

                            if (scaleData.heightScale <= 0.0) {
                                JOptionPane.showMessageDialog(TabJobSubmission.this,
                                        "The Height Scale field on the Run Parameters tab is less than or equal to 0.0.\n" +
                                                "Unable to submit the job.",
                                        "Job Submission Error",
                                        JOptionPane.ERROR_MESSAGE);
                                return;
                            }
                        } catch (NumberFormatException ex) {
                            JOptionPane.showMessageDialog(TabJobSubmission.this,
                                    "The Height Scale field on the Run Parameters tab is invalid.\n" +
                                            "Unable to submit the job.",
                                    "Job Submission Error",
                                    JOptionPane.ERROR_MESSAGE);
                            return;
                        }
                    }
                }

                // Stat Props
                if (mainData.getValue(TitanConstants.RUN_ID).compareTo("") != 0) {
                    try {
                        data.runID = Integer.parseInt(mainData.getValue(TitanConstants.RUN_ID));
                    } catch (NumberFormatException ex) {
                        JOptionPane.showMessageDialog(TabJobSubmission.this,
                                "The Stat Props Run ID field on the Run Parameters tab is invalid.\n" +
                                        "Unable to submit the job.",
                                "Job Submission Error",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                }

                if (mainData.getValue(TitanConstants.FLOW_OUTLINE_HGT).compareTo("") != 0) {
                    try {
                        data.flowOutlineHeight = Float.parseFloat(mainData.getValue(TitanConstants.FLOW_OUTLINE_HGT));
                    } catch (NumberFormatException ex) {
                        JOptionPane.showMessageDialog(TabJobSubmission.this,
                                "The Stat Props Height used to define flow outline (>0) field on the Run Parameters tab is invalid.\n" +
                                        "Unable to submit the job.",
                                "Job Submission Error",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                }

                if (mainData.getValue(TitanConstants.TEST_FLOW_HEIGHT_MIN).compareTo("") != 0) {
                    try {
                        data.flowHeight = Float.parseFloat(mainData.getValue(TitanConstants.TEST_FLOW_HEIGHT_MIN));
                    } catch (NumberFormatException ex) {
                        JOptionPane.showMessageDialog(TabJobSubmission.this,
                                "The Stat Props Test if flow Reaches Height field on the Run Parameters tab is invalid.\n" +
                                        "Unable to submit the job.",
                                "Job Submission Error",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                }

                if (mainData.getValue(TitanConstants.TEST_FLOW_X_LOC).compareTo("") != 0) {
                    try {
                        data.testX = Float.parseFloat(mainData.getValue(TitanConstants.TEST_FLOW_X_LOC));
                    } catch (NumberFormatException ex) {
                        JOptionPane.showMessageDialog(TabJobSubmission.this,
                                "The Stat Props Flow Height Test Point X Location field on the Run Parameters tab is invalid.\n" +
                                        "Unable to submit the job.",
                                "Job Submission Error",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                }

                if (mainData.getValue(TitanConstants.TEST_FLOW_Y_LOC).compareTo("") != 0) {
                    try {
                        data.testY = Float.parseFloat(mainData.getValue(TitanConstants.TEST_FLOW_Y_LOC));
                    } catch (NumberFormatException ex) {
                        JOptionPane.showMessageDialog(TabJobSubmission.this,
                                "The Stat Props Flow Height Test Point Y Location field on the Run Parameters tab is invalid.\n" +
                                        "Unable to submit the job.",
                                "Job Submission Error",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                }

                data.statPropsPrefix = mainData.getValue(TitanConstants.STAT_PROPS_PREFIX);

                // Outline Props
                data.outlinePropsEnabled = Boolean.valueOf(mainData.getValue(TitanConstants.OUTLINE_PROPS_ENABLED)).booleanValue();

                if (mainData.getValue(TitanConstants.MAX_LINEAR_SIZE).compareTo("") != 0) {
                    try {
                        data.maxLinearSize = Float.parseFloat(mainData.getValue(TitanConstants.MAX_LINEAR_SIZE));
                    } catch (NumberFormatException ex) {
                        JOptionPane.showMessageDialog(TabJobSubmission.this,
                                "The Outline Props Max Linear Size field on the Run Parameters tab is invalid.\n" +
                                        "Unable to submit the job.",
                                "Job Submission Error",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                }

                data.initSize = mainData.getValue(TitanConstants.INIT_SIZE);

                data.outlinePropsPrefix = mainData.getValue(TitanConstants.OUTLINE_PROPS_PREFIX);

                // Output Options Tab

                data.overwriteOutput = Boolean.valueOf(mainData.getValue(TitanConstants.OVERWRITE_OUTPUT)).booleanValue();

                data.restartOutputEnabled = Boolean.valueOf(mainData.getValue(TitanConstants.RESTART_OUTPUT_ENABLED)).booleanValue();

                if (mainData.getValue(TitanConstants.RESULT_OUTPUT_TIME_DELTA1).compareTo("") != 0) {
                    try {
                        data.resultOutputDelta1 = Float.parseFloat(mainData.getValue(TitanConstants.RESULT_OUTPUT_TIME_DELTA1));
                    } catch (NumberFormatException ex) {
                        JOptionPane.showMessageDialog(TabJobSubmission.this,
                                "The Restart Time Between Results Output field on the Output Options tab is invalid.\n" +
                                        "Unable to submit the job.",
                                "Job Submission Error",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                }

                if (mainData.getValue(TitanConstants.SAVE_TIME_DELTA1).compareTo("") != 0) {
                    try {
                        data.saveDelta1 = Integer.parseInt(mainData.getValue(TitanConstants.SAVE_TIME_DELTA1));
                    } catch (NumberFormatException ex) {
                        JOptionPane.showMessageDialog(TabJobSubmission.this,
                                "The Restart Time Between Saves field on the Output Options tab is invalid.\n" +
                                        "Unable to submit the job.",
                                "Job Submission Error",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                }

                data.keepAll = Boolean.valueOf(mainData.getValue(TitanConstants.KEEP_ALL)).booleanValue();
                data.keepRedundantData = Boolean.valueOf(mainData.getValue(TitanConstants.KEEP_REDUNDANT_DATA)).booleanValue();
                data.outputPrefix1 = mainData.getValue(TitanConstants.OUTPUT_PREFIX1);

                // Output types
                String visType = mainData.getValue(TitanConstants.VIS_OUTPUT);
                if (visType.indexOf(TitanConstants.VizTypes[TitanConstants.VIS_OUTPUT_TECPLOT]) != -1)
                    data.tecplot = true;
                if (visType.indexOf(TitanConstants.VizTypes[TitanConstants.VIS_OUTPUT_MSHPLOT]) != -1)
                    data.mshplot = true;
                if (visType.indexOf(TitanConstants.VizTypes[TitanConstants.VIS_OUTPUT_XDMF]) != -1)
                    data.xdmf = true;
                if (visType.indexOf(TitanConstants.VizTypes[TitanConstants.VIS_OUTPUT_GRASS]) != -1)
                    data.grass = true;
                if (visType.indexOf(TitanConstants.VizTypes[TitanConstants.VIS_OUTPUT_WEBVIZ]) != -1)
                    data.webviz = true;
                if (visType.indexOf(TitanConstants.VizTypes[TitanConstants.VIS_OUTPUT_GMFGVIZ]) != -1)
                    data.gmfg = true;

                if (mainData.getValue(TitanConstants.RESULT_OUTPUT_TIME_DELTA2).compareTo("") != 0) {
                    try {
                        data.resultOutputDelta2 = Float.parseFloat(mainData.getValue(TitanConstants.RESULT_OUTPUT_TIME_DELTA2));
                    } catch (NumberFormatException ex) {
                        JOptionPane.showMessageDialog(TabJobSubmission.this,
                                "The Time Series Time Between Results Output field on the Output Options tab is invalid.\n" +
                                        "Unable to submit the job.",
                                "Job Submission Error",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                }

                if (mainData.getValue(TitanConstants.SAVE_TIME_DELTA2).compareTo("") != 0) {
                    try {
                        data.saveDelta2 = Integer.parseInt(mainData.getValue(TitanConstants.SAVE_TIME_DELTA2));
                    } catch (NumberFormatException ex) {
                        JOptionPane.showMessageDialog(TabJobSubmission.this,
                                "The Time Series Time Between Saves field on the Output Options tab is invalid.\n" +
                                        "Unable to submit the job.",
                                "Job Submission Error",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                }

                data.outputPrefix2 = mainData.getValue(TitanConstants.OUTPUT_PREFIX2);

                // Material Model and Map Tab

                // Material map material model datum
                data.physicsModel = (mainData.getValue(TitanConstants.PHYSICS_MODEL));

                // Send pile vol_fract field for the TwoPhases-Pitman-Le model only
                if (TitanConstants.PhysicsModels[TitanConstants.PHYSICS_MODEL_TWOPHASES_PITMAN_LE].compareTo(data.physicsModel) == 0) {
                    TWOPHASES_PITMAN_LE_MODEL = true;
                }

                data.useMaterialMap = Boolean.valueOf(mainData.getValue(TitanConstants.USE_GIS_MAT_MAP)).booleanValue();

                // Check and write the physics model's parameters to data
                int numParms = 0;
                int physicsModelIndex = 0;
                boolean useMaterialMapCheck;
                float floatVal;

                try {
                    numParms = Integer.parseInt(mainData.getValue(TitanConstants.NUM_PHYSICS_MODEL_PARAMETERS));
                    data.numPhysicsModelParameters = numParms;
                } catch (NumberFormatException ex) {
                    JOptionPane.showMessageDialog(TabJobSubmission.this,
                            "The derived number of material model parameters is invalid.\n" +
                                    "Unable to submit the job.",
                            "Job Submission Error",
                            JOptionPane.ERROR_MESSAGE);
                    return;
                }

                for (int i = 0; i < TitanConstants.PhysicsModels.length; i++) {
                    if (data.physicsModel.compareTo(TitanConstants.PhysicsModels[i]) == 0) {
                        physicsModelIndex = i;
                        break;
                    }
                }

                if (data.useMaterialMap == true) {

                    useMaterialMapCheck = matModelMap.setUseMaterialMap(physicsModelIndex, false);

                    if (useMaterialMapCheck == false) {
                        int selectionOption = JOptionPane.showConfirmDialog(TabJobSubmission.this,
                                "Material Map Error detected.\n" +
                                        "Select the OK option to continue with Use Material Map set to False.",
                                "Job Submission Option",
                                JOptionPane.OK_CANCEL_OPTION);

                        if (selectionOption == JOptionPane.CANCEL_OPTION) {

                            return;
                        }
                    }
                }

                String[] matParmNames = TitanConstants.matParmNames[physicsModelIndex];
                String matParmName;
                String[] matParmUnits = TitanConstants.matParmUnits[physicsModelIndex];
                String matParmUnit;
                String label;

                int index;
                for (int i = 0; i < numParms; i++) {
                    if (physicsModelIndex == TitanConstants.PHYSICS_MODEL_COULOMB &&
                            data.useMaterialMap == true) {
                        if (i == TitanConstants.INT_FRICT_INDEX) {
                            matParmName = matParmNames[TitanConstants.INT_FRICT_INDEX];
                            matParmUnit = matParmUnits[TitanConstants.INT_FRICT_INDEX];
                            label = matParmName + " " + matParmUnit;
                        } else {
                            matParmName = matParmNames[TitanConstants.BED_FRICT_INDEX];
                            matParmUnit = matParmUnits[TitanConstants.BED_FRICT_INDEX];
                            label = TitanConstants.MATERIAL_MAP_NAMES[i - 1] + " " + matParmName + " " + matParmUnit;
                        }
                    } else {
                        matParmName = matParmNames[i];
                        matParmUnit = matParmUnits[i];
                        label = matParmName + " " + matParmUnit;
                    }
                    try {
                        // Check if valid float
                        floatVal = Float.parseFloat(mainData.getValue("MATERIAL_MODEL_TEXT_INPUT_" + i));
                        data.matParmsVals[i] = floatVal;
                    } catch (NumberFormatException ex) {
                        JOptionPane.showMessageDialog(TabJobSubmission.this,
                                "The " + label + " field on the Material Model and Map tab is invalid.\n" +
                                        "Unable to submit the job.",
                                "Job Submission Error",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                }

                data.stoppingCriteria = mainData.getValue(TitanConstants.STOPPING_CRITERIA);

                // Piles Tab

                // Check for presence of material sources
                if (pileData.length == 0 && fluxData.length == 0) {
                    JOptionPane.showMessageDialog(TabJobSubmission.this,
                            "There are no piles or flux sources defined.\n" +
                                    "Unable to submit the job.",
                            "Job Submission Error",
                            JOptionPane.ERROR_MESSAGE);
                    return;
                }

                data.pile = new TitanSimulationPile[pileData.length];
                for (int i = 0; i < pileData.length; i++) {
                    try {
                        data.pile[i] = simWriter.new TitanSimulationPile();
                        data.pile[i].type = Float.parseFloat(pileData[i].getValue(TitanConstants.PILE_TYPE));
                        data.pile[i].maxThickness = Float.parseFloat(pileData[i].getValue(TitanConstants.PILE_MAX_INIT_THICKNESS));
                        data.pile[i].centerVolumeX = Float.parseFloat(pileData[i].getValue(TitanConstants.PILE_CENTER_INIT_VOLUME_XC));
                        data.pile[i].centerVolumeY = Float.parseFloat(pileData[i].getValue(TitanConstants.PILE_CENTER_INIT_VOLUME_YC));

                        // Allow a blank field for volume fraction when physics model is not TwoPhases-Pitman-Le
                        if (TWOPHASES_PITMAN_LE_MODEL == true) {
                            data.pile[i].volumeFraction = Float.parseFloat(pileData[i].getValue(TitanConstants.PILE_VOLUME_FRACTION));
                        }

                        data.pile[i].majorExtent = Float.parseFloat(pileData[i].getValue(TitanConstants.PILE_MAJOR_EXTENT));
                        data.pile[i].minorExtent = Float.parseFloat(pileData[i].getValue(TitanConstants.PILE_MINOR_EXTENT));
                        data.pile[i].orientation = Float.parseFloat(pileData[i].getValue(TitanConstants.PILE_ORIENTATION_ANGLE));
                        data.pile[i].speed = Float.parseFloat(pileData[i].getValue(TitanConstants.PILE_INIT_SPEED));
                        data.pile[i].direction = Float.parseFloat(pileData[i].getValue(TitanConstants.PILE_INIT_DIRECTION));
                    } catch (NumberFormatException ex) {
                        JOptionPane.showMessageDialog(TabJobSubmission.this,
                                "There is a problem with the piles data.\n" +
                                        "Please verify no blank fields on the Piles tab;\n" +
                                        "when modifying a field, click enter <Enter> or another field.\n" +
                                        "Unable to submit the job.",
                                "Job Submission Error",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                }

                // Flux Sources Tab

                data.fluxSource = new TitanSimulationFluxSource[fluxData.length];
                for (int i = 0; i < fluxData.length; i++) {
                    try {
                        data.fluxSource[i] = simWriter.new TitanSimulationFluxSource();
                        data.fluxSource[i].centerX = Float.parseFloat(fluxData[i].getValue(TitanConstants.SRC_CENTER_XC));
                        data.fluxSource[i].centerY = Float.parseFloat(fluxData[i].getValue(TitanConstants.SRC_CENTER_YC));
                        data.fluxSource[i].majorExtent = Float.parseFloat(fluxData[i].getValue(TitanConstants.SRC_MAJOR_EXTENT));
                        data.fluxSource[i].minorExtent = Float.parseFloat(fluxData[i].getValue(TitanConstants.SRC_MINOR_EXTENT));
                        data.fluxSource[i].orientationAngle = Float.parseFloat(fluxData[i].getValue(TitanConstants.SRC_ORIENTATION_ANGLE));
                        data.fluxSource[i].initSpeed = Float.parseFloat(fluxData[i].getValue(TitanConstants.SRC_INIT_SPEED));
                        data.fluxSource[i].initDirection = Float.parseFloat(fluxData[i].getValue(TitanConstants.SRC_INIT_DIRECTION));
                        data.fluxSource[i].extrusionFluxRate = Float.parseFloat(fluxData[i].getValue(TitanConstants.SRC_EXTRUSION_FLUX_RATE));
                        data.fluxSource[i].activeTimeStart = Float.parseFloat(fluxData[i].getValue(TitanConstants.SRC_ACTIVE_TIME_START));
                        data.fluxSource[i].activeTimeEnd = Float.parseFloat(fluxData[i].getValue(TitanConstants.SRC_ACTIVE_TIME_END));
                    } catch (NumberFormatException ex) {
                        JOptionPane.showMessageDialog(TabJobSubmission.this,
                                "There is a problem with the flux sources data.\n" +
                                        "Please verify no blank fields on the Flux Sources tab;\n" +
                                        "when modifying a field, click enter <Enter> or another field.\n" +
                                        "Unable to submit the job.",
                                "Job Submission Error",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                }

                if (pileData.length == 0) {
                    boolean foundTimeStartZero = false;
                    for (int i = 0; i < fluxData.length; i++) {
                        if (data.fluxSource[i].activeTimeStart == 0.0f) {
                            foundTimeStartZero = true;
                        }
                    }
                    if (foundTimeStartZero == false) {
                        JOptionPane.showMessageDialog(TabJobSubmission.this,
                                "There is a problem with the flux sources data.\n" +
                                        "At least one flux source active time start must be 0.0 on the Flux Sources tab when no piles are defined.\n" +
                                        "Unable to submit the job.",
                                "Job Submission Error",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                }

                // Discharge Planes Tab

                data.plane = new TitanSimulationPlane[planeData.length];
                for (int i = 0; i < planeData.length; i++) {
                    try {
                        data.plane[i] = simWriter.new TitanSimulationPlane();
                        data.plane[i].a_e = Float.parseFloat(planeData[i].getValue(TitanConstants.PLANE_A_E));
                        data.plane[i].a_n = Float.parseFloat(planeData[i].getValue(TitanConstants.PLANE_A_N));
                        data.plane[i].b_e = Float.parseFloat(planeData[i].getValue(TitanConstants.PLANE_B_E));
                        data.plane[i].b_n = Float.parseFloat(planeData[i].getValue(TitanConstants.PLANE_B_N));
                    } catch (NumberFormatException ex) {
                        JOptionPane.showMessageDialog(TabJobSubmission.this,
                                "There is a problem with the discharge planes data.\n" +
                                        "Please verify no blank fields on the Discharge Planes tab;\n" +
                                        "when modifying a field, click enter <Enter> or another field.\n" +
                                        "Unable to submit the job.",
                                "Job Submission Error",
                                JOptionPane.ERROR_MESSAGE);
                        return;
                    }
                }

                // Create the titan python input file with the verified data and scaleData.
                // Note, on restarts do not need to create this file
                simWriter.writeSimulationData(data, scaleData);

                // Create files required for creating KML files when the KML button is clicked
                // on the JobDetailsDialog tab

                // Create the zone.txt file
                int zone = 0;
                boolean zoneFound = false;
                String cellhdmapfile;
                if (data.format.compareTo(TitanConstants.GisFormats[TitanConstants.GIS_FORMAT_GIS_GRASS]) == 0) {
                    cellhdmapfile =
                            mainData.getValue(TitanConstants.GIS_INFO_DIRECTORY) + File.separator +
                                    mainData.getValue(TitanConstants.GIS_SUBDIR) + File.separator +
                                    mainData.getValue(TitanConstants.GIS_MAPSET) + File.separator +
                                    "cellhd" + File.separator + mainData.getValue(TitanConstants.GIS_MAP);
                } else {
                    // For GDAL, GIS Map is a full path name
                    cellhdmapfile = mainData.getValue(TitanConstants.GIS_MAP);
                }
                BufferedWriter writer;
                String fn;

                if (data.zoneOverride == 0) {
                    // Read the zone from the GIS GRASS map file in the cellhd directory

                    try {
                        BufferedReader reader = new BufferedReader(new FileReader(cellhdmapfile));
                        String line;

                        // Look for zone:
                        Pattern p = Pattern.compile("zone:");

                        while ((line = reader.readLine()) != null) {
                            //System.out.println("line: " + line);
                            Matcher m = p.matcher(line);
                            if (m.find()) {
                                //System.out.println("zone: found");
                                try {
                                    int startIndex = line.lastIndexOf(":") + 1;
                                    //System.out.println("startIndex: " + startIndex);
                                    int stopIndex = line.length();
                                    //System.out.println("stopIndex: " + stopIndex);
                                    String sub = line.substring(startIndex, stopIndex);
                                    //System.out.println("sub.trim(): " + sub.trim());
                                    zone = Integer.parseInt(sub.trim());
                                    //System.out.println("zone: " + zone);
                                    zoneFound = true;
                                } catch (NumberFormatException ex) {
                                    zoneFound = false;
                                }
                                break;
                            }
                        }

                        reader.close();
                    } catch (IOException ex) {
                        System.out.println("Error reading: " + cellhdmapfile);
                        zoneFound = false;
                    }
                } else {
                    zone = data.zoneOverride;
                    zoneFound = true;
                }

                if (zoneFound == true) {

                    fn = new String(parmDir + File.separator + "zone.txt");
                    try {
                        writer = new BufferedWriter(new FileWriter(fn));
                        writer.write("" + zone + "\n");
                        writer.write("" + data.hemisphere + "\n");
                        writer.close();

                    } catch (IOException ex) {
                        int selectionOption = JOptionPane.showConfirmDialog(TabJobSubmission.this,
                                "Unable to create the zone.txt file required for creating KML files.\n\n" +
                                        "Select the OK Option to continue job submission.",
                                "Job Submission Option",
                                JOptionPane.OK_CANCEL_OPTION);

                        if (selectionOption != JOptionPane.OK_OPTION) {

                            return;

                        }
                    }
                } else {

                    int selectionOption = JOptionPane.showConfirmDialog(TabJobSubmission.this,
                            "Unable to create the zone.txt file required for creating KML files.\n" +
                                    "Please verify information on the GIS tab.\n" +
                                    "Please verify access to " + cellhdmapfile +".\n\n" +
                                    "Select the OK Option to continue job submission.",
                            "Job Submission Option",
                            JOptionPane.OK_CANCEL_OPTION);

                    if (selectionOption != JOptionPane.OK_OPTION) {

                        return;

                    }
                }

                // Previous versions of titan created a scale.data file containing the height scale.
                // Python files for creating the KML files need the height scale
                fn = new String(parmDir + File.separator + "height_scale_for_KML.data");

                try {
                    writer = new BufferedWriter(new FileWriter(fn));
                    writer.write("" + scaleData.heightScale + "\n");
                    writer.close();
                } catch (IOException ex) {
                    int selectionOption = JOptionPane.showConfirmDialog(TabJobSubmission.this,
                            "Unable to create the height_scale_for_KML.data file required for creating KML files.\n\n" +
                                    "Select the OK Option to continue job submission.",
                            "Job Submission Option",
                            JOptionPane.OK_CANCEL_OPTION);

                    if (selectionOption != JOptionPane.OK_OPTION) {

                        return;

                    }
                }
            }

            // Busy but all functionality needs to be done in the event handler

            // Set up for running the job

            // Create the run method

            RunMethod runMethod = null;

            submitContainer.createRunMethod(runStyle.getValue(), restartEnabled, parmDir, dateString);
            runMethod = submitContainer.getRunMethod();

            // Create the job monitor

            TitanJobMonitor jobMon = dba.new TitanJobMonitor();

            jobMon.submitMethod = new String(runStyle.getValue());

            JobSubmissionParameters jobParms = new JobSubmissionParameters();
            jobParms.addJobParameter(JobSubmissionParameters.QUEUE, queueNamePBS.getValue());
            jobParms.addJobParameter(JobSubmissionParameters.MAX_RUN_TIME, maxRuntimePBS.getValue());
            jobParms.addJobParameter(JobSubmissionParameters.RUN_MODE, runMode.getValue());
            jobParms.addJobParameter(JobSubmissionParameters.RUN_ARGS, runArgs.getValue());

            int numNodes = 1;
            int numCPUs = 1;
            int numNodesSubmit = 1;
            int numCPUsSubmit = 1;
            if (runMethod.allowMultiProcessors()) {
                try {
                    //check if SUBMIT is selected
                    if (runStyle.getValue().compareTo(JobSubmissionContainer.RUN_STYLE_SUBMIT) == 0) {
                        numNodesSubmit = Integer.parseInt(nodesSubmit.getValue());
                        numCPUsSubmit = Integer.parseInt(cpusSubmit.getValue());
                        if ((numNodesSubmit > 64) || (numNodesSubmit <= 0)) {
                            JOptionPane.showMessageDialog(TabJobSubmission.this,
                                    "The Number of Nodes is invalid.\n" +
                                            "Unable to submit the job.",
                                    "Job Submission Error",
                                    JOptionPane.ERROR_MESSAGE);
                            return;
                        }
                        if ((numCPUsSubmit <= 0) || (numCPUsSubmit > 8)) {
                            JOptionPane.showMessageDialog(TabJobSubmission.this,
                                    "The Number of CPUs is invalid.\n" +
                                            "Unable to submit the job.",
                                    "Job Submission Error",
                                    JOptionPane.ERROR_MESSAGE);
                            return;
                        }
                        if (!maxRuntimeSubmit.getValue().matches("\\d\\d\\:\\d\\d\\:\\d\\d")) {
                            JOptionPane.showMessageDialog(TabJobSubmission.this,
                                    "Format of Max time is invalid\n" +
                                            "hh:mm:ss",
                                    "Job Submission Error",
                                    JOptionPane.ERROR_MESSAGE);
                            return;
                        }
                        String tmpRunTime = maxRuntimeSubmit.getValue();
                        if (tmpRunTime.charAt(0) >= '7') {
                            if (tmpRunTime.charAt(1) > '2') {
                                JOptionPane.showMessageDialog(TabJobSubmission.this,
                                        "Max time is greater than 72 hrs.",
                                        "Job Submission Error",
                                        JOptionPane.ERROR_MESSAGE);
                                return;
                            }
                        }
                    } else {
                        numNodes = Integer.parseInt(nodesLocal.getValue());
                        numCPUs = Integer.parseInt(cpusLocal.getValue());
                    }
                } catch (NumberFormatException ex) {
                    JOptionPane.showMessageDialog(TabJobSubmission.this,
                            "The Number or Nodes or CPUs is invalid.\n" +
                                    "Unable to submit the job.",
                            "Job Submission Error",
                            JOptionPane.ERROR_MESSAGE);
                    return;
                }
            } // end if runMethod.allowMultiProcessors()

            jobParms.addJobParameter(JobSubmissionParameters.MAX_RUN_TIME_SUBMIT, maxRuntimeSubmit.getValue());
            jobParms.addJobParameter(JobSubmissionParameters.NODES, Integer.toString(numNodes));
            jobParms.addJobParameter(JobSubmissionParameters.CPUS, Integer.toString(numCPUs));
            jobParms.addJobParameter(JobSubmissionParameters.NODES_SUBMIT, Integer.toString(numNodesSubmit));
            jobParms.addJobParameter(JobSubmissionParameters.CPUS_SUBMIT, Integer.toString(numCPUsSubmit));
            jobParms.addJobParameter(JobSubmissionParameters.REQUIREMENTS, requirementsCondor.getValue());

            // Check for irods

            String useIrods = mainData.getValue("IRODS_OUTPUT");
            if (useIrods.indexOf("USE_IRODS") != -1) {
                jobParms.addJobParameter(JobSubmissionParameters.OUTPUT_IRODS, "true");
                jobParms.addJobParameter(JobSubmissionParameters.DATE, dateString);
                jobParms.addJobParameter(JobSubmissionParameters.SAVE_DIR, mainData.getValue(TitanConstants.RUN_NAME));
            } else {
                jobParms.addJobParameter(JobSubmissionParameters.OUTPUT_IRODS, "false");
            }

            RunMethod.GisParms gParms = runMethod.new GisParms();

            gParms.NUM_CELLS_ACROSS = mainData.getValue(TitanConstants.NUM_CELLS_ACROSS);
            gParms.GIS_FORMAT = mainData.getValue(TitanConstants.GIS_FORMAT);
            gParms.GIS_INFO_DIRECTORY = mainData.getValue(TitanConstants.GIS_INFO_DIRECTORY);
            gParms.GIS_SUBDIR = mainData.getValue(TitanConstants.GIS_SUBDIR);
            gParms.GIS_MAPSET = mainData.getValue(TitanConstants.GIS_MAPSET);
            gParms.GIS_MAP = mainData.getValue(TitanConstants.GIS_MAP);
            gParms.MIN_X_LOC = mainData.getValue(TitanConstants.MIN_X_LOC);
            gParms.MIN_Y_LOC = mainData.getValue(TitanConstants.MIN_Y_LOC);
            gParms.MAX_X_LOC = mainData.getValue(TitanConstants.MAX_X_LOC);
            gParms.MAX_Y_LOC = mainData.getValue(TitanConstants.MAX_Y_LOC);
            gParms.GIS_VECTOR = mainData.getValue(TitanConstants.GIS_VECTOR);
            gParms.HEMISPHERE = mainData.getValue(TitanConstants.HEMISPHERE);
            gParms.ZONEOVERRIDE = mainData.getValue(TitanConstants.ZONEOVERRIDE);

            runMethod.initialize(
                    gParms,
                    optionsManager.getBinaryMode(),
                    optionsManager.getBinaryDirectory(),
                    additional.getValue(),
                    jobParms,
                    restartFile,
                    restartMaxNumberTimeSteps,
                    restartMaxTime);

            // Run the job

            JobDetails jobDetails = runMethod.runJob();

            jobMon.dateSubmitted =
                    (jobDetails.getSubmissionDate().get(Calendar.MONTH) + 1) + "/" +
                            jobDetails.getSubmissionDate().get(Calendar.DAY_OF_MONTH) + "/" +
                            jobDetails.getSubmissionDate().get(Calendar.YEAR) + "  " +
                            jobDetails.getSubmissionDate().get(Calendar.HOUR_OF_DAY) + ":" +
                            jobDetails.getSubmissionDate().get(Calendar.MINUTE);
            jobMon.jobName = jobDetails.getJobName();
            jobMon.jobID = jobDetails.getJobID();
            jobMon.outputDirectory = jobDetails.getOutputDirectory();

            try {
                jobMon.submitHost = InetAddress.getLocalHost().getHostName();
            } catch (UnknownHostException ex) {
                jobMon.submitHost = new String("Unknown");
            }

            submitSummary.append("Job Submitted :\n");
            if (jobDetails.getErrorMessage() == null) {
                submitSummary.append("  Job Name    : " + jobDetails.getJobName() + "\n");
                submitSummary.append("  Job ID      : " + jobDetails.getJobID() + "\n");
                submitSummary.append("  Output Directory : " + jobDetails.getOutputDirectory() + "\n\n");

                //only do this if submit is called
                if (runStyle.getValue().compareTo(JobSubmissionContainer.RUN_STYLE_SUBMIT) == 0) {
                    Process job = jobDetails.getProcess();
                    Thread thread = this.new SubmitOutput(parmDir, jobDetails.getJobName(), job, submitSummary);
                    thread.start();
                } else {
                    submitSummary.append("----------------------------------------\n");
                }

            } else {
                submitSummary.append("  ERROR : " + jobDetails.getErrorMessage() + "\n");
                submitSummary.append("----------------------------------------\n");
            }

            //submitSummary.append("----------------------------------------\n");

            try {
                dba.addJob(jobMon);
            } catch (SQLException ex) {
                JOptionPane.showMessageDialog(TabJobSubmission.this,
                        ex.getMessage(),
                        "Titan Database Error",
                        JOptionPane.WARNING_MESSAGE);
            }
        }

        //------
        private class SubmitOutput extends Thread {

            public Process job;
            public JTextArea text;
            public String baseDir;
            public String jobName;
            private BufferedWriter writer;

            public SubmitOutput(String parmDir, String jobDetailsName, Process p, JTextArea t) {
                job = p;
                text = t;
                baseDir = new String(parmDir);
                jobName = jobDetailsName;
            }

            public void run() {

                try {
                    BufferedReader br = new BufferedReader(new InputStreamReader(job.getInputStream()));
                    String fn = new String(baseDir + File.separator + JobSubmissionContainer.HUB_SUBMIT_INFO);
                    writer = new BufferedWriter(new FileWriter(fn));

                    String line = "";
                    boolean getSlurmJobID = true;
                    int indexOfLeftParen;
                    int indexOfRightParen;
                    String slurmJobID;

                    text.append("Submit Output Stream :\n");
                    writer.write("Submit Output Stream :\n");
                    writer.write("Job Name: " + jobName + "\n");

                    do {
                        try {
                            line = (String) br.readLine();
                        } catch (IOException ex) {
                            text.append("IO Exception on Buffered Read");
                            writer.write("IO Exception on Buffered Read");
                        }

                        if (line != null) {
                            if (getSlurmJobID) {
                                getSlurmJobID = false;
                                // The Slurm Job ID is returned in the message as (Slurm Job ID) message
                                indexOfLeftParen = line.indexOf("(", 0);
                                indexOfRightParen = line.indexOf(")", indexOfLeftParen);
                                slurmJobID = line.substring(indexOfLeftParen + 1, indexOfRightParen);
                                text.append("Slurm Job ID: " + slurmJobID + "\n");
                                writer.write("Slurm Job ID: " + slurmJobID + "\n");
                            }
                            text.append("  " + line + "\n");
                            writer.write("  " + line + "\n");
                        }

                    } while (line != null);

                    writer.close();

                    submitSummary.append("----------------------------------------\n");
                } catch (IOException ex) {
                    ex.printStackTrace();
                }
            }
        }
    }
}
