package titan.gui;


import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.FileNotFoundException;
import java.net.URL;
import java.sql.SQLException;

import javax.help.CSH;
import javax.help.HelpBroker;
import javax.help.HelpSet;
import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;

import titan.io.INameValuePair;
import titan.io.NameValuePairGroup;
import titan.io.NameValuePairSimple;
import titan.io.ParameterInput;
import titan.io.ParameterListInput;
import titan.io.ParameterListOutput;
import titan.io.ParameterOutput;
import titan.io.TitanDBAccess;
import titan.jobs.JobSubmissionContainer;
import titan.options.OptionsDialog;
import titan.options.OptionsManager;

public class Titan extends JPanel {

    private static String APP_NAME = "Titan";

    private static TabFileIO fileIO;
    private static TabGIS gis;
    private static TabRunParameters runParameters;
    private static TabOutputOptions outputOptions;
    private static TabMaterialModelMap matModelMap;
    private static TabPiles piles;
    private static TabFluxSources flux;
    private static TabDischargePlanes discharge;
    private static TabJobSubmission jobSubmit;
    private static TabJobMonitor jobMonitor;
    private static DefaultTitanData defaultData;
    private static TitanDBAccess dba;
    private static OptionsManager optionsManager;
    private static JobSubmissionContainer submitContainer;

    public Titan() {

        // Note: on vhub only 1 node and 1 cpu are allowed per job.
        // Default, assume running titan2d on VHUB
        boolean eVHUB = true;
        String eInputDir = "";

        setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));

        String user = System.getProperty("user.home");

        if (user == null) {
            JOptionPane.showMessageDialog(this, "Unable to determine users home directory.\n" +
                    "Cannot create job monitoring database.");
        }
        File titanDir = new File(user + File.separator + ".titan");
        if (titanDir.exists() == false) {
            titanDir.mkdir();
        }

        optionsManager = new OptionsManager();

        createToolbar();

        try {

            try {

                String eVHUBStr = System.getenv("E_VHUB");
                if (eVHUBStr != null) {

                    try {
                        eVHUB = Boolean.valueOf(eVHUBStr).booleanValue();
                        System.out.println("E_VHUB environment variable eVHUBStr: " + eVHUBStr + " => eVHUB:" + eVHUB);

                    } catch (NumberFormatException ex) {
                        eVHUB = true;
                        JOptionPane.showMessageDialog(this, "The E_VHUB environment variable is defined but not set to \"true\" or \"false\"." +
                                "Setting eVHUB=true.");
                    }
                } else {
                    eVHUB = true;
                    //JOptionPane.showMessageDialog(this, "The E_VHUB environment variable is not defined." +
                    //"Setting eVHUB=true.");
                }
            } catch (Exception e) {

                eVHUB = true;
                JOptionPane.showMessageDialog(this, "The E_VHUB environment variable is not accessible using System.getenv()." +
                        "Setting eVHUB=true.");
            }

            try {
                // Allow user to set the start directory for the directory selectors
                String eInputDirStr = System.getenv("E_INPUTDIR");
                if (eInputDirStr != null) {
                    eInputDir = eInputDirStr;
                    System.out.println("E_INPUTDIR environment variable eInputDirStr: " + eInputDirStr + " => eInputDir:" + eInputDir);
                }
            } catch (Exception e) {

                 JOptionPane.showMessageDialog(this, "The E_INPUTDIR environment variable is not accessible using System.getenv().");
            }

            dba = new TitanDBAccess(titanDir.getPath() + File.separator + "titan_db");
            dba.open();

            submitContainer = new JobSubmissionContainer();

            JTabbedPane tp = new JTabbedPane();
            SaveTitanRun save = new SaveTitanRun();
            LoadTitanRun load = new LoadTitanRun();
            FetchData fetchData = new FetchData();
            fileIO = new TabFileIO(eInputDir, save, load);
            gis = new TabGIS();
            runParameters = new TabRunParameters(eInputDir);
            outputOptions = new TabOutputOptions();
            matModelMap = new TabMaterialModelMap(fetchData);
            piles = new TabPiles(false);
            flux = new TabFluxSources();
            discharge = new TabDischargePlanes();
            jobSubmit = new TabJobSubmission(fetchData, dba, optionsManager, save, matModelMap, submitContainer, eVHUB);
            jobMonitor = new TabJobMonitor(dba, submitContainer);

            tp.addTab("Load/Save", fileIO);
            tp.addTab("GIS", gis);
            tp.addTab("Run Parameters", runParameters);
            tp.addTab("Output Options", outputOptions);
            tp.addTab("Material Model and Map", matModelMap);
            tp.addTab("Piles", piles);
            tp.addTab("Flux Sources", flux);
            tp.addTab("Discharge Planes", discharge);
            tp.addTab("Job Submission", jobSubmit);
            tp.addTab("Job Monitor", jobMonitor);

            tp.setAlignmentX(Component.LEFT_ALIGNMENT);
            add(tp);

            // Set the defaults

            defaultData = new DefaultTitanData();
            gis.setData(defaultData.getDefaults());
            runParameters.setData(defaultData.getDefaults());
            outputOptions.setData(defaultData.getDefaults());
            matModelMap.setData(defaultData.getDefaults());
            piles.setData(new INameValuePair[0][0]);
            flux.setData(new INameValuePair[0][0]);
            discharge.setData(new INameValuePair[0][0]);

        } catch (ClassNotFoundException ex) {
            JOptionPane.showMessageDialog(this,
                    ex.getMessage(),
                    "Titan Initialization Error",
                    JOptionPane.ERROR_MESSAGE);
            Titan.Exit(this, dba, -1);
        } catch (SQLException ex) {
            JOptionPane.showMessageDialog(this,
                    ex.getMessage(),
                    "Titan Initialization Error",
                    JOptionPane.ERROR_MESSAGE);
            Titan.Exit(this, dba, -1);
        }
    }

    public static void updatePileTab(boolean TWOPHASES_COULOMB_MODEL) {
        piles.setHeadings(TWOPHASES_COULOMB_MODEL);
    }

    public static void Exit(Component parent, TitanDBAccess db, int status) {
        try {
            db.close();
        } catch (SQLException ex) {
            JOptionPane.showMessageDialog(parent,
                    ex.getMessage(),
                    "Database Error",
                    JOptionPane.ERROR_MESSAGE);
        }
        System.exit(status);
    }

    public TitanDBAccess getDB() {
        return dba;
    }

    protected void finalize() throws Throwable {
        System.out.println("Finalize!!!!");
        System.out.flush();
    }

    private void createToolbar() {
        JMenuBar menuBar = new JMenuBar();

        // File menu
        JMenu fileMenu = new JMenu("File");
        JMenuItem exitItem = new JMenuItem("Exit");
        exitItem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ev) {
                Titan.Exit(Titan.this, Titan.this.dba, 0);
            }
        });
        fileMenu.add(exitItem);
        menuBar.add(fileMenu);

        // Options menu
        JMenu optionsMenu = new JMenu("Options");
        JMenuItem runStyleItem = new JMenuItem("Run Style");
        runStyleItem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ev) {
                OptionsDialog od = new OptionsDialog(Titan.this, optionsManager);
                od.setVisible(true);
            }
        });
        optionsMenu.add(runStyleItem);
        menuBar.add(optionsMenu);

        HelpSet hs = null;
        ClassLoader cl = this.getClass().getClassLoader();

        try {

            // Help menu
            URL hsURL = HelpSet.findHelpSet(cl, "Titan.hs");
            hs = new HelpSet(null, hsURL);

            System.out.println("HelpSet: " + "Titan.hs" + " found");

            JMenu helpMenu = new JMenu("Help");
            JMenuItem helpItem = new JMenuItem("Help Me");

            HelpBroker hb = hs.createHelpBroker();
            CSH.setHelpIDString(helpItem, "welcome");
            helpItem.addActionListener(new CSH.DisplayHelpFromSource(hb));

            helpMenu.add(helpItem);

            menuBar.add(helpMenu);

        } catch (Exception ee) {
            System.out.println("HelpSet: " + ee.getMessage());
            System.out.println("HelpSet: " + "Titan.hs" + " not found");

            int selectionOption = JOptionPane.showConfirmDialog(Titan.this,
                    "HelpSet: Titan.hs not found. Continue without the Help menu?\n",
                  "Titan Initialization Question",
                    JOptionPane.OK_CANCEL_OPTION);

            if (selectionOption != JOptionPane.OK_OPTION) {
                Titan.Exit(this, dba, -1);
             }
        }

        menuBar.setAlignmentX(Component.LEFT_ALIGNMENT);
        add(menuBar);
    }

    public class SaveTitanRun {
        public boolean Save(String baseDir, String runName) {
            String directory = baseDir + File.separator + runName;
            String ascFilename = directory + File.separator +
                    runName + "." + TitanConstants.PARM_FILE_EXT;

            File dirFile = new File(directory);

            File testFile = dirFile;

            while (testFile != null) {
                if (testFile.exists()) {
                    if (testFile.canWrite() == false) {
                        JOptionPane.showMessageDialog(Titan.this,
                                "The Base Save Directory and/or Run Name fields specified on the Load/Save tab cannot be written to.\n",
                                "Save Error",
                                JOptionPane.ERROR_MESSAGE);
                        return false;
                    }
                    break;
                } else {
                    testFile = new File(testFile.getParent());
                }
            }

            // Set the load parameters to the same as the run just saved
            INameValuePair[] ionvp = new INameValuePair[1];
            ionvp[0] = new NameValuePairSimple(TitanConstants.LOAD_DIRECTORY, directory);
            fileIO.setData(ionvp);

            if (!dirFile.exists()) {
                dirFile.mkdir();
            }

            ParameterOutput output = new ParameterOutput(new File(ascFilename));
            INameValuePair[] data = gis.getData();
            output.writeData(data);
            data = runParameters.getData();
            output.writeData(data);
            data = outputOptions.getData();
            output.writeData(data);
            data = matModelMap.getData();  //*************
            output.writeData(data);
            output.close();

            ascFilename = directory + File.separator + runName +
                    "." + TitanConstants.PILE_FILE_EXT;
            ParameterListOutput outputList = new ParameterListOutput(new File(ascFilename));
            INameValuePair[][] dataList;
            dataList = piles.getData();
            outputList.writeData(dataList);
            outputList.close();

            ascFilename = directory + File.separator + runName +
                    "." + TitanConstants.FLUX_FILE_EXT;
            outputList = new ParameterListOutput(new File(ascFilename));
            dataList = flux.getData();
            outputList.writeData(dataList);
            outputList.close();

            ascFilename = directory + File.separator + runName +
                    "." + TitanConstants.PLANE_FILE_EXT;
            outputList = new ParameterListOutput(new File(ascFilename));
            dataList = discharge.getData();
            outputList.writeData(dataList);
            outputList.close();

            return true;
        }
    }

    public class LoadTitanRun {
        public void Load(String inputDirectory) {
            ParameterInput input = null;
            try {
                // ascprm file
                String ascprmFilename;

                if (inputDirectory.lastIndexOf(File.separator) == -1) {
                    ascprmFilename = inputDirectory + inputDirectory +
                            "." + TitanConstants.PARM_FILE_EXT;
                } else {
                    ascprmFilename = inputDirectory +
                            inputDirectory.substring(inputDirectory.lastIndexOf(File.separator)) +
                            "." + TitanConstants.PARM_FILE_EXT;
                }
                input = new ParameterInput(new File(ascprmFilename));

                INameValuePair nvp[];
                nvp = input.readData();

                // Set the save parameters to the same as the run just loaded
                INameValuePair[] ionvp = new INameValuePair[2];
                ionvp[0] = new NameValuePairSimple(TitanConstants.SAVE_DIRECTORY, inputDirectory.substring(0, inputDirectory.lastIndexOf(File.separator)));
                ionvp[1] = new NameValuePairSimple(TitanConstants.RUN_NAME, inputDirectory.substring(inputDirectory.lastIndexOf(File.separator) + 1));
                fileIO.setData(ionvp);

                // Set parameters on GIS, General and Material Map tabs
                gis.setData(nvp);
                runParameters.setData(nvp);
                outputOptions.setData(nvp);
                matModelMap.setData(nvp);

                // Coordinate directory selectors
                runParameters.setInputDirectory(fileIO.getInputDirectory());

            } catch (FileNotFoundException ex) {
                JOptionPane.showMessageDialog(Titan.this, "The directory chosen does not appear to be a Titan run directory.\n" +
                                "Titan Input File not found.\n" + ex.getMessage(),
                        "Load Error",
                        JOptionPane.ERROR_MESSAGE);
                gis.setData(defaultData.getDefaults());
                runParameters.setData(defaultData.getDefaults());
                outputOptions.setData(defaultData.getDefaults());
                matModelMap.setData(defaultData.getDefaults());
                piles.setData(new INameValuePair[0][0]);
                flux.setData(new INameValuePair[0][0]);
                discharge.setData(new INameValuePair[0][0]);
                // Failure to find the ASCPRM file is a terminal failure.
                // Do not attempt to read any of the other input files.
                return;
            }

            input.close();

            ParameterListInput inputList = null;
            INameValuePair nvpList[][];

            // ascflux file
            String ascflxFilename = inputDirectory +
                    inputDirectory.substring(inputDirectory.lastIndexOf(File.separator)) +
                    "." + TitanConstants.FLUX_FILE_EXT;
            try {
                inputList = new ParameterListInput(new File(ascflxFilename));
                nvpList = inputList.readData();
                flux.setData(nvpList);
            } catch (FileNotFoundException ex) {
                flux.setData(new INameValuePair[0][0]);
            }

            inputList.close();

            // ascpile
            String ascpileFilename = inputDirectory +
                    inputDirectory.substring(inputDirectory.lastIndexOf(File.separator)) +
                    "." + TitanConstants.PILE_FILE_EXT;

            String physicsModel = matModelMap.getPhysicsModel();
            //System.out.println("Saved material model: " + physicsModel);

            // If the physics model is TwoPhases_Coulomb, need to read in volume fraction
            if (TitanConstants.PhysicsModels[TitanConstants.PHYSICS_MODEL_TWOPHASES_COULOMB].compareTo(physicsModel) == 0) {
                piles.setHeadings(true);
            }

            try {
                inputList = new ParameterListInput(new File(ascpileFilename));
                nvpList = inputList.readData();
                piles.setData(nvpList);
            } catch (FileNotFoundException ex) {
                piles.setData(new INameValuePair[0][0]);
            }

            input.close();

            // ascplane
            String ascplaneFilename = inputDirectory +
                    inputDirectory.substring(inputDirectory.lastIndexOf(File.separator)) +
                    "." + TitanConstants.PLANE_FILE_EXT;
            try {
                inputList = new ParameterListInput(new File(ascplaneFilename));
                nvpList = inputList.readData();
                discharge.setData(nvpList);
            } catch (FileNotFoundException ex) {
                discharge.setData(new INameValuePair[0][0]);
            }

            input.close();

            JOptionPane.showMessageDialog(Titan.this,
                    "Successfully loaded the saved Titan run parameters",
                    "Titan Load",
                    JOptionPane.WARNING_MESSAGE);

        }
    }

    public class FetchData {

        public NameValuePairGroup getData() {

            NameValuePairGroup nvpg = new NameValuePairGroup();
            nvpg.addNameValuePair(fileIO.getData());
            nvpg.addNameValuePair(gis.getData());
            nvpg.addNameValuePair(runParameters.getData());
            nvpg.addNameValuePair(outputOptions.getData());
            nvpg.addNameValuePair(matModelMap.getData());
            return nvpg;
        }

        public NameValuePairGroup[] getPileData() {
            NameValuePairGroup nvpg[] = new NameValuePairGroup[piles.getData().length];
            for (int i = 0; i < piles.getData().length; i++) {
                nvpg[i] = new NameValuePairGroup();
                nvpg[i].addNameValuePair(piles.getData()[i]);
            }
            return nvpg;
        }

        public NameValuePairGroup[] getFluxData() {
            NameValuePairGroup nvpg[] = new NameValuePairGroup[flux.getData().length];
            for (int i = 0; i < flux.getData().length; i++) {
                nvpg[i] = new NameValuePairGroup();
                nvpg[i].addNameValuePair(flux.getData()[i]);
            }
            return nvpg;
        }

        public NameValuePairGroup[] getPlaneData() {
            NameValuePairGroup nvpg[] = new NameValuePairGroup[discharge.getData().length];
            for (int i = 0; i < discharge.getData().length; i++) {
                nvpg[i] = new NameValuePairGroup();
                nvpg[i].addNameValuePair(discharge.getData()[i]);
            }
            return nvpg;
        }
    }

    // This is the main class
    public static void main(String[] args) {

        Titan titan = new Titan();

        JFrame f = new JFrame();
        f.setSize(1200, 800);
        f.getContentPane().add(titan);
        //f.pack();
        f.setTitle(APP_NAME);
        f.setLocation(200, 200);
        f.setVisible(true);
        f.addWindowListener(titan.new TitanExitAdapter(titan));
    }

    private class TitanExitAdapter extends WindowAdapter {
        private Titan titan;

        public TitanExitAdapter(Titan t) {
            super();
            titan = t;
        }

        public void windowClosing(WindowEvent e) {
            Titan.Exit(titan, titan.getDB(), 0);
        }
    }

}
