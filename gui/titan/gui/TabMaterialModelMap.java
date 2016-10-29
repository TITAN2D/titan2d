package titan.gui;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.*;
import javax.swing.border.TitledBorder;

import titan.graphics.RadioButtonGroup;
import titan.graphics.TextInput;
import titan.io.INameValuePair;
import titan.io.NameValuePairComponent;
import titan.io.NameValuePairGroup;

import java.io.*;

import java.util.regex.Pattern;
import java.util.regex.Matcher;
import java.util.Vector;

import titan.gui.Titan;
import titan.gui.Titan.FetchData;

public class TabMaterialModelMap extends TitanTab {

    private RadioButtonGroup physicsModel;
    private RadioButtonGroup useMaterialMap;
    private RadioButtonGroup stoppingCriteria;
    private TextInput numPhysicsModelParameters;
    private FetchData dataFetch;
    public Vector<String> materialNamesVector;
    public Vector<TextInput> materialModelVector;
    public int numMatMapMaterials;

    public JPanel materialModelParametersPanel;

    public TabMaterialModelMap(FetchData create) {

        dataFetch = create;

        materialModelVector = new Vector<TextInput>();

        // Number of values stored to the test .ascprm file
        values = new NameValuePairComponent[4 + TitanConstants.MAX_PHYSICS_MODEL_PARAMETERS];
        int ascprmIndex = 0;

        setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));

        // Physics model
        physicsModel = new RadioButtonGroup("Material Model", TitanConstants.PhysicsModels, RadioButtonGroup.SINGLE_SELECTION);
        values[ascprmIndex++] = new NameValuePairComponent(TitanConstants.PHYSICS_MODEL, physicsModel);
        add(physicsModel.getPanel());

        // Use GIS material map
        useMaterialMap = new RadioButtonGroup("Use Material Map", TitanConstants.TrueFalse);
        values[ascprmIndex++] = new NameValuePairComponent(TitanConstants.USE_GIS_MAT_MAP, useMaterialMap);
        add(useMaterialMap.getPanel());

        // Stopping criteria
        stoppingCriteria = new RadioButtonGroup("Stopping Criteria", TitanConstants.StoppingCriterias);
        values[ascprmIndex++] = new NameValuePairComponent(TitanConstants.STOPPING_CRITERIA, stoppingCriteria);

        materialModelParametersPanel = new JPanel();
        materialModelParametersPanel.setLayout(new BoxLayout(materialModelParametersPanel, BoxLayout.Y_AXIS));
        materialModelParametersPanel.setBorder(new TitledBorder("Material Model Parameters"));

        // Add component to main panel
        add(materialModelParametersPanel);

        // Number of material model parameters (depends on selected material model and useMaterialMap)
        numPhysicsModelParameters = new TextInput(" ");
        values[ascprmIndex++] = new NameValuePairComponent(TitanConstants.NUM_PHYSICS_MODEL_PARAMETERS, numPhysicsModelParameters);
        materialModelParametersPanel.add(numPhysicsModelParameters.getPanel());
        numPhysicsModelParameters.setEditable(false);
        numPhysicsModelParameters.setValueVisible(false);

        for (int i = 0; i < TitanConstants.MAX_PHYSICS_MODEL_PARAMETERS; i++) {
            materialModelVector.add(new TextInput(""));
            values[ascprmIndex++] = new NameValuePairComponent("MATERIAL_MODEL_TEXT_INPUT_" + i, materialModelVector.get(i));
            materialModelVector.get(i).setValueVisible(false);
            materialModelParametersPanel.add(materialModelVector.get(i).getPanel());
        }

        add(stoppingCriteria.getPanel());

        // Add the action listeners

        physicsModel.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ev) {
                useMaterialMap.setEditable(true);
                stoppingCriteria.setEditable(true);
                setUseMaterialMap(TitanConstants.PHYSICS_MODEL_COULOMB, true);
                titan.gui.Titan.updatePileTab(false);
            }
        }, TitanConstants.PhysicsModels[TitanConstants.PHYSICS_MODEL_COULOMB]);

        physicsModel.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ev) {
                useMaterialMap.setEditable(false);
                stoppingCriteria.setEditable(false);
                setUseMaterialMap(TitanConstants.PHYSICS_MODEL_TWOPHASES_PITMAN_LE, true);
                titan.gui.Titan.updatePileTab(true);
            }
        }, TitanConstants.PhysicsModels[TitanConstants.PHYSICS_MODEL_TWOPHASES_PITMAN_LE]);

        physicsModel.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ev) {
                useMaterialMap.setEditable(false);
                stoppingCriteria.setEditable(false);
                setUseMaterialMap(TitanConstants.PHYSICS_MODEL_VOELLMY_SALM, true);
                titan.gui.Titan.updatePileTab(false);
            }
        }, TitanConstants.PhysicsModels[TitanConstants.PHYSICS_MODEL_VOELLMY_SALM]);

        physicsModel.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ev) {
                useMaterialMap.setEditable(false);
                stoppingCriteria.setEditable(false);
                setUseMaterialMap(TitanConstants.PHYSICS_MODEL_POULIQUEN_FORTERRE, true);
                titan.gui.Titan.updatePileTab(false);
            }
        }, TitanConstants.PhysicsModels[TitanConstants.PHYSICS_MODEL_POULIQUEN_FORTERRE]);

        // useMaterialMap ActionListeners valid for Coulomb model only
        useMaterialMap.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ev) {
                for (int i = 0; i < TitanConstants.PhysicsModels.length; i++) {
                    if (physicsModel.getValue().compareTo(TitanConstants.PhysicsModels[i]) == 0) {
                        setUseMaterialMap(i, true);
                        break;
                    }
                }
            }
        }, TitanConstants.TRUE);

        useMaterialMap.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ev) {
                for (int i = 0; i < TitanConstants.PhysicsModels.length; i++) {
                    if (physicsModel.getValue().compareTo(TitanConstants.PhysicsModels[i]) == 0) {
                        setUseMaterialMap(i, true);
                        break;
                    }
                }
            }
        }, TitanConstants.FALSE);
    }

    public void setData(INameValuePair[] data) {
        super.setData(data);

        // Physics Model
        if (physicsModel.getValue().compareTo(TitanConstants.PhysicsModels[TitanConstants.PHYSICS_MODEL_COULOMB]) == 0) {
            stoppingCriteria.setEditable(true);
            useMaterialMap.setEditable(true);
        } else {
            stoppingCriteria.setEditable(false);
            useMaterialMap.setEditable(false);
        }

        for (int i = 0; i < TitanConstants.PhysicsModels.length; i++) {
            if (physicsModel.getValue().compareTo(TitanConstants.PhysicsModels[i]) == 0) {
                setUseMaterialMap(i, true);
                break;
            }
        }
    }

    public boolean setUseMaterialMap(int physicsModelIndex, boolean set) {

        // set == true to verify the GIS info and setup for Material Map processing;
        // set == false to verify the GIS info
        if (physicsModelIndex == TitanConstants.PHYSICS_MODEL_COULOMB &&
                useMaterialMap.getValue() == TitanConstants.TRUE) {

            NameValuePairGroup mainData = dataFetch.getData();

            String mainDirectory = mainData.getValue(TitanConstants.GIS_INFO_DIRECTORY);
            String subDirectory = mainData.getValue(TitanConstants.GIS_SUBDIR);
            String mapset = mainData.getValue(TitanConstants.GIS_MAPSET);
            String map = mainData.getValue(TitanConstants.GIS_MAP);
            if (mainDirectory.compareTo("") == 0 ||
                    subDirectory.compareTo("") == 0 ||
                    mapset.compareTo("") == 0 ||
                    map.compareTo("") == 0) {

                JOptionPane.showMessageDialog(TabMaterialModelMap.this,
                        "A field on the GIS tab is blank.\n" +
                                "Setting Use Material Map to False.",
                        "Material Map Error",
                        JOptionPane.ERROR_MESSAGE);
                useMaterialMap.setValue(TitanConstants.FALSE);
                setMaterialModel(physicsModelIndex);
                return false;
            }

            String matPath = mainDirectory + File.separator + subDirectory + File.separator + mapset + "/cats/" + map + "_Mat";
            System.out.println(matPath);

            // vector to hold material names
            materialNamesVector = new Vector<String>();

            // Try opening material map file (file name ends in '_Mat')
            try {

                FileInputStream fstream = new FileInputStream(matPath);

                DataInputStream in = new DataInputStream(fstream);
                BufferedReader br = new BufferedReader(new InputStreamReader(in));

                // Regular expression
                Pattern pattern = Pattern.compile(": ");

                String strLine;
                //Read File Line By Line
                while ((strLine = br.readLine()) != null) {
                    Matcher matcher = pattern.matcher(strLine);
                    if (matcher.find()) {
                        if (!materialNamesVector.contains(strLine.substring(matcher.end()))) {
                            // if its not already in the vector, then add it
                            materialNamesVector.add(strLine.substring(matcher.end()));
                        }
                    }
                }

                //Close the input stream
                in.close();

            } catch (Exception e) {
                JOptionPane.showMessageDialog(TabMaterialModelMap.this,
                        "Exception thrown:\n" +
                                e.getMessage() + ".\n" +
                                "Setting Use Material Map to False.",
                        "Material Map Error",
                        JOptionPane.ERROR_MESSAGE);
                useMaterialMap.setValue(TitanConstants.FALSE);
                setMaterialModel(physicsModelIndex);
                return false;
            }

            if (materialNamesVector.size() == 0) {

                JOptionPane.showMessageDialog(TabMaterialModelMap.this,
                        "No material names found in " + matPath + ".\n" +
                                "Setting Use Material Map to False.",
                        "Material Map Error",
                        JOptionPane.ERROR_MESSAGE);
                useMaterialMap.setValue(TitanConstants.FALSE);
                setMaterialModel(physicsModelIndex);
                return false;
            }

            if (materialNamesVector.size() > TitanConstants.MAX_MAT_MAP_MATERIALS) {

                JOptionPane.showMessageDialog(TabMaterialModelMap.this,
                        "The number of material names in " + matPath + " is greater than " +
                                TitanConstants.MAX_MAT_MAP_MATERIALS + ".\n" +
                                "Setting Use Material Map to False.",
                        "Material Map Error",
                        JOptionPane.ERROR_MESSAGE);
                useMaterialMap.setValue(TitanConstants.FALSE);
                setMaterialModel(physicsModelIndex);
                return false;
            }

            //debug
            //for ( String materialName : materialNamesVector ) {
            //    System.out.println(materialName);
            //}

            //System.out.println("materialNamesVector.size(): " + materialNamesVector.size());
            numMatMapMaterials = materialNamesVector.size();

            if (set == true) {

                String materialName;

                for (int i = 0; i < materialNamesVector.size(); i++) {
                    materialName = materialNamesVector.get(i);
                    // Read by this and the TitanRunInput class
                    TitanConstants.MATERIAL_MAP_NAMES[i] = materialName;
                    //System.out.println("TitanConstants.MATERIAL_MODEL_NAMES[" + i + "]: " + TitanConstants.MATERIAL_MAP_NAMES[i]);
                }
            }
        }

        if (set == true)
            setMaterialModel(physicsModelIndex);

        return true;
    }

    public void setMaterialModel(int physicsModelIndex) {

        String[] matParmNames = TitanConstants.matParmNames[physicsModelIndex];
        String matParmName;
        String[] matParmUnits = TitanConstants.matParmUnits[physicsModelIndex];
        String matParmUnit;

        int index = 0;
        if (physicsModelIndex == TitanConstants.PHYSICS_MODEL_COULOMB &&
                useMaterialMap.getValue() == TitanConstants.TRUE) {

            int numParms = matParmNames.length - 1 + numMatMapMaterials;
            //System.out.println("numParms: " + numParms);

            numPhysicsModelParameters.setValue("" + numParms);

            // int_frict
            matParmName = matParmNames[TitanConstants.INT_FRICT_INDEX];
            matParmUnit = matParmUnits[TitanConstants.INT_FRICT_INDEX];
            materialModelVector.get(TitanConstants.INT_FRICT_INDEX).setLabelText(matParmName + " " + matParmUnit);
            materialModelVector.get(TitanConstants.INT_FRICT_INDEX).setValueVisible(true);
            index++;

            //bed_frict
            matParmName = matParmNames[TitanConstants.BED_FRICT_INDEX];
            matParmUnit = matParmUnits[TitanConstants.BED_FRICT_INDEX];
            for (int i = 1; i < numParms; i++) {
                materialModelVector.get(i).setLabelText(TitanConstants.MATERIAL_MAP_NAMES[i - 1] + " " + matParmName + " " + matParmUnit);
                materialModelVector.get(i).setValueVisible(true);
                index++;
            }
        } else {

            numPhysicsModelParameters.setValue("" + matParmNames.length);

            for (int i = 0; i < matParmNames.length; i++) {
                matParmName = matParmNames[i];
                matParmUnit = matParmUnits[i];
                materialModelVector.get(i).setLabelText(matParmName + " " + matParmUnit);
                materialModelVector.get(i).setValueVisible(true);
                index++;
            }
        }

        for (int i = index; i < TitanConstants.MAX_PHYSICS_MODEL_PARAMETERS; i++) {
            materialModelVector.get(i).setLabelText(" ");
            materialModelVector.get(i).setValueVisible(false);
        }
    }

    public String getPhysicsModel() {

        // Return the currently selected physics model
        String retVal = physicsModel.getValue();

        return retVal;
    }
}

