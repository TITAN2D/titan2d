package titan.gui;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.*;
import javax.swing.border.TitledBorder;

import titan.graphics.DirectorySelector;
import titan.graphics.RadioButtonGroup;
import titan.graphics.TextInput;
import titan.io.INameValuePair;
import titan.io.NameValuePairComponent;

public class TabRunParameters extends TitanTab {

	private RadioButtonGroup restartEnabled;
	private DirectorySelector restartFile;
	private TextInput maxTimeSteps;
	private TextInput maxTime;
	private TextInput cellsAcross;
	private RadioButtonGroup AMR;
	private RadioButtonGroup orderMethod;
	private RadioButtonGroup interfaceCapturingType;
	private RadioButtonGroup scaleSim;
	private TextInput lengthScale;
	private TextInput gravityScale;
	private TextInput heightScale;
	private TextInput runID;
	private TextInput flowOutlineHgt;
	private TextInput testFlowHgt;
	private TextInput testX;
	private TextInput testY;
	private TextInput statPropsPrefix;
	private RadioButtonGroup outlinePropsEnabled;
	private TextInput maxLinearSize;
	private RadioButtonGroup initSize;
	private TextInput outlinePropsPrefix;

	public TabRunParameters(String inputDir) {

		// General Parameters
		setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));

		values = new NameValuePairComponent[22];

		// Restart Enabled
		restartEnabled = new RadioButtonGroup("Restart", TitanConstants.TrueFalse);
		values[0] = new NameValuePairComponent(TitanConstants.RESTART_ENABLED, restartEnabled);

		// Restart file
		restartFile = new DirectorySelector("Restart File");
		restartFile.setChooserStyle(DirectorySelector.FILES_AND_DIRECTORIES);
		restartFile.setInputDirectory(inputDir);
		values[1] = new NameValuePairComponent(TitanConstants.RESTART_FILE, restartFile);

		// max number of time steps
		maxTimeSteps = new TextInput("Maximum Number of Time Steps");
		values[2] = new NameValuePairComponent(TitanConstants.MAX_NUM_TIME_STEPS, maxTimeSteps);

		// Maximum time
		maxTime = new TextInput("Maximum Time [s]");
		values[3] = new NameValuePairComponent(TitanConstants.MAX_TIME, maxTime);

		// Number of cells across smallest pile/flux-source diameter
		cellsAcross = new TextInput("Number of Computational Cells Across Smallest Pile/Flux-Source Diameter");
		values[4] = new NameValuePairComponent(TitanConstants.NUM_CELLS_ACROSS, cellsAcross);

		// AMR - Adaptive Mesh Refinement
		AMR = new RadioButtonGroup("AMR", TitanConstants.TrueFalse);
		values[5] = new NameValuePairComponent(TitanConstants.AMR, AMR);

		// Order method, First or Second
		orderMethod = new RadioButtonGroup("Order Method", TitanConstants.OrderMethod);
		values[6] = new NameValuePairComponent(TitanConstants.ORDER_METHOD, orderMethod);

		// Interface capturing type, Heuristic or LevelSet
		interfaceCapturingType = new RadioButtonGroup("Interface Capturing Type", TitanConstants.InterfaceCapturingType);
		values[7] = new NameValuePairComponent(TitanConstants.INTERFACE_CAPTURING_TYPE, interfaceCapturingType);

		// Scale Simulation and Scale Length
		JPanel scalePanel = new JPanel();
		scalePanel.setLayout(new BoxLayout(scalePanel, BoxLayout.Y_AXIS));
		scalePanel.setBorder(new TitledBorder("Scale Parameters"));

		scaleSim = new RadioButtonGroup("Scale Simulation", TitanConstants.TrueFalse);
		lengthScale = new TextInput("Length Scale [m]");
		gravityScale = new TextInput("Gravity Scale [ms^-2]");
		heightScale = new TextInput("Height Scale [m]");
		scalePanel.add(scaleSim.getPanel());
		scalePanel.add(lengthScale.getPanel());
		scalePanel.add(gravityScale.getPanel());
		scalePanel.add(heightScale.getPanel());
		values[8] = new NameValuePairComponent(TitanConstants.SCALE_SIM, scaleSim);
		values[9] = new NameValuePairComponent(TitanConstants.LENGTH_SCALE, lengthScale);
		values[10] = new NameValuePairComponent(TitanConstants.GRAVITY_SCALE, gravityScale);
		values[11] = new NameValuePairComponent(TitanConstants.HEIGHT_SCALE, heightScale);

		// Stat Props
		JPanel statPropsPanel = new JPanel();
		statPropsPanel.setLayout(new BoxLayout(statPropsPanel, BoxLayout.Y_AXIS));
		statPropsPanel.setBorder(new TitledBorder("Stat Props"));

		runID = new TextInput("Run ID");
		values[12] = new NameValuePairComponent(TitanConstants.RUN_ID, runID);

		flowOutlineHgt = new TextInput("Height used to define flow outline (>0) [m]");
		values[13] = new NameValuePairComponent(TitanConstants.FLOW_OUTLINE_HGT, flowOutlineHgt);

		// test if flow reaches height
		testFlowHgt = new TextInput("Test if flow Reaches Height [m]");
		values[14] = new NameValuePairComponent(TitanConstants.TEST_FLOW_HEIGHT_MIN, testFlowHgt);
		
		// at test point x
		testX = new TextInput("Flow Height Test Point X Location [UTM E]");
		values[15] = new NameValuePairComponent(TitanConstants.TEST_FLOW_X_LOC, testX);

		// at test point y
		testY = new TextInput("Flow Height Test Point Y Location [UTM N]");
		values[16] = new NameValuePairComponent(TitanConstants.TEST_FLOW_Y_LOC, testY);

		statPropsPrefix = new TextInput("Stat Props Prefix");
		values[17] = new NameValuePairComponent(TitanConstants.STAT_PROPS_PREFIX, statPropsPrefix);

		statPropsPanel.add(runID.getPanel());
		statPropsPanel.add(flowOutlineHgt.getPanel());
		statPropsPanel.add(testFlowHgt.getPanel());
		statPropsPanel.add(testX.getPanel());
		statPropsPanel.add(testY.getPanel());
		statPropsPanel.add(statPropsPrefix.getPanel());

		// Outline Props
		JPanel outlinePropsPanel = new JPanel();
		outlinePropsPanel.setLayout(new BoxLayout(outlinePropsPanel, BoxLayout.Y_AXIS));
		outlinePropsPanel.setBorder(new TitledBorder("Outline Props"));

	    // Enabled?
		outlinePropsEnabled = new RadioButtonGroup("Outline Props Enabled", TitanConstants.TrueFalse);
		values[18] = new NameValuePairComponent(TitanConstants.OUTLINE_PROPS_ENABLED, outlinePropsEnabled);

		// Max linear size
		maxLinearSize = new TextInput("Max Linear Size");
		values[19] = new NameValuePairComponent(TitanConstants.MAX_LINEAR_SIZE, maxLinearSize);

		// Init Size
		initSize = new RadioButtonGroup("Init Size", TitanConstants.InitSizes, RadioButtonGroup.SINGLE_SELECTION);
		values[20] = new NameValuePairComponent(TitanConstants.INIT_SIZE, initSize);

		outlinePropsPrefix = new TextInput("Outline Props Prefix");
		values[21] = new NameValuePairComponent(TitanConstants.OUTLINE_PROPS_PREFIX, outlinePropsPrefix);

		outlinePropsPanel.add(outlinePropsEnabled.getPanel());
		outlinePropsPanel.add(maxLinearSize.getPanel());
		outlinePropsPanel.add(initSize.getPanel());
		outlinePropsPanel.add(outlinePropsPrefix.getPanel());

		// Add components to main panel
		add(restartEnabled.getPanel());
		add(restartFile.getPanel());
		add(maxTimeSteps.getPanel());
		add(maxTime.getPanel());
		add(cellsAcross.getPanel());
		add(AMR.getPanel());
		add(orderMethod.getPanel());
		add(interfaceCapturingType.getPanel());
		add(scalePanel);
		add(statPropsPanel);
		add(outlinePropsPanel);

		restartEnabled.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent ev) {
				restartFile.setEditable(false);
				maxTimeSteps.setLabelText("Maximum Number of Time Steps");
				maxTime.setLabelText("Maximum Time [s]");
			}
		}, TitanConstants.FALSE);

		restartEnabled.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent ev) {
				restartFile.setEditable(true);
				maxTimeSteps.setLabelText("Additional Maximum Number of Time Steps");
				maxTime.setLabelText("Additional Maximum Time [s]");
			}
		}, TitanConstants.TRUE);

		scaleSim.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent ev) {
				lengthScale.setEditable(true);
				gravityScale.setEditable(true);
				heightScale.setEditable(true);
			}
		}, TitanConstants.TRUE);
		scaleSim.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent ev) {
				lengthScale.setEditable(false);
				gravityScale.setEditable(false);
				heightScale.setEditable(false);
			}
		}, TitanConstants.FALSE);
	}

	// Coordinate with the file io tab input file directory selector
	public void setInputDirectory (String inputDir) {
		restartFile.setInputDirectory(inputDir);
	}

	public void setData(INameValuePair[] data) {
        super.setData(data);

		if(scaleSim.getValue().compareTo(TitanConstants.FALSE) == 0) {
			lengthScale.setEditable(false);
			gravityScale.setEditable(false);
			heightScale.setEditable(false);
		} else {
			lengthScale.setEditable(true);
			gravityScale.setEditable(true);
			heightScale.setEditable(true);
		}

		if(restartEnabled.getValue().compareTo(TitanConstants.FALSE) == 0) {
			restartFile.setEditable(false);
			maxTimeSteps.setLabelText("Maximum Number of Time Steps");
			maxTime.setLabelText("Maximum Time [s]");
		} else {
			restartFile.setEditable(true);
		    maxTimeSteps.setLabelText("Additional Maximum Number of Time Steps");
		    maxTime.setLabelText("Additional Maximum Time [s]");
		}
	}
}
