package titan.gui;

import titan.graphics.CheckBoxGroup;
import titan.graphics.DirectorySelector;
import titan.graphics.RadioButtonGroup;
import titan.graphics.TextInput;
import titan.io.INameValuePair;
import titan.io.NameValuePairComponent;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

// This tab is display as the Time tab so rename
public class TabOutputOptions extends TitanTab {

    private RadioButtonGroup overwriteOutput;
    private RadioButtonGroup restartOutputEnabled;
    private TextInput resultDelta1;
    private TextInput saveDelta1;
    private RadioButtonGroup keepAll;
    private RadioButtonGroup keepRedundantData;
    private TextInput outputPrefix1;
    private RadioButtonGroup vizOutput;
    private TextInput resultDelta2;
    private TextInput saveDelta2;
    private TextInput outputPrefix2;

    public TabOutputOptions() {

        // General Parameters
        setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));

        values = new NameValuePairComponent[11];

        // Restart Output
        JPanel retartOutputPanel = new JPanel();
        retartOutputPanel.setLayout(new BoxLayout(retartOutputPanel, BoxLayout.Y_AXIS));
        retartOutputPanel.setBorder(new TitledBorder("Restart Output"));

        // Time Series Output
        JPanel timeSeriesOutputPanel = new JPanel();
        timeSeriesOutputPanel.setLayout(new BoxLayout(timeSeriesOutputPanel, BoxLayout.Y_AXIS));
        timeSeriesOutputPanel.setBorder(new TitledBorder("Time Series Output"));

        // Overwrite output
        overwriteOutput = new RadioButtonGroup("Overwrite Output", TitanConstants.TrueFalse);
        values[0] = new NameValuePairComponent(TitanConstants.OVERWRITE_OUTPUT, overwriteOutput);

        // Restart Output Enabled
        restartOutputEnabled = new RadioButtonGroup("Restart Output Enabled", TitanConstants.TrueFalse);
        values[1] = new NameValuePairComponent(TitanConstants.RESTART_OUTPUT_ENABLED, restartOutputEnabled);

        // Iterations between saves
        saveDelta1 = new TextInput("Iterations Between Restart File Output");
        values[2] = new NameValuePairComponent(TitanConstants.SAVE_TIME_DELTA1, saveDelta1);

        // Time between results output
        resultDelta1 = new TextInput("Time Between Restart File Output [s]");
        values[3] = new NameValuePairComponent(TitanConstants.RESULT_OUTPUT_TIME_DELTA1, resultDelta1);

        keepAll = new RadioButtonGroup("Keep All", TitanConstants.TrueFalse);
        values[4] = new NameValuePairComponent(TitanConstants.KEEP_ALL, keepAll);

        keepRedundantData = new RadioButtonGroup("Keep Redundant Data", TitanConstants.TrueFalse);
        values[5] = new NameValuePairComponent(TitanConstants.KEEP_REDUNDANT_DATA, keepRedundantData);

        // Output prefix 1
        outputPrefix1 = new TextInput("Restart Files Storage Directory");
        values[6] = new NameValuePairComponent(TitanConstants.OUTPUT_PREFIX1, outputPrefix1);
        outputPrefix1.setEditable(false);

        retartOutputPanel.add(restartOutputEnabled.getPanel());
        retartOutputPanel.add(saveDelta1.getPanel());
        retartOutputPanel.add(resultDelta1.getPanel());
        retartOutputPanel.add(keepAll.getPanel());
        retartOutputPanel.add(keepRedundantData.getPanel());
        retartOutputPanel.add(outputPrefix1.getPanel());

        vizOutput = new RadioButtonGroup("Visualization Output Type(s)", TitanConstants.VizTypes, RadioButtonGroup.MULTIPLE_SELECTION);
        vizOutput.setEditable(TitanConstants.VizTypes[TitanConstants.VIS_OUTPUT_WEBVIZ], false);
        vizOutput.setEditable(TitanConstants.VizTypes[TitanConstants.VIS_OUTPUT_GMFGVIZ], false);
        values[7] = new NameValuePairComponent(TitanConstants.VIS_OUTPUT, vizOutput);

        // Iterations between saves
        saveDelta2 = new TextInput("Iterations Between Snapshot Files Output");
        values[8] = new NameValuePairComponent(TitanConstants.SAVE_TIME_DELTA2, saveDelta2);

        // Time between results output
        resultDelta2 = new TextInput("Time Between Snapshot Files Output [s]");
        values[9] = new NameValuePairComponent(TitanConstants.RESULT_OUTPUT_TIME_DELTA2, resultDelta2);

        // Output prefix 2
        outputPrefix2 = new TextInput("Snapshot Files Storage Directory");
        values[10] = new NameValuePairComponent(TitanConstants.OUTPUT_PREFIX2, outputPrefix2);
        outputPrefix2.setEditable(false);

        timeSeriesOutputPanel.add(vizOutput.getPanel());
        timeSeriesOutputPanel.add(saveDelta2.getPanel());
        timeSeriesOutputPanel.add(resultDelta2.getPanel());
        timeSeriesOutputPanel.add(outputPrefix2.getPanel());

        // Add components to main panel
        add(overwriteOutput.getPanel());
        add(retartOutputPanel);
        add(timeSeriesOutputPanel);

        restartOutputEnabled.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ev) {
                saveDelta1.setEditable(false);
                resultDelta1.setEditable(false);
                keepAll.setEditable(false);
                keepRedundantData.setEditable(false);
                outputPrefix1.setEditable(false);
            }
        }, TitanConstants.FALSE);

        restartOutputEnabled.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent ev) {
                saveDelta1.setEditable(true);
                resultDelta1.setEditable(true);
                keepAll.setEditable(true);
                keepRedundantData.setEditable(true);
                outputPrefix1.setEditable(true);
            }
        }, TitanConstants.TRUE);
    }
    public void setData(INameValuePair[] data) {
        super.setData(data);

        if (restartOutputEnabled.getValue().compareTo(TitanConstants.FALSE) == 0) {
            saveDelta1.setEditable(false);
            resultDelta1.setEditable(false);
            keepAll.setEditable(false);
            keepRedundantData.setEditable(false);
            outputPrefix1.setEditable(false);
        } else {
            saveDelta1.setEditable(true);
            resultDelta1.setEditable(true);
            keepAll.setEditable(true);
            keepRedundantData.setEditable(true);
            outputPrefix1.setEditable(true);
        }
    }
}
