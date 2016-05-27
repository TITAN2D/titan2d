package titan.gui;

import java.awt.GridBagConstraints;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.Checkbox;
import java.io.File;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JPanel;
import javax.swing.border.TitledBorder;
import javax.swing.JTextArea;
import javax.swing.JScrollPane;
import javax.swing.JOptionPane;

import titan.graphics.CheckBoxGroup;
import titan.graphics.DirectorySelector;
import titan.graphics.TextInput;
import titan.io.INameValuePair;
import titan.io.NameValuePairComponent;
import titan.io.NameValuePairGroup;

public class TabFileIO extends TitanTab {

	private Titan.SaveTitanRun save;
	private Titan.LoadTitanRun load;
	private DirectorySelector inputFile;
	private DirectorySelector saveDir;
	private TextInput runName;
	private TextInput email;

	public TabFileIO(String inputDir, Titan.SaveTitanRun sv, Titan.LoadTitanRun ld) {

		setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));

		save = sv;
		load = ld;

		values = new INameValuePair[5];

		// General Parameters		
		JPanel loadPanel = new JPanel();
		loadPanel.setLayout(new BoxLayout(loadPanel, BoxLayout.Y_AXIS));
		loadPanel.setBorder(new TitledBorder("Load Titan Run"));
		
		// Load input file
		inputFile = new DirectorySelector("Input Directory");
		inputFile.setChooserStyle(DirectorySelector.DIRECTORIES_ONLY);
		inputFile.setInputDirectory(inputDir);
		values[0] = new NameValuePairComponent(TitanConstants.LOAD_DIRECTORY, inputFile);

		JButton loadBtn = new JButton("Load");
		loadBtn.setToolTipText("Load previously saved Titan run variables");
		loadBtn.addActionListener(new LoadAction());
		
		GridBagConstraints gbc1 = new GridBagConstraints();
		gbc1.gridx = 0;
		gbc1.gridy = 0;
		loadPanel.add(inputFile.getPanel(), gbc1);
		GridBagConstraints gbc2 = new GridBagConstraints();
		gbc2.gridx = 0;
		gbc2.gridy = 1;
		loadPanel.add(loadBtn, gbc2);

		// Save input
		JPanel savePanel = new JPanel();
		savePanel.setLayout(new BoxLayout(savePanel, BoxLayout.Y_AXIS));
		savePanel.setBorder(new TitledBorder("Save Titan Run"));
		
		saveDir = new DirectorySelector("Base Save Directory (REQUIRED)");
		saveDir.setChooserStyle(DirectorySelector.DIRECTORIES_ONLY);
		values[1] = new NameValuePairComponent(TitanConstants.SAVE_DIRECTORY, saveDir);
		
		runName = new TextInput("Run Name");
		values[2] = new NameValuePairComponent(TitanConstants.RUN_NAME, runName);

		JButton saveBtn = new JButton("Save");
		saveBtn.setToolTipText("Save current run variables for later use");
		saveBtn.addActionListener(new SaveAction());

		savePanel.add(saveDir.getPanel());
		savePanel.add(runName.getPanel());
		savePanel.add(saveBtn);
		
		//iRods checkbox
		//JPanel irodsPanel = new JPanel();
		//Checkbox irodsCheckbox = new Checkbox("Save output to iRODS?");
		String[] irodsArray = new String[1];
		irodsArray[0] = "USE_IRODS";
		CheckBoxGroup irodsCheckbox = new CheckBoxGroup("iRODS output:", irodsArray);
		values[3] = new NameValuePairComponent("IRODS_OUTPUT", irodsCheckbox);
		////values[3] = new NameValuePairComponent(TitanConstants.RUN_NAME, runName);
		////values[3] = new NameValuePairComponent("IRODS_SAVE", irodsCheckbox);
		//irodsPanel.add(irodsCheckbox.getPanel());
		//savePanel.add(irodsPanel);
		savePanel.add(irodsCheckbox.getPanel());

		// Make Load and Save directory selector buttons same size
		// for esthetic purposes
		int bSize = inputFile.getButtonWidth() > saveDir.getButtonWidth() ?
				inputFile.getButtonWidth() : saveDir.getButtonWidth();
		inputFile.setButtonWidth(bSize);
		saveDir.setButtonWidth(bSize);
		
		
		// Help box
		JPanel helpPanel = new JPanel();
		helpPanel.setLayout(new BoxLayout(helpPanel, BoxLayout.Y_AXIS));
		helpPanel.setBorder(new TitledBorder("Load/Save Instructions"));
		
		JTextArea textArea = new JTextArea();
		textArea.setText("The Load/Save tab helps you to load a set of previously " + 
		"saved Titan parameters, or save the current set, which can save you a lot " + 
		"of time/effort for frequently used sets of run parameters." + 
		"\n\nSAVE\nIn order to save " + 
		"all of the input variables to a new Titan run, you need to specify a save " + 
		"directory called a 'Base Save Directory' in the above 'Save Titan Run' box.  " +
		"To specify this directory, click on the 'Base Save Directory' button and choose " +
		"a directory where you would like to save your Titan input variables for this " +
		"run (NOTE: This directory is also where your Titan run results are placed).  After " + 
		"selecting a save directory, a name for the current run must be " +
		"selected.  This can be any name you want and it is actually a folder inside " +
		"the previous sub-directory where all the files with the parameters and results " +
		"are saved.  You need to press the 'Save' button in order to create the folder " +
		"(if it does not exist already) and to save the parameters that you have entered." +
		"\n\nLOAD\nTo load a previous saved Titan run, click on the 'Input Directory' " + 
		"button in the 'Load Titan Run' box and select the folder that contains the run " +
		"parameter information.  Once the 'Input Directory' text field has been populated, " +
		"click the 'Load' button to load the variables from the previous run directory that " +
		"you selected.\n\nNOTE: The Load steps are not necessary to run Titan.  " +
		"It is only for your convenience so that you can save time by " +
		"loading previous Titan runs.  The only text field that needs to be populated " +
		"in order for Titan to run, is the 'Base Save Directory' field." +
		"\n\nTITAN OUTPUT DIRECTORY\n" +
		"After running the Titan simulation, the output files are put into a directory that is " +
		"named by the time stamp of when it was run.  You can find this output directory in the " +
		"directory you specified in the 'Base Save Directory' field.  If you also specified " +
		"a Run Name, then your Titan output will be in base_save_dir/run_name/time_stamp." +
		"\n\nIRODS OUTPUT\n" +
		"By checking the 'USE_IRODS' box, your output from a submit Titan run will be placed " +
		"into your iRODS remote storage. Your output files will be placed in your home directory " +
		"in iRODS in the folder TitanVhubResults." + 
		"\n\nFor more information visit the following web pages:\n" +
		"http://vhub.org/tools/titan2d\nhttps://vhub.org/resources/3225/download/Titan_tutorial_VHub.pdf");
		
		textArea.setLineWrap(true);
		textArea.setWrapStyleWord(true);
		textArea.setEditable(false);
		
		JScrollPane scrollingArea = new JScrollPane(textArea);
		helpPanel.add(scrollingArea);

		// email address
		email = new TextInput("Email Address");
		values[4] = new NameValuePairComponent(TitanConstants.EMAIL_ADDRESS, email);

		// add JPanel's
		add(loadPanel);
		add(savePanel);
		//add(irodsPanel);
		add(helpPanel);
		add(email.getPanel());
	}
	
	private class LoadAction implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			load.Load(inputFile.getValue());
		}
	}
	
	private class SaveAction implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			boolean b = save.Save(saveDir.getValue(), runName.getValue());
			if (b) {
			    JOptionPane.showMessageDialog(TabFileIO.this,
						                   "Successfully Saved the current Titan run parameters",
						                   "Titan Save",
						                   JOptionPane.WARNING_MESSAGE);
			}
		}
	}

	// Coordinate with the run parameters tab restart file directory selector
	public String getInputDirectory () {
		return inputFile.getValue();
	}
}
