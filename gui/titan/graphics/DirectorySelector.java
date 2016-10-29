package titan.graphics;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.GroupLayout;
import javax.swing.JButton;
import javax.swing.filechooser.FileFilter;
import javax.swing.JFileChooser;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.GroupLayout.SequentialGroup;
import javax.swing.filechooser.FileNameExtensionFilter;

public class DirectorySelector implements IValueComponent {

	public static int FULL_PATH = 0;
	public static int DIR_ONLY = 1;

	public static int DIRECTORIES_ONLY = 0;
	public static int FILES_AND_DIRECTORIES = 1;

	private JPanel panel;
//	private String directory = new String("");
    private JTextField directory;
    private JButton select;
    private int outputStyle = FULL_PATH;
    private int showStyle = FILES_AND_DIRECTORIES;
	private String attribute = "";
	// For the Input File directory selector on the File IO tab and
	// the Restart File directory selector on the Run Parameters tab,
	// the user specified initial input directory.
	// Also, for the Restart File directory selector,
	// currently loaded input directory from the file io tab
	private String inputDirectory = "";

	public DirectorySelector(String attr) {

		attribute = attr;
		panel = new JPanel();
		GroupLayout layout = new GroupLayout(panel);
//		GridLayout layout = new GridLayout();
		panel.setLayout(layout);
		layout.setAutoCreateGaps(true);
		layout.setAutoCreateContainerGaps(true);
		
		select = new JButton(attr);
		
		ActionListener al = new ChooseDirectory();
		select.addActionListener(al);
		panel.add(select);

		directory = new JTextField(50);
		panel.add(directory);

		SequentialGroup horizGroup = layout.createSequentialGroup()
		.addComponent(select)
		.addComponent(directory);
		
		layout.setHorizontalGroup(horizGroup);
		
		layout.setVerticalGroup(
				layout.createSequentialGroup()
				.addGroup(layout.createParallelGroup(GroupLayout.Alignment.BASELINE)	
				.addComponent(select)
				.addComponent(directory)));
	}

	public Component getPanel() {
		return panel;
	}

	public String getValue() {
		return directory.getText();
	}

	public void setValue(String val) {
		directory.setText(val);
	}

	public void setOutputStyle(int style) {
		outputStyle = style;
	}

	public int getButtonWidth() {
		return select.getMaximumSize().width;
	}
	
	// Sets the preferred and maximum size of the button, which
	// should be enough to make the button appear that size
	// with the layout manager used.  Do not set the minimum size
	// in case the main app is resized and made smaller
	public void setButtonWidth(int size) {
		select.setMaximumSize(new Dimension(size, select.getMaximumSize().height));
		select.setPreferredSize(new Dimension(size, select.getMaximumSize().height));
	}
	
	/**
	 * Set the chooser to files and directories or directories only.
	 * @param style Should be : DIRECTORIES_ONLY or  FILES_AND_DIRECTORIES
	 */
	public void setChooserStyle(int style) {
		showStyle = style;
	}

	public void setInputDirectory (String inputDir) {
		// Currently loaded input directory from the file io tab
	    inputDirectory = inputDir;
	}

	public void setEditable(boolean state) {
		directory.setEditable(state);
		select.setEnabled(state);
	}
	
	private class ChooseDirectory implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			JFileChooser ds = new JFileChooser();

			if(showStyle == DIRECTORIES_ONLY) {
			    ds.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
			}
			else if(showStyle == JFileChooser.FILES_AND_DIRECTORIES) {
				ds.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
			}

			if(directory.getText().compareTo("") != 0){
				ds.setCurrentDirectory(new File(directory.getText()));
			} else if (inputDirectory.compareTo("") != 0) {
				ds.setCurrentDirectory(new File(inputDirectory));
			}

			if (attribute.compareTo("Restart File") == 0) {
				FileFilter filter = new FileNameExtensionFilter("Hierarchical Data Format (HDF5)","h5");
				ds.setFileFilter(filter);
			}

			int returnVal = ds.showOpenDialog(panel);
			if(returnVal == JFileChooser.APPROVE_OPTION) {
				File file = ds.getSelectedFile();
				if(outputStyle == FULL_PATH) {
				    directory.setText(file.getPath());
				}
				else if(outputStyle == DIR_ONLY) {
					directory.setText(file.getName());
				}
			}
		}
	}
}
