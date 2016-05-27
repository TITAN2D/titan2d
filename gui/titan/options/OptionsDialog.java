package titan.options;

import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JPanel;
import javax.swing.border.TitledBorder;

import titan.graphics.DirectorySelector;
import titan.graphics.RadioButtonGroup;
import titan.gui.TitanConstants;

public class OptionsDialog extends JDialog {

	private static final String binarySelectionString = new String("Use External Titan Binaries");
	private OptionsManager optionsManager;
	private RadioButtonGroup rb;
	private DirectorySelector ds;
	
	public OptionsDialog(Component owner, OptionsManager manager) {
		
		optionsManager = manager;
		
		getContentPane().setLayout(new BoxLayout(getContentPane(), BoxLayout.Y_AXIS));
		setTitle("Options");
		setResizable(false);
		this.setLocation(100, 100);
		setModal(true);
		this.setDefaultCloseOperation(JDialog.HIDE_ON_CLOSE);

		JPanel optionsPanel = new JPanel();
		optionsPanel.setLayout(new BoxLayout(optionsPanel, BoxLayout.Y_AXIS));
		
		JPanel binaryPanel = new JPanel();
		binaryPanel.setLayout(new BoxLayout(binaryPanel, BoxLayout.Y_AXIS));
		binaryPanel.setBorder(new TitledBorder("Titan Application Version"));
		
		rb = new RadioButtonGroup(binarySelectionString, TitanConstants.TrueFalse);
		
		String binMode = null;
		if(optionsManager.getBinaryMode() == true)
			binMode = new String(TitanConstants.TRUE);
		else
			binMode = new String(TitanConstants.FALSE);
		rb.setValue(binMode);

		binaryPanel.add(rb.getPanel());
		
		ds = new DirectorySelector("Titan Binary Directory");
		ds.setValue(optionsManager.getBinaryDirectory());
		ds.setChooserStyle(DirectorySelector.DIRECTORIES_ONLY);
		if(optionsManager.getBinaryMode() == false) {
			ds.setEditable(false);
		}
		binaryPanel.add(ds.getPanel());
		
		rb.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				ds.setEditable(true);
			}
		}, TitanConstants.TRUE);
		
		rb.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				ds.setEditable(false);
			}
		}, TitanConstants.FALSE);
		
		optionsPanel.add(binaryPanel);
		getContentPane().add(optionsPanel);
		
		JButton closeButton = new JButton("Close");
		closeButton.addActionListener(new CloseButtonListener());
		
		getContentPane().add(closeButton);

        pack();
        
        this.addWindowListener(new WindowAdapter() {
        	public void windowClosing(WindowEvent e) {
        		OptionsDialog.this.dispose();
        	}
	    });
	}
	
    private class CloseButtonListener implements ActionListener {
    	public void actionPerformed(ActionEvent ev) {
    		optionsManager.setBinaryMode(Boolean.parseBoolean(rb.getValue()));
    		if(rb.getValue().compareTo(TitanConstants.TRUE) == 0) {
    			optionsManager.setBinaryDirectory(ds.getValue());
    		}
    		else {
    		    optionsManager.setBinaryDirectory("");
    		}
    		dispose();
    	}
    }
	
}
