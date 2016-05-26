package titan.graphics;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridLayout;

import javax.swing.JLabel;
import javax.swing.JPanel;

public class LabelDisplay {

	JPanel panel;

    public LabelDisplay(String lbl, String val) {
		
		panel = new JPanel(new GridLayout());

		
//panel.setBackground(Color.RED);

		JLabel title = new JLabel(lbl);
//title.setOpaque(true);
//title.setBackground(Color.GREEN);
		panel.add(title);

		JLabel value = new JLabel(val);
//value.setOpaque(true);
//value.setBackground(Color.YELLOW);
		panel.add(value);

		// Do not allow height of object to be resized past
		// the preferred height
		Dimension myMaxSize = panel.getMaximumSize();
		myMaxSize.height = panel.getPreferredSize().height;
		panel.setMaximumSize(myMaxSize);
	}
	
	public Component getPanel() {
		return panel;
	}


}
