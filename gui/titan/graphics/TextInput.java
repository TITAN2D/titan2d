package titan.graphics;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridLayout;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

public class TextInput implements IValueComponent {

	JPanel panel;
	JLabel title;
	JTextField value;

	public TextInput(String lbl) {
		panel = new JPanel(new GridLayout());

		title = new JLabel(lbl);
		panel.add(title);

		value = new JTextField(20);
		panel.add(value);

		// Do not allow height of object to be resized past
		// the preferred height
		Dimension myMaxSize = panel.getMaximumSize();
		myMaxSize.height = panel.getPreferredSize().height;
		panel.setMaximumSize(myMaxSize);
	}
	
	public void setEditable(boolean b) {
		value.setEditable(b);
	}

	public Component getPanel() {
		return panel;
	}

	public String getLabelText() { return title.getText(); }

	public void setLabelText(String text) { title.setText(text); }

	public String getValue() {
		return value.getText();
	}

	public void setValue(String val) {
		value.setText(val);
	}

	public void setValueVisible(boolean b) {
		value.setVisible(b);
	}

	public void setToolTipText(String tip) {
	    value.setToolTipText(tip);
	}
}
