package titan.graphics;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridLayout;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

public class TextArea implements IValueComponent {

	JPanel panel;
	JTextArea value;
	
	public TextArea(String lbl) {
		panel = new JPanel(new GridLayout());

		JLabel title = new JLabel(lbl);
		panel.add(title);

		value = new JTextArea();
		value.setRows(4);
		
		JScrollPane textPain = new JScrollPane(value);
		textPain.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
		
		panel.add(textPain);

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

	public String getValue() {
		return value.getText();
	}

	public void setValue(String val) {
		value.setText(val);
	}

}
