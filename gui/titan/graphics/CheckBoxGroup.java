package titan.graphics;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.ActionListener;
import java.util.ArrayList;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.border.EtchedBorder;

public class CheckBoxGroup implements IValueComponent {
	
	private JPanel panel;
	private String[] values;
	JCheckBox[] btns;
	
	public CheckBoxGroup(String lbl, String[] buttons) {
		panel = new JPanel(new GridLayout());
		panel.setBorder(BorderFactory.createEtchedBorder(EtchedBorder.LOWERED));

		values = buttons;
		btns = new JCheckBox[values.length];

		JLabel title = new JLabel(lbl);
		panel.add(title);
		
		JPanel selectionPanel = new JPanel();
		selectionPanel.setLayout(new BoxLayout(selectionPanel, BoxLayout.Y_AXIS));

		for(int i = 0; i < buttons.length; i++) {
			btns[i] = new JCheckBox(buttons[i]);
			selectionPanel.add(btns[i]);
		}

		panel.add(selectionPanel);
		
		// Do not allow height of object to be resized past
		// the preferred height
		Dimension myMaxSize = panel.getMaximumSize();
		myMaxSize.height = panel.getPreferredSize().height;
		panel.setMaximumSize(myMaxSize);
	}
	
	public Component getPanel() {
		return panel;
	}
	
	public void setValue(String v) {
		// List could be comma separated for multi selection
		ArrayList<String> valueList = new ArrayList<String>();
		String vals = v;
		int idx;
		while((idx = vals.indexOf(",")) >= 0) {
            valueList.add(vals.substring(0, idx));
			vals = vals.substring(idx+1, vals.length());
		}
        valueList.add(vals);

		for(int i = 0; i < btns.length; i++) {
			btns[i].setSelected(false);
		}
		
        for(int i = 0; i < valueList.size(); i++) {
			for(int j = 0; j < values.length; j++) {
				if(values[j].compareTo((String)(valueList.get(i))) == 0) {
					btns[j].setSelected(true);
					break;
				}
			}
        }

	}
	
	public String getValue() {

		String selected = new String();

	    for(int i = 0; i < btns.length; i++) {
	    	if(btns[i].isSelected()) {
	        	if(selected.length() == 0) {
	        	    selected = btns[i].getText();
	        	}
	        	else {
	        		selected = selected + "," + btns[i].getText();
	        	}
	    	}
	    }
        return selected;
	}
	
	public void addActionListener(ActionListener al, String btnValue) {
		for(int i = 0; i < values.length; i++) {
			if(values[i].compareTo(btnValue) == 0) {
				btns[i].addActionListener(al);
				break;
			}
		}
	}
	

}
