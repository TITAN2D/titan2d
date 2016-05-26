package titan.gui;

import java.awt.Component;
import java.awt.Dimension;

import javax.swing.Box;
import javax.swing.JPanel;

import titan.io.INameValuePair;


public abstract class TitanTab extends JPanel {

	protected INameValuePair[] values;
	
    public INameValuePair[] getData() {
    	return values;
    }

    public void setData(INameValuePair[] data) {
		try {
    	    for(int i = 0; i < data.length; i++) {
    		    for(int j = 0; j < values.length; j++) {
    			    if(data[i].getName().compareTo(values[j].getName()) == 0) {
    				    values[j].setValue(data[i].getValue());
    				    break;
    			}
    		}
    	}
		} catch (Exception ee) {
			// This will occur when the number of allowable fields for a parameter has decreased.
		}
    }

    public String getValue(String key) {
		String value = null;
		for(int i = 0; i < values.length; i++) {
		    if(values[i].getName() == key) {
		    	value = new String(values[i].getValue());
		    	break;
		    }
		}
		return value;
    }
    
    public Component add(Component comp) {
    	Component c = super.add(comp);
		super.add(Box.createRigidArea(new Dimension(0,5)));
    	return c;
    }
}
