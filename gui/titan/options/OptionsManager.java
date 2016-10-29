package titan.options;

import titan.gui.TitanConstants;
import titan.io.NameValuePairGroup;
import titan.io.NameValuePairSimple;

public class OptionsManager {

	private NameValuePairGroup myValues;
	public static final String TITAN_BINARY_MODE = "TITAN_BINARY_MODE";
	public static final String TITAN_BINARY_DIR  = "TITAN_BINARY_DIR";
	
	public OptionsManager() {
		myValues = new NameValuePairGroup();
		myValues.addNameValuePair(new NameValuePairSimple[] {
				new NameValuePairSimple(TITAN_BINARY_MODE, TitanConstants.FALSE),
				new NameValuePairSimple(TITAN_BINARY_DIR, "")});
	}
	
	public void setBinaryMode(boolean mode) {
		if(mode == true)
		    myValues.setValue(TITAN_BINARY_MODE, TitanConstants.TRUE);
		else
			myValues.setValue(TITAN_BINARY_MODE, TitanConstants.FALSE);
	}
	
	public boolean getBinaryMode() {
		return Boolean.parseBoolean(myValues.getValue(TITAN_BINARY_MODE));
	}

	public void setBinaryDirectory(String dir) {
		myValues.setValue(TITAN_BINARY_DIR, dir);
	}
	
	public String getBinaryDirectory() {
		return myValues.getValue(TITAN_BINARY_DIR);
	}
}
