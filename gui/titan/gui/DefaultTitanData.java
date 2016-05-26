package titan.gui;

import titan.io.INameValuePair;
import titan.io.NameValuePairSimple;

public class DefaultTitanData {

	INameValuePair[] defaultData;
	
	public DefaultTitanData() {
		defaultData = new INameValuePair[TitanConstants.parameterKeyList.length];
		for(int i = 0; i < TitanConstants.parameterKeyList.length; i++) {
			defaultData[i] = new NameValuePairSimple(TitanConstants.parameterKeyList[i],
					                                 TitanConstants.parameterValueDefaults[i]);
		}

	}
	
	public INameValuePair[] getDefaults() {
		return defaultData;
	}
	
}
