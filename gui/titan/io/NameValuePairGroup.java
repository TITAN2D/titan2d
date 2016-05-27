package titan.io;

import java.util.ArrayList;

public class NameValuePairGroup {

	ArrayList<INameValuePair> myData;
	
	public NameValuePairGroup() {
		myData = new ArrayList<INameValuePair>();
	}
	
	public void addNameValuePair(INameValuePair nvp) {
		myData.add(nvp);
	}
	
	public void addNameValuePair(INameValuePair[] nvp) {
		for(int i = 0; i < nvp.length; i++) {
			myData.add(nvp[i]);
		}
	}
	
	public String getValue(String key) {
//System.out.println("  nvp len : " + myData.size());
		for(int i = 0; i < myData.size(); i++) {
//System.out.println("    key : " + ((INameValuePair)myData.get(i)).getName());
			if(((INameValuePair)myData.get(i)).getName().compareTo(key) == 0) {
				return ((INameValuePair)myData.get(i)).getValue();
			}
		}
		return null;
	}
	
	/**
	 * Changes the value of a name value pair held by the group.
	 * If the key does not exist in the group, no action is taken.
	 * @param key
	 * @param value
	 */
	public void setValue(String key, String value) {
		for(int i = 0; i < myData.size(); i++) {
		    if((myData.get(i)).getName().compareTo(key) == 0) {
		    	myData.get(i).setValue(value);
			    break;
		    }
	    }
	}
	
}
