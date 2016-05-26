package titan.io;

import titan.graphics.IValueComponent;

public class NameValuePairSimple implements INameValuePair {

	private String name;
	private String value;

	public NameValuePairSimple(String n, String c) {
		name = new String(n);
		value = new String(c);
	}
	
	public String getName() {
		return name;
	}

	public String getValue() {
		return value;
	}

	public void setValue(String v) {
		value = new String(v);
	}

}
