package titan.io;

import titan.graphics.IValueComponent;

public class NameValuePairComponent implements INameValuePair {

	private String name;
	private IValueComponent value;

	public NameValuePairComponent(String n, IValueComponent obj) {
		name = new String(n);
		value = obj;
	}
	
    public String getName() {
    	return name;
    }
    
    public String getValue() {
    	return value.getValue();
    }

    public void setValue(String v) {
    	value.setValue(v);
    }
}
