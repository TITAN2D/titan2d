package titan.io;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;

public class ParameterInput extends InputReader {
	
	public ParameterInput(File inputFile) throws FileNotFoundException {
        super(inputFile);
	}
	
	public INameValuePair[] readData() {
		String input;
		ArrayList nvpList = new ArrayList();
		try {
			while((input = stream.readLine()) != null) {
				removeSurroundingWhitespace(input);
				int equalIndex = input.indexOf("=");
				if(equalIndex != -1) {
					String name = removeSurroundingWhitespace(input.substring(0, equalIndex));
					String value = removeSurroundingWhitespace(input.substring(equalIndex+1, input.length()));
					nvpList.add(new NameValuePairSimple(name, value));
				}

			}
		}
		catch (IOException e) {}
		
		// Would rather typecast here, but not sure if I can
		INameValuePair[] nvp = new NameValuePairSimple[nvpList.size()];
		for(int i = 0; i < nvpList.size(); i++) {
			nvp[i] = (INameValuePair)(nvpList.get(i));
		}

		return nvp;
	}
	

}
