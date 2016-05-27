package titan.io;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;

public class ParameterListInput extends InputReader {
	
	public ParameterListInput(File inputFile) throws FileNotFoundException {
		super(inputFile);
	}
	
	public INameValuePair[][] readData() {
		String input;
		int maxNumParms = 0;
		int numParms = 0;
		int equalIndex;
		ArrayList list = new ArrayList();
		ArrayList nvpList = null;

		try {
			while((input = stream.readLine()) != null) {
				removeSurroundingWhitespace(input);
				if(input.indexOf("/=") != -1) {
					nvpList = new ArrayList();
					list.add(nvpList);
					numParms = 0;
				}
				else if((equalIndex = input.indexOf("=")) != -1) {
					String name = removeSurroundingWhitespace(input.substring(0, equalIndex));
					String value = removeSurroundingWhitespace(input.substring(equalIndex+1, input.length()));
					nvpList.add(new NameValuePairSimple(name, value));
					numParms++;
					if(numParms > maxNumParms) {
						maxNumParms = numParms;
					}
				}
			}
		}
		catch (IOException e) {}
		
		// Would rather typecast here, but not sure if I can
		INameValuePair[][] nvp = new NameValuePairSimple[list.size()][maxNumParms];
		for(int i = 0; i < list.size(); i++) {
			ArrayList al = (ArrayList)(list.get(i));
			for(int j = 0; j < al.size(); j++) {
			  nvp[i][j] = (INameValuePair)(al.get(j));
			}
		}

		return nvp;
	}
	
}
