package titan.io;

import java.io.File;
import java.io.IOException;

public class ParameterListOutput extends OutputWriter {

	public ParameterListOutput(File outputFile) {
        super(outputFile);
	}
	
	public void writeData(INameValuePair[][] data) {
		String output;
		try {
			for(int i = 0; i < data.length; i++) {
				output = new String("/=" + i + "\n");
				stream.write(output);
				for(int j = 0; j < data[i].length; j++) {
				    output = new String(data[i][j].getName() + "=" + data[i][j].getValue() + "\n");
				    stream.write(output);
				}
			}
		}
		catch(IOException e) { }
	}
}
