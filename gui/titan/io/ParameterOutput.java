package titan.io;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;

public class ParameterOutput extends OutputWriter {
	
	public ParameterOutput(File outputFile) {
        super(outputFile);
	}

	public void writeData(INameValuePair[] data) {
		String output;
		for(int i = 0; i < data.length; i++) {
			output = new String(data[i].getName() + "=" + data[i].getValue() + "\n");
			try {
			    stream.write(output);
			}
			catch(IOException e) {}
		}
	}

}
