package titan.io;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public abstract class OutputWriter {

	BufferedWriter stream;
	
	public OutputWriter(File outputFile) {
		try {
		    stream = new BufferedWriter(new FileWriter(outputFile));
		}
		catch (IOException e) {
			System.out.println("Problem with Output File : " + outputFile.getName());
			e.printStackTrace();
		}
	}

	public void close() {
		
		try {
			stream.close();
		}
		catch (IOException e) {
			System.out.println("Problem Closing stream.");
			e.printStackTrace();
		}
	}

}
