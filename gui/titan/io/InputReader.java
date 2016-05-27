package titan.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import javax.swing.JOptionPane;

public abstract class InputReader {

	BufferedReader stream;
	
	public InputReader(File inputFile) throws FileNotFoundException {
//		try {
		    stream = new BufferedReader(new FileReader(inputFile));
//		}
//		catch (FileNotFoundException e) { 
//			throw(e);
//		}
	}

	protected String removeSurroundingWhitespace(String str) {
		// Remove leading white space
        str = str.replaceAll("\\s+", "");
        // Remove trailing whitespace
        str = str.replaceAll("\\s+$", "");
		return str;
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
