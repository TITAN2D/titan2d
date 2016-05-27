package titan.io;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

public class GenericReader extends InputReader {

	public GenericReader(File inputFile) throws FileNotFoundException {
		super(inputFile);
	}
	
	public String read() {
		String contents = new String("");
        String line;
        try {
			while((line = stream.readLine()) != null) {
				contents = contents.concat(line + "\n");
			}
        }
        catch (IOException e) { }
        	
		return contents;
		
	}
	
}
