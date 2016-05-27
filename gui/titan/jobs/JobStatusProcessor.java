package titan.jobs;

import java.util.ArrayList;

public abstract class JobStatusProcessor {

    protected String removeEdgeWhitespace(String line) {
        String removed = line.replaceAll("^\\s+", "");
        removed = removed.replaceAll("\\s+$", "");
        return removed;
    }

    /**
     * Breaks the given string down into individual strings that are separated
     * by a blank space.
     * @param line
     * @return
     */
    protected String[] parseBlankSeparatedList(String line) {
        ArrayList myList = new ArrayList();

        String myLine = new String(line);

        while(myLine.contains(" ")) {
                String item = myLine.substring(0, myLine.indexOf(" "));
                myList.add(item);
                myLine = myLine.substring((myLine.indexOf(" ") + 1));
                myLine = removeEdgeWhitespace(myLine);
        }
        myLine = removeEdgeWhitespace(myLine);
        myList.add(myLine);

        String[] list = new String[myList.size()];
        for(int i = 0; i < myList.size(); i++) {
                list[i] = new String((String)myList.get(i));
        }

        return list;
    }

}
