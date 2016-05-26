package titan.jobs.condor;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

public class CondorHistoryProcessor extends CondorJobStatusProcessor {
	
	private static final int MY_ARG_LEN = 8;
	
	public CondorQData[] getJobs() {
		CondorQData[] jobData = null;
		try {
		    Process titanJob = Runtime.getRuntime().exec("condor_history");
		    BufferedReader rs = new BufferedReader(new InputStreamReader(titanJob.getInputStream()));
		    String line;
		    ArrayList<String> jobArray = new ArrayList<String>();
		    while((line = rs.readLine()) != null) {
		    	line = removeEdgeWhitespace(line);
	            if(parseBlankSeparatedList(line).length == MY_ARG_LEN)
		    	    jobArray.add(line);
		    }
		    
		    jobData = new CondorQData[jobArray.size()];
		    for(int job = 0; job < jobArray.size(); job++) {
		    	jobData[job] = new CondorQData();
		    	String[] items = parseBlankSeparatedList(jobArray.get(job));
		    	if(items.length == MY_ARG_LEN) {
		    		for(int i = 0; i < items.length; i++)
		    		    items[i] = removeEdgeWhitespace(items[i]);
			    	jobData[job].jobID = items[0];
			    	jobData[job].user = items[1];
			    	jobData[job].useTime = items[4];
			    	jobData[job].status = items[5].charAt(0);
		    	}
		    }
	
		}
		catch (IOException ex) { ex.printStackTrace(); }
	
		return jobData;
	}
	
	public CondorQData getJob(String jobID) {
		CondorQData jobData = null;
		try {
		    Process titanJob = Runtime.getRuntime().exec("condor_history " + jobID);
		    BufferedReader rs = new BufferedReader(new InputStreamReader(titanJob.getInputStream()));
		    String line;
		    ArrayList<String> jobArray = new ArrayList<String>();
		    while((line = rs.readLine()) != null) {
		    	line = removeEdgeWhitespace(line);
	            if(parseBlankSeparatedList(line).length == MY_ARG_LEN)
		    	    jobArray.add(line);
		    }
		    
		    if(jobArray.size() == 1) {
		    	jobData = new CondorQData();
		    	String[] items = parseBlankSeparatedList(jobArray.get(0));
		    	if(items.length == MY_ARG_LEN) {
		    		for(int i = 0; i < items.length; i++)
		    		    items[i] = removeEdgeWhitespace(items[i]);
			    	jobData.jobID = items[0];
			    	jobData.user = items[1];
			    	jobData.useTime = items[4];
			    	jobData.status = items[5].charAt(0);
		    	}
		    }
	
		}
		catch (IOException ex) { ex.printStackTrace(); }
	
		return jobData;
	}
	
}
