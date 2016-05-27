package titan.jobs.condor;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

import titan.jobs.JobStatusProcessor;

public abstract class CondorJobStatusProcessor extends JobStatusProcessor {
	

	
	public class CondorQData {
		public String jobID;
		public String user;
		public String useTime;
		public char status;
	}
	
}
