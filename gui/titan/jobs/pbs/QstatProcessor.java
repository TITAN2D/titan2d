package titan.jobs.pbs;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

import titan.jobs.JobStatusProcessor;
/**
 * If	job status is being displayed in the default  format  and
 * the  -f  option  is  not  specified, the following	items are
 * displayed	on  a  single  line,  in  the  specified   order,
 *  separated by white	space:
 *
 *  -  the job identifier	assigned by PBS.
 *  -  the job name given	by the submitter.
 *  -  the job owner
 *  -  the CPU time used
 *  -  the job state:
 *       E -  Job	is exiting after having	run.
 *       H -  Job	is held.
 *       Q -  job	is queued, eligable to run or routed.
 *       R -  job	is running.
 *       T -  job	is being moved to new location.
 *       W -  job	is waiting for its execution time
 *	    (-a	option)	to be reached.
 *       S -  (Unicos only) job is suspend.
 *
 *  -  the queue in which	the job	resides
 * 
 */
public class QstatProcessor extends JobStatusProcessor {

	/**
	 * Get all jobs in the PBS job manager
	 * @return Array of objects containing information describing
	 *         job in queue of the PBS job manager.
	 */
	public QstatData[] getJobs() {
		QstatData[] jobData = null;
		try {
		    Process titanJob = Runtime.getRuntime().exec("qstat -a");
		    BufferedReader rs = new BufferedReader(new InputStreamReader(titanJob.getInputStream()));
		    String line;
		    ArrayList<String> jobArray = new ArrayList<String>();
		    while((line = rs.readLine()) != null) {
		    	line = removeEdgeWhitespace(line);
//System.out.println("QSTAT -> " + line);
                if(parseBlankSeparatedList(line).length == 11)
		    	    jobArray.add(line);
		    }
		    
		    jobData = new QstatData[jobArray.size()];
		    for(int job = 0; job < jobArray.size(); job++) {
		    	jobData[job] = new QstatData();
		    	String[] items = parseBlankSeparatedList(jobArray.get(job));
		    	if(items.length == 11) {
		    		for(int i = 0; i < items.length; i++)
		    		    items[i] = removeEdgeWhitespace(items[i]);
			    	jobData[job].jobID = items[0];
			    	jobData[job].jobName = items[3];
			    	jobData[job].user = items[1];
			    	jobData[job].useTime = items[10];
			    	jobData[job].status = items[9].charAt(0);
			    	jobData[job].queue = items[2];
		    	}
		    }

		}
		catch (IOException ex) { ex.printStackTrace(); }

		return jobData;
	}
	
	/**
	 * 
	 * @param jobID
	 * @return Returns information describing job in queue for the specified
	 *         job.
	 */
	public QstatData getJob(String jobID) {
		QstatData jobData = null;
		try {
		    Process titanJob = Runtime.getRuntime().exec("qstat -a " + jobID);
		    BufferedReader rs = new BufferedReader(new InputStreamReader(titanJob.getInputStream()));
		    String line;
		    ArrayList<String> jobArray = new ArrayList<String>();
		    while((line = rs.readLine()) != null) {
		    	line = removeEdgeWhitespace(line);
                // Skip any header information
                if(parseBlankSeparatedList(line).length == 11 &&
                   jobID.startsWith(parseBlankSeparatedList(line)[0]))
		    	    jobArray.add(line);
		    }
		    
		    if(jobArray.size() == 1) {
			    	jobData = new QstatData();
			    	String[] items = parseBlankSeparatedList(jobArray.get(0));
			    	if(items.length == 11) {
			    		for(int i = 0; i < items.length; i++)
			    		    items[i] = removeEdgeWhitespace(items[i]);
				    	jobData.jobID = items[0];
				    	jobData.jobName = items[3];
				    	jobData.user = items[1];
				    	jobData.useTime = items[10];
				    	jobData.status = items[9].charAt(0);
				    	jobData.queue = items[2];
			    	}

		    }

		}
		catch (IOException ex) { ex.printStackTrace(); }

		return jobData;
	}
	
	public class QstatData {
		public String jobID;
		public String jobName;
		public String user;
		public String useTime;
		public char status;
		public String queue;
	}
}
