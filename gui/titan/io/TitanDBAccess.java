package titan.io;

import java.io.File;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

public class TitanDBAccess {

	private static final String driver = "org.apache.derby.jdbc.EmbeddedDriver";
	private String dbName;
	private String tableName;
	Connection conn;
	
	public TitanDBAccess(String db) throws ClassNotFoundException {
		dbName = new String(db);
		tableName = new String("TitanJobMonitor");

		try {
	        Class.forName(driver);
		}
		catch (ClassNotFoundException ex) {
			throw new ClassNotFoundException(
					"Titan is unable to find the job details database jar " +
	                " file.\nNeed to add the derby.jar file to the CLASSPATH.\n" +
	                "Exiting...");
		}

	}
	
	/**
	 * Open a connection to the job monitor database.
	 */
	public void open() throws SQLException {

			File dbFile = new File(dbName);
			String connectionURL;
			boolean dbExists = false;
			if(dbFile.exists()) dbExists = true;
			
			if(dbExists)
				connectionURL = new String("jdbc:derby:" + dbName + ";");
			else
			    connectionURL = new String("jdbc:derby:" + dbName + ";create=true");

			try {
			conn = DriverManager.getConnection(connectionURL);
			
			// If the database did not previously exists, create the job
			// monitor table
			if(dbExists == false) {
				Statement createTable = conn.createStatement();
				String createString = new String("CREATE table " + tableName +
						" (Job_Name VARCHAR(128) NOT NULL, Job_ID VARCHAR(128), " +
						" Submit_Method VARCHAR(128), Submit_Host VARCHAR(128)," +
						" Date_Submitted VARCHAR(128), Output_Directory VARCHAR(128))");
				createTable.execute(createString);
			}

			// The database has been successfully opened or created.
			// Save a backup of the database
		}
		catch (SQLException ex) {
			throw new SQLException(
					"There was a problem establishing a connection to the" +
					" job details database: " + dbName + ".\n" +
					"Two instances of Titan2D cannot access the job details database simultaneously.\n" +
					"Please verify an instance of Titan2D is not already running.");
		}
	}
	
	/**
	 * Add a job to the job monitor database.
	 * @param jm Contains information specific to the job to add to
	 *           the job monitor database.
	 */
	public void addJob(TitanJobMonitor jm) throws SQLException {
		String insertString = "";
		try {
		    Statement addCommand = conn.createStatement();
		    insertString = new String("insert into " + tableName +
		    		" values('" + jm.jobName +"', '" + jm.jobID + "', '" +
		    				jm.submitMethod + "', '" + jm.submitHost + "', '" +
		    				jm.dateSubmitted + "', '" + jm.outputDirectory + "')");
//System.out.println(insertString);
		    addCommand.execute(insertString);
		}
		catch (SQLException ex) {
			throw new SQLException(
					"Error while attempting to add an entry in the job " +
					"details database.\nThe error is not terminal, but " +
					"the job will not be displayed on the job monitoring " +
					"list.\nInsert String :\n" + insertString);
		}
		
	}
	/**
	 * 
	 * @param jm The jobName and jobID fields are used to delete any
	 *           matching entries in the job monitor database.  Both
	 *           fields must match.
	 */
	public void removeJob(TitanJobMonitor jm) throws SQLException {
		try {
		    Statement delCommand = conn.createStatement();
		    String deleteString = new String("delete from " + tableName +
		    		" where Job_Name='" + jm.jobName +"' and Job_ID='" + jm.jobID + "'");
//System.out.println(deleteString);
		    delCommand.execute(deleteString);
		}
		catch (SQLException ex) {
			throw new SQLException(
					"Error while attempting to remove an entry in the job " +
					"details database.\nThe error is not terminal, but " +
					"the job will continue to be displayed on the job\n" +
					"monitoring list.");
		}
	}

	/**
	 * 
	 * @return Returns an array containing all jobs for the current
	 *         user that are logged in the job monitor database.  A job may
	 *         not be in the database if the user has deleted the job
	 *         entry in the database.
	 */
	public TitanJobMonitor[] getAllJobs() throws SQLException {
		TitanJobMonitor[] allJobs = null;
		ResultSet results;
		try {
		    Statement queryCommand = conn.createStatement(ResultSet.TYPE_SCROLL_SENSITIVE, ResultSet.CONCUR_READ_ONLY);
		    String selectString = new String("select * from " + tableName + " order by Job_Name");
		    results = queryCommand.executeQuery(selectString);
		    
			results.last();
		    allJobs = new TitanJobMonitor[results.getRow()];
//System.out.println("Number of jobs : " + results.getRow());
		    results.beforeFirst();
		    
		    int idx = 0;
		    while(results.next()) {
		    	allJobs[idx] = new TitanJobMonitor();
		    	allJobs[idx].jobName = new String(results.getString(1));
		    	allJobs[idx].jobID = new String(results.getString(2));
		    	allJobs[idx].submitMethod = new String(results.getString(3));
		    	allJobs[idx].submitHost = new String(results.getString(4));
		    	allJobs[idx].dateSubmitted = new String(results.getString(5));
		    	allJobs[idx].outputDirectory = new String(results.getString(6));
//System.out.println("  " + allJobs[idx].jobName + "  " + allJobs[idx].jobID);
		    	idx++;
		    }
		}
		catch (SQLException ex) {
			throw new SQLException(
					"Error while attempting to retrieve job detail information " +
					"from the database.\nThe error is not terminal, but " +
					"job information will not to be displayed on the job\n" +
					"monitoring list.");
		}


		return allJobs;
	}
	
	/**
	 * 
	 * @return Returns an array containing all jobs for the current
	 *         user that are logged in the job monitor database.  A job may
	 *         not be in the database if the user has deleted the job
	 *         entry in the database.
	 */
	public TitanJobMonitor getJobByName(String jobName) throws SQLException {
		TitanJobMonitor job = new TitanJobMonitor();
		ResultSet results;
		try {
			if(conn == null) return null;
		    Statement queryCommand = conn.createStatement(ResultSet.TYPE_SCROLL_SENSITIVE, ResultSet.CONCUR_READ_ONLY);
		    String selectString = new String("select * from " + tableName +
		    		" where Job_Name='" + jobName + "'");
//System.out.println("Query command : " + selectString);
		    results = queryCommand.executeQuery(selectString);
		    
			results.last();
			if(results.getRow() != 1) {
				job = null;
			}
			else {
//System.out.println("Number of jobs : " + results.getRow());
			    results.first();
			    
		    	job.jobName = new String(results.getString(1));
		    	job.jobID = new String(results.getString(2));
		    	job.submitMethod = new String(results.getString(3));
		    	job.submitHost = new String(results.getString(4));
		    	job.dateSubmitted = new String(results.getString(5));
		    	job.outputDirectory = new String(results.getString(6));
			}
//System.out.println("  " + allJobs[idx].jobName + "  " + allJobs[idx].jobID);

		}
		catch (SQLException ex) {
			throw new SQLException(
					"Error while attempting to retrieve job detail information " +
					"for job " + jobName + ".\nThe error is not terminal, but " +
					"job information for the job will not to be displayed on the\n" +
					"job monitoring list.");
		}


		return job;
	}
	/**
	 * Close the connection to the job monitor database.
	 */
	public void close() throws SQLException {
		String connectionURL = new String("jdbc:derby:;shutdown=true");
		boolean gotSQLExc = false;
		try {
		    conn = DriverManager.getConnection(connectionURL);
		}
		catch (SQLException ex) {
			
			if(ex.getSQLState().equals("XJ015")) {
				gotSQLExc = true;
			}
		}
		if(!gotSQLExc) {
			throw new SQLException(
					"Unable to close the job details database.\n" +
					"Shutdown String : " + connectionURL);
		}
	}
	
	public class TitanJobMonitor {
		public String jobName;
		public String jobID;
		public String submitMethod;
		public String submitHost;
		public String dateSubmitted;
		public String status;
		public String outputDirectory;
	}
	
}
