package titan.gui;

import java.awt.Cursor;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.sql.SQLException;
import java.util.TimerTask;

import javax.swing.JButton;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JLabel;
import javax.swing.ListSelectionModel;
import javax.swing.Timer;
import javax.swing.table.TableColumn;
import javax.swing.border.TitledBorder;
import javax.swing.JTextArea;
import javax.swing.JScrollPane;
import javax.swing.BoxLayout;

import titan.io.INameValuePair;
import titan.io.NameValuePairSimple;
import titan.io.TitanDBAccess;
import titan.io.TitanDBAccess.TitanJobMonitor;
import titan.jobs.JobSubmissionContainer;


public class TabJobMonitor extends TitanTabList {

	private TitanDBAccess dba;
	private JobSubmissionContainer submitContainer;
	private RefreshJobList refreshList;
	private RefreshJobMonitor refreshAction;
	private JPopupMenu tableItemMenu;

	public TabJobMonitor(TitanDBAccess db, JobSubmissionContainer sc) {
        super(false, ListSelectionModel.SINGLE_SELECTION);

	//Note: to enable selection of multiple rows use:
	//super(false, ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
        //In this case, job details will show details for the first selected job only

		dba = db;
		submitContainer = sc;
	    
	        refreshList = new RefreshJobList();
	    
		addComponentListener(new OnDisplay());
		
		buttonPanel.remove(addButton);
		buttonPanel.remove(deleteButton);
		remove(buttonPanel);

		JButton refreshButton = new JButton("Refresh");
		refreshAction = new RefreshJobMonitor();
		refreshButton.addActionListener(refreshAction);
		
		JButton removeFromListButton = new JButton("Remove From List");
		removeFromListButton.addActionListener(new RemoveJob());
		JButton jobDetailsButton = new JButton("Job Details");
		jobDetailsButton.addActionListener(new ShowJobDetails());
		JButton deleteJobButton = new JButton("Delete Job");
		deleteJobButton.addActionListener(new DeleteJob());

		JPanel myButtonPanel = new JPanel();
		myButtonPanel.add(refreshButton);
		myButtonPanel.add(removeFromListButton);
		myButtonPanel.add(jobDetailsButton);
		myButtonPanel.add(deleteJobButton);
		add(myButtonPanel);
		
		
		// Help box
		JPanel helpPanel = new JPanel();
		helpPanel.setLayout(new BoxLayout(helpPanel, BoxLayout.Y_AXIS));
		helpPanel.setBorder(new TitledBorder("Button Explanations"));
		
		JTextArea textArea = new JTextArea();
		textArea.setText("First select a row in the above table by clicking on it, " +
		"and then use the buttons to preform actions on that specific Titan run.\n\n" + 
		"REFRESH - Refreshes the table of jobs\n\n" + 
		"REMOVE FROM LIST - Removes the selected job from the table " + 
		"(does not delete contents of job run)\n\n" +
		"JOB DETAILS - Shows the associated files of the job selected in the table\n\n" +
		"DELETE JOB - Only works for a job actively running");
		
		textArea.setLineWrap(true);
		textArea.setWrapStyleWord(true);
		textArea.setEditable(false);
		
		JScrollPane scrollingArea = new JScrollPane(textArea);
		helpPanel.add(scrollingArea);
		
		// add JPanel's
		add(helpPanel);
		
		String[] headings = new String[] { "Job Name",
                                           "Job ID",
                                           "Date Submitted",
                                           "Job Run Method",
                                           "Host Submitted From",
                                           "Status" };

		// Convoluted scheme used to determine column widths
		// arrived at by trial and error.
		int[] colSize = new int[] { 100, 100, 50, 10, 50, 20 };
        blankRow = new String[headings.length];
        for(int i = 0; i < headings.length; i++) {
        	blankRow[i] = "";
        }

		ids = new String[] { "JOB_NAME", "JOB_ID", "DATE_SUBMITTED", "RUN_METHOD",
				             "SUBMIT_HOST", "JOB_STATUS" };
		tm.setColumnIdentifiers(headings);

	    for(int i = 0; i < tm.getColumnCount(); i++) {
		    TableColumn column = table.getColumnModel().getColumn(i);
		    column.setPreferredWidth(colSize[i]);
	    }

	    this.addComponentListener(new DisplayListener());

	    tableItemMenu = new JPopupMenu();
	    JMenuItem removeItem = new JMenuItem("Remove From List");
	    JMenuItem detailsItem = new JMenuItem("Job Details");
	    JMenuItem deleteItem = new JMenuItem("Delete Job");
	    tableItemMenu.add(removeItem);
	    tableItemMenu.add(detailsItem);
	    tableItemMenu.add(deleteItem);
	    removeItem.addActionListener(new RemoveJob());
	    detailsItem.addActionListener(new ShowJobDetails());
	    deleteItem.addActionListener(new DeleteJob());
	    table.addMouseListener(new TableMouseListener());
	    
	}
	
	/**
	 * Listens for right mouse button down press events and selects the
	 * table row where the right mouse button down press occurs on and
	 * displays the table row item menu popup.
	 * @author rrassler
	 *
	 */
	private class TableMouseListener extends MouseAdapter {

		public void mousePressed(MouseEvent e) {
			if(e.getButton() == MouseEvent.BUTTON3) {
                table.setRowSelectionInterval(table.rowAtPoint(e.getPoint()), table.rowAtPoint(e.getPoint()));
				tableItemMenu.show(table, e.getPoint().x, e.getPoint().y);
			}
		}

	}
	
	private class DisplayListener extends ComponentAdapter {
		Timer timer = null;

		public void componentHidden(ComponentEvent e) {
			if(timer != null) timer.stop();
		}
		public void componentShown(ComponentEvent e) {
			// When the timer is canceled, a new timer and timer
			// task need to be created.  The garbage collector
			// should clean up memory for old timer/task.
			// (ms)
			timer = new Timer(300000, refreshAction);
			timer.setInitialDelay(60000);
			timer.setRepeats(true);
			timer.start();
		}
	}
	
	private void refreshJobList() {
		new Thread(refreshList).start();
	}
	
	private class RefreshJobMonitor implements ActionListener {
		
		// this actionPerformed might want to do its work in another thread
		// if possible
		public void actionPerformed(ActionEvent e) {
            refreshJobList();
		}
	}
	
	private class RemoveJob implements ActionListener {
		
		public void actionPerformed(ActionEvent e) {

			boolean jobRemoved = false;
			if(table.getSelectedRowCount() > 0) {
				int[] selectedRows = table.getSelectedRows();
				for(int i = selectedRows.length-1; i >=0 ; i--) {

				    // Remove entry in database
				    TitanJobMonitor jm = dba.new TitanJobMonitor();
				    jm.jobName = new String((String)tm.getValueAt(selectedRows[i], 0));
				    jm.jobID = new String((String)tm.getValueAt(selectedRows[i], 1));
				    try {
	                    dba.removeJob(jm);
	                    jobRemoved = true;
				    }
					catch (SQLException ex) {
						JOptionPane.showMessageDialog(TabJobMonitor.this,
								                   ex.getMessage(),
								                   "Titan Database Error",
								                   JOptionPane.WARNING_MESSAGE);
					}
		            // Remove row from table
//		            tm.removeRow(selectedRows[i]);
				}
				if(jobRemoved) refreshJobList();
		    }
		}
	}
	
	private class ShowJobDetails implements ActionListener {
		
		// this actionPerformed might want to do its work in another thread
		// if possible
		public void actionPerformed(ActionEvent e) {
			if(table.getSelectedRow() != -1) {
				// Remove entry in database
				try {
				    TitanJobMonitor jm = dba.getJobByName(
				    		(String)tm.getValueAt(table.getSelectedRow(), 0));

					if(jm == null) {
						JOptionPane optionPane = new JOptionPane();
						optionPane.showMessageDialog(TabJobMonitor.this,
								"Problem retrieving status for job " +
								(String)tm.getValueAt(table.getSelectedRow(), 0),
								"Job Details Error",
		                        JOptionPane.ERROR_MESSAGE);
					}
					else {
				        JobDetailsDialog details = new JobDetailsDialog(TabJobMonitor.this, jm, submitContainer);
				        details.setVisible(true);
				    }
				}
				catch (SQLException ex) {
					JOptionPane.showMessageDialog(TabJobMonitor.this,
							ex.getMessage(),
	                        "Job Details Error",
	                        JOptionPane.ERROR_MESSAGE);
				}
		    }
		}
	}
	
	private class DeleteJob implements ActionListener {
		
		public void actionPerformed(ActionEvent e) {

			boolean jobDeleted = false;
			if(table.getSelectedRowCount() > 0) {
				int[] selectedRows = table.getSelectedRows();
				for(int i = selectedRows.length-1; i >=0 ; i--) {

				    // Remove entry in database
					try {
					    TitanJobMonitor jm = dba.getJobByName(
					    		(String)tm.getValueAt(selectedRows[i], 0));
	                    if(submitContainer.deleteJob(jm.submitMethod, jm.jobID, jm.outputDirectory))
	                    	jobDeleted = true;
	                    else
	    				    JOptionPane.showMessageDialog(TabJobMonitor.this,
	    						"The job " + jm.jobID + " is not in a state\n" +
	                            "where it can be removed.",
	                            "Job Deletion Error",
	                            JOptionPane.INFORMATION_MESSAGE);
					}
					catch (SQLException ex) {
						JOptionPane.showMessageDialog(TabJobMonitor.this,
								                   ex.getMessage(),
								                   "Titan Database Error",
								                   JOptionPane.WARNING_MESSAGE);
					}
				}
				if(jobDeleted) refreshJobList();
		    }
		}
	}
	
	private class OnDisplay extends ComponentAdapter {
		public void componentShown(ComponentEvent e) {
			refreshJobList();
		}
	}
	
	private class RefreshJobList extends TimerTask {
		
		public synchronized void run() {

			Cursor origCursor = TabJobMonitor.this.getCursor();
			TabJobMonitor.this.setCursor(new Cursor(Cursor.WAIT_CURSOR));

	        try {
				TitanJobMonitor[] jobData = dba.getAllJobs();
				INameValuePair[][] data = new INameValuePair[jobData.length][ids.length];
				for(int job = 0; job < jobData.length; job++) {
					data[job][0] = new NameValuePairSimple(ids[0], jobData[job].jobName);
					data[job][1] = new NameValuePairSimple(ids[1], jobData[job].jobID);
					data[job][2] = new NameValuePairSimple(ids[2], jobData[job].dateSubmitted);
					data[job][3] = new NameValuePairSimple(ids[3], jobData[job].submitMethod);
					data[job][4] = new NameValuePairSimple(ids[4], jobData[job].submitHost);
					
				    data[job][5] = new NameValuePairSimple(ids[5],
				    		submitContainer.getJobStatus(jobData[job].submitMethod,
				    				jobData[job].jobID, jobData[job].outputDirectory));
				}
				setData(data);
	        }
	        catch (SQLException ex) {
	        	setCursor(origCursor);
				JOptionPane.showMessageDialog(TabJobMonitor.this,
						ex.getMessage(),
	                    "Job Monitor Error",
	                    JOptionPane.ERROR_MESSAGE);
	        }
	        setCursor(origCursor);
		}
	}
}
