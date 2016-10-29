package titan.gui;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.table.DefaultTableModel;

import titan.io.INameValuePair;
import titan.io.NameValuePairSimple;

public abstract class TitanTabList extends JPanel {

	protected JTable table;
	protected DefaultTableModel tm;
	protected String[] ids;
	protected JPanel buttonPanel;
	protected JButton addButton;
	protected JButton deleteButton;
	protected String[] blankRow;
	
	public TitanTabList(boolean isEditable) {
		this(isEditable, ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
	}
	
	/**
	 * 
	 * @param isEditable
	 * @param selectionMode selection mode (per ListSelectionModel.setSelectionMode).
	 */
	public TitanTabList(boolean isEditable, int selectionMode) {

		setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
		
		if(isEditable) {
			tm = new DefaultTableModel();
		}
		else {
		    tm = new DefaultTableModel() {
			    public boolean isCellEditable(int rowIndex, int columnIndex) {
				    return false;
			    }
		    }; 
		}

		table = new JTable(tm);
		table.setSelectionMode(selectionMode);
		
		JScrollPane tablePane = new JScrollPane(table);
		
		addButton = new JButton("Add Row");
		deleteButton = new JButton("Delete Row");
		buttonPanel = new JPanel();
		buttonPanel.add(addButton);
		buttonPanel.add(deleteButton);
		
		addButton.addActionListener(new AddEntry());
		deleteButton.addActionListener(new DeleteEntry());

		add(tablePane);
		add(buttonPanel);
//		add(addButton);
//		add(deleteButton);
    }

	public void setData(INameValuePair[][] nvp) {

		Object[] data = new Object[ids.length];
		tm.setRowCount(0);

		for(int row = 0; row < nvp.length; row++) {
			for(int col = 0; col < nvp[row].length; col++) {
				for(int id = 0; id < ids.length; id++) {
					if(nvp[row][col].getName().compareTo(ids[id]) == 0) {
						data[id] = nvp[row][col].getValue();
					}
				}
			}
			tm.addRow(data);
		}

	}
	
	public INameValuePair[][] getData() {
		
		INameValuePair[][] data = new INameValuePair[tm.getRowCount()][ids.length];
		
		for(int i = 0; i < tm.getRowCount(); i++) {
			for(int j = 0; j < ids.length; j++) {
				data[i][j] = new NameValuePairSimple(ids[j], (String)(tm.getValueAt(i,j)));
			}
		}
		return data;
	}
	
	private class AddEntry implements ActionListener {

		public void actionPerformed(ActionEvent e) {
			tm.addRow(blankRow);
		}
	}

	private class DeleteEntry implements ActionListener {

		public void actionPerformed(ActionEvent e) {
			if(table.getSelectedRow() != -1)
			    tm.removeRow(table.getSelectedRow());
		}
	}

}
