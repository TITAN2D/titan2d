package titan.gui;

import javax.swing.ListSelectionModel;

public class TabDischargePlanes extends TitanTabList {

	public TabDischargePlanes() {
		super(true, ListSelectionModel.SINGLE_SELECTION);
		
		addButton.setText("Add Discharge Plane");
		deleteButton.setText("Delete Discharge Plane");

		String[] headings = new String[] { "Point A UTM E",
                                           "Point A UTM N",
                                           "Point B UTM E",
                                           "Point b UTM N" };
        blankRow = new String[headings.length];
        for(int i = 0; i < headings.length; i++) {
        	blankRow[i] = "";
        }
		
		ids = new String[] { TitanConstants.PLANE_A_E,
				             TitanConstants.PLANE_A_N,
				             TitanConstants.PLANE_B_E,
				             TitanConstants.PLANE_B_N };
		tm.setColumnIdentifiers(headings);

	}

}
