package titan.gui;

import javax.swing.ListSelectionModel;

public class TabFluxSources extends TitanTabList {

	public TabFluxSources() {
		super(true, ListSelectionModel.SINGLE_SELECTION);
				
		addButton.setText("Add Flux Source");
		deleteButton.setText("Delete Flux Source");

		String[] headings = new String[] { "Extrusion Rate",
                                           "Active Time Start",
                                           "Active Time End",
                                           "Center of Source X",
                                           "Center of Source Y",
                                           "Major Extent",
                                           "Minor Extent",
                                           "Orientation Angle",
                                           "Init Speed",
                                           "Init Direction" };

        blankRow = new String[headings.length];
        for(int i = 0; i < headings.length; i++) {
        	blankRow[i] = "";
        }
		
		ids = new String[] { TitanConstants.SRC_EXTRUSION_FLUX_RATE,
				             TitanConstants.SRC_ACTIVE_TIME_START,
                             TitanConstants.SRC_ACTIVE_TIME_END,
				             TitanConstants.SRC_CENTER_XC,
			                 TitanConstants.SRC_CENTER_YC,
			                 TitanConstants.SRC_MAJOR_EXTENT,
			                 TitanConstants.SRC_MINOR_EXTENT,
			                 TitanConstants.SRC_ORIENTATION_ANGLE,
			                 TitanConstants.SRC_INIT_SPEED,
			                 TitanConstants.SRC_INIT_DIRECTION };
		tm.setColumnIdentifiers(headings);

	}

}
