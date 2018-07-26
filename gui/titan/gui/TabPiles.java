package titan.gui;

import javax.swing.ListSelectionModel;
import javax.swing.JButton;
import javax.swing.JOptionPane;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

public class TabPiles extends TitanTabList {

    public TabPiles(boolean TWOPHASES_COULOMB_MODEL) {
        super(true, ListSelectionModel.SINGLE_SELECTION);

        addButton.setText("Add Pile");
        deleteButton.setText("Delete Pile");

        JButton VolumeButton = new JButton("Calculate Volume");
        buttonPanel.add(VolumeButton);
        VolumeButton.addActionListener(new CalculateVolume());

        setHeadings(TWOPHASES_COULOMB_MODEL);
    }

    public void setHeadings(boolean TWOPHASES_COULOMB_MODEL) {

        if (TWOPHASES_COULOMB_MODEL == true) {
            String[] headings = new String[]{"Type (1:P 2:C)",
                    "Maximum Initial Thickness",
                    "Center of Initial Volume X",
                    "Center of Initial Volume Y",
                    "Major Extent",
                    "Minor Extent",
                    "Orientation Angle",
                    "Initial Speed",
                    "Initial Direction",
                    "Volume Fraction"};

            blankRow = new String[headings.length];
            for (int i = 0; i < headings.length; i++) {
                blankRow[i] = "";
            }

            ids = new String[]{
                    TitanConstants.PILE_TYPE,
                    TitanConstants.PILE_MAX_INIT_THICKNESS,
                    TitanConstants.PILE_CENTER_INIT_VOLUME_XC,
                    TitanConstants.PILE_CENTER_INIT_VOLUME_YC,
                    TitanConstants.PILE_MAJOR_EXTENT,
                    TitanConstants.PILE_MINOR_EXTENT,
                    TitanConstants.PILE_ORIENTATION_ANGLE,
                    TitanConstants.PILE_INIT_SPEED,
                    TitanConstants.PILE_INIT_DIRECTION,
                    TitanConstants.PILE_VOLUME_FRACTION};
            tm.setColumnIdentifiers(headings);

        } else {

            String[] headings = new String[]{"Type (1:P 2:C)",
                    "Maximum Initial Thickness",
                    "Center of Initial Volume X",
                    "Center of Initial Volume Y",
                    "Major Extent",
                    "Minor Extent",
                    "Orientation Angle",
                    "Initial Speed",
                    "Initial Direction"};

            blankRow = new String[headings.length];
            for (int i = 0; i < headings.length; i++) {
                blankRow[i] = "";
            }

            ids = new String[]{
                    TitanConstants.PILE_TYPE,
                    TitanConstants.PILE_MAX_INIT_THICKNESS,
                    TitanConstants.PILE_CENTER_INIT_VOLUME_XC,
                    TitanConstants.PILE_CENTER_INIT_VOLUME_YC,
                    TitanConstants.PILE_MAJOR_EXTENT,
                    TitanConstants.PILE_MINOR_EXTENT,
                    TitanConstants.PILE_ORIENTATION_ANGLE,
                    TitanConstants.PILE_INIT_SPEED,
                    TitanConstants.PILE_INIT_DIRECTION};
            tm.setColumnIdentifiers(headings);
        }
    }

    private class CalculateVolume implements ActionListener {
        public void actionPerformed(ActionEvent e) {

            if (table.getSelectedRowCount() != 1) {
                JOptionPane.showMessageDialog(TabPiles.this,
                        "No row selected",
                        "Calculate Volume Error",
                        JOptionPane.ERROR_MESSAGE);
                return;
            }

            int index = table.getSelectedRow();

            if (((String) (tm.getValueAt(index, 1))).isEmpty()) {
                JOptionPane.showMessageDialog(TabPiles.this,
                        "No Max Initial Thickness for selected row!",
                        "Calculate Volume Error",
                        JOptionPane.ERROR_MESSAGE);
                return;
            }

            if (((String) (tm.getValueAt(index, 4))).isEmpty()) {
                JOptionPane.showMessageDialog(TabPiles.this,
                        "No Major Extent for selected row!",
                        "Calculate Volume Error",
                        JOptionPane.ERROR_MESSAGE);
                return;
            }

            if (((String) (tm.getValueAt(index, 5))).isEmpty()) {
                JOptionPane.showMessageDialog(TabPiles.this,
                        "No Minor Extent for selected row!",
                        "Calculate Volume Error",
                        JOptionPane.ERROR_MESSAGE);
                return;
            }

            float maxInitThickness = Float.parseFloat((String) (tm.getValueAt(index, 1)));
            float majExtent = Float.parseFloat((String) (tm.getValueAt(index, 4)));
            float minExtent = Float.parseFloat((String) (tm.getValueAt(index, 5)));

            float vol = ((float) Math.PI) * maxInitThickness * majExtent * minExtent * ((float) .5);

            JOptionPane.showMessageDialog(TabPiles.this,
                    "Volume: " + vol + "  [m^3]",
                    "Calculated Volume",
                    JOptionPane.INFORMATION_MESSAGE);
            return;

        }
    }
}

