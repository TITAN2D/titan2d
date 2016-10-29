package titan.gui;

import javax.swing.BoxLayout;

import titan.graphics.DirectorySelector;
import titan.graphics.TextInput;
import titan.graphics.RadioButtonGroup;
import titan.io.INameValuePair;
import titan.io.NameValuePairComponent;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;


public class TabGIS extends TitanTab {

	private RadioButtonGroup format;
	private DirectorySelector mainDir;
	private DirectorySelector subDir;
	private DirectorySelector mapSet;
	private DirectorySelector map;
	private TextInput vector;
	private TextInput zoneOverride;
	private RadioButtonGroup hemisphere;
	private TextInput minX;
	private TextInput minY;
	private TextInput maxX;
	private TextInput maxY;

	public TabGIS() {

		setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));

		values = new NameValuePairComponent[12];

		// GIS Format
		format = new RadioButtonGroup("GIS Format", TitanConstants.GisFormats);
		values[0] = new NameValuePairComponent(TitanConstants.GIS_FORMAT, format);

		// GIS Main Directory
		mainDir = new DirectorySelector("GIS Information Main Directory");
		mainDir.setChooserStyle(DirectorySelector.DIRECTORIES_ONLY);
		values[1] = new NameValuePairComponent(TitanConstants.GIS_INFO_DIRECTORY, mainDir);

		// GIS Sub-directory
		subDir = new DirectorySelector("GIS Sub-Directory");
		subDir.setOutputStyle(DirectorySelector.DIR_ONLY);
		subDir.setChooserStyle(DirectorySelector.DIRECTORIES_ONLY);
		values[2] = new NameValuePairComponent(TitanConstants.GIS_SUBDIR, subDir);

		// GIS Map Set
		mapSet = new DirectorySelector("GIS Map Set");
		mapSet.setOutputStyle(DirectorySelector.DIR_ONLY);
		mapSet.setChooserStyle(DirectorySelector.DIRECTORIES_ONLY);
		values[3] = new NameValuePairComponent(TitanConstants.GIS_MAPSET, mapSet);
		
		// GIS Map
		map = new DirectorySelector("GIS Map");
		map.setOutputStyle(DirectorySelector.DIR_ONLY);
		map.setChooserStyle(DirectorySelector.DIRECTORIES_ONLY);
		values[4] = new NameValuePairComponent(TitanConstants.GIS_MAP, map);
		
		// GIS Vector
		vector = new TextInput("GIS Vector");
		values[5] = new NameValuePairComponent(TitanConstants.GIS_VECTOR, vector);

		int bSize = 0;
		bSize = bSize > mainDir.getButtonWidth() ? bSize : mainDir.getButtonWidth();
		bSize = bSize > subDir.getButtonWidth() ? bSize : subDir.getButtonWidth();
		bSize = bSize > mapSet.getButtonWidth() ? bSize : mapSet.getButtonWidth();
		bSize = bSize > map.getButtonWidth() ? bSize : map.getButtonWidth();
		mainDir.setButtonWidth(bSize);
		subDir.setButtonWidth(bSize);
		mapSet.setButtonWidth(bSize);
		map.setButtonWidth(bSize);

		// Zone Override and Hemisphere are used to create the zone.txt file required
		// for creating KML files when the KML button is clicked on the JobDetailsDialog tab.
		// Please see the Job Submission tab for more details

		// Zone
		zoneOverride = new TextInput("Zone Override");
		values[6] = new NameValuePairComponent(TitanConstants.ZONEOVERRIDE, zoneOverride);

		// Hemisphere
		hemisphere = new RadioButtonGroup("Hemisphere", TitanConstants.HemisphereHalves);
		values[7] = new NameValuePairComponent(TitanConstants.HEMISPHERE, hemisphere);

		// min x location
		minX = new TextInput("Minimum X Location [UTM E]");
		values[8] = new NameValuePairComponent(TitanConstants.MIN_X_LOC, minX);

		// min y location
		minY = new TextInput("Minimum Y Location [UTM N]");
		values[9] = new NameValuePairComponent(TitanConstants.MIN_Y_LOC, minY);

		// max x location
		maxX = new TextInput("Maximum X Location [UTM E]");
		values[10] = new NameValuePairComponent(TitanConstants.MAX_X_LOC, maxX);

		// max y location
		maxY = new TextInput("Maximum Y Location [UTM N]");
		values[11] = new NameValuePairComponent(TitanConstants.MAX_Y_LOC, maxY);

		// Add components to main panel
		add(format.getPanel());
		add(mainDir.getPanel());
		add(subDir.getPanel());
		add(mapSet.getPanel());
		add(map.getPanel());
		add(vector.getPanel());
		add(zoneOverride.getPanel());
		add(hemisphere.getPanel());
		add(minX.getPanel());
		add(minY.getPanel());
		add(maxX.getPanel());
		add(maxY.getPanel());

		format.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent ev) {
				mainDir.setEditable(true);
				subDir.setEditable(true);
				mapSet.setEditable(true);
				vector.setEditable(true);
			}
		}, TitanConstants.GisFormats[TitanConstants.GIS_FORMAT_GIS_GRASS]);

		format.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent ev) {
				mainDir.setEditable(false);
				subDir.setEditable(false);
				mapSet.setEditable(false);
				vector.setEditable(false);
			}
		}, TitanConstants.GisFormats[TitanConstants.GIS_FORMAT_GDAL]);
	}

	public void setData(INameValuePair[] data) {
		super.setData(data);

		if(format.getValue().compareTo(TitanConstants.GisFormats[TitanConstants.GIS_FORMAT_GIS_GRASS]) == 0) {
			mainDir.setEditable(true);
			subDir.setEditable(true);
			mapSet.setEditable(true);
			vector.setEditable(true);
		} else {
			mainDir.setEditable(false);
			subDir.setEditable(false);
			mapSet.setEditable(false);
			vector.setEditable(false);
		}
	}
}
