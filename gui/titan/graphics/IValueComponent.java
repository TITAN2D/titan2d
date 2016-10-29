package titan.graphics;

import java.awt.Component;

public interface IValueComponent {

	public Component getPanel();

	public String getValue();

	public void setValue(String val);
}
