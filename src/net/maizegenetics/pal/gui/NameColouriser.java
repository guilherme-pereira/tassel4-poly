// NameColouriser.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

package net.maizegenetics.pal.gui;

/**
 * Title:        NameColouriser.java
 * Description:  A means for mapping names to colours
 * @author      Matthew Goode
 * @version 1.0
 * @note
 *     <ul>
 *       <li> 14 August 2003 - Changed to use Label Displayers. Will need a name change at a later point
 *     </ul>
 */
import net.maizegenetics.pal.ids.Identifier;

import java.awt.*;
import java.util.Hashtable;
public final class NameColouriser implements java.io.Serializable {

	private final Hashtable displayMappings_ = new Hashtable();

	public NameColouriser() { }

	public NameColouriser(String name, Color colour) {
		addMapping(name,colour);
	}
	public NameColouriser(String[] names, Color colour) {
		for(int i = 0 ; i < names.length ; i++) {
			addMapping(names[i],colour);
		}
	}
	public NameColouriser(Identifier name, Color colour) {
		addMapping(name,colour);
	}

	public void addMapping(String name, Color colour) {
		displayMappings_.put(name,LabelDisplayer.Utils.buildDisplay(colour));
	}
	public void addMapping(String name, Color colour, int fontStyle) {
		displayMappings_.put(name,LabelDisplayer.Utils.buildDisplay(colour,fontStyle));
	}
	public void addMapping(String name, LabelDisplayer display) {
		displayMappings_.put(name, display);
	}

	public void addMapping(Identifier id, Color colour) {
		if(id!=null&&id.getName()!=null) {
			displayMappings_.put(id.getName(),LabelDisplayer.Utils.buildDisplay(colour));
		}
	}

	/**
	 * @param names Names
	 * @param colours associated colours
	 * @note assumes parallel arrays
	 */
	public void addMappings(String[] names, Color[] colours) {
		for(int i = 0 ; i < names.length ; i++) {
			displayMappings_.put(names[i],LabelDisplayer.Utils.buildDisplay(colours[i]));
		}
	}

	public LabelDisplayer getDisplay(String name, LabelDisplayer defaultDisplay) {
		if(name==null||!displayMappings_.containsKey(name)) {
			return defaultDisplay;
		}
		return (LabelDisplayer)displayMappings_.get(name);
	}

	public LabelDisplayer getDisplay(Identifier id, LabelDisplayer defaultDisplay) {
		if(id==null) {
			return defaultDisplay;
		}
		return getDisplay(id.getName(),defaultDisplay);
	}

}
