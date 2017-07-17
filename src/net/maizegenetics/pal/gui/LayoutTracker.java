// LayoutTracker.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

package net.maizegenetics.pal.gui;

/**
 * Title:        LayoutTracer
 * Description:  A means for tracking the layout of labels
 * @author			 Matthew Goode
 * @version 1.0
 */
import net.maizegenetics.pal.ids.Identifier;

import java.awt.*;
import java.util.Hashtable;

public class LayoutTracker {
	Hashtable layoutMappings_ = new Hashtable();
	public LayoutTracker() { }
	public void addMapping(String name, Rectangle bounds) {
		layoutMappings_.put(name,bounds);
	}
	public void addMapping(Identifier id, Rectangle bounds) {
		if(id!=null&&id.getName()!=null) {
			layoutMappings_.put(id.getName(),bounds);
		}
	}

	public Rectangle getBounds(String name) {
		if(name==null||!layoutMappings_.containsKey(name)) {
			return null;
		}
		return (Rectangle)layoutMappings_.get(name);
	}
	public Rectangle getBounds(Identifier id) {
		if(id==null) {
			return null;
		}
		return getBounds(id.getName());
	}

	public void reset() {
		layoutMappings_.clear();
	}

}