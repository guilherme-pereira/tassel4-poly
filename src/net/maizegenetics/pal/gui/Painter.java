// Painter.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

package net.maizegenetics.pal.gui;

import java.awt.*;


/**
 * An interface for objects which are used for painting themselves in a
 * general Graphics object
 *
 * @author Matthew Goode
 */

public interface Painter {
	void paint(Graphics g, int displayWidth, int displayHeight);
	void paint(Graphics g, int displayWidth, int displayHeight, LayoutTracker lt);
	Dimension getPreferredSize();

}
