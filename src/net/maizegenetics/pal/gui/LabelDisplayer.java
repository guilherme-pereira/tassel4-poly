// LabelDisplay.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

package net.maizegenetics.pal.gui;
import java.awt.*;

/**
 * A Label display object displays a label at a particular location on a graphics object.
 * The class may change attributes such as font, colour, etc... but should leave the graphics state unchanged (ie, revert back to original colour)
 * @author Matthew Goode
 * @note
 *     <ul>
 *       <li> 14 August 2003 - Created to allow greater flexibility in label display on trees
 *     </ul>
 */

public interface LabelDisplayer {

	public void display(Graphics g, String label, int x, int y );

	//=--=-=-=-=-=-==--=-=-=-=-==--==-=-=--==--==--=-==--=-==--==-=-=-=-=--=-==--=
	public static final class Utils {
		public static final LabelDisplayer buildDisplay(Color c) {
			return new ColourDisplay(c);
		}
		public static final LabelDisplayer buildDisplay(Color c, int style) {
			return new ColourAndFontStyleDisplay(c,style);
		}
		private static final class ColourDisplay implements LabelDisplayer {
			private final Color c_;
			public ColourDisplay(Color c) { this.c_ = c;	}
			public void display(Graphics g, String text, int x, int y ) {
				Color old = g.getColor();
				g.setColor(c_);
				g.drawString(text,x,y);
				g.setColor(old);
			}
		}
		private static final class ColourAndFontStyleDisplay implements LabelDisplayer {
			private final Color c_;
			private final int fontStyle_;
			public ColourAndFontStyleDisplay(Color c, int fontStyle) {
				this.c_ = c;
				this.fontStyle_ = fontStyle;
			}
			public void display(Graphics g, String text, int x, int y ) {
				Color old = g.getColor();
				Font oldFont = g.getFont();
				g.setColor(c_);
				System.out.println("Making bold:");
				Font newFont = oldFont.deriveFont(fontStyle_);
				g.setFont(newFont);
				g.drawString(text,x,y);
				g.setColor(old);
				g.setFont(oldFont);
			}
		}
	} //End of class Utils

}