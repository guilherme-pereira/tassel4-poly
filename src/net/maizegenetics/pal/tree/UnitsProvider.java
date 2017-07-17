// UnitsProvider.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

package net.maizegenetics.pal.tree;



/**
 * interface for classes that can provide the related Units used, (as
 *
 * @version $Id: UnitsProvider.java,v 1.1 2007/01/12 03:26:17 tcasstevens Exp $
 *
 * @author Matthew Goode
 */
public interface UnitsProvider extends Units {

	/**
	 * @return the units relating to this object.
	 */
	int getUnits();
}
