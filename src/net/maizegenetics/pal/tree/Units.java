// Units.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

package net.maizegenetics.pal.tree;

/**
 * interface holding unit constants
 *
 * @version $Id: Units.java,v 1.1 2007/01/12 03:26:17 tcasstevens Exp $
 *
 * @author Alexei Drummond
 * @author Matthew Goode
 */
public interface Units
{
		int EXPECTED_SUBSTITUTIONS = 0;
		int GENERATIONS = 1;
		int DAYS = 2;
		int MONTHS = 3;
		int YEARS = 4;
		int UNKNOWN = 5;

		String[] UNIT_NAMES = {"Expected Substitutions per Site", "Generations", "Days", "Months", "Years", "Unknown"};
		String[] SHORT_UNIT_NAMES = {"Expected Substitutions", "Generations", "Days", "Months", "Years", "Unknown"};



}
