// IdGenerator.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

package net.maizegenetics.pal.ids;

import net.maizegenetics.pal.io.FormattedOutput;

/**
 * Generates IdGroup objects given certain parameters. 
 * 
 * @version $Id: IdGenerator.java,v 1.1 2007/01/12 03:26:14 tcasstevens Exp $
 *
 * @author Alexei Drummond
 */
public class IdGenerator {

	/**
	 * generates a group of unique identifiers numbered from zero.
	 */
	public static IdGroup createIdGroup(int size) {
	
		int width = (int)Math.ceil(Math.log(size) / Math.log(10.0));
	
		IdGroup idGroup = new SimpleIdGroup(size);
	
		String name;
		for (int i = 0; i < size; i++) {
			name = (new Integer(i)).toString();
			name = FormattedOutput.space(width - name.length(), '0') + name;
			idGroup.setIdentifier(i, new Identifier(name));
		}
	
		return idGroup;
	}	
}

