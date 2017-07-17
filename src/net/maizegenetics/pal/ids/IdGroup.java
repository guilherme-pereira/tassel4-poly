// IdGroup.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.pal.ids;

/**
 * An indexed group of identifiers. For example of group of taxa
 * related by a phylogenetic tree.
 * <BR><B>NOTE:</B> Was called Taxa but not general enough.
 *
 * @version $Id: IdGroup.java,v 1.4 2009/07/07 16:19:37 tcasstevens Exp $
 *
 * @author Alexei Drummond
 */
public interface IdGroup extends java.io.Serializable {

    /**
     * Returns the number of identifiers in this group
     */
    public int getIdCount();

    /**
     * Returns the ith identifier.
     */
    public Identifier getIdentifier(int i);

    /**
     * Sets the ith identifier.
     */
    public void setIdentifier(int i, Identifier id);

    /**
     * Returns the index of the identifier with the given name.
     * Name should be the identifiers 'full name'.
     *
     * @param name full name of identifier
     *
     * @return index (-1 if not found)
     */
    public int whichIdNumber(String name);

    /**
     * Returns the index of the identifier.
     *
     * @param id identifier
     * @return index (-1 if not found)
     */
    public int whichIdNumber(Identifier id);
}
