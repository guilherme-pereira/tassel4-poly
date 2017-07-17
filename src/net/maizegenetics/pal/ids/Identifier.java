// Identifier.java
//
// (c) 1999-2000 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.pal.ids;

import java.io.Serializable;

import net.maizegenetics.prefs.TasselPrefs;

/**
 * An identifier for some sampled data. This will most often be for example, the
 * accession number of a DNA sequence, or the taxonomic name that the sequence
 * represents, et cetera.
 *
 * @author terry
 */
public class Identifier implements Serializable, Comparable {

    private static final long serialVersionUID = -7873729831795750538L;
    public static final String DELIMITER = ":";
    private final String myName;
    private final String[] myNameTokens;
    public static Identifier ANONYMOUS = new Identifier("");
    private final int hashCode;

    public Identifier(String name) {
        myName = name;
        myNameTokens = name.split(DELIMITER);
        hashCode = myName.hashCode();
    }

    public static Identifier getMergedInstance(Identifier id1, Identifier id2) {
        String[] first = id1.getFullNameTokens();
        String[] second = id2.getFullNameTokens();
        int count = Math.min(first.length, second.length);
        for (int i = 0; i < count; i++) {
            if (!first[i].equals(second[i])) {
                StringBuilder builder = new StringBuilder();
                for (int x = 0; x < i; x++) {
                    if (x != 0) {
                        builder.append(DELIMITER);
                    }
                    builder.append(first[x]);
                    return new Identifier(builder.toString());
                }
            }
        }
        return id1;
    }

    public String toString() {
        return getName();
    }

    // implements Comparable interface
    @Override
    public int compareTo(Object c) {
        if (this == c) {
            return 0;
        } else if (c instanceof Identifier) {
            return compareTo(((Identifier) c).getFullNameTokens());
        } else {
            throw new ClassCastException();
        }
    }

    @Override
    public boolean equals(Object c) {

        if (this == c) {
            return true;
        }

        TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES type = TasselPrefs.getIDJoinStrict();

        if (c instanceof Identifier) {
            if (type == TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES.Strict) {
                return getFullName().equals(((Identifier) c).getFullName());
            } else {
                return compareTo(c) == 0;
            }
        } else if (c instanceof String) {
            if (type == TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES.Strict) {
                return getFullName().equals(((String) c));
            } else {
                return compareTo((String) c) == 0;
            }
        } else {
            return false;
        }

    }

    public int compareTo(String c) {
        return compareTo(c.split(DELIMITER));
    }

    public int compareTo(String[] fullNameTokens) {

        TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES type = TasselPrefs.getIDJoinStrict();

        int count = Math.min(myNameTokens.length, fullNameTokens.length);
        if (type == TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES.NumLevels) {
            count = Math.min(count, TasselPrefs.getIDJoinNumLevels());
        }
        for (int i = 0; i < count; i++) {
            int current = myNameTokens[i].compareTo(fullNameTokens[i]);
            if (current != 0) {
                return current;
            }
        }

        if (type == TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES.Strict) {
            if (myNameTokens.length < fullNameTokens.length) {
                return -1;
            } else if (fullNameTokens.length < myNameTokens.length) {
                return 1;
            } else {
                return 0;
            }
        } else {
            return 0;
        }

    }

    public String getName() {
        return getNameLevel(0);
    }

    public String getFullName() {
        return myName;
    }

    public String getFullName(String delimiter) {
        if (delimiter.equals(DELIMITER)) {
            return myName;
        }
        return myName.replaceAll(DELIMITER, delimiter);
    }

    public String[] getFullNameTokens() {
        return myNameTokens;
    }

    /**
     * Returns requested level of name starting at index 0. 0 will generally be
     * most specific classification.
     *
     * @param index
     * @return Specified level.
     */
    public String getNameLevel(int index) {
        if (index < myNameTokens.length) {
            return myNameTokens[index];
        }
        return null;
    }

    /**
     * Returns name up to specified level (not including specified level. Levels
     * start at index 0.
     *
     * @param index
     * @return name up to specified level exclusive.
     */
    public String getNameToLevel(int index) {
        return getNameToLevel(index, DELIMITER);
    }

    public String getNameToLevel(int index, String delimiter) {

        int upto = 0;
        if (index > myNameTokens.length) {
            upto = myNameTokens.length;
        } else {
            upto = index;
        }
        if (upto == 0) {
            return null;
        }

        StringBuilder builder = new StringBuilder();
        for (int i = 0; i < upto; i++) {
            if (i != 0) {
                builder.append(delimiter);
            }
            builder.append(myNameTokens[i]);
        }

        return builder.toString();
    }

    /**
     * Returns number of name levels.
     *
     * @return number of name levels.
     */
    public int getNumNameLevels() {
        return myNameTokens.length;
    }

    @Override
    public int hashCode() {
        return hashCode;
    }
}
