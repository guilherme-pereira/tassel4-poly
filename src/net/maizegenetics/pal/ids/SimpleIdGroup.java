// SimpleIdGroup.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.pal.ids;

import java.io.Serializable;

import java.util.HashMap;
import java.util.List;

/**
 * SimpleIdGroup
 */
public class SimpleIdGroup implements IdGroup, Serializable {

    private Identifier[] ids;
    private HashMap<String, Integer> myIndices;
    //
    // Serialization code
    //
    private static final long serialVersionUID = -4266575329980153075L;

    /**
     * Constructor taking the size of the group.
     */
    public SimpleIdGroup(int size) {
        this(size, false);
    }

    /**
     * Constructor taking an array of strings.
     */
    public SimpleIdGroup(String[] labels) {
        this(labels.length);
        for (int i = 0; i < labels.length; i++) {
            setIdentifier(i, new Identifier(labels[i]));
        }
    }

    /**
     * Constructor taking the size of the group.
     *
     * @param size - the number of ids
     * @param createIDs - if true creates default Identifiers. Otherwise leaves
     * blank (for user to fill in)
     */
    public SimpleIdGroup(int size, boolean createIDs) {

        ids = new Identifier[size];
        myIndices = new HashMap<String, Integer>(size);
        if (createIDs) {
            for (int i = 0; i < size; i++) {
                setIdentifier(i, new Identifier("Taxa" + i));
            }
        }
    }

    /**
     * Constructor taking an array of identifiers.
     */
    public SimpleIdGroup(Identifier[] id) {
        this(id.length);
        for (int i = 0; i < id.length; i++) {
            setIdentifier(i, id[i]);
        }
    }

    public SimpleIdGroup(List<Identifier> ids) {
        this(ids.size());
        for (int i = 0, n = ids.size(); i < n; i++) {
            setIdentifier(i, ids.get(i));
        }
    }

    /**
     * Constructor taking two separate id groups and merging them.
     */
    public SimpleIdGroup(IdGroup a, IdGroup b) {
        this(a.getIdCount() + b.getIdCount());

        for (int i = 0; i < a.getIdCount(); i++) {
            setIdentifier(i, a.getIdentifier(i));
        }
        for (int i = 0; i < b.getIdCount(); i++) {
            setIdentifier(i + a.getIdCount(), b.getIdentifier(i));
        }
    }

    /**
     * Impersonating Constructor.
     */
    public static SimpleIdGroup getInstance(IdGroup a) {
        Identifier[] ids = new Identifier[a.getIdCount()];
        for (int i = 0; i < a.getIdCount(); i++) {
            ids[i] = a.getIdentifier(i);
        }
        return new SimpleIdGroup(ids);
    }

    /**
     * Impersonating Constructor.
     *
     * @param toIgnore - will ignore the identifier at the index specified by
     * toIgnore
     */
    public SimpleIdGroup(IdGroup a, int toIgnore) {
        this((toIgnore < 0 || toIgnore > a.getIdCount() ? a.getIdCount() : a.getIdCount() - 1));
        int index = 0;
        for (int i = 0; i < a.getIdCount(); i++) {
            if (i != toIgnore) {
                setIdentifier(index++, a.getIdentifier(i));
            }
        }
    }

    /**
     * Returns the number of identifiers in this group
     */
    @Override
    public int getIdCount() {
        return ids.length;
    }

    /**
     * Returns the ith identifier.
     */
    @Override
    public Identifier getIdentifier(int i) {
        return ids[i];
    }

    /**
     * Convenience method to return the name of identifier i
     */
    public final String getName(int i) {
        return ids[i].getName();
    }

    /**
     * Sets the ith identifier.
     */
    @Override
    public final void setIdentifier(int i, Identifier id) {
        ids[i] = id;
        myIndices.put(id.getFullName(), new Integer(i));
    }

    /**
     * Return index of identifier with name or -1 if not found
     */
    @Override
    public int whichIdNumber(String name) {

        Integer index = myIndices.get(name);
        if (index != null) {
            return index.intValue();
        }

        for (int i = 0, n = ids.length; i < n; i++) {
            if (ids[i].equals(name)) {
                return i;
            }
        }

        return -1;

    }

    /**
     * Returns a string representation of this IdGroup in the form of a
     * bracketed list.
     */
    @Override
    public String toString() {

        StringBuffer sb = new StringBuffer();
        sb.append("[ ");
        for (int i = 0; i < getIdCount(); i++) {
            sb.append(getIdentifier(i)).append(" ");
        }
        sb.append("]");
        return new String(sb);
    }

    @Override
    public int whichIdNumber(Identifier id) {
        return whichIdNumber(id.getFullName());
    }
}
