/*
 * IdGroupUtils
 */
package net.maizegenetics.pal.ids;

import java.util.ArrayList;
import java.util.List;
import java.util.TreeSet;

/**
 *
 * @author terry
 */
public final class IdGroupUtils {

    private IdGroupUtils() {
        // To prevent instantiation
    }

    /**
     * @return true if <i>sub</i> IdGroup completely contained within <i>full</i>, false otherwise
     */
    public static boolean isContainedWithin(IdGroup sub, IdGroup full) {
        for (int i = 0; i < sub.getIdCount(); i++) {
            boolean found = false;
            Identifier subID = sub.getIdentifier(i);
            for (int j = 0; j < full.getIdCount(); j++) {
                Identifier fullID = full.getIdentifier(j);
                if (fullID.equals(subID)) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                return false;
            }
        }
        return true;
    }

    /**
     * @return true if <i>id1</i> and <i>id2</i> share exactly the same identifiers (.equals() based, not reference base). The order is not important.
     */
    public static boolean isEqualIgnoringOrder(IdGroup id1, IdGroup id2) {
        return (isContainedWithin(id1, id2) && isContainedWithin(id2, id1));
    }

    /**
     * A convenience implementation of whichIdNumber that can be used by
     * IdGroup implementations
     * @return -1 if <i>s</i> not in <i>group</i>
     */
    public static int whichIdNumber(IdGroup group, String s) {
        for (int i = 0; i < group.getIdCount(); i++) {
            if (s.equals(group.getIdentifier(i).getFullName())) {
                return i;
            }
        }
        return -1;
    }

    /**
     * @param group1 an IdGroup
     * @param group2 another IdGroup
     * @return the Ids in the intersection of groups 1 and 2, sorted in ascending order
     */
    public static IdGroup getCommonIds(IdGroup group1, IdGroup group2) {
        return getCommonIds(new IdGroup[]{group1, group2});
    }

    /**
     * Intersect joins the specified groups.
     *
     * @param groups groups to join.
     * @return The ids from the intersect join, sorted in ascending order
     */
    public static IdGroup getCommonIds(IdGroup[] groups) {

        if ((groups == null) || (groups.length == 0)) {
            return null;
        }

        TreeSet<Identifier> intersectIds = new TreeSet<Identifier>();
        for (int x = 0; x < groups[0].getIdCount(); x++) {
            intersectIds.add(groups[0].getIdentifier(x));
        }
        for (int i = 1; i < groups.length; i++) {
            List temp = new ArrayList();
            for (int j = 0; j < groups[i].getIdCount(); j++) {
                temp.add(groups[i].getIdentifier(j));
            }
            intersectIds.retainAll(temp);
        }

        Identifier[] ids = new Identifier[intersectIds.size()];
        intersectIds.toArray(ids);

        return new SimpleIdGroup(ids);

    }

    /**
     * @param group1	an IdGroup
     * @param group2	another IdGroup
     * @return	the Ids in the union of groups 1 and 2, sorted in ascending order
     */
    public static IdGroup getAllIds(IdGroup group1, IdGroup group2) {
        return getAllIds(new IdGroup[]{group1, group2});
    }

    /**
     * Union joins the specified groups.
     *
     * @param groups groups to join.
     * @return The ids from the union join, sorted in ascending order
     */
    public static IdGroup getAllIds(IdGroup[] groups) {

        if ((groups == null) || (groups.length == 0)) {
            return null;
        }

        TreeSet<Identifier> allIds = new TreeSet<Identifier>();
        for (int i = 0; i < groups.length; i++) {
            int n = groups[i].getIdCount();
            for (int j = 0; j < n; j++) {
                allIds.add(groups[i].getIdentifier(j));
            }
        }

        Identifier[] ids = new Identifier[allIds.size()];
        allIds.toArray(ids);

        return new SimpleIdGroup(ids);

    }

    /**
     * @param original an IdGroup
     * @param include a boolean array with the same number of elements as the IdGroup original.
     * If the value of include is true for the ith element then that Identifier will be included
     * in the new subset. Otherwise, it will not.
     * @return a new IdGroup that is a subset of the original
     * @throws IllegalArgumentException if the number of Identifiers in original is not equal to the length of include.
     */
    public static IdGroup idGroupSubset(IdGroup original, boolean[] include) {
        int nOld = original.getIdCount();
        if (nOld != include.length) {
            throw new IllegalArgumentException("Size of IdGroup and include array are different.");
        }
        ArrayList<Identifier> newIds = new ArrayList<Identifier>();
        for (int i = 0; i < nOld; i++) {
            if (include[i]) {
                newIds.add(original.getIdentifier(i));
            }
        }
        return new SimpleIdGroup(newIds.toArray(new Identifier[newIds.size()]));
    }

    /**
     * Translates an array of identifiers into an array of strings
     */
    public static String[] getNames(Identifier[] ids) {
        String[] names = new String[ids.length];
        for (int i = 0; i < names.length; i++) {
            names[i] = ids[i].getName();
        }
        return names;
    }

    /**
     * Translates an array of identifiers into an array of strings, with optional removal of particular identifier
     * @param toIgnoreIndex the index of an idetifier to ignore, if <0 no element is ignored
     */
    public static String[] getNames(Identifier[] ids, int toIgnore) {
        if (toIgnore < 0 || toIgnore >= ids.length) {
            return getNames(ids);
        }
        String[] names = new String[ids.length - 1];
        int index = 0;
        for (int i = 0; i < names.length; i++) {
            if (i != toIgnore) {
                names[index] = ids[i].getName();
                index++;
            }
        }
        return names;
    }

    /**
     * Translates an an array of strings into an array of identifiers
     */
    public static Identifier[] getIdentifiers(String[] names) {
        Identifier[] ids = new Identifier[names.length];
        for (int i = 0; i < names.length; i++) {
            ids[i] = new Identifier(names[i]);
        }
        return ids;
    }

    /**
     * Translates an IdGroup into an array of identifiers
     */
    public static Identifier[] getIdentifiers(IdGroup idGroup) {
        Identifier[] ids = new Identifier[idGroup.getIdCount()];
        for (int i = 0; i < ids.length; i++) {
            ids[i] = idGroup.getIdentifier(i);
        }
        return ids;
    }

    /**
     * Translates an IdGroup into an array of strings
     */
    public static String[] getNames(IdGroup ids) {
        String[] names = new String[ids.getIdCount()];
        for (int i = 0; i < names.length; i++) {
            names[i] = ids.getIdentifier(i).getName();
        }
        return names;
    }

    /**
     * Translates an IDgroup into an array of strings, with optional removal of particular identifier
     * @param toIgnoreIndex the index of an idetifier to ignore, if <0 no element is ignored
     */
    public static String[] getNames(IdGroup ids, int toIgnore) {
        if (toIgnore < 0 || toIgnore >= ids.getIdCount()) {
            return getNames(ids);
        }
        int numberOfIDS = ids.getIdCount();
        String[] names = new String[numberOfIDS - 1];
        int index = 0;
        for (int i = 0; i < numberOfIDS; i++) {
            if (i != toIgnore) {
                names[index] = ids.getIdentifier(i).getName();
                index++;
            }
        }
        return names;
    }
}
