/*
 * DataSet.java
 *
 */
package net.maizegenetics.plugindef;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

/**
 * This is a set of Datum.
 */
public class DataSet implements Serializable {

    private static final long serialVersionUID = -5197800047652332969L;

    private final List<Datum> myList;
    private final Plugin myCreator;

    /**
     * Creates a new instance of Datum
     *
     * @param list list of data elements.
     * @param creator creating plugin.
     *
     */
    public DataSet(List<Datum> list, Plugin creator) {

        if (list == null) {
            myList = new ArrayList<Datum>();
        } else {
            myList = new ArrayList<Datum>(list);
        }
        myCreator = creator;

    }

    /**
     * Convenience constructor for a single Datum
     *
     * @param theDatum a single Datum
     * @param creator creating plugin.
     *
     */
    public DataSet(Datum theDatum, Plugin creator) {
        if (theDatum == null) {
            myList = new ArrayList<Datum>();
        } else {
            myList = new ArrayList<Datum>();
            myList.add(theDatum);
        }
        myCreator = creator;
    }

    /**
     * Creates a new instance of Datum
     *
     * @param list list of data elements.
     * @param creator creating plugin.
     *
     */
    public DataSet(Datum[] list, Plugin creator) {

        if (list == null) {
            myList = new ArrayList<Datum>();
        } else {
            myList = new ArrayList<Datum>(Arrays.asList(list));
        }
        myCreator = creator;

    }

    /**
     * Combines multiple data sets.
     *
     * @param list list of data sets.
     * @param creator creating plugin.
     *
     */
    public static DataSet getDataSet(List<DataSet> list, Plugin creator) {

        List temp = new ArrayList();

        if (list != null) {
            Iterator itr = list.iterator();
            while (itr.hasNext()) {
                DataSet current = (DataSet) itr.next();
                if (current != null) {
                    temp.addAll(current.getDataSet());
                }
            }
        }

        return new DataSet(temp, creator);

    }

    public List getDataSet() {
        return new ArrayList(myList);
    }

    public String toString() {

        StringBuilder buffer = new StringBuilder();
        buffer.append("DataSet...\n");
        if (getCreator() != null) {
            buffer.append("Creator: ");
            buffer.append(getCreator().getClass().getName());
            buffer.append("\n");
        }

        Iterator itr = myList.iterator();
        while (itr.hasNext()) {
            Datum current = (Datum) itr.next();
            buffer.append("name: ");
            buffer.append(current.getName());
            buffer.append("  type: ");
            buffer.append(current.getDataType());
            buffer.append("\n");

        }

        return buffer.toString();

    }

    public Datum getData(int i) {
        if (i < myList.size()) {
            return myList.get(i);
        } else {
            return null;
        }
    }

    public int getSize() {
        return myList.size();
    }

    public Plugin getCreator() {
        return myCreator;
    }

    public List<Datum> getDataOfType(Class theClass) {

        if (theClass == null) {
            return null;
        }

        List<Datum> list = new ArrayList<Datum>();
        Iterator itr = myList.iterator();
        while (itr.hasNext()) {
            Datum current = (Datum) itr.next();

            if (theClass.isInstance(current.getData())) {
                list.add(current);
            }
        }

        return list;

    }

    public List<Datum> getDataOfType(Class[] classes) {

        if ((classes == null) || (classes.length == 0)) {
            return new ArrayList(myList);
        }

        List<Datum> result = null;

        for (int i = 0, n = classes.length; i < n; i++) {

            List current = getDataOfType(classes[i]);
            if (result == null) {
                result = current;
            } else {
                result.addAll(current);
            }

        }

        return result;

    }

    public List<Datum> getDataWithName(String name) {

        if (name == null) {
            return null;
        }

        List<Datum> list = new ArrayList<Datum>();
        Iterator itr = myList.iterator();
        while (itr.hasNext()) {
            Datum current = (Datum) itr.next();

            if (name.equals(current.getName())) {
                list.add(current);
            }
        }

        return list;

    }

    public List<Datum> getDataWithName(String[] names) {

        if ((names == null) || (names.length == 0)) {
            return new ArrayList(myList);
        }

        List<Datum> result = null;

        for (int i = 0, n = names.length; i < n; i++) {

            List current = getDataWithName(names[i]);
            if (result == null) {
                result = current;
            } else {
                result.addAll(current);
            }

        }

        return result;

    }

    public List<Datum> getDataOfTypeWithName(Class[] classes, String[] names) {
        return (new DataSet(getDataOfType(classes), null)).getDataWithName(names);
    }
}
