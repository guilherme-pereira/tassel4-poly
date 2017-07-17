/*
 * TableRowHeaderListModel.java
 *
 * Created on February 11, 2009, 9:09 PM
 *
 */
package net.maizegenetics.gui;

import java.util.List;

import javax.swing.*;

/**
 *
 * @author  terryc
 */
public class TableRowHeaderListModel extends AbstractListModel {

    private List myList = null;

    /**
     * TableRowHeaderListModel Constructor.
     *
     * @param list list of row headers
     */
    public TableRowHeaderListModel(List list) {
        super();

        myList = list;
    }

    /** Returns the value at the specified index.
     *
     * @param index the requested index
     *
     * @return the value at <code>index</code>
     *
     */
    public Object getElementAt(int index) {
        if ((myList == null) || myList.size() == 0) {
            return null;
        }

        return myList.get(index);
    }

    /**
     * Returns the length of the list.
     *
     * @return the length of the list
     *
     */
    public int getSize() {
        if ((myList == null) || myList.size() == 0) {
            return 0;
        }

        return myList.size();
    }
}

