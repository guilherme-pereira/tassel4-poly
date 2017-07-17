/*
 * TableReportNoPagingTableModel.java
 *
 */
package net.maizegenetics.gui;


import javax.swing.table.AbstractTableModel;

import net.maizegenetics.pal.report.TableReport;

/**
 *
 * @author  terryc
 */
public class TableReportNoPagingTableModel extends AbstractTableModel {

    private TableReport myTable = null;
    private Object[] myColumnHeadings = null;

    public TableReportNoPagingTableModel(TableReport table) {
        if (table == null) {
            throw new IllegalArgumentException("TableReportNoPagingTableModel: init: table can not be null.");
        }

        myTable = table;
        myColumnHeadings = myTable.getTableColumnNames();

    }

    // Return values appropriate for the visible table part
    public int getRowCount() {
        return myTable.getRowCount();
    }

    public int getColumnCount() {
        return myTable.getColumnCount();
    }

    public Object getValueAt(int row, int col) {
        return myTable.getValueAt(row, col);
    }

    public String getColumnName(int col) {
        return myColumnHeadings[col].toString();
    }

    /**
     * Resets the table backing this matrix table model to
     * an empty table.
     */
    public void resetTable() {
    }

    public Object getColumnObject(int columnIndex) {
        return myColumnHeadings[columnIndex];
    }

    /**
     * Always returns false.
     */
    public boolean isCellEditable(int rowIndex, int columnIndex) {
        return false;
    }

    /**
     * No operation.
     */
    public void setValueAt(Object aValue, int rowIndex, int columnIndex) {
        // NO OPERATION
    }

    public void fireTableChanged() {
        fireTableStructureChanged();
    }

    public Class<?> getColumnClass(int columnIndex) {
        return getValueAt(0, columnIndex).getClass();
    }
}
