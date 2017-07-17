package net.maizegenetics.pal.report;

import java.io.Serializable;

import java.util.Iterator;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: ed
 * Date: Sep 28, 2006
 * Time: 9:37:46 PM
 * To change this template use File | Settings | File Templates.
 */
public class SimpleTableReport extends AbstractTableReport implements Serializable, TableReport {

    Object[][] theData;
    Object[] theColumnNames;
    Object[] theRowNames = null;
    String theName;

    //Implementation of TableReport Interface
    public SimpleTableReport(String theName, Object[] columnNames, Object[][] theData) {
        this.theData = theData;
        this.theColumnNames = columnNames;
        this.theName = theName;
    }

    public SimpleTableReport(TableReport tr) {
        this.theData = tr.getTableData();
        this.theColumnNames = tr.getTableColumnNames();
        this.theName = tr.getTableTitle();
    }

    /**
     * This method concatenates to rows together.  If the output has the same number of columns,
     * but the rows total the sum of the all the rows.
     */
    public static SimpleTableReport getInstance(List<SimpleTableReport> list) {

        if ((list == null) || (list.size() == 0)) {
            return null;
        }

        if (list.size() == 1) {
            return list.get(0);
        }

        int colNum = list.get(0).getTableColumnNames().length;
        int rowTotal = 0;
        Iterator<SimpleTableReport> itr = list.iterator();
        while (itr.hasNext()) {
            TableReport tr = itr.next();
            if (colNum != tr.getTableColumnNames().length) {
                return null;
            }
            rowTotal += tr.getTableData().length;
        }

        Object[][] fused = new Object[rowTotal][colNum];
        int currRow = 0;
        itr = list.iterator();
        while (itr.hasNext()) {
            Object[][] f = itr.next().getTableData();
            for (int i = 0; i < f.length; i++) {
                for (int j = 0; j < colNum; j++) {
                    fused[currRow][j] = f[i][j];
                }
                currRow++;
            }
        }

        return new SimpleTableReport(list.get(0).getTableTitle(), list.get(0).getTableColumnNames(), fused);

    }

    /**
     * Return column names for the table
     */
    public Object[] getTableColumnNames() {
        return theColumnNames;
    }

    /**
     * Return data for the table
     */
    public Object[][] getTableData() {
        return theData;
    }

    /**
     * Returns specified row.
     *
     * @param row row number
     *
     * @return row
     */
    public Object[] getRow(int row) {
        return theData[row];
    }

    /**
     * Return the name for the title of the ANOVA
     */
    public String getTableTitle() {
        return theName;
    }

    public int getRowCount() {
        return theData.length;
    }

    public int getElementCount() {
        return getRowCount() * getColumnCount();
    }

    public int getColumnCount() {
        return theColumnNames.length;
    }

    public void setRowNames(Object[] rowNames){
    	theRowNames = rowNames;
    }

    public Object getValueAt(int row, int col) {
        return theData[row][col];
    }
}
