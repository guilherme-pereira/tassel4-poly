/*
 * AbstractTableReport
 */
package net.maizegenetics.pal.report;

/**
 *
 * @author terry
 */
public abstract class AbstractTableReport implements TableReport {
	
	private int currentRowNumber = -1;
	private Object[] currentRow = null;
	
    /**
     * get the data elements
     *
     * @return the data elements
     */
    public Object[][] getTableData() {
        return getTableData(0, getRowCount() - 1);
    }

    /**
     * Get Table Data in specified range inclusive.
     * Indices start at 0.
     *
     * @param start start position
     * @param end end position
     * @return
     */
    public Object[][] getTableData(int start, int end) {

        if ((start < 0) || (end >= getRowCount())) {
            throw new IndexOutOfBoundsException("getTableData: start: " + start + "  end: " + end);
        }

        if (end < start) {
            return null;
        }

        Object[][] temp = new Object[end - start + 1][];

        for (int i = start; i <= end; i++) {
            temp[i] = getRow(i);
        }

        return temp;

    }

    public Object getValueAt(int row, int col) {
    	if (row != currentRowNumber) {
            currentRowNumber = row;
            currentRow = getRow(row);
        }
        return currentRow[col];
    }

}
