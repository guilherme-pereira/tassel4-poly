// TableReport.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.pal.report;

/**
 * interface for classes with data that can be presented in tables
 *
 * @author Ed Buckler
 */
public interface TableReport {

    /**
     * get the names of the columns
     *
     * @return columns names
     */
    public Object[] getTableColumnNames();

    /**
     * get the data elements
     *
     * @return the data elements
     */
    public Object[][] getTableData();

    /**
     * get the title of the table
     *
     * @return a String title
     */
    public String getTableTitle();

    /**
     * get the number of the columns
     *
     * @return columns names
     */
    public int getColumnCount();

    /**
     * get the number of columns
     *
     * @return columns names
     */
    public int getRowCount();

    /**
     * get the total number of elements in the dataset. Elements=rowCount * columnCount;
     *
     * @return columns names
     */
    public int getElementCount();

    /**
     * Returns specified row.
     *
     * @param row row number
     * 
     * @return row
     */
    public Object[] getRow(int row);

    /**
     * Get Table Data in specified range inclusive.
     * Indices start at 0.
     *
     * @param start start position
     * @param end end position
     * @return
     */
    public Object[][] getTableData(int start, int end);

    /**
     * Returns value at given row and column.
     *
     * @param row row number
     * @param col column number
     * @return data
     */
    public Object getValueAt(int row, int col);
    
}
