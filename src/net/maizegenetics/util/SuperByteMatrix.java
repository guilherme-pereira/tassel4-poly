/*
 *  SuperByteMatrix
 */
package net.maizegenetics.util;

/**
 *
 * @author Terry Casstevens
 */
public interface SuperByteMatrix {

    /**
     * Return number of rows.
     *
     * @return number of rows
     */
    public int getNumRows();

    /**
     * Return number of columns.
     *
     * @return number of columns
     */
    public int getNumColumns();

    /**
     * Sets value at give row and column.
     *
     * @param row row
     * @param column column
     * @param value value
     */
    public void set(int row, int column, byte value);

    /**
     * Gets value at given row and column.
     *
     * @param row row
     * @param column column
     *
     * @return value
     */
    public byte get(int row, int column);

    /**
     * Get all values for given row.
     *
     * @param row row
     *
     * @return values
     */
    public byte[] getAllColumns(int row);

    /**
     * Get values for given row from start column (inclusive) to end column
     * (exclusive).
     *
     * @param row row
     * @param start start
     * @param end end
     *
     * @return values
     */
    public byte[] getColumnRange(int row, int start, int end);

    /**
     * Get all values for give column.
     *
     * @param column column
     *
     * @return values
     */
    public byte[] getAllRows(int column);

    /**
     * Returns true if the matrix stored for better performance when column loop
     * inside row loop. False if matrix stored for better performance when
     * row loop inside column loop.
     *
     * @return true if the matrix stored for better performance when column loop
     * inside row loop. False if matrix stored for better performance when
     * row loop inside column loop.
     */
    public boolean isColumnInnerLoop();
}
