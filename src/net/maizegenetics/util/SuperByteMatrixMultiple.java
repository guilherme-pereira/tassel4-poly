/*
 *  SuperByteMatrixMultiple
 */
package net.maizegenetics.util;

/**
 *
 * @author Terry Casstevens
 */
public class SuperByteMatrixMultiple implements SuperByteMatrix {

    private final byte[][] myData;
    private final int myNumRows;
    private final int myNumColumns;
    private final int myNumRowsPerSingleDimArray;

    SuperByteMatrixMultiple(int rows, int columns) {

        myNumRows = rows;
        myNumColumns = columns;

        long numElements = (long) myNumRows * (long) myNumColumns;
        myNumRowsPerSingleDimArray = Integer.MAX_VALUE / myNumColumns;
        int numElementsPerSingleDimArray = myNumRowsPerSingleDimArray * myNumColumns;
        int numSingleDimArrays = (int) (numElements / (long) numElementsPerSingleDimArray);
        int numRemaining = (int) (numElements % (long) numElementsPerSingleDimArray);
        if (numRemaining != 0) {
            myData = new byte[numSingleDimArrays + 1][];
            for (int i = 0; i < numSingleDimArrays; i++) {
                myData[i] = new byte[numElementsPerSingleDimArray];
            }
            myData[numSingleDimArrays] = new byte[numRemaining];
        } else {
            myData = new byte[numSingleDimArrays][];
            for (int i = 0; i < numSingleDimArrays; i++) {
                myData[i] = new byte[numElementsPerSingleDimArray];
            }
        }

    }

    @Override
    public void set(int row, int column, byte value) {
        myData[getFirstIndex(row)][getSecondIndex(row, column)] = value;
    }

    @Override
    public byte get(int row, int column) {
        return myData[getFirstIndex(row)][getSecondIndex(row, column)];
    }

    @Override
    public byte[] getAllColumns(int row) {

        if ((row < 0) || (row >= myNumRows)) {
            throw new IndexOutOfBoundsException("SuperByteMatrixMultiple: getAllColumns: row: " + row);
        }

        int start = getSecondIndex(row, 0);
        byte[] result = new byte[myNumColumns];
        System.arraycopy(myData[getFirstIndex(row)], start, result, 0, myNumColumns);
        return result;

    }

    @Override
    public byte[] getColumnRange(int row, int start, int end) {

        if ((row < 0) || (row >= myNumRows)) {
            throw new IndexOutOfBoundsException("SuperByteMatrixMultiple: getColumnRange: row: " + row);
        }

        if ((start < 0) || (start >= myNumColumns)) {
            throw new IndexOutOfBoundsException("SuperByteMatrixMultiple: getColumnRange: start: " + start);
        }

        if ((end < 0) || (end >= myNumColumns)) {
            throw new IndexOutOfBoundsException("SuperByteMatrixMultiple: getColumnRange: end: " + end);
        }

        if (end < start) {
            throw new IllegalArgumentException("SuperByteMatrixMultiple: getColumnRange: end: " + end + " less than start: " + start);
        }

        int startIndex = getSecondIndex(row, start);
        int numElements = end - start;
        byte[] result = new byte[numElements];
        System.arraycopy(myData[getFirstIndex(row)], startIndex, result, 0, numElements);
        return result;

    }

    @Override
    public byte[] getAllRows(int column) {

        if ((column < 0) || (column >= myNumColumns)) {
            throw new IndexOutOfBoundsException("SuperByteMatrixMultiple: getAllRows: column: " + column);
        }

        byte[] result = new byte[myNumRows];
        for (int row = 0; row < myNumRows; row++) {
            result[row] = get(row, column);
        }
        return result;

    }

    private int getFirstIndex(int row) {
        return row / myNumRowsPerSingleDimArray;
    }

    private int getSecondIndex(int row, int column) {
        return (row % myNumRowsPerSingleDimArray) * myNumColumns + column;
    }

    @Override
    public int getNumRows() {
        return myNumRows;
    }

    @Override
    public int getNumColumns() {
        return myNumColumns;
    }

    @Override
    public boolean isColumnInnerLoop() {
        return true;
    }
}
