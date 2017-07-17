/*
 *  SuperByteMatrixSingleTranspose
 */
package net.maizegenetics.util;

/**
 *
 * @author Terry Casstevens
 */
public class SuperByteMatrixTranspose implements SuperByteMatrix {

    private final SuperByteMatrix myMatrix;

    SuperByteMatrixTranspose(int rows, int columns) {
        myMatrix = SuperByteMatrixBuilder.getInstance(columns, rows);
    }

    @Override
    public int getNumRows() {
        return myMatrix.getNumColumns();
    }

    @Override
    public int getNumColumns() {
        return myMatrix.getNumRows();
    }

    @Override
    public void set(int row, int column, byte value) {
        myMatrix.set(column, row, value);
    }

    @Override
    public byte get(int row, int column) {
        return myMatrix.get(column, row);
    }

    @Override
    public byte[] getAllColumns(int row) {
        return myMatrix.getAllRows(row);
    }

    @Override
    public byte[] getColumnRange(int row, int start, int end) {
        int length = end - start;
        byte[] result = new byte[length];
        for (int i = 0; i < length; i++) {
            result[i] = get(row, i);
        }
        return result;
    }

    @Override
    public byte[] getAllRows(int column) {
        return myMatrix.getAllColumns(column);
    }

    @Override
    public boolean isColumnInnerLoop() {
        return false;
    }
}
