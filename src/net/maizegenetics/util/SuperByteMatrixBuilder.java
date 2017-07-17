/*
 *  SuperByteMatrixBuilder
 */
package net.maizegenetics.util;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 *
 * @author Terry Casstevens
 */
public class SuperByteMatrixBuilder {

    private static final int TRANSPOSE_BLOCK_SIZE = 64;

    private SuperByteMatrixBuilder() {
    }

    /**
     * This returns a SuperByteMatrix designed for better performance when
     * column iteration loop inside row iteration loop.
     *
     * @param numRows number of rows
     * @param numColumns number of columns
     *
     * @return SuperByteMatrix (double dimension byte array)
     */
    public static SuperByteMatrix getInstance(int numRows, int numColumns) {
        long numElements = (long) numRows * (long) numColumns;
        if (numElements > (long) (Integer.MAX_VALUE - 10)) {
            return new SuperByteMatrixMultiple(numRows, numColumns);
        } else {
            return new SuperByteMatrixSingle(numRows, numColumns);
        }
    }

    public static SuperByteMatrix getInstanceCopy(SuperByteMatrix matrix) {
        int numRows = matrix.getNumRows();
        int numColumns = matrix.getNumColumns();
        SuperByteMatrix result = getInstance(numRows, numColumns);
        for (int r = 0; r < numRows; r++) {
            for (int c = 0; c < numColumns; c++) {
                result.set(r, c, matrix.get(r, c));
            }
        }
        return result;
    }

    /**
     * This returns a SuperByteMatrix designed for better performance when row
     * iteration loop inside column iteration loop.
     *
     * @param numRows number of rows
     * @param numColumns number of columns
     *
     * @return SuperByteMatrix (double dimension byte array)
     */
    public static SuperByteMatrix getInstanceTranspose(int numRows, int numColumns) {
        return new SuperByteMatrixTranspose(numRows, numColumns);
    }

    /**
     * This returns a SuperByteMatrix that performs better in the reverse row /
     * column iteration nesting than the given matrix.
     *
     * @param matrix input matrix
     *
     * @return SuperByteMatrix (double dimension byte array)
     */
    public static SuperByteMatrix getInstanceTranspose(SuperByteMatrix matrix) {


        int numRows = matrix.getNumRows();
        int numColumns = matrix.getNumColumns();
        SuperByteMatrix result;
        int numThreads = Runtime.getRuntime().availableProcessors();
        ExecutorService pool = Executors.newFixedThreadPool(numThreads);

        if ((matrix instanceof SuperByteMatrixSingle) || (matrix instanceof SuperByteMatrixMultiple)) {
            result = getInstanceTranspose(numRows, numColumns);
            int rowBlockSize = TRANSPOSE_BLOCK_SIZE;
            for (int rowOffset = 0; rowOffset < numRows; rowOffset += TRANSPOSE_BLOCK_SIZE) {
                if (numRows - rowOffset < TRANSPOSE_BLOCK_SIZE) {
                    rowBlockSize = numRows - rowOffset;
                }

                pool.execute(new TransposeColumnToRowProcess(matrix, result, numColumns, rowBlockSize, rowOffset));
            }

            try {
                pool.shutdown();
                if (!pool.awaitTermination(60, TimeUnit.SECONDS)) {
                    throw new IllegalStateException("SuperByteMatrixBuilder: getInstanceTranspose: processing threads timed out.");
                }
            } catch (Exception e) {
                e.printStackTrace();
            }

        } else if (matrix instanceof SuperByteMatrixTranspose) {
            result = getInstance(numRows, numColumns);
            int columnBlockSize = TRANSPOSE_BLOCK_SIZE;
            for (int columnOffset = 0; columnOffset < numColumns; columnOffset += TRANSPOSE_BLOCK_SIZE) {
                if (numColumns - columnOffset < TRANSPOSE_BLOCK_SIZE) {
                    columnBlockSize = numColumns - columnOffset;
                }

                pool.execute(new TransposeRowToColumnProcess(matrix, result, numRows, columnBlockSize, columnOffset));
            }

            try {
                pool.shutdown();
                if (!pool.awaitTermination(60, TimeUnit.SECONDS)) {
                    throw new IllegalStateException("SuperByteMatrixBuilder: getInstanceTranspose: processing threads timed out.");
                }
            } catch (Exception e) {
                e.printStackTrace();
            }

        } else {
            throw new IllegalArgumentException("SuperByteMatrixBuilder: getInstanceTranspose: Don't Know how to Transpose: " + matrix.getClass().getName());
        }

        return result;
    }

    private static class TransposeRowToColumnProcess implements Runnable {

        private final int myColumnBlockSize;
        private final int myColumnOffset;
        private final SuperByteMatrix mySourceMatrix;
        private final SuperByteMatrix myDestMatrix;
        private final int myNumRows;

        public TransposeRowToColumnProcess(SuperByteMatrix srcMatrix, SuperByteMatrix destMatrix, int numRows, int columnBlockSize, int columnOffset) {
            myColumnBlockSize = columnBlockSize;
            myColumnOffset = columnOffset;
            mySourceMatrix = srcMatrix;
            myDestMatrix = destMatrix;
            myNumRows = numRows;
        }

        @Override
        public void run() {
            byte[][] temp = new byte[TRANSPOSE_BLOCK_SIZE][myColumnBlockSize];

            int rowBlockSize = TRANSPOSE_BLOCK_SIZE;
            for (int rowOffset = 0; rowOffset < myNumRows; rowOffset += TRANSPOSE_BLOCK_SIZE) {
                if (myNumRows - rowOffset < TRANSPOSE_BLOCK_SIZE) {
                    rowBlockSize = myNumRows - rowOffset;
                }

                for (int c = 0; c < myColumnBlockSize; c++) {
                    for (int r = 0; r < rowBlockSize; r++) {
                        temp[r][c] = mySourceMatrix.get(r + rowOffset, c + myColumnOffset);
                    }
                }

                for (int r = 0; r < rowBlockSize; r++) {
                    for (int c = 0; c < myColumnBlockSize; c++) {
                        myDestMatrix.set(r + rowOffset, c + myColumnOffset, temp[r][c]);
                    }
                }
            }
        }
    }

    private static class TransposeColumnToRowProcess implements Runnable {

        private final int myRowBlockSize;
        private final int myRowOffset;
        private final SuperByteMatrix mySourceMatrix;
        private final SuperByteMatrix myDestMatrix;
        private final int myNumColumns;

        public TransposeColumnToRowProcess(SuperByteMatrix srcMatrix, SuperByteMatrix destMatrix, int numColumns, int rowBlockSize, int rowOffset) {
            myRowBlockSize = rowBlockSize;
            myRowOffset = rowOffset;
            mySourceMatrix = srcMatrix;
            myDestMatrix = destMatrix;
            myNumColumns = numColumns;
        }

        @Override
        public void run() {
            byte[][] temp = new byte[myRowBlockSize][TRANSPOSE_BLOCK_SIZE];

            int columnBlockSize = TRANSPOSE_BLOCK_SIZE;
            for (int columnOffset = 0; columnOffset < myNumColumns; columnOffset += TRANSPOSE_BLOCK_SIZE) {
                if (myNumColumns - columnOffset < TRANSPOSE_BLOCK_SIZE) {
                    columnBlockSize = myNumColumns - columnOffset;
                }

                for (int r = 0; r < myRowBlockSize; r++) {
                    for (int c = 0; c < columnBlockSize; c++) {
                        temp[r][c] = mySourceMatrix.get(r + myRowOffset, c + columnOffset);
                    }
                }

                for (int c = 0; c < columnBlockSize; c++) {
                    for (int r = 0; r < myRowBlockSize; r++) {
                        myDestMatrix.set(r + myRowOffset, c + columnOffset, temp[r][c]);
                    }
                }

            }
        }
    }
}
