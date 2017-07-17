package net.maizegenetics.matrixalgebra.Matrix;

import net.maizegenetics.matrixalgebra.decomposition.EigenvalueDecomposition;
import net.maizegenetics.matrixalgebra.decomposition.QRDecomposition;
import net.maizegenetics.matrixalgebra.decomposition.SingularValueDecomposition;

public interface DoubleMatrix {
	/**
	 * The row and column coordinates are not checked to make sure they are in the matrix
	 * @param row the zero-based row index
	 * @param col the zero-based column index
	 * @return the value at row, col
	 */
	double get(int row, int col);
	
	/**
	 * The row and column coordinates are checked to make sure they fall in the matrix. If not an error is generated.
	 * @param row the zero-based row index
	 * @param col the zero-based column index
	 * @return the value at row, col
	 */
	double getChecked(int row, int col);
	
	/**
	 * Sets the matrix value at row, col. The coordinates are not checked to make sure they fall in the matrix.
	 * @param row the zero-based row index
	 * @param col the zero-based column index
	 * @param value the value to be set at row, col
	 */
	void set(int row, int col, double value);
	
	/**
	 * Sets the matrix value at row, col. The coordinates are checked to make sure they fall in the matrix.
	 * @param row the zero-based row index
	 * @param col the zero-based column index
	 * @param value the value to be set at row, col
	 */
	void setChecked(int row, int col, double value);
	
	/**
	 * @return the transpose of this matrix
	 */
	DoubleMatrix transpose();
	
	/**
	 * Multiply this matrix times another.
	 * @param dm	a double matrix
	 * @param transpose	if true, this matrix will be transposed before multiplying
	 * @param transposedm	if true, dm will be transposed before multiplying
	 */
	DoubleMatrix mult(DoubleMatrix dm, boolean transpose, boolean transposedm);
	
	/**
	 * Using this function for combining multiplication and addition allows the implementing library to optimize the operations.
	 * Performs: alpha * XA + beta*B, where X is this matrix
	 * @param A	the matrix to be multiplied
	 * @param B	the matrix to be added, can be null
	 * @param alpha	scalar multiplier for A
	 * @param beta	scalar multiplier for B
	 * @param transpose	if true, X is transposed 
	 * @param transposeA if true, A is transposed
	 * @param transposeB if true, B is transposed
	 * @return the matrix resulting from the specified operations
	 */
	DoubleMatrix multadd(DoubleMatrix A, DoubleMatrix B, double alpha, double beta, boolean transpose, boolean transposeA);
	
	/**
	 * @param dm	a double matrix
	 * @return this matrix times dm
	 */
	DoubleMatrix mult(DoubleMatrix dm);
	
	/**
	 * @return X'X, where X is this matrix
	 */
	DoubleMatrix crossproduct();
	
	/**
	 * @param dm	a double matrix
	 * @return	X'dm, where X is this matrix
	 */
	DoubleMatrix crossproduct(DoubleMatrix dm);
	
	/**
	 * @return XX', where X is this matrix
	 */
	DoubleMatrix tcrossproduct();
	
	/**
	 * @param dm a double matrix
	 * @return Xdm', where X is this matrix
	 */
	DoubleMatrix tcrossproduct(DoubleMatrix dm);
	
	/**
	 * @param dm a DoubleMatrix
	 * @param rows true if rows are to concatenated, false if columns are to be concatenated
	 * @return returns X with dm appended, where X is this matrix
	 */
	DoubleMatrix concatenate(DoubleMatrix dm, boolean rows);
	
	/**
	 * This returns the inverse or a square matrix, without modifying the original matrix.
	 * @return the inverse of a square matrix if it is non-singular, null otherwise
	 */
	DoubleMatrix inverse();
	
	/**
	 * This returns the inverse or a square matrix, replacing the original with the inverse.
	 * @return the inverse of a square matrix if it is non-singular, null otherwise
	 */
	void invert();
	
	/**
	 * @return the generalized inverse of a square matrix
	 */
	DoubleMatrix generalizedInverse();
	
	/**
	 * inverts the matrix and returns the rank as the first element in rank[]. The original matrix is not modified.
	 * @return a generalized inverse of this matrix
	 */
	DoubleMatrix generalizedInverseWithRank(int[] rank);
	
	/**
	 * @param Y a DoubleMatrix 
	 * @return the least squares solutions for B, where XB = Y and X is this matrix
	 */
	DoubleMatrix solve(DoubleMatrix Y);
	
	/**
	 * @return the number of rows in this matrix
	 */
	int numberOfRows();
	
	/**
	 * @return the number of columns in this matrix
	 */
	int numberOfColumns();
	
	/**
	 * @param i a row index
	 * @return the ith row of this matrix as a column vector
	 */
	DoubleMatrix row(int i);
	
	/**
	 * @param j a column index
	 * @return the jth column of this matrix as a column vector
	 */
	DoubleMatrix column(int j);
	
	/**
	 * @return an array of three DoubleMatrix's. Where X is this matrix, the first is X'X, 
	 * the second is the inverse of X'X, and the third is I - XGX'.
	 */
	DoubleMatrix[] getXtXGM();
	
	/**
	 * @return a copy of this matrix
	 */
	DoubleMatrix copy();
	
	/**
	 * @return an Eigenvalue Decomposition of this matrix
	 */
	EigenvalueDecomposition getEigenvalueDecomposition();
	
	/**
	 * @return a Singular Value Decomposition of this matrix
	 */
	SingularValueDecomposition getSingularValueDecomposition();
	
	/**
	 * @return a QR Decomposition of this matrix
	 */
	QRDecomposition getQRDecomposition();
	
	/**
	 * @param dm	a DoubleMatrix
	 * @return a new DoubleMatrix created by subtracting dm from this matrix
	 */
	DoubleMatrix minus(DoubleMatrix dm);
	
	/**
	 * This function subtracts dm, modifying the original matrix
	 * @param dm	a DoubleMatrix
	 */
	void minusEquals(DoubleMatrix dm);
	
	/**
	 * @param dm	a DoubleMatrix
	 * @return a new DoubleMatrix created by adding dm to this matrix
	 */
	DoubleMatrix plus(DoubleMatrix dm);
	
	/**
	 * This function adds dm, modifying the original matrix
	 * @param dm	a DoubleMatrix
	 */
	void plusEquals(DoubleMatrix dm);
	
	/**
	 * Adds a scalar to this matrix and returns a new matrix. The original is not modified.
	 * @param s	a scalar
	 * @return	the sum of this matrix and a scalar s
	 */
	DoubleMatrix scalarAdd(double s);
	
	/**
	 * Adds a scalar s to this matrix, replacing the original matrix with the result.
	 * @param s	a scalar
	 */
	void scalarAddEquals(double s);
	
	/**
	 * Multiplies this matrix times a scalar and returns a new matrix. The original is not modified.
	 * @param s	a scalar
	 * @return	the product of this matrix and a scalar s
	 */
	DoubleMatrix scalarMult(double s);
	
	/**
	 * Multiplies this matrix times a scalar s, replacing the original matrix with the result.
	 * @param s	a scalar
	 */
	void scalarMultEquals(double s);
	
	/**
	 * Creates a new matrix consisting or the rows and columns of this matrix in the order specified. 
	 * If rows or columns is null then all rows or columns, respectively, will be included
	 * @param rows	the rows to be included in the new matrix
	 * @param	columns the columns to be included in the new matrix
	 * @return	a new matrix consisting of the specified rows and columns 
	 */
	DoubleMatrix getSelection(int[] rows, int[] columns);
	
	/**
	 * @param row
	 * @return	the sum of elements in this row
	 */
	double rowSum(int row);
	
	/**
	 * @param column
	 * @return	the sum of elements in this column
	 */
	double columnSum(int column);
	
	/**
	 * @return the column rank of this matrix
	 */
	int columnRank();
	
}
