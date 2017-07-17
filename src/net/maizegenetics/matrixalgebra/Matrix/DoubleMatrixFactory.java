package net.maizegenetics.matrixalgebra.Matrix;

import org.apache.log4j.Logger;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import org.ejml.ops.SpecializedOps;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;

public class DoubleMatrixFactory {
	private static Logger myLogger = Logger.getLogger(DoubleMatrixFactory.class);
	public enum FactoryType {ejml, jblas, colt, blas};
	private FactoryType myType;
	public static DoubleMatrixFactory DEFAULT;
	
	static{
		try {
			System.loadLibrary("TasselBlas");
			DEFAULT = new DoubleMatrixFactory(FactoryType.blas);
			myLogger.info("Using BLAS/LAPACK for DoubleMatrix operations");
		} catch (UnsatisfiedLinkError err) {
			//err.printStackTrace();
			DEFAULT = new DoubleMatrixFactory(FactoryType.ejml);
			myLogger.info("TasselBlas library for BLAS/LAPACK not found. Using EJML for DoubleMatrix operations.");
		} catch (SecurityException se) {
			//se.printStackTrace();
			DEFAULT = new DoubleMatrixFactory(FactoryType.ejml);
			myLogger.info("No permission to load blasDoubleMatrix library. Using EJML for DoubleMatrix operations.");
		}
	}
	
	public DoubleMatrixFactory(FactoryType type) {
		myType = type;
	}
	
	public static void setDefault(FactoryType type) {
		DoubleMatrixFactory.DEFAULT = new DoubleMatrixFactory(type);
	}
	
	public FactoryType getType() { return myType; }
	
	/**
	 * @param row	the number of rows in the matrix
	 * @param col	the number of columns in the matrix
	 * @return	a matrix with all elements equal to 0
	 */
	public DoubleMatrix make(int row, int col) {
		if (myType == FactoryType.ejml) return new EJMLDoubleMatrix(row, col);
		if (myType == FactoryType.colt) return new ColtDoubleMatrix(row, col);
		if (myType == FactoryType.blas) return new BlasDoubleMatrix(row, col);
		if (myType == FactoryType.jblas) return null;
		return null;
	}
	
	public DoubleMatrix make(int row, int col, double[] values) {
		return this.make(row, col, values, false);
	}

	public DoubleMatrix make(int row, int col, double[] values, boolean columnMajor) {
		if (columnMajor) {
			if (myType == FactoryType.ejml){
				DoubleMatrix temp = new EJMLDoubleMatrix(col, row, values);
				return temp.transpose();
			} 
			else if (myType == FactoryType.colt) return new ColtDoubleMatrix(row, col, values);
			else if (myType == FactoryType.blas) return BlasDoubleMatrix.getInstance(row, col, values, true);
		} else {
			if (myType == FactoryType.ejml) return new EJMLDoubleMatrix(row, col, values);
			if (myType == FactoryType.colt) {
				DoubleMatrix temp =  new ColtDoubleMatrix(col, row, values);
				return temp.transpose();
			}
			else if (myType == FactoryType.blas) return BlasDoubleMatrix.getInstance(row, col, values, false);
		}
		if (myType == FactoryType.jblas) return null;
		return null;
	}

	public DoubleMatrix make(double[][] values) {
		if (myType == FactoryType.ejml) return new EJMLDoubleMatrix(values);
		if (myType == FactoryType.colt) return new ColtDoubleMatrix(values);
		if (myType == FactoryType.blas) return new BlasDoubleMatrix(values);
		if (myType == FactoryType.jblas) return null;
		return null;
	}
	
	public DoubleMatrix make(int row, int col, double val) {
		if (myType == FactoryType.ejml) return new EJMLDoubleMatrix(row, col, val);
		if (myType == FactoryType.colt) return new ColtDoubleMatrix(row, col, val);
		if (myType == FactoryType.blas) return BlasDoubleMatrix.getInstance(row, col, val);
		if (myType == FactoryType.jblas) return null;
		return null;
	}
	
	public DoubleMatrix identity(int n) {
		if (myType == FactoryType.ejml) return new EJMLDoubleMatrix(n);
		if (myType == FactoryType.colt) return new ColtDoubleMatrix(n);
		if (myType == FactoryType.blas) return BlasDoubleMatrix.getIdentityMatrix(n);
		if (myType == FactoryType.jblas) return null;
		return null;
	}
	
	public DoubleMatrix diagonal(double[] diag) {
		if (myType == FactoryType.ejml) return new EJMLDoubleMatrix(diag);
		if (myType == FactoryType.colt) return new ColtDoubleMatrix(diag);
		if (myType == FactoryType.blas) return BlasDoubleMatrix.getDiagonalMatrix(diag);
		if (myType == FactoryType.jblas) return null;
		return null;
	}
	
	public DoubleMatrix compose(DoubleMatrix[][] components) {
		
		if (myType == FactoryType.ejml) {
			int totalRows = 0;
			int totalCols = 0;
			int nRows = components.length;
			int nCols = components[0].length;
			for (int i = 0; i < nRows; i++) totalRows += components[i][0].numberOfRows();
			for (int i = 0; i < nCols; i++) totalCols += components[0][i].numberOfColumns();
			
			DenseMatrix64F result = new DenseMatrix64F(totalRows, totalCols);
			int startRow = 0;
			for (int r = 0; r < nRows; r++) {
				int startCol = 0;
				for (int c = 0; c < nCols; c++) {
					DenseMatrix64F dm = ((EJMLDoubleMatrix) components[r][c]).myMatrix;
					CommonOps.insert(dm, result, startRow, startCol);
					startCol += dm.numCols;
				}
				startRow += components[r][0].numberOfRows();
			}
			return new EJMLDoubleMatrix(result);
		}
		
		if (myType == FactoryType.colt) {
			int nRows = components.length;
			int nCols = components[0].length;
			DoubleMatrix2D[][] coltComponents = new DoubleMatrix2D[nRows][nCols];
			for (int r = 0; r < nRows; r++) {
				for (int c = 0; c < nCols; c++) {
					coltComponents[r][c] = ((ColtDoubleMatrix) components[r][c]).myMatrix;
				}
			}
			return new ColtDoubleMatrix(DoubleFactory2D.dense.compose(coltComponents));
		}
		
		if (myType == FactoryType.blas) {
			return BlasDoubleMatrix.compose(components);
		}
		
		if (myType == FactoryType.jblas) return null;
		return null;
	}
}
