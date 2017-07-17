package net.maizegenetics.matrixalgebra.decomposition;

import org.ejml.data.DenseMatrix64F;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;

public class ColtEigenvalueDecomposition implements EigenvalueDecomposition {
	cern.colt.matrix.linalg.EigenvalueDecomposition myDecomposition;
	
	public ColtEigenvalueDecomposition(DoubleMatrix dm) {
		int nrows = dm.numberOfRows();
		int ncols = dm.numberOfColumns();
		DoubleMatrix2D matrix = DoubleFactory2D.dense.make(nrows, ncols);
		for (int r = 0; r < nrows; r++) {
			for (int c = 0; c < ncols; c++) {
				matrix.setQuick(r, c, dm.get(r,c));
			}
		}
		myDecomposition = new cern.colt.matrix.linalg.EigenvalueDecomposition(matrix);
	}
	
	public ColtEigenvalueDecomposition(DenseMatrix64F dm) {
		int nrows = dm.numRows;
		DoubleMatrix2D matrix = DoubleFactory2D.dense.make(dm.data, nrows);
		myDecomposition = new cern.colt.matrix.linalg.EigenvalueDecomposition(matrix);
	}

	public ColtEigenvalueDecomposition(DoubleMatrix2D matrix) {
		myDecomposition = new cern.colt.matrix.linalg.EigenvalueDecomposition(matrix);
	}
	
	@Override
	public double[] getEigenvalues() {
		return myDecomposition.getRealEigenvalues().toArray();
	}

	@Override
	public double getEigenvalue(int i) {
		return myDecomposition.getRealEigenvalues().get(i);
	}

	@Override
	public DoubleMatrix getEigenvectors() {
		DoubleMatrix2D V = myDecomposition.getV();
		int nrows = V.rows();
		int ncols = V.columns();
		DoubleMatrix dm = DoubleMatrixFactory.DEFAULT.make(nrows, ncols);
		
		for (int r = 0; r < nrows; r++) {
			for (int c = 0; c < ncols; c++) {
				dm.set(r, c, V.getQuick(r,c));
			}
		}
		
		return dm;
	}

	@Override
	public DoubleMatrix getEigenvalueMatrix() {
		double[] eigenvalues = getEigenvalues();
		return DoubleMatrixFactory.DEFAULT.diagonal(eigenvalues);
	}

}
