package net.maizegenetics.matrixalgebra.decomposition;

import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;

public interface EigenvalueDecomposition {
	double[] getEigenvalues();
	double getEigenvalue(int i);
	DoubleMatrix getEigenvectors();
	DoubleMatrix getEigenvalueMatrix();
}
