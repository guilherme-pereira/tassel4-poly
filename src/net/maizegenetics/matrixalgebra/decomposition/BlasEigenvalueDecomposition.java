package net.maizegenetics.matrixalgebra.decomposition;

import java.util.Arrays;

import net.maizegenetics.matrixalgebra.Matrix.BlasDoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;

public class BlasEigenvalueDecomposition implements EigenvalueDecomposition {
	
	protected double[] eigenvalues;
	protected double[] eigenvectors;
	protected BlasDoubleMatrix bdm;
	protected int info;
	
	// use public static native int eigenValueSymmetricDecomposition(int order, double[] A, double[] eigval) from BlasDoubleMatrix
	public BlasEigenvalueDecomposition(DoubleMatrix dm) {
		bdm = (BlasDoubleMatrix) dm;
		init();
	}
	
	protected void init() {
		double[] originalMatrix = Arrays.copyOf(bdm.getMatrix(), bdm.getSize());
		int order = bdm.numberOfRows();
		eigenvalues = new double[order];
		eigenvectors = new double[order * order];
		info = BlasDoubleMatrix.eigenValueSymmetricDecomposition(order, originalMatrix, eigenvalues, eigenvectors);
	}
	
	@Override
	public double[] getEigenvalues() {
		return eigenvalues;
	}

	@Override
	public double getEigenvalue(int i) {
		return eigenvalues[i];
	}

	@Override
	public DoubleMatrix getEigenvectors() {
		int nrows = bdm.numberOfRows();
		return BlasDoubleMatrix.getInstance(nrows, nrows, eigenvectors, true);
	}  

	@Override
	public DoubleMatrix getEigenvalueMatrix() {
		return BlasDoubleMatrix.getDiagonalMatrix(eigenvalues);
	}

}
