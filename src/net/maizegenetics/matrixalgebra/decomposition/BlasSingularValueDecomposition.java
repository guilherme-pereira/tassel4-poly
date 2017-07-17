package net.maizegenetics.matrixalgebra.decomposition;

import java.util.Arrays;

import org.apache.log4j.Logger;

import net.maizegenetics.matrixalgebra.Matrix.BlasDoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;

public class BlasSingularValueDecomposition implements SingularValueDecomposition {
	private BlasDoubleMatrix bdS = null;
	private BlasDoubleMatrix bdU = null;
	private BlasDoubleMatrix bdVT = null;
	public boolean successful = false;
	private static final Logger myLogger = Logger.getLogger(BlasSingularValueDecomposition.class);
	private double tol = 1e-12;
	
	public BlasSingularValueDecomposition(BlasDoubleMatrix bdm, char jobz) {
		int nrows = bdm.numberOfRows();
		int ncols = bdm.numberOfColumns();
		int ns = Math.min(nrows, ncols);
		double[] A = bdm.getMatrixCopy();
		double[] S = new double[ns];
		int usize = nrows * nrows;
		int lda, ldu, ldvt;
		lda = nrows;
		
		double[] U;
		if (jobz == 'N' || (jobz == 'O'& nrows >= ncols)) {
			U = new double[]{0};
			ldu = 1;
		} else {
			U = new double[usize];
			ldu = nrows;
		}
		
		double[] VT;
		int vtsize = ncols * ncols;
		if (jobz == 'N' || (jobz == 'O'& nrows < ncols)) {
			VT = new double[]{0};
			ldvt = 1;
		} else {
			VT = new double[vtsize];
			ldvt = ncols;
		}
		
		int result = BlasDoubleMatrix.singularValueDecompositionDgesdd(jobz, nrows, ncols, A, lda, S, U, ldu, VT, ldvt);
		if (result == 0) {
			successful = true;
			bdS = BlasDoubleMatrix.getInstance(ns, 1, S, true);
			if (jobz == 'A') {
				bdU = BlasDoubleMatrix.getInstance(nrows, nrows, U, true);
				bdVT = BlasDoubleMatrix.getInstance(ncols, ncols, VT, true);
			} else if (jobz == 'S') {
				if (nrows < ncols) {
					bdU = BlasDoubleMatrix.getInstance(nrows, nrows, U, true);
					int[] selection = new int[nrows];
					for (int i = 0; i < nrows; i++) selection[i] = i;
					double[] subvt = BlasDoubleMatrix.getSelectionFromDoubleArray(VT, ncols, ncols, selection, null);
					bdVT = BlasDoubleMatrix.getInstance(nrows, ncols, subvt, true);
				} else if (nrows > ncols) {
					int[] selection = new int[ncols];
					for (int i = 0; i < ncols; i++) selection[i] = i;
					double[] subu = BlasDoubleMatrix.getSelectionFromDoubleArray(U, nrows, nrows, null, selection);
					bdU = BlasDoubleMatrix.getInstance(nrows, ncols, subu, true);
					bdVT = BlasDoubleMatrix.getInstance(ncols, ncols, VT, true);
				} else {
					bdU = BlasDoubleMatrix.getInstance(nrows, nrows, U, true);
					bdVT = BlasDoubleMatrix.getInstance(ncols, ncols, VT, true);
				}
			} else if (jobz == 'O') {
				if (nrows >= ncols) {
					bdVT = BlasDoubleMatrix.getInstance(ncols, ncols, VT, true);
					int[] selection = new int[ncols];
					for (int i = 0; i < ncols; i++) selection[i] = i;
					U = BlasDoubleMatrix.getSelectionFromDoubleArray(A, nrows, ncols, null, selection);
					bdU = BlasDoubleMatrix.getInstance(nrows, ncols, U, true);
				} else {
					bdU = BlasDoubleMatrix.getInstance(nrows, nrows, U, true);
					int[] selection = new int[nrows];
					for (int i = 0; i < nrows; i++) selection[i] = i;
					VT = BlasDoubleMatrix.getSelectionFromDoubleArray(A, nrows, ncols, selection, null);
					bdVT = BlasDoubleMatrix.getInstance(nrows, ncols, VT, true);
				}
			}
			
		} else {
			myLogger.error("BlasSVD failed with a return value of " + result);
		}
	}
	
	public BlasSingularValueDecomposition(BlasDoubleMatrix bdm) {
		this(bdm, 'O');
	}

	@Override
	public DoubleMatrix getU(boolean transpose) {
		if (transpose) return bdU.transpose();
		else return bdU;
	}

	@Override
	public DoubleMatrix getV(boolean transpose) {
		if (transpose) return bdVT;
		else return bdVT.transpose();
	}

	@Override
	public DoubleMatrix getS() {
		return bdS;
	}

	@Override
	public double[] getSingularValues() {
		return bdS.getMatrix();
	}

	@Override
	public int getRank() {
		int rank = 0;
		double[] sv = bdS.getMatrix();
		int n = sv.length;
		for (int i = 0; i < n; i++) if (sv[i] > tol) rank++;
		return rank;
	}

	public double getTolerance() { return tol; }
	
	public void setTolerance(double tolerance) { tol = tolerance; }
	
	public boolean wasSuccessful() { return successful; }
}
