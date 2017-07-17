package net.maizegenetics.jGLiM;

import org.apache.log4j.Logger;

import net.maizegenetics.jGLiM.dm.SweepFast;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;


public class SymmetricMatrixInverterDM {

    private static final Logger myLogger = Logger.getLogger(SymmetricMatrixInverterDM.class);
	private SweepFast originalSweep;
	boolean[] included;
	
	
	public SymmetricMatrixInverterDM(DoubleMatrix symmetricMatrix) {
		originalSweep = new SweepFast();
		originalSweep.setUTA(SweepFast.UTAFromDoubleMatrix(symmetricMatrix));
		originalSweep.fullSweepSetDmin();
	}
	
	/**
	 * Calculates a new inverse when rows/columns of the original matrix have been removed. 
	 * Each element of exclude corresponds to a row/column of the matrix. If exclude is true, 
	 * then that row/column in removed.
	 * @param exclude	a boolean array with number of elements equal to the dimension of the original matrix.
	 * @return the inverse of the original matrix as updated by removing excluded rows and columns
	 */
	public DoubleMatrix getInverse(boolean[] exclude) {
		int n = exclude.length;
		myLogger.debug("exclude.length = " + n + ", dim(A) = " + originalSweep.getDimensionOfA());
		DoubleMatrix inverse;
		if (n == 0 || exclude == null) {
			inverse = originalSweep.getA();
			//zero out any singular rows and columns
			for (int i = 0; i < n; i++) {
				if (originalSweep.isSingular(i)) {
					for (int j = 0; j < n; j++) {
						inverse.set(i, j, 0);
						inverse.set(j, i, 0);
					}
				}
			}
			return inverse;
		} else {
			SweepFast sf = originalSweep.copy();
			for (int i = 0; i < n; i++) if (exclude[i]) sf.revg2sweep(i);
			sf.sweepSingularColumns();
			inverse = sf.getSubsetOfA(exclude);
			//zero out any singular rows and columns
			int count = 0;
			for (int i = 0; i < n; i++) {
				if (!exclude[i]) {
					if (sf.isSingular(i)) {
						for (int j = 0; j < n; j++) {
							inverse.set(count, j, 0);
							inverse.set(j, count, 0);
						}
					}
					count++;
				}
			}
		}

		return inverse;
	}
		
}
