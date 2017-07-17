package net.maizegenetics.matrixalgebra.decomposition;

import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;

/**
 * @author Peter Bradbury
 * created 7/22/2010
 */
public interface SingularValueDecomposition {
	/**
	 * for the decomposition of A, A = USV'
	 * @return U (orthogonal)
	 */
	DoubleMatrix getU(boolean transpose);
	
	/**
	 * for the decomposition of A, A = USV'
	 * @return V (orthogonal)
	 */
	DoubleMatrix getV(boolean transpose);
	
	/**
	 * for the decomposition of A, A = USV'
	 * @return S, the diagonal matrix of singular values
	 */
	DoubleMatrix getS();
	
	/**
	 * for the decomposition of A, A = USV'
	 * @return the singular values equal to the diagonal of S
	 */
	double[] getSingularValues();
	
	/**
	 * @return the rank of the matrix that was decomposed
	 */
	int getRank();
}
