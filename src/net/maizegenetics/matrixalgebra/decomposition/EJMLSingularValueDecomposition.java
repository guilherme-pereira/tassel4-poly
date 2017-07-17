package net.maizegenetics.matrixalgebra.decomposition;

import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.ops.SingularOps;

import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;
import net.maizegenetics.matrixalgebra.Matrix.EJMLDoubleMatrix;

public class EJMLSingularValueDecomposition implements
		SingularValueDecomposition {
	
	org.ejml.factory.SingularValueDecomposition<DenseMatrix64F> myDecomposition;
	boolean successful;
	int rank = -1;

	public EJMLSingularValueDecomposition(DenseMatrix64F matrix) {
		myDecomposition = DecompositionFactory.svd(matrix.numRows, matrix.numCols, true, true, true);
		successful = myDecomposition.decompose(matrix);
	}
	
	@Override
	public DoubleMatrix getS() {
		double[] singularValues = myDecomposition.getSingularValues();
		return DoubleMatrixFactory.DEFAULT.diagonal(singularValues);
	}

	@Override
	public double[] getSingularValues() {
		return myDecomposition.getSingularValues();
	}

	@Override
	public DoubleMatrix getU(boolean transpose) {
		return new EJMLDoubleMatrix(myDecomposition.getU(null, transpose));
	}

	@Override
	public DoubleMatrix getV(boolean transpose) {
		return new EJMLDoubleMatrix(myDecomposition.getV(null, transpose));
	}
	
	boolean wasSuccessful() {
		return successful;
	}

	@Override
	public int getRank() {
		rank = SingularOps.rank(myDecomposition, 1e-12);
		return rank;
	}

}
