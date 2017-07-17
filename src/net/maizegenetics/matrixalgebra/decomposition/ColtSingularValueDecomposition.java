package net.maizegenetics.matrixalgebra.decomposition;

import net.maizegenetics.matrixalgebra.Matrix.ColtDoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;

public class ColtSingularValueDecomposition implements
		SingularValueDecomposition {

	cern.colt.matrix.linalg.SingularValueDecomposition myDecomposition;
	
	public ColtSingularValueDecomposition(DoubleMatrix dm) {
		int rows = dm.numberOfRows();
		int cols = dm.numberOfColumns();
		DoubleMatrix2D matrix = DoubleFactory2D.dense.make(rows, cols);
		for (int r = 0; r < rows; r++) {
			for (int c = 0; c < cols; c++) {
				matrix.setQuick(r, c, dm.get(r, c));
			}
		}
		myDecomposition = new cern.colt.matrix.linalg.SingularValueDecomposition(matrix);
		
	}
	
	public ColtSingularValueDecomposition(DoubleMatrix2D matrix) {
		myDecomposition = new cern.colt.matrix.linalg.SingularValueDecomposition(matrix);
	}
	
	@Override
	public DoubleMatrix getU(boolean transpose) {
		if (transpose) {
			return new ColtDoubleMatrix(myDecomposition.getU().viewDice().copy());
		} else {
			return new ColtDoubleMatrix(myDecomposition.getU());
		}
	}

	@Override
	public DoubleMatrix getV(boolean transpose) {
		if (transpose) {
			return new ColtDoubleMatrix(myDecomposition.getV().viewDice().copy());
		} else {
			return new ColtDoubleMatrix(myDecomposition.getV());
		}
	}

	@Override
	public DoubleMatrix getS() {
		return new ColtDoubleMatrix(myDecomposition.getS());
	}

	@Override
	public double[] getSingularValues() {
		return myDecomposition.getSingularValues();
	}

	@Override
	public int getRank() {
		return myDecomposition.rank();
	}

}
