package net.maizegenetics.matrixalgebra.Matrix;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.math.Functions;
import net.maizegenetics.matrixalgebra.decomposition.ColtEigenvalueDecomposition;
import net.maizegenetics.matrixalgebra.decomposition.ColtSingularValueDecomposition;
import net.maizegenetics.matrixalgebra.decomposition.EigenvalueDecomposition;
import net.maizegenetics.matrixalgebra.decomposition.QRDecomposition;
import net.maizegenetics.matrixalgebra.decomposition.SingularValueDecomposition;

public class ColtDoubleMatrix implements DoubleMatrix {
	DoubleMatrix2D myMatrix;
	
	public ColtDoubleMatrix(DoubleMatrix2D aMatrix) {
		myMatrix = aMatrix;
	}
	
	public ColtDoubleMatrix(int row, int col) {
		myMatrix = DoubleFactory2D.dense.make(row, col);
	}
	
	public ColtDoubleMatrix(int row, int col, double[] values) {
		myMatrix = DoubleFactory2D.dense.make(values, row);
	}
	
	public ColtDoubleMatrix(int row, int col, double value) {
		myMatrix = DoubleFactory2D.dense.make(row, col, value);
	}
	
	public ColtDoubleMatrix(double[][] values) {
		myMatrix = DoubleFactory2D.dense.make(values);
	}
	
	public ColtDoubleMatrix(int size) {
		myMatrix = DoubleFactory2D.dense.identity(size);
	}
	
	public ColtDoubleMatrix(double[] diagonal) {
		myMatrix = DoubleFactory2D.dense.diagonal(DoubleFactory1D.dense.make(diagonal));
	}
	
	
	@Override
	public double get(int row, int col) {
		return myMatrix.getQuick(row, col);
	}

	@Override
	public double getChecked(int row, int col) {
		return myMatrix.get(row, col);
	}

	@Override
	public void set(int row, int col, double value) {
		myMatrix.setQuick(row, col, value);
	}

	@Override
	public void setChecked(int row, int col, double value) {
		myMatrix.set(row, col, value);
	}

	@Override
	public DoubleMatrix transpose() {
		return new ColtDoubleMatrix(myMatrix.viewDice().copy());
	}

	@Override
	public DoubleMatrix mult(DoubleMatrix dm, boolean transpose,
			boolean transposedm) {
		DoubleMatrix2D B = ((ColtDoubleMatrix) dm).myMatrix;
		return new ColtDoubleMatrix(myMatrix.zMult(B, null, 1.0, 0.0, transpose, transpose));
	}

	@Override
	public DoubleMatrix multadd(DoubleMatrix A, DoubleMatrix B, double alpha,
			double beta, boolean transpose, boolean transposeA) {
		DoubleMatrix2D C = ((ColtDoubleMatrix) A).myMatrix;
		DoubleMatrix2D D;
		if (B == null) D = null;
		else D = ((ColtDoubleMatrix) B).myMatrix;
		DoubleMatrix2D result = myMatrix.zMult(C, D, alpha, beta, transpose, transposeA);
		return new ColtDoubleMatrix(result);
	}

	@Override
	public DoubleMatrix mult(DoubleMatrix dm) {
		DoubleMatrix2D B = ((ColtDoubleMatrix) dm).myMatrix;
		return new ColtDoubleMatrix(myMatrix.zMult(B, null));
	}

	@Override
	public DoubleMatrix crossproduct() {
		return new ColtDoubleMatrix(myMatrix.viewDice().zMult(myMatrix, null));
	}

	@Override
	public DoubleMatrix crossproduct(DoubleMatrix dm) {
		DoubleMatrix2D B = ((ColtDoubleMatrix) dm).myMatrix;
		return new ColtDoubleMatrix(myMatrix.viewDice().zMult(B, null));
	}

	@Override
	public DoubleMatrix tcrossproduct() {
		return new ColtDoubleMatrix(myMatrix.zMult(myMatrix.viewDice(), null));
	}

	@Override
	public DoubleMatrix tcrossproduct(DoubleMatrix dm) {
		DoubleMatrix2D B = ((ColtDoubleMatrix) dm).myMatrix;
		return new ColtDoubleMatrix(myMatrix.zMult(B.viewDice(), null));
	}

	@Override
	public DoubleMatrix concatenate(DoubleMatrix dm, boolean rows) {
		DoubleMatrix2D B = ((ColtDoubleMatrix) dm).myMatrix;
		if (rows) {
			return new ColtDoubleMatrix(DoubleFactory2D.dense.appendRows(myMatrix, B));
		} else {
			return new ColtDoubleMatrix(DoubleFactory2D.dense.appendColumns(myMatrix, B));
		}
	}

	@Override
	public DoubleMatrix inverse() {
		return new ColtDoubleMatrix(Algebra.DEFAULT.inverse(myMatrix));
	}

	@Override
	public void invert() {
		myMatrix = Algebra.DEFAULT.inverse(myMatrix);
	}

	@Override
	public DoubleMatrix generalizedInverse() {
        Algebra A = new Algebra();
        boolean transposeMatrix = false;
        cern.colt.matrix.linalg.SingularValueDecomposition svd;
        if (myMatrix.rows() < myMatrix.columns()) {
            transposeMatrix = true;
            svd = new cern.colt.matrix.linalg.SingularValueDecomposition(myMatrix.viewDice());
        } else {
            svd = new cern.colt.matrix.linalg.SingularValueDecomposition(myMatrix);
        }
        
        DoubleMatrix2D invS = svd.getS();

        //calculate the inverse of S, a diagonal matrix with rank(aMatrix) non-zero elements
        int size = invS.rows();
        int r = 0;
        for (int i = 0; i < size; i++) {
            if (Math.abs(invS.get(i, i)) > 1E-10) {
                invS.set(i, i, 1 / invS.get(i, i));
                r++;
            }
            else
                invS.set(i, i, 0);
        }
        DoubleMatrix2D minv = A.mult(A.mult(svd.getV(), invS), A.transpose(svd.getU()));
        if (transposeMatrix) return new ColtDoubleMatrix(minv.viewDice().copy());
        return new ColtDoubleMatrix(minv);
	}

	@Override
	public DoubleMatrix solve(DoubleMatrix Y) {
		return new ColtDoubleMatrix(Algebra.DEFAULT.solve(myMatrix, ((ColtDoubleMatrix) Y).myMatrix));
	}

	@Override
	public int numberOfRows() {
		return myMatrix.rows();
	}

	@Override
	public int numberOfColumns() {
		return myMatrix.columns();
	}

	@Override
	public DoubleMatrix row(int i) {
		int[] rowIndex = new int[]{i};
		return new ColtDoubleMatrix(myMatrix.viewSelection(rowIndex,null).copy());
	}

	@Override
	public DoubleMatrix column(int j) {
		int[] colIndex = new int[]{j};
		return new ColtDoubleMatrix(myMatrix.viewSelection(null, colIndex).copy());
	}

	@Override
	public DoubleMatrix[] getXtXGM() {
		DoubleMatrix[] result = new DoubleMatrix[3];
		DoubleMatrix2D xtx = myMatrix.viewDice().zMult(myMatrix, null);
		DoubleMatrix2D g = Algebra.DEFAULT.inverse(xtx);
		DoubleMatrix2D m = DoubleFactory2D.dense.identity(myMatrix.rows());
		DoubleMatrix2D xgxt = myMatrix.zMult(g.zMult(myMatrix.viewDice(), null), null);
		m.assign(xgxt, Functions.minus);
		result[0] = new ColtDoubleMatrix(xtx);
		result[1] = new ColtDoubleMatrix(g);
		result[2] = new ColtDoubleMatrix(m);
		return result;
	}

	@Override
	public DoubleMatrix copy() {
		return new ColtDoubleMatrix(myMatrix.copy());
	}

	@Override
	public EigenvalueDecomposition getEigenvalueDecomposition() {
		return new ColtEigenvalueDecomposition(myMatrix);
	}

	@Override
	public SingularValueDecomposition getSingularValueDecomposition() {
		return new ColtSingularValueDecomposition(myMatrix);
	}

	@Override
	public QRDecomposition getQRDecomposition() {
		//TODO implement QRDecomposition
		return null;
	}

	@Override
	public DoubleMatrix minus(DoubleMatrix dm) {
		DoubleMatrix2D B = ((ColtDoubleMatrix)dm).myMatrix;
		return new ColtDoubleMatrix(myMatrix.copy().assign(B, Functions.minus));
	}

	@Override
	public void minusEquals(DoubleMatrix dm) {
		DoubleMatrix2D B = ((ColtDoubleMatrix)dm).myMatrix;
		myMatrix.assign(B, Functions.minus);
	}

	@Override
	public DoubleMatrix plus(DoubleMatrix dm) {
		DoubleMatrix2D B = ((ColtDoubleMatrix)dm).myMatrix;
		return new ColtDoubleMatrix(myMatrix.copy().assign(B, Functions.plus));
	}

	@Override
	public void plusEquals(DoubleMatrix dm) {
		DoubleMatrix2D B = ((ColtDoubleMatrix)dm).myMatrix;
		myMatrix.assign(B, Functions.plus);
	}

	@Override
	public DoubleMatrix scalarAdd(double s) {
		DoubleMatrix2D S = DoubleFactory2D.dense.make(myMatrix.rows(), myMatrix.columns(), s);
		return new ColtDoubleMatrix(S.assign(myMatrix, Functions.plus));
	}

	@Override
	public void scalarAddEquals(double s) {
		myMatrix.assign(Functions.plus(s));
	}

	@Override
	public DoubleMatrix scalarMult(double s) {
		return new ColtDoubleMatrix(myMatrix.copy().assign(Functions.mult(s)));
	}

	@Override
	public void scalarMultEquals(double s) {
		myMatrix.assign(Functions.mult(s));
	}

	@Override
	public DoubleMatrix getSelection(int[] rows, int[] columns) {
		return new ColtDoubleMatrix(myMatrix.viewSelection(rows,columns).copy());
	}

	@Override
	public double rowSum(int row) {
		return myMatrix.viewRow(row).aggregate(Functions.plus, Functions.identity);
	}

	@Override
	public double columnSum(int column) {
		return myMatrix.viewColumn(column).aggregate(Functions.plus, Functions.identity);
	}

	@Override
	public int columnRank() {
		return Algebra.DEFAULT.rank(myMatrix);
	}

	@Override
	public DoubleMatrix generalizedInverseWithRank(int[] rank) {
		double tol = 1e-10;
        Algebra A = new Algebra();
        boolean transposeMatrix = false;
        if (myMatrix.rows() < myMatrix.columns()) {
            transposeMatrix = true;
        }
        
        cern.colt.matrix.linalg.SingularValueDecomposition svd;
        if (transposeMatrix) {
        	svd = new cern.colt.matrix.linalg.SingularValueDecomposition(myMatrix.viewDice());
        } else {
        	svd = new cern.colt.matrix.linalg.SingularValueDecomposition(myMatrix);
        }
        DoubleMatrix2D invS = svd.getS();

        //calculate the inverse of S, a diagonal matrix with rank(aMatrix) non-zero elements
        int size = invS.rows();
        int r = 0;
        for (int i = 0; i < size; i++) {
            if (Math.abs(invS.get(i, i)) > tol) {
                invS.set(i, i, 1 / invS.get(i, i));
                r++;
            }
            else
                invS.set(i, i, 0);
        }
        rank[0] = r;
        DoubleMatrix2D minv = A.mult(A.mult(svd.getV(), invS), A.transpose(svd.getU()));
        if (transposeMatrix) minv = A.transpose(minv);
        return new ColtDoubleMatrix(minv);
	}

	@Override
	public String toString() {
		return myMatrix.toString();
	}

}
