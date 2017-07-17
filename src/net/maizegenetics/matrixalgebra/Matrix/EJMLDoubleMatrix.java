package net.maizegenetics.matrixalgebra.Matrix;

import org.ejml.factory.DecompositionFactory;
import org.ejml.factory.LinearSolver;
import org.ejml.factory.LinearSolverFactory;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import org.ejml.ops.SpecializedOps;

import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.decomposition.ColtEigenvalueDecomposition;
import net.maizegenetics.matrixalgebra.decomposition.EJMLEigenvalueDecomposition;
import net.maizegenetics.matrixalgebra.decomposition.EJMLSingularValueDecomposition;
import net.maizegenetics.matrixalgebra.decomposition.EigenvalueDecomposition;
import net.maizegenetics.matrixalgebra.decomposition.QRDecomposition;
import net.maizegenetics.matrixalgebra.decomposition.SingularValueDecomposition;

public class EJMLDoubleMatrix implements DoubleMatrix {

	public final DenseMatrix64F myMatrix;
	
	public EJMLDoubleMatrix(DenseMatrix64F aMatrix) {
		myMatrix = aMatrix;
	}
	
	public EJMLDoubleMatrix(int row, int col) {
		myMatrix = new DenseMatrix64F(row, col);
	}
	
	public EJMLDoubleMatrix(int row, int col, double[] values) {
		myMatrix = new DenseMatrix64F(row, col, true, values);
	}

	public EJMLDoubleMatrix(int row, int col, double value) {
		myMatrix = new DenseMatrix64F(row, col);
		if (value != 0) {
			for (int r = 0; r < row; r++) {
				for (int c = 0; c < col; c++) {
					myMatrix.set(r,c,value);
				}
					
			}
		}
	}

	public EJMLDoubleMatrix(double[][] values) {
		myMatrix = new DenseMatrix64F(values);
	}
	
	public EJMLDoubleMatrix(int size) {
		myMatrix = CommonOps.identity(size);
	}
	
	public EJMLDoubleMatrix(double[] diagonal) {
		myMatrix = CommonOps.diag(diagonal);
	}
	
	@Override
	public DoubleMatrix column(int j) {
		int n = myMatrix.numRows;
		DenseMatrix64F result = new DenseMatrix64F(n,1);
		SpecializedOps.subvector(myMatrix, 0, j, n, false, 0, result);
		return new EJMLDoubleMatrix(result);
	}

	@Override
	public int columnRank() {
		if (myMatrix.numCols == 1) return 1;
		EJMLSingularValueDecomposition svd = new EJMLSingularValueDecomposition(myMatrix);
		return svd.getRank();
	}

	@Override
	public DoubleMatrix concatenate(DoubleMatrix dm, boolean rows) {
		DenseMatrix64F result = null;
		DenseMatrix64F otherMatrix = ((EJMLDoubleMatrix) dm).myMatrix;
		int myNumberOfRows = myMatrix.numRows;
		int myNumberOfCols = myMatrix.numCols;
		int dmNumberOfRows = otherMatrix.numRows;
		int dmNumberOfCols = otherMatrix.numCols;
		if (rows) {
			if (myNumberOfCols != dmNumberOfCols) {
				StringBuilder sb = new StringBuilder("Non-conformable matrices in concatenate rows: ");
				sb.append(myNumberOfRows).append(" x ").append(myNumberOfCols);
				sb.append(", ").append(dmNumberOfRows).append(" x ").append(dmNumberOfCols);
				throw new IllegalArgumentException(sb.toString());
			}
			int totalRows = myNumberOfRows + dmNumberOfRows;
			result = new DenseMatrix64F(totalRows, myNumberOfCols);
			CommonOps.insert(myMatrix, result,0,0);
			CommonOps.insert(otherMatrix,result, dmNumberOfRows,0);
		} else {
			if (myNumberOfRows != dmNumberOfRows) {
				StringBuilder sb = new StringBuilder("Non-conformable matrices in concatenate columns: ");
				sb.append(myNumberOfRows).append(" x ").append(myNumberOfCols);
				sb.append(", ").append(dmNumberOfRows).append(" x ").append(dmNumberOfCols);
				throw new IllegalArgumentException(sb.toString());
			}
			int totalCol = myNumberOfCols + dmNumberOfCols;
			result = new DenseMatrix64F(myNumberOfRows, totalCol);
			CommonOps.insert(myMatrix, result, 0, 0);
			CommonOps.insert(otherMatrix, result, 0, myNumberOfCols);
		}
		return new EJMLDoubleMatrix(result);
	}

	@Override
	public DoubleMatrix copy() {
		return new EJMLDoubleMatrix(myMatrix.copy());
	}

	@Override
	public DoubleMatrix crossproduct() {
		int n = myMatrix.numCols;
		DenseMatrix64F result = new DenseMatrix64F(n,n);
		CommonOps.multTransA(myMatrix, myMatrix, result);
		return new EJMLDoubleMatrix(result);
	}

	@Override
	public DoubleMatrix crossproduct(DoubleMatrix dm) {
		int nrow = myMatrix.numCols;
		int ncol = dm.numberOfColumns();
		DenseMatrix64F otherMatrix = ((EJMLDoubleMatrix) dm).myMatrix;
		DenseMatrix64F result = new DenseMatrix64F(nrow,ncol);
		CommonOps.multTransA(myMatrix, otherMatrix, result);
		return new EJMLDoubleMatrix(result);
	}

	@Override
	public DoubleMatrix generalizedInverse() {
		DenseMatrix64F result = new DenseMatrix64F(myMatrix.numRows, myMatrix.numCols);
		CommonOps.pinv(myMatrix, result);
		return new EJMLDoubleMatrix(result);
	}

	@Override
	public DoubleMatrix generalizedInverseWithRank(int[] rank) {
		org.ejml.factory.SingularValueDecomposition<DenseMatrix64F> myDecomposition = DecompositionFactory.svd(myMatrix.numRows, myMatrix.numCols, true, true, false);
		myDecomposition.decompose(myMatrix);
		double tol = 1e-10;
		rank[0] = 0;
		DenseMatrix64F W = myDecomposition.getW(null);
		int n = W.getNumRows();
		for (int i = 0; i < n; i++) {
			double val = W.get(i,i);
			if (val < tol) {
				val = 0;
			} else {
				val = 1/val;
				rank[0]++;
			}
			W.set(i, i, val);
		}
		
		DenseMatrix64F V = myDecomposition.getV(null, false);
		DenseMatrix64F UT  = myDecomposition.getU(null, true);
		
		int nrows = V.getNumRows();
		int ncols = W.getNumCols();
		
		DenseMatrix64F VW = new DenseMatrix64F(nrows, ncols);
		CommonOps.mult(V, W, VW);
		
		ncols = UT.getNumCols();
		
		DenseMatrix64F inv = new DenseMatrix64F(nrows, ncols);
		CommonOps.mult(VW, UT, inv);
		return new EJMLDoubleMatrix(inv);
	}

	@Override
	public double get(int row, int col) {
		return myMatrix.get(row, col);
	}

	@Override
	public double getChecked(int row, int col) {
		return myMatrix.get(row,col);
	}

	@Override
	public EigenvalueDecomposition getEigenvalueDecomposition() {
		EJMLEigenvalueDecomposition decomposition = new EJMLEigenvalueDecomposition(myMatrix);
		if (decomposition.wasSuccessful()) return decomposition;
		return new ColtEigenvalueDecomposition(myMatrix);
	}

	@Override
	public DoubleMatrix[] getXtXGM() {
		DoubleMatrix[] dmarray = new DoubleMatrix[3];
		int ncol = myMatrix.numCols;
		int nrow = myMatrix.numRows;
		DenseMatrix64F result1 = new DenseMatrix64F(ncol, ncol);
		CommonOps.multTransA(myMatrix, myMatrix, result1);
		dmarray[0] = new EJMLDoubleMatrix(result1);
		
		DenseMatrix64F inverse = new DenseMatrix64F(ncol, ncol);
		CommonOps.invert(result1, inverse);
		dmarray[1] = new EJMLDoubleMatrix(inverse);
		
		DenseMatrix64F result2 = new DenseMatrix64F(ncol, nrow);
		DenseMatrix64F ident = CommonOps.identity(nrow);
		CommonOps.multTransB(inverse, myMatrix, result2);
		CommonOps.multAdd(-1, myMatrix, result2, ident);
		dmarray[2] = new EJMLDoubleMatrix(ident);
		return dmarray;
	}

	@Override
	public QRDecomposition getQRDecomposition() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public SingularValueDecomposition getSingularValueDecomposition() {
		return new EJMLSingularValueDecomposition(myMatrix);
	}

	@Override
	public DoubleMatrix inverse() {
		DenseMatrix64F inverse = new DenseMatrix64F(myMatrix.numRows, myMatrix.numCols);
		CommonOps.invert(myMatrix, inverse);
		return new EJMLDoubleMatrix(inverse);
	}

	@Override
	public void invert() {
		CommonOps.invert(myMatrix);
	}

	@Override
	public DoubleMatrix minus(DoubleMatrix dm) {
		DenseMatrix64F result = new DenseMatrix64F(myMatrix.numRows, myMatrix.numCols);
		CommonOps.sub(myMatrix, ((EJMLDoubleMatrix) dm).myMatrix, result);
		return new EJMLDoubleMatrix(result);
	}

	@Override
	public void minusEquals(DoubleMatrix dm) {
		CommonOps.subEquals(myMatrix, ((EJMLDoubleMatrix) dm).myMatrix);
	}

	@Override
	public DoubleMatrix mult(DoubleMatrix dm, boolean transpose,
			boolean transposedm) {
		DenseMatrix64F b = ((EJMLDoubleMatrix) dm).myMatrix;
		DenseMatrix64F result;
		
		int nrow, ncol;
		if (transpose) {
			nrow = myMatrix.numCols;
			if (transposedm) {
				ncol = b.numRows;
				result = new DenseMatrix64F(nrow, ncol);
				CommonOps.multTransAB(myMatrix, b, result);
			} else {
				ncol = b.numCols;
				result = new DenseMatrix64F(nrow, ncol);
				CommonOps.multTransA(myMatrix, b, result);
			}
		} else {
			nrow = myMatrix.numRows;
			if (transposedm) {
				ncol = b.numRows;
				result = new DenseMatrix64F(nrow, ncol);
				CommonOps.multTransB(myMatrix, b, result);
			} else {
				ncol = b.numCols;
				result = new DenseMatrix64F(nrow, ncol);
				CommonOps.mult(myMatrix, b, result);
			}
		}
		
		return new EJMLDoubleMatrix(result);
	}

	@Override
	public DoubleMatrix mult(DoubleMatrix dm) {
		DenseMatrix64F b = ((EJMLDoubleMatrix) dm).myMatrix;
		int nrow = myMatrix.numRows;
		int ncol = b.numCols;
		DenseMatrix64F result = new DenseMatrix64F(nrow, ncol);
		CommonOps.mult(myMatrix, b, result);
		return new EJMLDoubleMatrix(result);
	}

	@Override
	public DoubleMatrix multadd(DoubleMatrix A, DoubleMatrix B, double alpha,
			double beta, boolean transpose, boolean transposeA) {
		
		DoubleMatrix result = mult(A, transpose, transposeA);
		if (alpha != 1) result.scalarMultEquals(alpha);
		if (B == null) return result;

		if (beta == 1) return result.plus(B);
		return result.plus(B.scalarMult(beta));
	}

	@Override
	public int numberOfColumns() {
		return myMatrix.numCols;
	}

	@Override
	public int numberOfRows() {
		return myMatrix.numRows;
	}

	@Override
	public DoubleMatrix plus(DoubleMatrix dm) {
		DenseMatrix64F result = new DenseMatrix64F(myMatrix.numRows, myMatrix.numCols);
		CommonOps.add(myMatrix, ((EJMLDoubleMatrix) dm).myMatrix, result);
		return new EJMLDoubleMatrix(result);
	}

	@Override
	public void plusEquals(DoubleMatrix dm) {
		CommonOps.addEquals(myMatrix, ((EJMLDoubleMatrix) dm).myMatrix);
	}

	@Override
	public DoubleMatrix row(int i) {
		int n = myMatrix.numCols;
		DenseMatrix64F result = new DenseMatrix64F(1,n);
		SpecializedOps.subvector(myMatrix, 0, i, n, true, 0, result);
		return new EJMLDoubleMatrix(result);
	}
	
	@Override
	public DoubleMatrix scalarAdd(double s) {
		DenseMatrix64F result = new DenseMatrix64F(myMatrix.numRows, myMatrix.numCols);
		CommonOps.add(myMatrix, s, result);
		return new EJMLDoubleMatrix(result);
	}

	@Override
	public void scalarAddEquals(double s) {
		CommonOps.add(myMatrix, s);
	}

	@Override
	public DoubleMatrix scalarMult(double s) {
		DenseMatrix64F result = new DenseMatrix64F(myMatrix.numRows, myMatrix.numCols);
		CommonOps.scale(s, myMatrix, result);
		return new EJMLDoubleMatrix(result);
	}

	@Override
	public void scalarMultEquals(double s) {
		CommonOps.scale(s, myMatrix);
	}

	@Override
	public void set(int row, int col, double value) {
		myMatrix.set(row, col, value);
	}

	@Override
	public void setChecked(int row, int col, double value) {
		myMatrix.set(row, col, value);
	}

	@Override
	public DoubleMatrix solve(DoubleMatrix Y) {
		DenseMatrix64F data = ((EJMLDoubleMatrix) Y).myMatrix;
		DenseMatrix64F result = new DenseMatrix64F(myMatrix.numCols, data.numCols);
		LinearSolver<DenseMatrix64F> solver = LinearSolverFactory.leastSquares(myMatrix.numCols, data.numCols);
		solver.setA(myMatrix);
		solver.solve(data, result);
		return new EJMLDoubleMatrix(result);
	}

	@Override
	public DoubleMatrix tcrossproduct() {
		int ncol = myMatrix.numRows;
		DenseMatrix64F result = new DenseMatrix64F(ncol, ncol);
		CommonOps.multTransB(myMatrix, myMatrix, result);
		return new EJMLDoubleMatrix(result);
	}

	@Override
	public DoubleMatrix tcrossproduct(DoubleMatrix dm) {
		DenseMatrix64F other = ((EJMLDoubleMatrix) dm).myMatrix;
		int nrow = myMatrix.numRows;
		int ncol = other.numRows;
		DenseMatrix64F result = new DenseMatrix64F(nrow, ncol);
		CommonOps.multTransB(myMatrix, other, result);
		return new EJMLDoubleMatrix(result);
	}

	@Override
	public DoubleMatrix transpose() {
		DenseMatrix64F result = new DenseMatrix64F(myMatrix.numCols, myMatrix.numRows);
		CommonOps.transpose(myMatrix, result);
		return new EJMLDoubleMatrix(result);
	}

	@Override
	public DoubleMatrix getSelection(int[] rows, int[] columns) {
		if (rows == null) {
			if (columns == null) return copy();
			int nrow = myMatrix.numRows;
			int ncol = columns.length;
			DoubleMatrix result = DoubleMatrixFactory.DEFAULT.make(nrow, ncol);
			for (int r = 0; r < nrow; r++) {
				for (int c = 0; c < ncol; c++) {
					result.set(r, c, myMatrix.get(r, columns[c]));
				}
			}
			return result;
		} else if (columns == null) {
			int nrow = rows.length;
			int ncol = myMatrix.numCols;
			DoubleMatrix result = DoubleMatrixFactory.DEFAULT.make(nrow, ncol);
			for (int r = 0; r < nrow; r++) {
				for (int c = 0; c < ncol; c++) {
					result.set(r, c, myMatrix.get(rows[r], c));
				}
			}
			return result;
		} else {
			int nrow = rows.length;
			int ncol = columns.length;
			DoubleMatrix result = DoubleMatrixFactory.DEFAULT.make(nrow, ncol);
			for (int r = 0; r < nrow; r++) {
				for (int c = 0; c < ncol; c++) {
					result.set(r, c, myMatrix.get(rows[r], columns[c]));
				}
			}
			return result;
		}
	}

	@Override
	public String toString() {
		int nrows = Math.min(25, myMatrix.numRows);
		int ncols = Math.min(25, myMatrix.numCols);
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < nrows; i++) {
			for (int j = 0; j < ncols; j++) {
				if (j > 0) sb.append(" ");
				sb.append(myMatrix.get(i,j));
			}
			sb.append("\n");
		}
		return sb.toString();
	}

	@Override
	public double columnSum(int column) {
		int n = myMatrix.numRows;
		DenseMatrix64F vector = new DenseMatrix64F(n, 1);
		SpecializedOps.subvector(myMatrix, 0, column, n, false, 0, vector);
		return CommonOps.elementSum(vector);
	}

	@Override
	public double rowSum(int row) {
		int n = myMatrix.numCols;
		DenseMatrix64F vector = new DenseMatrix64F(n, 1);
		SpecializedOps.subvector(myMatrix, 0, row, n, false, 0, vector);
		return CommonOps.elementSum(vector);
	}

}
