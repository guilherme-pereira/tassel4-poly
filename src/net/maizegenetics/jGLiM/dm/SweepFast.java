package net.maizegenetics.jGLiM.dm;

import java.util.ArrayList;

import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;

public class SweepFast {
	//A is a double array representing a symmetric matrix. It holds the upper triangular part of the matrix,
	//including the diagonal, in row major order.
    private double[] A;
    private int dimA;
    private boolean[] singular;
    private double[] V;
    private double[] Dmin;
    
    public static final double TOL = 1e-12;

    public SweepFast(DoubleMatrix X, DoubleMatrix y) {
    	this(X.crossproduct(), X.crossproduct(y), y.crossproduct().get(0, 0));
    }
    
    public SweepFast(DoubleMatrix XTX, DoubleMatrix XTy, double yty) {
    	int dimXTX = XTX.numberOfRows();
    	dimA = dimXTX + 1;
    	int n = dimA * (dimA + 1) / 2;
    	A = new double[n];
    	
    	int count = 0;
    	for (int r = 0; r < dimXTX; r++) {
    		for (int c = r; c < dimXTX; c++) {
    			A[count++] = XTX.get(r, c);
    		}
    		count++;
    	}
    	
    	count = dimXTX;
    	int incr = dimXTX;
    	for (int r = 0; r < dimXTX; r++) {
    		A[count] = XTy.get(r, 0);
    		count += incr--;
    	}
    	
    	A[n-1] = yty;
    	
    	singular = new boolean[dimA];
    	V = new double[dimA];
    	for (int i = 0; i < dimA; i++) {
    		singular[i] = false;
    		V[i] = 1;
    	}
    	
    	initializeDmin();
    }
   
    public SweepFast(DoubleMatrix[][] xtxMatrices, DoubleMatrix[] xtyMatrices, double yty) {
    	init(xtxMatrices, xtyMatrices, yty);
    }
    
    public SweepFast() {}
    
    private void init(DoubleMatrix[][] xtxMatrices, DoubleMatrix[] xtyMatrices, double yty) {
    	int nx = xtyMatrices.length;
    	dimA = 0;
    	int[] dimX = new int[nx];
    	for (int i = 0; i < nx; i++) {
    		int num = xtxMatrices[0][i].numberOfColumns();
    		dimX[i] = num;
    		dimA += num;
    	}
    	dimA++; 
    	int n = dimA * (dimA + 1) / 2;
    	A = new double[n];

    	int count = 0;
    	for (int imain = 0; imain < nx; imain++) {
    		for (int isub = 0; isub < dimX[imain]; isub++) {
    			for (int j = isub; j < dimX[imain]; j++)
    				A[count++] = xtxMatrices[imain][imain].get(isub, j);
    			for (int k = imain + 1; k < nx; k++) {
    				for (int j = 0; j < dimX[k]; j++) 
    					A[count++] = xtxMatrices[imain][k].get(isub, j);
    			}
    			A[count++] = xtyMatrices[imain].get(isub, 0);
    		}
    	}
    	A[count] = yty;

    	singular = new boolean[dimA];
    	V = new double[dimA];
    	for (int i = 0; i < dimA; i++) {
    		singular[i] = false;
    		V[i] = 1;
    	}
    	
    	initializeDmin();
    }
    
    private void initializeDmin() {
    	Dmin = new double[dimA];
    	for (int i = 0; i < dimA; i++) Dmin[i] = TOL;
    }
    
    //returns false if the column was not swept (linearly combination of columns already in the model, i.e. singular)
    public boolean revg2sweep(int column) {
        if (singular[column]) {
            singular[column] = false;
            return false;
        }
        int diagndx = diagonal(column);
        double D = A[diagonal(column)];
        if (Math.abs(D) < Dmin[column]) {
            singular[column] = true;
            return false;
        }
        
        //get the starting values for this column 
        double[] originalColumnValues = new double[dimA];
        int ndx = column;
        int incr = dimA - 1;
        for (int i = 0; i <= column; i++) {
        	originalColumnValues[i] = A[ndx];
        	ndx += incr--;
        }
        ndx = diagndx + 1;
        for (int i = column + 1; i < dimA; i++) {
        	originalColumnValues[i] = A[ndx];
        	ndx++;
        }
        
        //now sweep the matrix for this column
        int count = 0;
        double vcol = V[column];
        for (int i = 0; i < column; i++) {
        	for (int j = i; j < dimA; j++) {
        		double B = originalColumnValues[i] / D;
        		if (column < j) A[count] -= B * originalColumnValues[j];
        		if (column > j) A[count] -= B * originalColumnValues[j] * V[j] * vcol;
        		count++;
        	}
        }
        count += dimA - column;
        for (int i = column + 1; i < dimA; i++) {
        	for (int j = i; j < dimA; j++) {
        		double B = V[i] *vcol * originalColumnValues[i] / D;
        		if (column < j) A[count] -= B * originalColumnValues[j];
        		if (column > j) A[count] -= B * originalColumnValues[j] * V[j] * vcol;
        		count++;
        	}
        }
        
        //divide the original row by D
        ndx = column;
        incr = dimA - 1;
        for (int i = 0; i <= column; i++) {
        	A[ndx] = -1 * originalColumnValues[i] / D;
        	ndx += incr--;
        }
        ndx = diagndx;
        for (int j = column; j < dimA; j++) {
        	A[ndx] = originalColumnValues[j] / D;
        	ndx++;
        }
        
        A[diagndx] = 1 / D;
        V[column] = V[column] * -1;
        
        return true;
    }
    
    private int diagonal(int d) {
    	if (d == 0) return 0;
    	return d * dimA - d * (d - 1) / 2;
    }
    
    //since A is upper triangular, col >= row
    private int coord(int row, int col) {
    	if (row == 0) return col;
    	return row * dimA - row * (row - 1) / 2 + col - row;
    }
    
    public int sweepSingularColumns() {
        int swept = 0;
        for (int i = 0; i < singular.length; i++) {
            if (singular[i]) {
                singular[i] = false;
                if (revg2sweep(i)) swept++; 
            }
        }
        return swept;
    }

    public void setcssDmin(double css[]) {
    	if (css == null) {
    		int n = dimA - 1;
    		Dmin = new double[n];
    		for (int i = 0; i < n; i++) Dmin[i] = TOL;
    	} else {
    		int n = css.length;
    		Dmin = new double[n];
    		for (int i = 0; i < n; i++) {
    			if (css[i] < 1) Dmin[i] = TOL;
    			else Dmin[i] = css[i] * TOL;
    		}
    	}
    }
    
    public void setDminFromA() {
    	int n = dimA;
    	Dmin = new double[n];
    	for (int i = 0; i < n; i++) Dmin[i] = Math.max(A[diagonal(i)] * TOL, TOL);
    }

    public double getResidualSS() {
    	return A[A.length - 1];
    }
    
    public double[] getBeta() {
    	int n = dimA - 1;
    	double[] beta = new double[n];
    	int incr = n;
    	int count = n;
    	for (int i = 0; i < n; i++) {
    		beta[i] = A[count];
    		count += incr--;
    	}
    	return beta;
    }
    
    public DoubleMatrix getXTXpart() {
    	int n = dimA - 1;
    	DoubleMatrix Amatrix = DoubleMatrixFactory.DEFAULT.make(n, n);
    	int count = 0;
    	for (int i = 0; i < n; i++) {
    		for (int j = i; j < n; j++) {
    			Amatrix.set(i, j, A[count++]);
    		}
    		count++;
    	}

    	for (int i = 0; i < n; i++) {
    		for (int j = i + 1; j < n; j++) {
    			Amatrix.set(j, i, Amatrix.get(i, j) * V[i] * V[j]);
    		}
    	}

    	return Amatrix;
    }
    
    public DoubleMatrix getA() {
    	DoubleMatrix Amatrix = DoubleMatrixFactory.DEFAULT.make(dimA, dimA);
    	int count = 0;
    	for (int i = 0; i < dimA; i++) {
    		for (int j = i; j < dimA; j++) {
    			Amatrix.set(i, j, A[count++]);
    		}
    	}
    	
    	count = 0;
    	for (int i = 0; i < dimA; i++) {
    		for (int j = i; j < dimA; j++) {
    			Amatrix.set(j, i, A[count++] * V[i] * V[j]);
    		}
    	}
    	return Amatrix;
    }
    
    public DoubleMatrix getSubsetOfA(boolean[] exclude) {
    	int n = 0;
    	for (boolean delete : exclude) if (!delete) n++;
    	DoubleMatrix subset = DoubleMatrixFactory.DEFAULT.make(n, n);
    	
    	int row = 0;
    	int col = 0;
    	int count = 0;
    	
    	for (int i = 0; i < n; i++) {
    		for (int j = i; j < n; j++) {
    			while (exclude[col] || exclude[row]) {
    				count++;
    				col++;
    				if (col == dimA) col = ++row;
    			}
    			subset.set(i, j, A[count++]);
				col++;
				if (col == dimA) col = ++row;
    		}
    	}

    	for (int i = 0; i < n; i++) {
    		for (int j = i + 1; j < n; j++) {
    			subset.set(j, i, subset.get(i, j));
    		}
    	}
    	
    	return subset;
    }
    
    public double[] getUTA() {
    	return A;
    }
    
    public void setUTA(double[] A) {
    	this.A = A;
    	dimA = (int) (Math.floor(Math.sqrt(2 * A.length))); 
    	
    	singular = new boolean[dimA];
    	V = new double[dimA];
    	for (int i = 0; i < dimA; i++) {
    		singular[i] = false;
    		V[i] = 1;
    	}
    }
    
    public void fullSweepSetDmin() {
    	Dmin = new double[dimA];
    	revg2sweep(0);
    	for (int i = 0; i < dimA; i++) Dmin[i] = Math.max(A[diagonal(i)] * TOL, TOL);

    	for (int k = 1; k < dimA; k++) {
    		revg2sweep(k);
    	}
    }
    
    public void XTXSweepSetDmin() {
    	Dmin = new double[dimA];
    	revg2sweep(0);
    	for (int i = 0; i < dimA; i++) Dmin[i] = Math.max(A[diagonal(i)] * TOL, TOL);

    	for (int k = 1; k < dimA - 1; k++) {
    		revg2sweep(k);
    	}
    }
    
    public SweepFast copy() {
    	SweepFast copy = new SweepFast();
    	
    	int n = A.length;
    	double[] newA = new double[n];
    	System.arraycopy(A, 0, newA, 0, n);
    	copy.A = newA;
    	
    	copy.dimA = dimA;
 
    	n = Dmin.length;
    	double[] newDmin = new double[n];
    	System.arraycopy(Dmin, 0, newDmin, 0, n);
    	copy.Dmin = newDmin;

    	n = singular.length;
    	boolean[] newsingular = new boolean[n];
    	System.arraycopy(singular, 0, newsingular, 0, n);
    	copy.singular = newsingular;
    	
    	n = V.length;
    	double[] newV = new double[n];
    	System.arraycopy(V, 0, newV, 0, n);
    	copy.V = newV;
    	
    	return copy;
    }
    
    public int getDimensionOfA() { return dimA; }
    
    public boolean isSingular(int column) { return singular[column]; }
    
    public static double[] UTAFromDoubleMatrix(DoubleMatrix A) {
    	int nrows = A.numberOfRows();
    	int n = nrows * (nrows + 1) / 2;
    	double[] uta = new double[n];
    	
    	int count = 0;
    	for (int i = 0; i < nrows; i ++) {
    		for (int j = i; j < nrows; j++) {
    			uta[count++] = A.get(i,j);
    		}
    	}
    	
    	return uta;
    }
    
}
