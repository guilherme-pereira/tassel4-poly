package net.maizegenetics.matrixalgebra.decomposition;

import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;

//modified from the EigenvalueDecomposition of Colt version 1.2 by converting A and V from 
//2-d arrays to 1-d arrays in order to improve performance

public class SymmetricEigenvalueDecomposition implements EigenvalueDecomposition {
	double[] V; //storage of eigenvectors
	int n;
	int size;
	double[] d,e; //storage of eigenvalues
	
	public SymmetricEigenvalueDecomposition(double[] matrix) {
		size = matrix.length;
		V = new double[size];
		for (int i = 0; i < size; i++) V[i] = matrix[i];
		n = (int) Math.round(Math.sqrt(size));
		d = new double[n];
		e = new double[n];
		
		// Tridiagonalize.
		tred2();
		
		// Diagonalize.
		tql2();

	}
	
	public double[] getRealEigenvalues(){
		return d;
	}
	
	public double[] getEigenvectorsAs1dArray() {
		return V;
	}
	
	/**
	Symmetric Householder reduction to tridiagonal form.
	*/
	private void tred2 () {
	   //  This is derived from the Algol procedures tred2 by
	   //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
	   //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
	   //  Fortran subroutine in EISPACK.

		  int start = size - n;
		  for (int j = 0; j < n; j++) {
			 d[j] = V[start + j] ; // d equals the last row of V
		  }

		  // Householder reduction to tridiagonal form.
		  int starti = size;
		  for (int i = n-1; i > 0; i--) {
	   
			 // Scale to avoid under/overflow.
			 starti -= n; 
			 double scale = 0.0;
			 double h = 0.0;
			 for (int k = 0; k < i; k++) {
				scale = scale + Math.abs(d[k]);
			 }
			 if (scale == 0.0) {
				e[i] = d[i-1];
				for (int j = 0; j < i; j++) {
				   d[j] = V[starti - n + j];
				   V[starti + j] = 0.0;
				   V[j * n + i] = 0.0;
				}
			 } else {
	   
				// Generate Householder vector.
	   
				for (int k = 0; k < i; k++) {
				   d[k] /= scale;
				   h += d[k] * d[k];
				}
				double f = d[i-1];
				double g = Math.sqrt(h);
				if (f > 0) {
				   g = -g;
				}
				e[i] = scale * g;
				h = h - f * g;
				d[i-1] = f - g;
				for (int j = 0; j < i; j++) {
				   e[j] = 0.0;
				}
	   
				// Apply similarity transformation to remaining columns.
				int startj = 0;
				for (int j = 0; j < i; j++) {
				   f = d[j];
				   V[startj + i] = f;
				   g = e[j] + V[startj + j] * f;
				   int cell = startj + n + j;
				   for (int k = j+1; k <= i-1; k++) {
					  g += V[cell] * d[k];
					  e[k] += V[cell] * f;
					  cell += n;
				   }
				   e[j] = g;
				   startj += n;
				}
				f = 0.0;
				for (int j = 0; j < i; j++) {
				   e[j] /= h;
				   f += e[j] * d[j];
				}
				double hh = f / (h + h);
				for (int j = 0; j < i; j++) {
				   e[j] -= hh * d[j];
				}
				for (int j = 0; j < i; j++) {
				   f = d[j];
				   g = e[j];
				   int cell = j * n + j;
				   for (int k = j; k <= i-1; k++) {
					  V[cell] -= (f * e[k] + g * d[k]);
					  cell += n;
				   }
				   cell = starti - n + j;
				   d[j] = V[cell];
				   cell += n;
				   V[cell] = 0.0;
				}
			 }
			 d[i] = h;
		  }
	   
		  // Accumulate transformations.
		  int rowi = 0;
		  int lastrow = size - n;
		  int lastcol = n - 1;
		  int cell;
		  for (int i = 0; i < n-1; i++) {
			 cell = rowi + i; 
			 V[lastrow + i] = V[cell];
			 V[cell] = 1.0;
			 double h = d[i+1];
			 if (h != 0.0) {
				int celliplus1 = i + 1; 
				for (int k = 0; k <= i; k++) {
				   d[k] = V[celliplus1] / h;
				   celliplus1 += n;
				}
				for (int j = 0; j <= i; j++) {
				   double g = 0.0;
				   celliplus1 = i + 1;
				   int cellj = j;
				   for (int k = 0; k <= i; k++) {
					  g += V[celliplus1] * V[cellj];
					  celliplus1 += n;
					  cellj += n;
				   }
				   cellj = j;
				   for (int k = 0; k <= i; k++) {
					  V[cellj] -= g * d[k];
					  cellj += n;
				   }
				}
			 }
			 int celliplus1 = i + 1;
			 for (int k = 0; k <= i; k++) {
				V[celliplus1] = 0.0;
				celliplus1 +=n;
			 }
			 rowi += n;
		  }
		  for (int j = 0; j < n; j++) {
			 cell = lastrow + j;
			 d[j] = V[cell];
			 V[cell] = 0.0;
		  }
		  V[lastrow + lastcol] = 1.0;
		  e[0] = 0.0;
	   }
	
	/**
	Symmetric tridiagonal QL algorithm.
	*/
	private void tql2 () {

		//  This is derived from the Algol procedures tql2, by
		//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
		//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
		//  Fortran subroutine in EISPACK.
	   
		  for (int i = 1; i < n; i++) {
			 e[i-1] = e[i];
		  }
		  e[n-1] = 0.0;
	   
		  double f = 0.0;
		  double tst1 = 0.0;
		  double eps = Math.pow(2.0,-52.0);
		  for (int l = 0; l < n; l++) {

			 // Find small subdiagonal element
	   
			 tst1 = Math.max(tst1,Math.abs(d[l]) + Math.abs(e[l]));
			 int m = l;
			 while (m < n) {
				if (Math.abs(e[m]) <= eps*tst1) {
				   break;
				}
				m++;
			 }
	   
			 // If m == l, d[l] is an eigenvalue,
			 // otherwise, iterate.
	   
			 if (m > l) {
				int iter = 0;
				do {
				   iter = iter + 1;  // (Could check iteration count here.)
	   
				   // Compute implicit shift
	   
				   double g = d[l];
				   double p = (d[l+1] - g) / (2.0 * e[l]);
				   double r = hypot(p,1.0);
				   if (p < 0) {
					  r = -r;
				   }
				   d[l] = e[l] / (p + r);
				   d[l+1] = e[l] * (p + r);
				   double dl1 = d[l+1];
				   double h = g - d[l];
				   for (int i = l+2; i < n; i++) {
					  d[i] -= h;
				   }
				   f = f + h;
	   
				   // Implicit QL transformation.
	   
				   p = d[m];
				   double c = 1.0;
				   double c2 = c;
				   double c3 = c;
				   double el1 = e[l+1];
				   double s = 0.0;
				   double s2 = 0.0;
				   for (int i = m-1; i >= l; i--) {
					  c3 = c2;
					  c2 = c;
					  s2 = s;
					  g = c * e[i];
					  h = c * p;
					  r = hypot(p,e[i]);
					  e[i+1] = s * r;
					  s = e[i] / r;
					  c = p / r;
					  p = c * d[i] - s * g;
					  d[i+1] = h + s * (c * g + s * d[i]);
	   
					  // Accumulate transformation.
					  int celliplus1 = i + 1;
					  int celli = 1;
					  for (int k = 0; k < n; k++) {
						 h = V[celliplus1];
						 V[celliplus1] = s * V[celli] + c * h;
						 V[celli] = c * V[celli] - s * h;
						 celliplus1 += n;
						 celli += n;
					  }
				   }
				   p = -s * s2 * c3 * el1 * e[l] / dl1;
				   e[l] = s * p;
				   d[l] = c * p;
	   
				   // Check for convergence.
	   
				} while (Math.abs(e[l]) > eps*tst1);
			 }
			 d[l] = d[l] + f;
			 e[l] = 0.0;
		  }
		 
		  // Sort eigenvalues and corresponding vectors.
	   
		  for (int i = 0; i < n-1; i++) {
			 int k = i;
			 double p = d[i];
			 for (int j = i+1; j < n; j++) {
				if (d[j] < p) {
				   k = j;
				   p = d[j];
				}
			 }
			 if (k != i) {
				d[k] = d[i];
				d[i] = p;
				int celli = i;
				int cellk = k;
				for (int j = 0; j < n; j++) {
				   p = V[celli];
				   V[celli] = V[cellk];
				   V[cellk] = p;
				   celli += n;
				   cellk += n;
				}
			 }
		  }
	   }
	
	//cern.colt.matrix.linalg.Algebra.hypot
	private double hypot(double a, double b) {
		double r;
		if (Math.abs(a) > Math.abs(b)) {
			r = b/a;
			r = Math.abs(a)*Math.sqrt(1+r*r);
		} else if (b != 0) {
			r = a/b;
			r = Math.abs(b)*Math.sqrt(1+r*r);
		} else {
			r = 0.0;
		}
		return r;
	}

	@Override
	public double getEigenvalue(int i) {
		return d[i];
	}

	@Override
	public DoubleMatrix getEigenvalueMatrix() {
		return DoubleMatrixFactory.DEFAULT.diagonal(d);
	}

	@Override
	public double[] getEigenvalues() {
		return d;
	}

	@Override
	public DoubleMatrix getEigenvectors() {
		return DoubleMatrixFactory.DEFAULT.make(n, n, V);
	}

	
}
