package net.maizegenetics.jGLiM.dm;

import java.util.Arrays;

import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;
import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.jet.math.Functions;

public class CovariateModelEffect implements ModelEffect {
    private final double[] covariate;
    private final int size;
    private final double sum;
    private final double sumsq;
    private Object id = null;

    public CovariateModelEffect(double[] covariate) {
    	this.covariate = covariate;
    	size = covariate.length;
    	double s = 0;
    	double ss = 0;
    	for (double cov : covariate) {
    		s += cov;
    		ss += cov * cov;
    	}
    	sum = s;
    	sumsq = ss;
    }
    
    public CovariateModelEffect(double[] covariate, Object id) {
    	this(covariate);
    	this.id = id;
    }
    
    private CovariateModelEffect(double[] covariate, int size, double sum, double sumsq, Object id){
    	this.covariate = Arrays.copyOf(covariate, covariate.length);
    	this.size = size;
    	this.sum = sum;
    	this.sumsq = sumsq;
    	this.id = id;
    }
    
	@Override
	public Object getID() {
		return id;
	}

	@Override
	public int getNumberOfLevels() {
		return 1;
	}

	@Override
	public void setID(Object id) {
		this.id = id;
	}

	@Override
	public int[] getLevelCounts() {
		return new int[]{size};
	}

	@Override
	public int getSize() {
		return covariate.length;
	}

	@Override
	public DoubleMatrix getX() {
		return DoubleMatrixFactory.DEFAULT.make(covariate.length, 1, covariate);
	}

	@Override
	public DoubleMatrix getXtX() {
		return DoubleMatrixFactory.DEFAULT.make(1, 1, sumsq);
	}

	@Override
	public DoubleMatrix getXty(double[] y) {
		double sumprod = 0;
		for (int i = 0; i < size; i++) sumprod += covariate[i] * y[i]; 
		return DoubleMatrixFactory.DEFAULT.make(1,1,sumprod);
	}

	@Override
	public DoubleMatrix getyhat(DoubleMatrix beta) {
		double scalar = beta.get(0, 0);
		DoubleMatrix yhat = DoubleMatrixFactory.DEFAULT.make(size, 1, covariate);
		yhat.scalarMultEquals(scalar);
		return yhat;
	}

	@Override
	public DoubleMatrix getyhat(double[] beta) {
		double scalar = beta[0];
		DoubleMatrix yhat = DoubleMatrixFactory.DEFAULT.make(size, 1, covariate);
		yhat.scalarMultEquals(scalar);
		return yhat;
	}

	public DoubleMatrix getXtX2(CovariateModelEffect cme) {
		double sumprod = 0;
		for (int i = 0; i < size; i++) sumprod += covariate[i] * cme.covariate[i];
		return DoubleMatrixFactory.DEFAULT.make(1, 1, sumprod);
	}
	
    public double[] getCovariate() { return covariate; }
    
    public double getSum() { return sum; }

    public double getSumSquares() { return sumsq; }

	@Override
	public ModelEffect getCopy() {
		return new CovariateModelEffect(covariate, size, sum, sumsq, id);
	}

	@Override
	public int getEffectSize() {
		return 1;
	}
    
}
