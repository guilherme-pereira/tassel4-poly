package net.maizegenetics.jGLiM.dm;

import java.util.Arrays;

import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;

public class NestedCovariateModelEffect implements ModelEffect {
	private Object id;
	private final int size;
	private final double[] covariate;
	private final FactorModelEffect fme;

	public NestedCovariateModelEffect(CovariateModelEffect cme, FactorModelEffect fme) {
		covariate = cme.getCovariate();
		size = covariate.length;
		this.fme = fme;
	}
	
	public NestedCovariateModelEffect(double[] covariate, FactorModelEffect fme) {
		this.covariate = covariate;
		size = covariate.length;
		this.fme = fme;
	}
	
	public NestedCovariateModelEffect(float[] covariate, FactorModelEffect fme) {
		int n = covariate.length;
		double[] cov = new double[n];
		for (int i = 0; i < n; i++) cov[i] = covariate[i];
		this.covariate = cov;
		size = covariate.length;
		this.fme = fme;
	}
	
	private NestedCovariateModelEffect(Object id, int size, double[] covariate, FactorModelEffect fme) {
		this.id = id;
		this.covariate = Arrays.copyOf(covariate, covariate.length);
		this.size = size;
		this.fme = (FactorModelEffect) fme.getCopy();
	}
	
	public static NestedCovariateModelEffect getInstance(CovariateModelEffect cme, ModelEffect me) {
		if (me instanceof FactorModelEffect) {
			return new NestedCovariateModelEffect(cme, (FactorModelEffect) me);
		} else {
			return null;
		}
	}
	
	public static NestedCovariateModelEffect getInstance(double[] covariate, ModelEffect me) {
		if (me instanceof FactorModelEffect) {
			return new NestedCovariateModelEffect(covariate, (FactorModelEffect) me);
		} else {
			return null;
		}
	}
	
	@Override
	public Object getID() {
		return id;
	}

	@Override
	public int[] getLevelCounts() {
		return fme.getLevelCounts();
	}

	@Override
	public int getNumberOfLevels() {
		return fme.getNumberOfLevels();
	}

	@Override
	public int getSize() {
		return size;
	}

	@Override
	public DoubleMatrix getX() {
		int nlevels = fme.getNumberOfLevels();
		int[] levels = fme.getLevels();
		DoubleMatrix x = DoubleMatrixFactory.DEFAULT.make(size, nlevels);
		for (int i = 0; i < size; i++) x.set(i, levels[i], covariate[i]);
		return x;
	}

	@Override
	public DoubleMatrix getXtX() {
		int nlevels = fme.getNumberOfLevels();
		int[] levels = fme.getLevels();
		DoubleMatrix xtx = DoubleMatrixFactory.DEFAULT.make(nlevels, nlevels, 0);
		for (int i = 0; i < size; i++) xtx.set(levels[i], levels[i], xtx.get(levels[i], levels[i]) + covariate[i] * covariate[i]);
		return xtx;
	}

	@Override
	public DoubleMatrix getXty(double[] y) {
		int nlevels = fme.getNumberOfLevels();
		int[] levels = fme.getLevels();
		DoubleMatrix xty = DoubleMatrixFactory.DEFAULT.make(nlevels, 1, 0);
		for (int i = 0; i < size; i++) xty.set(levels[i], 0, xty.get(levels[i], 0) + covariate[i] * y[i]);
		return xty;
	}

	@Override
	public DoubleMatrix getyhat(DoubleMatrix beta) {
		int[] levels = fme.getLevels();
		DoubleMatrix yhat = DoubleMatrixFactory.DEFAULT.make(size,1);
		for (int i = 0; i < size; i++) {
			yhat.set(i, 0, beta.get(levels[i], 0));
		}
		return yhat;
	}

	@Override
	public DoubleMatrix getyhat(double[] beta) {
		int[] levels = fme.getLevels();
		DoubleMatrix yhat = DoubleMatrixFactory.DEFAULT.make(size,1);
		for (int i = 0; i < size; i++) {
			yhat.set(i, 0, beta[levels[i]]);
		}
		return yhat;
	}

	@Override
	public void setID(Object id) {
		this.id = id;
	}

	public DoubleMatrix getXtX2(NestedCovariateModelEffect ncme) {
		FactorModelEffect fme2 = ncme.getFactorModelEffect();
		int nrows = fme.getNumberOfLevels();
		int ncols = fme2.getNumberOfLevels();
		int[] levels = fme.getLevels();
		int[] levels2 = fme2.getLevels();
		DoubleMatrix xtx2 = DoubleMatrixFactory.DEFAULT.make(nrows, ncols, 0);
		for (int i = 0; i < size; i++) {
			xtx2.set(levels[i], levels2[i], xtx2.get(levels[i], levels2[i]) + covariate[i] * ncme.covariate[i]);
		}
		return xtx2;
	}
	
	public FactorModelEffect getFactorModelEffect() { return fme; }

	@Override
	public ModelEffect getCopy() {
		return new NestedCovariateModelEffect(id, size, covariate, fme);
	}

	@Override
	public int getEffectSize() {
		return getNumberOfLevels();
	}
	
}
