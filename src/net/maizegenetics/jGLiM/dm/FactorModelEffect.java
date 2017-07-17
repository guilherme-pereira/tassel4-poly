package net.maizegenetics.jGLiM.dm;

import java.util.Arrays;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;

public class FactorModelEffect implements ModelEffect {
	private Object id;
	private final int[] levels;
	private final int[] levelCounts;
	private final int size;
	private final int numberOfLevels;
	private boolean restricted = false;
	
	public FactorModelEffect(int[] levels, boolean restricted) {
		this.levels = levels;
		size = levels.length;
		this.restricted = restricted;
		int maxLevel = 0;
		for (int level : levels) {
			maxLevel = Math.max(maxLevel, level);
		}
		numberOfLevels = maxLevel + 1;
		int[] counts = new int[numberOfLevels];
		for (int level : levels) {
			counts[level]++;
		}
		levelCounts = counts;
	}
	
	public FactorModelEffect(int[] levels, boolean restricted, Object id) {
		this(levels, restricted);
		this.id = id;
	}
	
	private FactorModelEffect(Object id, int[] levels, int[] levelCounts, int size, int numberOfLevels, boolean restricted) {
		this.id = id;
		this.levels = Arrays.copyOf(levels, levels.length);
		this.levelCounts = Arrays.copyOf(levelCounts, levelCounts.length);
		this.size = size;
		this.numberOfLevels = numberOfLevels;
		this.restricted = restricted;
	}
	
	@Override
	public Object getID() {
		return id;
	}

	@Override
	public int getSize() {
		return size;
	}

	@Override
	public DoubleMatrix getX() {
		if (restricted) return getX_restricted();
		return getX_unrestricted();
	}

	@Override
	public DoubleMatrix getXtX() {
		if (restricted) return getXtX_restricted();
		return getXtX_unrestricted();
	}

	@Override
	public DoubleMatrix getXty(double[] y) {
		if (restricted) return getXty_restricted(y);
		return getXty_unrestricted(y);
	}

	@Override
	public void setID(Object id) {
		this.id = id;
	}

	/**
	 * @param fme a FactorModelEffect
	 * @return the matrix product of this ModelEffect transposed and fme
	 */
	public DoubleMatrix getXtX2(FactorModelEffect fme) {
		int n1 = numberOfLevels;
		if (restricted) n1--;
		int n2 = fme.numberOfLevels;
		if (fme.restricted) n2--;
		
		DoubleMatrix xtx2 = DoubleMatrixFactory.DEFAULT.make(n1, n2, 0);
		for (int i = 0; i < size; i++) {
			int row = levels[i];
			int col = fme.levels[i];
			if (row < n1 && col < n2) xtx2.set(row, col, xtx2.get(row, col) + 1);
		}
		return xtx2;

	}

	@Override
	public DoubleMatrix getyhat(DoubleMatrix beta) {
		if (restricted) return getyhat_restricted(beta);
		return getyhat_unrestricted(beta);
	}

	@Override
	public DoubleMatrix getyhat(double[] beta) {
		if (restricted) return getyhat_restricted(beta);
		return getyhat_unrestricted(beta);
	}

	/**
	 * @return the level of each observation for this factor
	 */
	public int[] getLevels() { return levels; }
	
	public boolean getRestricted() { return restricted; }
	
	public void setRestricted(boolean restricted) { this.restricted = restricted; }
	
	public int getNumberOfLevels() {
		return numberOfLevels; 
	}
	
	@Override
	public int[] getLevelCounts() {
		return levelCounts;
	}
	
	//private methods
	
	private DoubleMatrix getX_unrestricted() {
		DoubleMatrix X = DoubleMatrixFactory.DEFAULT.make(size, numberOfLevels, 0);
		for (int i = 0; i < size; i++) {
			X.set(i, levels[i], 1);
		}
		return X;
	}
	
	private DoubleMatrix getX_restricted() {
		DoubleMatrix X = DoubleMatrixFactory.DEFAULT.make(size, numberOfLevels - 1, 0);
		for (int i = 0; i < size; i++) {
			if (levels[i] < numberOfLevels - 1) X.set(i, levels[i], 1);
		}
		return X;
	}
	
	private DoubleMatrix getXtX_unrestricted() {
		double[] counts = new double[numberOfLevels];
		for (int i = 0; i < numberOfLevels; i++) counts[i] = (int) levelCounts[i];
		DoubleMatrix xtx = DoubleMatrixFactory.DEFAULT.diagonal(counts);
		return xtx;
	}
	
	private DoubleMatrix getXtX_restricted() {
		double[] counts = new double[numberOfLevels - 1];
		for (int i = 0; i < numberOfLevels - 1; i++) counts[i] = (int) levelCounts[i];
		DoubleMatrix xtx = DoubleMatrixFactory.DEFAULT.diagonal(counts);
		return xtx;
	}
	
	private DoubleMatrix getXty_unrestricted(double[] y) {
		double[] sums = new double[numberOfLevels];
		for (int i = 0; i < size; i++) {
			sums[levels[i]] += y[i];
		}
		return DoubleMatrixFactory.DEFAULT.make(numberOfLevels, 1, sums);
	}

	private DoubleMatrix getXty_restricted(double[] y) {
		double[] sums = new double[numberOfLevels - 1];
		for (int i = 0; i < size; i++) {
			if (levels[i] < numberOfLevels - 1) sums[levels[i]] += y[i];
		}
		return DoubleMatrixFactory.DEFAULT.make(numberOfLevels - 1, 1, sums);
	}
	
	private DoubleMatrix getyhat_unrestricted(DoubleMatrix beta) {
		if (beta.numberOfRows() != numberOfLevels) throw new IllegalArgumentException("Number of factor levels is not the same as the size of beta.");
		double[] yhat = new double[size];
		for (int i = 0; i < size; i++) yhat[i] = beta.get(levels[i], 0);
		return DoubleMatrixFactory.DEFAULT.make(size, 1, yhat);
	}

	private DoubleMatrix getyhat_restricted(DoubleMatrix beta) {
		if (beta.numberOfRows() != numberOfLevels - 1) throw new IllegalArgumentException("Number of factor levels is not the same as the size of beta.");
		double[] yhat = new double[size];
		for (int i = 0; i < size; i++) {
			if (levels[i] < numberOfLevels - 1) yhat[i] = beta.get(levels[i], 0);	
		}
		return DoubleMatrixFactory.DEFAULT.make(size, 1, yhat);
	}

	private DoubleMatrix getyhat_unrestricted(double[] beta) {
		if (beta.length != numberOfLevels) throw new IllegalArgumentException("Number of factor levels is not the same as the size of beta.");
		double[] yhat = new double[size];
		for (int i = 0; i < size; i++) yhat[i] = beta[levels[i]];
		return DoubleMatrixFactory.DEFAULT.make(size, 1, yhat);
	}

	private DoubleMatrix getyhat_restricted(double[] beta) {
		if (beta.length != numberOfLevels - 1) throw new IllegalArgumentException("Number of factor levels is not the same as the size of beta.");
		double[] yhat = new double[size];
		for (int i = 0; i < size; i++) {
			if (levels[i] < numberOfLevels - 1) yhat[i] = beta[levels[i]];	
		}
		return DoubleMatrixFactory.DEFAULT.make(size, 1, yhat);
	}

	@Override
	public ModelEffect getCopy() {
		
		return new FactorModelEffect(id, 
				Arrays.copyOf(levels, levels.length),
				Arrays.copyOf(levelCounts,levelCounts.length),
				size,
				numberOfLevels,
				restricted);
	}

	@Override
	public int getEffectSize() {
		if (restricted) return numberOfLevels - 1;
		else return numberOfLevels;
	}

}
