package net.maizegenetics.jGLiM.dm;

import java.util.ArrayList;
import java.util.Arrays;

import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;

public class SweepFastLinearModel {
	SweepFast sf;
	ArrayList<ModelEffect> modelEffects;
	int[] dimX;
	double[] y;
	ArrayList<Double> incrementalss;
	ArrayList<Double> incrementaldf;
	double totalSS;
	double modeldf;
	
	public SweepFastLinearModel(ArrayList<ModelEffect> effects, double[] data) {
		modelEffects = effects;
		y = data;
		init();
	}
	
	public SweepFastLinearModel(ArrayList<ModelEffect> effects, DoubleMatrix[][] xtx, DoubleMatrix[] xty, double[] data) {
		modelEffects = effects;
		y = data;
		init(xtx, xty);
	}
	
	private void init() {
    	int neffects = modelEffects.size();
    	DoubleMatrix[][] xtx = new DoubleMatrix[neffects][neffects];
    	DoubleMatrix[] xty = new DoubleMatrix[neffects];
    	
    	for (int i = 0; i < neffects; i++) {
    		xty[i] = modelEffects.get(i).getXty(y);
    		xtx[i][i] = modelEffects.get(i).getXtX();
    		for (int j = i + 1; j < neffects; j++) {
    			xtx[i][j] = ModelEffectUtils.getXtY(modelEffects.get(i), modelEffects.get(j));
    		}
    	}
    	
    	init(xtx, xty);
	}
	
	private void init(DoubleMatrix[][] xtx, DoubleMatrix[] xty) {
		incrementalss = new ArrayList<Double>();
		incrementaldf = new ArrayList<Double>();
		
    	int neffects = xtx.length;
    	dimX = new int[neffects];
    	for (int i = 0; i < neffects; i++) {
    		dimX[i] = xtx[i][i].numberOfColumns();
    	}
    	
    	double sum = 0;
    	double sumsq = 0;
    	for (double d:y) {
    		sum += d;
    		sumsq += d * d;
    	}
    	totalSS = sumsq;
    	sf = new SweepFast(xtx, xty, totalSS);
    	
    	modeldf = 0;
    	double thisdf = 0;
    	int count = 0;
    	for (int i = 0; i < dimX[0]; i++) {
    		if (sf.revg2sweep(count++)) {
    			thisdf++;
    			modeldf++;
    		}
    	}
    	
    	double previousErrorSS = sf.getResidualSS();
    	incrementalss.add(new Double(totalSS - previousErrorSS));
    	incrementaldf.add(new Double(thisdf));
    	sf.setDminFromA();
    	
    	for (int i = 1; i < neffects; i++) {
    		thisdf = 0;
    		for (int j = 0; j < dimX[i]; j++) {
    			if (sf.revg2sweep(count++)) {
    		    	thisdf++;
    				modeldf++;
    			}
    		}
    		double newErrorSS = sf.getResidualSS();
    		incrementalss.add(new Double(previousErrorSS - newErrorSS));
    		incrementaldf.add(new Double(thisdf));
    		previousErrorSS = newErrorSS;
    	}
	}

	public double[] getMarginalSSdf(int effect) {
		int start = 0;
		for (int i = 0; i < effect; i++) start += dimX[i];
		int end = start + dimX[effect];
		
		double df = 0;
		
		for (int i = start; i < end; i++) sf.revg2sweep(i);
		sf.sweepSingularColumns();
		
		double reducedModelSS = sf.getResidualSS();
		for (int i = start; i < end; i++) {
			if (sf.revg2sweep(i)) df++;
		}

		return new double[]{reducedModelSS - sf.getResidualSS(), df};
	}
	
	public double[] getIncrementalSSdf(int effect) {
		return new double[]{incrementalss.get(effect), incrementaldf.get(effect)};
	}
	
	public double[] getResidualSSdf() {
		double ss = sf.getResidualSS();
		double df = y.length - modeldf;
		return new double[]{ss, df};
	}
	
	public double[] getFullModelSSdf() {
		double ss = totalSS - sf.getResidualSS(); 
		double df = 0;
		for (Double effectdf:incrementaldf) df += effectdf;
		return new double[]{ss, df};
	}
	
	public double[] getModelcfmSSdf() {
		double ss = totalSS - incrementalss.get(0) - sf.getResidualSS(); 
		double df = 0;
		int n = incrementaldf.size();
		for (int i = 1; i < n; i++) df += incrementaldf.get(i);
		return new double[]{ss, df};
	}
	
	public double[] getBeta() {
		return sf.getBeta();
	}
	
	public double[] getBlues(int effect, boolean restricted) {
		//only valid for factor model effects
		if (modelEffects == null) return null;
		if (!(modelEffects.get(effect) instanceof FactorModelEffect)) return null;
		double[] blues;
		double mean = 0;
		int neffects = modelEffects.size();
		double[] beta = getBeta();
		
		int startIndex = 0;
		int startIndexForEffect = 0;
		for (int i = 0; i < neffects; i++) {
			double effectAverage = 0;
			if (i != effect && (modelEffects.get(effect) instanceof FactorModelEffect) ) {
				FactorModelEffect fme = (FactorModelEffect) modelEffects.get(effect);
				for (int j = 0; j < dimX[i]; j++) {
					effectAverage += beta[startIndex + j];
				}
				if (restricted && effect > 0) effectAverage /= dimX[i] + 1; //effect 0 is the mean and is always not restricted
				else effectAverage /= dimX[i];
			} else if (i == effect) {
				startIndexForEffect = startIndex;
			}
			startIndex += dimX[effect];
			mean += effectAverage;
		}
		
		if (restricted) {
			blues = new double[dimX[effect] + 1];
		} else {
			blues = new double[dimX[effect]];
		}
		
		for (int i = 0; i < dimX[effect]; i++) {
			blues[i] = mean + beta[startIndexForEffect + i];
		}
		
		if (restricted) blues[blues.length - 1] = mean;
		
		return blues;
	}
	
	public DoubleMatrix getInverseOfXtX() {
		DoubleMatrix A = sf.getA();
		int n = A.numberOfRows();
		int[] ndx = new int[n - 1];
		for (int i = 0; i < n - 1; i++) ndx[i] = i;
		return A.getSelection(ndx, ndx);
	}

	public DoubleMatrix getPredictedValues() {
	
		if (modelEffects == null) return null;
		int numberOfEffects = modelEffects.size();
		if (numberOfEffects == 0) return null;
		double[] beta = getBeta();
		int start = 0;
		ModelEffect me = modelEffects.get(0);
		int effectSize = me.getEffectSize();
		double[] thisbeta = Arrays.copyOfRange(beta, start, start + effectSize);
		DoubleMatrix p = me.getyhat(thisbeta);
		start += effectSize;
		for (int i = 1; i < numberOfEffects; i++) {
			me = modelEffects.get(i);
			effectSize = me.getEffectSize();
			thisbeta = Arrays.copyOfRange(beta, start, start + effectSize);
			p.plusEquals(me.getyhat(thisbeta));
			start += effectSize;
		}
		return p;
	}
	
	public DoubleMatrix getResiduals() {
		DoubleMatrix r = DoubleMatrixFactory.DEFAULT.make(y.length, 1, y);
		r.minusEquals(getPredictedValues());
		return r;
	}
}
