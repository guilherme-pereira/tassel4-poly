package net.maizegenetics.jGLiM.dm;

import java.util.ArrayList;

import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;

public class SweepFastNestedModel {
	double markerSS;
	double markerdf;
	double errorSS;
	double errordf;
	double taxaSS;
	double taxadf;
	double modelcfmSS;
	double modelcfmdf;
	double[] beta;
	DoubleMatrix G; //the inverse of X'X
	
	/**
	 * The constructor for this class. It will calculate the SS and df for marker after all other effects in the model except taxa nested in marker.
	 * Taxa in marker will then be added to the model and its SS and df calculated. Error SS and df will be the final residual.
	 * @param marker	the ModelEffect for the marker being tested
	 * @param taxaInMarker	the ModelEffect for taxa nested in marker
	 * @param otherEffects	any other effects in the model. The mean is always the first of these.
	 * @param data	the dependent variable
	 */
	public SweepFastNestedModel(ArrayList<ModelEffect> modelEffects, double[] data) {
		
		//order to add effects is other, marker, taxaInMarker
		int numberOfEffects = modelEffects.size();
		
		DoubleMatrix[][] xtx = new DoubleMatrix[numberOfEffects][numberOfEffects];
		DoubleMatrix[] xty = new DoubleMatrix[numberOfEffects];
		double yty;
		int[] dimX = new int[numberOfEffects];
		
		for (int i = 0; i < numberOfEffects; i++) {
			xtx[i][i] = modelEffects.get(i).getXtX();
			dimX[i] = xtx[i][i].numberOfColumns();
			xty[i] = modelEffects.get(i).getXty(data);
			for (int j = i + 1; j < numberOfEffects; j++) {
				xtx[i][j] = ModelEffectUtils.getXtY(modelEffects.get(i), modelEffects.get(j));
			}
		}
		
		
		yty = 0;
		int nobs = data.length;
		for (int i = 0; i < nobs; i++) {
			yty += data[i] * data[i];
		}
		
		SweepFast sf = new SweepFast(xtx, xty, yty);

		//assume that the first effect is the mean
		sf.revg2sweep(0);
		double residualAfterMean = sf.getResidualSS();
		sf.setDminFromA();
		
		//now sweep the rest of the other effects up to marker
		double dfother = 0;
		int col = 1;
		for (int i = 1; i < numberOfEffects - 2; i++) {
			for (int j = 0; j < dimX[i]; j++) if (sf.revg2sweep(col++)) dfother++;
		}
		
		double residualAfterOther = sf.getResidualSS();
		
		//sweep marker
		markerdf = 0;
		for (int i = 0; i < dimX[numberOfEffects - 2]; i++) {
			if (sf.revg2sweep(col++)) markerdf++;
		}
		double residualAfterMarker = sf.getResidualSS();
		
		//sweep taxa in marker
		taxadf = 0;
		for (int i = 0; i < dimX[numberOfEffects - 1]; i++) {
			if (sf.revg2sweep(col++)) taxadf++;
		}
		
		//calculate SS and df needed for output
		markerSS = residualAfterOther - residualAfterMarker;
		taxaSS = residualAfterMarker - sf.getResidualSS();
		modelcfmSS = residualAfterMean - residualAfterMarker;
		modelcfmdf = dfother + markerdf;
		errordf = data.length - 1 - modelcfmdf;
		errorSS = residualAfterMarker;
		
		beta = sf.getBeta();
		G = sf.getXTXpart();
	}
	
	public double[] getMarkerSSdf() {
		return new double[]{markerSS, markerdf};
	}
	
	public double[] getErrorSSdf() {
		return new double[] {errorSS, errordf};
	}
	
	public double[] getTaxaInMarkerSSdf() {
		return new double[] {taxaSS, taxadf};
	}
	
	public double[] getModelcfmSSdf() {
		return new double[] {modelcfmSS, modelcfmdf};
	}
	
	public double[] getBeta() {
		return beta;
	}
	
	public DoubleMatrix getInverseOfXtX() { return G; }
}
