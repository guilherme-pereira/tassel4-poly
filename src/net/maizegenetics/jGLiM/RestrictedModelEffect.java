package net.maizegenetics.jGLiM;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

public class RestrictedModelEffect extends ModelEffect {
    
    public RestrictedModelEffect(int[] values, int numberOfLevels) {
        super(values, numberOfLevels);
        restrictLevels();
    }
    public RestrictedModelEffect(int[] values) {
        super(values);
        restrictLevels();
    }

    public void restrictLevels() {
        if (numberOfLevels > 1) {
            restrictedLevels = new int[numberOfLevels - 1];
            for (int lev = 0; lev < numberOfLevels - 1; lev++) restrictedLevels[lev] = lev;
        }
    }

    @Override
    public DoubleMatrix2D getXTCov(double[][] covariates) {
        return super.getXTCov(covariates).viewSelection(restrictedLevels, null).copy();
    }
    @Override
    public DoubleMatrix2D getXTX() {
        return super.getXTX().viewSelection(restrictedLevels, restrictedLevels).copy();
    }
    @Override
    public DoubleMatrix1D getXTy(double[] y) {
        return super.getXTy(y).viewSelection(restrictedLevels).copy();
    }
    @Override
    public int getNumberOfLevels() {
        return numberOfLevels - 1;
    }
	@Override
	public DoubleMatrix2D getX() {
		return super.getX().viewSelection(null,restrictedLevels).copy();
	}
    

}
