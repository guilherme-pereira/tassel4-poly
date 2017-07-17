package net.maizegenetics.jGLiM;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.jet.math.Functions;

public class CovariateModelEffect extends ModelEffect {
    double[] covariate;
    
    public CovariateModelEffect(double[] covariate) {
        super();
        this.covariate = covariate;
        size = covariate.length;
        numberOfLevels = 1;
    }

    public double getSumSquares() {
        double ss = 0;
        for (int i = 0; i < size; i++) {
            ss += covariate[i] * covariate[i];
        }
        return ss;
    }
    
    public double getSumProducts(double[] other) {
        double sp = 0;
        for (int i = 0; i < size; i++) {
            sp += covariate[i] * other[i];
        }
        return sp;
    }
    
    @Override
    public DoubleMatrix2D getXTCov(double[][] covariates) {
        int n = covariates.length;
        DoubleMatrix2D XTCov = DoubleFactory2D.dense.make(1,n,0);
        for (int i = 0; i < n; i++) {
            XTCov.setQuick(0,n,getSumProducts(covariates[i]));
        }
        
        return XTCov;
    }

    @Override
    public DoubleMatrix2D getXTX() {
        return DoubleFactory2D.dense.make(1,1,getSumSquares());
    }

    @Override
    public DoubleMatrix1D getXTy(double[] y) {
        return DoubleFactory1D.dense.make(1,getSumProducts(y));
    }

    @Override
    public DoubleMatrix1D getyhat(DoubleMatrix1D beta) {
        double b = beta.get(0);
        return DoubleFactory1D.dense.make(covariate).assign(Functions.mult(b));
    }

    public double[] getCovariate() {
        return covariate;
    }
    
    
    
}
