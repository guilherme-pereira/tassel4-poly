package net.maizegenetics.jGLiM;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

public class NestedCovariateModelEffect extends ModelEffect {
    CovariateModelEffect cme;
    ModelEffect withinME;
    
    public NestedCovariateModelEffect(CovariateModelEffect cme, ModelEffect me) {
        this.cme = cme;
        withinME = me;
    }

    @Override
    public int getNumberOfLevels() {
        return withinME.numberOfLevels;
    }

    public DoubleMatrix2D getX1TX2nested(NestedCovariateModelEffect me) {
        DoubleMatrix2D x1tx2 = DoubleFactory2D.dense.make(getNumberOfLevels(), me.getNumberOfLevels(), 0);
        int[] levels1 = withinME.values;
        int[] levels2 = me.withinME.values;
        
        for (int i = 0; i < withinME.size; i++) {
            double prod = cme.covariate[i] * me.cme.covariate[i];
            int lev1 = levels1[i];
            int lev2 = levels2[i];
            x1tx2.setQuick(lev1, lev2, x1tx2.getQuick(lev1,lev2) + prod);
        }
        
        return x1tx2;
    }
    
    public DoubleMatrix2D getX1TX2covariate(CovariateModelEffect othercme) {
        DoubleMatrix2D x1tx2 = DoubleFactory2D.dense.make(getNumberOfLevels(), 1, 0);
        int[] levels = withinME.values;
        
        for (int i = 0; i < withinME.size; i++) {
            double prod = cme.covariate[i] * othercme.covariate[i];
            int lev = levels[i];
            x1tx2.setQuick(lev, 0, x1tx2.getQuick(lev,0) + prod);
        }
        
        return x1tx2;
    }
    
    public DoubleMatrix2D getX1TX2modeleffect(ModelEffect me) {
        //for both regular and restricted model effects
        int melevels = me.getNumberOfLevels();
        DoubleMatrix2D x1tx2 = DoubleFactory2D.dense.make(getNumberOfLevels(), melevels, 0);
        int[] levels1 = withinME.values;
        int[] levels2 = me.values;
        
        for(int i = 0; i < withinME.size; i++) {
            int lev1 = levels1[i];
            int lev2 = levels2[i];
            if (lev2 < melevels)  
                x1tx2.setQuick(lev1,lev2,x1tx2.getQuick(lev1,lev2) + cme.covariate[i]);
        }
        return x1tx2;
    }
    
    public DoubleMatrix2D getX1TX2(ModelEffect me) {
        DoubleMatrix2D x1tx2 = null;
        if (me instanceof NestedCovariateModelEffect) {
            return getX1TX2nested((NestedCovariateModelEffect) me);
        }
        else if (me instanceof CovariateModelEffect) {
            return getX1TX2covariate((CovariateModelEffect) me);
        }
        else return getX1TX2modeleffect(me);
    }
    
    @Override
    public DoubleMatrix2D getXTX() {
        DoubleMatrix2D xtx = DoubleFactory2D.dense.make(withinME.numberOfLevels, withinME.numberOfLevels, 0);
        int[] levels = withinME.values;
        
        for (int i = 0; i < withinME.size; i++) {
            double cov = cme.covariate[i];
            int lev = levels[i];
            xtx.setQuick(lev, lev, xtx.getQuick(lev,lev) + cov*cov);
        }
        
        return xtx;
    }

    @Override
    public DoubleMatrix1D getXTy(double[] y) {
        DoubleMatrix1D xty = DoubleFactory1D.dense.make(withinME.numberOfLevels, 0);
        int[] levels = withinME.values;
        
        for (int i = 0; i < withinME.size; i++) {
            double prod = cme.covariate[i] * y[i];
            int lev = levels[i];
            xty.setQuick(lev, xty.getQuick(lev) + prod);
        }
        
        return xty;
    }

    @Override
    public DoubleMatrix2D getXTCov(double[][] covariates) {
        int rows = withinME.numberOfLevels;
        int cols = covariates.length;
        DoubleMatrix2D xtc = DoubleFactory2D.dense.make(rows, cols, 0);
        int[] levels = withinME.values;
        for (int i = 0; i < withinME.size; i++) {
            int lev = levels[i];
            for (int c = 0; c < cols; c++) {
                double prod = cme.covariate[i] * covariates[c][i];
                xtc.setQuick(lev, c, xtc.getQuick(lev,c) + prod);
            }
        }
        
        return xtc;
    }

    @Override
    public DoubleMatrix1D getyhat(DoubleMatrix1D beta) {
        DoubleMatrix1D yhat = DoubleFactory1D.dense.make(withinME.size, 0);
        int[] levels = withinME.values;
        for (int i = 0; i < withinME.size; i++) {
            yhat.setQuick(i, beta.getQuick(levels[i])*cme.covariate[i]);
        }
        
        return yhat;
    }
    
}
