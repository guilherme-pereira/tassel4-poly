package net.maizegenetics.jGLiM;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleFactory1D;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

public class ModelEffect {
    protected Object id;
    protected int[] values = null;
    protected int numberOfLevels = 0;
    protected int size = 0;
    protected boolean isDiscrete = true;
    protected int[] restrictedLevels = null;
    
    public ModelEffect(int[] values) {
        this.values = values;
        size = values.length;
        numberOfLevels = 0;
        for (int i = 0; i < size; i++) if (values[i] > numberOfLevels) numberOfLevels = values[i];
        numberOfLevels++;
    }

    public ModelEffect(int[] values, int numberOfLevels) {
        this.values = values;
        size = values.length;
        this.numberOfLevels = numberOfLevels;
    }

    public ModelEffect(){}

    public int[] getIntegerLevelCounts(){
        int[] count = new int[numberOfLevels];
        for (int i = 0; i < numberOfLevels; i++) count[i] = 0;
        for (int i = 0; i < size; i++) count[values[i]]++;
        return count;
    }

    public double[] getDoubleLevelCounts() {
        double[] count = new double[numberOfLevels];
        for (int i = 0; i < numberOfLevels; i++) count[i] = 0;
        for (int i = 0; i < size; i++) count[values[i]]++;
        return count;
    }

    public int[][] getIntegerLevelCounts(ModelEffect me) {
        int[][] count = new int[numberOfLevels][me.numberOfLevels];
        for (int i = 0; i < numberOfLevels; i++) {
            for (int j = 0; j < me.numberOfLevels; j++) {
                count[i][j] =0;
            }
        }
        for (int i = 0; i < size; i++) {
            count[values[i]][me.values[i]]++;
        }
        return count;
    }

    public double[][] getDoubleLevelCounts(ModelEffect me) {
        double[][] count = new double[numberOfLevels][me.numberOfLevels];
        for (int i = 0; i < numberOfLevels; i++) {
            for (int j = 0; j < me.numberOfLevels; j++) {
                count[i][j] =0;
            }
        }
        for (int i = 0; i < size; i++) {
            count[values[i]][me.values[i]]++;
        }
        return count;

    }

    public double[] getLevelSums(double[] y) {
        if (y.length != size) return null;
        double sum[] = new double[numberOfLevels];
        for (int i = 0; i < numberOfLevels; i++) sum[i] =0;
        for (int i = 0; i < size; i++) sum[values[i]] += y[i];
        return sum;
    }

    public DoubleMatrix2D getXTX() {
        return DoubleFactory2D.dense.diagonal(DoubleFactory1D.dense.make(getDoubleLevelCounts()));
    }

    public DoubleMatrix2D getX1TX2(ModelEffect me) {
        if (me instanceof RestrictedModelEffect) {
            return DoubleFactory2D.dense.make(getDoubleLevelCounts(me)).viewSelection(null, me.restrictedLevels);
        }
        else if (me instanceof CovariateModelEffect) {
            CovariateModelEffect cme = (CovariateModelEffect) me;
            return DoubleFactory2D.dense.make(getLevelSums(cme.covariate), numberOfLevels);
        }
        return DoubleFactory2D.dense.make(getDoubleLevelCounts(me)); 
    }

    public DoubleMatrix1D getXTy(double[] y) {
        return DoubleFactory1D.dense.make(getLevelSums(y));
    }

    //the higher dimension of the array is the number of covariates
    //the lower dimension is the number of data points, assumed equal to the number of values in this effect
    public DoubleMatrix2D getXTCov(double[][] covariates) {
        if (covariates[0].length != values.length) return null;
        DoubleMatrix2D xtc = DoubleFactory2D.dense.make(this.getNumberOfLevels(), covariates.length);
        for (int i = 0; i < covariates.length; i++) {
            xtc.viewColumn(i).assign(getXTy(covariates[i]));
        }
        return xtc;
    }

    public DoubleMatrix2D getX() {
    	DoubleMatrix2D X = DoubleFactory2D.dense.make(size, numberOfLevels, 0);
    	for (int r = 0; r < size; r++) X.setQuick(r, values[r], 1);
    	return X;
    }
    
    public DoubleMatrix1D getyhat(DoubleMatrix1D beta) {
        int nlevels = beta.size();
        DoubleMatrix1D yhat = DoubleFactory1D.dense.make(size,0);
        for (int i = 0; i < size; i++) {
            if (values[i] < nlevels) {
                yhat.setQuick(i, beta.getQuick (values[i]));
            }
        }
        return yhat;
    }

    public int[] getRestrictedLevels() {
        return restrictedLevels;
    }

    public int getNumberOfLevels() {
        return numberOfLevels;
    }
    
    public int getValue(int ndx) {return values[ndx];}

    public static DoubleMatrix2D getCovTCov(double[][] covariates) {
        DoubleMatrix2D cov = DoubleFactory2D.dense.make(covariates);
        return cov.viewDice().zMult(cov,null);
    }
    
    public static DoubleMatrix1D getyhat(DoubleMatrix1D beta, double[][] covariates) {
        int size = covariates[0].length;
        DoubleMatrix1D yhat = DoubleFactory1D.dense.make(size,0);
        for (int i = 0; i < covariates.length; i++) {
            double b = beta.getQuick(i);
            for (int j = 0; j < size; j++) {
                yhat.setQuick(j, yhat.getQuick(j) + b * covariates[i][j]);
            }
        }
        return yhat;
    }
    
    public static DoubleMatrix2D getX1TX2(ModelEffect me1, ModelEffect me2) {
        if (me1 instanceof CovariateModelEffect) {
            CovariateModelEffect cme = (CovariateModelEffect) me1;
            DoubleMatrix1D xty =  me2.getXTy(cme.getCovariate());
            DoubleMatrix2D x1tx2 = DoubleFactory2D.dense.make(1, xty.size());
            for (int i = 0; i < xty.size(); i++) x1tx2.setQuick(0,i,xty.getQuick(i));
            return x1tx2;
        }
        else if (me2 instanceof CovariateModelEffect) {
            CovariateModelEffect cme = (CovariateModelEffect) me2;
            DoubleMatrix1D xty =  me1.getXTy(cme.getCovariate());
            DoubleMatrix2D x1tx2 = DoubleFactory2D.dense.make(xty.size(),1);
            for (int i = 0; i < xty.size(); i++) x1tx2.setQuick(i,0,xty.getQuick(i));
            return x1tx2;
        }
        else if (me1 instanceof NestedCovariateModelEffect) {
            return me1.getX1TX2(me2);
        }
        else if (me2 instanceof NestedCovariateModelEffect) {
            return me2.getX1TX2(me1).viewDice().copy();
        }
        else if (me1 instanceof RestrictedModelEffect || me2 instanceof RestrictedModelEffect) {
            return DoubleFactory2D.dense.make(me1.getDoubleLevelCounts(me2)).viewSelection(me1.restrictedLevels, me2.restrictedLevels).copy();
        }
        
        return DoubleFactory2D.dense.make(me1.getDoubleLevelCounts(me2));
    }

    public static int[] getIntegerLevels(ArrayList alist) {
        ArrayList theLevels = new ArrayList();
        int[] intLevels = new int[alist.size()];
        for (int i = 0; i <  alist.size(); i++) {
            int ndx = theLevels.indexOf(alist.get(i));
            if (ndx == -1) {
                intLevels[i] = theLevels.size();
                theLevels.add(alist.get(i));
            }
            else intLevels[i] = ndx;
        }
        return intLevels;
    }

    public static int[] getIntegerLevels(Object[] originalLevels) {
    	int nLevels = originalLevels.length;
    	int[] intLevels = new int[nLevels];
    	HashMap<Object, Integer> levelMap = new HashMap<Object, Integer>();
    	for (int i = 0; i < nLevels; i++) {
    		Integer ndx = levelMap.get(originalLevels[i]);
    		if (ndx == null) {
    			ndx = new Integer(levelMap.size());
    			levelMap.put(originalLevels[i], ndx);
    		}
    		intLevels[i] = ndx.intValue();
    	}
    	return intLevels;
    }
    
    public static int[] getIntegerLevels(Object[] originalLevels, ArrayList<Object> levelIds) {
    	int nLevels = originalLevels.length;
    	int[] intLevels = new int[nLevels];
    	HashMap<Object, Integer> levelMap = new HashMap<Object, Integer>();
    	for (int i = 0; i < nLevels; i++) {
    		Integer ndx = levelMap.get(originalLevels[i]);
    		if (ndx == null) {
    			ndx = new Integer(levelMap.size());
    			levelMap.put(originalLevels[i], ndx);
    			if (levelIds != null) levelIds.add(originalLevels[i]);
    		}
    		intLevels[i] = ndx.intValue();
    	}
    	return intLevels;
    }
    
    public static <T> int[] getIntegerLevels(ArrayList<T> originalLevels, ArrayList<T> levelIds) {
    	int nLevels = originalLevels.size();
    	int[] intLevels = new int[nLevels];
    	HashMap<Object, Integer> levelMap = new HashMap<Object, Integer>();
    	for (int i = 0; i < nLevels; i++) {
    		Integer ndx = levelMap.get(originalLevels.get(i));
    		if (ndx == null) {
    			ndx = new Integer(levelMap.size());
    			levelMap.put(originalLevels.get(i), ndx);
    			if (levelIds != null) levelIds.add(originalLevels.get(i));
    		}
    		intLevels[i] = ndx.intValue();
    	}
    	return intLevels;
    }
    
    public Object getId() {
        return id;
    }

    public void setId(Object id) {
        this.id = id;
    }
    
    public int getSize() {return size;}
}
