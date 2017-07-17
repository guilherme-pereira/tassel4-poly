package net.maizegenetics.jGLiM;

import java.util.ArrayList;

import net.maizegenetics.pal.report.SimpleTableReport;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

public class LinearModelforStepwiseRegression {
    ArrayList<ModelEffect> modelEffects;
    int numberOfRequiredEffects = 1;
    DoubleMatrix2D[][] xtxmatrices;
    DoubleMatrix2D[][] xtymatrices;
    double[] data;
    double enterLimit = 1e-3;
    double exitLimit = 1e-3;
    
    LinearModelWithSweep lm;
    PartitionedLinearModel plm;
    
    public LinearModelforStepwiseRegression(ArrayList<ModelEffect> requiredEffects, double[] data) {
        modelEffects = requiredEffects;
        numberOfRequiredEffects = modelEffects.size();
        this.data = data;
        xtxmatrices = new DoubleMatrix2D[numberOfRequiredEffects][numberOfRequiredEffects];
        xtymatrices = new DoubleMatrix2D[numberOfRequiredEffects][1];
        for (int i = 0; i < numberOfRequiredEffects; i++) {
            ModelEffect me = requiredEffects.get(i);
            xtxmatrices[i][i] = me.getXTX();
            xtymatrices[i][0] = DoubleFactory2D.dense.make(me.getXTy(data).toArray(), me.getNumberOfLevels());
            for (int j = 0; j < i; j++) {
                xtxmatrices[i][j] = ModelEffect.getX1TX2(me, requiredEffects.get(j));
                xtxmatrices[j][i] = xtxmatrices[i][j].viewDice();
            }
        }
        
        lm = new LinearModelWithSweep(xtxmatrices, xtymatrices, data);
        plm = new PartitionedLinearModel(modelEffects, lm);
    }
    
    public void addEffect(ModelEffect me) {
        modelEffects.add(me);
        int newdim = xtxmatrices.length + 1;
        DoubleMatrix2D[][] oldxtx = xtxmatrices;
        DoubleMatrix2D[][] oldxty = xtymatrices;
        xtxmatrices = new DoubleMatrix2D[newdim][newdim];
        xtymatrices = new DoubleMatrix2D[newdim][1];
        for (int i = 0; i < newdim - 1; i++) {
            xtymatrices[i][0] = oldxty[i][0];
            for (int j = 0; j < newdim - 1; j++) {
                xtxmatrices[i][j] = oldxtx[i][j];
            }
        }
        xtxmatrices[newdim - 1][newdim - 1] = me.getXTX();
        xtymatrices[newdim - 1][0] = 
            DoubleFactory2D.dense.make(me.getXTy(data).toArray(), me.getNumberOfLevels());
        for (int i = 0; i < newdim - 1; i++) {
            xtxmatrices[i][newdim - 1] = ModelEffect.getX1TX2(modelEffects.get(i), me);
            xtxmatrices[newdim - 1][i] = xtxmatrices[i][newdim - 1].viewDice();
        }
        
        lm = new LinearModelWithSweep(xtxmatrices, xtymatrices, data);
        plm = new PartitionedLinearModel(modelEffects, lm);
    }
    
    public double[] testNewEffect(ModelEffect me) {
        //returns F, p-value for test
        //see if there is a term that can be added using a partitioned model
        
        plm.testNewModelEffect(me);
        double F = plm.getF();
        double p = plm.getp();
        return new double[]{F,p};
    }
    
    
    public ModelEffect backwardStep() {
        int numberOfModelEffects = modelEffects.size();
        if (numberOfModelEffects - numberOfRequiredEffects > 1) {
            double maxp = -1;
            double errordf = lm.getErrordf();
            double errorss = lm.getErrorSS();
            double errorms = errorss / errordf;
            int maxEffectnumber = -1;
            for (int i = numberOfRequiredEffects; i < numberOfModelEffects; i++) {
                double[] ssdf = lm.marginalEffectSSdf(i);
                double F = ssdf[0] / ssdf[1] / errorms;
                double p = -1;
                try { p = AbstractLinearModel.Ftest(F, ssdf[1], errordf); }
                catch(Exception e) {
                	System.err.println("Error calculating p value at effect = " + i);
                }
                if (p > maxp) {
                    maxp = p;
                    maxEffectnumber = i;
                }
            }

            if (maxp > exitLimit) {
                return removeTerm(maxEffectnumber);
            }
        }
        return null;
    }
    public ModelEffect removeTerm(int termNumber) {
        //remove the effect from the matrices
        int olddim = xtxmatrices.length;
        int newdim = olddim -1;
        DoubleMatrix2D[][] oldxtx = xtxmatrices;
        DoubleMatrix2D[][] oldxty = xtymatrices;
        xtxmatrices = new DoubleMatrix2D[newdim][newdim];
        xtymatrices = new DoubleMatrix2D[newdim][1];
        for (int i = 0; i < newdim; i++) {
            int ii = i;
            if (i >= termNumber) ii++;
            xtymatrices[i][0] = oldxty[ii][0];
            for (int j = 0; j < newdim; j++) {
                int jj = j;
                if (j >= termNumber) jj++;
                xtxmatrices[i][j] = oldxtx[ii][jj];
            }
        }
        //recalculate the model
        lm = new LinearModelWithSweep(xtxmatrices, xtymatrices, data);
        plm = new PartitionedLinearModel(modelEffects, lm);

        return modelEffects.remove(termNumber);
    }
    
    public double[] getyhat() {
        DoubleMatrix1D beta = lm.getBeta();
        double[] yhat = null;
        int numberOfEffects = modelEffects.size();
        int start = 0;
        for (int i = 0; i < numberOfEffects; i++) {
            ModelEffect me = modelEffects.get(i);
            DoubleMatrix1D effectyhat = me.getyhat(beta.viewPart(start, me.getNumberOfLevels()));
            if (i == 0) yhat = effectyhat.toArray();
            else for (int j = 0; j < yhat.length; j++) yhat[j] += effectyhat.getQuick(j);
            start += me.getNumberOfLevels();
        }
        return yhat;
    }
    
    public void changeData(double[] newdata) {
        data = newdata;
        int numberOfEffects = modelEffects.size();
        for (int i = 0; i < numberOfEffects; i++) {
            ModelEffect me = modelEffects.get(i);
            xtymatrices[i][0] = DoubleFactory2D.dense.make(me.getXTy(data).toArray(), me.getNumberOfLevels());
        }

        lm = new LinearModelWithSweep(xtxmatrices, xtymatrices, data);
        plm = new PartitionedLinearModel(modelEffects, lm);
    }
    
    public SimpleTableReport outputResults(String title, String traitname) {
        String[] heads = new String[]{"Trait","Term","SS","df", "MS", "F", "p", "Rsq"};
        int numberOfEffects = modelEffects.size();
        Object[][] results = new Object[numberOfEffects + 1][];
        double errordf = lm.getErrordf();
        double errorss = lm.getErrorSS();
        double totalss = lm.getTotalSS();
        double modelss = lm.getModelSS();
        double modeldf = lm.getModeldf();
        
        Object[] result;
        for (int i = 1; i < numberOfEffects; i++) {
            result = new Object[heads.length];
            result[0] = modelEffects.get(i).getId();
            double[] ssdf = lm.marginalEffectSSdf(i);
            int col = 0;
            result[col++] = traitname;
            result[col++] = modelEffects.get(i).getId();
            result[col++] = ssdf[0];
            result[col++] = ssdf[1];
            result[col++] = ssdf[0]/ssdf[1];
            double F = ssdf[0] / ssdf[1] / errorss * errordf;
            result[col++] = F;
            
            try {result[col++] = AbstractLinearModel.Ftest(F, ssdf[1], errordf);}
            catch(Exception e) {result[col++] = Double.NaN;}
            
            result[col++] = ssdf[0] / totalss;
            results[i-1] = result;
        }
        
        result = new Object[heads.length];
        int col = 0;
        result[col++] = traitname;
        result[col++] = "Model";
        result[col++] = modelss;
        result[col++] = modeldf;
        result[col++] = modelss / modeldf;
        double F = modelss / modeldf / errorss * errordf;
        result[col++] = F;
        
        try {result[col] = AbstractLinearModel.Ftest(F, modeldf, errordf);}
        catch(Exception e) {result[col] = Double.NaN;}
        col++;
        
        result[col++] = modelss / totalss;
        results[numberOfEffects - 1] = result;
        
        result = new Object[heads.length];
        col = 0;
        result[col++] = traitname;
        result[col++] = "Error";
        result[col++] = errorss;
        result[col++] = errordf;
        result[col++] = errorss / errordf;
        result[col++] = " ";
        result[col++] = " ";
        result[col++] = " ";
        results[numberOfEffects] = result;
        return new SimpleTableReport(title, heads, results);
    }
    
    public LinearModelWithSweep getLinearModel() {
        return lm;
    }

    public ArrayList<ModelEffect> getModelEffects() {
        return modelEffects;
    }

    public double getEnterLimit() {
        return enterLimit;
    }

    public void setEnterLimit(double enterLimit) {
        this.enterLimit = enterLimit;
    }

    public double getExitLimit() {
        return exitLimit;
    }

    public void setExitLimit(double exitLimit) {
        this.exitLimit = exitLimit;
    }
    
    
}
