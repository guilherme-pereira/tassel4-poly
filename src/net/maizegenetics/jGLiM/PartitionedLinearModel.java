package net.maizegenetics.jGLiM;

import java.util.ArrayList;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.jet.math.Functions;
import net.maizegenetics.jGLiM.ModelEffect;

public class PartitionedLinearModel {
    ArrayList<ModelEffect> modelEffectList;     //the model effects in the first partition
    LinearModelWithSweep lm = null;                    //the linear model representing the first partition
    DoubleMatrix2D G1;                          //the inverse for the first partition
    DoubleMatrix2D X1Ty;                        //XTy matrix for the first partition
    double modelss;
    double errorss;
    double modeldf;
    double errordf;
    
    public PartitionedLinearModel(ArrayList<ModelEffect> baseModel, LinearModelWithSweep lm) {
        modelEffectList = baseModel;
        this.lm = lm;
        X1Ty = DoubleFactory2D.dense.compose(lm.xtymatrices);
    }
    
    public LinearModelWithSweep getLinearModel() {return lm;}
    
    public void testNewModelEffect(ModelEffect me) {
        int numberOfTerms = modelEffectList.size();
        DoubleMatrix2D[][] x2tx1matrices = new DoubleMatrix2D[1][numberOfTerms];
        for (int t = 0; t < numberOfTerms; t++) {
            x2tx1matrices[0][t] = ModelEffect.getX1TX2(me, modelEffectList.get(t));
        }
        
        //from Searle. 1987. Linear Models for Unbalanced Data. p.264
        //for the added term, ssmodel = y'M1X2BX2'M1y, where
        //B = inverse(X2'M1X2)
        //X2'M1X2 = X2'X2 - X2'X1G1X1'X2
        //use A = X2'X1G1 and X2'M1X2 = AX1'X2
        //X2'M1y = X2'y - X2'X1G1X1'y = X2Ty - AX1Ty
        
        DoubleMatrix2D G1 = lm.getInverse();
        DoubleMatrix2D X2TX1 = DoubleFactory2D.dense.compose(x2tx1matrices);
        DoubleMatrix2D X2TX2 = me.getXTX();
        DoubleMatrix2D A = X2TX1.zMult(G1,null);
        DoubleMatrix2D X2TM1X2 = X2TX2.copy().assign(A.zMult(X2TX1.viewDice(), null), Functions.minus);
        
        int[] rank = new int[]{0};
        DoubleMatrix2D B = AbstractLinearModel.geninv(X2TM1X2, rank);
        
        DoubleMatrix1D X2Ty = me.getXTy(lm.data);
        DoubleMatrix1D X2TM1y = X2Ty.assign(A.zMult(X1Ty, null).viewColumn(0), Functions.minus);
        modelss = X2TM1y.zDotProduct(B.zMult(X2TM1y, null)); 
        modeldf = rank[0];
        errorss = lm.getErrorSS() - modelss;
        errordf = lm.getErrordf() - modeldf;
    }
    
    public double getModelSS() {return modelss;}
    public double getModeldf() {return modeldf;}
    public double getErrorSS() {return errorss;}
    public double getErrordf() {return errordf;}

    public double getF() {
        return modelss / errorss / modeldf * errordf;  
    }

    public double getp() {
        double p = 1;
        try{
            p = AbstractLinearModel.Ftest(getF(), modeldf, errordf);
        } catch(Exception e) {
//            System.err.println("error calculating p");
        }
        return p;
    }

}
