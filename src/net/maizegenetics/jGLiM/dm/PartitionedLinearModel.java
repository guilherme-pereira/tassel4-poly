package net.maizegenetics.jGLiM.dm;

import java.util.ArrayList;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.jet.math.Functions;

import net.maizegenetics.jGLiM.AbstractLinearModel;
import net.maizegenetics.jGLiM.LinearModelUtils;
import net.maizegenetics.jGLiM.dm.SweepFastLinearModel;
import net.maizegenetics.jGLiM.dm.ModelEffect;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;


public class PartitionedLinearModel {
	//	 ArrayList<ModelEffect> modelEffectList;     //the model effects in the first partition
	SweepFastLinearModel lm = null;                    //the linear model representing the first partition
	ArrayList<ModelEffect> baseModel;
	DoubleMatrix G1;                          //the inverse for the first partition
	DoubleMatrix[][] x2tx1matrices;
	DoubleMatrix X1Ty;
	double modelss;
	double errorss;
	double modeldf;
	double errordf;
	double[] y;

	public PartitionedLinearModel(ArrayList<ModelEffect> baseModel, SweepFastLinearModel lm) {
		this.lm = lm;
		y = lm.y;
		this.baseModel = baseModel;
		G1 = lm.getInverseOfXtX();
		x2tx1matrices = new DoubleMatrix[1][baseModel.size()];
		int nx1 = baseModel.size();
		DoubleMatrix[][] x1tymatrices = new DoubleMatrix[nx1][1];
		for (int i = 0; i < nx1; i++) {
			x1tymatrices[i][0] = baseModel.get(i).getXty(y);
		}
		X1Ty = DoubleMatrixFactory.DEFAULT.compose(x1tymatrices);
	}

	public SweepFastLinearModel getLinearModel() {return lm;}

	public void testNewModelEffect(ModelEffect me) {
		int nx1 = baseModel.size();
		for (int i = 0; i < nx1; i++) {
			x2tx1matrices[0][i] = ModelEffectUtils.getXtY(me, baseModel.get(i));
		}
        
        //from Searle. 1987. Linear Models for Unbalanced Data. p.264
        //for the added term, ssmodel = y'M1X2BX2'M1y, where
        //B = inverse(X2'M1X2)
        //X2'M1X2 = X2'X2 - X2'X1G1X1'X2
        //use A = X2'X1G1 and X2'M1X2 = X2'X2 - AX1'X2
        //X2'M1y = X2'y - X2'X1G1X1'y = X2Ty - AX1Ty
        

		DoubleMatrix X2TX1 = DoubleMatrixFactory.DEFAULT.compose(x2tx1matrices);
		DoubleMatrix X2TX2 = me.getXtX();
		DoubleMatrix A = X2TX1.mult(G1);
		DoubleMatrix X2TM1X2 = X2TX2.minus(A.tcrossproduct(X2TX1));
		double[] ssdf = lm.getResidualSSdf();
		DoubleMatrix X2Ty = me.getXty(y);
		DoubleMatrix X2TM1y = X2Ty.minus(A.mult(X1Ty));
		
		int[] rank = new int[]{0};
		DoubleMatrix B = X2TM1X2.generalizedInverseWithRank(rank);
		modelss = X2TM1y.crossproduct(B).mult(X2TM1y).get(0,0);
		errorss = ssdf[0] - modelss;
		modeldf = rank[0];
		errordf = ssdf[1] - modeldf;
	}

	public double testNewModelEffect(double[] covariate) {
		//from Searle. 1987. Linear Models for Unbalanced Data. p.264
		//for the added term, ssmodel = y'M1X2BX2'M1y, where
		//B = inverse(X2'M1X2)
		//X2'M1X2 = X2'X2 - X2'X1G1X1'X2
		//use A = X2'X1G1 and X2'X1G1X1'X2 = AX1'X2
		//X2'M1y = X2'y - X2'X1G1X1'y = X2Ty - AX1Ty
		
		int nx1 = baseModel.size();
		for (int i = 0; i < nx1; i++) {
			x2tx1matrices[0][i] = baseModel.get(i).getXty(covariate).transpose();
		}
		DoubleMatrix X2TX1 = DoubleMatrixFactory.DEFAULT.compose(x2tx1matrices);
		DoubleMatrix A = X2TX1.mult(G1);
		double ax1ty = A.mult(X1Ty).get(0, 0);
		double ax1tx2 = A.tcrossproduct(X2TX1).get(0, 0);
		
		double x2tx2 = 0;
		double x2ty = 0;
		int count = 0;
		for (double d:covariate) {
			x2tx2 += d * d;
			x2ty += d * y[count++];
		}
		
		double x2tm1x2 = x2tx2 - ax1tx2;
		double x2tm1y = x2ty - ax1ty;
		if (x2tm1x2 <1e-12) {
			return 0;
		}
		return x2tm1y * x2tm1y / x2tm1x2;
	}

	//only appropriate for a covariate which has one df
	public void setModelSS(double modelss) {
		this.modelss = modelss;
		double[] ssdf = lm.getResidualSSdf();
		errorss = ssdf[0] - modelss;
		modeldf = 1;
		errordf = ssdf[1] - modeldf;
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
		if (modeldf == 0) return 1;
		try{
			p = LinearModelUtils.Ftest(getF(), modeldf, errordf);
		} catch(Exception e) {
			//	            System.err.println("error calculating p");
		}
		return p;
	}
	
	public double[] getFp() {
		double F = modelss / errorss / modeldf * errordf;
		double p = Double.NaN;
		if (modeldf > 0) {
			try{
				p = LinearModelUtils.Ftest(getF(), modeldf, errordf);
			} catch(Exception e) {
				//	            System.err.println("error calculating p");
			}
		}
		return new double[]{F, p};
		
	}
}
