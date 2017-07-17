package net.maizegenetics.stats.PCA;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import cern.colt.matrix.linalg.SingularValueDecomposition;
import cern.jet.math.Functions;


public class PrincipalComponents {
    private DoubleMatrix2D genotypes;
    private SingularValueDecomposition svd;
    private boolean transposed = false;
    
    public static final int TYPE_CENTER = 1;
    public static final int TYPE_CENTER_AND_SCALE = 2;
    public static final int TYPE_NONE = 3;
    
    
    public PrincipalComponents(DoubleMatrix2D genotypes, int type) {
        this.genotypes = genotypes;
        if (type == TYPE_CENTER) centerGenotypesByMarker();
        else if (type == TYPE_CENTER_AND_SCALE) centerAndScaleGenotypes();
        int nrows = genotypes.rows();
        int ncols = genotypes.columns();
        if(nrows < ncols) {
            transposed = true;
            svd = new SingularValueDecomposition(genotypes.viewDice());
        }
        else svd = new SingularValueDecomposition(genotypes);
        
        System.out.println("finished");
    }
    
    public void centerGenotypesByMarker() {
        int nrows = genotypes.rows();
        int ncols = genotypes.columns();
        
        for (int c = 0; c < ncols; c++) {
            DoubleMatrix1D marker = genotypes.viewColumn(c);
            double mean = marker.zSum() / nrows;
            marker.assign(Functions.minus(mean));
        }
    }
    
    public void centerAndScaleGenotypes() {
        int nrows = genotypes.rows();
        int ncols = genotypes.columns();
        
        for (int c = 0; c < ncols; c++) {
            DoubleMatrix1D marker = genotypes.viewColumn(c);
            double mean = marker.aggregate(Functions.plus, Functions.identity) / nrows;
            marker.assign(Functions.minus(mean));
            double sd = marker.aggregate(Functions.plus, Functions.square);
            sd /= marker.size() - 1;
            sd = Math.sqrt(sd);
            marker.assign(Functions.div(sd));
        }
    }
    
    public DoubleMatrix2D getUS(int n) {
        if (transposed) {
            int ncol = svd.getV().columns();
            return svd.getV().viewPart(0, 0, ncol, n).zMult(svd.getS().viewPart(0, 0, n, n), null);
        }
        else {
            int ncol = svd.getU().columns();
            return svd.getU().viewPart(0, 0, ncol, n).zMult(svd.getS().viewPart(0, 0, n, n), null);
        }
    }

    public DoubleMatrix2D getSV(int n) {
        if (transposed) {
            int ncol = svd.getU().columns();
            return svd.getU().viewPart(0, 0, ncol, n).zMult(svd.getS().viewPart(0, 0, n, n), null).viewDice().copy();
        }
        else {
            int ncol = svd.getV().columns();
            return svd.getV().viewPart(0, 0, ncol, n).zMult(svd.getS().viewPart(0, 0, n, n), null).viewDice().copy();
        }
    }
    
    public DoubleMatrix2D getU() {
    	if (transposed) return svd.getV();
    	else return svd.getU();
    }
    
    public DoubleMatrix2D getV() {
    	if (transposed) return svd.getU();
    	else return svd.getV();
    }
    
    public DoubleMatrix1D getSingularValues() {
    	return DoubleFactory1D.dense.make(svd.getSingularValues());
    }
    
    public DoubleMatrix2D getEigenVectors() {
    	if (transposed) return svd.getU();
    	return svd.getV();
    }

    public DoubleMatrix2D getPrincipalComponents() {
    	if (transposed) return svd.getV().zMult(svd.getS(), null).viewDice().copy();
    	else return svd.getU().zMult(svd.getS(), null).viewDice().copy();
    }
    
    public DoubleMatrix1D getEigenValues() {
    	double[] s = svd.getSingularValues();
        int n = s.length;
        int size=genotypes.rows();

    	double[] eigenvalues = new double[n];
    	for (int i = 0; i < n; i++) eigenvalues[i] = s[i] * s[i]/(size - 1);
    	return DoubleFactory1D.dense.make(eigenvalues);
    }
}

