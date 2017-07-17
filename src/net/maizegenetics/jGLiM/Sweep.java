package net.maizegenetics.jGLiM;

import java.util.ArrayList;
import java.util.Iterator;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.jet.math.Functions;

public class Sweep {
    public DoubleMatrix2D XTX;
    public DoubleMatrix2D XTy;
    public DoubleMatrix2D yTy;
    public DoubleMatrix2D A;
    public boolean[] singular;
    double[] css; //corrected ss of the columns of X, equals 0 if the column has only two distinct values
    double[] Dmin;
    ArrayList<Integer> columnsSwept;
    
    static final double TOL = 1e-8;
    
    public static void main(String[] args) {
        Sweep s = new Sweep();
        s.makeData();
        System.out.println(s.A);
        s.sweep();
    }
    
    public void adjust() {
        for (int k = 0; k < XTX.rows(); k++) {
            double B = A.getQuick(k,k);
            A.viewRow(k).assign(Functions.div(B));
            for (int i=0; i< A.rows(); i++) if (i != k) {
                B = A.getQuick(i,k);
                DoubleMatrix1D sub = A.viewRow(k).copy();
                sub.assign(Functions.mult(B));
                A.viewRow(i).assign(sub,Functions.minus);
            }
        }
    }
    
    public void sweep() {
//        for (int k = 0; k < XTX.rows(); k++) {
//            sweep(k);
//        }
//        System.out.println(A);
        makeData();
        for (int k = 0; k < XTX.rows(); k++) {
            revg2sweep(k);
        }
        System.out.println(A);
        revg2sweep(2);
        System.out.println(A);
        revg2sweep(4);
        System.out.println(A);
        
    }
    
    public void sweep(int column) {
        double D = A.getQuick(column,column);
        A.viewRow(column).assign(Functions.div(D));
        for (int i=0; i< A.rows(); i++) if (i != column) {
            double B = A.getQuick(i,column);
            DoubleMatrix1D sub = A.viewRow(column).copy();
            sub.assign(Functions.mult(B));
            A.viewRow(i).assign(sub,Functions.minus);
            A.setQuick(i,column,-1*B/D);
        }
        A.setQuick(column,column,1/D);
    }
    
    //if the diagonal element is close to 0, that row and column is set to zero and the function returns false
    public boolean g2sweep(int column) {
        double D = A.getQuick(column,column);
        if (D < Dmin[column]) {
            for (int r = 0; r < A.rows(); r++) A.setQuick(r,column,0);
            for (int c = 0; c < A.rows(); c++) A.setQuick(column,c,0);
            return false;
        }
        A.viewRow(column).assign(Functions.div(D));
        for (int i=0; i< A.rows(); i++) if (i != column) {
            double B = A.getQuick(i,column);
            DoubleMatrix1D sub = A.viewRow(column).copy();
            sub.assign(Functions.mult(B));
            A.viewRow(i).assign(sub,Functions.minus);
            A.setQuick(i,column,-1*B/D);
        }
        A.setQuick(column,column,1/D);
        return true;
    }
    
    //returns false if the column was not swept (linearly combination of columns already in the model, i.e. singular)
    public boolean revg2sweep(int column) {
        if (singular[column]) {
            singular[column] = false;
            return false;
        }
        double D = A.getQuick(column,column);
        if (D < Dmin[column]) {
            singular[column] = true;
            return false;
        }
        A.viewRow(column).assign(Functions.div(D));
        for (int i=0; i< A.rows(); i++) if (i != column) {
            double B = A.getQuick(i,column);
            DoubleMatrix1D sub = A.viewRow(column).copy();
            sub.assign(Functions.mult(B));
            A.viewRow(i).assign(sub,Functions.minus);
            A.setQuick(i,column,-1*B/D);
        }
        A.setQuick(column,column,1/D);
        return true;
    }

    public int sweepSingularColumns() {
        int swept = 0;
        for (int i = 0; i < singular.length; i++) {
            if (singular[i]) {
                singular[i] = false;
                if (revg2sweep(i)) swept++; 
            }
        }
        return swept;
    }

    public void makeData() {

        DoubleMatrix2D X = DoubleFactory2D.dense.make(new double[]{1,1,1,1,1,1,1,1,1,1,
                0.652,0.356,0.301,0.757,0.037,0.246,0.058,0.069,0.079,0.249,
                0.802,.832,0.255,0.082,0.856,0.076,0.377,0.119,0.475,0.723,
                .1, .2, .3, .2, .1, .5, .6, .7, .2, .8,
                0.802,.832,0.255,0.082,0.856,0.076,0.377,0.119,0.475,0.723}, 10);
        DoubleMatrix2D y = DoubleFactory2D.dense.make(new double[]{13.605,11.952,5.235,5.456,9.038,0.657,3.816,1.303,5.4,10.222}, 10);
        XTX = X.viewDice().zMult(X,null);
        XTy = X.viewDice().zMult(y,null);
        yTy = y.viewDice().zMult(y,null);

        DoubleMatrix2D[][] composite = new DoubleMatrix2D[][]{{XTX, XTy},{XTy.viewDice(), yTy}};
        A = DoubleFactory2D.dense.compose(composite);
        
        double[] thisCSS = new double[X.columns()];
        double N = XTX.get(0,0);
        for (int i = 0; i < thisCSS.length; i++) {
            thisCSS[i] = XTX.get(i,i) - XTX.get(i,0) * XTX.get(i,0) / N;
        }
        
        setcssDmin(thisCSS);
    }
    
    public void makeA() {
        A = DoubleFactory2D.dense.compose(new DoubleMatrix2D[][]{{XTX, XTy},{XTy.viewDice(), yTy}});
        singular = new boolean[XTX.rows()];
        for (int i = 0; i < singular.length; i++) singular[i] = false;
        setcssDmin(null);
    }
    
    public void setcssDmin(double css[]) {
    	if (css == null) {
    		if (XTX == null) {
    			Dmin = new double[A.columns()];
    		}
    		else {
    			Dmin = new double[XTX.columns()];
    		}
			for (int i = 0; i < Dmin.length; i++) Dmin[i] = TOL;
    	} else {

    		this.css = css;
    		Dmin = new double[css.length];
    		for (int i = 0; i < css.length; i++) {
    			if (css[i] < TOL) Dmin[i] = TOL;
    			else Dmin[i] = css[i] * TOL;
    		}
    	}
    }
    
    public void setDminFromA() {
       int n = XTX.columns();
       Dmin = new double[n];
       for (int i = 0; i < n; i++) Dmin[i] = A.getQuick(i, i) * TOL;
    }
    
    //these two helper functions are for testing effects (after all others)
    //the routines assume that all columns have been swept in, in order, using revg2
    //sweepColumnsOut returns the columns that were swept out. That int array should be used to 
    //call sweepCoumnsIn to return the A matrix to its original state
    
    //returns df = columns in effect swept minus other columns swept
    public int sweepColumnsOutReversible(int[] cols) {
        columnsSwept = new ArrayList<Integer>();
        int n = cols.length;
        int df = 0;
        for (int i = 0; i < n; i++) {
            if (revg2sweep(cols[i])) {
                df++;
                columnsSwept.add(cols[i]);
            }
        }
        int m = A.columns();
        for (int i = cols[n - 1] + 1; i < m - 1; i++) {
            if (singular[i]) {
                singular[i] = false;
                if (revg2sweep(i)) {
                    df--;
                    columnsSwept.add(i);
                }
            }
        }
        return df;
    }
    
    public void sweepColumnsInReversible() {
        Iterator<Integer> it = columnsSwept.iterator();
        while (it.hasNext()) {
            revg2sweep(it.next());
        }
    }
}
