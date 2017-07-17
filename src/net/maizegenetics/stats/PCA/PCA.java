package net.maizegenetics.stats.PCA;


import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.EigenvalueDecomposition;

/**
 * Created by IntelliJ IDEA.
 * Author: Zhiwu Zhang
 * Date: Apr 26, 2005
 * Time: 10:17:57 AM
  */
public class PCA{

//fields
    private DoubleMatrix2D data;
    private boolean isPop;
    private boolean isCorrelation = true;
    PrincipalComponents prcomp;

//constructors: with data and its type (population or not)
public PCA(DoubleMatrix2D data) {
    this.data = data;
    this.isPop = false;
    build();
    }

//constructors: with data (not a sample as default) only
    public PCA(DoubleMatrix2D data, boolean isPop) {
    this.data = data;
    this.isPop = isPop;
    build();
    }

    public PCA(DoubleMatrix2D data, boolean isPop, boolean isCorrelation) {
    this.isPop = isPop;
    this.isCorrelation = isCorrelation;
    this.data = data;
    
    double sum = data.zSum();
    
    build();
    }

    public DoubleMatrix2D getEigenVectors() {
        System.out.println("Now is extracting EigenVectors...");
        return prcomp.getEigenVectors();
    }
    public DoubleMatrix1D getEigenValues() {
        System.out.println("Now is extracting Eigenvalues...");
        return prcomp.getEigenValues();
    }
    public DoubleMatrix2D getPC() {
        System.out.println("Now is calculating PC...");
        return prcomp.getPrincipalComponents();
    }


//methods: private
    private void build() {
        System.out.println("Doing singularValueDecomposition..., please wait.");

         if (isCorrelation) {
        	prcomp = new PrincipalComponents(data, PrincipalComponents.TYPE_CENTER_AND_SCALE);
        } else {
        	prcomp = new PrincipalComponents(data, PrincipalComponents.TYPE_CENTER);
        }

        System.out.println("PCA was built successfuly");
       
    }


}

