package net.maizegenetics.jGLiM;
/*
 * jGLiM: Java for General Linear Models
 * for more information: http://www.maizegenetics.net
 *
 * Copyright (C) 2005 Peter Bradbury
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 */

//package net.maizegenetics.jGLiM;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.SingularValueDecomposition;

/**
 * Created using IntelliJ IDEA.
 * Author: Peter Bradbury
 * Date: Jan 31, 2005
 * Time: 9:42:26 AM
 */
public abstract class AbstractLinearModel implements LinearModel {
    //fields
    protected double ssModel;
    protected double ssError;
    protected double Rsq;
    protected int dfModel;
    protected int dfError;
    protected DoubleMatrix1D beta = null;
    protected int analysisType;
    protected int numberOfPermutations = 0;
    protected double permute_p = -1;
    protected double[] Fvalues = null;
    protected LinearModel modelForError = null;
    ProgressIndicator progress = null;

    //type 1 fits the model sequentially
    public static final int ANALYSIS_TYPE_1 = 1;

    //type 2 fits each factor after fitting the other factors
    public static final int ANALYSIS_TYPE_2 = 2;

    //type 3 is type 2 with sigma restrictions
    //In general, type 3 provides the best tests of the null hypothesis that a factor has all zero effects
    //when the data is unbalanced.
    public static final int ANALYSIS_TYPE_3 = 3;

    //methods
    public Statistic getStatistic(int statisticType) {

        switch (statisticType) {
            case Statistic.STATISTIC_SS_MODEL:
                return new BasicStatistic(new Double(ssModel), Statistic.STATISTIC_SS_MODEL);
            case Statistic.STATISTIC_SS_ERROR:
                return new BasicStatistic(new Double(ssError), Statistic.STATISTIC_SS_ERROR);
            case Statistic.STATISTIC_DF_MODEL:
                return new BasicStatistic(new Double(dfModel), Statistic.STATISTIC_DF_MODEL);
            case Statistic.STATISTIC_DF_ERROR:
                return new BasicStatistic(new Double(dfError), Statistic.STATISTIC_DF_ERROR);
            case Statistic.STATISTIC_MS_MODEL:
                return new BasicStatistic(new Double(getMSModel()), Statistic.STATISTIC_MS_MODEL);
            case Statistic.STATISTIC_MS_ERROR:
                return new BasicStatistic(new Double(getMSError()), Statistic.STATISTIC_MS_ERROR);
            case Statistic.STATISTIC_F:
                return new BasicStatistic(new Double(getF()), Statistic.STATISTIC_F);
            case Statistic.STATISTIC_P:
                return new BasicStatistic(new Double(getp()), Statistic.STATISTIC_P);
            case Statistic.STATISTIC_LSESTIMATES:
                if (beta == null) return new BasicStatistic("NULL", Statistic.STATISTIC_LSESTIMATES);
                return new BasicStatistic(beta, Statistic.STATISTIC_LSESTIMATES);
            case Statistic.STATISTIC_RSQ:
                return new BasicStatistic(new Double(Rsq), Statistic.STATISTIC_RSQ);
            case Statistic.STATISTIC_ANALYSIS_TYPE:
                return new BasicStatistic(getNameForAnalysisType(analysisType),
                        Statistic.STATISTIC_ANALYSIS_TYPE);
            case Statistic.STATISTIC_PERMUTATION_NUMBER:
                return new BasicStatistic(new Integer(getNumberOfPermutations()), Statistic.STATISTIC_PERMUTATION_NUMBER);
            case Statistic.STATISTIC_NULL_FDIST:
                if (Fvalues == null) return new BasicStatistic("NULL", Statistic.STATISTIC_NULL_FDIST);
                return new BasicStatistic(Fvalues, Statistic.STATISTIC_NULL_FDIST);
            case Statistic.STATISTIC_PERMUTE_P:
                return new BasicStatistic(new Double(getPermutationp()), Statistic.STATISTIC_PERMUTE_P);
            default:
                return null;
        }
    }
    
    
    public int getdfError() {
        return dfError;
    }

    public int getdfModel() {
        return dfModel;
    }

    public double getF() {
        if (modelForError != null) {
            return (ssModel/dfModel)/(modelForError.getMSError());
        }
        return ssModel * dfError / ssError / dfModel;
    }

    public double getMSError() {
        return ssError / dfError;
    }

    public double getMSModel() {
        return ssModel / dfModel;
    }

    public double getp() {
        double dferr;
        double F = getF();
        if (modelForError != null) dferr = modelForError.getdfError();
        else dferr = dfError;
        if (F < 1E-10 || dfModel < 1 || dferr < 1) return Double.NaN;
        return Ftest(getF(), dfModel, dferr);
    }

    public double getSSError() {
        return ssError;
    }

    public double getSSModel() {
        return ssModel;
    }

    public DoubleMatrix1D getBeta() {
        return beta;
    }

    public double getRsq() {
        return Rsq;
    }

    public void setdfError(int dfError) {
        this.dfError = dfError;
    }

    public void setSSError(double SSError) {
        this.ssError = SSError;
    }

    public int getAnalysisType() {
        return analysisType;
    }

    public void setAnalysisType(int analysisType) {
        this.analysisType = analysisType;
    }

    public int getNumberOfPermutations() {
        return numberOfPermutations;
    }

    public void setNumberOfPermutations(int numberOfPermutations) {
        this.numberOfPermutations = numberOfPermutations;
    }

    public double[] getNullFDistribution() {
        return Fvalues;
    }

    public double getPermutationp() {
        return permute_p;
    }

    public LinearModel getModelForError() {
        return modelForError;
    }

    public void setModelForError(LinearModel lm) {
        modelForError = lm;
    }

    public ProgressIndicator getProgressIndicator() {
        return progress;
    }

    public void setProgressIndicator(ProgressIndicator progress) {
        this.progress = progress;
    }
    

    //public static methods begin here

    /**
     * Calculates the general inverse of a matrix using singular value decomposition
     *
     * @param aMatrix - a rectangular, Colt 2-dimensional matrix of doubles
     * @param rank    - an int array which will hold the rank of aMatrix as rank[0]. If null, it will not return rank.
     * @return a DoubleMatrix2D ojbect that is a general inverse of aMatrix
     */
    public static DoubleMatrix2D geninv(DoubleMatrix2D aMatrix, int[] rank) {
        Algebra A = new Algebra();
        boolean transposeMatrix = false;
        if (aMatrix.rows() < aMatrix.columns()) {
            transposeMatrix = true;
            aMatrix = aMatrix.viewDice();
        }
        SingularValueDecomposition svd = new SingularValueDecomposition(aMatrix);
        DoubleMatrix2D invS = svd.getS();
//        if (rank != null) rank[0] = svd.rank();

        //calculate the inverse of S, a diagonal matrix with rank(aMatrix) non-zero elements
        int size = invS.rows();
        int r = 0;
        for (int i = 0; i < size; i++) {
            if (Math.abs(invS.get(i, i)) > 1E-10) {
                invS.set(i, i, 1 / invS.get(i, i));
                r++;
            }
            else
                invS.set(i, i, 0);
        }
        if (rank != null) rank[0] = r;
        DoubleMatrix2D minv = A.mult(A.mult(svd.getV(), invS), A.transpose(svd.getU()));
        if (transposeMatrix) return minv.viewDice();
        return minv;
    }

    /**
     * Calculates the p-value associated with an F statistic.  The returned p-value is the probability
     * that a greater F is drawn from the F-distribution.
     *
     * @param F             - value of the F statistic
     * @param numeratordf   - degreees of freedom in the numerator of F
     * @param denominatordf - degreees of freedom in the denominator of F
     * @return p-value
     */
    public static double Ftest(double F, double numeratordf, double denominatordf) {
        double k = denominatordf / (denominatordf + numeratordf * F);
        return cern.jet.stat.Gamma.incompleteBeta(denominatordf / 2, numeratordf / 2, k);
    }

    /**
     * Calculates the rank of a rectangular matrix.
     *
     * @param A - a Colt 2-dimensional matrix of doubles
     * @return the rank of A
     */
    public static int rankOf(DoubleMatrix2D A) {
        SingularValueDecomposition svd;
        double tol = 1E-10;
        int r = 0;

        if (A.rows() >= A.columns()) {
            svd = new SingularValueDecomposition(A);
        } else {
            svd = new SingularValueDecomposition(A.viewDice());
        }

        double[] s = svd.getSingularValues();
        for (int i = 0; i < s.length; i++) {
            if (s[i] > tol) {
                r++;
            }
        }
        return r;
    }
    
    public static int significantDigits(int numerator) {
    		return (int) Math.log10(numerator) + 1;
    }

    public void permute(int numberOfPermutations) {
        this.numberOfPermutations = numberOfPermutations;
        permute();
    }
    
    public int[] getDefaultStatistics() {
        return new int[] {
                Statistic.STATISTIC_DF_MODEL,
                Statistic.STATISTIC_DF_ERROR,
                Statistic.STATISTIC_MS_MODEL,
                Statistic.STATISTIC_MS_ERROR,
                Statistic.STATISTIC_F,
                Statistic.STATISTIC_P,
                Statistic.STATISTIC_RSQ
        };
    }

    public static String getNameForAnalysisType(int analysisType) {
        switch(analysisType) {
            case ANALYSIS_TYPE_1:
                return "1 (sequential)";
            case ANALYSIS_TYPE_2:
                return "2 (each factor fitted after others)";
            case ANALYSIS_TYPE_3:
                return "3 (restricted model)";
            default:
                return "";
        }
    }

}
