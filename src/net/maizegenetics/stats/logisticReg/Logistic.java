/*
 * TASSEL - Trait Analysis by a aSSociation Evolution & Linkage
 * Copyright (C) 2003 Ed Buckler
 *
 * This software evaluates linkage disequilibrium nucletide diversity and 
 * associations. For more information visit http://www.maizegenetics.net
 *
 * This software is distributed under GNU general public license and without
 * any warranty ot technical support.
 *
 * You can redistribute and/or modify it under the terms of GNU General 
 * public license. 
 *
 */
/*
 *    Logistic.java
 *    Copyright (C) 1999 Len Trigg, Eibe Frank
 *
 *    This program is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */
package net.maizegenetics.stats.logisticReg;

import net.maizegenetics.pal.report.TableReport;
import net.maizegenetics.pal.report.AbstractTableReport;

import java.io.Serializable;
import java.util.Vector;


/**
 * Class for building and using a two-class logistic regression model
 * with a ridge estimator.  <p>
 *
 * Reference: le Cessie, S. and van Houwelingen, J.C. (1997). <i>
 * Ridge Estimators in Logistic Regression.</i> Applied Statistics,
 * Vol. 41, No. 1, pp. 191-201. <p>
 *
 * Missing values are replaced using a ReplaceMissingValuesFilter, and
 * nominal attributes are transformed into numeric attributes using a
 * NominalToBinaryFilter.<p>
 *
 * Valid options are:<p>
 *
 * -D <br>
 * Turn on debugging output.<p>
 *
 * @author Len Trigg (trigg@cs.waikato.ac.nz)
 * @author Eibe Frank (eibe@cs.waikato.ac.nz)
 * @version $Revision: 1.5 $
 */
public class Logistic extends AbstractTableReport implements Serializable, TableReport {

    Vector theLogisticResultsVector = new Vector();
    Matrix theYMatrix, theXMatrix;
    boolean includesMissing;
    double missingValue = Double.NaN;
    String[][] annotation = null;
    String[] annotationLabel = null;
    /** The log-likelihood of the built model */
    protected double m_LL;
    /** The log-likelihood of the null model */
    protected double m_LLn;
    /** The coefficients of the model */
    protected double[] m_Par;
    /** The number of attributes in the model */
    protected int m_NumPredictors;
    /** The index of the class attribute */
    protected int m_ClassIndex;
    /** The ridge parameter. */
    protected double m_Ridge = 1e-8;
    /** The filter used to make attributes numeric. */
//  private NominalToBinaryFilter m_NominalToBinary;
    /** The filter used to get rid of missing values. */
//  private ReplaceMissingValuesFilter m_ReplaceMissingValues;
    /** Debugging output */
    protected boolean m_Debug = false;

    public Logistic(boolean includesMissing, double missingValue, String[] annotationLabel) {
        this.includesMissing = includesMissing;
        this.missingValue = missingValue;
        this.annotationLabel = annotationLabel;
    }

    public Logistic(String[] annotationLabel) {
        this.annotationLabel = annotationLabel;
    }

    /**
     * This sets the data for the analysis.  The Y matrix are the dependent variables and should only have two states per column
     * theX is predictor variables that are used in all the analyses.  The annotation contains the comments for each of the columns states
     */
    public void setData(double[][] theY, double[][] theX, String[][] annotation) {
        this.annotation = annotation;
        theXMatrix = new Matrix(theX);
        theYMatrix = new Matrix(theY);
    }

    /**
     * This sets the data for the analysis.  The Y matrix are the dependent variables and should only have two states per column
     * theX is predictor variables that are used in all the analyses.  The annotation contains the comments for each of the columns states
     */
    public void setData(int[][] theY, double[][] theX, String[][] annotation) {
        this.annotation = annotation;
        theXMatrix = new Matrix(theX);
        theYMatrix = new Matrix(theY);
    }

    /**
     * Returns probability.
     */
    protected static double Norm(double z) {
        return net.maizegenetics.pal.statistics.ChiSquareDistribution.cdf(z * z, 1);
    //Statistics.chiSquaredProbability(z * z, 1);
    }

    /**
     * Evaluate the probability for this point using the current coefficients
     *
     * @param instDat the instance data
     * @return the probability for this instance
     */
    protected double evaluateProbability(double[] instDat) {

        double v = m_Par[0];
        for (int k = 1; k <= m_NumPredictors; k++) {
            v += m_Par[k] * instDat[k];
        }
        v = 1 / (1 + Math.exp(-v));
        return v;
    }

    /**
     * Calculates the log likelihood of the current set of
     * coefficients (stored in m_Par), given the data.
     *
     * @param X the instance data
     * @param Y the class values for each instance
     * @param jacobian the matrix which will contain the jacobian matrix after
     * the method returns
     * @param deltas an array which will contain the parameter adjustments after
     * the method returns
     * @return the log likelihood of the data.
     */
    protected double calculateLogLikelihood(double[][] X, double[] Y,
            Matrix jacobian, double[] deltas) {

        double LL = 0;
        double[][] Arr = new double[jacobian.rows][jacobian.rows];

        for (int j = 0; j < Arr.length; j++) {
            for (int k = 0; k < Arr.length; k++) {
                Arr[j][k] = 0;
            }
            deltas[j] = 0;
        }

        // For each data point
        for (int i = 0; i < X.length; i++) {

            // Evaluate the probability for this point using the current coefficients
            double p = evaluateProbability(X[i]);

            // Update the log-likelihood of this set of coefficients
            if (Y[i] == 1) {
                LL = LL - 2 * Math.log(p);
            } else {
                LL = LL - 2 * Math.log(1 - p);
            }

            double w = p * (1 - p);     // Weight
            double z = (Y[i] - p);      // The error of this prediction

            for (int j = 0; j < Arr.length; j++) {
                double xij = X[i][j];
                deltas[j] += xij * z;
                for (int k = j; k < Arr.length; k++) {
                    Arr[j][k] += xij * X[i][k] * w;
                }
            }
        }

        // Add ridge adjustment to first derivative
        for (int j = 0; j < m_Par.length; j++) {
            deltas[j] -= 2 * m_Ridge * m_Par[j];
        }

        // Add ridge adjustment to second derivative
        for (int j = 0; j < Arr.length; j++) {
            Arr[j][j] += 2 * m_Ridge;
        }

        // Fill out the rest of the array
        for (int j = 1; j < Arr.length; j++) {
            for (int k = 0; k < j; k++) {
                Arr[j][k] = Arr[k][j];
            }
        }

        for (int j = 0; j < Arr.length; j++) {
            jacobian.setRow(j, Arr[j]);
        }

        return LL;
    }

    public LogisticResults buildClassifier(double[] preY, double[][] preX) throws Exception {
        Matrix theMatrix = new Matrix(preY, preX);
        return buildClassifier(theMatrix, null);
    }

    /**
     * Builds the classifier
     *
     * @param theMatrix has the dependent variable (0 or 1) in the first column, while
     * the rest of the matrix are the independent variables for the model
     * @exception Exception if the classifier could not be built successfully
     */
    public LogisticResults buildClassifier(Matrix theMatrix, String[] annotation) throws Exception {
        //this routine will remove any row with the supplied missing value or Double.NaN
        Matrix analMat = theMatrix.subRemoveMissing(this.missingValue);

//    int nR = m_NumPredictors = train.numAttributes() - 1;
        int nR = m_NumPredictors = analMat.columns - 1;
//    int nC = train.numInstances();
        int nC = analMat.rows;

        double[][] X = new double[nC][nR + 1];       // Data values
        double[] Y = new double[nC];               // Class values
        double[] xMean = new double[nR + 1];           // Attribute means
        double[] xSD = new double[nR + 1];           // Attribute stddev's
        double sY0 = 0;                                 // Number of class 0
        double sY1 = 0;                                 // Number of class 1

        if (m_Debug) {
            System.out.println("Extracting data...");
        }
        double firstY = analMat.element[0][0];
        for (int i = 0; i < X.length; i++) {
//      Instance current = train.instance(i);
            X[i][0] = 1;
            int j = 1;
            for (int k = 0; k < nR; k++) {
//	if (k != m_ClassIndex) {
                double x = analMat.element[i][j];
                /*if ((i == j) && (i <= nR)) {
                x += 1e-8;
                }*/
                X[i][j] = x;
                xMean[j] = xMean[j] + x;
                xSD[j] = xSD[j] + x * x;
                j++;
//	}
            }
            Y[i] = analMat.element[i][0];
            //Response variables are collapsed into two classes of either 0 or 1
            if (Y[i] == firstY) {
                Y[i] = 0;
            } else {
                Y[i] = 1;
            }
            if (Y[i] == 0) {
                sY0 = sY0 + 1;
            } else {
                sY1 = sY1 + 1;
            }
        }
        xMean[0] = 0;
        xSD[0] = 1;
        for (int j = 1; j <= nR; j++) {
            xMean[j] = xMean[j] / nC;
            xSD[j] = xSD[j] / nC;
            xSD[j] = Math.sqrt(Math.abs(xSD[j] - xMean[j] * xMean[j]));
        }
        if (m_Debug) {

            // Output stats about input data
            System.out.println("Descriptives...");
            System.out.println("" + sY0 + " cases have Y=0; " + sY1 + " cases have Y=1.");
            System.out.println("\n Variable     Avg       SD    ");
            for (int j = 1; j <= nR; j++) {
                System.out.println(j + "  " + xMean[j] + "   " + xSD[j]);
            }
        }

        // Normalise input data and remove ignored attributes
        for (int i = 0; i < nC; i++) {
            for (int j = 0; j <= nR; j++) {
                if (xSD[j] != 0) {
                    X[i][j] = (X[i][j] - xMean[j]) / xSD[j];
                }
            }
        }

        if (m_Debug) {
            System.out.println("\nIteration History...");
        }
        m_Par = new double[nR + 1];    // Coefficients
        double LLp = 2e+10;             // Log-likelihood on previous iteration
        m_LL = 1e+10;                  // Log-likelihood of current iteration
        m_LLn = 0;                      // Log-likelihood of null hypothesis

        double[] deltas = new double[nR + 1];
        Matrix jacobian = new Matrix(nR + 1, nR + 1);

        // Set up parameters for null model
        m_Par[0] = Math.log((sY1 + 1) / (sY0 + 1));
        for (int j = 1; j < m_Par.length; j++) {
            m_Par[j] = 0;
        }

        // While the log-likelihood is changing (i.e. no maxima found)
        while (Math.abs(LLp - m_LL) > 0.00001) {
            LLp = m_LL;
            m_LL = calculateLogLikelihood(X, Y, jacobian, deltas);
            if (LLp == 1e+10) {
                m_LLn = m_LL;
            }
            if (m_Debug) {
                System.out.println("-2 Log Likelihood = " + m_LL + ((m_LLn == m_LL) ? " (Null Model)" : ""));
            }

            jacobian.lubksb(jacobian.ludcmp(), deltas);

            // Adjust the coefficients
            for (int j = 0; j < deltas.length; j++) {
                m_Par[j] += deltas[j];
            }
        }
        if (m_Debug) {
            System.out.println(" (Converged)");
        }

        // Convert coefficients back to non-normalized attribute units
        for (int j = 1; j <= nR; j++) {
            if (xSD[j] != 0) {
                m_Par[j] = m_Par[j] / xSD[j];
                m_Par[0] = m_Par[0] - m_Par[j] * xMean[j];
            }
        }

        /*    try{*/
        LogisticResults theLR = new LogisticResults(annotation, m_LL, m_LLn, m_Par, m_NumPredictors);
//    System.out.println("LogisticResults chiSquare Error: m_LL"+m_LL+"  "+m_LLn+"  "+m_Par+"  "+m_NumPredictors);
//    System.out.println(analMat.toString());
        return theLR;
    /*    }
    catch(IllegalArgumentException ee)
    {
    System.out.println("LogisticResults chiSquare Error: m_LL"+m_LL+"  "+m_LLn+"  "+m_Par+"  "+m_NumPredictors);
    System.out.println(analMat.toString());
    return null;
    }
     */

    }


    //Implementation of TableReport Interface
    public Object[] getTableColumnNames() {
        String[] theLabels;
        String[] basicLabels = {"Col", "ModelLogL", "NullLogL", "d.f.", "ChiSqr", "P"};
        if (annotationLabel == null) {
            return basicLabels;
        } else {
            theLabels = new String[annotationLabel.length + basicLabels.length];
            for (int i = 0; i < annotationLabel.length; i++) {
                theLabels[i] = annotationLabel[i];
            }
            for (int i = 0; i < basicLabels.length; i++) {
                theLabels[i + annotationLabel.length] = basicLabels[i];
            }
            return theLabels;
        }
    }

    /**
     * Returns specified row.
     *
     * @param row row number
     *
     * @return row
     */
    public Object[] getRow(int row) {

        Object[] data;
        java.text.NumberFormat nf = new java.text.DecimalFormat();
        nf.setMaximumFractionDigits(8);
        int basicCols = 8, labelOffset;
        if (annotationLabel == null) {
            data = new String[basicCols];
        } else {
            data = new String[annotationLabel.length + basicCols];
            labelOffset = annotationLabel.length;
        }
        LogisticResults theLogisticResults;
        theLogisticResults = (LogisticResults) theLogisticResultsVector.get(row);
        labelOffset = 0;
        if (annotationLabel != null) {
            for (int j = 0; j < annotationLabel.length; j++) {
                data[labelOffset++] = theLogisticResults.annotation[j];
            }
        }
        data[labelOffset++] = "" + row;
        data[labelOffset++] = "" + theLogisticResults.m_LL;
        data[labelOffset++] = "" + theLogisticResults.m_LLn;
        data[labelOffset++] = "" + theLogisticResults.df;
        data[labelOffset++] = "" + theLogisticResults.CSq;
        data[labelOffset++] = "" + theLogisticResults.P;
        return data;

    }

    public String getTableTitle() {
        return "Logistic Regression with pop structure";
    }

    public String toString() {
        if (theLogisticResultsVector.size() == 0) {
            return "Values need to be computed";
        }
        String cs = "";
        return "Needs to be implemented";
    }

    public int getRowCount() {
        return theLogisticResultsVector.size();
    }

    public int getElementCount() {
        throw new UnsupportedOperationException();
    }

    public int getColumnCount() {
        throw new UnsupportedOperationException();
    }
}

class LogisticResults implements Serializable, Cloneable {

    String[] annotation;
    /** The log-likelihood of the built model */
    protected double m_LL;
    /** The log-likelihood of the null model */
    protected double m_LLn;
    /** The coefficients of the model */
    protected double[] m_Par;
    /** The number of attributes in the model */
    protected int df;
    protected double CSq;
    protected double P;

    public LogisticResults(double m_LL, double m_LLn, double[] m_ParTemp, int df) {
        this.m_LL = m_LL;
        this.m_LLn = m_LLn;
        this.df = df;
        this.m_Par = new double[m_ParTemp.length];
        for (int i = 0; i < m_ParTemp.length; i++) {
            m_Par[i] = m_ParTemp[i];
        }
        CSq = m_LLn - m_LL;
//    try{
        if (Double.isInfinite(m_LL)) {
            P = Double.NaN;
        } else {
            P = 1.0 - net.maizegenetics.pal.statistics.ChiSquareDistribution.cdf(CSq, df);
        }
//    }
//    catch(IllegalArgumentException ee)
//      {System.out.println("LogisticResults chiSquare Error: m_LL"+m_LL+"  "+m_LLn+"  "+m_ParTemp+"  "+df);}
    }

    public LogisticResults(String[] annotationM, double m_LL, double m_LLn, double[] m_ParTemp, int df) {
        this(m_LL, m_LLn, m_ParTemp, df);
        annotation = new String[annotationM.length];
        for (int i = 0; i < annotationM.length; i++) {
            annotation[i] = annotationM[i];
        }
    }

    public Object clone() {
        try {
            LogisticResults lr = (LogisticResults) super.clone();
            lr.m_Par = (double[]) this.m_Par.clone();
            lr.annotation = (String[]) this.annotation.clone();
            return lr;
        } catch (CloneNotSupportedException e) {
            throw new InternalError(e.toString());
        }
    }

    public String toString() {
        String result = "Logistic Regression (2 classes)";
        if (m_Par == null) {
            return result + ": No model built yet.";
        }
        result += "\n\nOverall Model Fit...\n" + "  Chi Square=" + CSq + ";  df=" + df + ";  p=" + net.maizegenetics.pal.statistics.ChiSquareDistribution.cdf(CSq, df) //      Statistics.chiSquaredProbability(CSq, df)
                + "\n";

        result += "\nCoefficients...\n" + "Variable      Coeff.\n";
        for (int j = 1; j <= df; j++) {
            result += j + "   " + m_Par[j] + "\n";
        }
        result += "Intercept " + m_Par[0] + "\n";

        result += "\nOdds Ratios...\n" + "Variable         O.R.\n";
        for (int j = 1; j <= df; j++) {
            double ORc = Math.exp(m_Par[j]);
            result += j + " " + ORc + "\n";
        }
        return result;
    }
}

