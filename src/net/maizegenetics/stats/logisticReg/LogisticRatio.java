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
package net.maizegenetics.stats.logisticReg;

import net.maizegenetics.util.RightJustifiedFormat;

import java.io.Serializable;
import java.util.Vector;

/**
 * Title:        TASSEL
 * Description:  A java program to deal with diversity
 * Copyright:    Copyright (c) 2000
 * Company:      USDA-ARS/NCSU
 * @author Ed Buckler
 * @version 1.0
 */

public class LogisticRatio extends Logistic {

    Vector theLogisticRatioResultsVector=new Vector(1);
    Matrix theX2Matrix;
/*    int[] permutations;
    int[] permutationsGreaterThenMaxDiffLL;
    double[] permutedMaxP;
    double[] locusMaxP;*/

    //variable to holds information for permutations
    private int indexOfOverall=-1000, indexOfBestDiffL, currRep;
    private double bestDiffLL=-9999, permuteBestDiffLL=-9999;


    public LogisticRatio(boolean includesMissing, double missingValue, String[] annotationLabel) {
      super(includesMissing, missingValue, annotationLabel);
          //put a empty object in the space for the overall p-value object
      theLogisticRatioResultsVector.add(0,new Object());
    }

    public LogisticRatio(String[] annotationLabel) {
      super(annotationLabel);
          //put a empty object in the space for the overall p-value object
      theLogisticRatioResultsVector.add(0,new Object());
    }


    /**
     * This sets the data for the analysis.  The Y matrix are the dependent variables and should only have two states per column
     * theX is predictor variables that are used in all the analyses.  theX2 is the alternative predictor hypothesis.
     * The annotation contains the comments for each of the columns states.
     */
    public void setData(double[][] theY, double[][] theX, double[][] theX2, String[][] annotation) {
      this.annotation=annotation;
//      theXMatrix=new Matrix(theX);
//      theX2Matrix=new Matrix(theX2);
      setAndFilterXMatrices(theX, theX2, theY);
    }

    /**
     * This sets the data for the analysis.  The Y matrix are the dependent variables and should only have two states per column
     * theX is predictor variables that are used in all the analyses.  theX2 is the alternative predictor hypothesis.
     * The annotation contains the comments for each of the columns states.
     */
    public void setData(int[][] theY, double[][] theX, double[][] theX2, String[][] annotation) {
      this.annotation=annotation;
//      theXMatrix=new Matrix(theX);
//      theX2Matrix=new Matrix(theX2);
      Matrix ty=new Matrix(theY);
      setAndFilterXMatrices(theX, theX2, ty.element);
    }

    /**
     * this method removes rows with data that are missing from either X1 or X2
     * this is important so that the null logistic models will the same for both
     * note that the Y matrix is also filtered but the status of its missing data is not a concern at this level
     * (these are handled later in for each separate logistic regression)
     */
    private void setAndFilterXMatrices(double[][] theX, double[][] theX2, double[][] theY) {
      if(theX.length!=theX2.length) {theXMatrix=theX2Matrix=null; return;}
      double[][] tempX, tempX2;
      int[] includeRows=new int[theX.length];
      int count=0;
      boolean missingInRow;
      for (int i = 0; i<theX.length; i++) {
          missingInRow=false;
          for (int j=0; j<theX[0].length; j++)
            {if((theX[i][j]==missingValue)||(Double.isNaN(theX[i][j])))
              {missingInRow=true; break;}
            }
          for (int j=0; j<theX2[0].length; j++)
            {if((theX2[i][j]==missingValue)||(Double.isNaN(theX2[i][j])))
              {missingInRow=true; break;}
            }
          if(!missingInRow) {includeRows[count]=i;  count++;}
        }
      theXMatrix = new Matrix(count, theX[0].length);
      theX2Matrix = new Matrix(count,theX2[0].length);
      theYMatrix = new Matrix(count,theY[0].length);
      for (int i = 0; i<count; i++) {
          for (int j=0; j<theX[0].length; j++)
            {theXMatrix.element[i][j]=theX[includeRows[i]][j];}
          for (int j=0; j<theX2[0].length; j++)
            {theX2Matrix.element[i][j]=theX2[includeRows[i]][j];}
          for (int j=0; j<theY[0].length; j++)
            {theYMatrix.element[i][j]=theY[includeRows[i]][j];}
        }
    }

    /**
     * This cycles through each of the Y variables, calculates logistic regression and ratio for each.
     */
    public Vector computeAllSingleFactors(int[] indices) {
      Vector resultVector=new Vector();
      for(int i=0; i<theYMatrix.columns; i++)
        {Matrix analMat=theYMatrix.sub(0, theYMatrix.rows-1,i,i);
        try{
        Matrix joinMat=analMat.appendCols(theXMatrix);
        LogisticResults theLogisticResults1=buildClassifier(joinMat, annotation[i]);
        Matrix joinMat2=analMat.appendCols(theX2Matrix);
        LogisticResults theLogisticResults2=buildClassifier(joinMat2, annotation[i]);
        LogisticRatioResults theLogisticRatioResults=new LogisticRatioResults(theLogisticResults1,theLogisticResults2);
        theLogisticRatioResults.setIndex(indices[i]);
        theLogisticRatioResultsVector.add(theLogisticRatioResults);
        resultVector.add(theLogisticRatioResults);
        if(theLogisticRatioResults.diffLL>bestDiffLL)
          {this.bestDiffLL=theLogisticRatioResults.diffLL;
          LogisticRatioResults theLogisticRatioClone=(LogisticRatioResults)theLogisticRatioResults.clone();
          theLogisticRatioClone.setIndex(this.indexOfOverall);
          theLogisticRatioResultsVector.setElementAt(theLogisticRatioClone,0);
          }
//        System.out.println(theLogisticRatioResults.toString());
        }
        catch(Exception e)
          {System.out.println("LogisticRatio threw an exception"+e);}
        }
      return resultVector;
    }


        /**
     * This cycles through each of the Y variables, calculates logistic regression and ratio for each.
     */
    public Vector computeAllSingleFactorsPermutation(int[] indices, int rep) {
      Vector resultVector=new Vector();
      if(rep!=this.currRep)
        {currRep=rep;
//        System.out.println("rep="+rep+"  permuteBestDiffLL="+permuteBestDiffLL);
        if(permuteBestDiffLL!=-9999)  //then it has already run through the first permutation run
          {LogisticRatioResults ar=(LogisticRatioResults)theLogisticRatioResultsVector.get(0);
          ar.addPermutation(this.permuteBestDiffLL);
          permuteBestDiffLL=-9999;
          }
        }

      for(int i=0; i<theYMatrix.columns; i++)
        {Matrix analMat=theYMatrix.sub(0, theYMatrix.rows-1,i,i);
        try{
        Matrix joinMat=analMat.appendCols(theXMatrix);
        LogisticResults theLogisticResults1=buildClassifier(joinMat, annotation[i]);
        Matrix joinMat2=analMat.appendCols(theX2Matrix);
        LogisticResults theLogisticResults2=buildClassifier(joinMat2, annotation[i]);
        LogisticRatioResults permute=new LogisticRatioResults(theLogisticResults1,theLogisticResults2);
        permute.setIndex(indices[i]);
        int j=theLogisticRatioResultsVector.indexOf(permute);
        LogisticRatioResults real=(LogisticRatioResults)theLogisticRatioResultsVector.get(j);
        real.addPermutation(permute.diffLL);
        resultVector.add(permute);
        if((permute.diffLL>=permuteBestDiffLL)&&(!Double.isInfinite(permute.diffLL)))
          {this.permuteBestDiffLL=permute.diffLL;}
//        System.out.println(theLogisticRatioResults.toString());
        }
        catch(Exception e)
          {System.out.println("LogisticRatioPermutation threw an exception"+e);}
        }
      return resultVector;
    }

         //Implementation of TableReport Interface
    public Object[] getTableColumnNames() {
      String[] theLabels;
      String[] basicLabels={"Col","DiffLL","PermuteP","ModelLogL1","NullLogL1","d.f.1","ChiSqr1","P2", "ModelLogL2","NullLogL2","d.f.2","ChiSqr2","P2"};
      if(annotationLabel==null)
        {return basicLabels;}
       else
        {theLabels=new String[annotationLabel.length+basicLabels.length];
        for(int i=0; i<annotationLabel.length; i++)
          {theLabels[i]=annotationLabel[i];}
        for(int i=0; i<basicLabels.length; i++)
          {theLabels[i+annotationLabel.length]=basicLabels[i];}
        return theLabels;
        }
      }

    public Object[][] getTableData() {
      Object[][] data;
      java.text.NumberFormat nf=new java.text.DecimalFormat();
      nf.setMaximumFractionDigits(8);
      int basicCols=13, labelOffset;
      if(annotationLabel==null)
        {data=new String[theLogisticResultsVector.size()][basicCols];}
      else
        {data=new String[theLogisticRatioResultsVector.size()][annotationLabel.length+basicCols];
        labelOffset=annotationLabel.length;
        }
      LogisticResults theLogisticResults;
      LogisticRatioResults theLogisticRatioResults;
      for(int i=0; i<theLogisticRatioResultsVector.size(); i++)
        {theLogisticRatioResults=(LogisticRatioResults)theLogisticRatioResultsVector.get(i);
        labelOffset=0;
        if(annotationLabel!=null)
          {for(int j=0; j<annotationLabel.length; j++)
            {data[i][labelOffset++]=theLogisticRatioResults.lrModel1.annotation[j];}
          }
        if(theLogisticRatioResults.index==indexOfOverall)
          {data[i][labelOffset++]="OverallSignificance"+theLogisticRatioResults.permutations;}
         else
          {data[i][labelOffset++]=""+i;}
        data[i][labelOffset++]=""+theLogisticRatioResults.diffLL;
        data[i][labelOffset++]=""+theLogisticRatioResults.permutedP;
        theLogisticResults=theLogisticRatioResults.lrModel1;
        data[i][labelOffset++]=""+theLogisticResults.m_LL;
        data[i][labelOffset++]=""+theLogisticResults.m_LLn;
        data[i][labelOffset++]=""+theLogisticResults.df;
        data[i][labelOffset++]=""+theLogisticResults.CSq;
        data[i][labelOffset++]=""+theLogisticResults.P;
        theLogisticResults=theLogisticRatioResults.lrModel2;
        data[i][labelOffset++]=""+theLogisticResults.m_LL;
        data[i][labelOffset++]=""+theLogisticResults.m_LLn;
        data[i][labelOffset++]=""+theLogisticResults.df;
        data[i][labelOffset++]=""+theLogisticResults.CSq;
        data[i][labelOffset++]=""+theLogisticResults.P;
        }
      return data;
      }

    public String getTableTitle() {
      return "Logistic Regression Ratios with pop structure";
    }

  public String toString() {
    if(theLogisticRatioResultsVector.size()==0) {return "Needs to be implemented";}
    StringBuffer cs=new StringBuffer();
    Object[] header=this.getTableColumnNames();
    for(int i=0; i<header.length; i++)
      {cs.append(RightJustifiedFormat.monoString(header[i]));}
    cs.append("\n");
    Object[][] data=this.getTableData();
    for(int i=0; i<data.length; i++)
      {for(int j=0; j<data[i].length; j++)
        {cs.append(RightJustifiedFormat.monoString(data[i][j]));}
      cs.append("\n");
      }
    return cs.toString();
  }

}

class LogisticRatioResults  implements  Serializable, Cloneable, StatisticResult {
  /** The difference in the log-likelihoods of the two models */
  protected double diffLL;

  protected double P;

  protected LogisticResults lrModel1, lrModel2;
  int permutations=0;
  int permutationsGreaterThenDiffLL;
  double permutedP=Double.NaN;
  int index;

  public LogisticRatioResults(LogisticResults lr1, LogisticResults lr2)
    {this.lrModel1=lr1;
    this.lrModel2=lr2;
    diffLL=lrModel1.m_LL-lrModel2.m_LL;
    }

  public void addPermutation(double pDiffLL) {
    if(Double.isInfinite(pDiffLL)) return;  //this handles infinite likelihoods that get caused by missing data, resulting in variation in the dependent state
    permutations++;
    if(pDiffLL>=diffLL) {permutationsGreaterThenDiffLL++;}
    permutedP=(double)permutationsGreaterThenDiffLL/(double)permutations;
  }

  public Object clone() {
          try {
           LogisticRatioResults lrr=(LogisticRatioResults)super.clone();
           return lrr;
          }
          catch (CloneNotSupportedException e) {
              throw new InternalError(e.toString());
          }
      }

  public boolean equals(Object anObject) {
    LogisticRatioResults x=(LogisticRatioResults)anObject;
    if(x.index==this.index) {return true;}
      else {return false;}
    }

  public void setIndex(int theIndex) {
    this.index=theIndex;
    }

  public String toString() {
    String result = "Logistic Ratio Regression (2 classes)\n";
    result+="Model1 LL ="+lrModel1.m_LL+"   "+"Model2 LL ="+lrModel2.m_LL+"\n";
    result+="Difference in models ="+diffLL+"\n";
    result+=lrModel1.toString();
    result+=lrModel1.toString();
    return result;
  }

    //StatisticResult interface methods

  public double testStatistic() {
    return diffLL;
    }

  public double testPValue() {
//    System.out.println("diffLL="+diffLL);
    if(Double.isInfinite(diffLL)||Double.isNaN(diffLL)) return Double.NaN;
    return 1.0-net.maizegenetics.pal.statistics.ChiSquareDistribution.cdf(Math.abs(diffLL),1);
    }

  public double testLikelihood()  {
    return diffLL;
    }

  public String[] getAttributes() {
    return lrModel1.annotation;
    }

  public String getAttribute(int a) {
    return lrModel1.annotation[a];
    }
}
