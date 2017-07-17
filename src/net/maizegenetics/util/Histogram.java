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

/**
 * Title:        QTP Analyzer<p>
 * Description:  Your description<p>
 * Copyright:    Copyright (c) 1997<p>
 * Company:      NCSU<p>
 * @author Ed Buckler
 * @version null
 */
package net.maizegenetics.util;

import java.text.DecimalFormat;

public class Histogram {
    double[] data, binPos, binFreq, binProp;
    int numDataPoints, numBins=10;
    double topDatum, bottomDatum, binWidth;
    double binFudgeFactor=0.0001; //

    public Histogram(double[][] sqrData) {
        //only use the lower triangle
        numDataPoints=(sqrData.length*(sqrData.length-1))/2;
        data=new double[numDataPoints];
        binFreq=new double[numBins];
        int c=0;
        for(int i=0; i<sqrData.length; i++)
            {for(int j=0; j<i; j++)
                {data[c]=sqrData[i][j];
                c++;
                }
            }
        sortData();
        topDatum=data[numDataPoints-1];
        bottomDatum=data[0];
        findBinPositions();
        calculateBinFreq(numDataPoints);
    }

    public Histogram(double[][] sqrData, double top, double bottom, int bins) {
        //only use the lower triangle
        this.numBins=bins;
        topDatum=top;
        bottomDatum=bottom;
        numDataPoints=(sqrData.length*(sqrData.length-1))/2;
        data=new double[numDataPoints];
        binFreq=new double[numBins];
        int c=0;
        for(int i=0; i<sqrData.length; i++)
            {for(int j=0; j<i; j++)
                {data[c]=sqrData[i][j];
                c++;
                }
            }
        sortData();
        findBinPositions();
        calculateBinFreq(numDataPoints);
    }

   public void addData(double[][] sqrData) {
        int currDataPoints=(sqrData.length*(sqrData.length-1))/2;
        data=new double[currDataPoints];
        int c=0;
        for(int i=0; i<sqrData.length; i++)
            {for(int j=0; j<i; j++)
                {data[c]=sqrData[i][j];
                c++;
                }
            }
        sortData();
        java.util.Arrays.sort(data);
        calculateBinFreq(currDataPoints);
        numDataPoints+=currDataPoints;
   }

    private void sortData() {
//     QuickSort qs=new QuickSort(new CompareNumber());
//     qs.sort(data);
    }

   private void findBinPositions() {
      binWidth=(topDatum-bottomDatum)/(double)numBins;
      binPos=new double[numBins];
      for(int i=0; i<numBins; i++)
         {binPos[i]=bottomDatum+((i+1)*binWidth)+binFudgeFactor;}
   }

   private void calculateBinFreq(int currDataPoints) {
     int currBin=0;
     for(int i=0; i<currDataPoints; i++)
         {while(data[i]>binPos[currBin]) {currBin++;}
         binFreq[currBin]+=1.0;
         }
    }

   public void convertFreqToProportion() {
/*      double totalBinFreq=0;
      for(int i=0; i<numBins; i++)
         {totalBinFreq+=binFreq[i];}*/
      binProp=new double[numBins];
      for(int i=0; i<numBins; i++)
         {binProp[i]=binFreq[i]/(double)numDataPoints;}
   }

   public String writeHistogram() {
      String s=new String();
      DecimalFormat dF=new DecimalFormat("0.000");
      double currBottom=bottomDatum;
      for(int i=0; i<numBins; i++)
         {s+=("Bin"+i+" ["+dF.format(currBottom)+" - "+dF.format(binPos[i])+"] = "+dF.format(binProp[i])+"\n");
         currBottom=binPos[i];
         }
      s+="\n";
      return s;
   }

   public double getFreqOfSampleBelowCutoff(double cutoff) {
      double cumFreq=0;
      for(int i=0; i<numBins; i++)
         {cumFreq+=binFreq[i];
         if(cumFreq>cutoff)
            return binPos[i];
         }
      return binPos[numBins-1];
   }

   public String kolmogorovSmirnovTest(Histogram expectedHist, double alpha) {
      //Approximate test for two freq. distributions with large sample sizes
      //Sokal & Rohlf Biometry 1995 p437-439
      if((this.numBins!=expectedHist.numBins)||(this.binWidth!=expectedHist.binWidth)||(this.bottomDatum!=expectedHist.bottomDatum))
         {return "These frequencies distributions were not setup in the same way and can't be compared.\n";}
      DecimalFormat dF=new DecimalFormat("0.0####");
      double obsF=0, expF=0, d, maxD=-1, kAlpha, dAlpha;
      double obsN=(double)numDataPoints, expN=(double)expectedHist.numDataPoints;
      String s="Kolmogorov-Smirnov Two Sample Test\n";
      for(int i=0; i<numBins; i++)
         {obsF+=binFreq[i];
         expF+=expectedHist.binFreq[i];
         d=Math.abs((obsF/obsN)-(expF/expN));
         if(d>maxD) maxD=d;
         }
      kAlpha=Math.sqrt(-0.5*Math.log(alpha/2));
      dAlpha=kAlpha*Math.sqrt((obsN+expN)/(obsN*expN));
      s+="?="+dF.format(alpha)+" K="+dF.format(kAlpha)+" critical D="+dF.format(dAlpha)+"\n";
      s+="Observed D="+dF.format(maxD)+"\n";
      if(maxD>dAlpha)
         {s+="D is significant at the "+dF.format(alpha)+" level.\n\n";}
       else
         {s+="D is NOT significant at the "+dF.format(alpha)+" level.\n\n";}
      return s;
   }
}