// ContigencyTable.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.pal.statistics;

import net.maizegenetics.pal.math.MersenneTwisterFast;

/**
 * Class for permuting contigency tables and determining the likelihood of the table.
 * If determining of the probability of a 2x2 table use FisherExact as it is much
 * faster.
 *
 * @version $Id: ContigencyTable.java,v 1
 *
 * @author Ed Buckler
 */
public class ContigencyTable {

MersenneTwisterFast intRand;
int[][] contig;
int csum, rows, cols;
int[] crow,ccol, rowDist, colDist;
float[][] expectation;
private double[] f; //this holds a large series of factorials
double mcF; //this has all the marginal factorials for calculating stats
int maxSize;


  /**
   * constructor for Contigency table
   *
   * @param maxSize is the maximum sum that will be encountered by contigency table
   */
  public ContigencyTable(int maxSize) {
    intRand=new MersenneTwisterFast();
    this.maxSize=2*maxSize;
    f=new double[this.maxSize+1];
    f[0]=0.0;
    for(int i=1; i<=this.maxSize; i++)
      {f[i]=f[i-1]+Math.log(i);}
  }


  /**
   * sets the data for the contigency table, must be set before other methods are called.
   * If tcontig has a greater count than maxSize, then the contig is set to null
   *
   * @param contig is the array of integers with observed states
   */
public void setMatrix(int[][] tcontig)
  {  //this permutes of rowDist to rapidly do the permutations
  int i,j,k,count=0;
  rows=tcontig.length;
  cols=tcontig[0].length;
  contig=tcontig;
  csum=0;
  crow=new int[rows];
  ccol=new int[cols];
  for(i=0; i<rows; i++)
	{for(j=0; j<cols; j++)
          {csum+=contig[i][j];
          crow[i]+=contig[i][j];
          ccol[j]+=contig[i][j];
          }
	}
  if(csum>1000) {
     // System.out.println("Holly Crap! Stop:"+csum);

  }
  if(csum>maxSize)
      {contig=null;
      return;}
  rowDist=new int[csum];     //This sets up the row distribution so that the random number requires only one call
  colDist=new int[csum];
  for(i=0; i<rows; i++) //calculate expected values only once
    {for(j=0; j<cols; j++)
      {for(k=0; k<contig[i][j]; k++)
          {rowDist[count]=i;
          colDist[count]=j;
          count++;
          }
        }
    }
  if(maxSize>0)    //this will calculate the multiple used in calculating the fisher exact tests
    {mcF=0;
    for(i=0; i<rows; i++) {mcF+=f[crow[i]];}
     for(j=0; j<cols; j++) {mcF+=f[ccol[j]];}
     mcF-=f[2*csum];
    }
}


final float calcchiSquare()
{
	int i,j;
	float chi=0, E;

	for(i=0; i<rows; i++)
          {
          for(j=0; j<cols; j++)
              {if (contig[i][j]>0)
                  {E=expectation[i][j];
                  chi+=(((float)contig[i][j]-E)*((float)contig[i][j]-E)/E);
                  }
              }
          }
//	printf("The chi=%8.4f\n",chi);
	return chi;
}

final double calcLnFisherExactP()    //the natural log of the fisher exact P
{
	int i,j;
	double lnFisherExactP=mcF;

	for(i=0; i<rows; i++)
            {for(j=0; j<cols; j++)
                {lnFisherExactP-=f[contig[i][j]];
                }
            }
//	printf("The chi=%8.4f\n",chi);
	return lnFisherExactP;
}


final void randomcontig()
{
	int i,j,temp, r,sum=0;

    for(i=0; i<rows; i++)
      {for(j=0; j<cols; j++)
        {contig[i][j]=0;}
      }
    for(i=0; i<csum; i++)
      {r=intRand.nextInt(csum);
      temp=rowDist[i];
      rowDist[i]=rowDist[r];
      rowDist[r]=temp;
      }
    for(i=0; i<csum; i++)
      {contig[rowDist[i]][colDist[i]]++;}
}

  /**
   * This calculates the probability in a rapid approach, using the Chi Square as the test statistic.
   * It runs for maxPermutations permutations unless it find 10 values that beat the observed, and
   * then it stops and calculates the p-value.  This slighly biases the P-values but makes it much more rapid.
   *
   * @param maxPermutations Number of permutations used to calculate the probability
   * @return P-value (NaN is returned if bad contig was set)
   */
public double calcRapidContigencyChiSquare(int maxPermutations)
{//this version is rapid because it has a cut off if 10 reps beat the observed
//then it quits, this is good for rapid analysis and when you are mostly interested in the low p-values
	int r,reps=maxPermutations;
	double X, chiprob=0.0f, Origchi=0.0f;

        if(contig==null) {return Double.NaN;}

	for(r=0; r<=reps; r++)
	  {if (r>0) {randomcontig();}
	  X=calcchiSquare();
	  if(r==0)
	    {Origchi=X;}
	  else
	    {if(X>=Origchi) {chiprob+=1.0;}}
          if(chiprob>9.0) {reps=r-1;}
	  }
  return chiprob/(double)r;
}

  /**
   * This calculates the probability in the normal approach, using the Chi Square as the test statistic.
   *
   * @param permutations Number of permutations used to calculate the probability
   * @return P-value (NaN is returned if bad contig was set)
   */
public double calcContigencyChiSquare (int permutations)
{//this version is slow, but does fixed number of reps; better if you want the distribution of P
	int r;
	double X, chiprob=0.0f, Origchi=0.0f;

        if(contig==null) {return Double.NaN;}

	for(r=0; r<=permutations; r++)
	  {if (r>0) {randomcontig();}
	  X=calcchiSquare();
	  if(r==0)
	    {Origchi=X;}
	  else
	    {if(X>=Origchi) {chiprob+=1.0;}}
	  }
  return chiprob/(double)r;
}

  /**
   * This calculates the probability in the rapid permutational approach, using the method described
   * by Weir, B. S. (1996) Genetic Data Analysis II (Sinauer, Sunderland, MA)
   * It runs for 1000 permutations unless it find 10 values that beat the observed, and
   * then it stops and calculates the p-value.  This slighly biases the P-values but makes it much more rapid.
   *
   * @param reps is the number of permutations used to the probability
   * @return P-value (NaN is returned if bad contig was set)
   */
public double calcRapidMonteCarloExactTest(int maxPermutations)
    //this is the weir based approach
    //this version is rapid because it has a cut off if 10 reps beat the observed
//then it quits, this is good for rapid analysis and when you are mostly interested in the low p-values
{
	int r,reps=maxPermutations;
	double fE, mCFEprob=0.0f, origMCFE=0.0f;

        if(contig==null) {return Double.NaN;}

	for(r=0; r<=reps; r++)
		{if (r>0) {randomcontig();}
		fE=calcLnFisherExactP();
	if(r==0)
	  {origMCFE=fE;}
		  else
		  	{if(fE<=origMCFE)	{mCFEprob+=1.0;}
		  	}
    if(mCFEprob>9.0)
      {reps=r-1;}
		}
  return mCFEprob/(double)r;
}

  /**
   * This calculates the probability in the normal permutation approach, using the method described
   * by Weir, B. S. (1996) Genetic Data Analysis II (Sinauer, Sunderland, MA).
   *
   * @param permutations Number of permutations used to the probability
   * @return P-value (NaN is returned if bad contig was set)
   */
public double calcMonteCarloExactTest(int permutations)
    //this is the weir based approach
    //this version is slow, but does fixed number of reps; better if you want the distribution of P
{
	int r;
	double fE, mCFEprob=0.0f, origMCFE=0.0f;

        if(contig==null) {return Double.NaN;}

	for(r=0; r<=permutations; r++)
	  {if (r>0) {randomcontig();}
	  fE=calcLnFisherExactP();
    	if(r==0)
    	  {origMCFE=fE;}
    	else
    	  {if(fE<=origMCFE)	{mCFEprob+=1.0;}}
	}
  return mCFEprob/(double)r;
}


  void writeMatrix() {
    for(int i=0; i<rows; i++)
		  {System.out.print("r"+i+"   ");
      for(int j=0; j<cols; j++)
			  {System.out.print(contig[i][j]+"  ");}
      System.out.println();
		  }
  }

}