package net.maizegenetics.baseplugins.numericaltransform;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.ArrayList;
/**
 * User: Ed Buckler
 * Date: Oct 10, 2006
 * Time: 3:54:16 PM
 */
public class KNN {

    public static double[][] impute(double[][] data, int kNeighbors, boolean isManhattenDist, boolean isUnweighted) {
        int rows=data.length;
        int cols=data[0].length;
        double[][] normData=Conversion.normalizeData(data);
        double[][] distMat=createDistWithNaN(normData,isManhattenDist);
        double[][] result=new double[rows][cols];
         for(int i=0; i<rows; i++) {
            for(int j=0; j<cols; j++) {
                if(!Double.isNaN(data[i][j])) {result[i][j]=data[i][j];}
                else result[i][j]=calcKNN(i,j, data, distMat, kNeighbors, isUnweighted);
            }  //end j - cols or traits
         } //end i - taxa
        return result;
    }

    public static double calcKNN(int row, int col, double[][] data, double[][] distMat, int kNeighbors, boolean isUnweighted) {
        int theNextNeighbor=-1, kCount=0;
        double num=0, denom=0, minDist;
        ArrayList usedTaxa=new ArrayList();
        usedTaxa.add(row);
        while(kCount<kNeighbors) {
            minDist=1000000000;
            //this loops finds the nearest neighbor, defined by min distance, not used before, and has good data
            for(int i=0; i<distMat.length; i++) {
                if((distMat[row][i]<minDist)&&(usedTaxa.contains(i)==false)&&(!Double.isNaN(data[i][col]))) {
                    minDist=distMat[row][i];
                    theNextNeighbor=i;
                }
            }
            usedTaxa.add(theNextNeighbor);
            if(isUnweighted==false) {
                num+=(data[theNextNeighbor][col]*distMat[theNextNeighbor][col]);
                denom+=distMat[theNextNeighbor][col];
            } else {
                num+=data[theNextNeighbor][col];
                denom+=1.0;
            }
            kCount++;
        }
        return (num/denom);
    }

    public static double[][] createDistWithNaN(double[][] data, boolean isManhattenDist) {
        int rows=data.length;
        int cols=data[0].length;
        double[][] result=new double[rows][rows];
        double r, count, diff;
        for(int i=0; i<rows; i++) {
            for(int j=0; j<=i; j++) {
                r=0;
                count=0;
                for(int k=0; k<cols; k++) {
                    if((!Double.isNaN(data[i][k]))&&(!Double.isNaN(data[j][k]))) {
                        diff=Math.abs(data[i][k]-data[j][k]);
                        count+=1.0;
                        if(isManhattenDist) {r+=diff;}
                        else {r+=diff*diff;}  //else do Euclidian
                    }   //end of not NaN
                }   //end k
                if(count>0) {result[i][j]=result[j][i]=nNumD(r/count);}
                    else {result[i][j]=result[j][i]=Double.NaN;}
            } //end j
        }  //end i
        return result;
    }
     private static double nNumD(double d)  {
       return (new BigDecimal(d, new MathContext(5))).doubleValue();
  }
}
