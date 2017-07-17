/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.pal.popgen;

import java.util.TreeMap;
import net.maizegenetics.pal.alignment.*;
import net.maizegenetics.pal.distance.IBSDistanceMatrix;
import net.maizegenetics.util.OpenBitSet;


/**
 *  Scores error for each site of taxon based on IBS.
 * 
 * @author edbuckler
 */
public class IBSErrorByTaxon {
    private int blocks=-1;
    private int sites=0;
    private int[] errorCnt, mjCorrCnt, mnCorrCnt;
    

 
    public IBSErrorByTaxon(int tx, Alignment a, int minMinorCnt, int minMajorCnt, 
            int minSites, double maxDivergence) {
        OpenBitSet mjTbs=new OpenBitSet(a.getAllelePresenceForAllSites(tx, 0));
        OpenBitSet mnTbs=new OpenBitSet(a.getAllelePresenceForAllSites(tx, 1));
        blocks=mjTbs.getNumWords();
        int donorTaxaCnt=a.getSequenceCount();
        sites=a.getSiteCount();
        errorCnt=new int[sites];
        mjCorrCnt=new int[sites];
        mnCorrCnt=new int[sites];
        int mjCorrect=0, mnCorrect=0, error=0;
        for (int startBlock = 0; startBlock < blocks; startBlock++) {
            int[] window=getBlockWithMinCount(mjTbs.getBits(), mnTbs.getBits(), 
                    startBlock, minMinorCnt, minMajorCnt);      
            TreeMap<Double,long[][]> bestDonors=new TreeMap<Double,long[][]>();
            long[] iMj=a.getAllelePresenceForSitesBlock(tx, 0, window[0], window[1]);
            long[] iMn=a.getAllelePresenceForSitesBlock(tx, 1, window[0], window[1]);
            for (int d1 = 0; d1 < donorTaxaCnt; d1++) {
                if(d1==tx) continue;
                long[][] jM=new long[2][];
                jM[0]=a.getAllelePresenceForSitesBlock(d1, 0, window[0], window[1]);
                jM[1]=a.getAllelePresenceForSitesBlock(d1, 1, window[0], window[1]);
                double dist=IBSDistanceMatrix.computeHetBitDistances(iMj,iMn,jM[0],jM[1],minSites)[0];
                if(Double.isNaN(dist)) continue;
                if(dist<maxDivergence) {
                    bestDonors.put(new Double(dist), jM);
                }
            }
 //           System.out.printf("t:%d %s hits:%d %n",tx,Arrays.toString(window),bestDonors.size());
            OpenBitSet[] rsx=findCorrectErrors(iMj, iMn, bestDonors);
            mjCorrect+=recordCnts(startBlock, mjCorrCnt, rsx[0]);
            mnCorrect+=recordCnts(startBlock, mnCorrCnt, rsx[1]);
            error+=recordCnts(startBlock, errorCnt, rsx[2]);
            startBlock=window[1];
        }
        System.out.printf("%d:%s MjC:%d MnC:%d Er:%d %n",tx, a.getTaxaName(tx),
                mjCorrect, mnCorrect, error);
    }
    
    public int[] getMajorCorrectCounts() {
        return mjCorrCnt;
    }
    
    public int[] getMinorCorrectCounts() {
        return mnCorrCnt;
    }
    
    public int[] getErrorCounts() {
        return errorCnt;
    }
    
    private int recordCnts(int startBlock, int[] arrayCnt, OpenBitSet rsx) {
        int offset=startBlock*64;
        int[] setB=rsx.getIndicesOfSetBits();
        for (int i : setB) {
            arrayCnt[i+offset]++;
        }
        return setB.length;
    }
    
    private OpenBitSet[] findCorrectErrors(long[] iMj, long[] iMn, TreeMap<Double,long[][]> bestDonors){
        long[] mjCorrect=new long[iMj.length];
        long[] mnCorrect=new long[iMj.length];
        long[] wrong=new long[iMj.length];
        long[] mask=new long[iMj.length];
        for (int x = 0; x < iMj.length; x++) {
            mask[x]=iMj[x]|iMn[x];
        }
//        System.out.printf("%s %s %n",BitUtil.toPadString(iMj[0]),
//                            BitUtil.toPadString(iMn[0]));
//        System.out.println(bestDonors.toString());
        for (long[][] j : bestDonors.values()) {
            for (int x = 0; x < iMj.length; x++) {
                if(mask[x]==0) continue;
                //mask[0]=~0;
                mjCorrect[x]=mjCorrect[x]|(mask[x]&j[0][x]&iMj[x]);
                mnCorrect[x]=mnCorrect[x]|(mask[x]&j[1][x]&iMn[x]);
                wrong[x]=wrong[x]|(mask[x]&j[1][x]&iMj[x])|(mask[x]&j[0][x]&iMn[x]);
                mask[x]=(mask[x]&(~(mjCorrect[x]|mnCorrect[x]|wrong[x])));
//                if(x==-1) {
//                    System.out.printf("%s %s %n",BitUtil.toPadString(j[0][x]),
//                            BitUtil.toPadString(j[1][x]));
//                    System.out.printf("%s %s %s %s %n",BitUtil.toPadString(mjCorrect[x]),
//                            BitUtil.toPadString(mnCorrect[x]),BitUtil.toPadString(wrong[x]),
//                            BitUtil.toPadString(mask[x]));
//                }
            }
        }
    //    long[][] result={mjCorrect,mnCorrect,wrong};
        OpenBitSet[] result={new OpenBitSet(mjCorrect, mjCorrect.length), new OpenBitSet(mnCorrect, mnCorrect.length),
        new OpenBitSet(wrong, wrong.length)};
        return result;
    }
    
    
    /**
     * Given a start 64 site block, it expands to the left and right until it hits
     * the minimum Minor Site count in the target taxon
     * @param mnT - minor allele bit presence in a series of longs
     * @param startBlock
     * @param minMinorCnt
     * @return arrays of blocks {startBlock, endBlock}
     */
    private int[] getBlockWithMinCount(long[] mjT, long[] mnT, int startBlock, int minMinorCnt, int minMajorCnt) {
        int minorCnt=Long.bitCount(mnT[startBlock]);
        int majorCnt=Long.bitCount(mnT[startBlock]);
        int endBlock=startBlock;
        while((minorCnt<minMinorCnt)&&(majorCnt<minMajorCnt)&&((endBlock+1)<mjT.length)) {
            endBlock++;
            minorCnt+=Long.bitCount(mnT[endBlock]); 
            majorCnt+=Long.bitCount(mjT[endBlock]); 
        }
        if((blocks-endBlock)<(endBlock-startBlock)) endBlock=blocks-1;
        int[] result={startBlock, endBlock};
        return result;
    }
    
  
 
    
}

