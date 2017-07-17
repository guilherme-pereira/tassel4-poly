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
package net.maizegenetics.pal.distance;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.BitUtil;
import net.maizegenetics.util.ProgressListener;

/**
 * This class calculates an identity by state matrix. It is scaled so only
 * non-missing comparison are used.  It conducts bit level calculations of IBS for genotypes.
 * Only the two most common alleles are used in the distance calculations.
 * <p>
 * Please note that when heterozygous genotypes are used, Het to Het distance is 0.5 NOT 0.0. The default
 * along the identity diagonal is 0 (isTrueIBS = false), but changing isTrueIBS = true will calculate
 * the identity.
 * <p>
 * The distance estimates become wildly inaccuate when two few sites are used to calculate
 * distance.  The minSiteComp parameter can be used to control the minimum number of sites
 * used for a calculation.  If there are insufficient sites in the estimate, then Double.NaN
 * is returned.
 * 
 * @author Ed Buckler
 * @version 1.0
 */
public class IBSDistanceMatrix extends DistanceMatrix {

    private ProgressListener myListener = null;
    private int numSeqs;
    private Alignment theTBA = null;
    /**
     * Holds the average numbers of sites in the comparisons
     */
    private double avgTotalSites;
    private int minSitesComp = 0;
    private boolean isTrueIBS = false;

    /**
     * Compute observed distances for all taxa. Missing sites are ignored.
     *
     * @param theAlignment Alignment used to computed distances
     */
    public IBSDistanceMatrix(Alignment theAlignment) {
        this(theAlignment, 0, null);
    }

    /**
     * Compute observed distances for all taxa. Missing sites are ignored.
     * 
     * @param theAlignment Alignment used to computed distances
     * @param listener Listener to track progress in calculations
     */
    public IBSDistanceMatrix(Alignment theAlignment, ProgressListener listener) {
        this(theAlignment, 0, listener);
    }

    /**
     * Compute observed distances for all taxa. Missing sites are ignored.
     * 
     * @param theAlignment Alignment used to computed distances
     * @param minSiteComp Minimum number of sites needed to estimate distance
     * @param listener Listener to track progress in calculations
     */
    public IBSDistanceMatrix(Alignment theAlignment, int minSiteComp, ProgressListener listener) {
        this(theAlignment, minSiteComp, false, listener);
    }

    /**
     * Compute observed distances for all taxa. Missing sites are ignored.
     * 
     * @param theAlignment Alignment used to computed distances
     * @param minSiteComp Minimum number of sites needed to estimate distance
     * @param trueIBS estimate diagonal distance based IBS (default = false, i=i=0.0)
     * @param listener Listener to track progress in calculations
     */
    public IBSDistanceMatrix(Alignment theAlignment, int minSiteComp, boolean trueIBS, ProgressListener listener) {
        super();
        this.minSitesComp = minSiteComp;
        isTrueIBS = trueIBS;
        myListener = listener;
        numSeqs = theAlignment.getSequenceCount();
        theTBA = AlignmentUtils.optimizeForTaxa(theAlignment, listener);
        //  this should have an option to only use the 2 or 3 most common alleles
        setIdGroup(theAlignment.getIdGroup());
        computeHetBitDistances();
    }

    /**
     * This is a cleanest, fastest and most accurate way to calculate distance.
     */
    private void computeHetBitDistances() {
        avgTotalSites = 0;
        int count = 0;
        double[][] distance = new double[numSeqs][numSeqs];
        for (int i = 0; i < numSeqs; i++) {
            long[] iMj = theTBA.getAllelePresenceForAllSites(i, 0).getBits();
            long[] iMn = theTBA.getAllelePresenceForAllSites(i, 1).getBits();
            for (int j = i; j < numSeqs; j++) {
                if (j == i && !isTrueIBS) {
                    distance[i][i] = 0;
                } else {
                    long[] jMj = theTBA.getAllelePresenceForAllSites(j, 0).getBits();
                    long[] jMn = theTBA.getAllelePresenceForAllSites(j, 1).getBits();
                    double[] result=computeHetBitDistances(iMj, iMn, jMj, jMn, minSitesComp);
                    distance[i][j] = distance[j][i] = result[0];
                    avgTotalSites += result[1];  //this assumes not hets
                    count++;
                }
            }
            fireProgress((int) (((double) (i + 1) / (double) numSeqs) * 100.0));
        }
        setDistances(distance);
        avgTotalSites /= (double) count;
    }

    /**
     * Compute distance for a pair of taxa.
     * @param theTBA input alignment
     * @param taxon1 index of taxon 1
     * @param taxon2 index of taxon 2
     * @return array of {distance, number of sites used in comparison}
     */
    public static double[] computeHetBitDistances(Alignment theTBA, int taxon1, int taxon2) {
        return computeHetBitDistances(theTBA, taxon1, taxon2, 0, false);
    }

    /**
     * Compute distance for a pair of taxa.
     * @param theTBA input alignment
     * @param taxon1 index of taxon 1
     * @param taxon2 index of taxon 2
     * @param minSitesCompared Minimum number of sites needed to estimate distance
     * @param trueIBS estimate diagonal distance based IBS (default = false, i=i=0.0)
     * 
     * @return array of {distance, number of sites used in comparison}
     */
    public static double[] computeHetBitDistances(Alignment theTBA, int taxon1, int taxon2, int minSitesCompared, boolean isTrueIBS) {
        if(theTBA.isTBitFriendly()==false) theTBA = AlignmentUtils.optimizeForTaxa(theTBA);
        long[] iMj = theTBA.getAllelePresenceForAllSites(taxon1, 0).getBits();
        long[] iMn = theTBA.getAllelePresenceForAllSites(taxon1, 1).getBits();
        long[] jMj = theTBA.getAllelePresenceForAllSites(taxon2, 0).getBits();
        long[] jMn = theTBA.getAllelePresenceForAllSites(taxon2, 1).getBits();
        return computeHetBitDistances(iMj, iMn, jMj, jMn, minSitesCompared, 0, iMj.length-1); 
    }
    
    
    /**
     * Compute distance for a pair of taxa.  Optimized for calculations sites within a certain 
     * range of underlying word (64 sites chunks) in the TBit array
     * @param theTBA input alignment
     * @param taxon1 index of taxon 1
     * @param taxon2 index of taxon 2
     * @param minSitesCompared Minimum number of sites needed to estimate distance
     * @param startWord starting word for calculating distance site=(startWord*64)
     * @param endWord ending word for calculating distance inclusive site=(endWord*64+63)
     * @param maskBadSet Optional mask for sites (those set to 1 are kept) 
     * 
     * @return array of {distance, number of sites used in comparison}
     */
    public static double[] computeHetBitDistances(Alignment theTBA, int taxon1, int taxon2, 
            int minSitesCompared, int startWord, int endWord, BitSet maskBadSet) {
        if(theTBA.isTBitFriendly()==false) theTBA = AlignmentUtils.optimizeForTaxa(theTBA);
        long[] iMj = theTBA.getAllelePresenceForAllSites(taxon1, 0).getBits();
        long[] iMn = theTBA.getAllelePresenceForAllSites(taxon1, 1).getBits();
        if(maskBadSet!=null) {
            long[] maskBad=maskBadSet.getBits();
            for (int i = 0; i < iMj.length; i++) {iMj[i]=iMj[i]& maskBad[i];}
            for (int i = 0; i < iMn.length; i++) {iMn[i]=iMn[i]& maskBad[i];}
        }
        long[] jMj = theTBA.getAllelePresenceForAllSites(taxon2, 0).getBits();
        long[] jMn = theTBA.getAllelePresenceForAllSites(taxon2, 1).getBits();
        return computeHetBitDistances(iMj, iMn, jMj, jMn, minSitesCompared, 0, iMj.length-1);
    }
    
    /**
     * Calculation of distance using the bit vector of major and minor alleles.
     * @param iMj Vector of major alleles for taxon i
     * @param iMn Vector of minor alleles for taxon i
     * @param jMj Vector of major alleles for taxon j
     * @param jMn Vector of minor alleles for taxon j
     * @param minSitesCompared Minimum number of sites needed to estimate distance
     * @return array of {distance, number of sites used in comparison}
     */
    public static double[] computeHetBitDistances(long[] iMj, long[] iMn, long[] jMj, long[] jMn, 
            int minSitesCompared) {
        return computeHetBitDistances(iMj, iMn, jMj, jMn, minSitesCompared, 0, iMj.length-1);
    }

     /**
     * Calculation of distance using the bit vector of major and minor alleles.
     * @param iMj Vector of major alleles for taxon i
     * @param iMn Vector of minor alleles for taxon i
     * @param jMj Vector of major alleles for taxon j
     * @param jMn Vector of minor alleles for taxon j
     * @param minSitesCompared Minimum number of sites needed to estimate distance
     * @param endWord ending word for calculating distance inclusive site=(endWord*64+63)
     * @param maskBadSet Optional mask for sites (those set to 1 are kept) 
     * @return array of {distance, number of sites used in comparison}
     */
    public static double[] computeHetBitDistances(long[] iMj, long[] iMn, long[] jMj, long[] jMn, 
            int minSitesCompared, int startWord, int endWord) {
        int sameCnt = 0, diffCnt = 0, hetCnt = 0;
        for (int x = startWord; x <= endWord; x++) {
            long same = (iMj[x] & jMj[x]) | (iMn[x] & jMn[x]);
            long diff = (iMj[x] & jMn[x]) | (iMn[x] & jMj[x]);
            long hets = same & diff;
            sameCnt += BitUtil.pop(same);
            diffCnt += BitUtil.pop(diff);
            hetCnt += BitUtil.pop(hets);
        }
        int sites = sameCnt + diffCnt - hetCnt;
        double identity = ((double) (sameCnt) - (double)(0.5*hetCnt)) / (double) (sites);
        double dist = 1 - identity;
        if (sites > minSitesCompared) {
            return new double[] {dist,(double)sites};
        } else {
            return new double[] {Double.NaN,(double)sites};
        }
    }
    
    /*
     * Average number of sites used in calculating the distance matrix
     */
    public double getAverageTotalSites() {
        return avgTotalSites;
    }

    public String toString(int d) {
        double[][] distance = this.getDistances();
        /*Return a string representation of this matrix with 'd'
         displayed digits*/
        String newln = System.getProperty("line.separator");
        String outPut = new String();
        String num = new String();
        int i, j;
        java.text.NumberFormat nf = new java.text.DecimalFormat();
        nf.setMaximumFractionDigits(5);
        for (i = 0; i < distance.length; i++) {
            for (j = 0; j < distance[i].length; j++) {

                //Numeric x = new Numeric(distance[i][j]);
                num = nf.format(d);
                //num = x.toString(d);
                //ed change that screws up formatting
                //num=""+this.element[i][j];
                outPut = outPut + num + (char) 9;
            }
            outPut = outPut + newln;
        }
        return outPut;
    }

    @Override
    public String toString() {
        return this.toString(6);
    }

    protected void fireProgress(int percent) {
        if (myListener != null) {
            myListener.progress(percent, null);
        }

    }

    /*
     * Returns whether true IBS is calculated for the diagonal 
     */
    public boolean isTrueIBS() {
        return isTrueIBS;
    }
}
