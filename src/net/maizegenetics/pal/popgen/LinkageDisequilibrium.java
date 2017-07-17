// LinkageDisequilibrium.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.pal.popgen;

import cern.colt.map.OpenLongObjectHashMap;
import java.io.Serializable;
import java.io.StringWriter;
import java.util.Arrays;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.BitAlignment;
import net.maizegenetics.pal.report.TableReport;
import net.maizegenetics.pal.statistics.FisherExact;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.ProgressListener;
import org.apache.log4j.Logger;

/**
 * This class calculates D' and r^2 estimates of linkage disequilibrium. It also
 * calculates the significance of the LD by either Fisher Exact or the
 * multinomial permutation test. This class can work with either normal
 * alignments of annotated alignments. The alignments should be stripped of
 * invariable numSites.
 * <p>
 * {@link testDesign} sets matrix design for LD calculation.  Either all by all, 
 * sliding window, site by all, or site list.
 * <p>
 * There are multiple approaches for dealing with heterozygous sites.  {@link HetTreatment} sets
 * the way these are treated.  Haplotype assumes fully 
 * phased heterozygous sites (any hets are double counted).  This is the best approach for speed
 * when things are fully phased.  Homozygous converted all hets to missing.  Genotype
  * does a 3x3 genotype analysis (to be implemented)
 *<p>
 * 2 state estimates of D' and r^2 can be found reviewed and discussed in Weir
 * 1996
 *<p>
 * Multi-state loci (>=3) require an averaging approach. In TASSEL 3 in 2010,
 * Buckler removed these approach as the relative magnitudes and meaningfulness
 * of these approaches has never been clear. Additionally with the moving away
 * from SSR to SNPs these methods are less relevant. Researchers should convert
 * to biallelic - either by ignoring rarer classes or collapsing rarer states.
 *<p>
 * TODO: Add 3x3 (genotype) mode.
 *
 * @version $Id: LinkageDisequilibrium.java,v 2
 *
 * @author Ed Buckler
 */
public class LinkageDisequilibrium extends Thread implements Serializable, TableReport {
    
    /**
     * Design of test matrix.
     */
    public static enum testDesign {
        /** Test All by All site*/ All, 
        /** Sliding of LD comparisons*/ SlidingWindow, 
        /**Test on site versus all others*/SiteByAll, 
        /**Test all sites with the site list with all others*/SiteList
    };
    /**
     * Approaches for dealing with heterozygous sites.  Haplotype assumes fully 
     * phased heterozygous sites (any hets are double counted).  This is the best approach for speed
     * when things are fully phased.  Homozygous converted all hets to missing.  Genotype
     * does a 3x3 genotype analysis (to be implemented)
     */
    public static enum HetTreatment {
        Haplotype, Homozygous, Genotype
    };
    private static final Logger myLogger = Logger.getLogger(LinkageDisequilibrium.class);
    private Alignment myAlignment;
    private Alignment mySBitAlignment;
    private int myMinTaxaForEstimate = 20;
    private int myWindowSize = 50;
    private int myTestSite = -1;  // this is only set when one versus all numSites is calculated.
    private long myTotalTests = 0;
    private testDesign myCurrDesign = testDesign.SlidingWindow;
    /**
     * HashMap of results Key = (site1*siteNum + site2), Value = float[rsqr,d', pvalue, sampleSize]
     */
    private OpenLongObjectHashMap  myMapResults;
    private ProgressListener myListener = null;
    private FisherExact myFisherExact;
    private boolean myIsAccumulativeReport = false;
    private int myNumAccumulativeBins = 100;
    private float myAccumulativeInterval;
    private int[] myAccumulativeRValueBins;
    private int[] mySiteList;
    private static String NotImplemented = "NotImplemented";
    private static String NA = "N/A";
    private static Integer IntegerTwo = Integer.valueOf(2);
    private HetTreatment myHetTreatment = HetTreatment.Homozygous;


    /**
     * Constructor for doing LD analysis
     * @param alignment Input alignment with segregating sites
     * @param windowSize Size of sliding window
     * @param LDType 
     * @param testSite
     * @param listener
     * @param isAccumulativeReport
     * @param numAccumulateIntervals
     * @param sitesList
     * @param hetTreatment 
     */
    public LinkageDisequilibrium(Alignment alignment, int windowSize, testDesign LDType, int testSite, ProgressListener listener, boolean isAccumulativeReport, int numAccumulateIntervals, int[] sitesList, HetTreatment hetTreatment) {
        myAlignment = alignment;
        mySBitAlignment = AlignmentUtils.optimizeForSites(myAlignment, listener);
        myFisherExact = new FisherExact((2*myAlignment.getSequenceCount()) + 10);
        myWindowSize = windowSize;
        myCurrDesign = LDType;
        myTestSite = testSite;
        myListener = listener;
        myIsAccumulativeReport = isAccumulativeReport;
        if (myIsAccumulativeReport) {
            myNumAccumulativeBins = numAccumulateIntervals;
        }
        mySiteList = sitesList;
        if (mySiteList != null) {
            Arrays.sort(mySiteList);
        }
        myHetTreatment = hetTreatment;
    }

    /**
     * starts the thread to calculate LD
     */
    public void run() {
        initMatrices();
        switch (myHetTreatment) {
            case Haplotype:
                calculateBitLDForHaplotype(false);
                break;
            case Homozygous:
                calculateBitLDForHaplotype(true);
                break;
            case Genotype:
                calculateBitLDWithHets();
                break;
            default:
                myLogger.error("Unknown LD analysis type selected for heterozygotes; skipping");
                break;
        }
    }

    private void initMatrices() {
        int numSites = myAlignment.getSiteCount();
        if (myCurrDesign == testDesign.All) {
            myTotalTests = numSites * (numSites - 1) / 2;
        } else if (myCurrDesign == testDesign.SlidingWindow) {
            long n = Math.min(numSites - 1, myWindowSize);
            myTotalTests = ((n * (n + 1)) / 2) + (numSites - n - 1) * n;
        } else if (myCurrDesign == testDesign.SiteByAll) {
            myTotalTests = numSites - 1;
        } else if (myCurrDesign == testDesign.SiteList) {
            long n = mySiteList.length;
            myTotalTests = ((n * (n + 1)) / 2) + (numSites - n - 1) * n;
        }
        if (myIsAccumulativeReport) {
            myAccumulativeInterval = 1.0f / (float) myNumAccumulativeBins;
            myAccumulativeRValueBins = new int[myNumAccumulativeBins + 1];
        } else {
           myMapResults=new OpenLongObjectHashMap(numSites); 
        }

    }
    
    private long getMapKey(int r, int c) {
        return (c<r)?(((long)c*mySBitAlignment.getSiteCount())+r):(((long)r*mySBitAlignment.getSiteCount())+c);
    }

    private void calculateBitLDForHaplotype(boolean ignoreHets) {
        //It will ignore hets, make a new Alignment and set all het calls to missing. Otherwise set the pointer to the old alignment
        Alignment workingAlignment;
        if (ignoreHets) {
            workingAlignment = BitAlignment.getHomozygousNucleotideInstance(myAlignment, true);
        } else {
            workingAlignment = mySBitAlignment;
        }
        for (long currTest = 0; currTest < myTotalTests; currTest++) {
            int r = getRowFromIndex(currTest);
            int c = getColFromIndex(currTest);
            int currentProgress = (int) ((double) 100.0 * ((double) currTest / (double) myTotalTests));
            fireProgress(currentProgress);
            BitSet rMj = workingAlignment.getAllelePresenceForAllTaxa(r, 0);
            BitSet rMn = workingAlignment.getAllelePresenceForAllTaxa(r, 1);
            BitSet cMj = workingAlignment.getAllelePresenceForAllTaxa(c, 0);
            BitSet cMn = workingAlignment.getAllelePresenceForAllTaxa(c, 1);
            LDResult ldr=getLDForSitePair(rMj, rMn, cMj, cMn, 2, myMinTaxaForEstimate, -1.0f, myFisherExact);
            if (myIsAccumulativeReport) {
                if (ldr.r2 == Float.NaN) {
                    myAccumulativeRValueBins[myNumAccumulativeBins]++;
                } else if (ldr.r2 == 1.0f) {
                    myAccumulativeRValueBins[myNumAccumulativeBins - 1]++;
                } else {
                    int index = (int) Math.floor(ldr.r2 / myAccumulativeInterval);
                    myAccumulativeRValueBins[index]++;
                }
            } else {           
                long key=getMapKey(r,c);
                myMapResults.put(key, ldr);               
            }

        } //end of currTest
        myMapResults.trimToSize();
    }
    

    private void calculateBitLDWithHets() {
        //Do nothing; not implemented yet
        myLogger.error("Calculating LD with hets as a third state is not implemented yet; skipping");
        throw new IllegalStateException("LinkageDisequilibrium: calculateBitLDWithHets: Treating hets as a third state is not yet implemented");
    }

    public static double calculateDPrime(int countAB, int countAb, int countaB, int countab, int minTaxaForEstimate) {
        //this is the normalized D' is Weir Genetic Data Analysis II 1986 p120
        double freqR, freqC, freq, countR, countC, nonmissingSampleSize;
        nonmissingSampleSize = countAB + countAb + countaB + countab;
        if (nonmissingSampleSize < minTaxaForEstimate) {
            return Double.NaN;
        }
        countR = countab + countAb;
        countC = countab + countaB;
        freqR = (nonmissingSampleSize - countR) / nonmissingSampleSize;
        freqC = (nonmissingSampleSize - countC) / nonmissingSampleSize;
        // if((freqR==0)||(freqC==0)||(freqR==1)||(freqC==1)) return -999;  //changed by ed 8-13-2004
        if ((freqR == 0) || (freqC == 0) || (freqR == 1) || (freqC == 1)) {
            return Double.NaN;
        }
        freq = ((double) countAB / nonmissingSampleSize) - (freqR * freqC);
        if (freq < 0) {
            return freq / Math.max(-freqR * freqC, -(1 - freqR) * (1 - freqC));
        } else {
            return freq / Math.min((1 - freqR) * freqC, (1 - freqC) * freqR);
        }  //check these equations
    }

    public static double calculateRSqr(int countAB, int countAb, int countaB, int countab, int minTaxaForEstimate) {
        //this is the Hill & Robertson measure as used in Awadella Science 1999 286:2524
        double freqA, freqB, rsqr, nonmissingSampleSize;
        nonmissingSampleSize = countAB + countAb + countaB + countab;
        if (nonmissingSampleSize < minTaxaForEstimate) {
            return Double.NaN;
        }
        freqA = (double) (countAB + countAb) / nonmissingSampleSize;
        freqB = (double) (countAB + countaB) / nonmissingSampleSize;

        //Through missing data & incomplete datasets some alleles can be fixed this returns missing value
        if ((freqA == 0) || (freqB == 0) || (freqA == 1) || (freqB == 1)) {
            return Double.NaN;
        }

        rsqr = ((double) countAB / nonmissingSampleSize) * ((double) countab / nonmissingSampleSize);
        rsqr -= ((double) countaB / nonmissingSampleSize) * ((double) countAb / nonmissingSampleSize);
        rsqr *= rsqr;
        rsqr /= freqA * (1 - freqA) * freqB * (1 - freqB);
        return rsqr;
    }

    public static LDResult getLDForSitePair(BitSet rMj, BitSet rMn, BitSet cMj, BitSet cMn,
            int minMinorCnt, int minCnt, float minR2, FisherExact myFisherExact) {
       // float[] results = {Float.NaN, Float.NaN, Float.NaN, Float.NaN};
        LDResult results=new LDResult();
        int n = 0;
        int[][] contig = new int[2][2];
        n += contig[1][1] = (int) OpenBitSet.intersectionCount(rMn, cMn);
        n += contig[1][0] = (int) OpenBitSet.intersectionCount(rMn, cMj);
        if (contig[1][0] + contig[1][1] < minMinorCnt) {
            return results;
        }
        n += contig[0][1] = (int) OpenBitSet.intersectionCount(rMj, cMn);
        if (contig[0][1] + contig[1][1] < minMinorCnt) {
            return results;
        }
        n += contig[0][0] = (int) OpenBitSet.intersectionCount(rMj, cMj);
        results.n = n;
        if (n < minCnt) {
            return results;
        }
        double rValue = LinkageDisequilibrium.calculateRSqr(contig[0][0], contig[1][0], contig[0][1], contig[1][1], minCnt);
        results.r2 = (float) rValue;
        if (Double.isNaN(rValue)) {
            return results;
        }
        results.dprime = (float) LinkageDisequilibrium.calculateDPrime(contig[0][0], contig[1][0], contig[0][1], contig[1][1], minCnt);;
        if (rValue < minR2) {
            return results;
        }
        double pValue = myFisherExact.getTwoTailedP(contig[0][0], contig[1][0], contig[0][1], contig[1][1]);
        results.p = (float) pValue;
        return results;
    }

    private int getRowFromIndex(long index) {

        int row = 0;
        int n = myAlignment.getSiteCount();
        int w = myWindowSize;

        if (myCurrDesign == testDesign.SlidingWindow && n > w + 1 && index >= w * (w + 1) / (double) 2) {
            row = (int) Math.ceil(((double) index + 1 - w * (w + 1) / 2 + w * w) / w);
        } else if (myCurrDesign == testDesign.SiteByAll) {
            if (index < myTestSite) {
                row = myTestSite;
            } else {
                row = (int) index + 1;
            }
        } else if (myCurrDesign == testDesign.SiteList) {

            int k = (int) Math.ceil((n - 1.5) - Math.sqrt((n - 1.5) * (n - 1.5) + 2 * (n - index - 2)));
            int m = n * (k + 1) - ((k + 1) * (k + 2) / 2) - 1;

            if (m - index > n - mySiteList[k] - 1) {
                row = mySiteList[k];
            } else {
                row = n - 1 - m + (int) index;
            }
        } else {
            row = (int) Math.ceil((Math.sqrt(8 * (index + 1) + 1) - 1) / 2);
        }

        return row;
    }

    private int getColFromIndex(long index) {

        int row = getRowFromIndex(index);
        int col = 0;
        int n = myAlignment.getSiteCount();
        int w = myWindowSize;

        if (myCurrDesign == testDesign.SlidingWindow && n > w + 1 && index >= w * (w + 1) / (double) 2) {
            col = (int) ((row - 1 - (double) w * (w + 1) / 2 - w * (row - w) + 1 + index));
        } else if (myCurrDesign == testDesign.SiteByAll) {
            if (index < myTestSite) {
                col = (int) index;
            } else {
                col = myTestSite;
            }
        } else if (myCurrDesign == testDesign.SiteList) {

            int k = (int) Math.ceil((n - 1.5) - Math.sqrt((n - 1.5) * (n - 1.5) + 2 * (n - index - 2)));
            int m = n * (k + 1) - ((k + 1) * (k + 2) / 2) - 1;

            if (row != mySiteList[k]) {
                col = mySiteList[k];
            } else {
                col = n - m + (int) index - 2;
                int yy = Arrays.binarySearch(mySiteList, row);
                int y = Arrays.binarySearch(mySiteList, col);
                while (yy + (y + 1) != 0) {
                    if (y < 0) {
                        y = -(y + 1);
                    }
                    col = col - (yy - y);
                    yy = y;
                    y = Arrays.binarySearch(mySiteList, col);

                }
            }
        } else {
            col = row - (row * (row + 1) / 2 - (int) index);
        }

        return col;
    }

    /**
     * Returns P-value estimate for a given pair of numSites. If there were only
     * 2 alleles at each locus, then the Fisher Exact P-value (one-tail) is
     * returned. If more states then the permuted Monte Carlo test is used.
     *
     * @param r is site 1
     * @param c is site 2
     * @return P-value
     */
    public double getPVal(int r, int c) {
        long key=getMapKey(r, c);
        LDResult result=(LDResult)myMapResults.get(key);
        if(result==null) return Float.NaN;
        return result.p;
    }

    /**
     * Get number of gametes included in LD calculations (after missing data was
     * excluded)
     *
     * @param r is site 1
     * @param c is site 2
     * @return number of gametes
     */
    public int getSampleSize(int r, int c) {
        long key=getMapKey(r, c);
        LDResult result=(LDResult)myMapResults.get(key);
        if(result==null) return 0;
        return result.n;
    }

    /**
     * Returns D' estimate for a given pair of numSites
     *
     * @param r is site 1
     * @param c is site 2
     * @return D'
     */
    public float getDPrime(int r, int c) {
        long key=getMapKey(r, c);
        LDResult result=(LDResult)myMapResults.get(key);
        if(result==null) return Float.NaN;
        return result.dprime;
    }

    /**
     * Returns r^2 estimate for a given pair of numSites
     *
     * @param r is site 1
     * @param c is site 2
     * @return r^2
     */
    public float getRSqr(int r, int c) {
        long key=getMapKey(r, c);
        LDResult result=(LDResult)myMapResults.get(key);
        if(result==null) return Float.NaN;
        return result.r2;
    }

    public int getX(int row) {
        return getColFromIndex(row);
    }

    public int getY(int row) {
        return getRowFromIndex(row);
    }

    /**
     * Returns the counts of the numSites in the alignment
     */
    public int getSiteCount() {
        return myAlignment.getSiteCount();
    }

    /**
     * Returns an annotated aligment if one was used for this LD this could be
     * used to access information of locus position
     */
    public Alignment getAlignment() {
        return myAlignment;
    }

    /**
     * Returns representation of the LD results as a string
     */
    public String toString() {
        String delimit = "\t";
        StringWriter sw = new StringWriter();
        Object[] colNames = getTableColumnNames();
        for (int j = 0; j < colNames.length; j++) {
            sw.write(colNames[j].toString());
            sw.write(delimit);
        }
        sw.write("\n");

        for (int r = 0; r < myTotalTests; r++) {
            Object[] theRow = getRow(r);
            for (int i = 0; i < theRow.length; i++) {
                sw.write(theRow[i].toString());
                sw.write(delimit);
            }
        }
        return sw.toString();
    }

    @Override
    public Object[] getTableColumnNames() {
        String[] annotatedLabels = null;
        if (myIsAccumulativeReport) {
            annotatedLabels = new String[]{"R2BinMin", "R2BinMax", "Count"};
        } else {
            annotatedLabels = new String[]{"Locus1", "Position1", "Site1",
                "NumberOfStates1", "States1", "Frequency1", "Locus2", "Position2",
                "Site2", "NumberOfStates2", "States2", "Frequency2", "Dist_bp", "R^2", "DPrime", "pDiseq", "N"};
        }
        return annotatedLabels;
    }

    @Override
    public Object[][] getTableData() {
        return getTableData(0, getRowCount() - 1);
    }

    @Override
    public Object[][] getTableData(int start, int end) {
        if ((start < 0) || (end >= getRowCount())) {
            throw new IndexOutOfBoundsException("getTableData: start: " + start + "  end: " + end);
        }
        if (end < start) {
            return null;
        }
        Object[][] temp = new Object[end - start + 1][];
        for (int i = start; i <= end; i++) {
            temp[i] = getRow(i);
        }

        return temp;

    }

    @Override
    public Object[] getRow(int row) {

        if (myIsAccumulativeReport) {
            Object[] data = new Object[3];
            if (row == myNumAccumulativeBins) {
                data[0] = Double.NaN;
                data[1] = Double.NaN;
                data[2] = Integer.valueOf(myAccumulativeRValueBins[(int) row]);
            } else {
                double start = myAccumulativeInterval * (double) row;
                data[0] = Double.valueOf(start);
                data[1] = Double.valueOf(start + myAccumulativeInterval);
                data[2] = Integer.valueOf(myAccumulativeRValueBins[(int) row]);
            }
            return data;
        } else {
            int labelOffset = 0;
            Object[] data = new Object[17];

            int r = getRowFromIndex(row);
            int c = getColFromIndex(row);

            String rState = myAlignment.getMajorAlleleAsString(r) + ":" + myAlignment.getMinorAlleleAsString(r);
            Integer rStr = Integer.valueOf(r);

            String cState = myAlignment.getMajorAlleleAsString(c) + ":" + myAlignment.getMinorAlleleAsString(c);
            Integer cStr = Integer.valueOf(c);

            data[labelOffset++] = myAlignment.getLocusName(r);
            data[labelOffset++] = Integer.valueOf(myAlignment.getPositionInLocus(r));
            data[labelOffset++] = rStr;

            data[labelOffset++] = IntegerTwo;
            data[labelOffset++] = rState;
            data[labelOffset++] = NotImplemented;
            data[labelOffset++] = myAlignment.getLocusName(c);
            data[labelOffset++] = Integer.valueOf(myAlignment.getPositionInLocus(c));
            data[labelOffset++] = cStr;

            data[labelOffset++] = IntegerTwo;
            data[labelOffset++] = cState;
            data[labelOffset++] = NotImplemented;
            if (myAlignment.getLocusName(r).equals(myAlignment.getLocusName(c))) {
                data[labelOffset++] = Integer.valueOf(Math.abs(myAlignment.getPositionInLocus(r) - myAlignment.getPositionInLocus(c)));
            } else {
                data[labelOffset++] = NA;
            }
            data[labelOffset++] = getRSqr(r, c);
            data[labelOffset++] = getDPrime(r, c);
            data[labelOffset++] = getPVal(r, c);
            data[labelOffset++] = getSampleSize(r, c);

            return data;
        }

    }

    @Override
    public String getTableTitle() {
        return "Linkage Disequilibrium";
    }

    @Override
    public int getColumnCount() {
        return getTableColumnNames().length;
    }
    @Override
    public int getRowCount() {
        if (myIsAccumulativeReport) {
            return myNumAccumulativeBins + 1;
        } else {
            return (int) myTotalTests;
        }
    }
   @Override
    public int getElementCount() {
        return getRowCount() * getColumnCount();
    }
    
    @Override
    public Object getValueAt(int row, int col) {
        return getRow(row)[col];
    }

    protected void fireProgress(int percent) {
        if (percent < 0) {
            percent = -percent;
        }
        if (myListener != null) {
            myListener.progress(percent, null);
        }
    }
}

/**
 * Container class for reporting LD results
 * 
 * @author edbuckler
 */
class LDResult {
    public float r2=Float.NaN;
    public float dprime=Float.NaN;
    public float p=Float.NaN;
    public int n=0;

    @Override
    public String toString() {
        return "LDResult{" + "r2=" + r2 + ", dprime=" + dprime + ", p=" + p + ", n=" + n + '}';
    }
}
