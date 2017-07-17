/*
 * BitAlignmentHDF5
 */
package net.maizegenetics.pal.alignment;

import ch.systemsx.cisd.base.mdarray.MDArray;
import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import java.util.LinkedHashMap;
import java.util.Map;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.ProgressListener;
import net.maizegenetics.util.UnmodifiableBitSet;
import org.apache.log4j.Logger;

/**
 * This data alignment reads from HDF5 files
 *
 * @author Ed & Terry
 */
public class BitAlignmentHDF5 extends AbstractAlignment {

    private static final long serialVersionUID = -5197800047652332969L;
    private static final Logger myLogger = Logger.getLogger(BitAlignmentHDF5.class);
    private static final int NUM_UNITS_TO_CACHE_ON_GET = 64;
    private static final int MAX_CACHE_SIZE = 200;
    private final IHDF5Reader myHDF5;
    private final int myNumDataRows;
    private final int myNumSBitWords;
    private final int myNumTBitWords;
    private static final String SBIT_HDF5_PATH = HapMapHDF5Constants.SBIT + "/";
    private static final String TBIT_HDF5_PATH = HapMapHDF5Constants.TBIT + "/";
    private final Map<Integer, OpenBitSet[]> myCachedSites = new LinkedHashMap<Integer, OpenBitSet[]>((MAX_CACHE_SIZE * 3) / 2) {
        @Override
        protected boolean removeEldestEntry(Map.Entry eldest) {
            return size() > MAX_CACHE_SIZE;
        }
    };
    private final Map<Integer, OpenBitSet[]> myCachedTaxa = new LinkedHashMap<Integer, OpenBitSet[]>((MAX_CACHE_SIZE * 3) / 2) {
        @Override
        protected boolean removeEldestEntry(Map.Entry eldest) {
            return size() > MAX_CACHE_SIZE;
        }
    };
    private boolean myOptimizeSiteIteration = true;

    protected BitAlignmentHDF5(IHDF5Reader hdf5, IdGroup idGroup, byte[][] alleles, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {
        super(alleles, idGroup, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);
        myHDF5 = hdf5;
        myNumSBitWords = myHDF5.getIntAttribute(HapMapHDF5Constants.DEFAULT_ATTRIBUTES_PATH, HapMapHDF5Constants.NUM_SBIT_WORDS);
        myNumTBitWords = myHDF5.getIntAttribute(HapMapHDF5Constants.DEFAULT_ATTRIBUTES_PATH, HapMapHDF5Constants.NUM_TBIT_WORDS);
        if (retainsRareAlleles()) {
            myNumDataRows = myMaxNumAlleles + 1;
        } else {
            myNumDataRows = myMaxNumAlleles;
        }
    }

    public static BitAlignmentHDF5 getInstance(String filename) {
        return getInstance(filename, true);
    }

    public static BitAlignmentHDF5 getInstance(String filename, boolean iterateSites) {
        IHDF5Reader reader = HDF5Factory.openForReading(filename);

        String[] taxa = reader.readStringArray(HapMapHDF5Constants.TAXA);
        IdGroup idgroup = new SimpleIdGroup(taxa);

        byte[][] alleles = reader.readByteMatrix(HapMapHDF5Constants.ALLELES);

        MDArray<String> alleleStatesMDArray = reader.readStringMDArray(HapMapHDF5Constants.ALLELE_STATES);
        int[] dimensions = alleleStatesMDArray.dimensions();
        int numEncodings = dimensions[0];
        int numStates = dimensions[1];
        String[][] alleleStates = new String[numEncodings][numStates];
        for (int e = 0; e < numEncodings; e++) {
            for (int s = 0; s < numStates; s++) {
                alleleStates[e][s] = alleleStatesMDArray.get(e, s);
            }
        }

        int[] variableSites = reader.readIntArray(HapMapHDF5Constants.POSITIONS);

        int maxNumAlleles = reader.getIntAttribute(HapMapHDF5Constants.DEFAULT_ATTRIBUTES_PATH, HapMapHDF5Constants.MAX_NUM_ALLELES);

        boolean retainRare = reader.getBooleanAttribute(HapMapHDF5Constants.DEFAULT_ATTRIBUTES_PATH, HapMapHDF5Constants.RETAIN_RARE_ALLELES);

        String[] lociStrings = reader.readStringArray(HapMapHDF5Constants.LOCI);
        int numLoci = lociStrings.length;
        Locus[] loci = new Locus[numLoci];
        for (int i = 0; i < numLoci; i++) {
            loci[i] = new Locus(lociStrings[i]);
        }

        int[] lociOffsets = reader.readIntArray(HapMapHDF5Constants.LOCUS_OFFSETS);

        String[] snpIds = reader.readStringArray(HapMapHDF5Constants.SNP_IDS);

        BitAlignmentHDF5 result = null;
        if (NucleotideAlignmentConstants.isNucleotideEncodings(alleleStates)) {
            result = new BitNucleotideAlignmentHDF5(reader, idgroup, alleles, null, null, variableSites, maxNumAlleles, loci, lociOffsets, snpIds, retainRare);
        } else if (alleleStates.length == 1) {
            result = new BitAlignmentHDF5(reader, idgroup, alleles, null, null, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIds, retainRare);
        } else {
            result = new BitTextAlignmentHDF5(reader, idgroup, alleles, null, null, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIds, retainRare);
        }

        result.setOptimizeSiteIteration(iterateSites);
        return result;

    }

    private OpenBitSet[] getCachedSite(int site) {

        OpenBitSet[] result = myCachedSites.get(site);
        if (result != null) {
            return result;
        }

        int startSite = (site / NUM_UNITS_TO_CACHE_ON_GET) * NUM_UNITS_TO_CACHE_ON_GET;
        //int startSite = site & 0xFFFFFFC0;

        int numToRetrieve = Math.min(NUM_UNITS_TO_CACHE_ON_GET, getSiteCount() - startSite);
        OpenBitSet[][] block = new OpenBitSet[numToRetrieve][myNumDataRows];

        for (int aNum = 0; aNum < myNumDataRows; aNum++) {

            long[][] temp;
            synchronized (myHDF5) {
                temp = myHDF5.readLongMatrixBlockWithOffset(SBIT_HDF5_PATH + aNum, numToRetrieve, myNumSBitWords, startSite, 0);
            }
            for (int i = 0; i < numToRetrieve; i++) {
                block[i][aNum] = new OpenBitSet(temp[i], myNumSBitWords);
            }

        }

        for (int i = 0; i < numToRetrieve; i++) {
            myCachedSites.put(startSite + i, block[i]);
        }

        new Thread(new LookAheadSiteCache(startSite + numToRetrieve)).start();

        return block[site - startSite];

    }
    private static final int NUM_UNITS_TO_CACHE_IN_ADVANCE = 64;

    private class LookAheadSiteCache implements Runnable {

        private final int mySite;

        public LookAheadSiteCache(int site) {
            mySite = site;
        }

        @Override
        public void run() {

            int numToRetrieve = Math.min(NUM_UNITS_TO_CACHE_IN_ADVANCE, getSiteCount() - mySite);
            if (numToRetrieve == 0) {
                return;
            }
            OpenBitSet[][] block = new OpenBitSet[numToRetrieve][myNumDataRows];

            for (int aNum = 0; aNum < myNumDataRows; aNum++) {

                long[][] temp;
                synchronized (myHDF5) {
                    temp = myHDF5.readLongMatrixBlockWithOffset(SBIT_HDF5_PATH + aNum, numToRetrieve, myNumSBitWords, mySite, 0);
                }
                for (int i = 0; i < numToRetrieve; i++) {
                    block[i][aNum] = new OpenBitSet(temp[i], myNumSBitWords);
                }

            }

            for (int i = 0; i < numToRetrieve; i++) {
                myCachedSites.put(mySite + i, block[i]);
            }

        }
    }

    protected OpenBitSet[] getCachedTaxon(int taxon) {

        OpenBitSet[] result = myCachedTaxa.get(taxon);
        if (result != null) {
            return result;
        }

        int startTaxon = (taxon / NUM_UNITS_TO_CACHE_ON_GET) * NUM_UNITS_TO_CACHE_ON_GET;
        //int startTaxon = taxon & 0xFFFFFFC0;

        int numToRetrieve = Math.min(NUM_UNITS_TO_CACHE_ON_GET, getSequenceCount() - startTaxon);
        OpenBitSet[][] block = new OpenBitSet[numToRetrieve][myNumDataRows];

        for (int aNum = 0; aNum < myNumDataRows; aNum++) {

            long[][] temp;
            synchronized (myHDF5) {
                temp = myHDF5.readLongMatrixBlockWithOffset(TBIT_HDF5_PATH + aNum, numToRetrieve, myNumTBitWords, startTaxon, 0);
            }
            for (int i = 0; i < numToRetrieve; i++) {
                block[i][aNum] = new OpenBitSet(temp[i], myNumTBitWords);
            }

        }

        for (int i = 0; i < numToRetrieve; i++) {
            myCachedTaxa.put(startTaxon + i, block[i]);
        }

        new Thread(new LookAheadTaxonCache(startTaxon + numToRetrieve)).start();

        return block[taxon - startTaxon];
    }

    private class LookAheadTaxonCache implements Runnable {

        private final int myTaxon;

        public LookAheadTaxonCache(int taxon) {
            myTaxon = taxon;
        }

        @Override
        public void run() {

            int numToRetrieve = Math.min(NUM_UNITS_TO_CACHE_IN_ADVANCE, getSequenceCount() - myTaxon);
            if (numToRetrieve == 0) {
                return;
            }
            OpenBitSet[][] block = new OpenBitSet[numToRetrieve][myNumDataRows];

            for (int aNum = 0; aNum < myNumDataRows; aNum++) {

                long[][] temp;
                synchronized (myHDF5) {
                    temp = myHDF5.readLongMatrixBlockWithOffset(TBIT_HDF5_PATH + aNum, numToRetrieve, myNumTBitWords, myTaxon, 0);
                }
                for (int i = 0; i < numToRetrieve; i++) {
                    block[i][aNum] = new OpenBitSet(temp[i], myNumTBitWords);
                }

            }

            for (int i = 0; i < numToRetrieve; i++) {
                myCachedTaxa.put(myTaxon + i, block[i]);
            }

        }
    }

    @Override
    public byte getBase(int taxon, int site) {
        byte[] temp = getBaseArray(taxon, site);
        return (byte) ((temp[0] << 4) | temp[1]);
    }

    @Override
    public byte[] getBaseArray(int taxon, int site) {
        if (myOptimizeSiteIteration) {
            return getBaseArrayS(taxon, site);
        } else {
            return getBaseArrayT(taxon, site);
        }
    }

    public byte[] getBaseArrayS(int taxon, int site) {
        OpenBitSet[] sBitData = getCachedSite(site);
        byte[] result = new byte[2];
        result[0] = Alignment.UNKNOWN_ALLELE;
        result[1] = Alignment.UNKNOWN_ALLELE;
        try {
            int count = 0;
            for (int i = 0; i < myMaxNumAlleles; i++) {
                if (sBitData[i].fastGet(taxon)) {
                    if (count == 0) {
                        result[1] = myAlleles[site][i];
                    }
                    result[count++] = myAlleles[site][i];
                }
            }

            // Check For Rare Allele
            if (retainsRareAlleles() && sBitData[myMaxNumAlleles].fastGet(taxon)) {
                if (count == 0) {
                    result[1] = Alignment.RARE_ALLELE;
                }
                result[count] = Alignment.RARE_ALLELE;
            }

        } catch (IndexOutOfBoundsException e) {
            throw new IllegalStateException("BitAlignmentHDF5: getBaseArrayS: bit sets indicate more than two alleles for taxon: " + taxon + "   site: " + site);
        }
        return result;
    }

    public byte[] getBaseArrayT(int taxon, int site) {
        OpenBitSet[] tBitData = getCachedTaxon(taxon);
        byte[] result = new byte[2];
        result[0] = Alignment.UNKNOWN_ALLELE;
        result[1] = Alignment.UNKNOWN_ALLELE;
        try {
            int count = 0;
            for (int i = 0; i < myMaxNumAlleles; i++) {
                if (tBitData[i].fastGet(site)) {
                    if (count == 0) {
                        result[1] = myAlleles[site][i];
                    }
                    result[count++] = myAlleles[site][i];
                }
            }

            // Check For Rare Allele
            if (retainsRareAlleles() && tBitData[myMaxNumAlleles].fastGet(site)) {
                if (count == 0) {
                    result[1] = Alignment.RARE_ALLELE;
                }
                result[count] = Alignment.RARE_ALLELE;
            }

        } catch (IndexOutOfBoundsException e) {
            throw new IllegalStateException("BitNucleotideAlignmentHDF5: getBaseArrayT: bit sets indicate more than two alleles for taxon: " + taxon + "   site: " + site);
        }
        return result;
    }

    @Override
    public BitSet getAllelePresenceForAllTaxa(int site, int alleleNumber) {
        OpenBitSet[] sBitData = getCachedSite(site);
        return UnmodifiableBitSet.getInstance(sBitData[alleleNumber]);
    }

    @Override
    public BitSet getAllelePresenceForAllSites(int taxon, int alleleNumber) {
        OpenBitSet[] tBitData = getCachedTaxon(taxon);
        return UnmodifiableBitSet.getInstance(tBitData[alleleNumber]);
    }

    @Override
    public long[] getAllelePresenceForSitesBlock(int taxon, int alleleNumber, int startBlock, int endBlock) {
        OpenBitSet[] tBitData = getCachedTaxon(taxon);
        long[] result = new long[endBlock - startBlock];
        System.arraycopy(tBitData[alleleNumber].getBits(), startBlock, result, 0, endBlock - startBlock);
        return result;
    }

    @Override
    public int getTotalGametesNotMissing(int site) {
        OpenBitSet[] sBitData = getCachedSite(site);
        OpenBitSet temp = new OpenBitSet(getSequenceCount());
        for (int i = 0; i < myNumDataRows; i++) {
            temp.or(sBitData[i]);
        }
        return ((int) temp.cardinality()) * 2;

    }

    @Override
    public int getTotalNotMissing(int site) {

        OpenBitSet[] sBitData = getCachedSite(site);

        OpenBitSet temp = new OpenBitSet(getSequenceCount());
        for (int i = 0; i < myNumDataRows; i++) {
            temp.or(sBitData[i]);
        }
        return (int) temp.cardinality();

    }

    @Override
    public int getMinorAlleleCount(int site) {
        OpenBitSet[] sBitData = getCachedSite(site);
        if ((myMaxNumAlleles < 2) || (myAlleles[site][1] == Alignment.UNKNOWN_ALLELE)) {
            return 0;
        }

        OpenBitSet temp = new OpenBitSet(getSequenceCount());
        for (int i = 0; i < myNumDataRows; i++) {
            if (i != 1) {
                temp.or(sBitData[i]);
            }
        }
        temp.flip(0, temp.size());
        temp.and(sBitData[1]);

        return (int) temp.cardinality() + (int) sBitData[1].cardinality();

    }

    @Override
    public double getMinorAlleleFrequency(int site) {
        int minorAlleleCount = getMinorAlleleCount(site);
        if (minorAlleleCount == 0) {
            return 0.0;
        }
        return (double) minorAlleleCount / (double) getTotalGametesNotMissing(site);
    }

    @Override
    public double getMajorAlleleFrequency(int site) {
        int majorAlleleCount = getMajorAlleleCount(site);
        if (majorAlleleCount == 0) {
            return 0.0;
        }
        return (double) majorAlleleCount / (double) getTotalGametesNotMissing(site);
    }

    @Override
    public int getMajorAlleleCount(int site) {
        OpenBitSet[] sBitData = getCachedSite(site);
        if (myAlleles[site][0] == Alignment.UNKNOWN_ALLELE) {
            return 0;
        }

        OpenBitSet temp = new OpenBitSet(getSequenceCount());
        for (int i = 1; i < myNumDataRows; i++) {
            temp.or(sBitData[i]);
        }
        temp.flip(0, temp.size());
        temp.and(sBitData[0]);

        return (int) temp.cardinality() + (int) sBitData[0].cardinality();
    }

    @Override
    public boolean isHeterozygous(int taxon, int site) {
        OpenBitSet[] sBitData = getCachedSite(site);
        int count = 0;
        for (int i = 0; i < myNumDataRows; i++) {
            if (sBitData[i].fastGet(taxon)) {
                count++;
                if (count == 2) {
                    return true;
                }
            }
        }
        return false;
    }

    @Override
    public int getHeterozygousCount(int site) {
        OpenBitSet[] sBitData = getCachedSite(site);
        int result = 0;
        for (int i = 0; i < myNumDataRows; i++) {
            for (int j = i + 1; j < myNumDataRows; j++) {
                result += (int) OpenBitSet.intersectionCount(sBitData[i], sBitData[j]);
            }
        }
        return result;
    }

    @Override
    public boolean isPolymorphic(int site) {
        OpenBitSet[] sBitData = getCachedSite(site);
        boolean nonZero = false;
        for (int i = 0; i < myNumDataRows; i++) {
            int numTaxa = (int) sBitData[i].cardinality();
            if (numTaxa != 0) {
                if (nonZero) {
                    return true;
                }
                nonZero = true;
            }
        }
        return false;
    }

    @Override
    public Object[][] getDiploidCounts() {

        if (myAlleleStates.length != 1) {
            return super.getDiploidCounts();
        }

        long[][] counts = new long[16][16];
        for (int site = 0; site < myNumSites; site++) {
            OpenBitSet[] sBitData = getCachedSite(site);
            for (int i = 0; i < myMaxNumAlleles; i++) {
                byte indexI = myAlleles[site][i];
                counts[indexI][indexI] += sBitData[i].cardinality();
                for (int j = i + 1; j < myMaxNumAlleles; j++) {
                    byte indexJ = myAlleles[site][j];
                    long ijHet = OpenBitSet.intersectionCount(sBitData[i], sBitData[j]);
                    if (indexI < indexJ) {
                        counts[indexI][indexJ] += ijHet;
                    } else {
                        counts[indexJ][indexI] += ijHet;
                    }
                    counts[indexI][indexI] -= ijHet;
                    counts[indexJ][indexJ] -= ijHet;
                }
            }
        }

        int numAlleles = 0;
        long unknownCount = (long) getSequenceCount() * (long) myNumSites;
        for (byte x = 0; x < 16; x++) {
            for (byte y = x; y < 16; y++) {
                if (counts[x][y] != 0) {
                    numAlleles++;
                    unknownCount -= counts[x][y];
                }
            }
        }

        if (unknownCount > 0) {
            numAlleles++;
        }

        Object[][] result = new Object[2][numAlleles];
        int nextResult = 0;
        for (byte x = 0; x < 16; x++) {
            for (byte y = x; y < 16; y++) {
                if (counts[x][y] != 0) {
                    byte value = (byte) ((x << 4) | y);
                    result[0][nextResult] = getDiploidAsString(0, value);
                    result[1][nextResult++] = counts[x][y];
                }
            }
        }

        if (unknownCount > 0) {
            result[0][nextResult] = getDiploidAsString(0, UNKNOWN_DIPLOID_ALLELE);
            result[1][nextResult] = unknownCount;
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0; k < numAlleles - 1; k++) {

                if ((Long) result[1][k] < (Long) result[1][k + 1]) {

                    Object temp = result[0][k];
                    result[0][k] = result[0][k + 1];
                    result[0][k + 1] = temp;

                    Object tempCount = result[1][k];
                    result[1][k] = result[1][k + 1];
                    result[1][k + 1] = tempCount;

                    change = true;
                }
            }

        }

        return result;
    }

    @Override
    public Object[][] getDiploidssSortedByFrequency(int site) {
        OpenBitSet[] sBitData = getCachedSite(site);
        if (myAlleleStates.length != 1) {
            return super.getDiploidssSortedByFrequency(site);
        }

        int[][] counts = new int[16][16];
        for (int i = 0; i < myMaxNumAlleles; i++) {
            byte indexI = myAlleles[site][i];
            counts[indexI][indexI] += (int) sBitData[i].cardinality();
            for (int j = i + 1; j < myMaxNumAlleles; j++) {
                byte indexJ = myAlleles[site][j];
                int ijHet = (int) OpenBitSet.intersectionCount(sBitData[i], sBitData[j]);
                if (indexI < indexJ) {
                    counts[indexI][indexJ] += ijHet;
                } else {
                    counts[indexJ][indexI] += ijHet;
                }
                counts[indexI][indexI] -= ijHet;
                counts[indexJ][indexJ] -= ijHet;
            }
        }

        int numAlleles = 0;
        int unknownCount = getSequenceCount();
        for (byte x = 0; x < 16; x++) {
            for (byte y = x; y < 16; y++) {
                if (counts[x][y] != 0) {
                    numAlleles++;
                    unknownCount -= counts[x][y];
                }
            }
        }

        if (unknownCount > 0) {
            numAlleles++;
        }

        Object[][] result = new Object[2][numAlleles];
        int nextResult = 0;
        for (byte x = 0; x < 16; x++) {
            for (byte y = x; y < 16; y++) {
                if (counts[x][y] != 0) {
                    byte value = (byte) ((x << 4) | y);
                    result[0][nextResult] = getDiploidAsString(0, value);
                    result[1][nextResult++] = counts[x][y];
                }
            }
        }

        if (unknownCount > 0) {
            result[0][nextResult] = getDiploidAsString(0, UNKNOWN_DIPLOID_ALLELE);
            result[1][nextResult] = unknownCount;
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0; k < numAlleles - 1; k++) {

                if ((Integer) result[1][k] < (Integer) result[1][k + 1]) {

                    Object temp = result[0][k];
                    result[0][k] = result[0][k + 1];
                    result[0][k + 1] = temp;

                    Object tempCount = result[1][k];
                    result[1][k] = result[1][k + 1];
                    result[1][k + 1] = tempCount;

                    change = true;
                }
            }

        }

        return result;

    }

    @Override
    public int[][] getAllelesSortedByFrequency(int site) {

        OpenBitSet[] sBitData = getCachedSite(site);
        int[] counts = new int[16];
        for (int i = 0; i < myNumDataRows; i++) {
            byte indexI;
            if ((retainsRareAlleles()) && (i == myMaxNumAlleles)) {
                indexI = Alignment.RARE_ALLELE;
            } else {
                indexI = myAlleles[site][i];
            }
            counts[indexI] += (int) sBitData[i].cardinality() * 2;
            for (int j = i + 1; j < myNumDataRows; j++) {
                byte indexJ;
                if ((retainsRareAlleles()) && (j == myMaxNumAlleles)) {
                    indexJ = Alignment.RARE_ALLELE;
                } else {
                    indexJ = myAlleles[site][j];
                }
                int ijHet = (int) OpenBitSet.intersectionCount(sBitData[i], sBitData[j]);
                counts[indexI] -= ijHet;
                counts[indexJ] -= ijHet;
            }
        }

        int numAlleles = 0;
        for (byte x = 0; x < Alignment.UNKNOWN_ALLELE; x++) {
            if (counts[x] != 0) {
                numAlleles++;
            }
        }

        int current = 0;
        int[][] result = new int[2][numAlleles];
        for (byte x = 0; x < Alignment.UNKNOWN_ALLELE; x++) {
            if (counts[x] != 0) {
                result[0][current] = x;
                result[1][current++] = counts[x];
            }
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0; k < numAlleles - 1; k++) {

                if (result[1][k] < result[1][k + 1]) {

                    int temp = result[0][k];
                    result[0][k] = result[0][k + 1];
                    result[0][k + 1] = temp;

                    int tempCount = result[1][k];
                    result[1][k] = result[1][k + 1];
                    result[1][k + 1] = tempCount;

                    change = true;
                }
            }

        }

        return result;

    }

    @Override
    public boolean isSBitFriendly() {
        return true;
    }

    @Override
    public boolean isTBitFriendly() {
        return true;
    }

    @Override
    public int getTotalNumAlleles() {
        return myNumDataRows;
    }

    @Override
    public int getTotalGametesNotMissingForTaxon(int taxon) {

        OpenBitSet[] tBitData = getCachedTaxon(taxon);

        OpenBitSet temp = new OpenBitSet(getSiteCount());
        for (int i = 0; i < myNumDataRows; i++) {
            temp.or(tBitData[i]);
        }
        return ((int) temp.cardinality()) * 2;

    }

    @Override
    public int getHeterozygousCountForTaxon(int taxon) {

        OpenBitSet[] tBitData = getCachedTaxon(taxon);

        int result = 0;
        for (int i = 0; i < myNumDataRows; i++) {
            for (int j = i + 1; j < myNumDataRows; j++) {
                result += (int) OpenBitSet.intersectionCount(tBitData[i], tBitData[j]);
            }
        }
        return result;

    }

    @Override
    public int getTotalNotMissingForTaxon(int taxon) {

        OpenBitSet[] tBitData = getCachedTaxon(taxon);

        OpenBitSet temp = new OpenBitSet(getSiteCount());
        for (int i = 0; i < myNumDataRows; i++) {
            temp.or(tBitData[i]);
        }
        return (int) temp.cardinality();

    }

    @Override
    public void optimizeForTaxa(ProgressListener listener) {
        // do nothing
    }

    @Override
    public void optimizeForSites(ProgressListener listener) {
        // do nothing
    }

    /**
     * Sets whether this BitAlignmentHDF5 caches bit sets to optimize for site
     * iteration or taxa iteration.
     *
     * @param optimizeSiteIteration true to optimize for site iteration or false
     * to optimize for taxa iteration
     */
    public void setOptimizeSiteIteration(boolean optimizeSiteIteration) {
        myOptimizeSiteIteration = optimizeSiteIteration;
    }
}
