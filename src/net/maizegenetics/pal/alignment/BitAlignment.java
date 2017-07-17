/*
 * BitAlignment
 */
package net.maizegenetics.pal.alignment;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.BitUtil;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.ProgressListener;
import net.maizegenetics.util.UnmodifiableBitSet;
import org.apache.log4j.Logger;

/**
 * Bit Alignment holds alleles or taxa in BitSets. TBit holds the data organized
 * by taxon, while SBit is organized by site. TBit is efficient for comparison
 * for large numbers of sites between taxa (e.g. Genetic distance or
 * imputation). Sbit is efficient for comparisons between sites (e.g LD, and
 * better for accessing for GWAS but GWAS efficiency is determined by the linear
 * algebra).
 *
 * TODO: Need to create cache system for a long at a time.
 *
 * @author terry
 */
public class BitAlignment extends AbstractAlignment {
    
    private static final long serialVersionUID = -5197800047652332969L;
    private static final Logger myLogger = Logger.getLogger(BitAlignment.class);
    private BitSet[][] mySBitData;
    private BitSet[][] myTBitData;
    private int myNumDataRows;
    
    protected BitAlignment(Alignment a, int maxNumAlleles, boolean retainRareAlleles, boolean isSBit) {
        super(a, maxNumAlleles, retainRareAlleles);
        if (isSBit) {
            loadSBitAlleles(a, null);
        } else {
            loadTBitAlleles(a, null);
        }
    }
    
    protected BitAlignment(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isSBit) {
        super(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);
        if (isSBit) {
            loadSBitAlleles(data);
        } else {
            loadTBitAlleles(data);
        }
    }
    
    protected BitAlignment(IdGroup idGroup, byte[][] alleles, BitSet[][] data, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isSBit) {
        super(alleles, idGroup, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);
        if (isSBit) {
            mySBitData = data;
        } else {
            myTBitData = data;
        }
        myNumDataRows = myMaxNumAlleles;
        if (retainsRareAlleles()) {
            myNumDataRows++;
        }
    }
    
    public static Alignment getInstance(Alignment a, boolean isSBit) {
        return BitAlignment.getInstance(a, a.getMaxNumAlleles(), a.retainsRareAlleles(), isSBit);
    }
    
    public static Alignment getInstance(Alignment a, int maxNumAlleles, boolean retainRareAlleles, boolean isSBit) {
        
        if ((a instanceof BitAlignment) && (a.getMaxNumAlleles() == maxNumAlleles) && (a.retainsRareAlleles() == retainRareAlleles)) {
            return (BitAlignment) a;
        }
        
        String[][] alleleStates = a.getAlleleEncodings();
        if ((alleleStates == null) || (alleleStates.length == 0)) {
            throw new IllegalStateException("BitAlignment: init: allele states should not be empty.");
        }
        
        boolean isNucleotide = false;
        if (alleleStates.length == 1) {
            isNucleotide = true;
            if (alleleStates[0].length == NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES[0].length) {
                for (int i = 0; i < alleleStates.length; i++) {
                    if (!alleleStates[0][i].equals(NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES[0][i])) {
                        isNucleotide = false;
                    }
                }
            }
            
        }
        
        if (isNucleotide) {
            return new BitNucleotideAlignment(a, maxNumAlleles, retainRareAlleles, isSBit);
        } else if (alleleStates.length == 1) {
            return new BitAlignment(a, maxNumAlleles, retainRareAlleles, isSBit);
        } else {
            return new BitTextAlignment(a, maxNumAlleles, retainRareAlleles, isSBit);
        }
        
    }
    
    public static Alignment getInstance(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isSBit) {
        if ((alleleStates == null) || (alleleStates.length == 0)) {
            throw new IllegalArgumentException("BitAlignment: init: allele states can not be empty.");
        }
        if (alleleStates.length == 1) {
            return new BitAlignment(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isSBit);
        } else {
            return new BitTextAlignment(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isSBit);
        }
    }
    
    public static Alignment getNucleotideInstance(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isSBit) {
        return new BitNucleotideAlignment(idGroup, data, map, reference, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isSBit);
    }
    
    public static Alignment getNucleotideInstance(IdGroup idGroup, byte[][] alleles, BitSet[][] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isSBit) {
        return new BitNucleotideAlignment(idGroup, alleles, data, map, reference, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isSBit);
    }
    
    public static Alignment getNucleotideInstance(IdGroup idGroup, String[][] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isSBit) {
        
        if ((maxNumAlleles < 1) || (maxNumAlleles > NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES)) {
            throw new IllegalArgumentException("BitAlignment: getNucleotideInstance: max number of alleles must be between 1 and " + NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES + " inclusive: " + maxNumAlleles);
        }
        
        if ((data == null) || (data.length == 0)) {
            throw new IllegalArgumentException("BitAlignment: getNucleotideInstance: data can not be empty.");
        }
        
        if (data.length != idGroup.getIdCount()) {
            throw new IllegalArgumentException("BitAlignment: getNucleotideInstance: data rows not equal to number of identifers.");
        }
        
        byte[][] dataBytes = AlignmentUtils.getDataBytes(data, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES, NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES);
        
        return BitAlignment.getNucleotideInstance(idGroup, dataBytes, map, reference, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isSBit);
        
    }
    
    public static Alignment getNucleotideInstance(IdGroup idGroup, String[] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isSBit) {
        
        if ((maxNumAlleles < 1) || (maxNumAlleles > NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES)) {
            throw new IllegalArgumentException("BitAlignment: getNucleotideInstance: max number of alleles must be between 1 and " + NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES + " inclusive: " + maxNumAlleles);
        }
        
        if ((data == null) || (data.length == 0)) {
            throw new IllegalArgumentException("BitAlignment: getNucleotideInstance: data can not be empty.");
        }
        
        if (data.length != idGroup.getIdCount()) {
            throw new IllegalArgumentException("BitAlignment: getNucleotideInstance: data rows not equal to number of identifers.");
        }
        
        byte[][] dataBytes = AlignmentUtils.getDataBytes(data);
        
        return BitAlignment.getNucleotideInstance(idGroup, dataBytes, map, reference, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isSBit);
        
    }
    
    public static Alignment getNucleotideInstance(String[] data, boolean isSBit) {
        return BitAlignment.getNucleotideInstance(new SimpleIdGroup(data.length, true), data, null, null, null,
                TasselPrefs.getAlignmentMaxAllelesToRetain(), new Locus[]{Locus.UNKNOWN}, new int[]{0},
                null, TasselPrefs.getAlignmentRetainRareAlleles(), isSBit);
    }
    
    public static Alignment getInstance(IdGroup idGroup, String[][] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isSBit) {
        
        if ((maxNumAlleles < 1) || (maxNumAlleles > 14)) {
            throw new IllegalArgumentException("BitAlignment: getInstance: max number of alleles must be between 1 and 14 inclusive: " + maxNumAlleles);
        }
        
        if ((data == null) || (data.length == 0)) {
            throw new IllegalArgumentException("BitAlignment: getInstance: data can not be empty.");
        }
        
        if (data.length != idGroup.getIdCount()) {
            throw new IllegalArgumentException("BitAlignment: getInstance: data rows not equal to number of identifers.");
        }
        
        String[][] alleleStates = AlignmentUtils.getAlleleStates(data, maxNumAlleles);
        
        byte[][] dataBytes = AlignmentUtils.getDataBytes(data, alleleStates, maxNumAlleles);
        
        return BitAlignment.getInstance(idGroup, dataBytes, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isSBit);
        
    }

    /**
     * Create new alignment from given alignment that only retains two alleles
     * per site and only homozygous values. Heterzygous values are changed to
     * Unknown.
     *
     * @param a alignment
     * @param isSBit whether new alignment optimized for sites
     *
     * @return homozygous alignment
     */
    public static Alignment getHomozygousNucleotideInstance(Alignment a, boolean isSBit) {
        
        int numSites = a.getSiteCount();
        int numTaxa = a.getSequenceCount();
        
        int numAlleles = 2;
        byte[][] alleles = new byte[numSites][numAlleles];
        for (int s = 0; s < numSites; s++) {
            byte[] current = a.getAlleles(s);
            alleles[s][0] = current[0];
            alleles[s][1] = current[1];
        }
        
        BitSet[][] theData;
        int numDataRows = a.getTotalNumAlleles();
        if (isSBit) {
            a = AlignmentUtils.optimizeForSites(a);
            theData = new OpenBitSet[numAlleles][numSites];
            for (int s = 0; s < numSites; s++) {
                BitSet homo = new OpenBitSet(numTaxa);
                for (int al = 0; al < numDataRows; al++) {
                    BitSet current = a.getAllelePresenceForAllTaxa(s, al);
                    homo.xor(current);
                }
                theData[0][s] = new OpenBitSet(a.getAllelePresenceForAllTaxa(s, 0));
                theData[0][s].and(homo);
                theData[1][s] = new OpenBitSet(a.getAllelePresenceForAllTaxa(s, 1));
                theData[1][s].and(homo);
            }
        } else {
            a = AlignmentUtils.optimizeForTaxa(a);
            theData = new OpenBitSet[numAlleles][numTaxa];
            for (int t = 0; t < numTaxa; t++) {
                BitSet homo = new OpenBitSet(numSites);
                for (int al = 0; al < numDataRows; al++) {
                    BitSet current = a.getAllelePresenceForAllSites(t, al);
                    homo.xor(current);
                }
                theData[0][t] = new OpenBitSet(a.getAllelePresenceForAllSites(t, 0));
                theData[0][t].and(homo);
                theData[1][t] = new OpenBitSet(a.getAllelePresenceForAllSites(t, 1));
                theData[1][t].and(homo);
            }
        }
        
        return getNucleotideInstance(a.getIdGroup(), alleles, theData, a.getGeneticMap(), a.getReference(), a.getPhysicalPositions(), numAlleles, a.getLoci(), a.getLociOffsets(), a.getSNPIDs(), false, isSBit);
        
    }
    
    private void loadSBitAlleles(byte[][] data) {
        
        myNumDataRows = myMaxNumAlleles;
        if (retainsRareAlleles()) {
            myNumDataRows++;
        }
        int numSeqs = getSequenceCount();
        mySBitData = new OpenBitSet[myNumDataRows][myNumSites];
        for (int al = 0; al < myNumDataRows; al++) {
            for (int s = 0; s < myNumSites; s++) {
                mySBitData[al][s] = new OpenBitSet(numSeqs);
            }
        }
        ExecutorService pool = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
        for (int s = 0; s < myNumSites; s++) {
            pool.execute(new ProcessSite(data, mySBitData, s));
        }
        
        try {
            pool.shutdown();
            if (!pool.awaitTermination(600, TimeUnit.SECONDS)) {
                throw new IllegalStateException("BitAlignment: loadSBitAlleles: processing threads timed out.");
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalStateException("BitAlignment: loadSBitAlleles: processing threads problem.");
        }
        
    }
    
    private class ProcessSite implements Runnable {
        
        private BitSet[][] myData;
        private byte[][] myOrigData;
        private int mySite;
        
        public ProcessSite(byte[][] origData, BitSet[][] data, int site) {
            myData = data;
            myOrigData = origData;
            mySite = site;
        }
        
        public void run() {
            int numSeqs = getSequenceCount();
            byte[] cb = new byte[2];
            for (int t = 0; t < numSeqs; t++) {
                cb[0] = (byte) ((myOrigData[t][mySite] >>> 4) & 0xf);
                cb[1] = (byte) (myOrigData[t][mySite] & 0xf);
                for (int i = 0; i < 2; i++) {
                    if (cb[i] != Alignment.UNKNOWN_ALLELE) {
                        boolean isRare = true;
                        for (int j = 0; j < myMaxNumAlleles; j++) {
                            if (cb[i] == myAlleles[mySite][j]) {
                                myData[j][mySite].fastSet(t);
                                isRare = false;
                                break;
                            }
                        }
                        if (isRare && retainsRareAlleles()) {
                            myData[myMaxNumAlleles][mySite].fastSet(t);
                        }
                    }
                }
            }
        }
    }
    
    private void loadSBitAlleles(Alignment a, ProgressListener listener) {
        
        if (mySBitData != null) {
            return;
        }
        
        myNumDataRows = myMaxNumAlleles;
        if (retainsRareAlleles()) {
            myNumDataRows++;
        }
        int numSeqs = getSequenceCount();
        BitSet[][] temp = new OpenBitSet[myNumDataRows][myNumSites];
        for (int al = 0; al < myNumDataRows; al++) {
            for (int s = 0; s < myNumSites; s++) {
                temp[al][s] = new OpenBitSet(numSeqs);
            }
        }
        
        
        ExecutorService pool = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
        for (int s = 0; s < myNumSites; s++) {
            pool.execute(new ProcessLoadBitAllelesSite(a, temp, s, listener));
        }
        
        try {
            pool.shutdown();
            if (!pool.awaitTermination(600, TimeUnit.SECONDS)) {
                throw new IllegalStateException("BitAlignment: loadSBitAlleles: processing threads timed out.");
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalStateException("BitAlignment: loadSBitAlleles: processing threads problem.");
        }
        
        mySBitData = temp;
        
    }
    
    private class ProcessLoadBitAllelesSite implements Runnable {
        
        private BitSet[][] myData;
        private int mySite;
        private Alignment mySourceAlignment;
        private ProgressListener myListener;
        
        public ProcessLoadBitAllelesSite(Alignment a, BitSet[][] data, int site, ProgressListener listener) {
            myData = data;
            mySite = site;
            mySourceAlignment = a;
            myListener = listener;
        }
        
        public void run() {
            int numSeqs = getSequenceCount();
            for (int t = 0; t < numSeqs; t++) {
                byte[] cb = mySourceAlignment.getBaseArray(t, mySite);
                for (int i = 0; i < 2; i++) {
                    if (cb[i] != Alignment.UNKNOWN_ALLELE) {
                        boolean isRare = true;
                        for (int j = 0; j < myMaxNumAlleles; j++) {
                            if (cb[i] == myAlleles[mySite][j]) {
                                myData[j][mySite].fastSet(t);
                                isRare = false;
                                break;
                            }
                        }
                        if (isRare && retainsRareAlleles()) {
                            myData[myMaxNumAlleles][mySite].fastSet(t);
                        }
                    }
                }
            }
            
            if (myListener != null) {
                myListener.progress((int) (((double) (mySite + 1) / (double) myNumSites) * 100.0), null);
            }
        }
    }
    
    private void loadTBitAlleles(byte[][] data) {
        
        myNumDataRows = myMaxNumAlleles;
        if (retainsRareAlleles()) {
            myNumDataRows++;
        }
        myTBitData = new OpenBitSet[myNumDataRows][getSequenceCount()];
        for (int al = 0; al < myNumDataRows; al++) {
            for (int t = 0; t < getSequenceCount(); t++) {
                myTBitData[al][t] = new OpenBitSet(myNumSites);
            }
        }
        byte[] cb = new byte[2];
        for (int s = 0; s < myNumSites; s++) {
            for (int t = 0, n = getSequenceCount(); t < n; t++) {
                cb[0] = (byte) ((data[t][s] >>> 4) & 0xf);
                cb[1] = (byte) (data[t][s] & 0xf);
                for (int i = 0; i < 2; i++) {
                    if (cb[i] != Alignment.UNKNOWN_ALLELE) {
                        boolean isRare = true;
                        for (int j = 0; j < myMaxNumAlleles; j++) {
                            if (cb[i] == myAlleles[s][j]) {
                                myTBitData[j][t].fastSet(s);
                                isRare = false;
                                break;
                            }
                        }
                        if (isRare && retainsRareAlleles()) {
                            myTBitData[myMaxNumAlleles][t].fastSet(s);
                        }
                    }
                }
            }
        }
        
    }
    
    private void loadTBitAlleles(Alignment a, ProgressListener listener) {
        
        if (myTBitData != null) {
            return;
        }
        
        myNumDataRows = myMaxNumAlleles;
        if (retainsRareAlleles()) {
            myNumDataRows++;
        }
        int numTaxa = getSequenceCount();
        BitSet[][] temp = new OpenBitSet[myNumDataRows][getSequenceCount()];
        for (int al = 0; al < myNumDataRows; al++) {
            for (int t = 0; t < numTaxa; t++) {
                temp[al][t] = new OpenBitSet(myNumSites);
            }
        }
        
        ExecutorService pool = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
        for (int t = 0; t < numTaxa; t++) {
            pool.execute(new ProcessLoadBitAllelesTaxon(a, temp, t, listener));
        }
        
        try {
            pool.shutdown();
            if (!pool.awaitTermination(600, TimeUnit.SECONDS)) {
                throw new IllegalStateException("BitAlignment: loadTBitAlleles: processing threads timed out.");
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalStateException("BitAlignment: loadTBitAlleles: processing threads problem.");
        }
        
        myTBitData = temp;
        
    }
    
    private class ProcessLoadBitAllelesTaxon implements Runnable {
        
        private BitSet[][] myData;
        private int myTaxon;
        private Alignment mySourceAlignment;
        private ProgressListener myListener;
        
        public ProcessLoadBitAllelesTaxon(Alignment a, BitSet[][] data, int taxon, ProgressListener listener) {
            myData = data;
            myTaxon = taxon;
            mySourceAlignment = a;
            myListener = listener;
        }
        
        public void run() {
            int numSites = getSiteCount();
            for (int s = 0; s < numSites; s++) {
                byte[] cb = mySourceAlignment.getBaseArray(myTaxon, s);
                for (int i = 0; i < 2; i++) {
                    if (cb[i] != Alignment.UNKNOWN_ALLELE) {
                        boolean isRare = true;
                        for (int j = 0; j < myMaxNumAlleles; j++) {
                            if (cb[i] == myAlleles[s][j]) {
                                myData[j][myTaxon].fastSet(s);
                                isRare = false;
                                break;
                            }
                        }
                        if (isRare && retainsRareAlleles()) {
                            myData[myMaxNumAlleles][myTaxon].fastSet(s);
                        }
                    }
                }
            }
            
            if (myListener != null) {
                myListener.progress((int) (((double) (myTaxon + 1) / (double) getSequenceCount()) * 100.0), null);
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
        if (mySBitData != null) {
            return getBaseArraySBit(taxon, site);
        } else {
            return getBaseArrayTBit(taxon, site);
        }
    }
    
    private byte[] getBaseArraySBit(int taxon, int site) {
        byte[] result = new byte[2];
        result[0] = Alignment.UNKNOWN_ALLELE;
        result[1] = Alignment.UNKNOWN_ALLELE;
        try {
            int count = 0;
            for (int i = 0; i < myMaxNumAlleles; i++) {
                if (mySBitData[i][site].fastGet(taxon)) {
                    if (count == 0) {
                        result[1] = myAlleles[site][i];
                    }
                    result[count++] = myAlleles[site][i];
                }
            }

            // Check For Rare Allele
            if (retainsRareAlleles() && mySBitData[myMaxNumAlleles][site].fastGet(taxon)) {
                if (count == 0) {
                    result[1] = Alignment.RARE_ALLELE;
                }
                result[count] = Alignment.RARE_ALLELE;
            }
            
        } catch (IndexOutOfBoundsException e) {
            throw new IllegalStateException("BitAlignment: getBaseArray: bit sets indicate more than two alleles for taxon: " + taxon + "   site: " + site);
        }
        return result;
    }
    
    private byte[] getBaseArrayTBit(int taxon, int site) {
        byte[] result = new byte[2];
        result[0] = Alignment.UNKNOWN_ALLELE;
        result[1] = Alignment.UNKNOWN_ALLELE;
        try {
            int count = 0;
            for (int i = 0; i < myMaxNumAlleles; i++) {
                if (myTBitData[i][taxon].fastGet(site)) {
                    if (count == 0) {
                        result[1] = myAlleles[site][i];
                    }
                    result[count++] = myAlleles[site][i];
                }
            }

            // Check For Rare Allele
            if (retainsRareAlleles() && myTBitData[myMaxNumAlleles][taxon].fastGet(site)) {
                if (count == 0) {
                    result[1] = Alignment.RARE_ALLELE;
                }
                result[count] = Alignment.RARE_ALLELE;
            }
            
        } catch (IndexOutOfBoundsException e) {
            e.printStackTrace();
            throw new IllegalStateException("TBitAlignment: getBaseArray: bit sets indicate more than two alleles for taxon: " + taxon + "   site: " + site);
        }
        return result;
    }
    
    @Override
    public BitSet getAllelePresenceForAllTaxa(int site, int alleleNumber) {
        if (mySBitData != null) {
            return UnmodifiableBitSet.getInstance(mySBitData[alleleNumber][site]);
        } else {
            throw new IllegalStateException("BitAlignment: getAllelePresenceForAllTaxa: This alignment hasn't been optimized for Site Operations.");
        }
    }
    
    @Override
    public BitSet getAllelePresenceForAllSites(int taxon, int alleleNumber) {
        if (myTBitData != null) {
            return UnmodifiableBitSet.getInstance(myTBitData[alleleNumber][taxon]);
        } else {
            throw new IllegalStateException("BitAlignment: getAllelePresenceForAllSites: This alignment hasn't been optimized for Taxa Operations.");
        }
    }
    
    @Override
    public long[] getAllelePresenceForSitesBlock(int taxon, int alleleNumber, int startBlock, int endBlock) {
        if (myTBitData != null) {
            long[] result = new long[endBlock - startBlock];
            System.arraycopy(myTBitData[alleleNumber][taxon].getBits(), startBlock, result, 0, endBlock - startBlock);
            return result;
        } else {
            throw new IllegalStateException("BitAlignment: getAllelePresenceForSitesBlock: This alignment hasn't been optimized for Taxa Operations.");
        }
    }
    
    @Override
    public int getTotalGametesNotMissing(int site) {
        
        if (mySBitData == null) {
            return super.getTotalGametesNotMissing(site);
        }
        
        OpenBitSet temp = new OpenBitSet(getSequenceCount());
        for (int i = 0; i < myNumDataRows; i++) {
            temp.or(mySBitData[i][site]);
        }
        return ((int) temp.cardinality()) * 2;
        
    }
    
    @Override
    public int getTotalNotMissing(int site) {
        
        if (mySBitData == null) {
            return super.getTotalNotMissing(site);
        }
        
        OpenBitSet temp = new OpenBitSet(getSequenceCount());
        for (int i = 0; i < myNumDataRows; i++) {
            temp.or(mySBitData[i][site]);
        }
        return (int) temp.cardinality();
        
    }
    
    @Override
    public int getTotalGametesNotMissingForTaxon(int taxon) {
        
        if (myTBitData == null) {
            return super.getTotalGametesNotMissingForTaxon(taxon);
        }
        
        OpenBitSet temp = new OpenBitSet(getSiteCount());
        for (int i = 0; i < myNumDataRows; i++) {
            temp.or(myTBitData[i][taxon]);
        }
        return ((int) temp.cardinality()) * 2;
        
    }
    
    @Override
    public int getTotalNotMissingForTaxon(int taxon) {
        
        if (myTBitData == null) {
            return super.getTotalNotMissingForTaxon(taxon);
        }
        
        OpenBitSet temp = new OpenBitSet(getSiteCount());
        for (int i = 0; i < myNumDataRows; i++) {
            temp.or(myTBitData[i][taxon]);
        }
        return (int) temp.cardinality();
        
    }
    
    @Override
    public int getMinorAlleleCount(int site) {
        
        if ((myMaxNumAlleles < 2) || (myAlleles[site][1] == Alignment.UNKNOWN_ALLELE)) {
            return 0;
        }
        
        if (mySBitData == null) {
            return super.getMinorAlleleCount(site);
        }
        
        OpenBitSet temp = new OpenBitSet(getSequenceCount());
        for (int i = 0; i < myNumDataRows; i++) {
            if (i != 1) {
                temp.or(mySBitData[i][site]);
            }
        }
        temp.flip(0, temp.size());
        temp.and(mySBitData[1][site]);
        
        return (int) temp.cardinality() + (int) mySBitData[1][site].cardinality();
        
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
        
        if (myAlleles[site][0] == Alignment.UNKNOWN_ALLELE) {
            return 0;
        }
        
        if (mySBitData == null) {
            return super.getMajorAlleleCount(site);
        }
        
        OpenBitSet temp = new OpenBitSet(getSequenceCount());
        for (int i = 1; i < myNumDataRows; i++) {
            temp.or(mySBitData[i][site]);
        }
        temp.flip(0, temp.size());
        temp.and(mySBitData[0][site]);
        
        return (int) temp.cardinality() + (int) mySBitData[0][site].cardinality();
        
    }
    
    @Override
    public boolean isHeterozygous(int taxon, int site) {
        if (mySBitData != null) {
            return isHeterozygousSBit(taxon, site);
        } else {
            return isHeterozygousTBit(taxon, site);
        }
    }
    
    private boolean isHeterozygousSBit(int taxon, int site) {
        int count = 0;
        for (int i = 0; i < myNumDataRows; i++) {
            if (mySBitData[i][site].fastGet(taxon)) {
                count++;
                if (count == 2) {
                    return true;
                }
            }
        }
        return false;
    }
    
    private boolean isHeterozygousTBit(int taxon, int site) {
        int count = 0;
        for (int i = 0; i < myNumDataRows; i++) {
            if (myTBitData[i][taxon].fastGet(site)) {
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
        
        if (mySBitData == null) {
            return super.getHeterozygousCount(site);
        }
        
        int result = 0;
        for (int i = 0; i < myNumDataRows; i++) {
            for (int j = i + 1; j < myNumDataRows; j++) {
                result += (int) OpenBitSet.intersectionCount(mySBitData[i][site], mySBitData[j][site]);
            }
        }
        return result;
        
    }
    
    @Override
    public int getHeterozygousCountForTaxon(int taxon) {
        
        if (myTBitData == null) {
            return super.getHeterozygousCountForTaxon(taxon);
        }
        
        int result = 0;
        for (int i = 0; i < myNumDataRows; i++) {
            for (int j = i + 1; j < myNumDataRows; j++) {
                result += (int) OpenBitSet.intersectionCount(myTBitData[i][taxon], myTBitData[j][taxon]);
            }
        }
        return result;
        
    }
    
    @Override
    public boolean isPolymorphic(int site) {
        
        if (mySBitData == null) {
            return super.isPolymorphic(site);
        }
        
        boolean nonZero = false;
        for (int i = 0; i < myNumDataRows; i++) {
            int numTaxa = (int) mySBitData[i][site].cardinality();
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
        
        if ((myAlleleStates.length != 1) || (mySBitData == null)) {
            return super.getDiploidCounts();
        }
        
        long[][] counts = new long[16][16];
        for (int site = 0; site < myNumSites; site++) {
            for (int i = 0; i < myMaxNumAlleles; i++) {
                byte indexI = myAlleles[site][i];
                counts[indexI][indexI] += mySBitData[i][site].cardinality();
                for (int j = i + 1; j < myMaxNumAlleles; j++) {
                    byte indexJ = myAlleles[site][j];
                    long ijHet = OpenBitSet.intersectionCount(mySBitData[i][site], mySBitData[j][site]);
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
        
        if ((myAlleleStates.length != 1) || (mySBitData == null)) {
            return super.getDiploidssSortedByFrequency(site);
        }
        
        int[][] counts = new int[16][16];
        for (int i = 0; i < myMaxNumAlleles; i++) {
            byte indexI = myAlleles[site][i];
            counts[indexI][indexI] += (int) mySBitData[i][site].cardinality();
            for (int j = i + 1; j < myMaxNumAlleles; j++) {
                byte indexJ = myAlleles[site][j];
                int ijHet = (int) OpenBitSet.intersectionCount(mySBitData[i][site], mySBitData[j][site]);
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
        
        if (mySBitData == null) {
            return super.getAllelesSortedByFrequency(site);
        }
        
        int[] counts = new int[16];
        for (int i = 0; i < myNumDataRows; i++) {
            byte indexI;
            if ((retainsRareAlleles()) && (i == myMaxNumAlleles)) {
                indexI = Alignment.RARE_ALLELE;
            } else {
                indexI = myAlleles[site][i];
            }
            counts[indexI] += (int) mySBitData[i][site].cardinality() * 2;
            for (int j = i + 1; j < myNumDataRows; j++) {
                byte indexJ;
                if ((retainsRareAlleles()) && (j == myMaxNumAlleles)) {
                    indexJ = Alignment.RARE_ALLELE;
                } else {
                    indexJ = myAlleles[site][j];
                }
                int ijHet = (int) OpenBitSet.intersectionCount(mySBitData[i][site], mySBitData[j][site]);
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
        if (mySBitData != null) {
            return true;
        } else {
            return false;
        }
    }
    
    @Override
    public boolean isTBitFriendly() {
        if (myTBitData != null) {
            return true;
        } else {
            return false;
        }
    }
    
    @Override
    public int getTotalNumAlleles() {
        return myNumDataRows;
    }
    
    @Override
    public void optimizeForTaxa(ProgressListener listener) {
        if (myTBitData != null) {
            myLogger.info("optimizeForTaxa: Already Optimized for Taxa.");
            return;
        }
        myTBitData = BitUtil.transpose(mySBitData, myNumDataRows, myNumSites, getSequenceCount(), listener);
    }
    
    @Override
    public void optimizeForSites(ProgressListener listener) {
        if (mySBitData != null) {
            myLogger.info("optimizeForSites: Already Optimized for Sites.");
            return;
        }
        mySBitData = BitUtil.transpose(myTBitData, myNumDataRows, getSequenceCount(), myNumSites, listener);
    }
    
    
}
