/*
 * BitPhasedAlignment
 */
package net.maizegenetics.pal.alignment;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.BitUtil;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.ProgressListener;
import net.maizegenetics.util.UnmodifiableBitSet;
import org.apache.log4j.Logger;

/**
 *
 * @author terry
 */
public class BitPhasedAlignment extends AbstractAlignment {

    private static final long serialVersionUID = -5197800047652332969L;
    private static final Logger myLogger = Logger.getLogger(BitPhasedAlignment.class);
    private BitSet[][] mySBitData0;
    private BitSet[][] mySBitData1;
    private BitSet[][] myTBitData0;
    private BitSet[][] myTBitData1;
    private int myNumDataRows;

    protected BitPhasedAlignment(Alignment a, int maxNumAlleles, boolean retainRareAlleles, boolean isSBit) {
        super(a, maxNumAlleles, retainRareAlleles);
        if (isSBit) {
            loadSBitAlleles(a, null);
        } else {
            loadTBitAlleles(a, null);
        }
    }

    protected BitPhasedAlignment(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isSBit) {
        super(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);
        if (isSBit) {
            loadSBitAlleles(data);
        } else {
            loadTBitAlleles(data);
        }
    }

    protected BitPhasedAlignment(IdGroup idGroup, byte[][] alleles, BitSet[][] data0, BitSet[][] data1, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isSBit) {
        super(alleles, idGroup, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);
        if (isSBit) {
            mySBitData0 = data0;
            mySBitData1 = data1;
        } else {
            myTBitData0 = data0;
            myTBitData1 = data1;
        }
        myNumDataRows = myMaxNumAlleles;
        if (retainsRareAlleles()) {
            myNumDataRows++;
        }
    }

    public static Alignment getInstance(Alignment a, boolean isSBit) {
        return BitPhasedAlignment.getInstance(a, a.getMaxNumAlleles(), a.retainsRareAlleles(), isSBit);
    }

    public static Alignment getInstance(Alignment a, int maxNumAlleles, boolean retainRareAlleles, boolean isSBit) {

        if ((a instanceof BitPhasedAlignment) && (a.getMaxNumAlleles() == maxNumAlleles) && (a.retainsRareAlleles() == retainRareAlleles)) {
            return (BitPhasedAlignment) a;
        }

        String[][] alleleStates = a.getAlleleEncodings();
        if ((alleleStates == null) || (alleleStates.length == 0)) {
            throw new IllegalStateException("BitPhasedAlignment: init: allele states should not be empty.");
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
            return new BitPhasedNucleotideAlignment(a, maxNumAlleles, retainRareAlleles, isSBit);
        } else if (alleleStates.length == 1) {
            return new BitPhasedAlignment(a, maxNumAlleles, retainRareAlleles, isSBit);
        } else {
            return new BitPhasedTextAlignment(a, maxNumAlleles, retainRareAlleles, isSBit);
        }

    }

    public static Alignment getInstance(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isSBit) {
        if ((alleleStates == null) || (alleleStates.length == 0)) {
            throw new IllegalArgumentException("BitPhasedAlignment: init: allele states can not be empty.");
        }
        if (alleleStates.length == 1) {
            return new BitPhasedAlignment(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isSBit);
        } else {
            return new BitPhasedTextAlignment(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isSBit);
        }
    }

    public static Alignment getNucleotideInstance(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isSBit) {
        return new BitPhasedNucleotideAlignment(idGroup, data, map, reference, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isSBit);
    }

    public static Alignment getNucleotideInstance(IdGroup idGroup, byte[][] alleles, BitSet[][] data0, BitSet[][] data1, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isSBit) {
        return new BitPhasedNucleotideAlignment(idGroup, alleles, data0, data1, map, reference, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isSBit);
    }

    public static Alignment getNucleotideInstance(IdGroup idGroup, String[][] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isSBit) {

        if ((maxNumAlleles < 1) || (maxNumAlleles > NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES)) {
            throw new IllegalArgumentException("BitPhasedAlignment: getNucleotideInstance: max number of alleles must be between 1 and " + NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES + " inclusive: " + maxNumAlleles);
        }

        if ((data == null) || (data.length == 0)) {
            throw new IllegalArgumentException("BitPhasedAlignment: getNucleotideInstance: data can not be empty.");
        }

        if (data.length != idGroup.getIdCount()) {
            throw new IllegalArgumentException("BitPhasedAlignment: getNucleotideInstance: data rows not equal to number of identifers.");
        }

        byte[][] dataBytes = AlignmentUtils.getDataBytes(data, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES, NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES);

        return BitPhasedAlignment.getNucleotideInstance(idGroup, dataBytes, map, reference, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isSBit);

    }

    public static Alignment getNucleotideInstance(IdGroup idGroup, String[] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isSBit) {

        if ((maxNumAlleles < 1) || (maxNumAlleles > NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES)) {
            throw new IllegalArgumentException("BitPhasedAlignment: getNucleotideInstance: max number of alleles must be between 1 and " + NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES + " inclusive: " + maxNumAlleles);
        }

        if ((data == null) || (data.length == 0)) {
            throw new IllegalArgumentException("BitPhasedAlignment: getNucleotideInstance: data can not be empty.");
        }

        if (data.length != idGroup.getIdCount()) {
            throw new IllegalArgumentException("BitPhasedAlignment: getNucleotideInstance: data rows not equal to number of identifers.");
        }

        byte[][] dataBytes = AlignmentUtils.getDataBytes(data);

        return BitPhasedAlignment.getNucleotideInstance(idGroup, dataBytes, map, reference, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isSBit);

    }

    public static Alignment getInstance(IdGroup idGroup, String[][] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isSBit) {

        if ((maxNumAlleles < 1) || (maxNumAlleles > 14)) {
            throw new IllegalArgumentException("BitPhasedAlignment: getInstance: max number of alleles must be between 1 and 14 inclusive: " + maxNumAlleles);
        }

        if ((data == null) || (data.length == 0)) {
            throw new IllegalArgumentException("BitPhasedAlignment: getInstance: data can not be empty.");
        }

        if (data.length != idGroup.getIdCount()) {
            throw new IllegalArgumentException("BitPhasedAlignment: getInstance: data rows not equal to number of identifers.");
        }

        String[][] alleleStates = AlignmentUtils.getAlleleStates(data, maxNumAlleles);

        byte[][] dataBytes = AlignmentUtils.getDataBytes(data, alleleStates, maxNumAlleles);

        return BitPhasedAlignment.getInstance(idGroup, dataBytes, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isSBit);

    }

    private void loadSBitAlleles(byte[][] data) {

        myNumDataRows = myMaxNumAlleles;
        if (retainsRareAlleles()) {
            myNumDataRows++;
        }
        int numSeqs = getSequenceCount();
        mySBitData0 = new OpenBitSet[myNumDataRows][myNumSites];
        mySBitData1 = new OpenBitSet[myNumDataRows][myNumSites];
        for (int al = 0; al < myNumDataRows; al++) {
            for (int s = 0; s < myNumSites; s++) {
                mySBitData0[al][s] = new OpenBitSet(numSeqs);
                mySBitData1[al][s] = new OpenBitSet(numSeqs);
            }
        }
        ExecutorService pool = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
        for (int s = 0; s < myNumSites; s++) {
            pool.execute(new BitPhasedAlignment.ProcessSite(data, mySBitData0, mySBitData1, s));
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

        private BitSet[][] myData0;
        private BitSet[][] myData1;
        private byte[][] myOrigData;
        private int mySite;

        public ProcessSite(byte[][] origData, BitSet[][] data0, BitSet[][] data1, int site) {
            myData0 = data0;
            myData1 = data1;
            myOrigData = origData;
            mySite = site;
        }

        public void run() {
            int numSeqs = getSequenceCount();
            byte[] cb = new byte[2];
            for (int t = 0; t < numSeqs; t++) {
                cb[0] = (byte) ((myOrigData[t][mySite] >>> 4) & 0xf);
                cb[1] = (byte) (myOrigData[t][mySite] & 0xf);
                if (cb[0] != Alignment.UNKNOWN_ALLELE) {
                    boolean isRare = true;
                    for (int j = 0; j < myMaxNumAlleles; j++) {
                        if (cb[0] == myAlleles[mySite][j]) {
                            myData0[j][mySite].fastSet(t);
                            isRare = false;
                            break;
                        }
                    }
                    if (isRare && retainsRareAlleles()) {
                        myData0[myMaxNumAlleles][mySite].fastSet(t);
                    }
                }
                if (cb[1] != Alignment.UNKNOWN_ALLELE) {
                    boolean isRare = true;
                    for (int j = 0; j < myMaxNumAlleles; j++) {
                        if (cb[1] == myAlleles[mySite][j]) {
                            myData1[j][mySite].fastSet(t);
                            isRare = false;
                            break;
                        }
                    }
                    if (isRare && retainsRareAlleles()) {
                        myData1[myMaxNumAlleles][mySite].fastSet(t);
                    }
                }
            }
        }
    }

    private void loadSBitAlleles(Alignment a, ProgressListener listener) {

        if (mySBitData0 != null) {
            return;
        }

        myNumDataRows = myMaxNumAlleles;
        if (retainsRareAlleles()) {
            myNumDataRows++;
        }
        int numSeqs = getSequenceCount();
        BitSet[][] temp0 = new OpenBitSet[myNumDataRows][myNumSites];
        BitSet[][] temp1 = new OpenBitSet[myNumDataRows][myNumSites];
        for (int al = 0; al < myNumDataRows; al++) {
            for (int s = 0; s < myNumSites; s++) {
                temp0[al][s] = new OpenBitSet(numSeqs);
                temp1[al][s] = new OpenBitSet(numSeqs);
            }
        }


        ExecutorService pool = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
        for (int s = 0; s < myNumSites; s++) {
            pool.execute(new BitPhasedAlignment.ProcessLoadBitAllelesSite(a, temp0, temp1, s, listener));
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

        mySBitData0 = temp0;
        mySBitData1 = temp1;

    }

    private class ProcessLoadBitAllelesSite implements Runnable {

        private BitSet[][] myData0;
        private BitSet[][] myData1;
        private int mySite;
        private Alignment mySourceAlignment;
        private ProgressListener myListener;

        public ProcessLoadBitAllelesSite(Alignment a, BitSet[][] data0, BitSet[][] data1, int site, ProgressListener listener) {
            myData0 = data0;
            myData1 = data1;
            mySite = site;
            mySourceAlignment = a;
            myListener = listener;
        }

        public void run() {
            int numSeqs = getSequenceCount();
            for (int t = 0; t < numSeqs; t++) {
                byte[] cb = mySourceAlignment.getBaseArray(t, mySite);
                if (cb[0] != Alignment.UNKNOWN_ALLELE) {
                    boolean isRare = true;
                    for (int j = 0; j < myMaxNumAlleles; j++) {
                        if (cb[0] == myAlleles[mySite][j]) {
                            myData0[j][mySite].fastSet(t);
                            isRare = false;
                            break;
                        }
                    }
                    if (isRare && retainsRareAlleles()) {
                        myData0[myMaxNumAlleles][mySite].fastSet(t);
                    }
                }
                if (cb[1] != Alignment.UNKNOWN_ALLELE) {
                    boolean isRare = true;
                    for (int j = 0; j < myMaxNumAlleles; j++) {
                        if (cb[1] == myAlleles[mySite][j]) {
                            myData1[j][mySite].fastSet(t);
                            isRare = false;
                            break;
                        }
                    }
                    if (isRare && retainsRareAlleles()) {
                        myData1[myMaxNumAlleles][mySite].fastSet(t);
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
        myTBitData0 = new OpenBitSet[myNumDataRows][getSequenceCount()];
        myTBitData1 = new OpenBitSet[myNumDataRows][getSequenceCount()];
        for (int al = 0; al < myNumDataRows; al++) {
            for (int t = 0; t < getSequenceCount(); t++) {
                myTBitData0[al][t] = new OpenBitSet(myNumSites);
                myTBitData1[al][t] = new OpenBitSet(myNumSites);
            }
        }
        byte[] cb = new byte[2];
        for (int s = 0; s < myNumSites; s++) {
            for (int t = 0, n = getSequenceCount(); t < n; t++) {
                cb[0] = (byte) ((data[t][s] >>> 4) & 0xf);
                cb[1] = (byte) (data[t][s] & 0xf);
                if (cb[0] != Alignment.UNKNOWN_ALLELE) {
                    boolean isRare = true;
                    for (int j = 0; j < myMaxNumAlleles; j++) {
                        if (cb[0] == myAlleles[s][j]) {
                            myTBitData0[j][t].fastSet(s);
                            isRare = false;
                            break;
                        }
                    }
                    if (isRare && retainsRareAlleles()) {
                        myTBitData0[myMaxNumAlleles][t].fastSet(s);
                    }
                }
                if (cb[1] != Alignment.UNKNOWN_ALLELE) {
                    boolean isRare = true;
                    for (int j = 0; j < myMaxNumAlleles; j++) {
                        if (cb[1] == myAlleles[s][j]) {
                            myTBitData1[j][t].fastSet(s);
                            isRare = false;
                            break;
                        }
                    }
                    if (isRare && retainsRareAlleles()) {
                        myTBitData1[myMaxNumAlleles][t].fastSet(s);
                    }
                }
            }
        }

    }

    private void loadTBitAlleles(Alignment a, ProgressListener listener) {

        if (myTBitData0 != null) {
            return;
        }

        myNumDataRows = myMaxNumAlleles;
        if (retainsRareAlleles()) {
            myNumDataRows++;
        }
        int numTaxa = getSequenceCount();
        BitSet[][] temp0 = new OpenBitSet[myNumDataRows][getSequenceCount()];
        BitSet[][] temp1 = new OpenBitSet[myNumDataRows][getSequenceCount()];
        for (int al = 0; al < myNumDataRows; al++) {
            for (int t = 0; t < numTaxa; t++) {
                temp0[al][t] = new OpenBitSet(myNumSites);
                temp1[al][t] = new OpenBitSet(myNumSites);
            }
        }

        ExecutorService pool = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
        for (int t = 0; t < numTaxa; t++) {
            pool.execute(new BitPhasedAlignment.ProcessLoadBitAllelesTaxon(a, temp0, temp1, t, listener));
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

        myTBitData0 = temp0;
        myTBitData1 = temp1;

    }

    private class ProcessLoadBitAllelesTaxon implements Runnable {

        private BitSet[][] myData0;
        private BitSet[][] myData1;
        private int myTaxon;
        private Alignment mySourceAlignment;
        private ProgressListener myListener;

        public ProcessLoadBitAllelesTaxon(Alignment a, BitSet[][] data0, BitSet[][] data1, int taxon, ProgressListener listener) {
            myData0 = data0;
            myData1 = data1;
            myTaxon = taxon;
            mySourceAlignment = a;
            myListener = listener;
        }

        public void run() {
            int numSites = getSiteCount();
            for (int s = 0; s < numSites; s++) {
                byte[] cb = mySourceAlignment.getBaseArray(myTaxon, s);
                if (cb[0] != Alignment.UNKNOWN_ALLELE) {
                    boolean isRare = true;
                    for (int j = 0; j < myMaxNumAlleles; j++) {
                        if (cb[0] == myAlleles[s][j]) {
                            myData0[j][myTaxon].fastSet(s);
                            isRare = false;
                            break;
                        }
                    }
                    if (isRare && retainsRareAlleles()) {
                        myData0[myMaxNumAlleles][myTaxon].fastSet(s);
                    }
                }
                if (cb[1] != Alignment.UNKNOWN_ALLELE) {
                    boolean isRare = true;
                    for (int j = 0; j < myMaxNumAlleles; j++) {
                        if (cb[1] == myAlleles[s][j]) {
                            myData1[j][myTaxon].fastSet(s);
                            isRare = false;
                            break;
                        }
                    }
                    if (isRare && retainsRareAlleles()) {
                        myData1[myMaxNumAlleles][myTaxon].fastSet(s);
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
        if (mySBitData0 != null) {
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
            for (int i = 0; i < myMaxNumAlleles; i++) {
                if (mySBitData0[i][site].fastGet(taxon)) {
                    result[0] = myAlleles[site][i];
                    break;
                }
            }

            for (int i = 0; i < myMaxNumAlleles; i++) {
                if (mySBitData1[i][site].fastGet(taxon)) {
                    result[1] = myAlleles[site][i];
                    break;
                }
            }

            // Check For Rare Allele
            if (retainsRareAlleles() && mySBitData0[myMaxNumAlleles][site].fastGet(taxon)) {
                result[0] = Alignment.RARE_ALLELE;
            }
            if (retainsRareAlleles() && mySBitData1[myMaxNumAlleles][site].fastGet(taxon)) {
                result[1] = Alignment.RARE_ALLELE;
            }

        } catch (IndexOutOfBoundsException e) {
            throw new IllegalStateException("BitPhasedAlignment: getBaseArraySBit: bit sets indicate more than two alleles for taxon: " + taxon + "   site: " + site);
        }
        return result;
    }

    private byte[] getBaseArrayTBit(int taxon, int site) {
        byte[] result = new byte[2];
        result[0] = Alignment.UNKNOWN_ALLELE;
        result[1] = Alignment.UNKNOWN_ALLELE;
        try {

            for (int i = 0; i < myMaxNumAlleles; i++) {
                if (myTBitData0[i][taxon].fastGet(site)) {
                    result[0] = myAlleles[site][i];
                    break;
                }
            }

            for (int i = 0; i < myMaxNumAlleles; i++) {
                if (myTBitData1[i][taxon].fastGet(site)) {
                    result[1] = myAlleles[site][i];
                    break;
                }
            }

            // Check For Rare Allele
            if (retainsRareAlleles() && myTBitData0[myMaxNumAlleles][taxon].fastGet(site)) {
                result[0] = Alignment.RARE_ALLELE;
            }
            if (retainsRareAlleles() && myTBitData1[myMaxNumAlleles][taxon].fastGet(site)) {
                result[1] = Alignment.RARE_ALLELE;
            }

        } catch (IndexOutOfBoundsException e) {
            e.printStackTrace();
            throw new IllegalStateException("BitPhasedAlignment: getBaseArrayTBit: bit sets indicate more than two alleles for taxon: " + taxon + "   site: " + site);
        }
        return result;
    }

    @Override
    public boolean isPhased() {
        return true;
    }

    @Override
    public boolean isSBitFriendly() {
        if (mySBitData0 != null) {
            return true;
        } else {
            return false;
        }
    }

    @Override
    public boolean isTBitFriendly() {
        if (myTBitData0 != null) {
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
    public BitSet getPhasedAllelePresenceForAllSites(int taxon, boolean firstParent, int alleleNumber) {
        if (myTBitData0 != null) {
            if (firstParent) {
                return UnmodifiableBitSet.getInstance(myTBitData0[alleleNumber][taxon]);
            } else {
                return UnmodifiableBitSet.getInstance(myTBitData1[alleleNumber][taxon]);
            }
        } else {
            throw new IllegalStateException("BitPhasedAlignment: getPhasedAllelePresenceForAllSites: This alignment hasn't been optimized for Taxa Operations.");
        }
    }

    @Override
    public BitSet getPhasedAllelePresenceForAllTaxa(int site, boolean firstParent, int alleleNumber) {
        if (mySBitData0 != null) {
            if (firstParent) {
                return UnmodifiableBitSet.getInstance(mySBitData0[alleleNumber][site]);
            } else {
                return UnmodifiableBitSet.getInstance(mySBitData1[alleleNumber][site]);
            }
        } else {
            throw new IllegalStateException("BitPhasedAlignment: getPhasedAllelePresenceForAllTaxa: This alignment hasn't been optimized for Site Operations.");
        }
    }

    @Override
    public long[] getPhasedAllelePresenceForSitesBlock(int taxon, boolean firstParent, int alleleNumber, int startBlock, int endBlock) {
        if (myTBitData0 != null) {
            BitSet temp = getAllelePresenceForAllSites(taxon, alleleNumber);
            long[] result = new long[endBlock - startBlock];
            System.arraycopy(temp.getBits(), startBlock, result, 0, endBlock - startBlock);
            return result;
        } else {
            throw new IllegalStateException("BitPhasedAlignment: getPhasedAllelePresenceForSitesBlock: This alignment hasn't been optimized for Taxa Operations.");
        }
    }

    @Override
    public BitSet getAllelePresenceForAllTaxa(int site, int alleleNumber) {
        if (mySBitData0 != null) {
            OpenBitSet temp = new OpenBitSet(getSequenceCount());
            temp.or(mySBitData0[alleleNumber][site]);
            temp.or(mySBitData1[alleleNumber][site]);
            return temp;
        } else {
            throw new IllegalStateException("BitPhasedAlignment: getAllelePresenceForAllTaxa: This alignment hasn't been optimized for Site Operations.");
        }
    }

    @Override
    public BitSet getAllelePresenceForAllSites(int taxon, int alleleNumber) {
        if (myTBitData0 != null) {
            OpenBitSet temp = new OpenBitSet(getSequenceCount());
            temp.or(myTBitData0[alleleNumber][taxon]);
            temp.or(myTBitData1[alleleNumber][taxon]);
            return temp;
        } else {
            throw new IllegalStateException("BitPhasedAlignment: getAllelePresenceForAllSites: This alignment hasn't been optimized for Taxa Operations.");
        }
    }

    @Override
    public long[] getAllelePresenceForSitesBlock(int taxon, int alleleNumber, int startBlock, int endBlock) {
        if (myTBitData0 != null) {
            BitSet temp = getAllelePresenceForAllSites(taxon, alleleNumber);
            long[] result = new long[endBlock - startBlock];
            System.arraycopy(temp.getBits(), startBlock, result, 0, endBlock - startBlock);
            return result;
        } else {
            throw new IllegalStateException("BitPhasedAlignment: getAllelePresenceForSitesBlock: This alignment hasn't been optimized for Taxa Operations.");
        }
    }

    @Override
    public int getTotalGametesNotMissing(int site) {

        if (mySBitData0 == null) {
            return super.getTotalGametesNotMissing(site);
        }

        OpenBitSet temp0 = new OpenBitSet(getSequenceCount());
        OpenBitSet temp1 = new OpenBitSet(getSequenceCount());
        for (int i = 0; i < myNumDataRows; i++) {
            temp0.or(mySBitData0[i][site]);
            temp1.or(mySBitData1[i][site]);
        }
        return (int) (temp0.cardinality() + temp1.cardinality());

    }

    @Override
    public int getTotalNotMissing(int site) {

        if (mySBitData0 == null) {
            return super.getTotalNotMissing(site);
        }

        OpenBitSet temp = new OpenBitSet(getSequenceCount());
        for (int i = 0; i < myNumDataRows; i++) {
            temp.or(mySBitData0[i][site]);
            temp.or(mySBitData1[i][site]);
        }
        return (int) temp.cardinality();

    }

    @Override
    public int getTotalGametesNotMissingForTaxon(int taxon) {

        if (myTBitData0 == null) {
            return super.getTotalGametesNotMissingForTaxon(taxon);
        }

        OpenBitSet temp0 = new OpenBitSet(getSiteCount());
        OpenBitSet temp1 = new OpenBitSet(getSiteCount());
        for (int i = 0; i < myNumDataRows; i++) {
            temp0.or(myTBitData0[i][taxon]);
            temp1.or(myTBitData1[i][taxon]);
        }
        return (int) (temp0.cardinality() + temp1.cardinality());

    }

    @Override
    public int getTotalNotMissingForTaxon(int taxon) {

        if (myTBitData0 == null) {
            return super.getTotalNotMissingForTaxon(taxon);
        }

        OpenBitSet temp = new OpenBitSet(getSiteCount());
        for (int i = 0; i < myNumDataRows; i++) {
            temp.or(myTBitData0[i][taxon]);
            temp.or(myTBitData1[i][taxon]);
        }
        return (int) temp.cardinality();

    }

    @Override
    public int getMinorAlleleCount(int site) {

        if ((myMaxNumAlleles < 2) || (myAlleles[site][1] == Alignment.UNKNOWN_ALLELE)) {
            return 0;
        }

        if (mySBitData0 == null) {
            return super.getMinorAlleleCount(site);
        }

        return (int) (mySBitData0[1][site].cardinality() + mySBitData1[1][site].cardinality());

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

        if (mySBitData0 == null) {
            return super.getMajorAlleleCount(site);
        }

        return (int) (mySBitData0[0][site].cardinality() + mySBitData1[0][site].cardinality());

    }

    @Override
    public byte getMajorAllele(int site) {
        return myAlleles[site][0];
    }

    @Override
    public byte getMinorAllele(int site) {

        if (myMaxNumAlleles > 1) {
            return myAlleles[site][1];
        } else {
            return Alignment.UNKNOWN_ALLELE;
        }

    }

    @Override
    public void optimizeForTaxa(ProgressListener listener) {
        if (myTBitData0 != null) {
            myLogger.info("optimizeForTaxa: Already Optimized for Taxa.");
            return;
        }
        myTBitData0 = BitUtil.transpose(mySBitData0, myNumDataRows, myNumSites, getSequenceCount(), listener);
        myTBitData1 = BitUtil.transpose(mySBitData1, myNumDataRows, myNumSites, getSequenceCount(), listener);
    }

    @Override
    public void optimizeForSites(ProgressListener listener) {
        if (mySBitData0 != null) {
            myLogger.info("optimizeForSites: Already Optimized for Sites.");
            return;
        }
        mySBitData0 = BitUtil.transpose(myTBitData0, myNumDataRows, getSequenceCount(), myNumSites, listener);
        mySBitData1 = BitUtil.transpose(myTBitData1, myNumDataRows, getSequenceCount(), myNumSites, listener);
    }
}
