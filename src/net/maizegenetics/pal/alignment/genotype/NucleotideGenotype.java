/*
 *  NucleotideGenotype
 */
package net.maizegenetics.pal.alignment.genotype;

import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Map;

import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.util.SuperByteMatrix;

/**
 *
 * @author Terry Casstevens
 */
public class NucleotideGenotype extends ByteGenotype {

    private static final int SHIFT_AMOUNT = 10;
    private static final int NUM_SITES_TO_CACHE = 1 << SHIFT_AMOUNT;
    public static final int SITE_BLOCK_MASK = ~(NUM_SITES_TO_CACHE - 1);
    private static final int MAX_CACHE_SIZE = 150;
    private final Map<Integer, int[][][]> myCachedInternal = new LinkedHashMap<Integer, int[][][]>((3 * MAX_CACHE_SIZE) / 2) {
        @Override
        protected boolean removeEldestEntry(Map.Entry eldest) {
            return size() > MAX_CACHE_SIZE;
        }
    };
    private final Map<Integer, int[][][]> myCachedAlleleFreqs = Collections.synchronizedMap(myCachedInternal);
    private final int myMaxNumThreads = Runtime.getRuntime().availableProcessors();
    private int myNumRunningThreads = 0;

    private static int getStartSite(int site) {
        return site & SITE_BLOCK_MASK;
    }

    NucleotideGenotype(SuperByteMatrix genotype, boolean phased) {
        super(genotype, phased, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES);
    }

    private int[][] getCachedAlleleFreq(int site) {
        int startSite = getStartSite(site);
        int[][][] result = myCachedAlleleFreqs.get(startSite);
        if (result == null) {
            if (myNumRunningThreads < myMaxNumThreads) {
                new Thread(new LookAheadSiteStats(startSite + NUM_SITES_TO_CACHE * 10)).start();
            }
            result = calculateAlleleFreq(startSite);
        }

        if ((site == startSite) && (myNumRunningThreads < myMaxNumThreads)) {
            new Thread(new LookAheadSiteStats(startSite + NUM_SITES_TO_CACHE * 20)).start();
        }
        return result[site - startSite];
    }

    @Override
    public int[][] getAllelesSortedByFrequency(int site) {
        return getCachedAlleleFreq(site);
    }

    private int[][][] calculateAlleleFreq(int site) {
        int startSite = getStartSite(site);
        int numSites = Math.min(NUM_SITES_TO_CACHE, getSiteCount() - startSite);
        int numTaxa = getTaxaCount();
        int[][] alleleFreq = new int[numSites][6];
        for (int taxon = 0; taxon < numTaxa; taxon++) {
            for (int s = 0; s < numSites; s++) {
                byte[] b = getBaseArray(taxon, s + startSite);
                if (b[0] < 6) {
                    alleleFreq[s][b[0]]++;
                }
                if (b[1] < 6) {
                    alleleFreq[s][b[1]]++;
                }
            }
        }

        int[][][] alleleCounts = new int[numSites][][];
        for (int s = 0; s < numSites; s++) {
            for (byte i = 0; i < 6; i++) {
                // size | allele (the 5-i is to get the sort right, so if case of ties A is first)
                alleleFreq[s][i] = (alleleFreq[s][i] << 4) | (5 - i);
            }

            int numAlleles = sort6(alleleFreq[s]);
            alleleCounts[s] = new int[2][numAlleles];
            for (int i = 0; i < numAlleles; i++) {
                alleleCounts[s][0][i] = (byte) (5 - (0xF & alleleFreq[s][i]));
                alleleCounts[s][1][i] = alleleFreq[s][i] >>> 4;
            }
        }
        myCachedAlleleFreqs.put(startSite, alleleCounts);

        return alleleCounts;
    }

    private static int sort6(int[] data) {
        int countNotZero = 0;
        for (int j = 0; j < 5; j++) {
            int imax = j;
            for (int i = j + 1; i < 6; i++) {
                if (data[i] > data[imax]) {
                    imax = i;
                }
            }
            if (data[imax] > 0xF) {
                int temp = data[j];
                data[j] = data[imax];
                data[imax] = temp;
                countNotZero++;
            } else {
                return countNotZero;
            }
        }
        if (data[5] > 0xF) {
            countNotZero++;
        }
        return countNotZero;
    }

    private class LookAheadSiteStats implements Runnable {

        private final int myStartSite;

        public LookAheadSiteStats(int site) {
            myNumRunningThreads++;
            myStartSite = getStartSite(site);
        }

        @Override
        public void run() {
            try {
                if (myStartSite >= mySiteCount) {
                    return;
                }
                if (myCachedAlleleFreqs.get(myStartSite) == null) {
                    calculateAlleleFreq(myStartSite);
                }
            } finally {
                myNumRunningThreads--;
            }
        }
    }
}
