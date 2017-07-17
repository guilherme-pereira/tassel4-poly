/*
 *  AbstractGenotype
 */
package net.maizegenetics.pal.alignment.genotype;

import net.maizegenetics.pal.alignment.AlignmentNew;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import org.apache.log4j.Logger;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 *
 * @author Terry Casstevens
 */
public abstract class AbstractGenotype implements Genotype {

    private static final Logger myLogger = Logger.getLogger(AbstractGenotype.class);
    protected final int myTaxaCount;
    protected final int mySiteCount;
    private final String[][] myAlleleEncodings;
    private final boolean myIsPhased;

    AbstractGenotype(int numTaxa, int numSites, boolean phased, String[][] alleleEncodings) {
        myTaxaCount = numTaxa;
        mySiteCount = numSites;
        myIsPhased = phased;
        myAlleleEncodings = alleleEncodings;
    }

    @Override
    public byte[] getBaseArray(int taxon, int site) {
        return AlignmentUtils.getDiploidValues(getBase(taxon, site));
    }

    @Override
    public byte[] getBaseRange(int taxon, int startSite, int endSite) {
        byte[] result = new byte[endSite - startSite];
        for (int i = startSite; i < endSite; i++) {
            result[i - startSite] = getBase(taxon, i);
        }
        return result;
    }

    @Override
    public byte[] getBaseRow(int taxon) {
        byte[] result = new byte[mySiteCount];
        for (int i = 0; i < mySiteCount; i++) {
            result[i] = getBase(taxon, i);
        }
        return result;
    }

    @Override
    public String getBaseAsString(int taxon, int site) {
        String[][] alleleStates = getAlleleEncodings();
        byte[] temp = getBaseArray(taxon, site);
        return alleleStates[0][temp[0]] + ":" + alleleStates[0][temp[1]];
    }

    @Override
    public String getBaseAsStringRange(int taxon, int startSite, int endSite) {
        StringBuilder builder = new StringBuilder();
        for (int i = startSite; i < endSite; i++) {
            if (i != startSite) {
                builder.append(";");
            }
            builder.append(getBaseAsString(taxon, i));
        }
        return builder.toString();
    }

    @Override
    public String getBaseAsStringRow(int taxon) {
        return getBaseAsStringRange(taxon, 0, mySiteCount);
    }

    @Override
    public String[] getBaseAsStringArray(int taxon, int site) {
        String[][] alleleStates = getAlleleEncodings();
        byte[] temp = getBaseArray(taxon, site);
        return new String[]{alleleStates[0][temp[0]], alleleStates[0][temp[1]]};
    }

    @Override
    public boolean isHeterozygous(int taxon, int site) {
        byte[] values = getBaseArray(taxon, site);
        if (values[0] == values[1]) {
            return false;
        } else {
            return true;
        }
    }

    @Override
    public int getHeterozygousCount(int site) {
        int result = 0;
        for (int i = 0, n = myTaxaCount; i < n; i++) {
            if (isHeterozygous(i, site)) {
                result++;
            }
        }
        return result;
    }

    @Override
    public boolean isPolymorphic(int site) {

        byte first = AlignmentNew.UNKNOWN_ALLELE;
        for (int i = 0, n = myTaxaCount; i < n; i++) {
            byte[] current = getBaseArray(i, site);
            if (current[0] != AlignmentNew.UNKNOWN_ALLELE) {
                if (first == AlignmentNew.UNKNOWN_ALLELE) {
                    first = current[0];
                } else if (first != current[0]) {
                    return true;
                }
            }
            if (current[1] != AlignmentNew.UNKNOWN_ALLELE) {
                if (first == AlignmentNew.UNKNOWN_ALLELE) {
                    first = current[1];
                } else if (first != current[1]) {
                    return true;
                }
            }
        }

        return false;

    }

    @Override
    public boolean isAllPolymorphic() {

        for (int i = 0, n = mySiteCount; i < n; i++) {
            if (!isPolymorphic(i)) {
                return false;
            }
        }

        return true;

    }

    @Override
    public boolean isPhased() {
        return myIsPhased;
    }

    @Override
    public boolean retainsRareAlleles() {
        return true;
    }

    @Override
    public String[][] getAlleleEncodings() {
        return myAlleleEncodings;
    }

    @Override
    public String[] getAlleleEncodings(int site) {
        if (myAlleleEncodings.length == 1) {
            return myAlleleEncodings[0];
        } else {
            return myAlleleEncodings[site];
        }
    }

    @Override
    public String getBaseAsString(int site, byte value) {
        return getAlleleEncodings(site)[value];
    }

    @Override
    public String getDiploidAsString(int site, byte value) {
        String[] alleleStates = getAlleleEncodings(site);
        return alleleStates[(value >>> 4) & 0xf] + ":" + alleleStates[value & 0xf];
    }

    @Override
    public int getMaxNumAlleles() {
        return NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES;
    }

    @Override
    public int getTotalNumAlleles() {
        int result = getMaxNumAlleles();
        if (retainsRareAlleles()) {
            result++;
        }
        if (isPhased()) {
            result++;
        }
        return result;
    }

    @Override
    public int getTotalGametesNotMissing(int site) {

        int result = 0;
        for (int i = 0, n = myTaxaCount; i < n; i++) {
            byte[] current = getBaseArray(i, site);
            if (current[0] != AlignmentNew.UNKNOWN_ALLELE) {
                result++;
            }
            if (current[1] != AlignmentNew.UNKNOWN_ALLELE) {
                result++;
            }
        }
        return result;

    }

    @Override
    public int getTotalNotMissing(int site) {

        int result = 0;
        for (int i = 0, n = myTaxaCount; i < n; i++) {
            byte current = getBase(i, site);
            if (current != AlignmentNew.UNKNOWN_DIPLOID_ALLELE) {
                result++;
            }
        }
        return result;

    }

    public byte[] getMajorAlleleForAllSites() {
        byte[] result = new byte[mySiteCount];
        for (int i = 0; i < mySiteCount; i++) {
            result[i] = getMajorAllele(i);
        }
        return result;
    }

    public byte[] getMinorAlleleForAllSites() {
        byte[] result = new byte[mySiteCount];
        for (int i = 0; i < mySiteCount; i++) {
            result[i] = getMinorAllele(i);
        }
        return result;
    }

    @Override
    public int getMinorAlleleCount(int site) {

        int[][] alleles = getAllelesSortedByFrequency(site);

        if (alleles[0].length >= 2) {
            return alleles[1][1];
        } else {
            return 0;
        }

    }

    @Override
    public byte getMinorAllele(int site) {
        int[][] alleles = getAllelesSortedByFrequency(site);

        if (alleles[0].length >= 2) {
            return (byte) alleles[0][1];
        } else {
            return AlignmentNew.UNKNOWN_ALLELE;
        }
    }

    @Override
    public String getMinorAlleleAsString(int site) {
        return getBaseAsString(site, getMinorAllele(site));
    }

    @Override
    public byte[] getMinorAlleles(int site) {
        int[][] alleles = getAllelesSortedByFrequency(site);
        int resultSize = alleles[0].length - 1;
        byte[] result = new byte[resultSize];
        for (int i = 0; i < resultSize; i++) {
            result[i] = (byte) alleles[0][i + 1];
        }
        return result;
    }

    @Override
    public int getMajorAlleleCount(int site) {

        int[][] alleles = getAllelesSortedByFrequency(site);

        if (alleles[0].length >= 1) {
            return alleles[1][0];
        } else {
            return 0;
        }

    }

    @Override
    public byte getMajorAllele(int site) {
        int[][] alleles = getAllelesSortedByFrequency(site);

        if (alleles[0].length >= 1) {
            return (byte) alleles[0][0];
        } else {
            return AlignmentNew.UNKNOWN_ALLELE;
        }
    }

    @Override
    public String getMajorAlleleAsString(int site) {
        return getBaseAsString(site, getMajorAllele(site));
    }

    @Override
    public double getMajorAlleleFrequency(int site) {

        int[][] alleles = getAllelesSortedByFrequency(site);

        int numAlleles = alleles[0].length;
        if (numAlleles >= 1) {
            int totalNonMissing = 0;
            for (int i = 0; i < numAlleles; i++) {
                totalNonMissing = totalNonMissing + alleles[1][i];
            }
            return (double) alleles[1][0] / (double) totalNonMissing;
        } else {
            return 0.0;
        }

    }

    @Override
    public double getMinorAlleleFrequency(int site) {

        int[][] alleles = getAllelesSortedByFrequency(site);

        int numAlleles = alleles[0].length;
        if (numAlleles >= 2) {
            int totalNonMissing = 0;
            for (int i = 0; i < numAlleles; i++) {
                totalNonMissing = totalNonMissing + alleles[1][i];
            }
            return (double) alleles[1][1] / (double) totalNonMissing;
        } else {
            return 0.0;
        }

    }

    @Override
    public int[][] getAllelesSortedByFrequency(int site) {

        int[] stateCnt = new int[16];
        for (int i = 0; i < myTaxaCount; i++) {
            byte[] dipB = getBaseArray(i, site);
            if (dipB[0] != AlignmentNew.UNKNOWN_ALLELE) {
                stateCnt[dipB[0]]++;
            }
            if (dipB[1] != AlignmentNew.UNKNOWN_ALLELE) {
                stateCnt[dipB[1]]++;
            }
        }

        int count = 0;
        for (int j = 0; j < 16; j++) {
            if (stateCnt[j] != 0) {
                count++;
            }
        }

        int result[][] = new int[2][count];
        int index = 0;
        for (int k = 0; k < 16; k++) {
            if (stateCnt[k] != 0) {
                result[0][index] = k;
                result[1][index] = stateCnt[k];
                index++;
            }
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0; k < count - 1; k++) {

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
    public Object[][] getDiploidsSortedByFrequency(int site) {

        Integer ONE_INTEGER = 1;

        Map<String, Integer> diploidValueCounts = new HashMap<String, Integer>();
        for (int r = 0; r < myTaxaCount; r++) {
            String current = getBaseAsString(r, site);
            Integer num = diploidValueCounts.get(current);
            if (num == null) {
                diploidValueCounts.put(current, ONE_INTEGER);
            } else {
                diploidValueCounts.put(current, ++num);
            }
        }

        Object[][] result = new Object[2][diploidValueCounts.size()];

        int i = 0;
        Iterator itr = diploidValueCounts.keySet().iterator();
        while (itr.hasNext()) {
            String key = (String) itr.next();
            Integer count = (Integer) diploidValueCounts.get(key);
            result[0][i] = key;
            result[1][i++] = count;
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0, n = diploidValueCounts.size() - 1; k < n; k++) {

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
    public Object[][] getDiploidCounts() {

        Map<String, Long> diploidValueCounts = new HashMap<String, Long>();
        for (int c = 0; c < mySiteCount; c++) {
            Object[][] diploids = getDiploidsSortedByFrequency(c);
            for (int i = 0; i < diploids[0].length; i++) {
                String current = (String) diploids[0][i];
                Long count = (long) ((Integer) diploids[1][i]).intValue();
                Long num = diploidValueCounts.get(current);
                if (num == null) {
                    diploidValueCounts.put(current, count);
                } else {
                    diploidValueCounts.put(current, (num + count));
                }
            }
        }

        Object[][] result = new Object[2][diploidValueCounts.size()];

        int i = 0;
        Iterator itr = diploidValueCounts.keySet().iterator();
        while (itr.hasNext()) {
            String key = (String) itr.next();
            Long count = diploidValueCounts.get(key);
            result[0][i] = key;
            result[1][i++] = count;
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0, n = diploidValueCounts.size() - 1; k < n; k++) {

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
    public Object[][] getMajorMinorCounts() {

        String[][] alleleStates = getAlleleEncodings();

        if (alleleStates.length != 1) {
            return new Object[0][0];
        }

        long[][] counts = new long[16][16];

        if (getMaxNumAlleles() >= 2) {
            for (int site = 0; site < mySiteCount; site++) {
                byte[] alleles = getAlleles(site);
                byte indexI = alleles[0];
                byte indexJ = alleles[1];
                if (indexJ == AlignmentNew.UNKNOWN_ALLELE) {
                    indexJ = indexI;
                }
                counts[indexI][indexJ]++;
            }
        } else {
            for (int site = 0; site < mySiteCount; site++) {
                byte[] alleles = getAlleles(site);
                byte indexI = alleles[0];
                counts[indexI][indexI]++;
            }
        }

        int numAlleles = 0;
        for (byte x = 0; x < 16; x++) {
            for (byte y = 0; y < 16; y++) {
                if (counts[x][y] != 0) {
                    numAlleles++;
                }
            }
        }

        Object[][] result = new Object[2][numAlleles];
        int nextResult = 0;
        for (byte x = 0; x < 16; x++) {
            for (byte y = 0; y < 16; y++) {
                if (counts[x][y] != 0) {
                    result[0][nextResult] = getBaseAsString(0, x) + ":" + getBaseAsString(0, y);
                    result[1][nextResult++] = counts[x][y];
                }
            }
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
    public int getTotalGametesNotMissingForTaxon(int taxon) {

        int result = 0;
        for (int i = 0, n = mySiteCount; i < n; i++) {
            byte[] current = getBaseArray(taxon, i);
            if (current[0] != AlignmentNew.UNKNOWN_ALLELE) {
                result++;
            }
            if (current[1] != AlignmentNew.UNKNOWN_ALLELE) {
                result++;
            }
        }
        return result;

    }

    @Override
    public int getHeterozygousCountForTaxon(int taxon) {
        int result = 0;
        for (int i = 0, n = mySiteCount; i < n; i++) {
            if (isHeterozygous(taxon, i)) {
                result++;
            }
        }
        return result;
    }

    @Override
    public int getTotalNotMissingForTaxon(int taxon) {

        int result = 0;
        for (int i = 0, n = mySiteCount; i < n; i++) {
            byte current = getBase(taxon, i);
            if (current != AlignmentNew.UNKNOWN_DIPLOID_ALLELE) {
                result++;
            }
        }
        return result;

    }

    @Override
    public byte[] getAlleles(int site) {
        int[][] alleles = getAllelesSortedByFrequency(site);
        int resultSize = alleles[0].length;
        int maxNumAlleles = getMaxNumAlleles();
        byte[] result = new byte[maxNumAlleles];
        for (int i = 0; i < maxNumAlleles; i++) {
            result[i] = (i < resultSize) ? (byte) alleles[0][i] : AlignmentNew.UNKNOWN_ALLELE;
        }
        return result;
    }

    @Override
    public int getSiteCount() {
        return mySiteCount;
    }

    @Override
    public int getTaxaCount() {
        return myTaxaCount;
    }

    @Override
    public byte[] getGenotypeForAllSites(int taxon) {
        int numSites = getSiteCount();
        byte[] result = new byte[numSites];
        for (int i = 0; i < numSites; i++) {
            result[i] = getBase(taxon, i);
        }
        return result;
    }

    @Override
    public byte[] getGenotypeForSiteRange(int taxon, int start, int end) {
        int numSites = end - start;
        byte[] result = new byte[numSites];
        for (int i = start; i < end; i++) {
            result[i] = getBase(taxon, i);
        }
        return result;
    }

    @Override
    public byte[] getGenotypeForAllTaxa(int site) {
        int numTaxa = getTaxaCount();
        byte[] result = new byte[numTaxa];
        for (int i = 0; i < numTaxa; i++) {
            result[i] = getBase(i, site);
        }
        return result;
    }
}
