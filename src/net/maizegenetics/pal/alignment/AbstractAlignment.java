/*
 * AbstractAlignment
 */
package net.maizegenetics.pal.alignment;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.ProgressListener;

/**
 *
 * @author terry, Gabriel Rodrigues Alves Margarido
 */
abstract public class AbstractAlignment implements Alignment {

    private static final long serialVersionUID = -5197800047652332969L;
    private IdGroup myIdGroup;
    private GeneticMap myGeneticMap;
    protected byte[] myReference;
    protected int myNumSites;
    protected String[][] myAlleleStates;
    private int[] myVariableSites;
    protected int myMaxNumAlleles = Alignment.DEFAULT_MAX_NUM_ALLELES;
    protected byte[][] myAlleles;
    protected Locus[] myLoci;
    /**
     * Loci offsets hold first site of each locus.
     */
    protected int[] myLociOffsets;
    private String[] mySNPIDs;
    private boolean myRetainRareAlleles = false;
    private Alignment myOriginalAlignment = null;

    public AbstractAlignment(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {
        if (idGroup.getIdCount() != data.length) {
            throw new IllegalArgumentException("AbstractAlignment: init: id group count: " + idGroup.getIdCount() + " doesn't equal number of data rows: " + data.length);
        }
        myNumSites = data[0].length;
        init(idGroup, map, reference, alleleStates, variableSites, maxNumAlleles, snpIDs, loci, lociOffsets, retainRareAlleles);
        initAlleles(data);
    }

    public AbstractAlignment(byte[][] alleles, IdGroup idGroup, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {
        myNumSites = alleles.length;
        init(idGroup, map, reference, alleleStates, variableSites, maxNumAlleles, snpIDs, loci, lociOffsets, retainRareAlleles);
        myAlleles = alleles;
    }

    public AbstractAlignment(Alignment a, int maxNumAlleles, boolean retainRareAlleles) {
        if (maxNumAlleles > a.getMaxNumAlleles()) {
            throw new IllegalArgumentException("AbstractAlignment: init: max number of alleles can't be larger than original alignment.");
        }
        myNumSites = a.getSiteCount();
        init(a.getIdGroup(), a.getGeneticMap(), a.getReference(), a.getAlleleEncodings(), a.getPhysicalPositions(), maxNumAlleles, a.getSNPIDs(), a.getLoci(), a.getLociOffsets(), retainRareAlleles);
        initAlleles(a);
        myOriginalAlignment = a;
    }

    /**
     * Constructor for FilterAlignment and CombineAlignment. Most attributes are
     * stored by wrapped alignments.
     *
     * @param idGroup id group
     */
    public AbstractAlignment(IdGroup idGroup) {
        myIdGroup = idGroup;
    }

    /**
     * Constructor for MutableAlignment.
     *
     * @param alleleStates
     */
    public AbstractAlignment(String[][] alleleStates) {
        init(null, null, null, alleleStates, null, 1, null, new Locus[]{new Locus("dummy", "0", 0, 0, null, null)}, new int[]{0}, false);
    }

    // TERRY - Need to check if this needed?
    public AbstractAlignment(IdGroup idGroup, String[][] alleleStates) {
        myIdGroup = idGroup;
        init(null, null, null, alleleStates, null, 1, null, new Locus[]{new Locus("dummy", "0", 0, 0, null, null)}, new int[]{0}, false);
    }

    private void init(IdGroup idGroup, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, String[] snpIDs, Locus[] loci, int[] lociOffsets, boolean retainRareAlleles) {

        if ((idGroup != null) && (idGroup.getIdCount() != 0)) {
            myIdGroup = SimpleIdGroup.getInstance(idGroup);
        }

        myGeneticMap = map;

        if ((reference == null) || (reference.length == 0)) {
            myReference = null;
        } else {
            myReference = reference;
            if (getSiteCount() != myReference.length) {
                throw new IllegalArgumentException("AbstractAlignment: init: reference differs in length with sequences.");
            }
        }

        if ((variableSites == null) || (variableSites.length == 0)) {
            myVariableSites = null;
        } else {
            myVariableSites = variableSites;
            if (getSiteCount() != myVariableSites.length) {
                throw new IllegalArgumentException("AbstractAlignment: init: variable sites differs in length with sequences.");
            }
        }

        if ((alleleStates == null) || (alleleStates.length == 0)) {
            throw new IllegalArgumentException("AbstractAlignment: init: allele states can't be empty.");
        } else if ((alleleStates.length != 1) && (alleleStates.length != getSiteCount())) {
            throw new IllegalArgumentException("AbstractAlignment: init: number of allele states must be either 1 or the number of sites.");
        }
        myAlleleStates = alleleStates;

        if ((maxNumAlleles < 1) || (maxNumAlleles > 14)) {
            throw new IllegalArgumentException("AbstractAlignment: init: max number of alleles must be between 1 and 14 inclusive: " + maxNumAlleles);
        }
        myMaxNumAlleles = maxNumAlleles;

        if ((snpIDs == null) || (snpIDs.length == 0)) {
            mySNPIDs = null;
        } else {
            mySNPIDs = snpIDs;
        }

        if ((loci == null) || (loci.length == 0)) {
            throw new IllegalArgumentException("AbstractAlignment: init: must have at least one locus.");
        } else if (loci.length != lociOffsets.length) {
            throw new IllegalArgumentException("AbstractAlignment: init: number of loci must match number of loci offsets.");
        } else {
            myLoci = loci;
            myLociOffsets = lociOffsets;
            for (int i = 0; i < lociOffsets.length; i++) {
                String name = myLoci[i].getName();
                int end = (i < myLoci.length - 1) ? myLociOffsets[i + 1] - 1 : myNumSites - 1;
                myLoci[i] = new Locus(name, name, myLociOffsets[i], end, null, null);
            }
        }

        myRetainRareAlleles = retainRareAlleles;

    }

    private void initAlleles(Alignment a) {
        myAlleles = new byte[getSiteCount()][myMaxNumAlleles];
        for (int i = 0, n = getSiteCount(); i < n; i++) {
            byte[] alleles = a.getAlleles(i);
            for (int j = 0; j < myMaxNumAlleles; j++) {
                myAlleles[i][j] = (j < alleles.length) ? alleles[j] : Alignment.UNKNOWN_ALLELE;
            }
        }
    }

    private void initAlleles(byte[][] data) {
        myAlleles = new byte[getSiteCount()][myMaxNumAlleles];
        for (int i = 0, n = getSiteCount(); i < n; i++) {
            byte[] alleles = AlignmentUtils.getAlleles(data, i);
            for (int j = 0; j < myMaxNumAlleles; j++) {
                myAlleles[i][j] = (j < alleles.length) ? alleles[j] : Alignment.UNKNOWN_ALLELE;
            }
        }
    }

    @Override
    public byte[] getBaseArray(int taxon, int site) {
        byte[] result = new byte[2];
        byte combinedBase = getBase(taxon, site);
        result[0] = (byte) ((combinedBase >>> 4) & 0xf);
        result[1] = (byte) (combinedBase & 0xf);
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
        return getBaseAsStringRange(taxon, 0, getSiteCount());
    }

    @Override
    public String[] getBaseAsStringArray(int taxon, int site) {
        String[][] alleleStates = getAlleleEncodings();
        byte[] temp = getBaseArray(taxon, site);
        return new String[]{alleleStates[0][temp[0]], alleleStates[0][temp[1]]};
    }

    @Override
    public byte[] getBaseRow(int taxon) {
        int numSites = getSiteCount();
        byte[] result = new byte[numSites];
        for (int i = 0; i < numSites; i++) {
            result[i] = getBase(taxon, i);
        }
        return result;
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
    public byte getBase(int taxon, Locus locus, int physicalPosition) {
        return getBase(taxon, getSiteOfPhysicalPosition(physicalPosition, locus));
    }

    @Override
    public BitSet getAllelePresenceForAllSites(int taxon, int alleleNumber) {
        throw new UnsupportedOperationException();
    }

    @Override
    public BitSet getAllelePresenceForAllTaxa(int site, int alleleNumber) {
        throw new UnsupportedOperationException();
    }

    @Override
    public long[] getAllelePresenceForSitesBlock(int taxon, int alleleNumber, int startBlock, int endBlock) {
        throw new UnsupportedOperationException();
    }

    @Override
    public BitSet getPhasedAllelePresenceForAllSites(int taxon, boolean firstParent, int alleleNumber) {
        throw new UnsupportedOperationException();
    }

    @Override
    public BitSet getPhasedAllelePresenceForAllTaxa(int site, boolean firstParent, int alleleNumber) {
        throw new UnsupportedOperationException();
    }

    @Override
    public long[] getPhasedAllelePresenceForSitesBlock(int taxon, boolean firstParent, int alleleNumber, int startBlock, int endBlock) {
        throw new UnsupportedOperationException();
    }

    @Override
    public IdGroup getIdGroup() {
        return myIdGroup;
    }

    @Override
    public String getTaxaName(int index) {
        return getIdGroup().getIdentifier(index).getName();
    }

    @Override
    public String getFullTaxaName(int index) {
        return getIdGroup().getIdentifier(index).getFullName();
    }

    @Override
    public int getSequenceCount() {
        return getIdGroup().getIdCount();
    }

    @Override
    public int getTaxaCount() {
        return getSequenceCount();
    }

    @Override
    public int[][] getAllelesSortedByFrequency(int site) {

        int[] stateCnt = new int[16];
        for (int i = 0; i < getSequenceCount(); i++) {
            byte[] dipB = getBaseArray(i, site);
            if (dipB[0] != Alignment.UNKNOWN_ALLELE) {
                stateCnt[dipB[0]]++;
            }
            if (dipB[1] != Alignment.UNKNOWN_ALLELE) {
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
    public byte getMajorAllele(int site) {
        if (myAlleles != null) {
            return myAlleles[site][0];
        } else {
            int[][] alleles = getAllelesSortedByFrequency(site);

            if (alleles[0].length >= 1) {
                return (byte) alleles[0][0];
            } else {
                return Alignment.UNKNOWN_ALLELE;
            }
        }
    }

    @Override
    public String getMajorAlleleAsString(int site) {
        return getBaseAsString(site, getMajorAllele(site));
    }

    @Override
    public byte getMinorAllele(int site) {
        if (myAlleles != null) {
            if (myMaxNumAlleles > 1) {
                return myAlleles[site][1];
            } else {
                return Alignment.UNKNOWN_ALLELE;
            }
        } else {
            int[][] alleles = getAllelesSortedByFrequency(site);

            if (alleles[0].length >= 2) {
                return (byte) alleles[0][1];
            } else {
                return Alignment.UNKNOWN_ALLELE;
            }
        }
    }

    @Override
    public String getMinorAlleleAsString(int site) {
        return getBaseAsString(site, getMinorAllele(site));
    }

    @Override
    public byte[] getAlleles(int site) {
        if (myAlleles != null) {
            return myAlleles[site];  //TODO this may not be correct for BitAlignment, as it returns UNKNOWN as one of the alleles
        } else {
            int[][] alleles = getAllelesSortedByFrequency(site);
            int resultSize = alleles[0].length;
            int maxNumAlleles = getMaxNumAlleles();
            byte[] result = new byte[maxNumAlleles];
            for (int i = 0; i < maxNumAlleles; i++) {
                result[i] = (i < resultSize) ? (byte) alleles[0][i] : Alignment.UNKNOWN_ALLELE;
            }
            return result;
        }
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
    public boolean isPolymorphic(int site) {

        byte first = Alignment.UNKNOWN_ALLELE;
        for (int i = 0, n = getSequenceCount(); i < n; i++) {
            byte[] current = getBaseArray(i, site);
            if (current[0] != Alignment.UNKNOWN_ALLELE) {
                if (first == Alignment.UNKNOWN_ALLELE) {
                    first = current[0];
                } else if (first != current[0]) {
                    return true;
                }
            }
            if (current[1] != Alignment.UNKNOWN_ALLELE) {
                if (first == Alignment.UNKNOWN_ALLELE) {
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

        for (int i = 0, n = getSiteCount(); i < n; i++) {
            if (!isPolymorphic(i)) {
                return false;
            }
        }

        return true;

    }

    @Override
    public int getIndelSize(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean isIndel(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public float getSiteScore(int seq, int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public float[][] getSiteScores() {
        if (hasSiteScores() == false) {
            return null;
        }
        float[][] f = new float[getSequenceCount()][getSiteCount()];
        for (int i = 0; i < getSequenceCount(); i++) {
            for (int j = 0; j < getSiteCount(); j++) {
                f[i][j] = getSiteScore(i, j);
            }
        }
        return f;
    }

    @Override
    public boolean hasSiteScores() {
        return false;
    }

    @Override
    public SITE_SCORE_TYPE getSiteScoreType() {
        return Alignment.SITE_SCORE_TYPE.None;
    }

    @Override
    public String getLocusName(int site) {
        return getLocus(site).getName();
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
        for (int i = 0, n = getSequenceCount(); i < n; i++) {
            if (isHeterozygous(i, site)) {
                result++;
            }
        }
        return result;
    }

    @Override
    public int getHeterozygousCountForTaxon(int taxon) {
        int result = 0;
        for (int i = 0, n = getSiteCount(); i < n; i++) {
            if (isHeterozygous(taxon, i)) {
                result++;
            }
        }
        return result;
    }

    @Override
    public boolean hasReference() {
        byte[] reference = getReference();
        if ((reference != null) && (reference.length != 0)) {
            return true;
        }
        return false;
    }

    @Override
    public byte getReferenceAllele(int site) {
        if (hasReference()) {
            return myReference[site];
        } else {
            return Alignment.UNKNOWN_DIPLOID_ALLELE;
        }
    }

    @Override
    public byte[] getReference() {
        return myReference;
    }

    @Override
    public byte[] getReference(int startSite, int endSite) {

        if (!hasReference()) {
            return null;
        }

        byte[] result = new byte[endSite - startSite];
        for (int i = startSite; i < endSite; i++) {
            result[i] = getReferenceAllele(i);
        }
        return result;

    }

    @Override
    public Alignment[] getAlignments() {
        return new Alignment[]{this};
    }

    @Override
    public boolean isPhased() {
        return false;
    }

    @Override
    public boolean isPositiveStrand(int site) {
        throw new UnsupportedOperationException("Not supported.");
    }

    @Override
    public String getGenomeAssembly() {
        throw new UnsupportedOperationException("Not supported.");
    }

    @Override
    public GeneticMap getGeneticMap() {
        return myGeneticMap;
    }

    @Override
    public int[] getPhysicalPositions() {
        if (myVariableSites == null) {
            return null;
        }
        return myVariableSites.clone();
    }

    @Override
    public int getSiteCount() {
        return myNumSites;
    }

    @Override
    public int getPositionInLocus(int site) {
        try {
            return myVariableSites[site];
        } catch (Exception e) {
            return site;
        }
    }

    @Override
    public Locus getLocus(int site) {
        for (int i = 1; i < myLociOffsets.length; i++) {
            if (myLociOffsets[i] > site) {
                return myLoci[i - 1];
            }
        }
        return myLoci[myLoci.length - 1];
    }

    @Override
    public Locus[] getLoci() {
        return myLoci;
    }

    @Override
    public Locus getLocus(String name) {

        name = name.trim();
        Locus[] temp = getLoci();
        for (int i = 0; i < temp.length; i++) {
            if (temp[i].getChromosomeName().equalsIgnoreCase(name)) {
                return temp[i];
            }
        }
        return null;

    }

    @Override
    public int getNumLoci() {
        return myLoci.length;
    }

    @Override
    public int[] getLociOffsets() {
        return myLociOffsets;
    }

    @Override
    public int getLocusSiteCount(Locus locus) {
        int[] startEnd = getStartAndEndOfLocus(locus);
        return startEnd[1] - startEnd[0];
    }

    @Override
    public int[] getStartAndEndOfLocus(Locus locus) {
        for (int i = 0; i < getNumLoci(); i++) {
            if (locus.equalName(myLoci[i])) {
                int end = 0;
                if (i == getNumLoci() - 1) {
                    end = getSiteCount();
                } else {
                    end = myLociOffsets[i + 1];
                }
                return new int[]{myLociOffsets[i], end};
            }
        }
        throw new IllegalArgumentException("AbstractAlignment: getStartAndEndOfLocus: this locus not defined: " + locus.getName());
    }

    @Override
    public String[] getSNPIDs() {
        return mySNPIDs;
    }

    @Override
    public String getSNPID(int site) {
        if (mySNPIDs == null) {
            return "S" + getLocus(site).getChromosomeName() + "_" + getPositionInLocus(site);
        } else {
            return mySNPIDs[site];
        }
    }

    @Override
    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus) {
        if (myVariableSites == null) {
            return physicalPosition;
        }
        if (locus == null) {
            locus = myLoci[0];
        }
        int[] startEnd = getStartAndEndOfLocus(locus);
        return Arrays.binarySearch(myVariableSites, startEnd[0], startEnd[1], physicalPosition);
    }

    @Override
    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus, String snpID) {
        if (myVariableSites == null) {
            return physicalPosition;
        }
        if (locus == null) {
            locus = myLoci[0];
        }
        int[] startEnd = getStartAndEndOfLocus(locus);
        int result = Arrays.binarySearch(myVariableSites, startEnd[0], startEnd[1], physicalPosition);
        if (result < 0) {
            return result;
        } else {
            if (snpID.equals(getSNPID(result))) {
                return result;
            } else {
                int index = result - 1;
                while ((index >= startEnd[0]) && (getPositionInLocus(index) == physicalPosition)) {
                    if (snpID.equals(getSNPID(index))) {
                        return index;
                    }
                    index--;
                }
                index = result + 1;
                while ((index < startEnd[1]) && (getPositionInLocus(index) == physicalPosition)) {
                    if (snpID.equals(getSNPID(index))) {
                        return index;
                    }
                    index++;
                }
                return -result - 1;
            }
        }
    }

    @Override
    public boolean retainsRareAlleles() {
        return myRetainRareAlleles;
    }

    @Override
    public String[][] getAlleleEncodings() {
        return myAlleleStates;
    }

    @Override
    public String[] getAlleleEncodings(int site) {
        if (myAlleleStates.length == 1) {
            return myAlleleStates[0];
        } else {
            return myAlleleStates[site];
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
        return myMaxNumAlleles;
    }

    @Override
    public int getTotalGametesNotMissing(int site) {

        int result = 0;
        for (int i = 0, n = getSequenceCount(); i < n; i++) {
            byte[] current = getBaseArray(i, site);
            if (current[0] != Alignment.UNKNOWN_ALLELE) {
                result++;
            }
            if (current[1] != Alignment.UNKNOWN_ALLELE) {
                result++;
            }
        }
        return result;

    }

    @Override
    public int getTotalNotMissing(int site) {

        int result = 0;
        for (int i = 0, n = getSequenceCount(); i < n; i++) {
            byte current = getBase(i, site);
            if (current != Alignment.UNKNOWN_DIPLOID_ALLELE) {
                result++;
            }
        }
        return result;

    }

    @Override
    public int getTotalGametesNotMissingForTaxon(int taxon) {

        int result = 0;
        for (int i = 0, n = getSiteCount(); i < n; i++) {
            byte[] current = getBaseArray(taxon, i);
            if (current[0] != Alignment.UNKNOWN_ALLELE) {
                result++;
            }
            if (current[1] != Alignment.UNKNOWN_ALLELE) {
                result++;
            }
        }
        return result;

    }

    @Override
    public int getTotalNotMissingForTaxon(int taxon) {

        int result = 0;
        for (int i = 0, n = getSiteCount(); i < n; i++) {
            byte current = getBase(taxon, i);
            if (current != Alignment.UNKNOWN_DIPLOID_ALLELE) {
                result++;
            }
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
    public int getMajorAlleleCount(int site) {

        int[][] alleles = getAllelesSortedByFrequency(site);

        if (alleles[0].length >= 1) {
            return alleles[1][0];
        } else {
            return 0;
        }

    }

    @Override
    public Object[][] getMajorMinorCounts() {

        String[][] alleleStates = getAlleleEncodings();

        if (alleleStates.length != 1) {
            return new Object[0][0];
        }

        int numSites = getSiteCount();
        long[][] counts = new long[16][16];

        if (getMaxNumAlleles() >= 2) {
            for (int site = 0; site < numSites; site++) {
                byte[] alleles = getAlleles(site);
                byte indexI = alleles[0];
                byte indexJ = alleles[1];
                if (indexJ == UNKNOWN_ALLELE) {
                    indexJ = indexI;
                }
                counts[indexI][indexJ]++;
            }
        } else {
            for (int site = 0; site < numSites; site++) {
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
    public Object[][] getDiploidCounts() {

        int numSites = getSiteCount();
        int numTaxa = getSequenceCount();

        Map<String, Long> diploidValueCounts = new HashMap<String, Long>();
        for (int c = 0; c < numSites; c++) {
            Object[][] diploids = getDiploidssSortedByFrequency(c);
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
    public Object[][] getDiploidssSortedByFrequency(int site) {

        Integer ONE_INTEGER = 1;
        int numTaxa = getSequenceCount();

        Map<String, Integer> diploidValueCounts = new HashMap<String, Integer>();
        for (int r = 0; r < numTaxa; r++) {
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
    public int getTotalNumAlleles() {
        throw new UnsupportedOperationException("Not supported.");
    }

    @Override
    public void optimizeForTaxa(ProgressListener listener) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void optimizeForSites(ProgressListener listener) {
        throw new UnsupportedOperationException();
    }

    @Override
    public short[] getDepthForAlleles(int taxon, int site) {
        return new short[]{1, 1};
    }

    @Override
    public byte[] getAllelesByScope(ALLELE_SCOPE_TYPE scope, int site) {
        if (scope == ALLELE_SCOPE_TYPE.Frequency) {
            return getAlleles(site);
        } else if (scope == ALLELE_SCOPE_TYPE.Reference) {
            return AlignmentUtils.getDiploidValues(getReferenceAllele(site));
        } else if (myOriginalAlignment != null) {
            return myOriginalAlignment.getAllelesByScope(scope, site);
        } else if (scope == ALLELE_SCOPE_TYPE.Global_Frequency) {
            return getAlleles(site);
        } else {
            throw new UnsupportedOperationException("AbstractAlignment: getAllelesByScope: This Alignment does not support scope: " + scope);
        }
    }

    @Override
    public BitSet getAllelePresenceForAllTaxaByScope(ALLELE_SCOPE_TYPE scope, int site, int alleleNumber) {
        if (scope == ALLELE_SCOPE_TYPE.Frequency) {
            return getAllelePresenceForAllTaxa(site, alleleNumber);
        } else if (scope == ALLELE_SCOPE_TYPE.Reference) {
            byte[] reference = AlignmentUtils.getDiploidValues(getReferenceAllele(site));
            int numTaxa = getSequenceCount();
            BitSet result = new OpenBitSet(numTaxa);
            for (int i = 0; i < numTaxa; i++) {
                byte[] current = getBaseArray(i, site);
                if ((current[0] == reference[alleleNumber]) || (current[1] == reference[alleleNumber])) {
                    result.fastSet(i);
                }
            }
            return result;
        } else if (myOriginalAlignment != null) {
            return myOriginalAlignment.getAllelePresenceForAllTaxaByScope(scope, site, alleleNumber);
        } else if (scope == ALLELE_SCOPE_TYPE.Global_Frequency) {
            return getAllelePresenceForAllTaxa(site, alleleNumber);
        } else {
            throw new UnsupportedOperationException("AbstractAlignment: getAllelePresenceForAllTaxaByScope: This Alignment does not support scope: " + scope);
        }
    }
}
