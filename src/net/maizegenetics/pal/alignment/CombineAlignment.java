/*
 * CombineAlignment
 */
package net.maizegenetics.pal.alignment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.IdGroupUtils;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.ProgressListener;

/**
 *
 * @author terry
 */
public class CombineAlignment extends AbstractAlignment {

    private static final long serialVersionUID = -5197800047652332969L;
    private final Alignment[] myAlignments;
    private final int[] mySiteOffsets;
    private final Map myLoci = new HashMap();
    private Locus[] myLociList;
    private int[] myLociOffsets;

    private CombineAlignment(IdGroup subIdGroup, Alignment[] alignments) {
        super(subIdGroup);
        myAlignments = alignments;
        mySiteOffsets = new int[alignments.length + 1];

        mySiteOffsets[0] = 0;
        int count = 0;
        for (int i = 0; i < alignments.length; i++) {
            count = alignments[i].getSiteCount() + count;
            mySiteOffsets[i + 1] = count;

            Locus[] loci = alignments[i].getLoci();
            for (int j = 0; j < loci.length; j++) {
                myLoci.put(loci[j], alignments[i]);
            }
        }

        initLoci();
    }

    /**
     * This factory method combines given alignments. If only one alignment,
     * then it is returned unchanged. Otherwise, this requires that each
     * alignment has the same Identifiers in the same order.
     *
     * @param alignments
     * @return
     */
    public static Alignment getInstance(Alignment[] alignments) {

        if ((alignments == null) || (alignments.length == 0)) {
            throw new IllegalArgumentException("CombineAlignment: getInstance: must provide alignments.");
        }

        if (alignments.length == 1) {
            return alignments[0];
        }

        IdGroup firstGroup = alignments[0].getIdGroup();
        for (int i = 1; i < alignments.length; i++) {
            if (!areIdGroupsEqual(firstGroup, alignments[i].getIdGroup())) {
                throw new IllegalArgumentException("CombineAlignment: getInstance: IdGroups do not match.");
            }
        }

        return new CombineAlignment(firstGroup, alignments);

    }

    /**
     * This factory method combines given alignments. If only one alignment,
     * then it is returned unchanged. If isUnion equals true, a union join of
     * the Identifiers will be used to construct the combination. Any alignment
     * not containing one of the Identifiers will return unknown value for those
     * locations. If isUnion equals false, a intersect join of the Identifiers
     * will be used.
     *
     * @param alignments alignments to combine
     * @param isUnion whether to union or intersect join
     * @return
     */
    public static Alignment getInstance(Alignment[] alignments, boolean isUnion) {

        if ((alignments == null) || (alignments.length == 0)) {
            throw new IllegalArgumentException("CombineAlignment: getInstance: must provide alignments.");
        }

        if (alignments.length == 1) {
            return alignments[0];
        }

        IdGroup[] groups = new IdGroup[alignments.length];
        for (int i = 0; i < alignments.length; i++) {
            groups[i] = alignments[i].getIdGroup();
        }
        IdGroup newTaxa = null;
        if (isUnion) {
            newTaxa = IdGroupUtils.getAllIds(groups);
        } else {
            newTaxa = IdGroupUtils.getCommonIds(groups);
        }

        Alignment[] newAlignments = new Alignment[alignments.length];
        for (int i = 0; i < alignments.length; i++) {
            newAlignments[i] = FilterAlignment.getInstance(alignments[i], newTaxa);
        }

        return new CombineAlignment(newTaxa, newAlignments);

    }

    private static boolean areIdGroupsEqual(IdGroup first, IdGroup second) {

        if (first.getIdCount() != second.getIdCount()) {
            return false;
        }

        for (int i = 0, n = first.getIdCount(); i < n; i++) {
            if (!first.getIdentifier(i).equals(second.getIdentifier(i))) {
                return false;
            }
        }

        return true;

    }

    private void initLoci() {

        List offsets = new ArrayList();
        List<Locus> loci = new ArrayList();
        for (int i = 0; i < myAlignments.length; i++) {
            loci.addAll(Arrays.asList(myAlignments[i].getLoci()));
            int[] tempOffsets = myAlignments[i].getLociOffsets();
            for (int j = 0; j < tempOffsets.length; j++) {
                offsets.add(tempOffsets[j] + mySiteOffsets[i]);
            }
        }

        myLociList = new Locus[loci.size()];
        myLociList = loci.toArray(myLociList);

        myLociOffsets = new int[offsets.size()];
        for (int i = 0; i < offsets.size(); i++) {
            myLociOffsets[i] = (Integer) offsets.get(i);
        }

        if (myLociOffsets.length != myLociList.length) {
            throw new IllegalStateException("CombineAlignment: initLoci: number loci offsets should equal number of loci.");
        }

    }

    public byte getBase(int taxon, int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getBase(taxon, site - mySiteOffsets[translate]);
    }

    @Override
    public byte[] getBaseRange(int taxon, int startSite, int endSite) {

        int siteCount = getSiteCount();
        if ((startSite == 0) && (endSite == siteCount)) {
            byte[] result = new byte[siteCount];
            int count = 0;
            for (int i = 0; i < myAlignments.length; i++) {
                int currentNumSites = myAlignments[i].getSiteCount();
                for (int j = 0; j < currentNumSites; j++) {
                    result[count++] = myAlignments[i].getBase(taxon, j);
                }
            }
            return result;
        } else {
            return super.getBaseRange(taxon, startSite, endSite);
        }

    }

    @Override
    public byte getBase(int taxon, Locus locus, int physicalPosition) {
        int site = getSiteOfPhysicalPosition(physicalPosition, locus);
        int translate = translateSite(site);
        return myAlignments[translate].getBase(taxon, site - mySiteOffsets[translate]);
    }

    /**
     * Returns which alignment to use.
     *
     * @param site
     * @return alignment index.
     */
    public int translateSite(int site) {

        for (int i = 1; i < mySiteOffsets.length; i++) {
            if (mySiteOffsets[i] > site) {
                return i - 1;
            }
        }
        throw new IndexOutOfBoundsException("CombineAlignment: translateSite: index out of range: " + site);

    }

    @Override
    public boolean hasReference() {

        for (int i = 0; i < myAlignments.length; i++) {
            if (!myAlignments[i].hasReference()) {
                return false;
            }
        }

        return true;
    }

    @Override
    public String[] getSNPIDs() {

        int numSites = getSiteCount();
        String[] result = new String[numSites];
        int count = 0;
        for (int i = 0; i < myAlignments.length; i++) {
            for (int j = 0, n = myAlignments[i].getSiteCount(); j < n; j++) {
                result[count++] = myAlignments[i].getSNPID(j);
            }
        }

        return result;

    }

    @Override
    public String getSNPID(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getSNPID(site - mySiteOffsets[translate]);
    }

    @Override
    public int getSiteCount() {
        return mySiteOffsets[mySiteOffsets.length - 1];
    }

    @Override
    public int getLocusSiteCount(Locus locus) {
        return ((Alignment) myLoci.get(locus)).getLocusSiteCount(locus);
    }

    @Override
    public int getPositionInLocus(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getPositionInLocus(site - mySiteOffsets[translate]);
    }

    @Override
    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus) {
        Alignment align = ((Alignment) myLoci.get(locus));
        int i = -1;
        for (int j = 0; j < myAlignments.length; j++) {
            if (myAlignments[j] == align) {
                i = j;
                break;
            }
        }
        if (i == -1) {
            return -1;
        }
        return mySiteOffsets[i] + align.getSiteOfPhysicalPosition(physicalPosition, locus);
    }

    @Override
    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus, String snpID) {
        Alignment align = ((Alignment) myLoci.get(locus));
        int i = -1;
        for (int j = 0; j < myAlignments.length; j++) {
            if (myAlignments[j] == align) {
                i = j;
                break;
            }
        }
        if (i == -1) {
            return -1;
        }
        return mySiteOffsets[i] + align.getSiteOfPhysicalPosition(physicalPosition, locus, snpID);
    }

    @Override
    public Locus getLocus(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getLocus(site - mySiteOffsets[translate]);
    }

    @Override
    public Locus[] getLoci() {
        return myLociList;
    }

    @Override
    public int getNumLoci() {
        if (myLociList == null) {
            return 0;
        } else {
            return myLociList.length;
        }
    }

    @Override
    public float[][] getSiteScores() {

        if (!hasSiteScores()) {
            return null;
        }

        int numSeqs = getSequenceCount();
        float[][] result = new float[numSeqs][getSiteCount()];
        for (int a = 0, n = myAlignments.length; a < n; a++) {
            if (myAlignments[a].hasSiteScores()) {
                for (int s = 0, m = myAlignments[a].getSiteCount(); s < m; s++) {
                    for (int t = 0; t < numSeqs; t++) {
                        result[t][mySiteOffsets[a] + s] = myAlignments[a].getSiteScore(t, s);
                    }
                }
            }
        }

        return result;

    }

    @Override
    public float getSiteScore(int seq, int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getSiteScore(seq, site - mySiteOffsets[translate]);
    }

    @Override
    public boolean hasSiteScores() {
        for (Alignment align : myAlignments) {
            if (align.hasSiteScores()) {
                return true;
            }
        }
        return false;
    }

    @Override
    public int getIndelSize(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getIndelSize(site - mySiteOffsets[translate]);
    }

    @Override
    public boolean isIndel(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].isIndel(site - mySiteOffsets[translate]);
    }

    @Override
    public byte getReferenceAllele(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getReferenceAllele(site - mySiteOffsets[translate]);
    }

    @Override
    public Alignment[] getAlignments() {
        return myAlignments;
    }

    @Override
    public byte getMajorAllele(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getMajorAllele(site - mySiteOffsets[translate]);
    }

    @Override
    public byte getMinorAllele(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getMinorAllele(site - mySiteOffsets[translate]);
    }

    @Override
    public byte[] getMinorAlleles(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getMinorAlleles(site - mySiteOffsets[translate]);
    }

    @Override
    public byte[] getAlleles(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getAlleles(site - mySiteOffsets[translate]);
    }

    @Override
    public double getMinorAlleleFrequency(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getMinorAlleleFrequency(site - mySiteOffsets[translate]);
    }

    @Override
    public int[][] getAllelesSortedByFrequency(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getAllelesSortedByFrequency(site - mySiteOffsets[translate]);
    }

    @Override
    public byte[] getBaseArray(int taxon, int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getBaseArray(taxon, site - mySiteOffsets[translate]);
    }

    @Override
    public byte[] getBaseRow(int taxon) {
        byte[] result = new byte[getSiteCount()];
        for (int i = 0; i < myAlignments.length; i++) {
            byte[] current = myAlignments[i].getBaseRow(taxon);
            System.arraycopy(current, 0, result, myLociOffsets[i], current.length);
        }
        return result;
    }

    @Override
    public BitSet getAllelePresenceForAllSites(int taxon, int alleleNumber) {
        throw new UnsupportedOperationException("CombineAlignment: getAllelePresenceForAllSites: This operation isn't possible as it spans multiple Alignments. It needs to be optimized for taxa first.");
    }

    @Override
    public BitSet getAllelePresenceForAllTaxa(int site, int alleleNumber) {
        int translate = translateSite(site);
        return myAlignments[translate].getAllelePresenceForAllTaxa(site - mySiteOffsets[translate], alleleNumber);
    }

    @Override
    public long[] getAllelePresenceForSitesBlock(int taxon, int alleleNumber, int startBlock, int endBlock) {
        throw new UnsupportedOperationException("CombineAlignment: getAllelePresenceForSitesBlock: This operation isn't possible as it spans multiple Alignments. It needs to be optimized for taxa first.");
    }

    @Override
    public String getBaseAsString(int taxon, int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getBaseAsString(taxon, site - mySiteOffsets[translate]);
    }

    @Override
    public String[] getBaseAsStringArray(int taxon, int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getBaseAsStringArray(taxon, site - mySiteOffsets[translate]);
    }

    @Override
    public byte[] getReference(int startSite, int endSite) {
        int numSites = endSite - startSite;
        byte[] result = new byte[numSites];
        for (int i = 0; i < numSites; i++) {
            result[i] = getReferenceAllele(startSite + i);
        }
        return result;
    }

    @Override
    public byte[] getReference() {

        for (int i = 0; i < myAlignments.length; i++) {
            if (!myAlignments[i].hasReference()) {
                return null;
            }
        }

        byte[] result = new byte[getSiteCount()];
        int count = 0;
        for (int i = 0; i < myAlignments.length; i++) {
            byte[] current = myAlignments[i].getReference();
            for (int j = 0; j < current.length; j++) {
                result[count++] = current[j];
            }
        }
        return result;

    }

    @Override
    public boolean isHeterozygous(int taxon, int site) {
        int translate = translateSite(site);
        return myAlignments[translate].isHeterozygous(taxon, site - mySiteOffsets[translate]);
    }

    @Override
    public int[] getPhysicalPositions() {

        boolean allNull = true;
        for (int i = 0; i < myAlignments.length; i++) {
            int[] current = myAlignments[0].getPhysicalPositions();
            if ((current != null) && (current.length != 0)) {
                allNull = false;
                break;
            }
        }

        if (allNull) {
            return null;
        } else {
            int[] result = new int[getSiteCount()];
            int count = 0;
            for (int i = 0; i < myAlignments.length; i++) {
                int[] current = myAlignments[i].getPhysicalPositions();
                for (int j = 0; j < current.length; j++) {
                    result[count++] = current[j];
                }
            }
            return result;
        }
    }

    @Override
    public String getLocusName(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getLocusName(site - mySiteOffsets[translate]);
    }

    @Override
    public int[] getLociOffsets() {
        return myLociOffsets;
    }

    @Override
    public SITE_SCORE_TYPE getSiteScoreType() {
        SITE_SCORE_TYPE first = myAlignments[0].getSiteScoreType();
        for (int i = 1; i < myAlignments.length; i++) {
            if (first != myAlignments[i].getSiteScoreType()) {
                return SITE_SCORE_TYPE.MixedScoreTypes;
            }
        }
        return first;
    }

    @Override
    public boolean isAllPolymorphic() {
        for (int i = 0; i < myAlignments.length; i++) {
            if (!myAlignments[i].isAllPolymorphic()) {
                return false;
            }
        }
        return true;
    }

    @Override
    public boolean isPolymorphic(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].isPolymorphic(site - mySiteOffsets[translate]);
    }

    @Override
    public double getMajorAlleleFrequency(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getMajorAlleleFrequency(site - mySiteOffsets[translate]);
    }

    @Override
    public String getGenomeAssembly() {
        String first = myAlignments[0].getGenomeAssembly();
        if (first == null) {
            return null;
        }
        for (int i = 1; i < myAlignments.length; i++) {
            String current = myAlignments[i].getGenomeAssembly();
            if ((current != null) && (!first.equals(current))) {
                return null;
            }
        }
        return first;
    }

    @Override
    public boolean isPositiveStrand(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].isPositiveStrand(site - mySiteOffsets[translate]);
    }

    @Override
    public boolean isPhased() {
        for (int i = 0; i < myAlignments.length; i++) {
            if (myAlignments[i].isPhased() == false) {
                return false;
            }
        }
        return true;
    }

    @Override
    public GeneticMap getGeneticMap() {
        GeneticMap result = myAlignments[0].getGeneticMap();
        for (int i = 1; i < myAlignments.length; i++) {
            GeneticMap current = myAlignments[i].getGeneticMap();
            if ((current == null) || (!current.equals(result))) {
                return null;
            }
        }
        return result;
    }

    @Override
    public boolean retainsRareAlleles() {
        for (int i = 0; i < myAlignments.length; i++) {
            if (myAlignments[i].retainsRareAlleles() == false) {
                return false;
            }
        }
        return true;
    }

    @Override
    public String[][] getAlleleEncodings() {

        if (myAlleleStates != null) {
            return myAlleleStates;
        }

        boolean allTheSame = true;
        String[][] encodings = myAlignments[0].getAlleleEncodings();
        if (encodings.length == 1) {
            for (int i = 1; i < myAlignments.length; i++) {
                String[][] current = myAlignments[i].getAlleleEncodings();
                if ((current.length == 1) && (encodings[0].length == current[0].length)) {
                    for (int j = 0; j < encodings[0].length; j++) {
                        if (!current[0][j].equals(encodings[0][j])) {
                            allTheSame = false;
                            break;
                        }
                    }
                } else {
                    allTheSame = false;
                    break;
                }

                if (!allTheSame) {
                    break;
                }
            }
        } else {
            allTheSame = false;
        }

        if (allTheSame) {
            myAlleleStates = encodings;
        } else {
            String[][] result = new String[getSiteCount()][];
            int count = 0;
            for (int i = 0; i < myAlignments.length; i++) {
                for (int j = 0, n = myAlignments[i].getSiteCount(); j < n; j++) {
                    result[count++] = myAlignments[i].getAlleleEncodings(j);
                }
            }
            myAlleleStates = result;
        }

        return myAlleleStates;

    }

    @Override
    public String[] getAlleleEncodings(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getAlleleEncodings(site - mySiteOffsets[translate]);
    }

    @Override
    public String getBaseAsString(int site, byte value) {
        int translate = translateSite(site);
        return myAlignments[translate].getBaseAsString(site - mySiteOffsets[translate], value);
    }

    @Override
    public int getMaxNumAlleles() {
        int result = 999999;
        for (int i = 0; i < myAlignments.length; i++) {
            if (myAlignments[i].getMaxNumAlleles() < result) {
                result = myAlignments[i].getMaxNumAlleles();
            }
        }
        return result;
    }

    @Override
    public int getTotalGametesNotMissing(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getTotalGametesNotMissing(site - mySiteOffsets[translate]);
    }

    @Override
    public int getHeterozygousCount(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getHeterozygousCount(site - mySiteOffsets[translate]);
    }

    @Override
    public int getMinorAlleleCount(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getMinorAlleleCount(site - mySiteOffsets[translate]);
    }

    @Override
    public int getMajorAlleleCount(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getMajorAlleleCount(site - mySiteOffsets[translate]);
    }

    @Override
    public Object[][] getDiploidssSortedByFrequency(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getDiploidssSortedByFrequency(site - mySiteOffsets[translate]);
    }

    @Override
    public boolean isSBitFriendly() {
        for (int i = 0; i < myAlignments.length; i++) {
            if (!myAlignments[i].isSBitFriendly()) {
                return false;
            }
        }
        return true;
    }

    @Override
    public boolean isTBitFriendly() {
        return false;
    }

    @Override
    public void optimizeForTaxa(ProgressListener listener) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void optimizeForSites(ProgressListener listener) {
        for (int i = 0; i < myAlignments.length; i++) {
            myAlignments[i].optimizeForSites(listener);
        }
    }

    @Override
    public byte[] getAllelesByScope(ALLELE_SCOPE_TYPE scope, int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getAllelesByScope(scope, site - mySiteOffsets[translate]);
    }

    @Override
    public BitSet getAllelePresenceForAllTaxaByScope(ALLELE_SCOPE_TYPE scope, int site, int alleleNumber) {
        int translate = translateSite(site);
        return myAlignments[translate].getAllelePresenceForAllTaxaByScope(scope, site - mySiteOffsets[translate], alleleNumber);
    }
}
