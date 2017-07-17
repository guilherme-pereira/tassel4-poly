/*
 * MutableSingleEncodeAlignment
 */
package net.maizegenetics.pal.alignment;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.TreeSet;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import org.apache.log4j.Logger;

/**
 *
 * @author terry, Gabriel Rodrigues Alves Margarido
 */
public class MutableSingleEncodeAlignment extends AbstractAlignment implements MutableAlignment {

    private static final Logger myLogger = Logger.getLogger(MutableSingleEncodeAlignment.class);
    private boolean myIsDirty = true;
    private byte[][] myData;
    private List<Identifier> myIdentifiers = new ArrayList<Identifier>();
    private final int myMaxTaxa;
    private final int myMaxNumSites;
    private int myNumSites = 0;
    private int myNumSitesStagedToRemove = 0;
    protected int[] myVariableSites;
    private List<Locus> myLocusToLociIndex = new ArrayList<Locus>();
    protected int[] myLocusIndices;
    private int[] myLocusOffsets = null;
    protected String[] mySNPIDs;

    protected MutableSingleEncodeAlignment(Alignment a, int maxNumTaxa, int maxNumSites) {
        super(a.getAlleleEncodings());
        myMaxNumAlleles = a.getMaxNumAlleles();

        if (a.getAlleleEncodings().length != 1) {
            throw new IllegalArgumentException("MutableSingleEncodeAlignment: init: must only have one allele encoding.");
        }

        if (a.getSiteCount() > maxNumSites) {
            throw new IllegalArgumentException("MutableSingleEncodeAlignment: init: initial number of sites can't be more than max number of sites.");
        }

        if (a.getSequenceCount() > maxNumTaxa) {
            throw new IllegalArgumentException("MutableSingleEncodeAlignment: init: initial number of taxa can't be more than max number of taxa.");
        }

        myMaxTaxa = maxNumTaxa;
        myMaxNumSites = maxNumSites;
        myNumSites = a.getSiteCount();
        myReference = new byte[myMaxNumSites];
        Arrays.fill(myReference, Alignment.UNKNOWN_DIPLOID_ALLELE);
        initData();
        initTaxa(a.getIdGroup());
        loadAlleles(a);
        loadLoci(a);
        System.arraycopy(a.getSNPIDs(), 0, mySNPIDs, 0, a.getSiteCount());
        System.arraycopy(a.getPhysicalPositions(), 0, myVariableSites, 0, a.getSiteCount());
    }

    public static MutableSingleEncodeAlignment getInstance(Alignment a, int maxTaxa, int maxNumSites) {
        if (a.getAlleleEncodings() == NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES) {
            return MutableNucleotideAlignment.getInstance(a, maxTaxa, maxNumSites);
        } else {
            return MutableSingleEncodeAlignment.getInstance(a, maxTaxa, maxNumSites);
        }
    }

    public static MutableSingleEncodeAlignment getInstance(Alignment a) {
        return getInstance(a, a.getSequenceCount(), a.getSiteCount());
    }

    protected MutableSingleEncodeAlignment(String[][] encodings, IdGroup idGroup, int initNumSites, int maxNumTaxa, int maxNumSites) {
        super(encodings);

        if (initNumSites > maxNumSites) {
            throw new IllegalArgumentException("MutableNucleotideAlignment: init: initial number of sites can't be more than max number of sites.");
        }

        if (idGroup.getIdCount() > maxNumTaxa) {
            throw new IllegalArgumentException("MutableNucleotideAlignment: init: initial number of taxa can't be more than max number of taxa.");
        }

        myMaxTaxa = maxNumTaxa;
        myMaxNumSites = maxNumSites;
        myNumSites = initNumSites;
        myReference = new byte[myMaxNumSites];
        Arrays.fill(myReference, Alignment.UNKNOWN_DIPLOID_ALLELE);

        initData();
        initTaxa(idGroup);
    }

    protected MutableSingleEncodeAlignment(String[][] encodings, List<Identifier> idGroup, int[] variableSites, List<Locus> locusToLociIndex, int[] locusIndices, String[] siteNames) {
        super(encodings);

        if ((variableSites.length != locusIndices.length) || (variableSites.length != siteNames.length)) {
            throw new IllegalArgumentException("MutableSingleEncodeAlignment: init: number variable sites, loci, and site names must be same.");
        }

        myMaxTaxa = idGroup.size();
        myMaxNumSites = siteNames.length;
        myNumSites = siteNames.length;

        myVariableSites = variableSites;
        myLocusToLociIndex = locusToLociIndex;
        myLocusIndices = locusIndices;
        mySNPIDs = siteNames;

        myData = new byte[myMaxTaxa][myMaxNumSites];
        for (int t = 0; t < myMaxTaxa; t++) {
            Arrays.fill(myData[t], Alignment.UNKNOWN_DIPLOID_ALLELE);
        }

        myReference = new byte[myMaxNumSites];
        Arrays.fill(myReference, Alignment.UNKNOWN_DIPLOID_ALLELE);

        myIdentifiers = idGroup;
    }

    public static MutableSingleEncodeAlignment getInstance(String[][] encodings, IdGroup idGroup, int initNumSites, int maxNumTaxa, int maxNumSites) {
        return new MutableSingleEncodeAlignment(encodings, idGroup, initNumSites, maxNumTaxa, maxNumSites);
    }

    public static MutableSingleEncodeAlignment getInstance(String[][] encodings, IdGroup idGroup, int initNumSites) {
        return MutableSingleEncodeAlignment.getInstance(encodings, idGroup, initNumSites, idGroup.getIdCount(), initNumSites);
    }

    public static MutableSingleEncodeAlignment getInstance(String[][] encodings, List<Identifier> idGroup, int[] variableSites, List<Locus> locusToLociIndex, int[] locusIndices, String[] siteNames) {
        return new MutableSingleEncodeAlignment(encodings, idGroup, variableSites, locusToLociIndex, locusIndices, siteNames);
    }

    public static MutableSingleEncodeAlignment getInstance(Alignment[] alignments) {

        if ((alignments == null) || (alignments.length == 0)) {
            return null;
        }

        String[][] resultEncodings = alignments[0].getAlleleEncodings();
        if (resultEncodings.length != 1) {
            throw new IllegalArgumentException("MutableSingleEncodeAlignment: getInstance: Alignments must have single allele encoding.");
        }

        String[][][] encodings = new String[alignments.length][][];
        for (int i = 0; i < alignments.length; i++) {
            encodings[i] = alignments[i].getAlleleEncodings();
        }
        if (!AlignmentUtils.areEncodingsEqual(encodings)) {
            throw new IllegalArgumentException("MutableSingleEncodeAlignment: getInstance: Alignments must have same allele encoding.");
        }

        TreeSet<Identifier> taxa = new TreeSet<Identifier>();
        List<Integer> physicalPositions = new ArrayList<Integer>();
        List<Locus> locusToLociIndex = new ArrayList<Locus>();
        List<Integer> locusIndices = new ArrayList<Integer>();
        LinkedHashMap<String, Integer> siteNameIndex = new LinkedHashMap<>();

        for (int i = 0; i < alignments.length; i++) {

            IdGroup currentIds = alignments[i].getIdGroup();
            for (int j = 0, n = currentIds.getIdCount(); j < n; j++) {
                Identifier current = currentIds.getIdentifier(j);
                if (taxa.contains(current)) {
                    Identifier match = taxa.floor(current);
                    Identifier merged = Identifier.getMergedInstance(match, current);
                    taxa.remove(match);
                    taxa.add(merged);
                } else {
                    taxa.add(current);
                }
            }

            for (int s = 0, m = alignments[i].getSiteCount(); s < m; s++) {
                String currentSiteName = alignments[i].getSNPID(s);
                int currentPhysicalPos = alignments[i].getPositionInLocus(s);
                Locus currentLocus = alignments[i].getLocus(s);
                int index = -1;
                Integer temp = siteNameIndex.get(currentSiteName);
                if (temp != null) {
                    index = temp;
                }
                if (index == -1) {
                    siteNameIndex.put(currentSiteName, siteNameIndex.size());
                    physicalPositions.add(currentPhysicalPos);
                    int locusIndex = -1;
                    for (int li = 0; li < locusToLociIndex.size(); li++) {
                        if (currentLocus.getChromosomeName().equals(locusToLociIndex.get(li).getChromosomeName())) {
                            locusIndex = li;
                            locusToLociIndex.set(li, Locus.getMergedInstance(currentLocus, locusToLociIndex.get(li)));
                            break;
                        }
                    }
                    if (locusIndex == -1) {
                        locusIndices.add(locusToLociIndex.size());
                        locusToLociIndex.add(currentLocus);
                    } else {
                        locusIndices.add(locusIndex);
                    }
                } else {
                    if (i == 0) {
                        throw new IllegalStateException("MutableSingleEncodeAlignment: getInstance: Duplicate site name in alignment: " + currentSiteName);
                    } else {
                        if (currentPhysicalPos != physicalPositions.get(index)) {
                            throw new IllegalStateException("MutableSingleEncodeAlignment: getInstance: Physical Positions do not match for site name: " + currentSiteName);
                        }
                        int locusIndex = -1;
                        for (int li = 0; li < locusToLociIndex.size(); li++) {
                            if (currentLocus.getChromosomeName().equals(locusToLociIndex.get(li).getChromosomeName())) {
                                locusIndex = li;
                                break;
                            }
                        }
                        if (locusIndices.get(index) != locusIndex) {
                            throw new IllegalStateException("MutableSingleEncodeAlignment: getInstance: Loci do not match for site name: " + currentSiteName + " expecting: " + locusToLociIndex.get(locusIndices.get(index)) + " but doesn't match: " + currentLocus);
                        }
                    }
                }
            }

        }

        int[] variableSites = new int[physicalPositions.size()];
        for (int i = 0, n = physicalPositions.size(); i < n; i++) {
            variableSites[i] = physicalPositions.get(i);
        }

        int[] locusIndicesArray = new int[locusIndices.size()];
        for (int i = 0, n = locusIndices.size(); i < n; i++) {
            locusIndicesArray[i] = (int) locusIndices.get(i);
        }

        String[] siteNamesArray = new String[siteNameIndex.size()];
        siteNameIndex.keySet().toArray(siteNamesArray);

        List taxaList = new ArrayList<Identifier>(taxa);
        MutableSingleEncodeAlignment result = null;

        encodings = new String[2][][];
        encodings[0] = NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES;
        encodings[1] = resultEncodings;
        if (AlignmentUtils.areEncodingsEqual(encodings)) {
            result = MutableNucleotideAlignment.getInstance(taxaList, variableSites, locusToLociIndex, locusIndicesArray, siteNamesArray);
        } else {
            result = getInstance(resultEncodings, taxaList, variableSites, locusToLociIndex, locusIndicesArray, siteNamesArray);
        }
        result.sortSitesByPhysicalPositionEmptyData();
        result.setClean();

        for (int i = 0; i < alignments.length; i++) {
            myLogger.info("Merging Alignment: " + (i + 1) + " of " + alignments.length);
            Alignment currentAlignment = alignments[i];
            IdGroup ids = currentAlignment.getIdGroup();
            int numSeqs = ids.getIdCount();
            int[] taxaIndices = new int[numSeqs];
            for (int t = 0; t < numSeqs; t++) {
                taxaIndices[t] = taxaList.indexOf(ids.getIdentifier(t));
            }
            for (int s = 0, n = currentAlignment.getSiteCount(); s < n; s++) {
                String siteName = currentAlignment.getSNPID(s);
                int physicalPosition = currentAlignment.getPositionInLocus(s);
                Locus locus = currentAlignment.getLocus(s);
                for (int li = 0; li < locusToLociIndex.size(); li++) {
                    if (locus.getChromosomeName().equals(locusToLociIndex.get(li).getChromosomeName())) {
                        locus = locusToLociIndex.get(li);
                        break;
                    }
                }
                int site = result.getSiteOfPhysicalPosition(physicalPosition, locus, siteName);
                if (site < 0) {
                    throw new IllegalStateException("MutableSingleEncodeAlignment: getInstance: physical position: " + physicalPosition + " in locus: " + locus.getName() + " not found.");
                }

                for (int t = 0; t < numSeqs; t++) {
                    result.setBase(taxaIndices[t], site, currentAlignment.getBase(t, s));
                }
            }
        }

        return result;
    }

    private void initData() {
        myData = new byte[myMaxTaxa][myMaxNumSites];
        for (int t = 0; t < myMaxTaxa; t++) {
            Arrays.fill(myData[t], Alignment.UNKNOWN_DIPLOID_ALLELE);
        }
        myLocusIndices = new int[myMaxNumSites];
        Arrays.fill(myLocusIndices, Integer.MAX_VALUE);
        myVariableSites = new int[myMaxNumSites];
        Arrays.fill(myVariableSites, -1);
        mySNPIDs = new String[myMaxNumSites];
        Arrays.fill(mySNPIDs, null);
    }

    private void initTaxa(IdGroup idGroup) {
        for (int i = 0, n = idGroup.getIdCount(); i < n; i++) {
            myIdentifiers.add(idGroup.getIdentifier(i));
        }
    }

    private void loadAlleles(Alignment a) {

        int numSites = a.getSiteCount();
        int numSeqs = a.getSequenceCount();

        for (int s = 0; s < numSites; s++) {
            for (int t = 0; t < numSeqs; t++) {
                myData[t][s] = a.getBase(t, s);
            }
        }

    }

    private void loadLoci(Alignment a) {

        Locus[] loci = a.getLoci();
        for (int i = 0; i < loci.length; i++) {
            myLocusToLociIndex.add(loci[i]);
        }

        int[] offsets = a.getLociOffsets();
        for (int i = 0; i < offsets.length - 1; i++) {
            for (int j = offsets[i]; j < offsets[i + 1]; j++) {
                myLocusIndices[j] = i;
            }
        }
        for (int j = offsets[offsets.length - 1], n = a.getSiteCount(); j < n; j++) {
            myLocusIndices[j] = offsets.length - 1;
        }
    }

    public byte getBase(int taxon, int site) {
        return myData[taxon][site];
    }

    public boolean isSBitFriendly() {
        return false;
    }

    public boolean isTBitFriendly() {
        return false;
    }

    @Override
    public int getSiteCount() {
        return myNumSites;
    }

    @Override
    public int getSequenceCount() {
        return myIdentifiers.size();
    }

    @Override
    public IdGroup getIdGroup() {
        Identifier[] ids = new Identifier[myIdentifiers.size()];
        myIdentifiers.toArray(ids);
        return new SimpleIdGroup(ids);
    }

    @Override
    public String getTaxaName(int index) {
        return myIdentifiers.get(index).getName();
    }

    @Override
    public String getFullTaxaName(int index) {
        return myIdentifiers.get(index).getFullName();
    }

    @Override
    public boolean isPhased() {
        return false;
    }

    @Override
    public boolean isPositiveStrand(int site) {
        return true;
    }

    @Override
    public String getGenomeAssembly() {
        return "AGPV2";
    }

    @Override
    public int[] getPhysicalPositions() {
        return myVariableSites.clone();
    }

    @Override
    public int getPositionInLocus(int site) {
        try {
            if (myVariableSites[site] < 0) {
                return site;
            }
            return myVariableSites[site];
        } catch (Exception e) {
            return site;
        }
    }

    @Override
    public Locus getLocus(int site) {
        return myLocusToLociIndex.get(myLocusIndices[site]);
    }

    @Override
    public Locus[] getLoci() {
        Locus[] result = new Locus[myLocusToLociIndex.size()];
        for (int i = 0; i < myLocusToLociIndex.size(); i++) {
            result[i] = myLocusToLociIndex.get(i);
        }
        return result;
    }

    @Override
    public int getNumLoci() {
        return myLocusToLociIndex.size();
    }

    @Override
    public int[] getLociOffsets() {

        if (isDirty()) {
            throw new IllegalStateException("MutableSingleEncodeAlignment: getLociOffsets: this alignment is dirty.");
        }

        if (myLocusOffsets != null) {
            return myLocusOffsets;
        }

        List<Integer> result = new ArrayList<Integer>();
        int current = myLocusIndices[0];
        result.add(0);
        for (int i = 0, n = getSiteCount(); i < n; i++) {
            if (myLocusIndices[i] != current) {
                result.add(i);
                current = myLocusIndices[i];
            }
        }
        myLocusOffsets = new int[result.size()];
        for (int i = 0, n = result.size(); i < n; i++) {
            myLocusOffsets[i] = result.get(i);
        }
        return myLocusOffsets;

    }

    @Override
    public int[] getStartAndEndOfLocus(Locus locus) {

        if (isDirty()) {
            throw new IllegalStateException("MutableSingleEncodeAlignment: getStartAndEndOfLocus: this alignment is dirty.");
        }

        Locus[] loci = getLoci();
        int[] lociOffsets = getLociOffsets();
        int numLoci = getNumLoci();
        for (int i = 0; i < numLoci; i++) {
            if (locus.equals(loci[i])) {
                int end = 0;
                if (i == numLoci - 1) {
                    end = getSiteCount();
                } else {
                    end = lociOffsets[i + 1];
                }
                return new int[]{lociOffsets[i], end};
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
        if ((mySNPIDs == null) || (mySNPIDs.length == 0) || (mySNPIDs[site] == null)) {
            return "S" + getLocus(site).getChromosomeName() + "_" + getPositionInLocus(site);
        }
        return mySNPIDs[site];
    }

    @Override
    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus) {

        if (isDirty()) {
            throw new IllegalStateException("MutableSingleEncodeAlignment: getSiteOfPhysicalPosition: this alignment is dirty.");
        }

        if (myVariableSites == null) {
            return physicalPosition;
        }
        try {
            if (locus == null) {
                locus = myLocusToLociIndex.get(0);
            }
            int[] startEnd = getStartAndEndOfLocus(locus);
            return Arrays.binarySearch(myVariableSites, startEnd[0], startEnd[1], physicalPosition);
        } catch (Exception e) {
            e.printStackTrace();
            return -1;
        }
    }

    @Override
    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus, String snpID) {

        if (isDirty()) {
            throw new IllegalStateException("MutableSingleEncodeAlignment: getSiteOfPhysicalPosition: this alignment is dirty.");
        }

        if (myVariableSites == null) {
            return physicalPosition;
        }
        if (locus == null) {
            locus = myLocusToLociIndex.get(0);
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
        return false;
    }

    @Override
    public byte[] getAlleles(int site) {
        int[][] alleles = getAllelesSortedByFrequency(site);
        int resultSize = alleles[0].length;
        byte[] result = new byte[resultSize];
        for (int i = 0; i < resultSize; i++) {
            result[i] = (byte) alleles[0][i];
        }
        return result;
    }

    // Mutable Methods...
    public void setBase(int taxon, int site, byte newBase) {
        myData[taxon][site] = newBase;
    }

    public void setBase(Identifier taxon, String siteName, Locus locus, int physicalPosition, byte newBase) {

        int taxonIndex = myIdentifiers.indexOf(taxon);
        if (taxonIndex == -1) {
            throw new IllegalArgumentException("MutableSingleEncodeAlignment: setBase: taxon not found.");
        }

        int site = getSiteOfPhysicalPosition(physicalPosition, locus, siteName);
        if (site < 0) {
            throw new IllegalStateException("MutableSingleEncodeAlignment: setBase: physical position: " + physicalPosition + " in locus: " + locus.getName() + " not found.");
        } else {
            if (!siteName.equals(getSNPID(site))) {
                throw new IllegalStateException("MutableSingleEncodeAlignment: setBase: site names at physical position: " + physicalPosition + " in locus: " + locus.getName() + " does not match: " + siteName);
            }
        }

        myData[taxonIndex][site] = newBase;

    }

    public void setBaseRange(int taxon, int startSite, byte[] newBases) {
        for (int i = 0; i < newBases.length; i++) {
            myData[taxon][startSite++] = newBases[i];
        }
    }

    public void setReferenceAllele(int site, byte diploidAllele) {
        myReference[site] = diploidAllele;
    }

    public void addSite(int site) {

        if (myMaxNumSites < myNumSites + 1) {
            throw new IllegalStateException("MutableSingleEncodeAlignment: addSite: this exceeds max num of sites: " + myMaxNumSites);
        }

        for (int t = 0, n = getSequenceCount(); t < n; t++) {
            for (int s = myNumSites; s > site; s--) {
                myData[t][s] = myData[t][s - 1];
            }
            myData[t][site] = Alignment.UNKNOWN_DIPLOID_ALLELE;
        }

        for (int s = myNumSites; s > site; s--) {
            myVariableSites[s] = myVariableSites[s - 1];
            myLocusIndices[s] = myLocusIndices[s - 1];
            mySNPIDs[s] = mySNPIDs[s - 1];
            myReference[s] = myReference[s - 1];
        }
        myVariableSites[site] = -1;
        myLocusIndices[site] = Integer.MAX_VALUE;
        mySNPIDs[site] = null;
        myReference[site] = Alignment.UNKNOWN_DIPLOID_ALLELE;

        myNumSites++;

        setDirty();

    }

    @Override
    public void setSNPID(int site, String name) {
        mySNPIDs[site] = name;
    }

    public void removeSite(int site) {

        myNumSites--;

        for (int t = 0, n = getSequenceCount(); t < n; t++) {
            for (int s = site; s < myNumSites; s++) {
                myData[t][s] = myData[t][s + 1];
            }
            myData[t][myNumSites] = Alignment.UNKNOWN_DIPLOID_ALLELE;
        }

        for (int s = site; s < myNumSites; s++) {
            myVariableSites[s] = myVariableSites[s + 1];
            myLocusIndices[s] = myLocusIndices[s + 1];
            mySNPIDs[s] = mySNPIDs[s + 1];
            myReference[s] = myReference[s + 1];
        }
        myVariableSites[myNumSites] = -1;
        myLocusIndices[myNumSites] = Integer.MAX_VALUE;
        mySNPIDs[myNumSites] = null;
        myReference[myNumSites] = Alignment.UNKNOWN_DIPLOID_ALLELE;

    }

    public void clearSiteForRemoval(int site) {

        myNumSitesStagedToRemove++;

        for (int t = 0, n = getSequenceCount(); t < n; t++) {
            myData[t][site] = Alignment.UNKNOWN_DIPLOID_ALLELE;
        }

        myVariableSites[site] = Integer.MAX_VALUE;
        myLocusIndices[site] = Integer.MAX_VALUE;
        mySNPIDs[site] = null;
        myReference[site] = Alignment.UNKNOWN_DIPLOID_ALLELE;

    }

    public void addTaxon(Identifier id) {
        if (getSequenceCount() + 1 > myMaxTaxa) {
            throw new IllegalStateException("MutableSingleEncodeAlignment: addTaxon: this exceeds max num of taxa: " + myMaxTaxa);
        }
        myIdentifiers.add(id);
    }

    public void setTaxonName(int taxon, Identifier id) {
        if (taxon >= myIdentifiers.size()) {
            throw new IllegalStateException("MutableSingleEncodeAlignment: setTaxonName: this taxa index does not exist: " + taxon);
        }
        myIdentifiers.set(taxon, id);
    }

    public void removeTaxon(int taxon) {

        myIdentifiers.remove(taxon);

        int numTaxa = getSequenceCount();
        for (int s = 0; s < myNumSites; s++) {
            for (int t = taxon; t < numTaxa; t--) {
                myData[t][s] = myData[t + 1][s];
            }
            myData[numTaxa][s] = Alignment.UNKNOWN_DIPLOID_ALLELE;
        }

    }

    public void clean() {
        sortSitesByPhysicalPosition();
        removeUnusedLoci();
        myIsDirty = false;
        myNumSites -= myNumSitesStagedToRemove;
        myNumSitesStagedToRemove = 0;
    }

    public boolean isDirty() {
        return myIsDirty;
    }

    private void setDirty() {
        myLocusOffsets = null;
        myIsDirty = true;
    }

    private void setClean() {
        myIsDirty = false;
    }

    private void removeUnusedLoci() {
        boolean[] isUsed = new boolean[myLocusToLociIndex.size()];
        Arrays.fill(isUsed, false);
        for (int i = 0; i < myLocusIndices.length; i++) {
            if ((myLocusIndices[i] >= 0) && (myLocusIndices[i] != Integer.MAX_VALUE)) {
                isUsed[myLocusIndices[i]] = true;
            }
        }
        int removed = 0;
        for (int i = 0; i < isUsed.length; i++) {
            if (!isUsed[i]) {
                myLocusToLociIndex.remove(i - removed);
                decrementAllHigherLociIndices(i - removed);
                removed++;
            }
        }

    }

    private void decrementAllHigherLociIndices(int starting) {
        for (int i = 0; i < myLocusIndices.length; i++) {
            if (myLocusIndices[i] > starting) {
                myLocusIndices[i]--;
            }
        }
    }

    protected void sortSitesByPhysicalPosition() {

        Swapper swapperPos = new Swapper() {
            public void swap(int a, int b) {
                int it;
                it = myLocusIndices[a];
                myLocusIndices[a] = myLocusIndices[b];
                myLocusIndices[b] = it;

                byte bt;
                for (int t = 0, n = getSequenceCount(); t < n; t++) {
                    bt = getBase(t, a);
                    setBase(t, a, getBase(t, b));
                    setBase(t, b, bt);
                }

                it = myVariableSites[a];
                myVariableSites[a] = myVariableSites[b];
                myVariableSites[b] = it;

                String st = mySNPIDs[a];
                mySNPIDs[a] = mySNPIDs[b];
                mySNPIDs[b] = st;

                bt = myReference[a];
                myReference[a] = myReference[b];
                myReference[b] = bt;
            }
        };
        IntComparator compPos = new IntComparator() {
            public int compare(int a, int b) {
                if (myLocusIndices[a] < myLocusIndices[b]) {
                    return -1;
                }
                if (myLocusIndices[a] > myLocusIndices[b]) {
                    return 1;
                }
                if (myVariableSites[a] < myVariableSites[b]) {
                    return -1;
                }
                if (myVariableSites[a] > myVariableSites[b]) {
                    return 1;
                }
                return 0;
            }
        };

        GenericSorting.quickSort(0, getSiteCount(), compPos, swapperPos);

    }

    private void sortSitesByPhysicalPositionEmptyData() {

        Swapper swapperPos = new Swapper() {
            public void swap(int a, int b) {
                int it;
                it = myLocusIndices[a];
                myLocusIndices[a] = myLocusIndices[b];
                myLocusIndices[b] = it;

                it = myVariableSites[a];
                myVariableSites[a] = myVariableSites[b];
                myVariableSites[b] = it;

                String st = mySNPIDs[a];
                mySNPIDs[a] = mySNPIDs[b];
                mySNPIDs[b] = st;

                byte bt = myReference[a];
                myReference[a] = myReference[b];
                myReference[b] = bt;
            }
        };
        IntComparator compPos = new IntComparator() {
            public int compare(int a, int b) {
                if (myLocusIndices[a] < myLocusIndices[b]) {
                    return -1;
                }
                if (myLocusIndices[a] > myLocusIndices[b]) {
                    return 1;
                }
                if (myVariableSites[a] < myVariableSites[b]) {
                    return -1;
                }
                if (myVariableSites[a] > myVariableSites[b]) {
                    return 1;
                }
                return 0;
            }
        };

        GenericSorting.quickSort(0, this.getSiteCount(), compPos, swapperPos);

    }

    public void setPositionOfSite(int site, int position) {
        if ((site < 0) || (site >= getSiteCount())) {
            throw new IllegalArgumentException("MutableSingleEncodeAlignment: setPositionOfSite: site outside of range: " + site);
        }
        myVariableSites[site] = position;

        setDirty();
    }

    public void setLocusOfSite(int site, Locus locus) {
        if ((site < 0) || (site >= getSiteCount())) {
            throw new IllegalArgumentException("MutableSingleEncodeAlignment: setLocusOfSite: site outside of range: " + site);
        }
        int index = getLocusIndex(locus);
        if (index < 0) {
            myLocusToLociIndex.add(locus);
            myLocusIndices[site] = myLocusToLociIndex.size() - 1;
        } else {
            myLocusIndices[site] = index;
        }

        setDirty();
    }

    private int getLocusIndex(Locus locus) {
        for (int i = 0; i < myLocusToLociIndex.size(); i++) {
            if (myLocusToLociIndex.get(i).equals(locus)) {
                return i;
            }
        }
        return -1;
    }

    public void setDepthForAlleles(int taxon, int site, short[] values) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public void setCommonAlleles(int site, byte[] values) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
