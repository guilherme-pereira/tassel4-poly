/*
 * CoreAlignment
 */
package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.alignment.score.SiteScore;
import net.maizegenetics.pal.alignment.genotype.Genotype;
import net.maizegenetics.pal.alignment.bit.BitStorage;
import net.maizegenetics.pal.alignment.bit.DynamicBitStorage;
import net.maizegenetics.pal.alignment.depth.AlleleDepth;
import net.maizegenetics.pal.site.PositionList;
import net.maizegenetics.pal.site.Chromosome;
import net.maizegenetics.pal.taxa.TaxaList;
import net.maizegenetics.util.BitSet;
import org.apache.log4j.Logger;

/**
 *
 * @author terry
 */
public class CoreAlignment implements AlignmentNew {

    private static final Logger myLogger = Logger.getLogger(CoreAlignment.class);
    private final Genotype myGenotype;
    private BitStorage myFreqBitStorage;
    private BitStorage myReferenceBitStorage;
    private final PositionList myPositionList;
    private final TaxaList myTaxaList;
    private final SiteScore mySiteScore;
    private final AlleleDepth myAlleleDepth;

    public CoreAlignment(Genotype genotype, PositionList positionList, TaxaList taxaList, SiteScore siteScore, AlleleDepth alleleDepth) {
        myGenotype = genotype;
        myPositionList=positionList;
        myTaxaList = taxaList;
        mySiteScore = siteScore;
        myAlleleDepth = alleleDepth;
    }

    @Override
    public byte getBase(int taxon, int site) {
        return myGenotype.getBase(taxon, site);
    }

    @Override
    public byte[] getBaseArray(int taxon, int site) {
        return myGenotype.getBaseArray(taxon, site);
    }

    @Override
    public byte getBase(int taxon, Chromosome chromosome, int physicalPosition) {
        return myGenotype.getBase(taxon, myPositionList.getSiteOfPhysicalPosition(physicalPosition, chromosome));
    }

    @Override
    public byte[] getBaseRange(int taxon, int startSite, int endSite) {
        return myGenotype.getBaseRange(taxon, startSite, endSite);
    }

    @Override
    public byte[] getBaseRow(int taxon) {
        return myGenotype.getBaseRow(taxon);
    }

    @Override
    public BitSet getAllelePresenceForAllSites(int taxon, int alleleNumber) {
        return getBitStorage(ALLELE_SCOPE_TYPE.Frequency).getAllelePresenceForAllSites(taxon, alleleNumber);
    }

    @Override
    public BitSet getAllelePresenceForAllTaxa(int site, int alleleNumber) {
        return getBitStorage(ALLELE_SCOPE_TYPE.Frequency).getAllelePresenceForAllTaxa(site, alleleNumber);
    }

    @Override
    public long[] getAllelePresenceForSitesBlock(int taxon, int alleleNumber, int startBlock, int endBlock) {
        return getBitStorage(ALLELE_SCOPE_TYPE.Frequency).getAllelePresenceForSitesBlock(taxon, alleleNumber, startBlock, endBlock);
    }

    @Override
    public BitSet getPhasedAllelePresenceForAllSites(int taxon, boolean firstParent, int alleleNumber) {
        return getBitStorage(ALLELE_SCOPE_TYPE.Frequency).getPhasedAllelePresenceForAllSites(taxon, firstParent, alleleNumber);
    }

    @Override
    public BitSet getPhasedAllelePresenceForAllTaxa(int site, boolean firstParent, int alleleNumber) {
        return getBitStorage(ALLELE_SCOPE_TYPE.Frequency).getPhasedAllelePresenceForAllTaxa(site, firstParent, alleleNumber);
    }

    @Override
    public long[] getPhasedAllelePresenceForSitesBlock(int taxon, boolean firstParent, int alleleNumber, int startBlock, int endBlock) {
        return getBitStorage(ALLELE_SCOPE_TYPE.Frequency).getPhasedAllelePresenceForSitesBlock(taxon, firstParent, alleleNumber, startBlock, endBlock);
    }

    @Override
    public String getBaseAsString(int taxon, int site) {
        return myGenotype.getBaseAsString(taxon, site);
    }

    @Override
    public String getBaseAsStringRange(int taxon, int startSite, int endSite) {
        return myGenotype.getBaseAsStringRange(taxon, startSite, endSite);
    }

    @Override
    public String getBaseAsStringRow(int taxon) {
        return myGenotype.getBaseAsStringRow(taxon);
    }

    @Override
    public String[] getBaseAsStringArray(int taxon, int site) {
        return myGenotype.getBaseAsStringArray(taxon, site);
    }

    @Override
    public byte getReferenceAllele(int site) {
        return myPositionList.getReferenceAllele(site);
    }

    @Override
    public byte[] getReference(int startSite, int endSite) {
        return myPositionList.getReference(startSite, endSite);
    }

    @Override
    public byte[] getReference() {
        return myPositionList.getReference();
    }

    @Override
    public boolean hasReference() {
        return myPositionList.hasReference();
    }

    @Override
    public boolean isHeterozygous(int taxon, int site) {
        return myGenotype.isHeterozygous(taxon, site);
    }

    @Override
    public int getHeterozygousCount(int site) {
        return myGenotype.getHeterozygousCount(site);
    }

    @Override
    public PositionList getPositionList() {
        return myPositionList;
    }

    @Override
    public String[] getSNPIDs() {
        return myPositionList.getSNPIDs();
    }

    @Override
    public String getSNPID(int site) {
        return myPositionList.getSNPID(site);
    }

    @Override
    public int getSiteCount() {
        return myPositionList.getSiteCount();
    }

    @Override
    public int getChromosomeSiteCount(Chromosome chromosome) {
        return myPositionList.getChromosomeSiteCount(chromosome);
    }

    @Override
    public int[] getStartAndEndOfChromosome(Chromosome chromosome) {
        return myPositionList.getStartAndEndOfChromosome(chromosome);
    }

    @Override
    public int getSequenceCount() {
        return myTaxaList.getTaxaCount();
    }

    @Override
    public int getTaxaCount() {
        return myTaxaList.getTaxaCount();
    }

    @Override
    public int getPositionInChromosome(int site) {
        return myPositionList.getPositionInChromosome(site);
    }

    @Override
    public int getSiteOfPhysicalPosition(int physicalPosition, Chromosome chromosome) {
        return myPositionList.getSiteOfPhysicalPosition(physicalPosition, chromosome);
    }

    @Override
    public int getSiteOfPhysicalPosition(int physicalPosition, Chromosome chromosome, String snpID) {
        return myPositionList.getSiteOfPhysicalPosition(physicalPosition, chromosome, snpID);
    }

    @Override
    public int[] getPhysicalPositions() {
        return myPositionList.getPhysicalPositions();
    }

    @Override
    public String getChromosomeName(int site) {
        return myPositionList.getChromosomeName(site);
    }

    @Override
    public Chromosome getChromosome(int site) {
        return myPositionList.getChromosome(site);
    }

    @Override
    public Chromosome getChromosome(String name) {
        return myPositionList.getChromosome(name);
    }

    @Override
    public Chromosome[] getChromosomes() {
        return myPositionList.getChromosomes();
    }

    @Override
    public int getNumChromosomes() {
        return myPositionList.getNumChromosomes();
    }

    @Override
    public int[] getChromosomesOffsets() {
        return myPositionList.getChromosomesOffsets();
    }

    @Override
    public float getSiteScore(int seq, int site) {
        if (mySiteScore == null) {
            throw new IllegalStateException("CoreAlignment: getSiteScore: This Alignment has no Site Scores.");
        }
        return mySiteScore.getSiteScore(seq, site);
    }

    @Override
    public float[][] getSiteScores() {
        if (mySiteScore == null) {
            throw new IllegalStateException("CoreAlignment: getSiteScores: This Alignment has no Site Scores.");
        }
        return mySiteScore.getSiteScores();
    }

    @Override
    public boolean hasSiteScores() {
        if (mySiteScore == null) {
            return false;
        } else {
            return true;
        }
    }

    @Override
    public SITE_SCORE_TYPE getSiteScoreType() {
        return mySiteScore.getSiteScoreType();
    }

    @Override
    public int getIndelSize(int site) {
        return myPositionList.getIndelSize(site);
    }

    @Override
    public boolean isIndel(int site) {
        return myPositionList.isIndel(site);
    }

    @Override
    public boolean isAllPolymorphic() {
        return myGenotype.isAllPolymorphic();
    }

    @Override
    public boolean isPolymorphic(int site) {
        return myGenotype.isPolymorphic(site);
    }

    @Override
    public byte getMajorAllele(int site) {
        return myGenotype.getMajorAllele(site);
    }

    @Override
    public String getMajorAlleleAsString(int site) {
        return myGenotype.getMajorAlleleAsString(site);
    }

    @Override
    public byte getMinorAllele(int site) {
        return myGenotype.getMinorAllele(site);
    }

    @Override
    public String getMinorAlleleAsString(int site) {
        return myGenotype.getMinorAlleleAsString(site);
    }

    @Override
    public byte[] getMinorAlleles(int site) {
        return myGenotype.getMinorAlleles(site);
    }

    @Override
    public byte[] getAlleles(int site) {
        return myGenotype.getAlleles(site);
    }

    @Override
    public double getMinorAlleleFrequency(int site) {
        return myGenotype.getMinorAlleleFrequency(site);
    }

    @Override
    public double getMajorAlleleFrequency(int site) {
        return myGenotype.getMajorAlleleFrequency(site);
    }

    @Override
    public TaxaList getTaxaList() {
        return myTaxaList;
    }

    @Override
    public String getTaxaName(int index) {
        return myTaxaList.getTaxaName(index);
    }

    @Override
    public String getFullTaxaName(int index) {
        return myTaxaList.getFullTaxaName(index);
    }

    @Override
    public String getGenomeAssembly() {
        return myPositionList.getGenomeAssembly();
    }

    @Override
    public boolean isPositiveStrand(int site) {
        return myPositionList.isPositiveStrand(site);
    }

    @Override
    public AlignmentNew[] getAlignments() {
        return new AlignmentNew[]{this};
    }

    @Override
    public int[][] getAllelesSortedByFrequency(int site) {
        return myGenotype.getAllelesSortedByFrequency(site);
    }

    @Override
    public Object[][] getDiploidsSortedByFrequency(int site) {
        return myGenotype.getDiploidsSortedByFrequency(site);
    }

    @Override
    public boolean isPhased() {
        return myGenotype.isPhased();
    }

    @Override
    public GeneticMap getGeneticMap() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean retainsRareAlleles() {
        return myGenotype.retainsRareAlleles();
    }

    @Override
    public String[][] getAlleleEncodings() {
        return myGenotype.getAlleleEncodings();
    }

    @Override
    public String[] getAlleleEncodings(int site) {
        return myGenotype.getAlleleEncodings(site);
    }

    @Override
    public String getBaseAsString(int site, byte value) {
        return myGenotype.getBaseAsString(site, value);
    }

    @Override
    public String getDiploidAsString(int site, byte value) {
        return myGenotype.getDiploidAsString(site, value);
    }

    @Override
    public int getMaxNumAlleles() {
        return myGenotype.getMaxNumAlleles();
    }

    @Override
    public int getTotalNumAlleles() {
        return myGenotype.getTotalNumAlleles();
    }

    @Override
    public int getTotalGametesNotMissing(int site) {
        return myGenotype.getTotalGametesNotMissing(site);
    }

    @Override
    public int getTotalNotMissing(int site) {
        return myGenotype.getTotalNotMissing(site);
    }

    @Override
    public int getMinorAlleleCount(int site) {
        return myGenotype.getMinorAlleleCount(site);
    }

    @Override
    public int getMajorAlleleCount(int site) {
        return myGenotype.getMajorAlleleCount(site);
    }

    @Override
    public Object[][] getDiploidCounts() {
        return myGenotype.getDiploidCounts();
    }

    @Override
    public Object[][] getMajorMinorCounts() {
        return myGenotype.getMajorMinorCounts();
    }

    @Override
    public int getTotalGametesNotMissingForTaxon(int taxon) {
        return myGenotype.getTotalGametesNotMissingForTaxon(taxon);
    }

    @Override
    public int getHeterozygousCountForTaxon(int taxon) {
        return myGenotype.getHeterozygousCountForTaxon(taxon);
    }

    @Override
    public int getTotalNotMissingForTaxon(int taxon) {
        return myGenotype.getTotalNotMissingForTaxon(taxon);
    }

    @Override
    public short[] getDepthForAlleles(int taxon, int site) {
        return myAlleleDepth.getDepthForAlleles(taxon, site);
    }

    @Override
    public byte[] getAllelesByScope(ALLELE_SCOPE_TYPE scope, int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public BitSet getAllelePresenceForAllTaxaByScope(ALLELE_SCOPE_TYPE scope, int site, int alleleNumber) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public BitStorage getBitStorage(ALLELE_SCOPE_TYPE scopeType) {
        switch (scopeType) {
            case Frequency:
                if (myFreqBitStorage == null) {
                    myFreqBitStorage = new DynamicBitStorage(myGenotype, scopeType, myGenotype.getMajorAlleleForAllSites(), myGenotype.getMinorAlleleForAllSites());
                }
                return myFreqBitStorage;
            case Reference:
                if (myReferenceBitStorage == null) {
                    myReferenceBitStorage = DynamicBitStorage.getInstance(myGenotype, scopeType, getReference());
                }
                return myReferenceBitStorage;
            default:
                myLogger.warn("getBitStorage: Unsupported type: " + scopeType);
                return null;
        }
    }
}
