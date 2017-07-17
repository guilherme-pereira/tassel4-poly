/*
 * Alignment
 */
package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.ProgressListener;

import java.io.Serializable;

/**
 * This supports heterozygous diploid alignments.
 *
 * @author terry, Gabriel Rodrigues Alves Margarido
 */
public interface Alignment extends Serializable {

    /**
     * This encoding is used to lump together allele values with frequencies too
     * low to be retained as one of the maximum number of alleles.
     */
    public static byte RARE_ALLELE = 0xE;
    public static byte RARE_DIPLOID_ALLELE = (byte) 0xEE;
    public static String RARE_ALLELE_STR = "Z";
    public static byte UNKNOWN_ALLELE = 0xF;
    public static byte UNKNOWN_DIPLOID_ALLELE = (byte) 0xFF;
    public static String UNKNOWN_ALLELE_STR = "N";
    public static String UNKNOWN_DIPLOID_ALLELE_STR = "N:N";
    public static char UNKNOWN_ALLELE_CHAR = 'N';
    /**
     * Default number of alleles to retain.
     */
    public static final int DEFAULT_MAX_NUM_ALLELES = 2;

    public static enum SITE_SCORE_TYPE {

        None, MixedScoreTypes, QualityScore, ImputedProbablity, Dosage
    };
    final static public byte POSITION_TYPE_ALL_GROUP = 0;
    // final static public byte POSITION_TYPE_SILENT_GROUP = 1;
    // final static public byte POSITION_TYPE_SYNONYMOUS_GROUP = 2;
    // final static public byte POSITION_TYPE_NONCODING_GROUP = 3;
    // final static public byte POSITION_TYPE_NONTRANSSCRIBED_GROUP = 4;
    // final static public byte POSITION_TYPE_INTRON_GROUP = 5;
    // final static public byte POSITION_TYPE_INDEL_GROUP = 6;
    // final static public byte POSITION_TYPE_NONCODINGINDEL_GROUP = 7;
    // final static public byte POSITION_TYPE_NONSYNONYMOUS_GROUP = 8;
    // final static public byte POSITION_TYPE_CODING_GROUP = 9;
    // final static public byte POSITION_TYPE_CODINGINDEL_GROUP = 10;
    // final static public byte POSITION_TYPE_TRANSCRIBED_GROUP = 11;
    final static public String[] POSITION_TYPE_GROUP_TEXT = {"All", "Silent", "Synonymous", "Noncoding",
        "Nontranscribed", "Intron", "Indel", "Noncoding Indel", "Nonsynonymous", "Coding",
        "Coding Indel", "Transcribed"};
    // final static public byte POSITION_TYPE_NONTRANSCRIBED_TYPE = 'N';
    // final static public byte POSITION_TYPE_ANON_CODING_TYPE = 'C';
    // final static public byte POSITION_TYPE_CODON1_TYPE = '1';
    // final static public byte POSITION_TYPE_CODON2_TYPE = '2';
    // final static public byte POSITION_TYPE_CODON3_TYPE = '3';
    // final static public byte POSITION_TYPE_INTRON_TYPE = 'I';

    /**
     * This defines the possible allele scope types.
     */
    public static enum ALLELE_SCOPE_TYPE {
        /**
         * This is the default where alleles are sorted by frequency. Same as
         * getAlleles().
         */
        Frequency,
        /**
         * This sorts alleles based on there depth value.
         */
        Depth,
        /**
         * This uses the allele frequency of a base/global Alignment determine
         * sort order of alleles. That Alignment is usually a superset.
         */
        Global_Frequency,
        /**
         * This sorts alleles based on the reference sequence.
         */
        Reference
    };

    /**
     * Returns diploid value (genotype) for a given taxon and site.
     *
     * @param taxon taxon
     * @param site site
     *
     * @return high four bits generally encode the more frequent allele and the
     * lower four bits encode the less frequent allele.
     */
    public byte getBase(int taxon, int site);

    /**
     * Returns diploid values for given taxon and site. Same values as
     * getBase(), except two values are already separated into two bytes.
     *
     * @param taxon taxon
     * @param site site
     *
     * @return first byte (index 0) holds first allele value in right-most four
     * bits. second byte (index 1) holds second allele value in right-most four
     * bits.
     */
    public byte[] getBaseArray(int taxon, int site);

    /**
     * Returns diploid values for given taxon, locus, and physical position. The
     * locus and physical position should map to an unique site.
     *
     * @param taxon taxon
     * @param locus locus
     * @param physicalPosition physical position
     *
     * @return first four bits are the first allele value and the second four
     * bits are the second allele value.
     */
    public byte getBase(int taxon, Locus locus, int physicalPosition);

    /**
     * Returns sequence of diploid allele values for given taxon in specified
     * range (end site excluded). Each value in array is what would be returned
     * by getBase().
     *
     * @param taxon taxon
     * @param startSite start site
     * @param endSite end site
     *
     * @return sequence of diploid allele values.
     */
    public byte[] getBaseRange(int taxon, int startSite, int endSite);

    /**
     * Returns sequence of diploid allele values for all sites for given taxon.
     * Each value in array is what would be returned by getBase().
     *
     * @param taxon taxon
     *
     * @return sequence of diploid allele values.
     */
    public byte[] getBaseRow(int taxon);

    /**
     * Returns sequence of true/false values indicating whether taxon at each
     * site matches a specific allele (based on frequency). Allele number of
     * value 0 would be the major allele. Allele number of value 1 would be the
     * minor allele. Allele number of value 2 would be the third most frequent
     * allele value and so on.
     *
     * @param taxon taxon
     * @param alleleNumber allele number
     *
     * @return sequence of true/false values.
     */
    public BitSet getAllelePresenceForAllSites(int taxon, int alleleNumber);

    /**
     * Returns sequence of true/false values indicating whether site at each
     * taxon matches a specific allele (based on frequency). Allele number of
     * value 0 would be the major allele. Allele number of value 1 would be the
     * minor allele. Allele number of value 2 would be the third most frequent
     * allele value and so on.
     *
     * @param site site
     * @param alleleNumber allele number
     *
     * @return sequence of true/false values.
     */
    public BitSet getAllelePresenceForAllTaxa(int site, int alleleNumber);

    /**
     * Returns sequence of true/false values indicating whether taxon at sites
     * (in given blocks, 64 sites per block including start block but excluding
     * end block) matches a specific allele (based on frequency). Allele number
     * of value 0 would be the major allele. Allele number of value 1 would be
     * the minor allele. Allele number of value 2 would be the third most
     * frequent allele value and so on.
     *
     * @param taxon taxon
     * @param alleleNumber allele number
     * @param startBlock starting block
     * @param endBlock end block
     *
     * @return sequence of true/false values.
     */
    public long[] getAllelePresenceForSitesBlock(int taxon, int alleleNumber, int startBlock, int endBlock);

    /**
     * Returns sequence of true/false values indicating whether taxon at each
     * site for given parent matches a specific allele (based on frequency).
     * Allele number of value 0 would be the major allele. Allele number of
     * value 1 would be the minor allele. Allele number of value 2 would be the
     * third most frequent allele value and so on.
     *
     * @param taxon taxon
     * @param firstParent true for first parent (false for second parent)
     * @param alleleNumber allele number
     *
     * @return sequence of true/false values.
     */
    public BitSet getPhasedAllelePresenceForAllSites(int taxon, boolean firstParent, int alleleNumber);

    /**
     * Returns sequence of true/false values indicating whether site at each
     * taxon for given parent matches a specific allele (based on frequency).
     * Allele number of value 0 would be the major allele. Allele number of
     * value 1 would be the minor allele. Allele number of value 2 would be the
     * third most frequent allele value and so on.
     *
     * @param site site
     * @param firstParent true for first parent (false for second parent)
     * @param alleleNumber allele number
     *
     * @return sequence of true/false values.
     */
    public BitSet getPhasedAllelePresenceForAllTaxa(int site, boolean firstParent, int alleleNumber);

    /**
     * Returns sequence of true/false values indicating whether taxon at sites
     * (in given blocks, 64 sites per block including start block but excluding
     * end block) for given parent matches a specific allele (based on
     * frequency). Allele number of value 0 would be the major allele. Allele
     * number of value 1 would be the minor allele. Allele number of value 2
     * would be the third most frequent allele value and so on.
     *
     * @param taxon taxon
     * @param firstParent true for first parent (false for second parent)
     * @param alleleNumber allele number
     * @param startBlock starting block
     * @param endBlock end block
     *
     * @return sequence of true/false values.
     */
    public long[] getPhasedAllelePresenceForSitesBlock(int taxon, boolean firstParent, int alleleNumber, int startBlock, int endBlock);

    /**
     * Returns string representation of diploid values returned by getBase() for
     * given taxon and site. The two allele values will be separated by a colon
     * (:) delimiter.
     *
     * @param taxon taxon
     * @param site site
     *
     * @return string representation of diploid values.
     */
    public String getBaseAsString(int taxon, int site);

    /**
     * Returns string representation of diploid alleles for given taxon in
     * specified range (end site excluded). Each value in string is what would
     * be returned by getBaseAsString().
     *
     * @param taxon taxon
     * @param startSite start site
     * @param endSite end site
     *
     * @return string representation of alleles in range
     */
    public String getBaseAsStringRange(int taxon, int startSite, int endSite);

    /**
     * Returns string representation of diploid alleles for given taxon for all
     * sites. Each value in string is what would be returned by
     * getBaseAsString().
     *
     * @param taxon taxon
     *
     * @return string representation of alleles
     */
    public String getBaseAsStringRow(int taxon);

    /**
     * Returns string representation of diploid values returned by
     * getBaseArray() for given taxon and site. Same two allele values as
     * getBaseAsString(), except already separated into two Strings.
     *
     * @param taxon taxon
     * @param site site
     *
     * @return string representations of diploid values.
     */
    public String[] getBaseAsStringArray(int taxon, int site);

    /**
     * Return reference diploid allele values at given site.
     *
     * @param site site
     *
     * @return first four bits are the first allele value and the second four
     * bits are the second allele value.
     */
    public byte getReferenceAllele(int site);

    /**
     * Returns reference sequence of diploid allele values for given taxon in
     * specified range (end site not included). Each value in array contains
     * both diploid values. First four bits holds the first allele, and the
     * second four bits holds the second allele.
     *
     * @param startSite start site
     * @param endSite end site
     *
     * @return reference sequence of diploid allele values.
     */
    public byte[] getReference(int startSite, int endSite);

    /**
     * Returns reference sequence of diploid allele values. Each value in array
     * contains both diploid values. First four bits holds the first allele, and
     * the second four bits holds the second allele.
     *
     * @return reference sequence of diploid allele values.
     */
    public byte[] getReference();

    /**
     * Return whether this alignment has defined reference sequence.
     *
     * @return true if this alignment has reference sequence.
     */
    public boolean hasReference();

    /**
     * Returns whether allele values at given taxon and site are heterozygous.
     * If two values returned by getBase() are different, this will return
     * false.
     *
     * @param taxon taxon
     * @param site site
     *
     * @return whether heterozygous
     */
    public boolean isHeterozygous(int taxon, int site);

    /**
     * Returns number of heterozygous taxa at given site.
     *
     * @param site site
     *
     * @return number of heterozygous taxa
     */
    public int getHeterozygousCount(int site);

    /**
     * Get SNP IDs.
     *
     * @return site names.
     */
    public String[] getSNPIDs();

    /**
     * Get SNP ID for specified site.
     *
     * @param site site
     * @return site name
     */
    public String getSNPID(int site);

    /**
     * Returns total number of sites of this alignment.
     *
     * @return number of sites
     */
    public int getSiteCount();

    /**
     * Return number of sites for given locus.
     *
     * @param locus locus
     *
     * @return number of sites
     */
    public int getLocusSiteCount(Locus locus);

    /**
     * Get the first (inclusive) and last (exclusive) site of the specified
     * locus in this alignment.
     *
     * @param locus locus
     *
     * @return first and last site
     */
    public int[] getStartAndEndOfLocus(Locus locus);

    /**
     * Returns number of sequences (taxa).
     *
     * @return number of sequences
     */
    public int getSequenceCount();
    
    /**
     * Returns number of taxa (same as getSequenceCount()
     * 
     * @return number of taxa
     */
    public int getTaxaCount();

    /**
     * Returns the physical position at given site.
     *
     * @param site site
     *
     * @return physical position
     */
    public int getPositionInLocus(int site);

    /**
     * Return site of given physical position in locus. If the physical position
     * doesn't exist, (-(insertion point) - 1) is returned. If locus is not
     * found, an exception is thrown.
     *
     * @param physicalPosition physical position
     * @param locus locus. if null, the first locus is used.
     *
     * @return index
     */
    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus);

    /**
     * Return site of given physical position / SNP ID in locus. If the physical
     * position doesn't exist, (-(insertion point) - 1) is returned. If locus is
     * not found, an exception is thrown. This is to support multiple sites with
     * the same physical position but different SNP IDs.
     *
     * @param physicalPosition physical position
     * @param locus locus. if null, the first locus is used.
     * @param snpID SNP ID
     *
     * @return index
     */
    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus, String snpID);

    /**
     * Returns all physical positions.
     *
     * @return physical positions.
     */
    public int[] getPhysicalPositions();

    /**
     * Return Locus Name for given site.
     *
     * @param site site
     *
     * @return Locus Name
     */
    public String getLocusName(int site);

    /**
     * Return Locus for given site.
     *
     * @param site site
     *
     * @return Locus
     */
    public Locus getLocus(int site);

    /**
     * Return Locus with matching name. First to match will be returned.
     *
     * @param name name
     *
     * @return Locus
     */
    public Locus getLocus(String name);

    /**
     * Return all loci.
     *
     * @return loci
     */
    public Locus[] getLoci();

    /**
     * Return number of loci.
     *
     * @return number of loci
     */
    public int getNumLoci();

    /**
     * Returns starting site for each locus.
     *
     * @return starting site for each locus.
     */
    public int[] getLociOffsets();

    /**
     * Returns the site score of the given sequence and site.
     *
     * @param seq sequence index
     * @param site site
     *
     * @return site score.
     */
    public float getSiteScore(int seq, int site);

    /**
     * Returns the site scores.
     *
     * @return site scores.
     */
    public float[][] getSiteScores();

    /**
     * Returns true if this alignment has site scores.
     *
     * @return true if this alignment has site scores.
     */
    public boolean hasSiteScores();

    /**
     * Return what type of site scores this alignment has.
     *
     * @return site score type.
     */
    public SITE_SCORE_TYPE getSiteScoreType();

    /**
     * Return size of indel at given site.
     *
     * @param site site
     *
     * @return indel size
     */
    public int getIndelSize(int site);

    /**
     * Returns whether give site is an indel.
     *
     * @param site site
     *
     * @return true if indel
     */
    public boolean isIndel(int site);

    /**
     * Returns whether all sites are polymorphic.
     *
     * @return true if all sites are polymorphic.
     */
    public boolean isAllPolymorphic();

    /**
     * Return whether given site is polymorphic.
     *
     * @param site site
     *
     * @return true if given site is polymorphic.
     */
    public boolean isPolymorphic(int site);

    /**
     * Return most common allele at given site. Gap is included as state.
     * Heterozygous count one for each allele value. Homozygous counts two for
     * the allele value.
     *
     * @param site site
     *
     * @return most common allele
     */
    public byte getMajorAllele(int site);

    /**
     * Return most common allele at given site. Gap is included as state.
     * Heterozygous count one for each allele value. Homozygous counts two for
     * the allele value.
     *
     * @param site site
     *
     * @return most common allele as String
     */
    public String getMajorAlleleAsString(int site);

    /**
     * Return most common minor allele at given site. Gap is included as state.
     * Heterozygous count one for each allele value. Homozygous counts two for
     * the allele value.
     *
     * @param site site
     *
     * @return most common minor allele
     */
    public byte getMinorAllele(int site);

    /**
     * Return most common minor allele at given site. Gap is included as state.
     * Heterozygous count one for each allele value. Homozygous counts two for
     * the allele value.
     *
     * @param site site
     *
     * @return most common minor allele as String
     */
    public String getMinorAlleleAsString(int site);

    /**
     * Return all minor alleles at given site. Gap is included as state.
     * Heterozygous count one for each allele value. Homozygous counts two for
     * the allele value.
     *
     * @param site site
     *
     * @return all minor alleles
     */
    public byte[] getMinorAlleles(int site);

    /**
     * Returns all alleles at given site in order of frequency. Gap is included
     * as state. Heterozygous count one for each allele value. Homozygous counts
     * two for the allele value.
     *
     * @param site site
     *
     * @return all alleles
     */
    public byte[] getAlleles(int site);

    /**
     * Return frequency for most common minor allele at given site. Gap is
     * included as state. Heterozygous count one for each allele value.
     * Homozygous counts two for the allele value.
     *
     * @param site site
     *
     * @return frequency
     */
    public double getMinorAlleleFrequency(int site);

    /**
     * Return frequency for major allele at given site. Gap is included as
     * state. Heterozygous count one for each allele value. Homozygous counts
     * two for the allele value.
     *
     * @param site site
     *
     * @return frequency
     */
    public double getMajorAlleleFrequency(int site);

    /**
     * Return id group of this alignment.
     *
     * @return id group.
     */
    public IdGroup getIdGroup();

    /**
     * Return taxa name at given index.
     *
     * @param index
     *
     * @return taxa name
     */
    public String getTaxaName(int index);

    /**
     * Return full taxa name at given index.
     *
     * @param index
     * @return full taxa name
     */
    public String getFullTaxaName(int index);

    /**
     * Gets the Genome Assembly.
     *
     * @return the genome assembly.
     */
    public String getGenomeAssembly();

    /**
     * Return whether is positive strand at given site.
     *
     * @param site site
     *
     * @return whether is positive strand.
     */
    public boolean isPositiveStrand(int site);

    /**
     * Returns individual alignments within this alignment.
     *
     * @return list of alignments.
     */
    public Alignment[] getAlignments();

    /**
     * Return sorted list of alleles from highest frequency to lowest at given
     * site in alignment. Resulting double dimension array holds alleles (bytes)
     * in result[0]. And the counts are in result[1]. Counts haploid values
     * twice and diploid values once. Higher ploids are not supported.
     *
     * @param site site
     *
     * @return sorted list of alleles and counts
     */
    public int[][] getAllelesSortedByFrequency(int site);

    /**
     * Return sorted list of diploid vales from highest frequency to lowest at
     * given site in alignment. Resulting double dimension array holds diploids
     * (Strings) in result[0]. And the counts are in result[1] (Integers).
     *
     * @param site site
     *
     * @return sorted list of diploids and counts
     */
    public Object[][] getDiploidssSortedByFrequency(int site);

    /**
     * Returns whether this alignment is phased.
     *
     * @return true if phased.
     */
    public boolean isPhased();

    /**
     * Returns Genetic Map of this Alignment if available.
     *
     * @return Genetic Map
     */
    public GeneticMap getGeneticMap();

    /**
     * Returns true if this Alignment retains rare alleles. If false, rare
     * alleles are recorded as unknown.
     *
     * @return whether rare alleles are retained.
     */
    public boolean retainsRareAlleles();

    /**
     * Returns allele values as strings for all sites. The first dimension of
     * the array indexes the sites. The second dimension indexes the allele
     * values for given site. The indices for the allele values are used as the
     * codes to store data. These codes (indices) are returned by the getBase()
     * methods. If only one array of allele values is returned, that is the
     * encoding for all sites.
     *
     * @return allele values for all sites.
     */
    public String[][] getAlleleEncodings();

    /**
     * Same as getAlleleEncodings() for only one site.
     *
     * @param site site
     *
     * @return allele values for given site.
     */
    public String[] getAlleleEncodings(int site);

    /**
     * Returns String representation of allele value at site.
     *
     * @param site site
     * @param value allele value
     *
     * @return String representation
     */
    public String getBaseAsString(int site, byte value);

    /**
     * Returns String representation of diploid allele value at site.
     *
     * @param site site
     * @param value diploid allele value
     *
     * @return String representation
     */
    public String getDiploidAsString(int site, byte value);

    /**
     * Return max number of alleles retained by this alignment.
     *
     * @return max number of alleles.
     */
    public int getMaxNumAlleles();

    /**
     * Returns max number of alleles retained by this alignment including rare
     * (if retained) and unknown (if phased).
     *
     * @return total number alleles
     */
    public int getTotalNumAlleles();

    /**
     * Returns total number of non-missing allele values for given site. This
     * can be twice the number of taxa, as diploid values are supported.
     *
     * @param site site
     * @return number of non-missing allele values.
     */
    public int getTotalGametesNotMissing(int site);

    /**
     * Returns total number of non-missing taxa for given site. Taxa are
     * considered missing only if both allele values are Unknown (N).
     *
     * @param site site
     *
     * @return number of non-missing taxa..
     */
    public int getTotalNotMissing(int site);

    /**
     * Returns the minor allele count for given site.
     *
     * @param site site
     * @return minor allele count
     */
    public int getMinorAlleleCount(int site);

    /**
     * Returns the major allele count for given site.
     *
     * @param site site
     * @return major allele count
     */
    public int getMajorAlleleCount(int site);

    /**
     * Returns counts of all diploid combinations from highest frequency to
     * lowest for whole alignment. Resulting double dimension array holds
     * diploids (Strings) in result[0]. And the counts are in result[1] (Longs).
     *
     * @return diploid counts.
     */
    public Object[][] getDiploidCounts();

    /**
     * Returns counts of all major/minor allele combinations from highest
     * frequency to lowest for whole alignment. Resulting double dimension array
     * holds major/minor allele (Strings) in result[0]. And the counts are in
     * result[1] (Longs).
     *
     * @return diploid counts.
     */
    public Object[][] getMajorMinorCounts();

    /**
     * Returns total number of non-missing allele values for given taxon. This
     * can be twice the number of sites, as diploid values are supported.
     *
     * @param taxon taxon
     *
     * @return number of non-missing allele values.
     */
    public int getTotalGametesNotMissingForTaxon(int taxon);

    /**
     * Returns number of heterozygous sites at given taxon.
     *
     * @param taxon taxon
     *
     * @return number of heterozygous sites
     */
    public int getHeterozygousCountForTaxon(int taxon);

    /**
     * Returns total number of non-missing sites for given taxon. Sites are
     * considered missing only if both allele values are Unknown (N).
     *
     * @param taxon taxon
     *
     * @return number of non-missing sites.
     */
    public int getTotalNotMissingForTaxon(int taxon);

    /**
     * Return whether alignment will execute quickly for site optimized
     * operations. SBitAlignment is obviously friendly. But so would
     * FilterAlignment is only sites have been filtered for example.
     *
     * @return whether optimized for site operations.
     */
    public boolean isSBitFriendly();

    /**
     * Return whether alignment will execute quickly for taxa optimized
     * operations. TBitAlignment is obviously friendly. But so would
     * FilterAlignment is only taxa have been filtered for example.
     *
     * @return whether optimized for taxa operations.
     */
    public boolean isTBitFriendly();

    /**
     * Optimizes this Alignment for Taxa based operations.
     */
    public void optimizeForTaxa(ProgressListener listener);

    /**
     * Optimizes this Alignment for Site based operations.
     */
    public void optimizeForSites(ProgressListener listener);

    /**
     * Returns depth count for each diploid allele at the given taxon and site.
     *
     * @param taxon taxon
     * @param site site
     *
     * @return two counts
     */
    public short[] getDepthForAlleles(int taxon, int site);

    /**
     * Returns all alleles at given site in order defined by scope.
     *
     * @param scope scope
     * @param site site
     *
     * @return alleles
     */
    public byte[] getAllelesByScope(ALLELE_SCOPE_TYPE scope, int site);

    /**
     * Returns sequence of true/false values indicating whether site at each
     * taxon matches a specific allele (based on scope).
     *
     * @param scope scope
     * @param site site
     * @param alleleNumber allele number
     *
     * @return sequence of true/false values.
     */
    public BitSet getAllelePresenceForAllTaxaByScope(ALLELE_SCOPE_TYPE scope, int site, int alleleNumber);
}
