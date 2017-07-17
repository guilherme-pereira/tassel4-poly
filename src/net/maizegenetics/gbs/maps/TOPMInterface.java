/*
 * TOPMInterface
 */
package net.maizegenetics.gbs.maps;

import java.io.File;
import net.maizegenetics.gbs.tagdist.Tags;
import net.maizegenetics.pal.alignment.Locus;

/**
 * Tags on Physical Map (TOPM) Interface. Define methods for the TOPM classes.
 * The TOPM class holds information relating tags to the physical map, variants
 * (SNPs) called by the tag, and quality of mapping information.
 * <p>
 * Physical map information is recorded with chromosome, strand, start and end positions.  If a tags maps to multiple
 * locations with equal quality it is recorded in multimaps with the count of multi-mapping.
 * <p>
 * Each SNP call for a tag is record in the variant fields.  The position by the variant offset,
 * which is added to the startPosition (negative values indicate before the startPosition).  Each variant is also
 * defined by the SNP call implied by the given tag.  The called are stored in a byte using the
 * NucleotideAlignmentConstants (A=00, C=01, G=02, T=03).
 *
 * @author Ed Buckler and Terry Casstevens
 */
public interface TOPMInterface extends Tags {

    public final static byte BYTE_MISSING = Byte.MIN_VALUE;
    public final static int INT_MISSING = Integer.MIN_VALUE;

    /**
     * Adds Variant definition.
     *
     * @param tagIndex tag index
     * @param offset offset
     * @param base base value one of: NucleotideAlignmentConstants.A_ALLELE,
     * NucleotideAlignmentConstants.C_ALLELE,
     * NucleotideAlignmentConstants.G_ALLELE,
     * NucleotideAlignmentConstants.T_ALLELE,
     * NucleotideAlignmentConstants.GAP_ALLELE,
     * NucleotideAlignmentConstants.INSERT_ALLELE, Alignment.UNKNOWN_ALLELE
     *
     * @return variant index or -1 if max number variants reached at this tag.
     */
    public int addVariant(int tagIndex, byte offset, byte base);

    /**
     * Compares tags at given indices.
     *
     * @param index1 first index
     * @param index2 second index
     *
     * @return -1 if tag at first index less than tag at second index. 1 if tag
     * at first index greater. 0 if equal.
     */
    public int compare(int index1, int index2);

    /**
     * Returns chromosome at give tag index.
     *
     * @param index tag index
     *
     * @return chromosome
     */
    public int getChromosome(int index);

    /**
     * Returns an array whose <i>values</i> are the distinct chromosomes in this
     * file, as stored in the chromosome[] array. The indices are arbitrary.
     *
     * @return list of distinct chromosomes
     */
    public int[] getChromosomes();

    /**
     * Round(Log2(P))
     *
     * @param index tag index
     *
     * @return DcoP
     */
    public byte getDcoP(int index);

    /**
     * Returns Divergence
     *
     * @param index tag index
     *
     * @return Divergence
     */
    public byte getDivergence(int index);

    /**
     * Returns End Position of Tag
     *
     * @param index tag index
     *
     * @return end position of tag
     */
    public int getEndPosition(int index);

    /**
     * Returns Loci created from getChromosomes().
     *
     * @return Loci
     */
    public Locus[] getLoci();

    /**
     * Get Locus representing chromosome for given tag.
     *
     * @param tagIndex tag index
     *
     * @return Locus
     */
    public Locus getLocus(int tagIndex);

    public byte getMapP(int index);

    /**
     * Returns maximum number of variants stored per tag.
     *
     * @return maximum number of variants
     */
    public int getMaxNumVariants();

    public byte getMultiMaps(int index);

    /**
     * Returns chromosome, strand, and start position for given tag.
     *
     * @param index tag index
     *
     * @return index 0 is chromosome, 1 is strand, and 2 is start position
     */
    public int[] getPositionArray(int index);

    /**
     * Consider removing this. -Terry
     */
    public int getReadIndexForPositionIndex(int posIndex);

    /**
     * Returns number of tags.
     *
     * @return number of tags
     */
    public int getSize();

    /**
     * Returns start position of given tag.
     *
     * @param index tag index
     *
     * @return start position
     */
    public int getStartPosition(int index);

    /**
     * Returns Strand for given tag.
     *
     * @param tagIndex tag index
     *
     * @return 1 = same sense as reference FASTA file. -1 = opposite sense.
     * unknown = Byte.MIN_VALUE
     */
    public byte getStrand(int tagIndex);

    /**
     * Returns variant definition at given tag and variant index.
     *
     * @param tagIndex tag index
     * @param variantIndex variant index
     *
     * @return variant definition (see addVariant())
     */
    public byte getVariantDef(int tagIndex, int variantIndex);

    /**
     * Returns the set of unique positions for the given chromosome
     *
     * @param chromosome
     * @return unique positions on chromosome
     */
    public int[] getUniquePositions(int chromosome);

    /**
     * Returns an array containing all variant definitions for given tag.
     *
     * @param tagIndex tag index
     *
     * @return variant definitions for tag
     */
    public byte[] getVariantDefArray(int tagIndex);

    /**
     * Returns variant position offset from start position at given tag and
     * variant index.
     *
     * @param tagIndex tag index
     * @param variantIndex variant index
     *
     * @return variant position offset
     */
    public byte getVariantPosOff(int tagIndex, int variantIndex);

    /**
     * Returns an array containing all variant position offsets for given tag.
     *
     * @param tagIndex tag index
     *
     * @return variant position offsets for tag
     */
    public byte[] getVariantPosOffArray(int tagIndex);

    /**
     * Returns variant position offsets for all tags. First index of result is
     * tag and second is variant.
     *
     * @return all variant position offsets
     */
    public byte[][] getVariantOff();

    /**
     * Returns variant definitions for all tags. First index of result is tag
     * and second is variant.
     *
     * @return all variant definitions
     */
    public byte[][] getVariantDef();

    /**
     * Sets chromosome, strand, start position, and end position for given tag.
     *
     * @param index tag index
     * @param chromosome chromosome
     * @param strand strand
     * @param positionMin start position
     * @param positionMax end position
     */
    public void setChromoPosition(int index, int chromosome, byte strand, int positionMin, int positionMax);

    /**
     * Set Divergence for given tag
     *
     * @param index tag index
     * @param divergence divergence
     */
    public void setDivergence(int index, byte divergence);

    /**
     * Set MapP for given tag
     *
     * @param index tag index
     * @param mapP MapP
     */
    public void setMapP(int index, byte mapP);

    public void setMapP(int index, double mapP);

    /**
     * Set variant definition at given tag and variant index.
     *
     * @param tagIndex tag index
     * @param variantIndex variant index
     * @param def definition (see addVariant())
     */
    public void setVariantDef(int tagIndex, int variantIndex, byte def);

    /**
     * Set variant position offset at given tag and variant index.
     *
     * @param tagIndex tag index
     * @param variantIndex variant index
     * @param offset position offset
     */
    public void setVariantPosOff(int tagIndex, int variantIndex, byte offset);

    /**
     * Clears all variant definitions and position offsets.
     */
    public void clearVariants();

    public void writeTextFile(File outfile);
}
