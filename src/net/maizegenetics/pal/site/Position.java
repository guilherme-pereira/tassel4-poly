package net.maizegenetics.pal.site;

import net.maizegenetics.pal.util.GeneralAnnotation;

/**
 * Defines a genomic positions and its known variants.  Includes attributes of chromosome, position, strand,
 * centiMorgans, name (or SNP ID), whether this position is a nucleotide, or includes an indel.
 *
 * @author Ed Buckler
 */
public interface Position extends Comparable<Position>, GeneralAnnotation {

    /**Return the locus (generally a chromosome) of a site*/
    Chromosome getChromosome();

    /**Return the physical position of a site*/
    int getPosition();

    /**Return the strand for a site definition*/
    byte getStrand();

    /**Return the strand for a site definition*/
    float getCM();

    /**Return the ID (name) for a site*/
    String getSNPID();

    /**Whether the position is a nucleotide position or another marker type (SSR, AFLP, RAPD, CNV, which are recoded
     * with text states)*/
    boolean isNucleotide();

    /**Whether the position includes indels, which would be defined in the variants*/
    boolean isIndel();

    /**Returns the nature of the polymorphism {"ACTAT","-"} or {"A","C","G"} or {"100","103","106"}
     */
    String[] getKnownVariants();

    /**Return the minor allele frequency in a global scope*/
    float getGlobalMAF();

    /**Returns the proportion of genotypes scored at a given site*/
    float getGlobalSiteCoverage();

    /**Return the allele specified by alleleType, if unknown Alignment.Unknown is return*/
    byte getAllele(Allele alleleType);


    /**
     * Allele types recorded in an annotated position.  If unknown,
     * Alignment.UNKNOWN_ALLELE is returned.
     */
    public enum Allele {  //The indices are used in effectively as map (EnumMap is not used as it requires 4X more memory)
        /**Reference Allele*/
        REF(0),
        /**Major (most frequent) allele from the globally defined alignment*/
        GLBMAJ(1),
        /**Minor (second most frequent) allele from the globally defined alignment*/
        GLBMIN(2),
        /**Ancestral allele defined by evolutionary comparison*/
        ANC(3),
        /**High depth allele as defined from DNA sequencing analysis*/
        HIDEP(4);
        private final int index;
        /**Count of the number of allele types*/
        public final static int COUNT=Allele.values().length;
        Allele(int index) {this.index=index;}
        /**Sequential index that can be use for primitive arrays*/
        public int index() {return index;}
    }
}
