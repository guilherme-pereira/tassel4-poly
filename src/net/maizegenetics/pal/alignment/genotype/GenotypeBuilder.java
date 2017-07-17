/*
 *  GenotypeBuilder
 */
package net.maizegenetics.pal.alignment.genotype;

import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.util.SuperByteMatrix;
import net.maizegenetics.util.SuperByteMatrixBuilder;

/**
 *
 * @author Terry Casstevens
 */
public class GenotypeBuilder {

    private SuperByteMatrix myGenotype;
    private boolean myIsPhased = false;
    private String[][] myAlleleEncodings = NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES;

    private GenotypeBuilder(SuperByteMatrix genotype) {
        myGenotype = genotype;
    }

    /**
     * Get Genotype Builder given number of taxa and sites. Performance
     * optimized for site loop inside taxon loop. Default is unphased and
     * NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES encoding.
     *
     * @param numTaxa number of taxa
     * @param numSites number of sites.
     *
     * @return Genotype Builder
     */
    public static GenotypeBuilder getInstance(int numTaxa, int numSites) {
        return getUnphasedNucleotideGenotypeBuilder(numTaxa, numSites);
    }

    /**
     * Get Genotype Builder given number of taxa and sites. Performance
     * optimized for taxon loop inside site loop. Default is unphased and
     * NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES encoding.
     *
     * @param numTaxa number of taxa
     * @param numSites number of sites.
     *
     * @return Genotype Builder
     */
    public static GenotypeBuilder getInstanceTranspose(int numTaxa, int numSites) {
        return new GenotypeBuilder(SuperByteMatrixBuilder.getInstanceTranspose(numTaxa, numSites));
    }

    /**
     * Get Genotype Builder given number of taxa and sites. Performance
     * optimized for site loop inside taxon loop. Default is unphased and
     * NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES encoding.
     *
     * @param numTaxa number of taxa
     * @param numSites number of sites.
     *
     * @return Genotype Builder
     */
    public static GenotypeBuilder getUnphasedNucleotideGenotypeBuilder(int numTaxa, int numSites) {
        return new GenotypeBuilder(SuperByteMatrixBuilder.getInstance(numTaxa, numSites));
    }

    public void setBase(int taxon, int site, byte value) {
        myGenotype.set(taxon, site, value);
    }

    public void setBaseRangeForTaxon(int taxon, int startSite, byte[] value) {
        //TODO this needs an array copy method, startSite was eliminated
        for (int i = 0; i < value.length; i++) {
            myGenotype.set(taxon, i+startSite, value[i]);
        }
    }

    public void isPhased(boolean isPhased) {
        myIsPhased = isPhased;
    }

    public void alleleEncodings(String[][] alleleEncodings) {
        myAlleleEncodings = alleleEncodings;
    }

    public Genotype build() {
        SuperByteMatrix temp = myGenotype;
        myGenotype = null;
        if (NucleotideAlignmentConstants.isNucleotideEncodings(myAlleleEncodings)) {
            return new NucleotideGenotype(temp, myIsPhased);
        } else {
            return new ByteGenotype(temp, myIsPhased, myAlleleEncodings);
        }
    }

    public Genotype buildHDF5(String filename) {
        SuperByteMatrix temp = myGenotype;
        myGenotype = null;
        throw new UnsupportedOperationException();
    }
}
