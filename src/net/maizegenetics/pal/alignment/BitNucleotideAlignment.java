/*
 * BitNucleotideAlignment
 */
package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.util.BitSet;

/**
 *
 * @author terry
 */
public class BitNucleotideAlignment extends BitAlignment {

    private static final long serialVersionUID = -5197800047652332969L;

    protected BitNucleotideAlignment(Alignment a, int maxNumAlleles, boolean retainRareAlleles, boolean isSBit) {
        super(a, maxNumAlleles, retainRareAlleles, isSBit);
    }

    protected BitNucleotideAlignment(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isSBit) {
        super(idGroup, data, map, reference, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isSBit);
    }

    protected BitNucleotideAlignment(IdGroup idGroup, byte[][] alleles, BitSet[][] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isSBit) {
        super(idGroup, alleles, data, map, reference, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isSBit);
    }

    @Override
    public String getBaseAsStringRange(int taxon, int startSite, int endSite) {
        StringBuilder builder = new StringBuilder();
        for (int i = startSite; i < endSite; i++) {
            builder.append(getBaseAsString(taxon, i));
        }
        return builder.toString();
    }
    
    @Override
    public String getBaseAsString(int taxon, int site) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(getBase(taxon, site));
    }

    @Override
    public String getDiploidAsString(int site, byte value) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(value);
    }

    @Override
    public boolean isIndel(int site) {
        int[][] alleles = getAllelesSortedByFrequency(site);
        int numAlleles = Math.min(alleles[0].length, 2);
        for (int i = 0; i < numAlleles; i++) {
            if ((alleles[0][i] == NucleotideAlignmentConstants.INSERT_ALLELE) || (alleles[0][i] == NucleotideAlignmentConstants.GAP_ALLELE)) {
                return true;
            }
        }
        return false;
    }
}
