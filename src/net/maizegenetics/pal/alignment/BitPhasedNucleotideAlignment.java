/*
 * BitPhasedNucleotideAlignment
 */
package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.util.BitSet;

/**
 *
 * @author terry
 */
public class BitPhasedNucleotideAlignment extends BitPhasedAlignment {

    private static final long serialVersionUID = -5197800047652332969L;

    protected BitPhasedNucleotideAlignment(Alignment a, int maxNumAlleles, boolean retainRareAlleles, boolean isSBit) {
        super(a, maxNumAlleles, retainRareAlleles, isSBit);
    }

    protected BitPhasedNucleotideAlignment(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isSBit) {
        super(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isSBit);
    }

    protected BitPhasedNucleotideAlignment(IdGroup idGroup, byte[][] alleles, BitSet[][] data0, BitSet[][] data1, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isSBit) {
        super(idGroup, alleles, data0, data1, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isSBit);
    }

    @Override
    public String getBaseAsString(int taxon, int site) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(getBase(taxon, site));
    }

    @Override
    public String[] getBaseAsStringArray(int taxon, int site) {
        return new String[]{getBaseAsString(taxon, site)};
    }
}
