/*
 * BitTextAlignment
 */
package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.ids.IdGroup;

/**
 *
 * @author terry
 */
public class BitTextAlignment extends BitAlignment {
    
    private static final long serialVersionUID = -5197800047652332969L;

    protected BitTextAlignment(Alignment a, int maxNumAlleles, boolean retainRareAlleles, boolean isSBit) {
        super(a, maxNumAlleles, retainRareAlleles, isSBit);
    }

    protected BitTextAlignment(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isSBit) {
        super(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isSBit);
    }

    @Override
    public String getBaseAsString(int taxon, int site) {
        byte[] temp = getBaseArray(taxon, site);
        return myAlleleStates[site][temp[0]] + ":" + myAlleleStates[site][temp[1]];
    }

    @Override
    public String[] getBaseAsStringArray(int taxon, int site) {
        byte[] temp = getBaseArray(taxon, site);
        return new String[]{myAlleleStates[site][temp[0]], myAlleleStates[site][temp[1]]};
    }

    @Override
    public String getBaseAsString(int site, byte value) {
        return myAlleleStates[site][value];
    }

    @Override
    public String getDiploidAsString(int site, byte value) {
        return myAlleleStates[site][(value >>> 4) & 0xf] + ":" + myAlleleStates[site][value & 0xf];
    }
}
