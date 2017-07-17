/*
 * BitTextAlignmentHDF5
 */
package net.maizegenetics.pal.alignment;

import ch.systemsx.cisd.hdf5.IHDF5Reader;
import net.maizegenetics.pal.ids.IdGroup;

/**
 *
 * @author terry
 */
public class BitTextAlignmentHDF5 extends BitAlignmentHDF5 {

    private static final long serialVersionUID = -5197800047652332969L;

    protected BitTextAlignmentHDF5(IHDF5Reader hdf5, IdGroup idGroup, byte[][] alleles, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {
        super(hdf5, idGroup, alleles, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);
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
