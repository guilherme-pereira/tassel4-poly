/*
 * MutableNucleotideAlignment
 */
package net.maizegenetics.pal.alignment;

import java.util.List;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;

/**
 *
 * @author terry
 */
public class MutableNucleotideAlignment extends MutableSingleEncodeAlignment implements MutableAlignment {

    protected MutableNucleotideAlignment(Alignment a, int maxNumTaxa, int maxNumSites) {
        super(a, maxNumTaxa, maxNumSites);
    }

    public static MutableNucleotideAlignment getInstance(Alignment a, int maxTaxa, int maxNumSites) {
        if (a.getAlleleEncodings() == NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES) {
            return new MutableNucleotideAlignment(a, maxTaxa, maxNumSites);
        } else {
            throw new IllegalArgumentException("MutableNucleotideAlignment: getInstance: alignment must be nucleotide data.");
        }
    }

    public static MutableNucleotideAlignment getInstance(Alignment a) {
        return getInstance(a, a.getSequenceCount(), a.getSiteCount());
    }

    protected MutableNucleotideAlignment(IdGroup idGroup, int initNumSites, int maxNumTaxa, int maxNumSites) {
        super(NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES, idGroup, initNumSites, maxNumTaxa, maxNumSites);
    }

    protected MutableNucleotideAlignment(List<Identifier> idGroup, int[] variableSites, List<Locus> locusToLociIndex, int[] locusIndices, String[] siteNames) {
        super(NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES, idGroup, variableSites, locusToLociIndex, locusIndices, siteNames);
    }

    public static MutableNucleotideAlignment getInstance(IdGroup idGroup, int initNumSites, int maxNumTaxa, int maxNumSites) {
        return new MutableNucleotideAlignment(idGroup, initNumSites, maxNumTaxa, maxNumSites);
    }

    public static MutableNucleotideAlignment getInstance(IdGroup idGroup, int maxNumSites) {
        return new MutableNucleotideAlignment(idGroup, 0, idGroup.getIdCount(), maxNumSites);
    }

    public static MutableSingleEncodeAlignment getInstance(List<Identifier> idGroup, int[] variableSites, List<Locus> locusToLociIndex, int[] locusIndices, String[] siteNames) {
        return new MutableNucleotideAlignment(idGroup, variableSites, locusToLociIndex, locusIndices, siteNames);
    }

    @Override
    public String getBaseAsString(int taxon, int site) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(getBase(taxon, site));
    }

    @Override
    public String getDiploidAsString(int site, byte value) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(value);
    }
}
