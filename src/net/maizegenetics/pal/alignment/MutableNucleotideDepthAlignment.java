/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.pal.alignment;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import java.util.List;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;

/**
 *
 * @author yz79 & glaubitz, Gabriel Rodrigues Alves Margarido
 */
public class MutableNucleotideDepthAlignment extends MutableNucleotideAlignment implements MutableAlignment {

    // allelic depth information, short[allele][taxa][site] = depth, max: 32767
    private short[][][] myAlleleDepth;
    private final static int MAX_NUM_ALLELES = 6;  // will hold all 6 alleles (A,C,G,T,+,-) [= (byte) 0,1,2,3,4,5] (Jeff G)

    private MutableNucleotideDepthAlignment(Alignment a, int maxNumTaxa, int maxNumSites) {
        super(a, maxNumTaxa, maxNumSites);
        initAllelicDepthArrays(maxNumTaxa, maxNumSites);
    }

    private MutableNucleotideDepthAlignment(IdGroup idGroup, int initNumSites, int maxNumTaxa, int maxNumSites) {
        super(idGroup, initNumSites, maxNumTaxa, maxNumSites);
        initAllelicDepthArrays(maxNumTaxa, maxNumSites);
    }

    private MutableNucleotideDepthAlignment(List<Identifier> idGroup, int[] variableSites, List<Locus> locusToLociIndex, int[] locusIndices, String[] siteNames) {
        super(idGroup, variableSites, locusToLociIndex, locusIndices, siteNames);
        initAllelicDepthArrays(idGroup.size(), siteNames.length);
    }

    public static MutableNucleotideDepthAlignment getInstance(Alignment a, int maxTaxa, int maxNumSites) {
        if (a.getAlleleEncodings() == NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES) {
            return new MutableNucleotideDepthAlignment(a, maxTaxa, maxNumSites);
        } else {
            throw new IllegalArgumentException("MutableNucleotideAlignment: getInstance: alignment must be nucleotide data.");
        }
    }

    public static MutableNucleotideDepthAlignment getInstance(IdGroup idGroup, int maxNumSites) {
        return new MutableNucleotideDepthAlignment(idGroup, 0, idGroup.getIdCount(), maxNumSites);
    }

    public static MutableNucleotideDepthAlignment getInstance(IdGroup idGroup, int maxNumSites, int NumberOfAllelesToKeep) {
        throw new UnsupportedOperationException("This constructor is not supported by MutableNucleotideDepthAlignment.");
    }

    public static MutableNucleotideDepthAlignment getInstance(IdGroup idGroup, int initNumSites, int maxNumTaxa, int maxNumSites) {
        return new MutableNucleotideDepthAlignment(idGroup, initNumSites, maxNumTaxa, maxNumSites);
    }

    public static MutableNucleotideDepthAlignment getInstance(IdGroup idGroup, int initNumSites, int maxNumTaxa, int maxNumSites, int NumberOfAllelesToKeep) {
        throw new UnsupportedOperationException("This constructor is not supported by MutableNucleotideDepthAlignment.");
    }

    public static MutableNucleotideDepthAlignment getInstance(List<Identifier> idGroup, int[] variableSites, List<Locus> locusToLociIndex, int[] locusIndices, String[] siteNames) {
        return new MutableNucleotideDepthAlignment(idGroup, variableSites, locusToLociIndex, locusIndices, siteNames);
    }

    public void initAllelicDepthArrays(int numMaxTaxa, int numMaxSites) {
        myAlleleDepth = new short[MAX_NUM_ALLELES][numMaxTaxa][numMaxSites];
        for (int i = 0; i < myAlleleDepth.length; i++) {
            for (int j = 0; j < myAlleleDepth[i].length; j++) {
                for (int k = 0; k < myAlleleDepth[i][j].length; k++) {
                    myAlleleDepth[i][j][k] = (short) 0;  // initialize depth to zero (Jeff Glaubitz)
                }
            }
        }
    }

    @Override 
    public byte[] getAllelesByScope(ALLELE_SCOPE_TYPE scope, int site) {
        if (scope == ALLELE_SCOPE_TYPE.Depth) {
            throw new UnsupportedOperationException("This method is not supported by MutableNucleotideDepthAlignment.");
        } else {
            return super.getAllelesByScope(scope, site);
        }
    }

    public short[][] getDepthsForAllSites(int taxon) {
        int nSites = this.getSiteCount();
        short[][] depthsForTaxon = new short[MAX_NUM_ALLELES][nSites];
        for (int allele = 0; allele < myAlleleDepth.length; allele++) {
            for (int site = 0; site < nSites; site++) {
                depthsForTaxon[allele][site] = myAlleleDepth[allele][taxon][site];
            }
        }
        return depthsForTaxon;
    }
    
    // Return depth for all 6 alleles, even if they have never been set (Jeff G)
    @Override 
    public short[] getDepthForAlleles(int taxon, int site) {
        short[] alleleDepths = new short[myAlleleDepth.length];
        for (int i = 0; i < myAlleleDepth.length; i++) {
            alleleDepths[i] = myAlleleDepth[i][taxon][site];
        }
        return alleleDepths;
    }

    public short getDepthForAllele(int taxon, int site, int allele) {
        return myAlleleDepth[allele][taxon][site];
    }
    
    @Override
    public void setCommonAlleles(int site, byte[] values) {
        throw new UnsupportedOperationException("This method is not supported by MutableNucleotideDepthAlignment.");
    }

    @Override
    public void setDepthForAlleles(int taxon, int site, short[] values) {
        for (int i = 0; i < values.length; i++) {
            myAlleleDepth[i][taxon][site] = values[i];
        }
    }

    public void setDepthForAllele(int taxon, int site, int allele, int value) {
        if (value < 0) {
            myAlleleDepth[allele][taxon][site] = 0;
        } else if (value > 32767) {
            myAlleleDepth[allele][taxon][site] = 32767;
        } else {
            myAlleleDepth[allele][taxon][site] = (short) value;
        }
    }
    
    public void incrementDepthForAllele(int taxon, int site, int allele) {
        int newDepth = getDepthForAllele(taxon, site, allele);
        newDepth++;
        setDepthForAllele(taxon, site, allele, newDepth);
    } 
    
    @Override
    protected void sortSitesByPhysicalPosition() {

        Swapper swapperPos = new Swapper() {
            public void swap(int a, int b) {
                int it;
                it = myLocusIndices[a];
                myLocusIndices[a] = myLocusIndices[b];
                myLocusIndices[b] = it;

                byte bt;
                for (int t = 0, n = getSequenceCount(); t < n; t++) {
                    bt = getBase(t, a);
                    setBase(t, a, getBase(t, b));
                    setBase(t, b, bt);
                }

                it = myVariableSites[a];
                myVariableSites[a] = myVariableSites[b];
                myVariableSites[b] = it;

                String st = mySNPIDs[a];
                mySNPIDs[a] = mySNPIDs[b];
                mySNPIDs[b] = st;

                bt = myReference[a];
                myReference[a] = myReference[b];
                myReference[b] = bt;

                short dp;
                for (int i = 0; i < MAX_NUM_ALLELES; i++) {
                    for (int taxa = 0; taxa < getSequenceCount(); taxa++) {
                        dp = myAlleleDepth[i][taxa][a];
                        myAlleleDepth[i][taxa][a] = myAlleleDepth[i][taxa][b];
                        myAlleleDepth[i][taxa][b] = dp;
                    }
                }
            }
        };
        IntComparator compPos = new IntComparator() {
            public int compare(int a, int b) {
                if (myLocusIndices[a] < myLocusIndices[b]) {
                    return -1;
                }
                if (myLocusIndices[a] > myLocusIndices[b]) {
                    return 1;
                }
                if (myVariableSites[a] < myVariableSites[b]) {
                    return -1;
                }
                if (myVariableSites[a] > myVariableSites[b]) {
                    return 1;
                }
                return 0;
            }
        };

        GenericSorting.quickSort(0, getSiteCount(), compPos, swapperPos);

    }
}
