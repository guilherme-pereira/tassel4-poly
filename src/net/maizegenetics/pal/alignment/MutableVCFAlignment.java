/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.pal.alignment;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import java.util.ArrayList;
import java.util.List;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;

/**
 *
 * @author yz79, Gabriel Rodrigues Alves Margarido
 */
public class MutableVCFAlignment extends MutableNucleotideAlignment implements MutableAlignment {

    // allelic depth information, short[allele][taxa][site] = depth, max: 32767
    private short[][][] myAlleleDepth;
    // possible alleles for each site
    private byte[][] myCommonAlleles;
    private static int myMaxNumAlleles = 3;

    private MutableVCFAlignment(Alignment a, int maxNumTaxa, int maxNumSites) {
        super(a, maxNumTaxa, maxNumSites);
        initAllelicDepthArrays(maxNumTaxa, maxNumSites);
    }

    private MutableVCFAlignment(IdGroup idGroup, int initNumSites, int maxNumTaxa, int maxNumSites) {
        super(idGroup, initNumSites, maxNumTaxa, maxNumSites);
        initAllelicDepthArrays(maxNumTaxa, maxNumSites);
    }

    private MutableVCFAlignment(List<Identifier> idGroup, int[] variableSites, List<Locus> locusToLociIndex, int[] locusIndices, String[] siteNames) {
        super(idGroup, variableSites, locusToLociIndex, locusIndices, siteNames);
        initAllelicDepthArrays(idGroup.size(), siteNames.length);
    }

    public static MutableVCFAlignment getInstance(Alignment a, int maxTaxa, int maxNumSites) {
        if (a.getAlleleEncodings() == NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES) {
            return new MutableVCFAlignment(a, maxTaxa, maxNumSites);
        } else {
            throw new IllegalArgumentException("MutableNucleotideAlignment: getInstance: alignment must be nucleotide data.");
        }
    }

    public static MutableVCFAlignment getInstance(IdGroup idGroup, int maxNumSites) {
        return new MutableVCFAlignment(idGroup, 0, idGroup.getIdCount(), maxNumSites);
    }

    public static MutableVCFAlignment getInstance(IdGroup idGroup, int maxNumSites, int NumberOfAllelesToKeep) {
        myMaxNumAlleles = NumberOfAllelesToKeep;
        return new MutableVCFAlignment(idGroup, 0, idGroup.getIdCount(), maxNumSites);
    }

    public static MutableVCFAlignment getInstance(IdGroup idGroup, int initNumSites, int maxNumTaxa, int maxNumSites) {
        return new MutableVCFAlignment(idGroup, initNumSites, maxNumTaxa, maxNumSites);
    }

    public static MutableVCFAlignment getInstance(IdGroup idGroup, int initNumSites, int maxNumTaxa, int maxNumSites, int NumberOfAllelesToKeep) {
        myMaxNumAlleles = NumberOfAllelesToKeep;
        return new MutableVCFAlignment(idGroup, initNumSites, maxNumTaxa, maxNumSites);
    }

    public static MutableVCFAlignment getInstance(List<Identifier> idGroup, int[] variableSites, List<Locus> locusToLociIndex, int[] locusIndices, String[] siteNames) {
        return new MutableVCFAlignment(idGroup, variableSites, locusToLociIndex, locusIndices, siteNames);
    }

    public void initAllelicDepthArrays(int numMaxTaxa, int numMaxSites) {
        myAlleleDepth = new short[myMaxNumAlleles][numMaxTaxa][numMaxSites];
        myCommonAlleles = new byte[myMaxNumAlleles][numMaxSites];
        for (int i = 0; i < myAlleleDepth.length; i++) {
            for (int j = 0; j < myAlleleDepth[i].length; j++) {
                for (int k = 0; k < myAlleleDepth[i][j].length; k++) {
                    myAlleleDepth[i][j][k] = (short) -1;
                }
            }
            for (int j = 0; j < myCommonAlleles[i].length; j++) {
                myCommonAlleles[i][j] = (byte) -1;
            }
        }
    }

    @Override
    public byte[] getAllelesByScope(ALLELE_SCOPE_TYPE scope, int site) {
        if (scope == ALLELE_SCOPE_TYPE.Depth) {
            ArrayList<Byte> outArray = new ArrayList<Byte>();
            for (int i = 0; i < myCommonAlleles.length; i++) {
                if (myCommonAlleles[i][site] != (byte) -1) {
                    outArray.add(myCommonAlleles[i][site]);
                } else {
                    break;
                }
            }
            byte[] out = new byte[outArray.size()];
            for (int i = 0; i < outArray.size(); i++) {
                out[i] = outArray.get(i).byteValue();
            }
            return out;
        } else {
            return super.getAllelesByScope(scope, site);
        }
    }

    @Override
    public short[] getDepthForAlleles(int taxon, int site) {
        ArrayList<Short> outArray = new ArrayList<Short>();
        for (int i = 0; i < myAlleleDepth.length; i++) {
            if (myAlleleDepth[i][taxon][site] != (short) -1) {
                outArray.add(myAlleleDepth[i][taxon][site]);
            } else {
                break;
            }
        }
        short[] out = new short[outArray.size()];
        for (int i = 0; i < outArray.size(); i++) {
            out[i] = outArray.get(i).shortValue();
        }
        return out;
    }

    @Override
    public void setCommonAlleles(int site, byte[] values) {
        for (int i = 0; i < values.length; i++) {
            myCommonAlleles[i][site] = values[i];
        }
    }

    @Override
    public void setDepthForAlleles(int taxon, int site, short[] values) {
        for (int i = 0; i < values.length; i++) {
            myAlleleDepth[i][taxon][site] = values[i];
        }
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

                byte al;
                for (int i = 0; i < myMaxNumAlleles; i++) {
                    al = myCommonAlleles[i][a];
                    myCommonAlleles[i][a] = myCommonAlleles[i][b];
                    myCommonAlleles[i][b] = al;
                }

                bt = myReference[a];
                myReference[a] = myReference[b];
                myReference[b] = bt;

                short dp;
                for (int i = 0; i < myMaxNumAlleles; i++) {
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
