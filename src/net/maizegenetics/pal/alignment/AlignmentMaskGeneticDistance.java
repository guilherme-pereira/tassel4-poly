/*
 * AlignmentMaskGeneticDistance
 */
package net.maizegenetics.pal.alignment;

import java.awt.Color;
import java.util.LinkedHashMap;
import java.util.Map;
import net.maizegenetics.pal.distance.IBSDistanceMatrix;
import net.maizegenetics.pal.ids.Identifier;

/**
 *
 * @author terry
 */
public class AlignmentMaskGeneticDistance extends AbstractAlignmentMask {

    private static final long serialVersionUID = -5197800047652332969L;
    private Map<Integer, Byte> myCache = new LinkedHashMap<Integer, Byte>() {
        protected boolean removeEldestEntry(Map.Entry eldest) {
            return size() > 100;
        }
    };
    private final int myTaxonReference;
    private final Alignment myAlignment;
    private Alignment myTBitAlignment = null;

    private AlignmentMaskGeneticDistance(Alignment align, int taxonReference, String name, Color color) {
        super(align, name, color, AlignmentMask.MaskType.reference);
        myTaxonReference = taxonReference;
        myAlignment = align;
    }

    public static AlignmentMaskGeneticDistance getInstanceCompareReference(Alignment align, Identifier id) {
        int index = align.getIdGroup().whichIdNumber(id);
        if (index == -1) {
            throw new IllegalArgumentException("AlignmentMaskGeneticDistance: getInstanceCompareReference: unknown id: " + id);
        }
        return getInstanceCompareReference(align, index);
    }

    public static AlignmentMaskGeneticDistance getInstanceCompareReference(Alignment align, String id) {
        int index = align.getIdGroup().whichIdNumber(id);
        if (index == -1) {
            throw new IllegalArgumentException("AlignmentMaskGeneticDistance: getInstanceCompareReference: unknown id: " + id);
        }
        return getInstanceCompareReference(align, index);
    }

    public static AlignmentMaskGeneticDistance getInstanceCompareReference(Alignment align, int index) {
        if ((index < 0) || (index >= align.getSequenceCount())) {
            throw new IllegalArgumentException("AlignmentMaskGeneticDistance: getInstanceCompareReference: unknown index: " + index);
        }
        String name = align.getTaxaName(index) + " Genetic Distance";
        return new AlignmentMaskGeneticDistance(align, index, name, null);
    }

    @Override
    public byte getMask(int taxon, int site) {

        Byte result = myCache.get(taxon);
        if (result != null) {
            return result;
        }

        if (myTBitAlignment == null) {
            myTBitAlignment = AlignmentUtils.optimizeForTaxa(myAlignment);
        }

        result = (byte) (IBSDistanceMatrix.computeHetBitDistances(myTBitAlignment, taxon, myTaxonReference)[0] * 255.0);
        myCache.put(taxon, result);
        return result;

    }
}
