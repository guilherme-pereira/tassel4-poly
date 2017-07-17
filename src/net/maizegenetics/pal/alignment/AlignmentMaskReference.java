/*
 * AlignmentMaskReference
 */
package net.maizegenetics.pal.alignment;

import java.awt.Color;
import net.maizegenetics.pal.ids.Identifier;

/**
 *
 * @author terry
 */
public class AlignmentMaskReference extends AbstractAlignmentMask {

    private static final long serialVersionUID = -5197800047652332969L;
    private final int myTaxonReference;
    private final Alignment myAlignment;

    private AlignmentMaskReference(Alignment align, int taxonReference, String name, Color color) {
        super(align, name, color, AlignmentMask.MaskType.reference);
        myTaxonReference = taxonReference;
        myAlignment = align;
    }

    public static AlignmentMaskReference getInstanceCompareReference(Alignment align) {
        return getInstanceCompareReference(align, -1);
    }

    public static AlignmentMaskReference getInstanceCompareReference(Alignment align, Identifier id) {
        int index = align.getIdGroup().whichIdNumber(id);
        if (index == -1) {
            throw new IllegalArgumentException("AlignmentMaskReference: getInstanceCompareReference: unknown id: " + id);
        }
        return getInstanceCompareReference(align, index);
    }

    public static AlignmentMaskReference getInstanceCompareReference(Alignment align, String id) {
        int index = align.getIdGroup().whichIdNumber(id);
        if (index == -1) {
            throw new IllegalArgumentException("AlignmentMaskReference: getInstanceCompareReference: unknown id: " + id);
        }
        return getInstanceCompareReference(align, index);
    }

    public static AlignmentMaskReference getInstanceCompareReference(Alignment align, int index) {
        if ((index < -1) || (index >= align.getSequenceCount())) {
            throw new IllegalArgumentException("AlignmentMaskReference: getInstanceCompareReference: unknown index: " + index);
        }
        String name;
        if (index == -1) {
            name = "Alignment Reference";
        } else {
            name = align.getTaxaName(index) + " Reference";
        }
        return new AlignmentMaskReference(align, index, name, getNextColor());
    }

    @Override
    public byte getMask(int taxon, int site) {
        if ((myTaxonReference == -1) && (AlignmentUtils.isEqualOrUnknown(myAlignment.getBase(taxon, site), myAlignment.getReferenceAllele(site)))) {
            return 0;
        } else if (AlignmentUtils.isEqualOrUnknown(myAlignment.getBaseArray(taxon, site), myAlignment.getBaseArray(myTaxonReference, site))) {
            return 0;
        } else {
            return 1;
        }
    }
}
