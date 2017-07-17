/*
 * AlignmentMaskBoolean
 */
package net.maizegenetics.pal.alignment;

import java.awt.Color;

import net.maizegenetics.pal.ids.Identifier;

/**
 *
 * @author terry
 */
public class AlignmentMaskBoolean extends AbstractAlignmentMask {

    private static final long serialVersionUID = -5197800047652332969L;
    private final byte[][] myMask;
    

    public AlignmentMaskBoolean(Alignment align, byte[][] mask, String name, MaskType type) {
        this(align, mask, name, getNextColor(), type);
    }

    public AlignmentMaskBoolean(Alignment align, byte[][] mask, String name, Color color, MaskType type) {

        super(align, name, color, type);

        if (mask.length != align.getSequenceCount()) {
            throw new IllegalArgumentException("AlignmentMask: init: number of mask rows should equal number of sequences.");
        }

        int numBytesNeeded = getNumMaskColumns(align.getSiteCount());
        if (numBytesNeeded != mask[0].length) {
            throw new IllegalArgumentException("AlignmentMask: init: incorrect number of mask columns: " + mask[0].length + "  should be: " + numBytesNeeded);
        }

        myMask = mask;

    }

    public static AlignmentMaskBoolean getInstanceCompareReference(Alignment align) {
        String name = align.getIdGroup().getIdentifier(0).getName() + " Reference";
        return getInstanceCompareReference(align, align.getBaseRange(0, 0, align.getSiteCount()), name);
    }

    public static AlignmentMaskBoolean getInstanceCompareReference(Alignment align, Identifier id) {
        int index = align.getIdGroup().whichIdNumber(id);
        if (index == -1) {
            throw new IllegalArgumentException("AlignmentMask: getInstanceCompareReference: unknown id: " + id);
        }
        String name = id.getName() + " Reference";
        return getInstanceCompareReference(align, align.getBaseRange(index, 0, align.getSiteCount()), name);
    }

    public static AlignmentMaskBoolean getInstanceCompareReference(Alignment align, int index) {
        if ((index < 0) || (index >= align.getSequenceCount())) {
            throw new IllegalArgumentException("AlignmentMask: getInstanceCompareReference: unknown index: " + index);
        }
        String name = align.getTaxaName(index) + " Reference";
        return getInstanceCompareReference(align, align.getBaseRange(index, 0, align.getSiteCount()), name);
    }

    public static AlignmentMaskBoolean getInstanceCompareReference(Alignment align, String id) {
        int index = align.getIdGroup().whichIdNumber(id);
        if (index == -1) {
            throw new IllegalArgumentException("AlignmentMask: getInstanceCompareReference: unknown id: " + id);
        }
        return getInstanceCompareReference(align, align.getBaseRange(index, 0, align.getSiteCount()), id + " Reference");
    }

    public static AlignmentMaskBoolean getInstanceCompareReference(Alignment align, byte[] ref, String name) {

        if ((align == null) || (ref == null)) {
            throw new IllegalArgumentException("AlignmentMask: getInstanceCompareReference: alignment or reference can not be null.");
        }

        if (align.getSiteCount() != ref.length) {
            throw new IllegalArgumentException("AlignmentMask: getInstanceCompareReference: ref length should equal alignment site count.");
        }

        int numMaskColumns = getNumMaskColumns(ref.length);
        byte[][] mask = new byte[align.getSequenceCount()][numMaskColumns];

        for (int c = 0, n = align.getSiteCount(); c < n; c++) {

            int currentByteCol = c / 8;
            byte currentColMask = (byte) (0x80 >>> (c % 8));

            for (int r = 0, m = align.getSequenceCount(); r < m; r++) {

                // TERRY - If allele order switch it should still probably match?
                if (align.getBase(r, c) != ref[c]) {
                    mask[r][currentByteCol] = (byte) (mask[r][currentByteCol] | currentColMask);
                }

            }

        }

        return new AlignmentMaskBoolean(align, mask, name, MaskType.reference);

    }

    public static AlignmentMaskBoolean getInstanceCompareAlignments(Alignment align1, Alignment align2, String name, MaskType type) {

        if ((align1.getSequenceCount() != align2.getSequenceCount()) ||
                (align1.getSiteCount() != align2.getSiteCount())) {
            throw new IllegalArgumentException("AlignmentMaskBoolean: getInstanceCompareAlignments: both alignments should have same number of sequences and sites.");
        }

        int numMaskColumns = getNumMaskColumns(align1.getSiteCount());
        byte[][] mask = new byte[align1.getSequenceCount()][numMaskColumns];

        for (int c = 0, n = align1.getSiteCount(); c < n; c++) {

            int currentByteCol = c / 8;
            byte currentColMask = (byte) (0x80 >>> (c % 8));

            for (int r = 0, m = align1.getSequenceCount(); r < m; r++) {

                // TERRY - If allele order switch it should still probably match?
                if (align1.getBase(r, c) != align2.getBase(r, c)) {
                    mask[r][currentByteCol] = (byte) (mask[r][currentByteCol] | currentColMask);
                }

            }

        }

        return new AlignmentMaskBoolean(align1, mask, name, type);

    }

    public static int getNumMaskColumns(int numSites) {
        int numMaskColumns = numSites / 8;
        if (numSites % 8 > 0) {
            numMaskColumns++;
        }
        return numMaskColumns;
    }

    public byte getMask(int taxon, int site) {
        int maskColumn = site / 8;
        int shift = 7 - (site % 8);
        return (byte) ((myMask[taxon][maskColumn] >>> shift) & 0x1);
    }
}
