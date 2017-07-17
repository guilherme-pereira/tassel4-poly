/*
 * AlignmentMask
 */
package net.maizegenetics.pal.alignment;

import java.awt.Color;

import java.io.Serializable;


/**
 *
 * @author terry
 */
public interface AlignmentMask extends Serializable {

    /**
     * This defines the type of mask.
     *
     * reference - Mask created using taxon as reference.
     * imputed - Mask created to identify imputed values.
     * compared = Mask created to identify differences between two alignments.
     */
    public enum MaskType {reference, imputed, compared};

    public byte getMask(int taxon, int site);

    public Color getColor();

    public void setColor(Color color);

    public Alignment getAlignment();

    public MaskType getMaskType();

}
