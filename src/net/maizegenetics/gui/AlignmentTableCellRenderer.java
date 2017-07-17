/*
 * AlignmentTableCellRenderer
 */
package net.maizegenetics.gui;

import java.awt.Color;
import java.awt.Component;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import javax.swing.JTable;
import javax.swing.SwingConstants;
import javax.swing.table.DefaultTableCellRenderer;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AlignmentMask;
import net.maizegenetics.pal.alignment.AlignmentMaskGeneticDistance;
import net.maizegenetics.pal.alignment.AlignmentMaskReference;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;

/**
 *
 * @author terry
 */
public class AlignmentTableCellRenderer extends DefaultTableCellRenderer {

    private static int[] NUCLEOTIDE_COLORS = new int[16];

    static {
        NUCLEOTIDE_COLORS[NucleotideAlignmentConstants.A_ALLELE] = 0xFA0000;
        NUCLEOTIDE_COLORS[NucleotideAlignmentConstants.C_ALLELE] = 0x9BCD9B;
        NUCLEOTIDE_COLORS[NucleotideAlignmentConstants.G_ALLELE] = 0x4876FF;
        NUCLEOTIDE_COLORS[NucleotideAlignmentConstants.T_ALLELE] = 0x777D7E;
        NUCLEOTIDE_COLORS[NucleotideAlignmentConstants.GAP_ALLELE] = 0xFFF68F;
        NUCLEOTIDE_COLORS[NucleotideAlignmentConstants.INSERT_ALLELE] = 0xFF8C00;
    }
    private static final Color MAJOR_ALLELE_COLOR = new Color(0xe6cf45);
    private static final Color MINOR_ALLELE_COLOR = new Color(0x45a1e6);
    private static final Color HETEROZYGOUS_COLOR = new Color(0xe64545);
    private static final Color MAJOR_MINOR_ALLELE_COLOR = new Color((MAJOR_ALLELE_COLOR.getRGB() + MINOR_ALLELE_COLOR.getRGB()) % 0xFFFFFF);
    private static final Color[] COLORS_256 = generateColors(256);
    private static final Map<String, Color> COLORS_NUCLEOTIDES = generateNucleotideColors();

    public static enum RENDERING_TYPE {

        Nucleotide, NucleotideHeterozygous, MajorAllele, MinorAllele, MajorMinorAllele, Heterozygous, ReferenceMasks, GeneticDistanceMasks, None
    };
    private final AlignmentTableModel myAlignmentTableModel;
    private final Alignment myAlignment;
    private AlignmentMask[] myMasks;
    private RENDERING_TYPE myRenderingType = RENDERING_TYPE.Nucleotide;
    private final Map<Integer, byte[]> myCachedAlleles = new LinkedHashMap<Integer, byte[]>() {
        protected boolean removeEldestEntry(Map.Entry eldest) {
            return size() > 100;
        }
    };

    public AlignmentTableCellRenderer(AlignmentTableModel model, Alignment alignment, AlignmentMask[] masks) {
        myAlignmentTableModel = model;
        myAlignment = alignment;
        myMasks = masks;
    }

    @Override
    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int col) {

        switch (myRenderingType) {
            case Nucleotide:
                return getNucleotideRendering(table, value, isSelected, hasFocus, row, col);
            case NucleotideHeterozygous:
                return getNucleotideHeterozygousRendering(table, value, isSelected, hasFocus, row, col);
            case MajorAllele:
                return getMajorAlleleRendering(table, value, isSelected, hasFocus, row, col);
            case MinorAllele:
                return getMinorAlleleRendering(table, value, isSelected, hasFocus, row, col);
            case Heterozygous:
                return getHeterozygousRendering(table, value, isSelected, hasFocus, row, col);
            case MajorMinorAllele:
                return getMajorMinorAlleleRendering(table, value, isSelected, hasFocus, row, col);
            case ReferenceMasks:
                return getReferenceMasksRendering(table, value, isSelected, hasFocus, row, col);
            case GeneticDistanceMasks:
                return getGeneticDistanceMasksRendering(table, value, isSelected, hasFocus, row, col);
            case None:
                return getDefaultRendering(table, value, isSelected, hasFocus, row, col);
            default:
                return getDefaultRendering(table, value, isSelected, hasFocus, row, col);
        }

    }

    private Component getDefaultRendering(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int col) {

        Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, col);

        setHorizontalAlignment(SwingConstants.CENTER);

        if (isSelected) {
            comp.setBackground(Color.DARK_GRAY);
        } else {
            comp.setBackground(null);
        }

        return comp;

    }

    private Component getNucleotideRendering(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int col) {

        Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, col);

        setHorizontalAlignment(SwingConstants.CENTER);

        String alleles = myAlignment.getBaseAsString(row, myAlignmentTableModel.getRealColumnIndex(col));

        if (isSelected) {
            comp.setBackground(Color.DARK_GRAY);
        } else if (alleles.equals(Alignment.UNKNOWN_ALLELE_STR)) {
            comp.setBackground(null);
        } else {
            comp.setBackground(COLORS_NUCLEOTIDES.get(alleles));
        }

        return comp;

    }

    private Component getNucleotideHeterozygousRendering(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int col) {

        Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, col);

        setHorizontalAlignment(SwingConstants.CENTER);

        int site = myAlignmentTableModel.getRealColumnIndex(col);

        if (isSelected) {
            comp.setBackground(Color.DARK_GRAY);
        } else if (myAlignment.isHeterozygous(row, site)) {
            String alleles = myAlignment.getBaseAsString(row, myAlignmentTableModel.getRealColumnIndex(col));
            comp.setBackground(COLORS_NUCLEOTIDES.get(alleles));
        } else {
            comp.setBackground(null);
        }

        return comp;

    }

    private Component getMajorAlleleRendering(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int col) {

        Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, col);

        setHorizontalAlignment(SwingConstants.CENTER);

        int site = myAlignmentTableModel.getRealColumnIndex(col);
        byte[] alleles = myCachedAlleles.get(site);
        if (alleles == null) {
            alleles = myAlignment.getAlleles(site);
            myCachedAlleles.put(site, alleles);
        }
        byte major = Alignment.UNKNOWN_ALLELE;
        if (alleles.length > 0) {
            major = alleles[0];
        }
        byte[] diploidValues = myAlignment.getBaseArray(row, site);

        if (isSelected) {
            comp.setBackground(Color.DARK_GRAY);
        } else if ((diploidValues[0] == major) || (diploidValues[1] == major)) {
            comp.setBackground(MAJOR_ALLELE_COLOR);
        } else {
            comp.setBackground(null);
        }

        return comp;

    }

    private Component getHeterozygousRendering(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int col) {

        Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, col);

        setHorizontalAlignment(SwingConstants.CENTER);

        int site = myAlignmentTableModel.getRealColumnIndex(col);

        if (isSelected) {
            comp.setBackground(Color.DARK_GRAY);
        } else if (myAlignment.isHeterozygous(row, site)) {
            comp.setBackground(HETEROZYGOUS_COLOR);
        } else {
            comp.setBackground(null);
        }

        return comp;

    }

    private Component getMinorAlleleRendering(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int col) {

        Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, col);

        setHorizontalAlignment(SwingConstants.CENTER);

        int site = myAlignmentTableModel.getRealColumnIndex(col);
        byte[] alleles = myCachedAlleles.get(site);
        if (alleles == null) {
            alleles = myAlignment.getAlleles(site);
            myCachedAlleles.put(site, alleles);
        }
        byte minor = Alignment.UNKNOWN_ALLELE;
        if (alleles.length > 1) {
            minor = alleles[1];
        }
        byte[] diploidValues = myAlignment.getBaseArray(row, site);

        if (isSelected) {
            comp.setBackground(Color.DARK_GRAY);
        } else if ((diploidValues[0] == minor) || (diploidValues[1] == minor)) {
            comp.setBackground(MINOR_ALLELE_COLOR);
        } else {
            comp.setBackground(null);
        }

        return comp;

    }

    private Component getMajorMinorAlleleRendering(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int col) {

        Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, col);

        setHorizontalAlignment(SwingConstants.CENTER);

        int site = myAlignmentTableModel.getRealColumnIndex(col);
        byte[] alleles = myCachedAlleles.get(site);
        if (alleles == null) {
            alleles = myAlignment.getAlleles(site);
            if (alleles.length > 1) {
                if (alleles[0] == Alignment.UNKNOWN_ALLELE) {
                    alleles = new byte[0];
                } else if (alleles[1] == Alignment.UNKNOWN_ALLELE) {
                    alleles = new byte[]{alleles[0]};
                }
            } else if (alleles.length == 1) {
                if (alleles[0] == Alignment.UNKNOWN_ALLELE) {
                    alleles = new byte[0];
                }
            }
            myCachedAlleles.put(site, alleles);
        }

        byte[] diploidValues = myAlignment.getBaseArray(row, site);
        if (alleles.length > 1) {
            byte major = alleles[0];
            byte minor = alleles[1];
            if (isSelected) {
                comp.setBackground(Color.DARK_GRAY);
            } else if (((diploidValues[0] == major) && (diploidValues[1] == minor))
                    || ((diploidValues[0] == minor) && (diploidValues[1] == major))) {
                comp.setBackground(MAJOR_MINOR_ALLELE_COLOR);
            } else if ((diploidValues[0] == major) || (diploidValues[1] == major)) {
                comp.setBackground(MAJOR_ALLELE_COLOR);
            } else if ((diploidValues[0] == minor) || (diploidValues[1] == minor)) {
                comp.setBackground(MINOR_ALLELE_COLOR);
            } else {
                comp.setBackground(null);
            }
        } else if (alleles.length == 1) {
            byte major = alleles[0];
            if (isSelected) {
                comp.setBackground(Color.DARK_GRAY);
            } else if ((diploidValues[0] == major) || (diploidValues[1] == major)) {
                comp.setBackground(MAJOR_ALLELE_COLOR);
            } else {
                comp.setBackground(null);
            }
        } else {
            if (isSelected) {
                comp.setBackground(Color.DARK_GRAY);
            } else {
                comp.setBackground(null);
            }
        }

        return comp;

    }

    public Component getReferenceMasksRendering(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int col) {

        Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, col);

        setHorizontalAlignment(SwingConstants.CENTER);

        AlignmentMask[] masks = getAlignmentMasksOfClass(AlignmentMaskReference.class);

        if (isSelected) {
            comp.setBackground(Color.DARK_GRAY);
        } else if ((masks == null) || (masks.length == 0)) {
            comp.setBackground(null);
        } else if (masks.length == 1) {

            if (masks[0].getMask(row, myAlignmentTableModel.getRealColumnIndex(col)) == 0x1) {
                comp.setBackground(masks[0].getColor());
            } else {
                comp.setBackground(null);
            }

        } else {

            int red = 0;
            int green = 0;
            int blue = 0;

            boolean changed = false;
            for (int i = 0; i < masks.length; i++) {
                if (masks[i].getMask(row, myAlignmentTableModel.getRealColumnIndex(col)) == 0x1) {
                    red = red + masks[i].getColor().getRed();
                    green = green + masks[i].getColor().getGreen();
                    blue = blue + masks[i].getColor().getBlue();
                    changed = true;
                }
            }

            if (changed) {
                red = red % 256;
                green = green % 256;
                blue = blue % 256;
                comp.setBackground(new Color(red, green, blue));
            } else {
                comp.setBackground(null);
            }

        }

        return comp;

    }

    public Component getGeneticDistanceMasksRendering(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int col) {

        Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, col);

        setHorizontalAlignment(SwingConstants.CENTER);

        AlignmentMask[] masks = getAlignmentMasksOfClass(AlignmentMaskGeneticDistance.class);

        if (isSelected) {
            comp.setBackground(Color.DARK_GRAY);
        } else if ((masks == null) || (masks.length == 0)) {
            comp.setBackground(null);
        } else {

            int aveDistance = 0;
            for (int i = 0; i < masks.length; i++) {
                aveDistance += 0xFF & masks[i].getMask(row, myAlignmentTableModel.getRealColumnIndex(col));
            }
            aveDistance /= masks.length;

            comp.setBackground(COLORS_256[aveDistance]);

        }

        return comp;

    }

    private static Color[] generateColors(int n) {
        Color[] cols = new Color[n];
        for (int i = 0; i < n; i++) {
            cols[i] = Color.getHSBColor((float) i / (float) n, 0.85f, 1.0f);
        }
        return cols;
    }

    private static Map<String, Color> generateNucleotideColors() {
        Map<String, Color> result = new HashMap<String, Color>();
        Iterator itr = NucleotideAlignmentConstants.NUCLEOTIDE_IUPAC_HASH.entrySet().iterator();
        while (itr.hasNext()) {
            Map.Entry current = (Map.Entry) itr.next();
            result.put((String) current.getValue(), Color.BLACK);
        }
        int numCodes = result.size();
        itr = result.entrySet().iterator();
        int count = 0;
        while (itr.hasNext()) {
            Map.Entry current = (Map.Entry) itr.next();
            Color color = Color.getHSBColor((float) count / (float) (numCodes - 1), 0.7f, 0.9f);
            result.put((String) current.getKey(), color);
            count++;
        }
        return result;
    }

    private AlignmentMask[] getAlignmentMasksOfClass(Class type) {
        if ((myMasks == null) || (myMasks.length == 0)) {
            return null;
        }
        List<AlignmentMask> masks = new ArrayList<AlignmentMask>();
        for (int i = 0; i < myMasks.length; i++) {
            if (type.isInstance(myMasks[i])) {
                masks.add(myMasks[i]);
            }
        }
        AlignmentMask[] result = new AlignmentMask[masks.size()];
        masks.toArray(result);
        return result;
    }

    public RENDERING_TYPE getRenderingType() {
        return myRenderingType;
    }

    public void setRenderingType(RENDERING_TYPE type) {
        myRenderingType = type;
    }

    public void setMasks(AlignmentMask[] masks) {
        myMasks = masks;
    }
}
