package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.IdGroupUtils;

import net.maizegenetics.pal.report.TableReport;

/**
 * @author terry
 */
public class MarkerPhenotype implements TableReport {

    private Alignment myAlignment;  //alignment to hold markers
    private Phenotype myPhenotype;  //alignment to hold characters

    /**
     * Constructor
     *
     * @param alignment alignment
     * @param phenotype phenotype
     */
    private MarkerPhenotype(Alignment alignment, Phenotype phenotype) {
        myAlignment = alignment;
        myPhenotype = phenotype;
    }

    public static MarkerPhenotype getInstance(Alignment aa, Phenotype ca, boolean union) {
        IdGroup idGroup = getIdGroup(aa.getIdGroup(), ca.getTaxa(), union);
        Alignment align = FilterAlignment.getInstance(aa, idGroup);
        Phenotype phenotype = FilterPhenotype.getInstance(ca, idGroup, null);
        return new MarkerPhenotype(align, phenotype);
    }

    public static MarkerPhenotype getInstance(MarkerPhenotype aac, IdGroup group) {
        Alignment aa = FilterAlignment.getInstance(aac.getAlignment(), group);
        Phenotype ca = FilterPhenotype.getInstance(aac.getPhenotype(), group, null);
        return new MarkerPhenotype(aa, ca);
    }

    private static IdGroup getIdGroup(IdGroup group1, IdGroup group2, boolean union) {
        if (union) {
            return IdGroupUtils.getAllIds(group1, group2);
        } else {
            return IdGroupUtils.getCommonIds(group1, group2);
        }
    }

    public String toString() {

        StringBuilder builder = new StringBuilder();

        builder.append("  ");
        builder.append(myAlignment.getSequenceCount());
        builder.append(" ");
        builder.append(myPhenotype.getNumberOfTraits());
        builder.append("Taxa/Trait\t");
        for (int j = 0; j < myPhenotype.getNumberOfTraits(); j++) {
            builder.append(myPhenotype.getTrait(j).getName());
            builder.append("\t");
        }
        builder.append("\n");
        for (int i = 0; i < myAlignment.getSequenceCount(); i++) {
            builder.append(myAlignment.getIdGroup().getIdentifier(i).getName());
            builder.append("\t");
            for (int j = 0; j < myPhenotype.getNumberOfTraits(); j++) {
                builder.append(myPhenotype.getData(i, j));
                builder.append("\t");
            }
            builder.append(myAlignment.getBaseAsStringRow(i));
        }
        return builder.toString();
    }

    public Alignment getAlignment() {
        return myAlignment;
    }

    public Phenotype getPhenotype() {
        return myPhenotype;
    }

    public Object[] getTableColumnNames() {
        String[] basicLabels = new String[myPhenotype.getNumberOfTraits() + 2];
        basicLabels[0] = "Taxa";
        for (int c = 0; c < myPhenotype.getNumberOfTraits(); c++) {
            basicLabels[c + 1] = myPhenotype.getTrait(c).getName();
        }
        basicLabels[myPhenotype.getNumberOfTraits() + 1] = "Haplotype";
        return basicLabels;
    }

    public Object[][] getTableData() {
        return getTableData(0, getRowCount() - 1);
    }

    public String getTableTitle() {
        return "Phenotypes and Genotypes";
    }

    public int getColumnCount() {
        return myPhenotype.getNumberOfTraits() + 2;
    }

    public int getRowCount() {
        return myAlignment.getSequenceCount();
    }

    public int getElementCount() {
        return getRowCount() * getColumnCount();
    }

    public Object[] getRow(int row) {

        Object[] data;
        data = new String[myPhenotype.getNumberOfTraits() + 2];
        data[0] = myAlignment.getIdGroup().getIdentifier(row).getName();
        for (int c = 0; c < myPhenotype.getNumberOfTraits(); c++) {
            data[c + 1] = "" + myPhenotype.getData(row, c);
        }
        int siteCount = Math.min(myAlignment.getSiteCount(), 10);
        StringBuilder builder = new StringBuilder();
        builder.append(myAlignment.getBaseAsStringRange(row, 0, siteCount));
        if (myAlignment.getSiteCount() > 10) {
            builder.append("...");
        }
        data[myPhenotype.getNumberOfTraits() + 1] = builder.toString();
        return data;

    }

    public Object[][] getTableData(int start, int end) {
        if ((start < 0) || (end >= getRowCount())) {
            throw new IndexOutOfBoundsException("getTableData: start: " + start + "  end: " + end);
        }

        if (end < start) {
            return null;
        }

        Object[][] temp = new Object[end - start + 1][];

        for (int i = start; i <= end; i++) {
            temp[i] = getRow(i);
        }

        return temp;
    }

    public Object getValueAt(int row, int col) {
        int haplotypeColumn = myPhenotype.getColumnCount();
        if (col == haplotypeColumn) {
            int siteCount = Math.min(myAlignment.getSiteCount(), 10);
            StringBuilder builder = new StringBuilder();
            builder.append(myAlignment.getBaseAsStringRange(row, 0, siteCount));
            if (myAlignment.getSiteCount() > 10) {
                builder.append("...");
            }
            return builder.toString();
        }
        return myPhenotype.getValueAt(row, col);
    }
}
