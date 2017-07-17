package net.maizegenetics.pal.popgen;

import net.maizegenetics.pal.report.TableReport;
import net.maizegenetics.pal.report.AbstractTableReport;

import java.io.Serializable;
import java.util.Vector;
import net.maizegenetics.pal.alignment.Alignment;

/**
 *This class provides the distribution of polymorphisms
 *
 * @author Ed Buckler
 * @version 1.0
 */
public class PolymorphismDistribution extends AbstractTableReport implements TableReport, Serializable {

    Vector polyDistResultsVector = new Vector();
    int maxSeqCount = -1;

    public PolymorphismDistribution() {
    }

    //consider whether to output pooled minority alleles or everything
    public void addDistribution(String label, Alignment theSP, boolean poolMinor) {
        maxSeqCount = theSP.getSequenceCount() * 2;
        int[] pdist = new int[maxSeqCount];
        for (int i = 0; i < theSP.getSiteCount(); i++) {
            if (theSP.isPolymorphic(i)) {
                int[][] alleleCounts = theSP.getAllelesSortedByFrequency(i);
                int numAlleles = alleleCounts[0].length;
                if ((poolMinor == false) || (numAlleles == 2)) {
                    pdist[alleleCounts[1][1]]++;
                } else {
                    int sum = 0;
                    for (int a = 1; a < numAlleles; a++) {
                        sum += alleleCounts[1][a];
                    }
                    pdist[sum]++;
                }
            } else {
                pdist[0]++;
            }
        }

        PolymorphismDistributionResults pdr = new PolymorphismDistributionResults(label, pdist, poolMinor);
        polyDistResultsVector.add(pdr);
    }

    //Implementation of TableReport Interface
    public Object[] getTableColumnNames() {
        String[] basicLabels = new String[1 + polyDistResultsVector.size()];
        basicLabels[0] = "Site_Freq";
        PolymorphismDistributionResults pdr;
        for (int i = 0; i < polyDistResultsVector.size(); i++) {
            pdr = (PolymorphismDistributionResults) polyDistResultsVector.get(i);
            basicLabels[i + 1] = pdr.label;
        }
        return basicLabels;
    }

    public Object[][] getTableData() {
        Object[][] data;
        int basicCols = 1, labelOffset;
        PolymorphismDistributionResults pdr;
        data = new String[maxSeqCount + 1][basicCols + polyDistResultsVector.size()];
        data[0][0] = "N";
        //label along side for frequency
        for (int i = 0; i < maxSeqCount; i++) {
            data[i + 1][0] = "" + i;
        }
        for (int i = 0; i < polyDistResultsVector.size(); i++) {
            pdr = (PolymorphismDistributionResults) polyDistResultsVector.get(i);
            //label of taxa sample across the top
            data[0][i + 1] = "" + pdr.polyDist.length;
            for (int j = 0; j < pdr.polyDist.length; j++) {
                data[j + 1][i + 1] = "" + pdr.polyDist[j];
            }
        }
        return data;
    }

    /**
     * Returns specified row.
     *
     * @param row row number
     *
     * @return row
     */
    public Object[] getRow(int row) {
        throw new UnsupportedOperationException();
    }

    public String getTableTitle() {
        return "Polymorphism Distribution";
    }

    public String toString() {
        if (polyDistResultsVector.size() == 0) {
            return "Needs to be run";
        }
        StringBuffer cs = new StringBuffer();
        Object[] header = this.getTableColumnNames();
        for (int i = 0; i < header.length; i++) {
            cs.append(header[i]);
        }
        cs.append("\n");
        Object[][] data = this.getTableData();
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {//System.out.println("i="+i+" j="+j+" data="+(String)data[i][j]);
                cs.append(data[i][j]);
            }
            cs.append("\n");
        }
        return cs.toString();
    }

    public int getRowCount() {
        return maxSeqCount + 1;
    }

    public int getElementCount() {
        throw new UnsupportedOperationException();
    }

    public int getColumnCount() {
        throw new UnsupportedOperationException();
    }
}

class PolymorphismDistributionResults implements Serializable {

    protected int[] polyDist;
    protected String label;
    protected boolean poolMinor;
    private int index;

    public PolymorphismDistributionResults(String label, int[] dist, boolean poolMinor) {
        this.label = label;
        this.poolMinor = poolMinor;
        this.polyDist = dist;
    }

    public boolean equals(Object anObject) {
        PolymorphismDistributionResults x = (PolymorphismDistributionResults) anObject;
        if (x.index == this.index) {
            return true;
        } else {
            return false;
        }
    }

    public void setIndex(int theIndex) {
        this.index = theIndex;
    }

    public String toString() {
        String result = "Polymorphism Distribution for " + label + "\n";
        result += "Pool Minor =" + poolMinor + "\n";
        result += "Dist =" + polyDist.toString() + "\n";
        return result;
    }
}
