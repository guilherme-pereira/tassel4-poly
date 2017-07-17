/*
 * GenotypeSummaryPlugin
 */
package net.maizegenetics.baseplugins;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.Frame;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import java.net.URL;

import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTabbedPane;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.report.SimpleTableReport;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;

import org.apache.log4j.Logger;

/**
 *
 * @author terry
 */
public class GenotypeSummaryPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(GenotypeSummaryPlugin.class);
    private static final String NA = "NA";
    private static final Double ZERO_DOUBLE = 0.0;
    private static final int ZERO_INT = 0;
    private long myNumGametesMissing = 0;
    private long myNumHeterozygous = 0;
    private boolean myIsOverview = true;
    private boolean myIsSiteSummary = true;
    private boolean myIsTaxaSummary = true;

    public GenotypeSummaryPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet performFunction(DataSet input) {

        try {

            myNumGametesMissing = 0;
            myNumHeterozygous = 0;

            List<Datum> alignInList = input.getDataOfType(Alignment.class);

            if (alignInList.size() != 1) {
                String gpMessage = "Invalid selection.  Please select one genotype alignment.";
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), gpMessage);
                } else {
                    myLogger.error(gpMessage);
                }
                return null;
            }

            if (isInteractive()) {
                GenotypeSummaryPluginDialog theDialog = new GenotypeSummaryPluginDialog();
                theDialog.setLocationRelativeTo(getParentFrame());
                theDialog.setVisible(true);
                if (theDialog.isCancel()) {
                    return null;
                }
                myIsOverview = theDialog.isOverview();
                myIsSiteSummary = theDialog.isSiteSummary();
                myIsTaxaSummary = theDialog.isTaxaSummary();
                theDialog.dispose();
            }

            Datum current = alignInList.get(0);
            Alignment alignment = (Alignment) current.getData();
            String name = current.getName();

            List<Datum> summaryTables = new ArrayList<Datum>();


            SimpleTableReport siteSummary = null;
            if (myIsSiteSummary) {
                siteSummary = getSiteSummary(alignment);
            }

            SimpleTableReport taxaSummary = null;
            if (myIsTaxaSummary) {
                taxaSummary = getTaxaSummary(alignment);
            }

            SimpleTableReport[] overallSummaries = null;
            if (myIsOverview) {
                overallSummaries = getOverallSummary(alignment);
            }


            if (myIsOverview) {
                summaryTables.add(new Datum(name + "_OverallSummary", overallSummaries[0], "Overall Summary of " + name));
                summaryTables.add(new Datum(name + "_AlleleSummary", overallSummaries[1], "Allele Summary of " + name));
            }

            if (siteSummary != null) {
                summaryTables.add(new Datum(name + "_SiteSummary", siteSummary, "Site Summary of " + name));
            }
            if (taxaSummary != null) {
                summaryTables.add(new Datum(name + "_TaxaSummary", taxaSummary, "Taxa Summary of " + name));
            }

            if (summaryTables.isEmpty()) {
                return null;
            }

            DataSet output = new DataSet(summaryTables, this);
            fireDataSetReturned(new PluginEvent(output, GenotypeSummaryPlugin.class));

            return output;

        } finally {
            fireProgress(100);
        }

    }

    private SimpleTableReport[] getOverallSummary(Alignment alignment) {

        Object[] firstColumnNames = new String[]{"Stat Type", "Value"};

        long numSites = alignment.getSiteCount();
        long numTaxa = alignment.getSequenceCount();

        Object[][] diploidValueCounts = alignment.getDiploidCounts();
        int numAlleles = diploidValueCounts[0].length;

        if (!myIsSiteSummary) {
            int totalGametes = (int) numTaxa * 2;
            for (int i = 0; i < numSites; i++) {
                int totalGametesNotMissing = alignment.getTotalGametesNotMissing(i);
                int totalGametesMissing = totalGametes - totalGametesNotMissing;
                myNumGametesMissing = myNumGametesMissing + (long) totalGametesMissing;
                int numHeterozygous = alignment.getHeterozygousCount(i);
                myNumHeterozygous = myNumHeterozygous + (long) numHeterozygous;
            }
        }

        long totalGametes = numSites * numTaxa * 2L;
        long totalGametesNotMissing = totalGametes - myNumGametesMissing;

        long numDiploidsMissing = 0;
        for (int j = 0; j < numAlleles; j++) {
            if ((diploidValueCounts[0][j].equals(Alignment.UNKNOWN_ALLELE_STR)) || (diploidValueCounts[0][j].equals(Alignment.UNKNOWN_DIPLOID_ALLELE_STR))) {
                numDiploidsMissing = (Long) diploidValueCounts[1][j];
                break;
            }
        }

        long totalDiploids = numSites * numTaxa;
        long totalDiploidsNotMissing = totalDiploids - numDiploidsMissing;
        int count = 0;

        Object[][] data = new Object[14][firstColumnNames.length];

        data[count][0] = "Number of Taxa";
        data[count++][1] = (double) numTaxa;

        data[count][0] = "Number of Sites";
        data[count++][1] = (double) numSites;

        data[count][0] = "Sites x Taxa";
        data[count++][1] = (double) totalDiploids;

        data[count][0] = "Number Not Missing";
        data[count++][1] = (double) totalDiploidsNotMissing;

        data[count][0] = "Proportion Not Missing";
        data[count++][1] = (double) totalDiploidsNotMissing / (double) totalDiploids;

        data[count][0] = "Number Missing";
        data[count++][1] = (double) numDiploidsMissing;

        data[count][0] = "Proportion Missing";
        data[count++][1] = (double) numDiploidsMissing / (double) totalDiploids;

        data[count][0] = "Number Gametes";
        data[count++][1] = (double) totalGametes;

        data[count][0] = "Gametes Not Missing";
        data[count++][1] = (double) totalGametesNotMissing;

        data[count][0] = "Proportion Gametes Not Missing";
        data[count++][1] = (double) totalGametesNotMissing / (double) totalGametes;

        data[count][0] = "Gametes Missing";
        data[count++][1] = (double) myNumGametesMissing;

        data[count][0] = "Proportion Gametes Missing";
        data[count++][1] = (double) myNumGametesMissing / (double) totalGametes;

        data[count][0] = "Number Heterozygous";
        data[count++][1] = (double) myNumHeterozygous;

        data[count][0] = "Proportion Heterozygous";
        data[count++][1] = (double) myNumHeterozygous / (double) totalDiploids;


        Object[][] majorMinorDiploidValueCounts = alignment.getMajorMinorCounts();
        int numMajorMinorAlleles = majorMinorDiploidValueCounts[0].length;

        Object[] alleleColumnNames = new String[]{"Alleles", "Number", "Proportion", "Frequency"};
        Object[][] data2 = new Object[numAlleles + numMajorMinorAlleles][alleleColumnNames.length];

        count = 0;
        for (int i = 0; i < numAlleles; i++) {
            String value = (String) diploidValueCounts[0][i];
            Long numValue = (Long) diploidValueCounts[1][i];
            data2[count][0] = value;
            data2[count][1] = numValue;
            data2[count][2] = numValue.doubleValue() / (double) totalDiploids;
            data2[count++][3] = numValue.doubleValue() / (double) totalDiploidsNotMissing;
        }

        numDiploidsMissing = 0;
        for (int j = 0; j < numMajorMinorAlleles; j++) {
            if ((majorMinorDiploidValueCounts[0][j].equals(Alignment.UNKNOWN_ALLELE_STR)) || (majorMinorDiploidValueCounts[0][j].equals(Alignment.UNKNOWN_DIPLOID_ALLELE_STR))) {
                numDiploidsMissing = (Long) majorMinorDiploidValueCounts[1][j];
                break;
            }
        }

        for (int i = 0; i < numMajorMinorAlleles; i++) {
            String value = (String) majorMinorDiploidValueCounts[0][i];
            Long numValue = (Long) majorMinorDiploidValueCounts[1][i];
            data2[count][0] = value;
            data2[count][1] = numValue;
            data2[count++][2] = numValue.doubleValue() / (double) numSites;
        }

        return new SimpleTableReport[]{new SimpleTableReport("Overall Summary", firstColumnNames, data),
            new SimpleTableReport("Allele Summary", alleleColumnNames, data2)};
    }

    private SimpleTableReport getSiteSummary(Alignment alignment) {

        String[] firstColumnNames = new String[]{"Site Number", "Site Name", "Physical Position", "Number of Taxa", "Major Allele", "Major Allele Gametes", "Major Allele Proportion", "Major Allele Frequency",
            "Minor Allele", "Minor Allele Gametes", "Minor Allele Proportion", "Minor Allele Frequency"};
        String[] lastColumnNames = new String[]{"Gametes Missing", "Proportion Missing", "Number Heterozygous", "Proportion Heterozygous",
            "Inbreeding Coefficient", "Inbreeding Coefficient Scaled by Missing"};

        List<String> columnNames = new ArrayList<String>(Arrays.asList(firstColumnNames));

        int maxAlleles = alignment.getMaxNumAlleles();
        if (alignment.retainsRareAlleles()) {
            maxAlleles++;
        }
        for (int i = 2; i < maxAlleles; i++) {
            String alleleHeading = "Allele " + (i + 1);
            columnNames.add(alleleHeading);
            columnNames.add(alleleHeading + " Gametes");
            columnNames.add(alleleHeading + " Proportion");
            columnNames.add(alleleHeading + " Frequency");
        }

        columnNames.addAll(Arrays.asList(lastColumnNames));

        int numSites = alignment.getSiteCount();
        int numTaxa = alignment.getSequenceCount();
        Object[][] data = new Object[numSites][columnNames.size()];
        int totalGametes = numTaxa * 2;

        String[] snpIDs = alignment.getSNPIDs();
        boolean hasSnpIDs = ((snpIDs != null) && (snpIDs.length != 0));

        int[] physicalPositions = alignment.getPhysicalPositions();
        boolean hasPhysicalPositions = ((physicalPositions != null) && (physicalPositions.length != 0));

        for (int i = 0; i < numSites; i++) {

            int totalNotMissing = alignment.getTotalNotMissing(i);
            int totalGametesNotMissing = alignment.getTotalGametesNotMissing(i);
            int count = 0;

            data[i][count++] = i;
            if (hasSnpIDs) {
                data[i][count++] = snpIDs[i];
            } else {
                count++;
            }
            if (hasPhysicalPositions) {
                data[i][count++] = physicalPositions[i];
            } else {
                count++;
            }
            data[i][count++] = numTaxa;

            int[][] alleles = alignment.getAllelesSortedByFrequency(i);
            int numAlleles = alleles[0].length;

            for (int a = 0; a < numAlleles; a++) {
                data[i][count++] = alignment.getBaseAsString(i, (byte) alleles[0][a]);
                data[i][count++] = alleles[1][a];
                data[i][count++] = (double) alleles[1][a] / (double) totalGametes;
                data[i][count++] = (double) alleles[1][a] / (double) totalGametesNotMissing;
            }

            for (int b = 0; b < (maxAlleles - numAlleles); b++) {
                data[i][count++] = NA;
                data[i][count++] = ZERO_INT;
                data[i][count++] = ZERO_DOUBLE;
                data[i][count++] = ZERO_DOUBLE;
            }

            int totalGametesMissing = totalGametes - totalGametesNotMissing;
            myNumGametesMissing = myNumGametesMissing + (long) totalGametesMissing;
            data[i][count++] = totalGametesMissing;
            data[i][count++] = (double) totalGametesMissing / (double) totalGametes;

            int numHeterozygous = alignment.getHeterozygousCount(i);
            myNumHeterozygous = myNumHeterozygous + (long) numHeterozygous;
            data[i][count++] = numHeterozygous;
            data[i][count++] = (double) numHeterozygous / (double) totalNotMissing;

            data[i][count++] = "TBD";
            data[i][count++] = "TBD";

        }

        String[] columnNameStrings = new String[columnNames.size()];
        columnNames.toArray(columnNameStrings);
        return new SimpleTableReport("Site Summary", columnNameStrings, data);

    }

    private SimpleTableReport getTaxaSummary(Alignment alignment) {

        Object[] columnNames = new String[]{"Taxa", "Taxa Name", "Number of Sites", "Gametes Missing", "Proportion Missing",
            "Number Heterozygous", "Proportion Heterozygous", "Inbreeding Coefficient",
            "Inbreeding Coefficient Scaled by Missing"};
        int numSites = alignment.getSiteCount();
        int numTaxa = alignment.getSequenceCount();
        Object[][] data = new Object[numTaxa][columnNames.length];

        int totalGametes = numSites * 2;
        for (int i = 0; i < numTaxa; i++) {

            int totalGametesNotMissing = alignment.getTotalGametesNotMissingForTaxon(i);
            int totalGametesMissing = totalGametes - totalGametesNotMissing;
            int numHeterozygous = alignment.getHeterozygousCountForTaxon(i);
            int totalSitesNotMissing = alignment.getTotalNotMissingForTaxon(i);

            int count = 0;
            data[i][count++] = i;
            data[i][count++] = alignment.getIdGroup().getIdentifier(i).getFullName();
            data[i][count++] = numSites;
            data[i][count++] = totalGametesMissing;
            data[i][count++] = (double) totalGametesMissing / (double) totalGametes;
            data[i][count++] = numHeterozygous;
            data[i][count++] = (double) numHeterozygous / (double) totalSitesNotMissing;
            data[i][count++] = "Inbreeding Coefficient";
            data[i][count++] = "ICSBM";
        }

        return new SimpleTableReport("Taxa Summary", columnNames, data);

    }

    public void setCaculateOverview(boolean calculate) {
        myIsOverview = calculate;
    }

    public boolean getCalculateOverview() {
        return myIsOverview;
    }

    public void setCalculateSiteSummary(boolean calculate) {
        myIsSiteSummary = calculate;
    }

    public boolean getCalculateSiteSummary() {
        return myIsSiteSummary;
    }

    public void setCalculateTaxaSummary(boolean calculate) {
        myIsTaxaSummary = calculate;
    }

    public boolean getCalculateTaxaSummary() {
        return myIsTaxaSummary;
    }

    public ImageIcon getIcon() {
        URL imageURL = GenotypeSummaryPlugin.class.getResource("images/summary.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    public String getButtonName() {
        return "Geno Summary";
    }

    public String getToolTipText() {
        return "Genotype Summary";
    }

    class GenotypeSummaryPluginDialog extends JDialog {

        private JTabbedPane myTabbedPane = new JTabbedPane();
        private JRadioButton myOverview = new JRadioButton();
        private JRadioButton mySiteSummary = new JRadioButton();
        private JRadioButton myTaxaSummary = new JRadioButton();
        private boolean myIsCancel = true;

        public GenotypeSummaryPluginDialog() {
            super((Frame) null, "Genotype Summary", true);


            JButton okButton = new JButton();
            okButton.setActionCommand("Ok");
            okButton.setText("Ok");
            okButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    myIsCancel = false;
                    setVisible(false);
                }
            });
            JButton closeButton = new JButton();
            closeButton.setText("Close");
            closeButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    myIsCancel = true;
                    setVisible(false);
                }
            });

            JPanel panel = new JPanel();
            panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));

            //Radio Buttons
            myOverview.setText("Genotype Overview");
            myOverview.setSelected(myIsOverview);
            mySiteSummary.setText("Site Summary");
            mySiteSummary.setSelected(myIsSiteSummary);
            myTaxaSummary.setText("Taxa Summary");
            myTaxaSummary.setSelected(myIsTaxaSummary);

            panel.add(myOverview);
            panel.add(mySiteSummary);
            panel.add(myTaxaSummary);

            myTabbedPane.add(panel, "Genotype Summary");

            JPanel pnlButtons = new JPanel();
            pnlButtons.setLayout(new FlowLayout());
            pnlButtons.add(okButton);
            pnlButtons.add(closeButton);
            getContentPane().add(myTabbedPane, BorderLayout.CENTER);
            getContentPane().add(pnlButtons, BorderLayout.SOUTH);

            pack();
        }

        public boolean isCancel() {
            return myIsCancel;
        }

        public boolean isOverview() {
            return myOverview.isSelected();
        }

        public boolean isSiteSummary() {
            return mySiteSummary.isSelected();
        }

        public boolean isTaxaSummary() {
            return myTaxaSummary.isSelected();
        }
    }
}
