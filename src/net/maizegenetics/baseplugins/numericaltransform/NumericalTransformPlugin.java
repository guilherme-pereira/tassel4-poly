/*
 * NumericalTransformPlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.baseplugins.numericaltransform;

import net.maizegenetics.baseplugins.NumericalGenotypePlugin;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.Phenotype;
import net.maizegenetics.pal.alignment.SimplePhenotype;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import java.net.URL;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import net.maizegenetics.pal.alignment.Trait;

/**
 *
 * @author Ed Buckler
 */
public class NumericalTransformPlugin extends AbstractPlugin {
    //todo this is very friendly to running non-interactive currently

    /** Creates a new instance of NumericalTransformPlugin */
    public NumericalTransformPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {

        try {

            List<Datum> phenotypeList = input.getDataOfType(Phenotype.class);
            List<Datum> alignmentList = input.getDataOfType(Alignment.class);

            if (phenotypeList.size() + alignmentList.size() < 1) {
                JOptionPane.showMessageDialog(getParentFrame(), "Invalid selection. Please select genotype or phenotype data.");
                return null;
            }

            List<Datum> outputList = new ArrayList<Datum>();
            Iterator<Datum> itr = phenotypeList.iterator();
            while (itr.hasNext()) {
                Datum current = itr.next();
                List<Datum> td = processDatum(current, isInteractive());
                if (td != null) {
                    outputList.addAll(td);
                }
            }

            Iterator<Datum> itr2 = alignmentList.iterator();
            while (itr2.hasNext()) {
                Datum current = itr2.next();
                List<Datum> td = processDatum(current, isInteractive());
                if (td != null) {
                    outputList.addAll(td);
                }
            }

            if (outputList.isEmpty()) {
                return null;
            }

            DataSet output = new DataSet(outputList, this);
            // I am setting the firing class as the metadata
            // so that the control panel know where the event is coming from
            fireDataSetReturned(new PluginEvent(output, NumericalTransformPlugin.class));

            return output;

        } finally {
            fireProgress(100);
        }

    }

    private List<Datum> processDatum(Datum inDatum, boolean isInteractive) {

        List<Datum> resultList = null;

        if (isInteractive) {

            NumTransformDialog dialog = null;
            try {
                dialog = new NumTransformDialog(getParentFrame(), inDatum);
                dialog.setLocationRelativeTo(getParentFrame());
                dialog.setVisible(true);
                resultList = dialog.getResults();
            } catch (Exception e) {
                e.printStackTrace();
            } finally {
                dialog.dispose();
            }

        }

        return resultList;

    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        URL imageURL = NumericalTransformPlugin.class.getResource("Transform.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    public String getButtonName() {
        return "Transform";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Transform Numerical Dataset";
    }
}

/**
 * @author
 * @version 1.0
 */
class NumTransformDialog extends JDialog {

    private Datum myInput;
    private List<Datum> myResults = new ArrayList<Datum>();
    private Phenotype myCharacterAlignment;
    private Alignment myGenotype;
    private JTabbedPane myTabbedPane = new JTabbedPane();
    private int myTransformTabIndex = -1;
    private int myImputationTabIndex = -1;
    private int myPcaTabIndex = -1;
    private int myNgTabIndex = -1;
    private ImputePanel myPnlImputation;
    private TransformationPanel myPnlTranformation;
    private PCAPanel myPnlPCA;
    private NumericalGenotypePanel myPnlNumericalGenotype;
    private JTable myTblTraits;
    private Frame myParentFrame = null;

    public NumTransformDialog(Frame frame, Datum theDatum) throws Exception {

        super(frame, true);
        myInput = theDatum;
        myParentFrame = frame;

        if (theDatum.getData() instanceof Phenotype) {
            myCharacterAlignment = (Phenotype) theDatum.getData();
        } else if (theDatum.getData() instanceof Alignment) {
            myGenotype = (Alignment) theDatum.getData();
        } else {
            throw new Exception("Must be Phenotype or Genotype");
        }

        try {
            jbInit();
            pack();
        } catch (Exception e) {
            System.out.println("Error in NumTransform Dialog");
        }

    }

    private void jbInit() throws Exception {
        if (myCharacterAlignment != null) {
            myTblTraits = createTraitTable(myCharacterAlignment);
            myTblTraits.selectAll();
        }
        myTabbedPane.setDebugGraphicsOptions(0);
        JButton btnCreateAlignment = new JButton();
        btnCreateAlignment.setActionCommand("Create Dataset");
        btnCreateAlignment.setText("Create Dataset");
        btnCreateAlignment.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                btnCreateAlignment_actionPerformed(e);
            }
        });
        JButton btnClose = new JButton();
        btnClose.setText("Close");
        btnClose.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                btnClose_actionPerformed(e);
            }
        });

        if (myCharacterAlignment != null) {
            myPnlTranformation = new TransformationPanel(myInput, myParentFrame);
            myPnlImputation = new ImputePanel(myInput);
            myPnlPCA = new PCAPanel(myInput, myParentFrame);

            myTabbedPane.add(myPnlTranformation, "Trans");
            myTransformTabIndex = myTabbedPane.getTabCount() - 1;
            myTabbedPane.add(myPnlImputation, "Impute");
            myImputationTabIndex = myTabbedPane.getTabCount() - 1;
            myTabbedPane.add(myPnlPCA, "PCA");
            myPcaTabIndex = myTabbedPane.getTabCount() - 1;
        } else if (myGenotype != null) {
            myPnlNumericalGenotype = new NumericalGenotypePanel(myInput);
            myTabbedPane.add(myPnlNumericalGenotype, "Numerical Genotype");
            myNgTabIndex = myTabbedPane.getTabCount() - 1;
        }

        JScrollPane scpTraits = new JScrollPane(myTblTraits);
        JPanel pnlButtons = new JPanel();
        pnlButtons.setLayout(new FlowLayout());
        pnlButtons.add(btnCreateAlignment);
        pnlButtons.add(btnClose);
        this.getContentPane().add(scpTraits, BorderLayout.CENTER);
        this.getContentPane().add(myTabbedPane, BorderLayout.EAST);
        this.getContentPane().add(pnlButtons, BorderLayout.SOUTH);

    }

    private JTable createTraitTable(Phenotype ca) {
        String[] tableColumnNames = {"Column", "Percent Missing Data"};
        Object[] missingData = Conversion.getPercentMissingData(ca, 2);
        Object[][] tableData = new Object[ca.getNumberOfTraits()][2];
        for (int i = 0; i < ca.getNumberOfTraits(); i++) {
            tableData[i][0] = ca.getTrait(i).getName() + "." + ca.getTrait(i).getProperty(Trait.FACTOR_ENV);
            tableData[i][1] = missingData[i];
        }
        JTable tblAvailableColumns = new JTable(tableData, tableColumnNames);
        tblAvailableColumns.setEnabled(true);
        tblAvailableColumns.setRowSelectionAllowed(true);
        tblAvailableColumns.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
        return tblAvailableColumns;
    }

    private void btnCreateAlignment_actionPerformed(ActionEvent e) {

        Datum result = null;

        int selectedTab = myTabbedPane.getSelectedIndex();

        if (selectedTab == myTransformTabIndex) {
            result = myPnlTranformation.createNormalizedData(myTblTraits);
        } else if (selectedTab == myImputationTabIndex) {
            result = myPnlImputation.createImputedData(myTblTraits);
        } else if (selectedTab == myPcaTabIndex) {
            List<Datum> resultL = myPnlPCA.createPCAData(myTblTraits);
            if (resultL != null) {
                myResults.addAll(resultL);
            }
        } else if (selectedTab == myNgTabIndex) {
            if (myPnlNumericalGenotype.isCollapseSelected()) {
                SimplePhenotype temp = NumericalGenotypePlugin.collapseTransform(myGenotype);
                result = new Datum(myInput.getName() + "_Collasped", temp, null);
            } else if (myPnlNumericalGenotype.isSeparateSelected()) {
                SimplePhenotype temp = NumericalGenotypePlugin.separatedTransform(myGenotype);
                result = new Datum(myInput.getName() + "_Separated", temp, null);
            }

        }

        if (result != null) {
            myResults.add(result);
        }

        setVisible(false);

    }

    List<Datum> getResults() {
        return myResults;
    }

    private void btnClose_actionPerformed(ActionEvent e) {
        this.setVisible(false);
    }
}
