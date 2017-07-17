package net.maizegenetics.baseplugins.numericaltransform;

import net.maizegenetics.pal.alignment.Phenotype;
import net.maizegenetics.pal.alignment.SimplePhenotype;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.alignment.Trait;
import net.maizegenetics.plugindef.Datum;

import javax.swing.*;
import java.awt.*;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: ed
 * Date: Oct 5, 2006
 * Time: 4:00:46 PM
 * To change this template use File | Settings | File Templates.
 */
public class ImputePanel extends JPanel {

    private Datum theDatum;
    private static final String DEFAULT_MIN_REQ_DATA = "0.80";
    //    private static final int        SCALE                       = 5;
    private double minRequiredData;
    //   private DataTreePanel           dtp;
    private Phenotype aCharacterAlignment;
    private JTextField tfdMinPercentRequired = new JTextField();
    private JSpinner spnK;
    private JRadioButton rdoUnweightedAverage = new JRadioButton("Unweighted Average");
    private JRadioButton rdoWeightedAverage = new JRadioButton("Weighted Average");
    private JRadioButton rdoManhattenDistance = new JRadioButton("Manhatten Distance");
    private JRadioButton rdoEuclidDistance = new JRadioButton("Euclid Distance");
    private ButtonGroup bgrpAverage = new ButtonGroup();
    private ButtonGroup bgrpDistance = new ButtonGroup();
    private JCheckBox cbxUseStandardizedData = new JCheckBox();

    public ImputePanel(Datum theDatum) {
        try {
            if (!(theDatum.getData() instanceof Phenotype)) {
                throw new Exception("Must be Character Alignment");
            }
            this.theDatum = theDatum;
            this.aCharacterAlignment = (Phenotype) theDatum.getData();
            jbInit();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void jbInit() throws Exception {
        // titledBorder1 = new TitledBorder(BorderFactory.createEtchedBorder(Color.white, new Color(178, 178, 178)), "Column(s) To Transform");
        JPanel fillerPanel = new JPanel();
        this.setLayout(new BorderLayout());
        JPanel optPanel = new JPanel(new GridLayout(6, 1, 5, 5));
        optPanel.add(rdoManhattenDistance);
        optPanel.add(rdoEuclidDistance);
        bgrpDistance.add(rdoManhattenDistance);
        bgrpDistance.add(rdoEuclidDistance);
        rdoManhattenDistance.setSelected(true);
        optPanel.add(rdoUnweightedAverage);
        optPanel.add(rdoWeightedAverage);
        bgrpAverage.add(rdoUnweightedAverage);
        bgrpAverage.add(rdoWeightedAverage);
        rdoUnweightedAverage.setSelected(true);

        // setup JSpinner to have correct values
        JPanel spnPanel = new JPanel(new BorderLayout());
        SpinnerNumberModel aSpinnerNumberModel = new SpinnerNumberModel(1, 1, aCharacterAlignment.getData().length, 1);
        spnK = new JSpinner(aSpinnerNumberModel);
        spnPanel.add(new JLabel("Number of Neighbors (K):"), BorderLayout.WEST);
        spnPanel.add(spnK, BorderLayout.CENTER);
        spnK.setOpaque(true);
        spnK.setValue(3);
        optPanel.add(spnPanel);

        JPanel minPPanel = new JPanel(new BorderLayout());
        tfdMinPercentRequired.setText(DEFAULT_MIN_REQ_DATA);
        minPPanel.add(new JLabel("Min. Freq. of Row Data:"), BorderLayout.WEST);
        minPPanel.add(tfdMinPercentRequired, BorderLayout.CENTER);
        minRequiredData = new Double(DEFAULT_MIN_REQ_DATA).doubleValue();
        optPanel.add(minPPanel);

        //       cbxUseStandardizedData.setText("Use Standardized Data");
        //       optPanel.add(cbxUseStandardizedData);

        this.add(optPanel, BorderLayout.NORTH);
        this.add(fillerPanel, BorderLayout.CENTER);
        tfdMinPercentRequired.addFocusListener(new FocusAdapter() {

            /**
             * Invoked when a component loses the keyboard focus.
             */
            public void focusLost(FocusEvent e) {
                boolean failed = true;
                String text = tfdMinPercentRequired.getText();
                BigDecimal bd = null;
                try {
                    bd = new BigDecimal(text, new MathContext(2, RoundingMode.HALF_UP));
                    failed = false;  // this statement will not be reached if Exception is thrown
                } catch (NumberFormatException nfe) {
                    // do nothing
                }
                if (failed) {
                    //              JOptionPane.showMessageDialog(thisDialog, "Please enter numbers only.");
                    tfdMinPercentRequired.setText(DEFAULT_MIN_REQ_DATA);
                } else {
                    tfdMinPercentRequired.setText(bd.toString());
                    minRequiredData = bd.doubleValue();
                }
            }
        });
    }

    public Datum createImputedData(JTable tblTraits) {
        int[] colsSelected = null;       // set of columns to be used to calculate distance (should be correlated columns)
        colsSelected = tblTraits.getSelectedRows();
        int colCount = colsSelected.length;
        int includedCount = 0;
        //find all the rows with enough data to keep
        int[] includedRowTemp = new int[aCharacterAlignment.getNumberOfTaxa()];
        for (int i = 0; i < aCharacterAlignment.getNumberOfTaxa(); i++) {
            double goodData = 0;
            for (int j = 0; j < colCount; j++) {
                if (!Double.isNaN(aCharacterAlignment.getData(i, colsSelected[j]))) {
                    goodData++;
                }
            }
            goodData = goodData / colCount;
            if (goodData >= minRequiredData) {
                includedRowTemp[includedCount++] = i;
            }
        }
        //rebuild the data set
        Identifier[] newIDs = new Identifier[includedCount];

        int traitCount = colsSelected.length;
        java.util.List<Trait> newtraits = new ArrayList<Trait>();
        for (int t = 0; t < traitCount; t++) {
            newtraits.add(Trait.getInstance(aCharacterAlignment.getTrait(colsSelected[t])));
        }

        double[][] tempData = new double[includedCount][colsSelected.length];
        for (int i = 0; i < includedCount; i++) {
            for (int j = 0; j < colCount; j++) {
                newIDs[i] = aCharacterAlignment.getTaxa().getIdentifier(includedRowTemp[i]);
                tempData[i][j] = aCharacterAlignment.getData(includedRowTemp[i], colsSelected[j]);
            }
        }
        //      for(int j = 0; j < colCount; j++){
        //              newTraits[j]=aCharacterAlignment.getTraitName(colsSelected[j]);
        //              newEnvs[j]=aCharacterAlignment.getEnvironmentName(colsSelected[j]);
        //          }
        int kNeighbors = Integer.parseInt(spnK.getValue().toString());
        double[][] theImputedData = KNN.impute(tempData, kNeighbors, rdoManhattenDistance.isSelected(), rdoUnweightedAverage.isSelected());
        //SimplePhenotype sca = new SimplePhenotype(new SimpleIdGroup(newIDs), theImputedData, aCharacterAlignment.getFactorNameCopy(), newtraits);
        SimplePhenotype sca = new SimplePhenotype(new SimpleIdGroup(newIDs), newtraits, theImputedData);
        StringWriter sw = new StringWriter();
        //sca.report(new PrintWriter(sw));
        String theComment = sw.toString() + "\nImputed Phenotypic Values." + "\nTaxa with insufficient data: " + (aCharacterAlignment.getNumberOfTaxa() - sca.getNumberOfTaxa()) + "\nK = " + kNeighbors + minRequiredData + "% cutoff):\n";
        String theName = theDatum.getName() + "_" + colCount + "_imputed";
        Datum result = new Datum(theName, sca, theComment);
        return result;
    }
}
