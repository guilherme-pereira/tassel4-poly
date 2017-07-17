package net.maizegenetics.baseplugins.numericaltransform;

import net.maizegenetics.pal.alignment.Phenotype;
import net.maizegenetics.pal.alignment.SimplePhenotype;
import net.maizegenetics.plugindef.Datum;

import javax.swing.*;
import java.awt.*;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.LinkedHashMap;

/**
 * Created by IntelliJ IDEA.
 * User: dallas/ed
 * Date: Oct 5, 2006
 * Time: 4:00:46 PM
 * To change this template use File | Settings | File Templates.
 */
public class TransformationPanel extends JPanel {

    private Datum theDatum;
    private static final String BASE_10 = "Base 10";
    private static final String BASE_2 = "Base 2";
    private static final String E_CONSTANT = "Natural";

    // private static final String     DEFAULT_MIN_REQ_DATA        = "80.00";
    private static final int SCALE = 5;

    // private double                  minRequiredData;
    private Phenotype aCharacterAlignment;
    // private String[]                columnNames;
    // private int[]                   imputedValueCount;          // count # of values imputed for the column; shares
    // normality panel components
    //   private JPanel                  pnlNormality                = new JPanel();
    private JRadioButton rdoPower = new JRadioButton();
    private JRadioButton rdoLog = new JRadioButton();
    private JCheckBox ckbxStandardize = new JCheckBox("Standardize");
    private JComboBox cbxPower = new JComboBox();
    private JComboBox cbxLog = new JComboBox();
    private LinkedHashMap powerMap = new LinkedHashMap();
    private String[] dataLog = {BASE_10, BASE_2, E_CONSTANT};
    private boolean rdoPowerSelected = true;
    private int colSelected[];        // column selected for transformation or
    // normality testing
    private double[] transformedColumnData[];
    private StringBuffer procedure;
    private StringBuffer procedureReport;
    private Frame myParentFrame = null;

    public TransformationPanel(Datum theDatum, Frame parentFrame) {
        myParentFrame = parentFrame;
        try {
            if (!(theDatum.getData() instanceof Phenotype)) {
                throw new Exception("Must be Character Alignment");
            }
            this.theDatum = theDatum;
            this.aCharacterAlignment = (Phenotype) theDatum.getData();
            powerMap.put("2", new Double(2));
            powerMap.put("3", new Double(3));
            powerMap.put("4", new Double(4));
            powerMap.put("10", new Double(10));
            powerMap.put("1/2", new Double(0.5));
            powerMap.put("1/3", new Double((double) 1 / 3));
            powerMap.put("1/4", new Double(0.25));
            powerMap.put("1/10", new Double(0.1));
            jbInit();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void jbInit() throws Exception {
        //      titledBorder1 = new TitledBorder(BorderFactory.createEtchedBorder(Color.white, new Color(178, 178, 178)), "Column(s) To Transform");
        JPanel fillerPanel = new JPanel();
        this.setLayout(new BorderLayout());
        JPanel optPanel = new JPanel(new GridLayout(3, 1, 5, 5));


        JPanel jpow = new JPanel(new BorderLayout());
        rdoPower.setText("Raise to Power");
        cbxPower = new JComboBox(powerMap.keySet().toArray());
        // cbxPower.setPreferredSize(new Dimension(50,30));
        jpow.add(rdoPower, BorderLayout.WEST);
        jpow.add(cbxPower, BorderLayout.CENTER);
        optPanel.add(jpow);

        JPanel jlog = new JPanel(new BorderLayout());
        rdoLog.setText("Take Log");
        rdoLog.setSelected(!rdoPowerSelected);
        cbxLog = new JComboBox(dataLog);
        // cbxLog.setPreferredSize(new Dimension(50,30));
        jlog.add(rdoLog, BorderLayout.WEST);
        jlog.add(cbxLog, BorderLayout.CENTER);
        optPanel.add(jlog);

        optPanel.add(ckbxStandardize);
        ButtonGroup rbg = new ButtonGroup();
        rbg.add(rdoPower);
        rbg.add(rdoLog);
        this.add(optPanel, BorderLayout.NORTH);
        this.add(fillerPanel, BorderLayout.CENTER);
    }

    public Datum createNormalizedData(JTable tblNormalityTraits) {
        // colSelected = lstColumnNames.getSelectedIndices();
        colSelected = tblNormalityTraits.getSelectedRows();

        double[][] tempData = null;

        // need to add one to account for the the taxa column at position 0
        //      for (int i = 0; i < colSelected.length; i++) {
        //          colSelected[i]++;
        //      }
        if (colSelected == null || colSelected.length == 0) {
            JOptionPane.showMessageDialog(myParentFrame, "Please select a column to test.");
        } else {
            tempData = Conversion.parseColumnData(aCharacterAlignment, colSelected);
        }
        if (tempData == null) {
            JOptionPane.showMessageDialog(myParentFrame, "There is no numeric data in the selected column.");
            return null;
        }
        procedure = new StringBuffer("(");
        procedureReport = new StringBuffer("\nProcedure(s) applied:");
        if (rdoPower.isSelected()) {
            double powerValue = 1;
            String power = (String) cbxPower.getSelectedItem();
            try {
                powerValue = ((Double) powerMap.get(power)).doubleValue();
            } catch (NumberFormatException nfe) {
                nfe.printStackTrace();
            }
            for (int i = 0; i < tempData.length; i++) {
                for (int j = 0; j < tempData[i].length; j++) {
                    if (!Double.isNaN(tempData[i][j])) {
                        BigDecimal bd = new BigDecimal(Math.pow(tempData[i][j], powerValue));
                        bd = bd.setScale(SCALE, RoundingMode.HALF_UP);
                        tempData[i][j] = bd.doubleValue();
                    }
                }
            }
            procedure.append(" power: " + powerValue);
            procedureReport.append("\n\tPower: " + powerValue);

        } else if (rdoLog.isSelected()) {
            for (int i = 0; i < tempData.length; i++) {
                for (int j = 0; j < tempData[j].length; j++) {
                    if (tempData[i][j] <= 0.000001) {
                        tempData[i][j] = Double.NaN;
                    }
                    if (!Double.isNaN(tempData[i][j])) {
                        if ((String) cbxLog.getSelectedItem() == E_CONSTANT) {
                            BigDecimal bd = new BigDecimal(Math.log(tempData[i][j]));
                            bd = bd.setScale(SCALE, RoundingMode.HALF_UP);
                            tempData[i][j] = bd.doubleValue();
                        } else if ((String) cbxLog.getSelectedItem() == BASE_2) {
                            BigDecimal bd = new BigDecimal(Math.log(tempData[i][j]) / Math.log(2));
                            bd = bd.setScale(SCALE, RoundingMode.HALF_UP);
                            tempData[i][j] = bd.doubleValue();
                        } else if ((String) cbxLog.getSelectedItem() == BASE_10) {
                            BigDecimal bd = new BigDecimal(Math.log(tempData[i][j]) / Math.log(10));
                            bd = bd.setScale(SCALE, RoundingMode.HALF_UP);
                            tempData[i][j] = bd.doubleValue();
                        }
                    }
                }
            }
            procedure.append(" log " + (String) cbxLog.getSelectedItem());
            procedureReport.append("\n\tLog " + (String) cbxLog.getSelectedItem());
        }
        if (ckbxStandardize.isSelected()) {
            tempData = Conversion.normalizeData(tempData);
            procedure.append(" normalized");
            procedureReport.append("\n\tNormalization");
        }
        transformedColumnData = tempData;

        // if there was nothing selected, do not put anything on the data tree
        if (procedure.toString().equals("(")) {
            JOptionPane.showMessageDialog(myParentFrame, "Please select a transformation to do.");
            return null;
        }

        SimplePhenotype sca = Conversion.reconstituteDataset(aCharacterAlignment, colSelected, transformedColumnData);

        StringWriter sw = new StringWriter();
        //sca.report(new PrintWriter(sw));

        String theComment = sw.toString() + "Converted Column:\n" //+ lstColumnNames.getSelectedValue()
                + procedureReport.toString();

        String theName = theDatum.getName() + ": " + colSelected.length + " col(s) " + procedure.append(")").toString();
        Datum result = new Datum(theName, sca, theComment);
        return result;
    }
}
