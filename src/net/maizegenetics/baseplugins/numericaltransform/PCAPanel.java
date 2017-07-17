package net.maizegenetics.baseplugins.numericaltransform;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

import net.maizegenetics.pal.alignment.Phenotype;
import net.maizegenetics.pal.alignment.SimplePhenotype;
import net.maizegenetics.pal.alignment.Trait;
import net.maizegenetics.pal.report.SimpleTableReport;
import net.maizegenetics.pal.report.TableReport;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.stats.PCA.PCA;
import net.maizegenetics.stats.PCA.PrincipalComponents;

import javax.swing.*;
import java.awt.*;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;

import java.io.StringWriter;
import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;

import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: Karin Holmberg and Michael Oak - with Dallas Kroon
 * Date: Jun 12, 2006
 * Time: 4:39:55 PM
 * To change this template use File | Settings | File Templates.
 */
public class PCAPanel extends JPanel {

    private Datum theDatum;
    // PCA Panel Components
    private JPanel pnlPCA = new JPanel();
    private String DEFAULT_EIGEN = "0";
    private String DEFAULT_CUM = "0.10";
    private String DEFAULT_COMP = "0";
    private JTextField tfdEigenvalue = new JTextField(5);
    private JTextField tfdCumulative = new JTextField(5);
    private JTextField tfdNumComp = new JTextField(5);
    private GridBagLayout gridBagLayout1 = new GridBagLayout();
    //  private JCheckBox               chbxMissingData             = new JCheckBox();
    private JSpinner spnComp;

    //    private int                     colSelected[];
    private StringBuffer procedure;
    private StringBuffer procedureReport;
    //Button Group
    private ButtonGroup bgrpPCA = new ButtonGroup();
    private ButtonGroup bgrpVar = new ButtonGroup();
    //    private TASSELMainFrame        parent;
    private Phenotype aCharacterAlignment;
    private String[] columnNames;
    // Filter Options
    private JRadioButton rdoEigen = new JRadioButton();
    private JRadioButton rdoCum = new JRadioButton();
    private JRadioButton rdoNum = new JRadioButton();
    private JRadioButton rdoCor = new JRadioButton();
    private JRadioButton rdoCovar = new JRadioButton();
    private JLabel lblMethod = new JLabel("Method");
    private JLabel lblOutput = new JLabel("Output");
    private double eigenvalue;
    private double userMinPropOfVar;
    private Frame myParentFrame = null;

    public PCAPanel(Datum theDatum, Frame parentFrame) throws Exception {
        myParentFrame = parentFrame;
        try {
            if (!(theDatum.getData() instanceof Phenotype)) {
                throw new Exception("Must be Character Alignment");
            }
            this.theDatum = theDatum;
            this.aCharacterAlignment = (Phenotype) theDatum.getData();
            Object[] trColumnNames = aCharacterAlignment.getTableColumnNames();
            String[] tempColNames = new String[trColumnNames.length];
            columnNames = new String[trColumnNames.length - 1];
            // start at 1 as the 0 column holds taxa names
            for (int i = 1; i < trColumnNames.length; i++) {
                if (trColumnNames instanceof String[]) {
                    tempColNames = (String[]) trColumnNames;
                }
                columnNames = new String[tempColNames.length - 1];
                //imputedValueCount = new int[tempColNames.length - 1];
                System.arraycopy(tempColNames, 1, columnNames, 0, columnNames.length);
            }
            //DEFAULT_CUM = "" + nNum(1.0 / (double) columnNames.length);
            DEFAULT_CUM = "" + nNum(1.0 / (double) aCharacterAlignment.getNumberOfTraits());
            jbInit();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void jbInit() throws Exception {
        // titledBorder1 = new TitledBorder(BorderFactory.createEtchedBorder(Color.white, new Color(178, 178, 178)), "Columns To Transform");

        pnlPCA.setLayout(gridBagLayout1);
        //Radio Buttons
        rdoEigen.setText("Eigenvalue \u2265");
        rdoCum.setText("Var Prop % \u2265");
        rdoNum.setText("Components = ");
        rdoCor.setText("Correlation");
        rdoCovar.setText("Covariance");
        tfdEigenvalue.setText(DEFAULT_EIGEN);
        tfdCumulative.setText(DEFAULT_CUM);
        tfdNumComp.setText(DEFAULT_COMP);
        //Button Group: Allows only one radiobutton to be selected at one time, default is to select Eigenvalue radiobutton
        bgrpPCA.add(rdoEigen);
        bgrpPCA.add(rdoCum);
        bgrpPCA.add(rdoNum);
        bgrpPCA.setSelected(rdoEigen.getModel(), true);

        bgrpVar.add(rdoCor);
        bgrpVar.add(rdoCovar);
        bgrpVar.setSelected(rdoCor.getModel(), true);

        SpinnerNumberModel aSpinnerNumberModel = new SpinnerNumberModel(columnNames.length, 1, columnNames.length, 1);
        spnComp = new JSpinner(aSpinnerNumberModel);
        spnComp.setMinimumSize(new Dimension(100, 100));
        spnComp.setOpaque(true);
        //Panel Setup
        pnlPCA.add(lblMethod, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0, GridBagConstraints.NORTH, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 7));
        pnlPCA.add(rdoCor, new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0, GridBagConstraints.NORTH, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 7));
        pnlPCA.add(rdoCovar, new GridBagConstraints(1, 2, 1, 1, 0.0, 0.0, GridBagConstraints.NORTH, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 7));
        pnlPCA.add(lblOutput, new GridBagConstraints(1, 3, 1, 1, 0.0, 0.0, GridBagConstraints.NORTH, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 7));

        pnlPCA.add(rdoEigen, new GridBagConstraints(1, 4, 1, 1, 0.0, 0.0, GridBagConstraints.NORTH, GridBagConstraints.NONE, new Insets(34, -10, 0, 50), 0, 7));
        pnlPCA.add(rdoCum, new GridBagConstraints(1, 5, 1, 1, 0.0, 0.0, GridBagConstraints.NORTH, GridBagConstraints.NONE, new Insets(0, 10, 0, 50), 7, 7));
        pnlPCA.add(rdoNum, new GridBagConstraints(1, 6, 1, 1, 0.0, 0.0, GridBagConstraints.NORTH, GridBagConstraints.NONE, new Insets(0, 0, 0, 50), 0, 7));
        pnlPCA.add(tfdEigenvalue, new GridBagConstraints(2, 4, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(40, -50, 5, 10), 1, 1));
        pnlPCA.add(tfdCumulative, new GridBagConstraints(2, 5, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(6, -50, 3, 10), 1, 1));
        pnlPCA.setBorder(BorderFactory.createEtchedBorder());
        pnlPCA.add(spnComp, new GridBagConstraints(2, 6, 15, 15, 1, 1, GridBagConstraints.NORTH, GridBagConstraints.NONE, new Insets(5, -50, 0, 10), 1, 1));
        this.add(pnlPCA);

        //Eigenvalue Errors
        eigenvalue = new BigDecimal(DEFAULT_EIGEN).doubleValue();
        tfdEigenvalue.addFocusListener(new FocusAdapter() {

            public void focusLost(FocusEvent e) {
                String minEigenvalue = tfdEigenvalue.getText();
                boolean failed = true;
                BigDecimal bd = null;
                try {
                    bd = new BigDecimal(minEigenvalue, new MathContext(2, RoundingMode.HALF_UP));
                    failed = false;
                } catch (NumberFormatException nfe) {
                }
                if (failed) {
                    JOptionPane.showMessageDialog(myParentFrame, "Please enter numbers only.");
                    tfdEigenvalue.setText(DEFAULT_EIGEN);
                } else {
                    tfdEigenvalue.setText(bd.toString());
                    eigenvalue = bd.doubleValue();
                    if (eigenvalue < 0) {
                        JOptionPane.showMessageDialog(myParentFrame, "Eigenvalue must be greater than or equal to 0.");
                        tfdEigenvalue.setText(DEFAULT_EIGEN);
                    }
                }
            }
        });

        //Cumulative Errors
        userMinPropOfVar = new BigDecimal(DEFAULT_CUM).doubleValue();
        tfdCumulative.addFocusListener(new FocusAdapter() {

            public void focusLost(FocusEvent e) {
                String minCumulative = tfdCumulative.getText();
                boolean failed = true;
                BigDecimal bd = null;
                try {
                    bd = new BigDecimal(minCumulative, new MathContext(2, RoundingMode.HALF_UP));
                    failed = false;
                } catch (NumberFormatException nfe) {
                }
                if (failed) {
                    JOptionPane.showMessageDialog(myParentFrame, "Please enter numbers only.");
                    tfdCumulative.setText(DEFAULT_CUM);
                } else {
                    tfdCumulative.setText(bd.toString());
                    userMinPropOfVar = bd.doubleValue();
                    if (userMinPropOfVar > 1 || userMinPropOfVar <= 0) {
                        JOptionPane.showMessageDialog(myParentFrame, "Cumulative proportion must be between 0 and 1.");
                        tfdCumulative.setText(DEFAULT_CUM);
                    }
                }
            }
        });

    }
    //Create Dataset Button

    public java.util.List<Datum> createPCAData(JTable tblTraits) {
        //if (selectedTab == PCATabIndex) {    ]
        int[] colSelected = tblTraits.getSelectedRows();
        if (colSelected == null || colSelected.length == 0) {
            JOptionPane.showMessageDialog(myParentFrame, "Please select a column to test.");
            return null;
        }
        PCA thePCA = null;
        try {
            thePCA = doPCAAnalysis(colSelected, this.aCharacterAlignment);
        } catch (Exception e) {
            JOptionPane.showMessageDialog(myParentFrame, "There is bad data in the PCA.  Impute missing data or remove non-variable columns.");
            return null;
        }
        if (thePCA != null) {
            return makeOutputFiles(thePCA);
        }
        return null;
    }

    //Create PCA Data
    private PCA doPCAAnalysis(int[] colSelected, Phenotype ca) {
        double[][] tempData = null;
        tempData = Conversion.parseColumnData(ca, colSelected);
        if (tempData == null) {
            return null;
        }
        DoubleMatrix2D values = new DenseDoubleMatrix2D(tempData);
//        values.assign(tempData);
        
        //Procedure Reports
        procedure = new StringBuffer("(");
        procedureReport = new StringBuffer("\nProcedure(s) applied:");
        // Filter by Eigenvalue
        procedure.append(" eigenvalues");
        procedureReport.append("\n\tPCA Eigenvalues: ");
        return new PCA(values, false, rdoCor.getModel().isSelected());

    }  //end create PCA

    private java.util.List<Datum> makeOutputFiles(PCA thePCAAnalysis) {
        java.util.List<Datum> resultList = new ArrayList<Datum>();
        DoubleMatrix2D principalComponents = thePCAAnalysis.getPC().viewDice();
        DoubleMatrix1D eigenvalues = thePCAAnalysis.getEigenValues();
        DoubleMatrix2D eigenvectors = thePCAAnalysis.getEigenVectors();
        double[] eigenvalueArray = eigenvalues.toArray();
        int numberOfPC = eigenvalueArray.length;
        int filterNumberOfPC = 0;
        int taxaNumber = principalComponents.rows();
        
        //Determine the number of columns and rows necessary in new matrix filtered by entered eigenvalue
        //int filterNumberOfPC = 0;
        String compString = spnComp.getValue().toString();
        int maxcomp = Integer.parseInt(compString);
        double rsum = 0;
        for (int n = 0; n < numberOfPC; n++) {
            rsum += eigenvalueArray[n];
        }
        for (int n = 0; n < numberOfPC; n++) {
            if (((rdoEigen.isSelected()) && (eigenvalueArray[n] >= eigenvalue)) ||
                    (((rdoNum.isSelected()) && (n < maxcomp))) ||
                    (((rdoCum.isSelected()) && ((eigenvalueArray[n] / rsum) > userMinPropOfVar)))) {
                filterNumberOfPC++;
            }
        }
        if (filterNumberOfPC == 0) {
            filterNumberOfPC = 1;
        }

        //Generate eigenFilter
        double[][] eigenFilter = principalComponents.viewPart(0, 0, taxaNumber, filterNumberOfPC).toArray();

        //Number Principal Components
        ArrayList<Trait> PCNames = new ArrayList<Trait>();
        for (int n = 1; n < filterNumberOfPC + 1; n++) {
            PCNames.add(new Trait("PC " + n, false, Trait.TYPE_COVARIATE));
        }
        
        // Creates Eigenvectors to be displayed in Data Tree
        double[][] eigenVectorsArray = eigenvectors.viewPart(0, 0, numberOfPC, filterNumberOfPC).toArray();

        if (eigenvalue >= 0) {
            //SimplePhenotype sca = new SimplePhenotype(aCharacterAlignment, eigenFilterT, PCNames);
            SimplePhenotype sca = new SimplePhenotype(aCharacterAlignment.getTaxa(), PCNames, eigenFilter);
            TableReport sca2 = createEigenVectorReport(eigenVectorsArray);
            TableReport sca3 = createCumulativeVarianceReport(eigenvalueArray);

            StringWriter sw = new StringWriter();
            String reportPanelMsg = new StringBuffer() + (sw.toString()) + ("Converted Column:\n") +
                    (procedureReport.toString()).toString();

            resultList.add(new Datum("PC for " + theDatum.getName(), sca, reportPanelMsg));
            resultList.add(new Datum("Eigenvectors for  " + theDatum.getName(), sca2, reportPanelMsg));
            resultList.add(new Datum("Eigenvalues for  " + theDatum.getName(), sca3, reportPanelMsg));

        }
        return resultList;
    }

    TableReport createEigenVectorReport(double[][] eigenVectorsArray) {
        // Creates column labels for eigenvectors
        Object[] evColNames = new String[eigenVectorsArray[0].length + 1];
        // EV = new String[eigenVectorsArray[0].length] ;
        evColNames[0] = "Trait.Env";
        for (int n = 0; n < eigenVectorsArray[0].length; n++) {
            evColNames[n + 1] = "EigenVector" + (n + 1);
        }
        Object[][] evOutput = new Object[eigenVectorsArray.length][eigenVectorsArray[0].length + 1];
        for (int i = 0; i < eigenVectorsArray.length; i++) {
            evOutput[i][0] = aCharacterAlignment.getTrait(i).getName() + "." + aCharacterAlignment.getTrait(i).getProperty(Trait.FACTOR_ENV);
            for (int n = 0; n < eigenVectorsArray[0].length; n++) {
                evOutput[i][n + 1] = nNum(eigenVectorsArray[i][n]);
            }
        }
        SimpleTableReport str = new SimpleTableReport("EigenVectors", evColNames, evOutput);
        return str;
    }

    TableReport createCumulativeVarianceReport(double[] eigenvalueArray) {
        // Creates Individual Percent Array
        double[] propVariancePC = new double[eigenvalueArray.length];
        double[] cumPropVariancePC = new double[propVariancePC.length];

        double sum = 0;
        for (int n = 0; n < eigenvalueArray.length; n++) {
            sum += eigenvalueArray[n];
        }
        double[] rpropVariancePC = new double[eigenvalueArray.length];
        for (int n = 0; n < eigenvalueArray.length; n++) {
            propVariancePC[n] = eigenvalueArray[n] / sum;
        }
        // Creates Labels for Percents Array
        Object[] cumlativeColNames = new String[4];
        cumlativeColNames[0] = "PC";
        cumlativeColNames[1] = "Eigenvalues";
        cumlativeColNames[2] = "Individual Proportion";
        cumlativeColNames[3] = "Cumulative Proportion";
        // Creates Cumulative Percent Array
        cumPropVariancePC[0] = propVariancePC[0];
        for (int n = 1; n < eigenvalueArray.length; n++) {
            cumPropVariancePC[n] = propVariancePC[n] + cumPropVariancePC[n - 1];
        }
        // Creates an Array with Individual and Cumulative Percents to be displayed
        Object[][] propOutput = new Object[eigenvalueArray.length][4];
        for (int n = 0; n < eigenvalueArray.length; n++) {
            //propOutput[n][0] = "PC" + (n + 1);
            propOutput[n][0] =  (n + 1);
            propOutput[n][1] = nNum(eigenvalueArray[n]);
            propOutput[n][2] = nNum(propVariancePC[n]);
            propOutput[n][3] = nNum(cumPropVariancePC[n]);
        }
        SimpleTableReport sca3 = new SimpleTableReport("Prop of Variance", cumlativeColNames, propOutput);
        return sca3;
    }

    private static String nNum(double d) {
        return (new BigDecimal(d, new MathContext(5))).toString();
    }

    private static double nNumD(double d) {
        return (new BigDecimal(d, new MathContext(5))).doubleValue();
    }
}


