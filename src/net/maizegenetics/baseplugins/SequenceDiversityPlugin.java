/*
 * SequenceDiversityPlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.baseplugins;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.popgen.DiversityAnalyses;
import net.maizegenetics.pal.popgen.PolymorphismDistribution;
import net.maizegenetics.pal.report.SimpleTableReport;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;

import javax.swing.*;
import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.FocusEvent;

import java.net.URL;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;

/**
 *
 * @author Ed Buckler
 */
public class SequenceDiversityPlugin extends AbstractPlugin {

    Vector typeOfSitesToAnalyze = new Vector();
    boolean isSlidingWindowAnalysis = false;
    int startSite = 0;
    int endSite = 0;
    int windowSize = 500;
    int stepSize = 100;
    //todo implement Tajima and Fu & Li Selection Tests

    /** Creates a new instance of SequenceDiversityPlugin */
    public SequenceDiversityPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
        typeOfSitesToAnalyze.add(new Integer(Alignment.POSITION_TYPE_ALL_GROUP));
    }

    public DataSet performFunction(DataSet input) {

        try {
            List<Datum> alignInList = input.getDataOfType(Alignment.class);
            if (alignInList.size() < 1) {
                JOptionPane.showMessageDialog(getParentFrame(), "Invalid selection.  Please select sequence or marker alignment.");
                return null;
            }

            List result = new ArrayList();
            Iterator<Datum> itr = alignInList.iterator();
            while (itr.hasNext()) {
                DataSet tds = null;
                Datum current = itr.next();
                Alignment aa = (Alignment) current.getData();
                if (isInteractive()) {
                    DiversityDialog myDialog = new DiversityDialog(aa);
                    myDialog.setLocationRelativeTo(getParentFrame());
                    myDialog.setVisible(true);
                    if (myDialog.isCancel()) {
                        return null;
                    }
                    isSlidingWindowAnalysis = myDialog.isSlidingWindowAnalysis();
                    //typeOfSitesToAnalyze = myDialog.getTypeOfSitesToAnalyze();
                    startSite = myDialog.getStartSite();
                    endSite = myDialog.getEndSite();
                    windowSize = myDialog.getWindowSize();
                    stepSize = myDialog.getStepSize();
                } else {
                    if ((startSite + 1) > aa.getSiteCount()) {
                        startSite = 0;
                    }
                    if ((endSite < 1) || ((endSite + 1) > aa.getSiteCount())) {
                        endSite = aa.getSiteCount() - 1;
                    }
                    if ((windowSize + 1) > aa.getSiteCount()) {
                        windowSize = aa.getSiteCount() - 1;
                    }
                    if ((stepSize + 1) > aa.getSiteCount()) {
                        stepSize = aa.getSiteCount() - 1;
                    }
                }
                tds = processDatum(current);
                if (tds != null) {
                    result.add(tds);
                    fireDataSetReturned(new PluginEvent(tds, SequenceDiversityPlugin.class));
                }
            }

            return DataSet.getDataSet(result, this);
        } finally {
            fireProgress(100);
        }
    }

    public DataSet processDatum(Datum input) {
        Alignment aa = (Alignment) input.getData();
        PolymorphismDistribution pda = new PolymorphismDistribution();
        DiversityAnalyses theDA = new DiversityAnalyses(aa, typeOfSitesToAnalyze, isSlidingWindowAnalysis,
                startSite, endSite, windowSize, stepSize, pda);
        List<Datum> results = new ArrayList<Datum>();
        results.add(new Datum("PolyDist:" + input.getName(), new SimpleTableReport(pda), "Polymorphism Distribution"));
        results.add(new Datum("Diversity:" + input.getName(), new SimpleTableReport(theDA), "Diversity Analysis"));
        DataSet tds = new DataSet(results, this);
        return tds;
    }

    public Vector getTypeOfSitesToAnalyze() {
        return typeOfSitesToAnalyze;
    }

    public void setTypeOfSitesToAnalyze(Vector typeOfSitesToAnalyze) {
        this.typeOfSitesToAnalyze = typeOfSitesToAnalyze;
    }

    public boolean isSlidingWindowAnalysis() {
        return isSlidingWindowAnalysis;
    }

    public void setSlidingWindowAnalysis(boolean slidingWindowAnalysis) {
        isSlidingWindowAnalysis = slidingWindowAnalysis;
    }

    public int getStartSite() {
        return startSite;
    }

    public void setStartSite(int startSite) {
        this.startSite = startSite;
    }

    public int getEndSite() {
        return endSite;
    }

    public void setEndSite(int endSite) {
        this.endSite = endSite;
    }

    public int getWindowSize() {
        return windowSize;
    }

    public void setWindowSize(int windowSize) {
        this.windowSize = windowSize;
    }

    public int getStepSize() {
        return stepSize;
    }

    public void setStepSize(int stepSize) {
        this.stepSize = stepSize;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        URL imageURL = SequenceDiversityPlugin.class.getResource("images/Diversity.gif");
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
        return "Diversity";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Basic description of diversity";
    }
}

/**
 *This dialog collects informations on which diversity tests to carry out.
 *
 * @author Ed Buckler
 * @version 1.0
 */
class DiversityDialog extends JDialog {

    int start = 0, end, step = 100, window = 500;
    boolean runAnalysis = false;
    Alignment theAlignment;
    JPanel panel1 = new JPanel();
    JButton runButton = new JButton();
    JPanel jPanel1 = new JPanel();
    JTextField stepTextField = new JTextField();
    JTextField endTextField = new JTextField();
    JTextField startTextField = new JTextField();
    JLabel jLabel3 = new JLabel();
    JLabel jLabel2 = new JLabel();
    JLabel jLabel1 = new JLabel();
    JCheckBox slidingCheckBox = new JCheckBox();
    JCheckBox silentCheckBox = new JCheckBox();
    JCheckBox noncodingCheckBox = new JCheckBox();
    JCheckBox intronCheckBox = new JCheckBox();
    JCheckBox nontranscribedCheckBox = new JCheckBox();
    JCheckBox synonymousCheckBox = new JCheckBox();
    JCheckBox indelNonCodingCheckBox = new JCheckBox();
    JCheckBox codingCheckBox = new JCheckBox();
    JCheckBox codingIndelsCheckBox = new JCheckBox();
    JCheckBox overallCheckBox = new JCheckBox();
    JCheckBox allIndelCheckBox = new JCheckBox();
    JTextField windowTextField = new JTextField();
    JLabel jLabel4 = new JLabel();
    JCheckBox nonSynCheckBox = new JCheckBox();
    JCheckBox transcribedCheckBox = new JCheckBox();
    JPanel jPanel2 = new JPanel();
    JPanel jPanel3 = new JPanel();
    JButton closeButton = new JButton();
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    GridBagLayout gridBagLayout2 = new GridBagLayout();
    GridBagLayout gridBagLayout3 = new GridBagLayout();

    public DiversityDialog(Alignment aa) {
        super((Frame) null, "Diversity Surveys", true);
        theAlignment = aa;
        try {
            end = theAlignment.getSiteCount() - 1;
            jbInit();
            turnOffOptionsIfNotAnnotated();
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    private void jbInit() throws Exception {
        panel1.setLayout(gridBagLayout1);
        runButton.setText("Run");
        runButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                runButton_actionPerformed(e);
            }
        });
        closeButton.setText("Close");
        closeButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                closeButton_actionPerformed(e);
            }
        });
        stepTextField.setText("" + step);
        stepTextField.addFocusListener(new java.awt.event.FocusAdapter() {

            public void focusLost(FocusEvent e) {
                stepTextField_focusLost(e);
            }
        });
        endTextField.setText("" + end);
        endTextField.addFocusListener(new java.awt.event.FocusAdapter() {

            public void focusLost(FocusEvent e) {
                endTextField_focusLost(e);
            }
        });
        startTextField.setText("" + start);
        startTextField.addFocusListener(new java.awt.event.FocusAdapter() {

            public void focusLost(FocusEvent e) {
                startTextField_focusLost(e);
            }
        });
        jLabel3.setText("Step");
        jLabel2.setText("End Base");
        jLabel1.setText("Start Base");
        slidingCheckBox.setText("Sliding Window");
        jPanel1.setLayout(gridBagLayout2);
        jPanel1.setBorder(BorderFactory.createEtchedBorder());
        silentCheckBox.setText("Silent");
        silentCheckBox.setSelected(false);
        noncodingCheckBox.setText("Noncoding");
        intronCheckBox.setText("Intron");
        nontranscribedCheckBox.setText("Non-transcribed");
        synonymousCheckBox.setText("Synonymous");
        indelNonCodingCheckBox.setText("Noncoding Indels");
        codingCheckBox.setText("Coding");
        codingIndelsCheckBox.setText("Coding Indels");
        overallCheckBox.setText("Overall");
        overallCheckBox.setSelected(true);
        allIndelCheckBox.setText("Indels");
        windowTextField.setText("" + window);
        windowTextField.addFocusListener(new java.awt.event.FocusAdapter() {

            public void focusLost(FocusEvent e) {
                windowTextField_focusLost(e);
            }
        });
        jLabel4.setText("Window");
        nonSynCheckBox.setText("Nonsynonymous");
        transcribedCheckBox.setText("Transcribed");
        jPanel2.setPreferredSize(new Dimension(200, 300));
        jPanel2.setLayout(gridBagLayout3);
        closeButton.setText("Close");
        closeButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                closeButton_actionPerformed(e);
            }
        });
        panel1.setMinimumSize(new Dimension(200, 200));
        panel1.setPreferredSize(new Dimension(200, 250));
        panel1.add(silentCheckBox, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        panel1.add(noncodingCheckBox, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        panel1.add(intronCheckBox, new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        panel1.add(synonymousCheckBox, new GridBagConstraints(0, 3, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        panel1.add(nontranscribedCheckBox, new GridBagConstraints(0, 4, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        panel1.add(indelNonCodingCheckBox, new GridBagConstraints(0, 5, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        panel1.add(transcribedCheckBox, new GridBagConstraints(0, 6, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        panel1.add(nonSynCheckBox, new GridBagConstraints(0, 7, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        panel1.add(codingCheckBox, new GridBagConstraints(0, 8, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        panel1.add(codingIndelsCheckBox, new GridBagConstraints(0, 9, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        panel1.add(overallCheckBox, new GridBagConstraints(0, 10, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        panel1.add(allIndelCheckBox, new GridBagConstraints(0, 11, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));



        jPanel1.add(slidingCheckBox, new GridBagConstraints(0, 0, 2, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        jPanel1.add(stepTextField, new GridBagConstraints(1, 1, 1, 1, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
        jPanel1.add(jLabel3, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 31, 7, 0), 0, 0));
        jPanel1.add(jLabel4, new GridBagConstraints(0, 2, 2, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(12, 12, 8, 62), 20, 0));
        jPanel1.add(windowTextField, new GridBagConstraints(1, 2, 1, 1, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(8, 15, 8, 15), 0, 0));
        this.getContentPane().add(panel1, BorderLayout.WEST);
        jPanel2.add(jLabel1, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(58, 36, 0, 0), 0, 0));
        jPanel2.add(startTextField, new GridBagConstraints(1, 0, 1, 1, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(58, 10, 0, 50), 41, 0));
        jPanel2.add(jLabel2, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(11, 36, 0, 0), 0, 0));
        jPanel2.add(endTextField, new GridBagConstraints(1, 1, 1, 1, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(9, 10, 0, 49), 42, 0));
        this.getContentPane().add(jPanel3, BorderLayout.SOUTH);
        jPanel3.add(runButton, null);
        jPanel3.add(closeButton, null);
        this.getContentPane().add(jPanel2, BorderLayout.EAST);
        jPanel2.add(jPanel1, new GridBagConstraints(0, 2, 2, 1, 1.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(27, 43, 43, 14), 0, 4));
    }

    private void endTextField_focusLost(FocusEvent e) {
        try {
            end = Integer.parseInt(endTextField.getText());
            if ((end <= start) || (end > (theAlignment.getSiteCount() - 1))) {
                end = this.theAlignment.getSiteCount() - 1;
            }
        } catch (Exception ee) {
            end = this.theAlignment.getSiteCount() - 1;
        }
        endTextField.setText(end + "");
    }

    private void turnOffOptionsIfNotAnnotated() {
        boolean annotatedSites = false;
        //for (int i = 0; i < theAlignment.getSiteCount(); i++) {
        //    if (theAlignment.getPositionType(i) != 0) {
        //        annotatedSites = true;
        //        continue;
        //    }
        //}
        if (annotatedSites == false) {
            codingCheckBox.setEnabled(false);
            noncodingCheckBox.setEnabled(false);
            silentCheckBox.setEnabled(false);
            intronCheckBox.setEnabled(false);
            nontranscribedCheckBox.setEnabled(false);
            codingIndelsCheckBox.setEnabled(false);
            synonymousCheckBox.setEnabled(false);
            indelNonCodingCheckBox.setEnabled(false);
            nonSynCheckBox.setEnabled(false);
            transcribedCheckBox.setEnabled(false);
        }
    }

    void startTextField_focusLost(FocusEvent e) {
        try {
            start = Integer.parseInt(startTextField.getText());
            if ((end <= start) || (start < 0)) {
                start = 0;
            }
        } catch (Exception ee) {
            start = 0;
        }
        startTextField.setText(start + "");
    }

    /** Returns whether the run button was chosen*/
    public boolean isRunAnalysis() {
        return runAnalysis;
    }

    /** Returns whether sliding windows of diversity should be calculated */
    public boolean isSlidingWindowAnalysis() {
        return slidingCheckBox.isSelected();
    }

    /** Returns whether sliding windows of diversity should be calculated */
    /*
    public Vector getTypeOfSitesToAnalyze() {
        Vector grp = new Vector();
        if (overallCheckBox.isSelected()) {
            grp.add(new Integer(Alignment.POSITION_TYPE_ALL_GROUP));
        }
        if (codingCheckBox.isSelected()) {
            grp.add(new Integer(Alignment.POSITION_TYPE_CODING_GROUP));
        }
        if (noncodingCheckBox.isSelected()) {
            grp.add(new Integer(Alignment.POSITION_TYPE_NONCODING_GROUP));
        }
        if (silentCheckBox.isSelected()) {
            grp.add(new Integer(Alignment.POSITION_TYPE_SILENT_GROUP));
        }
        if (intronCheckBox.isSelected()) {
            grp.add(new Integer(Alignment.POSITION_TYPE_INTRON_GROUP));
        }
        if (nontranscribedCheckBox.isSelected()) {
            grp.add(new Integer(Alignment.POSITION_TYPE_NONTRANSSCRIBED_GROUP));
        }
        if (codingIndelsCheckBox.isSelected()) {
            grp.add(new Integer(Alignment.POSITION_TYPE_CODINGINDEL_GROUP));
        }
        if (synonymousCheckBox.isSelected()) {
            grp.add(new Integer(Alignment.POSITION_TYPE_SYNONYMOUS_GROUP));
        }
        if (indelNonCodingCheckBox.isSelected()) {
            grp.add(new Integer(Alignment.POSITION_TYPE_NONCODINGINDEL_GROUP));
        }
        if (allIndelCheckBox.isSelected()) {
            grp.add(new Integer(Alignment.POSITION_TYPE_INDEL_GROUP));
        }
        if (nonSynCheckBox.isSelected()) {
            grp.add(new Integer(Alignment.POSITION_TYPE_NONSYNONYMOUS_GROUP));
        }
        if (transcribedCheckBox.isSelected()) {
            grp.add(new Integer(Alignment.POSITION_TYPE_TRANSCRIBED_GROUP));
        }
        return grp;
    }
    */

    /** Return the number of number of bases to step in each sliding window*/
    public int getStepSize() {
        return step;
    }

    /** Return the number of number of bases in each sliding window*/
    public int getWindowSize() {
        return window;
    }

    /** Return the last site to use in the analysis*/
    public int getEndSite() {
        return end;
    }

    /** Return the first site to use in the analysis*/
    public int getStartSite() {
        return start;
    }

    void runButton_actionPerformed(ActionEvent e) {
        runAnalysis = true;
        setVisible(false);
    }

    void closeButton_actionPerformed(ActionEvent e) {
        runAnalysis = false;
        setVisible(false);
    }

    boolean isCancel() {
        return !runAnalysis;
    }

    void stepTextField_focusLost(FocusEvent e) {
        try {
            step = Integer.parseInt(stepTextField.getText());
            if ((step <= 0) || (step > (theAlignment.getSiteCount() - 1))) {
                step = this.theAlignment.getSiteCount() - 1;
            }
        } catch (Exception ee) {
            step = 100;
        }
        stepTextField.setText(step + "");
    }

    void windowTextField_focusLost(FocusEvent e) {
        try {
            window = Integer.parseInt(windowTextField.getText());
            if ((window <= 0) || (window > (theAlignment.getSiteCount() - 1))) {
                window = this.theAlignment.getSiteCount() - 1;
            }
        } catch (Exception ee) {
            window = 400;
        }
        windowTextField.setText(window + "");
    }
}
