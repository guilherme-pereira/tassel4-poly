/*
 * LinkageDisequilibriumPlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.baseplugins;

import java.awt.event.KeyEvent;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.FilterAlignment;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.pal.popgen.LinkageDisequilibrium;

import javax.swing.*;

import java.awt.*;
import java.awt.event.ActionEvent;

import java.net.URL;

import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;

/**
 *
 * @author Ed Buckler
 */
public class LinkageDisequilibriumPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(LinkageDisequilibriumPlugin.class);
    private boolean myIsRapidAnalysis = true;
    private int myPermutationNumber = 1000;
    private int myWindowSize = 50;
    private LinkageDisequilibrium.testDesign myLDType = LinkageDisequilibrium.testDesign.SlidingWindow;
    private int myTestSite = -1;
    private int myNumAccumulateIntervals = 100;
    private boolean myIsAccumulateResults = false;
    private FilterAlignment myPossibleAlignmentForSiteList;
    private String myPossibleAlignmentName;
    private int[] myPossibleSiteList;
    private LinkageDisequilibrium.HetTreatment myHetTreatment = LinkageDisequilibrium.HetTreatment.Haplotype;

    /**
     * Creates a new instance of LinkageDisequilibriumPlugin
     */
    public LinkageDisequilibriumPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {

        try {
            List<Datum> alignInList = input.getDataOfType(Alignment.class);
            if (alignInList.size() < 1) {
                JOptionPane.showMessageDialog(getParentFrame(), "Invalid selection.  Please select sequence or marker alignment.");
                return null;
            }

            Datum current = alignInList.get(0);

            if (alignInList.size() > 1) {
                try {
                    FilterAlignment temp = (FilterAlignment) alignInList.get(1).getData();
                    if (temp.getBaseAlignment() == current.getData()) {
                        myPossibleAlignmentForSiteList = temp;
                        myPossibleAlignmentName = alignInList.get(1).getName();
                    }
                } catch (Exception e) {
                    // do nothing
                }
            }

            if (isInteractive()) {
                LinkageDiseqDialog myDialog = new LinkageDiseqDialog(((Alignment) current.getData()).getSiteCount(), myPossibleAlignmentName);
                myDialog.setLocationRelativeTo(getParentFrame());
                myDialog.setVisible(true);
                if (myDialog.isCancel()) {
                    return null;
                }
                myLDType = myDialog.getLDType();

                if (myLDType == LinkageDisequilibrium.testDesign.SlidingWindow) {
                    myWindowSize = myDialog.getWindowSize();
                }

                myIsAccumulateResults = myDialog.isAccumulateResults();
                if (myIsAccumulateResults) {
                    myNumAccumulateIntervals = myDialog.getNumAccumulateIntervals();
                }

                if (myLDType == LinkageDisequilibrium.testDesign.SiteByAll) {
                    myTestSite = myDialog.getTestSite();

                }

                if (myLDType == LinkageDisequilibrium.testDesign.SiteList) {
                    myPossibleSiteList = myPossibleAlignmentForSiteList.getBaseSitesShown();

                }

                myHetTreatment = myDialog.getHetTreatment();
            }

            List result = new ArrayList();
            DataSet tds = null;

            tds = processDatum(current);
            if (tds != null) {
                result.add(tds);
                fireDataSetReturned(new PluginEvent(tds, LinkageDisequilibriumPlugin.class));
            }

            return DataSet.getDataSet(result, this);
        } finally {
            fireProgress(100);
        }

    }

    private DataSet processDatum(Datum input) {
        Alignment aa = (Alignment) input.getData();
        LinkageDisequilibrium theLD = new LinkageDisequilibrium(aa, myWindowSize, myLDType, myTestSite, this, myIsAccumulateResults, myNumAccumulateIntervals, myPossibleSiteList, myHetTreatment);
        try {
            theLD.run();
            Datum td = new Datum("LD:" + input.getName(), theLD, "LD Analysis");
            DataSet tds = new DataSet(td, this);
            return tds;
        } catch (Exception e) {
            e.printStackTrace();
            StringBuilder builder = new StringBuilder();
            builder.append("Unable to run Linkage Disequilibrium analysis ");
            builder.append(e.getMessage());
            String str = builder.toString();
            if (isInteractive()) {
                JOptionPane.showMessageDialog(getParentFrame(), str);
            } else {
                myLogger.error("processDatum: " + str);
            }
        }
        return null;
    }

    public boolean isRapidAnalysis() {
        return myIsRapidAnalysis;
    }

    public void setRapidAnalysis(boolean rapidAnalysis) {
        myIsRapidAnalysis = rapidAnalysis;
    }

    public int isPermutationNumber() {
        return myPermutationNumber;
    }

    public void setPermutationNumber(int permutationNumber) {
        myPermutationNumber = permutationNumber;
    }

    public void setLDType(LinkageDisequilibrium.testDesign type) {
        myLDType = type;
    }

    public LinkageDisequilibrium.testDesign getLDType() {
        return myLDType;
    }

    public void setWinSize(int winSize) {
        myWindowSize = winSize;
    }

    public int getWinSize() {
        return myWindowSize;
    }

    public void setNumAccumulateIntervals(int numIntervals) {
        myNumAccumulateIntervals = numIntervals;
    }

    public int getNumAccumulateIntervals() {
        return myNumAccumulateIntervals;
    }

    public void setIsAccumulateResults(boolean accumulate) {
        myIsAccumulateResults = accumulate;
    }

    public boolean getIsAccumulateResults() {
        return myIsAccumulateResults;
    }

    public void setTestSite(int site) {
        myTestSite = site;
    }

    public int getTestSite() {
        return myTestSite;
    }

    public void setHetTreatment(LinkageDisequilibrium.HetTreatment treatment) {
        myHetTreatment = treatment;
    }

    public LinkageDisequilibrium.HetTreatment getHetTreatment() {
        return myHetTreatment;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    @Override
    public ImageIcon getIcon() {
        URL imageURL = LinkageDisequilibriumPlugin.class.getResource("images/LDPlot.gif");
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
    @Override
    public String getButtonName() {
        return "Linkage Disequilibrium";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    @Override
    public String getToolTipText() {
        return "Linkage Disequilibrium";
    }
}

/**
 * Title: TASSEL Description: A java program to deal with diversity Copyright:
 * Copyright (c) 2000 Company: USDA-ARS/NCSU
 *
 * @author Ed Buckler
 * @version 1.0
 */
class LinkageDiseqDialog extends JDialog {

    // private int numberPermutations = 1000;
    private boolean myRunAnalysis = false;
    // private JTextField permNumberTextField = new JTextField();
    // private JLabel jLabel1 = new JLabel();
    private int myTestSite = 0;
    private String myAlignmentForSiteList;
    private int myNumSites;
    private JPanel myPanel = new JPanel();
    private JPanel myLDSelectionPanel = new JPanel();
    private JLabel myLDTypeLabel = new JLabel("Select LD type: ");
    private JComboBox myLDType;
    private JPanel myLDOptionsPanel = new JPanel();
    private JLabel myFullMatrixLabel = new JLabel();
    private JTextField myWindowSizeTextField = new JTextField();
    private JLabel myWindowSizeLabel = new JLabel("LD Window Size: ");
    private JLabel myWindowSizeCountLabel = new JLabel();
    private int myWindowSize = 50;
    private JTextField mySiteByAllTextField = new JTextField();
    private JLabel mySiteByAllLabel = new JLabel("Site: ");
    private JLabel mySiteByAllCountLabel = new JLabel();
    private JLabel mySiteListLabel = new JLabel();
    private JPanel myAccumulateOptionsPanel = new JPanel();
    private JCheckBox myAccumulativeResultsBox = new JCheckBox("Accumulate R2 Results");
    private JTextField myAccumulativeResultsTextField = new JTextField();
    private JLabel myAccumulativeResultsLabel = new JLabel("Number Of Intervals: ");
    private int myNumAccumulativeInterval = 100;
    private JPanel myHetTreatmentPanel = new JPanel();
    private JLabel myHetTreatmentLabel = new JLabel("How to treat heterozygous calls:");
    private JComboBox myHetTreatment;
    private JPanel myButtonsPanel = new JPanel();
    private JButton myRunButton = new JButton("Run");
    private JButton myCloseButton = new JButton("Close");

    public LinkageDiseqDialog(int numSites, String alignmentForSiteList) {
        super((Frame) null, "Linkage Disequilibrium", true);
        myAlignmentForSiteList = alignmentForSiteList;
        myNumSites = numSites;
        // numberPermutations = numberPermutations;
        try {
            jbInit();
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    private void jbInit() throws Exception {

        String[] ldTypes;
        if (myAlignmentForSiteList != null) {
            ldTypes = new String[4];
            ldTypes[0] = "Full Matrix";
            ldTypes[1] = "Sliding Window";
            ldTypes[2] = "Site by All";
            ldTypes[3] = "Site List";
        } else {
            ldTypes = new String[3];
            ldTypes[0] = "Full Matrix";
            ldTypes[1] = "Sliding Window";
            ldTypes[2] = "Site by All";
        }

        myLDType = new JComboBox(ldTypes);
        myLDType.setSelectedIndex(1);

        myFullMatrixLabel.setText("Full LD with " + myNumSites * ((myNumSites - 1) / 2) + " comparisons.");

        long n = Math.min(myNumSites - 1, myWindowSize);

        myWindowSizeCountLabel.setText("Sliding Window LD with " + (((n * (n + 1)) / 2) + (myNumSites - n - 1) * n) + " comparisons.");
        mySiteByAllCountLabel.setText("Site by All LD with " + (myNumSites - 1) + " comparisons.");

        if (myAlignmentForSiteList != null) {
            mySiteListLabel.setText("Sites From: " + myAlignmentForSiteList);
        }

        myLDSelectionPanel.setLayout(new GridBagLayout());
        myLDSelectionPanel.add(myLDTypeLabel, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0, GridBagConstraints.LINE_END, GridBagConstraints.NONE, new Insets(10, 7, 0, 0), 0, 0));
        myLDSelectionPanel.add(myLDType, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0, GridBagConstraints.LINE_START, GridBagConstraints.NONE, new Insets(10, 0, 0, 0), 0, 0));
//        myLDSelectionPanel.add(myAccumulativeResultsBox, new GridBagConstraints(0, 1, 2, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));

        myLDType.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                ldType_actionPerformed(e);
            }
        });

        myAccumulativeResultsBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                accumulativeResultsBox_actionPerformed(e);
            }
        });

        myLDOptionsPanel.setLayout(new GridBagLayout());
        myLDOptionsPanel.add(myFullMatrixLabel, new GridBagConstraints(0, 0, 2, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(10, 0, 0, 0), 0, 0));
        myFullMatrixLabel.setVisible(false);

        myWindowSizeTextField.setText("" + myWindowSize);
        myLDOptionsPanel.add(myWindowSizeLabel, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0, GridBagConstraints.LINE_END, GridBagConstraints.NONE, new Insets(10, 7, 0, 0), 0, 0));
        myLDOptionsPanel.add(myWindowSizeTextField, new GridBagConstraints(1, 0, 1, 1, 1.0, 0.0, GridBagConstraints.LINE_START, GridBagConstraints.HORIZONTAL, new Insets(10, 0, 0, 7), 0, 0));
        myLDOptionsPanel.add(myWindowSizeCountLabel, new GridBagConstraints(0, 1, 2, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(7, 7, 0, 7), 0, 0));

        myWindowSizeTextField.addKeyListener(new java.awt.event.KeyListener() {
            public void keyTyped(KeyEvent e) {
//                windowSizeTextField_actionPerformed(e);
            }

            public void keyPressed(KeyEvent ke) {
//                throw new UnsupportedOperationException("Not supported yet.");
            }

            public void keyReleased(KeyEvent e) {
                windowSizeTextField_actionPerformed(e);
            }
        });

        myLDOptionsPanel.add(mySiteByAllLabel, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0, GridBagConstraints.LINE_END, GridBagConstraints.NONE, new Insets(10, 7, 0, 0), 0, 0));
        myLDOptionsPanel.add(mySiteByAllTextField, new GridBagConstraints(1, 0, 1, 1, 1.0, 0.0, GridBagConstraints.LINE_START, GridBagConstraints.HORIZONTAL, new Insets(10, 0, 0, 7), 0, 0));
        myLDOptionsPanel.add(mySiteByAllCountLabel, new GridBagConstraints(0, 1, 2, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(7, 7, 0, 7), 0, 0));
        mySiteByAllLabel.setVisible(false);
        mySiteByAllTextField.setVisible(false);
        mySiteByAllCountLabel.setVisible(false);

        myLDOptionsPanel.add(mySiteListLabel, new GridBagConstraints(0, 0, 2, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(10, 0, 0, 0), 0, 0));
        mySiteListLabel.setVisible(false);

        //Heterozygous treatment options
        String[] hetTypes = {"Ignore (inbred lines only)", "Set to missing", "Treat as third state"};
        myHetTreatment = new JComboBox(hetTypes);
        myHetTreatment.setSelectedIndex(1);
        /*myHetTreatment.addActionListener(new java.awt.event.ActionListener() {	//Nothing reall happens when switch option, only at running, so no action listener needed
         public void actionPerformed(ActionEvent e) {
         hetTreatment_actionPerformed(e);
         }
         });*/
        myHetTreatmentPanel.setLayout(new GridBagLayout());
        myHetTreatmentPanel.add(myHetTreatmentLabel, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0, GridBagConstraints.LINE_START, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        myHetTreatmentPanel.add(myHetTreatment, new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));

        //Accumulate R2 results panel
        myAccumulativeResultsTextField.setText(myNumAccumulativeInterval + "");
        myAccumulateOptionsPanel.setLayout(new GridBagLayout());
        myAccumulateOptionsPanel.add(myAccumulativeResultsBox, new GridBagConstraints(0, 0, 2, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(25, 7, 0, 7), 0, 0));
        myAccumulateOptionsPanel.add(myAccumulativeResultsLabel, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0, GridBagConstraints.LINE_END, GridBagConstraints.NONE, new Insets(7, 7, 0, 0), 0, 0));
        myAccumulateOptionsPanel.add(myAccumulativeResultsTextField, new GridBagConstraints(1, 1, 1, 1, 1.0, 0.0, GridBagConstraints.LINE_START, GridBagConstraints.HORIZONTAL, new Insets(7, 0, 0, 7), 0, 0));
        myAccumulativeResultsLabel.setVisible(false);
        myAccumulativeResultsTextField.setVisible(false);

        myButtonsPanel.setLayout(new GridBagLayout());
        myButtonsPanel.add(myRunButton, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0, GridBagConstraints.PAGE_END, GridBagConstraints.NONE, new Insets(0, 150, -5, 0), 0, 0));
        myButtonsPanel.add(myCloseButton, new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0, GridBagConstraints.PAGE_END, GridBagConstraints.NONE, new Insets(0, 0, -5, 0), 0, 0));

        myPanel.setLayout(new GridLayout(5, 1));
        myPanel.add(myLDSelectionPanel);
        myPanel.add(myLDOptionsPanel);
        myPanel.add(myHetTreatmentPanel);
        myPanel.add(myAccumulateOptionsPanel);
        myPanel.add(myButtonsPanel);

        myRunButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                runButton_actionPerformed(e);
            }
        });

        myCloseButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                closeButton_actionPerformed(e);
            }
        });
        getContentPane().add(myPanel);
        getContentPane().setPreferredSize(new Dimension(320, 375));
    }

    private void hideOptions() {
        myWindowSizeTextField.setVisible(false);
        myWindowSizeLabel.setVisible(false);
        mySiteByAllTextField.setVisible(false);
        mySiteByAllLabel.setVisible(false);
    }

    public boolean isAccumulateResults() {
        return myAccumulativeResultsBox.isSelected();
    }

    public LinkageDisequilibrium.testDesign getLDType() {
        if (myLDType.getSelectedIndex() == 0) {
            return LinkageDisequilibrium.testDesign.All;
        } else if (myLDType.getSelectedIndex() == 1) {
            return LinkageDisequilibrium.testDesign.SlidingWindow;
        } else if (myLDType.getSelectedIndex() == 2) {
            return LinkageDisequilibrium.testDesign.SiteByAll;
        } else if (myLDType.getSelectedIndex() == 3) {
            return LinkageDisequilibrium.testDesign.SiteList;
        } else {
            throw new IllegalStateException("LinkageDisequilibriumPlugin: getLDType: No known LD Type selected.");
        }
    }

    /**
     * Returns whether the run button was chosen
     */
    public boolean isRunAnalysis() {
        return myRunAnalysis;
    }

    /**
     * Returns whether the run button was chosen
     */
    public boolean isCancel() {
        return !myRunAnalysis;
    }

    /**
     * Return the window size
     */
    public int getWindowSize() {
        return myWindowSize;
    }

    public int getTestSite() {
        return myTestSite;
    }

    public int getNumAccumulateIntervals() {
        return myNumAccumulativeInterval;
    }

    public LinkageDisequilibrium.HetTreatment getHetTreatment() {
        if (myHetTreatment.getSelectedIndex() == 0) {
            return LinkageDisequilibrium.HetTreatment.Haplotype;
        } else if (myHetTreatment.getSelectedIndex() == 1) {
            return LinkageDisequilibrium.HetTreatment.Homozygous;
        } else if (myHetTreatment.getSelectedIndex() == 2) {
            return LinkageDisequilibrium.HetTreatment.Genotype;
        } else {
            throw new IllegalStateException("LinkageDisequilibriumPlugin: getLDType: No known LD Type selected.");
        }
    }

    void ldType_actionPerformed(ActionEvent e) {
        if (myLDType.getSelectedIndex() == 0) {
            myFullMatrixLabel.setVisible(true);
            myWindowSizeLabel.setVisible(false);
            myWindowSizeCountLabel.setVisible(false);
            myWindowSizeTextField.setVisible(false);
            mySiteByAllLabel.setVisible(false);
            mySiteByAllCountLabel.setVisible(false);
            mySiteByAllTextField.setVisible(false);
            mySiteListLabel.setVisible(false);
        } else if (myLDType.getSelectedIndex() == 1) {
            myFullMatrixLabel.setVisible(false);
            myWindowSizeLabel.setVisible(true);
            myWindowSizeCountLabel.setVisible(true);
            myWindowSizeTextField.setVisible(true);
            mySiteByAllLabel.setVisible(false);
            mySiteByAllCountLabel.setVisible(false);
            mySiteByAllTextField.setVisible(false);
            mySiteListLabel.setVisible(false);
        } else if (myLDType.getSelectedIndex() == 2) {
            myFullMatrixLabel.setVisible(false);
            myWindowSizeLabel.setVisible(false);
            myWindowSizeCountLabel.setVisible(false);
            myWindowSizeTextField.setVisible(false);
            mySiteByAllLabel.setVisible(true);
            mySiteByAllCountLabel.setVisible(true);
            mySiteByAllTextField.setVisible(true);
            mySiteListLabel.setVisible(false);
        } else if (myLDType.getSelectedIndex() == 3) {
            myFullMatrixLabel.setVisible(false);
            myWindowSizeLabel.setVisible(false);
            myWindowSizeCountLabel.setVisible(false);
            myWindowSizeTextField.setVisible(false);
            mySiteByAllLabel.setVisible(false);
            mySiteByAllCountLabel.setVisible(false);
            mySiteByAllTextField.setVisible(false);
            mySiteListLabel.setVisible(true);
        }
    }

    void windowSizeTextField_actionPerformed(KeyEvent e) {
        try {

            int newSize = Integer.parseInt(myWindowSizeTextField.getText());
            long n = Math.min(myNumSites - 1, newSize);
            myWindowSizeCountLabel.setText("Sliding Window LD with " + (((n * (n + 1)) / 2) + (myNumSites - n - 1) * n) + " comparisons.");

        } catch (Exception err) {
            System.out.println("error!");
        }
    }

    void accumulativeResultsBox_actionPerformed(ActionEvent e) {
        if (myAccumulativeResultsBox.isSelected()) {
            myAccumulativeResultsLabel.setVisible(true);
            myAccumulativeResultsTextField.setVisible(true);
        } else {
            myAccumulativeResultsLabel.setVisible(false);
            myAccumulativeResultsTextField.setVisible(false);
        }
    }

    void runButton_actionPerformed(ActionEvent e) {

        //        try {
        //            numberPermutations = Integer.parseInt(permNumberTextField.getText());
        //        } catch (Exception ee) {
        //            permNumberTextField.setText("Set Integer");
        //            return;
        //        }

        if (getLDType() == LinkageDisequilibrium.testDesign.SlidingWindow) {
            try {
                myWindowSize = Integer.parseInt(myWindowSizeTextField.getText());
            } catch (Exception ee) {
                myWindowSizeTextField.setText("Set Integer");
                return;
            }
        }

        if (getLDType() == LinkageDisequilibrium.testDesign.SiteByAll) {
            try {
                myTestSite = Integer.parseInt(mySiteByAllTextField.getText());
            } catch (Exception ee) {
                mySiteByAllTextField.setText("Set Integer");
                return;
            }
        }

        if (isAccumulateResults()) {
            try {
                myNumAccumulativeInterval = Integer.parseInt(myAccumulativeResultsTextField.getText());
            } catch (Exception ee) {
                myAccumulativeResultsTextField.setText("Set Integer");
                return;
            }
        }
        myRunAnalysis = true;
        setVisible(false);
    }

    void closeButton_actionPerformed(ActionEvent e) {
        myRunAnalysis = false;
        setVisible(false);
    }
}