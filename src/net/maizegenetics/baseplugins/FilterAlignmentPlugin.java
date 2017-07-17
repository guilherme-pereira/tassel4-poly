/*
 * FilterAlignmentPlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.baseplugins;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.alignment.CombineAlignment;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.pal.alignment.AlignmentUtils;

import javax.swing.*;
import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;

import java.net.URL;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import net.maizegenetics.pal.alignment.BitAlignment;

import org.apache.log4j.Logger;

/**
 *
 * @author Ed Buckler, Terry, Jon
 */
public class FilterAlignmentPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(FilterAlignmentPlugin.class);
    private int myStart = 0;
    private int myEnd = -1;
    private int myStartPos = -1;
    private int myEndPos = -1;
    private Locus myLocus = null;
    private String myLocusStr = null;
    private int myMinCount = 1;
    private double myMinFreq = 0.01;
    private double myMaxFreq = 1.0;
    private boolean myExtractIndels = false;
    private boolean myFilterMinorSNPs = false;
    private boolean myIsUseAllSiteTypes = true;
    private boolean myDoSlidingHaps = false;
    private int myWinSize = 3;
    private int myStepSize = 3;

    /**
     * Creates a new instance of FilterAlignmentPlugin
     */
    public FilterAlignmentPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet performFunction(DataSet input) {

        try {

            List<Datum> alignInList = input.getDataOfType(Alignment.class);

            if (alignInList.size() < 1) {
                String gpMessage = "Invalid selection.  Please select genotype alignment.";
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), gpMessage);
                } else {
                    myLogger.error(gpMessage);
                }
                return null;
            }
            List<Datum> alignOutList = new ArrayList<Datum>();
            Iterator<Datum> itr = alignInList.iterator();
            while (itr.hasNext()) {
                Datum current = itr.next();
                Datum result = processDatum(current, isInteractive());
                if (result != null) {
                    alignOutList.add(result);
                }
            }

            if (alignOutList.isEmpty()) {
                return null;
            }

            DataSet output = new DataSet(alignOutList, this);
            fireDataSetReturned(new PluginEvent(output, FilterAlignmentPlugin.class));

            return output;

        } finally {
            fireProgress(100);
        }

    }

    private Datum processDatum(Datum inDatum, boolean isInteractive) {
        Alignment aa = (Alignment) inDatum.getData();

        if (myEnd == -1) {
            myEnd = aa.getSiteCount() - 1;
        }

        Locus theLocus = myLocus;
        if ((theLocus == null) && (myLocusStr != null)) {
            theLocus = aa.getLocus(myLocusStr);
            if (theLocus == null) {
                throw new IllegalStateException("FilterAlignmentPlugin: processDatum: Alignment doesn't contain locus: " + myLocusStr);
            }
        }

        if (myStartPos != -1) {
            myStart = aa.getSiteOfPhysicalPosition(myStartPos, theLocus);
            if (myStart < 0) {
                myStart = -(myStart + 1);
            }
        }

        if (myEndPos != -1) {
            myEnd = aa.getSiteOfPhysicalPosition(myEndPos, theLocus);
            if (myEnd < 0) {
                myEnd = -(myEnd + 2);
            }
        }

        if (myEnd < myStart) {
            throw new IllegalStateException("FilterAlignmentPlugin: processDatum: Start Site can't be after End Site.");
        }

        if (isInteractive) {
            DataFilterAlignmentDialog theDialog = new DataFilterAlignmentDialog(aa, getParentFrame());
            theDialog.setLocationRelativeTo(getParentFrame());
            theDialog.setVisible(true);
            if (theDialog.isCanceled()) {
                return null;
            }
            aa = theDialog.getChromFilteredAlignment();
            myStart = theDialog.getStart();
            myEnd = theDialog.getEnd();
            myMinCount = theDialog.getMinimumCount();
            myMinFreq = theDialog.getMinimumFrequency();
            myMaxFreq = theDialog.getMaximumFrequency();
            myExtractIndels = theDialog.isExtractIndels();
            myFilterMinorSNPs = theDialog.isRemoveMinorSNPs();
            myIsUseAllSiteTypes = theDialog.isAllSiteIncluded();
            // char[] siteType = theDialog.getIncludedPositionTypes();
            myDoSlidingHaps = theDialog.isUseSlidingWindow();
            myWinSize = theDialog.getWindowSize();
            myStepSize = theDialog.getStepSize();

            theDialog.dispose();
        }

        if (myStart >= aa.getSiteCount()) {
            throw new IllegalArgumentException("FilterAlignmentPlugin: starting site can't be past end of alignment.");
        }
        if ((myEnd < 0) || (myEnd < myStart)) {
            throw new IllegalArgumentException("FilterAlignmentPlugin: end site can't be less than zero or less that starting site.");
        }

        if (myStart < 0) {
            myStart = 0;
        }
        if (myEnd >= aa.getSiteCount()) {
            myEnd = aa.getSiteCount() - 1;
        }

        Alignment naa = aa;

        if (myFilterMinorSNPs) {
            naa = BitAlignment.getInstance(naa, 2, false, true);
        }

        if ((myStart != 0) || (myEnd < (naa.getSiteCount() - 1))) {
            naa = AlignmentUtils.removeSitesOutsideRange(naa, myStart, myEnd);
        }
        if (myExtractIndels) {
            //naa = AnnotatedAlignmentUtils.removeSitesBasedOnFreqIgnoreMissing(naa, myMinFreq, myMinCount);
            throw new UnsupportedOperationException();
        } else {
            //naa = AnnotatedAlignmentUtils.removeSitesBasedOnFreqIgnoreGapsMissing(naa, myMinFreq, myMinCount);
        }
        naa = AlignmentUtils.removeSitesBasedOnFreqIgnoreMissing(naa, myMinFreq, myMaxFreq, myMinCount);
        if (myDoSlidingHaps) {
            //naa = AnnotatedAlignmentUtils.extractSlidingHaplotypes(naa, myWinSize, myStepSize);
            throw new UnsupportedOperationException();
        }
        String theComment;
        StringBuilder builder = new StringBuilder();
        Locus[] loci = naa.getLoci();
        builder.append(inDatum.getName());
        builder.append("_");
        if ((loci != null) && (loci.length != 0)) {
            boolean first = true;
            for (int i = 0; i < loci.length; i++) {
                String name = loci[i].getName();
                if ((name != null) && (name.length() != 0)) {
                    if (first) {
                        builder.append("chr");
                        first = false;
                    }
                    builder.append(name);
                    builder.append("_");
                }
            }
        }
        if (naa.getSiteCount() > 1) {
            builder.append(naa.getPositionInLocus(0));
            builder.append("-");
            builder.append(naa.getPositionInLocus(naa.getSiteCount() - 1));
        }
        String theName = builder.toString();
        if (myDoSlidingHaps) {
            //theName = "Sliding_Haps_" + inDatum.getName();
            theComment = "Sliding Haplotypes.\n";
        } else if (myExtractIndels) {
            //theName = "Indels_" + inDatum.getName();
            theComment = "Indels\n";
        } else if (myFilterMinorSNPs) {
            //theName = "Point_" + inDatum.getName();
            theComment = "Point Poly.\nMinor SNPs Removed\n";
        } else {
            theComment = "Point Poly.\n";
        }
        if (naa.getSiteCount() != 0) {
            myLogger.info("Resulting Number Sites: " + naa.getSiteCount());
            return new Datum(theName, naa, theComment);
        } else {
            if (isInteractive()) {
                JOptionPane.showMessageDialog(getParentFrame(), "No available SNPs given filter parameters.");
            } else {
                myLogger.warn("No available SNPs given filter parameters.");

            }
            return null;
        }
    }

    public int getStart() {
        return myStart;
    }

    public void setStart(int start) {
        myStart = start;
    }

    public int getEnd() {
        return myEnd;
    }

    public void setEnd(int end) {
        myEnd = end;
    }

    public int getStartPos() {
        return myStartPos;
    }

    public void setStartPos(int start) {
        myStartPos = start;
    }

    public int getEndPos() {
        return myEndPos;
    }

    public void setEndPos(int end) {
        myEndPos = end;
    }

    public String getLocusStr() {
        return myLocusStr;
    }

    public void setLocusStr(String name) {
        myLocus = null;
        myLocusStr = name;
    }

    public Locus getLocus() {
        return myLocus;
    }

    public void setLocus(Locus locus) {
        myLocusStr = null;
        myLocus = locus;
    }

    public int getMinCount() {
        return myMinCount;
    }

    public void setMinCount(int minCount) {
        myMinCount = minCount;
    }

    public double getMinFreq() {
        return myMinFreq;
    }

    public void setMinFreq(double minFreq) {
        myMinFreq = minFreq;
    }

    public double getMaxFreq() {
        return myMaxFreq;
    }

    public void setMaxFreq(double maxFreq) {
        myMaxFreq = maxFreq;
    }

    public boolean isExtractIndels() {
        return myExtractIndels;
    }

    public void setExtractIndels(boolean extractIndels) {
        myExtractIndels = extractIndels;
    }

    public boolean isFilterMinorSNPs() {
        return myFilterMinorSNPs;
    }

    public void setFilterMinorSNPs(boolean filterMinorSNPs) {
        myFilterMinorSNPs = filterMinorSNPs;
    }

    public boolean isUseAllSiteTypes() {
        return myIsUseAllSiteTypes;
    }

    public void setUseAllSiteTypes(boolean useAllSiteTypes) {
        myIsUseAllSiteTypes = useAllSiteTypes;
    }

    public boolean isDoSlidingHaps() {
        return myDoSlidingHaps;
    }

    public void setDoSlidingHaps(boolean doSlidingHaps) {
        myDoSlidingHaps = doSlidingHaps;
    }

    public int getWinSize() {
        return myWinSize;
    }

    public void setWinSize(int winSize) {
        myWinSize = winSize;
    }

    public int getStepSize() {
        return myStepSize;
    }

    public void setStepSize(int stepSize) {
        myStepSize = stepSize;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        URL imageURL = FilterAlignmentPlugin.class.getResource("images/Filter.gif");
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
        return "Sites";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Filter and Process Alignment";
    }
}

class DataFilterAlignmentDialog extends JDialog {

    Alignment theAlignment;
    Alignment chromFilteredAlignment;
    int start = 0, end, startPos, endPos, minCount = 0, totalSeq, siteCount = 0;
    String[] chromsAvailable;
    String[] chromsSelected;
    double minFreq = 0.01;
    double maxFreq = 1.0;
    double minPercentage = 0.5;
    boolean isCanceled = true;
    private final static int TEXT_FIELD_WIDTH = 8;
    private final static int INVALID_VALUE = -999;   // what is returned when inappropriate to return any value
    private JPanel mainPanel = new JPanel();
    private JButton filterButton = new JButton();
    private JButton cancelButton = new JButton();
    private JButton chromSelectButton = new JButton();
    private JLabel lblFilterAlignment = new JLabel();
    private JLabel lblTotalSequences = new JLabel();
    private JLabel lblMinCount = new JLabel();
    private JLabel lblMinFreq = new JLabel();
    private JLabel lblMaxFreq = new JLabel();
    private JLabel lblStartSite = new JLabel();
    private JLabel lblDistanceFromEndSite = new JLabel();
    private JLabel lblEndSite = new JLabel();
    private JLabel lblSeqLength = new JLabel();
    private JLabel lblWinSize = new JLabel();
    private JLabel lblStepSize = new JLabel();
    private JLabel lblPosType = new JLabel();
    private JLabel lblSiteIndex = new JLabel();
    private JLabel lblSitePos = new JLabel();
    private JTextField countTextField = new JTextField();
    private JTextField endTextField = new JTextField();
    private JTextField startTextField = new JTextField();
    private JTextField endPosTextField = new JTextField();
    private JTextField startPosTextField = new JTextField();
    private JTextField freqTextField = new JTextField();
    private JTextField maxFreqTextField = new JTextField();
    private JTextField winSizeTextField = new JTextField(TEXT_FIELD_WIDTH);
    private JTextField stepSizeTextField = new JTextField(TEXT_FIELD_WIDTH);
    private JPanel checkBoxPanel = new JPanel();
    //private JCheckBox indelCheckBox = new JCheckBox();
    private JCheckBox removeMinorCheckBox = new JCheckBox();
    private JCheckBox slidingHapCheckBox = new JCheckBox();
    private GridBagLayout gridBagLayout2 = new GridBagLayout();
    private String lblEndString = "End Position:";
    private boolean doBatchAnalysis = false;  // for batch analysis logic
    private boolean isStartTextFieldNumeric = true;
    private boolean isEndTextFieldNumeric = true;
    private boolean isStartPosTextFieldNumeric = true;
    private boolean isEndPosTextFieldNumeric = true;
    private boolean isChromSelectionValid = true;
    private ChromosomeFilterDialog myChromFilter;

    public DataFilterAlignmentDialog(Alignment a, Frame f) {
        super(f, "Filter Alignment", true);
        theAlignment = a;
        chromFilteredAlignment = theAlignment;
        chromsAvailable = new String[theAlignment.getNumLoci()];
        for (int i = 0; i < chromsAvailable.length; i++) {
            chromsAvailable[i] = theAlignment.getLoci()[i].getName().trim();
        }
        totalSeq = theAlignment.getSequenceCount();
        siteCount = theAlignment.getSiteCount();
        lblSeqLength.setText(" of " + (siteCount - 1) + " sites");
        lblMinCount.setText("Minimum PERCENTAGE:");
        start = 0;
        end = siteCount - 1;
        startPos = theAlignment.getPositionInLocus(0);
        endPos = theAlignment.getPositionInLocus(siteCount - 1);
        minCount = TasselPrefs.getFilterAlignPluginMinCount();
        minFreq = TasselPrefs.getFilterAlignPluginMinFreq();
        maxFreq = TasselPrefs.getFilterAlignPluginMaxFreq();
        myChromFilter = new ChromosomeFilterDialog(chromsAvailable, f);

        try {
            initUI();
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    void initUI() throws Exception {

        // mainPanel.setBackground(SystemColor.menu);
        mainPanel.setMinimumSize(new Dimension(640, 480));
        mainPanel.setPreferredSize(new Dimension(640, 480));
        mainPanel.setLayout(gridBagLayout2);
        filterButton.setText("Filter");
        filterButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                filterButton_actionPerformed(e);
            }
        });
        cancelButton.setMaximumSize(new Dimension(63, 27));
        cancelButton.setMinimumSize(new Dimension(63, 27));
        cancelButton.setText("Cancel");
        cancelButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                cancelButton_actionPerformed(e);
            }
        });
        chromSelectButton.setText("Select Chromosomes...");
        chromSelectButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                chromSelectButton_actionPerformed(e);
            }
        });

        JPanel buttonPanel = new JPanel();
        buttonPanel.setLayout(new GridLayout(1, 3));
        buttonPanel.add(filterButton);
        if (chromsAvailable.length != 1) {
            buttonPanel.add(chromSelectButton);
        }
        buttonPanel.add(cancelButton);

        lblFilterAlignment.setFont(new java.awt.Font("Dialog", 1, 16));
        lblFilterAlignment.setText("Filter Alignment");
        // lblMinCount.setText("Minimum Count:");
        lblMinFreq.setText("Minimum Frequency:");
        lblMaxFreq.setText("Maximum Frequency:");
        lblStartSite.setText("Start Position:");
        lblPosType.setText("Position Type:");
        lblSiteIndex.setText(" Position index");
        lblSitePos.setText(" Physical Position (AGP)");
        lblEndSite.setText(lblEndString);

        countTextField.setMinimumSize(new Dimension(40, 25));
        countTextField.setPreferredSize(new Dimension(63, 25));

        if (doBatchAnalysis) {
            countTextField.setText(minPercentage + "");
        } else {
            countTextField.setText(minCount + "");
        }

        countTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(FocusEvent e) {
                countTextField_focusLost(e);
            }
        });
        // distanceFromEndTextField.setText(siteCount - myEnd+"");
        //        this.setDistanceFromEndTextField(myEnd);
        //        distanceFromEndTextField.addFocusListener(new java.awt.event.FocusAdapter() {
        //
        //            public void focusLost(FocusEvent e) {
        //                distanceFromEndTextField_focusLost(e);
        //            }
        //        });
        //        distanceFromEndTextField.setPreferredSize(new Dimension(63, 25));
        //        distanceFromEndTextField.setMinimumSize(new Dimension(40, 25));

        setEndTextField();
        endTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(FocusEvent e) {
                endTextField_focusLost(e);
            }
        });
        endTextField.setPreferredSize(new Dimension(63, 25));
        endTextField.setMinimumSize(new Dimension(40, 25));

        startTextField.setText(start + "");
        startTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(FocusEvent e) {
                startTextField_focusLost(e);
            }
        });
        startTextField.setPreferredSize(new Dimension(63, 25));
        startTextField.setMinimumSize(new Dimension(40, 25));

        endPosTextField.setText(endPos + "");
        endPosTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(FocusEvent e) {
                endPosTextField_focusLost(e);
            }
        });
        endPosTextField.setPreferredSize(new Dimension(63, 25));
        endPosTextField.setMinimumSize(new Dimension(40, 25));

        startPosTextField.setText(startPos + "");
        startPosTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(FocusEvent e) {
                startPosTextField_focusLost(e);
            }
        });
        startPosTextField.setPreferredSize(new Dimension(63, 25));
        startPosTextField.setMinimumSize(new Dimension(40, 25));

        try {
            minFreq = TasselPrefs.getFilterAlignPluginMinFreq();
        } catch (NullPointerException npe) {
            npe.printStackTrace();
            System.err.println(" There is an issue with the settings: " + npe);
        }

        freqTextField.setText(minFreq + "");
        freqTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(FocusEvent e) {
                freqTextField_focusLost(e);
            }
        });
        freqTextField.setPreferredSize(new Dimension(63, 25));
        freqTextField.setMinimumSize(new Dimension(40, 25));

        maxFreqTextField.setText(maxFreq + "");
        maxFreqTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(FocusEvent e) {
                maxFreqTextField_focusLost(e);
            }
        });
        maxFreqTextField.setPreferredSize(new Dimension(63, 25));
        maxFreqTextField.setMinimumSize(new Dimension(40, 25));


        if (!doBatchAnalysis) {
            lblMinCount.setText("Minimum Count:");
            lblTotalSequences.setText(" out of " + totalSeq + " sequences");
        }
        /*        siteGroupPanel.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(), "Site Type")   );
         siteGroupPanel.setOpaque(false);
        
         siteGroupPanel.setLayout(gridBagLayout1);
        
         transcribedRadioButton.setText("Transcribed");
         transcribedRadioButton.setOpaque(false);
         noncodingRadioButton.setText("Noncoding");
         noncodingRadioButton.setOpaque(false);
         codingRadioButton.setText("Coding");
         codingRadioButton.setOpaque(false);
         allRadioButton.setText("All");
         allRadioButton.setOpaque(false);
         allRadioButton.setSelected(true);
         */
        //indelCheckBox.setText("Extract Indels");
        //indelCheckBox.setOpaque(false);
        // aminoCheckBox.setText("Convert To Amino Acid");
        // aminoCheckBox.setOpaque(false);

        removeMinorCheckBox.setText("Remove minor SNP states");
        removeMinorCheckBox.setOpaque(false);
        slidingHapCheckBox.setText("Generate haplotypes via sliding window");
        slidingHapCheckBox.setOpaque(false);
        slidingHapCheckBox.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                boolean enable = slidingHapCheckBox.isSelected();
                if (!enable) {
                    winSizeTextField.setText("");
                    stepSizeTextField.setText("");
                }
                winSizeTextField.setEnabled(enable);
                stepSizeTextField.setEnabled(enable);
                lblWinSize.setEnabled(enable);
                lblStepSize.setEnabled(enable);
            }
        });   //.addActionListener(new JCheckBoxListener());
        lblStepSize.setText("Step Length");
        lblWinSize.setText("Haplotype Length");
        lblStepSize.setLabelFor(stepSizeTextField);
        lblWinSize.setLabelFor(winSizeTextField);
        winSizeTextField.setEnabled(false);
        stepSizeTextField.setEnabled(false);
        lblWinSize.setEnabled(false);
        lblStepSize.setEnabled(false);

        JPanel winSizePanel = new JPanel();
        FlowLayout fl = new FlowLayout(FlowLayout.RIGHT);
        winSizePanel.setLayout(fl);
        winSizePanel.add(lblWinSize);
        winSizePanel.add(winSizeTextField);

        JPanel stepSizePanel = new JPanel();
        stepSizePanel.setLayout(fl);
        stepSizePanel.add(lblStepSize);
        stepSizePanel.add(stepSizeTextField);

        checkBoxPanel.setLayout(new GridLayout(6, 1));
        //checkBoxPanel.add(indelCheckBox);
        checkBoxPanel.add(removeMinorCheckBox);
        checkBoxPanel.add(slidingHapCheckBox);
        checkBoxPanel.add(winSizePanel);
        checkBoxPanel.add(stepSizePanel);
        mainPanel.add(lblFilterAlignment, new GridBagConstraints(1, 0, 4, 1, 1.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(5, 0, 0, 5), 0, 12));

        mainPanel.add(lblMinCount, new GridBagConstraints(0, 1, 2, 1, 1.0, 0.0, GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 37, 14));
        mainPanel.add(lblMinFreq, new GridBagConstraints(0, 2, 1, 1, 1.0, 0.0, GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 37, 14));
        mainPanel.add(lblMaxFreq, new GridBagConstraints(0, 3, 1, 1, 1.0, 0.0, GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 37, 14));
        mainPanel.add(lblPosType, new GridBagConstraints(0, 4, 1, 1, 1.0, 0.0, GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 37, 14));
        mainPanel.add(lblStartSite, new GridBagConstraints(0, 5, 1, 1, 1.0, 1.0, GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 37, 14));
        mainPanel.add(lblDistanceFromEndSite, new GridBagConstraints(0, 6, 1, 1, 1.0, 1.0, GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 37, 14));

        mainPanel.add(countTextField, new GridBagConstraints(2, 1, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 12, 0));
        mainPanel.add(lblTotalSequences, new GridBagConstraints(3, 1, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 42, 3));
        mainPanel.add(freqTextField, new GridBagConstraints(2, 2, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 12, 0));
        mainPanel.add(maxFreqTextField, new GridBagConstraints(2, 3, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 12, 0));
        mainPanel.add(lblSiteIndex, new GridBagConstraints(2, 4, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 42, 3));
        mainPanel.add(startTextField, new GridBagConstraints(2, 5, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 12, 0));

        mainPanel.add(lblSitePos, new GridBagConstraints(3, 4, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 42, 3));
        mainPanel.add(startPosTextField, new GridBagConstraints(3, 5, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 12, 0));

        lblSitePos.setVisible(theAlignment.getNumLoci() == 1 && startPos >= 0);
        startPosTextField.setVisible(theAlignment.getNumLoci() == 1 && startPos >= 0);

        if (!doBatchAnalysis) {
            mainPanel.add(lblEndSite, new GridBagConstraints(0, 6, 1, 1, 1.0, 1.0, GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 37, 14));
            mainPanel.add(endTextField, new GridBagConstraints(2, 6, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 12, 0));
            mainPanel.add(endPosTextField, new GridBagConstraints(3, 6, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 12, 0));
            endPosTextField.setVisible(theAlignment.getNumLoci() == 1 && endPos >= 0);
            mainPanel.add(lblSeqLength, new GridBagConstraints(2, 7, 1, 1, 1.0, 0.6, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 3));
        }

        mainPanel.add(checkBoxPanel, new GridBagConstraints(0, 8, 2, 1, 1.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(15, 2, 11, 19), 0, 15));
        mainPanel.add(buttonPanel, new GridBagConstraints(0, 12, 4, 1, 1.0, 1.0, GridBagConstraints.SOUTH, GridBagConstraints.NONE, new Insets(20, 50, 14, 12), 25, 12));

        this.add(mainPanel, BorderLayout.CENTER);
    }

    private void filterButton_actionPerformed(ActionEvent e) {
        if (start < 0) {
            JOptionPane.showMessageDialog(this.getParent(), "Start site must be non negative.");
        } else if (start >= siteCount) {
            JOptionPane.showMessageDialog(this.getParent(), "Start site must be less than " + siteCount + ".");
        } else if (isStartTextFieldNumeric == false) {
            JOptionPane.showMessageDialog(this.getParent(), "Start site must be a number between 0 and " + (siteCount - 1) + " inclusive.");
        } else if (end < 0) {
            JOptionPane.showMessageDialog(this.getParent(), "End site must be non negative.");
        } else if (end >= siteCount) {
            JOptionPane.showMessageDialog(this.getParent(), "End site must be less than " + siteCount + ".");
        } else if (isEndTextFieldNumeric == false) {
            JOptionPane.showMessageDialog(this.getParent(), "End site must be a number between 0 and " + (siteCount - 1) + " inclusive.");
        } else if (start > end) {
            JOptionPane.showMessageDialog(this.getParent(), "Start site must be less than the end site.");
        } else if (startPos < 0 && endPos >= 0 && theAlignment.getNumLoci() == 1) {
            JOptionPane.showMessageDialog(this.getParent(), "Start position must be non negative.");
        } else if (startPos > theAlignment.getPositionInLocus(siteCount - 1) && theAlignment.getNumLoci() == 1) {
            JOptionPane.showMessageDialog(this.getParent(), "No available SNPs with positions greater than " + theAlignment.getPositionInLocus(siteCount - 1) + ".");
        } else if (isStartPosTextFieldNumeric == false) {
            JOptionPane.showMessageDialog(this.getParent(), "Start position must be a number between 0 and " + theAlignment.getPositionInLocus(siteCount - 1) + " inclusive.");
        } else if (endPos < 0 && startPos >= 0 && theAlignment.getNumLoci() == 1) {
            JOptionPane.showMessageDialog(this.getParent(), "End Position must be non negative.");
        } else if (endPos < theAlignment.getPositionInLocus(0) && theAlignment.getNumLoci() == 1) {
            JOptionPane.showMessageDialog(this.getParent(), "No available SNPs with positions less than " + theAlignment.getPositionInLocus(0) + ".");
        } else if (isEndPosTextFieldNumeric == false) {
            JOptionPane.showMessageDialog(this.getParent(), "End position must be a number greater than " + theAlignment.getPositionInLocus(0) + ".");
        } else if (startPos > endPos && theAlignment.getNumLoci() == 1) {
            JOptionPane.showMessageDialog(this.getParent(), "Start position must be less than the end position.");
        } else if (!isChromSelectionValid) {
            JOptionPane.showMessageDialog(this.getParent(), "Invalid chromosome selection");
        } else {
            isCanceled = false;
            setVisible(false);
        }
    }

    public boolean isAllSiteIncluded() {
        return true;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public double getMinimumFrequency() {
        return minFreq;
    }

    public double getMaximumFrequency() {
        return maxFreq;
    }

    /**
     * For use when not doing batch... it is useful to be able to specify the
     * exact number of sequences which must be present.
     *
     * @see public double getMinPercentage()
     * @return
     */
    public int getMinimumCount() {
        if (doBatchAnalysis) {
            throw new RuntimeException("This method should not be called when using an Minimum Sequence Percentage");
        }
        return minCount;
    }

    /**
     * For use when doing batch... given that one cannot know the number of
     * sequences in advance of actual execution, one must use a relational
     * number to specify the number of sequences which must be present.
     *
     * @see public int getMinimumCount()
     * @return
     */
    public double getMinPercentage() {
        if (!doBatchAnalysis) {
            throw new RuntimeException("This method should not be called when using an absolute Minimum Sequence Count");
        }
        return minPercentage;
    }

    public boolean isExtractIndels() {
        //return indelCheckBox.isSelected();
        return false;
    }

    public boolean isRemoveMinorSNPs() {
        return removeMinorCheckBox.isSelected();
    }

    public boolean isUseSlidingWindow() {
        return slidingHapCheckBox.isSelected();
    }

    public int getWindowSize() {

        int windowSize = INVALID_VALUE;
        if (slidingHapCheckBox.isSelected()) {
            try {
                windowSize = Integer.parseInt(winSizeTextField.getText().trim());
            } catch (NumberFormatException nfe) {
                JOptionPane.showMessageDialog(this.getParent(), "Please enter an integer.");
            }
        }
        return windowSize;
    }

    public int getStepSize() {

        int stepSize = INVALID_VALUE;
        if (slidingHapCheckBox.isSelected()) {
            try {
                stepSize = Integer.parseInt(stepSizeTextField.getText().trim());
            } catch (NumberFormatException nfe) {
                JOptionPane.showMessageDialog(this.getParent(), "Please enter an integer.");
            }
        }
        return stepSize;
    }

    private void cancelButton_actionPerformed(ActionEvent e) {
        isCanceled = true;
        setVisible(false);
    }

    private void chromSelectButton_actionPerformed(ActionEvent e) {
        myChromFilter.setLocationRelativeTo(this);
        myChromFilter.setVisible(true);
        if (!myChromFilter.isCanceled()) {
            chromsSelected = myChromFilter.getChromsSelected();
            Alignment[] selectedAlignments = new Alignment[chromsSelected.length];
            List<Datum> availableAlignments = SeparatePlugin.separateAlignmentIntoLoci(theAlignment, null);
            for (int i = 0; i < chromsSelected.length; i++) {
                for (int j = 0; j < availableAlignments.size(); j++) {
                    Alignment current = (Alignment) availableAlignments.get(j).getData();
                    if (current.getLoci().length == 1) {
                        if (chromsSelected[i].equals(current.getLocusName(0))) {
                            selectedAlignments[i] = current;
                        }
                    }
                }
            }
            chromFilteredAlignment = CombineAlignment.getInstance(selectedAlignments);
            lblSitePos.setVisible(chromFilteredAlignment.getNumLoci() == 1 && startPos >= 0);
            startPosTextField.setVisible(chromFilteredAlignment.getNumLoci() == 1 && startPos >= 0);
            endPosTextField.setVisible(chromFilteredAlignment.getNumLoci() == 1 && endPos >= 0);
            totalSeq = chromFilteredAlignment.getSequenceCount();
            siteCount = chromFilteredAlignment.getSiteCount();
            lblSeqLength.setText(" of " + (siteCount - 1) + " sites");
            lblMinCount.setText("Minimum Count:");
            start = 0;
            end = siteCount - 1;
            startPos = chromFilteredAlignment.getPositionInLocus(0);
            endPos = chromFilteredAlignment.getPositionInLocus(siteCount - 1);
            if (doBatchAnalysis) {
                countTextField.setText(minPercentage + "");
            } else {
                countTextField.setText(minCount + "");
            }
            setEndTextField();
            startTextField.setText(start + "");
            endPosTextField.setText(endPos + "");
            startPosTextField.setText(startPos + "");
        }
    }

    public String[] getChromsSelected() {
        return chromsSelected;
    }

    public Alignment getChromFilteredAlignment() {
        return chromFilteredAlignment;
    }

    public boolean isCanceled() {
        return isCanceled;
    }

    private void setEndTextField() {
        this.endTextField.setText("" + end);
    }

    private void endTextField_focusLost(FocusEvent e) {

        try {
            int endFromField = Integer.parseInt(endTextField.getText().trim());
            if (endFromField < 0) {
                throw new IllegalArgumentException("End Site Can't be Negative.");
            }

            if (endFromField >= siteCount) {
                throw new IllegalArgumentException("End Site Can't be Greater Than: " + (siteCount - 1));
            }

            end = endFromField;
            endPos = theAlignment.getPositionInLocus(end);
            endPosTextField.setText(String.valueOf(endPos));
        } catch (Exception ee) {
            StringBuilder builder = new StringBuilder();
            builder.append("\nProblem with End Site: ");
            builder.append(endTextField.getText().trim());
            builder.append("\n");
            builder.append("Number Should be a Positive Integer between:\n0 and ");
            builder.append(siteCount - 1);
            builder.append("\n");
            builder.append(ee.getMessage());
            builder.append("\n");
            JOptionPane.showMessageDialog(this.getParent(), builder.toString(), "Error", JOptionPane.ERROR_MESSAGE);
            try {
                end = theAlignment.getSiteOfPhysicalPosition(endPos, null);
                endTextField.setText(String.valueOf(end));
            } catch (Exception ex) {
                ex.printStackTrace();
                // do nothing
            }
        }

    }

    private void startTextField_focusLost(FocusEvent e) {

        try {
            int startFromField = Integer.parseInt(startTextField.getText().trim());
            if (startFromField < 0) {
                throw new IllegalArgumentException("Start Site Can't be Negative.");
            }

            if (startFromField >= siteCount) {
                throw new IllegalArgumentException("Start Site Can't be Greater Than: " + (siteCount - 1));
            }

            start = startFromField;
            startPos = theAlignment.getPositionInLocus(start);
            startPosTextField.setText(String.valueOf(startPos));
        } catch (Exception ee) {
            StringBuilder builder = new StringBuilder();
            builder.append("\nProblem with Start Site: ");
            builder.append(startTextField.getText().trim());
            builder.append("\n");
            builder.append("Number Should be a Positive Integer between:\n0 and ");
            builder.append(siteCount - 1);
            builder.append("\n");
            builder.append(ee.getMessage());
            builder.append("\n");
            JOptionPane.showMessageDialog(this.getParent(), builder.toString(), "Error", JOptionPane.ERROR_MESSAGE);
            try {
                start = theAlignment.getSiteOfPhysicalPosition(startPos, null);
                startTextField.setText(String.valueOf(start));
            } catch (Exception ex) {
                ex.printStackTrace();
                // do nothing
            }
        }

    }

    private void endPosTextField_focusLost(FocusEvent e) {

        try {
            int endPosFromField = Integer.parseInt(endPosTextField.getText().trim());
            if (endPosFromField < 0) {
                throw new IllegalArgumentException("End Position Can't be Negative.");
            }

            int endSite = theAlignment.getSiteOfPhysicalPosition(endPosFromField, null);

            if (endSite < 0) {
                endSite = -(endSite + 1);
            }

            if (endSite >= siteCount) {
                endSite = siteCount - 1;
            }

            end = endSite;
            endPos = theAlignment.getPositionInLocus(end);
            endTextField.setText(String.valueOf(end));
        } catch (Exception ee) {
            ee.printStackTrace();
            StringBuilder builder = new StringBuilder();
            builder.append("\nProblem with End Physical Position: ");
            builder.append(endPosTextField.getText().trim());
            builder.append("\n");
            builder.append("Number Should be a Positive Integer between: \n");
            builder.append(theAlignment.getPositionInLocus(0));
            builder.append(" and ");
            builder.append(theAlignment.getPositionInLocus(siteCount - 1));
            builder.append("\n");
            builder.append(ee.getMessage());
            builder.append("\n");
            JOptionPane.showMessageDialog(this.getParent(), builder.toString(), "Error", JOptionPane.ERROR_MESSAGE);
            try {
                endPos = theAlignment.getPositionInLocus(end);
                endPosTextField.setText(String.valueOf(endPos));
            } catch (Exception ex) {
                ex.printStackTrace();
                // do nothing
            }
        }

    }

    private void startPosTextField_focusLost(FocusEvent e) {

        try {
            int startPosFromField = Integer.parseInt(startPosTextField.getText().trim());
            if (startPosFromField < 0) {
                throw new IllegalArgumentException("Start Position Can't be Negative.");
            }

            int startSite = theAlignment.getSiteOfPhysicalPosition(startPosFromField, null);

            if (startSite < 0) {
                startSite = -(startSite + 1);
            }

            if (startSite >= siteCount) {
                startSite = siteCount - 1;
            }

            start = startSite;
            startPos = theAlignment.getPositionInLocus(start);
            startTextField.setText(String.valueOf(start));
        } catch (Exception ee) {
            StringBuilder builder = new StringBuilder();
            builder.append("\nProblem with Start Physical Position: ");
            builder.append(startPosTextField.getText().trim());
            builder.append("\n");
            builder.append("Number Should be a Positive Integer between: \n");
            builder.append(theAlignment.getPositionInLocus(0));
            builder.append(" and ");
            builder.append(theAlignment.getPositionInLocus(siteCount - 1));
            builder.append("\n");
            builder.append(ee.getMessage());
            builder.append("\n");
            JOptionPane.showMessageDialog(this.getParent(), builder.toString(), "Error", JOptionPane.ERROR_MESSAGE);
            try {
                startPos = theAlignment.getPositionInLocus(start);
                startPosTextField.setText(String.valueOf(startPos));
            } catch (Exception ex) {
                ex.printStackTrace();
                // do nothing
            }
        }

    }

    private void freqTextField_focusLost(FocusEvent e) {
        try {
            String input = freqTextField.getText().trim();
            double tmpMinFreq = -0.1;
            if (input != null) {
                tmpMinFreq = Double.parseDouble(input);
            }
            if ((tmpMinFreq > 1.0) || (tmpMinFreq < 0.0)) {
                tmpMinFreq = minFreq;
            }
            minFreq = tmpMinFreq;

        } catch (NumberFormatException nfe) {
            JOptionPane.showMessageDialog(this.getParent(), "Could not parse \"Minimum Frequency\".  "
                    + "Please enter a value between 0.0 and 1.0");
        } catch (Exception ee) {
            ee.printStackTrace();
        }
        freqTextField.setText(minFreq + "");
        // save this value to the users settings so it shows up the next time
        TasselPrefs.putFilterAlignPluginMinFreq(minFreq);
    }

    private void maxFreqTextField_focusLost(FocusEvent e) {
        try {
            String input = maxFreqTextField.getText().trim();
            double tmpMaxFreq = -0.1;
            if (input != null) {
                tmpMaxFreq = Double.parseDouble(input);
            }
            if ((tmpMaxFreq > 1.0) || (tmpMaxFreq < 0.0)) {
                tmpMaxFreq = maxFreq;
            }
            maxFreq = tmpMaxFreq;

        } catch (NumberFormatException nfe) {
            JOptionPane.showMessageDialog(this.getParent(), "Could not parse \"Maximum Frequency\".  "
                    + "Please enter a value between 0.0 and 1.0");
        } catch (Exception ee) {
            ee.printStackTrace();
        }
        maxFreqTextField.setText(maxFreq + "");
        // save this value to the users settings so it shows up the next time
        TasselPrefs.putFilterAlignPluginMaxFreq(maxFreq);
    }

    private void countTextField_focusLost(FocusEvent e) {
        if (doBatchAnalysis) {
            double minPercentageOriginal = minPercentage;
            try {
                minPercentage = Double.parseDouble(countTextField.getText().trim());
                if (minPercentage > 1 || minPercentage < 0) {
                    minPercentage = minPercentageOriginal;
                }
            } catch (Exception ee) {
                ee.printStackTrace();
                minPercentage = minPercentageOriginal;
            }
            countTextField.setText(minPercentage + "");
        } else {
            int minCountOriginal = minCount;
            try {
                minCount = Integer.parseInt(countTextField.getText().trim());
                if ((minCount > theAlignment.getSequenceCount()) || (minCount < 0)) {
                    minCount = minCountOriginal;
                }
            } catch (Exception ee) {
                ee.printStackTrace();
                minCount = minCountOriginal;
            }
            countTextField.setText(minCount + "");
            TasselPrefs.putFilterAlignPluginMinCount(minCount);
        }
    }

    public Dimension getMinimumSize() {
        return new Dimension(600, 600);
    }
}

class ChromosomeFilterDialog extends JDialog {

    private int numChromsSelected;
    private boolean isCanceled = true;
    private JPanel mainPanel = new JPanel();
    private JPanel checkBoxPanel = new JPanel();
    private JLabel lblChromSelect = new JLabel();
    private JButton okayButton = new JButton();
    private JButton cancelButton = new JButton();
    private JCheckBox selectAllCheckBox = new JCheckBox();
    private JCheckBox[] selectChromsCheckBoxes;
    private GridBagLayout gridBagLayout = new GridBagLayout();

    public ChromosomeFilterDialog(String[] chromNames, Frame f) {
        super(f, "Filter Alignment", true);
        selectChromsCheckBoxes = new JCheckBox[chromNames.length];
        selectAllCheckBox.setText("Select/Deselect All");
        selectAllCheckBox.setSelected(true);
        for (int i = 0; i < chromNames.length; i++) {
            selectChromsCheckBoxes[i] = new JCheckBox(chromNames[i]);
            selectChromsCheckBoxes[i].setSelected(true);
        }
        numChromsSelected = chromNames.length;
        try {
            initUI();
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    void initUI() throws Exception {

        lblChromSelect.setFont(new java.awt.Font("Dialog", 1, 16));
        lblChromSelect.setText("Select Chromosomes to Filter");

        mainPanel.setMinimumSize(new Dimension(480, 480));
        mainPanel.setPreferredSize(new Dimension(480, 480));
        mainPanel.setLayout(gridBagLayout);

        okayButton.setMaximumSize(new Dimension(63, 27));
        okayButton.setMinimumSize(new Dimension(63, 27));
        okayButton.setText("Select");
        okayButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                okayButton_actionPerformed(e);
            }
        });

        cancelButton.setMaximumSize(new Dimension(63, 27));
        cancelButton.setMinimumSize(new Dimension(63, 27));
        cancelButton.setText("Cancel");
        cancelButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                cancelButton_actionPerformed(e);
            }
        });

        JPanel buttonPanel = new JPanel();
        buttonPanel.add(okayButton);
        buttonPanel.add(cancelButton);

        selectAllCheckBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                selectAllCheckBox_actionPerformed(e);
            }
        });

        for (int i = 0; i < selectChromsCheckBoxes.length; i++) {
            selectChromsCheckBoxes[i].addActionListener(new java.awt.event.ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    checkBox_actionPerformed(e);
                }
            });
        }
        checkBoxPanel.setLayout(new GridLayout(0, 2));
        for (int i = 0; i < selectChromsCheckBoxes.length; i++) {
            checkBoxPanel.add(selectChromsCheckBoxes[i]);
        }
        checkBoxPanel.add(selectAllCheckBox);

        mainPanel.add(lblChromSelect, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 0, 5, 0), 0, 0));
        mainPanel.add(checkBoxPanel, new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(20, 40, 5, 0), 15, 0));
        mainPanel.add(buttonPanel, new GridBagConstraints(0, 2, 1, 1, 1.0, 1.0, GridBagConstraints.PAGE_END, GridBagConstraints.NONE, new Insets(5, 0, 5, 0), 0, 0));

        this.add(mainPanel, BorderLayout.CENTER);
    }

    public String[] getChromsSelected() {
        String[] chroms = new String[numChromsSelected];
        int j = 0;
        for (int i = 0; i < selectChromsCheckBoxes.length; i++) {
            if (selectChromsCheckBoxes[i].isSelected()) {
                chroms[j] = selectChromsCheckBoxes[i].getText();
                j++;
            }
        }
        return chroms;
    }

    private void okayButton_actionPerformed(ActionEvent e) {
        if (numChromsSelected != 0) {
            isCanceled = false;
            setVisible(false);
        } else {
            JOptionPane.showMessageDialog(this.getParent(), "Please select at least one chromosome.");
        }
    }

    private void cancelButton_actionPerformed(ActionEvent e) {
        isCanceled = true;
        setVisible(false);
    }

    private void selectAllCheckBox_actionPerformed(ActionEvent e) {
        if (((JCheckBox) e.getSource()).isSelected()) {
            for (int i = 0; i < selectChromsCheckBoxes.length; i++) {
                selectChromsCheckBoxes[i].setSelected(true);
            }
            numChromsSelected = selectChromsCheckBoxes.length;
        } else {
            for (int i = 0; i < selectChromsCheckBoxes.length; i++) {
                selectChromsCheckBoxes[i].setSelected(false);
            }
            numChromsSelected = 0;
        }
    }

    private void checkBox_actionPerformed(ActionEvent e) {
        if (!((JCheckBox) e.getSource()).isSelected()) {
            selectAllCheckBox.setSelected(false);
            numChromsSelected--;
        } else {
            numChromsSelected++;
            if (numChromsSelected == selectChromsCheckBoxes.length) {
                selectAllCheckBox.setSelected(true);
            }
        }
    }

    public boolean isCanceled() {
        return isCanceled;
    }
}
