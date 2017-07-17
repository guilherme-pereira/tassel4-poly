/*
 * FilterTaxaPropertiesPlugin.java
 *
 * Created on May 9, 2013
 *
 */
package net.maizegenetics.baseplugins;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;

import javax.swing.*;
import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;

import java.net.URL;

import java.util.ArrayList;
import java.util.List;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.prefs.TasselPrefs;

import org.apache.log4j.Logger;

/**
 *
 * @author Terry
 */
public class FilterTaxaPropertiesPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(FilterTaxaPropertiesPlugin.class);
    private double myMinNotMissing = TasselPrefs.FILTER_TAXA_PROPS_PLUGIN_MIN_NOT_MISSING_DEFAULT;
    private double myMinHeterozygous = TasselPrefs.FILTER_TAXA_PROPS_PLUGIN_MIN_HET_DEFAULT;
    private double myMaxHeterozygous = TasselPrefs.FILTER_TAXA_PROPS_PLUGIN_MAX_HET_DEFAULT;

    /**
     * Creates a new instance of FilterTaxaPropertiesPlugin
     */
    public FilterTaxaPropertiesPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet performFunction(DataSet input) {

        try {

            List<Datum> alignInList = input.getDataOfType(Alignment.class);

            if (alignInList.size() != 1) {
                String gpMessage = "Invalid selection.  Please Select One Genotype Alignment.";
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), gpMessage);
                } else {
                    myLogger.error(gpMessage);
                }
                return null;
            }
            Datum current = alignInList.get(0);
            Datum result = processDatum(current, isInteractive());

            if (result != null) {
                DataSet output = new DataSet(result, this);
                fireDataSetReturned(new PluginEvent(output, FilterTaxaPropertiesPlugin.class));
                return output;
            } else {
                return null;
            }

        } finally {
            fireProgress(100);
        }

    }

    private Datum processDatum(Datum inDatum, boolean isInteractive) {
        Alignment aa = (Alignment) inDatum.getData();

        if (isInteractive) {
            FilterTaxaPropertiesDialog theDialog = new FilterTaxaPropertiesDialog();
            theDialog.setLocationRelativeTo(getParentFrame());
            theDialog.setVisible(true);
            if (theDialog.isCancel()) {
                return null;
            }

            myMinNotMissing = theDialog.getMinNotMissingProportion();
            myMinHeterozygous = theDialog.getMinHeterozygousProportion();
            myMaxHeterozygous = theDialog.getMaxHeterozygousProportion();

            theDialog.dispose();
        }

        Alignment result = getFilteredAlignment(aa);

        StringBuilder builder = new StringBuilder();
        builder.append("Filter Alignment by Taxa Properties...\n");
        builder.append("   Min. Proportion of Sites Present: ");
        builder.append(myMinNotMissing);
        builder.append("\n");
        builder.append("   Heterozygous Proportion: ");
        builder.append(myMinHeterozygous);
        builder.append(" - ");
        builder.append(myMaxHeterozygous);
        builder.append("\n");
        String theComment = builder.toString();

        if (result == aa) {
            if (isInteractive()) {
                JOptionPane.showMessageDialog(getParentFrame(), "The Alignment is Unchanged.");
            } else {
                myLogger.warn("The Alignment is Unchanged: " + inDatum.getName());
            }
        }

        if (result.getSequenceCount() != 0) {
            String theName = inDatum.getName() + "_" + result.getSequenceCount() + "Taxa";
            myLogger.info("Resulting Number Sequences: " + result.getSequenceCount());
            return new Datum(theName, result, theComment);
        } else {
            if (isInteractive()) {
                JOptionPane.showMessageDialog(getParentFrame(), "No remaining Taxa given filter parameters.");
            } else {
                myLogger.warn("No remaining Taxa given filter parameters.");

            }
            return null;
        }
    }

    private Alignment getFilteredAlignment(Alignment alignment) {
        int numSites = alignment.getSiteCount();
        int numTaxa = alignment.getSequenceCount();
        IdGroup ids = alignment.getIdGroup();

        List<Identifier> keepTaxaList = new ArrayList<Identifier>();
        for (int t = 0; t < numTaxa; t++) {

            progress((int) ((double) t / (double) numTaxa * 100.0), null);

            if (myMinNotMissing != 0.0) {
                int totalNotMissing = alignment.getTotalNotMissingForTaxon(t);
                double percentNotMissing = (double) totalNotMissing / (double) numSites;
                if (percentNotMissing < myMinNotMissing) {
                    continue;
                }
            }

            if ((myMinHeterozygous != 0.0) || (myMaxHeterozygous != 1.0)) {
                int numHeterozygous = alignment.getHeterozygousCountForTaxon(t);
                int totalSitesNotMissing = alignment.getTotalNotMissingForTaxon(t);
                double percentHets = (double) numHeterozygous / (double) totalSitesNotMissing;
                if ((percentHets < myMinHeterozygous) || (percentHets > myMaxHeterozygous)) {
                    continue;
                }
            }

            keepTaxaList.add(ids.getIdentifier(t));
        }

        IdGroup taxa = new SimpleIdGroup(keepTaxaList);
        return FilterAlignment.getInstance(alignment, taxa, false);
    }

    public double getMinNotMissing() {
        return myMinNotMissing;
    }

    public void setMinNotMissing(int minNotMissing) {
        if ((minNotMissing < 0.0) || (minNotMissing > 1.0)) {
            throw new IllegalArgumentException("FilterTaxaPropertiesPlugin: setMinNotMissing: Value must be between 0.0 and 1.0: " + minNotMissing);
        }
        myMinNotMissing = minNotMissing;
    }

    public double getMinHeterozygous() {
        return myMinHeterozygous;
    }

    public void setMinHeterozygous(double minHeterozygous) {
        if ((minHeterozygous < 0.0) || (minHeterozygous > 1.0)) {
            throw new IllegalArgumentException("FilterTaxaPropertiesPlugin: setMinHeterozygous: Value must be between 0.0 and 1.0: " + minHeterozygous);
        }
        myMinHeterozygous = minHeterozygous;
    }

    public double getMaxHeterozygous() {
        return myMaxHeterozygous;
    }

    public void setMaxHeterozygous(double maxHeterozygous) {
        if ((maxHeterozygous < 0.0) || (maxHeterozygous > 1.0)) {
            throw new IllegalArgumentException("FilterTaxaPropertiesPlugin: setMaxHeterozygous: Value must be between 0.0 and 1.0: " + maxHeterozygous);
        }
        myMaxHeterozygous = maxHeterozygous;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        URL imageURL = FilterTaxaPropertiesPlugin.class.getResource("images/Filter_horizontal.gif");
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
        return "Taxa";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Filter Alignment Based Taxa Properties";
    }
}

class FilterTaxaPropertiesDialog extends JDialog {

    private static final int TEXT_FIELD_WIDTH = 10;
    private JTabbedPane myTabbedPane = new JTabbedPane();
    private JTextField myMinNotMissingField = new JTextField(TEXT_FIELD_WIDTH);
    private JTextField myMinHeterozygousField = new JTextField(TEXT_FIELD_WIDTH);
    private JTextField myMaxHeterozygousField = new JTextField(TEXT_FIELD_WIDTH);
    private double myMinNotMissing = TasselPrefs.getFilterTaxaPropsMinNotMissingFreq();
    private double myMinHeterozygous = TasselPrefs.getFilterTaxaPropsMinHetFreq();
    private double myMaxHeterozygous = TasselPrefs.getFilterTaxaPropsMaxHetFreq();
    private boolean myIsCancel = true;
    private boolean myIsOkToContinue = true;

    public FilterTaxaPropertiesDialog() {
        super((Frame) null, null, true);

        JButton okButton = new JButton();
        okButton.setActionCommand("Ok");
        okButton.setText("Ok");
        okButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                if (myIsOkToContinue) {
                    myIsCancel = false;
                    setVisible(false);
                } else {
                    myIsOkToContinue = true;
                }
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

        panel.add(getLine("Min Proportion of Sites Present", myMinNotMissingField));
        panel.add(getLine("Min Heterozygous Proportion", myMinHeterozygousField));
        panel.add(getLine("Max Heterozygous Proportion", myMaxHeterozygousField));

        myMinNotMissingField.setText(String.valueOf(myMinNotMissing));
        myMinNotMissingField.addFocusListener(new FocusAdapter() {
            @Override
            public void focusLost(FocusEvent e) {
                myMinNotMissing = focusLostParseFreq(myMinNotMissingField, myMinNotMissing);
                TasselPrefs.putFilterTaxaPropsMinNotMissingFreq(myMinNotMissing);
            }
        });

        myMinHeterozygousField.setText(String.valueOf(myMinHeterozygous));
        myMinHeterozygousField.addFocusListener(new FocusAdapter() {
            @Override
            public void focusLost(FocusEvent e) {
                myMinHeterozygous = focusLostParseFreq(myMinHeterozygousField, myMinHeterozygous);
                TasselPrefs.putFilterTaxaPropsMinHetFreq(myMinHeterozygous);
            }
        });

        myMaxHeterozygousField.setText(String.valueOf(myMaxHeterozygous));
        myMaxHeterozygousField.addFocusListener(new FocusAdapter() {
            @Override
            public void focusLost(FocusEvent e) {
                myMaxHeterozygous = focusLostParseFreq(myMaxHeterozygousField, myMaxHeterozygous);
                TasselPrefs.putFilterTaxaPropsMaxHetFreq(myMaxHeterozygous);
            }
        });

        myTabbedPane.add(panel, "Filter Taxa by Properties");

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

    public double getMinNotMissingProportion() {
        return myMinNotMissing;
    }

    public double getMinHeterozygousProportion() {
        return myMinHeterozygous;
    }

    public double getMaxHeterozygousProportion() {
        return myMaxHeterozygous;
    }

    private double focusLostParseFreq(JTextField field, double lastValue) {

        myIsOkToContinue = true;

        double freq = -0.1;
        try {
            String input = field.getText().trim();

            if (input != null) {
                freq = Double.parseDouble(input);
            }
            if ((freq > 1.0) || (freq < 0.0)) {
                myIsOkToContinue = false;
                JOptionPane.showMessageDialog(this, "Please enter a value between 0.0 and 1.0");
                freq = lastValue;
            }

        } catch (NumberFormatException nfe) {
            myIsOkToContinue = false;
            JOptionPane.showMessageDialog(this, "Please enter a value between 0.0 and 1.0");
            freq = lastValue;
        }

        field.setText(String.valueOf(freq));
        return freq;
    }

    private JPanel getLine(String label, JTextField ref) {

        JPanel result = new JPanel(new FlowLayout(FlowLayout.RIGHT));

        result.add(new JLabel(label));
        ref.setEditable(true);
        ref.setHorizontalAlignment(JTextField.LEFT);
        ref.setAlignmentX(JTextField.CENTER_ALIGNMENT);
        ref.setAlignmentY(JTextField.CENTER_ALIGNMENT);
        ref.setMaximumSize(ref.getPreferredSize());
        result.add(ref);

        return result;

    }
}
