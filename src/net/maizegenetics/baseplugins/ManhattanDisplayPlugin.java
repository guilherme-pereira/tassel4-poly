/*
 * ManhattanDisplayPlugin
 */
package net.maizegenetics.baseplugins;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.net.URL;
import java.util.ArrayList;
import javax.swing.ImageIcon;
import net.maizegenetics.plugindef.DataSet;
import java.util.List;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import net.maizegenetics.baseplugins.chart.XYScatterMultipleYPanel;
import net.maizegenetics.pal.report.TableReport;
import net.maizegenetics.plugindef.Datum;

/**
 *
 * @author yz79
 */
public class ManhattanDisplayPlugin extends AbstractDisplayPlugin {

    /**
     * Creates a new instance of ManhattanDisplayPlugin
     */
    public ManhattanDisplayPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {
        try {
            List<Datum> tableInList = input.getDataOfType(TableReport.class);
            if (tableInList.size() != 1) {
                String message = "Invalid selection.  Please select one table result.";
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), message);
                } else {
                    System.out.println(message);
                }
                return null;
            }
            TableReport myTableReport = (TableReport) tableInList.get(0).getData();
            if (isInteractive()) {
                try {
                    ArrayList<Integer> indexes = splitTable(myTableReport);
                    String[] traits = getTraits(myTableReport);
                    if (traits.length > 1) {
                        PlotOptionsDialog myOptions = new PlotOptionsDialog(this.getParentFrame(), getTraits(myTableReport));
                        myOptions.setLocationRelativeTo(getParentFrame());
                        myOptions.setVisible(true);
                        if (myOptions.isCanceled() == false) {
                            int index = myOptions.getTraitIndex();
                            ManhattanDisplayPluginDialog myDialog = new ManhattanDisplayPluginDialog(this.getParentFrame(), this, myTableReport, indexes.get((index - 1) * 2), indexes.get((index - 1) * 2 + 1));
                            myDialog.setLocationRelativeTo(getParentFrame());
                            myDialog.setVisible(true);
                        }
                    } else if (traits.length == 1) {
                        ManhattanDisplayPluginDialog myDialog = new ManhattanDisplayPluginDialog(this.getParentFrame(), this, myTableReport, indexes.get(0), indexes.get(1));
                        myDialog.setLocationRelativeTo(getParentFrame());
                        myDialog.setVisible(true);
                    }
                } catch (Exception ex) {
                    ex.printStackTrace();
                    JOptionPane.showMessageDialog(this.getParentFrame(), "Unable to create Manhattan plot " + ex);
                } catch (Error er) {
                    er.printStackTrace();
                    JOptionPane.showMessageDialog(this.getParentFrame(), "Unable to create Manhattan plot " + er);
                }
            }
            return null;
        } finally {
            fireProgress(100);
        }
    }

    private ArrayList<Integer> splitTable(TableReport table) {
        ArrayList<Integer> indexes = new ArrayList<Integer>();
        int numRows = table.getRowCount();
        String previousTrait = "";
        for (int i = 0; i < numRows; i++) {
            if (!previousTrait.equals((String) table.getValueAt(i, 0))) {
                if (!((String) table.getValueAt(i, 1)).equals("None")) {
                    indexes.add(new Integer(i));
                    previousTrait = (String) table.getValueAt(i, 0);
                    if (i > 1) {
                        indexes.add(new Integer(i));
                    }
                } else if (i != 0) {
                    indexes.add(new Integer(i));
                    indexes.add(new Integer(i + 1));
                    previousTrait = (String) table.getValueAt(i + 1, 0);
                }
            }
        }
        indexes.add(new Integer(numRows));
        return indexes;
    }

    private String[] getTraits(TableReport table) {
        ArrayList<String> traitArray = new ArrayList<String>();
        int numRows = table.getRowCount();
        String previousTrait = "";
        for (int i = 0; i < numRows; i++) {
            if (!previousTrait.equals((String) table.getValueAt(i, 0))) {
                previousTrait = (String) table.getValueAt(i, 0);
                traitArray.add(previousTrait);
            }
        }

        String[] traits = new String[traitArray.size()];
        for (int i = 0; i < traitArray.size(); i++) {
            traits[i] = traitArray.get(i);
        }
        return traits;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    @Override
    public ImageIcon getIcon() {
        URL imageURL = QQDisplayPlugin.class.getResource("images/ManhattanPlot.gif");
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
        return "Manhattan Plot";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    @Override
    public String getToolTipText() {
        return "Display Manhattan Plot";
    }
}

class ManhattanDisplayPluginDialog extends JDialog {

    XYScatterMultipleYPanel myManhattanPlot;
    TableReport myTableReport;
    JButton myCloseButton = new JButton();
    JPanel myMainPanel;
    JPanel myOptionPanel = new JPanel();

    public ManhattanDisplayPluginDialog(Frame f, ManhattanDisplayPlugin plugin, TableReport theTableReport, int start, int end) {
        super(f, "Manhattan Plot", false);
        myTableReport = theTableReport;
        try {
            jbInit();
            myManhattanPlot = new XYScatterMultipleYPanel(plugin, theTableReport, start, end);
//            myQQFigurePanel = new QQComponent(theTableReport);
            getContentPane().add(myManhattanPlot, BorderLayout.CENTER);
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
            JOptionPane.showMessageDialog(this.getParent(), "Unable to create Manhattan plot " + ex);
        }
        repaint();
    }

    void jbInit() throws Exception {
        myCloseButton.setText("Close");
        myCloseButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                closeButton_actionPerformed(e);
            }
        });

        myOptionPanel.add(myCloseButton, new GridBagConstraints(0, 0, 0, 0, 0.0, 0.0, GridBagConstraints.SOUTHEAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        getContentPane().add(myOptionPanel);
    }

    void closeButton_actionPerformed(ActionEvent e) {
        dispose();
    }
}

class PlotOptionsDialog extends JDialog {

    boolean isCanceled = true;
    private final static int TEXT_FIELD_WIDTH = 8;
    private final static int INVALID_VALUE = -999;   // what is returned when inappropriate to return any value
    private JPanel mainPanel = new JPanel();
    private JButton okayButton = new JButton();
    private JButton cancelButton = new JButton();
    private JComboBox traitList;
    private JPanel checkBoxPanel = new JPanel();
    private GridBagLayout gridBagLayout2 = new GridBagLayout();

    public PlotOptionsDialog(Frame f, String[] traits) {
        super(f, "Manhattan Plot Options", true);

        String[] dropDownList = new String[traits.length + 1];
        dropDownList[0] = "Select trait";
        for (int i = 1; i < dropDownList.length; i++) {
            dropDownList[i] = traits[i - 1];
        }
        traitList = new JComboBox(dropDownList);

        try {
            initUI();
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    void initUI() throws Exception {

        // mainPanel.setBackground(SystemColor.menu);
        mainPanel.setMinimumSize(new Dimension(300, 80));
        mainPanel.setPreferredSize(new Dimension(300, 80));
        mainPanel.setLayout(gridBagLayout2);

        GridBagConstraints c = new GridBagConstraints();

        traitList.setSelectedIndex(0);
        traitList.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                traitSelector_actionPerformed(e);
            }
        });
        c.gridx = 0;
        c.gridy = 0;
        c.fill = GridBagConstraints.HORIZONTAL;
        c.gridwidth = 2;
        mainPanel.add(traitList, c);

        okayButton.setMaximumSize(new Dimension(63, 27));
        okayButton.setMinimumSize(new Dimension(63, 27));
        okayButton.setText("Okay");
        okayButton.setEnabled(false);
        okayButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                okayButton_actionPerformed(e);
            }
        });
        c.gridx = 0;
        c.gridy = 1;
        c.gridwidth = 1;
        mainPanel.add(okayButton, c);

        cancelButton.setMaximumSize(new Dimension(63, 27));
        cancelButton.setMinimumSize(new Dimension(63, 27));
        cancelButton.setText("Cancel");
        cancelButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                cancelButton_actionPerformed(e);
            }
        });
        c.gridx = 1;
        c.gridy = 1;
        mainPanel.add(cancelButton, c);

        this.add(mainPanel, BorderLayout.CENTER);
    }

    private void okayButton_actionPerformed(ActionEvent e) {
        isCanceled = false;
        setVisible(false);
    }

    private void cancelButton_actionPerformed(ActionEvent e) {
        isCanceled = true;
        setVisible(false);
    }

    private void traitSelector_actionPerformed(ActionEvent e) {
        if (traitList.getSelectedIndex() != 0) {
            okayButton.setEnabled(true);
        } else {
            okayButton.setEnabled(false);
        }
    }

    public String getTrait() {
        return (String) traitList.getSelectedItem();
    }

    public int getTraitIndex() {
        return traitList.getSelectedIndex();
    }

    public boolean isCanceled() {
        return isCanceled;
    }

    public Dimension getMinimumSize() {
        return new Dimension(600, 600);
    }
}