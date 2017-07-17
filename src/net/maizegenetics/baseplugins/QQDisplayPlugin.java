/*
 * QQDisplayPlugin.java
 *
 * Created on December 13, 2010
 *
 */
package net.maizegenetics.baseplugins;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.net.URL;
import java.util.ArrayList;
import javax.swing.ImageIcon;
import net.maizegenetics.plugindef.DataSet;
import java.util.List;
import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.JTextField;
import javax.swing.ListModel;
import javax.swing.ListSelectionModel;
import javax.swing.event.ChangeEvent;
import net.maizegenetics.baseplugins.chart.XYScatterAndLinePanel;
import net.maizegenetics.pal.report.TableReport;
import net.maizegenetics.plugindef.Datum;

/**
 *
 * @author yz79
 */
public class QQDisplayPlugin extends AbstractDisplayPlugin {

    /**
     * Creates a new instance of QQDisplayPlugin
     */
    public QQDisplayPlugin(Frame parentFrame, boolean isInteractive) {
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
                    PlotOptionsQQDialog myOptions = new PlotOptionsQQDialog(this.getParentFrame(), getTraits(myTableReport), splitTable(myTableReport));
                    myOptions.setLocationRelativeTo(getParentFrame());
                    myOptions.setVisible(true);
                    if (myOptions.isCanceled() == false) {
                        QQDisplayPluginDialog myDialog = new QQDisplayPluginDialog(this.getParentFrame(), this, myTableReport, myOptions.getSliderValue(), splitTable(myTableReport), myOptions.getTraitIndices());
                        myDialog.setLocationRelativeTo(getParentFrame());
                        myDialog.setVisible(true);
                    }
                } catch (Exception ex) {
                    ex.printStackTrace();
                    JOptionPane.showMessageDialog(this.getParentFrame(), "Unable to create QQ plot " + ex);
                } catch (Error er) {
                    er.printStackTrace();
                    JOptionPane.showMessageDialog(this.getParentFrame(), "Unable to create QQ plot " + er);
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
        URL imageURL = QQDisplayPlugin.class.getResource("images/QQPlot.gif");
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
        return "QQ Plot";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    @Override
    public String getToolTipText() {
        return "Display QQ Plot";
    }
}

class QQDisplayPluginDialog extends JDialog {

    XYScatterAndLinePanel myQQPlot;
    TableReport myTableReport;
    JButton myCloseButton = new JButton();
    JPanel myMainPanel;
    JPanel myOptionPanel = new JPanel();

    public QQDisplayPluginDialog(Frame f, QQDisplayPlugin plugin, TableReport theTableReport, int countToDisplay, ArrayList<Integer> tableIndices, int[] indices) {
        super(f, "QQ Plot", false);
        myTableReport = theTableReport;
        try {
            jbInit();
            myQQPlot = new XYScatterAndLinePanel(plugin, theTableReport, countToDisplay, tableIndices, indices);
//            myQQFigurePanel = new QQComponent(theTableReport);
            getContentPane().add(myQQPlot, BorderLayout.CENTER);
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
            JOptionPane.showMessageDialog(this.getParent(), "Unable to create QQ plot " + ex);
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

class PlotOptionsQQDialog extends JDialog {

    boolean isCanceled = true;
    String[] myTraits;
    private JPanel mainPanel = new JPanel();
    private JButton okayButton = new JButton();
    private JButton cancelButton = new JButton();
    private JLabel sliderLabel1 = new JLabel();
    private JSlider slider = new JSlider();
    private JLabel countLabel1 = new JLabel();
    private JLabel listLabel1 = new JLabel();
    private JLabel listLabel2 = new JLabel();
    private JTextField countTextField = new JTextField();
    private GridBagLayout gridBagLayout2 = new GridBagLayout();
    private JList list1 = new JList();
    private JList list2 = new JList();
    private JButton addAllButton = new JButton();
    private JButton addOneButton = new JButton();
    private JButton removeAllButton = new JButton();
    private JButton removeOneButton = new JButton();

    public PlotOptionsQQDialog(Frame f, String[] traits, ArrayList<Integer> indexes) {
        super(f, "QQ Plot Options", true);

        myTraits = traits;

        int numSites = indexes.get(1) - indexes.get(0);
        slider.setMinimum(1);
        slider.setMaximum(numSites);
        slider.setValue((int) (numSites * 0.01));
        countTextField.setText("" + (int) (numSites * 0.01));

        try {
            initUI();
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    void initUI() throws Exception {

        mainPanel.setMinimumSize(new Dimension(400, 230));
        mainPanel.setPreferredSize(new Dimension(400, 230));
        mainPanel.setLayout(gridBagLayout2);

        listLabel1.setText("Available Traits");
        listLabel2.setText("Traits to Plot");

        list1.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
        list1.setLayoutOrientation(JList.VERTICAL);
        list1.setVisibleRowCount(5);
        list1.setMinimumSize(new Dimension(100, 104));
        list1.setBackground(Color.white);
        list1.setBorder(BorderFactory.createLineBorder(Color.black));

        list2.setListData(myTraits);
        list2.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
        list2.setLayoutOrientation(JList.VERTICAL);
        list2.setVisibleRowCount(5);
        list2.setMinimumSize(new Dimension(100, 104));
        list2.setBorder(BorderFactory.createLineBorder(Color.black));

        addAllButton.setText(">>");
        addAllButton.setMaximumSize(new Dimension(63, 27));
        addAllButton.setMinimumSize(new Dimension(63, 27));
        addAllButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                addAllButton_actionPerformed(e);
            }
        });

        addOneButton.setText(">");
        addOneButton.setMaximumSize(new Dimension(63, 27));
        addOneButton.setMinimumSize(new Dimension(63, 27));
        addOneButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                addOneButton_actionPerformed(e);
            }
        });

        removeOneButton.setText("<");
        removeOneButton.setMaximumSize(new Dimension(63, 27));
        removeOneButton.setMinimumSize(new Dimension(63, 27));
        removeOneButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                removeOneButton_actionPerformed(e);
            }
        });

        removeAllButton.setText("<<");
        removeAllButton.setMaximumSize(new Dimension(63, 27));
        removeAllButton.setMinimumSize(new Dimension(63, 27));
        removeAllButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                removeAllButton_actionPerformed(e);
            }
        });

        sliderLabel1.setText("Plot Density:");
        sliderLabel1.setToolTipText("Plot Density is then number of Significant Points plotted without removing data. The rest of the data is plotted at a regular interval to maintain the trend without drawing all the points.");

        slider.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(ChangeEvent ce) {
                slider_actionPerformed(ce);
            }
        });

        if (myTraits.length == 1) {
            mainPanel.setMinimumSize(new Dimension(400, 100));
            mainPanel.setPreferredSize(new Dimension(400, 100));
            listLabel1.setVisible(false);
            listLabel2.setVisible(false);
            list1.setVisible(false);
            list2.setVisible(false);
            addAllButton.setVisible(false);
            addOneButton.setVisible(false);
            removeOneButton.setVisible(false);
            removeAllButton.setVisible(false);
        }

        countTextField.addKeyListener(new java.awt.event.KeyListener() {
            public void keyTyped(KeyEvent ke) {
//                countTextField_keyTyped(ke);
            }

            public void keyPressed(KeyEvent ke) {
//                throw new UnsupportedOperationException("Not supported yet.");
            }

            public void keyReleased(KeyEvent ke) {
                countTextField_keyTyped(ke);
            }
        });

        countLabel1.setText("of " + slider.getMaximum() + " points per trait.");

        okayButton.setMaximumSize(new Dimension(95, 27));
        okayButton.setMinimumSize(new Dimension(95, 27));
        okayButton.setText("Okay");
        okayButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                okayButton_actionPerformed(e);
            }
        });

        cancelButton.setMaximumSize(new Dimension(95, 27));
        cancelButton.setMinimumSize(new Dimension(95, 27));
        cancelButton.setText("Cancel");
        cancelButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                cancelButton_actionPerformed(e);
            }
        });

        mainPanel.add(listLabel1, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(10, 0, 5, 0), 0, 0));
        mainPanel.add(listLabel2, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(10, 0, 5, 0), 0, 0));
        mainPanel.add(list1, new GridBagConstraints(0, 1, 1, 4, 1.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 5, 0, 0), 0, 0));
        mainPanel.add(list2, new GridBagConstraints(2, 1, 1, 4, 1.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 5), 0, 0));
        mainPanel.add(addAllButton, new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
        mainPanel.add(addOneButton, new GridBagConstraints(1, 2, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
        mainPanel.add(removeOneButton, new GridBagConstraints(1, 3, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
        mainPanel.add(removeAllButton, new GridBagConstraints(1, 4, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
        mainPanel.add(sliderLabel1, new GridBagConstraints(0, 5, 1, 2, 0.0, 0.0, GridBagConstraints.LINE_END, GridBagConstraints.NONE, new Insets(0, 0, 5, 0), 0, 0));
        mainPanel.add(slider, new GridBagConstraints(1, 5, 2, 1, 1.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 0, 0), 0, 0));
        mainPanel.add(countTextField, new GridBagConstraints(1, 6, 1, 1, 0.7, 0.7, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
        mainPanel.add(countLabel1, new GridBagConstraints(2, 6, 1, 1, 0.0, 0.0, GridBagConstraints.LINE_START, GridBagConstraints.NONE, new Insets(0, 5, 0, 0), 0, 0));
        mainPanel.add(okayButton, new GridBagConstraints(1, 7, 1, 1, 0.0, 0.0, GridBagConstraints.LINE_END, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        mainPanel.add(cancelButton, new GridBagConstraints(2, 7, 1, 1, 0.0, 0.0, GridBagConstraints.LINE_START, GridBagConstraints.NONE, new Insets(0, 5, 0, 0), 0, 0));

        this.add(mainPanel, BorderLayout.CENTER);
    }

    private void addAllButton_actionPerformed(ActionEvent e) {
        list1.setListData(new Object[0]);
        list2.setListData(myTraits);
        okayButton.setEnabled(true);
    }

    private void addOneButton_actionPerformed(ActionEvent e) {
        Object[] selected = list1.getSelectedValues();
        ListModel listModel = list2.getModel();
        String[] toList = new String[listModel.getSize() + selected.length];
        String[] fromList = new String[myTraits.length - toList.length];
        int j = 0;
        int k = 0;
        int l = 0;
        int m = 0;
        for (int i = 0; i < myTraits.length; i++) {
            if (j < selected.length && myTraits[i].equals(selected[j])) {
                toList[l] = myTraits[i];
                j++;
                l++;
            } else if (k < listModel.getSize() && myTraits[i].equals(listModel.getElementAt(k))) {
                toList[l] = myTraits[i];
                k++;
                l++;
            } else {
                fromList[m] = myTraits[i];
                m++;
            }
        }
        list1.setListData(fromList);
        list2.setListData(toList);
        if (toList.length > 0) {
            okayButton.setEnabled(true);
        }
    }

    private void removeOneButton_actionPerformed(ActionEvent e) {
        Object[] selected = list2.getSelectedValues();
        ListModel listModel = list1.getModel();
        String[] toList = new String[listModel.getSize() + selected.length];
        String[] fromList = new String[myTraits.length - toList.length];
        int j = 0;
        int k = 0;
        int l = 0;
        int m = 0;
        for (int i = 0; i < myTraits.length; i++) {
            if (j < selected.length && myTraits[i].equals(selected[j])) {
                toList[l] = myTraits[i];
                j++;
                l++;
            } else if (k < listModel.getSize() && myTraits[i].equals(listModel.getElementAt(k))) {
                toList[l] = myTraits[i];
                k++;
                l++;
            } else {
                fromList[m] = myTraits[i];
                m++;
            }
        }
        list2.setListData(fromList);
        list1.setListData(toList);
        if (fromList.length < 1) {
            okayButton.setEnabled(false);
        }
    }

    private void removeAllButton_actionPerformed(ActionEvent e) {
        list2.setListData(new Object[0]);
        list1.setListData(myTraits);
        okayButton.setEnabled(false);
    }

    private void okayButton_actionPerformed(ActionEvent e) {
        isCanceled = false;
        setVisible(false);
    }

    private void cancelButton_actionPerformed(ActionEvent e) {
        isCanceled = true;
        setVisible(false);
    }

    private void slider_actionPerformed(ChangeEvent cd) {
        countTextField.setText("" + slider.getValue());
    }

    private void countTextField_keyTyped(KeyEvent e) {
        try {
            if (!countTextField.getText().equals("")) {
                int value = Integer.valueOf(countTextField.getText());
                if (value >= slider.getMinimum() && value <= slider.getMaximum()) {
                    slider.setValue(value);
                } else if (value <= slider.getMinimum()) {
                    slider.setValue(slider.getMinimum());
                    countTextField.setText("" + slider.getMinimum());
                } else if (value >= slider.getMaximum()) {
                    slider.setValue(slider.getMaximum());
                    countTextField.setText("" + slider.getMaximum());
                }
            }
        } catch (NumberFormatException nfe) {
            countTextField.setText("" + slider.getValue());
        }
    }

    public boolean isCanceled() {
        return isCanceled;
    }

    public int getSliderValue() {
        return slider.getValue();
    }

    public Dimension getMinimumSize() {
        return new Dimension(600, 600);
    }

    public String[] getSelectedTraits() {
        ListModel list = list2.getModel();
        String[] traits = new String[list.getSize()];
        for (int i = 0; i < traits.length; i++) {
            traits[i] = list.getElementAt(i).toString();
        }
        return traits;
    }

    public int[] getTraitIndices() {
        ListModel list = list2.getModel();
        int[] indices = new int[list.getSize()];
        int j = 0;
        for (int i = 0; i < myTraits.length; i++) {
            if (j < list.getSize() && myTraits[i].equals(list.getElementAt(j))) {
                indices[j] = i;
                j++;
            }
        }
        return indices;
    }
}