/*
 * SynonymizerPlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.baseplugins;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.Phenotype;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.IdentifierSynonymizer;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;

import javax.swing.*;
import javax.swing.table.TableModel;
import javax.swing.table.DefaultTableModel;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.net.URL;
import java.text.Collator;
import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;

/**
 *
 * @author Ed Buckler
 */
public class SynonymizerPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(SynonymizerPlugin.class);

    /**
     * Creates a new instance of SynonymizerPlugin
     */
    public SynonymizerPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {

        try {

            List<Datum> data = new ArrayList<Datum>();
            for (int i = 0, n = input.getSize(); i < n; i++) {
                Datum current = input.getData(i);
                Object currentData = current.getData();
                if (currentData instanceof Alignment) {
                    IdGroup idGroup = ((Alignment) currentData).getIdGroup();
                    Datum idGroupDatum = new Datum(current.getName(), idGroup, current.getComment());
                    data.add(idGroupDatum);
                } else if (currentData instanceof Phenotype) {
                    IdGroup idGroup = ((Phenotype) currentData).getTaxa();
                    Datum idGroupDatum = new Datum(current.getName(), idGroup, current.getComment());
                    data.add(idGroupDatum);
                } else {
                    data.add(current);
                }
            }
            DataSet newInput = new DataSet(data, this);

            int alignCnt = newInput.getDataOfType(IdGroup.class).size();
            int synCnt = newInput.getDataOfType(IdentifierSynonymizer.class).size();
            if ((synCnt == 0) && (alignCnt > 1)) {  //create a new synonymizer
                Datum td = createSynonymizer(newInput);
                DataSet output = new DataSet(td, this);
                fireDataSetReturned(new PluginEvent(output, SynonymizerPlugin.class));
                return output;
            } else if ((synCnt == 1) && (alignCnt > 0)) {   //apply synonymizer to alignments
                applySynonymsToIdGroups(newInput);
            } else if ((synCnt == 1) && (alignCnt == 0)) {
                if (isInteractive()) {
                    Datum inputDatum = newInput.getDataOfType(IdentifierSynonymizer.class).get(0);
                    IdentifierSynonymizer is = (IdentifierSynonymizer) inputDatum.getData();
                    SynonymizerDialog theSD = new SynonymizerDialog(is, getParentFrame());
                    theSD.setLocationRelativeTo(getParentFrame());
                    theSD.setVisible(true);
                }
            } else {
                String msg = "To create a synonym list:\n Please first select the reference taxa names and then the synonym taxa names (use Ctrl key)\n"
                        + "To apply a synonym list to a dataset:\n Select a synonym list and then the taxa names to be changed (use Ctrl key)";
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), msg);
                } else {
                    myLogger.error(msg);
                }
            }

            return null;

        } finally {
            fireProgress(100);
        }
    }

    private Datum createSynonymizer(DataSet input) {
        Datum td = null;
        StringBuilder synonymSets = new StringBuilder();
        for (int i = 1; i < input.getSize(); i++) {
            synonymSets.append(input.getData(i).getName());
            synonymSets.append("\n");
        }
        boolean performFunction = true;
        String msg = "You have selected to apply synonym list " + input.getData(0).getName() + " to the following dataset:\n"
                + synonymSets.toString();
        if (isInteractive()) {
            int response = JOptionPane.showOptionDialog(getParentFrame(), msg, "Verify Selection",
                    JOptionPane.OK_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null, null, null);
            if (response == JOptionPane.CANCEL_OPTION) {
                performFunction = false;
            }
        } else {
            myLogger.info(msg);
        }
        if (performFunction) {
            List<Datum> idList = input.getDataOfType(IdGroup.class);
            IdGroup[] aa = new IdGroup[idList.size() - 1];
            for (int i = 1; i < idList.size(); i++) {
                aa[i - 1] = (IdGroup) idList.get(i).getData();
            }
            IdentifierSynonymizer ts = new IdentifierSynonymizer((IdGroup) idList.get(0).getData(), aa);
            StringWriter sw = new StringWriter();
            ts.report(new PrintWriter(sw));
            td = new Datum(input.getData(0).getName() + " Synonyms", ts, "Taxa synonyms\n" + sw.toString());
        }
        return td;
    }

    private void applySynonymsToIdGroups(DataSet input) {
        StringBuilder synonymSets = new StringBuilder();
        for (int i = 1; i < input.getSize(); i++) {
            synonymSets.append(input.getData(i).getName());
            synonymSets.append("\n");
        }
        boolean performFunction = true;
        String msg = "You have selected " + input.getData(0).getName() + " as the reference name dataset.\n"
                + "The synonyms will be extracted from the following: \n" + synonymSets.toString();
        if (isInteractive()) {
            int response = JOptionPane.showOptionDialog(getParentFrame(), msg, "Verify Selection",
                    JOptionPane.OK_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null, null, null);
            if (response == JOptionPane.CANCEL_OPTION) {
                performFunction = false;
            }
        } else {
            myLogger.info(msg);
        }
        if (performFunction) {
            IdentifierSynonymizer is = (IdentifierSynonymizer) input.getDataOfType(IdentifierSynonymizer.class).get(0).getData();
            List<Datum> idList = input.getDataOfType(IdGroup.class);
            IdGroup[] aa = new IdGroup[idList.size()];
            for (int i = 0; i < idList.size(); i++) {
                aa[i] = (IdGroup) idList.get(i).getData();
            }
            is.changeAlignmentIdentifiers(aa);
        }
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        URL imageURL = SynonymizerPlugin.class.getResource("images/Synonymizer.gif");
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
        return "Synonymizer";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Unify Taxa Names";
    }
}

/**
 */
class SynonymizerDialog extends JDialog {

    private JPanel jPanel1 = new JPanel();
    private JTextField ThresholdTextField = new JTextField();
    private JButton setThresholdButton = new JButton();
    private JButton CancelButton = new JButton();
    private Frame theFrame;
    private double threshold = 1.0;
    boolean isCanceled;
    JButton okButton = new JButton();
    JList matchList = new JList();
    JScrollPane newRealNameScrollPane1 = new JScrollPane();
    JTable synTable = new JTable();
    JButton selectSynButton = new JButton();
    JButton setNoSynButton = new JButton();
    JLabel jLabel1 = new JLabel();
    JPanel jPanel2 = new JPanel();
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    JScrollPane theATP;
    JTable theNameTable;
    IdentifierSynonymizer theTS;
    JCheckBox cbxSortAlphabetically = new JCheckBox();

    public SynonymizerDialog(IdentifierSynonymizer ts, Frame theFrame) {
        super((Frame) theFrame, true);
        this.theFrame = theFrame;
        this.theTS = ts;
        try {
            theNameTable = new JTable(new DefaultTableModel(theTS.getTableData(), theTS.getTableColumnNames()));
            theNameTable.setAutoCreateRowSorter(true);
            theNameTable.setCellEditor(null);
            theNameTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
            matchList.setAutoscrolls(true);
            theATP = new JScrollPane(theNameTable);
            jbInit();
            pack();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public SynonymizerDialog() {
        try {
            jbInit();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void jbInit() throws Exception {
        jPanel1.setLayout(gridBagLayout1);
        setThresholdButton.setText("Apply threshold");
        setThresholdButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                setThresholdButton_actionPerformed(e);
            }
        });
        CancelButton.setText("Cancel");
        CancelButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                CancelButton_actionPerformed(e);
            }
        });
        okButton.setText("OK");
        okButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                okButton_actionPerformed(e);
            }
        });
        ThresholdTextField.setPreferredSize(new Dimension(30, 30));
        selectSynButton.setFont(new java.awt.Font("Dialog", Font.BOLD, 14));
        selectSynButton.setToolTipText("Set synonym to selected taxon");
        selectSynButton.setText("<");
        selectSynButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                selectSynButton_actionPerformed(e);
            }
        });
        // setNoSynButton.setFont(new java.awt.Font("Dialog", Font.BOLD, 14));
        setNoSynButton.setToolTipText("Set selected taxon to no synonym");
        setNoSynButton.setText("No Synonym");
        setNoSynButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                setNoSynButton_actionPerformed(e);
            }
        });
        theNameTable.addMouseListener(new MouseAdapter() {
            public void mouseClicked(MouseEvent e) {
                synTable_mouseClicked(e);
            }
        });
        jLabel1.setFont(new java.awt.Font("Dialog", Font.BOLD, 14));
        jLabel1.setText("Synonymizer");
        newRealNameScrollPane1.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
        cbxSortAlphabetically.setActionCommand("jCheckBox1");
        cbxSortAlphabetically.setHorizontalAlignment(SwingConstants.CENTER);
        cbxSortAlphabetically.setHorizontalTextPosition(SwingConstants.TRAILING);
        cbxSortAlphabetically.setSelectedIcon(null);
        cbxSortAlphabetically.setText("Sort Alphabetically");
        cbxSortAlphabetically.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                cbxSortAlphabetically_actionPerformed(e);
            }
        });

        this.getContentPane().add(jPanel1, BorderLayout.CENTER);

        jPanel2.add(theATP);
        newRealNameScrollPane1.getViewport().add(matchList);
        jPanel1.add(ThresholdTextField, new GridBagConstraints(1, 2, 1, 2, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 50, 0));
        jPanel1.add(setThresholdButton, new GridBagConstraints(2, 3, 3, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 37, 7));
        jPanel1.add(newRealNameScrollPane1, new GridBagConstraints(2, 1, 2, 1, 1.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 100, 300));
        jPanel1.add(selectSynButton, new GridBagConstraints(1, 1, 1, 1, 0.5, 0.5, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 10, 4));
        jPanel1.add(setNoSynButton, new GridBagConstraints(0, 2, 1, 2, 0.5, 0.5, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 10, 4));
        jPanel1.add(CancelButton, new GridBagConstraints(2, 4, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(19, 0, 10, 25), 15, 7));
        jPanel1.add(jLabel1, new GridBagConstraints(0, 0, 3, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(8, 20, 0, 0), 184, 5));
        jPanel1.add(theATP, new GridBagConstraints(0, 1, 1, 1, 2.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 200, 300));
        jPanel1.add(okButton, new GridBagConstraints(0, 4, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(19, 69, 10, 49), 33, 10));
        jPanel1.add(cbxSortAlphabetically, new GridBagConstraints(3, 2, 1, 1, 0.0, 0.0, GridBagConstraints.SOUTHWEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        this.setSize(800, 600);
        this.setTitle(" Threshold for synonymizer ");
    }

    public IdentifierSynonymizer getIdentifierSynonymizer() {
        return theTS;
    }

    void deleteByThreshold(double threshold) {
        TableModel dm = theNameTable.getModel();
        String synName, realName;
        double score;
        for (int i = 0; i < dm.getRowCount(); i++) {
            synName = (String) dm.getValueAt(i, 0);
            realName = (String) dm.getValueAt(i, 1);
            score = IdentifierSynonymizer.scoreMatch(synName, realName, true, false, false);
            if (score < threshold) {
                dm.setValueAt("", i, 1);
                dm.setValueAt("-1", i, 2);
            }
        }
    }

    void setThresholdButton_actionPerformed(ActionEvent e) {
        threshold = getMatchThreshold();
        deleteByThreshold(threshold);
    }

    void CancelButton_actionPerformed(ActionEvent e) {
        isCanceled = true;
        this.setVisible(false);
    }

    double getMatchThreshold() {
        double th = threshold;
        try {
            th = Double.parseDouble(ThresholdTextField.getText().trim());
            if ((th < 0) || (th > 1)) {
                throw new NumberFormatException();
            }
        } catch (NumberFormatException nfe) {
            JOptionPane.showMessageDialog(theFrame, "Please enter an double between 0 and 1.");
        }
        return th;
    }

    public boolean isCanceled() {
        return isCanceled;
    }

    void selectSynButton_actionPerformed(ActionEvent e) {
        try {
            String newRealName = (String) matchList.getSelectedValue();
            int theRow = theNameTable.getSelectedRow();
            theNameTable.getModel().setValueAt(newRealName, theRow, 1);
            theNameTable.getModel().setValueAt("" + theTS.getPreferredIndex(newRealName), theRow, 2);
        } catch (Exception ex) {
            System.out.println("Make sure both a row and a new name are selected");
        }
    }

    void setNoSynButton_actionPerformed(ActionEvent e) {
        try {
            int theRow = theNameTable.getSelectedRow();
            theNameTable.getModel().setValueAt("", theRow, 1);
            theNameTable.getModel().setValueAt("-1", theRow, 2);
        } catch (Exception ex) {
            System.out.println("Make sure a row is selected");
        }
    }

    void synTable_mouseClicked(MouseEvent e) {
        if (cbxSortAlphabetically.isSelected()) {
            sortListAlphabetically();
        } else {
            sortListByMatchScore();
        }
    }

    void okButton_actionPerformed(ActionEvent e) {
        TableModel dm = theNameTable.getModel();
        String synName, realName;
        int newID;
        for (int i = 0; i < dm.getRowCount(); i++) {
            synName = (String) dm.getValueAt(i, 0);
            newID = Integer.parseInt((String) dm.getValueAt(i, 2));
            if (theTS.getPreferredIndex(synName) != newID) {
                System.out.println("synName=" + synName + "  " + theTS.getPreferredName(synName));
                theTS.setRealID(synName, newID);
                System.out.println("synName=" + synName + "  " + theTS.getPreferredName(synName));
            }
        }
        isCanceled = false;
        this.setVisible(false);
    }

    void sortListByMatchScore() {
        Object theSynonym = theNameTable.getModel().getValueAt(theNameTable.getSelectedRow(), 0);
        ArrayList findOrderedMatches = theTS.findOrderedMatches((String) theSynonym, 4);
        DefaultListModel dlm = new DefaultListModel();
        Object[] a = findOrderedMatches.toArray();
        for (int i = 0; i < a.length; i++) {
            dlm.insertElementAt(a[i], i);
        }
        matchList.setModel(dlm);
    }

    private void cbxSortAlphabetically_actionPerformed(ActionEvent e) {

        if (cbxSortAlphabetically.isSelected()) {
            sortListAlphabetically();
        } else {
            sortListByMatchScore();
        }
    }

    private void sortListAlphabetically() {
        DefaultListModel listModel = (DefaultListModel) matchList.getModel();

        int itemCount = listModel.getSize();
        String[] a = new String[itemCount];

        listModel.copyInto(a);

        sortArray(Collator.getInstance(), a);

        for (int i = 0; i < itemCount; i++) {
            listModel.setElementAt(a[i], i);
        }
    }

    private void sortArray(Collator collator, String[] strArray) {
        String tmp;
        if (strArray.length == 1) {
            return;
        }
        for (int i = 0; i < strArray.length; i++) {
            for (int j = i + 1; j < strArray.length; j++) {
                if (collator.compare(strArray[i], strArray[j]) > 0) {
                    tmp = strArray[i];
                    strArray[i] = strArray[j];
                    strArray[j] = tmp;
                }
            }
        }
    }
}
