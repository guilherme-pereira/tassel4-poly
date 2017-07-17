/*
 * CreateTreePlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.baseplugins;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.pal.distance.IBSDistanceMatrix;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;

import java.net.URL;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 *
 * @author Ed Buckler
 */
public class CreateTreePlugin extends AbstractPlugin {
    //  Settings theSettings;

    private boolean isNeighborJoining = true;
    private boolean isReturnDistanceMatrix = true;
    //private boolean isReturnGroupTable = false;

    /** Creates a new instance of CreateTreePlugin */
    public CreateTreePlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {

        try {

            List<Datum> alignInList = input.getDataOfType(Alignment.class);
            if (alignInList.size() < 1) {
                String message = "Invalid selection.  Please select sequence or marker alignment.";
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), message);
                } else {
                    System.out.println(message);
                }
                return null;
            }

            if (isInteractive()) {
                CreateTreeDialog myDialog = new CreateTreeDialog();
                myDialog.setLocationRelativeTo(getParentFrame());
                myDialog.setVisible(true);
                if (myDialog.isCancel()) {
                    return null;
                }
                isNeighborJoining = myDialog.isNeighborJoiningTree();
                isReturnDistanceMatrix = myDialog.isSaveMatrix();
                //isReturnGroupTable = myDialog.isSaveGroups();
                myDialog.dispose();
            }

            List result = new ArrayList();
            Iterator<Datum> itr = alignInList.iterator();
            while (itr.hasNext()) {
                Datum current = itr.next();
                DataSet tds = processDatum(current, isNeighborJoining, isReturnDistanceMatrix);
                result.add(tds);
                if (tds != null) {
                    fireDataSetReturned(new PluginEvent(tds, CreateTreePlugin.class));
                }
            }

            return DataSet.getDataSet(result, this);

        } finally {
            fireProgress(100);
        }

    }

    public DataSet processDatum(Datum input, boolean isNJ, boolean isSaveMatrix) {
        Alignment aa = (Alignment) input.getData();
        // SitePattern sp = new SitePattern(aa);
        IBSDistanceMatrix adm = new IBSDistanceMatrix(aa, this);
        // adm.recompute(sp);
        net.maizegenetics.pal.tree.Tree theTree;
        List<Datum> results = new ArrayList<Datum>();
        if (isNJ) {
            theTree = new net.maizegenetics.pal.tree.NeighborJoiningTree(adm);
            results.add(new Datum("Tree:" + input.getName(), theTree, "NJ Tree"));
        } else {
            theTree = new net.maizegenetics.pal.tree.UPGMATree(adm);
            results.add(new Datum("Tree:" + input.getName(), theTree, "UPGMA Tree"));
        }
        if (isSaveMatrix) {
            results.add(new Datum("Matrix:" + input.getName(), adm, "Distance Matrix"));
        }

        //TreeCut tc = new TreeCut(theTree);
        /*
        // display grouping results;
        if (isReturnGroupTable) {
        }
         */

        DataSet tds = new DataSet(results, this);
        return tds;
    }

    public boolean isNeighborJoining() {
        return isNeighborJoining;
    }

    public void setNeighborJoining(boolean neighborJoining) {
        isNeighborJoining = neighborJoining;
    }

    public boolean isReturnDistanceMatrix() {
        return isReturnDistanceMatrix;
    }

    public void setReturnDistanceMatrix(boolean returnDistanceMatrix) {
        isReturnDistanceMatrix = returnDistanceMatrix;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        URL imageURL = CreateTreePlugin.class.getResource("images/Tree.gif");
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
        return "Cladogram";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Create a cladogram";
    }
}

/**
 * Title:        TASSEL
 * Description:  A java program to deal with diversity
 * Copyright:    Copyright (c) 2000
 * Company:      USDA-ARS/NCSU
 * @author Ed Buckler
 * @version 1.0
 */
class CreateTreeDialog extends JDialog {

    boolean saveMatrix = true;
    //boolean saveGroups = false;
    boolean runAnalysis = false;
    JPanel panel1 = new JPanel();
    JCheckBox saveMatrixCheckBox = new JCheckBox();
    //JCheckBox saveGroupsCheckBox = new JCheckBox();
    JButton runButton = new JButton();
    JButton closeButton = new JButton();
    //  String[] subModels = {"None","Max Likelihood"};
    //  JComboBox subModelComboBox = new JComboBox(subModels);
    //  JLabel jLabel2 = new JLabel();
    String[] clustringModels = {"Neighbor Joining", "UPGMA"};
    JComboBox clusteringComboBox = new JComboBox(clustringModels);
    JLabel jLabel3 = new JLabel();
    GridBagLayout gridBagLayout1 = new GridBagLayout();

    public CreateTreeDialog() {
        super((Frame) null, "Create Tree", true);
        try {
            jbInit();
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    void jbInit() throws Exception {
        panel1.setLayout(gridBagLayout1);
        saveMatrixCheckBox.setSelected(saveMatrix);
        saveMatrixCheckBox.setText("Save distance matrix");
        //saveGroupsCheckBox.setSelected(saveGroups);
        //saveGroupsCheckBox.setText("Save taxa groups");
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
        //    jLabel2.setText("Substitution Model");
        jLabel3.setText("Clustering Method");
        getContentPane().add(panel1);
        panel1.add(runButton, new GridBagConstraints(0, 6, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(76, 55, 26, 0), 32, 5));
        panel1.add(closeButton, new GridBagConstraints(1, 6, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(77, 20, 26, 156), 23, 5));
        //    panel1.add(subModelComboBox,  new GridBagConstraints(0, 1, 2, 1, 1.0, 0.0
        //            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 105, 0, 152), 50, -1));
        //    panel1.add(jLabel2,  new GridBagConstraints(0, 0, 2, 1, 0.0, 0.0
        //            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(36, 107, 0, 171), 35, 9));
        panel1.add(saveMatrixCheckBox, new GridBagConstraints(1, 2, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(6, 67, 0, 28), 38, 0));
        //panel1.add(saveGroupsCheckBox, new GridBagConstraints(1, 3, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(6, 67, 0, 28), 38, 0));
        panel1.add(clusteringComboBox, new GridBagConstraints(0, 5, 2, 1, 1.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 111, 0, 146), 41, -1));
        panel1.add(jLabel3, new GridBagConstraints(0, 4, 2, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(6, 113, 0, 165), 36, 9));
    }

    /** Returns whether the run button was chosen*/
    public boolean isRunAnalysis() {
        return runAnalysis;
    }

    public boolean isCancel() {
        return !runAnalysis;
    }

    /** Returns whether the matrix should be saved */
    public boolean isSaveMatrix() {
        return saveMatrix;
    }

    //public boolean isSaveGroups() {
    //    return saveGroups;
    //}

    /** Returns whether the matrix should be saved */
    public boolean isNeighborJoiningTree() {
        String s = (String) clusteringComboBox.getSelectedItem();
        if (s.equals(clustringModels[0])) {
            return true;
        }
        return false;
    }

    void runButton_actionPerformed(ActionEvent e) {
        saveMatrix = saveMatrixCheckBox.isSelected();
        //saveGroups = saveGroupsCheckBox.isSelected();
        runAnalysis = true;
        setVisible(false);
    }

    void closeButton_actionPerformed(ActionEvent e) {
        runAnalysis = false;
        setVisible(false);
    }
}
