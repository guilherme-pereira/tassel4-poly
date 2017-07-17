/*
 * TreeDisplayPlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.baseplugins;

import net.maizegenetics.pal.tree.Tree;
import net.maizegenetics.pal.alignment.Phenotype;
import net.maizegenetics.pal.alignment.Trait;
import net.maizegenetics.pal.gui.QuantitativeLegendComponent;
import net.maizegenetics.pal.gui.TreeComponent;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.net.URL;
import java.util.List;

/**
 *
 * @author Ed Buckler
 */
public class TreeDisplayPlugin extends AbstractDisplayPlugin {

    int treeMode = TreeComponent.NORMAL_BW;

    /** Creates a new instance of TreeDisplayPlugin */
    public TreeDisplayPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {

        try {
            List<Datum> treeInList = input.getDataOfType(Tree.class);
            List<Datum> caInList = input.getDataOfType(Phenotype.class);
            if (treeInList.size() != 1) {
                String message = "Invalid selection.  Please select a 1 tree and upto 1 numerical.";
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), message);
                } else {
                    System.out.println(message);
                }
                return null;
            }
            Phenotype ca = null;
            if (caInList.size() > 0) {
                ca = (Phenotype) caInList.get(0).getData();
            }
            Tree theTree = (Tree) treeInList.get(0).getData();
            if (isInteractive()) {
                TreePluginDialog myDialog = new TreePluginDialog(this, theTree, ca);
                myDialog.setLocationRelativeTo(getParentFrame());
                myDialog.setVisible(true);
            } else if (getSaveFile() != null) {
                TreeComponent tc = new TreeComponent(theTree, false);
                tc.setSize(getImageWidth(), getImageHeight());
                tc.setMode(getTreeMode());
                saveDataToFile(tc, getSaveFile());
            }

            return null;
        } finally {
            fireProgress(100);
        }

    }

    public int getTreeMode() {
        return treeMode;
    }

    public void setTreeMode(int treeMode) {
        this.treeMode = treeMode;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        URL imageURL = TreeDisplayPlugin.class.getResource("images/Tree.gif");
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
        return "Tree Plot";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Display Cladogram of tree";
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
class TreePluginDialog extends JDialog {

    JPanel panel1 = new JPanel();
    BorderLayout borderLayout1 = new BorderLayout();
    Tree theTree;
    Phenotype pheno;
    TreeDisplayPlugin theTreeDisplayPlugin = null;
    QuantitativeLegendComponent theQuantitativeLegendComponent;
    TreeComponent theQTreeComponent;
    JPanel jPanel2 = new JPanel();
    JSpinner fontSpinner = new JSpinner(new SpinnerNumberModel(15, 1, 40, 1));
    JButton normalButton = new JButton();
    JPanel jPanel1 = new JPanel();
    JButton saveButton = new JButton();
    JButton circularButton = new JButton();
    FlowLayout flowLayout1 = new FlowLayout();
    BorderLayout borderLayout2 = new BorderLayout();
    JPanel jPanel3 = new JPanel();
    JButton quantitativeButton = new JButton();
    JComboBox traitComboBox = new JComboBox();
    JPanel printPanel = new JPanel(true);
    BorderLayout borderLayout3 = new BorderLayout();
//  JComboBox formatComboBox = new JComboBox();

    public TreePluginDialog(TreeDisplayPlugin tdp, Tree theTree, Phenotype pheno) {
        super(tdp.getParentFrame(), "Tree", false);
        //   if(frame instanceof TASSELMainFrame) {theTASSELMainFrame=(TASSELMainFrame)frame;}
        theTreeDisplayPlugin = tdp;
        this.theTree = theTree;
        this.pheno = pheno;
        try {
            theQTreeComponent = new TreeComponent(theTree, "");
            theQTreeComponent.setSize(500, 500);
            theQTreeComponent.setMode(TreeComponent.NORMAL_BW);
            jbInit();
            if (pheno == null) {
                jPanel2.remove(jPanel3);
                printPanel.remove(theQuantitativeLegendComponent);
            } else {
                String[] s = new String[pheno.getNumberOfTraits()];
                for (int i = 0; i < pheno.getNumberOfTraits(); i++) {
                    traitComboBox.addItem(pheno.getTrait(i).getName() + "." + pheno.getTrait(i).getProperty(Trait.FACTOR_ENV));
                }
            }
            pack();
            this.setSize(500, 500);
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    void jbInit() throws Exception {
        panel1.setLayout(borderLayout1);
        normalButton.setText("Normal");
        normalButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                normalButton_actionPerformed(e);
            }
        });
        saveButton.setText("Save");
        saveButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                saveButton_actionPerformed(e);
            }
        });
        circularButton.setText("Circular");
        circularButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                circularButton_actionPerformed(e);
            }
        });
        jPanel1.setLayout(flowLayout1);
        jPanel2.setLayout(borderLayout2);
        quantitativeButton.setToolTipText("Display quantitative traits");
        quantitativeButton.setText("Quantiative");
        quantitativeButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                quantitativeButton_actionPerformed(e);
            }
        });
//    traitComboBox.setItems(new String[] {"Item 1", "Item 2"});
        traitComboBox.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                traitComboBox_actionPerformed(e);
            }
        });
        fontSpinner.addChangeListener(new javax.swing.event.ChangeListener() {

            public void stateChanged(ChangeEvent e) {
                fontSpinner_stateChanged(e);
            }
        });
        printPanel.setLayout(borderLayout3);
//    formatComboBox.setText("Save SVG");
//    String[] s=AbstractDisplayPlugin.getPossibleGraphicOutFormats();
//    formatComboBox=new JComboBox(s);
        panel1.setBorder(BorderFactory.createLineBorder(Color.GREEN));
        panel1.add(jPanel2, BorderLayout.SOUTH);
        jPanel2.add(jPanel1, BorderLayout.CENTER);
        jPanel1.add(normalButton, null);
        jPanel1.add(circularButton, null);
        jPanel1.add(new JLabel(" Font Size:"), null);
        jPanel1.add(fontSpinner, null);
        //   jPanel1.add(printButton, null);
        jPanel1.add(saveButton, null);
//    jPanel1.add(formatComboBox, null);
        jPanel2.add(jPanel3, BorderLayout.SOUTH);
        jPanel3.add(quantitativeButton, null);
        jPanel3.add(traitComboBox, null);
        panel1.add(printPanel, BorderLayout.CENTER);
        printPanel.add(theQTreeComponent, BorderLayout.CENTER);
        theQuantitativeLegendComponent = new QuantitativeLegendComponent(0, 0);
        printPanel.add(theQuantitativeLegendComponent, BorderLayout.EAST);
        getContentPane().add(panel1);
    }

    void normalButton_actionPerformed(ActionEvent e) {
        theQTreeComponent.setMode(TreeComponent.NORMAL_BW);
        theQTreeComponent.repaint();
    }

    void circularButton_actionPerformed(ActionEvent e) {
        theQTreeComponent.setMode(TreeComponent.CIRCULAR_BW);
        theQTreeComponent.repaint();
    }

    void traitComboBox_actionPerformed(ActionEvent e) {
        int trait = traitComboBox.getSelectedIndex();
        theQTreeComponent.setTrait(pheno, trait);
        theQTreeComponent.repaint();
        theQuantitativeLegendComponent.setTrait(pheno, trait);
        theQuantitativeLegendComponent.repaint();
    }

    void quantitativeButton_actionPerformed(ActionEvent e) {
        theQTreeComponent.setMode(TreeComponent.QUANTITATIVE);
        traitComboBox_actionPerformed(null);
    }

    void fontSpinner_stateChanged(ChangeEvent e) {
        Font f = theQTreeComponent.getFont();
        int size = (new Integer(fontSpinner.getValue().toString())).intValue();
        theQTreeComponent.setLabelFontSize(size);
        theQTreeComponent.repaint();
    }

    public JPanel getPrintPanel() {
        return printPanel;
    }

    void saveButton_actionPerformed(ActionEvent e) {
//    String s=formatComboBox.getSelectedItem().toString();
        theTreeDisplayPlugin.saveDataToFile(printPanel);
    }
}
