/*
 * GenotypeImputationPlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.baseplugins;

import javax.swing.*;
import java.awt.*;
import java.awt.event.FocusEvent;

import java.net.URL;

import java.util.ArrayList;
import java.util.List;

import net.maizegenetics.baseplugins.GenotypeImputationPlugin.ImpMethod;
import net.maizegenetics.gui.DialogUtils;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AlignmentMask;
import net.maizegenetics.pal.alignment.AlignmentMask.MaskType;
import net.maizegenetics.pal.alignment.AlignmentMaskBoolean;
import net.maizegenetics.pal.popgen.BasicImputation;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.util.ExceptionUtils;
import net.maizegenetics.util.Utils;

import org.apache.log4j.Logger;

/**
 * Plugin supports the imputation of genotypic data using a variety of methods.
 * It initially will only support imputation by length
 *
 * @author Ed Buckler
 */
public class GenotypeImputationPlugin extends AbstractDisplayPlugin {
    
    private static final Logger myLogger = Logger.getLogger(GenotypeImputationPlugin.class);
    
    public enum ImpMethod {
        
        Length, MajorAllele, SimilarWindow, IBDProb
    };
    private int minLength = 30;  //below the length threshold impute with the common allele
    private int maxMismatch = 1;
    private double minProb = 0.001;
    private ImpMethod currMethod = ImpMethod.Length;

    /**
     * Creates a new instance of GenotypeImputationPlugin
     */
    public GenotypeImputationPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }
    
    public DataSet performFunction(DataSet input) {
        
        StringBuilder builder = new StringBuilder();
        builder.append("\nThere are many possible algorithms for Imputation.\n");
        builder.append("After consideration, we may make one or more available\n");
        builder.append("via TASSEL.  But currently, please impute your data by\n");
        builder.append("another software package.");
        if (isInteractive()) {
            JOptionPane.showMessageDialog(getParentFrame(), builder.toString(), "Info", JOptionPane.INFORMATION_MESSAGE);
        } else {
            myLogger.info("GenotypeImputationPlugin: performFunction: " + builder.toString());
        }
        return null;
    }
    
    public DataSet performFunctionOLD(DataSet input) {
        
        List<Datum> alignInList = input.getDataOfType(Alignment.class);
        List result = new ArrayList();
        
        try {
            
            if (alignInList.size() == 0) {
                throw new IllegalArgumentException("Invalid selection.  Please select 1 or more alignments.");
            }
            if (isInteractive()) {
                GenotypeImputationPluginDialog myDialog = new GenotypeImputationPluginDialog(null);
                myDialog.setLocationRelativeTo(getParentFrame());
                myDialog.setVisible(true);
                if (myDialog.isRunAnalysis()) {
                    minLength = myDialog.getMinLength();
                    maxMismatch = myDialog.getMaxMismatch();
                    minProb = myDialog.getMinProb();
                    currMethod = myDialog.getChosenImpMethod();
                } else {
                    return null;
                }
            }
            for (Datum d : alignInList) {
                DataSet ds = processDatum(d);
                if (ds != null) {
                    result.addAll(ds.getDataSet());
                    fireDataSetReturned(new PluginEvent(ds, GenotypeImputationPlugin.class));
                }
            }
            
        } catch (Exception e) {
            e.printStackTrace();
            StringBuilder builder = new StringBuilder();
            builder.append(Utils.shortenStrLineLen(ExceptionUtils.getExceptionCauses(e), 50));
            String str = builder.toString();
            if (isInteractive()) {
                DialogUtils.showError(str, getParentFrame());
            } else {
                myLogger.error(str);
            }
            
            return null;
        } finally {
            fireProgress(100);
        }
        
        DataSet resultSet = new DataSet(result, this);
        return resultSet;
        
    }
    
    public DataSet processDatum(Datum inDatum) {
        Alignment align = (Alignment) inDatum.getData();
        Alignment impP1A = BasicImputation.imputeBySite(align, minLength, maxMismatch);
        
        String theName, theComment;
        theName = inDatum.getName() + "_Imp";
        theComment = "Imputed Data Of\n" + inDatum.getComment();
        Datum outDatum = new Datum(theName, impP1A, theComment);
        
        AlignmentMask mask = AlignmentMaskBoolean.getInstanceCompareAlignments(impP1A, align, "Imputed", MaskType.imputed);
        Datum maskDatum = new Datum(mask.toString(), mask, null);
        
        DataSet result = new DataSet(new Datum[]{outDatum, maskDatum}, this);
        
        return result;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        URL imageURL = GenotypeImputationPlugin.class.getResource("images/ImputeSNP.gif");
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
        return "Impute SNPs";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Impute genotypic data";
    }
    
    public void setMinLength(int length) {
        minLength = length;
    }
    
    public void setMaxMisMatch(int max) {
        maxMismatch = max;
    }
    
    public void setMinProb(double num) {
        minProb = num;
    }
    
    public void setMethod(ImpMethod method) {
        currMethod = method;
    }
}

/**
 * Title: TASSEL Description: A java program to deal with diversity Copyright:
 * Copyright (c) 2000 Company: USDA-ARS/NCSU
 *
 * @author Ed Buckler
 * @version 1.0
 */
class GenotypeImputationPluginDialog extends JDialog {
    
    boolean runAnalysis = false;
    private GenotypeImputationPlugin.ImpMethod chosenImpMethod = GenotypeImputationPlugin.ImpMethod.Length;
    private int minLength = 31;
    private int maxMismatch = 1;
    private double minProb = 0.001;
    
    public GenotypeImputationPluginDialog(Frame f) {
        super((Frame) f, "Imputation", true);
        try {
            initComponents();
            this.lengthImputeRadioButton.setSelected(true);
            this.majorAlleleRadioButton.setEnabled(false);
            this.ibdProbRadioButton.setEnabled(false);
            this.similarWindowRadioButton.setEnabled(false);
            //pack();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
    
    public boolean isRunAnalysis() {
        return runAnalysis;
    }
    
    public int getMaxMismatch() {
        return maxMismatch;
    }
    
    public int getMinLength() {
        return minLength;
    }
    
    public double getMinProb() {
        return minProb;
    }
    
    public ImpMethod getChosenImpMethod() {
        return chosenImpMethod;
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">
    private void initComponents() {
        
        buttonGroup1 = new javax.swing.ButtonGroup();
        jPanel1 = new javax.swing.JPanel();
        jLabel1 = new javax.swing.JLabel("Imputation of Genotypic Data");
        lengthImputeRadioButton = new javax.swing.JRadioButton("Maximum shared length");
        lengthTextField = new javax.swing.JTextField("" + minLength);
        maxMismatchTextField = new javax.swing.JTextField("" + maxMismatch);
        jLabel2 = new javax.swing.JLabel("Min Length");
        jLabel4 = new javax.swing.JLabel("Maximum Mismatch");
        majorAlleleRadioButton = new javax.swing.JRadioButton("Major Allele");
        similarWindowRadioButton = new javax.swing.JRadioButton("Simlarity by window");
        ibdProbRadioButton = new javax.swing.JRadioButton("IBD Probabilty");
        jLabel3 = new javax.swing.JLabel("Min Probability");
        minProbTextField = new javax.swing.JTextField("" + minProb);
        runButton = new javax.swing.JButton("Run");
        cancelButton = new javax.swing.JButton("Close");
        buttonGroup1.add(lengthImputeRadioButton);
        buttonGroup1.add(majorAlleleRadioButton);
        buttonGroup1.add(similarWindowRadioButton);
        buttonGroup1.add(ibdProbRadioButton);
        jLabel1.setFont(new java.awt.Font("Tahoma", 0, 18));
        //  lengthImputeRadioButton.setText("Maximum shared length");
        lengthTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(FocusEvent e) {
                lengthTextFieldFocusLost(e);
            }
        });
        maxMismatchTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(FocusEvent e) {
                maxMismatchTextFieldFocusLost(e);
            }
        });
        minProbTextField.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(FocusEvent e) {
                minProbTextFieldFocusLost(e);
            }
        });
        runButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                runButtonActionPerformed(evt);
            }
        });
        cancelButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                cancelButtonActionPerformed(evt);
            }
        });
        javax.swing.GroupLayout jPanel1Layout = new javax.swing.GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
                jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING).addGroup(jPanel1Layout.createSequentialGroup().addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE).addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING).addComponent(lengthImputeRadioButton).addComponent(jLabel1).addGroup(jPanel1Layout.createSequentialGroup().addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING).addGroup(jPanel1Layout.createSequentialGroup().addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING).addGroup(jPanel1Layout.createSequentialGroup().addComponent(majorAlleleRadioButton).addGap(58, 58, 58)).addGroup(jPanel1Layout.createSequentialGroup().addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING).addComponent(jLabel2).addComponent(jLabel4)).addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED))).addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING).addComponent(maxMismatchTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 44, javax.swing.GroupLayout.PREFERRED_SIZE).addComponent(lengthTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 44, javax.swing.GroupLayout.PREFERRED_SIZE))).addGroup(javax.swing.GroupLayout.Alignment.LEADING, jPanel1Layout.createSequentialGroup().addComponent(similarWindowRadioButton).addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)).addComponent(ibdProbRadioButton, javax.swing.GroupLayout.Alignment.LEADING).addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING).addComponent(runButton, javax.swing.GroupLayout.PREFERRED_SIZE, 64, javax.swing.GroupLayout.PREFERRED_SIZE).addGroup(jPanel1Layout.createSequentialGroup().addComponent(jLabel3).addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED).addComponent(minProbTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 64, javax.swing.GroupLayout.PREFERRED_SIZE)))).addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED).addComponent(cancelButton))).addGap(156, 156, 156)));
        jPanel1Layout.setVerticalGroup(
                jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING).addGroup(jPanel1Layout.createSequentialGroup().addContainerGap().addComponent(jLabel1).addGap(29, 29, 29).addComponent(lengthImputeRadioButton).addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED).addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE).addComponent(lengthTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE).addComponent(jLabel2)).addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED).addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE).addComponent(maxMismatchTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE).addComponent(jLabel4)).addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED).addComponent(majorAlleleRadioButton).addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED).addComponent(similarWindowRadioButton).addGap(5, 5, 5).addComponent(ibdProbRadioButton).addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED).addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING).addComponent(minProbTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE).addComponent(jLabel3)).addGap(18, 18, 18).addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING).addComponent(runButton).addComponent(cancelButton)).addContainerGap(28, Short.MAX_VALUE)));
        
        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
                layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING).addGroup(layout.createSequentialGroup().addComponent(jPanel1, javax.swing.GroupLayout.PREFERRED_SIZE, 340, javax.swing.GroupLayout.PREFERRED_SIZE).addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)));
        layout.setVerticalGroup(
                layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING).addComponent(jPanel1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE));
        
        pack();
        
    }// </editor-fold>

    private void lengthTextFieldFocusLost(FocusEvent evt) {
        try {
            int v = Integer.parseInt(lengthTextField.getText().trim());
            if (v < 0) {
                v = minLength;
            }
            minLength = v;
        } catch (Exception ee) {
            System.err.println(ee);
        }
    }
    
    private void maxMismatchTextFieldFocusLost(FocusEvent evt) {
        try {
            int v = Integer.parseInt(maxMismatchTextField.getText().trim());
            if (v >= 0) {
                v = this.maxMismatch;
            }
            maxMismatch = v;
        } catch (Exception ee) {
            System.err.println(ee);
        }
    }
    
    private void minProbTextFieldFocusLost(FocusEvent evt) {
        try {
            double v = Double.parseDouble(minProbTextField.getText().trim());
            if ((v < 0) || (v > 1)) {
                v = this.minProb;
            }
            minProb = v;
        } catch (Exception ee) {
            System.err.println(ee);
        }
    }
    
    private void runButtonActionPerformed(java.awt.event.ActionEvent evt) {
        // TODO add your handling code here:
        if (lengthImputeRadioButton.isSelected()) {
            chosenImpMethod = ImpMethod.Length;
        } else if (ibdProbRadioButton.isSelected()) {
            chosenImpMethod = ImpMethod.IBDProb;
        } else if (majorAlleleRadioButton.isSelected()) {
            chosenImpMethod = ImpMethod.MajorAllele;
        } else if (similarWindowRadioButton.isSelected()) {
            chosenImpMethod = ImpMethod.SimilarWindow;
        }
        lengthTextFieldFocusLost(null);
        maxMismatchTextFieldFocusLost(null);
        minProbTextFieldFocusLost(null);
        runAnalysis = true;
        setVisible(false);
    }
    
    private void cancelButtonActionPerformed(java.awt.event.ActionEvent evt) {
        // TODO add your handling code here:
        runAnalysis = false;
        setVisible(false);
    }
    // Variables declaration - do not modify
    private javax.swing.ButtonGroup buttonGroup1;
    private javax.swing.JButton cancelButton;
    private javax.swing.JRadioButton ibdProbRadioButton;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JRadioButton lengthImputeRadioButton;
    private javax.swing.JTextField lengthTextField;
    private javax.swing.JTextField minProbTextField;
    private javax.swing.JTextField maxMismatchTextField;
    private javax.swing.JRadioButton majorAlleleRadioButton;
    private javax.swing.JButton runButton;
    private javax.swing.JRadioButton similarWindowRadioButton;
    // End of variables declaration
}
