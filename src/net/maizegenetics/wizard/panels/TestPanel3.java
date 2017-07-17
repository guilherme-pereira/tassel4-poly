package net.maizegenetics.wizard.panels;

import javax.swing.*;
import javax.swing.border.*;
import java.awt.*;
import java.awt.event.*;

import net.maizegenetics.wizard.*;

public class TestPanel3 extends JPanel {
 
    private JLabel anotherBlankSpace;
    private JLabel blankSpace;
    private ButtonGroup connectorGroup;
    private JLabel jLabel1;
    private JPanel jPanel1;
//    private JLabel progressDescription;
//    private JProgressBar progressSent;
    private JLabel welcomeTitle;
    private JLabel yetAnotherBlankSpace1;
    private JLabel yetAnotherBlankSpace2;
    
    private JPanel contentPanel;
    private JLabel iconLabel;
    private JSeparator separator;
    private JLabel textLabel;
    private JPanel titlePanel;

    private JCheckBox ldCheckBox = new JCheckBox("LD Analysis");
    private JCheckBox glmCheckBox = new JCheckBox("GLM");
    private JCheckBox mlmCheckBox = new JCheckBox("MLM");
    private JCheckBox imputeCheckBox = new JCheckBox("Impute");

//    private JCheckBox ldPlotCheckBox = new JCheckBox("LD plot");

        
    public TestPanel3() {
        
        super();
                
        contentPanel = getContentPanel();
        ImageIcon icon = getImageIcon();
        
        titlePanel = new javax.swing.JPanel();
        textLabel = new javax.swing.JLabel();
        iconLabel = new javax.swing.JLabel();
        separator = new javax.swing.JSeparator();

        setLayout(new java.awt.BorderLayout());

        titlePanel.setLayout(new java.awt.BorderLayout());
//        titlePanel.setBackground(Color.gray);
//
//        textLabel.setBackground(Color.gray);
        textLabel.setFont(new Font("MS Sans Serif", Font.BOLD, 14));
        textLabel.setText("Select Analysis Types");
        textLabel.setBorder(new EmptyBorder(new Insets(10, 10, 10, 10)));
        textLabel.setOpaque(true);

//        iconLabel.setBackground(Color.gray);
//        if (icon != null)
//            iconLabel.setIcon(icon);
        
        titlePanel.add(textLabel, BorderLayout.CENTER);
        titlePanel.add(iconLabel, BorderLayout.EAST);
        titlePanel.add(separator, BorderLayout.SOUTH);

        add(titlePanel, BorderLayout.NORTH);
        JPanel secondaryPanel = new JPanel();
        secondaryPanel.add(contentPanel, BorderLayout.NORTH);
        add(secondaryPanel, BorderLayout.WEST);
        
    }
    
    public boolean isLDChecked() {
        return ldCheckBox.isSelected();
    }

    public boolean isGLMChecked() {
        return glmCheckBox.isSelected();
    }

    public boolean isMLMChecked() {
        return mlmCheckBox.isSelected();
    }

    public boolean isImputeChecked() {
        return imputeCheckBox.isSelected();
    }

//    public boolean isLDPlotChecked() {
//        return ldPlotCheckBox.isSelected();
//    }
    
//    public void setProgressText(String s) {
//        progressDescription.setText(s);
//    }
    
//    public void setProgressValue(int i) {
//        progressSent.setValue(i);
//    }
    
    private JPanel getContentPanel() {            
        
        JPanel contentPanel1 = new JPanel();
        
        connectorGroup = new ButtonGroup();
        welcomeTitle = new JLabel();
        jPanel1 = new JPanel();
        blankSpace = new JLabel();
//        progressSent = new JProgressBar();
//        progressDescription = new JLabel();
        anotherBlankSpace = new JLabel();
        yetAnotherBlankSpace1 = new JLabel();
        yetAnotherBlankSpace2 = new JLabel();
        jLabel1 = new JLabel();

        contentPanel1.setLayout(new java.awt.BorderLayout());

        welcomeTitle.setText("Please select the types of analysis you wish to perform");
        contentPanel1.add(welcomeTitle, java.awt.BorderLayout.NORTH);

        jPanel1.setLayout(new java.awt.GridLayout(0, 1));

        jPanel1.add(blankSpace);
        jPanel1.add(ldCheckBox);
        //jPanel1.add(glmCheckBox);
        // Delete line when functionality is added
        glmCheckBox.setVisible(false);
        //jPanel1.add(mlmCheckBox);
        // Delete line when functionality is added
        mlmCheckBox.setVisible(false);
        jPanel1.add(imputeCheckBox);

        //jPanel1.add(blankSpace);

//        jPanel1.add(ldPlotCheckBox);

        jPanel1.add(yetAnotherBlankSpace2);

//        progressSent.setStringPainted(true);
//        jPanel1.add(progressSent);

//        progressDescription.setFont(new java.awt.Font("MS Sans Serif", 1, 11));
//        progressDescription.setText("Connecting to Server...");
//        jPanel1.add(progressDescription);

        jPanel1.add(anotherBlankSpace);

        jPanel1.add(yetAnotherBlankSpace1);

        contentPanel1.add(jPanel1, java.awt.BorderLayout.CENTER);

        jLabel1.setText("Hit \"Next\" after desired analysis types have been selected.");
        contentPanel1.add(jLabel1, java.awt.BorderLayout.SOUTH);
        
        return contentPanel1;
    }
    
    private ImageIcon getImageIcon() {        
        return null;
    }

    public void addCheckBoxActionListener(ActionListener l) {
        ldCheckBox.addActionListener(l);
        glmCheckBox.addActionListener(l);
        mlmCheckBox.addActionListener(l);
        imputeCheckBox.addActionListener(l);
    }

    public void addLDCheckBoxActionListener(ActionListener l) {
        ldCheckBox.addActionListener(l);
    }

    public void addGLMCheckBoxActionListener(ActionListener l) {
        glmCheckBox.addActionListener(l);
    }

    public void addMLMCheckBoxActionListener(ActionListener l) {
        mlmCheckBox.addActionListener(l);
    }

    public void addImputeCheckBoxActionListener(ActionListener l) {
        imputeCheckBox.addActionListener(l);
    }
}
