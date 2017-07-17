/*
 * TASSEL - Trait Analysis by a aSSociation Evolution & Linkage
 * Copyright (C) 2003 Ed Buckler
 *
 * This software evaluates linkage disequilibrium nucletide diversity and 
 * associations. For more information visit http://www.maizegenetics.net
 *
 * This software is distributed under GNU general public license and without
 * any warranty ot technical support.
 *
 * You can redistribute and/or modify it under the terms of GNU General 
 * public license. 
 *
 */
/**
 * Title:        QEX (Quantitative Analysis of Expression)<p>
 * Description:  <p>
 * Copyright:    Copyright (c) Edward Buckler<p>
 * Company:      <p>
 * @author Edward Buckler
 * @version 1.0
 */
/**
 * Title:        QEX (Quantitative Analysis of Expression)<p>
 * Description:  <p>
 * Copyright:    Copyright (c) Edward Buckler<p>
 * Company:      <p>
 * @author Edward Buckler
 * @version 1.0
 */
package net.maizegenetics.tassel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowEvent;
import java.net.URL;

public class AboutBox extends JDialog implements ActionListener {

    JPanel panel1 = new JPanel();
    JPanel panel2 = new JPanel();
    JPanel insetsPanel1 = new JPanel();
    JPanel insetsPanel2 = new JPanel();
    JPanel textPanel = new JPanel();
    JButton button1 = new JButton();
    JLabel imageControl1 = new JLabel();
    ImageIcon imageIcon;
    BorderLayout borderLayout1 = new BorderLayout();
    BorderLayout borderLayout2 = new BorderLayout();
    FlowLayout flowLayout1 = new FlowLayout();
    FlowLayout flowLayout2 = new FlowLayout();
    GridLayout gridLayout1 = new GridLayout();
    private String product = "TASSEL";
    private String version;
    private String comments = "Trait Analysis by Association, Evolution and Linkage";
    private String contact = "Contacts...";
    private String contact1 = "General Questions / Pipeline - Terry Casstevens (tmc46@cornell.edu)";
    private String contact2 = "Analysis Tools / MLM / GLM - Peter Bradbury (pjb39@cornell.edu), Zhiwu Zhang (zz19@cornell.edu)";
    private String contact3 = "Overall Project Lead - Ed Buckler (esb33@cornell.edu)";
    private String programmers = "Produced by Edward Buckler, Peter Bradbury, Dallas Kroon, Yogesh Ramdoss, AJ Fink, Zhiwu Zhang, and Terry Casstevens";
    private String acknowledgements = "Libraries used: PAL, Batik, COLT, JFreeChart, and others";
    private String copyright = "Copyright 2004 (c) by Edward Buckler";

    public AboutBox(TASSELMainFrame parent) {
        super(parent);

        version = "Version " + parent.version + " on " + parent.versionDate;

        enableEvents(AWTEvent.WINDOW_EVENT_MASK);
        try {
            jbInit();
        } catch (Exception e) {
            e.printStackTrace();
        }
        if (imageIcon != null) {
            imageControl1.setIcon(imageIcon);
        }
        pack();
    }

    private void jbInit() throws Exception {
        URL imageURL = getClass().getResource("images/AboutPicture.GIF");

        if (imageURL != null) {
            imageIcon = new ImageIcon(imageURL);
        }
        this.setTitle("About");
        setResizable(false);
        panel1.setLayout(borderLayout1);
        panel2.setLayout(borderLayout2);
        insetsPanel1.setLayout(flowLayout1);
        insetsPanel2.setLayout(flowLayout1);
        insetsPanel2.setBorder(new EmptyBorder(10, 10, 10, 10));
        gridLayout1.setRows(9);
        gridLayout1.setColumns(1);
        JLabel label1 = new JLabel(version);
        JLabel label2 = new JLabel(comments);
        JLabel contactLabel = new JLabel(contact);
        JLabel contact1Label = new JLabel(contact1);
        JLabel contact2Label = new JLabel(contact2);
        JLabel contact3Label = new JLabel(contact3);
        JLabel label3 = new JLabel(programmers);
        JLabel label4 = new JLabel(acknowledgements);
        JLabel label5 = new JLabel(copyright);

        textPanel.setLayout(gridLayout1);
        textPanel.setBorder(new EmptyBorder(10, 60, 10, 10));
        button1.setText("OK");
        button1.addActionListener(this);
        insetsPanel2.add(imageControl1, null);
        panel2.add(insetsPanel2, BorderLayout.WEST);
        this.getContentPane().add(panel1, null);
        textPanel.add(label1, null);
        textPanel.add(label2, null);
        textPanel.add(contactLabel, null);
        textPanel.add(contact1Label, null);
        textPanel.add(contact2Label, null);
        textPanel.add(contact3Label, null);
        textPanel.add(label3, null);
        textPanel.add(label4, null);
        textPanel.add(label5, null);
        panel2.add(textPanel, BorderLayout.CENTER);
        insetsPanel1.add(button1, null);
        panel1.add(insetsPanel1, BorderLayout.SOUTH);
        panel1.add(panel2, BorderLayout.NORTH);
    }

    protected void processWindowEvent(WindowEvent e) {
        if (e.getID() == WindowEvent.WINDOW_CLOSING) {
            cancel();
        }
        super.processWindowEvent(e);
    }

    void cancel() {
        dispose();
    }

    public void actionPerformed(ActionEvent e) {
        if (e.getSource() == button1) {
            cancel();
        }
    }
}
