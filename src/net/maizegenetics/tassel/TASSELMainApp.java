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
// Title:      TASSELMainApp
// Version:
// Copyright:  Copyright (c) 1998
// Author:     Ed Buckler
package net.maizegenetics.tassel;

import javax.swing.*;
import java.awt.*;

import net.maizegenetics.pipeline.TasselPipeline;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.util.ExceptionUtils;

import org.apache.log4j.PropertyConfigurator;

public class TASSELMainApp {

    private final TASSELMainFrame frame;

    //Construct the application
    public TASSELMainApp() {

        java.util.Properties props = new java.util.Properties();
        props.setProperty("log4j.logger.net.maizegenetics", "INFO, stdout");
        props.setProperty("log4j.appender.stdout",
                "org.apache.log4j.ConsoleAppender");
        props.setProperty("log4j.appender.stdout.layout",
                "org.apache.log4j.TTCCLayout");
        PropertyConfigurator.configure(props);

        try {
            UIManager.setLookAndFeel(new com.sun.java.swing.plaf.windows.WindowsLookAndFeel());
        } catch (Exception e) {
        }

        TasselPrefs.setPersistPreferences(true);

        frame = new TASSELMainFrame();
        frame.validate();

        //Center the window
        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
        Dimension frameSize = frame.getSize();
        if (frameSize.height > screenSize.height) {
            frameSize.height = screenSize.height;
        }
        if (frameSize.width > screenSize.width) {
            frameSize.width = screenSize.width;
        }
        frame.setLocation((screenSize.width - frameSize.width) / 2, (screenSize.height - frameSize.height) / 2);
        frame.setVisible(true);
    }

    public TASSELMainFrame getTASSELMainFrame() {
        return frame;
    }

    //Main method
    static public void main(String[] args) {
        try {
            TASSELMainApp mainApp = new TASSELMainApp();
            if (args.length > 0) {
                new TasselPipeline(args, mainApp.getTASSELMainFrame());
            }
        } catch (Throwable e) {
            String userMessage = "TASSEL has experienced an error.  "
                    + ExceptionUtils.getExceptionCauses(e);
            if (e instanceof java.lang.OutOfMemoryError) {
                userMessage = "You have used up all of the memory allocated to the Java Virtual Machine.  "
                        + "It is recommneded that you adjust your heap settings and possibly add more memory to the computer.  "
                        + "Additionally, some operations are not recommended on a full dataset, i.e., select SNPs *before* determining LD";
            }
            JOptionPane.showMessageDialog(null, userMessage, "Fatal Error", JOptionPane.ERROR_MESSAGE);
            e.printStackTrace();
        }
    }
}
