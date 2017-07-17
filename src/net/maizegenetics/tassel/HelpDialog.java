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
package net.maizegenetics.tassel;

import javax.swing.*;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;
import java.awt.*;
import java.io.IOException;
import java.net.URL;

//Designed and programmed by Dr. Edward Buckler and his Bioinformatics Team: Peter Bradbury, Terry Casstevens, Chunguang Du, Dallas Kroon, Jack Liu, David Remington, Jeff Thornsberry, and Zhiwu Zhang.
public class HelpDialog extends JDialog {

    private JEditorPane htmlPane;

    public HelpDialog(Frame frame) {
        super(frame, "TASSEL Help", false);

        //Create the HTML viewing pane.
        htmlPane = new JEditorPane();
        htmlPane.setEditable(false);
        initHelp();
        JScrollPane htmlView = new JScrollPane(htmlPane);
        Dimension minimumSize = new Dimension(400, 400);
        htmlView.setMinimumSize(minimumSize);
        htmlView.setPreferredSize(new Dimension(600, 400));

        //Add the split pane to this frame.
        getContentPane().add(htmlView, BorderLayout.CENTER);

        pack();
    }

    private void initHelp() {
        String s = "Home.html";
        try {
            displayURL(s);
        } catch (Exception e) {
            System.err.println("Couldn't create net.maizegenetics.help URL: " + s);
        }
    }

    private void displayURL(String url) {
        try {
            htmlPane.setPage(HelpDialog.class.getResource(url));
            htmlPane.addHyperlinkListener(new HyperlinkListener() {
                public void hyperlinkUpdate(HyperlinkEvent hyperlinkEvent) {
                    HyperlinkEvent.EventType type = hyperlinkEvent.getEventType();
                    final URL url = hyperlinkEvent.getURL();
                    if (type == HyperlinkEvent.EventType.ENTERED) {
                        System.out.println("URL: " + url);
                    } else if (type == HyperlinkEvent.EventType.ACTIVATED) {
                        System.out.println("Activated");

                        //do some thing here
                        try {
                            Runtime.getRuntime().exec("rundll32 url.dll,FileProtocolHandler " + url);
                            //Runtime.getRuntime().exec("rundll32 url.dll,FileProtocolHandler " + "www.cnn.com");
                            //"C:\\Maize\\tassel\\src\\net\\maizegenetics\\help\\LoadSSR.gif");
                            //"C:\\Maize\\tassel\\src\\net\\maizegenetics\\help\\Overview.html");
                        } catch (Exception er) {
                            er.printStackTrace();
                        }

                    }
                }
            });

        } catch (IOException e) {
            System.err.println("Attempted to read a bad URL: " + url);
        }
    }
}
