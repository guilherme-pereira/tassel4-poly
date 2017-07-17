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
package net.maizegenetics.gui;

import java.awt.*;
import java.io.EOFException;
import java.io.LineNumberReader;
import java.io.StringReader;

// Print a file into the text area.
public class PrintTextArea {

    Frame theFrame;

    public PrintTextArea(Frame frame) {
        theFrame = frame;
    }

    public void printThis(String s) {
        PrintJob pjob = theFrame.getToolkit().getPrintJob(theFrame, "Cool Stuff", null);
        // PrintJob pjob = new PrintJob();
        if (pjob != null) {
            Graphics pg = pjob.getGraphics();
            if (pg != null) {
                //String s = textArea.getText();
                printLongString(pjob, pg, s);
                pg.dispose();
            }
            pjob.end();
        }
    }

    // Print string to graphics via printjob
    // Does not deal with word wrap or tabs
    void printLongString(PrintJob pjob, Graphics pg, String s) {
        int pageNum = 1;
        int linesForThisPage = 0;
        int linesForThisJob = 0;
        int topMargin = 20, leftMargin = 20;
        // Note: String is immutable so won't change while printing.
        if (!(pg instanceof PrintGraphics)) {
            throw new IllegalArgumentException("Graphics context not PrintGraphics");
        }
        StringReader sr = new StringReader(s);
        LineNumberReader lnr = new LineNumberReader(sr);
        String nextLine;
        int pageHeight = pjob.getPageDimension().height - (topMargin * 2);
        Font helv = new Font("Helvetica", Font.PLAIN, 10);
        //have to set the font to get any output
        pg.setFont(helv);
        FontMetrics fm = pg.getFontMetrics(helv);
        int fontHeight = fm.getHeight();
        int fontDescent = fm.getDescent();
        int curHeight = topMargin;
        try {
            do {
                nextLine = lnr.readLine();
                if (nextLine != null) {
                    if ((curHeight + fontHeight) > pageHeight) {
                        // New Page
                        System.out.println("" + linesForThisPage + " lines printed for page " + pageNum);
                        pageNum++;
                        linesForThisPage = 0;
                        pg.dispose();
                        pg = pjob.getGraphics();
                        if (pg != null) {
                            pg.setFont(helv);
                        }
                        curHeight = topMargin;
                    }
                    curHeight += fontHeight;
                    if (pg != null) {
                        pg.drawString(nextLine, leftMargin, curHeight - fontDescent);
                        linesForThisPage++;
                        linesForThisJob++;
                    } else {
                        System.out.println("pg null");
                    }
                }
            } while (nextLine != null);
        } catch (EOFException eof) {
            // Fine, ignore
        } catch (Throwable t) { // Anything else
            t.printStackTrace();
        }
        System.out.println("" + linesForThisPage + " lines printed for page " + pageNum);
        System.out.println("pages printed: " + pageNum);
        System.out.println("total lines printed: "
                + linesForThisJob);
    }
}
