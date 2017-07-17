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
// QuantitativeLegendComponent.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.pal.gui;

import net.maizegenetics.pal.alignment.Phenotype;

import java.awt.*;

import java.text.DecimalFormat;

/**
 * An AWT Component for displaying a quantitative legend
 *
 *
 *
 * @author Ed Buckler
 * @version $Id:
 */
public class QuantitativeLegendComponent extends Component {

    Color theColor = new Color(0, 0, 0);
    int hoff = 70, h2off = 70, voff = 20;
    double maximum, minimum;
    //hoff is on the left side for site labels
    //h2off is on the right side for legends

    public QuantitativeLegendComponent(double maximum, double minimum) {
        this.maximum = maximum;
        this.minimum = minimum;
        Dimension d = this.getSize();
//    System.out.println(d);
        try {
            jbInit();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        setVisible(true);
        d = this.getSize();
//    System.out.println(d);
    }

    /**
     * This determines what is displayed in the lower left corner.
     * Options are: P_VALUE, DPRIME, and RSQUARE
     */
    public void setRange(double max, double min) {
        if (max < min) {
            return;
        }
        maximum = max;
        minimum = min;
    }

    public void setTrait(Phenotype ca, int trait) {
        maximum = -9999999;
        minimum = 9999999;
        double ct;
        for (int i = 0; i < ca.getNumberOfTaxa(); i++) {
            ct = ca.getData(i, trait);
            if (ct == Phenotype.MISSING) {
                continue;
            }
            if (ct < minimum) {
                minimum = ct;
            }
            if (ct > maximum) {
                maximum = ct;
            }
        }
    }

    private void jbInit() throws Exception {
        this.setBackground(Color.red);
        this.setSize(50, 400);
//    this.setPreferredSize(new Dimension(400, 400));
//    this.setLayout(borderLayout1);
    }

    private Color getColor(double value) {
        double val = (value - minimum) / (maximum - minimum);
        //return Color.getHSBColor((float)val, 1.0f, 1.0f);  //color scale
        return Color.getHSBColor(0.0f, 0.0f, (float) val);  //black and white
    }

    protected void paintComponent(Graphics g) {
        addLegend(g);
    }

    public void paint(Graphics g) {
        paintComponent(g);
    }

    /**
     * Returns the preferred size for drawing
     * (that is the size that will show everything nicely)
     */
    public Dimension getPreferredSize() {
        return new Dimension(80, 400);
    }

    private void addLegend(Graphics g) {
        Dimension d = this.getSize();
//    System.out.println(d);
//    setSize(50,300);
        int localX = 5;
        int mid = (int) (d.height * 0.9);
        g.setColor(Color.white);
        g.fillRect(0, 0, d.width, d.height);
        g.setColor(Color.black);
        g.drawString("Values ", localX, 10);
        addLegendGraph(g, false, 10, 20, mid - 10);
    }

    private void addLegendGraph(Graphics g, boolean prob, int xStart, int yStart, int yEnd) {
        DecimalFormat dF; //=new DecimalFormat("0.0000");
        int yInc, currY = yStart;
        int barWidth = 10;
        if (prob) {
            yInc = (yEnd - yStart) / 4;
            g.setColor(Color.white);
            g.fillRect(xStart, currY, barWidth, yInc);
            g.setColor(Color.black);
            g.drawRect(xStart, currY, barWidth, yInc);
            g.drawString(">0.01", xStart + barWidth + 5, currY + 10);
            currY += yInc;
            g.setColor(Color.blue);
            g.fillRect(xStart, currY, barWidth, yInc);
            g.setColor(Color.black);
            g.drawRect(xStart, currY, barWidth, yInc);
            g.drawString("<0.01", xStart + barWidth + 5, currY + 10);
            currY += yInc;
            g.setColor(Color.green);
            g.fillRect(xStart, currY, barWidth, yInc);
            g.setColor(Color.black);
            g.drawRect(xStart, currY, barWidth, yInc);
            g.drawString("<0.001", xStart + barWidth + 5, currY + 10);
            currY += yInc;
            g.setColor(Color.red);
            g.fillRect(xStart, currY, barWidth, yInc);
            g.setColor(Color.black);
            g.drawRect(xStart, currY, barWidth, yInc);
            g.drawString("<0.0001", xStart + barWidth + 5, currY + 10);
        } else {
            yInc = (yEnd - yStart) / 11;
            dF = new DecimalFormat("0.00");
            double inc = (maximum - minimum) / 11;
            for (double d = maximum; d >= minimum; d -= inc) {
                g.setColor(getColor(d));
                g.fillRect(xStart, currY, barWidth, yInc);
                g.setColor(Color.black);
                g.drawRect(xStart, currY, barWidth, yInc);
                g.drawString(dF.format(d), xStart + barWidth + 5, currY + 10);
                currY += yInc;
            }
        }
    }
}