/**
 * VerticalLabelUI
 */
package net.maizegenetics.gui;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import javax.swing.Icon;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.plaf.basic.BasicLabelUI;

public class VerticalLabelUI extends BasicLabelUI {

    private static final VerticalLabelUI INSTANCE = new VerticalLabelUI();

    private VerticalLabelUI() {
    }

    public static VerticalLabelUI getInstance() {
        return INSTANCE;
    }

    @Override
    public int getBaseline(JComponent c, int width, int height) {
        super.getBaseline(c, width, height);
        return -1;
    }

    @Override
    public Component.BaselineResizeBehavior getBaselineResizeBehavior(JComponent c) {
        super.getBaselineResizeBehavior(c);
        return Component.BaselineResizeBehavior.OTHER;
    }

    @Override
    protected String layoutCL(JLabel label, FontMetrics fontMetrics,
            String text, Icon icon, Rectangle viewR, Rectangle iconR,
            Rectangle textR) {

        Rectangle tempViewR = transpose(viewR, new Rectangle());
        Rectangle tempIconR = transpose(iconR, new Rectangle());
        Rectangle tempTextR = transpose(textR, new Rectangle());

        text = super.layoutCL(label, fontMetrics, text, icon, tempViewR, tempIconR, tempTextR);

        viewR = copy(tempViewR, viewR);
        iconR = copy(tempIconR, iconR);
        textR = copy(tempTextR, textR);
        return text;
    }

    @Override
    public void paint(Graphics g, JComponent c) {
        Graphics2D g2d = (Graphics2D) g.create();
        g2d.rotate(-Math.PI / 2, c.getSize().height / 2, c.getSize().height / 2);
        super.paint(g2d, c);
    }

    @Override
    public Dimension getPreferredSize(JComponent c) {
        return transposeDimension(super.getPreferredSize(c));
    }

    @Override
    public Dimension getMaximumSize(JComponent c) {
        return transposeDimension(super.getMaximumSize(c));
    }

    @Override
    public Dimension getMinimumSize(JComponent c) {
        return transposeDimension(super.getMinimumSize(c));
    }

    private Dimension transposeDimension(Dimension from) {
        return new Dimension(from.height, from.width);
    }

    private Rectangle transpose(Rectangle from, Rectangle to) {
        to.x = from.y;
        to.y = from.x;
        to.width = from.height;
        to.height = from.width;
        return to;
    }

    private Rectangle copy(Rectangle from, Rectangle to) {
        to.x = from.x;
        to.y = from.y;
        to.setSize(from.getSize());
        return to;
    }
}