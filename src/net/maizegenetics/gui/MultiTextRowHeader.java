package net.maizegenetics.gui;

import javax.swing.*;
import javax.swing.table.TableCellRenderer;
import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.awt.image.VolatileImage;

/**
 * Created by IntelliJ IDEA.
 * User: ed
 * Date: Nov 22, 2004
 * Time: 11:01:21 AM
 * To change this template use File | Settings | File Templates.
 */
public class MultiTextRowHeader extends JPanel implements TableCellRenderer {

    int textRows = 1;
    char delimiter = '.';
    JLabel[] theLabels;
    JLabel theLabel;
    VolatileImage offscreenImage;
    AffineTransform at;
    FontMetrics fm;
    Rectangle2D r2d;

    public MultiTextRowHeader(int textRows, char delimiter) {
        super();
        this.textRows = textRows;
        this.delimiter = delimiter;
        this.setLayout(new GridLayout(textRows, 1));
        if (textRows > 1) {
            theLabels = new JLabel[textRows];
        }
    }

    public Component getTableCellRendererComponent(JTable table, Object value,
            boolean isSelected, boolean hasFocus, int rowIndex, int vColIndex) {
        // 'value' is column header value of column 'vColIndex'
        // rowIndex is always -1
        // isSelected is always false
        // hasFocus is always false

        // Configure the component with the specified value
        String theText = value.toString();
        if (this.getComponentCount() == 0) {
            init(theText);
        }
        // Set tool tip if desired
        setToolTipText((String) value);

        // Since the renderer is a component, return itself
        return this;
    }

    private void init(String theText) {
        String theRemainder;
        removeAll();
        int begin = 0, end = 1000;
        if (textRows == 1) {
            theLabel = new JLabel(makeIcon(theText));
            this.add(theLabel);
        } else {
            for (int i = 0; i < textRows; i++) {
                if (end > 0) {
                    end = theText.indexOf(delimiter, begin);
                    if (end > 0) {
                        theRemainder = theText.substring(begin, end);
                    } else {
                        theRemainder = theText.substring(begin);
                    }
                    theLabels[i] = new JLabel(theRemainder);
                    begin = end + 1;
                } else {
                    theLabels[i] = new JLabel("");
                }
                this.add(theLabels[i]);
            }
        }
    }

    public ImageIcon makeIcon(String text) {
        int textLen = Math.min(text.length(), 20);
        text = text.substring(0, textLen);
        //BufferedImage takes lots of time and space, it has been replaced with volatile image
        if (offscreenImage == null) {
            BufferedImage testImage = new BufferedImage(60, 200, BufferedImage.TYPE_INT_RGB);    //60,200
            Graphics2D g2d = (Graphics2D) testImage.getGraphics();
            fm = g2d.getFontMetrics(new Font("Monospaced", Font.PLAIN, 12));
            Rectangle2D r2d = fm.getStringBounds("XXXXXXXXXXXXXXXXXXXXXXXX".substring(0, textLen + 1), g2d);
            GraphicsEnvironment ge = GraphicsEnvironment.getLocalGraphicsEnvironment();
            GraphicsConfiguration gc = ge.getDefaultScreenDevice().getDefaultConfiguration();
            offscreenImage = gc.createCompatibleVolatileImage((int) (r2d.getHeight() * 1.2), (int) (r2d.getWidth() * 1.2));
            at = new AffineTransform();
            // counter clockwise 90 degrees
            at.setToRotation(-Math.PI / 2);
        }
        Graphics2D g2d = (Graphics2D) offscreenImage.getGraphics();
        g2d.setColor(Color.WHITE);
        g2d.fillRect(0, 0, offscreenImage.getWidth(), offscreenImage.getHeight());
        g2d.setColor(Color.BLACK);
        g2d.setTransform(at);
        g2d.drawString(text, -(int) (offscreenImage.getHeight() * 0.9), (int) (offscreenImage.getWidth() * 0.9));
        ImageIcon small = new ImageIcon(offscreenImage);
        return small;
    }
}
