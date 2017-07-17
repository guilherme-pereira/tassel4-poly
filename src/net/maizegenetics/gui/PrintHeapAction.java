/*
 * PrintHeapAction.java
 *
 * Created on June 18, 2010
 */
package net.maizegenetics.gui;

import java.awt.Frame;
import java.awt.event.ActionEvent;

import java.net.URL;
import java.util.*;

import javax.swing.*;

import java.text.NumberFormat;

import net.maizegenetics.util.Sizeof;
import net.maizegenetics.util.Utils;

/**
 *
 * @author terry
 */
public class PrintHeapAction extends AbstractAction implements Runnable {

    private static final Map myInstances = new HashMap(1);
    public static final String TEXT = "Print Heap";
    private final Frame myParentFrame;
    private long myLastRecordedHeap = 0;
    private static final NumberFormat myFormat = NumberFormat.getInstance();

    static {
        myFormat.setMaximumFractionDigits(2);
    }

    private PrintHeapAction(Frame parentFrame) {
        super(TEXT, getIcon());
        myParentFrame = parentFrame;
    }

    public static PrintHeapAction getInstance(Frame parentFrame) {
        PrintHeapAction result = (PrintHeapAction) myInstances.get(parentFrame);

        if (result == null) {
            result = new PrintHeapAction(parentFrame);
            myInstances.put(parentFrame, result);
        }

        return result;
    }

    public void actionPerformed(ActionEvent e) {
        Thread thread = new Thread(this);
        thread.start();
    }

    /**
     * Method executed when this thread starts that removes all items from the
     * associated tree.
     */
    public void run() {

        StringBuffer buffer = new StringBuffer();
        try {
            buffer.append("Current Heap Size: ");
            long current = Sizeof.getMemoryUse();
            String currentStr = myFormat.format(current / 1048576l);
            buffer.append(currentStr);
            buffer.append(" MB");
            buffer.append("\nMax Available Heap: ");
            buffer.append(Utils.getMaxHeapSizeMB());
            buffer.append(" MB");
            buffer.append("\nDelta Since Last Record: ");
            String delta = myFormat.format((current - myLastRecordedHeap) / 1048576l);
            buffer.append(delta);
            buffer.append(" MB\n");
            myLastRecordedHeap = current;
        } catch (Exception ex) {
            buffer.append("\n");
            buffer.append(ex.getMessage());
        }

        String str = buffer.toString();

        JOptionPane dialog = new JOptionPane(str, JOptionPane.INFORMATION_MESSAGE);
        dialog.showMessageDialog(myParentFrame, str, "Heap Used", JOptionPane.INFORMATION_MESSAGE);

    }

    private static ImageIcon getIcon() {
        URL imageURL = PrintHeapAction.class.getResource("MemoryUse16.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }
}
