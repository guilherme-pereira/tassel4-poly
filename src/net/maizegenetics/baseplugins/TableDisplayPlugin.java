/*
 * TableDisplayPlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.baseplugins;


import net.maizegenetics.pal.report.TableReport;
import net.maizegenetics.pal.report.TableReportUtils;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.baseplugins.AbstractDisplayPlugin.FileFormat;
import net.maizegenetics.gui.TableReportNoPagingTableModel;

import javax.swing.*;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.print.PageFormat;
import java.awt.print.Printable;
import java.awt.print.PrinterException;
import java.awt.print.PrinterJob;

import java.net.URL;

import java.util.Iterator;
import java.util.List;


/**
 *
 * @author Ed Buckler
 */
public class TableDisplayPlugin extends AbstractDisplayPlugin {

    private String myDelimiter;
    private TablePluginDialog myDialog = null;

    /** Creates a new instance of TableDisplayPlugin */
    public TableDisplayPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {

        try {
            List<Datum> alignInList = input.getDataOfType(TableReport.class);
            if (alignInList.size() < 1) {
                String message = "Invalid selection.  Please select Table Report.";
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), message);
                } else {
                    System.out.println(message);
                }
                return null;
            }
            Iterator<Datum> itr = alignInList.iterator();
            while (itr.hasNext()) {
                Datum current = itr.next();
                processDatum(current);
            }

            return null;
        } finally {
            fireProgress(100);
        }

    }

    private void processDatum(Datum input) {
        TableReport tr = (TableReport) input.getData();

        if (isInteractive()) {
            myDialog = new TablePluginDialog(this, tr);
            myDialog.setLocationRelativeTo(getParentFrame());
            myDialog.setVisible(true);
        } else if (getSaveFile() != null) {
            saveDataToFile(tr, myDelimiter);
        }
    }

    public void saveDataToFile(TableReport tr, String delimit, FileFormat[] formats) {
        TableReportUtils.saveDelimitedTableReport(tr, delimit, getSaveFileByChooser(formats, myDialog));
    }

    public void saveDataToFile(TableReport tr, String delimit) {
        TableReportUtils.saveDelimitedTableReport(tr, delimit, getSaveFile());
    }

    public String getDelimiter() {
        return myDelimiter;
    }

    public void setDelimiter(String theDelimiter) {
        myDelimiter = theDelimiter;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        URL imageURL = TableDisplayPlugin.class.getResource("images/Table.gif");
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
        return "Table";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Present data in table";
    }

}

/**
 * Title:        TASSEL
 * Description:  A java program to deal with diversity
 * Copyright:    Copyright (c) 2000
 * Company:      USDA-ARS/NCSU
 * @author Ed Buckler
 * @version 1.0
 */
class TablePluginDialog extends JDialog implements Printable {

    net.maizegenetics.pal.report.TableReport theTableSource;
    TableDisplayPlugin theTableDisplayPlugin;
    JPanel panel1 = new JPanel();
    BorderLayout borderLayout1 = new BorderLayout();
    JScrollPane jScrollPane1 = new JScrollPane();
    JTable jTable;
    JPanel jPanel1 = new JPanel();
    JButton saveTabButton = new JButton();
    JButton saveCommaButton = new JButton();
    JButton printButton = new JButton();

    public TablePluginDialog(TableDisplayPlugin plugin, TableReport theTableSource) {
        super(plugin.getParentFrame(), theTableSource.getTableTitle(), false);
        theTableDisplayPlugin = plugin;
        this.theTableSource = theTableSource;
        jTable = new JTable(new TableReportNoPagingTableModel(theTableSource));
        jTable.setAutoCreateRowSorter(true);



        //Set up tool tips for column headers.
        jTable.getTableHeader().setToolTipText("Click to specify sorting; Control-Click to specify secondary sorting");

        try {
            jbInit();
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    void jbInit() throws Exception {
        panel1.setLayout(borderLayout1);
        jScrollPane1.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        jScrollPane1.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
        saveTabButton.setText("Export (Tab)");
        saveTabButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                saveTabButton_actionPerformed(e);
            }
        });
        saveCommaButton.setText("Export (CSV)");
        saveCommaButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                saveCommaButton_actionPerformed(e);
            }
        });
        printButton.setText("Print");
        printButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                printButton_actionPerformed(e);
            }
        });
        getContentPane().add(panel1);
        panel1.add(jScrollPane1, BorderLayout.CENTER);
        jScrollPane1.getViewport().add(jTable, null);
        panel1.add(jPanel1, BorderLayout.SOUTH);
        jPanel1.add(printButton, null);
        jPanel1.add(saveCommaButton, null);
        jPanel1.add(saveTabButton, null);
    }

    void saveTabButton_actionPerformed(ActionEvent e) {
        theTableDisplayPlugin.saveDataToFile(theTableSource, "\t", new FileFormat[]{FileFormat.txt});
    }

    void saveCommaButton_actionPerformed(ActionEvent e) {
        theTableDisplayPlugin.saveDataToFile(theTableSource, ",", new FileFormat[]{FileFormat.csv});
    }

    void printButton_actionPerformed(ActionEvent e) {
        sendToPrinter();
    }

    void sendToPrinter() {
        PrinterJob printJob = PrinterJob.getPrinterJob();
        printJob.setPrintable(this);
        if (printJob.printDialog()) {
            try {
                printJob.print();
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }
    }

    public int print(Graphics g, PageFormat pageFormat,
            int pageIndex) throws PrinterException {
        Graphics2D g2 = (Graphics2D) g;
        g2.setColor(Color.black);
        int fontHeight = g2.getFontMetrics().getHeight();
        int fontDesent = g2.getFontMetrics().getDescent();

        //leave room for page number
        double pageHeight = pageFormat.getImageableHeight() - fontHeight;
        double pageWidth = pageFormat.getImageableWidth();
        double tableWidth = (double) jTable.getColumnModel().getTotalColumnWidth();
        double scale = 1;
        if (tableWidth >= pageWidth) {
            scale = pageWidth / tableWidth;
        }

        double headerHeightOnPage =
                jTable.getTableHeader().getHeight() * scale;
        double tableWidthOnPage = tableWidth * scale;

        double oneRowHeight = (jTable.getRowHeight()
                + jTable.getRowMargin()) * scale;
        int numRowsOnAPage =
                (int) ((pageHeight - headerHeightOnPage) / oneRowHeight);
        double pageHeightForTable = oneRowHeight * numRowsOnAPage;
        int totalNumPages = (int) Math.ceil(((double) jTable.getRowCount()) / numRowsOnAPage);
        if (pageIndex >= totalNumPages) {
            return NO_SUCH_PAGE;
        }

        g2.translate(pageFormat.getImageableX(),
                pageFormat.getImageableY());
        g2.drawString("Page: " + (pageIndex + 1), (int) pageWidth / 2 - 35,
                (int) (pageHeight + fontHeight - fontDesent));//bottom center

        g2.translate(0f, headerHeightOnPage);
        g2.translate(0f, -pageIndex * pageHeightForTable);

        //If this piece of the table is smaller than the size available,
        //clip to the appropriate bounds.
        if (pageIndex + 1 == totalNumPages) {
            int lastRowPrinted = numRowsOnAPage * pageIndex;
            int numRowsLeft = jTable.getRowCount() - lastRowPrinted;
            g2.setClip(0, (int) (pageHeightForTable * pageIndex),
                    (int) Math.ceil(tableWidthOnPage),
                    (int) Math.ceil(oneRowHeight * numRowsLeft));
        } //else clip to the entire area available.
        else {
            g2.setClip(0, (int) (pageHeightForTable * pageIndex),
                    (int) Math.ceil(tableWidthOnPage),
                    (int) Math.ceil(pageHeightForTable));
        }

        g2.scale(scale, scale);
        jTable.paint(g2);
        g2.scale(1 / scale, 1 / scale);
        g2.translate(0f, pageIndex * pageHeightForTable);
        g2.translate(0f, -headerHeightOnPage);
        g2.setClip(0, 0, (int) Math.ceil(tableWidthOnPage),
                (int) Math.ceil(headerHeightOnPage));
        g2.scale(scale, scale);
        jTable.getTableHeader().paint(g2);//paint header at top

        return Printable.PAGE_EXISTS;
    }
}
