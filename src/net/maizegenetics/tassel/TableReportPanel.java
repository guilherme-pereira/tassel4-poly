package net.maizegenetics.tassel;

import net.maizegenetics.gui.TableReportNoPagingTableModel;
import net.maizegenetics.pal.report.TableReport;

import net.maizegenetics.util.DoubleFormat;

import javax.swing.*;
import javax.swing.table.*;
import java.awt.*;

/**
 */
public class TableReportPanel extends JPanel {

    private JTable myDataTable;
    private int myTaxaColumnIndex = -1;

    public TableReportPanel(TableReport theTableSource) {

        TableModel theModel = null;

        if (theTableSource instanceof TableModel) {
            theModel = (TableModel) theTableSource;
        } else {
            theModel = new TableReportNoPagingTableModel(theTableSource);
        }

        myDataTable = new JTable(theModel);

        myDataTable.setDefaultRenderer(Double.class, new DefaultTableCellRenderer() {
            public void setValue(Object value) {
                setText(DoubleFormat.format((Double) value));
            }
        });

        ((DefaultTableCellRenderer) myDataTable.getDefaultRenderer(String.class)).setHorizontalAlignment(JLabel.RIGHT);
        ((DefaultTableCellRenderer) myDataTable.getDefaultRenderer(Double.class)).setHorizontalAlignment(JLabel.RIGHT);
        ((DefaultTableCellRenderer) myDataTable.getDefaultRenderer(Integer.class)).setHorizontalAlignment(JLabel.RIGHT);

        myDataTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        TableRowSorter<TableModel> theSorter = new TableRowSorter<TableModel>(theModel);
        myDataTable.setRowSorter(theSorter);
        JScrollPane jsp = new JScrollPane(myDataTable);

        String[][] rowtable = getRowHeadings(theTableSource);
        if (myTaxaColumnIndex != -1) {

            JTable rowNameTable = new JTable(rowtable, new String[]{"Taxa"}) {
                public boolean isCellEditable(int row, int column) {
                    return false;
                }
            };
            rowNameTable.setRowSorter(theSorter);
            jsp.setRowHeaderView(rowNameTable);
            int numTaxaToCheck = Math.min(500, rowtable.length);
            int maxChars = 0;
            for (int i = 0; i < numTaxaToCheck; i++) {
                if (rowtable[i][0].length() > maxChars) {
                    maxChars = rowtable[i][0].length();
                }
            }
            int preferedWidth = Math.max(100, maxChars * 9);
            jsp.getRowHeader().setPreferredSize(new Dimension(preferedWidth, myDataTable.getHeight()));
            jsp.setCorner(JScrollPane.UPPER_LEFT_CORNER, rowNameTable.getTableHeader());

            TableColumn taxaColumn = myDataTable.getColumnModel().getColumn(myTaxaColumnIndex);
            myDataTable.getColumnModel().removeColumn(taxaColumn);
        }

        setLayout(new BorderLayout());
        add(jsp, BorderLayout.CENTER);

        setVisible(true);

    }

    private String[][] getRowHeadings(TableReport report) {

        Object[] headings = report.getTableColumnNames();
        for (int i = 0, n = headings.length; i < n; i++) {
            if (headings[i].toString().equalsIgnoreCase("taxa")) {
                myTaxaColumnIndex = i;
                break;
            }
        }

        if (myTaxaColumnIndex == -1) {
            return null;
        }

        String[][] result = new String[report.getRowCount()][1];
        for (int i = 0, n = report.getRowCount(); i < n; i++) {
            result[i][0] = report.getValueAt(i, myTaxaColumnIndex).toString();
        }

        return result;

    }
}
