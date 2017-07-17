/*
 * Grid2dDisplayPlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.baseplugins;

import net.maizegenetics.pal.report.TableReport;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableColumnModel;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableColumn;
import javax.swing.border.Border;
import javax.swing.border.LineBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.*;
import java.net.URL;
import java.util.*;
import java.util.List;

import ext.swing.ZToolBar;
import ext.swing.ZTableModel;
import ext.swing.ZDialog;
import ext.util.ZDoubleUtil;

/**
 *
 * @author Ed Buckler
 */
public class Grid2dDisplayPlugin extends AbstractDisplayPlugin {

    String defaultRow = "Site";
    String defaultCol = "Environment";
    String defaultValue = "PermuteP";

    /** Creates a new instance of Grid2dDisplayPlugin */
    public Grid2dDisplayPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {
        List<Datum> tableInList = input.getDataOfType(TableReport.class);
        for (int i = 0; i < tableInList.size(); i++) {
            TableReport theTableReport = (TableReport) tableInList.get(i).getData();
            if (isInteractive()) {
                Grid2DDialog myDialog = new Grid2DDialog(this, theTableReport);
                myDialog.setRowComboBox(this.defaultRow);
                myDialog.setColumnComboBox(this.defaultCol);
                myDialog.setValueComboBox(this.defaultValue);
                myDialog.setLocationRelativeTo(getParentFrame());
                myDialog.setVisible(true);
            } else if (getSaveFile() != null) {
                System.out.println("Grid2dDisplayPlugin not fully functional");
                Grid2DDialog myDialog = new Grid2DDialog(this, theTableReport);
                myDialog.setRowComboBox(this.defaultRow);
                myDialog.setColumnComboBox(this.defaultCol);
                myDialog.setValueComboBox(this.defaultValue);
                myDialog.setLocationRelativeTo(getParentFrame());
                myDialog.setVisible(true);
                Component c = myDialog.getChart();
                c.setSize(getImageWidth(), getImageHeight());
                File tempFile = getSaveFile();
                if (tableInList.size() > 1) {
                    tempFile = new File(getSaveFile().getParent(), i + "_" + getSaveFile().getName());
                }
                saveDataToFile(c, tempFile);
            }
        }

        return null;
    }

    public String getDefaultRow() {
        return defaultRow;
    }

    public void setDefaultRow(String defaultRow) {
        this.defaultRow = defaultRow;
    }

    public String getDefaultCol() {
        return defaultCol;
    }

    public void setDefaultCol(String defaultCol) {
        this.defaultCol = defaultCol;
    }

    public String getDefaultValue() {
        return defaultValue;
    }

    public void setDefaultValue(String defaultValue) {
        this.defaultValue = defaultValue;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        URL imageURL = Grid2dDisplayPlugin.class.getResource("images/2DPlot.gif");
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
        return "2D Plot";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Display a 2-D Plot of Table Data";
    }
}

/**
 * 2-dimensional plot for table data
 * by Jack Liu
 */
class Grid2DDialog extends JDialog {

    private int selectedRow = 0, selectedColumn = 0, selectedZ = 0;
//      TASSELMainFrame theTASSELMainFrame;
    Grid2dDisplayPlugin theGrid2dDisplayPlugin;
    private Object tableHeader[], tableData[][];
    private boolean rowLabels = true, colLabels = true;
    private Vector headerVector;
//      private double[][] numericData = null;
    private String[] numericHeader = null;
    private double[][] matrix = null;
    private String[][] modelMatrix = null;
    private String[] rowHeader = null;
    private String[] columnHeader = null;
    private int curCellSize = 20;
    private Color fgColor = new Color(255, 255, 204);
    private Color bgColor = Color.blue;
    private Color color1 = new Color(0, 0, 0);
    private Color color2 = new Color(128, 128, 128);
    private Color color3 = new Color(196, 196, 196);
    private Color color4 = new Color(255, 255, 255);
    private float min = 0;
    private float max = 0;
    private float range = 0;
    private float cutoff = 0;
    private Grid2DDialog.MyRenderer rdr = new Grid2DDialog.MyRenderer();
    private JPanel pHelp = new Grid2DDialog.pValueHelp();
    // JBuilder components
    private JFileChooser fc = new JFileChooser();
    private JPanel panel = new JPanel();
    private BorderLayout borderLayout1 = new BorderLayout();
    private JScrollPane LP_chart = new JScrollPane();
    private JLabel statusBar = new JLabel();
    private JTable chart = new JTable();
    private JPanel P_toolBar = new JPanel();
    private ZToolBar toolBar = new ZToolBar();
    private JPanel rightToolPanel = new JPanel();
    private JComboBox CB_cellSize = new JComboBox();
    private JPanel cutoffBar = new JPanel();
    private JLabel L_cutoff = new JLabel();
    private JLabel L_cellSize = new JLabel();
    private JColorChooser cc = new JColorChooser();
    private Border cellBorder;
    private JTextField TF_cutoff = new JTextField();
    private JSlider S_cutoff = new JSlider();
    private JCheckBox C_pValueData = new JCheckBox();
    private JCheckBox C_upperTriangle = new JCheckBox();
    BorderLayout borderLayout2 = new BorderLayout();
    JPanel dataPanel = new JPanel();
    JComboBox columnComboBox = new JComboBox();
    JLabel jLabel1 = new JLabel();
    JComboBox rowComboBox = new JComboBox();
    JLabel jLabel2 = new JLabel();
    JComboBox valueComboBox = new JComboBox();
    JLabel jLabel3 = new JLabel();
    JCheckBox rowLabelCheckBox = new JCheckBox();
    JCheckBox colLabelCheckBox = new JCheckBox();
    JButton svgButton = new JButton();
    GridBagLayout gridBagLayout1 = new GridBagLayout();
//        JComboBox formatComboBox;

    /** My renderer for cell */
    private class MyRenderer extends DefaultTableCellRenderer {

        private Border border1 =
                new LineBorder(new Color(214, 214, 214), 1);
        private Border border2 = new LineBorder(Color.white, 1);

        public Component getTableCellRendererComponent(JTable table,
                Object value, boolean isSelected, boolean hasFocus, int row,
                int column) {
            this.setBorder(border1);
            setToolTipText(rowHeader[row] + ":"
                    + columnHeader[column] + "("
                    + value.toString() + ")");
            if (rowLabels && (column == 0)) {
                setText(value.toString());
                setBackground(Color.white);
                return this;
            } else {
                setText("");
            }


            if ((C_upperTriangle.isSelected() == true)
                    && (row > column)) {
                this.setBorder(border2);
                setBackground(Color.white);

                return this;
            }
            if (Double.isNaN(matrix[row][column])) {
                setBackground(Color.lightGray);
                return this;
            }
            if (C_pValueData.isSelected() == false) {
                if (matrix[row][column] >= cutoff) {
                    setBackground(fgColor);
                } else {
                    setBackground(bgColor);
                }
            } else {
                if (matrix[row][column] <= 0.001) {
                    setBackground(color1);
                } else if (matrix[row][column] <= 0.01) {
                    setBackground(color2);
                } else if (matrix[row][column] <= 0.05) {
                    setBackground(color3);
                } else {
                    setBackground(color4);
                }
            }

            return this;
        }
    }

    /** Help panel */
    private class pValueHelp extends JPanel {

        public pValueHelp() {
            setLayout(new GridLayout(4, 1));

            JPanel p1 = new JPanel();
            JLabel label1 = new DefaultTableCellRenderer();

            label1.setPreferredSize(new Dimension(18, 18));
            label1.setBackground(color1);
            p1.add(label1);
            p1.add(new JLabel(" < 0.001"));
            p1.setLayout(new FlowLayout(0));
            add(p1);

            JPanel p2 = new JPanel();
            JLabel label2 = new DefaultTableCellRenderer();

            label2.setPreferredSize(new Dimension(18, 18));
            label2.setBackground(color2);
            p2.add(label2);
            p2.add(new JLabel(" 0.001 - 0.01"));
            p2.setLayout(new FlowLayout(0));
            add(p2);

            JPanel p3 = new JPanel();
            JLabel label3 = new DefaultTableCellRenderer();

            label3.setPreferredSize(new Dimension(18, 18));
            label3.setBackground(color3);
            p3.add(label3);
            p3.add(new JLabel(" 0.01 - 0.05"));
            p3.setLayout(new FlowLayout(0));
            add(p3);

            JPanel p4 = new JPanel();
            JLabel label4 = new DefaultTableCellRenderer();

            label4.setPreferredSize(new Dimension(18, 18));
            label4.setBackground(color4);
            p4.add(label4);
            p4.add(new JLabel(" > 0.05"));
            p4.setLayout(new FlowLayout(0));
            add(p4);
        }
    }

    /**
     * Simple constructor
     */
    public Grid2DDialog(Grid2dDisplayPlugin gdp, TableReport theTable) {
        super(gdp.getParentFrame(), "Plot Stuff", false);
        theGrid2dDisplayPlugin = gdp;
        try {
            jbInit();
            tableHeader = theTable.getTableColumnNames();
            headerVector = new Vector();
            tableData = theTable.getTableData();
            //this.setModal(true);
            this.setSize(new Dimension(600, 510));
            String[] newHeader = new String[theTable.getTableColumnNames().length];

            for (int i = 0; i < newHeader.length; i++) {
                newHeader[i] = (String) tableHeader[i];
                headerVector.add(tableHeader[i]);
            }
            setComboBoxes(newHeader, theTable.getTableData());
            // processAndSetData(theTable.getTableColumnNames(), theTable.getTableData());
            showAtCenter();
            //         setVisible(true);

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * This method finds the columns with numeric data
     */
    private boolean isNumericData(String[] header, Object[][] data) {
        if (header.length <= 1) {
            return false;
        }
        if (data.length < 1) {
            return false;
        }
        for (int i = 0; i < data.length; i++) {
            if (data[i].length != header.length) {
                return false;
            }
        }
        String[] tempColumnHeader = new String[header.length];
        //rowHeader = new String[data.length];
        double[][] tempColumnData = new double[header.length][data.length];
        int[] badDataCount = new int[header.length];
        //rowData = new double[rowHeader.length][columnHeader.length];
        try {
            for (int i = 0; i < data.length; i++) {
                for (int j = 0; j < header.length; j++) {
                    tempColumnHeader[j] = header[j];
                    try {
                        double d = Double.valueOf(data[i][j].toString()).doubleValue();
                        tempColumnData[j][i] = d;
                    } catch (NumberFormatException ex) {
                        tempColumnData[j][i] = Double.NaN;
                        badDataCount[j]++;
                        //System.out.println("data["+i+"]["+j+"].toString()="+data[i][j].toString());
                    }
                }
            }
        } catch (NumberFormatException ex) {
            return false;
        }
        int goodCol = 0;
        for (int j = 0; j < header.length; j++) {
            if (badDataCount[j] < (data.length * 0.9)) {
                goodCol++;
            }
        }
        numericHeader = new String[goodCol];
        //    rowHeader = new String[data.length];
        //     numericData = new double[goodCol][data.length];
        goodCol = 0;
        for (int j = 0; j < header.length; j++) {
            if (badDataCount[j] < (data.length / 2)) {
                numericHeader[goodCol] = tempColumnHeader[j];
//            numericData[goodCol]=tempColumnData[j];
                goodCol++;
            }
        }
        return true;
    }

    /** Set data, throw exception if invalid data */
    public void setComboBoxes(String[] header, Object[][] data) throws Exception {
        if (isNumericData(header, data) == false) {
            throw new Exception("Invalid data for 1-D plotting");
        }
        rowComboBox.setEnabled(false);
        rowComboBox.removeAllItems();
        for (int i = 0; i < header.length; i++) {
            rowComboBox.addItem(header[i]);
        }
        rowComboBox.setEnabled(true);
        rowComboBox.setSelectedIndex(0);
        columnComboBox.setEnabled(false);
        columnComboBox.removeAllItems();
        for (int i = 0; i < header.length; i++) {
            columnComboBox.addItem(header[i]);
        }
        columnComboBox.setEnabled(true);
        columnComboBox.setSelectedIndex(0);
        valueComboBox.setEnabled(false);
        valueComboBox.removeAllItems();
        for (int i = 0; i < numericHeader.length; i++) {
            valueComboBox.addItem(numericHeader[i]);
        }
        valueComboBox.setEnabled(true);
        valueComboBox.setSelectedIndex(0);
    }

    public boolean processAndSetData(Object[] header, Object[][] data) {

        if ((selectedRow == selectedZ) || (selectedColumn == selectedZ) || (selectedColumn == selectedRow)) {
            return false;
        }
        Vector rowVec = new Vector();
        Vector colVec = new Vector();

//The following codes carries out the pivot of the data
        for (int i = 0; i < data.length; i++) {
            if (!rowVec.contains(data[i][selectedRow])) {
                rowVec.add(data[i][selectedRow]);
            }
            if (!colVec.contains(data[i][selectedColumn])) {
                colVec.add(data[i][selectedColumn]);
            }
        }
        Collections.sort(rowVec, new net.maizegenetics.util.StringNumberComparator());
        Collections.sort(colVec, new net.maizegenetics.util.StringNumberComparator());
        int extraCol = (rowLabels) ? 1 : 0;
        columnHeader = new String[colVec.size() + extraCol];
        rowHeader = new String[rowVec.size()];
        matrix = new double[rowHeader.length][columnHeader.length];
        modelMatrix = new String[rowHeader.length][columnHeader.length];
//          newData=new Object[rowVec.size()][colVec.size()];
        for (int i = 0; i < rowVec.size(); i++) {
            for (int j = 0; j < colVec.size() + extraCol; j++) {
                modelMatrix[i][j] = "0.0";
            }
        }
        for (int i = 0; i < rowVec.size(); i++) {
            rowHeader[i] = (String) rowVec.get(i);
        }
        if (rowLabels) {
            columnHeader[0] = (String) tableHeader[selectedZ];
            for (int i = 0; i < rowVec.size(); i++) {
                modelMatrix[i][0] = (String) rowVec.get(i);
            }
        }
        for (int j = 0; j < colVec.size(); j++) {
            columnHeader[j + extraCol] = (String) colVec.get(j);
        }
        for (int i = 0; i < data.length; i++) {
            double d;
            if (data[i][selectedZ].toString().equals("NaN")) {
                d = Double.NaN;
            } else {
                d = Double.valueOf(data[i][selectedZ].toString()).doubleValue();
            }
            matrix[rowVec.indexOf(data[i][selectedRow])][colVec.indexOf(data[i][selectedColumn]) + extraCol] = d;
            modelMatrix[rowVec.indexOf(data[i][selectedRow])][colVec.indexOf(data[i][selectedColumn]) + extraCol] = data[i][selectedZ].toString();
        } //end of rows

        min = (float) ZDoubleUtil.min(matrix);
        max = (float) ZDoubleUtil.max(matrix);
        range = max - min;

        S_cutoff.setValue(5000);
        LP_chart.getViewport().remove(chart);

        chart.setModel(new ZTableModel(modelMatrix, columnHeader));
        TableColumnModel cm = chart.getColumnModel();
        if (colLabels) {
            JTableHeader th = new JTableHeader(cm);
            chart.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);
            chart.setTableHeader(th);
        } else {
            chart.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
            chart.setTableHeader(null);
        }


        for (int i = 0; i < columnHeader.length; i++) {
            TableColumn c = cm.getColumn(i);

            c.setCellRenderer(rdr);
            if (colLabels) {
                c.setResizable(true);
            } else {
                c.setResizable(false);
            }
            c.setMinWidth(1);
        }
        LP_chart.getViewport().add(chart);
        CB_cellSize.setSelectedIndex(5);
        setStatusBar();
        return true;
    }

    private void setCellSize(int cSize) {
        curCellSize = cSize;

        redraw();
    }

    private void redraw() {
        try {
            chart.setRowHeight(curCellSize);

            cutoff =
                    Float.valueOf(TF_cutoff.getText()).floatValue();

            TableColumnModel cm = chart.getColumnModel();

            for (int i = 0; i < columnHeader.length; i++) {
                TableColumn c = cm.getColumn(i);

                c.setPreferredWidth(curCellSize);
                if (rowLabels && (i == 0)) {
                    c.setPreferredWidth(80);  //ed add on
                }
            }

            chart.repaint();
        } catch (NumberFormatException ex) {
        }
    }

    private void setStatusBar() {
        float minimum = (float) ZDoubleUtil.min(matrix);
        float maximum = (float) ZDoubleUtil.max(matrix);
        float mean = (float) ZDoubleUtil.mean(matrix);
        float sd = (float) ZDoubleUtil.sd(matrix);

        statusBar.setText("Statistics - Min: " + minimum + "  Max: "
                + maximum + "  Mean: " + mean + " SD: " + sd);
    }

    /**
     * Override
     */
    private void chooseBackground() {
        Color c = cc.showDialog(this, "Choose backgroud",
                chart.getBackground());

        if (c == null) {
            return;
        }

        if (!c.equals(chart.getBackground())) {
            bgColor = c;

            redraw();
        }
    }

    private void chooseForeground() {
        Color c = cc.showDialog(this, "Choose foregroud",
                chart.getForeground());

        if (c == null) {
            return;
        }

        if (!c.equals(chart.getForeground())) {
            fgColor = c;

            redraw();
        }
    }

    /**
     * Show and center the dialog
     */
    public void showAtCenter() {
        Dimension screenSize =
                Toolkit.getDefaultToolkit().getScreenSize();
        Dimension thisSize = this.getSize();

        if (thisSize.height > screenSize.height) {
            thisSize.height = screenSize.height;
        }

        if (thisSize.width > screenSize.width) {
            thisSize.width = screenSize.width;
        }

        this.setLocation((screenSize.width - thisSize.width) / 2,
                (screenSize.height - thisSize.height) / 2);
        this.setVisible(true);
        redraw();
    }

    private void jbInit() throws Exception {
        setTitle("2-D chart");

        cellBorder = BorderFactory.createEmptyBorder(15, 15, 15, 15);

        panel.setLayout(borderLayout1);
        statusBar.setBorder(BorderFactory.createLoweredBevelBorder());
        statusBar.setToolTipText("");
        statusBar.setText("    Statistics - ");
        P_toolBar.setLayout(borderLayout2);
        toolBar.setAlignmentX((float) 0.0);
        cutoffBar.setAlignmentX((float) 0.0);
        cutoffBar.setLayout(gridBagLayout1);
        L_cutoff.setText(" Cutoff  ");
        L_cellSize.setText(" Cell size  ");
        chart.setBackground(Color.white);
        chart.setToolTipText("");
//               chart.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        chart.setIntercellSpacing(new Dimension(0, 0));
//             chart.setTableHeader(null);
        TF_cutoff.setMaximumSize(new Dimension(80, 21));
        TF_cutoff.setMinimumSize(new Dimension(80, 21));
        TF_cutoff.setPreferredSize(new Dimension(80, 21));
        TF_cutoff.setToolTipText("");
        TF_cutoff.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                TF_cutoff_actionPerformed(e);
            }
        });
        S_cutoff.setMinorTickSpacing(2000);
        S_cutoff.setMajorTickSpacing(10000);
        S_cutoff.setPaintTicks(true);
        S_cutoff.setMaximum(100000);
        S_cutoff.setToolTipText("Drag the slide to set cutoff!");
        S_cutoff.addChangeListener(new javax.swing.event.ChangeListener() {

            public void stateChanged(ChangeEvent e) {
                S_cutoff_stateChanged(e);
            }
        });
        C_pValueData.setText("P-Value");
        C_pValueData.setToolTipText("");
        C_pValueData.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                C_pValueData_actionPerformed(e);
            }
        });
        C_upperTriangle.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                C_upperTriangle_actionPerformed(e);
            }
        });
        C_upperTriangle.setToolTipText("");
        C_upperTriangle.setText("only upper triangle");
        jLabel1.setText("Column:");
        jLabel2.setText("Row:");
        jLabel3.setText("Value:");
        rowComboBox.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                rowComboBox_actionPerformed(e);
            }
        });
        columnComboBox.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                columnComboBox_actionPerformed(e);
            }
        });
        valueComboBox.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                valueComboBox_actionPerformed(e);
            }
        });
        rowLabelCheckBox.setToolTipText("Show Row Labels");
        rowLabelCheckBox.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                rowLabelCheckBox_actionPerformed(e);
            }
        });
        colLabelCheckBox.setToolTipText("Show Column Labels");
        colLabelCheckBox.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                colLabelCheckBox_actionPerformed(e);
            }
        });
        svgButton.setMargin(new Insets(2, 2, 2, 2));
        svgButton.setText("Save");
        svgButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                svgButton_actionPerformed(e);
            }
        });
        this.getContentPane().add(panel, BorderLayout.CENTER);
        panel.add(LP_chart, BorderLayout.CENTER);
        LP_chart.getViewport().add(chart);
        panel.add(statusBar, BorderLayout.SOUTH);
        panel.add(P_toolBar, BorderLayout.NORTH);
        P_toolBar.add(toolBar, BorderLayout.NORTH);
        P_toolBar.add(cutoffBar, BorderLayout.SOUTH);
        toolBar.addButton("exit", "", "exit");
        toolBar.getButton("exit").addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                dispose();
            }
        });
        toolBar.addButton("foreground", "", "choose foreground");
        toolBar.getButton("foreground").addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                chooseForeground();
            }
        });
        toolBar.addButton("background", "", "choose background");
        toolBar.getButton("background").addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                chooseBackground();
            }
        });
        toolBar.addSep();
        toolBar.add(L_cellSize, null);

        for (int i = 2; i < 20; i++) {
            CB_cellSize.addItem(new Integer(i * 2));
        }

        CB_cellSize.setPreferredSize(new Dimension(40, 21));
        CB_cellSize.setMinimumSize(new Dimension(40, 21));
        CB_cellSize.setMaximumSize(new Dimension(40, 21));
        CB_cellSize.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                Object obj = CB_cellSize.getSelectedItem();

                if (obj == null) {
                    return;
                }

                Integer iobj = (Integer) obj;

                setCellSize(iobj.intValue());
            }
        });
        toolBar.addSep();
        toolBar.add(CB_cellSize);
//        toolBar.add(C_upperTriangle, null);
//        toolBar.addSep();
//        C_pValueData.setSize(50,20);
//        toolBar.add(C_pValueData, null);
//        String[] s=AbstractDisplayPlugin.getPossibleGraphicOutFormats();
//        formatComboBox=new JComboBox(s);
        rightToolPanel.add(C_upperTriangle);
        rightToolPanel.add(C_pValueData);
        rightToolPanel.add(svgButton);
//        rightToolPanel.add(formatComboBox);
        toolBar.addButton("help", "", "cutoff for pValue data");
        toolBar.getButton("help").addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                showHelp();
            }
        });
        toolBar.add(rightToolPanel);
        toolBar.getButton("help").setEnabled(false);
        cutoffBar.add(L_cutoff, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(9, 0, 8, 0), 0, 0));
        cutoffBar.add(TF_cutoff, new GridBagConstraints(1, 0, 1, 1, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(6, 0, 0, 0), 0, 0));
        cutoffBar.add(S_cutoff, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 20, 0));
//        cutoffBar.add(svgButton,  new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
//                ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        P_toolBar.add(dataPanel, BorderLayout.CENTER);
        dataPanel.add(jLabel2, null);
        dataPanel.add(rowLabelCheckBox, null);
        rowLabelCheckBox.setSelected(true);
        dataPanel.add(rowComboBox, null);
        dataPanel.add(jLabel1, null);
        dataPanel.add(colLabelCheckBox, null);
        colLabelCheckBox.setSelected(true);
        dataPanel.add(columnComboBox, null);
        dataPanel.add(jLabel3, null);
        dataPanel.add(valueComboBox, null);
        this.C_pValueData.setSelected(true);
        C_pValueData_actionPerformed(null);
    }

    private void S_cutoff_stateChanged(ChangeEvent e) {
        redraw();

        float currentValue = S_cutoff.getValue();
        String str = String.valueOf(min
                + currentValue * range / 100000);

        TF_cutoff.setText(str);
    }

    private void TF_cutoff_actionPerformed(ActionEvent e) {
        String currentValue = TF_cutoff.getText();

        try {
            float f = Float.valueOf(currentValue).floatValue();

            if (f < min) {
                f = min;
            }

            if (f > max) {
                f = max;
            }

            S_cutoff.setValue((int) ((f - min) * 100000 / range));
        } catch (NumberFormatException ex) {
            float c = S_cutoff.getValue();

            TF_cutoff.setText(String.valueOf(min
                    + c * range / 100000));
        }
    }

    private void C_pValueData_actionPerformed(ActionEvent e) {
        if (C_pValueData.isSelected() == true) {
            S_cutoff.setEnabled(false);
            TF_cutoff.setEnabled(false);
            toolBar.getButton("help").setEnabled(true);
        } else {
            S_cutoff.setEnabled(true);
            TF_cutoff.setEnabled(true);
            toolBar.getButton("help").setEnabled(false);
        }

        redraw();
    }

    private void C_upperTriangle_actionPerformed(ActionEvent e) {
        redraw();
    }

    private void showHelp() {
        ZDialog d = new ZDialog();

        d.getContentPane().add(pHelp);
        d.pack();
        d.showAtCenter();
    }

    void setRowComboBox(String s) {
        rowComboBox.setSelectedItem(s);
        rowComboBox_actionPerformed(null);
    }

    void rowComboBox_actionPerformed(ActionEvent e) {
        if (rowComboBox.getSelectedIndex() == -1) {
            return;
        }
        selectedRow = headerVector.indexOf(rowComboBox.getSelectedItem());
        processAndSetData(tableHeader, tableData);
    }

    void setColumnComboBox(String s) {
        columnComboBox.setSelectedItem(s);
        columnComboBox_actionPerformed(null);
    }

    void columnComboBox_actionPerformed(ActionEvent e) {
        if (columnComboBox.getSelectedIndex() == -1) {
            return;
        }
        selectedColumn = headerVector.indexOf(columnComboBox.getSelectedItem());
        processAndSetData(tableHeader, tableData);
    }

    void setValueComboBox(String s) {
        valueComboBox.setSelectedItem(s);
        valueComboBox_actionPerformed(null);
    }

    void valueComboBox_actionPerformed(ActionEvent e) {
        if (valueComboBox.getSelectedIndex() == -1) {
            return;
        }
        selectedZ = headerVector.indexOf(valueComboBox.getSelectedItem());
        processAndSetData(tableHeader, tableData);
    }

    void rowLabelCheckBox_actionPerformed(ActionEvent e) {
        rowLabels = rowLabelCheckBox.isSelected();
        processAndSetData(tableHeader, tableData);
    }

    void colLabelCheckBox_actionPerformed(ActionEvent e) {
        colLabels = colLabelCheckBox.isSelected();
        processAndSetData(tableHeader, tableData);
    }

    void svgButton_actionPerformed(ActionEvent e) {
        //     LP_chart.paint(svgGenerator);
//      String s=formatComboBox.getSelectedItem().toString();
        theGrid2dDisplayPlugin.saveDataToFile(LP_chart);
        // ZJIMIUtil.saveAsJPG(ZImageUtil.getComponentImage(chart),fileName);
    }

    Component getChart() {
        redraw();
        return LP_chart;
    }
}
