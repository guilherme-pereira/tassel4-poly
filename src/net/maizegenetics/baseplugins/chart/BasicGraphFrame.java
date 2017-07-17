package net.maizegenetics.baseplugins.chart;

import net.maizegenetics.pal.report.TableReport;

import org.jfree.chart.editor.ChartEditor;
import org.jfree.chart.editor.ChartEditorManager;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.WindowEvent;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2004</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */
public class BasicGraphFrame extends JFrame {
    // String[] graphStrings = { "Histogram", "XY Plot", "Bar Chart", "Pie Chart"};

    ChartDisplayPlugin theChartDisplayPlugin;
    BorderLayout borderLayout1 = new BorderLayout();
    TableReport theTable;
    JComboBox graphTypeComboBox = new JComboBox(ChartDisplayPlugin.getPossibleCharts());
    JLabel jLabel1 = new JLabel();
    JPanel headerPanel = new JPanel();
    JButton saveButton = new JButton();
    JButton propertyButton = new JButton();
    BasicChartPanel chartPanel;
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    JFileChooser theFileChooser;

    public BasicGraphFrame(ChartDisplayPlugin cdp, TableReport theTable) {
        this.theTable = theTable;
        theChartDisplayPlugin = cdp;
        try {
            jbInit();
            chartPanel = new HistogramPanel(theTable);
            this.getContentPane().add(chartPanel, BorderLayout.CENTER);
            this.pack();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    void jbInit() throws Exception {
        headerPanel.setMinimumSize(new Dimension(396, 29));
        headerPanel.setLayout(gridBagLayout1);
        propertyButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                propertyButton_actionPerformed(e);
            }
        });
        graphTypeComboBox.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                graphTypeComboBox_actionPerformed(e);
            }
        });
        this.addWindowListener(new java.awt.event.WindowAdapter() {

            public void windowClosed(WindowEvent e) {
                this_windowClosed(e);
            }
        });
        saveButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                saveButton_actionPerformed(e);
            }
        });
        headerPanel.add(jLabel1, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 144), 26, 11));
        headerPanel.add(graphTypeComboBox, new GridBagConstraints(0, 0, 1, 1, 1.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 105, 0, 0), 132, 0));
        headerPanel.add(saveButton, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 7, 0, 0), 2, 0));
        headerPanel.add(propertyButton, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 4), 2, 1));
        this.getContentPane().setLayout(borderLayout1);
        jLabel1.setText("Graph Type:");
        saveButton.setText("Save");
        propertyButton.setText("Properties");
        this.getContentPane().add(headerPanel, BorderLayout.NORTH);

    }

    void graphTypeComboBox_actionPerformed(ActionEvent e) {
        ChartDisplayPlugin.ChartType ct = ChartDisplayPlugin.ChartType.valueOf((String) graphTypeComboBox.getSelectedItem());
        this.getContentPane().remove(chartPanel);
        switch (ct) {
            case Histogram:
                chartPanel = new HistogramPanel(theTable);
                break;
            case XYScatter:
                chartPanel = new XYScatterPanel(theTable);
                break;
            case BarChart:
                chartPanel = new BarChartPanel(theTable);
                break;
            case PieChart:
                chartPanel = new PieChartPanel(theTable);
                break;
        }
        this.getContentPane().add(chartPanel, BorderLayout.CENTER);
        pack();
    }

    void propertyButton_actionPerformed(ActionEvent e) {
        ChartEditor panel = ChartEditorManager.getChartEditor(chartPanel.getChart());
        int result = JOptionPane.showConfirmDialog(this, panel, "Chart Properties",
                JOptionPane.OK_CANCEL_OPTION,
                JOptionPane.PLAIN_MESSAGE);
        if (result == JOptionPane.OK_OPTION) {
            panel.updateChart(chartPanel.getChart());
            //     panel.updateChartProperties(chartPanel.getChart());
        }

    }

    void this_windowClosed(WindowEvent e) {
        System.exit(0);
    }

    void saveButton_actionPerformed(ActionEvent e) {
        //JFreeChart aChart=chartPanel.chart;
        theChartDisplayPlugin.saveDataToFile(chartPanel.getMainComponent());
    }
}
