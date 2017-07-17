package net.maizegenetics.baseplugins.chart;

import net.maizegenetics.pal.report.TableReport;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.data.general.PieDataset;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;

/**
 * <p>Title: PieChartPanel</p>
 * <p>Description: A Panel for selecting and designing a pie chart</p>
 * <p>Copyright: Copyright (c) 2004</p>
 * <p>Company: USDA-ARS</p>
 * @author Ed Buckler
 * @version 1.0
 */


public class PieChartPanel extends BasicChartPanel {
  BorderLayout borderLayout1 = new BorderLayout();
  PieDataset dataset;
  ChartPanel chartPanel;
  TableReport theTable;
  JPanel controlPanel = new JPanel();
  JLabel jLabel2 = new JLabel();

  String[] columnNames;
  int bins=5;
  JComboBox categoryComboBox;
  JCheckBox threeDCheckBox = new JCheckBox();
  GridBagLayout gridBagLayout1 = new GridBagLayout();

  public PieChartPanel(TableReport theTable) {
    this.theTable=theTable;
    try {
      Object[] colNames = theTable.getTableColumnNames();
      columnNames = new String[colNames.length+1];
      columnNames[0] = "None";
      for (int i = 1; i < columnNames.length; i++) {
        columnNames[i] = (String) colNames[i - 1];
      }

      dataset = null;
      chart = createChart(dataset,false);
      chartPanel = new ChartPanel(chart);
      chartPanel.setPreferredSize(new java.awt.Dimension(500, 270));
      chartPanel.setMouseZoomable(true, false);
      jbInit();
    }
    catch(Exception ex) {
      ex.printStackTrace();
    }
  }

  void jbInit() throws Exception {
    categoryComboBox = new JComboBox(columnNames);
    this.setLayout(borderLayout1);
    controlPanel.setLayout(gridBagLayout1);
    jLabel2.setText("Category");
    controlPanel.setMinimumSize(new Dimension(394, 50));
    controlPanel.setPreferredSize(new Dimension(394, 50));
    categoryComboBox.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(ActionEvent e) {
        categoryComboBox_actionPerformed(e);
      }
    });
    threeDCheckBox.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(ActionEvent e) {
        threeDCheckBox_actionPerformed(e);
      }
    });
    threeDCheckBox.setOpaque(true);
    threeDCheckBox.setText("3-D");
    this.add(chartPanel, BorderLayout.CENTER);
    this.add(controlPanel, BorderLayout.NORTH);
    controlPanel.add(categoryComboBox,  new GridBagConstraints(1, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(1, 0, 28, 0), 120, 0));
    controlPanel.add(jLabel2,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 1, 28, 0), 7, -1));
    controlPanel.add(threeDCheckBox,  new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(1, 59, 28, 55), 45, -7));
  }


  /**
 * Creates a sample {@link org.jfree.data.statistics.HistogramDataset}.
 *
 * @return The dataset.
 */
PieDataset createDataset(int categoryColumn) {
    TableReportPieDataset dataset=new TableReportPieDataset(theTable, categoryColumn);
    return dataset;
}


  /**
   * Creates a chart.
   *
   * @param dataset  a dataset.
   *
   * @return The chart.
   */
  JFreeChart createChart(PieDataset dataset, boolean is3D) {
    String chartName="TEST";
    String categoryName="CHOOSE CATEGORY";
    try{
      categoryName=(String)categoryComboBox.getSelectedItem();
      chartName="Frequency of "+categoryName;
    }
    catch (Exception ex) {
      System.out.println("Pie chart labels not ready");
      chartName = "Please choose a categorical variable";
      categoryName = "Unknown";
    }
    if(is3D) {
      chart = ChartFactory.createPieChart3D(chartName, dataset, true, true, false);
    }
    else {
      chart = ChartFactory.createPieChart(chartName, dataset, true, true, false);
      }
    return chart;
  }

  void updateChart() {
    int categorySeries;
    boolean is3D;
    categorySeries = categoryComboBox.getSelectedIndex()-1;
    is3D=threeDCheckBox.isSelected();
    if(categorySeries>-1){
      dataset = createDataset(categorySeries);
      chart = createChart(dataset, is3D);
      chartPanel.setChart(chart);
    }
  }

  void categoryComboBox_actionPerformed(ActionEvent e) {
    updateChart();
  }
  void threeDCheckBox_actionPerformed(ActionEvent e) {
    updateChart();
  }

  public JComponent getMainComponent() {
        return chartPanel;
    }
}

