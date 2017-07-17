package net.maizegenetics.baseplugins.chart;

import net.maizegenetics.pal.report.TableReport;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.renderer.category.BoxAndWhiskerRenderer;
import org.jfree.chart.renderer.category.CategoryItemRenderer;
import org.jfree.chart.renderer.category.StatisticalBarRenderer;
import org.jfree.data.category.CategoryDataset;
import org.jfree.data.statistics.HistogramDataset;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;





/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2004</p>
 * <p>Company: </p>
 * @author Ed Buckler
 * @version 1.0
 */

/**
 * todo Need to make conconcated alignment report sequence haplotypes into table
 */

public class BarChartPanel extends BasicChartPanel {
  boolean isBoxWhisker=false, isBar=false, errBar=true;
  BorderLayout borderLayout1 = new BorderLayout();
//  StatisticalCategoryDataset dataset;
  CategoryDataset dataset;
//  JFreeChart chart;
  ChartPanel chartPanel;
  TableReport theTable;
  JPanel controlPanel = new JPanel();
  JComboBox series1ComboBox;// = new JComboBox();
  JComboBox series2ComboBox;// = new JComboBox();
  JLabel jLabel1 = new JLabel();
  JLabel jLabel2 = new JLabel();

  String[] columnNames;
  int bins=5;
  JComboBox categoryComboBox;
  JLabel jLabel3 = new JLabel();
  ButtonGroup errorButtonGroup = new ButtonGroup();
  JRadioButton errorRadioButton = new JRadioButton();
  JRadioButton deviationRadioButton = new JRadioButton();
  JRadioButton noErrRadioButton = new JRadioButton();
  JRadioButton boxWhiskerRadioButton = new JRadioButton();
  Component component1;
  GridBagLayout gridBagLayout1 = new GridBagLayout();

  public BarChartPanel(TableReport theTable) {
    this.theTable=theTable;
    try {
      Object[] colNames = theTable.getTableColumnNames();
      columnNames = new String[colNames.length+1];
      columnNames[0] = "None";
      for (int i = 1; i < columnNames.length; i++) {
        columnNames[i] = (String) colNames[i - 1];
      }

      dataset = null;
      chart = createChart(dataset,false, false);
      chartPanel = new ChartPanel(chart);
      chartPanel.setPreferredSize(new java.awt.Dimension(500, 270));
      chartPanel.setMouseZoomable(true, false);
      jbInit();
    }
    catch(Exception ex) {
      ex.printStackTrace();
    }
  }

  public BarChartPanel(TableReport theTable, int seriesCat, int seriesY1, int seriesY2, boolean isBoxWhisker,
                       boolean isBar, boolean errBar) {
    this.theTable=theTable;
    this.isBoxWhisker=isBoxWhisker;
    this.isBar=isBar;
    this.errBar=errBar;
    dataset= this.createDataset(seriesCat, seriesY1, seriesY2, isBoxWhisker, errBar);
    chart = createChart(dataset,isBoxWhisker, isBar);
    chartPanel = new ChartPanel(chart);
    chartPanel.setPreferredSize(new java.awt.Dimension(500, 270));
  }

  void jbInit() throws Exception {
    categoryComboBox = new JComboBox(columnNames);
    component1 = Box.createHorizontalStrut(8);
    this.setLayout(borderLayout1);
    controlPanel.setLayout(gridBagLayout1);
    jLabel1.setText("Y2");
    jLabel2.setText("Category");
    series1ComboBox = new JComboBox(columnNames);
    series2ComboBox = new JComboBox(columnNames);
    series1ComboBox.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(ActionEvent e) {
        series1ComboBox_actionPerformed(e);
      }
    });
    series2ComboBox.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(ActionEvent e) {
        series2ComboBox_actionPerformed(e);
      }
    });
    controlPanel.setMinimumSize(new Dimension(394, 60));
    controlPanel.setPreferredSize(new Dimension(394, 60));
    categoryComboBox.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(ActionEvent e) {
        categoryComboBox_actionPerformed(e);
      }
    });
    jLabel3.setText("Y1");
    errorRadioButton.setToolTipText("Standard error of the mean");
    errorRadioButton.setText("Std Err");
    errorRadioButton.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(ActionEvent e) {
        errorRadioButton_actionPerformed(e);
      }
    });
    deviationRadioButton.setOpaque(true);
    deviationRadioButton.setToolTipText("Stardard deviation");
    deviationRadioButton.setText("Std Dev");
    deviationRadioButton.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(ActionEvent e) {
        deviationRadioButton_actionPerformed(e);
      }
    });
    noErrRadioButton.setSelected(true);
    noErrRadioButton.setText("None");
    noErrRadioButton.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(ActionEvent e) {
        noErrRadioButton_actionPerformed(e);
      }
    });
    boxWhiskerRadioButton.setToolTipText("Box and Whisker Plot");
    boxWhiskerRadioButton.setText("Box");
    boxWhiskerRadioButton.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(ActionEvent e) {
        boxWhiskerRadioButton_actionPerformed(e);
      }
    });
    this.add(chartPanel, BorderLayout.CENTER);
    this.add(controlPanel, BorderLayout.NORTH);
    errorButtonGroup.add(errorRadioButton);
    errorButtonGroup.add(deviationRadioButton);
    controlPanel.add(series1ComboBox,  new GridBagConstraints(0, 0, 2, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 30, 0, 12), 101, 3));
    controlPanel.add(series2ComboBox,  new GridBagConstraints(2, 0, 2, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 29, 0, 0), 90, 3));
    controlPanel.add(jLabel1,  new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 42), 14, 7));
    controlPanel.add(jLabel3,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 10, 0, 19), 15, 7));
    controlPanel.add(boxWhiskerRadioButton,  new GridBagConstraints(4, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 15, 0, 2), 9, -2));
    controlPanel.add(categoryComboBox,  new GridBagConstraints(1, 1, 1, 2, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 3, 0), 83, 3));
    controlPanel.add(jLabel2,  new GridBagConstraints(0, 1, 1, 2, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 3, 0), 8, 7));
    controlPanel.add(noErrRadioButton,  new GridBagConstraints(2, 1, 1, 2, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 13, 3, 0), 5, -2));
    controlPanel.add(errorRadioButton,  new GridBagConstraints(3, 1, 1, 2, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 3, 0), 1, -2));
    controlPanel.add(deviationRadioButton,  new GridBagConstraints(4, 1, 1, 2, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 3, 2), 6, -2));
    controlPanel.add(component1,  new GridBagConstraints(0, 1, 5, GridBagConstraints.REMAINDER, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 2), 382, 2));
    errorButtonGroup.add(noErrRadioButton);
    errorButtonGroup.add(boxWhiskerRadioButton);
  }


  /**
 * Creates a sample {@link HistogramDataset}.
 *
 * @return The dataset.
 */
CategoryDataset createDataset(int categoryColumn, int series1Column,
   int series2Column, boolean boxwhisker, boolean stderr) {
    CategoryDataset dataset;
    int[] seriesY;
    if(series2Column>-1) {seriesY=new int[2]; seriesY[0]=series1Column; seriesY[1]=series2Column;}
    else {seriesY=new int[1]; seriesY[0]=series1Column;}
    if(boxwhisker)
    {dataset=new TableReportBoxWhiskerCatDataset(theTable, categoryColumn, seriesY);}
    else {dataset=new TableReportStatCategoryDataset(theTable, categoryColumn, seriesY, stderr);}
//    final BasicXYDataset dataset = new BasicXYDataset(theTable,series1Column, series2Column);
    return dataset;
}


  /**
   * Creates a chart.
   *
   * @param dataset  a dataset.
   *
   * @return The chart.
   */
  JFreeChart createChart(CategoryDataset dataset, boolean isBoxWhisker, boolean isBar) {
    String chartName;
    String categoryName;
    try{
      System.out.println("RowKey"+dataset.getRowKey(0));
      System.out.println("ColKey"+dataset.getColumnKey(0));
      chartName="Mean "+dataset.getRowKey(0);
      if(dataset.getRowCount()>1) {chartName+=" and "+dataset.getRowKey(1);}
      categoryName=(String)categoryComboBox.getSelectedItem();
      chartName+=" by "+categoryName;
    }
    catch(Exception ex) {
      System.out.println("Bar chart labels not ready");
      chartName="Please choose a categorical variable and a numeric variable";
      categoryName="Unknown";
    }
    if (isBoxWhisker) {
      final CategoryAxis xAxis = new CategoryAxis(categoryName);
      final NumberAxis yAxis = new NumberAxis("Mean");
      yAxis.setAutoRangeIncludesZero(false);
      final BoxAndWhiskerRenderer renderer = new BoxAndWhiskerRenderer();
      renderer.setFillBox(false);
      final CategoryPlot plot = new CategoryPlot((TableReportBoxWhiskerCatDataset)dataset, xAxis, yAxis,
                                                 renderer);
      chart = new JFreeChart(chartName, new Font("Helvetica", Font.BOLD, 14),
                                  plot, true);

    }
    else if(isBar) {
      final CategoryAxis xAxis = new CategoryAxis(categoryName);
      final ValueAxis yAxis = new NumberAxis("Mean");
      final CategoryItemRenderer renderer = new StatisticalBarRenderer();
      final CategoryPlot plot = new CategoryPlot((TableReportStatCategoryDataset)dataset, xAxis, yAxis, renderer);
      chart = new JFreeChart(chartName, new Font("Helvetica", Font.BOLD, 14),
                                  plot, true);
    }
    else {
      chart = ChartFactory.createBarChart(chartName,categoryName,"Mean", dataset, PlotOrientation.VERTICAL,
        true, true, false);
      }
    return chart;
  }

  void updateChart() {
    int series1, series2, categorySeries;
    series1 = series1ComboBox.getSelectedIndex() - 1;
    series2 = series2ComboBox.getSelectedIndex() - 1;
    categorySeries = categoryComboBox.getSelectedIndex()-1;
    if((series1>-1)&&(categorySeries>-1)){
      dataset = createDataset(categorySeries,series1, series2, isBoxWhisker, errBar);
      chart = createChart(dataset, isBoxWhisker, isBar);
      chartPanel.setChart(chart);
    }
  }

  void series1ComboBox_actionPerformed(ActionEvent e) {
    updateChart();
  }

  void series2ComboBox_actionPerformed(ActionEvent e) {
    updateChart();
  }

  void categoryComboBox_actionPerformed(ActionEvent e) {
    updateChart();
  }
  void updateBarStyle() {
    if(this.boxWhiskerRadioButton.isSelected()==true) {isBoxWhisker=true; isBar=false;}
    else if(noErrRadioButton.isSelected()==true) {isBoxWhisker=false; isBar=false; errBar=false;}
    else if(errorRadioButton.isSelected()==true) {isBoxWhisker=false; isBar=true; errBar=true;}
    else {isBoxWhisker=false; isBar=true; errBar=false;}
    updateChart();
  }

  void noErrRadioButton_actionPerformed(ActionEvent e) {
    updateBarStyle();
  }

  void errorRadioButton_actionPerformed(ActionEvent e) {
    updateBarStyle();
  }

  void deviationRadioButton_actionPerformed(ActionEvent e) {
    updateBarStyle();
  }

  void boxWhiskerRadioButton_actionPerformed(ActionEvent e) {
     updateBarStyle();
  }


    public JComponent getMainComponent() {
        return chartPanel;
    }
}

