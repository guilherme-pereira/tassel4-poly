package net.maizegenetics.baseplugins.chart;

import net.maizegenetics.pal.report.TableReport;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StandardXYItemRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.DefaultXYItemRenderer;
import org.jfree.data.function.Function2D;
import org.jfree.data.function.LineFunction2D;
import org.jfree.data.general.DatasetUtilities;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.statistics.Regression;

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


public class XYScatterPanel extends BasicChartPanel {
  BorderLayout borderLayout1 = new BorderLayout();
  TableReportXYDataset dataset, dataset2;
//  JFreeChart chart;
  ChartPanel chartPanel;
  TableReport theTable;
  JPanel controlPanel = new JPanel();
  JComboBox seriesXComboBox;// = new JComboBox();
  JComboBox seriesY1ComboBox;// = new JComboBox();
  JLabel jLabel1 = new JLabel();
  JLabel jLabel2 = new JLabel();

  String[] columnNames;
  JCheckBox lineChartCheckBox = new JCheckBox();
  JComboBox seriesY2ComboBox;
  JLabel jLabel3 = new JLabel();
  JCheckBox regressionCheckBox = new JCheckBox();
  JCheckBox multiYAxisCheckBox = new JCheckBox();
  GridBagLayout gridBagLayout1 = new GridBagLayout();

  public XYScatterPanel(TableReport theTable) {
    this.theTable=theTable;
    try {
      Object[] colNames = theTable.getTableColumnNames();
      columnNames = new String[colNames.length+1];
      columnNames[0] = "None";
      for (int i = 1; i < columnNames.length; i++) {
        columnNames[i] = (String) colNames[i - 1];
      }

      dataset = null;
      chart = createChart(dataset,null,false,regressionCheckBox.isSelected());
      chartPanel = new ChartPanel(chart);
      chartPanel.setPreferredSize(new java.awt.Dimension(500, 270));
      chartPanel.setMouseZoomable(true, false);
      jbInit();
    }
    catch(Exception ex) {
      ex.printStackTrace();
    }
  }

  public XYScatterPanel(TableReport theTable, int seriesX, int seriesY1, int seriesY2, boolean isRegression) {
    //this constructor is for the non-interactive version
    this.theTable=theTable;
    try {
        dataset = createDataset(seriesX, seriesY1, seriesY2);
        dataset2=null;
        chart = createChart(dataset, dataset2, false, isRegression);
        chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new java.awt.Dimension(500, 270));
        this.setVisible(true);
    }
    catch(Exception ex) {
      ex.printStackTrace();
    }
  }

  void jbInit() throws Exception {
    seriesY2ComboBox = new JComboBox(columnNames);
    this.setLayout(borderLayout1);
    controlPanel.setLayout(gridBagLayout1);
    jLabel1.setText("Y1");
    jLabel2.setText("X");
    seriesXComboBox = new JComboBox(columnNames);
    seriesY1ComboBox = new JComboBox(columnNames);
    seriesXComboBox.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(ActionEvent e) {
        seriesXComboBox_actionPerformed(e);
      }
    });
    seriesY1ComboBox.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(ActionEvent e) {
        seriesY1ComboBox_actionPerformed(e);
      }
    });
    lineChartCheckBox.setToolTipText("Connect scatter with a line");
    lineChartCheckBox.setText("Line");
    lineChartCheckBox.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(ActionEvent e) {
        lineChartCheckBox_actionPerformed(e);
      }
    });
    controlPanel.setMinimumSize(new Dimension(393, 50));
    controlPanel.setPreferredSize(new Dimension(393, 50));
    seriesY2ComboBox.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(ActionEvent e) {
        seriesY2ComboBox_actionPerformed(e);
      }
    });
    jLabel3.setText("Y2");
    regressionCheckBox.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(ActionEvent e) {
        regressionCheckBox_actionPerformed(e);
      }
    });
    regressionCheckBox.setToolTipText("Plot a linear regression line");
    regressionCheckBox.setText("Regression");
    seriesXComboBox.setToolTipText("Numerical variable for x-axis");
    seriesY1ComboBox.setToolTipText("Numerical variable for first Y");
    multiYAxisCheckBox.setText("2 Y Axes");
    multiYAxisCheckBox.setToolTipText("Plot Y2 on a second axis scaling");
    multiYAxisCheckBox.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(ActionEvent e) {
        multiYAxisCheckBox_actionPerformed(e);
      }
    });
    seriesY2ComboBox.setToolTipText("Numerical variable for optional second Y");
    this.add(chartPanel, BorderLayout.CENTER);
    this.add(controlPanel, BorderLayout.NORTH);
    controlPanel.add(jLabel2,  new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(7, 12, 1, 111), 26, -4));
    controlPanel.add(lineChartCheckBox,  new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(6, 20, 1, 0), 0, -11));
    controlPanel.add(seriesXComboBox,  new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 24, 1, 0), 104, 0));
    controlPanel.add(regressionCheckBox,  new GridBagConstraints(2, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(6, 0, 1, 0), 2, -11));
    controlPanel.add(multiYAxisCheckBox,  new GridBagConstraints(3, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(6, 0, 1, 0), 4, -11));
    controlPanel.add(seriesY1ComboBox,  new GridBagConstraints(0, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(2, 25, 0, 0), 104, 0));
    controlPanel.add(jLabel1,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 6, 0, 99), 35, -4));
    controlPanel.add(jLabel3,  new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 11, 0, 34), 12, 0));
    controlPanel.add(seriesY2ComboBox,  new GridBagConstraints(1, 0, 2, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(2, 29, 0, 0), 104, 0));
  }


  /**
 * Creates a sample {@link HistogramDataset}.
 *
 * @return The dataset.
 */
TableReportXYDataset createDataset(int seriesXColumn, int seriesY1Column, int seriesY2Column) {
  TableReportXYDataset dataset;
  if (seriesY2Column < 0) {
    dataset = new TableReportXYDataset(theTable, seriesXColumn, seriesY1Column);
  }
  else
    {dataset = new TableReportXYDataset(theTable, seriesXColumn, seriesY1Column, seriesY2Column);}
    return dataset;
}


  /**
   * Creates a chart.
   *
   * @param dataset  a dataset.
   *
   * @return The chart.
   */
  JFreeChart createChart(TableReportXYDataset dataset, TableReportXYDataset dataset2, boolean lineChart, boolean isRegression) {
    String name="Please select numeric variables";
    String xName="X";
    String y1Name="Y";
    String y2Name="Y2";
    if(dataset!=null) {
      xName=dataset.getXName();
      y1Name=dataset.getSeriesName(0);
      name=xName+" vs. "+y1Name;
      if(dataset.getSeriesCount()==2)
        {y2Name=dataset.getSeriesName(1);
          name+=" and "+y2Name;}
      else if(dataset2!=null)
          {y2Name=dataset2.getSeriesName(0);
            name+=" and "+y2Name;}

    }
    if(lineChart) {
      chart = ChartFactory.createXYLineChart(
          name,
          xName, y1Name,
          dataset,
          PlotOrientation.VERTICAL,
          true,
          true,
          false
          );

    }
    else {
      chart = ChartFactory.createScatterPlot(
        name,
        xName,y1Name,
        dataset,
        PlotOrientation.VERTICAL,
        true,
        true,
        false
        );
      }
    chart.getXYPlot().setForegroundAlpha(0.75f);
    if(dataset2!=null) {
      create2ndAxis(chart.getXYPlot(),dataset2,0,Color.BLUE, lineChart);
    }
    if(isRegression)
      {createRegression(chart.getXYPlot(), dataset, 0, Color.RED);
      if(dataset.getSeriesCount()==2) createRegression(chart.getXYPlot(), dataset, 1, Color.BLUE);
       if(dataset2!=null) createRegression(chart.getXYPlot(),  dataset2, 0, Color.BLUE);
      }
    return chart;
  }
  private void create2ndAxis(XYPlot plot,  XYDataset data, int series,Color theColor, boolean isLine) {
      //final NumberAxis axis2 = new NumberAxis(data.getSeriesName(series));
    final NumberAxis axis2 = new NumberAxis((String)data.getSeriesKey(series));
    axis2.setAutoRangeIncludesZero(false);
    plot.setRangeAxis(1, axis2);
    plot.setDataset(1, data);
    plot.mapDatasetToRangeAxis(1, 1);
    StandardXYItemRenderer renderer2; //= new StandardXYItemRenderer();
 //   renderer2.setSeriesPaint(0, theColor);
//    if(isLine) {renderer2.setPlotShapes(false); renderer2.setPlotLines(true);}
//    else {renderer2.setPlotShapes(true); renderer2.setPlotLines(false);}
    if(isLine) {renderer2 = new StandardXYItemRenderer(StandardXYItemRenderer.LINES);}
    else {renderer2 = new StandardXYItemRenderer(StandardXYItemRenderer.SHAPES);}
        renderer2.setSeriesPaint(0, theColor);
    plot.setRenderer(1, renderer2);
  }
  private void createRegression(XYPlot plot,  XYDataset data, int series, Color theColor) {
    // calculate the regression and create subplot 2...
     double[] coefficients = Regression.getOLSRegression(data, series);
     double max=DatasetUtilities.findMaximumDomainValue(data).doubleValue();
     double min=DatasetUtilities.findMinimumDomainValue(data).doubleValue();
     Function2D curve = new LineFunction2D(coefficients[0], coefficients[1]);
     XYDataset regressionData = DatasetUtilities.sampleFunction2D(
         curve,
         min,max, 100,
         data.getSeriesKey(series)+" Fitted Reg. Line"
     );
     int datasetCount=plot.getDatasetCount();
     plot.setDataset(datasetCount, regressionData);
     XYItemRenderer renderer2 = new StandardXYItemRenderer(StandardXYItemRenderer.LINES);
     renderer2.setSeriesPaint(0, theColor);
     plot.setRenderer(datasetCount, renderer2);

  }

  void updateChart() {
    int seriesX, seriesY1, seriesY2;
    boolean lineChart;
    seriesX = seriesXComboBox.getSelectedIndex() - 1;
    seriesY1 = seriesY1ComboBox.getSelectedIndex() - 1;
    seriesY2=seriesY2ComboBox.getSelectedIndex() - 1;
    lineChart=lineChartCheckBox.isSelected();
    if((seriesX>-1)&(seriesY1>-1)){
      if(multiYAxisCheckBox.isSelected()) {
        dataset = createDataset(seriesX, seriesY1, -1);
        dataset2 = createDataset(seriesX, seriesY2, -1);
      } else
      {dataset = createDataset(seriesX, seriesY1, seriesY2); dataset2=null;}
      chart = createChart(dataset, dataset2,lineChart, regressionCheckBox.isSelected());
//      chart.getXYPlot().getDomainAxis().setLabel((String)seriesXComboBox.getSelectedItem());
//      chart.getXYPlot().getRangeAxis().setLabel((String)seriesY1ComboBox.getSelectedItem());
      chartPanel.setChart(chart);
    }
  }

  void seriesXComboBox_actionPerformed(ActionEvent e) {
    updateChart();
  }

  void seriesY1ComboBox_actionPerformed(ActionEvent e) {
    updateChart();
  }

  void lineChartCheckBox_actionPerformed(ActionEvent e) {
    updateChart();
  }
  void seriesY2ComboBox_actionPerformed(ActionEvent e) {
    updateChart();
  }
  void regressionCheckBox_actionPerformed(ActionEvent e) {
    updateChart();
  }
  void multiYAxisCheckBox_actionPerformed(ActionEvent e) {
     updateChart();
  }

  public JComponent getMainComponent() {
        return chartPanel;
    }

}

