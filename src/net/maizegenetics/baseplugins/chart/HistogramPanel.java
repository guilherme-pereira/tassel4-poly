package net.maizegenetics.baseplugins.chart;

import net.maizegenetics.pal.report.TableReport;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.IntervalXYDataset;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.statistics.HistogramType;

import javax.swing.*;
import javax.swing.event.CaretEvent;
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


//todo catch errors on bins empty issue



public class HistogramPanel extends BasicChartPanel {
  BorderLayout borderLayout1 = new BorderLayout();
  IntervalXYDataset dataset;
  ChartPanel chartPanel;
  TableReport theTable;
  JPanel controlPanel = new JPanel();
  JComboBox series1ComboBox;
  JComboBox series2ComboBox;
  JLabel jLabel1 = new JLabel();
  JLabel jLabel2 = new JLabel();
  JTextField binsTextField = new JTextField();
  JLabel jLabel3 = new JLabel();

  String[] columnNames;
  int bins=5;
  GridBagLayout gridBagLayout1 = new GridBagLayout();

  public HistogramPanel(TableReport theTable) {
    this.theTable=theTable;
    try {
      Object[] colNames = theTable.getTableColumnNames();
      columnNames = new String[colNames.length+1];
      columnNames[0] = "None";
      for (int i = 1; i < columnNames.length; i++) {
        columnNames[i] = (String) colNames[i - 1];
      }

      dataset = null;
      chart = createChart(dataset);
      chartPanel = new ChartPanel(chart);
      chartPanel.setPreferredSize(new java.awt.Dimension(500, 270));
      chartPanel.setMouseZoomable(true, false);
      jbInit();
    }
    catch(Exception ex) {
      ex.printStackTrace();
    }
  }

  public HistogramPanel(TableReport theTable, int series1, int series2, int bins) {
    //this constructor is for the non-interactive version
    this.theTable=theTable;
    try {
        this.bins=bins;
        dataset = createDataset(series1, series2);
        chart = createChart(dataset);
        chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new java.awt.Dimension(500, 270));
        this.setVisible(true);
    }
    catch(Exception ex) {
      ex.printStackTrace();
    }
  }

  void jbInit() throws Exception {
    this.setLayout(borderLayout1);
    controlPanel.setLayout(gridBagLayout1);
    jLabel1.setText("Series 2");
    jLabel2.setText("Series 1");
    series1ComboBox = new JComboBox(columnNames);
    series2ComboBox = new JComboBox(columnNames);
    binsTextField.setText(""+bins);
    binsTextField.addCaretListener(new javax.swing.event.CaretListener() {
      public void caretUpdate(CaretEvent e) {
        binsTextField_caretUpdate(e);
      }
    });
    binsTextField.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(ActionEvent e) {
        binsTextField_actionPerformed(e);
      }
    });
    jLabel3.setText("Bins");
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
    this.add(chartPanel, BorderLayout.CENTER);
    this.add(controlPanel, BorderLayout.NORTH);
    controlPanel.add(jLabel2,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 3, -4));
    controlPanel.add(binsTextField,   new GridBagConstraints(5, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 4), 23, 0));
    controlPanel.add(jLabel3,  new GridBagConstraints(4, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 6, -1));
    controlPanel.add(series1ComboBox,  new GridBagConstraints(1, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 84, 0));
    controlPanel.add(series2ComboBox,  new GridBagConstraints(3, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 80, 0));
    controlPanel.add(jLabel1,  new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 3, -4));
  }

  /**
 * Creates a sample {@link HistogramDataset}.
 *
 * @return The dataset.
 */
IntervalXYDataset createDataset(int series1Column, int series2Column) {
    final HistogramDataset dataset = new HistogramDataset();
    dataset.setType(HistogramType.FREQUENCY);
    double[] series1Data, series2Data;
    Object[] theNames = theTable.getTableColumnNames();
    if(series1Column>=0) {
      series1Data = getFilteredNumericData(series1Column);
      if(series1Data.length>0) dataset.addSeries( (String) theNames[series1Column], series1Data, bins);
    }
    if(series2Column>=0) {
      series2Data = getFilteredNumericData(series2Column);
      if(series2Data.length>0) dataset.addSeries( (String) theNames[series2Column], series2Data, bins);
    }
    return dataset;
}

  double[] getFilteredNumericData(int column) {
    Object[][] theRawData=theTable.getTableData();
    int countGood=0;
    double[] tempData=new double[theRawData.length];
    for(int i=0; i<theRawData.length; i++) {
      try{
        double d = Double.valueOf(theRawData[i][column].toString()).doubleValue();
        if (!Double.isNaN(d)) {
           tempData[countGood] = d;
           countGood++;
        }
      }
      catch (NumberFormatException ex) {}
    }
    double[] theGoodData=new double[countGood];
    for(int i=0; i<countGood; i++) {
      theGoodData[i]=tempData[i];
    }
    return theGoodData;
  }
  /**
   * Creates a chart.
   *
   * @param dataset  a dataset.
   *
   * @return The chart.
   */
  JFreeChart createChart(final IntervalXYDataset dataset) {
      String name="Non-numeric choice";
      if(dataset==null) {name="Non-numeric choice";}
      else if(dataset.getSeriesCount()==1) {name=dataset.getSeriesKey(0)+" Distribution";}
      else if(dataset.getSeriesCount()==2) {name=dataset.getSeriesKey(0)+" & "+dataset.getSeriesKey(1)+" Distribution";};
      final JFreeChart chart = ChartFactory.createHistogram(
          name,
          null,
          null,
          dataset,
          PlotOrientation.VERTICAL,
          true,
          false,
          false
      );
      chart.getXYPlot().setForegroundAlpha(0.75f);
      return chart;
  }

  void updateChart() {
    int series1, series2;
    series1 = series1ComboBox.getSelectedIndex() - 1;
    series2 = series2ComboBox.getSelectedIndex() - 1;
    dataset = createDataset(series1, series2);
    chart = createChart(dataset);
    chartPanel.setChart(chart);
  }

  void series1ComboBox_actionPerformed(ActionEvent e) {
    updateChart();
  }

  void series2ComboBox_actionPerformed(ActionEvent e) {
    updateChart();
  }

  void binsTextField_actionPerformed(ActionEvent e) {
    int newbin;
    try{newbin=Integer.parseInt(this.binsTextField.getText());
      bins=newbin;
      updateChart();
      }
    catch (NumberFormatException ex) {
      binsTextField.setText(""+bins);
    }
  }

  void binsTextField_caretUpdate(CaretEvent e) {
     binsTextField_actionPerformed(null);
  }

  public JComponent getMainComponent() {
        return chartPanel;
    }
}

