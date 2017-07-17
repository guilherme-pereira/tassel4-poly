package net.maizegenetics.baseplugins.chart;

import org.jfree.chart.JFreeChart;

import javax.swing.*;
/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2004</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

public abstract class BasicChartPanel extends JPanel {
  JFreeChart chart;


  public BasicChartPanel() {
  }

  public JFreeChart getChart() {
    return chart;
  }

  public abstract JComponent getMainComponent();

}
