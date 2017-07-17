/*
 * ChartDisplayPlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.baseplugins.chart;

import net.maizegenetics.pal.report.TableReport;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.baseplugins.AbstractDisplayPlugin;

import javax.swing.*;
import java.awt.*;
import java.net.URL;
import java.util.List;

/**
 *
 * @author Ed Buckler
 */
public class ChartDisplayPlugin extends AbstractDisplayPlugin {

    public enum ChartType {

        Histogram, XYScatter, BarChart, PieChart
    };
    ChartType chartMode = ChartType.Histogram;
    int series1 = -1;   //X for XY, category for BarChart
    int series2 = -1;   //Y1 for XY,
    int series3 = -1;   //Y2 for XY, not used for histogram
    int bins = 5;
    boolean isRegression = false;  //XY setting
    boolean isBoxWhisker = false, isBar = true, errBar = false;  //Barchart settings

    /**
     * Creates a new instance of ChartDisplayPlugin
     */
    public ChartDisplayPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {

        try {

            List<Datum> tableInList = input.getDataOfType(TableReport.class);
            if (tableInList.size() != 1) {
                String message = "Invalid selection.  Please select one table result.";
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), message);
                } else {
                    System.out.println(message);
                }
                return null;
            }
            TableReport theTR = (TableReport) tableInList.get(0).getData();
            if (isInteractive()) {
                BasicGraphFrame myDialog = new BasicGraphFrame(this, theTR);
                myDialog.setLocationRelativeTo(getParentFrame());
                myDialog.setVisible(true);
            } else if (getSaveFile() != null) {
                BasicChartPanel chartPanel = null;
                switch (chartMode) {
                    case Histogram:
                        chartPanel = new HistogramPanel(theTR, series1, series2, bins);
                        break;
                    case XYScatter:
                        chartPanel = new XYScatterPanel(theTR, series1, series2, series3, isRegression);
                        break;
                    case BarChart:
                        chartPanel = new BarChartPanel(theTR, series1, series2, series3,
                                isBoxWhisker, isBar, errBar);
                        break;
                    case PieChart:
                        chartPanel = new PieChartPanel(theTR);
                        break;
                }
                saveDataToFile(chartPanel.getMainComponent(), getSaveFile());
            }

            return null;

        } finally {
            fireProgress(100);
        }

    }

    public ChartType getChartMode() {
        return chartMode;
    }

    public void setChartMode(ChartType chartMode) {
        this.chartMode = chartMode;
    }

    public int getSeries1() {
        return series1;
    }

    public void setSeries1(int series1) {
        this.series1 = series1;
    }

    public int getSeries2() {
        return series2;
    }

    public void setSeries2(int series2) {
        this.series2 = series2;
    }

    public int getSeries3() {
        return series3;
    }

    public void setSeries3(int series3) {
        this.series3 = series3;
    }

    public int getBins() {
        return bins;
    }

    public void setBins(int bins) {
        this.bins = bins;
    }

    public boolean isRegression() {
        return isRegression;
    }

    public void setRegression(boolean regression) {
        isRegression = regression;
    }

    public boolean isBoxWhisker() {
        return isBoxWhisker;
    }

    public void setBoxWhisker(boolean boxWhisker) {
        isBoxWhisker = boxWhisker;
    }

    public boolean isBar() {
        return isBar;
    }

    public void setBar(boolean bar) {
        isBar = bar;
    }

    public boolean isErrBar() {
        return errBar;
    }

    public void setErrBar(boolean errBar) {
        this.errBar = errBar;
    }

    public static String[] getPossibleCharts() {
        String[] s = new String[ChartType.values().length];
        int i = 0;
        for (ChartType p : ChartType.values()) {
            s[i] = "" + p;
            i++;
        }
        return s;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        URL imageURL = ChartDisplayPlugin.class.getResource("BarChart.gif");
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
        return "Chart";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Charting Tools";
    }
}
