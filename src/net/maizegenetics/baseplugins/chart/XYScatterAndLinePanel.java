/*
 * XYScatterAndLinePanel
 */
package net.maizegenetics.baseplugins.chart;

import java.awt.Color;
import java.awt.Component;
import java.awt.event.ActionEvent;
import java.util.ArrayList;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComponent;
import net.maizegenetics.baseplugins.QQDisplayPlugin;
import net.maizegenetics.pal.report.TableReport;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StandardXYItemRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.function.Function2D;
import org.jfree.data.function.LineFunction2D;
import org.jfree.data.general.DatasetUtilities;
import org.jfree.data.xy.XYDataset;

/**
 *
 * @author yz79
 */
public class XYScatterAndLinePanel extends BasicChartPanel {

    QQDisplayPlugin myQQDisplayPlugin;
    ChartPanel myChartPanel;
    JButton saveButton = new JButton("save...");
    TableReportQQDataset[] datasets;
    TableReport myTableReport;

    public XYScatterAndLinePanel(QQDisplayPlugin plugin, TableReport table, int countToDisplay, ArrayList<Integer> tableIndices, int[] indices) {
        myQQDisplayPlugin = plugin;
        myTableReport = table;
//        ArrayList<Integer> indexes = splitTable(table);
        datasets = new TableReportQQDataset[indices.length];
        for (int i = 0; i < datasets.length; i++) {
            datasets[i] = new TableReportQQDataset(table, tableIndices.get(indices[i] * 2).intValue(), tableIndices.get(indices[i] * 2 + 1).intValue(), countToDisplay);
        }
        try {
            chart = createChart(datasets[0]);
            myChartPanel = new ChartPanel(chart);
            myChartPanel.setPreferredSize(new java.awt.Dimension(900, 500));
            for (int i = 1; i < datasets.length; i++) {
                addSeries(chart.getXYPlot(), datasets[i], i);
            }
            createLine(chart.getXYPlot(), datasets[0], 0, Color.BLACK);
            jbInit();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    private void jbInit() throws Exception {
        this.setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
        this.add(myChartPanel);
        saveButton.setAlignmentX(Component.RIGHT_ALIGNMENT);
        saveButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                saveButton_actionPerformed(e);
            }
        });
        this.add(saveButton);
    }

    private void saveButton_actionPerformed(ActionEvent e) {
        myQQDisplayPlugin.saveDataToFile(myChartPanel);
    }

    public JFreeChart createChart(TableReportQQDataset dataset) {
        String name = "Please select numeric variables";
        String xName = "X";
        String y1Name = "Y";
        String y2Name = "Y2";
        if (dataset != null) {
            xName = dataset.getXName();
            y1Name = "-Log(P-Value)";
            name = xName + " vs. " + y1Name;
            if (dataset.getSeriesCount() == 2) {
                y2Name = dataset.getSeriesName(1);
                name += " and " + y2Name;
            }
        }
        chart = ChartFactory.createScatterPlot(
                name,
                xName, y1Name,
                dataset,
                PlotOrientation.VERTICAL,
                true,
                true,
                false);
        chart.getXYPlot().setForegroundAlpha(0.75f);
        chart.getXYPlot().getRenderer().setToolTipGenerator(new XYAndLineToolTipGenerator());
        return chart;
    }

    private void createLine(XYPlot plot, XYDataset data, int series, Color theColor) {
        Function2D curve = new LineFunction2D(0, 1);
        double max = DatasetUtilities.findMaximumDomainValue(data).doubleValue();
        double min = DatasetUtilities.findMinimumDomainValue(data).doubleValue();
        XYDataset regressionData = DatasetUtilities.sampleFunction2D(
                curve,
                min, max, 2,
                "Expected Values");
        int datasetCount = plot.getDatasetCount();
        plot.setDataset(datasetCount, regressionData);
        XYItemRenderer renderer2 = new StandardXYItemRenderer(StandardXYItemRenderer.LINES);
        renderer2.setSeriesPaint(0, theColor);
        plot.setRenderer(datasetCount, renderer2);
        setAxis(plot, max, min);
    }

    private void addSeries(XYPlot plot, XYDataset data, int index) {
        plot.setDataset(index, data);
        XYLineAndShapeRenderer renderer2 = new XYLineAndShapeRenderer();
//        renderer2.setSeriesPaint(0);
        renderer2.setBaseShapesVisible(true);
        renderer2.setBaseLinesVisible(false);
        renderer2.setToolTipGenerator(new XYAndLineToolTipGenerator());
        plot.setRenderer(index, renderer2);
    }

    private void setAxis(XYPlot plot, double domainMax, double domainMin) {
        ValueAxis xAxis = plot.getDomainAxis();
        ValueAxis yAxis = plot.getRangeAxis();
        double xMax = domainMax;
        double xMin = domainMin;
        double yMax = yAxis.getUpperBound();
        double yMin = yAxis.getLowerBound();
        xAxis.setRange(0, xMax);
        yAxis.setRange(0, yMax);
        plot.setDomainAxis(xAxis);
        plot.setRangeAxis(yAxis);
    }

    public JComponent getMainComponent() {
        return myChartPanel;
    }

    public ArrayList<Integer> splitTable(TableReport table) {
        ArrayList<Integer> indexes = new ArrayList<Integer>();
        int numRows = table.getRowCount();
        String previousTrait = "";
        for (int i = 0; i < numRows; i++) {
            if (!previousTrait.equals((String) table.getValueAt(i, 0))) {
                if (!((String) table.getValueAt(i, 1)).equals("None")) {
                    indexes.add(new Integer(i));
                    previousTrait = (String) table.getValueAt(i, 0);
                    if (i > 1) {
                        indexes.add(new Integer(i));
                    }
                } else if (i != 0) {
                    indexes.add(new Integer(i));
                    indexes.add(new Integer(i + 1));
                    previousTrait = (String) table.getValueAt(i + 1, 0);
                }
            }
        }
        indexes.add(new Integer(numRows));
        return indexes;
    }
}