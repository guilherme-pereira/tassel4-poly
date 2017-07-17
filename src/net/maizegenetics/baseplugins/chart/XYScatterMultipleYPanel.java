/*
 * XYScatterMultipleYPanel
 */
package net.maizegenetics.baseplugins.chart;

import java.awt.Color;
import java.awt.Component;
import java.awt.event.ActionEvent;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComponent;
import net.maizegenetics.baseplugins.ManhattanDisplayPlugin;
import net.maizegenetics.pal.report.TableReport;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StandardXYItemRenderer;
import org.jfree.data.xy.XYDataset;

/**
 *
 * @author yz79
 */
public class XYScatterMultipleYPanel extends BasicChartPanel {

    ManhattanDisplayPlugin myManhattanDisplayPlugin;
    ChartPanel myChartPanel;
    JButton saveButton = new JButton("save...");
    TableReportManhattanDataset dataset;
    TableReport myTableReport;

    public XYScatterMultipleYPanel(ManhattanDisplayPlugin plugin, TableReport theTable, int start, int end) {
        myManhattanDisplayPlugin = plugin;
        myTableReport = theTable;
        try {
            dataset = new TableReportManhattanDataset(theTable, start, end);
            chart = createChart(dataset);
            myChartPanel = new ChartPanel(chart);
            myChartPanel.setPreferredSize(new java.awt.Dimension(900, 500));
            myTableReport = theTable;
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
        myManhattanDisplayPlugin.saveDataToFile(myChartPanel);
    }

    public JFreeChart createChart(TableReportManhattanDataset dataset) {
        String name = "Please select numeric variables";
        String xName = "X";
        String y1Name = "Y";
        String y2Name = "Y2";
        if (dataset != null) {
            xName = dataset.getXName();
            y1Name = "-Log(P-Value)";
            name = "P-Values by Chromosome for " + dataset.getTrait();
//            if(dataset.getSeriesCount()==2) {
//                y2Name=dataset.getSeriesName(1);
//                name+=" and "+y2Name;
//            }
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
        chart.getXYPlot().getRenderer().setToolTipGenerator(new XYMultipleYToolTipGenerator());
        return chart;
    }

//    private void createLine(XYPlot plot, XYDataset data, int series, Color theColor) {
//        Function2D curve = new LineFunction2D(0, 1);
//        double max=DatasetUtilities.findMaximumDomainValue(data).doubleValue();
//        double min=DatasetUtilities.findMinimumDomainValue(data).doubleValue();
//        XYDataset regressionData = DatasetUtilities.sampleFunction2D(
//            curve,
//            min,max, 2,
//            data.getSeriesKey(series) + "Test"
//        );
//        int datasetCount=plot.getDatasetCount();
//        plot.setDataset(datasetCount, regressionData);
//        XYItemRenderer renderer2 = new StandardXYItemRenderer(StandardXYItemRenderer.LINES);
//        renderer2.setSeriesPaint(0, theColor);
//        plot.setRenderer(datasetCount, renderer2);
//        setAxis(plot, max);
//    }
    private void create2ndAxis(XYPlot plot, XYDataset data, int series, Color theColor, boolean isLine) {
        //final NumberAxis axis2 = new NumberAxis(data.getSeriesName(series));
        final NumberAxis axis2 = new NumberAxis((String) data.getSeriesKey(series));
        axis2.setAutoRangeIncludesZero(false);
        plot.setRangeAxis(1, axis2);
        plot.setDataset(1, data);
        plot.mapDatasetToRangeAxis(1, 1);
        StandardXYItemRenderer renderer2; //= new StandardXYItemRenderer();
        //   renderer2.setSeriesPaint(0, theColor);
        //    if(isLine) {renderer2.setPlotShapes(false); renderer2.setPlotLines(true);}
        //    else {renderer2.setPlotShapes(true); renderer2.setPlotLines(false);}
        if (isLine) {
            renderer2 = new StandardXYItemRenderer(StandardXYItemRenderer.LINES);
        } else {
            renderer2 = new StandardXYItemRenderer(StandardXYItemRenderer.SHAPES);
        }
        renderer2.setSeriesPaint(0, theColor);
        plot.setRenderer(1, renderer2);
    }

//    private void setAxis(XYPlot plot, double domainMax) {
//        ValueAxis xAxis = plot.getDomainAxis();
//        ValueAxis yAxis = plot.getRangeAxis();
//        double xMax = domainMax;
//        double yMax = yAxis.getUpperBound();
//        xAxis.setRange(0, xMax);
//        yAxis.setRange(0, yMax);
//        plot.setDomainAxis(xAxis);
//        plot.setRangeAxis(yAxis);
//    }
    public JComponent getMainComponent() {
        return myChartPanel;
    }
}