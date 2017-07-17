/*
 * XYAndLineToolTipGenerator
 */
package net.maizegenetics.baseplugins.chart;

import java.text.DecimalFormat;
import java.util.Map;
import org.jfree.chart.labels.AbstractXYItemLabelGenerator;
import org.jfree.chart.labels.XYToolTipGenerator;
import org.jfree.data.xy.XYDataset;

/**
 *
 * @author yz79
 */
public class XYAndLineToolTipGenerator extends AbstractXYItemLabelGenerator implements XYToolTipGenerator {

    public void XYAndLineToolTipGenerator() {
    }

    @Override
    public String generateToolTip(XYDataset dataset, int series, int item) {
        TableReportQQDataset myDataset = (TableReportQQDataset) dataset;
        int[] positions = myDataset.getPositions();
        String[] markers = myDataset.getMarkers();
        Map myTable = myDataset.getLookupTable();
        int index = (Integer) myTable.get(myDataset.getYValue(series, item));
        DecimalFormat df = new DecimalFormat("#0.000");
        StringBuilder sb = new StringBuilder("SNP ID: ");
        sb.append(markers[index]);
        sb.append(", -Log(P-Value): ");
        sb.append(df.format(myDataset.getYValue(series, item)));
        sb.append(", Expected: ");
        sb.append(df.format(myDataset.getXValue(series, item)));
        sb.append(", Position: ");
        sb.append(positions[index]);
        return sb.toString();
    }
}
