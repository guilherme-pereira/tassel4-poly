/*
 * DistanceMatrixRangesPlugin
 */
package net.maizegenetics.baseplugins;

import java.awt.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import javax.swing.*;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.distance.IBSDistanceMatrix;
import net.maizegenetics.pal.report.SimpleTableReport;
import net.maizegenetics.pal.report.TableReport;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;
import org.apache.log4j.Logger;

/**
 *
 * @author terryc
 */
public class DistanceMatrixRangesPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(DistanceMatrixRangesPlugin.class);
    private int[] myPhysicalPositions = null;
    private String myLocus = null;
    private String myTaxon = null;

    public DistanceMatrixRangesPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {

        try {

            List<Datum> alignInList = input.getDataOfType(Alignment.class);
            if (alignInList.size() < 1) {
                String message = "Invalid selection.  Please select sequence alignment.";
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), message);
                } else {
                    myLogger.error(message);
                }
                return null;
            }

            List result = new ArrayList();
            Iterator<Datum> itr = alignInList.iterator();
            while (itr.hasNext()) {
                Datum current = itr.next();
                DataSet tds = processDatum(current);
                result.add(tds);
                if (tds != null) {
                    fireDataSetReturned(new PluginEvent(tds, DistanceMatrixRangesPlugin.class));
                }
            }

            return DataSet.getDataSet(result, this);

        } finally {
            fireProgress(100);
        }

    }

    public DataSet processDatum(Datum input) {

        Alignment aa = (Alignment) input.getData();
        int numTaxa = aa.getSequenceCount();
        int interestedTaxa = aa.getIdGroup().whichIdNumber(myTaxon);
        Object[][] theData = new Object[numTaxa][myPhysicalPositions.length];
        for (int t = 0; t < numTaxa; t++) {
            theData[t][0] = aa.getTaxaName(t);
        }
        String[] columnNames = new String[myPhysicalPositions.length];
        columnNames[0] = "Taxa";
        for (int i = 1, n = myPhysicalPositions.length; i < n; i++) {
            Alignment alignment = FilterAlignment.getInstance(aa, myLocus, myPhysicalPositions[i - 1], myPhysicalPositions[i]);
            columnNames[i] = String.valueOf(myPhysicalPositions[i - 1] + "-" + String.valueOf(myPhysicalPositions[i]));
            if (alignment == null) {
                for (int t = 0; t < numTaxa; t++) {
                    theData[t][i] = Double.NaN;
                }
            } else {
                for (int t = 0; t < numTaxa; t++) {
                    theData[t][i] = Double.valueOf(IBSDistanceMatrix.computeHetBitDistances(alignment, interestedTaxa, t)[0]);
                }
            }
        }

        TableReport result = new SimpleTableReport("Distance Matrix Ranges", columnNames, theData);

        return new DataSet(new Datum("MatrixRanges:" + input.getName(), result, "Distance Matrix Ranges"), this);

    }

    public void setPhysicalPositions(String[] positions) {
        int[] result = new int[positions.length];
        for (int i = 0; i < positions.length; i++) {
            result[i] = Integer.valueOf(positions[i]);
        }
        setPhysicalPositions(result);
    }

    public void setPhysicalPositions(int[] positions) {
        myPhysicalPositions = positions;
    }

    public int[] getPhysicalPositions() {
        return myPhysicalPositions;
    }

    public void setLocus(String locus) {
        myLocus = locus;
    }

    public String getLocus() {
        return myLocus;
    }

    public void setTaxon(String taxon) {
        myTaxon = taxon;
    }

    public String getTaxon() {
        return myTaxon;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        return null;
    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    public String getButtonName() {
        return "Distance Matrix";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Create a distance matrix";
    }
}
