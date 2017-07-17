/*
 * DistanceMatrixPlugin
 */
package net.maizegenetics.baseplugins;

import java.awt.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import javax.swing.*;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.distance.DistanceMatrix;
import net.maizegenetics.pal.distance.IBSDistanceMatrix;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.util.ProgressListener;
import org.apache.log4j.Logger;

/**
 *
 * @author terryc
 */
public class DistanceMatrixPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(DistanceMatrixPlugin.class);

    public DistanceMatrixPlugin(Frame parentFrame, boolean isInteractive) {
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
                    fireDataSetReturned(new PluginEvent(tds, DistanceMatrixPlugin.class));
                }
            }

            return DataSet.getDataSet(result, this);

        } finally {
            fireProgress(100);
        }

    }

    public DataSet processDatum(Datum input) {
        Alignment aa = (Alignment) input.getData();
        DistanceMatrix adm = getDistanceMatrix(aa, this);
        return new DataSet(new Datum("Matrix:" + input.getName(), adm, "Distance Matrix"), this);
    }

    public static DistanceMatrix getDistanceMatrix(Alignment alignment) {
        return new IBSDistanceMatrix(alignment);
    }

    public static DistanceMatrix getDistanceMatrix(Alignment alignment, ProgressListener listener) {
        return new IBSDistanceMatrix(alignment, listener);
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
