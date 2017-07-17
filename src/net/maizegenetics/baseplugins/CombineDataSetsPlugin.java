/*
 * CombineDataSetsPlugin.java
 *
 * Created on January 5, 2007, 2:25 AM
 *
 */
package net.maizegenetics.baseplugins;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Plugin;
import net.maizegenetics.plugindef.PluginEvent;

import javax.swing.*;

import java.util.Iterator;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author Terry Casstevens
 */
public class CombineDataSetsPlugin extends AbstractPlugin {

    private Map myDataSets = new LinkedHashMap();
    private Map myOnceDataSets = new LinkedHashMap();

    /**
     * Creates a new instance of CombineDataSetsPlugin
     */
    public CombineDataSetsPlugin() {
        super(null, false);
    }

    /**
     * Returns combined data set if all inputs have been received.
     *
     * @param dataSet Not used.
     */
    public DataSet performFunction(DataSet dataSet) {

        try {

            List dataSets = null;
            synchronized (myDataSets) {
                if ((myDataSets.containsValue(null)) || myOnceDataSets.containsValue(null)) {
                    return null;
                }

                dataSets = new ArrayList();

                dataSets.addAll(myDataSets.values());
                dataSets.addAll(myOnceDataSets.values());

                reset();
            }

            DataSet result = DataSet.getDataSet(dataSets, this);
            fireDataSetReturned(result);

            return result;

        } finally {
            fireProgress(100);
        }

    }

    /**
     * Same as performFunction except this doesn't wait to receive all inputs.
     */
    public void flush() {

        List dataSets = new ArrayList();

        dataSets.addAll(myDataSets.values());
        dataSets.addAll(myOnceDataSets.values());

        reset();

        if (dataSets.size() > 0) {
            fireDataSetReturned(DataSet.getDataSet(dataSets, this));
        }

    }

    public void reset() {

        // Clear only values.
        // Method dataSetReturned knows what inputs
        // to expect based on keys stored here.
        Set keys = myDataSets.keySet();
        for (Iterator itr = keys.iterator(); itr.hasNext();) {
            myDataSets.put(itr.next(), null);
        }

    }

    public String getToolTipText() {
        return "";
    }

    public ImageIcon getIcon() {
        return null;
    }

    public String getButtonName() {
        return "Combine";
    }

    public void dataSetReturned(PluginEvent event) {

        DataSet input = (DataSet) event.getSource();
        Plugin creator = input.getCreator();

        if (myOnceDataSets.containsKey(creator)) {
            Object value = myOnceDataSets.get(creator);
            if (value != null) {
                throw new IllegalStateException("CombineDataSetsPlugin: dataSetReturned: this plugin should only return data once: " + creator);
            } else {
                myOnceDataSets.put(creator, input);
            }
        } else if (myDataSets.containsKey(creator)) {
            Object value = myDataSets.get(creator);
            if (value != null) {
                throw new IllegalStateException("CombineDataSetsPlugin: dataSetReturned: this plugin should only return data once per iteration: " + creator);
            } else {
                myDataSets.put(creator, input);
            }
        } else {
            throw new IllegalStateException("CombineDataSetsPlugin: dataSetReturned: can not receive data from unknown plugin: " + creator);
        }

        performFunction(null);

    }

    /**
     * Add given plugin as source to receive data sets only once and use that
     * data set in every resulting output.
     *
     * @param plugin plugin
     */
    public void receiveDataSetOnceFrom(Plugin plugin) {
        super.receiveInput(plugin);
        myOnceDataSets.put(plugin, null);
    }

    /**
     * Add given plugin as source to receive data sets iteratively.
     *
     * @param plugin plugin
     */
    public void receiveDataSetFrom(Plugin plugin) {
        super.receiveInput(plugin);
        myDataSets.put(plugin, null);
    }

    /**
     * Sets up this plugin to receive input from another plugin.
     *
     * @param input input
     */
    public void receiveInput(Plugin input) {
        receiveDataSetFrom(input);
    }

    public String toString() {

        StringBuilder str = new StringBuilder();

        Iterator itr = myDataSets.values().iterator();
        while (itr.hasNext()) {
            DataSet current = (DataSet) itr.next();
            if (current != null) {
                str.append(current.toString());
            }
        }

        return str.toString();

    }
}
