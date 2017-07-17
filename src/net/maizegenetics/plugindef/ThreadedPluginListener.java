/*
 * ThreadedPluginListener
 */
package net.maizegenetics.plugindef;

import net.maizegenetics.gui.DialogUtils;
import net.maizegenetics.util.Utils;

import org.apache.log4j.Logger;

/**
 *
 * @author terry
 */
public class ThreadedPluginListener extends Thread {

    private static final Logger myLogger = Logger.getLogger(ThreadedPluginListener.class);
    private final PluginListener myPluginListener;
    private final PluginEvent myEvent;

    public ThreadedPluginListener(PluginListener pluginListener, PluginEvent event) {
        myPluginListener = pluginListener;
        myEvent = event;
    }

    @Override
    public void run() {
        try {
            myPluginListener.dataSetReturned(myEvent);
        } catch (OutOfMemoryError e) {
            e.printStackTrace();
            StringBuilder builder = new StringBuilder();
            builder.append("Out of Memory: ");

            Plugin plugin = null;
            try {
                plugin = (Plugin) myPluginListener;
                builder.append(Utils.getBasename(plugin.getClass().getName()));
                builder.append(" could not complete task: ");
            } catch (Exception exp) {
                // do nothing
            }

            builder.append("\n");

            long heapMaxSize = Runtime.getRuntime().maxMemory() / 1048576l;
            builder.append("Current Max Heap Size: ");
            builder.append(heapMaxSize);
            builder.append(" Mb\n");
            builder.append("Use -Xmx option in start_tassel.pl or start_tassel.bat\n");
            builder.append("to increase heap size.");
            builder.append(" Included with tassel standalone zip.");

            String str = builder.toString();

            try {
                if (plugin.isInteractive()) {
                    DialogUtils.showError(str, plugin.getParentFrame());
                }
            } catch (Exception ex) {
                // do nothing
            }
            myLogger.error(str);
        }
    }

    public PluginListener getPluginListener() {
        return myPluginListener;
    }
}
