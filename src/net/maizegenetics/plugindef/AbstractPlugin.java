/*
 * AbstractPlugin.java
 *
 * Created on December 22, 2006, 5:03 PM
 *
 */
package net.maizegenetics.plugindef;

import java.awt.Frame;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import javax.swing.JMenu;
import javax.swing.JPanel;

/**
 *
 * @author terryc
 */
abstract public class AbstractPlugin implements Plugin {

    private final List myListeners = new ArrayList();
    private final List<Plugin> myInputs = new ArrayList();
    private final Frame myParentFrame;
    private final boolean myIsInteractive;
    private boolean myTrace = false;
    private boolean myThreaded = false;

    /** Creates a new instance of AbstractPlugin */
    public AbstractPlugin() {
        this(null, true);
    }

    /** Creates a new instance of AbstractPlugin */
    public AbstractPlugin(Frame parentFrame, boolean isInteractive) {
        myParentFrame = parentFrame;
        myIsInteractive = isInteractive;
    }

    /**
     * Returns menu that can be added to main menu bar.
     *
     * @return menu
     */
    public JMenu getMenu() {
        return null;
    }

    /**
     * Sets up this plugin to receive input from another plugin.
     *
     * @param input input
     */
    public void receiveInput(Plugin input) {

        if (input == null) {
            throw new IllegalArgumentException("AbstractPlugin: receiveInput: input can not be null.");
        }

        if (!myInputs.contains(input)) {
            myInputs.add(input);
        }

        input.addListener(this);

    }

    /**
     * GUI Panel for this plugin.
     *
     * @return panel
     */
    public JPanel getPanel() {
        return null;
    }

    /**
     * If interactive = true, the plugin will create dialogs and panels to interacts with the user
     *
     * @return boolean
     */
    public boolean isInteractive() {
        return myIsInteractive;
    }

    /**
     * Parent Frame for this plugin.  Can be null.
     *
     * @return frame
     */
    public Frame getParentFrame() {
        return myParentFrame;
    }

    /**
     *  Adds listener to this plugin.
     *
     * @param listener listener to add
     */
    public void addListener(PluginListener listener) {

        synchronized (myListeners) {
            if ((listener != null) && (!myListeners.contains(listener))) {
                myListeners.add(listener);
            }
        }

    }

    protected List getListeners() {
        return myListeners;
    }

    public List<Plugin> getInputs() {
        return myInputs;
    }

    /**
     * Returns data set after complete.
     *
     * @param event event
     */
    protected void fireDataSetReturned(PluginEvent event) {

        synchronized (myListeners) {
            Iterator itr = myListeners.iterator();
            while (itr.hasNext()) {
                try {
                    if (myThreaded) {
                        PluginListener current = (PluginListener) itr.next();
                        ThreadedPluginListener thread = new ThreadedPluginListener(current, event);
                        thread.start();
                    } else {
                        PluginListener current = (PluginListener) itr.next();
                        current.dataSetReturned(event);
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }

    }

    /**
     * Returns data set after complete.
     *
     * @param data data set
     */
    protected void fireDataSetReturned(DataSet data) {
        fireDataSetReturned(new PluginEvent(data));
    }

    /**
     * Returns progress of execution.
     *
     * @param event event
     */
    protected void fireProgress(PluginEvent event) {

        synchronized (myListeners) {
            Iterator itr = myListeners.iterator();
            while (itr.hasNext()) {
                PluginListener current = (PluginListener) itr.next();
                current.progress(event);
            }
        }

    }

    /**
     * Returns progress of execution.
     *
     * @param percent percentage between 0 and 100 inclusive.
     */
    protected void fireProgress(Integer percent) {

        if ((percent < 0) || (percent > 100)) {
            throw new IllegalArgumentException("AbstractPlugin: fireProgress: percent must be between 0 and 100 inclusive.  arg: " + percent);
        }

        Datum percentage = new Datum("Percent", percent, null);
        fireProgress(new PluginEvent(new DataSet(percentage, this)));

    }

    //
    // Methods for PluginListener.
    //
    /**
     * Returns data set after complete.
     *
     * @param event event
     */
    public void dataSetReturned(PluginEvent event) {

        DataSet input = (DataSet) event.getSource();

        performFunction(input);

    }

    /**
     * No operation for this abstract class.
     */
    public void progress(PluginEvent event) {
        // The default action of a plugin is to do
        // nothing when another plugin reports its
        // progress.   This is intended to be implemented
        // by GUI applications to show the user the
        // progress of an interactive action.
    }

    public void reverseTrace(int indent) {

        if (myTrace) {
            return;
        }

        indent(indent);
        System.out.println(getClass().getName());

        Iterator itr = myInputs.iterator();
        while (itr.hasNext()) {
            try {
                AbstractPlugin current = (AbstractPlugin) itr.next();
                current.reverseTrace(indent + 3);
            } catch (Exception e) {
                // do nothing
            }
        }

        myTrace = true;

    }

    public void trace(int indent) {

        if (myTrace) {
            return;
        }

        indent(indent);
        System.out.println(getClass().getName());

        Iterator itr = myListeners.iterator();
        while (itr.hasNext()) {
            try {
                AbstractPlugin current = (AbstractPlugin) itr.next();
                current.trace(indent + 3);
            } catch (Exception e) {
                // do nothing
            }
        }

        myTrace = true;

    }

    private void indent(int indent) {

        for (int i = 0; i < indent; i++) {
            System.out.print(" ");
        }

    }

    public void setThreaded(boolean threaded) {
        myThreaded = threaded;
    }

    public boolean cancel() {
        return false;
    }

    public void run() {
        performFunction(null);
    }

    @Override
    public void progress(int percent, Object meta) {
        fireProgress(percent);
    }

    @Override
    public void setParameters(String [] args) {
        throw new UnsupportedOperationException();
    }
}
