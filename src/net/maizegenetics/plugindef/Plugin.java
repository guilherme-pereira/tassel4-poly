/*
 * Plugin.java
 *
 */
package net.maizegenetics.plugindef;

import net.maizegenetics.util.ProgressListener;

import javax.swing.*;
import java.awt.*;


/**
 *
 * @author terryc
 */
public interface Plugin extends PluginListener, ProgressListener, Runnable {

    /**
     * Returns menu that can be added to main menu bar.
     *
     * @return menu
     */
    public JMenu getMenu();

    /**
     * Performs function of this plugin.
     *
     * @param input input
     *
     * @return resulting data set or null.
     */
    public DataSet performFunction(DataSet input);

    /**
     * Sets up this plugin to receive input from another plugin.
     *
     * @param input input
     */
    public void receiveInput(Plugin input);

    /**
     * GUI Panel for this plugin.
     *
     * @return panel
     */
    public JPanel getPanel();

    /**
     * If interactive = true, the plugin will create dialogs and panels to interacts with the user
     *
     * @return boolean
     */
    public boolean isInteractive();

    /**
     * Parent Frame for this plugin.  Can be null.
     *
     * @return frame
     */
    public Frame getParentFrame();

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon();

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    public String getButtonName();

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText();

    /**
     *  Adds listener to this plugin.
     *
     * @param listener listener to add
     */
    public void addListener(PluginListener listener);

    /**
     * Set whether this plugin is threaded.
     *
     * @param threaded whether to be threaded.
     */
    public void setThreaded(boolean threaded);

    /**
     * Attempt to cancel processing.
     *
     * @return true if plugin will cancel itself.
     */
    public boolean cancel();

    /**
     * Allows self-describing Plugins to use args
     * to set parameters specific to itself.
     *
     * @param args arguments
     */
    public void setParameters(String [] args);
}
