/*
 * PluginAction.java
 */

package net.maizegenetics.plugindef;

import javax.swing.*;


/**
 *
 * @author terryc
 */
public interface PluginAction extends Action{
    
    
    /**
     * Return associated plugin.
     *
     * @return tassel plugin
     */
    public Plugin getPlugin();
    
    
    /**
     * GUI Panel for associated plugin.
     *
     * @return panel
     */
    public JPanel getPanel();
    
}
