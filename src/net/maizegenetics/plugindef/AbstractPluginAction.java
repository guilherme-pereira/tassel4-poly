/*
 * AbstractPluginAction.java
 *
 */

package net.maizegenetics.plugindef;

import javax.swing.*;


/**
 *
 * @author terryc
 */
abstract public class AbstractPluginAction extends AbstractAction implements PluginAction {
    
    private final Plugin myPlugin;
    
    
    /** Creates a new instance of AbstractPluginAction */
    public AbstractPluginAction(Plugin plugin, String name, Icon icon) {
        
        super(name, icon);
        
        if (plugin == null) {
            throw new IllegalArgumentException("AbstractPluginAction: init: plugin can not be null.");
        }
        
        myPlugin = plugin;
        
    }
    
    
    /**
     * Return associated plugin.
     *
     * @return tassel plugin
     */
    public Plugin getPlugin() {
        return myPlugin;
    }
    
    
    /**
     * GUI Panel for associated plugin.
     *
     * @return panel
     */
    public JPanel getPanel() {
        return myPlugin.getPanel();
    }
    
}
