/*
 * PluginListener.java
 */

package net.maizegenetics.plugindef;



/**
 *
 * @author terryc
 */
public interface PluginListener {
    
    /**
     * Returns Tassel data set after complete.
     *
     * @param event event
     */
    public void dataSetReturned (PluginEvent event);
    
    
    /**
     * Returns progress of execution.
     *
     * @param event event
     */
    public void progress (PluginEvent event);
    
}
