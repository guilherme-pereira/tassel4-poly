/*
 * PluginEvent.java
 *
 */

package net.maizegenetics.plugindef;


import java.util.EventObject;


/**
 *
 * @author terryc
 */
public class PluginEvent extends EventObject {
    
    private final Object myMetaData;
    
    
    /** Creates a new instance of PluginEvent */
    public PluginEvent(Object source) {
        super(source);
        myMetaData = null;
    }
    
    
    /** Creates a new instance of PluginEvent */
    public PluginEvent(Object source, Object metaData) {
        super(source);
        myMetaData = metaData;
    }
    
    
    public Object getMetaData() {
        return myMetaData;
    }
    
}
