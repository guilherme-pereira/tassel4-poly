/*
 * FilterDataSetPlugin.java
 *
 * Created on January 31, 2007, 2:55 PM
 *
 */

package net.maizegenetics.baseplugins;

import javax.swing.ImageIcon;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;

/**
 *
 * @author terryc
 */
public class FilterDataSetPlugin extends AbstractPlugin{
    
    private final String [] myNames;
    private final Class [] myTypes;
    
    
    /**
     * Creates a new instance of FilterDataSetPlugin
     */
    public FilterDataSetPlugin(String [] names, Class [] types) {
        super(null, false);
        myNames = names;
        myTypes = types;
    }
    
    public DataSet performFunction(DataSet input) {
        DataSet result = new DataSet(input.getDataOfTypeWithName(myTypes, myNames), this);
        fireDataSetReturned(result);
        return result;
    }
    
    public String getToolTipText() {
        return "";
    }
    
    public ImageIcon getIcon() {
        return null;
    }
    
    public String getButtonName() {
        return "";
    }
    
}
