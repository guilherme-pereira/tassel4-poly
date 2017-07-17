/*
 * IntersectionAlignmentPlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.baseplugins;

import java.awt.Frame;

import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;

import javax.swing.*;
import java.net.URL;

/**
 *
 * @author Ed Buckler
 */
public class IntersectionAlignmentPlugin extends UnionAlignmentPlugin {

    /**
     * Creates a new instance of IntersectionAlignmentPlugin
     */
    public IntersectionAlignmentPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {
        Datum joinedDatum = processData(input, false);
        if (joinedDatum == null) {
            return null;
        }
        DataSet output = new DataSet(joinedDatum, this);
        //I am setting the firing class as the metadata - so that the control panel know where the event is coming from
        fireDataSetReturned(new PluginEvent(output, IntersectionAlignmentPlugin.class));
        return output;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    @Override
    public ImageIcon getIcon() {
        URL imageURL = IntersectionAlignmentPlugin.class.getResource("images/IntersectJoin.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    @Override
    public String getButtonName() {
        return "Intersect Join";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    @Override
    public String getToolTipText() {
        return "Join Datasets by Intersecting Taxa";
    }
}
