package net.maizegenetics.baseplugins;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.pal.alignment.*;
import net.maizegenetics.pal.distance.Kinship;

import javax.swing.*;
import java.net.URL;
import java.awt.Container;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Author: Zhiwu Zhang
 * Date: Apr 29, 2007
 */
public class KinshipPlugin extends AbstractPlugin {

    public KinshipPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);

    }

    public DataSet performFunction(DataSet input) {

        try {

            List<Datum> alignInList = input.getDataSet();

            if (alignInList.isEmpty()) {
                String message = "Nothing selected. Please select pedigree data.";
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), message);
                } else {
                    System.out.println(message);
                }
                return null;
            }

            List result = new ArrayList();
            Iterator itr = alignInList.iterator();
            while (itr.hasNext()) {

                Datum current = (Datum) itr.next();
                String datasetName = current.getName();
                Kinship kin = null;

                try {

                    if (current.getData() instanceof Alignment) {
                        //this section implements additional options for calculating kinship
                        Alignment theAlignment = (Alignment) current.getData();
                        kin = new Kinship(theAlignment);

                    } else if (current.getData() instanceof SimplePhenotype) { //pedigree data
                        SimplePhenotype ped = (SimplePhenotype) current.getData();
                        kin = new Kinship(ped);
                    } else {
                        String message = "Invalid selection. Can't create kinship matrix from: " + datasetName;
                        if (isInteractive()) {
                            JOptionPane.showMessageDialog(getParentFrame(), message);
                        } else {
                            System.out.println(message);
                        }
                    }

                } catch (Exception e) {
                    e.printStackTrace();
                    String message = "Problem creating kinship matrix from: " + datasetName + "\n" + e.getClass().getName() + ": " + e.getMessage();
                    if (isInteractive()) {
                        JOptionPane.showMessageDialog(getParentFrame(), message);
                    } else {
                        System.out.println(message);
                        e.printStackTrace();
                    }
                }

                if (kin != null) {
                    //add kin to datatree;
                    DataSet ds = new DataSet(new Datum("kin_" + datasetName, kin.getDm(), "kinship matrix created from " + datasetName), this);
                    result.add(ds);
                    fireDataSetReturned(new PluginEvent(ds, KinshipPlugin.class));
                }

            }

            return DataSet.getDataSet(result, this);

        } finally {
            fireProgress(100);
        }

    }

    public ImageIcon getIcon() {
        URL imageURL = KinshipPlugin.class.getResource("images/Kin.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    public String getButtonName() {
        return "Kinship";
    }

    public String getToolTipText() {
        return "Calculate kinship from marker data";
    }

}
