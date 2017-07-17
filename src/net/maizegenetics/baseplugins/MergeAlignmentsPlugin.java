/*
 * MergeAlignmentsPlugin
 */
package net.maizegenetics.baseplugins;

import java.awt.Frame;
import java.net.URL;

import java.util.List;

import javax.swing.ImageIcon;
import javax.swing.JOptionPane;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.MutableSingleEncodeAlignment;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;

import org.apache.log4j.Logger;

/**
 *
 * @author terry
 */
public class MergeAlignmentsPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(MergeAlignmentsPlugin.class);

    public MergeAlignmentsPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {

        List<Datum> inputs = input.getDataOfType(Alignment.class);

        if ((inputs == null) || (inputs.size() < 2)) {
            if (isInteractive()) {
                JOptionPane.showMessageDialog(getParentFrame(), "Must select at least two alignments.");
            } else {
                myLogger.warn("performFunction: Must select at least two alignments.");
            }
            return null;
        }

        try {
            Alignment[] alignments = new Alignment[inputs.size()];
            for (int i = 0; i < inputs.size(); i++) {
                alignments[i] = (Alignment) ((Datum) inputs.get(i)).getData();
            }

            Alignment alignment = MutableSingleEncodeAlignment.getInstance(alignments);
            DataSet result = new DataSet(new Datum("Merged Alignment", alignment, null), this);

            fireDataSetReturned(new PluginEvent(result, MergeAlignmentsPlugin.class));

            return result;
        } finally {
            fireProgress(100);
        }

    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = SeparatePlugin.class.getResource("images/Merge.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Merge Alignments";
    }

    @Override
    public String getToolTipText() {
        return "Merge Alignments";
    }
}
