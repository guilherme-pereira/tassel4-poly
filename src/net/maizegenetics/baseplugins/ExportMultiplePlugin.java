/*
 * ExportMultiplePlugin.java
 *
 * Created on December 21, 2010
 *
 */
package net.maizegenetics.baseplugins;

import java.awt.Frame;

import java.util.List;

import javax.swing.ImageIcon;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

import org.apache.log4j.Logger;

/**
 *
 * @author terry
 */
public class ExportMultiplePlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(ExportMultiplePlugin.class);
    private FileLoadPlugin.TasselFileType[] myFileTypes = null;
    private String[] mySaveFiles = null;
    private final ExportPlugin myExportPlugin;

    /** Creates a new instance of ExportMultiplePlugin */
    public ExportMultiplePlugin(Frame parentFrame) {
        super(parentFrame, false);
        myExportPlugin = new ExportPlugin(parentFrame, false);
    }

    public DataSet performFunction(DataSet input) {

        List data = input.getDataSet();

        int numSaveFiles = 0;
        if (mySaveFiles != null) {
            numSaveFiles = mySaveFiles.length;
        }

        if ((numSaveFiles != 0) && (numSaveFiles != 1) && (numSaveFiles != data.size())) {
            throw new IllegalStateException("ExportMultiplePlugin: performFunction: number of save files should be either 0, 1 or number of input data sets.");
        }

        if ((myFileTypes != null) && (myFileTypes.length != 0)) {
            if ((myFileTypes.length != 1) && (myFileTypes.length != data.size())) {
                throw new IllegalStateException("ExportMultiplePlugin: performFunction: number of files types should be either 0, 1 or number of input data sets.");
            }
        }

        for (int i = 0, n = data.size(); i < n; i++) {

            Datum datum = (Datum) data.get(i);
            DataSet current = new DataSet(datum, input.getCreator());

            if (numSaveFiles == 0) {
                myExportPlugin.setSaveFile(datum.getName());
            } else if (numSaveFiles == 1) {
                String temp;
                if (data.size() == 1) {
                    temp = mySaveFiles[0];
                } else {
                    temp = mySaveFiles[0].replaceFirst("\\.", (i + 1) + ".");
                    if (temp.length() == mySaveFiles[0].length()) {
                        temp = mySaveFiles[0] + (i + 1);
                    }
                }
                myExportPlugin.setSaveFile(temp);
            } else {
                myExportPlugin.setSaveFile(mySaveFiles[i]);
            }

            if ((myFileTypes == null) || (myFileTypes.length == 0)) {
                // do nothing
            } else if (myFileTypes.length == 1) {
                myExportPlugin.setAlignmentFileType(myFileTypes[0]);
            } else {
                myExportPlugin.setAlignmentFileType(myFileTypes[i]);
            }

            DataSet filename = myExportPlugin.performFunction(current);
            myLogger.info("Datum: " + datum.getName() + " Written to: " + filename.getDataOfType(String.class).get(0).getData());

        }

        return null;

    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        return null;
    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    public String getButtonName() {
        return "Export Multiple";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Export multiple data sets to files on your computer.";
    }

    public String[] getSaveFiles() {
        return mySaveFiles;
    }

    public void setSaveFiles(String[] saveFiles) {
        mySaveFiles = saveFiles;
    }

    public void setAlignmentFileTypes(FileLoadPlugin.TasselFileType[] types) {
        myFileTypes = types;
    }

    public void setAlignmentFileType(FileLoadPlugin.TasselFileType type) {
        myFileTypes = new FileLoadPlugin.TasselFileType[]{type};
    }
}
