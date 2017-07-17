/*
 * TextToBinaryPlugin
 */
package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;
import java.io.File;

import javax.swing.ImageIcon;

import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.tagdist.TagsByTaxaBit;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByte;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.ArgsEngine;

import org.apache.log4j.Logger;

/**
 *
 * @author terry
 */
public class TextToBinaryPlugin extends AbstractPlugin {

    private Logger myLogger = Logger.getLogger(TextToBinaryPlugin.class);
    private ArgsEngine myEngine = null;
    private String myInput;
    private String myOutput;
    private BinaryToTextPlugin.FILE_TYPES myType;

    public TextToBinaryPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    @Override
    public DataSet performFunction(DataSet input) {

        switch (getType()) {
            case TOPM:
                TagsOnPhysicalMap topm = new TagsOnPhysicalMap(myInput, false);
                topm.writeBinaryFile(new File(myOutput));
                break;
            case TagCounts:
                TagCounts tc = new TagCounts(myInput, FilePacking.Text);
                tc.writeTagCountFile(myOutput, FilePacking.Bit, 0);
                break;
            case TBTBit:
                TagsByTaxaBit tbtbit = new TagsByTaxaBit(myInput, FilePacking.Text);
                tbtbit.writeDistFile(new File(myOutput), FilePacking.Bit, 0);
                break;
            case TBTByte:
                TagsByTaxaByte tbtbyte = new TagsByTaxaByte(myInput, FilePacking.Text);
                tbtbyte.writeDistFile(new File(myOutput), FilePacking.Byte, 0);
                break;
        }

        return null;

    }

    @Override
    public void setParameters(String[] args) {
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        if (myEngine == null) {
            myEngine = new ArgsEngine();
            myEngine.add("-i", "--input-file", true);
            myEngine.add("-o", "--output-file", true);
            myEngine.add("-t", "--file-type", true);
        }

        myEngine.parse(args);

        if (myEngine.getBoolean("-i")) {
            myInput = myEngine.getString("-i");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify the input file name.");
        }

        if (myEngine.getBoolean("-o")) {
            myOutput = myEngine.getString("-o");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify the output file name.");
        }

        if (myEngine.getBoolean("-t")) {
            String temp = myEngine.getString("-t");
            if (temp.equalsIgnoreCase(BinaryToTextPlugin.FILE_TYPES.TOPM.toString())) {
                setType(BinaryToTextPlugin.FILE_TYPES.TOPM);
            } else if (temp.equalsIgnoreCase(BinaryToTextPlugin.FILE_TYPES.TagCounts.toString())) {
                setType(BinaryToTextPlugin.FILE_TYPES.TagCounts);
            } else if (temp.equalsIgnoreCase(BinaryToTextPlugin.FILE_TYPES.TBTBit.toString())) {
                setType(BinaryToTextPlugin.FILE_TYPES.TBTBit);
            } else if (temp.equalsIgnoreCase(BinaryToTextPlugin.FILE_TYPES.TBTByte.toString())) {
                setType(BinaryToTextPlugin.FILE_TYPES.TBTByte);
            }
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify the file type.");
        }

    }

    private void printUsage() {
        myLogger.info(
                "\nUsage:\n"
                + "TextToBinaryPlugin <options>\n"
                + " -i  Input File Name\n"
                + " -o  Output File Name\n"
                + " -t  File Type (TOPM, TagCounts, TBTBit, TBTByte)\n");
    }

    public void setInput(String filename) {
        myInput = filename;
    }

    public String getInput() {
        return myInput;
    }

    public void setOutput(String filename) {
        myOutput = filename;
    }

    public String getOutput() {
        return myOutput;
    }

    public void setType(BinaryToTextPlugin.FILE_TYPES type) {
        myType = type;
    }

    public BinaryToTextPlugin.FILE_TYPES getType() {
        return myType;
    }

    @Override
    public ImageIcon getIcon() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getButtonName() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getToolTipText() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
