/*
 * BinaryToText
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
public class BinaryToTextPlugin extends AbstractPlugin {

    private Logger myLogger = Logger.getLogger(BinaryToTextPlugin.class);

    public static enum FILE_TYPES {

        TOPM, TagCounts, TBTBit, TBTByte
    };
    private ArgsEngine myEngine = null;
    private String myInput;
    private String myOutput;
    private FILE_TYPES myType;

    public BinaryToTextPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    @Override
    public DataSet performFunction(DataSet input) {

        switch (getType()) {
            case TOPM:
                TagsOnPhysicalMap topm = new TagsOnPhysicalMap(myInput, true);
                topm.writeTextFile(new File(myOutput));
                break;
            case TagCounts:
                TagCounts tc = new TagCounts(myInput, FilePacking.Bit);
                tc.writeTagCountFile(myOutput, FilePacking.Text, 0);
                break;
            case TBTBit:
                TagsByTaxaBit tbtbit = new TagsByTaxaBit(myInput, FilePacking.Bit);
                tbtbit.writeDistFile(new File(myOutput), FilePacking.Text, 0);
                break;
            case TBTByte:
                TagsByTaxaByte tbtbyte = new TagsByTaxaByte(myInput, FilePacking.Byte);
                tbtbyte.writeDistFile(new File(myOutput), FilePacking.Text, 0);
                break;
        }

        return null;

    }

    @Override
    public void setParameters(String[] args) {
        if (args == null || args.length == 0) {
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
            if (temp.equalsIgnoreCase(FILE_TYPES.TOPM.toString())) {
                setType(FILE_TYPES.TOPM);
            } else if (temp.equalsIgnoreCase(FILE_TYPES.TagCounts.toString())) {
                setType(FILE_TYPES.TagCounts);
            } else if (temp.equalsIgnoreCase(FILE_TYPES.TBTBit.toString())) {
                setType(FILE_TYPES.TBTBit);
            } else if (temp.equalsIgnoreCase(FILE_TYPES.TBTByte.toString())) {
                setType(FILE_TYPES.TBTByte);
            }
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify the file type.");
        }

    }

    private void printUsage() {
        myLogger.info(
                "\nUsage:\n"
                + "BinaryToTextPlugin <options>\n"
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

    public void setType(FILE_TYPES type) {
        myType = type;
    }

    public FILE_TYPES getType() {
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
