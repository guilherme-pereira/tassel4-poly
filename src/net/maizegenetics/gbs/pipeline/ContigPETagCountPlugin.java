/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;
import java.io.File;
import javax.swing.ImageIcon;
import net.maizegenetics.gbs.tagdist.PETagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.ArgsEngine;
import org.apache.log4j.Logger;

/**
 *
 * @author Fei Lu
 */
public class ContigPETagCountPlugin extends AbstractPlugin {

    static long timePoint1;
    private ArgsEngine engine = null;
    private Logger logger = Logger.getLogger(MergePETagCountPlugin.class);
    private String inputFileS = null;
    private String outputFileS = null;
    
    public ContigPETagCountPlugin() {
        super(null, false);
    }

    public ContigPETagCountPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        logger.info(
                "\n\nUsage is as follows:\n"
                + " -i  input PETagCount filename\n"
                + " -o  output PETagCount filename\n");
    }

    @Override
    public DataSet performFunction(DataSet input) {
        contigPETagCount();
        return null;
    }
    
    /**
     * Make contigs of PETagCounts from forward and backward tags
     */
    public void contigPETagCount () {
        PETagCounts p = new PETagCounts(inputFileS, TagsByTaxa.FilePacking.Bit).getCollapsedPETagCounts();
        p.contigPETags();
        p.writeDistFile(outputFileS, TagsByTaxa.FilePacking.Bit, 0);
    }
    
    @Override
    public void setParameters(String[] args) {
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        if (engine == null) {
            engine = new ArgsEngine();
            engine.add("-i", "--input-file", true);
            engine.add("-o", "--output-file", true);
            engine.parse(args);
        }

        if (engine.getBoolean("-i")) {
            inputFileS = engine.getString("-i");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify the input directory of PETagCounts files");
        }

        if (engine.getBoolean("-o")) {
            outputFileS = engine.getString("-o");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify the merged PETagCounts file");
        }
        
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
