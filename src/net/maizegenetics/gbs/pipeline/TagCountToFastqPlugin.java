/*
 * TagCountToFastqPlugin
 */
package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

import javax.swing.ImageIcon;

import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.gbs.util.BaseEncoder;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.Utils;

import org.apache.log4j.Logger;

/*
 * Converts a TagCounts binary (*.cnt) file (presumably a master tag list) to a fastq file that can be used as input
 * for BWA or bowtie2 (and possibly additional aligners).  The same function can be performed with
 * MergeMultipleTagCountPlugin using the -t option and a single Master Tag List file in the input directory, but
 * having a separate plugin to do this reduces confusion and eliminates the risk of merging the master tag list back on
 * itself.
 *
 * @author glaubitz (jcg233) (modified from MergeMultipleTagCountPlugin)
 */
public class TagCountToFastqPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(TagCountToFastqPlugin.class);
    private ArgsEngine myArgsEngine = null;
    private String inFileName;
    private DataInputStream inStream;
    private String outputFileName = null;
    private DataOutputStream outStream;
    private int nTags, tagLengthInLong;
    private long[] tag = new long[2]; // [indexOfTagLong], emptyLong=0, fileFinish=Long.MAX
    private int tagCount; // tag count
    private byte tagLength;
    private int tagsRead = 0;
    private int minCount = 0;

    public TagCountToFastqPlugin() {
        super(null, false);
    }

    public TagCountToFastqPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        myLogger.info(
                "\n\n\nThe options for the TagCountToFastqPlugin are:\n"
                + "-i Input binary tag count (*.cnt) file\n"
                + "-o Output fastq file to use as input for BWA or bowtie2\n"
                + "-c Minimum count of reads for a tag to be output (default: 1)\n\n\n");
    }

    @Override
    public void setParameters(String[] args) {
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        if (myArgsEngine == null) {
            myArgsEngine = new ArgsEngine();
            myArgsEngine.add("-i", "--input_file", true);
            myArgsEngine.add("-o", "--output_file", true);
            myArgsEngine.add("-c", "--min_count", true);
        }
        myArgsEngine.parse(args);
        if (myArgsEngine.getBoolean("-i")) {
            inFileName = myArgsEngine.getString("-i");
            File inFile = new File(inFileName);
            if (!inFile.isFile()) {
                printUsage();
                throw new IllegalArgumentException("Problem opening the input file: " + myArgsEngine.getString("-i"));
            }
            inFile = null;
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify an input file.\n");
        }
        if (myArgsEngine.getBoolean("-o")) {
            outputFileName = myArgsEngine.getString("-o");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify an output file.");
        }
        if (myArgsEngine.getBoolean("-c")) {
            minCount = Integer.parseInt(myArgsEngine.getString("-c"));
        } else {
            minCount = 1;
        }
    }

    @Override
    public DataSet performFunction(DataSet input) {
        try {
            inStream = new DataInputStream(new BufferedInputStream(new FileInputStream(inFileName), 655360));
            nTags = inStream.readInt();
            tagLengthInLong = inStream.readInt();
            myLogger.info("Opened the input file: " + inFileName + "  nTags=" + nTags);
            //outStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFileName), 655360));
            outStream = Utils.getDataOutputStream(outputFileName, 655360);
            int tagsWritten = 0;
            while (inStream.available()!=0) {
                readNextTag();
                if (tagCount >= minCount) {
                    writeFASTQ();
                    tagsWritten++;
                }
                if (tagsRead % 500000 == 1) {
                    System.out.printf("tagsRead=%d tagsWritten=%d %n", tagsRead, tagsWritten);
                    myLogger.info(BaseEncoder.getSequenceFromLong(tag));
                }
            }
            outStream.flush();
            outStream.close();
            System.out.println("Finished converting binary tag count file to fastq."
                    + "\nTotal number of tags read: " + tagsRead
                    + "\nTotal number of tags written: " + tagsWritten + " (above minCount of " + minCount + ")"
                    + "\nOuput fastq file: " + outputFileName + "\n\n");
        } catch (Exception e) {
            myLogger.info("Catch in reading TagCount file e=" + e);
            e.printStackTrace();
        }
        return null;
    }

    private void readNextTag() {
        try {
            for (int j = 0; j < tagLengthInLong; j++) {
                tag[j] = inStream.readLong();
            }
            tagLength = inStream.readByte();
            tagCount = inStream.readInt();
            tagsRead++;
        } catch (IOException eof) {
            try {
                myLogger.info("Finished reading input file.");
                inStream.close();
                inStream = null;
            } catch (IOException eof2) {
                myLogger.info("Catch closing" + eof2);
                inStream = null;
            }
        }
    }

    private void writeFASTQ() {
        try {
            outStream.writeBytes("@length=" + tagLength + "count=" + tagCount + "\n");   //Length & count header
            String tagSequence = BaseEncoder.getSequenceFromLong(tag);
            tagSequence = tagSequence.substring(0, tagLength);  //Remove any poly-A padding
            outStream.writeBytes(tagSequence + "\n+\n");    //Sequence and "+" symbol
            for (int i = 0; i < tagLength; i++) {
                outStream.writeBytes("f");
            }           //Bogus quality string
            outStream.writeBytes("\n");
        } catch (IOException e) {
            myLogger.info("Catch in writing TagCount file e=" + e);
            e.printStackTrace();
        }
    }

    @Override
    public String getToolTipText() {
        return null;
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return null;
    }
}
