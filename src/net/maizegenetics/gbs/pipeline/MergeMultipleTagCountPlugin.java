/*
 * MergeMultipleTagCountPlugin
 */
package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

import javax.swing.ImageIcon;

import net.maizegenetics.gbs.tagdist.AbstractTags;
import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.gbs.util.BaseEncoder;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.DirectoryCrawler;

import org.apache.log4j.Logger;

/*
 * Implements an external mergesort to combine multiple tag-count files.
 * @author edbuckler
 */
public class MergeMultipleTagCountPlugin extends AbstractPlugin {

    private static boolean textOutput = false; //Specifies text output format
    private static long[][] ctags;  //array of current tags, [indexOfFile][indexOfTagLong], emptyLong=0, fileFinish=Long.MAX
    private static int[] ctagCnt; //array of tag counts [indexOfFile]
    private static byte[] ctagLength;  //array of tag length [indexOfFile]
    private static int[] chunkTagSizes;  // array of number of tags in each file  [indexOfFile]
//    private int[] ctagRemaining; //How many tags still remain to be read from each chunk file
    private static int tagLengthInLong = 2;
    private static int tagsRead = 0, outCnt = 0;
    private static int rwOpen = 0;
    private static DataInputStream[] rw;
    private static DataOutputStream outStream;
    private static ArgsEngine myArgsEngine = null;
    static File inputDirectory = null;
    static String[] inputFileNames = null;
    static String outputFileName = null;
    static int minCount = 0;
    private static final Logger myLogger = Logger.getLogger(MergeMultipleTagCountPlugin.class);

    public MergeMultipleTagCountPlugin() {
        super(null, false);
    }

    public MergeMultipleTagCountPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    @Override
    public void setParameters(String[] args) {
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        if (myArgsEngine == null) {
            myArgsEngine = new ArgsEngine();
            myArgsEngine.add("-i", "--input-directory", true);
            myArgsEngine.add("-o", "--output_file", true);
            myArgsEngine.add("-c", "--min-count", true);
            myArgsEngine.add("-t", "--fastq");
        }

        myArgsEngine.parse(args);
        if (myArgsEngine.getBoolean("-o")) {
            outputFileName = myArgsEngine.getString("-o");
            if (new File(outputFileName).isDirectory()) {
                printUsage();
                throw new IllegalArgumentException("The output name you supplied is a directory, not a file.");
            }
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify an output file.");
        }

        if (myArgsEngine.getBoolean("-i")) {
            inputDirectory = new File(myArgsEngine.getString("-i"));
            if (!inputDirectory.isDirectory()) {
                printUsage();
                throw new IllegalArgumentException("The input name you supplied is not a directory.");
            }
            inputFileNames = DirectoryCrawler.listFileNames(".*\\.cnt", inputDirectory.getAbsolutePath());
            if (inputFileNames == null || inputFileNames.length == 0) {
                printUsage();
                throw new IllegalArgumentException("Couldn't find any files ending in \".cnt\" in the directory you specified.");
            }
            myLogger.info("Merging the following .cnt files...");
            for (String filename : inputFileNames) {
                myLogger.info(filename);
            }
            myLogger.info("...to \"" + outputFileName + "\".");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify an input file.\n");

        }
        if (myArgsEngine.getBoolean("-t")) {
            textOutput = true;
        }

        if (myArgsEngine.getBoolean("-c")) {
            minCount = Integer.parseInt(myArgsEngine.getString("-c"));
        } else {
            minCount = 1;
        }

    }

    private void printUsage() {
        myLogger.info(
                "\nUsage is as follows:\n"
                + "-i  Input directory containing .cnt files\n"
                + "-o  Output file name\n"
                + "-c Minimum count of reads to be output (default 1)\n"
                + "-t Specifies that reads should be output in FASTQ text format.");
    }

    @Override
    public DataSet performFunction(DataSet input) {
        mergeChunks(inputFileNames, inputDirectory, outputFileName, minCount);
        return null;
    }

    public static void mergeChunks(String[] chunkFileNames, File inputDirectory, String outputFileName, int minCount) {
        rw = new DataInputStream[chunkFileNames.length];
        chunkTagSizes = new int[rw.length];
        ctagCnt = new int[rw.length];
        ctagLength = new byte[rw.length];
        rwOpen = rw.length;
        try {
            for (int f = 0; f < rw.length; f++) {
                String infile = chunkFileNames[f];
                rw[f] = new DataInputStream(new BufferedInputStream(new FileInputStream(infile), 4000000));
                chunkTagSizes[f] = rw[f].readInt();
                tagLengthInLong = rw[f].readInt();
                myLogger.info("Opened :" + infile + " tags=" + chunkTagSizes[f]);
            }
            outStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFileName + ".fq"), 655360));
            ctags = new long[rw.length][tagLengthInLong];
            outStack("B:");
            int t = 0;
            while (rwOpen > 0) {
                long[] minTag = updateCurrentTags();
                //     myLogger.info(BaseEncoder.getSequenceFromLong(minTag));
                //     outStack("U:"+t);
                if (textOutput) {
                    writeFASTQ(minTag, minCount);
                } else {
                    writeTags(minTag, minCount);
                }
                //     outStack("A:"+t);
                if (t % 1000000 == 0) {
                    System.out.printf("t=%d tagsRead=%d outCnt=%d rwOpen=%d %n", t, tagsRead, outCnt, rwOpen);
                    outStack("A:" + t);
                    myLogger.info(BaseEncoder.getSequenceFromLong(minTag));
                }
                t++;
            }
            outStream.flush();
            outStream.close();
        } catch (Exception e) {
            myLogger.info("Catch in reading TagCount file e=" + e);
            e.printStackTrace();
        }
        //Binary files need a header
        if (!textOutput) {
            prependHeader(outputFileName, outCnt, tagLengthInLong);
        }
    }

    private static void outStack(String prefix) {
        System.out.print(prefix + ":");
        for (int f = 0; f < rw.length; f++) {
            System.out.print(ctags[f][0] + ":");
        }
        myLogger.info("");
    }

    private static void writeTags(long[] writeTag, int minCount) {
        int count = 0;
        int tagLength = -1;

        //Loop through input buffers, compare current records, write smallest to output buffer
        for (int f = 0; f < rw.length; f++) {
            if (AbstractTags.compareTags(ctags[f], writeTag) == 0) {
                count += ctagCnt[f];
                tagLength = ctagLength[f];
                for (int j = 0; j < tagLengthInLong; j++) {
                    ctags[f][j] = 0;
                }
            }
        }
        if (count >= minCount) {
            outCnt++;
        } else {
            return;
        }
        //write to file here.
        try {
            long test1;
            for (int i = 0; i < tagLengthInLong; i++) {
                test1 = writeTag[i];
                outStream.writeLong(test1);
            }
            byte test2 = (byte) tagLength;
            int test3 = count;
            outStream.writeByte(test2);
            outStream.writeInt(test3);
            if (count == 0) {
                myLogger.info("");
            }
        } catch (IOException e) {
            myLogger.info("Catch in writing TagCount file e=" + e);
            e.printStackTrace();
        }
    }

    private static void writeFASTQ(long[] writeTag, int minCount) {
        int count = 0;
        int tagLength = -1;
        String tagSequence = "";

        //Loop through input buffers, compare current records, write smallest to output buffer
        for (int f = 0; f < rw.length; f++) {
            if (AbstractTags.compareTags(ctags[f], writeTag) == 0) {
                count += ctagCnt[f];
                tagLength = ctagLength[f];
                for (int j = 0; j < tagLengthInLong; j++) {
                    ctags[f][j] = 0;
                }
            }
        }
        if (count >= minCount) {
            outCnt++;
        } else {
            return;
        }
        //write to file here.
        try {
            outStream.writeBytes("@length=" + tagLength + "count=" + count + "\n");   //Length & count header
            tagSequence = BaseEncoder.getSequenceFromLong(writeTag);
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

    private static void prependHeader(String fileName, int tagCount, int tagLengthInLong) {
        myLogger.info("Adding header to " + fileName + ".");
        File inputFile = new File(fileName + ".fq");
        File outputFile = new File(fileName);
        long[] testLong = new long[2];
        int tagsWritten = 0;
        try {
            DataInputStream input = new DataInputStream(new BufferedInputStream(new FileInputStream(inputFile), 4000000));
            DataOutputStream output = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile), 65536));

            output.writeInt(tagCount);
            output.writeInt(tagLengthInLong);

            for (int i = 0; i < tagCount; i++) {
                for (int j = 0; j < tagLengthInLong; j++) {
                    output.writeLong(input.readLong());
                }
                output.writeByte(input.readByte());
                output.writeInt(input.readInt());
                tagsWritten++;
                if (tagsWritten % 1000000 == 0) {
                    myLogger.info("Wrote " + tagsWritten + " records.");
                }
            }
            input.close();
            output.close();

            if (!inputFile.delete()) {
                myLogger.info("WARNING: Failure to delete file:\n\t" + inputFile.getCanonicalPath());  //Delete old file
            }
        } catch (Exception e) {
            myLogger.info("Caught exception while prepending read count to read count file: " + e);
            e.printStackTrace();
        }

    }

    private static long[] updateCurrentTags() {
        long[] minTag = new long[tagLengthInLong];
        minTag[0] = Long.MAX_VALUE;
        for (int f = 0; f < rw.length; f++) {
            if (ctags[f][0] == 0) {
                readNextTag(f);
            }
            if (AbstractTags.compareTags(ctags[f], minTag) < 0) {
                minTag = ctags[f].clone();
            }
        }
        return minTag;
    }

    private static void readNextTag(int f) {
        if (rw[f] == null) {
            return;
        }
        try {
            for (int j = 0; j < tagLengthInLong; j++) {
                ctags[f][j] = rw[f].readLong();
            }
            ctagLength[f] = rw[f].readByte();
            ctagCnt[f] = rw[f].readInt();
            tagsRead++;
        } catch (IOException eof) {
            try {
                myLogger.info("Finished reading file " + f + ".");
                rw[f].close();
                rw[f] = null;
                for (int i = 0; i < tagLengthInLong; i++) {

                    ctags[f][i] = Long.MAX_VALUE;
                }
                rwOpen--;
            } catch (IOException eof2) {
                myLogger.info("Catch closing" + eof2);
                rw[f] = null;
            }
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
