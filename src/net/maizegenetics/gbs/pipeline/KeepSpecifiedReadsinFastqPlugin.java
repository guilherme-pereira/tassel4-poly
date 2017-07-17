/*
 * KeepSpecifiedReadsinFastqPlugin
 */
package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;
import java.io.BufferedWriter;
import java.io.File;
import javax.swing.ImageIcon;
import net.maizegenetics.gbs.homology.ParseBarcodeRead;
import net.maizegenetics.gbs.homology.ReadBarcodeResult;
import net.maizegenetics.gbs.maps.TOPMInterface;
import net.maizegenetics.gbs.maps.TOPMUtils;
import net.maizegenetics.gbs.util.FastqReader;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

/**
 *
 * @author terry
 */
public class KeepSpecifiedReadsinFastqPlugin extends AbstractPlugin {

    private final Logger myLogger = Logger.getLogger(KeepSpecifiedReadsinFastqPlugin.class);
    private static String FASTQ_FILENAME_REGEX = "(?i).*\\.fq$|.*\\.fq\\.gz$|.*\\.fastq$|.*_fastq\\.txt$|.*_fastq\\.gz$|.*_fastq\\.txt\\.gz$|.*_sequence\\.txt$|.*_sequence\\.txt\\.gz$";
    private ArgsEngine myArgsEngine = null;
    private String[] myInputFastqFileNames = null;
    private String myOutputDir = null;
    private String myTOPMFilename = null;
    private String myKeyFilename = null;
    private String myEnzyme = null;
    private TOPMInterface myTOPM = null;
    private int[] myChrs;
    private int myStartPos;
    private int myEndPos;

    public KeepSpecifiedReadsinFastqPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    @Override
    public DataSet performFunction(DataSet input) {
        BufferedWriter writer = null;
        try {
            myTOPM = TOPMUtils.readTOPM(myTOPMFilename);
            for (int i = 0; i < myInputFastqFileNames.length; i++) {
                myLogger.info("Processing Fastq: " + myInputFastqFileNames[i]);
                String outputFile = myOutputDir + "/" + Utils.getFilename(myInputFastqFileNames[i]);
                myLogger.info("Output Fastq: " + outputFile);
                writer = Utils.getBufferedWriter(outputFile);
                ParseBarcodeRead thePBR = getParseBarcodeRead(myInputFastqFileNames[i]);
                FastqReader reader = new FastqReader(myInputFastqFileNames[i]);
                String[] current = reader.getNext();
                while (current != null) {
                    ReadBarcodeResult rr = thePBR.parseReadIntoTagAndTaxa(current[1], current[3], true, 0);
                    if (rr != null) {
                        int index = myTOPM.getTagIndex(rr.getRead());
                        if (index >= 0) {
                            int chr = myTOPM.getChromosome(index);
                            int start = myTOPM.getStartPosition(index);
                            int end = myTOPM.getEndPosition(index);
                            if (isMatchingChr(chr)) {
                                if (((start >= myStartPos) && (start <= myEndPos))
                                        || ((end >= myStartPos) && (end <= myEndPos))) {
                                    writer.write(current[0]);
                                    writer.write("\n");
                                    writer.write(current[1]);
                                    writer.write("\n");
                                    writer.write(current[2]);
                                    writer.write("\n");
                                    writer.write(current[3]);
                                    writer.write("\n");
                                }
                            }
                        }
                    }
                    current = reader.getNext();
                }
                writer.close();
                reader.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            try {
                writer.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        return null;
    }

    private void printUsage() {
        myLogger.info(
                "\nThe options for the KeepSpecifiedReadsinFastqPlugin:\n"
                + "-input            Input directory containing fastq files\n"
                + "-topm             TOPM\n"
                + "-k                Key File\n"
                + "-e                Enzyme\n"
                + "-output           Output directory for resulting fastq files\n"
                + "-chrs             Comma separated list of chromosomes. (no spaces)\n"
                + "-startPosition    Start Position\n"
                + "-endPosition      End Position\n");
    }

    @Override
    public void setParameters(String[] args) {
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        if (myArgsEngine == null) {
            myArgsEngine = new ArgsEngine();
            myArgsEngine.add("-input", "-input", true);
            myArgsEngine.add("-topm", "-topm", true);
            myArgsEngine.add("-k", "-k", true);
            myArgsEngine.add("-e", "-e", true);
            myArgsEngine.add("-output", "-output", true);
            myArgsEngine.add("-chrs", "-chrs", true);
            myArgsEngine.add("-startPosition", "-startPosition", true);
            myArgsEngine.add("-endPosition", "-endPosition", true);
        }
        myArgsEngine.parse(args);

        String tempDirectory = myArgsEngine.getString("-input");
        if ((tempDirectory != null) && tempDirectory.length() != 0) {
            File inputDirectory = new File(tempDirectory);
            if (!inputDirectory.isDirectory()) {
                printUsage();
                throw new IllegalArgumentException("KeepSpecifiedReadsinFastqPlugin: setParameters: The input name you supplied is not a directory: " + tempDirectory);
            }
            myInputFastqFileNames = DirectoryCrawler.listFileNames(FASTQ_FILENAME_REGEX, inputDirectory.getAbsolutePath());
            if (myInputFastqFileNames.length == 0 || myInputFastqFileNames == null) {
                printUsage();
                throw new IllegalArgumentException("KeepSpecifiedReadsinFastqPlugin: setParameters: No fastq files in: " + tempDirectory);
            } else {
                myLogger.info("setParameters: Using these fastq files:");
                for (String filename : myInputFastqFileNames) {
                    myLogger.info("setParameters: found fastq: " + filename);
                }
            }
        } else {
            printUsage();
            throw new IllegalArgumentException("KeepSpecifiedReadsinFastqPlugin: setParameters: No input directory specified.");
        }

        myTOPMFilename = myArgsEngine.getString("-topm");
        if ((myTOPMFilename == null) || (myTOPMFilename.length() == 0)) {
            printUsage();
            throw new IllegalArgumentException("KeepSpecifiedReadsinFastqPlugin: setParameters: Must define TOPM");
        }
        File tempFile = new File(myTOPMFilename);
        if (!tempFile.exists()) {
            printUsage();
            throw new IllegalArgumentException("KeepSpecifiedReadsinFastqPlugin: setParameters: The TOPM file doesn't exist: " + myTOPMFilename);
        }

        myKeyFilename = myArgsEngine.getString("-k");
        if ((myKeyFilename == null) || (myKeyFilename.length() == 0)) {
            printUsage();
            throw new IllegalArgumentException("KeepSpecifiedReadsinFastqPlugin: setParameters: Must define Key File");
        }
        tempFile = new File(myKeyFilename);
        if (!tempFile.exists()) {
            printUsage();
            throw new IllegalArgumentException("KeepSpecifiedReadsinFastqPlugin: setParameters: The Key file doesn't exist: " + myKeyFilename);
        }

        myEnzyme = myArgsEngine.getString("-e");
        if ((myEnzyme == null) || (myEnzyme.length() == 0)) {
            printUsage();
            throw new IllegalArgumentException("KeepSpecifiedReadsinFastqPlugin: setParameters: Must define TOPM");
        }

        myOutputDir = myArgsEngine.getString("-output");
        if ((myOutputDir != null) && myOutputDir.length() != 0) {
            File outputDirectory = new File(myOutputDir);
            if (!outputDirectory.isDirectory()) {
                printUsage();
                throw new IllegalArgumentException("KeepSpecifiedReadsinFastqPlugin: setParameters: The output name you supplied is not a directory: " + myOutputDir);
            }
        } else {
            printUsage();
            throw new IllegalArgumentException("KeepSpecifiedReadsinFastqPlugin: setParameters: No output directory specified.");
        }

        if (myOutputDir.equals(tempDirectory)) {
            throw new IllegalArgumentException("KeepSpecifiedReadsinFastqPlugin: setParameters: Output Directory should be different from Input Directory");
        }

        String tempChrs = myArgsEngine.getString("-chrs");
        if ((tempChrs == null) || tempChrs.length() == 0) {
            throw new IllegalArgumentException("KeepSpecifiedReadsinFastqPlugin: setParameters: No chromosomes specified.");
        }
        String[] tokenChrs = tempChrs.split(",");
        myChrs = new int[tokenChrs.length];
        for (int i = 0; i < tokenChrs.length; i++) {
            try {
                myChrs[i] = Integer.parseInt(tokenChrs[i]);
            } catch (NumberFormatException e) {
                throw new IllegalArgumentException("KeepSpecifiedReadsinFastqPlugin: setParameters: chromosome not a number: " + tokenChrs[i]);
            }
        }

        String tempStart = myArgsEngine.getString("-startPosition");
        if ((tempStart == null) || tempStart.length() == 0) {
            throw new IllegalArgumentException("KeepSpecifiedReadsinFastqPlugin: setParameters: No Start Position specified.");
        }
        try {
            myStartPos = Integer.parseInt(tempStart);
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("KeepSpecifiedReadsinFastqPlugin: setParameters: Start Position not a number: " + tempStart);
        }

        String tempEnd = myArgsEngine.getString("-endPosition");
        if ((tempEnd == null) || tempEnd.length() == 0) {
            throw new IllegalArgumentException("KeepSpecifiedReadsinFastqPlugin: setParameters: No End Position specified.");
        }
        try {
            myEndPos = Integer.parseInt(tempEnd);
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("KeepSpecifiedReadsinFastqPlugin: setParameters: End Position not a number: " + tempEnd);
        }

        if (myStartPos > myEndPos) {
            throw new IllegalArgumentException("KeepSpecifiedReadsinFastqPlugin: setParameters: Start Position: " + myStartPos + " can't be larger than End Position: " + myEndPos);
        }

    }

    private boolean isMatchingChr(int chr) {
        for (int i = 0; i < myChrs.length; i++) {
            if (myChrs[i] == chr) {
                return true;
            }
        }
        return false;
    }

    private ParseBarcodeRead getParseBarcodeRead(String fastqFile) {

        String[] np = Utils.getFilename(fastqFile).split("_");

        ParseBarcodeRead result;
        if (np.length == 3) {
            result = new ParseBarcodeRead(myKeyFilename, myEnzyme, np[0], np[1]);
        } else if (np.length == 5) {
            result = new ParseBarcodeRead(myKeyFilename, myEnzyme, np[1], np[3]);
        } else if (np.length == 4) {
            result = new ParseBarcodeRead(myKeyFilename, myEnzyme, np[0], np[2]);
        } else if (np.length == 6) {
            result = new ParseBarcodeRead(myKeyFilename, myEnzyme, np[1], np[3]);
        } else {
            myLogger.error("Error in parsing file name:");
            myLogger.error("   The filename does not contain either 3 or 5 underscore-delimited values.");
            myLogger.error("   Expect: flowcell_lane_fastq.txt OR code_flowcell_s_lane_fastq.txt");
            myLogger.error("   Filename: " + fastqFile);
            return null;
        }

        myLogger.info("Total barcodes found in lane:" + result.getBarCodeCount());
        return result;

    }

    @Override
    public ImageIcon getIcon() {
        return null;
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
