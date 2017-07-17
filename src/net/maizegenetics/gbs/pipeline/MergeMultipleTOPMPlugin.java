/*
 * MergeMultipleTOPMPlugin
 */
package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;
import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Arrays;
import javax.swing.ImageIcon;
import net.maizegenetics.gbs.maps.TOPMInterface;
import net.maizegenetics.gbs.maps.TOPMUtils;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.util.DirectoryCrawler;
import org.apache.log4j.Logger;

/**
 *
 * @author terry
 */
public class MergeMultipleTOPMPlugin extends AbstractPlugin {

    private final Logger myLogger = Logger.getLogger(MergeMultipleTOPMPlugin.class);
    private static String TOPM_FILENAME_REGEX = "(?i).*\\.topm$|.*\\.topm\\.bin";
    private ArgsEngine myArgsEngine = null;
    private String[] myTOPMFileNames = null;
    private String myOutputFilename = null;
    private String myOrigFilename = null;
    private TOPMInterface myOrigTOPM = null;
    private int myOrigTagCount = 0;
    private byte[][] myOrigVariantOff = null;
    private byte[][] myOrigVariantDef = null;
    private boolean[] myChangedRows = null;
    private int[] myChromosomeChangedCounts = new int[20];

    public MergeMultipleTOPMPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    @Override
    public DataSet performFunction(DataSet input) {

        myOrigTOPM = TOPMUtils.readTOPM(myOrigFilename);
        myOrigTagCount = myOrigTOPM.getTagCount();
        myLogger.info("performFunction: Number of Original Tags: " + myOrigTagCount);
        myOrigVariantOff = myOrigTOPM.getVariantOff();
        myOrigVariantDef = myOrigTOPM.getVariantDef();
        myChangedRows = new boolean[myOrigTagCount];
        Arrays.fill(myChangedRows, false);

        for (int i = 0; i < myTOPMFileNames.length; i++) {
            if (!myTOPMFileNames[i].equals(myOrigFilename)) {
                processTOPM(myTOPMFileNames[i]);
            }
        }

        for (int x = 0; x < myChromosomeChangedCounts.length; x++) {
            if (myChromosomeChangedCounts[x] != 0) {
                myLogger.info("performFunction: chromosome: " + x + " changed: " + myChromosomeChangedCounts[x]);
            }
        }

        TOPMUtils.writeTOPM(myOrigTOPM, myOutputFilename);

        return null;
    }

    private void printUsage() {
        myLogger.info(
                "\nThe options for the MergeMultipleTOPMPlugin:\n"
                + "-input  Input directory containing TOPM files\n"
                + "-orig Original TOPM"
                + "-result  TOPM Output Filename\n");
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
            myArgsEngine.add("-orig", "-orig", true);
            myArgsEngine.add("-result", "-result", true);
        }
        myArgsEngine.parse(args);

        String tempDirectory = myArgsEngine.getString("-input");
        if ((tempDirectory != null) && tempDirectory.length() != 0) {
            File topmDirectory = new File(tempDirectory);
            if (!topmDirectory.isDirectory()) {
                printUsage();
                throw new IllegalArgumentException("MergeMultipleTOPMPlugin: setParameters: The input name you supplied is not a directory: " + tempDirectory);
            }
            myTOPMFileNames = DirectoryCrawler.listFileNames(TOPM_FILENAME_REGEX, topmDirectory.getAbsolutePath());
            if (myTOPMFileNames.length == 0 || myTOPMFileNames == null) {
                printUsage();
                throw new IllegalArgumentException("MergeMultipleTOPMPlugin: setParameters: No TOPM files in: " + tempDirectory);
            } else {
                myLogger.info("setParameters: Using these TOPM files:");
                for (String filename : myTOPMFileNames) {
                    myLogger.info(filename);
                }
            }
        }

        myOrigFilename = myArgsEngine.getString("-orig");
        if ((myOrigFilename == null) || (myOrigFilename.length() == 0)) {
            printUsage();
            throw new IllegalArgumentException("MergeMultipleTOPMPlugin: setParameters: Must define original file");
        }
        File origFile = new File(myOrigFilename);
        if (!origFile.exists()) {
            printUsage();
            throw new IllegalArgumentException("MergeMultipleTOPMPlugin: setParameters: The original file doesn't exist: " + myOrigFilename);
        }

        myOutputFilename = myArgsEngine.getString("-result");
        if ((myOutputFilename == null) || (myOutputFilename.length() == 0)) {
            printUsage();
            throw new IllegalArgumentException("MergeMultipleTOPMPlugin: setParameters: Must define result file");
        }
        File outputFile = new File(myOutputFilename);
        if (outputFile.exists()) {
            printUsage();
            throw new IllegalArgumentException("MergeMultipleTOPMPlugin: setParameters: The output file already exists: " + myOutputFilename);
        }

    }

    private void processTOPM(String filename) {

        myLogger.info("processTOPM: " + filename);
        DataInputStream dis = null;
        int tagsInput = 0;
        try {
            dis = new DataInputStream(new BufferedInputStream(new FileInputStream(filename), 65536));
            int tagNum = dis.readInt();
            int tagLengthInLong = dis.readInt();
            int maxVariants = dis.readInt();
            for (int row = 0; row < tagNum; row++) {
                tagsInput++;
                processTag(dis, row, tagLengthInLong, maxVariants);
                if (row % 1000000 == 0) {
                    myLogger.info("processTOPM: Tags Read: " + row);
                }
            }
            myLogger.info("processTOPM: Number of Tags: " + tagsInput);
        } catch (Exception e) {
            myLogger.error("processTOPM: Error Reading Tag: " + tagsInput);
            e.printStackTrace();
            throw new IllegalStateException("MergeMultipleTOPMPlugin: processTOPM: Problem processing: " + filename);
        } finally {
            try {
                dis.close();
            } catch (Exception e) {
                // do nothing
            }
        }

    }

    private void processTag(DataInputStream dis, int row, int tagLengthInLong, int maxVariants) throws IOException {

        long[] tags = new long[tagLengthInLong];
        for (int j = 0; j < tagLengthInLong; j++) {
            tags[j] = dis.readLong();
        }
        byte tagLength = dis.readByte();
        byte multimaps = dis.readByte();
        int chromosome = dis.readInt();
        byte strand = dis.readByte();
        int startPosition = dis.readInt();
        int endPosition = dis.readInt();
        byte divergence = dis.readByte();
        byte[] variantPosOff = new byte[maxVariants];
        byte[] variantDef = new byte[maxVariants];
        for (int j = 0; j < maxVariants; j++) {
            variantPosOff[j] = dis.readByte();
            variantDef[j] = dis.readByte();
        }
        byte dcoP = dis.readByte();
        byte mapP = dis.readByte();

        if (myOrigTOPM.getTagLength(row) != tagLength) {
            myLogger.error("processTag: " + row + " Tag Length: " + tagLength + " doesn't match Original: " + myOrigTOPM.getTagLength(row));
        }

        if (myOrigTOPM.getMultiMaps(row) != multimaps) {
            myLogger.error("processTag: " + row + " Multi Maps: " + multimaps + " doesn't match Original: " + myOrigTOPM.getMultiMaps(row));
        }

        if (myOrigTOPM.getChromosome(row) != chromosome) {
            myLogger.error("processTag: " + row + " Chromosome: " + chromosome + " doesn't match Original: " + myOrigTOPM.getChromosome(row));
        }

        if (myOrigTOPM.getStrand(row) != strand) {
            myLogger.error("processTag: " + row + " Strand: " + strand + " doesn't match Original: " + myOrigTOPM.getStrand(row));
        }

        if (myOrigTOPM.getStartPosition(row) != startPosition) {
            myLogger.error("processTag: " + row + " Start Position: " + startPosition + " doesn't match Original: " + myOrigTOPM.getStartPosition(row));
        }

        if (myOrigTOPM.getEndPosition(row) != endPosition) {
            myLogger.error("processTag: " + row + " End Position: " + endPosition + " doesn't match Original: " + myOrigTOPM.getEndPosition(row));
        }

        if (myOrigTOPM.getDivergence(row) != divergence) {
            myLogger.error("processTag: " + row + " Divergence: " + divergence + " doesn't match Original: " + myOrigTOPM.getDivergence(row));
        }

        if (myOrigTOPM.getDcoP(row) != dcoP) {
            myLogger.error("processTag: " + row + " DcoP: " + dcoP + " doesn't match Original: " + myOrigTOPM.getDcoP(row));
        }

        if (myOrigTOPM.getMapP(row) != mapP) {
            myLogger.error("processTag: " + row + " MapP: " + mapP + " doesn't match Original: " + myOrigTOPM.getMapP(row));
        }

        boolean variantsEqual = true;

        for (int i = 0; i < maxVariants; i++) {
            if (myOrigVariantDef[row][i] != variantDef[i]) {
                variantsEqual = false;
                break;
            }
            if (myOrigVariantOff[row][i] != variantPosOff[i]) {
                variantsEqual = false;
                break;
            }
        }

        if (!variantsEqual) {

            if (myChangedRows[row]) {
                throw new IllegalStateException("MergeMultipleTOPMPlugin: processTag: " + row + " has already been merged.");
            } else {
                myChangedRows[row] = true;
                if (chromosome < myChromosomeChangedCounts.length) {
                    myChromosomeChangedCounts[chromosome]++;
                }
            }

            for (int i = 0; i < maxVariants; i++) {
                myOrigTOPM.setVariantDef(row, i, variantDef[i]);
                myOrigTOPM.setVariantPosOff(row, i, variantPosOff[i]);
            }

        }

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
