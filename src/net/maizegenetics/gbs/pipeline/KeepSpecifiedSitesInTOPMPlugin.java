/*
 * KeepSpecifiedSitesInTOPMPlugin
 */
package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;
import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Pattern;
import javax.swing.ImageIcon;
import net.maizegenetics.gbs.maps.TOPMInterface;
import net.maizegenetics.gbs.maps.TOPMUtils;
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
public class KeepSpecifiedSitesInTOPMPlugin extends AbstractPlugin {

    private static final Pattern WHITESPACE_PATTERN = Pattern.compile("\\s");
    private final Logger myLogger = Logger.getLogger(KeepSpecifiedSitesInTOPMPlugin.class);
    private static String SITE_LIST_FILENAME_REGEX = "(?i).*\\.txt$";
    private static int PAD_POSITION = 300;
    private ArgsEngine myArgsEngine = null;
    private String[] mySiteListFileNames = null;
    private String myOutputFilename = null;
    private String myOrigFilename = null;
    private TOPMInterface myOrigTOPM = null;
    private int myOrigTagCount = 0;
    private byte[][] myOrigVariantOff = null;
    private byte[][] myOrigVariantDef = null;
    private int[] myNumVariantsKeptPerChrom = new int[20];
    private int[] myTagsWithVariants = new int[20];

    public KeepSpecifiedSitesInTOPMPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    @Override
    public DataSet performFunction(DataSet input) {

        myOrigTOPM = TOPMUtils.readTOPM(myOrigFilename);
        myOrigTagCount = myOrigTOPM.getTagCount();
        myLogger.info("performFunction: Number of Original Tags: " + myOrigTagCount);
        myOrigVariantOff = myOrigTOPM.getVariantOff();
        myOrigVariantDef = myOrigTOPM.getVariantDef();
        myOrigTOPM.clearVariants();

        for (int i = 0; i < mySiteListFileNames.length; i++) {
            if (!mySiteListFileNames[i].equals(myOrigFilename)) {
                processSiteList(mySiteListFileNames[i]);
            }
        }

        for (int x = 0; x < myNumVariantsKeptPerChrom.length; x++) {
            if (myNumVariantsKeptPerChrom[x] != 0) {
                myLogger.info("performFunction: chromosome: " + x + " variants kept: " + myNumVariantsKeptPerChrom[x]);
            }
        }

        for (int x = 0; x < myTagsWithVariants.length; x++) {
            if (myTagsWithVariants[x] != 0) {
                myLogger.info("performFunction: Chromosome: " + x + " Number Tags with Variants Defined: " + myTagsWithVariants[x]);
            }
        }

        TOPMUtils.writeTOPM(myOrigTOPM, myOutputFilename);

        return null;
    }

    private void printUsage() {
        myLogger.info(
                "\n\n\nThe options for the KeepSpecifiedSitesInTOPMPlugin are:\n"
                + "   -input   Input directory containing Site List files\n"
                + "   -orig    Original TOPM\n"
                + "   -result  Output, site-filtered TOPM\n\n\n");
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
                throw new IllegalArgumentException("KeepSpecifiedSitesInTOPMPlugin: setParameters: The input name you supplied is not a directory: " + tempDirectory);
            }
            mySiteListFileNames = DirectoryCrawler.listFileNames(SITE_LIST_FILENAME_REGEX, topmDirectory.getAbsolutePath());
            if (mySiteListFileNames.length == 0 || mySiteListFileNames == null) {
                printUsage();
                throw new IllegalArgumentException("KeepSpecifiedSitesInTOPMPlugin: setParameters: No Site List files in: " + tempDirectory);
            } else {
                myLogger.info("setParameters: Using these Site List files:");
                for (String filename : mySiteListFileNames) {
                    myLogger.info("setParameters: found site list: " + filename);
                }
            }
        }

        myOrigFilename = myArgsEngine.getString("-orig");
        if ((myOrigFilename == null) || (myOrigFilename.length() == 0)) {
            printUsage();
            throw new IllegalArgumentException("KeepSpecifiedSitesInTOPMPlugin: setParameters: Must define original file");
        }
        File origFile = new File(myOrigFilename);
        if (!origFile.exists()) {
            printUsage();
            throw new IllegalArgumentException("KeepSpecifiedSitesInTOPMPlugin: setParameters: The original file doesn't exist: " + myOrigFilename);
        }

        myOutputFilename = myArgsEngine.getString("-result");
        if ((myOutputFilename == null) || (myOutputFilename.length() == 0)) {
            printUsage();
            throw new IllegalArgumentException("KeepSpecifiedSitesInTOPMPlugin: setParameters: Must define result file");
        }
        File outputFile = new File(myOutputFilename);
        if (outputFile.exists()) {
            printUsage();
            throw new IllegalArgumentException("KeepSpecifiedSitesInTOPMPlugin: setParameters: The output file already exists: " + myOutputFilename);
        }

    }

    private void processSiteList(String filename) {

        myLogger.info("processSiteList: " + filename);
        BufferedReader reader = Utils.getBufferedReader(filename);
        try {

            List<Integer> positions = new ArrayList<Integer>();
            String line = reader.readLine();
            String chr = WHITESPACE_PATTERN.split(line)[0];
            while (line != null) {
                String[] tokens = WHITESPACE_PATTERN.split(line);
                if (tokens.length != 2) {
                    throw new IllegalArgumentException("KeepSpecifiedSitesInTOPMPlugin: processSiteList: file not correctly formatted: " + filename);
                }
                if (!tokens[0].equals(chr)) {
                    throw new IllegalArgumentException("KeepSpecifiedSitesInTOPMPlugin: processSiteList: all positions must be from same chromosome: " + filename);
                }
                positions.add(Integer.valueOf(tokens[1]));

                line = reader.readLine();
            }
            reader.close();

            int numPositions = positions.size();
            int[] orderedPositions = new int[numPositions];
            for (int j = 0; j < numPositions; j++) {
                orderedPositions[j] = positions.get(j);
            }
            Arrays.sort(orderedPositions);

            int[] numTagsWithUnfoundSitesAndMaxVariants = new int[numPositions];

            int numVariants = myOrigVariantOff[0].length;
            int chrInt = Integer.valueOf(chr);
            int tagCount = myOrigTOPM.getTagCount();
            for (int i = 0; i < tagCount; i++) {

                if (myOrigTOPM.getChromosome(i) == chrInt) {

                    int startPos = myOrigTOPM.getStartPosition(i);
                    int endPos = myOrigTOPM.getEndPosition(i);
                    byte strand = myOrigTOPM.getStrand(i);

                    if (strand == -1) {

                        if (endPos > startPos) {
                            throw new IllegalStateException("KeepSpecifiedSitesInTOPMPlugin: processSiteList: tag: " + i + " strand: " + strand + " end pos: " + endPos + " is greater than start pos: " + startPos);
                        }
                        int posIndex = Arrays.binarySearch(orderedPositions, Math.max(endPos - PAD_POSITION, 0));
                        if (posIndex < 0) {
                            posIndex = -(posIndex + 1);
                        }
                        int variantAdded = 0;
                        while ((posIndex < numPositions) && (orderedPositions[posIndex] <= startPos + PAD_POSITION)) {
                            boolean found = false;
                            int currentPosition = orderedPositions[posIndex];
                            for (int x = 0; x < numVariants; x++) {
                                if ((myOrigVariantOff[i][x] != Byte.MIN_VALUE) && (myOrigVariantDef[i][x] != Byte.MIN_VALUE)) {
                                    int tagPosition = myOrigVariantOff[i][x] + startPos;
                                    if (tagPosition == currentPosition) {
                                        found = true;
                                        myOrigTOPM.addVariant(i, myOrigVariantOff[i][x], myOrigVariantDef[i][x]);
                                        variantAdded = 1;
                                        if (chrInt < myNumVariantsKeptPerChrom.length) {
                                            myNumVariantsKeptPerChrom[chrInt]++;
                                        }
                                    }
                                }
                            }
                            if (found) {
                                numTagsWithUnfoundSitesAndMaxVariants[posIndex] = -1;
                            } else if ((!found) && (numVariants == myOrigTOPM.getMaxNumVariants()) && (numTagsWithUnfoundSitesAndMaxVariants[posIndex] != -1)
                                    && (currentPosition <= startPos) && (currentPosition >= endPos)) {
                                numTagsWithUnfoundSitesAndMaxVariants[posIndex]++;
                            }
                            posIndex++;
                        }
                        myTagsWithVariants[chrInt] += variantAdded;

                    } else if (strand == 1) {

                        if (startPos > endPos) {
                            throw new IllegalStateException("KeepSpecifiedSitesInTOPMPlugin: processSiteList: tag: " + i + " strand: " + strand + " start pos: " + startPos + " is greater than end pos: " + endPos);
                        }
                        int posIndex = Arrays.binarySearch(orderedPositions, Math.max(startPos - PAD_POSITION, 0));
                        if (posIndex < 0) {
                            posIndex = -(posIndex + 1);
                        }
                        int variantAdded = 0;
                        while ((posIndex < numPositions) && (orderedPositions[posIndex] <= endPos + PAD_POSITION)) {
                            boolean found = false;
                            int currentPosition = orderedPositions[posIndex];
                            for (int x = 0; x < numVariants; x++) {
                                if ((myOrigVariantOff[i][x] != Byte.MIN_VALUE) && (myOrigVariantDef[i][x] != Byte.MIN_VALUE)) {
                                    int tagPosition = myOrigVariantOff[i][x] + startPos;
                                    if (tagPosition == currentPosition) {
                                        found = true;
                                        myOrigTOPM.addVariant(i, myOrigVariantOff[i][x], myOrigVariantDef[i][x]);
                                        variantAdded = 1;
                                        if (chrInt < myNumVariantsKeptPerChrom.length) {
                                            myNumVariantsKeptPerChrom[chrInt]++;
                                        }
                                    }
                                }
                            }
                            if (found) {
                                numTagsWithUnfoundSitesAndMaxVariants[posIndex] = -1;
                            } else if ((!found) && (numVariants == myOrigTOPM.getMaxNumVariants()) && (numTagsWithUnfoundSitesAndMaxVariants[posIndex] != -1)
                                    && (currentPosition >= startPos) && (currentPosition <= endPos)) {
                                numTagsWithUnfoundSitesAndMaxVariants[posIndex]++;
                            }
                            posIndex++;
                        }
                        myTagsWithVariants[chrInt] += variantAdded;

                    } else {
                        throw new IllegalStateException("KeepSpecifiedSitesInTOPMPlugin: processSiteList: tag: " + i + " unknown strand: " + strand);
                    }

                }
            }

            for (int i = 0; i < numPositions; i++) {
                if (numTagsWithUnfoundSitesAndMaxVariants[i] > 0) {
                    myLogger.info("chromosome: " + chrInt + " position: " + orderedPositions[i] + " tags with no variant info: " + numTagsWithUnfoundSitesAndMaxVariants[i]);
                }
            }

        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalStateException("KeepSpecifiedSitesInTOPMPlugin: processSiteList: Problem processing: " + filename);
        } finally {
            try {
                reader.close();
            } catch (Exception e) {
                // do nothing
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
