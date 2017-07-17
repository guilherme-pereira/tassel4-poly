/*
 * QseqToTBTPlugin
 */
package net.maizegenetics.gbs.pipeline;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;

import net.maizegenetics.util.MultiMemberGZIPInputStream;
import net.maizegenetics.gbs.homology.ParseBarcodeRead;
import net.maizegenetics.gbs.homology.ReadBarcodeResult;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.Tags;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.tagdist.TagsByTaxaBit;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByte;
import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;

import java.awt.Frame;

import javax.swing.ImageIcon;

import org.apache.log4j.Logger;

/**
 * This pipeline converts a series of qseq files to TagsByTaxa files (one per qseq file).
 * It requires a list of existing tags (Tags object), which may come from a TagCounts file or TOPM file.
 *
 * @author james
 */
public class QseqToTBTPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(QseqToTBTPlugin.class);
    private ArgsEngine myArgsEngine = null;
    private String[] myQseqFileS = null;
    private String myKeyFile = null;
    private String myEnzyme = null;
    private String myOutputDir = null;
    private int myMinCount = 1;
    private Tags myMasterTags = null;
    private boolean useTBTByte = false;

    public QseqToTBTPlugin() {
        super(null, false);
    }

    public QseqToTBTPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    @Override
    public DataSet performFunction(DataSet input) {
        matchTagsToTaxa(myQseqFileS, myKeyFile, myEnzyme, myMasterTags, myOutputDir, myMinCount, useTBTByte);
        return null;
    }

    private void printUsage() {
        myLogger.info(
                "\nUsage is as follows:\n"
                + "-i  Input directory containing .qseq files\n"
                + "-k  Barcode key file\n"
                + "-e  Enzyme used to create the GBS library, if it differs from the one listed in the key file.\n"
                + "-o  Output directory\n"
                + "-c  Minimum taxa count within a qseq file for a tag to be output (default 1)\n" // Nb: using TagsByTaxaBit, so max count PER TAXON = 1
                + "-y  Output to tagsByTaxaByte (tag counts per taxon from 0 to 127) instead of tagsByTaxaBit (0 or 1)\n"
                + "One of either:\n"
                + "    -t  Tag count file, OR A\n"
                + "    -m  Physical map file containing alignments\n");
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
            myArgsEngine.add("-k", "--key-file", true);
            myArgsEngine.add("-e", "--enzyme", true);
            myArgsEngine.add("-o", "--output-directory", true);
            myArgsEngine.add("-c", "--min-count", true);
            myArgsEngine.add("-y", "--TBTbyte", false);
            myArgsEngine.add("-t", "--tag-count", true);
            myArgsEngine.add("-m", "--physical-map", true);
        }
        myArgsEngine.parse(args);

        String tempDirectory = myArgsEngine.getString("-i");
        if (tempDirectory != null) {
            File qseqDirectory = new File(tempDirectory);
            if (!qseqDirectory.isDirectory()) {
                printUsage();
                throw new IllegalArgumentException("setParameters: The input name you supplied is not a directory: " + tempDirectory);
            }
            myQseqFileS = DirectoryCrawler.listFileNames(".*_qseq\\.txt$|.*_qseq\\.txt\\.gz$", qseqDirectory.getAbsolutePath());
            if (myQseqFileS.length == 0 || myQseqFileS == null) {
                printUsage();
                throw new IllegalArgumentException("Couldn't find any files that end with \"_qseq.txt\" or \"_qseq.txt.gz\" in the supplied directory: " + tempDirectory);
            } else {
                myLogger.info("QseqToTBTPlugin: setParameters: Using the following .qseq files:");
                for (String filename : myQseqFileS) {
                    myLogger.info(filename);
                }
            }
        }
        if (myArgsEngine.getBoolean("-k")) {
            myKeyFile = myArgsEngine.getString("-k");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify a key file (option -k).");
        }
        if (myArgsEngine.getBoolean("-e")) {
            myEnzyme = myArgsEngine.getString("-e");
        } else {
            System.out.println("No enzyme specified.  Using enzyme listed in key file.");
        }
        if (myArgsEngine.getBoolean("-o")) {
            myOutputDir = myArgsEngine.getString("-o");
            File outDirectory = new File(myOutputDir);
            if (!outDirectory.isDirectory()) {
                printUsage();
                throw new IllegalArgumentException("The output name you supplied (option -o) is not a directory: " + myOutputDir);
            }
            outDirectory = null;
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify an output directory (option -o).");
        }
        if (myArgsEngine.getBoolean("-c")) {
            myMinCount = Integer.parseInt(myArgsEngine.getString("-c"));
        } else {
            myMinCount = 1;
        }
        if (myArgsEngine.getBoolean("-y")) {
            useTBTByte = true;
        }

        // Create Tags object from tag count file with option -t, or from TOPM file with option -m
        if (myArgsEngine.getBoolean("-t")) {
            if (myArgsEngine.getBoolean("-m")) {
                printUsage();
                throw new IllegalArgumentException("Options -t and -m are mutually exclusive.");
            }
            myMasterTags = new TagCounts(myArgsEngine.getString("-t"), FilePacking.Bit);
        } else if (myArgsEngine.getBoolean("-m")) {
            if (myArgsEngine.getBoolean("-t")) {
                printUsage();
                throw new IllegalArgumentException("Options -t and -m are mutually exclusive.");
            }
            myMasterTags = new TagsOnPhysicalMap(myArgsEngine.getString("-m"), true);
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify a tagCounts file (-t) *OR* a TagsOnPhysicalMap file (-m)");
        }
    }

    /**
     * Uses an existing Tags object to create one TagsByTaxa file for each qseq file in the input directory.
     *
     * Output TBT files written to the outputDir, using qseq file names with extension changed to .tbt.bin (or .tbt.txt)
     *
     * @param qseqFileS      Array of qseq file names (Illumina-created files with raw read sequence, quality score, machine name, etc.)
     * @param keyFileS       A key file (list of taxa by barcode, lane & flow cell, including plate maps)
     * @param enzyme         The enzyme used to make the library (currently ApeKI or PstI)
     * @param theMasterTags  A Tags object: list of tags to be included in the final TBT
     * @param outputDir      String containing the path of the output directory to contain tags-by-taxa files
     * @param minCount       The minimum number of times a tag must show up in a qseq file before it is included in the corresponding TBT file
     */
    public static void matchTagsToTaxa(String[] qseqFileS, String keyFileS, String enzyme, Tags theMasterTags, String outputDir, int minCount, boolean useTBTByte) {
        for (int laneNum = 0; laneNum < qseqFileS.length; laneNum++) {
            System.out.println("\nWorking on qseq file: " + qseqFileS[laneNum]);
            TagsByTaxa theTBT = null;
            System.gc();

            File outfile;
            FilePacking outFormat = useTBTByte ? FilePacking.Byte : FilePacking.Bit;
            String outFileS = outputDir + qseqFileS[laneNum].substring(qseqFileS[laneNum].lastIndexOf(File.separator));
            String replaceS = (outFormat == FilePacking.Text) ? ".tbt.txt" : ((outFormat == FilePacking.Byte) ? ".tbt.byte" : ".tbt.bin");
            outfile = new File(outFileS.replaceAll("_qseq\\.txt$|_qseq\\.txt\\.gz$", replaceS));


            //Skip input file if a corresponding output file has already been written.
            if (outfile.isFile()) {
                System.out.println(
                        "An output file " + outfile.getName() + "\n"
                        + " already exists in the output directory for file " + qseqFileS[laneNum] + ".  Skipping.");
                continue;
            }

            int goodBarcodedReads = 0, allReads = 0, goodMatched = 0;
            File qseqFile = new File(qseqFileS[laneNum]);
            String[] np = qseqFile.getName().split("_");

            //Create a new object to hold barcoded tags.  The constructor can optionally process a group of fastq
            //files.  A minimum quality score for inclusion of a read can also be provided.
            ParseBarcodeRead thePBR;
            if (np.length == 3) {
                thePBR = new ParseBarcodeRead(keyFileS, enzyme, np[0], np[1]);
            } else if (np.length == 4) {
                thePBR = new ParseBarcodeRead(keyFileS, enzyme, np[0], np[2]);
            } else if (np.length == 5) {
                thePBR = new ParseBarcodeRead(keyFileS, enzyme, np[1], np[3]);
            } else {
                System.out.println("Error in parsing file name:");
                System.out.println("   The filename does not contain either 3 or 5 underscore-delimited values.");
                System.out.println("   Expect: flowcell_lane_qseq.txt OR code_flowcell_s_lane_qseq.txt");
                System.out.println("   Filename: " + qseqFileS[laneNum]);
                return;
            }
            System.out.println("Total barcodes found in lane:" + thePBR.getBarCodeCount());
            if (thePBR.getBarCodeCount() == 0) {
                System.out.println("No barcodes found.  Skipping this flowcell lane.");
                continue;
            }

            //Fill an array with taxon names.
            String[] taxaNames = new String[thePBR.getBarCodeCount()];
            for (int i = 0; i < taxaNames.length; i++) {
                taxaNames[i] = thePBR.getTheBarcodes(i).getTaxaName();
            }


            if (useTBTByte) {
                theTBT = new TagsByTaxaByte(taxaNames, theMasterTags);
            } else {
                theTBT = new TagsByTaxaBit(taxaNames, theMasterTags);
            }

            // Read the qseq file and assign reads to tags and taxa
            String temp = "";
            goodBarcodedReads = 0;
            allReads = 0;
            goodMatched = 0;
            try {
                BufferedReader br;
                //Read in qseq file as a gzipped text stream if its name ends in ".gz", otherwise read as text
                if (qseqFileS[laneNum].endsWith(".gz")) {
                    br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(qseqFileS[laneNum]))));
                } else {
                    br = new BufferedReader(new FileReader(qseqFileS[laneNum]), 65536);
                }
                String sl, qualS = "";
                while ((temp = br.readLine()) != null) {
                    String[] jj = temp.split("\\s");
                    allReads++;
                    if (allReads % 1000000 == 0) {
                        System.out.println("Total Reads:" + allReads + " goodReads:" + goodBarcodedReads + " goodMatched:" + goodMatched);
                    }
                    sl = jj[8];
                    qualS = jj[9];
                    ReadBarcodeResult rr = thePBR.parseReadIntoTagAndTaxa(sl, qualS, false, 0);
                    if (rr != null) {
                        goodBarcodedReads++;
                        int t = theTBT.getIndexOfTaxaName(rr.getTaxonName());
                        int h = theTBT.getTagIndex(rr.getRead());
                        if (h > -1) {
                            theTBT.addReadsToTagTaxon(h, t, 1);
                            goodMatched++;
                        }
                    }
                }
                br.close();
            } catch (Exception e) {
                System.out.println("Catch testBasicPipeline c=" + goodBarcodedReads + " e=" + e);
                System.out.println(temp);
                e.printStackTrace();
            }
            System.out.println("Timing process (writing TagsByTaxa file)...");
            long timePoint1 = System.currentTimeMillis();
            theTBT.writeDistFile(outfile, outFormat, minCount);
            System.out.println("...process (writing TagsByTaxa file) took " + (System.currentTimeMillis() - timePoint1) + " milliseconds.");
            System.out.println("Total number of reads in lane=" + allReads);
            System.out.println("Total number of good, barcoded reads=" + goodBarcodedReads);
            int filesDone = laneNum + 1;
            System.out.println("Finished reading " + filesDone + " of " + qseqFileS.length + " sequence files: " + qseqFileS[laneNum] + "\n");
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
