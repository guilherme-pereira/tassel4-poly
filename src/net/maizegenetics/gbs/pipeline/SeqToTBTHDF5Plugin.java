/*
 * SeqToTBTHDF5Plugin
 */
package net.maizegenetics.gbs.pipeline;

import cern.colt.list.IntArrayList;

import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.IOException;

import net.maizegenetics.gbs.homology.ParseBarcodeRead;
import net.maizegenetics.gbs.homology.ReadBarcodeResult;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.Tags;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByteHDF5TaxaGroups;
import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.util.MultiMemberGZIPInputStream;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;

import java.awt.Frame;

import java.util.HashMap;

import javax.swing.ImageIcon;

import org.apache.log4j.Logger;

/**
 * This pipeline converts a series of fastq or qseq files to a single TagsByTaxa HDF5 file.
 * It requires a list of existing tags (Tags object), which may come from a TagCounts file or TOPM file.
 *
 * @author ed and james
 */
public class SeqToTBTHDF5Plugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(SeqToTBTHDF5Plugin.class);
    private ArgsEngine myArgsEngine = null;
    private String[] myFastqFileS = null;
    private String myKeyFile = null;
    private String myEnzyme = null;
    private String myOutputTBTHDF5 = null;
    private String myOutputLogFile = null;
    private Tags myMasterTags = null;
    private static int maxGoodReads = 500000000; // maximum number of good barcoded reads expected in a fastq file
    IntArrayList[] taxaReads;
    int[] readsPerSample, mappedReadsPerSample;
    int goodBarcodedReads = 0, allReads = 0, goodMatched = 0;
    HashMap<String, Integer> taxaNameToIndices;

    public SeqToTBTHDF5Plugin() {
        super(null, false);
    }

    public SeqToTBTHDF5Plugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    @Override
    public DataSet performFunction(DataSet input) {
        File possibleFile = new File(myOutputTBTHDF5);
        if (possibleFile.exists() == false) {
            TagsByTaxaByteHDF5TaxaGroups theTBT = new TagsByTaxaByteHDF5TaxaGroups(myMasterTags, myOutputTBTHDF5);
        }
        matchTagsToTaxa(myFastqFileS, myKeyFile, myEnzyme, myMasterTags, myOutputTBTHDF5, myOutputLogFile);
        return null;
    }

    private void printUsage() {
        myLogger.info(
                "\nUsage is as follows:\n"
                + "-i  Input directory containing .fastq files\n"
                + "-k  Barcode key file\n"
                + "-e  Enzyme used to create the GBS library, if it differs from the one listed in the key file.\n"
                + "-o  Output HDF5 file\n"
                + "-s  Max good reads per lane. (Optional. Default is 500,000,000).\n"
                + "-L  Output log file \n"
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
            myArgsEngine.add("-o", "--output-HDF5", true);
            myArgsEngine.add("-L", "--outputlogfile", true);
            myArgsEngine.add("-s", "--max-reads", true);
            myArgsEngine.add("-t", "--tag-count", true);
            myArgsEngine.add("-m", "--physical-map", true);
        }
        myArgsEngine.parse(args);

        String tempDirectory = myArgsEngine.getString("-i");

        if (myArgsEngine.getBoolean("-s")) {
            maxGoodReads = Integer.parseInt(myArgsEngine.getString("-s"));
        }

        if (tempDirectory != null) {
            File fastqDirectory = new File(tempDirectory);
            if (!fastqDirectory.isDirectory()) {
                printUsage();
                throw new IllegalArgumentException("setParameters: The input name you supplied is not a directory: " + tempDirectory);
            }
            myFastqFileS = DirectoryCrawler.listFileNames("(?i).*\\.fq$|.*\\.fq\\.gz$|.*\\.fastq$|.*_fastq\\.txt$|.*_fastq\\.gz$|.*_fastq\\.txt\\.gz$|.*_sequence\\.txt$|.*_sequence\\.txt\\.gz$|.*_qseq\\.txt$|.*_qseq\\.txt\\.gz$", fastqDirectory.getAbsolutePath());
            //    (?i) denotes case insensitive;                 \\. denotes escape . so it doesn't mean 'any char'
            // NOTE: If you add addtional file naming conventions here, you must also add them to the "outfile = new File(outFileS.replaceAll()" list below (near the bottom)
            if (myFastqFileS.length == 0 || myFastqFileS == null) {
                printUsage();
                throw new IllegalArgumentException(
                        "Couldn't find any files that end with \".fq\", \".fq.gz\", \".fastq\", \"_fastq.txt\", \"_fastq.gz\", \"_fastq.txt.gz\", \"_sequence.txt\", or \"_sequence.txt.gz\" in the supplied directory: "
                        + tempDirectory);
            } else {
                myLogger.info("FastqToTBTPlugin: setParameters: Using the following fastq files:");
                for (String filename : myFastqFileS) {
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
            myLogger.warn("No enzyme specified.  Using enzyme listed in key file.");
        }
        if (myArgsEngine.getBoolean("-o")) {
            myOutputTBTHDF5 = myArgsEngine.getString("-o");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify an output directory (option -o).");
        }
        if (myArgsEngine.getBoolean("-L")) {
            myOutputLogFile = myArgsEngine.getString("-L");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify a log file (option -L).");
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
     * Uses an existing Tags object to create one TagsByTaxa file for each fastq file in the input directory.
     *
     * Output TBT files written to the outputDir, using fastq file names with extension changed to .tbt.bin (or .tbt.txt)
     *
     * @param fastqFileS      Array of fastq file names (Illumina-created files with raw read sequence, quality score, machine name, etc.)
     * @param keyFileS       A key file (list of taxa by barcode, lane & flow cell, including plate maps)
     * @param enzyme         The enzyme used to make the library (currently ApeKI or PstI)
     * @param theMasterTags  A Tags object: list of tags to be included in the final TBT
     * @param outputDir      String containing the path of the output directory to contain tags-by-taxa files
     * @param minCount       The minimum number of times a tag must show up in a fastq file before it is included in the corresponding TBT file
     */
    public void matchTagsToTaxa(String[] fastqFileS, String keyFileS, String enzyme, Tags theMasterTags, String outputTBT, String outputLog) {
        for (int laneNum = 0; laneNum < fastqFileS.length; laneNum++) {
            //            IntArrayList[] taxaReads;
            //            int[] readsPerSample, mappedReadsPerSample;
            //            int goodBarcodedReads = 0, allReads = 0, goodMatched = 0;

            goodBarcodedReads = allReads = goodMatched = 0;

            boolean isFastQ = true;
            if (fastqFileS[laneNum].substring(fastqFileS[laneNum].lastIndexOf(File.separator)).contains("qseq")) {
                isFastQ = false;
            }

            System.out.println("\nWorking on fastq file: " + fastqFileS[laneNum]);

            File fastqFile = new File(fastqFileS[laneNum]);

            ParseBarcodeRead thePBR = initParseBarcodeRead(keyFileS, enzyme, fastqFile);
            if (thePBR == null) {
                continue;  //if no barcodes found move on
            }
            //Fill an array with taxon names.
            String[] taxaNames = new String[thePBR.getBarCodeCount()];
            taxaReads = new IntArrayList[thePBR.getBarCodeCount()];
            mappedReadsPerSample = new int[thePBR.getBarCodeCount()];
            readsPerSample = new int[thePBR.getBarCodeCount()];
            taxaNameToIndices = new HashMap<String, Integer>();

            for (int i = 0; i < taxaNames.length; i++) {
                taxaNames[i] = thePBR.getTheBarcodes(i).getTaxaName();
                taxaReads[i] = new IntArrayList(200000000 / taxaNames.length);
                taxaNameToIndices.put(taxaNames[i], i);
            }
            int currLine = 0;
            goodBarcodedReads = 0;
            allReads = 0;
            goodMatched = 0;
            try {
                BufferedReader br = getBufferedReader(fastqFileS[laneNum]);
                String[] seqAndQual;
                while (((seqAndQual = getNextSeq(br, isFastQ)) != null) && (goodBarcodedReads < maxGoodReads)) {
                    allReads++;
                    currLine++;
                    //After quality score is read, decode barcode using the current sequence & quality  score
                    ReadBarcodeResult rr = thePBR.parseReadIntoTagAndTaxa(seqAndQual[0], seqAndQual[1], true, 0);
                    if (rr != null) {  //If read is barcoded, note this fact and find tag in the TOPM file
                        goodBarcodedReads++;
                        int t = taxaNameToIndices.get(rr.getTaxonName());
                        int h = theMasterTags.getTagIndex(rr.getRead());
                        readsPerSample[t]++;
                        if (h > -1) {
                            goodMatched++;
                            mappedReadsPerSample[t]++;
                            taxaReads[t].add(h);
                        }
                    }
                    if (allReads % 1000000 == 0) {
                        System.out.println("Total Reads:" + allReads + " goodReads:" + goodBarcodedReads + " goodMatched:" + goodMatched);
                    }
                }
                br.close();
            } catch (Exception e) {
                System.out.println("Catch testBasicPipeline c=" + goodBarcodedReads + " e=" + e);
                e.printStackTrace();
            }
            System.out.println("Timing process (writing TagsByTaxa file)...");
            long timePoint1 = System.currentTimeMillis();
            writeReport(fastqFileS[laneNum], outputLog);
            writeTBT(outputTBT);
            System.out.println("...process (writing TagsByTaxa file) took " + (System.currentTimeMillis() - timePoint1) + " milliseconds.");
            System.out.println("Total number of reads in lane=" + allReads);
            System.out.println("Total number of good, barcoded reads=" + goodBarcodedReads);
            int filesDone = laneNum + 1;
            System.out.println("Finished reading " + filesDone + " of " + fastqFileS.length + " sequence files: " + fastqFileS[laneNum] + "\n");
//            //Write report to a file named after the TBT file, with the extension ".log"
//            float goodPct = ((float) goodBarcodedReads) / ((float) allReads);
//            float goodMappedPct = ((float) goodMatched) / ((float) allReads);
//            try {
//                DataOutputStream report = new DataOutputStream(new FileOutputStream(outputLog, true));
//                report.writeBytes(
//                        "File: " + fastqFileS[laneNum] + '\n'
//                        + "Total reads: " + allReads + '\n'
//                        + "Accepted reads (with barcode and cut site): " + goodBarcodedReads + "(" + goodPct + " of total)" + '\n'
//                        + "Accepted reads found in TOPM: " + goodMatched + "(" + goodMappedPct + " of total)" + '\n'
//                        + "name" + '\t' + "read count" + '\t' + "fraction of total" + '\t' + "mapped read count" + '\t' + "fraction mapped of total" + '\n');
//               // TagsByTaxaByteHDF5TaxaGroups theTBT = new TagsByTaxaByteHDF5TaxaGroups(outputTBT);
//                for (String name : taxaNameToIndices.keySet()) {
//                    int t=taxaNameToIndices.get(name);
//                    int count = readsPerSample[t];
//                    int mappedCount = mappedReadsPerSample[t];
//                    float pct = ((float) count) / ((float) allReads);
//                    float mappedPct = ((float) mappedCount) / ((float) count);
//                    byte[] tagByteDist=getTagsDistribution(theMasterTags.getTagCount(),taxaReads[t]);
//                    if(TagsByTaxaByteHDF5TaxaGroups.addTaxon(outputTBT, name, tagByteDist)) {
//                        System.out.printf("Taxon %s written to %s %n", name, outputTBT);
//                    } else {
//                        System.err.printf("Error Writing Taxon %s written to %s %n", name, outputTBT);
//                    }
//                    report.writeBytes(name + '\t' + count + '\t' + pct + '\t' + mappedCount + '\t' + mappedPct + '\n');
//                }
//                report.close();
//             //   theTBT.getFileReadyForClosing();
//            } catch (Exception e) {
//                myLogger.warn("Caught exception while writing report file for lane " + laneNum + ": " + e);
//            }
        }
    }

    private synchronized boolean writeReport(String seqFile, String outputLog) {
        //Write report to a file named after the TBT file, with the extension ".log"
        float goodPct = ((float) goodBarcodedReads) / ((float) allReads);
        float goodMappedPct = ((float) goodMatched) / ((float) allReads);
        try {
            DataOutputStream report = new DataOutputStream(new FileOutputStream(outputLog, true));
            report.writeBytes(
                    "File: " + seqFile + '\n'
                    + "Total reads: " + allReads + '\n'
                    + "Accepted reads (with barcode and cut site): " + goodBarcodedReads + "(" + goodPct + " of total)" + '\n'
                    + "Accepted reads found in TOPM: " + goodMatched + "(" + goodMappedPct + " of total)" + '\n'
                    + "name" + '\t' + "read count" + '\t' + "fraction of total" + '\t' + "mapped read count" + '\t' + "fraction mapped of total" + '\n');
            // TagsByTaxaByteHDF5TaxaGroups theTBT = new TagsByTaxaByteHDF5TaxaGroups(outputTBT);
            for (String name : taxaNameToIndices.keySet()) {
                int t = taxaNameToIndices.get(name);
                int count = readsPerSample[t];
                int mappedCount = mappedReadsPerSample[t];
                float pct = ((float) count) / ((float) allReads);
                float mappedPct = ((float) mappedCount) / ((float) count);
                report.writeBytes(name + '\t' + count + '\t' + pct + '\t' + mappedCount + '\t' + mappedPct + '\n');
            }
            report.close();
            //   theTBT.getFileReadyForClosing();
        } catch (Exception e) {
            myLogger.warn("Caught exception while writing report file for file " + seqFile + ": " + e);
            return false;
        }
        return true;
    }

    private synchronized boolean writeTBT(String outputTBT) {
        TagsByTaxaByteHDF5TaxaGroups theTBT = new TagsByTaxaByteHDF5TaxaGroups(outputTBT);
        for (String name : taxaNameToIndices.keySet()) {
            int t = taxaNameToIndices.get(name);
            byte[] tagByteDist = getTagsDistribution(theTBT.getTagCount(), taxaReads[t]);
            theTBT.addTaxon(name, tagByteDist);
            System.out.printf("Taxon %s written to %s %n", name, outputTBT);
        }
        theTBT.getFileReadyForClosing();
        theTBT = null;
        try {
            System.out.println("Sleeping thread" + Thread.currentThread().getName());
            Thread.sleep(5000);
        } catch (InterruptedException ie) {
            System.err.println("Sleeping failed");
            return false;
        }
        return true;
    }

    private static String[] getNextSeq(BufferedReader br, boolean isFastQ) {
        String[] s = new String[2];
        try {
            if (isFastQ) {
                String temp = br.readLine();  //read sequence identifier
                if (temp == null) {
                    return null;
                }
                while (temp.charAt(0) != '@') {
                    temp = br.readLine();
                    if (temp == null) {
                        return null;
                    }
                }
                s[0] = br.readLine(); //read sequence
                br.readLine();  //read quality identifier
                s[1] = br.readLine();  //read quality
            } else {
                String temp = br.readLine();  //read sequence identifier
                if (temp == null) {
                    return null;
                }
                String[] jj = temp.split("\\s");
                s[0] = jj[8];
                s[1] = jj[9];
            }
        } catch (IOException e) {
            System.out.println("File closing");
            return null;
        }
        return s;
    }

    private static BufferedReader getBufferedReader(String fileName) {
        BufferedReader br = null;
        try {
            if (fileName.endsWith(".gz")) {
                br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(fileName))));
            } else {
                br = new BufferedReader(new FileReader(fileName), 65536);
            }
        } catch (IOException e) {
            System.out.println("Failed to open file:" + fileName);

        }
        return br;
    }

    private static byte[] getTagsDistribution(int tagNumber, IntArrayList tagDist) {
        byte[] result = new byte[tagNumber];
        for (int i : tagDist.elements()) {
            if (result[i] < Byte.MAX_VALUE) {
                result[i]++;
            }
        }
        return result;
    }

    private static ParseBarcodeRead initParseBarcodeRead(String keyFileS, String enzyme, File fastqFile) {
        ParseBarcodeRead thePBR = null;
        String[] np = fastqFile.getName().split("_");

        //Create a new object to hold barcoded tags.  The constructor can optionally process a group of fastq
        //files.  A minimum quality score for inclusion of a read can also be provided.

        if (np.length == 3) {
            thePBR = new ParseBarcodeRead(keyFileS, enzyme, np[0], np[1]);
        } else if (np.length == 5) {
            thePBR = new ParseBarcodeRead(keyFileS, enzyme, np[1], np[3]);
        } else if (np.length == 4) {
            thePBR = new ParseBarcodeRead(keyFileS, enzyme, np[0], np[2]);
        } else if (np.length == 6) {
            thePBR = new ParseBarcodeRead(keyFileS, enzyme, np[1], np[3]);
        } else {
            System.out.println("Error in parsing file name:");
            System.out.println("   The filename does not contain either 3 or 5 underscore-delimited values.");
            System.out.println("   Expect: flowcell_lane_fastq.txt OR code_flowcell_s_lane_fastq.txt");
            System.out.println("   Filename: " + fastqFile.getName());
            return null;
        }
        System.out.println("Total barcodes found in lane:" + thePBR.getBarCodeCount());
        if (thePBR.getBarCodeCount() == 0) {
            System.out.println("No barcodes found.  Skipping this flowcell.");
            return null;
        }
        return thePBR;
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
