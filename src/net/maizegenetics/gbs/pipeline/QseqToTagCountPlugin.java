package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;

import net.maizegenetics.util.MultiMemberGZIPInputStream;
import javax.swing.ImageIcon;
import net.maizegenetics.gbs.homology.ParseBarcodeRead;
import net.maizegenetics.gbs.homology.ReadBarcodeResult;
import net.maizegenetics.gbs.tagdist.TagCountMutable;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.DirectoryCrawler;

import org.apache.log4j.Logger;

/** 
 * Derives a tagCount list for each qseq file in the qseqDirectory.
 *
 * Keeps only good reads having a barcode and a cut site and no N's in the
 * useful part of the sequence.  Trims off the barcodes and truncates sequences
 * that (1) have a second cut site, or (2) read into the common adapter.
 * 
 */
public class QseqToTagCountPlugin extends AbstractPlugin {

    static long timePoint1;
    private ArgsEngine engine = null;
    private Logger logger = Logger.getLogger(QseqToTagCountPlugin.class);
    String directoryName = null;
    String keyfile = null;
    String enzyme = null;
    int maxGoodReads = 300000000;
    int minCount = 1;
    String outputDir = null;

    public QseqToTagCountPlugin() {
        super(null, false);
    }

    public QseqToTagCountPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        logger.info(
                "\n\nUsage is as follows:\n"
                + " -i  Input directory containing qseq files as text or gzipped text.\n"
                + "     NOTE: Directory will be searched recursively and should\n"
                + "     be written WITHOUT a slash after its name.\n\n"
                + " -k  Key file listing barcodes for each sample\n"
                + " -e  Enzyme used to create the GBS library, if it differs from the one listed in the key file.\n"
                + " -s  Max good reads per lane. (Optional. Default is 300,000,000).\n"
                + " -c  Minimum tag count (default is 1).\n"
                + " -o  Output directory to contain .cnt files (one per qseq file, defaults to input directory).\n\n");
    }

    public DataSet performFunction(DataSet input) {
        File qseqDirectory = new File(directoryName);
        if (!qseqDirectory.isDirectory()) {
            throw new IllegalStateException("The input name you supplied is not a directory.");
        }
        countTags(keyfile, enzyme, directoryName, outputDir, maxGoodReads, minCount);  //TODO change to perform function
        return null;
    }

    @Override
    public void setParameters(String[] args) {
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

//        try{
        if (engine == null) {
            engine = new ArgsEngine();
            engine.add("-i", "--input-directory", true);
            engine.add("-k", "--key-file", true);
            engine.add("-e", "--enzyme", true);
            engine.add("-s", "--max-reads", true);
            engine.add("-c", "--min-count", true);
            engine.add("-o", "--output-file", true);
            engine.parse(args);
        }

        if (engine.getBoolean("-i")) {
            directoryName = engine.getString("-i");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify the location of your .qseq files.");
        }

        if (engine.getBoolean("-k")) {
            keyfile = engine.getString("-k");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify a barcode key file.");
        }

        if (engine.getBoolean("-e")) {
            enzyme = engine.getString("-e");
        } else {
            System.out.println("No enzyme specified.  Using enzyme listed in key file.");
//                    printUsage(); throw new IllegalArgumentException("Please specify the enzyme used to create the GBS library.");
        }

        if (engine.getBoolean("-s")) {
            maxGoodReads = Integer.parseInt(engine.getString("-s"));
        }

        if (engine.getBoolean("-c")) {
            minCount = Integer.parseInt(engine.getString("-c"));
        }

        if (engine.getBoolean("-o")) {
            outputDir = engine.getString("-o");
        } else {
            outputDir = directoryName;
        }
//        }catch (Exception e){
//            System.out.println("Caught exception while setting parameters of "+this.getClass()+": "+e);
//        }
    }

    /**
     * Derives a tagCount list for each qseq file in the qseqDirectory.
     *
     * @param keyFileS  A key file (a sample key by barcode, with a plate map included).
     * @param enzyme  The enzyme used to create the library (currently ApeKI or PstI).
     * @param qseqDirectory  Directory containing the qseq files (will be recursively searched).
     * @param outputDir  Directory to which the tagCounts files (one per qseq file) will be written.
     * @param maxGoodReads  The maximum number of barcoded reads expected in a qseq file
     * @param minCount  The minimum number of occurrences of a tag in a qseq file for it to be included in the output tagCounts file
     */
    public static void countTags(String keyFileS, String enzyme, String qseqDirectory, String outputDir, int maxGoodReads, int minCount) {
        BufferedReader br;
        String[] countFileNames = null;

        File inputDirectory = new File(qseqDirectory);
        File[] qseqFiles = DirectoryCrawler.listFiles(".*_qseq\\.txt$|.*_qseq\\.txt\\.gz$", inputDirectory.getAbsolutePath());

        if (qseqFiles.length == 0 || qseqFiles == null) {
            System.out.println("Couldn't find any files that end with \"_qseq.txt\" or \"_qseq.txt.gz\" in the supplied directory.");
        } else {
            System.out.println("Using the following .qseq files:");
            countFileNames = new String[qseqFiles.length];
            for (int i = 0; i < qseqFiles.length; i++) {
                countFileNames[i] = qseqFiles[i].getName().replaceAll("_qseq\\.txt$|_qseq\\.txt\\.gz$", ".cnt");
                System.out.println(qseqFiles[i].getAbsolutePath());
            }
        }

        int allReads = 0, goodBarcodedReads = 0;
        for (int laneNum = 0; laneNum < qseqFiles.length; laneNum++) {
            File outputFile = new File(outputDir + File.separator + countFileNames[laneNum]);
            if (outputFile.isFile()) {
                System.out.println(
                        "An output file " + countFileNames[laneNum] + "\n"
                        + " already exists in the output directory for file " + qseqFiles[laneNum] + ".  Skipping.");
                continue;
            }

            TagCountMutable theTC = null;
            System.out.println("Reading qseq file: " + qseqFiles[laneNum]);
            String[] filenameField = qseqFiles[laneNum].getName().split("_");
            ParseBarcodeRead thePBR;  // this reads the key file and store the expected barcodes for this lane
            if (filenameField.length == 3) {
                thePBR = new ParseBarcodeRead(keyFileS, enzyme, filenameField[0], filenameField[1]);
            } else if (filenameField.length == 4) {
                thePBR = new ParseBarcodeRead(keyFileS, enzyme, filenameField[0], filenameField[2]);
            } else if (filenameField.length == 5) {
                thePBR = new ParseBarcodeRead(keyFileS, enzyme, filenameField[1], filenameField[3]);
            } else {
                System.out.println("Error in parsing file name:");
                System.out.println("   The filename does not contain either 3 or 5 underscore-delimited values.");
                System.out.println("   Expect: flowcell_lane_qseq.txt OR code_flowcell_s_lane_qseq.txt");
                System.out.println("   Filename: " + qseqFiles[laneNum]);
                return;
            }
            System.out.println("Total barcodes found in lane:" + thePBR.getBarCodeCount());
            if (thePBR.getBarCodeCount() == 0) {
                System.out.println("No barcodes found.  Skipping this flowcell lane.");
                continue;
            }
            String[] taxaNames = new String[thePBR.getBarCodeCount()];
            for (int i = 0; i < taxaNames.length; i++) {
                taxaNames[i] = thePBR.getTheBarcodes(i).getTaxaName();
            }

            try {
                //Read in qseq file as a gzipped text stream if its name ends in ".gz", otherwise read as text
                if (qseqFiles[laneNum].getName().endsWith(".gz")) {
                    br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(qseqFiles[laneNum]))));
                } else {
                    br = new BufferedReader(new FileReader(qseqFiles[laneNum]), 65536);
                }
                String sequence = "", qualityScore = "";
                String temp;

                try {
                    theTC = new TagCountMutable(2, maxGoodReads);
                } catch (OutOfMemoryError e) {
                    System.out.println(
                            "Your system doesn't have enough memory to store the number of sequences"
                            + "you specified.  Try using a smaller value for the maximum number of good reads (-s option).");
                }
                allReads = 0;
                goodBarcodedReads = 0;
                while (((temp = br.readLine()) != null) && goodBarcodedReads < maxGoodReads) {
                    String[] jj = temp.split("\t");
                    allReads++;
                    try {
                        sequence = jj[8];
                        qualityScore = jj[9];
                    } catch (NullPointerException e) {
                        System.out.println("Read a line that lacks a sequence and "
                                + "quality score in fields 9 and 10.  Your file may have been corrupted.");
                        System.exit(0);
                    }
                    ReadBarcodeResult rr = thePBR.parseReadIntoTagAndTaxa(sequence, qualityScore, false, 0);
                    if (rr != null) {
                        goodBarcodedReads++;
                        theTC.addReadCount(rr.getRead(), rr.getLength(), 1);
                    }
                    if (allReads % 1000000 == 0) {
                        System.out.println("Total Reads:" + allReads + " Reads with barcode and cut site overhang:" + goodBarcodedReads);
                    }
                }
                System.out.println("Total number of reads in lane=" + allReads);
                System.out.println("Total number of good barcoded reads=" + goodBarcodedReads);
                System.out.println("Timing process (sorting, collapsing, and writing TagCount to file).");
                timePoint1 = System.currentTimeMillis();
                theTC.collapseCounts();
                theTC.writeTagCountFile(outputDir + File.separator + countFileNames[laneNum], FilePacking.Bit, minCount);
                System.out.println("Process took " + (System.currentTimeMillis() - timePoint1) + " milliseconds.");
                br.close();

            } catch (Exception e) {
                System.out.println("Catch testBasicPipeline c=" + goodBarcodedReads + " e=" + e);
                e.printStackTrace();
            }
            System.out.println("Finished reading " + (laneNum + 1) + " of " + qseqFiles.length + " sequence files.");
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
