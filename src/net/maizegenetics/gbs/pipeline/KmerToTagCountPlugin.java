package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;

import javax.swing.ImageIcon;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;

import net.maizegenetics.util.MultiMemberGZIPInputStream;
import net.maizegenetics.gbs.util.BaseEncoder;
import net.maizegenetics.gbs.tagdist.TagCountMutable;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.ArgsEngine;

import org.apache.log4j.Logger;

/**
 * Derives a tagCount list for each fastq file in the input directory.
 *
 * Keeps only good reads having a barcode and a cut site and no N's in the
 * useful part of the sequence. Trims off the barcodes and truncates sequences
 * that (1) have a second cut site, or (2) read into the common adapter.
 *
 */
public class KmerToTagCountPlugin extends AbstractPlugin {

    static long timePoint1;
    private ArgsEngine engine = null;
    private Logger logger = Logger.getLogger(KmerToTagCountPlugin.class);
    String inputFile = null;
    int maxGoodReads = 200000000;
    int minCount = 1;
    int maxCount = 100;
    String outputFile = null;
    private static String nullS = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";

    public KmerToTagCountPlugin() {
        super(null, false);
    }

    public KmerToTagCountPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        logger.info(
                "\n\nUsage is as follows:\n"
                + " -i  Input file containing kmer files in fasta or gzipped fasta file.\n"
                + "     NOTE: Directory will be searched recursively and should\n"
                + "     be written WITHOUT a slash after its name.\n\n"
                + " -s  Max good reads per file. (Optional. Default is 200,000,000).\n"
                + " -cm  Minimum kmer depth (default is 1).\n"
                + " -cM  Maximum kmer depth (default is 100).\n"
                + " -o  Output .cnt files.\n\n");
    }

    public DataSet performFunction(DataSet input) {
        countTags(inputFile, outputFile, maxGoodReads, minCount, maxCount);
        return null;
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
            engine.add("-s", "--max-reads", true);
            engine.add("-cM", "--max-count", true);
            engine.add("-cm", "--min-count", true);
            engine.add("-o", "--output-file", true);
            engine.parse(args);
        }

        if (engine.getBoolean("-i")) {
            inputFile = engine.getString("-i");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify the input kmer files.");
        }


        if (engine.getBoolean("-s")) {
            maxGoodReads = Integer.parseInt(engine.getString("-s"));
        }

        if (engine.getBoolean("-cm")) {
            minCount = Integer.parseInt(engine.getString("-cm"));
        }
        if (engine.getBoolean("-cM")) {
            maxCount = Integer.parseInt(engine.getString("-cM"));
        }

        if (engine.getBoolean("-o")) {
            outputFile = engine.getString("-o");
        } else {
            outputFile = inputFile.substring(0, inputFile.indexOf(".")) + ".cnt";
        }

    }

    /**
     * Derives a tagCount list for each fastq file in the fastqDirectory.
     *
     *
     * @param myInput input file
     * @param myOutput outputfile.
     * @param maxGoodReads The maximum number of reads expected in a fastq file
     * @param minCount The minimum number of kmer depth
     * @param minCount The maximum number of kmer depth
     */
    public static void countTags(String myInput, String myOutput, int maxGoodReads, int minCount, int maxCount) {
        BufferedReader br;
        String[] countFileNames = null;
        int allReads = 0, goodBarcodedReads = 0;
        File myoutputFile = new File(myOutput);
        File myinputFile = new File(myInput);

        if (!myinputFile.isFile()) {
            System.out.println(
                    "Input file " + myInput + "\n"
                    + " does not exists");
            return;
        }
        if (myoutputFile.isFile()) {
            System.out.println(
                    "An output file " + myOutput + "\n"
                    + " already exists");
            return;
        }

        TagCountMutable theTC = null;
        System.out.println("Reading kmer file: " + myInput);


        try {
            //Read in qseq file as a gzipped text stream if its name ends in ".gz", otherwise read as text
            if (myInput.endsWith(".gz")) {
                br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(myinputFile))));
            } else {
                br = new BufferedReader(new FileReader(myinputFile), 65536);
            }
            String temp;

            try {
                theTC = new TagCountMutable(2, maxGoodReads);
            } catch (OutOfMemoryError e) {
                System.out.println(
                        "Your system doesn't have enough memory to store the number of sequences"
                        + "you specified.  Try using a smaller value for the minimum number of reads.");
            }
            allReads = 0;
            int kmercount;
            String sl;
            long[] binary_s1 = new long[2];
            while (((temp = br.readLine()) != null) && (allReads < maxGoodReads)) {
                allReads++;
                kmercount = Integer.parseInt(temp.replaceAll("[^\\d]", ""));
                sl = br.readLine().trim();
                sl = sl + nullS;
                sl = sl.substring(0, 64);
                if ((kmercount >= minCount) && (kmercount <= maxCount)) {
                    try {
                        //The quality score is every 4th line; the sequence is every 4th line starting from the 2nd.
                        binary_s1 = BaseEncoder.getLongArrayFromSeq(sl);
                        theTC.addReadCount(binary_s1, 31, kmercount);
                        if (allReads % 1000000 == 0) {
                            System.out.println("Total Reads:" + allReads + " Reads with barcode and cut site overhang:" + goodBarcodedReads);
                        }

                    } catch (NullPointerException e) {
                        System.out.println("Unable to correctly parse the sequence and "
                                + "quality score from fastq file.  Your fastq file may have been corrupted.");
                        System.exit(0);
                    }
                }
            }
            System.out.println("Total number of reads in lane=" + allReads);
            System.out.println("Timing process (sorting, collapsing, and writing TagCount to file).");
            timePoint1 = System.currentTimeMillis();
            theTC.collapseCounts();
            theTC.writeTagCountFile(myOutput, FilePacking.Bit, minCount);
            System.out.println("Process took " + (System.currentTimeMillis() - timePoint1) + " milliseconds.");
            br.close();
        } catch (Exception e) {
            System.out.println("Catch testBasicPipeline c=" + goodBarcodedReads + " e=" + e);
            e.printStackTrace();
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
