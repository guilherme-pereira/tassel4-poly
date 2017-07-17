/*
 * MergeTagsByTaxaFilesPlugin
 */
package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;

import java.io.File;

import java.util.TreeSet;

import javax.swing.ImageIcon;

import net.maizegenetics.gbs.tagdist.TagCountMutable;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxaBitFileMap;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByteFileMap;
import net.maizegenetics.gbs.tagdist.TagsByTaxaShortFileMap;
import net.maizegenetics.gbs.tagdist.TagsByTaxaUtils;
import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.gbs.util.ReadsByTaxa;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.DirectoryCrawler;

import org.apache.log4j.ConsoleAppender;
import org.apache.log4j.Logger;
import org.apache.log4j.SimpleLayout;

/**
 * Merges multiple TagsByTaxa files that are too large to fit in memory when
 * combined.  Currently, the output file stores only presence or absence of a
 * tag in a taxon, so the merged count  is the boolean OR of the individual counts.
 *
 * The program loops over files to determine their size, creates a RandomAccessFile object
 * large enough to hold the merged data, then loops over each file again to
 * fill the RandomAccessFile.
 * @author edbuckler
 */
public class MergeTagsByTaxaFilesPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(MergeTagsByTaxaFilesPlugin.class);
    static File inputDirectory = null, hapmapFile = null;
    static String[] infiles = null;
    static String outfilename = null;
    static ArgsEngine myArgsEngine = null;
    static boolean combineSynonymousTaxa = false;
    static int maxTags = 200000000;

    public MergeTagsByTaxaFilesPlugin() {
        super(null, false);
    }

    public MergeTagsByTaxaFilesPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    /**Determines the dimensions of the merged matrix, and creates a RandomAccessFile on disk large enough
    to hold the data.*/
    public static final void createMergeOutfile(String[] infiles, String outfile) {
        TagCountMutable theTCM = new TagCountMutable(2, maxTags);  // 50,000,000 works on a Windows7 machine with 4GB of RAM - JCG. 
        myLogger.info("Mutable Tag Count Created.  MaxTags:" + theTCM.getSize());
        myLogger.info("CurrentSize:" + theTCM.getCurrentSize());
        TreeSet<String> allTaxa = new TreeSet();
        int filesRead = 0;
        for (String inName : infiles) {
            filesRead++;
            myLogger.info("Scanning " + inName + " (file " + filesRead + " of " + infiles.length + ").");
            TagsByTaxa tbtFM = newTBT(inName);

            if (combineSynonymousTaxa) {
                tbtFM.truncateTaxonNames();
            }
            theTCM.addReadCounts(tbtFM, 1);
            for (String name : tbtFM.getTaxaNames()) {
                allTaxa.add(name);
            }
            myLogger.info("CurrentSize:" + theTCM.getCurrentSize());
            myLogger.info("Current Taxa:" + allTaxa.size());
            theTCM.collapseCounts();
            myLogger.info("CurrentSize:" + theTCM.getCurrentSize());
            myLogger.info("Size:" + theTCM.getSize());
        }
        theTCM.shrinkToCurrentRows();
        myLogger.info("Size:" + theTCM.getSize());
        String[] tn = allTaxa.toArray(new String[0]);
        TagsByTaxa tbtOut = newTBT(outfile, tn, theTCM);
    }

    /**Inserts tag count values into the RandomAccessFile created by @link{createMergeOutfile}.*/
    public static final void fillMergeOutfile(String[] infiles, String outfile) {
        TagsByTaxa tbtOut = newTBT(outfile);
        tbtOut.setMethodByRows(true);
        int count = 0;
        int filesRead = 0;
        for (String inName : infiles) {
            filesRead++;
            myLogger.info("Scanning " + inName + " (file " + filesRead + " of " + infiles.length + ").");
            TagsByTaxa tbtFM = newTBT(inName);
            if (combineSynonymousTaxa) {
                tbtFM.truncateTaxonNames();
            }
            int[] theTR = taxaRedirect(tbtFM.getTaxaNames(), tbtOut.getTaxaNames());
            for (int i = 0; i < tbtFM.getTagCount(); i++) {
                int toTag = tbtOut.getTagIndex(tbtFM.getTag(i));
                if (toTag < 0) {
                    continue;
                }
                for (int t = 0; t < tbtFM.getTaxaCount(); t++) {
                    if (theTR[t] < 0) {
                        continue;
                    }
                    int tagCount = tbtOut.getReadCountForTagTaxon(toTag, theTR[t]) + tbtFM.getReadCountForTagTaxon(i, t);
                    if (tagCount > 0) {
                        tbtOut.setReadCountForTagTaxon(toTag, theTR[t], tagCount);
                        count++;
                    }
                }
                if (count % 100000 == 0) {
                    System.out.printf("Tag:%d BitSet:%d %n", i, count);
                }
            }
        }
        tbtOut.getFileReadyForClosing();
    }

    public static int[] taxaRedirect(String[] fromNames, String[] toNames) {
        int[] theRedirect = new int[fromNames.length];
        for (int t = 0; t < fromNames.length; t++) {
            theRedirect[t] = -1;
            for (int i = 0; i < toNames.length; i++) {  //this is required if names are not sorted
                if (fromNames[t].equals(toNames[i])) {
                    theRedirect[t] = i;
                    break;
                }
            }
        }
        return theRedirect;
    }

    @Override
    public DataSet performFunction(DataSet input) {
    	createMergeOutfile(infiles, outfilename);
        fillMergeOutfile(infiles, outfilename);
        if (hapmapFile != null) {
            myLogger.info("Calling SNPs in good reads.");
            ReadsByTaxa rbt = new ReadsByTaxa();
            rbt.readTBTFile(new File(outfilename));
            Clusters cls = new Clusters(rbt);
            cls.networkFilter(); //Really powerful, which generates 85% single locus SNPs in 282 maize lines
            cls.alleleFrequencyFileter(rbt, (double) 0.20, (double) 0.30); //Only for linkage pop which has allele frequency peaks
            cls.heteozygoteFilter(rbt); //Seems not useful, need to be refined or removed
            cls.writeHapMap(rbt, hapmapFile.getAbsolutePath(), (float) 0.9);
        }
        return null;
    }

    public void setParameters(String[] args) {
        myLogger.addAppender(new ConsoleAppender(new SimpleLayout()));
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        if (myArgsEngine == null) {
            myArgsEngine = new ArgsEngine();
            myArgsEngine.add("-i", "--input-directory", true);
            myArgsEngine.add("-o", "--output_file", true);
            myArgsEngine.add("-s", "--max_tags", true);
            myArgsEngine.add("-x", "--combine-synonymous-taxa");
            myArgsEngine.add("-h", "--write-hapmap", true);
        }

        myArgsEngine.parse(args);

        if (myArgsEngine.getBoolean("-x")) {
            combineSynonymousTaxa = true;
        }
        if (myArgsEngine.getBoolean("-h")) {
            hapmapFile = new File(myArgsEngine.getString("-h"));
        }
        if (myArgsEngine.getBoolean("-s")) {
            maxTags = Integer.parseInt(myArgsEngine.getString("-s"));
        }
        if (myArgsEngine.getBoolean("-o")) {
            outfilename = myArgsEngine.getString("-o");
            File outfile = new File(outfilename);
            if (outfile.isDirectory()) {
                printUsage();
                throw new IllegalArgumentException("The output filename you provided is a directory, not a file.");
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
            } else {
                infiles = DirectoryCrawler.listFileNames(".*\\.tbt\\.bin|.*\\.tbt\\.byte|.*\\.tbt\\.shrt", inputDirectory.getAbsolutePath());
                myLogger.info("Merging the following .tbt.bin or .tbt.byte or .tbt.shrt files...");
                for (String filename : infiles) {
                    if (!TagsByTaxaUtils.format(filename).equals(TagsByTaxaUtils.format(outfilename))) {
                        myLogger.warn("Input file extension does not match output file extension.");
                    }
                    myLogger.info(filename);
                }
                myLogger.info("...to \"" + outfilename + "\".");
            }
        } else {
            printUsage();
            throw new IllegalArgumentException("You forgot to provide an input directory name.");
        }
    }

    private static void printUsage() {
        myLogger.info(
                "\n\n\nUsage is as follows:\n"
                + "-i  Input directory containing .tbt.bin or .tbt.byte files\n"
                + "-o  Output file name\n"
                + "-s  Maximum number of tags the TBT can hold while merging (default: " + maxTags + ")\n"
                + "-x  Merge tag counts of taxa with identical names (default: false)\n"
                + "-h  Call snps in output and write to HapMap file with the provided name\n\n\n",
                new IllegalArgumentException());
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "FilterErrorForBiparental";
    }

    @Override
    public String getToolTipText() {
        return "FilterErrorForBiparental";
    }

    private static TagsByTaxa newTBT(String filename, String[] taxonNames, TagCountMutable tcm) {
        if (TagsByTaxaUtils.format(filename).equals(TagsByTaxa.FilePacking.Bit)) {
            return new TagsByTaxaBitFileMap(filename, taxonNames, tcm);
        } else if (TagsByTaxaUtils.format(filename).equals(TagsByTaxa.FilePacking.Byte)) {
            return new TagsByTaxaByteFileMap(filename, taxonNames, tcm);
        } else if (TagsByTaxaUtils.format(filename).equals(TagsByTaxa.FilePacking.Short)) {
            return new TagsByTaxaShortFileMap(filename, taxonNames, tcm);
        } else {
            return null;
        }
    }

    private static TagsByTaxa newTBT(String filename) {
        if (TagsByTaxaUtils.format(filename).equals(TagsByTaxa.FilePacking.Bit)) {
            return new TagsByTaxaBitFileMap(filename);
        } else if (TagsByTaxaUtils.format(filename).equals(TagsByTaxa.FilePacking.Byte)) {
            return new TagsByTaxaByteFileMap(filename);
        } else if (TagsByTaxaUtils.format(filename).equals(TagsByTaxa.FilePacking.Short)) {
            return new TagsByTaxaShortFileMap(filename);
        } else {
            return null;
        }
    }
}
