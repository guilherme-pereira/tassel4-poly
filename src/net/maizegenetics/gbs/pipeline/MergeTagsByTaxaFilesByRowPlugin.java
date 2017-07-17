/*
 * MergeTagsByTaxaFilesByRowPlugin
 */
package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;

import java.util.ArrayList;

import javax.swing.ImageIcon;

import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.tagdist.TagCountMutable;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.tagdist.TagsByTaxaBitFileMap;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByteFileMap;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByteFileMapSeq;
import net.maizegenetics.gbs.tagdist.TagsByTaxaUtils;
import net.maizegenetics.util.ArgsEngine;
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
public class MergeTagsByTaxaFilesByRowPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(MergeTagsByTaxaFilesByRowPlugin.class);
    static File inputDirectory = null, hapmapFile = null;
    static String topmFileName = "/Volumes/nextgen/Zea/build20120110/topm/allZea_mappedonly_chr5-10_20120115.topm";
    static String topmFileName2 = "/Volumes/LaCie/zea20120110c510.topm";
    static String[] infiles = null;
    static String outfilename = null;
    static ArgsEngine myArgsEngine = null;
    static boolean combineSynonymousTaxa = false;
    static int maxTags = 200000000;

    public MergeTagsByTaxaFilesByRowPlugin() {
        super(null, false);
    }

    public MergeTagsByTaxaFilesByRowPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    /**Determines the dimensions of the merged matrix, and creates a RandomAccessFile on disk large enough
    to hold the data.*/
    public static final void createMergeOutfile(String[] infiles, String outfile) {
        TagsOnPhysicalMap theTOPM = new TagsOnPhysicalMap(topmFileName2, true);
        ArrayList<String> allTaxa = new ArrayList();
        int filesRead = 0;
        int tagLengthInLong = theTOPM.getTagSizeInLong();
        TagsByTaxaByteFileMapSeq[] inTBT = new TagsByTaxaByteFileMapSeq[infiles.length];
        for (int i = 0; i < infiles.length; i++) {
            String inName = infiles[i];
            filesRead++;
            myLogger.info("Scanning " + inName + " (file " + filesRead + " of " + infiles.length + ").");
            inTBT[i] = new TagsByTaxaByteFileMapSeq(inName);
            for (String name : inTBT[i].getTaxaNames()) {
                allTaxa.add(name);
            }
        }
        System.out.println("Total Taxa:" + allTaxa.size());
        System.out.println(allTaxa.toString());
        try {
            DataOutputStream fw = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfilename), 65536 * 1024));
            fw.writeInt(theTOPM.getTagCount());
            fw.writeInt(theTOPM.getTagSizeInLong());
            fw.writeInt(allTaxa.size());
            for (int t = 0; t < allTaxa.size(); t++) {
                fw.writeUTF(allTaxa.get(t));
            }
            for (int i = 0; i < theTOPM.getTagCount(); i++) {
                long[] currtag = theTOPM.getTag(i);
                for (int j = 0; j < tagLengthInLong; j++) {
                    fw.writeLong(currtag[j]);
                }
                fw.writeByte(theTOPM.getTagLength(i));
                for (int f = 0; f < infiles.length; f++) {
                    byte[] ob = inTBT[f].advanceToTagDist(currtag);
                    fw.write(ob);
                }
                if (i % 10000 == 0) {
                    System.out.println("TagNumber out:" + i);
                }
            }
            fw.close();
            for (int f = 0; f < infiles.length; f++) {
                inTBT[f].getFileReadyForClosing();
            }
        } catch (Exception e) {
            System.out.println("Catch in writeTextDistFile writing output file e=" + e);
            e.printStackTrace();
        }

//        TagsByTaxaByte combTBT=new TagsByTaxaByte(outfilename,FilePacking.Byte);
//        combTBT.writeDistFile(new File(outfilename+".txt"), TagsByTaxa.FilePacking.Text, -1);
//        for (int f=0; f<infiles.length; f++) {
//               TagsByTaxaByte ainTBT=new TagsByTaxaByte(infiles[f],FilePacking.Byte);
//               ainTBT.writeDistFile(new File(infiles[f]+".txt"), TagsByTaxa.FilePacking.Text, -1);
//           }

    }

    public void creatMergeTBTByTagCount(String infileDirS, String outfileS, String tagCountFileS) {
        File[] tbtFiles = new File(infileDirS).listFiles();
        String[] tbtFileSs = new String[tbtFiles.length];
        for (int i = 0; i < tbtFiles.length; i++) {
            tbtFileSs[i] = tbtFiles[i].getAbsolutePath();
        }
        TagCounts tc = new TagCounts(tagCountFileS, FilePacking.Byte);
        this.createMergeTBTByTagCount(tbtFileSs, outfileS, tc);
    }

    private void createMergeTBTByTagCount(String[] infiles, String outfileS, TagCounts tc) {
        ArrayList<String> allTaxa = new ArrayList();
        int filesRead = 0;
        int tagLengthInLong = tc.getTagSizeInLong();
        TagsByTaxaByteFileMapSeq[] inTBT = new TagsByTaxaByteFileMapSeq[infiles.length];
        for (int i = 0; i < infiles.length; i++) {
            String inName = infiles[i];
            filesRead++;
            //myLogger.info("Scanning "+inName+" (file "+filesRead+" of "+ infiles.length+").");
            inTBT[i] = new TagsByTaxaByteFileMapSeq(inName);
            for (String name : inTBT[i].getTaxaNames()) {
                allTaxa.add(name);
            }
        }
        System.out.println("Total Taxa:" + allTaxa.size());
        System.out.println(allTaxa.toString());
        try {
            DataOutputStream fw = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS), 65536 * 1024));
            fw.writeInt(tc.getTagCount());
            fw.writeInt(tc.getTagSizeInLong());
            fw.writeInt(allTaxa.size());
            for (int t = 0; t < allTaxa.size(); t++) {
                fw.writeUTF(allTaxa.get(t));
            }
            for (int i = 0; i < tc.getTagCount(); i++) {
                long[] currtag = tc.getTag(i);
                for (int j = 0; j < tagLengthInLong; j++) {
                    fw.writeLong(currtag[j]);
                }
                fw.writeByte(tc.getTagLength(i));
                for (int f = 0; f < infiles.length; f++) {
                    byte[] ob = inTBT[f].advanceToTagDist(currtag);
                    fw.write(ob);
                }
                if (i % 10000 == 0) {
                    System.out.println("TagNumber out:" + i);
                }
            }
            fw.close();
            for (int f = 0; f < infiles.length; f++) {
                inTBT[f].getFileReadyForClosing();
            }
        } catch (Exception e) {
            System.out.println("Catch in writeTextDistFile writing output file e=" + e);
            e.printStackTrace();
        }

//        TagsByTaxaByte combTBT=new TagsByTaxaByte(outfilename,FilePacking.Byte);
//        combTBT.writeDistFile(new File(outfilename+".txt"), TagsByTaxa.FilePacking.Text, -1);
//        for (int f=0; f<infiles.length; f++) {
//               TagsByTaxaByte ainTBT=new TagsByTaxaByte(infiles[f],FilePacking.Byte);
//               ainTBT.writeDistFile(new File(infiles[f]+".txt"), TagsByTaxa.FilePacking.Text, -1);
//           }

    }

    @Override
    public DataSet performFunction(DataSet input) {
        createMergeOutfile(infiles, outfilename);

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
                infiles = DirectoryCrawler.listFileNames(".*\\.tbt\\.bin|.*\\.tbt\\.byte", inputDirectory.getAbsolutePath());
                myLogger.info("Merging the following .tbt.bin files...");
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
                + "-i  Input directory containing .tbt.bin files\n"
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
        } else {
            return null;
        }
    }

    private static TagsByTaxa newTBT(String filename) {
        if (TagsByTaxaUtils.format(filename).equals(TagsByTaxa.FilePacking.Bit)) {
            return new TagsByTaxaBitFileMap(filename);
        } else if (TagsByTaxaUtils.format(filename).equals(TagsByTaxa.FilePacking.Byte)) {
            return new TagsByTaxaByteFileMap(filename);
        } else {
            return null;
        }
    }

    public static void main(String[] args) {
//         args = new String[] {
//            "-i", "/Users/edbuckler/SolexaAnal/GBS/build20120110/tbt/",
//            "-o", "/Users/edbuckler/SolexaAnal/GBS/build20120110/test.tbt.byte"
//        };
//         
        args = new String[]{
            "-i", "/Volumes/nextgen/Zea/build20120110/tbt/",
            "-o", "/Volumes/LaCie/zea20120110c510b.tbt.byte",};

        MergeTagsByTaxaFilesByRowPlugin testClass = new MergeTagsByTaxaFilesByRowPlugin();
        testClass.setParameters(args);
        testClass.performFunction(null);
    }
}